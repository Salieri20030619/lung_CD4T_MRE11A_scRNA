import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1')
from pathlib import Path
import sys,json,hashlib,subprocess,importlib.util,concurrent.futures
import numpy as np,pandas as pd
from scipy import sparse
from ctxcore.rnkdb import FeatherRankingDatabase
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/death_priority_20261018_v2';N=O/'joint_network';B=P/'11_reports/three_cohort_four_axis_20261017_v16';R=P/'00_raw/SCENIC_resources_v10/Mouse';N.mkdir(exist_ok=True)
for k in ['ROOT','PREV1','PREV2','PREV7']:os.environ[k]=str(P)
os.environ['OUT']=str(N);os.environ['CACHE']=str(R.parent);os.environ['WORKERS']='8'
spec=importlib.util.spec_from_file_location('original_pipeline',str(P/'02_analysis/scenic1_20260929_103513_xRoEnr/pipeline.py'));old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
dbs=[R/'mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather',R/'mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'];anno=R/'motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl';tflist=R/'allTFs_mm.txt'
alias={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
def canon(v):return alias.get(str(v).upper(),str(v).upper())
universe=set(FeatherRankingDatabase(str(dbs[0]),name='a').genes)&set(FeatherRankingDatabase(str(dbs[1]),name='b').genes);mapping={}
for g in sorted(universe):mapping.setdefault(canon(g),[]).append(g)
xs=[];ms=[];gs=[]
for ds in ['GSE314943','GSE183682']:
 stem=ds+'__Core';xs.append(sparse.load_npz(B/'inputs'/(stem+'.cellcounts.npz')));g=pd.read_csv(B/'inputs'/(stem+'.genes.txt'),header=None)[0].to_numpy();gs.append(g);m=pd.read_csv(B/'inputs'/(stem+'.cells.tsv.gz'),sep='\t');m['series']=ds;m['dataset']=ds;m['cell']=ds+'|'+m.cell;m['core']=True;m['bg']=m.background_fraction;ms.append(m)
common=sorted(set(gs[0])&set(gs[1]));x=sparse.vstack([v[:,pd.Index(g).get_indexer(common)] for v,g in zip(xs,gs)],format='csr');m=pd.concat(ms,ignore_index=True);genes=np.array([mapping[g][0] if g in mapping and len(mapping[g])==1 else 'UNMAPPED_'+g for g in common])
assert m['sample'].nunique()==10 and m.cell.is_unique and len(m)==3724
pd.DataFrame({'canonical':common,'motif_symbol':genes,'mapping_status':['UNIQUE' if not g.startswith('UNMAPPED_') else 'UNMAPPED_OR_AMBIGUOUS' for g in genes]}).to_csv(N/'mapping.tsv',sep='\t',index=False)
protocol={'stage':'joint_10_sample_network_reconstruction','cohort':'GSE314943 6 + GSE183682 4','arms':['All','BG1'],'seeds':[17,29,43],'sample_cap':250,'original_method':'SCENIC original pipeline prepare and CLI parameters reused','workers_per_job':8,'concurrent_jobs':6,'genes':'common measured, unique canonical-to-motif mapping, detection>=max(10,1%cells); full library logCP10k; no outcome-selected genes','stability':'motif edges retained in>=2of3 seeds; resampling stability not independent validation','downstream':'remove candidate/program/self genes before stable-regulon AUCell; score all retained cells; sample-level joint effects and associations','old_network_projection':'preserved separate comparator; this stage completes actual joint network inference'}
(N/'protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
resources=[]
for p in dbs+[anno,tflist]:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for ch in iter(lambda:f.read(8*1024*1024),b''):h.update(ch)
 resources.append(dict(path=str(p),sha256=h.hexdigest(),bytes=p.stat().st_size))
(N/'resource_manifest.json').write_text(json.dumps(resources,indent=2),encoding='utf-8')
arms=[]
for arm in ['All','BG1']:
 d=N/'Mouse'/arm
 if not (d/'training_audit.tsv').exists():d=old.prepare('Mouse',arm,x,genes,m,dbs,tflist)
 assert pd.read_csv(d/'cells.tsv',sep='\t')['sample'].nunique()==10;arms.append(d)
print('JOINT_NETWORK_INPUTS_READY',flush=True)
def fit_one(d,seed):
 sd=d/f'seed{seed}';py=sys.executable;cli=Path(py).parent/'pyscenic'
 commands=[('grn',[py,'-m','pyscenic.cli.arboreto_with_multiprocessing',str(sd/'train.loom'),str(sd/'tfs.txt'),'--method','grnboost2','--num_workers','8','--seed',str(seed),'-o',str(sd/'adj.tsv')]),('ctx',[str(cli),'ctx',str(sd/'adj.tsv'),*[str(p) for p in dbs],'--annotations_fname',str(anno),'--expression_mtx_fname',str(sd/'train.loom'),'--mode','custom_multiprocessing','--num_workers','8','--output',str(sd/'regulons.csv')])]
 for stage,cmd in commands:
  print('START',d.name,seed,stage,flush=True)
  with (sd/(stage+'.log')).open('wb') as log:p=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
  (sd/(stage+'_exit.json')).write_text(json.dumps({'returncode':p.returncode,'command':cmd}),encoding='utf-8')
  if p.returncode:raise RuntimeError(f'{d.name} seed{seed} {stage} failed exit{p.returncode}')
  print('DONE',d.name,seed,stage,flush=True)
 return {'arm':d.name,'seed':seed,'status':'success'}
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
  futures=[ex.submit(fit_one,d,seed) for d in arms for seed in [17,29,43]]
  rows=[f.result() for f in concurrent.futures.as_completed(futures)]
 (N/'network_inference_complete.json').write_text(json.dumps({'status':'success','runs':rows},indent=2),encoding='utf-8');print('ALL_JOINT_NETWORKS_COMPLETE',flush=True)
