import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1',TMPDIR='/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/tmp/death_priority_v3')
from pathlib import Path
import sys,json,subprocess,concurrent.futures,time
import numpy as np,pandas as pd,loompy
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
V=P/'11_reports/death_priority_20261018_v2';O=P/'11_reports/death_priority_20261018_v3';N=O/'network_loo';R=P/'00_raw/SCENIC_resources_v10/Mouse'
for s in ['code','tables','logs','external','figures','network_loo']: (O/s).mkdir(parents=True,exist_ok=True)
dbs=[R/'mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather',R/'mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather']
anno=R/'motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl'
jobs=[]
for arm in ['All','BG1']:
 d=V/'joint_network/Mouse'/arm/'seed17'
 m=pd.read_csv(d/'train_cells.tsv',sep='\t')
 with loompy.connect(str(d/'train.loom'),'r') as ds: x=ds[:,:];g=np.asarray(ds.ra.Gene).astype(str);cells=np.asarray(ds.ca.CellID).astype(str)
 assert np.array_equal(cells,m.cell.to_numpy()) and m['sample'].nunique()==10
 for sam in sorted(m['sample'].unique()):
  out=N/arm/sam;out.mkdir(parents=True,exist_ok=True);sel=m['sample'].ne(sam).to_numpy()
  xx=x[:,sel];keep=xx.var(axis=1)>1e-8;gg=g[keep];tfs=[t for t in (d/'tfs.txt').read_text().splitlines() if t in set(gg)]
  if not (out/'train.loom').exists():loompy.create(str(out/'train.loom'),xx[keep],{'Gene':gg},{'CellID':cells[sel]})
  (out/'tfs.txt').write_text('\n'.join(tfs)+'\n',encoding='utf-8')
  m[sel].to_csv(out/'train_cells.tsv',sep='\t',index=False)
  jobs.append((arm,sam,out))
(O/'protocol_v3.json').write_text(json.dumps({'network_loo':{'arms':['All','BG1'],'seed':17,'folds':20,'comparison':'paired to original seed17; retain original cell draw and gene universe, remove one biological donor and zero-variance genes per fold','workers_per_job':10,'concurrent_jobs':8,'interpretation':'donor sensitivity; not independent biological validation'},'external':'real perturbation datasets, verify biological units before inference','report':'all portrait; no old network projection'},indent=2),encoding='utf-8')
def fit(j):
 arm,sam,d=j;py=sys.executable
 for stage,cmd in [('grn',[py,'-m','pyscenic.cli.arboreto_with_multiprocessing',str(d/'train.loom'),str(d/'tfs.txt'),'--method','grnboost2','--num_workers','10','--seed','17','-o',str(d/'adj.tsv')]),('ctx',[str(Path(py).parent/'pyscenic'),'ctx',str(d/'adj.tsv'),*map(str,dbs),'--annotations_fname',str(anno),'--expression_mtx_fname',str(d/'train.loom'),'--mode','custom_multiprocessing','--num_workers','10','--output',str(d/'regulons.csv')])]:
  marker=d/(stage+'_exit.json')
  if marker.exists() and json.loads(marker.read_text())['returncode']==0:continue
  print('START',arm,sam,stage,flush=True)
  with (d/(stage+'.log')).open('wb') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
  marker.write_text(json.dumps({'returncode':p.returncode,'command':cmd,'finished':time.time()}),encoding='utf-8')
  if p.returncode:raise RuntimeError(str((arm,sam,stage,p.returncode)))
  print('DONE',arm,sam,stage,flush=True)
 return {'arm':arm,'left_out':sam,'status':'success'}
if __name__=='__main__':
 (N/'supervisor_pid.json').write_text(json.dumps({'pid':os.getpid(),'workers_per_job':10,'parallel_jobs':8}),encoding='utf-8')
 with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex: results=list(ex.map(fit,jobs))
 (N/'complete.json').write_text(json.dumps({'status':'success','folds':results},indent=2),encoding='utf-8')
 print('NETWORK_LOO_COMPLETE',flush=True)
