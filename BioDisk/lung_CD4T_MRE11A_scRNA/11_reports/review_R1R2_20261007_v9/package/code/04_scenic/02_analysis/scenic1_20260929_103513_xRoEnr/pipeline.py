import os, sys, json, hashlib, subprocess, tarfile, re
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

ROOT=Path(os.environ['ROOT']); OUT=Path(os.environ['OUT']); CACHE=Path(os.environ['CACHE'])
P1=Path(os.environ['PREV1']); P2=Path(os.environ['PREV2']); P7=Path(os.environ['PREV7'])
SEEDS=(17,29,43); WORKERS=int(os.environ.get('WORKERS','8')); CAP=250
BG=set('SFTPA1 SFTPA2 SFTPB SFTPC SFTPD SCGB1A1 SCGB3A1 SCGB3A2 ABCA3 SLC34A2 NAPSA LYZ LYZ1 LYZ2 S100A8 S100A9 CSF1R LST1 TYROBP FCER1G'.split())
def bgmask(gs):
 return np.array([g.upper() in BG or bool(re.match(r'^(IGH[VDJMAEG]|IGK[VJC]|IGL[VJC])',g.upper())) for g in gs])
def sha(p,algo='sha256'):
 h=hashlib.new(algo)
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def wr(d,p):
 pd.DataFrame(d).to_csv(p,sep='\t',index=False)
def run(cmd):
 print('RUN='+' '.join(map(str,cmd)),flush=True)
 subprocess.run(list(map(str,cmd)),check=True)
def download(url,p):
 if p.is_file() and p.stat().st_size>0:return
 part=Path(str(p)+'.part')
 run(['curl','-L','--fail','--show-error','--silent','--retry','12','--retry-all-errors','--retry-delay','10','--connect-timeout','30','--speed-time','180','--speed-limit','1024','--continue-at','-','-o',part,url])
 assert part.stat().st_size>0
 part.replace(p)
def resource(species):
 from ctxcore.rnkdb import FeatherRankingDatabase
 assembly,latin,anno,tf=('mm10','mus_musculus','mgi','allTFs_mm.txt') if species=='Mouse' else ('hg38','homo_sapiens','hgnc','allTFs_hg38.txt')
 base='https://resources.aertslab.org/cistarget/'
 d=CACHE/species;d.mkdir(exist_ok=True,parents=True);rows=[];dbs=[]
 for win in ['10kbp_up_10kbp_down','500bp_up_100bp_down']:
  name=f'{assembly}_{win}_full_tx_v10_clust.genes_vs_motifs.rankings.feather';p=d/name
  url=base+f'databases/{latin}/{assembly}/refseq_r80/mc_v10_clust/gene_based/'+name
  download(url+'.sha1sum.txt',Path(str(p)+'.sha1sum.txt'));download(url,p)
  expect=Path(str(p)+'.sha1sum.txt').read_text().split()[0]
  assert len(expect)==40 and sha(p,'sha1')==expect, f'CHECKSUM_MISMATCH: {p}; preserve file for diagnosis'
  db=FeatherRankingDatabase(str(p),name=win);assert db.total_genes>10000
  dbs.append(p);rows.append(dict(file=str(p),url=url,sha256=sha(p),official_sha1=expect,genes=db.total_genes))
 for name,pathpart in [(f'motifs-v10nr_clust-nr.{anno}-m0.001-o0.0.tbl','motif2tf/'),(tf,'tf_lists/')]:
  p=d/name;url=base+pathpart+name;download(url,p);rows.append(dict(file=str(p),url=url,sha256=sha(p)))
 wr(rows,OUT/'tables'/f'{species}_resources.tsv')
 return dbs,d/f'motifs-v10nr_clust-nr.{anno}-m0.001-o0.0.tbl',d/tf

def load_input(species):
 if species=='Mouse':
  pref=P1/'comm_input/CD4';x=mmread(str(pref)+'_counts.mtx.gz').tocsr().T.tocsr()
  g=pd.read_csv(str(pref)+'_genes.tsv',sep='\t');m=pd.read_csv(str(pref)+'_cells.tsv',sep='\t')
  names=g.mouse_symbol.astype(str).to_numpy();canon=g.canonical.astype(str).to_numpy()
  assert len(m)==2311 and m['sample'].nunique()==6 and set(m.group)=={'Saline','Silica'}
  m['series']='GSE314943';m['core']=True
 else:
  x=mmread(P2/'input/human_counts.mtx').tocsr().T.tocsr()
  g=pd.read_csv(P2/'input/human_genes.tsv',sep='\t');m=pd.read_csv(P2/'input/human_cells.tsv',sep='\t')
  names=g.symbol.astype(str).to_numpy();canon=np.char.upper(names.astype(str))
  c=m.core.astype(str).str.upper();assert c.isin(['TRUE','FALSE']).all();m['core']=c.eq('TRUE')
  assert len(m)==6523
 assert x.shape==(len(m),len(g)) and m.cell.is_unique and m[['cell','sample','group','subtype','library']].notna().all().all()
 assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.rint(x.data))
 assert (m.library>0).all() and np.all(np.asarray(x.sum(axis=1)).ravel()<=m.library.to_numpy()+1e-6)
 # Calculate background on the frozen canonical gene universe before symbol remapping.
 m['bg']=np.asarray(x[:,bgmask(canon)].sum(axis=1)).ravel()/m.library.to_numpy()
 valid=np.array([bool(v.strip()) and v not in ('nan','NA') and not v.upper().startswith('ENS') for v in names])
 u,inv=np.unique(names[valid],return_inverse=True)
 proj=sparse.csr_matrix((np.ones(valid.sum()),(np.flatnonzero(valid),inv)),shape=(x.shape[1],len(u)))
 x=(x@proj).tocsr();wr(dict(original=names,retained=valid),OUT/'tables'/f'{species}_symbol_audit.tsv')
 assert not pd.Index(u).has_duplicates
 return x,u,m

def select_cells(species,m,arm):
 sel=m.core.to_numpy().copy()
 if arm=='BG1':sel &= m.bg.to_numpy()<=.01
 if species=='Human':
  aud=pd.read_csv(P7/'tables/human_sample_audit.tsv',sep='\t')
  a=aud[(aud.scope=='Core')&(aud['filter']==arm)]
  assert not a['sample'].duplicated().any()
  sel &= m['sample'].isin(a['sample']).to_numpy()
  cnt=m.loc[sel].groupby('sample',observed=True).size()
  assert set(cnt.index)==set(a['sample'])
  assert all(cnt[r['sample']]==r['cells'] for _,r in a.iterrows()),'Frozen human cell count mismatch'
 else:
  cnt=m.loc[sel].groupby('sample',observed=True).size();ok=cnt[cnt>=20].index
  sel &= m['sample'].isin(ok).to_numpy()
 z=m.loc[sel].reset_index(drop=True)
 assert z['sample'].nunique()>=6 and z.group.nunique()==2
 assert (z.groupby('sample',observed=True).size()>=20).all()
 assert not z[['sample','group','series']].drop_duplicates()['sample'].duplicated().any()
 return sel,z

def capped_indices(m,seed):
 rng=np.random.default_rng(seed);ids=[]
 for ss,v in m.groupby('sample',observed=True,sort=True):
  ids.extend(rng.choice(v.index.to_numpy(),size=min(CAP,len(v)),replace=False).tolist())
 return np.array(sorted(ids),dtype=int)

def prepare(species,arm,x,genes,m,dbs,tflist):
 import loompy
 from ctxcore.rnkdb import FeatherRankingDatabase
 sel,z=select_cells(species,m,arm);xx=x[sel]
 d=OUT/species/arm;d.mkdir(parents=True,exist_ok=True)
 z.to_csv(d/'cells.tsv',sep='\t',index=False)
 universe=set(FeatherRankingDatabase(str(dbs[0]),name='db').genes)
 universe &= set(FeatherRankingDatabase(str(dbs[1]),name='db').genes)
 det=np.asarray((xx>0).sum(axis=0)).ravel();threshold=max(10,int(np.ceil(.01*len(z))))
 keep=(det>=threshold)&np.array([g in universe for g in genes]);gg=genes[keep];xx=xx[:,keep]
 tf=[t.strip() for t in tflist.read_text().splitlines() if t.strip() in set(gg)]
 assert len(gg)>=1000 and len(tf)>=100, f'Insufficient genes/TF overlap {species} {arm}: {len(gg)} {len(tf)}'
 wr(dict(gene=genes,detected_cells=det,in_database=[g in universe for g in genes],retained=keep),d/'gene_filter.tsv')
 (d/'tfs.txt').write_text('\n'.join(tf)+'\n')
 # CP10k using original full library; no integrated assay, imputation, or outcome-selected genes.
 vals=xx.multiply(1e4/z.library.to_numpy()[:,None]).tocsr().astype('float32');vals.data=np.log1p(vals.data)
 rows=[]
 for seed in SEEDS:
  idx=capped_indices(z,seed);sd=d/f'seed{seed}';sd.mkdir(exist_ok=True)
  tr=vals[idx];var=np.asarray(tr.power(2).mean(axis=0)).ravel()-np.asarray(tr.mean(axis=0)).ravel()**2
  k=var>1e-8;tg=gg[k];tfs=[t for t in tf if t in set(tg)];assert len(tfs)>=100
  (sd/'tfs.txt').write_text('\n'.join(tfs)+'\n');z.iloc[idx].to_csv(sd/'train_cells.tsv',sep='\t',index=False)
  loompy.create(str(sd/'train.loom'),tr[:,k].T.tocsc(),{'Gene':tg.astype(str)},{'CellID':z.cell.iloc[idx].astype(str).to_numpy()})
  for ss,n in z.iloc[idx].groupby('sample',observed=True).size().items():rows.append(dict(seed=seed,sample=ss,cells=int(n),genes=int(k.sum()),tfs=len(tfs)))
 # All retained cells are scored; cap is used only for network discovery.
 loompy.create(str(d/'score.loom'),vals.T.tocsc(),{'Gene':gg.astype(str)},{'CellID':z.cell.astype(str).to_numpy()})
 wr(rows,d/'training_audit.tsv')
 print(f'INPUT_PASS={species} {arm} cells={len(z)} samples={z["sample"].nunique()} genes={len(gg)}',flush=True)
 return d

def fit(d,dbs,anno):
 py=sys.executable;cli=Path(py).parent/'pyscenic'
 for seed in SEEDS:
  sd=d/f'seed{seed}'
  print(f'SCENIC_RUN_START={d.parent.name}/{d.name}/seed{seed}',flush=True)
  run([py,'-m','pyscenic.cli.arboreto_with_multiprocessing',sd/'train.loom',sd/'tfs.txt','--method','grnboost2','--num_workers',WORKERS,'--seed',seed,'-o',sd/'adj.tsv'])
  run([cli,'ctx',sd/'adj.tsv',*dbs,'--annotations_fname',anno,'--expression_mtx_fname',sd/'train.loom','--mode','custom_multiprocessing','--num_workers',WORKERS,'--output',sd/'regulons.csv'])
  run([cli,'aucell',d/'score.loom',sd/'regulons.csv','--num_workers',WORKERS,'--seed',seed,'--output',sd/'auc.csv'])
  print(f'SCENIC_RUN_DONE={d.parent.name}/{d.name}/seed{seed}',flush=True)

def summarize(d):
 from pyscenic.cli.utils import load_signatures
 m=pd.read_csv(d/'cells.tsv',sep='\t');edges=[];acts=[];stat=[]
 for seed in SEEDS:
  sd=d/f'seed{seed}';regs=load_signatures(str(sd/'regulons.csv'))
  assert len(regs)>0,'No motif-supported regulons'
  for r in regs:
   for g,w in r.gene2weight.items():edges.append(dict(seed=seed,regulon=r.name,tf=r.transcription_factor,target=g,weight=w))
  a=pd.read_csv(sd/'auc.csv',index_col=0);a.index=a.index.astype(str)
  assert a.index.is_unique and set(a.index)==set(m.cell) and a.columns.is_unique
  a=a.loc[m.cell];assert np.isfinite(a.to_numpy()).all()
  for scope in ['Core']+sorted(m.subtype.unique()):
   md=m if scope=='Core' else m[m.subtype==scope]
   for ss,v in md.groupby('sample',observed=True):
    if len(v)<20:continue
    mu=a.loc[v.cell].mean()
    for reg,value in mu.items():acts.append(dict(seed=seed,scope=scope,sample=ss,group=v.group.iloc[0],series=v.series.iloc[0],cells=len(v),regulon=reg,mean_auc=value))
  stat.append(dict(seed=seed,regulons=len(regs),cells=len(a)))
 e=pd.DataFrame(edges);wr(e,d/'regulon_edges.tsv')
 support=e.groupby(['regulon','tf','target'],observed=True).agg(seed_support=('seed','nunique'),mean_weight=('weight','mean')).reset_index()
 support['stable_2of3']=support.seed_support>=2;wr(support,d/'edge_stability.tsv');wr(acts,d/'sample_activity.tsv');wr(stat,d/'run_summary.tsv')
 # These are gene symbols, not orthology calls; cross-species validation is performed downstream.
 targets={'IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3','MRE11','MRE11A'}
 wr(support[support.target.str.upper().isin(targets)],d/'candidate_edges.tsv')

def main():
 for p,stage in [(P1,'mechanism_step1'),(P2,'mechanism_step2_human_and_driver_audit'),(P7,'mechanism_step7_human_TF')]:
  s=json.loads((p/'status.json').read_text());assert s['status']=='PASS' and s['stage']==stage
 OUT.joinpath('tables').mkdir(exist_ok=True)
 prepared=[]
 # Validate both inputs and all arms before expensive GRN inference.
 for species in ['Mouse','Human']:
  dbs,anno,tfs=resource(species);x,g,m=load_input(species)
  for arm in ['All','BG1']:prepared.append((prepare(species,arm,x,g,m,dbs,tfs),dbs,anno))
 print('ANALYSIS_LOGIC_PASS INPUT_COUNTS_SCOPE_SYMBOLS_AND_TRAINING',flush=True)
 for d,dbs,anno in prepared:fit(d,dbs,anno);summarize(d)
 methods='''SCENIC discovery: mouse GSE314943 Core CD4 and human frozen Step7 Core CD4 analyzed separately. All and BG1 (background fraction <=1%) are separate network runs. Human subjects/cell counts match Step7 per filter. Sample training cap 250, three seeds 17/29/43, all cells retained for AUC scoring. The cap limits large sample dominance; unequal small samples are not oversampled. Gene detection >=max(10,1% retained cells), motif database intersection; no HVG or differential gene preselection. Original-library CP10k log1p expression. GRNBoost2 multiprocessing, cisTarget v10 motif pruning using two windows, AUCell. Edge stability is repeatability across random subsampling/seeds, not independent replication. Same subtype-conditioned downstream comparisons, group/source/depth adjustment, shared-target removal, and human/mouse orthology comparisons remain downstream. No causal or direct binding claims are generated here. Step8 macrophage results remain unchanged.\n'''
 (OUT/'methods.txt').write_text(methods)
 status=dict(status='PASS',stage='lung_scenic1_network',output=str(OUT),arms=4,seeds=list(SEEDS),runs=12,pending=['sample_level_effects_and_source_adjustment','module_overlap_exclusion','cross_species_and_prior_comparison','macrophage_linkage'])
 (OUT/'status.json').write_text(json.dumps(status,indent=2))
 ret=ROOT/'tmp/lung_scenic1_return.tar.gz';part=Path(str(ret)+'.part')
 with tarfile.open(part,'w:gz') as t:
  for p in sorted(OUT.rglob('*')):
   if not p.is_file() or p.suffix in ('.loom','.pyc') or p.name in ('adj.tsv','auc.csv'):continue
   t.add(p,arcname=str(p.relative_to(OUT)))
 with tarfile.open(part) as t:
  for member in t:
   if member.isfile():
    with t.extractfile(member) as f:
     for b in iter(lambda:f.read(2**20),b''):pass
 part.replace(ret)
 print('RETURN='+str(ret),flush=True);print('RETURN_BYTES='+str(ret.stat().st_size),flush=True);print('RETURN_SHA256='+sha(ret),flush=True)
 print('LUNG_SCENIC1_PASS',flush=True)
if __name__=='__main__':main()

