from pathlib import Path
import os,json
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import pandas as pd,numpy as np
from scipy.io import mmread
from threadpoolctl import threadpool_limits
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');A=P/'02_analysis';O=P/'11_reports/human_gap_diagnostics_20261005_v1';R=O/'results';V=P/'11_reports/human_alignment_20261005_v1/results'
def rd(p):return pd.read_csv(p,sep='\t')
def wr(d,n):pd.DataFrame(d).to_csv(R/n,sep='\t',index=False,na_rep='NA')
gs='IL15 IL15RA IL2RB IL2RG IKZF2 STAT1 FOXP3 CTLA4 TIGIT CD4 CD8A CD8B TRAC NKG7 GNLY MS4A1 LYZ'.split();allcells=[];audit=[]
for lin,pre,col in [('CD4',A/'mechanism_step2_20260927_154919_EsN4gk/input/human_','symbol'),('Mac',A/'scenic3_20261001_123452_Kvp9h5/input/human_mac_','gene')]:
 m=rd(str(pre)+'cells.tsv');g=rd(str(pre)+'genes.tsv')[col].astype(str);print('READ_RAW',lin,flush=True)
 with threadpool_limits(limits=2):x=mmread(str(pre)+'counts.mtx').tocsr()
 assert x.shape==(len(g),len(m)) and np.all(x.data>=0) and np.all(x.data==np.floor(x.data))
 lib=np.asarray(x.sum(axis=0)).ravel();assert np.allclose(lib,m.library);m['n_genes']=np.asarray((x>0).sum(axis=0)).ravel();m['raw_mt_fraction']=np.asarray(x[g.str.startswith('MT-').to_numpy()].sum(axis=0)).ravel()/np.maximum(lib,1)
 for gene in gs:
  ix=np.flatnonzero(g.eq(gene));audit.append(dict(lineage=lin,gene=gene,raw_features=len(ix),symbols=';'.join(g.iloc[ix])));m['count_'+gene]=np.asarray(x[ix].sum(axis=0)).ravel() if len(ix) else np.nan
 del x
 if lin=='Mac':
  rev=rd(P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release/prior_review_repaired/scenic3_review/results/reviewed_cells.tsv.gz');m=m.merge(rev[['cell','review_scope','label','r4']],on='cell',validate='one_to_one')
 else:
  cs=rd(A/'mechanism_step2_20260927_154919_EsN4gk/objects/human_cell_scores.tsv.gz');print('CELL_SCORE_COLS',list(cs),flush=True);extra=[c for c in cs if c.startswith('score_') or c in ['bg','mt']];m=m.merge(cs[['cell']+extra],on='cell',validate='one_to_one')
 m['lineage']=lin;allcells.append(m);print('RAW_DONE',lin,len(m),flush=True)
cells=pd.concat(allcells,ignore_index=True);wr(cells,'measured_selected_genes_cells_SERVER_ONLY.tsv.gz');wr(audit,'raw_feature_mapping_audit.tsv')
des=rd(V/'stage1_donor_design.tsv');scopes=['Primary','Core','CD4-CCR7','CD4-GZMA','CD4-IL7R','Treg','CD4state'];agg=[];bins=[];clusters=[]
for lin in ['CD4','Mac']:
 m=cells[(cells.lineage==lin)&cells['sample'].isin(des['sample'])]
 for scope in (scopes if lin=='CD4' else ['MacAM','MonoMac','MacWide']):
  if lin=='CD4':z=m if scope=='Primary' else m[m.core.eq(True)] if scope=='Core' else m[m.subtype==scope];chosen=['IL2RB','IL2RG','IKZF2','STAT1']
  else:z=m[m.review_scope.isin(['MacAM','MonoMac'])] if scope=='MacWide' else m[m.review_scope==scope];chosen=['IL15','IL15RA']
  for sample,v in z.groupby('sample'):
   tag=dict(lineage=lin,scope=scope,sample=sample,group=v.group.iloc[0],series=v.series.iloc[0],cells=len(v),mean_depth=v.library.mean(),median_depth=v.library.median(),n_genes_mean=v.n_genes.mean(),mt_mean=v.raw_mt_fraction.mean())
   for gene in chosen:
    cnt=v['count_'+gene].to_numpy();agg.append(dict(tag,gene=gene,counts=cnt.sum(),detected=int((cnt>0).sum()),pct=(cnt>0).mean(),mean_CP10k=np.mean(cnt/v.library*1e4),mean_log1pCP10k=np.mean(np.log1p(cnt/v.library*1e4))))
   q=pd.qcut(v.library.rank(method='first'),4,labels=False,duplicates='drop') if len(v)>=4 else pd.Series(np.zeros(len(v)),index=v.index)
   for quart,t in v.groupby(q):
    for gene in chosen:bins.append(dict(tag,quartile=int(quart)+1,bin_cells=len(t),gene=gene,bin_mean_depth=t.library.mean(),detected=int((t['count_'+gene]>0).sum()),pct=(t['count_'+gene]>0).mean()))
  if lin=='Mac':
   for (label,group,sample),v in z.groupby(['label','group','sample']):
    for gene in chosen:clusters.append(dict(scope=scope,label=label,group=group,sample=sample,gene=gene,cells=len(v),detected=int((v['count_'+gene]>0).sum()),counts=v['count_'+gene].sum()))
agg=pd.DataFrame(agg);wr(agg,'raw_donor_gene_measurements.tsv.gz');wr(bins,'within_donor_depth_quartiles.tsv.gz');wr(clusters,'macrophage_cluster_gene_coverage.tsv.gz')
old=rd(V/'alignment_gene_donors.tsv.gz');old=old[old['filter']=='All'];check=agg.merge(old,on=['lineage','scope','sample','gene'],suffixes=('_raw','_old'));assert len(check)==len(agg);assert np.allclose(check.pct_raw,check.pct_old,atol=1e-12);assert np.allclose(check.mean_CP10k_raw,check.mean_CP10k_old,atol=1e-10)
wr(check[['lineage','scope','sample','gene','pct_raw','pct_old','mean_CP10k_raw','mean_CP10k_old']],'raw_vs_previous_mapping_validation.tsv.gz')
# Read only existing target-TF cell AUC columns, no AUCell/network re-running.
md=rd(A/'scenic1_20260929_103513_xRoEnr/Human/All/cells.tsv');cv=cells[cells.lineage=='CD4'];auc=[]
for seed in [17,29,43]:
 for variant in ['Original','DropOverlap']:
  f=A/f'scenic1_20260929_103513_xRoEnr/Human/All/seed{seed}/auc.csv' if variant=='Original' else A/f'scenic2_20261001_092655_ZKNWD6/objects/Human_All_{seed}_DropOverlap_auc.csv.gz'
  header=pd.read_csv(f,nrows=0);cols=[header.columns[0]]+[c for c in header if c.startswith(('IKZF2(','STAT1('))];d=pd.read_csv(f,usecols=cols).rename(columns={cols[0]:'cell'});d=d.melt(id_vars='cell',var_name='regulon',value_name='auc');d['seed']=seed;d['variant']=variant;d['tf']=d.regulon.str.split('(',regex=False).str[0];auc.append(d)
auc=pd.concat(auc,ignore_index=True);auc=auc.merge(cv,on='cell',validate='many_to_one');wr(auc,'existing_TF_auc_selected_cell_context_SERVER_ONLY.tsv.gz')
(O/'provenance/raw_complete.json').write_text(json.dumps({'raw_cells':len(cells),'mapping_comparisons':len(check),'max_pct_error':float(abs(check.pct_raw-check.pct_old).max()),'auc_cell_rows':len(auc)},indent=2));print('GAP_RAW_COMPLETE',len(cells),len(auc),flush=True)
