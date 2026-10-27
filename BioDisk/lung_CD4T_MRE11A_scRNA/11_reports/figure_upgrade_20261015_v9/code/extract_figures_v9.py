from pathlib import Path
import sys,json,hashlib
import pandas as pd,numpy as np,anndata as ad
from scipy.io import mmread
O=Path(sys.argv[1]);P=O.parents[1];M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF';V=P/'11_reports/integrated_macrophage_CD4_20261015_v4';manifest=[]
def use(p):p=Path(p);manifest.append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()});return p
G=['CD3D','CD3E','TRAC','CCR7','IL7R','FOXP3','IL2RA','TNFRSF9','LAG3','CXCR6','IKZF2','STAT1','GSDMD','C1QA','C1QB','FOLR2','SPP1','MRC1','FABP4','MARCO','CXCL16','C3','IL15','IL18','TNF','IL27','EBI3','IFNG']
source={'TNF':['TNF'],'IL15':['IL15','IL15RA'],'IL18':['IL18'],'CXCL16':['CXCL16'],'C3':['C3'],'IL27':['IL27','EBI3'],'IFNG':['IFNG']};rec={'TNF':['TNFRSF1A','TNFRSF1B'],'IL15':['IL2RB','IL2RG'],'IL18':['IL18R1','IL18RAP'],'CXCL16':['CXCR6'],'C3':['C3AR1'],'IL27':['IL27RA','IL6ST'],'IFNG':['IFNGR1','IFNGR2']};joint=[];markers=[]
for lineage,obj,annot in [('CD4','stage5_20260905_222457_a2dApT/objects/CD4_core.h5ad','Subtype'),('Mac','stage11_20260906_100314_pALwAB/objects/macrophages_core.h5ad','MacType')]:
 md=pd.read_csv(use(M/'comm_input'/f'{lineage}_cells.tsv'),sep='\t');gn=pd.read_csv(use(M/'comm_input'/f'{lineage}_genes.tsv'),sep='\t');x=mmread(use(M/'comm_input'/f'{lineage}_counts.mtx.gz')).T.tocsr();gi=pd.Index(gn.canonical);h=ad.read_h5ad(use(P/'02_analysis/GSE314943'/obj),backed='r');assert set(md.cell)==set(h.obs_names);ix=h.obs_names.get_indexer(md.cell);assert np.array_equal(h.obs.iloc[ix][annot].astype(str).to_numpy(),md.subtype.astype(str).to_numpy());um=np.asarray(h.obsm['X_umap_Harmony'])[ix];h.file.close();out=md.copy();out['UMAP1']=um[:,0];out['UMAP2']=um[:,1]
 for g in G:
  j=gi.get_indexer([g])[0];out[g]=np.log1p(x[:,j].toarray().ravel()*1e4/md.library.to_numpy()) if j>=0 else np.nan
 if lineage=='CD4':
  a=pd.read_csv(use(P/'02_analysis/scenic2_20261001_092655_ZKNWD6/objects/Mouse_All_17_DropOverlap_auc.csv.gz'),index_col=0).loc[md.cell];out['IKZF2_AUC']=a['Ikzf2(+)'].to_numpy();out['STAT1_AUC']=a['Stat1(+)'].to_numpy()
 out.to_csv(O/'inputs'/f'{lineage}_embedding_expression.tsv.gz',sep='\t',index=False,na_rep='NA')
 for (state,sample),part in md.groupby(['subtype','sample']):
  mask=(md.subtype.eq(state)&md['sample'].eq(sample)).to_numpy();n=int(mask.sum());group=part.group.iloc[0]
  for g in G:
   j=gi.get_indexer([g])[0];z=x[mask,j].toarray().ravel() if j>=0 else None;markers.append(dict(lineage=lineage,state=state,sample=sample,group=group,gene=g,cells=n,mean_logCP10k=np.mean(np.log1p(z*1e4/md.loc[mask,'library'])) if z is not None else np.nan,pct=float((z>0).mean()*100) if z is not None else np.nan))
  for lig,gs in (source if lineage=='Mac' else rec).items():
   ii=gi.get_indexer(gs);ok=(ii>=0).all();z=x[mask][:,ii].toarray()>0 if ok else None;rule='any' if lineage=='CD4' and lig=='TNF' else 'all';pos=(z.any(axis=1) if rule=='any' else z.all(axis=1)) if ok else None;joint.append(dict(lineage=lineage,state=state,sample=sample,group=group,ligand=lig,components=';'.join(gs),rule=rule,cells=n,positive_cells=int(pos.sum()) if ok else np.nan,pct=float(pos.mean()*100) if ok else np.nan))
pd.DataFrame(markers).to_csv(O/'inputs/marker_sample_summary.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(joint).to_csv(O/'inputs/communication_complex_expression.tsv',sep='\t',index=False,na_rep='NA')
# Reconstruct exactly the overlap-excluded SCENIC graph; verify retained target counts.
E=P/'02_analysis/scenic1_20260929_103513_xRoEnr/Mouse/All/regulon_edges.tsv';e=pd.read_csv(use(E),sep='\t');e=e[e.tf.str.upper().isin(['IKZF2','STAT1'])].copy();e['target_upper']=e.target.str.upper().replace({'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'});mods=pd.read_csv(use(M/'tables/module_definition.tsv'),sep='\t');exclude=set(mods.gene)|{'IL1R2','PENK','TNFRSF9','LAG3','TTC39C','FOXP1','MRE11'};e['excluded']=e.target_upper.isin(exclude)|(e.target_upper==e.tf.str.upper());coverage=pd.read_csv(use(V/'03_receiver/TF_regulon_coverage.tsv'),sep='\t');checks=[]
for (seed,reg),b in e.groupby(['seed','regulon']):
 c=coverage[(coverage.species=='Mouse')&(coverage['filter']=='All')&(coverage.seed==seed)&(coverage.variant=='DropOverlap')&(coverage.regulon==reg)].iloc[0];kept=b.loc[~b.excluded,'target'].nunique();assert kept==c.retained_genes,(seed,reg,kept,c.retained_genes);checks.append(dict(seed=int(seed),regulon=reg,retained_genes=kept))
e.to_csv(O/'inputs/SCENIC_focal_edges_all_seeds.tsv',sep='\t',index=False);keep=e[~e.excluded];support=keep.groupby(['regulon','target_upper']).seed.nunique().rename('seed_support').reset_index();keep=keep.merge(support,on=['regulon','target_upper']);keep.to_csv(O/'inputs/SCENIC_focal_edges_retained.tsv',sep='\t',index=False);print('REGULON_PARITY',checks)
(O/'provenance/input_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8');(O/'provenance/extraction_checks.json').write_text(json.dumps({'cell_roster_and_labels':'PASS','coordinates':'existing Harmony UMAP; no new embedding','SCENIC_retained_target_count_parity':checks}),encoding='utf-8');print('EXTRACTION_COMPLETE')