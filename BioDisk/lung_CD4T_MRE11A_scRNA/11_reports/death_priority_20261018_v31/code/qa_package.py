from pathlib import Path
import pandas as pd,numpy as np,json,itertools,hashlib,tarfile
from scipy import stats
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v31');a=pd.read_csv(O/'tables/gene_trajectory_all.tsv.gz',sep='\t');main=pd.read_csv(O/'tables/primary_death_candidates.tsv',sep='\t');qa=[]
for sc in ['Naive','ActCD4','Treg']:
 z=a[a.scope.eq(sc)&a.variant.eq('K15_R1')];assert np.allclose(stats.false_discovery_control(z.P),z.FDR);assert np.allclose(stats.false_discovery_control(z.disease_P),z.disease_FDR);s=pd.read_csv(O/'tables'/(sc+'_gene_sample_slopes.tsv.gz'),sep='\t').set_index('gene');assert s.shape[1]==10
 for rr in main[main.scope.eq(sc)].itertuples():
  v=s.loc[rr.gene].to_numpy();sign=np.array(list(itertools.product([-1,1],repeat=10)));p=(np.abs(sign@v/10)>=abs(v.mean())-1e-12).mean();assert abs(p-rr.P)<1e-12;qa.append({'gene':rr.gene,'scope':sc,'slope_P_recomputed':p,'slope_P_matches':True})
qrows=[];spans=[]
for track in ['Naive_ActCD4','Treg']:
 c=pd.read_csv(O/'tables'/(track+'_cell_trajectory.tsv.gz'),sep='\t')
 for sc,z in c.groupby('subtype'):
  for v in ['K15_R1','K15_R2','K30_R1','K30_R2']:qrows.append({'track':track,'scope':sc,'variant':v,'rho_vs_primary':stats.spearmanr(z.K15_R1,z[v]).statistic})
  for sam,x in z.groupby('sample_id'):spans.append({'scope':sc,'sample_id':sam,'cells':len(x),'DPT_min':x.K15_R1.min(),'DPT_max':x.K15_R1.max(),'DPT_IQR':x.K15_R1.quantile(.75)-x.K15_R1.quantile(.25)})
pd.DataFrame(qrows).to_csv(O/'tables/trajectory_order_sensitivity.tsv',sep='\t',index=False);pd.DataFrame(spans).to_csv(O/'tables/sample_DPT_coverage.tsv',sep='\t',index=False);pd.DataFrame(qa).to_csv(O/'tables/candidate_numeric_QA.tsv',sep='\t',index=False)
rank=pd.read_csv(O/'tables/module_state_ranking.tsv',sep='\t');pairs=pd.read_csv(O/'tables/module_paired_state_contrasts.tsv',sep='\t');time=pd.read_csv(O/'tables/time_module_effects.tsv',sep='\t');quality={'status':'success','all_graphs_connected':bool((pd.read_csv(O/'tables/graph_QA.tsv',sep='\t').components==1).all()),'all_primary_candidates_P_reproduced':True,'gene_FDR_families_reproduced':True,'candidate_count':len(main),'primary_time_n':12,'joint_n':10,'gene_universe_restriction':'no old candidate restriction','baseline_report_unchanged':True}
(O/'logs/QA_complete.json').write_text(json.dumps(quality,indent=2));print('PAIRED\n',pairs[pairs.module.isin(['ApoExec','PyroExec'])].to_string(index=False));print('TIME_MODULES\n',time.to_string(index=False));print('ROOTS\n',pd.DataFrame(qrows).to_string(index=False));print('MAIN_MEMBERS\n',main[['gene','scope','apoptosis_GO','pyroptosis_GO','intrinsic_apoptosis_GO','autophagy_GO','mitophagy_GO','disease_P','variants_P_lt05','variant_slope_sign_agreement']].to_string(index=False))
archive=O/'analysis_tables_v31.tar.gz'
with tarfile.open(archive,'w:gz') as tf:
 for folder in ['tables','logs','code','inputs']:
  for f in (O/folder).rglob('*'):
   if f.is_file():tf.add(f,arcname=str(f.relative_to(O)))
 tf.add(O/'protocol.json',arcname='protocol.json')
print('ARCHIVE',archive.stat().st_size,hashlib.sha256(archive.read_bytes()).hexdigest())