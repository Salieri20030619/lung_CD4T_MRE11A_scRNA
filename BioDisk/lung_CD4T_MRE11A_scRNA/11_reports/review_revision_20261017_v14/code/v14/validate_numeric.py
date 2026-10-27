from pathlib import Path
import pandas as pd,numpy as np,json,itertools,hashlib
from scipy.stats import rankdata,false_discovery_control
W=Path('C:/Users/31937/lung_review_20261017_v14');O=W/'revision_v14';Q=O/'supplement_v14';S=['Saline_1','Saline_2','Saline_3','Silica_1','Silica_2','Silica_3']
rd=lambda n:pd.read_csv(Q/n,sep='\t')
a=rd('S21_complete_TF_fixed.tsv');checks=rd('S21_BH_recalculation_audit.tsv').to_dict('records')
for fam,b in a.groupby('family'):
 ok=b.P.notna();err=float(np.max(abs(false_discovery_control(b.loc[ok,'P'].to_numpy())-b.loc[ok,'BH'].to_numpy())));checks.append(dict(table='TF_fixed',family=fam,rows=len(b),valid_tests=int(ok.sum()),max_abs_BH_error=err,status='PASS' if err<1e-10 else 'MISMATCH'))
pd.DataFrame(checks).to_csv(Q/'S21_BH_recalculation_audit.tsv',sep='\t',index=False)
# Recompute four key correlations independently from packaged six-mouse inputs.
act=rd('S22_TF_all_regulon_sample_activity.tsv.gz');gene=rd('S22_gene_sample_expression.tsv');xf=rd('S22_Core_fixed_composition_TF_inputs.tsv').set_index('sample').loc[S,'Stat1(+)'].to_numpy();yf=rd('S22_fixed_composition_gene_inputs.tsv').set_index('sample').loc[S,'GSDMD'].to_numpy();x=act[(act.scope=='Core')&(act.seed==17)&(act.tf=='STAT1')].set_index('sample').loc[S,'mean_auc'].to_numpy();y=gene[(gene.lineage=='CD4')&(gene.scope=='Core')&(gene.gene=='GSDMD')].set_index('sample').loc[S,'CP10k'].to_numpy()
rows=[]
for composition,xx,yy,fn in [('Original',x,y,'S21_complete_STAT1_gene.tsv'),('Fixed',xf,yf,'S21_complete_STAT1_fixed.tsv')]:
 for adj in [False,True]:
  rx=rankdata(xx);ry=rankdata(yy)
  for inds in ([np.arange(3),np.arange(3,6)] if adj else [np.arange(6)]):rx[inds]-=rx[inds].mean();ry[inds]-=ry[inds].mean()
  perms=(list(a)+list(b) for a in itertools.permutations(range(3)) for b in itertools.permutations(range(3,6))) if adj else itertools.permutations(range(6));den=np.linalg.norm(rx)*np.linalg.norm(ry);r=rx@ry/den;null=[rx[list(ix)]@ry/den for ix in perms];p=sum(abs(v)>=abs(r)-1e-12 for v in null)/len(null);orig=rd(fn);orig=orig[(orig.y=='GSDMD')&(orig.seed==17)&(orig.model==('ConditionAdjusted' if adj else 'Pooled'))&orig.scope.isin(['Core','Core_fixed'])].iloc[0];assert abs(r-orig.rho)<1e-10 and abs(p-orig.P)<1e-10
  rows.append(dict(composition=composition,model=orig.model,rho=r,P=p,BH=orig.BH,samples=';'.join(S),allocations=len(null),status='PASS'))
pd.DataFrame(rows).to_csv(Q/'S21_key_exact_P_recalculation.tsv',sep='\t',index=False)
vals=np.arange(6);n=[]
for ix in itertools.combinations(range(6),3):
 mask=np.zeros(6,bool);mask[list(ix)]=True;n.append(abs(vals[mask].mean()-vals[~mask].mean()))
assert sum(np.isclose(n,max(n)))/len(n)==.1
# GO not used in the original seven-program TF family; report overlap explicitly.
e=rd('S25_reconstructed_focal_regulon_members.tsv');go=set(pd.read_csv(O/'supplement/S20_GO_autophagy_measured_members.tsv',sep='\t').gene);ctrl=set(pd.read_csv(O/'supplement/S20_GO_autophagy_control_members.tsv',sep='\t').gene);out=[]
for (seed,tf),b in e.groupby(['seed','tf_canonical']):
 tg=set(b.target_canonical);out.append(dict(seed=seed,tf=tf,GO_member_overlap=len(tg&go),GO_control_overlap=len(tg&ctrl),GO_member_genes=';'.join(sorted(tg&go)),GO_control_genes=';'.join(sorted(tg&ctrl)),used_in_original_TF_program_tests=False))
pd.DataFrame(out).to_csv(Q/'S25_GO_overlap_not_in_TF_test_family.tsv',sep='\t',index=False)
a=rd('S23_targeted_proximal_associations.tsv');a['status']=np.where(a.tf.isin(['STAT5A','STAT5B']),'NOT_EVALUABLE_primary_seed17_regulon_absent',a.status);a.loc[(a.sender=='AMlike')&(a.ligand=='IFNB1'),'status']='NOT_EVALUABLE_sender_constant';a.to_csv(Q/'S23_targeted_proximal_associations.tsv',sep='\t',index=False,na_rep='NA')
# Baseline naming map preserves hashes; six remote names have encoding damage.
r=json.loads((W/'analysis/provenance/baseline_remote_sha256.json').read_text(encoding='utf-8'));local={hashlib.sha256(f.read_bytes()).hexdigest():f.relative_to(W/'baseline_v13').as_posix() for f in (W/'baseline_v13').rglob('*') if f.is_file()};mapping=[dict(remote_path=p,local_path=local[h],sha256=h) for p,h in r.items() if p!=local[h]];pd.DataFrame(mapping).to_csv(O/'provenance/v14/baseline_filename_map.tsv',sep='\t',index=False,encoding='utf-8')
status={'BH_families_checked':len(checks),'BH_all_pass':all(z['status']=='PASS' for z in checks),'key_exact_tests_reproduced':len(rows),'baseline_files_content_identical':len(r),'baseline_filename_mapping':len(mapping),'minimum_group_exact_P':.1,'inference_scope':'BH and four focal exact tests; no full upstream re-run'}
(O/'provenance/v14/local_numeric_validation.json').write_text(json.dumps(status,indent=2),encoding='utf-8');print(json.dumps(status))