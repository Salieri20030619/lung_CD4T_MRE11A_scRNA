from pathlib import Path
import os,json,itertools,math
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np
from scipy import stats
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/phenotype_candidate_discovery_20261005_v1';R=O/'results';H=P/'11_reports/human_first_20261004_v1/results';src=(O/'code/discover_compute.py').read_text();exec(src[src.index('def wr('):src.index('scores=rd(')]);rd=lambda p:pd.read_csv(p,sep='\t',low_memory=False)
(O/'provenance/DISCOVERY_AUDIT_ADDENDUM.txt').write_bytes((O/'code/DISCOVERY_AUDIT_ADDENDUM.txt').read_bytes())
models=rd(R/'complete_direct_phenotype_discovery_models.tsv.gz');extreme=models[models.rho.abs().ge(1-1e-12)&models.P.notna()];wr(extreme,'all_exact_rank_correlation_cautions.tsv');inp=rd(R/'direct_phenotype_donor_inputs_SERVER_ONLY.tsv.gz');d=inp[(inp.ligand=='TXLNA')&(inp.sender=='MacAM')&(inp.receiver=='Treg')&(inp.metric=='score_OXPHOS')].reset_index(drop=True);assert len(d)==6
rows=[]
for step in ['Unadjusted','Group']:
 baseline=corr(d,step)['rho'];count=total=0
 if step=='Unadjusted':perms=itertools.permutations(range(len(d)))
 else:
  groups=[list(v.index) for _,v in d.groupby('group')];pp=[list(itertools.permutations(ii)) for ii in groups]
  def gen():
   for z in itertools.product(*pp):
    order=np.arange(len(d))
    for ii,jj in zip(groups,z):order[ii]=jj
    yield order
  perms=gen()
 for perm in perms:
  q=d.copy();q['x']=d.x.to_numpy()[list(perm)];r=corr(q,step)['rho'];count+=abs(r)>=abs(baseline)-1e-12;total+=1
 rows.append(dict(ligand='TXLNA',receiver='Treg',metric='score_OXPHOS',sender='MacAM',step=step,n=6,baseline_rho=baseline,permutations=total,as_or_more_extreme=count,exact_P=count/total,interpretation='exhaustive diagnostic; postselection; conditional exchangeability assumption'))
wr(rows,'TXLNA_exact_small_sample_diagnostic.tsv')
ss=rd(R/'SERPINA1_fixed_source_identity_SERVER_ONLY.tsv');score=rd(H/'stage1_reused_all_program_inputs.tsv');score=score[(score.scope=='Subtype:CD4-IL7R')&(score['filter']=='All')][['sample','cells','library','score_Mitophagy']].rename(columns={'cells':'receiver_cells','library':'receiver_depth','score_Mitophagy':'y'});d=ss.merge(score,on='sample');d=d[d.sender_cells.ge(20)&d.receiver_cells.ge(20)&np.isfinite(d[['x','y','sender_depth','receiver_depth']]).all(axis=1)];wr(d,'SERPINA1_source_paired_donors_SERVER_ONLY.tsv');rows=[]
for sender,z in d.groupby('sender'):
 for step in ['Group','Source','SourceBothDepth']:rows.append(dict(ligand='SERPINA1',receiver='CD4-IL7R',metric='score_Mitophagy',sender=sender,step=step,**corr(z,step)))
q=bh(pd.DataFrame(rows),['step']);wr(q,'SERPINA1_source_identity_associations.tsv')
orig=models[(models.ligand=='SERPINA1')&(models.receiver=='CD4-IL7R')&(models.metric=='score_Mitophagy')];wr(orig,'SERPINA1_original_discovery_context.tsv');check=q[q.sender.isin(['MonoMac','MacAM'])].merge(orig,on=['ligand','receiver','metric','sender','step'],suffixes=('_new','_old'));assert np.allclose(check.rho_new,check.rho_old,equal_nan=True,atol=1e-10);assert np.allclose(check.P_new,check.P_old,equal_nan=True,rtol=1e-8)
summary=ss.groupby(['sender','group']).agg(donors=('sample','size'),cells=('sender_cells','sum'),positive=('positive','sum'),mature_marker2=('mature_marker2','sum'),mono_marker2=('mono_marker2','sum')).reset_index();wr(summary,'SERPINA1_source_identity_summary.tsv')
e=rd(H/'stage2_evidence_all_candidates.tsv.gz');wr(e[e.ligand=='SERPINA1'],'SERPINA1_all_existing_candidates.tsv');lr=rd(H/'stage2_all_LR_supported_and_unsupported.tsv.gz');wr(lr[(lr['from']=='SERPINA1')&lr.receiver.eq('CD4-IL7R')],'SERPINA1_existing_LR_context.tsv')
print('EXACT',pd.DataFrame(rows).to_string(index=False));print('SERPINA_MODELS',q[['sender','step','n','n_case','n_control','rho','P','FDR','status']].to_string(index=False));print('IDENTITY',summary.to_string(index=False));print('TXLNA',rd(R/'TXLNA_exact_small_sample_diagnostic.tsv').to_string(index=False));(O/'provenance/audit_diagnostics_complete.json').write_text(json.dumps({'exact_rank_correlation_rows':len(extreme),'SERPINA1_original_reproduced':True,'no_new_endpoint':True},indent=2))
