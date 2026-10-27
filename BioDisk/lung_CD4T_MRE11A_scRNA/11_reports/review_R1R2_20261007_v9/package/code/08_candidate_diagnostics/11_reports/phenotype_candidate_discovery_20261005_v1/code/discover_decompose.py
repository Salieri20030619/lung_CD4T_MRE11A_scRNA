from pathlib import Path
import os,sys,json
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import pandas as pd,numpy as np
from scipy.io import mmread
from scipy import stats
from threadpoolctl import threadpool_limits
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/phenotype_candidate_discovery_20261005_v1';R=O/'results';H=P/'11_reports/human_first_20261004_v1/results';A=P/'02_analysis/mechanism_step2_20260927_154919_EsN4gk'
src=(O/'code/discover_compute.py').read_text();exec(src[src.index('def wr('):src.index('scores=rd(')])
rd=lambda p:pd.read_csv(p,sep='\t',low_memory=False)
sys.path.insert(0,str(A));from shared import canonical
ctrl=rd(A/'tables/human_score_controls.tsv.gz');members=rd(H/'stage1_complete_program_members.tsv');m=rd(A/'input/human_cells.tsv');g=rd(A/'input/human_genes.tsv').symbol.astype(str);canon=np.array([canonical(v) for v in g]);old=rd(A/'objects/human_cell_scores.tsv.gz');m=m.merge(old[['cell','bg','score_Mitophagy','score_TReg']],on='cell',validate='one_to_one')
with threadpool_limits(limits=2):x=mmread(A/'input/human_counts.mtx').tocsr()
out=[];memberrows=[];checks=[]
for prog,rec in [('Mitophagy','CD4-IL7R'),('TReg','CD4-GZMA')]:
 tg=members[members.program==prog].gene.tolist();cg=ctrl[ctrl.module==prog].gene.tolist();assert set(tg+cg)<=set(canon)
 logs={s:np.log1p(np.asarray(x[canon==s].sum(axis=0)).ravel()/m.library*1e4) for s in set(tg+cg)}
 target=np.mean([logs[s] for s in tg],axis=0);background=np.mean([logs[s] for s in cg],axis=0);err=float(np.max(abs(target-background-m['score_'+prog])));assert err<1e-10;checks.append(dict(program=prog,target_n=len(tg),control_n=len(cg),max_cell_reconstruction_error=err))
 for role,genes in [('target',tg),('original_fixed_control',cg)]:
  for gene in genes:memberrows.append(dict(program=prog,role=role,gene=gene))
 mm=m.assign(target_raw=target,background_raw=background,original_score=target-background)
 for filt in ['All','BG1']:
  zz=mm[(mm.subtype==rec)&((mm.bg<=.01) if filt=='BG1' else True)]
  for ss,z in zz.groupby('sample'):
   for component in ['target_raw','background_raw','original_score']:out.append(dict(receiver=rec,program=prog,filter=filt,sample=ss,group=z.group.iloc[0],series=z.series.iloc[0],cells=len(z),library=z.library.mean(),component=component,value=z[component].mean()))
wr(memberrows,'original_target_and_control_complete_members.tsv.gz');wr(out,'program_components_donor_SERVER_ONLY.tsv.gz');wr(checks,'original_program_reconstruction.tsv')
d=pd.DataFrame(out);dis=[];cors=[];inp=rd(R/'selected_fixed_donor_inputs_SERVER_ONLY.tsv.gz');inp=inp[inp.followup_roster=='All_eligible']
for (rec,prog,filt,component),z in d.groupby(['receiver','program','filter','component']):
 z=z[z.cells>=20]
 for model in ['Base','Source','SourceDep']:dis.append(dict(receiver=rec,program=prog,filter=filt,component=component,model=model,**ols(z,'value',model)))
 if filt!='All':continue
 for lig,v in inp[inp.receiver==rec].groupby('ligand'):
  q=v.drop(columns='y').merge(z[['sample','value']],on='sample').rename(columns={'value':'y'})
  for step in ['Unadjusted','Group','Source','SourceBothDepth']:cors.append(dict(receiver=rec,program=prog,ligand=lig,component=component,step=step,**corr(q,step)))
dis=bh(pd.DataFrame(dis),['filter','model']);cs=bh(pd.DataFrame(cors),['step']);wr(dis,'program_component_disease_models.tsv');wr(cs,'program_component_selected_associations.tsv');(O/'provenance/decomposition_complete.json').write_text(json.dumps({'cell_score_reproduction':checks,'disease_tests':len(dis),'association_tests':len(cs),'interpretation':'fixed original target vs control decomposition; no new score selection'},indent=2));print('COMPONENT_DISEASE',dis[(dis['filter']=='All')&(dis.model=='Source')][['receiver','component','n','effect','P','FDR']].to_string(index=False));print('COMPONENT_ASSOC',cs[cs.step=='Source'][['ligand','component','n','rho','P','FDR']].to_string(index=False))
