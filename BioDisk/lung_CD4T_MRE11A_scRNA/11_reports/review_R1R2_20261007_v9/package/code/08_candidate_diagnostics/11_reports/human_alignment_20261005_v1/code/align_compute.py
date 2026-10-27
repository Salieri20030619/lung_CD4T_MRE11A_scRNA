from pathlib import Path
import os,json,hashlib
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import pandas as pd,numpy as np
from scipy import stats
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/human_alignment_20261005_v1';R=O/'results';H=P/'11_reports/human_first_20261004_v1/results';T=P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release';A=P/'02_analysis'
L='IL15 TNF CLSTN1 ADAM15 PLAU DSC2 CXCL2 ICAM2 VTN VEGFB GSTP1 RTN4'.split();TF=['IKZF2','STAT1'];SC=['Core','CD4-CCR7','CD4-GZMA','CD4-IL7R','Treg'];K=['receiver','filter','model','test','direction']
def rd(f):return pd.read_csv(f,sep='\t',low_memory=False)
def wr(d,n):pd.DataFrame(d).to_csv(R/n,sep='\t',index=False,na_rep='NA')
def bh(d,keys):
 d['FDR']=np.nan;d['BH_family']=''
 for key,ix in d.groupby(keys,dropna=False).groups.items():
  ix=[i for i in ix if np.isfinite(d.loc[i,'P'])];d.loc[ix,'FDR']=stats.false_discovery_control(d.loc[ix,'P']) if ix else [];d.loc[ix,'BH_family']=str(key)+';m='+str(len(ix))
 return d
def corr(d,x='x',y='y',two=False):
 depths=['receiver_depth']+(['sender_depth'] if two else []);need=[x,y]+depths
 q=d[np.isfinite(d[need]).all(axis=1)&(d.receiver_cells>=20)].copy()
 if two:q=q[q.sender_cells>=20]
 n=len(q);out=dict(n=n,n_case=int(q.group.eq('CWP').sum()),n_control=int(q.group.eq('Control').sum()),samples=';'.join(q['sample']),df=np.nan,rho=np.nan,P=np.nan,status='INSUFFICIENT')
 if min(out['n_case'],out['n_control'])<2:return out
 X=np.column_stack([np.ones(n),q.group.eq('CWP'),pd.get_dummies(q.series,drop_first=True,dtype=float)]+[stats.rankdata(q[c]) for c in (['sender_depth','receiver_depth'] if two else depths)])
 rank=np.linalg.matrix_rank(X);df=n-rank-1;out['df']=df
 if rank<X.shape[1] or df<3:out['status']='NON_IDENTIFIABLE_OR_DF_LT3';return out
 a=stats.rankdata(q[x]);b=stats.rankdata(q[y]);a-=X@np.linalg.lstsq(X,a,rcond=None)[0];b-=X@np.linalg.lstsq(X,b,rcond=None)[0]
 if min(np.std(a),np.std(b))<1e-10:out['status']='NO_RESIDUAL_VARIATION';return out
 rho=np.clip(np.corrcoef(a,b)[0,1],-1,1);p=2*stats.t.sf(abs(rho)*np.sqrt(df/max(1e-15,1-rho*rho)),df);out.update(rho=rho,P=p,status='EXPLORATORY_SAME_DATA');return out
def ols(d,model):
 q=d[np.isfinite(d.value)&(d.cells>=20)].copy();n=len(q);out=dict(n=n,n_case=int(q.group.eq('CWP').sum()),n_control=int(q.group.eq('Control').sum()),samples=';'.join(q['sample']),df=np.nan,effect=np.nan,P=np.nan,CI_low=np.nan,CI_high=np.nan,status='INSUFFICIENT')
 if min(out['n_case'],out['n_control'])<2:return out
 cols=[np.ones(n),q.group.eq('CWP')]
 if model!='Base':cols.append(pd.get_dummies(q.series,drop_first=True,dtype=float))
 if model=='SourceDep':cols.append(np.log10(q.mean_depth))
 X=np.column_stack(cols);rank=np.linalg.matrix_rank(X);df=n-rank;out['df']=df
 if rank<X.shape[1] or df<2:out['status']='NON_IDENTIFIABLE';return out
 b=np.linalg.lstsq(X,q.value,rcond=None)[0];res=q.value-X@b;se=np.sqrt(np.sum(res**2)/df*np.linalg.inv(X.T@X)[1,1]);p=2*stats.t.sf(abs(b[1]/se),df) if se>0 else np.nan;ci=stats.t.ppf(.975,df)*se
 out.update(effect=b[1],P=p,CI_low=b[1]-ci,CI_high=b[1]+ci,status='FIT');return out
assert 'EXPORT_COMPLETE' in (O/'provenance/export.log').read_text()
spec=O/'code/ALIGNMENT_SPEC.txt';(O/'provenance/ALIGNMENT_SPEC.txt').write_bytes(spec.read_bytes())
design=rd(H/'stage1_donor_design.tsv');meta=rd(H/'stage2_scope_sample_counts.tsv');genes=rd(R/'alignment_gene_donors.tsv.gz');progs=rd(R/'TF_excluded_program_donors.tsv.gz')
full=rd(H/'stage2_evidence_all_candidates.tsv.gz');cand=full[full.ligand.isin(L)].copy();wr(cand,'all_fixed_candidates_full_context.tsv.gz')
de=rd(H/'stage2_sender_disease_DE.tsv.gz');wr(de[de.gene.isin(L+['IL15RA'])],'sender_disease_fixed_candidates.tsv')
lr=rd(H/'stage2_all_LR_supported_and_unsupported.tsv.gz');wr(lr[lr['from'].isin(L)],'all_fixed_candidate_LR_pairs.tsv.gz')
wr(rd(H/'stage2_receptor_complex_audit.tsv.gz').query("ligand == 'IL15'"),'IL15_receiver_complex.tsv.gz')
for n in ['stage2_targetset_members.tsv.gz','stage2_background_members.tsv.gz','stage2_targetset_status.tsv','stage1_all_finite_models.tsv','stage1_complete_program_members.tsv','stage1_donor_design.tsv','stage2_scope_sample_counts.tsv']:(R/n).write_bytes((H/n).read_bytes())
lt=[]
for f in H.glob('stage2_LT_*.tsv.gz'):
 d=rd(f);lt.append(d[d.ligand.isin(L)])
lt=pd.concat(lt,ignore_index=True);wr(lt,'fixed_ligand_target_complete_members.tsv.gz')
net=rd(T/'data/TF_network_all_edges.tsv');net=net[(net.species=='Human')&net.tf.isin(TF)].copy();wr(net,'existing_TF_edges_complete.tsv.gz')
program=rd(H/'stage1_complete_program_members.tsv');overlap=set(program.gene);old7=set('IL1R2 PENK TTC39C FOXP1 TNFRSF9 LAG3 MRE11'.split())
members=rd(H/'stage2_targetset_members.tsv.gz');bg=rd(H/'stage2_background_members.tsv.gz');ov=[];complete=[]
for key,v in members.groupby(K):
 b=bg
 for k,val in zip(K[:3],key[:3]):b=b[b[k]==val]
 universe=set(b.gene);target=set(v.gene)&universe
 for tf in TF:
  for seed in [17,29,43]:
   e=net[(net['filter']==key[1])&(net.tf==tf)&(net.seed==seed)]
   for variant in ['Original','DropOverlap']:
    excluded={tf}|(old7|overlap if variant=='DropOverlap' else set());reg=set(e.target)-excluded;regbg=reg&universe;hit=target&regbg
    n=len(universe);p=stats.hypergeom.sf(len(hit)-1,n,len(regbg),len(target)) if n and len(target) else np.nan
    rec=dict(zip(K,key));rec.update(tf=tf,seed=seed,variant=variant,background_n=n,target_n=len(target),regulon_background_n=len(regbg),overlap_n=len(hit),overlap_genes=';'.join(sorted(hit)),tf_in_target=tf in target,P=p);ov.append(rec)
    for g in sorted(hit):complete.append(dict(rec,gene=g))
ov=bh(pd.DataFrame(ov),['filter','variant','seed']);wr(ov,'TF_disease_target_overlap.tsv');wr(complete,'TF_disease_overlap_complete_members.tsv.gz')
# Donor-level TF readouts; existing AUC roster/cells and depth remain explicit.
rna=genes[(genes.lineage=='CD4')&genes.gene.isin(TF)].copy();rna=rna.rename(columns={'gene':'tf','mean_log1pCP10k':'value'});rna['readout']='RNA';rna['variant']='RNA';rna['seed']=0
auc=rd(T/'data/S2_TF_sample_activity.tsv');auc=auc[(auc.species=='Human')&auc.scope.isin(SC)&auc.variant.isin(['Original','DropOverlap'])].copy();auc=auc.rename(columns={'mean_auc':'value'});auc['readout']='AUC'
auc=auc.drop(columns=['series'],errors='ignore').merge(design[['sample','series']],on='sample',validate='many_to_one')
ac=rd(A/'scenic1_20260929_103513_xRoEnr/Human/All/cells.tsv');depth=[]
for filt in ['All','BG1']:
 c=ac if filt=='All' else ac[ac.bg<=.01]
 for scope in SC:
  cc=c if scope=='Core' else c[c.subtype==scope]
  for s,g in cc.groupby('sample'):depth.append(dict(filter=filt,scope=scope,sample=s,mean_depth=g.library.mean()))
auc=auc.merge(pd.DataFrame(depth),on=['filter','scope','sample'],validate='many_to_one')
readouts=pd.concat([rna,auc],ignore_index=True);wr(readouts,'TF_readout_donors.tsv.gz')
audit=readouts.merge(design[['sample','group','series']],on='sample',suffixes=('_readout','_design'));assert audit.group_readout.eq(audit.group_design).all() and audit.series_readout.eq(audit.series_design).all()
dr=[];rk=['readout','variant','seed','filter','scope','tf']
for key,d in readouts.groupby(rk):
 for model in ['Base','Source','SourceDep']:dr.append(dict(zip(rk,key),model=model,**ols(d,model)))
dr=bh(pd.DataFrame(dr),['readout','variant','seed','filter','model']);wr(dr,'TF_disease_models.tsv')
# Fixed bridge contexts, TF-excluded program readout.
contexts=[(s,'Nominal05') for s in SC]+[('CD4-GZMA','FDR05')];rows=[];datasets={};idx=0
for scope,test in contexts:
 for tf in TF:
  pg=progs[(progs.receiver==scope)&(progs['filter']=='All')&(progs.model=='Source')&(progs.test==test)&(progs.direction=='Up')&(progs.excluded_tf==tf)][['sample','target_score','genes_n']]
  r=readouts[(readouts.scope==scope)&(readouts['filter']=='All')&(readouts.tf==tf)]
  for (readout,variant,seed),dd in r.groupby(['readout','variant','seed']):
   base=dd[['sample','group','series','cells','mean_depth','value']].rename(columns={'cells':'receiver_cells','mean_depth':'receiver_depth'})
   d=base.merge(pg,on='sample').rename(columns={'value':'x','target_score':'y'})
   tag=dict(receiver=scope,test=test,tf=tf,readout=readout,variant=variant,seed=seed,model='Source',filter='All',direction='Up')
   out=dict(tag,link='TF_to_program',ligand='NONE',sender='NONE',**corr(d));out['id']=idx;rows.append(out);datasets[idx]=(d,False);idx+=1
   for sender in ([] if test=='FDR05' else ['MacAM','MonoMac','MacWide']):
    for lig in L:
     sx=genes[(genes.key==sender)&(genes.gene==lig)][['sample','cells','mean_depth','mean_CP10k']].rename(columns={'cells':'sender_cells','mean_depth':'sender_depth','mean_CP10k':'x'})
     d=base.rename(columns={'value':'y'}).merge(sx,on='sample');out=dict(tag,link='ligand_to_TF',ligand=lig,sender=sender,**corr(d,two=True));out['id']=idx;rows.append(out);datasets[idx]=(d,True);idx+=1
rows=bh(pd.DataFrame(rows),['readout','variant','seed','link']);wr(rows,'fixed_TF_bridge_associations.tsv.gz')
# Frozen old BH leads, reproduce exactly and donor deletion sensitivity.
old=rd(H/'stage2_same_donor_associations.tsv.gz');old=old[old.FDR<.05].copy();expr=rd(H/'stage2_LR_expression_all_donors.tsv.gz');scores=rd(H/'stage2_disease_program_donor_scores.tsv.gz');sens=[];check=[]
def influence(d,two,tag):
 original=corr(d,two=two);ss=[]
 eligible=d[d.receiver_cells.ge(20)]
 if two:eligible=eligible[eligible.sender_cells.ge(20)]
 for s in eligible['sample'].unique():
  z=corr(d[d['sample']!=s],two=two);ss.append(dict(tag,omitted=s,**z,sign_flip=bool(np.isfinite(z['rho']) and np.isfinite(original['rho']) and np.sign(z['rho'])!=np.sign(original['rho']))))
 return ss
for oi,v in old.iterrows():
 pg=scores
 for k in K:pg=pg[pg[k]==v[k]]
 rm=meta[(meta.lineage=='CD4')&(meta.scope==v.receiver)&(meta['filter']==v['filter'])][['sample','cells','mean_depth']].rename(columns={'cells':'receiver_cells','mean_depth':'receiver_depth'})
 sm=meta[(meta.lineage=='Mac')&(meta.scope==v.sender)][['sample','cells','mean_depth']].rename(columns={'cells':'sender_cells','mean_depth':'sender_depth'})
 sx=expr[(expr.key==v.sender)&(expr.gene==v.ligand)][['sample','mean_CP10k']].rename(columns={'mean_CP10k':'x'})
 d=pg[['sample','target_score']].rename(columns={'target_score':'y'}).merge(design[['sample','group','series']],on='sample').merge(rm,on='sample').merge(sm,on='sample').merge(sx,on='sample')
 tag={k:v[k] for k in K+['ligand','sender']};tag.update(kind='previous_BH_lead',old_rho=v.rho,old_P=v.P,old_FDR=v.FDR)
 z=corr(d,two=True);check.append(dict(tag,**z));assert np.isclose(z['rho'],v.rho,atol=1e-10);assert np.isclose(z['P'],v.P,rtol=1e-8)
 sens+=influence(d,True,tag)
for _,v in rows[rows.FDR<.05].iterrows():
 d,two=datasets[v.id];tag={k:v[k] for k in ['receiver','test','tf','readout','variant','seed','link','ligand','sender','rho','P','FDR']};tag['kind']='new_BH_bridge';sens+=influence(d,two,tag)
wr(check,'previous_nine_associations_reproduced.tsv');wr(sens,'donor_deletion_influence.tsv.gz')
ss=pd.DataFrame(sens);su=[]
for key,d in ss[ss.kind=='previous_BH_lead'].groupby(K+['ligand','sender']):su.append(dict(zip(K+['ligand','sender'],key),deletions=len(d),computed=int(d.P.notna().sum()),nonidentifiable=int(d.P.isna().sum()),rho_min=d.rho.min(),rho_max=d.rho.max(),sign_flips=int(d.sign_flip.sum()),max_P=d.P.max(),most_changed_donor=d.loc[(d.rho-d.old_rho).abs().idxmax(),'omitted'] if d.rho.notna().any() else 'NA'))
wr(su,'previous_nine_influence_summary.tsv')
summary={'candidate_rows':len(cand),'ligand_target_rows':len(lt),'disease_TF_models':len(dr),'disease_TF_BH':int((dr.FDR<.05).sum()),'bridge_rows':len(rows),'bridge_computed':int(rows.P.notna().sum()),'bridge_BH':int((rows.FDR<.05).sum()),'old_nine_reproduced':len(check),'spec_sha256':hashlib.sha256(spec.read_bytes()).hexdigest()}
(O/'provenance/alignment_complete.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary));print('TF_DISEASE_BH',dr[(dr.FDR<.05)&(dr['filter']=='All')].to_string(index=False));print('BRIDGE_BH',rows[rows.FDR<.05].to_string(index=False));print('INFLUENCE',pd.DataFrame(su).to_string(index=False))
