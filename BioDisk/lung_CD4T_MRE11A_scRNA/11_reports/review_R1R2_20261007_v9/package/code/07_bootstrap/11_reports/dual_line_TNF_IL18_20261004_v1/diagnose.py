"""Frozen finite descriptive diagnostics. No old model/P/FDR/LOO refits."""
from pathlib import Path
import pandas as pd, numpy as np
from scipy.stats import rankdata
import hashlib,json,itertools,shutil

B=Path(__file__).resolve().parent; I=B/'inputs'; O=B/'results'; O.mkdir(exist_ok=True)
OLD=I/'previous_results'
if not OLD.exists():
 OLD.mkdir()
 source=B.parent/'recovery_delivery_20261004_v1/breakpoint_TNF_IL18_20261004_v1/results'
 for name in ['aligned_chain_models.tsv','available_per_donor_influence.tsv.gz','input_manifest.tsv','influence_summary.tsv','chain_detection_by_donor.tsv','donor_roster_reasons.tsv']:
  shutil.copyfile(source/name,OLD/name)
def rd(p):return pd.read_csv(p,sep='\t')
def wr(d,n):pd.DataFrame(d).to_csv(O/n,sep='\t',index=False,na_rep='NA')
def rho(a,b):
 a=np.asarray(a,float);b=np.asarray(b,float)
 if np.std(a)<1e-12 or np.std(b)<1e-12:return np.nan
 return float(np.corrcoef(a,b)[0,1])
def spr(a,b):return rho(rankdata(a),rankdata(b))
def residual(a,D):return a-D@np.linalg.lstsq(D,a,rcond=None)[0]
def design(z,extra=None):
 cols=[np.ones(len(z)),z.group.eq('CWP').to_numpy(float)];names=['Intercept','G']
 for s in sorted(z.series.unique())[1:]:cols.append(z.series.eq(s).to_numpy(float));names.append('S_'+s)
 for c in ['receiver_library','sender_library']+(extra or []):
  a=rankdata(z[c]);cols.append((a-a.mean())/a.std());names.append(c)
 return np.column_stack(cols),names
def zscore(a):return (a-a.mean())/a.std()
def selection(chain,modality):
 q=m[(m.chain==chain)&(m.modality==modality)&(m.policy=='aligned_steps')]
 ids=q.iloc[0].sample_ids.split(';');sender,scope,_=chain.split('/')
 return pair[(pair.sender==sender)&(pair.scope==scope)].set_index('sample').loc[ids].reset_index()
def oldrow(chain,modality,edge,policy):
 q=m[(m.chain==chain)&(m.modality==modality)&(m.edge==edge)&(m.policy==policy)]
 if edge!='ligand_TF':q=q[q.y=='State_Effector']
 if policy=='aligned_steps':q=q[q.covariates==base]
 assert len(q)==1
 return q.iloc[0]
pair=rd(I/'paired_inputs.tsv');m=rd(OLD/'aligned_chain_models.tsv');loo=rd(OLD/'available_per_donor_influence.tsv.gz')
base='G+S+receiver_library+sender_library'
manifest=[]
for path in list(I.iterdir())+list(OLD.iterdir()):
 if path.is_file():manifest.append({'path':str(path.resolve()),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size})
oldmanifest=rd(OLD/'input_manifest.tsv');expected=oldmanifest[oldmanifest.path.str.endswith('front_half_20261004_v1/results/paired_inputs.tsv')].iloc[0].sha256
assert hashlib.sha256((I/'paired_inputs.tsv').read_bytes()).hexdigest()==expected
tnf='MacAM/Core/TNF-STAT1';il='MacAM/CD4-IL7R/IL18-IKZF2'
decomp=[];vifs=[];explained=[];resrows=[];loomatched=[];sources=[];corrs=[]
for mod in ['RNA','AUC17']:
 z=selection(tnf,mod); D,names=design(z);T=rankdata(z.Treg_fraction);t=residual(T,D);tf='STAT1_RNA' if mod=='RNA' else 'STAT1_AUC_17'
 for title,ix in [('G_Source',list(range(4))),('depth_only',[0,4,5]),('G_Source_receiver_depth',list(range(5))),('base',list(range(6)))]:
  rr=residual(T,D[:,ix]);explained.append({'modality':mod,'n':len(z),'design':title,'Treg_rank_R2':1-np.sum(rr**2)/np.sum((T-T.mean())**2)})
 for withT in [False,True]:
  DD,nn=design(z,['Treg_fraction'] if withT else [])
  for j in range(1,len(nn)):
   r=residual(DD[:,j],np.delete(DD,j,axis=1));vifs.append({'modality':mod,'adjustment':'base+Treg' if withT else 'base','covariate':nn[j],'VIF':np.sum((DD[:,j]-DD[:,j].mean())**2)/np.sum(r*r)})
 vars=['L_TNF',tf,'State_Effector','Treg_fraction','receiver_library','sender_library']
 for a,b in itertools.combinations(vars,2):corrs.append({'axis':'TNF','modality':mod,'x':a,'y':b,'raw_spearman':spr(z[a],z[b])})
 for edge,x,y in [('ligand_TF','L_TNF',tf),('TF_phenotype',tf,'State_Effector'),('ligand_phenotype','L_TNF','State_Effector')]:
  a=residual(rankdata(z[x]),D);b=residual(rankdata(z[y]),D)
  rxy=rho(a,b);rxt=rho(a,t);ryt=rho(b,t);product=rxt*ryt;den=np.sqrt((1-rxt*rxt)*(1-ryt*ryt));adjusted=(rxy-product)/den
  oldb=oldrow(tnf,mod,edge,'aligned_steps');olda=oldrow(tnf,mod,edge,'Treg_composition')
  assert np.isclose(rxy,oldb.rho,atol=1e-12) and np.isclose(adjusted,olda.rho,atol=1e-12)
  decomp.append({'modality':mod,'n':len(z),'edge':edge,'base_spec':oldb.spec_id,'Treg_spec':olda.spec_id,'base_rho':oldb.rho,'Treg_rho':olda.rho,'delta_rho':olda.rho-oldb.rho,'r_xT_given_base':rxt,'r_yT_given_base':ryt,'removed_covariance_product':product,'normalization_denominator':den,'numerator_after_remove_T':rxy-product,'partial_formula_rho':adjusted,'standardized_beta_x_given_T':(rxy-product)/(1-rxt*rxt),'negative_covariance_suppression':product<0})
  for j,ss in enumerate(z['sample']):resrows.append({'modality':mod,'edge':edge,'sample':ss,'group':z.group.iloc[j],'series':z.series.iloc[j],'x_rank_resid':a[j],'y_rank_resid':b[j],'Treg_rank_resid':t[j],'x_rank_resid_after_T':a[j]-np.dot(a,t)/np.dot(t,t)*t[j],'y_rank_resid_after_T':b[j]-np.dot(b,t)/np.dot(t,t)*t[j]})
  ll=loo[loo.signature.eq(oldb.signature)][['omitted_sample','rho']].drop_duplicates('omitted_sample').rename(columns={'rho':'base_LOO_rho'}).merge(loo[loo.signature.eq(olda.signature)][['omitted_sample','rho']].drop_duplicates('omitted_sample').rename(columns={'rho':'Treg_LOO_rho'}),on='omitted_sample',validate='one_to_one')
  assert len(ll)==len(z)
  ll['delta_rho']=ll.Treg_LOO_rho-ll.base_LOO_rho;ll['modality']=mod;ll['edge']=edge;ll['origin']='existing saved omissions only';loomatched.extend(ll.to_dict('records'))
 for s,g in z.groupby('series'):
  sources.append({'modality':mod,'source':s,'n':len(g),'cases':int(g.group.eq('CWP').sum()),'Treg_min':g.Treg_fraction.min(),'Treg_median':g.Treg_fraction.median(),'Treg_max':g.Treg_fraction.max()})
wr(decomp,'TNF_Treg_covariance_decomposition.tsv');wr(vifs,'TNF_covariate_VIF.tsv');wr(explained,'TNF_Treg_explained_variance.tsv');wr(resrows,'TNF_rank_residuals_by_donor.tsv');wr(loomatched,'TNF_existing_LOO_adjustment_changes.tsv');wr(sources,'TNF_Treg_source_distribution.tsv')
# Frozen Treg numerator/denominator: no new annotation or cell identity selection.
counts=rd(I/'receiver_subtype_counts.tsv');fr=[]
for ss,q in counts.groupby('sample'):
 total=q.cells.sum();tr=q.loc[q.subtype.str.contains('Treg',case=False),'cells'].sum();p=pair[(pair.sender=='MacAM')&(pair.scope=='Core')&(pair['sample']==ss)].iloc[0]
 assert total==p.receiver_cells and np.isclose(tr/total,p.Treg_fraction)
 fr.append({'sample':ss,'Treg_cells':tr,'Core_cells':total,'Treg_fraction':tr/total,'State_Effector':p.State_Effector,'group':p.group,'series':p.series})
wr(fr,'TNF_Treg_numerator_denominator_check.tsv')
# IL18: frozen targets only; unweighted AUCell means inferred edge weights are descriptive metadata.
def canon(g):return {'ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','MRE11A':'MRE11'}.get(str(g).upper(),str(g).upper())
edges=rd(I/'regulon_edges.tsv');edges=edges[edges.tf.str.upper().eq('IKZF2')].copy();edges['target']=edges.target.map(canon)
excluded=set(rd(I/'module_definition.tsv').gene.map(canon))|set('IL1R2 PENK TTC39C FOXP1 TNFRSF9 LAG3 MRE11 IKZF2'.split())
edges['retained_DropOverlap']=~edges.target.isin(excluded);wr(edges,'IKZF2_frozen_targets_and_weights.tsv')
coverage=rd(I/'regulon_exclusion_coverage.tsv');sizes=[];over=[]
for seed in [17,29,43]:
 q=edges[edges.seed==seed];actual=coverage[(coverage.species=='Human')&(coverage['filter']=='All')&(coverage.seed==seed)&(coverage.variant=='DropOverlap')&(coverage.regulon=='IKZF2(+)')].iloc[0]
 assert q.target.nunique()==actual.original_genes and q[q.retained_DropOverlap].target.nunique()==actual.retained_genes
 sizes.append({'seed':seed,'original_genes':int(actual.original_genes),'DropOverlap_genes':int(actual.retained_genes),'removed':';'.join(sorted(set(q.loc[~q.retained_DropOverlap,'target']))),'weighted_AUC':False})
for a,b in itertools.combinations([17,29,43],2):
 for variant in ['original','DropOverlap']:
  q=edges if variant=='original' else edges[edges.retained_DropOverlap]
  ga=set(q.loc[q.seed==a,'target']);gb=set(q.loc[q.seed==b,'target']);both=sorted(ga&gb)
  w1=q[q.seed==a].set_index('target').weight.reindex(both);w2=q[q.seed==b].set_index('target').weight.reindex(both)
  over.append({'seed1':a,'seed2':b,'variant':variant,'n1':len(ga),'n2':len(gb),'intersection':len(both),'union':len(ga|gb),'Jaccard':len(both)/len(ga|gb),'only_seed1':';'.join(sorted(ga-gb)),'only_seed2':';'.join(sorted(gb-ga)),'shared_weight_spearman':spr(w1,w2) if len(both)>2 else np.nan,'weights_used_for_AUC':False})
wr(sizes,'IKZF2_target_sizes.tsv');wr(over,'IKZF2_target_overlap.tsv')
z=selection(il,'AUC17');D,names=design(z);aucs={s:z['IKZF2_AUC_'+str(s)].to_numpy() for s in [17,29,43]};ranks={s:rankdata(v) for s,v in aucs.items()};res={s:residual(v,D) for s,v in ranks.items()};seedcorr=[];diffs=[];contexts=[]
for a,b in itertools.combinations([17,29,43],2):
 seedcorr.append({'seed1':a,'seed2':b,'n':len(z),'AUC_pearson':rho(aucs[a],aucs[b]),'AUC_spearman':rho(ranks[a],ranks[b]),'rank_residual_correlation':rho(res[a],res[b]),'max_absolute_rank_change':np.max(np.abs(ranks[b]-ranks[a])),'median_absolute_rank_change':np.median(np.abs(ranks[b]-ranks[a]))})
 delta=ranks[b]-ranks[a];dres=zscore(res[b])-zscore(res[a]);L=zscore(residual(rankdata(z.L_IL18),D));contribution=L*dres/len(z)
 for j,ss in enumerate(z['sample']):diffs.append({'seed1':a,'seed2':b,'sample':ss,'series':z.series.iloc[j],'group':z.group.iloc[j],'auc1':aucs[a][j],'auc2':aucs[b][j],'rank1':ranks[a][j],'rank2':ranks[b][j],'rank_change':delta[j],'standardized_residual_difference':dres[j],'additive_contribution_to_ligand_rho_change':contribution[j],'not_a_LOO_or_causal_effect':True})
 assert np.isclose(contribution.sum(),rho(res[b],L)-rho(res[a],L))
 for c in ['L_IL18','IKZF2_RNA','State_Effector','State_Activation','receiver_library','sender_library','receiver_cells','sender_cells','R_IL18']:
  contexts.append({'seed1':a,'seed2':b,'variable':c,'rank_change_spearman':spr(delta,z[c]),'n':len(z),'inference':'descriptive; no P/FDR'})
for seed in [17,29,43]:
 for c in ['L_IL18','IKZF2_RNA','State_Effector','State_Activation','receiver_library','sender_library','receiver_cells','R_IL18']:
  row={'axis':'IL18','modality':'AUC_seed'+str(seed),'x':'IKZF2_AUC_'+str(seed),'y':c,'raw_spearman':spr(aucs[seed],z[c])}
  if c not in ['receiver_library','sender_library']:row['rank_residual_correlation']=rho(res[seed],residual(rankdata(z[c]),D))
  corrs.append(row)
wr(seedcorr,'IKZF2_seed_activity_agreement.tsv');wr(diffs,'IKZF2_seed_donor_rank_changes.tsv');wr(contexts,'IKZF2_seed_change_context.tsv');wr(corrs,'fixed_roster_descriptive_correlations.tsv')
# Directly copy relevant inferential results; no statistical recalculation.
f=m[m.chain.isin([tnf,il])].copy();wr(f,'existing_two_axis_models_REUSED.tsv')
infl=rd(OLD/'influence_summary.tsv');wr(infl[infl.spec_id.isin(f.spec_id)],'existing_two_axis_influence_REUSED.tsv')
det=rd(OLD/'chain_detection_by_donor.tsv');det=det[(det.chain==il)&(det.roster=='AUC17')].merge(z[['sample','receiver_library','sender_library']],on='sample',validate='one_to_one')
det['receptor_detected']=det.receptor_pct>0;wr(det,'IL18_receptor_depth_by_donor.tsv')
dr=[]
for present,q in det.groupby('receptor_detected'):
 dr.append({'receptor_detected':bool(present),'n':len(q),'receiver_cells_median':q.receiver_cells.median(),'receiver_cells_min':q.receiver_cells.min(),'receiver_cells_max':q.receiver_cells.max(),'receiver_library_median':q.receiver_library.median(),'sender_library_median':q.sender_library.median(),'receptor_positive_cells':q.receptor_positive_cells.sum(),'sample_ids':';'.join(q['sample'])})
wr(dr,'IL18_receptor_depth_summary.tsv')
roster=rd(OLD/'donor_roster_reasons.tsv');wr(roster[roster.chain.isin([tnf,il])],'existing_roster_reasons_REUSED.tsv')
wr(manifest,'input_manifest.tsv')
checks={'status':'PASS','plan_frozen_before_diagnostics':True,'old_model_or_LOO_runs':0,'new_P_or_FDR':0,'TNF_decomposition_rows':len(decomp),'old_omission_rows_reused':len(loomatched),'IKZF2_target_counts_match_original_AUC_audit':True,'Treg_fraction_matches_all_19_numerators_denominators':True,'old_results_unchanged':True,'limitations':['No new same-identity cell-level effect estimates: existing summaries do not separate Treg/non-Treg per-gene/state means.','Seed changes combine frozen target-set and AUCell ranking random-seed effects; no crossed rerun attempted.']}
for entry in manifest:
 assert hashlib.sha256(Path(entry['path']).read_bytes()).hexdigest()==entry['sha256']
(B/'verification.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
print(json.dumps(checks,indent=2));print(pd.DataFrame(decomp).to_string(index=False));print(pd.DataFrame(sizes).to_string(index=False));print(pd.DataFrame(seedcorr).to_string(index=False))
