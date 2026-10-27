# Independent focal numerical verification of new leave-one-mouse-out results.
auc=pd.read_csv(T/'network_LOO_sample_AUC.tsv',sep='\t')
sc=pd.read_csv(P/'analysis/tables/sample_scores.tsv',sep='\t')
tests=pd.read_csv(T/'network_LOO_candidate_associations.tsv',sep='\t')
tests=tests[(tests.gene=='IL2RA')&(tests.tf=='NFKB2')&(tests.status=='OK')]
verified=[]
for r in tests.itertuples():
 aa=auc[(auc.arm==r.arm)&(auc.left_out==r.left_out)&(auc.scope==r.scope)&(auc.tf==r.tf)]
 z=sc[(sc.scope==r.scope)&(sc.variant==r.variant)&sc.dataset.ne('GSE274372')&sc['sample'].ne(r.left_out)].merge(aa[['sample','dataset','AUC']],on=['sample','dataset']).sort_values(['dataset','sample']).reset_index(drop=True)
 ok('LOO independent paired sample count '+r.arm+' '+r.left_out+' '+r.model,len(z)==9 and z['sample'].is_unique)
 labels=z.dataset.astype(str)+'_'+z.day.astype(str)
 C=pd.get_dummies(labels,drop_first=False,dtype=float).to_numpy()
 if 'Condition' in r.model:C=np.column_stack([C,z.group.eq('Silica').to_numpy(float)])
 M=np.eye(len(z))-C@np.linalg.pinv(C);x=stats.rankdata(z.AUC);y=stats.rankdata(z.G_IL2RA);rx=M@x;ry=M@y;rho=rx@ry/(np.linalg.norm(rx)*np.linalg.norm(ry))
 if 'Condition' in r.model:labels=labels+'_'+z.group
 blocks=[np.flatnonzero(labels.to_numpy()==v) for v in sorted(labels.unique())]
 ge=n=0
 for part in itertools.product(*[list(itertools.permutations(b)) for b in blocks]):
  ix=np.arange(len(z))
  for b,k in zip(blocks,part):ix[b]=k
  yp=M@y[ix];den=np.linalg.norm(rx)*np.linalg.norm(yp)
  if den<=1e-12:continue
  rr=rx@yp/den;ge+=abs(rr)>=abs(rho)-1e-12;n+=1
 pv=ge/n
 ok('LOO independent rho '+r.arm+' '+r.left_out+' '+r.model,abs(rho-r.rho)<1e-10)
 ok('LOO independent exact P '+r.arm+' '+r.left_out+' '+r.model,abs(pv-r.P)<1e-10)
 verified.append(dict(arm=r.arm,left_out=r.left_out,model=r.model,n=9,rho=rho,P=pv,permutations=n))
pd.DataFrame(verified).to_csv(T/'independent_network_LOO_QA.tsv',sep='\t',index=False,encoding='utf-8')
