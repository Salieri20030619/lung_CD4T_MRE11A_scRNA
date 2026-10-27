refs=pd.read_csv(T/'network_seed17_reference_associations.tsv',sep='\t')
for fam,z in refs.groupby('family'):
 z=z[z.P.notna()]
 if len(z):ok('FullSeed17 BH '+fam,np.allclose(bh(z.P),z.FDR,atol=1e-10))
ra=pd.read_csv(T/'network_seed17_reference_AUC.tsv',sep='\t');rs=pd.read_csv(P/'analysis/tables/sample_scores.tsv',sep='\t')
refnums=[]
for r in refs[refs.gene.eq('IL2RA')&refs.tf.eq('NFKB2')&refs.status.eq('OK')].itertuples():
 aa=ra[(ra.arm==r.arm)&(ra.scope==r.scope)&ra.tf.eq('NFKB2')]
 z=rs[(rs.variant==r.variant)&(rs.scope==r.scope)&rs.dataset.ne('GSE274372')].merge(aa[['dataset','sample','AUC']],on=['dataset','sample']).sort_values(['dataset','sample']).reset_index(drop=True)
 labels=z.dataset.astype(str)+'_'+z.day.astype(str);C=pd.get_dummies(labels,dtype=float).to_numpy()
 if 'Condition' in r.model:C=np.column_stack([C,z.group.eq('Silica').to_numpy(float)])
 M=np.eye(len(z))-C@np.linalg.pinv(C);rx=M@stats.rankdata(z.AUC);y=stats.rankdata(z.G_IL2RA);ry=M@y;rho=rx@ry/(np.linalg.norm(rx)*np.linalg.norm(ry))
 if 'Condition' in r.model:labels=labels+'_'+z.group
 blocks=[np.flatnonzero(labels.to_numpy()==a) for a in sorted(labels.unique())];ge=nn=0
 for part in itertools.product(*[list(itertools.permutations(b)) for b in blocks]):
  ix=np.arange(len(z))
  for b,pp in zip(blocks,part):ix[b]=pp
  yp=M@y[ix];rp=rx@yp/(np.linalg.norm(rx)*np.linalg.norm(yp));ge+=abs(rp)>=abs(rho)-1e-12;nn+=1
 pv=ge/nn
 ok('FullSeed17 n10 '+r.arm+' '+r.model,len(z)==10 and z['sample'].is_unique)
 ok('FullSeed17 rho '+r.arm+' '+r.model,abs(rho-r.rho)<1e-10)
 ok('FullSeed17 exact P '+r.arm+' '+r.model,abs(pv-r.P)<1e-10)
 refnums.append(dict(arm=r.arm,model=r.model,n=10,rho=rho,P=pv,permutations=nn))
pd.DataFrame(refnums).to_csv(T/'independent_seed17_reference_QA.tsv',sep='\t',index=False,encoding='utf-8')
