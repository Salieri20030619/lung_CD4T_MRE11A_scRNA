from pathlib import Path
import os,json,hashlib
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
import pandas as pd,numpy as np,scipy.sparse as sp
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');B=P/'11_reports/three_cohort_four_axis_20261017_v16';O=P/'11_reports/append_analysis_20261017_v18';plan=json.loads((O/'protocol_frozen.json').read_text());cand=sorted(set(g for g,s in plan['focus_pairs']));mods=pd.read_csv(O/'inputs/fixed_modules.tsv',sep='\t');gly=pd.read_csv(O/'inputs/glycolysis_GO_members.tsv',sep='\t');mods=pd.concat([mods,gly[['module','gene']]]).drop_duplicates();mods['removed_candidate']=mods.gene.isin(cand);mods.to_csv(O/'tables/module_member_audit.tsv',sep='\t',index=False);sets={m:sorted(set(z.gene)-set(cand)) for m,z in mods.groupby('module')};ligs=sorted(set(z[2] for z in plan['routes']));lr=pd.read_csv(B/'tables/expressed_LR_with_decoy_literature.tsv',sep='\t');selected=lr[lr.ligand_parts.isin(ligs)];selected.to_csv(O/'tables/frozen_route_LR_coverage.tsv',sep='\t',index=False);rec=set(';'.join(selected.receptor_parts).split(';'));need=sorted(set(cand)|set(mods.gene)|set(ligs)|rec|{'ENO1'});rows=[];expression=[];coverage=[];manifest=[]
for ds,core in [('GSE314943','Core'),('GSE183682','Core'),('GSE274372','TLike')]:
 stem=ds+'__'+core;files=[B/'inputs'/f'{stem}{s}' for s in ['.cellcounts.npz','.cells.tsv.gz','.genes.txt']]
 for f in files:manifest.append(dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
 x=sp.load_npz(files[0]);m=pd.read_csv(files[1],sep='\t');g=pd.read_csv(files[2],header=None)[0].tolist();ix={v:i for i,v in enumerate(g)}
 for mod,genes in sets.items():
  mapped=[v for v in genes if v in ix];coverage.append(dict(dataset=ds,module=mod,defined_genes=len(genes),mapped_genes=len(mapped),mapped_fraction=len(mapped)/len(genes),mapped=';'.join(mapped),status='SCORABLE' if len(mapped)>=3 and len(mapped)>=len(genes)*.5 else 'INSUFFICIENT_COVERAGE'))
 for variant in ['Raw','BG1']:
  for scope in ([core,'Naive','ActCD4','Treg'] if core=='Core' else [core]):
   take=np.ones(len(m),bool) if scope==core else m.subtype.eq(scope).to_numpy()
   if variant=='BG1':take &= m.background_fraction.to_numpy()<=.01
   for sample,mm in m.loc[take].groupby('sample'):
    ids=mm.index.to_numpy();counts=np.asarray(x[ids].sum(0)).ravel();library=mm.library.sum();cp=counts/library*1e4;log=np.log1p(cp);r=dict(dataset=ds,variant=variant,scope=scope,sample=sample,group=mm.group.iloc[0],day=int(mm.day.iloc[0]),cells=len(ids),library=library)
    for gene in cand+['ENO1']:r['G_'+gene]=log[ix[gene]] if gene in ix else np.nan;r['CP_'+gene]=cp[ix[gene]] if gene in ix else np.nan
    for mod,genes in sets.items():
     mapped=[v for v in genes if v in ix];r['M_'+mod]=float(log[[ix[v] for v in mapped]].mean()) if len(mapped)>=3 and len(mapped)>=len(genes)*.5 else np.nan
    rows.append(r)
    for gene in need:
     if gene in ix:v=x[ids,ix[gene]].toarray().ravel();umi=v.sum();positive=(v>0).sum();multi=(v>=2).sum()
     else:umi=positive=multi=np.nan
     expression.append(dict(dataset=ds,variant=variant,scope=scope,sample=sample,group=r['group'],gene=gene,cells=len(ids),UMI=umi,positive_cells=positive,multiUMI_cells=multi,CP10k=umi/library*1e4 if gene in ix else np.nan))
df=pd.DataFrame(rows);fixed=[]
for (ds,var,sample),z in df[df.scope.isin(['Naive','ActCD4','Treg'])].groupby(['dataset','variant','sample']):
 if set(z.scope)!=set(['Naive','ActCD4','Treg']):continue
 r=z.iloc[0].to_dict();r['scope']='Core';r['variant']='FixedComp' if var=='Raw' else 'BG1FixedComp';r['cells']=int(z.cells.sum());r['library']=float(z.library.sum())
 for c in df.columns:
  if c.startswith('M_') or c.startswith('CP_'):r[c]=float(z[c].mean())
  elif c.startswith('G_'):r[c]=float(np.log1p(z['CP_'+c[2:]].mean()))
 fixed.append(r)
df=pd.concat([df,pd.DataFrame(fixed)],ignore_index=True);df.to_csv(O/'tables/sample_gene_and_program_scores.tsv',sep='\t',index=False);pd.DataFrame(expression).to_csv(O/'tables/sample_selected_gene_expression.tsv.gz',sep='\t',index=False);pd.DataFrame(coverage).to_csv(O/'tables/module_coverage.tsv',sep='\t',index=False)
sender=[]
for scope in ['AllMac','AMlike']:
 stem='GSE314943__'+scope;c=pd.read_csv(B/'inputs'/f'{stem}.counts.tsv.gz',sep='\t',index_col=0);f=pd.read_csv(B/'inputs'/f'{stem}.fraction.tsv.gz',sep='\t',index_col=0);m=pd.read_csv(B/'inputs'/f'{stem}.samples.tsv',sep='\t')
 for r in m.itertuples():
  for lig in ligs:sender.append(dict(sender=scope,sample=r.sample,group=r.group,cells=r.cells,ligand=lig,CP10k=float(c.at[lig,r.sample]/r.library*1e4) if lig in c.index else np.nan,fraction=float(f.at[lig,r.sample]) if lig in f.index else np.nan))
pd.DataFrame(sender).to_csv(O/'tables/sender_ligand_expression.tsv',sep='\t',index=False)
whole=pd.read_csv(B/'tables/whole_lung_source_expression.tsv.gz',sep='\t');whole[whole.gene.isin(ligs)].to_csv(O/'tables/whole_lung_competing_sources.tsv',sep='\t',index=False)
# Keep the original DE family; subset only for convenient display, not recomputing FDR.
de=pd.read_csv(B/'tables/DE_all_models.tsv',sep='\t',low_memory=False);de[de.gene.isin(cand+ligs)].to_csv(O/'tables/original_DE_context.tsv',sep='\t',index=False)
(O/'inputs/input_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8');(O/'logs/export_complete.txt').write_text('success',encoding='utf-8');print('EXPORTED',df.shape,'expressions',len(expression),'module coverage',len(coverage))