from pathlib import Path
import pandas as pd,numpy as np,hashlib,json,gzip,shutil,zipfile
P=Path('technical_delivery')
original=P/'prior_review_original/scenic3_review'
repair=P/'prior_review_repaired/scenic3_review'
assert not repair.exists()
shutil.copytree(original,repair)
source=P/'source_tables/S3_human_refined_cells_original.tsv.gz'
expected='a5790136023a07c6529789620923d8268493106a59dedc4cba1f33af4f45d9aa'
assert hashlib.sha256(source.read_bytes()).hexdigest()==expected
old=original/'results/reviewed_cells.tsv.gz'
try:
 with gzip.open(old,'rb') as f: f.read()
 damaged=False
except (EOFError,OSError) as e:
 damaged=True
 old_error=str(e)
assert damaged
c=pd.read_csv(source,sep='\t')
ann=pd.read_csv(original/'results/reviewed_annotation.tsv',sep='\t')
assert c.cell.is_unique and len(c)==31120
assert set(c.r4)==set(ann.cluster)
ann=ann.set_index('cluster')
c['label']=c.r4.map(ann.label);c['review_scope']=c.r4.map(ann.scope)
assert c[['label','review_scope']].notna().all().all()
cols=['cell','sample','group','series','parent','r2','r4','r8','label','review_scope']
c[cols].to_csv(repair/'results/reviewed_cells.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
rr=pd.read_csv(repair/'results/reviewed_cells.tsv.gz',sep='\t')
pd.testing.assert_frame_equal(rr,c[cols],check_dtype=False)
with gzip.open(repair/'results/reviewed_cells.tsv.gz','rb') as f: f.read()
counts=c.groupby(['r4','label','sample','group','series'],observed=True).size().reset_index(name='cells')
before=pd.read_csv(original/'results/annotation_sample_counts.tsv',sep='\t')
keys=['r4','label','sample','group','series']
pd.testing.assert_frame_equal(counts.sort_values(keys).reset_index(drop=True),before.sort_values(keys).reset_index(drop=True),check_dtype=False)
unchanged=[]
for src in original.rglob('*'):
 if src.is_file() and src!=old:
  dst=repair/src.relative_to(original)
  assert src.read_bytes()==dst.read_bytes()
  unchanged.append(str(src.relative_to(original)))
audit=dict(operation='technical_gzip_repair_only',source_sha256=expected,original_library_file_id='libfile_6306ca7bb35081919d568f9766de9031',original_zip_sha256='94e330fbedf5c262fa9185e455568fa2ddd6e16048202efe685a82e29aa6c474',original_gzip_error=old_error,cells=len(c),sample_count_table_exact_match=True,unchanged_other_files=len(unchanged),models_refit=False,labels_changed=False,counts=c.groupby('review_scope').size().to_dict())
(repair/'technical_repair_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
(P/'provenance/prior_review_library_identity.json').write_text(json.dumps(dict(library_file_id=audit['original_library_file_id'],file_id='file_00000000748c81fdb493c54a6d8ddf57',file_name='lung_SCENIC3_review_complete.zip',sha256=audit['original_zip_sha256'],transfer='Current official Library helper on user project Linux, followed by authorized SSH copy to desktop; Linux identity metadata succeeded; Windows xattrs unsupported'),indent=2),encoding='utf-8')
review=pd.read_csv(P/'data/myeloid_working_annotation_audit.tsv',sep='\t')
review=review.drop(columns=['working_label_from_prior_review','review_action']).merge(ann.reset_index()[['cluster','label','scope','evidence']],on='cluster',validate='one_to_one')
review['review_action']='Reuse original delivered working mapping; no new relabeling. Confirm status if adopting as final atlas label.'
review.to_csv(P/'data/myeloid_working_annotation_audit.tsv',sep='\t',index=False)
models=pd.read_csv(original/'results/reviewed_sender_TF_models.tsv',sep='\t')
assert len(models)==12960 and models.FDR.notna().sum()==8208
assert (models.FDR<.05).sum()==0
models['tf_gene']=models.tf.str.split('(').str[0]
focus=models[models.tf_gene.isin(['IKZF2','STAT1'])].copy()
focus['nominal_p_lt_005']=focus.p.lt(.05) & focus.p.notna()
focus['fdr_lt_005']=focus.FDR.lt(.05) & focus.FDR.notna()
focus['acceptance_label']=np.where(focus.p.isna(),'NO_VALID_P',np.where(focus.p<.05,'NOMINAL_ONLY','NOT_NOMINAL'))
focus.to_csv(P/'data/prior_sender_TF_focused_with_acceptance.tsv',sep='\t',index=False)
main=focus[(focus.tf_gene=='IKZF2')&(focus.sender=='MacAM')&(focus.scope=='CD4-IL7R')&(focus.ligand=='IL18')]
main.to_csv(P/'data/fixed_IL18_IKZF2_all_sensitivities.tsv',sep='\t',index=False)
summary=focus.groupby(['tf_gene','model']).agg(rows=('p','size'),valid_P=('p','count'),nominal=('nominal_p_lt_005','sum'),FDR_positive=('fdr_lt_005','sum')).reset_index()
summary.to_csv(P/'data/sender_model_acceptance_inventory.tsv',sep='\t',index=False)
with zipfile.ZipFile('lung_SCENIC3_review_technical_repair_v1.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in sorted(repair.rglob('*')):
  if f.is_file():z.write(f,Path('scenic3_review')/f.relative_to(repair))
print(json.dumps(audit,indent=2))
print('FOCUSED_MODEL_ACCEPTANCE_INVENTORY');print(summary.to_string(index=False))
print('MAINLINE')
print(main[(main['filter']=='All')&(main.sender_filter=='All')].to_string(index=False))
for f in original.rglob('*.gz'):
 try:
  with gzip.open(f,'rb') as g:g.read()
  print('ORIGINAL_GZIP_OK',f.name)
 except (EOFError,OSError):print('ORIGINAL_GZIP_DAMAGED',f.name)

