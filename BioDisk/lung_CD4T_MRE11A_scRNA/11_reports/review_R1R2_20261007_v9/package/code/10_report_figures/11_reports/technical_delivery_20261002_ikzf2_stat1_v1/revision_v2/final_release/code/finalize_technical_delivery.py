from pathlib import Path
import json,hashlib,shutil,zipfile
import pandas as pd
import pymupdf

P=Path('technical_delivery')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
(P/'code').mkdir(exist_ok=True)
for name in ['extract_technical_delivery.py','repair_and_audit_review.py','build_technical_figures.py','build_delivery_report.py','finalize_technical_delivery.py','apply_revision.py','verify_saved_pdf.py']:
 shutil.copy2(name,P/'code'/name)
original=Path('library_sources/lung_SCENIC3_review_complete.zip')
assert sha(original)=='94e330fbedf5c262fa9185e455568fa2ddd6e16048202efe685a82e29aa6c474'
(P/'provenance/library_originals').mkdir(exist_ok=True)
shutil.copy2(original,P/'provenance/library_originals'/original.name)
(P/'provenance/library_originals/README.txt').write_text('Original Library archive retained byte-for-byte. Its reviewed_cells.tsv.gz is truncated. Use prior_review_repaired/scenic3_review/results/reviewed_cells.tsv.gz for the technically repaired cell table; all other original files and model results are unchanged.\n',encoding='utf-8')
caps=json.loads((P/'figure_manifest.json').read_text(encoding='utf-8'))
extra=' 人Core两TF在Source/SourceDep所有展示规格均未达名义P<0.05；STAT1负方向仅为描述，不能称显著反向。'
for info in caps:
 if info['figure']=='F05_regulon_activity_effects' and extra not in info['caption_zh']:info['caption_zh']+=extra
(P/'figure_manifest.json').write_text(json.dumps(caps,ensure_ascii=False,indent=2),encoding='utf-8')
(P/'图注与缩写说明.txt').write_text('\n\n'.join(i['figure']+' | '+i['title']+'\n'+i['caption_zh']+'\n来源：'+i['source'] for i in caps),encoding='utf-8')
with (P/'图注与缩写说明.txt').open('a',encoding='utf-8') as fh:
 fh.write('\n\n通用缩写：All=原不过滤臂；BG1=固定背景比例≤1%敏感性臂；s17/s29/s43=随机种子，不是独立生物重复；CP10k=按原全基因库归一化的每万计数；AUCell=regulon活性积分；LOO=逐一移除受试者/文库；N=名义P<0.05；F=原BH FDR<0.05；NoDom=原规则移除主导对照CTL068/CTL070；NoParent16=原父群16排除敏感性；CWP=原人肺病例标签，不改写为新诊断；SourceDep=组别/来源/接收深度调整；BothEpi=再加发送深度和发送上皮比例。\n')
repaired=P/'prior_review_repaired/scenic3_review'
old=P/'prior_review_original/scenic3_review'
unchanged=[]
for f in old.rglob('*'):
 if f.is_file() and f.name!='reviewed_cells.tsv.gz':
  target=repaired/f.relative_to(old)
  assert sha(f)==sha(target),str(f)
  unchanged.append(str(f.relative_to(old)))
cells=pd.read_csv(repaired/'results/reviewed_cells.tsv.gz',sep='\t')
models=pd.read_csv(repaired/'results/reviewed_sender_TF_models.tsv',sep='\t')
assert len(cells)==31120 and len(models)==12960 and len(unchanged)==22
audit=pd.read_csv(P/'data/nominal_acceptance_audit.tsv',sep='\t')
expected=audit.apply(lambda r:'NO_VALID_P' if pd.isna(r.p) else 'NOT_NOMINAL' if r.p>=.05 else 'FDR_PASS' if pd.notna(r.FDR) and r.FDR<.05 else 'NOMINAL_ONLY',axis=1)
assert (audit.acceptance==expected).all()
assert len(audit)==12391
new=audit[audit.evidence_type=='HUMAN_MODULE_GROUP_EFFECT'].sort_values('source_row_0based')
src=pd.read_csv(P/'source_tables/M2_human_source_models.tsv',sep='\t')
src=src[src.metric.str.startswith('score_')]
assert len(new)==504
for target,source_column in [('p','p_nominal'),('FDR','FDR'),('estimate','coefficient'),('se','se'),('CI_low','CI_low'),('CI_high','CI_high'),('df','residual_df'),('status','status'),('FDR_status','FDR_status'),('n_control','n_control'),('n_case','n_case')]:
 pd.testing.assert_series_equal(new[target].reset_index(drop=True),src[source_column].reset_index(drop=True),check_names=False,check_dtype=False)
original_audit=pd.read_csv(Path('../technical_delivery/data/nominal_acceptance_audit.tsv'),sep='\t')
pd.testing.assert_frame_equal(audit.iloc[:len(original_audit)][original_audit.columns].reset_index(drop=True),original_audit,check_dtype=False)
for file in (P/'figures').glob('*'):
 assert sha(file)==sha(Path('../technical_delivery/figures')/file.name),str(file)
assert ((models.FDR<.05)&models.FDR.notna()).sum()==0
doc=pymupdf.open(P/'output/pdf/IKZF2_STAT1_technical_delivery_zh.pdf')
assert len(doc)==14
outliers=[]
for idx,page in enumerate(doc):
 for b in page.get_text('dict')['blocks']:
  for line in b.get('lines',[]):
   for span in line['spans']:
    r=pymupdf.Rect(span['bbox'])
    if r.x0<0 or r.y0<0 or r.x1>page.rect.width+1 or r.y1>page.rect.height+1:outliers.append(idx+1)
assert not outliers
qa={'pdf_pages':14,'new_figures':8,'visual_review':'All 14 final PDF pages rendered with Poppler and visually inspected individually. No layout defect was established in v1; this revision adds receptor limitations, M2 audit coverage, and scope/direction interpretation.','off_page_text_spans':0,'acceptance_rows':len(audit),'acceptance_flags_match_original_p_fdr':True,'reviewed_cells_rows':len(cells),'original_model_rows':len(models),'original_files_byte_unchanged':unchanged,'original_library_zip_sha256':sha(original),'new_statistical_models_fit':0,'scientific_thresholds_or_groups_changed':False,'working_labels_changed':False,'source_hash_count':116,'excluded_from_release':'Unpacked damaged original review working directory; exact original ZIP retained with repair warning.','patient_fdr_family':'Original complete family including ETV7; no candidate-only BH recalculation.'}
(P/'provenance/final_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
revision={'version':'v2','added_original_model_rows':504,'total_acceptance_rows':12391,'M2_valid_p':396,'M2_no_valid_p':108,'M2_nominal':18,'M2_FDR_status':'NOT_COMPUTED_GROUP_N_LT10 for all 504 original rows','original_11887_audit_rows_unchanged':True,'original_16_figure_files_byte_unchanged':True,'receptor_scope':'Human All CD4-IL7R, >=20 cells per sample','receptor_samples':16,'IL18R1_IL18RAP_codetection_median_pct':0,'IL18R1_IL18RAP_codetection_range_pct':[0,2.127659574468085],'model_sample_denominator':'Fixed sender-TF All/Matched n=11; differs from receptor-only n=16','interpretation':'Core IKZF2-PyroExec/ApoExec inverse associations are not a positive causal chain; CD4-IL7R sender association and Core module associations cannot be joined as one tested chain.','layout':'No v1 layout defect established; preserved figure layouts. All revised PDF pages visually checked; human heatmap also checked with a second renderer and embedded image pixels match original PNG.'}
(P/'provenance/revision_v2_audit.json').write_text(json.dumps(revision,indent=2),encoding='utf-8')
for note in P.glob('*.txt'):
 if ' v1' in note.read_text(encoding='utf-8'):note.write_text(note.read_text(encoding='utf-8').replace(' v1',' v2').replace('所有保留模型','本轮纳入审计的既有模型'),encoding='utf-8')
status=json.loads((P/'technical_delivery_status.json').read_text())
status.update({'verification':'PASS','original_models_unchanged':True,'release_version':'v2'})
(P/'technical_delivery_status.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
files=[p for p in P.rglob('*') if p.is_file() and 'prior_review_original' not in p.parts and p.name!='SHA256SUMS.tsv']
manifest='path\tbytes\tsha256\n'+''.join(f'{f.relative_to(P).as_posix()}\t{f.stat().st_size}\t{sha(f)}\n' for f in sorted(files))
(P/'SHA256SUMS.tsv').write_text(manifest,encoding='utf-8')
files.append(P/'SHA256SUMS.tsv')
zpath=Path('lung_IKZF2_STAT1_technical_delivery_v1.zip')
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in sorted(files):z.write(f,f.relative_to(P).as_posix())
with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
deliverables=[zpath,P/'output/pdf/IKZF2_STAT1_technical_delivery_zh.pdf',Path('lung_SCENIC3_review_technical_repair_v1.zip')]
release={'archive_members':len(files),'files':[{'path':str(p.resolve()),'bytes':p.stat().st_size,'sha256':sha(p)} for p in deliverables]}
Path('technical_release_verification.json').write_text(json.dumps(release,indent=2),encoding='utf-8')
print(json.dumps(release,indent=2))
