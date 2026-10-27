from pathlib import Path
import json,hashlib,gzip,shutil,datetime
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/GSE296682_patient_localization_20261005_v1'
for d in ['code','results','provenance','figures','deliverables','public_cache']:(O/d).mkdir(parents=True,exist_ok=True)
spec='''Frozen before expression values are read. Patient-side localization only in two metadata Silicosis samples B2_RZ725 and B6_RZ725; shared RZ725_Pool4/10xv3 is not independent batches. Preserve official annotations and review marker support separately, no new whole-library clustering/ligand discovery, no n=2 correlations, no disease DE, no new NicheNet/SCENIC/virtual KO.
Fixed candidates IL15 TNF CLSTN1 ADAM15 PLAU DSC2 CXCL2 ICAM2 VTN VEGFB GSTP1 RTN4. Fixed components/TF IL15RA IL2RB IL2RG TNFRSF1A TNFRSF1B IKZF2 STAT1. Review panel T: CD3D CD3E TRAC CD4 IL7R LTB CCR7 SELL TCF7; CD8/NK: CD8A CD8B NKG7 GNLY KLRD1 PRF1; Treg: FOXP3 IL2RA CTLA4 TIGIT; myeloid: LYZ LST1 TYROBP AIF1 FCER1G; mature Mac: C1QA C1QB C1QC APOE APOC1 MRC1 MSR1 CD68 CD163 CSF1R MERTK TREM2 GPNMB LGMN SPP1 FOLR2 LYVE1; AM: FABP4 PPARG MARCO; Mono: FCN1 VCAN S100A8 S100A9; competitor DC/epithelial/platelet: CD1C FCER1A CLEC10A EPCAM KRT8 KRT18 SCGB1A1 PPBP PF4. Stable ID/symbol audit; do not guess unresolvable aliases.
Use author's retained cells and original labels as primary; separately record lineage mixtures based on concordant multigene evidence, not candidate expression. No automatic relabel monocyte-like macrophages as mature or monocyte-derived. Treg is not assumed CD4 from metadata; CD4 helper n2/3 and Treg n14/3 descriptive, no strong quantitative conclusions. No forced main lung GZMA/IL7R correspondence.
For each donor and author major/minor label report n/N positive, continuous count and/or expression mean/median/quantiles, feature/layer provenance, exact denominator. Prefer bona fide count layer if available; if transformed/adjusted, label it explicitly and do not invent raw UMIs. No pooled donor statistical tests or P/FDR. Display actual expression first; fixed 5%/10% coverage sensitivity afterward, no outcome-driven thresholds. Co-detection same cell for sender IL15+IL15RA, receiver IL2RB+IL2RG, TNFRSF1A+TNFRSF1B, and TCR+CD4/CD8/NK context. Marginal or joint RNA expression does not establish protein trans-presentation, signaling or causality. Zero detection does not prove biological absence.
Only minimum official data needed, server cache first; inspect checksums and object layers. No FASTQ, duplicate RDS/H5AD downloads, local matrix staging, public upload of private data, or access-control bypass. Only Chinese report and necessary aggregate figures may transfer; matrices, donor/cell detail and ZIP remain server.
'''
f=O/'provenance/FROZEN_LOCALIZATION_SPEC.txt'
if not f.exists():f.write_text(spec)
else:assert f.read_text()==spec
print('FROZEN_SPEC',str(f),'UTC',datetime.datetime.now(datetime.timezone.utc).isoformat())
def hs(p):
 m=hashlib.md5();s=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):m.update(b);s.update(b)
 return dict(path=str(p),bytes=p.stat().st_size,md5=m.hexdigest(),sha256=s.hexdigest())
rows=[]
for p in (P/'00_raw/GSE296682').iterdir():
 if p.is_file():rows.append(hs(p))
(O/'provenance/existing_cache_hashes.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
for p in [P/'02_analysis/screen6_structure/GSE296682_structure.json',P/'02_analysis/screen6_prescreen_20260909_153209_hPfqvr/GSE296682_sample_metadata.json']:
 print('PRIOR_STRUCTURE',str(p),p.read_text()[:18000]);(O/'provenance'/p.name).write_bytes(p.read_bytes())
for p in (P/'logs').glob('*GSE296682*'):
 print('DOWNLOAD_PROVENANCE',p.name,p.read_text(errors='replace')[-2500:]);(O/'provenance'/p.name).write_bytes(p.read_bytes())
q=P/'00_raw/GSE296682/GSE296682_BAL_FINAL_metadata.txt.gz'
with gzip.open(q,'rb') as f:data=f.read()
print('METADATA_EXPANDED',len(data),hashlib.md5(data).hexdigest(),hashlib.sha256(data).hexdigest(),data[:1300].decode(errors='replace'))
