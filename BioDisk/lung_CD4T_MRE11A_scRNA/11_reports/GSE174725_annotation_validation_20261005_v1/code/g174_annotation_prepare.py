from pathlib import Path
import json,importlib.util,importlib.metadata as im,collections
import pandas as pd
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/GSE174725_annotation_validation_20261005_v1'
for n in ['code','results','provenance','figures','deliverables']:(O/n).mkdir(parents=True,exist_ok=True)
print('PACKAGES',{p:(im.version(p) if importlib.util.find_spec(p) else None) for p in ['scanpy','anndata','scrublet','skimage','igraph','leidenalg'] if p!='skimage'})
for n in ['02_analysis/screen2_Trefine_nuzVQu/GSE174725/marker_detection.tsv','02_analysis/screen2_Trefine_nuzVQu/GSE174725/top40.tsv']:
 d=pd.read_csv(P/n,sep='\t');d=d[(d.resolution=='F08')&d.cluster.astype(str).isin(['0','1','12','2','3'])]
 if 'gene' in d and 'pct' in d and 'mean_count' in d:d=d[d.gene.isin('CD3D CD3E TRAC CD4 CD8A CD8B IL7R FOXP3 IL2RA CTLA4 NKG7 GNLY LYZ LST1 TYROBP FCER1G MS4A1 CD79A EPCAM'.split())]
 print('CD4_REVIEW',n);print(d.to_string(index=False)[:19000])
for folder in [P/'11_reports/human_gap_diagnostics_20261005_v1/results',P/'11_reports/human_first_20261004_v1/results']:
 for p in folder.iterdir():
  if any(w in p.name for w in ['members','target','background']) and p.suffix in ['.gz','.tsv']:
   d=pd.read_csv(p,sep='\t',nrows=2);print('PROGRAM_FILE',p.name,list(d.columns),d.to_dict('records'))
spec='''Stage A fixed annotation audit. No human main-cohort rerun. Reuse CSV recovery/mapping and old Core/Expanded CD4 labels. No existing formal BALF macrophage annotation found. Build one broad cohort graph at Leiden0.6 seed42 on common unambiguous feature space, all cells passing existing full-CSV QC (200features,300counts,MT<20); retain five donors. Exclude 12 fixed candidates and IKZF2/STAT1 plus MT/ribosomal from HVGs. Review broad positive/negative lineage panels, cluster profiles, per-donor dominance and mixed modules; one conditional myeloid-only graph at Leiden0.4 if broad macrophage/monocyte mixing requires it, not parameter search. Existing T F08 identities remain separate; no forced mapping to primary-cohort IL7R/GZMA states. Existing raw->Core count audit reused. Use installed per-donor Scrublet with fixed expected rate0.06 and seed42 if available, no threshold tuning; computational flags are uncertainty not biological doublet proof. Without empty droplets/no BAM, ambient correction or genotype-based contamination confirmation is unavailable. Retain original and flagged exclusions as separate sensitivity. Freeze biological direction/program test scope after annotation and before examining candidate group differences.\n'''
(O/'provenance/ANNOTATION_SPEC.txt').write_text(spec)
