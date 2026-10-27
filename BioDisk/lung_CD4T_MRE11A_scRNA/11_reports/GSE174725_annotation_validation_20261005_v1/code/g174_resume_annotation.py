from pathlib import Path
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/GSE174725_annotation_validation_20261005_v1');t=(O/'code/g174_build_annotation.py').read_text();exec(t[:t.index('mp=pd.read_csv')])
checkpoint=R/'graph_checkpoint_SERVER_ONLY.h5ad'
if checkpoint.exists():
 a=ad.read_h5ad(checkpoint);doublets=pd.read_csv(R/'scrublet_per_donor_audit.tsv',sep='\t').to_dict('records');exec(t[t.index("C=a.layers['counts'];X=a.X"):])
else:
 a=ad.read_h5ad(R/'common_counts_QC_SERVER_ONLY.h5ad');exec(t[t.index('# Scanpy includes'):])
