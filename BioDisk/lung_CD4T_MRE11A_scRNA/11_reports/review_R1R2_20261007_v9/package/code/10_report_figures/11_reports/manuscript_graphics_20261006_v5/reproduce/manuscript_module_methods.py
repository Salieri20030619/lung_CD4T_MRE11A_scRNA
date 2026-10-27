from pathlib import Path
import re
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');T=P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release'
for name in ['M1_methods.txt','M2_methods.txt']:
 f=T/'provenance'/name
 if f.exists():print('FILE',str(f.relative_to(P)));print(f.read_text())
roots=[P/'02_analysis/CD4_DEG_collect_20260914_100836_XJWbue/source/scripts',P/'02_analysis/GSE314943/stage7_20260905_234200_aiQYL0',P/'02_analysis/GSE314943/stage8_20260906_003153_6mFKmg']
rx=re.compile(r'gseGO|gseKEGG|fgsea|rank|apopt|pyro|mito|singscore|GSVA|msig|CP10|log1p|log2|design|Wilcox|paired|lm\(|set.seed|nPerm|pAdjust|org\.Mm|GO:|KEGG|mean|term',re.I)
for root in roots:
 for f in root.rglob('*'):
  if f.is_file() and f.suffix.lower() in ['.r','.py','.txt'] and (f.suffix.lower() in ['.r','.py']) and (root.name!='scripts' or re.match(r'1[4-8]_',f.name)):
   ls=f.read_text(errors='replace').splitlines();ix=set()
   for i,l in enumerate(ls):
    if rx.search(l):ix.update(range(max(0,i-1),min(len(ls),i+2)))
   print('CODE',str(f.relative_to(P)));print('\n'.join(f'{i+1}: {ls[i]}' for i in sorted(ix))[:22000])
