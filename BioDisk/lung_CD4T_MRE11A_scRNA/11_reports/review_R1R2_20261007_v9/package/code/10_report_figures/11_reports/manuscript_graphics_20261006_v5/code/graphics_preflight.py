from pathlib import Path
import subprocess,json
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
for env in ['r_4.4','r_seurat','r_seurat_new','r_analysis','R_4.5.2']:
 exe=Path('/mnt/BioDisk/wangmingyuan/miniconda3/envs')/env/'bin/Rscript'
 if exe.exists():
  r=subprocess.run([str(exe),'-e','p<-c("SeuratObject","Seurat","enrichplot","DOSE","ggplot2","jsonlite","cowplot");cat(R.version.string,"\\n"); for(k in p)cat(k,if(requireNamespace(k,quietly=TRUE))as.character(packageVersion(k))else"missing","\\n")'],capture_output=True,text=True,timeout=90)
  print('R_ENV',env,r.stdout,r.stderr[-1200:])
for name,starts in [('scripts/16_GSE288887_CD4_pathway_analysis.R',[(450,482),(634,707),(1595,1613)]),('02_analysis/GSE314943/stage7_20260905_234200_aiQYL0/analysis.R',[(1,48),(83,110)])]:
 ls=(P/name).read_text().splitlines();print('SOURCE',name)
 for a,b in starts:print('\n'.join(f'{i+1}: {ls[i]}' for i in range(a-1,min(b,len(ls)))))
for f in (P/'scripts').glob('*.R'):
 if 'lineage_prelim' in f.read_text(errors='replace'):
  ls=f.read_text().splitlines();print('LINEAGE_SOURCE',f.name);print('\n'.join(f'{i+1}: {l}' for i,l in enumerate(ls) if any(s in l for s in ['lineage_prelim','readRDS','saveRDS','scored.rds'])))
