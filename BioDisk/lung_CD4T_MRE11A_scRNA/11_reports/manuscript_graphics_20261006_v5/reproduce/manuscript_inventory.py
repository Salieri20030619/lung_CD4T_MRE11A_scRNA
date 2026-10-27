from pathlib import Path
import re,json,subprocess
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
pattern=re.compile(r'gsea|fgsea|apopt|pyro|mito|oxphos|geneset|module|pathway',re.I)
roots=[P/'02_analysis',P/'11_reports'];found=[]
for root in roots:
 for f in root.rglob('*'):
  if f.is_file() and pattern.search(f.name) and f.suffix.lower() in ['.tsv','.csv','.r','.py','.txt','.json','.png','.pdf','.gmt','.md','.html']:
   found.append({'path':str(f.relative_to(P)),'bytes':f.stat().st_size})
print('MATCHING_FILES',json.dumps(found,ensure_ascii=False))
print('PROCESSES',subprocess.run(['ps','-eo','pid,etime,args'],capture_output=True,text=True).stdout[:0])
proc=subprocess.run(['ps','-eo','pid,etime,args'],capture_output=True,text=True).stdout.splitlines()
print('\n'.join(x for x in proc if str(P) in x and not any(k in x for k in ['base64','ps -eo','sshd:'])))
