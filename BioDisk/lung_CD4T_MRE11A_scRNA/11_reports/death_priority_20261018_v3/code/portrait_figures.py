import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MPLCONFIGDIR='/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/tmp/matplotlib_v3')
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');V=P/'11_reports/death_priority_20261018_v2';O=P/'11_reports/death_priority_20261018_v3';F=O/'figures';F.mkdir(exist_ok=True)
original_subplots=plt.subplots
def portrait(*args,**kwargs):
 args=list(args)
 nr= args[0] if len(args)>0 else kwargs.get('nrows',1)
 nc= args[1] if len(args)>1 else kwargs.get('ncols',1)
 old=kwargs.get('figsize',(8,5))
 if nr==1 and nc>1:
  if len(args)>1:args[0],args[1]=nc,1
  else:kwargs.update(nrows=nc,ncols=1)
  kwargs['figsize']=(7.6,min(11,3*nc))
 elif nr>1 and nc>1:kwargs['figsize']=(8,10)
 else:kwargs['figsize']=(8,max(4.3,min(9,old[1])))
 return original_subplots(*args,**kwargs)
plt.subplots=portrait
for name in []:
 source=(V/'code'/name).read_text(encoding='utf-8')
 source=source.replace("F=O/'figures'","F=Path('"+str(F)+"')")
 exec(compile(source,name,'exec'),{})
 print('PORTRAIT',name,flush=True)
# Draw only the final network's existing result figures; do not rerun or modify numerical QA.
source=(V/'code/joint_network_finish.py').read_text(encoding='utf-8')
source=source[source.index('fig,axs=plt.subplots(1,2'):source.index("status={'status'")]
import pandas as pd,numpy as np,json,re
O=V;coverage=pd.read_csv(V/'tables/joint_network_TF_coverage.tsv',sep='\t');runs=pd.read_csv(V/'tables/joint_network_run_summary.tsv',sep='\t');effects=pd.read_csv(V/'tables/joint_network_TF_disease_effects.tsv',sep='\t');edges=pd.read_csv(V/'tables/joint_network_candidate_edges.tsv',sep='\t');aa=pd.read_csv(V/'tables/joint_network_TF_associations.tsv',sep='\t');figs=[]
def save(name,title,source):
 plt.savefig(F/(name+'.png'),bbox_inches='tight',dpi=200,facecolor='white');plt.close();figs.append(dict(id=name,title=title,source=source,file=str(F/(name+'.png'))))
source=source.replace("sel.to_csv(O/'tables/joint_network_display_edges.tsv',sep='\\t',index=False)","None")
# Exclude the existing display-table write regardless of whitespace.
source='\n'.join(line.replace("sel.to_csv(O/'tables/joint_network_display_edges.tsv',sep='\\t',index=False)","None") for line in source.splitlines())
source=re.sub(r';sel.to_csv[^\n]*','',source)
run=runs
exec(compile(source,'network_figures','exec'))
print('PORTRAIT_FIGURES_COMPLETE',flush=True)
