from pathlib import Path
import zipfile,hashlib,pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
W=Path('C:/Users/31937/lung_four_axes_20261017_v15');f=W/'v15_final_delta.zip';assert hashlib.sha256(f.read_bytes()).hexdigest()=='d0974570caa1d629cb2dab40a8b5470cc251af385cdc9ec461a13d2a48e1880c'
with zipfile.ZipFile(f) as z:assert z.testzip() is None;z.extractall(W/'tables')
T=W/'tables';F=W/'figures';a=pd.read_csv(T/'06_novelty/decontx_Treg_program_samples.tsv',sep='\t');ret=pd.read_csv(T/'06_novelty/decontx_Treg_gene_retention.tsv',sep='\t');ids=['Saline_1','Saline_2','Saline_3','Silica_1','Silica_2','Silica_3'];col=['#3f91b6']*3+['#ce7357']*3
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180});fig,axs=plt.subplots(1,3,figsize=(13,4.6))
for ax,variant,title in zip(axs[:2],['Raw','Corrected_new_library'],['Original data\nrho=-0.943; P=0.0167','Model-corrected sensitivity\nrho=-0.657; P=0.175']):
 b=a[a.variant==variant].set_index('sample').loc[ids];ax.scatter(b.DUSP5,b.PyroExec,c=col,s=40)
 for s,x,y,c in zip(ids,b.DUSP5,b.PyroExec,col):ax.annotate(s.replace('Saline_','C').replace('Silica_','S'),(x,y),xytext=(4,4),textcoords='offset points',fontsize=8,color=c)
 ax.set_xlabel('Treg DUSP5 CP10k');ax.set_ylabel('PyroExec score');ax.set_title(title)
keep=['DUSP5','PYCARD','CASP1','GSDMD','IL1B','IL18'];q=ret.groupby('gene')[['raw_UMI','corrected_UMI']].sum().loc[keep];pct=100*q.corrected_UMI/q.raw_UMI;axs[2].bar(np.arange(len(q)),pct,color=['#287e91']+['#91b3bf']*3+['#c66c4d','#91b3bf']);axs[2].set_xticks(np.arange(len(q)),q.index,rotation=45,ha='right');axs[2].set_ylim(0,115);axs[2].set_ylabel('Corrected / original UMI (%)');axs[2].set_title('Treg transcript retention\nPooled counts, descriptive only')
for j,v in enumerate(pct):axs[2].text(j,v+2,f'{v:.1f}',ha='center',fontsize=8)
fig.tight_layout();fig.savefig(F/'F03_model_diagnostic.png');fig.savefig(F/'F03_model_diagnostic.pdf');plt.close(fig)
print('DELTA_AND_DIAGNOSTIC_FIGURE_READY')