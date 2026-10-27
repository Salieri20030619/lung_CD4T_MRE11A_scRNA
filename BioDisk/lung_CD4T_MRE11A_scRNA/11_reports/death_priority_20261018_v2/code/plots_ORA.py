from pathlib import Path
import pandas as pd,numpy as np,json,textwrap,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v2');F=O/'figures'
d=pd.read_csv(O/'tables/joint_ORA_all.tsv',sep='\t');cols=[(sc,di) for sc in ['Naive','ActCD4','Treg'] for di in ['up','down']];manifest=[];shown=[]
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':190})
for db,fid in [('GO_BP','ORA01'),('KEGG','ORA02')]:
 z=d[d.database==db];ids=[]
 for sc,di in cols:
  q=z[(z.scope==sc)&(z.direction==di)&(z.P<.05)].sort_values(['P','pathway']).head(2)
  ids.extend(q.pathway)
 ids=list(dict.fromkeys(ids));names=z.drop_duplicates('pathway').set_index('pathway').description
 fig,ax=plt.subplots(figsize=(9,7));labels=[]
 for i,term in enumerate(ids):
  name=names[term].replace(' - Mus musculus (house mouse)','');labels.append('\n'.join(textwrap.wrap(name,38)))
  for j,(sc,di) in enumerate(cols):
   q=z[(z.pathway==term)&(z.scope==sc)&(z.direction==di)]
   if not len(q):continue
   r=q.iloc[0]
   if r.P>=.05:continue
   ax.scatter(j,i,s=25+15*np.sqrt(r.hits),c=[min(-np.log10(max(r.FDR,1e-30)),5)],cmap='YlOrRd',vmin=0,vmax=5,edgecolors='black' if r.FDR<.05 else '#aaaaaa',lw=.7);shown.append({'figure':fid,**r.to_dict()})
 ax.set(xticks=range(6),xticklabels=[sc+'\n'+di for sc,di in cols],yticks=range(len(ids)),yticklabels=labels,title=db+' | top 2 nominal terms per state/direction')
 ax.invert_yaxis();ax.grid(alpha=.15);ax.tick_params(axis='y',labelsize=9);fig.colorbar(ScalarMappable(norm=Normalize(0,5),cmap='YlOrRd'),ax=ax,label='-log10 FDR (capped at 5)',shrink=.65)
 handles=[ax.scatter([],[],s=25+15*np.sqrt(v),c='#aaaaaa',label=str(v)) for v in [5,20,50]];ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(0,-.14),ncol=3,frameon=False,title='Overlapping genes')
 plt.tight_layout();plt.savefig(F/(fid+'.png'),bbox_inches='tight');plt.close();manifest.append({'id':fid,'title':db+'联合疾病差异ORA','source':'joint_ORA_all.tsv','file':str(F/(fid+'.png'))})
pd.DataFrame(shown).to_csv(O/'tables/ORA_display_values.tsv',sep='\t',index=False)
summ=d.groupby(['scope','database','direction']).agg(tests=('P','size'),nominal=('P',lambda s:int((s<.05).sum())),FDR=('FDR',lambda s:int((s<.05).sum()))).reset_index();summ.to_csv(O/'tables/ORA_count_summary.tsv',sep='\t',index=False)
(F/'manifest_ORA.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print(summ.to_string(index=False))
