from pathlib import Path
import pandas as pd,numpy as np,json,itertools,textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.stats import false_discovery_control,hypergeom
O=Path(__file__).parent;H=O.parent;T=O/'tables';U=O/'server_tables';F=O/'figures';F.mkdir(exist_ok=True);PD=O/'plot_data';PD.mkdir(exist_ok=True)
V=H/'lung_append_analysis_20261017_v18';B=H/'lung_three_cohort_20261017_v16/data';W=H/'lung_feedback_revision_20261017_v19';TR=O/'sources/11_reports/customer_sequence_20261015_v2/03A_state_pseudotime'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
colors={'Saline':'#3278a8','Silica':'#d46746'};states=['Naive','ActCD4','Treg'];statecolors={'Naive':'#4c9e9c','ActCD4':'#dc8451','Treg':'#8b74aa'};manifest=[];qa=[]
def rd(p):return pd.read_csv(p,sep='\t',low_memory=False)
def save(fig,id,module,title,caption,data=None):
 fig.tight_layout();fig.savefig(F/(id+'.png'),dpi=220,bbox_inches='tight');fig.savefig(F/(id+'.pdf'),bbox_inches='tight');plt.close(fig)
 if data is not None:data.to_csv(PD/(id+'.tsv'),sep='\t',index=False,encoding='utf-8')
 manifest.append(dict(id=id,module=module,title=title,caption=caption,file=str(F/(id+'.png'))))
def dots(ax,z,y):
 for k,gr in enumerate(['Saline','Silica']):
  a=z[z.group==gr];ax.scatter(k+np.linspace(-.09,.09,len(a)),a[y],c=colors[gr],s=35,zorder=3);ax.hlines(a[y].mean(),k-.16,k+.16,color='black',lw=1.2)
 ax.set_xticks([0,1],['Saline','Silica'])
# 1 Atlas
atlas=rd(U/'atlas_cells.tsv.gz');ct=atlas.Broad.value_counts();fig,ax=plt.subplots(figsize=(11,7));palette=plt.cm.tab20(np.linspace(0,1,len(ct)))
for color,label in zip(palette,ct.index):
 z=atlas[atlas.Broad==label];ax.scatter(z.UMAP1,z.UMAP2,s=2,color=color,label=f'{label} ({len(z):,})',rasterized=True,alpha=.65)
ax.set(xlabel='UMAP 1',ylabel='UMAP 2',title='GSE314943 | frozen whole-lung annotation | 37,155 cells');ax.legend(bbox_to_anchor=(1.01,1),loc='upper left',markerscale=4,fontsize=9.5,frameon=False)
save(fig,'F01','免疫全景','全肺图谱定位主要免疫群及其组织环境','主队列3个对照、3个矽尘样本，共37,155个质控后细胞。使用冻结注释及UMAP坐标；Like/Unresolved保留其身份不确定性。细胞数是图谱规模，不是统计样本量。',atlas[['sample','group','Broad','UMAP1','UMAP2']])
mk=rd(U/'atlas_markers.tsv');genes=['CD3D','TRAC','MS4A1','CD79A','NKG7','LYZ2','CSF1R','S100A8','FCER1A','SFTPC','PECAM1','COL1A1'];a=mk[mk.gene.isin(genes)].groupby(['state','gene'],as_index=False)[['fraction','mean_logCP10k']].mean();order=ct.index.tolist();fig,ax=plt.subplots(figsize=(11,7));norm=Normalize(0,a.mean_logCP10k.max())
for r in a.itertuples():ax.scatter(genes.index(r.gene),order.index(r.state),s=5+180*r.fraction,c=[r.mean_logCP10k],cmap='viridis',norm=norm)
ax.set_xticks(range(len(genes)),genes,rotation=55,ha='right');ax.set_yticks(range(len(order)),order);ax.invert_yaxis();ax.set_title('Lineage markers | sample-balanced expression and detection');fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='viridis'),ax=ax,label='Mean log1p(CP10k), sample-balanced')
save(fig,'F02','免疫全景','标志基因支持主要谱系的身份区分','点大小为样本内阳性比例再跨样本平均，颜色为逐细胞log1p(CP10k)的样本均值再平均。用于核查冻结谱系注释，不以单一标志基因重新定义细胞。',a)
# 2 CD4 states
cd=rd(U/'cd4_cells.tsv.gz');fig,axs=plt.subplots(1,2,figsize=(11,4.7))
for st in states:
 z=cd[cd.Subtype==st];axs[0].scatter(z.UMAP1,z.UMAP2,s=7,color=statecolors[st],label=f'{st} (n={len(z)})',alpha=.6,rasterized=True)
for gr in colors:
 z=cd[cd.group==gr];axs[1].scatter(z.UMAP1,z.UMAP2,s=7,color=colors[gr],label=gr,alpha=.5,rasterized=True)
for ax in axs:ax.set(xlabel='UMAP 1',ylabel='UMAP 2');ax.legend(frameon=False,fontsize=9)
axs[0].set_title('Core CD4 states (2,311 cells)');axs[1].set_title('Sample group overlay')
save(fig,'F03','CD4亚群','核心CD4群分为Naive、ActCD4和Treg状态','使用主队列冻结Harmony UMAP。Naive 933、ActCD4 801、Treg 577；嵌入用于状态展示，不把空间距离当作细胞转化或疾病效应检验。',cd[['sample','group','Subtype','UMAP1','UMAP2']])
cnt=cd.groupby(['sample','group','Subtype'],observed=True).size().rename('cells').reset_index();cnt['fraction']=cnt.cells/cnt.groupby('sample').cells.transform('sum');samp=sorted(cnt['sample'].unique());wide=cnt.pivot(index='sample',columns='Subtype',values='fraction').reindex(samp);fig,axs=plt.subplots(1,2,figsize=(12,4.6));bottom=np.zeros(len(samp))
for st in states:axs[0].bar(range(len(samp)),wide[st],bottom=bottom,color=statecolors[st],label=st);bottom+=wide[st].to_numpy()
axs[0].set_xticks(range(len(samp)),samp,rotation=45,ha='right');axs[0].set_ylabel('Fraction within core CD4');axs[0].legend(frameon=False);comp=[]
for j,st in enumerate(states):
 z=cnt[cnt.Subtype==st].sort_values('sample');x=z.fraction.to_numpy();obs=x[3:].mean()-x[:3].mean();null=[]
 for ids in itertools.combinations(range(6),3):m=np.zeros(6,bool);m[list(ids)]=True;null.append(x[m].mean()-x[~m].mean())
 pv=np.mean(np.abs(null)>=abs(obs)-1e-12);comp.append(dict(state=st,effect=obs,P=pv,n=6,family='CD4_state_proportions'))
 for off,gr in [(-.13,'Saline'),(.13,'Silica')]:
  zz=z[z.group==gr];axs[1].scatter(j+off+np.linspace(-.04,.04,len(zz)),zz.fraction,c=colors[gr],s=35)
 axs[1].text(j,1.02,f'P={pv:.2g}',ha='center',transform=axs[1].get_xaxis_transform(),fontsize=9)
cc=pd.DataFrame(comp);cc['FDR']=false_discovery_control(cc.P);cc.to_csv(T/'CD4_proportion_tests.tsv',sep='\t',index=False);axs[1].set_xticks(range(3),states);axs[1].set_ylabel('Independent-sample fraction');axs[1].set_ylim(0,1)
save(fig,'F04','CD4亚群','疾病组CD4状态组成在独立样本层面呈现变化','每点为独立样本，比例分母为该样本Core CD4。双侧精确组标签排列20种，FDR在3个亚群家族内计算；3+3样本按探索性方向与效应解释。比例不等于组织绝对细胞数量。',cnt)
# 3 DE
DE=rd(B/'tables/DE_all_models.tsv');fig,axs=plt.subplots(1,2,figsize=(11,4.8));focus={'Treg':['IL1R2','SKIL','FOXP1'],'ActCD4':['TTC39C','MDFIC','DDIT4','GSN']};deplot=[]
for ax,sc in zip(axs,['Treg','ActCD4']):
 z=DE[(DE.model=='joint_stratified')&(DE.scope==sc)].copy();z['neglogP']=-np.log10(z.PValue.clip(lower=1e-300));c=np.where(z.FDR<.05,'#bc4d42',np.where(z.PValue<.05,'#d6a342','#b9c1c7'));ax.scatter(z.logFC,z.neglogP,s=9,c=c,alpha=.6,rasterized=True);ax.axhline(-np.log10(.05),color='#666666',ls='--',lw=.7);ax.axvline(0,color='#888888',lw=.6)
 for gi,g in enumerate(focus[sc]):
  r=z[z.gene==g].iloc[0];ax.scatter([r.logFC],[r.neglogP],color='#172e45',s=24);ax.annotate(g,(r.logFC,r.neglogP),xytext=(8,18+gi*13 if g in ['TTC39C','MDFIC'] else 5),textcoords='offset points',fontsize=9,arrowprops=dict(arrowstyle='-',color='#555555',lw=.5))
 ax.set(title=sc+' | joint stratified (n=10)',xlabel='log2 fold change (Silica / Saline)',ylabel='-log10 raw P');deplot.append(z)
save(fig,'F05','差异分析','Treg候选的疾病差异与ActCD4程序候选呈不同证据特征','联合edgeR模型控制队列/时间层，含原6样本及183682的4样本。红为FDR<0.05，黄为仅原始P<0.05，灰为其余；完整原家族FDR保留。联合结果不算独立复现。',pd.concat(deplot))
f=rd(W/'tables/R05_current_focus_conclusions.tsv');sel=f[((f.gene.isin(['IL1R2','SKIL','FOXP1']))&(f.scope=='Treg'))|((f.gene.isin(['TTC39C','MDFIC','DDIT4','GSN']))&(f.scope=='ActCD4'))].copy();fig,ax=plt.subplots(figsize=(11,5));yy=np.arange(len(sel))
for off,prefix,color,label in [(-.12,'314_disease','#3278a8','Original six'),(.12,'joint_stratified','#d46746','Joint ten (includes six)')]:
 ax.scatter(sel[prefix+'_logFC'],yy+off,color=color,label=label,s=50)
 for k,r in enumerate(sel.itertuples(index=False)):
  rr=sel.iloc[k];ax.text(4.1,k+off,f"P={rr[prefix+'_PValue']:.3g}; FDR={rr[prefix+'_FDR']:.3g}",fontsize=8,va='center')
ax.set_yticks(yy,[g+' / '+s for g,s in zip(sel.gene,sel.scope)]);ax.invert_yaxis();ax.axvline(0,c='grey',lw=.7);ax.set_xlim(-2.2,7.6);ax.set_xlabel('log2 fold change (right-side text: original P/FDR)');ax.legend(loc='upper center',bbox_to_anchor=(.6,1.15),ncol=2,frameon=False)
save(fig,'F06','差异分析','IL1R2与SKIL具有疾病差异依据，TTC39C优势来自程序关联','并列原六样本与联合十样本效应，二者样本重叠。TTC39C联合DE P=0.163、FDR=0.803，不能用其凋亡程序关联替代疾病差异证据。',sel)
# 4 ORA / 5 GSEA; current v16 joint DE enrichment
ora=rd(U/'ORA_all.tsv');gsea=rd(U/'GSEA_all.tsv');members=rd(U/'enrichment_members.tsv')
for family,z in ora.groupby('family'):
 err=float(np.max(np.abs(false_discovery_control(z.P)-z.FDR)));assert err<1e-10;assert np.allclose(hypergeom.sf(z.k-1,z.N,z.K,z.n),z.P,atol=1e-10);qa.append(dict(check='ORA_hypergeom_BH',family=family,error=err))
for family,z in gsea.groupby('family'):
 ok=z.pval.notna();err=float(np.max(np.abs(false_discovery_control(z.loc[ok,'pval'])-z.loc[ok,'padj'])));assert err<1e-10;qa.append(dict(check='GSEA_BH',family=family,error=err))
for num,db in [(7,'GO_BP'),(8,'KEGG')]:
 fig,axs=plt.subplots(1,2,figsize=(13,6));selected=[]
 for ax,sc in zip(axs,['Treg','ActCD4']):
  z=ora[(ora.database==db)&(ora.scope==sc)&(ora.k>0)].sort_values(['P','pathway','direction']).head(10).copy();z['plot_label']=z.description.str.replace(' - Mus musculus (house mouse)','',regex=False);y=np.arange(len(z));ax.scatter(-np.log10(z.P.clip(lower=1e-300)),y,s=20+z.k*5,c=z.direction.map({'up':'#d46746','down':'#3278a8'}));ax.set_yticks(y,[textwrap.fill(s,32)+' ['+d+']' for s,d in zip(z.plot_label,z.direction)],fontsize=9.5);ax.invert_yaxis();ax.set(xlabel='-log10 raw ORA P',title=sc+' | '+db)
  for i,rr in enumerate(z.itertuples()):ax.annotate(f'q={rr.FDR:.2g}',(-np.log10(max(rr.P,1e-300)),i),xytext=(5,0),textcoords='offset points',fontsize=9,va='center')
  ax.margins(x=.38);selected.append(z)
 save(fig,f'F{num:02d}','GO/KEGG',('GO生物过程' if db=='GO_BP' else 'KEGG通路')+'揭示疾病差异涉及的功能集合','基于当前联合模型原始P<0.05的差异基因，按上/下调分别作超几何检验。背景为同模型通过表达过滤且可映射的受检基因；集合10-500成员；每个状态/数据库/方向完整家族BH校正。每状态按原P展示前10条，不以候选主题筛选；点面积反映命中基因数。',pd.concat(selected))
fig,axs=plt.subplots(1,2,figsize=(13,6));gsselect=[]
for ax,sc in zip(axs,['Treg','ActCD4']):
 z=gsea[gsea.scope==sc].sort_values(['pval','pathway']).head(12).copy();labs=[textwrap.fill(s.replace(' - Mus musculus (house mouse)',''),35) for s in z.description];ax.barh(np.arange(len(z)),z.NES,color=np.where(z.NES>0,'#d46746','#3278a8'));ax.set_yticks(np.arange(len(z)),labs,fontsize=9.5);ax.invert_yaxis();ax.axvline(0,c='grey',lw=.7);ax.set(xlabel='Normalized enrichment score (NES)',title=sc+' | ranked all-tested genes');gsselect.append(z)
save(fig,'F09','GSEA','全受检基因排序补充差异阈值之外的通路变化','排名为sign(logFC)×sqrt(edgeR QL F)。fgseaMultilevel，10-500成员，seed=20261017，P估计下限1e-8；GO与KEGG分家族校正。每状态按P展示前12项，正NES偏向疾病上调端。完整P、FDR和log2err在附表；达到下限的P解读为≤1e-8。',pd.concat(gsselect))
fig,axs=plt.subplots(1,2,figsize=(12,4.5));curve_records=[]
for ax,sc in zip(axs,['Treg','ActCD4']):
 r=gsea[(gsea.scope==sc)&(gsea.database=='KEGG')].sort_values(['pval','pathway']).iloc[0];ranks=rd(U/('GSEA_ranks_'+sc+'.tsv'));gsset=set(members[(members.database=='KEGG')&(members.pathway==r.pathway)].gene);hit=ranks.gene.isin(gsset).to_numpy();weight=abs(ranks.rank_stat.to_numpy());increments=np.where(hit,weight/max(weight[hit].sum(),1e-99),-1/(~hit).sum());es=np.cumsum(increments);peak=es.max() if es.max()>-es.min() else es.min();assert abs(peak-r.ES)<1e-8;qa.append(dict(check='GSEA_running_ES',family=sc+'|'+r.pathway,error=abs(peak-r.ES)));ax.plot(np.arange(len(es))+1,es,c='#3278a8',lw=1.4);ax.axhline(0,c='grey',lw=.6);ax.vlines(np.flatnonzero(hit)+1,es.min()-.04,es.min()-.025,color='#777777',lw=.35);ax.set(title=sc+' | '+textwrap.fill(r.description.replace(' - Mus musculus (house mouse)',''),38)+f'\nNES={r.NES:.2f}; P={r.pval:.2g}; FDR={r.padj:.2g}',xlabel='Rank position (Silica-up to Silica-down)',ylabel='Running enrichment score');curve_records.extend([dict(scope=sc,pathway=r.pathway,rank=i+1,gene=g,hit=bool(hh),running_ES=float(e)) for i,(g,hh,e) in enumerate(zip(ranks.gene,hit,es))])
save(fig,'F10','GSEA','代表通路的运行富集曲线展示方向与成员位置','每状态选择原P最小的KEGG通路，来自图9的同一完整排序及同一统计结果。运行ES由成员及排名权重独立重建，并与fgsea原ES一致；曲线不产生额外检验。',pd.DataFrame(curve_records))
# 6 communication: qualified expression counts, NOT inferred functional strength
lr=rd(B/'tables/expressed_LR_with_decoy_literature.tsv');lr=lr[(lr.dataset=='GSE314943')&(~lr.ligand_parts.eq('PROS1'))].copy();nn=lr.groupby(['sender','scope','group']).interaction.nunique().reset_index(name='eligible_interactions');fig,axs=plt.subplots(1,2,figsize=(11,4.5))
for ax,gr in zip(axs,['Saline','Silica']):
 z=nn[nn.group==gr];pv=z.pivot(index='sender',columns='scope',values='eligible_interactions').reindex(columns=['Naive','ActCD4','Treg','Core']);im=ax.imshow(pv,cmap='Blues',vmin=0,vmax=nn.eligible_interactions.max());ax.set_xticks(range(4),pv.columns);ax.set_yticks(range(len(pv)),pv.index);ax.set_title(gr+' | expression-compatible LR count')
 for (i,j),v in np.ndenumerate(pv.to_numpy()):ax.text(j,i,str(int(v)),ha='center',va='center',color='white' if v>nn.eligible_interactions.max()*.55 else '#183a54')
 fig.colorbar(im,ax=ax,shrink=.7)
save(fig,'F11','细胞通讯','巨噬来源与CD4状态存在表达相容的配体受体入口','计数为同一数据库相互作用ID的去重数，须满足冻结发送/接收覆盖规则；全部PROS1支持在此展示中撤回。AllMac与AMlike为重叠来源定义，不能相加。图中数量不是通讯强度，也不是功能激活或来源特异因果证据。',nn)
lig=rd(V/'tables/sender_ligand_expression.tsv');show=['CCL2','IL1A','IL1B','TNF','TGFB1','IL16'];show=[g for g in show if g in set(lig.ligand)];pv=lig[(lig.sender=='AMlike')&lig.ligand.isin(show)].pivot(index='ligand',columns='sample',values='CP10k').reindex(index=show,columns=sorted(lig['sample'].unique()));fig,ax=plt.subplots(figsize=(10,4.8));im=ax.imshow(np.log1p(pv),cmap='YlOrRd',aspect='auto');ax.set_xticks(range(len(pv.columns)),pv.columns,rotation=35,ha='right');ax.set_yticks(range(len(pv)),pv.index);ax.set_title('AMlike sender | frozen-route ligand RNA');fig.colorbar(im,ax=ax,label='log1p(pseudobulk CP10k)')
save(fig,'F12','细胞通讯','上游配体表达按样本及发送细胞范围展示','选择与优先候选讨论相关且位于冻结路线的配体。每列为独立样本，颜色仅代表RNA表达量；AMlike证据不与Broad Mac来源优势拼接。受体必要组分及来源竞争详见完整路线表。',lig[(lig.sender=='AMlike')&lig.ligand.isin(show)])
# 7 TF
auc=rd(V/'tables/clean_TF_sample_AUC.tsv');aa=auc[(auc.variant=='Raw')&auc.scope.isin(states)].copy();aa['column']=aa.scope+' / '+aa['sample'];pv=aa.pivot(index='tf',columns='column',values='AUC');pv=pv.reindex(columns=[st+' / '+s for st in states for s in samp]);z=(pv.sub(pv.mean(axis=1),axis=0)).div(pv.std(axis=1).replace(0,np.nan),axis=0).fillna(0);fig,ax=plt.subplots(figsize=(12,5));im=ax.imshow(z,cmap='RdBu_r',aspect='auto',vmin=-2.5,vmax=2.5);ax.set_xticks(range(len(z.columns)),[s.replace(' / ','\n') for s in z.columns],rotation=60,ha='right',fontsize=9);ax.set_yticks(range(len(z)),z.index);ax.set_title('Candidate/program-overlap-cleaned regulon AUCell');fig.colorbar(im,ax=ax,label='Row z score (display only)')
save(fig,'F13','转录因子网络','去重叠regulon活性刻画CD4状态调控环境','12个冻结TF的样本平均AUCell；成员去除候选/程序重叠规则沿用v18，图中行标准化仅用于展示。活性推断不等于TF蛋白活性实测；RNA与AUCell分开解释。',aa)
reg=rd(V/'tables/clean_regulon_members.tsv');reg=reg[reg.retained.astype(str).str.lower().eq('true')];top=reg.sort_values('weight',ascending=False).groupby('tf').head(2);tfs=sorted(top.tf.unique());targets=sorted(top.canonical.unique());fig,ax=plt.subplots(figsize=(12,7));yp={g:i/max(1,len(tfs)-1) for i,g in enumerate(tfs)};yt={g:i/max(1,len(targets)-1) for i,g in enumerate(targets)}
for r in top.itertuples():ax.plot([0,1],[yp[r.tf],yt[r.canonical]],color='#9aaebd',alpha=.6,lw=.8)
for g,y in yp.items():ax.scatter(0,y,s=120,color='#3278a8');ax.text(-.035,y,g,ha='right',va='center',fontsize=9)
for g,y in yt.items():ax.scatter(1,y,s=70,color='#d5a55c');ax.text(1.035,y,g,ha='left',va='center',fontsize=10)
ax.set_xlim(-.2,1.3);ax.set_ylim(-.06,1.06);ax.axis('off');ax.set_title('Representative retained TF-target network | top 2 weights per TF')
save(fig,'F14','转录因子网络','代表TF靶网络呈现调控先验的具体成员','每TF按原模型权重展示前2个保留靶基因，完整成员随统计表保留。连线为SCENIC推断/筛选的调控关系，不表示正负调控方向或已证实因果；该展示子集不用于重新排序候选。',top)
# 8 scores
scores=rd(V/'tables/sample_gene_and_program_scores.tsv');eff=rd(V/'tables/programme_group_effects.tsv')
for num,mods in [(15,['PyroExec','Inflam']),(16,['ApoExec','MitoISR'])]:
 fig,axs=plt.subplots(3,2,figsize=(10,9));data=[]
 for i,sc in enumerate(states):
  zz=scores[(scores.dataset=='GSE314943')&(scores.variant=='Raw')&(scores.scope==sc)]
  for j,mod in enumerate(mods):
   ax=axs[i,j];dots(ax,zz,'M_'+mod);r=eff[(eff.cohort.astype(str)=='314')&(eff.variant=='Raw')&(eff.scope==sc)&(eff.program==mod)].iloc[0];ax.set_title(f'{sc} | {mod}\ndelta={r.effect:+.3f}; P={r.P:.3g}; FDR={r.FDR:.3g}',fontsize=10);ax.set_ylabel('Mean gene log1p(CP10k)');data.append(zz[['sample','group','scope','M_'+mod]].assign(program=mod).rename(columns={'M_'+mod:'score'}))
 save(fig,f'F{num}','焦亡/凋亡评分',('焦亡执行与炎症程序分开评价' if num==15 else '凋亡执行与线粒体应激程序提供不同终点线索'),'采用v18样本汇总计数评分，9候选移除规则与完整成员覆盖沿用。每点独立样本，3+3；双侧组标签排列20种，FDR保留原52程序-状态检验家族。程序表达不等于死亡实测，Inflam、PyroExec、ApoExec与MitoISR不可互换。',pd.concat(data))
# 9 trajectory
pt=rd(TR/'pseudotime_cell_membership.tsv.gz').rename(columns={'index':'cell'});tm=rd(U/'trajectory_cells.tsv.gz');merged=tm[['cell','UMAP1','UMAP2']].merge(pt,on='cell',validate='one_to_one');fig,axs=plt.subplots(1,2,figsize=(11,4.6));im=axs[0].scatter(merged.UMAP1,merged.UMAP2,c=merged.K15_R1,s=6,cmap='viridis',rasterized=True);fig.colorbar(im,ax=axs[0],label='Frozen DPT K15_R1')
for st in ['Naive','ActCD4']:
 z=merged[merged.subtype==st];axs[1].scatter(z.UMAP1,z.UMAP2,c=statecolors[st],s=6,label=st,alpha=.6,rasterized=True)
for ax in axs:ax.set(xlabel='UMAP 1',ylabel='UMAP 2')
axs[0].set_title('Naive / ActCD4 state ordering');axs[1].set_title('Same cells, frozen state labels');axs[1].legend(frameon=False)
save(fig,'F17','状态轨迹','Naive与ActCD4可作状态排序，不能据此推断Treg谱系','复用1,734细胞的冻结DPT，不重选根。4套邻居/根设置的排序相关rho=0.904-0.999。该轨迹不包含Treg，也不等同真实时间；图不支持Naive必然转化为损伤细胞的因果判断。',merged)
cur=rd(TR/'program_sample_quantile_curves.tsv');tests=rd(TR/'program_trajectory_tests.tsv');sl=rd(TR/'sample_level_slopes.tsv.gz');fig,axs=plt.subplots(2,2,figsize=(11,7));chosen=[]
for i,st in enumerate(['Naive','ActCD4']):
 for j,mod in enumerate(['PyroExec','ApoExec']):
  ax=axs[i,j];z=cur[(cur.state==st)&(cur.module==mod)]
  for sample,b in z.groupby('sample'):ax.plot(b.mean_pt,b.score,'o-',c=colors['Saline' if sample.startswith('Saline') else 'Silica'],alpha=.65,lw=1,ms=4)
  r=tests[(tests.state==st)&(tests.run=='K15_R1')&(tests.feature==mod)].iloc[0];ax.set(title=f'{st} | {mod}\nslope={r.mean_slope:+.3f}; P={r.P:.3g}; FDR={r.BH:.3g}',xlabel='Within-state rescaled DPT',ylabel='Matched-control program score');chosen.append(z)
  ss=sl[(sl.state==st)&(sl.run=='K15_R1')&(sl.feature==mod)].slope.to_numpy();signs=np.array(list(itertools.product([-1,1],repeat=len(ss))));pv=np.mean(abs(signs@ss/len(ss))>=abs(ss.mean())-1e-12);assert abs(pv-r.P)<1e-12;qa.append(dict(check='trajectory_signflip',family=st+'|'+mod,error=abs(pv-r.P)))
save(fig,'F18','状态轨迹','样本级程序曲线显示状态排序中的表达变化','每条线是一个独立样本的4个分位区间；区间和细胞不是重复。P来自6个样本的文库深度调整斜率、64种符号翻转，FDR在同一排序的全部状态-程序家族内。此处为历史匹配对照评分，不能与图15-16汇总计数评分直接比值。',pd.concat(chosen))
# 10 NicheNet
ni=rd(B/'tables/NicheNet_activity_all.tsv');ni=ni[(ni.model=='joint_stratified')&(ni.direction=='up')&ni.scope.isin(states)&(ni.test_ligand!='PROS1')];toplig=ni.sort_values('activity_rank').groupby('scope').head(6).test_ligand.unique();pv=ni[ni.test_ligand.isin(toplig)].pivot_table(index='test_ligand',columns='scope',values='pearson').reindex(columns=states);fig,ax=plt.subplots(figsize=(10,6));im=ax.imshow(pv,cmap='YlGnBu',aspect='auto');ax.set_xticks(range(3),states);ax.set_yticks(range(len(pv)),pv.index);ax.set_title('NicheNet ligand activity | joint disease-up target sets');fig.colorbar(im,ax=ax,label='Ligand-target prior prediction Pearson')
save(fig,'F19','NicheNet','配体活性模型为巨噬来源信号提供排序线索','每状态取原活动排名前6个配体的并集展示，PROS1不作为当前功能支持。Pearson是NicheNet先验预测指标，不是独立样本相关P值；分状态背景与目标集不同，不能直接当作来源特异作用强度。空白为该状态无对应活动结果，不是数值0。',ni[ni.test_ligand.isin(toplig)])
routes=rd(B/'tables/candidate_complete_routes.tsv');rt=routes[(routes.ligand!='PROS1')&routes.gene.isin(['IL1R2','TTC39C','SKIL'])&routes.family.str.startswith('joint_stratified')].copy();rt=rt[((rt.gene=='TTC39C')&(rt.scope=='ActCD4'))|((rt.gene.isin(['IL1R2','SKIL']))&(rt.scope=='Treg'))];rt=rt.sort_values(['activity_rank','target_rank'],na_position='last').drop_duplicates(['gene','ligand','receptor']).groupby('gene').head(6);fig,ax=plt.subplots(figsize=(12,5.5));genes=['IL1R2','TTC39C','SKIL'];ligs=sorted(rt.ligand.unique());posg={g:i for i,g in enumerate(genes)};posl={g:i/max(1,len(ligs)-1)*2 for i,g in enumerate(ligs)}
for r in rt.itertuples():ax.plot([0,1],[posl[r.ligand],posg[r.gene]],c='#9caeb8',lw=1,alpha=.65)
for g,y in posl.items():ax.scatter(0,y,c='#3278a8',s=90);ax.text(-.03,y,g,ha='right',va='center',fontsize=9)
for g,y in posg.items():ax.scatter(1,y,c='#d46746',s=180);ax.text(1.04,y,g,ha='left',va='center',fontsize=11)
ax.set_xlim(-.25,1.4);ax.set_ylim(-.2,2.2);ax.axis('off');ax.set_title('Prior-supported candidate routes (representative subset; no PROS1)')
save(fig,'F20','NicheNet','候选与上游配体的模型连接用于安排后续验证','展示当前联合DE对应的既有完整路线子集，每候选至多6条，按原活动/靶排序选取。完整受体组分保留在源表；IL1R2相关诱饵受体不等于胞内信号传导。若某候选缺少路线，空缺保持，不补造连接。本图TTC39C在当前联合名义DE对应路线中无连接，原因包括联合DE未达P<0.05，不表示所有历史先验均不存在。',rt)
# 11 current candidates
exp=rd(T/'candidate_sample_detection.tsv');fig,axs=plt.subplots(1,3,figsize=(12,4.4));chosen=[]
for ax,(g,sc) in zip(axs,[('IL1R2','Treg'),('TTC39C','ActCD4'),('SKIL','Treg')]):
 z=exp[(exp.dataset=='GSE314943')&(exp.variant=='Raw')&(exp.gene==g)&(exp.scope==sc)];dots(ax,z,'positive_fraction');ax.set(title=g+' / '+sc,ylabel='RNA-positive cell fraction');ax.set_ylim(0,max(.35,z.positive_fraction.max()*1.2));chosen.append(z)
save(fig,'F21','候选与验证','优先候选在目标亚群内具有可追溯的样本检出','每点为独立样本阳性率；IL1R2 5/6样本检出，TTC39C/SKIL 6/6检出。均需进一步验证蛋白和干预效率；RNA阳性率不是预实验成功率。',pd.concat(chosen))
fig,axs=plt.subplots(1,3,figsize=(12,4.4));assoc=rd(V/'tables/all_associations.tsv');chosen=[]
for ax,(g,sc,mod) in zip(axs,[('IL1R2','Treg','Inhib'),('TTC39C','ActCD4','ApoExec'),('SKIL','Treg','TAct')]):
 z=scores[scores.dataset.isin(['GSE314943','GSE183682'])&(scores.variant=='Raw')&(scores.scope==sc)]
 for ds,m in [('GSE314943','o'),('GSE183682','^')]:
  zz=z[z.dataset==ds];ax.scatter(zz['G_'+g],zz['M_'+mod],c=zz.group.map(colors),marker=m,s=45,label=ds)
 r=assoc[(assoc.cohort=='Joint')&(assoc.variant=='Raw')&(assoc.scope==sc)&(assoc.model=='JointStratum')&(assoc.x=='G_'+g)&(assoc.y=='M_'+mod)].iloc[0];ax.set(title=f'{g} / {sc}\nadjusted rho={r.rho:.3f}; P={r.P:.3g}\nFDR={r.FDR:.3g}',xlabel=g+' log1p(CP10k)',ylabel=mod+' score');chosen.append(z)
axs[0].legend(fontsize=7,frameon=False)
save(fig,'F22','候选与验证','IL1R2、TTC39C和SKIL分别对应状态、凋亡及活化问题','点展示原始样本评分；图题rho/P来自控制队列/时间层的秩关联模型，并非图中原始点的简单相关。圆为主队列，三角为第二队列；颜色与前图一致。三条路线的主要实验终点不同，不以一个总分混排。',pd.concat(chosen))
pd.DataFrame(qa).to_csv(T/'v20_numeric_QA.tsv',sep='\t',index=False);(O/'figure_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8');print('FIGURES',len(manifest),'QA',len(qa));print(cc.to_string(index=False))