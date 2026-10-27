from pathlib import Path
import pandas as pd,numpy as np,json,math
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/three_cohort_four_axis_20261017_v16'
alias={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
def c(x):x=str(x).upper();return alias.get(x,x)
g=pd.read_csv(O/'resources/GO_mouse_membership.tsv',sep='\t').dropna(subset=['SYMBOL']);rows=[]
for r in g.itertuples():rows.append(dict(gene=c(r.SYMBOL),source='GO_mouse',term_id=r.GOALL,term=r.TERM,evidence=r.EVIDENCEALL,mechanism='pyroptosis' if 'pyropt' in r.TERM else 'intrinsic_apoptosis',role='negative_regulation' if 'negative regulation' in r.TERM else 'positive_regulation' if 'positive regulation' in r.TERM else 'regulation' if 'regulation' in r.TERM else 'pathway_member',context='mouse_annotation_not_CD4_specific'))
for pid,mech in [('R-MMU-109606','intrinsic_apoptosis'),('R-MMU-5620971','pyroptosis')]:
 for a in json.loads((O/'resources'/f'{pid}.json').read_text()):
  if a.get('geneName'):rows.append(dict(gene=c(a['geneName'][0]),source='Reactome_mouse',term_id=pid,term=mech,evidence=a.get('identifier'),mechanism=mech,role='pathway_participant',context='mouse_pathway_not_CD4_specific'))
f=pd.DataFrame(rows).drop_duplicates();f.to_csv(O/'tables/function_membership.tsv',sep='\t',index=False);print('FUNCTION_GENES',f.gene.nunique())
db=pd.read_csv(O/'resources/CellChat_interaction.tsv',sep='\t',keep_default_na=False);cx=pd.read_csv(O/'resources/CellChat_complex.tsv',sep='\t',keep_default_na=False).set_index('complex')
def parts(v):return [c(x) for x in cx.loc[v].tolist() if x] if v in cx.index else [c(v)]
db=db[db.annotation.isin(['Secreted Signaling','ECM-Receptor','Cell-Cell Contact'])].copy();db['ligand_parts']=db.ligand.map(lambda x:';'.join(parts(x)));db['receptor_parts']=db.receptor.map(lambda x:';'.join(parts(x)));db.to_csv(O/'resources/CellChat_protein_LR_expanded.tsv',sep='\t',index=False)
def load(ds,sc):
 stem=ds+'__'+sc;m=pd.read_csv(O/'inputs'/f'{stem}.samples.tsv',sep='\t');p=pd.read_csv(O/'inputs'/f'{stem}.fraction.tsv.gz',sep='\t',index_col=0);return m,p
coverage=[];eligible={}
for ds,sc in [('GSE314943',s) for s in ['AllMac','AMlike','Core','Naive','ActCD4','Treg']]+[('GSE274372',s) for s in ['AllMac','TLike']]:
 m,p=load(ds,sc)
 for group,mm in m.groupby('group'):
  ss=mm.loc[mm.cells>=10,'sample'].tolist();need=math.ceil(.5*len(ss));ok=(p[ss]>=.05).sum(1)>=need if ss else pd.Series(False,index=p.index);eligible[ds,sc,group]=set(p.index[ok]);
  coverage.append(dict(dataset=ds,scope=sc,group=group,total_samples=len(mm),eligible_samples=len(ss),required_samples=need,expressed_genes=int(ok.sum())))
pd.DataFrame(coverage).to_csv(O/'tables/communication_sample_coverage.tsv',sep='\t',index=False)
lr=[];bg=[]
for ds,sc in [('GSE314943',s) for s in ['Core','Naive','ActCD4','Treg']]+[('GSE274372','TLike')]:
 groups=sorted(k[2] for k in eligible if k[:2]==(ds,sc));background=set().union(*(eligible[ds,sc,gr] for gr in groups));bg.extend(dict(dataset=ds,scope=sc,gene=v) for v in sorted(background))
 for sender in (['AllMac','AMlike'] if ds=='GSE314943' else ['AllMac']):
  for gr in groups:
   se=eligible[ds,sender,gr];re=eligible[ds,sc,gr]
   for r in db.itertuples(index=False):
    if set(r.ligand_parts.split(';'))<=se and set(r.receptor_parts.split(';'))<=re:
     lr.append(dict(dataset=ds,scope=sc,sender=sender,group=gr,interaction=r.interaction_name,ligand=r.ligand,ligand_parts=r.ligand_parts,receptor=r.receptor,receptor_parts=r.receptor_parts,annotation=r.annotation,evidence=r.evidence))
pd.DataFrame(lr).to_csv(O/'tables/expressed_curated_LR.tsv',sep='\t',index=False);pd.DataFrame(bg).to_csv(O/'inputs/receiver_background.tsv',sep='\t',index=False)
print('LR',len(lr));print(pd.DataFrame(lr).groupby(['dataset','scope']).size().to_string())