from pathlib import Path
import json,sys
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.io import mmread
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
O=Path('/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/18_macrophage_supplement_20261019')
V=P/'11_reports/death_priority_20261018_v2'
sys.path.insert(0,str(V/'code'));import statslib as st
st.O=O
def wr(x,n):pd.DataFrame(x).to_csv(O/'tables'/n,sep='\t',index=False)
panels={
 'Inflammatory':['IL1B','TNF','NFKBIA','PTGS2','CCL3','CCL4'],
 'IFN_response':['STAT1','IRF1','ISG15','IFIT1','IFIT3','CXCL10'],
 'Antigen_presentation':['CD74','H2-AA','H2-AB1','H2-EB1','CIITA'],
 'Alveolar_lipid':['FABP4','PPARG','MARCO','ABCG1','LPL'],
 'SPP1_remodeling':['SPP1','GPNMB','LGALS3','APOE','CTSB','CTSD'],
 'C1Q_phagocytosis':['C1QA','C1QB','C1QC','MERTK','MSR1']}
humanpanels={**panels,'Antigen_presentation':['CD74','HLA-DRA','HLA-DRB1','HLA-DPA1','CIITA']}
ligands=['IL1B','TNF','IL6','IL18','TGFB1','SPP1','CXCL16','LGALS9','F11R','SIRPB1C','IGF1']
endpoints=['TAct','TReg','Inhib','Naive','Cycle','MitoROS','MitoISR','ApoExec','PyroExec']
protocol=dict(date='2026-10-19',mouse_panels=panels,human_panels=humanpanels,ligands=ligands,endpoints=endpoints,
 panel_status='Prespecified analyst-defined descriptive panels; not validated subtype classifiers. No M1/M2 assignment.',
 score='mean log1p donor pseudobulk CP10k; full library normalization; >=3 measured panel genes; >=10 cells per donor/source',
 association='rank partial correlation; permutations within dataset/day; sensitivity adds condition and permutes within dataset/day/condition; leave-one-donor sign',
 families='separate cohort/sender/receiver variant/model/state-or-ligand; all three scopes and nine endpoints included in each family',
 human='paired donors; primary same-source G227; exploratory pooled source-adjusted and condition-adjusted; retain CWP label',
 nichenet='disease nominal P<0.05 up/down intersect frozen T-cell modules; background all eligible receiver genes; >=5 response genes; activity has no P')
(O/'protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
aliases={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
def canon(x):return np.array([aliases.get(str(g).upper(),str(g).upper()) for g in x])
coverage=[];members=[]
def scores(x,g,m,cohort,panel):
 g=canon(g);u,inv=np.unique(g,return_inverse=True)
 x=(x@sparse.csr_matrix((np.ones(len(g)),(np.arange(len(g)),inv)),shape=(len(g),len(u)))).tocsr();g=u
 assert x.shape[0]==len(m) and not m.cell.duplicated().any()
 assert np.all(np.asarray(x.sum(1)).ravel()<=m.library.to_numpy()+1e-5)
 print('INPUT',cohort,x.shape,'full_library_ratio',float(np.median(np.asarray(x.sum(1)).ravel()/m.library.to_numpy())),flush=True)
 strict=np.asarray((x[:,np.isin(g,['C1QA','C1QB','C1QC','ADGRE1','MERTK','MSR1','MARCO'])]>0).sum(1)).ravel()>=2
 out=[]
 for source,keep in [('AllMac',np.ones(len(m),bool)),('MacStrict',strict)]:
  for sam,idx0 in m.groupby('sample').groups.items():
   idx=np.array(list(idx0));idx=idx[keep[idx]];z=m.loc[idx];r=dict(cohort=cohort,source=source,sample=sam,cells=len(idx))
   if not len(idx):coverage.append(r);continue
   r.update({c:z[c].iloc[0] for c in ['dataset','day','group']});coverage.append(r.copy())
   if len(idx)<10:continue
   lib=float(z.library.sum());cp=np.asarray(x[idx].sum(0)).ravel()/lib*1e4
   for name,gs in panel.items():
    avail=sorted(set(gs)&set(g));r['S_'+name]=float(np.log1p(cp[np.isin(g,avail)]).mean()) if len(avail)>=3 else np.nan
    if sam==m['sample'].iloc[0]:members.extend(dict(cohort=cohort,source=source,panel=name,gene=a,measured=a in g) for a in gs)
   for a in ligands:
    j=np.flatnonzero(g==a);r['L_'+a]=float(np.log1p(cp[j[0]])) if len(j) else np.nan
   out.append(r)
 return pd.DataFrame(out)
M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF/comm_input'
m=pd.read_csv(M/'Mac_cells.tsv',sep='\t');m['dataset']='GSE314943';m['day']=140
x=mmread(M/'Mac_counts.mtx.gz').T.tocsr();g=pd.read_csv(M/'Mac_genes.tsv',sep='\t').canonical
mouse=[scores(x,g,m,'Mouse',panels)];del x
a=ad.read_h5ad(V/'objects/GSE183682_whole_clustered.h5ad');ann=json.loads((V/'inputs/GSE183682_annotation.json').read_text())
print('ANNOTATION',ann,flush=True)
keep=a.obs.cluster.astype(str).map(ann).eq('Mac').to_numpy()
m=a.obs.loc[keep].copy().reset_index(names='cell');m['dataset']='GSE183682'
x=a.layers['counts'][keep].tocsr();mouse.append(scores(x,a.var_names,m,'Mouse',panels));del a,x
ms=pd.concat(mouse,ignore_index=True);wr(ms,'mouse_mac_scores.tsv')
tc=pd.read_csv(V/'tables/sample_scores.tsv',sep='\t')
lr=pd.read_csv(V/'tables/joint_LR_eligibility.tsv',sep='\t');eligible=lr[lr.eligible.astype(str).str.lower().eq('true')]
for source,mac in ms.groupby('source'):
 for col in [c for c in mac if c.startswith(('S_','L_'))]:
  st.group_effect(mac,col,'Joint',source,'Mac')
 for variant in ['Raw','BG1']:
  for scope in ['Naive','ActCD4','Treg']:
   t=tc[(tc.variant==variant)&(tc.scope==scope)]
   d=t.merge(mac.drop(columns=['cohort','cells','dataset','day','group','source']),on='sample',validate='one_to_one')
   assert len(d)==10
   for model in ['JointStratum','JointStratumCondition']:
    for kind,prefix in [('MacState','S_'),('MacLigand','L_')]:
     for col in [c for c in mac if c.startswith(prefix)]:
      for ep in endpoints:
       er=eligible[(eligible.source==source)&(eligible.scope==scope)&eligible.ligand_parts.eq(col[2:])]
       st.assoc(d,col,'M_'+ep,'Mouse',variant,kind,scope,model,source,extra=dict(LR_eligible=bool(len(er)) if kind=='MacLigand' else False,LR_interactions=';'.join(er.interaction) if kind=='MacLigand' else '',family=f'Mouse|{source}|{variant}|{model}|{kind}'))
wr(coverage,'mac_coverage.tsv');wr(members,'mac_panel_members.tsv')
st.finalize(st.results,'mouse_associations.tsv');st.finalize(st.effectrows,'mouse_mac_disease_effects.tsv')
print('MOUSE_DONE',len(st.results),flush=True)
H=P/'02_analysis/scenic3_20261001_123452_Kvp9h5/input'
m=pd.read_csv(H/'human_mac_cells.tsv',sep='\t');m['dataset']=m.series;m['day']=0
x=mmread(H/'human_mac_counts.mtx').T.tocsr();g=pd.read_csv(H/'human_mac_genes.tsv',sep='\t').gene
hs=scores(x,g,m,'Human',humanpanels);wr(hs,'human_mac_scores.tsv');del x
ht=pd.read_csv(P/'11_reports/death_priority_20261018_v3/tables/human_sample_scores.tsv',sep='\t')
st.results=[];st.testids=set();pair=[]
for cohort in ['G227','AllHuman']:
 for source,mac in hs.groupby('source'):
  for sub in ['CD4-CCR7','CD4-GZMA','Treg']:
   t=ht[ht.subtype.eq(sub)].copy()
   if cohort=='G227':t=t[t.series=='G227']
   d=t.merge(mac.drop(columns=['cohort','cells','day','group','source']),on='sample',validate='one_to_one');d['day']=0
   assert (d.series==d.dataset).all()
   d['group_original']=d.group;d['group']=d.group.map({'CWP':'Silica','Control':'Saline'})
   assert d.group.notna().all()
   pair.extend(dict(cohort=cohort,source=source,scope=sub,sample=r['sample'],dataset=r['dataset'],group=r['group_original'],T_cells=r['cells']) for r in d.to_dict('records'))
   for model in ['JointStratum','JointStratumCondition']:
    for kind,prefix in [('MacState','S_'),('MacLigand','L_')]:
     for col in [c for c in mac if c.startswith(prefix) and c!='L_SIRPB1C']:
      for ep in endpoints:
       st.assoc(d,col,ep,cohort,'HumanFrozen',kind,sub,model,source,extra=dict(family=f'{cohort}|{source}|{model}|{kind}',group_mapping='CWP case indicator; not silicosis diagnosis'))
st.finalize(st.results,'human_associations.tsv');wr(pair,'human_pairing.tsv');wr(coverage,'mac_coverage.tsv');wr(members,'mac_panel_members.tsv')
(O/'logs/analysis_complete.json').write_text(json.dumps(dict(status='success',mouse_rows=len(ms),human_rows=len(hs),human_tests=len(st.results))),encoding='utf-8')
print('ALL_COMPLETE',len(st.results),flush=True)
