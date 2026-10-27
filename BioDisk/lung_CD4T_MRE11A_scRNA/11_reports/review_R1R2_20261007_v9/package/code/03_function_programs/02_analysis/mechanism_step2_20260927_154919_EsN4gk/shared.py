import os, re, json, gzip, tarfile, hashlib, platform, itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata
from scipy.io import mmwrite

TARGETS = ['IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3','MRE11']
SETS = {
 'OXPHOS':'NDUFA1 NDUFA2 NDUFA9 NDUFB8 NDUFS1 NDUFS2 NDUFS3 NDUFV1 SDHA SDHB SDHC UQCRC1 UQCRC2 UQCRFS1 CYC1 COX4I1 COX5A COX5B COX6A1 ATP5F1A ATP5F1B ATP5F1C ATP5MC1 ATP5PO',
 'MitoISR':'ATF4 ATF5 DDIT3 DDIT4 CHAC1 ASNS TRIB3 PPP1R15A SESN2 DELE1 EIF2AK1 HSPA9 HSPD1 HSPE1 LONP1 CLPP OMA1',
 'MitoROS':'SOD2 PRDX3 PRDX5 GPX4 TXN2 TXNRD2 HMOX1 NQO1 GCLC GCLM SRXN1 FTH1',
 'Mitophagy':'PINK1 PRKN BNIP3 BNIP3L FUNDC1 SQSTM1 OPTN TBK1 ULK1 BECN1 MAP1LC3B',
 'Inflam':'NLRP3 AIM2 NLRC4 PYCARD CASP1',
 'PyroExec':'PYCARD CASP1 GSDMD IL1B IL18',
 'ApoExec':'BAX BAK1 APAF1 CASP9 CASP3 CASP7 CYCS',
 'TAct':'CD69 IL2RA ICOS TNFRSF4 BATF IRF4 NR4A1 NR4A2 EGR1',
 'TReg':'FOXP3 IL2RA CTLA4 IKZF2 IKZF4 LRRC32 ENTPD1 TIGIT',
 'Inhib':'PDCD1 CTLA4 TIGIT HAVCR2 TOX TOX2 ENTPD1',
 'Naive':'CCR7 SELL TCF7 LEF1 KLF2 IL7R MAL',
 'Cycle':'MKI67 TOP2A PCNA MCM2 MCM3 MCM4 MCM5 MCM6 MCM7 UBE2C'
}
SETS={k:v.split() for k,v in SETS.items()}
EPI='SFTPA1 SFTPA2 SFTPB SFTPC SFTPD SCGB1A1 SCGB3A1 SCGB3A2 ABCA3 SLC34A2 NAPSA'.split()
MYE='LYZ LYZ1 LYZ2 S100A8 S100A9 CSF1R LST1 TYROBP FCER1G'.split()
SOURCE_GENES='IL1A IL1B IL18 CXCL16 TNF IL6 IL10 TGFB1 IL12A IL12B IL23A CCL3 CCL4 CCL5 SPP1 LGALS9 CD80 CD86 CD274 PDCD1LG2 TNFSF9'.split()
RECEPTORS='IL1R1 IL1RAP IL1R2 IL18R1 IL18RAP CXCR6 TNFRSF1A TNFRSF1B IL6R IL6ST IL10RA IL10RB TGFBR1 TGFBR2 IL12RB1 IL12RB2 IL23R CCR5 CD44 ITGAV ITGB1 ITGB3 HAVCR2 CD28 CTLA4 PDCD1 TNFRSF9'.split()
SPECS=[
 ('GSE314943','CD4','GSE314943/stage5_20260905_222457_a2dApT/objects/CD4_core.h5ad',2311),
 ('GSE314943','Mac','GSE314943/stage11_20260906_100314_pALwAB/objects/macrophages_core.h5ad',1796),
 ('GSE183682','CD4','screen2_CD4_target_wz1lfx/GSE183682/CD4_core.h5ad',1413),
 ('GSE274372','TLike','GSE274372/stage9_20260908_140721_Yj56oq/objects/CD4_author_reclustered.h5ad',888),
 ('GSE174725','CD4','screen2_Trefine_nuzVQu/GSE174725/T_refined.h5ad',10286)]

def canonical(s):
    s=str(s).strip().upper()
    return {'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}.get(s,s)

def rho(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float);ok=np.isfinite(a)&np.isfinite(b);a=a[ok];b=b[ok]
    if len(a)<3 or np.ptp(a)==0 or np.ptp(b)==0:return np.nan
    return float(np.clip(np.corrcoef(rankdata(a),rankdata(b))[0,1],-1,1))

def partial(a,b,z):
    a=np.asarray(a,float);b=np.asarray(b,float);z=np.asarray(z,float)
    if z.ndim==1:z=z[:,None]
    ok=np.isfinite(a)&np.isfinite(b)&np.isfinite(z).all(axis=1);a=a[ok];b=b[ok];z=z[ok]
    if len(a)<6 or np.ptp(a)==0 or np.ptp(b)==0:return np.nan
    d=np.column_stack([np.ones(len(a))]+[rankdata(z[:,j]) for j in range(z.shape[1])])
    if len(a)-np.linalg.matrix_rank(d)<3:return np.nan
    u=rankdata(a);v=rankdata(b);u-=d@np.linalg.lstsq(d,u,rcond=None)[0];v-=d@np.linalg.lstsq(d,v,rcond=None)[0]
    if np.std(u)<1e-10 or np.std(v)<1e-10:return np.nan
    return float(np.clip(np.corrcoef(u,v)[0,1],-1,1))

def effect(b,c):
    b=np.asarray(b,float);c=np.asarray(c,float);b=b[np.isfinite(b)];c=c[np.isfinite(c)]
    d=float(c.mean()-b.mean()) if len(b) and len(c) else np.nan
    loo=[]
    if min(len(b),len(c))>=2:
        loo=[c.mean()-np.delete(b,i).mean() for i in range(len(b))]+[np.delete(c,i).mean()-b.mean() for i in range(len(c))]
    return dict(n_control=len(b),n_case=len(c),control_mean=b.mean() if len(b) else np.nan,case_mean=c.mean() if len(c) else np.nan,
      difference=d,loo_same=np.mean(np.sign(loo)==np.sign(d)) if loo and d!=0 else np.nan,
      loo_min=min(loo) if loo else np.nan,loo_max=max(loo) if loo else np.nan)

def bgmask(symbols):
    return np.array([bool(re.match(r'^(IGH[VDJMAEG]|IGK[VJC]|IGL[VJC])',s)) or s in EPI or s in MYE for s in symbols])

def scores(y, x, symbols, samples):
    # Equal-sample reference expression; controls fixed before group comparisons.
    reference=np.mean([np.asarray(y[np.asarray(samples)==s].mean(axis=0)).ravel() for s in sorted(set(samples))],axis=0)
    excluded=set(TARGETS)|set(itertools.chain.from_iterable(SETS.values()))|set(EPI)|set(MYE)
    pool=np.array([s not in excluded and not re.match(r'^(MT-|RPL|RPS|MRPL|MRPS|IGH|IGK|IGL)',s) for s in symbols])&(reference>0)
    edges=np.unique(np.quantile(reference[reference>0],np.linspace(0,1,26)))
    bins=np.digitize(reference,edges[1:-1]);coverage=[];controls=[];values={};members={}
    for name,expected in SETS.items():
        j=np.flatnonzero(np.isin(symbols,expected));observed=j[np.asarray(x[:,j].sum(axis=0)).ravel()>0]
        present=list(symbols[j]);eligible=len(j)>=3 and len(j)/len(expected)>=.6 and len(observed)>=3
        selected=set()
        for ix in j:
            same=np.flatnonzero(pool & (bins==bins[ix]))
            options=same if len(same)>=20 else np.flatnonzero(pool)
            order=np.argsort(np.abs(reference[options]-reference[ix]),kind='stable')
            selected.update(options[order[:20]].tolist())
        cc=np.array(sorted(selected),dtype=int)
        eligible=eligible and len(cc)>=20
        coverage.append(dict(module=name,expected=len(expected),present=len(j),expressed=len(observed),coverage=len(j)/len(expected),eligible=eligible,
          genes=';'.join(present),missing=';'.join(sorted(set(expected)-set(present))),controls=len(cc)))
        for ix in cc:controls.append(dict(module=name,gene=symbols[ix],reference_mean=reference[ix],bin=int(bins[ix])))
        mean=np.asarray(y[:,j].mean(axis=1)).ravel() if len(j) else np.full(y.shape[0],np.nan)
        values[name]=mean-np.asarray(y[:,cc].mean(axis=1)).ravel() if eligible else np.full(y.shape[0],np.nan)
        members[name]=np.asarray((x[:,j]>0).sum(axis=1)).ravel() if len(j) else np.zeros(y.shape[0])
    return values,members,coverage,controls


