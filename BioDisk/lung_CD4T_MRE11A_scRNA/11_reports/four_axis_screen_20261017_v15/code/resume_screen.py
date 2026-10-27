import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2')
from pathlib import Path
import sys,json,hashlib,itertools,re,subprocess,time,shutil
import numpy as np,pandas as pd
from scipy.io import mmread
from scipy.stats import rankdata,false_discovery_control
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/four_axis_screen_20261017_v15';V=P/'11_reports/integrated_macrophage_CD4_20261015_v4';M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF';BASE=P/'11_reports/methods_complete_20261016_v13';S=['Saline_1','Saline_2','Saline_3','Silica_1','Silica_2','Silica_3'];SC=['Naive','ActCD4','Treg','Core'];PR=['PyroExec','MitoROS','MitoISR'];CTX=['Inflam','ApoExec','OXPHOS','Mitophagy'];manifest=[];ALIASES={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
def rd(f):
 f=Path(f);manifest.append(dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()));return pd.read_csv(f,sep='\t')
def wr(a,n):a.to_csv(O/n,sep='\t',index=False,encoding='utf-8',na_rep='NA')
def bh(p):
 p=np.asarray(p,float);q=np.full(len(p),np.nan);ok=np.isfinite(p)
 if ok.any():q[ok]=false_discovery_control(p[ok])
 return q
PERMS=np.array(list(itertools.permutations(range(6))));WITHIN=np.array([a+b for a in itertools.permutations(range(3)) for b in itertools.permutations(range(3,6))])
def ranks(a,adj,grp=np.array([0,0,0,1,1,1])):
 z=rankdata(a,axis=0).astype(float)
 for g in (np.unique(grp) if adj else [None]):
  ix=np.ones(len(grp),bool) if g is None else grp==g;z[ix]-=z[ix].mean(axis=0)
 return z
# paired columns permit a different leave-one-candidate-out module for each candidate.
def corr_pairs(X,Y,adj,exact=True):
 if X.ndim==1:X=X[:,None]
 if Y.ndim==1:Y=Y[:,None]
 if X.shape[1]==1 and Y.shape[1]>1:X=np.broadcast_to(X,Y.shape)
 if Y.shape[1]==1 and X.shape[1]>1:Y=np.broadcast_to(Y,X.shape)
 rx=ranks(X,adj);ry=ranks(Y,adj);den=np.sqrt((rx*rx).sum(0)*(ry*ry).sum(0));ok=np.isfinite(X).all(0)&np.isfinite(Y).all(0)&(den>1e-12);r=np.divide((rx*ry).sum(0),den,out=np.full(len(den),np.nan),where=ok);pv=np.full(len(r),np.nan)
 if exact:
  perms=WITHIN if adj else PERMS
  for start in range(0,len(r),512):
   end=min(start+512,len(r));null=np.einsum('psg,sg->pg',rx[perms,start:end],ry[:,start:end],optimize=True)/np.where(den[start:end]>1e-12,den[start:end],np.nan);pp=(np.abs(null)>=np.abs(r[start:end])[None,:]-1e-12).mean(0);pv[start:end]=np.where(ok[start:end],pp,np.nan)
 loo=[]
 for k in range(6):
  keep=np.arange(6)!=k;xx=ranks(X[keep],adj,np.array([0,0,0,1,1,1])[keep]);yy=ranks(Y[keep],adj,np.array([0,0,0,1,1,1])[keep]);dd=np.sqrt((xx*xx).sum(0)*(yy*yy).sum(0));rr=np.divide((xx*yy).sum(0),dd,out=np.full(len(dd),np.nan),where=dd>1e-12);loo.append(rr)
 loo=np.array(loo);valid=np.isfinite(loo);agree=((np.sign(loo)==np.sign(r)[None,:])&valid).sum(0);return r,pv,agree,valid.sum(0)

print('RESUME_FROM_SAVED_QUALITY_DAMAGE',flush=True)
gn=rd(M/'comm_input/CD4_genes.tsv');genes=pd.Index(gn.canonical);quality=rd(O/'01_quality/all_genes_quality.tsv.gz');deall=rd(O/'02_DE/all_frozen_edgeR_results.tsv.gz');damage=rd(O/'03_damage/all_gene_program_tests.tsv.gz');agg={}
for sc in SC:
 a=rd(O/f'01_quality/{sc}_all_gene_sample_values.tsv.gz')
 for branch in ['All','BG1']:
  b=a[a.branch==branch].pivot(index='sample',columns='gene',values='CP10k').reindex(index=S,columns=genes);agg[(sc,branch)]=b.to_numpy()
for f in [M/'comm_input/CD4_counts.mtx.gz',M/'comm_input/CD4_cells.tsv']:
 manifest.append(dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))

r='''args<-commandArgs(TRUE);O<-args[1];P<-dirname(dirname(O));canon<-function(x){x<-toupper(trimws(x));a<-c(MRE11A='MRE11',ATP5A1='ATP5F1A',ATP5B='ATP5F1B',ATP5C1='ATP5F1C',ATP5G1='ATP5MC1',ATP5O='ATP5PO',PARK2='PRKN');j<-x%in%names(a);x[j]<-a[x[j]];x};lt<-readRDS(file.path(P,'00_raw/nichenet_v2_mouse/ligand_target_matrix_nsga2r_final_mouse.rds'));rownames(lt)<-canon(rownames(lt));colnames(lt)<-canon(colnames(lt));stopifnot(!anyDuplicated(rownames(lt)),!anyDuplicated(colnames(lt)));res<-lapply(seq_len(ncol(lt)),function(i){w<-lt[,i];ord<-order(-w,rownames(lt));ord<-head(ord[w[ord]>0],200);if(length(ord)==0)return(NULL);data.frame(ligand=colnames(lt)[i],gene=rownames(lt)[ord],weight=w[ord],prior_rank=seq_along(ord))});write.table(do.call(rbind,res),file.path(O,'inputs/top200_ligand_target_prior.tsv'),sep="\\t",quote=FALSE,row.names=FALSE);writeLines(capture.output(sessionInfo()),file.path(O,'provenance/prior_R_session.txt'))'''
(O/'code/export_prior.R').write_text(r,encoding='utf-8');res=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(O/'code/export_prior.R'),str(O)],capture_output=True,text=True);(O/'provenance/prior_export.log').write_text(res.stdout+res.stderr,encoding='utf-8');assert res.returncode==0,res.stderr
prior=rd(O/'inputs/top200_ligand_target_prior.tsv');lr=rd(V/'inputs/LR_prior.tsv');send=rd(V/'02_sender/all_ligand_sample_expression.tsv.gz');rec=rd(V/'03_receiver/receptor_candidate_sample_expression.tsv.gz');elig=rd(V/'02_sender/ligand_scope_eligibility.tsv');macs=[];edgesupport=[]
for source in ['AMlike','AllMac']:
 sd=rd(V/f'02_sender/{source}_DE_all_genes.tsv.gz');assert sd.canonical.is_unique;sd=sd.set_index('canonical');sender=send[send.scope==source]
 for sc in SC:
  de=deall[deall.scope==sc].set_index('canonical');gg=set(de.index)&set(genes);R=rec[(rec.scope==sc)&(rec.group=='Silica')].groupby('gene').UMI.apply(lambda a:int((a>0).sum()));L=sender[sender.group=='Silica'].groupby('gene').UMI.apply(lambda a:int((a>0).sum()));frames=[]
  for lig,b in prior.groupby('ligand',sort=True):
   if L.get(lig,0)<2:continue
   receptors=sorted(set(lr.loc[lr.ligand==lig,'receptor']));detrec=[r for r in receptors if R.get(r,0)>=2]
   if not detrec:continue
   b=b[b.gene.isin(gg)].copy()
   if b.empty:continue
   receptor_sig=[r for r in detrec if r in de.index and de.loc[r,'PValue']<.05];lp=float(sd.loc[lig,'PValue']) if lig in sd.index else np.nan;lf=float(sd.loc[lig,'logFC']) if lig in sd.index else np.nan;lb=float(sd.loc[lig,'FDR']) if lig in sd.index else np.nan;support=(lp<.05) or bool(receptor_sig)
   G=b.gene.to_numpy();ix=genes.get_indexer(G);xx=sender[sender.gene==lig].set_index('sample').reindex(S).CP10k.to_numpy();Y=agg[(sc,'All')][:,ix];Ybg=agg[(sc,'BG1')][:,ix]
   for adj in [False,True]:
    r,p,lc,lv=corr_pairs(xx,Y,adj);rb,pb,lcb,lvb=corr_pairs(xx,Ybg,adj,exact=False);f=b.copy();f['scope']=sc;f['sender']=source;f['model']='Adjusted' if adj else 'Pooled';f['rho']=r;f['P']=p;f['LOO_agree']=lc;f['LOO_valid']=lv;f['BG1_rho']=rb;f['ligand_DE_P']=lp;f['ligand_DE_BH']=lb;f['ligand_logFC']=lf;f['detected_receptors']=';'.join(detrec);f['DE_receptors']=';'.join(receptor_sig);f['disease_support']=support;f['coherent']=np.sign(lf)*np.sign(r)==np.sign(de.loc[G,'logFC'].to_numpy()) if lp<.05 else bool(receptor_sig);f['n_samples']=6;frames.append(f)
  a=pd.concat(frames,ignore_index=True) if frames else pd.DataFrame();a['family']=source+'_'+sc+'_'+a.model+'_all_prior_links';a['BH']=a.groupby('family').P.transform(bh);macs.append(a);wr(a,'04_macrophage/'+source+'_'+sc+'_all_links.tsv.gz');print('MAC',source,sc,len(a),flush=True)
mac=pd.concat(macs,ignore_index=True)
# Consolidate each evidence axis before intersection. No cross-state/source stitching.
d=damage[(damage.branch=='All')&(damage.program_family=='primary')&(damage.model=='Pooled')].copy();adj=damage[(damage.branch=='All')&(damage.program_family=='primary')&(damage.model=='Adjusted')][['scope','gene','module','rho','P','BH','LOO_agree','LOO_valid']].rename(columns={c:'adj_'+c for c in ['rho','P','BH','LOO_agree','LOO_valid']});bgd=damage[(damage.branch=='BG1')&(damage.program_family=='primary')&(damage.model=='Pooled')][['scope','gene','module','rho','P','BH']].rename(columns={c:'BG1_'+c for c in ['rho','P','BH']});d=d.merge(adj,on=['scope','gene','module'],validate='1:1').merge(bgd,on=['scope','gene','module'],validate='1:1').merge(quality[['scope','canonical','logFC']],left_on=['scope','gene'],right_on=['scope','canonical'],validate='m:1');d['pooled_pass']=(d.P<.05)&(np.sign(d.rho)*np.sign(d.logFC)>0)&(d.LOO_agree>=5)&(np.sign(d.BG1_rho)==np.sign(d.rho));d['primary_pass']=d.pooled_pass&(d.adj_P<.05)&(np.sign(d.adj_rho)==np.sign(d.rho));wr(d,'03_damage/primary_program_evidence.tsv.gz')
a=mac[mac.model=='Pooled'].copy();b=mac[mac.model=='Adjusted'][['sender','scope','ligand','gene','rho','P','BH','LOO_agree','LOO_valid']].rename(columns={c:'adj_'+c for c in ['rho','P','BH','LOO_agree','LOO_valid']});a=a.merge(b,on=['sender','scope','ligand','gene'],validate='1:1');a['pooled_pass']=(a.P<.05)&a.disease_support&a.coherent&(a.LOO_agree>=5)&(np.sign(a.BG1_rho)==np.sign(a.rho));a['primary_pass']=a.pooled_pass&(a.adj_P<.05)&(np.sign(a.adj_rho)==np.sign(a.rho));wr(a,'04_macrophage/consolidated_all_links.tsv.gz')
summary=quality.rename(columns={'canonical':'gene'}).copy()
for col,df in [('damage_pooled',d[d.pooled_pass]),('damage_primary',d[d.primary_pass]),('mac_pooled',a[a.pooled_pass]),('mac_primary',a[a.primary_pass])]:
 cnt=df.groupby(['scope','gene']).size();summary[col]=[int(cnt.get((s,g),0)) for s,g in zip(summary.scope,summary.gene)]
summary['data_primary']=summary.quality_pass&summary.DE_pass&(summary.damage_primary>0)&(summary.mac_primary>0)&summary.scope.ne('Core');summary['data_exploratory']=summary.quality_pass&summary.DE_pass&(summary.damage_pooled>0)&(summary.mac_pooled>0)&summary.scope.ne('Core');summary['missing_axes']=summary.apply(lambda r:';'.join(k for k,b in [('quality',r.quality_pass),('DE',r.DE_pass),('damage_primary',r.damage_primary>0),('mac_primary',r.mac_primary>0),('state_specific',r.scope!='Core')] if not b),axis=1);wr(summary,'05_intersection/all_gene_state_evidence.tsv.gz');wr(summary[summary.data_primary],'05_intersection/primary_data_candidates.tsv');wr(summary[summary.data_exploratory],'05_intersection/exploratory_data_candidates.tsv');wr(d[d.gene.isin(summary.loc[summary.data_primary|summary.data_exploratory,'gene'])],'05_intersection/shortlist_damage_links.tsv');wr(a[a.gene.isin(summary.loc[summary.data_primary|summary.data_exploratory,'gene'])],'05_intersection/shortlist_macrophage_links.tsv')
counts=[]
for sc,b in summary.groupby('scope'):
 counts.append(dict(scope=sc,tested=len(b),DE=int(b.DE_pass.sum()),Q_DE=int((b.DE_pass&b.quality_pass).sum()),damage_pooled=int(((b.damage_pooled>0)&b.DE_pass&b.quality_pass).sum()),damage_primary=int(((b.damage_primary>0)&b.DE_pass&b.quality_pass).sum()),mac_primary=int(((b.mac_primary>0)&b.DE_pass&b.quality_pass).sum()),data_primary=int(b.data_primary.sum()),data_exploratory=int(b.data_exploratory.sum())))
wr(pd.DataFrame(counts),'05_intersection/stage_counts.tsv');wr(pd.DataFrame(manifest).drop_duplicates('path'),'provenance/input_manifest.tsv');status={'stage':'computational_screen','status':'COMPLETE','primary_gene_states':int(summary.data_primary.sum()),'exploratory_gene_states':int(summary.data_exploratory.sum()),'primary_genes':sorted(summary.loc[summary.data_primary,'gene'].unique()),'exploratory_genes':sorted(summary.loc[summary.data_exploratory,'gene'].unique()),'novelty_pending':True};(O/'provenance/computation_complete.json').write_text(json.dumps(status,indent=2),encoding='utf-8');print(pd.DataFrame(counts).to_string(index=False));print(json.dumps(status),flush=True)