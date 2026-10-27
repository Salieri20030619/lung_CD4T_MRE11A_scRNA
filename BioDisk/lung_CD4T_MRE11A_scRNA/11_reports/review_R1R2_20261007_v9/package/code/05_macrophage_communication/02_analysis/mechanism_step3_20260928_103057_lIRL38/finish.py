import os,json,tarfile,hashlib
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(os.environ['OUT']);root=Path(os.environ['ROOT'])
required=['DE_status.tsv','geneset_status.tsv','ligand_activities.tsv','candidate_prior_coverage.tsv','candidate_ligand_prior.tsv','supported_LR_edges.tsv','human_composition_correlations.tsv','human_correlation_LOO.tsv','human_source_models_BH.tsv','receptor_complex_codetection.tsv','paired_ligand_target_correlations.tsv']
for f in required:
    assert (out/'tables'/f).is_file(),f
st=pd.read_csv(out/'tables/geneset_status.tsv',sep='\t')
assert len(st)==16,st
ac=pd.read_csv(out/'tables/ligand_activities.tsv',sep='\t')
assert len(st[st.status=='PASS'])==(ac[['receiver','filter','test','direction']].drop_duplicates().shape[0] if len(ac) else 0)
captions=[]
if len(ac):
    keys=['receiver','test','direction','test_ligand']
    merged=ac[ac['filter']=='All'][keys+['rank','aupr_corrected']].merge(ac[ac['filter']=='BG1'][keys+['rank','aupr_corrected']],on=keys,how='outer',suffixes=('_All','_BG1'))
    merged['top20_both']=(merged['rank_All']<=20)&(merged['rank_BG1']<=20)
    merged.to_csv(out/'tables/ligand_rank_sensitivity.tsv',sep='\t',index=False,na_rep='NA')
else:
    pd.DataFrame(columns=['receiver','test','direction','test_ligand','rank_All','rank_BG1','top20_both']).to_csv(out/'tables/ligand_rank_sensitivity.tsv',sep='\t',index=False)
if len(ac):
    for (receiver,direction),v in ac[(ac['filter']=='All')&(ac.test=='FDR05')].groupby(['receiver','direction']):
        v=v.sort_values(['rank','test_ligand']).head(15).iloc[::-1]
        labels=[x if len(x)<=10 else 'L'+str(i+1) for i,x in enumerate(v.test_ligand)]
        fig,ax=plt.subplots(figsize=(6.8,5.4));ax.barh(labels,v.aupr_corrected,color='#3979a8')
        ax.set_xlabel('AUPR adj');ax.set_title(receiver+' / '+direction);ax.axvline(0,color='gray',linewidth=.5)
        fig.tight_layout();stem=f'NicheNet_{receiver}_{direction}';fig.savefig(out/'figures'/f'{stem}.png',dpi=300);fig.savefig(out/'figures'/f'{stem}.pdf');plt.close(fig)
        captions.append(stem+': FDR-defined receiver gene set; corrected AUPR ligand ranking. '+ '; '.join(f'{a}={b}' for a,b in zip(labels,v.test_ligand)))
(out/'captions.txt').write_text('\n'.join(captions))
summary=dict(status='PASS',stage='mechanism_step3_nichenet',output=str(out),previous=os.environ['PREV'],previous2=os.environ['PREV2'],nichenet_sets_run=int((st.status=='PASS').sum()),nichenet_sets_skipped=int((st.status!='PASS').sum()),note='PASS means workflow completed, not proof of any communication axis or cell death.')
(out/'status.json').write_text(json.dumps(summary,indent=2))
archive=root/'tmp/lung_mechanism_step3_return.tar.gz';part=archive.with_suffix('.tar.gz.part')
files=[p for p in out.rglob('*') if p.is_file() and 'objects' not in p.relative_to(out).parts]
with tarfile.open(part,'w:gz') as tar:
    for p in sorted(files):tar.add(p,arcname=str(p.relative_to(out)))
with tarfile.open(part,'r:gz') as tar:
    for m in tar:
        if m.isfile():
            h=tar.extractfile(m)
            while h.read(1024*1024):pass
part.replace(archive)
print('RETURN='+str(archive),flush=True)
print('RETURN_BYTES='+str(archive.stat().st_size),flush=True)
print('RETURN_SHA256='+hashlib.sha256(archive.read_bytes()).hexdigest(),flush=True)
print('NICHE_NET_SETS_RUN='+str(summary['nichenet_sets_run']),flush=True)
print('LUNG_MECHANISM_STEP3_PASS',flush=True)

