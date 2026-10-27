import os, ast, json, hashlib, tarfile
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def labels_text(s):
    return s.astype('string').fillna('Unknown').astype(str)

def main():
    import anndata as ad
    from PIL import Image
    out=Path(os.environ['OUT']);root=Path(os.environ['ROOT'])
    names=['mouse_marker_sample.tsv','mouse_TF_identity_adjustment.tsv','mouse_FOXP3_strata.tsv','mouse_identity_sample.tsv','mouse_identity_cells.tsv.gz','CD4_receptor_detection.tsv','CD4_receptor_codetection.tsv','human_cluster_markers.tsv','human_cluster_scorecard.tsv','human_cluster_sample_composition.tsv','human_marker_coverage.tsv','human_cluster_marker_expression.tsv','human_sender_cluster_expression.tsv','human_sender_parent_sensitivity.tsv','human_sender_sample_covariates.tsv','human_refined_cells.tsv.gz']
    for name in names:
        p=out/'tables'/name
        assert p.is_file() and p.stat().st_size>0, f'MISSING_TABLE={p}'
        assert not pd.read_csv(p,sep='\t',nrows=1).empty, f'EMPTY_TABLE={p}'
    p=out/'objects/human_myeloid_refined.h5ad';assert p.is_file(),f'MISSING_CHECKPOINT={p}'
    a=ad.read_h5ad(p,backed='r')
    try:
        xy=np.asarray(a.obsm['X_umap']);obs=a.obs.copy()
        assert xy.shape==(len(obs),2) and np.isfinite(xy).all() and obs.index.is_unique
        md=pd.read_csv(out/'tables/human_refined_cells.tsv.gz',sep='\t',dtype={'cell':str})
        assert np.array_equal(md.cell.to_numpy(),obs.index.astype(str).to_numpy())
        for key in ['r4','parent','group','series']:
            assert np.array_equal(labels_text(md[key]).to_numpy(),labels_text(obs[key]).to_numpy()),key
        algorithm=str(a.uns.get('clustering_algorithm','recorded_in_original_run'))
    finally:
        a.file.close()
    print('ANALYSIS_LOGIC_PASS SAVED_TABLES_CELL_IDS_UMAP',flush=True)
    mapping=[]
    for key,name in [('r4','MacClust'),('parent','MacParent'),('group','MacGroup'),('series','MacSource')]:
        labs=labels_text(obs[key]);levels=sorted(labs.unique());fig,ax=plt.subplots(figsize=(8,6))
        for j,lab in enumerate(levels):
            keep=labs.eq(lab).to_numpy();label=lab if len(lab)<=10 else 'L'+str(j)
            ax.scatter(xy[keep,0],xy[keep,1],s=2,c=[plt.get_cmap('tab20')(j%20)],label=label,rasterized=True)
            mapping.append(dict(figure=name,label=label,original=lab))
        ax.set(xlabel='UMAP1',ylabel='UMAP2',title=name)
        ax.legend(markerscale=4,bbox_to_anchor=(1.01,1),loc='upper left',frameon=False,fontsize=7);fig.tight_layout()
        for ext in ['png','pdf']:
            dest=out/'figures'/f'{name}.{ext}';part=Path(str(dest)+'.part')
            fig.savefig(part,format=ext,dpi=300)
            if ext=='png':
                with Image.open(part) as im:im.verify()
            else:
                with part.open('rb') as f:assert f.read(5)==b'%PDF-'
            part.replace(dest)
        plt.close(fig);print('FIGURE_PASS='+name,flush=True)
    pd.DataFrame(mapping).to_csv(out/'tables/figure_label_map.tsv',sep='\t',index=False)
    (out/'figures/captions.txt').write_text('MacClust: human candidate myeloid clusters at resolution 0.4; final identities await marker review. MacParent: original broad clusters 3/6/16/22. MacGroup: Control/CWP/Unknown. MacSource: original dataset series. UMAP1/UMAP2: unchanged saved embedding. Full labels: tables/figure_label_map.tsv. Algorithm: '+algorithm+'. No reclustering or statistical refitting performed during resume.\n')
    # Recover the exact original methods text without importing or running analyze.py.
    tree=ast.parse((out/'analyze.py').read_text());texts=[]
    for n in ast.walk(tree):
        if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='write_text':
            v=n.func.value
            if isinstance(v,ast.BinOp) and isinstance(v.right,ast.Constant) and v.right.value=='methods.txt':texts.append(ast.literal_eval(n.args[0]))
    assert len(texts)==1
    (out/'methods.txt').write_text(texts[0]+'Resume: categorical labels converted to nullable strings before filling missing values; saved UMAP, cluster labels, matrices, and statistical results unchanged.\n')
    status=dict(status='PASS',stage='lung_scenic3_identity_sender_refinement',output=str(out),scenic1=str(root/'02_analysis/scenic1_20260929_103513_xRoEnr'),scenic2=str(root/'02_analysis/scenic2_20261001_092655_ZKNWD6'),resumed_from='saved_human_myeloid_refined.h5ad',pending=['review_human_myeloid_subtype_assignments','final_sender_TF_models_using_reviewed_labels'])
    (out/'status.json').write_text(json.dumps(status,indent=2))
    ret=root/'tmp/lung_scenic3_return.tar.gz';part=Path(str(ret)+'.part')
    with tarfile.open(part,'w:gz') as t:
        for f in sorted(out.rglob('*')):
            if f.is_file() and not {'input','objects','__pycache__'}&set(f.relative_to(out).parts) and not f.name.endswith('.part'):t.add(f,arcname=str(f.relative_to(out)))
    with tarfile.open(part) as t:
        for member in t:
            if member.isfile():
                with t.extractfile(member) as f:
                    while f.read(2**20):pass
    part.replace(ret);print('RETURN='+str(ret),flush=True);print('RETURN_SHA256='+hashlib.sha256(ret.read_bytes()).hexdigest(),flush=True)
    print('LUNG_SCENIC3_PASS',flush=True)
if __name__=='__main__':main()

