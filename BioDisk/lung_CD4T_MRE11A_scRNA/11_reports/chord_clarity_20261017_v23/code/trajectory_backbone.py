from pathlib import Path
import json,numpy as np,pandas as pd,h5py
from anndata.io import read_elem
from scipy.sparse.csgraph import connected_components,minimum_spanning_tree
from scipy import sparse
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/chord_clarity_20261017_v23';T=O/'plot_data';T.mkdir(exist_ok=True)
f=P/'02_analysis/GSE314943/stage9_20260906_080538_wN3sjO/objects/K15_trajectory.h5ad'
with h5py.File(f,'r') as h:
 obs=read_elem(h['obs']);uns=read_elem(h['uns']);print('GRAPH_KEYS',list(h['obsp'].keys()),'NEIGHBORS',uns.get('neighbors',{}),flush=True)
 key=uns.get('neighbors',{}).get('connectivities_key','connectivities');G=read_elem(h['obsp'][key]).tocsr();U=h['obsm']['X_umap_Harmony'][:]
pt=pd.read_csv(P/'11_reports/customer_sequence_20261015_v2/03A_state_pseudotime/pseudotime_cell_membership.tsv.gz',sep='\t').rename(columns={'index':'cell'}).set_index('cell');pt=pt.loc[obs.index];assert len(pt)==1734
G=G.maximum(G.T);ncomp,lab=connected_components(G,directed=False);counts=np.bincount(lab);comp=int(np.argmax(counts));ids=np.flatnonzero(lab==comp);assert len(ids)>=len(obs)*.8
q=np.asarray(pd.qcut(pd.Series(pt.K15_R1.to_numpy()[ids]).rank(method='first'),12,labels=False),int);H=sparse.csr_matrix((np.ones(len(ids)),(np.arange(len(ids)),q)),shape=(len(ids),12));ng=np.asarray(H.sum(0)).ravel();mass=np.asarray((H.T@G[ids][:,ids]@H).todense());np.fill_diagonal(mass,0);strength=mass/np.sqrt(ng[:,None]*ng[None,:]);tree=minimum_spanning_tree(sparse.csr_matrix(-strength)).tocoo();nodes=[]
for b in range(12):
 jj=ids[q==b];xy=U[jj];med=np.median(xy,axis=0);which=jj[np.argmin(np.sum((xy-med)**2,axis=1))];nodes.append(dict(node=b,cells=len(jj),UMAP1=U[which,0],UMAP2=U[which,1],median_DPT=float(pt.K15_R1.iloc[jj].median()),representative_cell=obs.index[which]))
nodes=pd.DataFrame(nodes);edges=[]
for i,j in zip(tree.row,tree.col):
 if nodes.loc[i,'median_DPT']>nodes.loc[j,'median_DPT']:i,j=j,i
 edges.append(dict(source=int(i),target=int(j),graph_connectivity=float(strength[i,j]),source_DPT=nodes.loc[i,'median_DPT'],target_DPT=nodes.loc[j,'median_DPT'],x1=nodes.loc[i,'UMAP1'],y1=nodes.loc[i,'UMAP2'],x2=nodes.loc[j,'UMAP1'],y2=nodes.loc[j,'UMAP2']))
assert len(edges)==11 and all(e['graph_connectivity']>0 and e['source_DPT']<e['target_DPT'] for e in edges)
nodes.to_csv(T/'trajectory_backbone_nodes.tsv',sep='\t',index=False);pd.DataFrame(edges).to_csv(T/'trajectory_backbone_edges.tsv',sep='\t',index=False)
qa={'source':str(f),'source_graph':key,'cells':len(obs),'components':ncomp,'component_sizes':counts.tolist(),'overlay_cells':len(ids),'bins':12,'edges':len(edges),'method':'12 DPT-quantile bins within largest connected component; sum frozen-neighbor connectivity /sqrt(bin sizes); maximum spanning tree; positions closest cell to bin median UMAP; arrows low to high existing K15_R1 DPT','interpretation':'descriptive neighborhood backbone, not independently inferred lineage or real-time path'};(O/'trajectory_backbone_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8');print(json.dumps(qa),flush=True)