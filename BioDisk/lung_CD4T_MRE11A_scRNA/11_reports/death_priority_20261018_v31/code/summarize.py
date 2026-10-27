from pathlib import Path
import pandas as pd,json
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v31')
r=pd.read_csv(O/'tables/module_state_ranking.tsv',sep='\t');print('RANKING\n',r.to_string(index=False));print('DISEASE\n',pd.read_csv(O/'tables/module_disease_effects.tsv',sep='\t').to_string(index=False));t=pd.read_csv(O/'tables/time_all12_D28_vs_D7.tsv',sep='\t');print('TIME20\n',t[t.PValue<.05][['gene','logFC','PValue','FDR']].sort_values('PValue').to_string(index=False))
a=pd.read_csv(O/'tables/gene_trajectory_all.tsv.gz',sep='\t');a=a[a.variant.eq('K15_R1')];g=pd.read_csv(O/'tables/GO_priority_members.tsv',sep='\t');death=set(g.loc[g.GOALL.isin(['GO:0006915','GO:0070269']),'gene']);tm=pd.read_csv(O/'tables/time_all12_time_omnibus.tsv',sep='\t');rows=[]
for sc,z in a.groupby('scope'):
 A=set(z.loc[z.P<.05,'gene']);B=set(t.loc[t.PValue<.05,'gene']);C=set(tm.loc[tm.PValue<.05,'gene']);rows.append({'scope':sc,'trajectory':len(A),'D28D7':len(B),'overlap':len(A&B),'death_primary':sorted(A&B&death),'secondary_death':sorted(A&C&death)})
print('INTERSECTIONS',json.dumps(rows));print('GRAPH_QA\n',pd.read_csv(O/'tables/graph_QA.tsv',sep='\t').to_string(index=False))