from pathlib import Path
import pandas as pd,json
P=Path(__file__).parent;T=P/'tables'
e=pd.read_csv(T/'network_LOO_all_edges.tsv.gz',sep='\t');b=pd.read_csv(T/'network_baseline_seed17_edges.tsv.gz',sep='\t')
a=pd.read_csv(T/'network_LOO_fit_audit.tsv',sep='\t')
rows=[]
for (arm,left),q in e.groupby(['arm','left_out']):
 full=b[b.arm.eq(arm)&b.tf.ne(b.target)];fold=q[q.tf.ne(q.target)]
 bs=set(zip(full.tf,full.target));fs=set(zip(fold.tf,fold.target))
 r=a[a.arm.eq(arm)&a.left_out.eq(left)].iloc[0].to_dict()
 r.update(edges=len(fs),baseline_seed17_edges=len(bs),shared=len(bs&fs),jaccard=len(bs&fs)/len(bs|fs),self_members_excluded=int(q.tf.eq(q.target).sum()))
 rows.append(r)
out=pd.DataFrame(rows);assert len(out)==20
out.to_csv(T/'network_LOO_nonself_fit_audit.tsv',sep='\t',index=False,encoding='utf-8')
(P/'structural_QA.json').write_text(json.dumps({'status':'success','issue':'STRUCT01: software TF-self membership is not interpreted as an autoregulatory edge','resolution':'nonself global edge metrics and display; raw memberships preserved; RNA versus own-TF activity tests unchanged','folds_checked':20,'self_members_removed':int(out.self_members_excluded.sum())},indent=2),encoding='utf-8')
f=P/'supplement_plots.py';s=f.read_text(encoding='utf-8')
s=s.replace("d.to_csv(T/'network_LOO_candidate_edges_annotated.tsv'","d['self_signature_member']=d.tf.eq(d.target);d.to_csv(T/'network_LOO_candidate_edges_annotated.tsv'")
s=s.replace("d=d[(d.arm=='All')&d.stable_2of3","d=d[(d.arm=='All')&d.tf.ne(d.target)&d.stable_2of3")
f.write_text(s,encoding='utf-8')
f=P/'supplement_v3_report.py';s=f.read_text(encoding='utf-8')
s=s.replace("la=xr('network_LOO_fit_audit.tsv')","la=xr('network_LOO_nonself_fit_audit.tsv')")
s=s.replace('全图边集合的Jaccard约为0.19至0.26','剔除TF自身成员后，全图边集合的Jaccard约为0.19至0.26')
s=s.replace('展示原All三种子稳定候选边在10折中的出现次数；','展示原All三种子稳定候选边在10折中的出现次数，剔除TF自身成员；')
s=s.replace('network_LOO_fit_audit.tsv；network_LOO_candidate_edges.tsv','network_LOO_nonself_fit_audit.tsv；network_LOO_candidate_edges.tsv')
f.write_text(s,encoding='utf-8')
print(out.groupby('arm').jaccard.agg(['min','max']).to_string())
