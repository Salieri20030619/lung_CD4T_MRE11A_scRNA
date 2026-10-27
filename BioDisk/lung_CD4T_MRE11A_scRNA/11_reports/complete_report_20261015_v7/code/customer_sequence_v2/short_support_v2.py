from pathlib import Path
import pandas as pd,sys
O=Path(sys.argv[1]);P=O.parents[1]
r=pd.read_csv(O/'03A_state_pseudotime/genomewide_trajectory_tests.tsv.gz',sep='\t')
for state,g in [('Naive','AHNAK'),('ActCD4','MDFIC'),('ActCD4','TTC39C'),('ActCD4','ST6GALNAC3')]:
 a=r[(r.state==state)&(r.feature==g)];print(state,g,[(x.run,round(x.mean_slope,3),x.P) for x in a.itertuples()])
for name in ['candidate_ligand_prior.tsv','candidate_prior_coverage.tsv','supported_LR_edges.tsv']:
 f=P/'02_analysis/mechanism_step3_20260928_103057_lIRL38/tables'/name;a=pd.read_csv(f,sep='\t');print(name,a.shape,list(a.columns));print(a.head(1).to_dict('records'))
for d in ['02_analysis/mechanism_step7_20260929_100321_6hm4dZ','02_analysis/mechanism_step6_20260928_163015_VoRQqA']:
 print(d,[(str(f.relative_to(P/d)),f.stat().st_size) for f in (P/d).rglob('*') if f.is_file() and any(w in f.name.lower() for w in ['prior','regulon','target'])][:15])