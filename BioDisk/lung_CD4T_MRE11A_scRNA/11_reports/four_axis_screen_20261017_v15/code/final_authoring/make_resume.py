from pathlib import Path
W=Path('C:/Users/31937/lung_four_axes_20261017_v15');s=(W/'run_screen.py').read_text(encoding='utf-8');prefix=s.split("print('START'")[0];tail=s.split('# Export all top-200 positive prior links without looking at damage results.')[1]
load="""
print('RESUME_FROM_SAVED_QUALITY_DAMAGE',flush=True)
gn=rd(M/'comm_input/CD4_genes.tsv');genes=pd.Index(gn.canonical);quality=rd(O/'01_quality/all_genes_quality.tsv.gz');deall=rd(O/'02_DE/all_frozen_edgeR_results.tsv.gz');damage=rd(O/'03_damage/all_gene_program_tests.tsv.gz');agg={}
for sc in SC:
 a=rd(O/f'01_quality/{sc}_all_gene_sample_values.tsv.gz')
 for branch in ['All','BG1']:
  b=a[a.branch==branch].pivot(index='sample',columns='gene',values='CP10k').reindex(index=S,columns=genes);agg[(sc,branch)]=b.to_numpy()
for f in [M/'comm_input/CD4_counts.mtx.gz',M/'comm_input/CD4_cells.tsv']:
 manifest.append(dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
"""
(W/'resume_screen.py').write_text(prefix+load+tail,encoding='utf-8')
p=W/'run_driver.py';s=p.read_text(encoding='utf-8').replace("str(O/'code/run_screen.py')","str(O/'code/resume_screen.py')");p.write_text(s,encoding='utf-8')