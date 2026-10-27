import os, json, importlib, importlib.util
from pathlib import Path
for variable,expected in [('PREV1','mechanism_step1'),('PREV2','mechanism_step2_human_and_driver_audit'),('PREV8','mechanism_step8_human_mac_input_audit'),('SCENIC1','lung_scenic1_network'),('SCENIC2','lung_scenic2_activity_sensitivity')]:
    p=Path(os.environ[variable]);s=json.loads((p/'status.json').read_text())
    assert s['status']=='PASS' and s['stage']==expected,(variable,s)
for variable,files in [
 ('PREV1',['comm_input/CD4_counts.mtx.gz','comm_input/CD4_genes.tsv','comm_input/CD4_cells.tsv','tables/module_definition.tsv']),
 ('PREV2',['input/human_counts.mtx','input/human_genes.tsv','input/human_cells.tsv']),
 ('PREV8',['objects/marker_audit_checkpoint.rds'])]:
    for name in files:
        p=Path(os.environ[variable])/name;assert p.is_file() and p.stat().st_size>0,p
for arm in ['All','BG1']:
    for seed in [17,29,43]:
        p=Path(os.environ['SCENIC2'])/'objects'/f'Mouse_{arm}_{seed}_DropOverlap_auc.csv.gz'
        assert p.is_file() and p.stat().st_size>0,p
    for sp in ['Mouse','Human']:assert (Path(os.environ['SCENIC1'])/sp/arm/'cells.tsv').is_file()
versions={}
for pkg in ['numpy','pandas','scipy','scanpy','anndata','matplotlib']:
    mod=importlib.import_module(pkg);versions[pkg]=getattr(mod,'__version__','NA')
assert importlib.util.find_spec('leidenalg') or importlib.util.find_spec('louvain'),'Clustering dependency unavailable'
(Path(os.environ['OUT'])/'python_versions.json').write_text(json.dumps(versions,indent=2))
print('ANALYSIS_LOGIC_PASS FROZEN_INPUTS_AND_DEPENDENCIES',flush=True)

