from pathlib import Path
import zipfile,hashlib,json,collections,shutil,pandas as pd
W=Path('C:/Users/31937/lung_review_20261017_v14');f=W/'review_analysis_bundle.zip';assert hashlib.sha256(f.read_bytes()).hexdigest()=='e1f8eecadd92fd4d0e2086fee50547ab26af70014c81c2a336951ca39f858c7b'
with zipfile.ZipFile(f) as z:z.extractall(W/'analysis')
remote=json.loads((W/'analysis/provenance/baseline_remote_sha256.json').read_text(encoding='utf-8'));local={str(f.relative_to(W/'baseline_v13')).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in (W/'baseline_v13').rglob('*') if f.is_file() and not f.name.endswith('.inspect.ndjson')}
assert collections.Counter(remote.values())==collections.Counter(local.values()),'baseline contents differ'
print('BASELINE_PASS',len(local),'filename differences',len(set(remote)-set(local)))
O=W/'revision_v14';shutil.copytree(W/'baseline_v13',O,dirs_exist_ok=True,ignore=shutil.ignore_patterns('*.inspect.ndjson'));shutil.copytree(W/'analysis/supplement_v14',O/'supplement_v14',dirs_exist_ok=True);shutil.copytree(W/'analysis/provenance',O/'provenance/v14',dirs_exist_ok=True)
for name in ['S23_targeted_proximal_associations.tsv','S23_proximal_AUC_disease_original_family.tsv']:
 a=pd.read_csv(O/'supplement_v14'/name,sep='\t');print(name)
 if 'ligand' in a:print(a[(a.scope=='ActCD4')&a.ligand.isin(['IL15','IFNB1','IL27','IFNG'])][['sender','ligand','tf','model','rho','P','status']].to_string(index=False))
 else:print(a[(a.scope=='ActCD4')][['regulon','delta_mean_auc','P','BH']].to_string(index=False))
print('FAMILYFIX',pd.read_csv(O/'supplement_v14/S21_complete_TF_fixed.tsv',sep='\t').columns.to_list())