"""Restore deduplicated source-code aliases into a NEW directory; never run them."""
from pathlib import Path
import argparse,csv,hashlib,json,shutil
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--package',type=Path,required=True)
parser.add_argument('--destination',type=Path,help='Must not already exist; omit for dry-run')
args=parser.parse_args();base=args.package.resolve()
rows=list(csv.DictReader((base/'provenance/code_manifest.tsv').open(encoding='utf-8'),delimiter='\t'))
prefix='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/'
plan=[]
for r in rows:
 assert r['original_path'].startswith(prefix)
 rel=Path(r['original_path'][len(prefix):]);assert not rel.is_absolute() and '..' not in rel.parts
 source=(base/r['package_path']).resolve();assert source.is_relative_to(base)
 assert hashlib.sha256(source.read_bytes()).hexdigest()==r['package_sha256']
 plan.append((source,rel))
if args.destination:
 dest=args.destination.resolve()
 assert not dest.exists(),'Destination must be new; nothing will be overwritten.'
 assert str(dest) not in [prefix.rstrip('/'),'/']
 assert not base.is_relative_to(dest) and not dest.is_relative_to(base),'Use a directory outside this package.'
 dest.mkdir(parents=True)
 for source,rel in plan:
  target=dest/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
 print(json.dumps({'restored_alias_paths':len(plan),'destination':str(dest),'analysis_executed':False}))
else:print(json.dumps({'verified_source_paths':len(plan),'mode':'dry_run','analysis_executed':False}))
