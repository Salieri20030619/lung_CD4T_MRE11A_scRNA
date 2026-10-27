from pathlib import Path
import hashlib,json,collections
W=Path('C:/Users/31937/lung_review_20261017_v14');r=json.loads((W/'analysis/provenance/baseline_remote_sha256.json').read_text(encoding='utf-8'));l={f.relative_to(W/'baseline_v13').as_posix():hashlib.sha256(f.read_bytes()).hexdigest() for f in (W/'baseline_v13').rglob('*') if f.is_file()};print('COUNT',len(r),len(l))
for p,h in r.items():
 if h not in l.values():print('REMOTE_ONLY',p,h)
for p,h in l.items():
 if h not in r.values():print('LOCAL_ONLY',p,h)