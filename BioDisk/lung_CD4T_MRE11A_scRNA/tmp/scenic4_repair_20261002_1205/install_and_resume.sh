#!/usr/bin/env bash
set -Eeuo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA
OUT="$ROOT/02_analysis/scenic4_20261002_091545_iazXf1"
AUDIT="$ROOT/tmp/scenic4_repair_20261002_1205"
PY="$ROOT/software/pyscenic_0.12.1_py310/bin/python"
export AUDIT
(
exec 9>"$ROOT/tmp/lung_scenic4.lock"
flock -n 9
while IFS=$'\t' read -r url name digest; do
  printf '%s  %s\n' "$digest" "$AUDIT/wheels/$name" | sha256sum -c -
done < "$AUDIT/official_wheels.tsv"
"$PY" -u -B -m pip --isolated --disable-pip-version-check install --dry-run --no-index --find-links "$AUDIT/wheels" --constraint "$AUDIT/environment_before.txt" --report "$AUDIT/offline_install_plan.json" matplotlib==3.7.5 cycler==0.12.1 fonttools==4.53.1 kiwisolver==1.4.7 pyparsing==3.1.4 > "$AUDIT/offline_install_plan.log" 2>&1
"$PY" -B - <<'CHECK_PLAN'
import os,json
from pathlib import Path
p=Path(os.environ['AUDIT'])
plan=json.loads((p/'offline_install_plan.json').read_text())
actual={i['metadata']['name'].lower():i['metadata']['version'] for i in plan['install']}
expected={'matplotlib':'3.7.5','cycler':'0.12.1','fonttools':'4.53.1','kiwisolver':'1.4.7','pyparsing':'3.1.4'}
assert actual==expected,actual
print('INSTALL_PLAN_VERIFIED='+json.dumps(actual),flush=True)
CHECK_PLAN
"$PY" -u -B -m pip --isolated --disable-pip-version-check install --no-index --find-links "$AUDIT/wheels" --constraint "$AUDIT/environment_before.txt" --report "$AUDIT/offline_install_actual.json" matplotlib==3.7.5 cycler==0.12.1 fonttools==4.53.1 kiwisolver==1.4.7 pyparsing==3.1.4 > "$AUDIT/install.log" 2>&1
"$PY" -B -m pip freeze > "$AUDIT/environment_after.txt"
"$PY" -B -m pip check > "$AUDIT/pip_check_after.txt" 2>&1
"$PY" -B - <<'CHECK_ENV'
import os,json
from pathlib import Path
p=Path(os.environ['AUDIT'])
def parse(name):
 return dict(x.split('==',1) for x in (p/name).read_text().splitlines() if '==' in x)
before=parse('environment_before.txt');after=parse('environment_after.txt')
assert all(after.get(k)==v for k,v in before.items())
added={k:v for k,v in after.items() if k not in before}
(p/'environment_change.json').write_text(json.dumps({'existing_versions_unchanged':True,'added':added},indent=2))
print('ENVIRONMENT_EXISTING_VERSIONS_UNCHANGED=PASS',flush=True)
CHECK_ENV
sha256sum -c "$AUDIT/original_sha256.txt"
)
mkdir "$OUT/repair_audit"
cp "$AUDIT"/environment_before.txt "$AUDIT"/environment_after.txt "$AUDIT"/environment_change.json "$AUDIT"/pip_check_before.txt "$AUDIT"/pip_check_after.txt "$AUDIT"/offline_install_plan.json "$AUDIT"/offline_install_actual.json "$AUDIT"/original_sha256.txt "$AUDIT"/official_wheels.tsv "$OUT/repair_audit/"
cat > "$OUT/repair_audit/repair_notes.txt" <<'NOTES'
2026-10-02: authorized technical recovery after missing matplotlib import.
Original shell script, validate.py, failed log and exported inputs preserved.
Original baseline SHA256 records were verified before recovery.
Added matplotlib 3.7.5, cycler 0.12.1, fonttools 4.53.1, kiwisolver 1.4.7, pyparsing 3.1.4 into project-local pyscenic Python 3.10 only.
All previously installed package versions held fixed; pip check passed before and after.
Wheels obtained from official PyPI/files.pythonhosted.org and verified against official SHA256 metadata.
Unused importlib_resources wheel was downloaded but NOT installed.
Initial online pip preview was slow; server wheel download required extended timeout and resume. Local download attempt timed out.
resume_run.sh reuses existing input exports and executes original unchanged validate.py.
Full runtime import preflight added to recovery wrapper, plus dedicated log and exclusive lock.
No scientific changes: frozen IKZF2/STAT1 regulons, GSE183682 days7/56, All/BG1 <=0.01, seeds17/29/43, >=20 cells, original DropOverlap exclusions, shared measured gene universe, coverage reported without new cutoff.
No inferential P or seed-as-biological-replicate tests. Step8 unchanged.
Raw recovery audit directory: /mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/tmp/scenic4_repair_20261002_1205
NOTES
LOG="$ROOT/logs/lung_scenic4_resume_20261002_1205.log"
test ! -e "$LOG"
nohup bash "$OUT/resume_run.sh" > "$LOG" 2>&1 < /dev/null &
pid=$!
printf '%s\n' "$pid" > "$AUDIT/resume.pid"
printf 'RESUME_LAUNCHED_PID=%s LOG=%s TIME=%s\n' "$pid" "$LOG" "$(date -Is)"
