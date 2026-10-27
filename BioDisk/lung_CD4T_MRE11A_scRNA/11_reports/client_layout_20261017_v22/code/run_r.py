from pathlib import Path
import os,subprocess,json
O=Path(__file__).parent;env=os.environ.copy()
for k in list(env):
 if k.startswith('LC_'):env.pop(k)
env.update(LANG='C',LC_ALL='C')
r=subprocess.run([r'D:\R-4.5.2\bin\x64\Rscript.exe','--vanilla',str(O/'reference_plots.R')],env=env,capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
(O/'reference_plots.stdout.log').write_bytes(r.stdout);(O/'reference_plots.stderr.log').write_bytes(r.stderr);(O/'R_exit.json').write_text(json.dumps({'exit_code':r.returncode,'completion_marker':(O/'R_plots_complete.txt').read_text().strip()}),encoding='utf-8');print('R_EXIT',r.returncode);print(r.stderr.decode('utf-8'));r.check_returncode()