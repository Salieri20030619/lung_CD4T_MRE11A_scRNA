from pathlib import Path
import subprocess,time,json
O=Path(__file__).parent;P='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports';remote=P+'/client_layout_20261017_v22'
s=(O/'reference_plots.R').read_text(encoding='utf-8');s=s.replace('C:/Users/31937/lung_client_layout_20261017_v22',remote).replace('C:/Users/31937/lung_client_rebuild_20261017_v20',P+'/client_rebuild_20261017_v20').replace('C:/Users/31937/lung_feedback_revision_20261017_v19',P+'/feedback_revision_20261017_v19').replace('server_tables/','tables/');(O/'reference_plots_server.R').write_text(s,encoding='utf-8',newline='')
def run(cmd):
 for i in range(4):
  r=subprocess.run(cmd,capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
  if r.returncode==0:return r
  if r.returncode!=255:break
  time.sleep(2)
 (O/'server_failure.log').write_bytes(r.stdout+r.stderr);r.check_returncode()
run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','mkdir -p '+remote+'/code '+remote+'/figures '+remote+'/plot_data'])
run(['scp','-q',str(O/'reference_plots_server.R'),'sgs10:'+remote+'/code/'])
r=run(['ssh','-T','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript --vanilla '+remote+'/code/reference_plots_server.R']);(O/'R_server.log').write_bytes(r.stdout+r.stderr)
run(['scp','-q','sgs10:'+remote+'/figures/*',str(O/'figures')+'/']);run(['scp','-q','sgs10:'+remote+'/plot_data/*',str(O/'plot_data')+'/']);run(['scp','-q','sgs10:'+remote+'/R_plots_complete.txt','sgs10:'+remote+'/R_plot_session.txt',str(O)+'/'])
(O/'R_server_exit.json').write_text(json.dumps({'exit_code':r.returncode,'completion_marker':(O/'R_plots_complete.txt').read_text().strip(),'local_windows_status':'graphics generated but R crashed on exit; replaced by server output'},indent=2),encoding='utf-8');print('SERVER_R_SUCCESS')