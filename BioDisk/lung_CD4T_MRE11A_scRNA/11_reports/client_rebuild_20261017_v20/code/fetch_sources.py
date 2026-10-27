from pathlib import Path
import json,subprocess,tarfile
O=Path(__file__).parent;P='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/';ps=json.loads((O/'server_file_inventory.json').read_text());sel=[]
for p in ps:
 rel=p[len(P):]
 ok=(rel.startswith('11_reports/customer_sequence_20261015_v2/03A_state_pseudotime/') or rel.startswith('11_reports/customer_sequence_20261015_v2/provenance/prospective_specification') or rel.startswith('11_reports/customer_sequence_20261015_v2/code/') or rel.startswith('11_reports/manuscript_graphics_20261006_v5/source_tables/mouse_pathways/') or rel.startswith('11_reports/manuscript_graphics_20261006_v5/provenance/plot_input_manifest') or rel.startswith('11_reports/manuscript_graphics_20261006_v5/code/graphics_inputs') or rel.startswith('11_reports/manuscript_graphics_20261006_v5/figures/fig1_'))
 if ok and p.endswith(('.tsv','.tsv.gz','.json','.py','.R','.png','.txt')):sel.append(rel)
(O/'fetch_manifest.txt').write_text('\n'.join(sel)+'\n',encoding='utf-8',newline='');kw=dict(check=True,creationflags=subprocess.CREATE_NO_WINDOW)
remote=P+'11_reports/client_rebuild_20261017_v20';subprocess.run(['ssh','-T','sgs10','mkdir -p '+remote+'/code '+remote+'/inputs '+remote+'/tables'],**kw)
subprocess.run(['scp','-q',str(O/'fetch_manifest.txt'),'sgs10:'+remote+'/inputs/fetch_manifest.txt'],**kw)
subprocess.run(['ssh','-T','sgs10','tar -czf '+remote+'/inputs/existing_sources.tar.gz -C '+P+' -T '+remote+'/inputs/fetch_manifest.txt'],**kw)
subprocess.run(['scp','-q','sgs10:'+remote+'/inputs/existing_sources.tar.gz',str(O/'existing_sources.tar.gz')],**kw)
with tarfile.open(O/'existing_sources.tar.gz') as t:t.extractall(O/'sources',filter='data')
print('downloaded',len(sel),'files', (O/'existing_sources.tar.gz').stat().st_size)