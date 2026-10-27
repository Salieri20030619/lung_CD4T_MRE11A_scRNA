from pathlib import Path
import json,subprocess,collections
O=Path(__file__).parent;ps=json.loads((O/'server_file_inventory.json').read_text());P='/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/'
for fragment in ['customer_sequence','manuscript_graphics_20261006_v5','broad_clustering_v1','CD4_final']:
 sel=[p[len(P):] for p in ps if fragment in p and '/package/' not in p and '/code/' not in p and '/source_tables/' not in p and '/supplement/' not in p]
 print('\n',fragment,'count',len(sel));print('\n'.join(sel[:50]))