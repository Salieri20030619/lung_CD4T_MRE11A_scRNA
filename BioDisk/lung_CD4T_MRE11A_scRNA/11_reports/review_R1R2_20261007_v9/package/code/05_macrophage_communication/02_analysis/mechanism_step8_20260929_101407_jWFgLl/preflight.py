import os,json
from pathlib import Path
import pandas as pd
root=Path(os.environ['ROOT']);p7=Path(os.environ['PREV7']);p2=Path(os.environ['PREV2']);out=Path(os.environ['OUT']);obj=Path(os.environ['HUMAN_OBJECT'])
for p,stage in [(p7,'mechanism_step7_human_TF'),(p2,'mechanism_step2_human_and_driver_audit')]:
 s=json.loads((p/'status.json').read_text());assert s['status']=='PASS' and s['stage']==stage
inventory=pd.read_csv(p7/'tables/human_mac_object_inventory.tsv',sep='\t');row=inventory[inventory.path==str(obj)]
assert len(row)==1 and obj.is_file() and obj.stat().st_size==int(row.bytes.iloc[0]),'Selected object absent from verified inventory or size changed'
assert (p2/'input/human_cells.tsv').is_file()
(out/'selected_input.json').write_text(json.dumps({'path':str(obj),'bytes':obj.stat().st_size,'inventory':str(p7/'tables/human_mac_object_inventory.tsv'),'identity_status':'TO_BE_AUDITED'},indent=2))
print('INPUT_INVENTORY_CHECK=PASS',flush=True)
mem=Path('/proc/meminfo')
if mem.exists():
 for line in mem.read_text().splitlines():
  if line.startswith('MemAvailable:'):print(line,flush=True)

