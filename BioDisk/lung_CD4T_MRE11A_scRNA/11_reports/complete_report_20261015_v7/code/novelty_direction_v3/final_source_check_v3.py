from pathlib import Path
import sys,json,xml.etree.ElementTree as ET
O=Path(sys.argv[1]);root=ET.fromstring((O/'literature/cache/PMC7290873.xml').read_bytes())
for row in root.iter('tr'):
 text=' '.join(row.itertext())
 if 'ST6GALNAC3' in text:print('ROW_LABEL',text[:80])
print('DIRECTION_QA',json.loads((O/'provenance/direction_complete.json').read_text(encoding='utf-8')))