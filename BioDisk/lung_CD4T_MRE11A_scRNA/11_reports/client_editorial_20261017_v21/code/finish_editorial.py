from pathlib import Path
import json,shutil
O=Path(__file__).parent;P=O/'editorial_replacements.json';r=json.loads(P.read_text(encoding='utf-8'))
r.insert(0,['尚不能将任一候选写成已确认的“巨噬细胞诱导CD4线粒体损伤或焦亡”机制。','现有证据尚不足以确认任一候选介导“巨噬细胞诱导CD4线粒体损伤或焦亡”的机制。'])
r.append(['整合多条件交集筛选获得的55个基因与研究初始阶段纳入的42个候选基因，两组重叠7个，形成90个基因的综合评价集合。','综合差异表达、功能注释、细胞通讯网络及候选相关证据，对90个基因进行评价。'])
r.append(['不能直接改写成促进凋亡','尚不足以证明其促进凋亡'])
P.write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8')
# Keep the full screening provenance in internal files instead of customer narrative.
shutil.copy2(O.parent/'lung_feedback_revision_20261017_v19/tables/R01_selection_flow.tsv',O/'internal_selection_flow.tsv')