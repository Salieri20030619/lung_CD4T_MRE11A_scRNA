from pathlib import Path
import json,ast
O=Path(r'C:\Users\31937\lung_selection_joint_revision_20261017_v30');p=O/'plots.py';s=p.read_text(encoding='utf-8').replace('top=.91,hspace=.5','top=.85,hspace=.5');p.write_text(s,encoding='utf-8')
p=O/'report.py';s=p.read_text(encoding='utf-8').replace("loc('焦亡/凋亡评分｜焦亡执行与炎症程序分开评价')","loc('焦亡/凋亡评分｜联合分析未支持焦亡相关程序的疾病增强')").replace('程序效应减弱反映结果差异','部分程序效应减弱反映结果差异')
# Short lines for long identifier lists; preserve every identifier.
for line in s.splitlines():
 if line.startswith("p('ARNT2、"):
  val=ast.literal_eval(line[2:-1]);genes=val.split('、');new='\n'.join('、'.join(genes[i:i+7]) for i in range(0,len(genes),7));s=s.replace(line,'p('+repr(new)+')')
p.write_text(s,encoding='utf-8')