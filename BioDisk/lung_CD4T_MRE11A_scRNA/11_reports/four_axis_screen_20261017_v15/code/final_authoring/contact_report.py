from pathlib import Path
from PIL import Image,ImageDraw
W=Path('C:/Users/31937/lung_four_axes_20261017_v15');files=sorted((W/'report_render').glob('page_*.png'))
for start in range(0,len(files),4):
 out=Image.new('RGB',(1684,1230),'#dbe3e8')
 for j,f in enumerate(files[start:start+4]):
  x=j%2*842;y=j//2*615;out.paste(Image.open(f).convert('RGB'),(x,y+20));ImageDraw.Draw(out).text((x+8,y+3),f.name,fill='black')
 out.save(W/f'qa_report_{start//4+1}.png')