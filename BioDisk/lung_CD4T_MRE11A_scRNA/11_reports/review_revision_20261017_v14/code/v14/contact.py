from pathlib import Path
from PIL import Image,ImageOps,ImageDraw
W=Path('C:/Users/31937/lung_review_20261017_v14')
for name,files in [('main',sorted((W/'render').glob('page_*.png'))),('response',sorted((W/'render_response').glob('response_*.png')))]:
 for start in range(0,len(files),6):
  canvas=Image.new('RGB',(1684,1845),'#dbe3e8')
  for i,f in enumerate(files[start:start+6]):
   im=Image.open(f).convert('RGB');x=(i%2)*842;y=(i//2)*615;canvas.paste(im,(x,y+20));ImageDraw.Draw(canvas).text((x+10,y+3),f.name,fill='black')
  canvas.save(W/f'qa_{name}_{start//6+1}.png')