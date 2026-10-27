from pathlib import Path
from PIL import Image,ImageDraw
O=Path(__file__).parent;R=O/'render';R.mkdir(exist_ok=True);fs=sorted((O/'figures').glob('*.png'))
for j in range(0,len(fs),8):
 sheet=Image.new('RGB',(1400,440*((min(8,len(fs)-j)+1)//2)),'#dfe5e8');dr=ImageDraw.Draw(sheet)
 for k,p in enumerate(fs[j:j+8]):
  im=Image.open(p).convert('RGB');im.thumbnail((690,410));x=k%2*700;y=k//2*440;sheet.paste(im,(x,y+25));dr.text((x+5,y+5),p.stem,fill='black')
 sheet.save(R/('figures_'+str(j//8+1)+'.png'))