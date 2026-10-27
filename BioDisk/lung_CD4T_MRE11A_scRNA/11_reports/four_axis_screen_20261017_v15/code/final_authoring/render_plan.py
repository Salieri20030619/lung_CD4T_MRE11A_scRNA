from pathlib import Path
import fitz
W=Path('C:/Users/31937/lung_four_axes_20261017_v15');T=W/'plan_render';T.mkdir(exist_ok=True);d=fitz.open(W/'deliverables/four_axis_frozen_plan_v15.pdf');print('PAGES',len(d))
for i,p in enumerate(d):
 assert '\ufffd' not in p.get_text();p.get_pixmap().save(T/f'page_{i+1}.png')