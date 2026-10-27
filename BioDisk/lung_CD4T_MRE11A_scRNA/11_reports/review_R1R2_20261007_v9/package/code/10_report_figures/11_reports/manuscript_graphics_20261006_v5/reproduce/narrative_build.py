from pathlib import Path
for part in ['narrative_document.py','narrative_results.py']:
 exec(compile(Path(part).read_text(encoding='utf8'),part,'exec'))
supp=Path('manuscript_supplement.py').read_text(encoding='utf8').replace('方法1.6','方法1.9').replace('不仅本表7程序','不仅所展示的7个程序')
exec(compile(supp,'manuscript_supplement.py','exec'))
exec(compile(Path('manuscript_finish.py').read_text(encoding='utf8'),'manuscript_finish.py','exec'))
