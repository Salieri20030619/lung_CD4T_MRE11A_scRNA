from pathlib import Path
for part in ['manuscript_document.py','manuscript_results.py','manuscript_supplement.py','manuscript_finish.py']:
 exec(compile(Path(part).read_text(encoding='utf8'),part,'exec'))
