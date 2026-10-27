from pathlib import Path
s=Path('manuscript_qa.py').read_text(encoding='utf8').replace("D=Path('manuscript_delivery')","D=Path('narrative_delivery')")
exec(compile(s,'manuscript_qa.py','exec'))
