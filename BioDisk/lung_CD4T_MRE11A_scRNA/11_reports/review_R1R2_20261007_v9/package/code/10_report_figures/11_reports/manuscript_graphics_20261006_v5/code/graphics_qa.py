from pathlib import Path
src=Path('manuscript_qa.py').read_text(encoding='utf8').replace("D=Path('manuscript_delivery')","D=Path('graphics_delivery')")
exec(compile(src,'graphics_qa','exec'))
src=Path('manuscript_numeric_audit.py').read_text(encoding='utf8').replace("D=Path('manuscript_delivery')","D=Path('graphics_delivery')")
exec(compile(src,'graphics_numeric_audit','exec'))
