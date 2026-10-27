from pathlib import Path
B=Path(__file__).parent
exec(compile((B/'revision_document.py').read_text(encoding='utf-8'),str(B/'revision_document.py'),'exec'))
exec(compile((B/'revision_document_tail.py').read_text(encoding='utf-8'),str(B/'revision_document_tail.py'),'exec'))
exec(compile((B/'revision_document_finish.py').read_text(encoding='utf-8'),str(B/'revision_document_finish.py'),'exec'))
