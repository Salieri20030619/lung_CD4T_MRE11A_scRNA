from pathlib import Path
import tokenize,io,ast
p=Path(r'C:\Users\31937\lung_selection_joint_revision_20261017_v30/report.py');s=p.read_text(encoding='utf-8');out=[]
for t in tokenize.generate_tokens(io.StringIO(s).readline):
 if t.type==tokenize.STRING:
  try:
   v=ast.literal_eval(t.string)
   if isinstance(v,str) and v.startswith('ARNT2、'):
    genes=v.split('、');v='\n'.join('、'.join(genes[i:i+6]) for i in range(0,len(genes),6));t=t._replace(string=repr(v))
  except (ValueError,SyntaxError):pass
 out.append(t)
p.write_text(tokenize.untokenize(out),encoding='utf-8')