from pathlib import Path
from bs4 import BeautifulSoup
import json,re
p=Path('C:/codex_test/tmp/silicosis_v11')
soup=BeautifulSoup((p/'shared_conversation.html').read_text(encoding='utf-8'),'html.parser')
raw=next(t.get_text() for t in soup.find_all('script') if 'streamController.enqueue(' in t.get_text())
offset=raw.index('enqueue(')+len('enqueue(')
payload,_=json.JSONDecoder().raw_decode(raw[offset:]);flat=json.loads(payload)
def key(k):
 return flat[int(k[1:])] if isinstance(k,str) and re.fullmatch(r'_\d+',k) else k
messages=[]
for i,v in enumerate(flat):
 if not isinstance(v,dict):continue
 d={key(k):x for k,x in v.items()}
 if 'content' not in d or 'author' not in d:continue
 content=flat[d['content']] if d['content']>=0 else None
 author=flat[d['author']] if d['author']>=0 else None
 if not isinstance(content,dict) or not isinstance(author,dict):continue
 content={key(k):x for k,x in content.items()}; author={key(k):x for k,x in author.items()}
 role=flat[author['role']] if 'role' in author else 'unknown'
 parts=flat[content['parts']] if 'parts' in content else []
 if not isinstance(parts,list):continue
 texts=[flat[x] for x in parts if isinstance(x,int) and x>=0 and isinstance(flat[x],str)]
 if not texts:continue
 time=flat[d['create_time']] if isinstance(d.get('create_time'),int) and d['create_time']>=0 else None
 mid=flat[d['id']] if isinstance(d.get('id'),int) and d['id']>=0 else str(i)
 messages.append({'flat_index':i,'id':mid,'role':role,'create_time':time,'text':'\n'.join(texts)})
seen=set();unique=[]
for m in messages:
 if m['id'] not in seen: unique.append(m);seen.add(m['id'])
unique.sort(key=lambda x:x['create_time'] or 0)
(p/'shared_messages.json').write_text(json.dumps(unique,ensure_ascii=False,indent=2),encoding='utf-8')
hits=[m for m in unique if any(g in m['text'].upper() for g in ['PENK','IL1R2','TTC39C','GM2663','TNFRSF8']) and any(k in m['text'] for k in ['筛选','创新','成熟','排名','文献','优先'])]
(p/'screening_history_hits.json').write_text(json.dumps(hits,ensure_ascii=False,indent=2),encoding='utf-8')
print('flat entries',len(flat),'messages',len(unique),'screening hits',len(hits))
for m in hits:print(m['id'],m['create_time'],m['role'],len(m['text']),m['text'][:180].replace('\n',' '))
