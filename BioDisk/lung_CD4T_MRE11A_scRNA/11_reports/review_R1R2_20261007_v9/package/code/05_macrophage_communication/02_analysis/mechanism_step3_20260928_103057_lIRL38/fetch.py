import os, json, hashlib, subprocess, time
from pathlib import Path
out=Path(os.environ['OUT']); cache=Path(os.environ['CACHE']); cache.mkdir(parents=True,exist_ok=True)
prev=Path(os.environ['PREV']); prev2=Path(os.environ['PREV2'])
for p,stage in [(prev,'mechanism_step1'),(prev2,'mechanism_step2_human_and_driver_audit')]:
    s=json.loads((p/'status.json').read_text()); assert s['status']=='PASS' and s['stage']==stage,(p,s)
for line in ['CD4','Mac']:
    for suf in ['counts.mtx.gz','genes.tsv','cells.tsv']:
        assert (prev/'comm_input'/f'{line}_{suf}').is_file()
base=['curl','-L','--fail','--show-error','--silent','--connect-timeout','30','--speed-time','300','--speed-limit','128']
if os.environ.get('PROXY_URL'):base+=['--proxy',os.environ['PROXY_URL'],'--noproxy','']
def get(url,path,resume=True):
    for attempt in range(1,7):
        args=base+(['--continue-at','-'] if resume else [])+['-o',str(path),url]
        rc=subprocess.run(args).returncode
        if rc==0:return
        print(f'DOWNLOAD_RETRY file={path.name} attempt={attempt} curl_exit={rc}',flush=True)
        time.sleep(10)
    raise RuntimeError('Download failed: '+path.name+'; check server proxy/tunnel, then rerun.')
meta=cache/'zenodo_7074291.json'
if not meta.exists():
    tmp=meta.with_suffix('.json.part');get('https://zenodo.org/api/records/7074291',tmp,False)
    j=json.loads(tmp.read_text());assert 'files' in j;tmp.replace(meta)
j=json.loads(meta.read_text()); records={v['key']:v for v in j['files']}
rows=[]
for name in ['lr_network_mouse_21122021.rds','ligand_target_matrix_nsga2r_final_mouse.rds']:
    spec=records[name]; expected=spec['checksum'];algorithm,digest=expected.split(':',1)
    path=cache/name
    def verify(p):
        h=hashlib.new(algorithm)
        with p.open('rb') as f:
            for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
        return h.hexdigest()==digest and p.stat().st_size==int(spec['size'])
    if path.exists():
        if not verify(path):raise RuntimeError('Cached file failed checksum: '+str(path))
    else:
        part=path.with_suffix('.rds.part')
        print('DOWNLOAD_START='+name,flush=True)
        get('https://zenodo.org/records/7074291/files/'+name+'?download=1',part)
        if not verify(part):raise RuntimeError('Downloaded file failed checksum: '+str(part))
        part.replace(path)
    rows.append(dict(file=name,bytes=path.stat().st_size,checksum=expected,status='PASS'))
    print('RESOURCE_PASS='+name,flush=True)
(out/'resource_manifest.json').write_text(json.dumps(rows,indent=2))
print('NICHE_NET_RESOURCE_CHECK=PASS',flush=True)

