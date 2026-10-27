import csv, gzip, os, pathlib, re, shutil, tarfile

root = pathlib.Path(os.environ['DEST'])
with (root / 'filelist.txt').open() as f:
    rows = list(csv.DictReader(f, delimiter='\t'))

expected = {
    r['Name']: int(r['Size'])
    for r in rows if r['#Archive/File'] == 'File'
}
archive = [r for r in rows if r['Name'] == 'GSE274372_RAW.tar']
assert len(archive) == 1 and len(expected) == 36, 'Unexpected GEO manifest'
assert (root / 'GSE274372_RAW.tar').stat().st_size == int(
    archive[0]['Size']
), 'TAR size mismatch'

out = root / 'extracted'
out.mkdir(exist_ok=True)

with tarfile.open(root / 'GSE274372_RAW.tar') as tf:
    members = tf.getmembers()
    files = [m for m in members if m.isfile()]
    assert len(files) == 36
    assert {m.name for m in files} == set(expected), 'TAR member mismatch'
    for m in members:
        assert m.isfile(), 'Unexpected archive entry'
        assert pathlib.PurePosixPath(m.name).name == m.name
        assert m.size == expected[m.name], m.name
        target = out / m.name
        temp = target.with_suffix(target.suffix + '.tmp')
        with tf.extractfile(m) as src, temp.open('wb') as dst:
            shutil.copyfileobj(src, dst)
        temp.replace(target)

manifest = []
for name, size in sorted(expected.items()):
    path = out / name
    assert path.stat().st_size == size, name
    with gzip.open(path, 'rb') as f:
        while f.read(1024 * 1024):
            pass
    match = re.fullmatch(
        r'(GSM\d+)_(\d+\.silicosis\.day(\d+)\.(\d+))\.(.+)\.gz',
        name
    )
    assert match, name
    manifest.append([
        match[1], match[2], match[3], match[4],
        match[5], size, str(path)
    ])

for name in ['GSE274372_cds_object.rds.gz', 'GSE274372_family.soft.gz']:
    with gzip.open(root / name, 'rb') as f:
        while f.read(1024 * 1024):
            pass

with gzip.open(root / 'GSE274372_family.soft.gz', 'rt') as f:
    soft = f.read()
assert '^SERIES = GSE274372' in soft
assert len(re.findall(r'^\^SAMPLE = GSM', soft, re.M)) == 12

with (root / 'download_manifest.tsv').open('w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow([
        'GSM', 'sample', 'day', 'replicate_label',
        'file_type', 'bytes', 'path'
    ])
    writer.writerows(manifest)

with (root / 'download_audit.txt').open('w') as f:
    f.write(
        'SAMPLES=12\nEXPECTED_FILES=36\nAUDITED_PASS=36\n'
        'GZIP_CHECK=PASS\nGSE274372_DOWNLOAD_PASS\n'
    )
print((root / 'download_audit.txt').read_text(), flush=True)
