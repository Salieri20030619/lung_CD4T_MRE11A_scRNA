import subprocess,sys,pathlib
p=subprocess.run(['ssh','-T','-o','BatchMode=yes','sgs10','/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python','-'],input=pathlib.Path(sys.argv[1]).read_bytes(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,creationflags=subprocess.CREATE_NO_WINDOW)
sys.stdout.buffer.write(p.stdout);sys.stderr.buffer.write(p.stderr);sys.exit(p.returncode)