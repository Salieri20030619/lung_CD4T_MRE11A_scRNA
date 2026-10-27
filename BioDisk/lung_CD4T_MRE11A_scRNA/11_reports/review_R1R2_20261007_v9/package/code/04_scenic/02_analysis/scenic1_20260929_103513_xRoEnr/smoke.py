import sys, subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from pyscenic.cli.utils import load_signatures
from ctxcore.rnkdb import FeatherRankingDatabase
import loompy
p=Path(sys.argv[1]);p.mkdir(exist_ok=True)
rng=np.random.default_rng(17);a=rng.poisson(2,size=(80,32)).astype(float)
a[:,10]=a[:,0]*2+rng.poisson(1,size=80)
names=[f'T{i}' for i in range(32)]
loompy.create(str(p/'smoke.loom'),a.T,{'Gene':np.array(names)},{'CellID':np.array([f'C{i}' for i in range(80)])})
(p/'tf.txt').write_text('\n'.join(names[:8])+'\n')
subprocess.run([sys.executable,'-m','pyscenic.cli.arboreto_with_multiprocessing',str(p/'smoke.loom'),str(p/'tf.txt'),'--num_workers','2','--seed','17','-o',str(p/'adj.tsv')],check=True)
d=pd.read_csv(p/'adj.tsv',sep='\t');assert {'TF','target','importance'}<=set(d.columns) and len(d)>0 and np.isfinite(d.importance).all()
for command in ['ctx','aucell']:subprocess.run([str(Path(sys.executable).parent/'pyscenic'),command,'--help'],stdout=subprocess.DEVNULL,check=True)
print('ENVIRONMENT_SMOKE_PASS',flush=True)

