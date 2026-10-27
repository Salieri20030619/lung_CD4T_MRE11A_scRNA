import inspect
import pyscenic.aucell as a
import ctxcore.recovery as r
print('CREATE_RANKINGS',inspect.getsource(a.create_rankings));print('AUC4R',inspect.getsource(a.aucell4r));print('AUCS',inspect.getsource(r.aucs));print('AUC_SIGNATURE',inspect.signature(a.aucell))
