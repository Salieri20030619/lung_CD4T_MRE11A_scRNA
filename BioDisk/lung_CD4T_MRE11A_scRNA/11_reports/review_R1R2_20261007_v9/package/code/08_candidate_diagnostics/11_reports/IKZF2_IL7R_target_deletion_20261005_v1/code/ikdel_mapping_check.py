exec(open(__file__.replace('ikdel_mapping_check.py','ikdel_compute.py')).read().split('all_models=[]')[0])
from pyscenic.cli.utils import load_signatures
for seed in [17,29,43]:
 q=net[net.seed==seed];reg=[r for r in load_signatures(str(N/f'seed{seed}/regulons.csv')) if r.name=='IKZF2(+)'][0]
 print('MAPPING',seed,'n',len(q),'unique',q.target.is_unique,'technical_not_in_universe',sorted(set(q.target)-set(genes)),'original_not_in_universe',sorted(set(reg.gene2weight)-set(genes)),'technical_not_in_original',sorted(set(q.target)-set(reg.gene2weight)),flush=True)
print('MAPPING_CHECK_DONE',flush=True)
