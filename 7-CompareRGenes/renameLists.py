keys = {}

import glob
for k in glob.glob('*Key*'):
    name = k.split('_')[0]
    keys[name] = {}
    for line in open(k):
        old, new = line.rstrip().split()
        keys[name][old] = new

for k in glob.glob('*lst'):
    if 'renamed' in k:
        continue
    with open(k.replace('.lst', '_renamed.lst'), 'w') as out:
        name = k.split('.')[0].split('_')[0]
        for line in open(k):
            ll = line.split()
            old = ll[0]
            try:
                ll[0] = keys[name][old]
            except KeyError:
                # these are filtered genes not in other cultivars
                continue
            out.write('\t'.join(ll) + '\n')

