import sys
import numpy as np
 
g_d={}
inds = int(sys.argv[2])
for l in open(sys.argv[1]):
    l_arr=l.rstrip().split("\t")
    if(l_arr[0]=="Individual"):
        genes=l_arr[1:]
        cnt=0
    else:
        g_d[cnt]=set()
        gl=l_arr[1:]
        for x in range(0, len(genes)):
            if(int(gl[x])==1):
                g_d[cnt].add(genes[x])
            x=x+1
        cnt=cnt+1
 
c_d={}
for i in range(2, inds):
    c_d[i]=set()
    for b in range(0, inds):
        c_d[i].add(",".join(str(k) for k in sorted(np.random.choice(inds, i, replace=False))))
 
 
print('Combinations\tcore\tpangenome\tvariable')
for y in c_d:
    for z in c_d[y]:
        sl=[]
        for q in z.split(","):
            sl.append(g_d[int(q)])
        combinations = y
        core = len(set.intersection(*sl))
        pangenome = len(set.union(*sl))
        variable = pangenome - core
        #print("CORE"+"\t"+str(y)+"\t"+str(len(set.intersection(*sl))))
        #print("PANGENOME"+"\t"+str(y)+"\t"+str(len(set.union(*sl))))
        ll = map(str, [combinations, core, pangenome, variable])
        print('\t'.join(ll))
 
