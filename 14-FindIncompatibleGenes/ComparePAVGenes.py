
import shelve

inds = []
variable_genes = set()
pav_dict = {}
with open('./NewPAV_Table_nice_names_intersection_of_individuals_Napus_Only.csv') as fh:
    header = fh.readline().split()
    genes = header[1:]
    for g in genes:
        pav_dict[g] = []
    for line in fh:
        ll = line.split()
        ind = ll[0]
        inds.append(ind)
        alleles = ll[1:]
        for a, b in zip(genes, alleles):
            pav_dict[a].append(int(b))
            if b == '0':
                variable_genes.add(a)

newpav_dict = { your_key: pav_dict[your_key] for your_key in variable_genes } 

import pandas as pd
#df = pd.DataFrame.from_dict(newpav_dict)

gene_order = []
keep = set()
for line in open('Darmor_v1.4_Gene_Order.txt'):
    ll = line.split('\t')
    names = ll[-1].split(';')[0].replace('ID=','')
    if names in variable_genes:# and ll[0] in set(['chrA01_RaGOO','chrC01_RaGOO']):
        keep.add(names)
        gene_order.append(names)

newpav_dict = { your_key: newpav_dict[your_key] for your_key in keep } 

#df = df[gene_order]
#df = df.iloc[0:1000,0:1000]

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#f = plt.figure(figsize=(19, 15))
#c = df.corr()
#c = np.corrcoef(df.transpose())
#print(c)
#c = pd.DataFrame(data = c)
#print(c)
#s = c.unstack()
#so = s.sort_values(kind="quicksort")

##print(so[:100])
#print('-------')
#print(so[-101:])
#plt.matshow(c, fignum=f.number)
#plt.xticks(range(df.shape[1]), df.columns, fontsize=14, rotation=45)
#plt.yticks(range(df.shape[1]), df.columns, fontsize=14)
#cb = plt.colorbar()
#cb.ax.tick_params(labelsize=14)
#plt.title('Correlation Matrix', fontsize=16);
#plt.savefig('hi.png')

#pair_different_dict= {} # for each pair, the score of their difference

#pair_different_dict = set()
#pair_different_dict = shelve.open('myset')

second_newpav_dict=  newpav_dict.copy()
for g in newpav_dict:

    #for other_g in newpav_dict:
    for other_g in second_newpav_dict:
        if g == other_g:
            continue
        this_pair = ' '.join(tuple(sorted([g, other_g])))
        g_alleles, other_g_alleles = newpav_dict[g], newpav_dict[other_g]
        shared_zero = 0
        shared_one = 0
        different = 0
        a_zeros = 0
        b_zeros = 0
        for a,b in zip(g_alleles, other_g_alleles):
            if a + b == 0:
                #shared 0
                shared_zero += 1
            if a + b == 2:
                # shared 1
                shared_one += 1
            if a + b == 1:
                different += 1
            if a == 0:
                a_zeros += 1
            if b == 0:
                b_zeros += 1
        #pair_different_dict[this_pair] = (float(diff_score) / total_score , float(same_score) / total_score)
        #pair_different_dict.add(this_pair)
        #pair_different_dict[this_pair] = 1
        print(this_pair, shared_zero, shared_one, different, a_zeros, b_zeros)
    del second_newpav_dict[g]
