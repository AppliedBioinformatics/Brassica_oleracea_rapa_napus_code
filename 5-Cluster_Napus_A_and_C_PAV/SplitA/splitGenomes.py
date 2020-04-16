fh = open('../../1-FilterPAVTables/NewPAV_Table_nice_names_intersection_of_individuals.filtered.csv')

header = fh.readline().split()
indexes_dict = {} #key: chrom, value: indexes (columns) for this chrom
for index, h in enumerate(header):
    if h == 'Individual':
        continue
    hh = h.split('.')
    chrom = hh[2]
    if 'jcf' in chrom or 'harsh' in chrom:
        chrom = 'pangenome'
    if chrom not in indexes_dict:
        indexes_dict[chrom] = [index]
    else:
        indexes_dict[chrom].append(index)

fh_dict = {}
for i in indexes_dict:
    outfh = open('NewPAV_Table_nice_names_intersection_of_individuals.filtered_Split_%s.csv'%i, 'w')
    fh_dict[i] = outfh
    this_indexes = indexes_dict[i]
    newheader = ['Individual']
    for x in this_indexes:
        newheader.append(header[x])
    outfh.write('\t'.join(newheader) + '\n')

for line in fh:
    ll = line.split()
    for i in indexes_dict:
        outfh = fh_dict[i]
        newll = [ll[0]]
        this_indexes = indexes_dict[i]
        for x in this_indexes:
            newll.append(ll[x])
        outfh.write('\t'.join(newll) + '\n')
