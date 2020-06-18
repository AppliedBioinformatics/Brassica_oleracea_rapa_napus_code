
from glob import glob

rapa_to_skip = set([x.rstrip().split()[0] for x in open('Rapa_Garbage_Genes_Remove_Us.txt')])

for g in glob('*.csv'):
    if '.filtered.' in g: continue
    if g == 'Oleracea_PAV_nice_names_intersection_of_individuals.csv':
        hits = 'Oleracea_EVM.without_RNASeq_FINAL.No_contaminants_no_small.gff3_prot.fa_with_Oleraceav21_Hits.txt'
    elif g == 'Rapa_PAV_nice_names_intersection_of_individuals.csv':
        hits = 'Rapa_EVM.without_RNASeq_FINAL.No_contaminants_no_small.gff3_prot_with_Rapa25_30_Hits.txt'
    elif g == 'NewPAV_Table_nice_names_intersection_of_individuals.csv':
        hits = 'all.newPan.prot.Names_With_Darmorv5_Hits.txt'
    to_keep = set()
    for line in open(hits):
        ll = line.split()
        assert ll[0] not in to_keep
        to_keep.add(ll[0])
    with open(g) as fh, open(g.replace('.csv','.filtered.csv'), 'w') as out:
        header = fh.readline().split()
        ind = header[0]
        genes = header[1:]
        indexes_to_skip = set()
        newheader = [ind]
        for index, ge in enumerate(genes):
            if ge not in to_keep and ('jcf' not in ge) and ('harsh' not in ge):
                indexes_to_skip.add(index)
            elif (g == 'Rapa_PAV_nice_names_intersection_of_individuals.csv' and ge in rapa_to_skip):
                indexes_to_skip.add(index)
                print('skipping %s'%ge)
            else:
                newheader.append(ge)
        out.write('\t'.join(newheader) + '\n')
        for line in fh:
            ll = line.split()
            ind, alleles = ll[0], ll[1:]
            newll = [ind]
            for index, a in enumerate(alleles):
                if index in indexes_to_skip:
                    continue
                newll.append(a)
            out.write('\t'.join(newll) + '\n')
