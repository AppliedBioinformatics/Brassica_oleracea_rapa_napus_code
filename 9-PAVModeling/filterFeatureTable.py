
to_keep = set([x.rstrip() for x in open('all.newPan.prot.Names_With_Darmorv5_Hits.txt')])
with open('Napus_Feature_Table.csv') as fh, open('Napus_Feature_Table.Filtered.csv','w') as out:
    out.write(fh.readline())
    for line in fh:
        ll = line.split()
        if ll[0] in to_keep:
            out.write(line)


to_keep = set([x.rstrip() for x in open('Oleracea_EVM.without_RNASeq_FINAL.No_contaminants_no_small.gff3_prot.fa_with_Oleraceav21_Hits.txt')])
with open('Oleracea_Feature_Table.csv') as fh, open('Oleracea_Feature_Table.Filtered.csv','w') as out:
    out.write(fh.readline())
    for line in fh:
        ll = line.split()
        if ll[0] in to_keep:
            out.write(line)
to_keep = set([x.rstrip() for x in open('Rapa_EVM.without_RNASeq_FINAL.No_contaminants_no_small.gff3_prot_with_Rapa25_30_Hits.txt')])
with open('Rapa_Feature_Table.csv') as fh, open('Rapa_Feature_Table.Filtered.csv','w') as out:
    out.write(fh.readline())
    for line in fh:
        ll = line.split()
        if ll[0] in to_keep:
            out.write(line)
