files = [['Rapa_EVM.without_RNASeq_FINAL.No_contaminants_no_small.gff3_prot.fa_swissprot_first_hit.Description.long.txt','../Rapa_RenameKey.txt'],
        ['all.newPan.prot.fa_swissprot_first_hit.Description.long.txt', '../Napus_RenameKey.txt'],
        ['Oleracea_EVM.without_RNASeq_FINAL.No_contaminants_no_small.gff3_prot.fa_swissprot_first_hit.Description.long.txt','../Oleracea_RenameKey.txt']]

for a,b in files:
    print(a,b)
    x = dict( line.split() for line in open(b))
    out = open(a.replace('.txt','.BrassicaNamingStandard.txt'), 'w')
    for line in open(a):
        ll = line.split('\t')
        try:
            ll[0] = x[ll[0]]
        except:
            continue
        out.write('\t'.join(ll))

