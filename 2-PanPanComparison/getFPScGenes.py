
from Bio import SeqIO
toget = set([x.rstrip().split()[1] for x in open('Rapa_Pangenome_OnlyFPSc_Variable_GeneNames.txt')])

for s in SeqIO.parse('Rapa_EVM.without_RNASeq_FINAL.No_contaminants_no_small.gff3_prot.fa', 'fasta'):
    if s.id in toget:
        print(s.format('fasta'))
