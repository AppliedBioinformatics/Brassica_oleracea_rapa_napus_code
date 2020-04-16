
library(VennDiagram)
#all_genes 132913
#all_three 123816
#lost_in_oler 8073
#            oler_chrC = 340
#            oler_chrA = 6330
#            oler_pan = 529
#            oler_unpl =874
#lost_in_rapa 995
#            rapa_chrC = 754
#            rapa_chrA = 15
#            rapa_pan = 211
#            rapa_unpl =15
#lost_in_napus 0
#only_napus 28
#only_oler 0
#only_rapa 0

png('Venn_napus.png', height=12, width=12, unit='cm', res=300)
source('Venn_out.txt')

venn.diagram(myvenn, filename = 'Venn_other_napus.png')
dev.off()
