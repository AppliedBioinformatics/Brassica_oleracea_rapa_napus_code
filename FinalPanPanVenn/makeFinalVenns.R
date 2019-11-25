library(tidyverse)
library(eulerr)
library(cowplot)
# Oleracea_PAV_nice_names_intersection_of_individuals.csv
vd <- venn(c("napus&oleracea&rapa"=65186, "napus&oleracea"=275, "napus&rapa"=0, "oleracea&rapa"=91, napus=0, rapa=0, oleracea=88))
oleracea <- plot(vd, labels=c(expression(italic('B. napus')), expression(italic('B. oleracea')), expression(italic('B. rapa'))))
oleracea

ggsave('Oleracea_pangenome.png', plot=oleracea)
# NewPAV_Table_nice_names_intersection_of_individuals.csv
vd <- venn(c("napus&oleracea&rapa"=123714, "napus&oleracea"=964, "napus&rapa"=8024, "oleracea&rapa"=102, napus=17, rapa=49, oleracea=31))
napus <- plot(vd, labels=c(expression(italic('B. napus')), expression(italic('B. oleracea')), expression(italic('B. rapa'))))
napus
ggsave('Napus_pangenome.png', plot=napus)
# Rapa_PAV_nice_names_intersection_of_individuals.csv

vd <- venn(c("napus&oleracea&rapa"=62160, "napus&oleracea"=0, "napus&rapa"=1287, "oleracea&rapa"=43, napus=0, rapa=54, oleracea=0))

rapa <- plot(vd, labels=c(expression(italic('B. napus')), expression(italic('B. oleracea')), expression(italic('B. rapa'))))
rapa
ggsave('Rapa_pangenome.png', plot = rapa)

top_row <- plot_grid(oleracea, rapa, labels=c('A','B'), ncol=2)
bottom_row <- plot_grid(NULL, napus, NULL,ncol=3,
                        rel_widths=c(0.25,0.5,0.25),
                        label_x = 1, labels='C')
grid <- plot_grid(top_row, bottom_row, ncol=1)
save_plot('Napus_oleracea_rapa_grid.png', plot=grid, base_height = 8)


vd <- venn(c("napus&oleracea&rapa"=123714, "napus&oleracea"=964, "napus&rapa"=8024, "oleracea&rapa"=102, napus=17, rapa=103, oleracea=104))
panpan <- plot(vd, labels=c(expression(italic('B. napus')), expression(italic('B. oleracea')), expression(italic('B. rapa'))))
panpan
ggsave('Pan_Pan_Venn.png', plot=panpan)
