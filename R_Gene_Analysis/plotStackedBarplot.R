library(tidyverse)
library(cowplot)
df <- readxl::read_xlsx('./R_gene_stats.xlsx')
df
df$NLR <- NULL
long_df <- df %>% gather(Type, Count, CN:RLP_lysm)
df
long_df$Genome <- factor(long_df$Genome, levels=c('B. rapa','B. napus A', 'B. oleracea', 'B. napus C', 'B. rapa new contigs','B. oleracea new contigs', 'B. napus new contigs'))

long_df


rlps <- long_df %>% filter(str_detect(Type, '^RLP'))
rlks <- long_df %>% filter(str_detect(Type, '^RLK'))
rlps_rlks <- long_df %>% filter(str_detect(Type, '^RLP|^RLK'))
nbs <- long_df %>% filter(! str_detect(Type, '^RLP|^RLK'))
rlps


a <- ggplot(nbs, aes(fill=Type,x=Genome, y =Count)) + geom_bar(stat='identity') +
  scale_fill_brewer(palette = "Dark2") + coord_flip() + 
  theme(axis.text.y = element_text(face = "italic")) +
  scale_x_discrete(limits = rev(levels(nbs$Genome)))
a
save_plot(plot =a, 'R_genes_NBS.png', base_height = 6, base_aspect_ratio = 1.3)


ggplot(rlps, aes(fill=Type,x=Genome, y =Count)) + geom_bar(stat='identity') +
  scale_fill_brewer(palette = "Dark2") + coord_flip() + 
  theme(axis.text.y = element_text(face = "italic")) +
  scale_x_discrete(limits = rev(levels(rlps$Genome)))

ggplot(rlks, aes(fill=Type,x=Genome, y =Count)) + geom_bar(stat='identity') +
  scale_fill_brewer(palette = "Dark2") + coord_flip() + 
  theme(axis.text.y = element_text(face = "italic")) +
  scale_x_discrete(limits = rev(levels(rlks$Genome)))


b <- ggplot(rlps_rlks, aes(fill=Type,x=Genome, y =Count)) + geom_bar(stat='identity') +
  scale_fill_brewer(palette = "Dark2") + coord_flip() + 
  theme(axis.text.y = element_text(face = "italic")) +
  scale_x_discrete(limits = rev(levels(rlps_rlks$Genome)))


both <- plot_grid(a, b, ncol=1, labels='AUTO', align='v')
save_plot(plot =both, 'R_genes.png', base_height = 6, base_aspect_ratio = 1.3)
