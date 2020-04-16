
#core_vs_variable <- data.frame(stringsAsFactors=TRUE,
#      Status = c("Variable", "Variable", "Variable", "Core", "Core", "Core", "Variable",
#                 "Variable", "Variable", "Core", "Core", "Core"),
#     Species = c("B. rapa", "B. napus\n(with syn-\nthetics)", "B. napus\n(no syn-\nthetics)", "B. rapa", "B. napus\n(with syn-\nthetics)", "B. napus\n(no syn-\nthetics)",
#                 "B. oleracea", "B. napus\n(with syn-\nthetics)", "B. napus\n(no syn-\nthetics)", "B. oleracea", "B. napus\n(with syn-\nthetics)", "B. napus\n(no syn-\nthetics)"),
#      Genome = c("A-genome", "A-genome", "A-genome", "A-genome", "A-genome", "A-genome", "C-genome", "C-genome", "C-genome",
#                 "C-genome", "C-genome", "C-genome"),
#       Count = c(21729, 16844, 10447, 42089, 32145, 38542, 14667, 28464, 18233, 50973, 41938, 52169),
#     Perc = c(21729/63818*100, 16844/48989*100, 10447/48989*100, 42089/63818*100, 32145/48989*100, 38542/48989*100,
#              14667/65640*100,28464/70402*100, 18233/70402*100, 50973/65640*100, 41938/70402*100, 52169/70402*100)
#              
#)
library(tidyverse)
#write_tsv(core_vs_variable, 'core_vs_variable.tsv')
core_vs_variable <- read_tsv('core_vs_variable.tsv')
core_vs_variable$Perc
View(core_vs_variable)
library(tidyverse)
totals <- core_vs_variable %>% group_by(Species, Genome) %>% summarise(total = sum(Count))
totals


core_vs_variable
library(ggplot2)
library(cowplot)
theme_set(theme_light())
p <- ggplot(data = core_vs_variable, aes(
  x = factor(Species, levels = c('B. rapa', 'B. oleracea', "B. napus\n(no syn-\nthetics)",'B. napus\n(with syn-\nthetics)')),
  y = Perc,
  group = Status
)) +
  ylab("Percentage core/variable") +
  xlab('Species') +
  #geom_point(size = 3) +
  #geom_line() +
  geom_col(aes(fill=Status)) +
  scale_color_brewer(palette='Dark2') +
  facet_grid( ~ Genome, space = 'free', scale = 'free') +
  theme(axis.text.x = element_text(face = "italic")) 
p

save_plot(p, filename = 'Core_Variable_Genes_Per_A_C_Genome.png', base_aspect_ratio = 1.3, base_height = 5)
