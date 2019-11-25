# from http://appliedbioinformatics.com.au/wiki/index.php/PAV_PCA_and_clustering

#####PCA using logisticPCA######
library(logisticPCA)
library(ggfortify)
library(calibrate)
library(data.table)
library(ggplot2)
library(cowplot)
library(tidyverse)
theme_set(theme_cowplot())

# load in the sample IDs of the special subgroups
synthetic_lines <- read_csv('synthetic_lines.txt', col_names=FALSE)
synthetic_lines
fpsc_lines <- read_csv('FPSc_samples.txt', col_names=FALSE)

##Load in entire genome
all_st <- fread('NewPAV_Table.tsv', sep='\t', header=T)

# Now subsample for the C-genome
# keep only napus/oleracea individuals
all_st.c <- all_st %>% filter(!grepl('Brassica_rapa', Individual))

names(all_st.c)[1:5]

# pull out column of individual names
sample_type <- all_st.c$Individual

# clean up column
sample_type <- sample_type %>% str_replace_all('/scratch/pawsey0149/pbayer/DOwnloadPanPanData/fastq/', '') %>% 
  str_replace_all('_MERGED.bam_mosdepth.per-base.bed.gz_vs_exons.bed', '') %>% 
  str_replace_all('_R1_001.fastq_newPan_sorted', '') %>% 
  str_replace_all('_L005', '') %>% 
  str_replace_all('_L006', '') %>% 
  str_replace_all('_L007', '') %>% 
  str_replace_all('_L008', '')
sample_type
# make data frame from column
sample_type <- data.frame(Individual=sample_type, stringsAsFactors=FALSE)

synth_check <- paste0(synthetic_lines$X1, collapse='|')
# detect species
sample_species <- sample_type %>% mutate(type = ifelse(
  str_detect(Individual, synth_check),
  'Synth',
  ifelse(
  str_detect(Individual, 'Brassica_rapa'),
  'B_rapa',
  ifelse(
    str_detect(Individual, 'Brassica_oleracea'),
    'B_oleracea',
    ifelse(
      str_detect(Individual, 'Brassica_napus'),
      'B_napus',
      ifelse(str_detect(Individual, 'Harsh'),
             'B_napus',
             'NOPE')
    )
  )
)
)) %>% select(type)

# double check species, there should be no NOPE species left
sample_species %>% table()

# pull out only C genome genes
all_st.c <- all_st.c %>%  select(starts_with('evm.model.chrC')) 

#convert to MATRIX, run SVD
all_st.c_t.m <- as.matrix(sapply(all_st.c, as.numeric))
allc_logsvd_model <- logisticSVD(all_st.c_t.m, k = 2, quiet = F, partial_decomp = TRUE)
sample_names
# PLOT
p <- plot(allc_logsvd_model, type = 'scores') +
  geom_point(aes(colour = factor(
    sample_species$type,
    labels = c('B. napus',
               'B. oleracea',
               'Synthetic\nB. napus')
  ))) +
  labs(colour = "Species") +
  theme(legend.text = element_text(face = 'italic')) +
  guides(colour=guide_legend(keywidth=0.1,
                             keyheight = 0.35,
                             default.unit='inch'))

save_plot('C_pav.pca.log.png', p,
          base_height = 5,
          base_aspect_ratio = 1.3) # make room for figure legend

# Then for only the C-chromosomes
# First, run the PCA for each subgroup in teh C genome
filenames <-
  dir('SplitC', pattern = "Napus_PAV_C_R_Format_Split.*csv", full.names =
        TRUE)

c_big_pc <- data.frame() # stores chromosome, pc1, pc2
all_pcs <- data.frame()

for (i in 1:length(filenames)) {
  thisfile <- filenames[i]
  
  c_t <- fread(thisfile, sep = '\t', header = T)
  c_t.m <- as.matrix(sapply(c_t, as.numeric))
  
  rownames(c_t.m) <- sample_type$Individual
  
  c_logsvd_model <- logisticSVD(c_t.m, k = 2, quiet = F, partial_decomp = TRUE)

  p <- plot(c_logsvd_model, type = 'scores') +
    geom_point(aes(colour = factor(
      sample_species$type,
      labels = c('B. napus',
                 'B. oleracea',
                 'Synthetic\nB. napus')
    ))) +
    labs(colour = "Species") +
    theme(legend.text = element_text(face = 'italic')) +
    guides(colour = guide_legend(
      keywidth = 0.1,
      keyheight = 0.35,
      default.unit = 'inch'
    ))
  
  save_plot(paste(thisfile, '.png', sep = ''), p,
            base_height = 5,
            base_aspect_ratio = 1.3) # make room for figure legend

  pc1 <- c_logsvd_model$A[,1]
  all_pcs <- rbind(all_pcs, as.data.frame(t(pc1)))
  
  thisframe <- as.data.frame(cbind(thisfile, c_logsvd_model$A, sample_species$type))
  c_big_pc <- bind_rows(big_pc, thisframe)
}

# now run animations!
# devtools::install_github('thomasp85/gganimate') 
library(gganimate)

names(c_big_pc) <- c('chrom','PC1','PC2','sample_type')


c_big_pc$PC1 <- as.double(c_big_pc$PC1)
c_big_pc$PC2 <- as.double(c_big_pc$PC2)
head(c_big_pc)
c_big_pc$chrom <- str_replace(c_big_pc$chrom, 'SplitA/Napus_PAV_A_R_Format_Split_','') %>% 
  str_replace('.csv','')

p <- ggplot(data=c_big_pc, aes(x=PC1, y=PC2, colour = factor(
  sample_type,
  labels = c('B. napus',
             'B. oleracea',
             'Synthetic\nB. napus')
))) +
  geom_point() +
  labs(colour = 'Species',
       title= 'Chrom: {closest_state}') +
  theme(legend.text = element_text(face = 'italic')) +
  guides(colour = guide_legend(
    keywidth = 0.1,
    keyheight = 0.35,
    default.unit = 'inch'
  )) +
  transition_states(chrom, 3, 3) +
  enter_fade()   +
  exit_shrink() + 
  ease_aes('sine-in-out') +
  theme_cowplot(12)

anim <- animate(p)
anim_save('C_genome_animation.gif', anim)

# Now, run a by-chromosome PCA
all_pcs <- as.matrix(all_pcs)

rownames(all_pcs) <- c('C1','C2','C3','C4','C5','C6','C7','C8','C9')

mine <- prcomp(all_pcs)
p <- autoplot(mine, shape=FALSE)
save_plot('PCA_of_PC1_of_C_chroms.png', p,
          base_height = 5,
          base_aspect_ratio = 1.3)


# Now plot for the A-genome only

# keep only napus/rapa individuals
all_st.a <- all_st %>% filter(!grepl('Brassica_oleracea', Individual))

names(all_st.a)[1:5]

# pull out column of individual names
sample_type <- all_st.a$Individual
sample_type
# clean up column
sample_type <- sample_type %>% str_replace_all('/scratch/pawsey0149/pbayer/DOwnloadPanPanData/fastq/', '') %>% 
  str_replace_all('_MERGED.bam_mosdepth.per-base.bed.gz_vs_exons.bed', '') %>% 
  str_replace_all('_R1_001.fastq_newPan_sorted', '') %>% 
  str_replace_all('_L005', '') %>% 
  str_replace_all('_L006', '') %>% 
  str_replace_all('_L007', '') %>% 
  str_replace_all('_L008', '')

# make data frame from column
sample_type <- data.frame(Individual=sample_type, stringsAsFactors=FALSE)

synth_check <- paste0(synthetic_lines$X1, collapse='|')
fpsc_check <- paste0(fpsc_lines$X1, collapse='|')

# detect species
sample_species <- sample_type %>% mutate(type = ifelse(
  str_detect(Individual, synth_check),
  'Synth',
  ifelse(
    str_detect(Individual, fpsc_check),
    'FPSc',
    ifelse(
    str_detect(Individual, 'Brassica_rapa'),
    'B_rapa',
    ifelse(
      str_detect(Individual, 'Brassica_oleracea'),
      'B_oleracea',
      ifelse(
        str_detect(Individual, 'Brassica_napus'),
        'B_napus',
        ifelse(str_detect(Individual, 'Harsh'),
               'B_napus',
               'NOPE')
      )
    )
  )
)
)) %>% select(type)

# double check species
sample_species %>% table()

# pull out only C genome genes
all_st.a <- all_st.a %>%  select(starts_with('evm.model.chrA')) 

#convert to MATRIX, run SVD
all_st.a_t.m <- as.matrix(sapply(all_st.a, as.numeric))
alla_logsvd_model <- logisticSVD(all_st.a_t.m, k = 2, quiet = F, partial_decomp = TRUE)

# PLOT
p <- plot(alla_logsvd_model, type = 'scores') +
  geom_point(aes(colour = 
    factor(sample_species$type,
    labels = c('B. napus',
               'B. rapa',
               'B. rapa FPSc',
               'Synthetic\nB. napus')
    ))) +
  labs(colour = "Species") +
  theme(legend.text = element_text(face = 'italic')) +
  guides(colour=guide_legend(keywidth=0.1,
                             keyheight = 0.35,
                             default.unit='inch'))
p

save_plot('A_pav.pca.log.png', p,
          base_height = 5,
          base_aspect_ratio = 1.3) # make room for figure legend

# now do the same for the A genome split
# First, run the PCA for each subgroup in teh C genome
filenames <-
  dir('SplitA', pattern = "Napus_PAV_A_R_Format_Split*", full.names = TRUE)
filenames
big_pc <- data.frame() # stores chromosome, pc1, pc2
all_pcs <- data.frame() # stores only pc1

for (i in 1:length(filenames)) {
  thisfile <- filenames[i]
  
  a_t <- fread(thisfile, sep = '\t', header = T)
  a_t.m <- as.matrix(sapply(a_t, as.numeric))
  
  rownames(a_t.m) <- sample_type$Individual
  
  a_logsvd_model <- logisticSVD(a_t.m, k = 2, quiet = F, partial_decomp = TRUE)
  
  p <- plot(a_logsvd_model, type = 'scores') +
    geom_point(aes(colour = factor(
      sample_species$type,
      labels = c('B. napus',
                 'B. rapa',
                 'B. rapa FPSc',
                 'Synthetic\nB. napus')
    ))) +
    labs(colour = "Species") +
    theme(legend.text = element_text(face = 'italic')) +
    guides(colour = guide_legend(
      keywidth = 0.1,
      keyheight = 0.35,
      default.unit = 'inch'
    ))
  
  save_plot(paste(thisfile, '.png', sep = ''), p,
            base_height = 5,
            base_aspect_ratio = 1.3) # make room for figure legend
  
  pc1 <- a_logsvd_model$A[,1]
  all_pcs <- rbind(all_pcs, as.data.frame(t(pc1)))
  
  thisframe <- as.data.frame(cbind(thisfile, a_logsvd_model$A, sample_species$type))
  big_pc <- bind_rows(big_pc, thisframe)
}
# now run animations!
# devtools::install_github('thomasp85/gganimate') 
library(gganimate)

names(big_pc) <- c('chrom','PC1','PC2','sample_type')


big_pc$PC1 <- as.double(big_pc$PC1)
big_pc$PC2 <- as.double(big_pc$PC2)
head(big_pc)
big_pc$chrom <- str_replace(big_pc$chrom, 'SplitA/Napus_PAV_A_R_Format_Split_','') %>% 
  str_replace('.csv','')
big_pc

p <- ggplot(data=big_pc, aes(x=PC1, y=PC2, colour = factor(
  sample_type,
  labels = c('B. napus',
             'B. rapa',
             'B. rapa FPSc',
             'Synthetic\nB. napus')
))) +
  geom_point() +
  labs(colour = 'Species',
       title= 'Chrom: {closest_state}') +
  theme(legend.text = element_text(face = 'italic')) +
  guides(colour = guide_legend(
    keywidth = 0.1,
    keyheight = 0.35,
    default.unit = 'inch'
  )) +
  transition_states(chrom, 3, 3) +
  enter_fade()   +
  exit_shrink() + 
  ease_aes('sine-in-out') +
  theme_cowplot(12)
  
anim <- animate(p)
anim_save('A_genome_animation.gif', anim)

# Now, run a by-chromosome PCA
all_pcs <- as.matrix(all_pcs)

rownames(all_pcs) <- c('A01','A02','A03','A04','A05','A06','A07','A08','A09', 'A10')
getwd()
mine <- prcomp(all_pcs)
p <- autoplot(mine, shape=FALSE)
save_plot('PCA_of_PC1_of_A_chroms.png', p,
          base_height = 5,
          base_aspect_ratio = 1.3)
p
