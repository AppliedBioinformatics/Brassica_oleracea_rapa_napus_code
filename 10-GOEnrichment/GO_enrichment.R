library("topGO")
library(wordcloud)
library(tidyverse)
writeTopGO <- function(background, foreground, GOterm, top, txt){
  # give properly formatted background in format: GO:0005838	GSBRNA2T00088508001;GSBRNA2T00088313001;GSBRNA2T00035842001 
  annAT <- readMappings(background, sep="\t", IDsep=";")
  allgenes <- unique(unlist(annAT))
  # give file with your genes of interest, one gene_id per line
  mygenes <-scan(foreground ,what="")
  geneList <- factor(as.integer(allgenes %in% mygenes))
  names(geneList) <- allgenes
  GOdata <-new ("topGOdata", ontology = GOterm, allGenes = geneList, nodeSize = top, annot=annFUN.GO2genes, GO2genes=annAT)
  
  # using ClassicCount:
  #test.stat <-new ("classicCount", testStatistic = GOFisherTest, name = "Fisher Test")
  #resultsFisherC <-getSigGroups (GOdata, test.stat)
  
  # using weight01:
  weight01.fisher <- runTest(GOdata, statistic = "fisher")
  
  # using ClassicCount:
  #allRes <- GenTable(GOdata, classicFisher= resultsFisherC, topNodes = 30)
  
  # using weight01:
  allRes <- GenTable(GOdata, classicFisher=weight01.fisher,topNodes=30)
  
  names(allRes)[length(allRes)] <- "p.value"
  
  write.table(allRes, file=txt, sep="\t",quote=FALSE,row.names=FALSE)
}

#different node sizes give slightly different results, 5 or 10 are recommened as most meaningful by the author
#Background_GO_Darmor_noAED1.txt
#Foreground_GO_Darmor_noAED1.txt

# napus
writeTopGO("all.newPan.prot.fa_swissprot_first_hit.mergedFinal.OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Napus_Pangenome_With_Synths_Variable_GeneNames.txt", "BP", 5, "Darmor_noAED1_5topBPweight01.txt")
writeTopGO("all.newPan.prot.fa_swissprot_first_hit.mergedFinal.OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Napus_Pangenome_With_Synths_Variable_GeneNames.txt", "MF", 5, "Darmor_noAED1_5topMFweight01.txt")
writeTopGO("all.newPan.prot.fa_swissprot_first_hit.mergedFinal.OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Napus_Pangenome_With_Synths_Variable_GeneNames.txt", "CC", 5, "Darmor_noAED1_5topCCweight01.txt")

writeTopGO("Oleracea_EVM_swissprot_hfirst_hit_OnlyPlantHits.txt", "Oleracea_Pangenome_Variable_GeneNames.txt", "BP", 5, "Oleracea_Pangenome_5topBPweight01.txt")
writeTopGO("Oleracea_EVM_swissprot_hfirst_hit_OnlyPlantHits.txt", "Oleracea_Pangenome_Variable_GeneNames.txt", "MF", 5, "Oleracea_Pangenome_5topMFweight01.txt")
writeTopGO("Oleracea_EVM_swissprot_hfirst_hit_OnlyPlantHits.txt", "Oleracea_Pangenome_Variable_GeneNames.txt", "CC", 5, "Oleracea_Pangenome_5topCCweight01.txt")

writeTopGO("Rapa_EVM_swissprot_frist_hit_OnlyPlantHits.txt", "Rapa_Pangenome_WithFPSc_Variable_GeneNames.txt", "BP", 5, "Rapa_Pangenome_5topBPweight01.txt")
writeTopGO("Rapa_EVM_swissprot_frist_hit_OnlyPlantHits.txt", "Rapa_Pangenome_WithFPSc_Variable_GeneNames.txt", "MF", 5, "Rapa_Pangenome_5topMFweight01.txt")
writeTopGO("Rapa_EVM_swissprot_frist_hit_OnlyPlantHits.txt", "Rapa_Pangenome_WithFPSc_Variable_GeneNames.txt", "CC", 5, "Rapa_Pangenome_5topCCweight01.txt")


# let's try 10 too

writeTopGO("all.newPan.prot.fa_swissprot_first_hit.mergedFinal.OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Napus_Pangenome_With_Synths_Variable_GeneNames.txt", "BP", 10, "Darmor_noAED1_10topBPweight01.txt")
writeTopGO("all.newPan.prot.fa_swissprot_first_hit.mergedFinal.OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Napus_Pangenome_With_Synths_Variable_GeneNames.txt", "MF", 10, "Darmor_noAED1_10topMFweight01.txt")
writeTopGO("all.newPan.prot.fa_swissprot_first_hit.mergedFinal.OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Napus_Pangenome_With_Synths_Variable_GeneNames.txt", "CC", 10, "Darmor_noAED1_10topCCweight01.txt")

writeTopGO("Oleracea_EVM_swissprot_hfirst_hit_OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Oleracea_Pangenome_Variable_GeneNames.txt", "BP", 10, "Oleracea_Pangenome_10topBPweight01.txt")
writeTopGO("Oleracea_EVM_swissprot_hfirst_hit_OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Oleracea_Pangenome_Variable_GeneNames.txt", "MF", 10, "Oleracea_Pangenome_10topMFweight01.txt")
writeTopGO("Oleracea_EVM_swissprot_hfirst_hit_OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Oleracea_Pangenome_Variable_GeneNames.txt", "CC", 10, "Oleracea_Pangenome_10topCCweight01.txt")

writeTopGO("Rapa_EVM_swissprot_frist_hit_OnlyPlantHits.txt", "Rapa_Pangenome_WithFPSc_Variable_GeneNames.txt", "BP", 10, "Rapa_Pangenome_10topBPweight01.txt")
writeTopGO("Rapa_EVM_swissprot_frist_hit_OnlyPlantHits.txt", "Rapa_Pangenome_WithFPSc_Variable_GeneNames.txt", "MF", 10, "Rapa_Pangenome_10topMFweight01.txt")
writeTopGO("Rapa_EVM_swissprot_frist_hit_OnlyPlantHits.txt", "Rapa_Pangenome_WithFPSc_Variable_GeneNames.txt", "CC", 10, "Rapa_Pangenome_10topCCweight01.txt")


# and let's try classicCOunt

writeTopGO <- function(background, foreground, GOterm, top, txt){
  # give properly formatted background in format: GO:0005838	GSBRNA2T00088508001;GSBRNA2T00088313001;GSBRNA2T00035842001 
  annAT <- readMappings(background, sep="\t", IDsep=";")
  allgenes <- unique(unlist(annAT))
  # give file with your genes of interest, one gene_id per line
  mygenes <-scan(foreground ,what="")
  geneList <- factor(as.integer(allgenes %in% mygenes))
  names(geneList) <- allgenes
  GOdata <-new ("topGOdata", ontology = GOterm, allGenes = geneList, nodeSize = top, annot=annFUN.GO2genes, GO2genes=annAT)
  
  # using ClassicCount:
  test.stat <-new ("classicCount", testStatistic = GOFisherTest, name = "Fisher Test")
  resultsFisherC <-getSigGroups (GOdata, test.stat)
  
  # using weight01:
  #weight01.fisher <- runTest(GOdata, statistic = "fisher")
  
  # using ClassicCount:
  allRes <- GenTable(GOdata, classicFisher= resultsFisherC, topNodes = 30)
  
  # using weight01:
  #allRes <- GenTable(GOdata, classicFisher=weight01.fisher,topNodes=30)
  
  names(allRes)[length(allRes)] <- "p.value"
  
  write.table(allRes, file=txt, sep="\t",quote=FALSE,row.names=FALSE)
}


writeTopGO("all.newPan.prot.fa_swissprot_first_hit.mergedFinal.OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Napus_Pangenome_With_Synths_Variable_GeneNames.txt", "BP", 5, "Darmor_noAED1_5topBPclassicCount01.txt")
writeTopGO("all.newPan.prot.fa_swissprot_first_hit.mergedFinal.OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Napus_Pangenome_With_Synths_Variable_GeneNames.txt", "MF", 5, "Darmor_noAED1_5topMFclassicCount01.txt")
writeTopGO("all.newPan.prot.fa_swissprot_first_hit.mergedFinal.OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Napus_Pangenome_With_Synths_Variable_GeneNames.txt", "CC", 5, "Darmor_noAED1_5topCCclassicCount01.txt")

writeTopGO("Oleracea_EVM_swissprot_hfirst_hit_OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Oleracea_Pangenome_Variable_GeneNames.txt", "BP", 5, "Oleracea_Pangenome_5topBPclassicCount01.txt")
writeTopGO("Oleracea_EVM_swissprot_hfirst_hit_OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Oleracea_Pangenome_Variable_GeneNames.txt", "MF", 5, "Oleracea_Pangenome_5topMFclassicCount01.txt")
writeTopGO("Oleracea_EVM_swissprot_hfirst_hit_OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Oleracea_Pangenome_Variable_GeneNames.txt", "CC", 5, "Oleracea_Pangenome_5topCCclassicCount01.txt")

writeTopGO("Rapa_EVM_swissprot_frist_hit_OnlyPlantHits.txt", "Rapa_Pangenome_WithFPSc_Variable_GeneNames.txt", "BP", 5, "Rapa_Pangenome_5topBPclassicCount01.txt")
writeTopGO("Rapa_EVM_swissprot_frist_hit_OnlyPlantHits.txt", "Rapa_Pangenome_WithFPSc_Variable_GeneNames.txt", "MF", 5, "Rapa_Pangenome_5topMFclassicCount01.txt")
writeTopGO("Rapa_EVM_swissprot_frist_hit_OnlyPlantHits.txt", "Rapa_Pangenome_WithFPSc_Variable_GeneNames.txt", "CC", 5, "Rapa_Pangenome_5topCCclassicCount01.txt")


# let's try 10 too

writeTopGO("all.newPan.prot.fa_swissprot_first_hit.mergedFinal.OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Napus_Pangenome_With_Synths_Variable_GeneNames.txt", "BP", 10, "Darmor_noAED1_10topBPclassicCount01.txt")
writeTopGO("all.newPan.prot.fa_swissprot_first_hit.mergedFinal.OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Napus_Pangenome_With_Synths_Variable_GeneNames.txt", "MF", 10, "Darmor_noAED1_10topMFclassicCount01.txt")
writeTopGO("all.newPan.prot.fa_swissprot_first_hit.mergedFinal.OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Napus_Pangenome_With_Synths_Variable_GeneNames.txt", "CC", 10, "Darmor_noAED1_10topCCclassicCount01.txt")

writeTopGO("Oleracea_EVM_swissprot_hfirst_hit_OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Oleracea_Pangenome_Variable_GeneNames.txt", "BP", 10, "Oleracea_Pangenome_10topBPclassicCount01.txt")
writeTopGO("Oleracea_EVM_swissprot_hfirst_hit_OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Oleracea_Pangenome_Variable_GeneNames.txt", "MF", 10, "Oleracea_Pangenome_10topMFclassicCount01.txt")
writeTopGO("Oleracea_EVM_swissprot_hfirst_hit_OnlyPlantHits.txt_out.mergedFinal.OnlyPlantHits.txt", "Oleracea_Pangenome_Variable_GeneNames.txt", "CC", 10, "Oleracea_Pangenome_10topCCclassicCount01.txt")

writeTopGO("Rapa_EVM_swissprot_frist_hit_OnlyPlantHits.txt", "Rapa_Pangenome_WithFPSc_Variable_GeneNames.txt", "BP", 10, "Rapa_Pangenome_10topBPclassicCount01.txt")
writeTopGO("Rapa_EVM_swissprot_frist_hit_OnlyPlantHits.txt", "Rapa_Pangenome_WithFPSc_Variable_GeneNames.txt", "MF", 10, "Rapa_Pangenome_10topMFclassicCount01.txt")
writeTopGO("Rapa_EVM_swissprot_frist_hit_OnlyPlantHits.txt", "Rapa_Pangenome_WithFPSc_Variable_GeneNames.txt", "CC", 10, "Rapa_Pangenome_10topCCclassicCount01.txt")

# Now comes the plotting
library("wordcloud")

files <- list.files(path = '.', pattern='*5topBP*weight*')

files
files <- files[!grepl('*png', files)]

pal3 <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#CAB2D6","#6A3D9A","#FFFF99","#B15928")


for (i in seq_along(files)) {
  t <- read_tsv(files[i])

  t$p.value <- as.numeric(str_replace(t$p.value, "< 1e-30", '1e-30'))
  png(paste(c(files[i], 'wordcloud.png'), collapse=''))
  wordcloud(t$Term, t$p.value, min.freq=1, random.order = FALSE, rot.per=0.0, colors=pal3)
  dev.off()
}

