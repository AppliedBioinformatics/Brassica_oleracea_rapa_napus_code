library("topGO")


background <- gzfile('Napus_vs_Swissprot_best_score_only.csv_go_terms.csv.gz', 'rt')
foreground <-
  "NewPAV_Table_with_Oleracea_Rapa_Only_genes_NapusIndsOnly_VariableGenes.txt"
GOterm <- 'BP'
top <- 5
txt <- 'Napus_5topBPweight01.txt'

annAT <- readMappings(background, sep = " ", IDsep = ";")
allgenes <- unique(unlist(annAT))

# give file with your genes of interest, one gen_eid per line
mygenes <- scan(foreground , what = "")
geneList <- factor(as.integer(allgenes %in% mygenes))
names(geneList) <- allgenes

GOdata <- new (
  "topGOdata",
  ontology = GOterm,
  
  allGenes = geneList,
  #nodeSize = top, # REVIGO DOES NOT USE THIS?
  annot = annFUN.GO2genes,
  GO2genes = annAT
)

test.stat <-
  methods::new(
    'weightCount',
    testStatistic = topGO::GOFisherTest,
    name = 'Fisher test',
    sigratio = 'ratio'
  )
resultWeight <- topGO::getSigGroups(GOdata, test.stat)
topGO::geneData(resultWeight)
allRes <- topGO::GenTable(
  GOdata,
  weightFisher = resultWeight,
  orderBy = 'weightFisher',
  ranksOf = 'weightFisher',
  topNodes = length(topGO::score(resultWeight))
)

utils::write.table(
  allRes,
  file = 'napus_weightFisher.csv',
  quote = F,
  row.names = F,
  col.names = T
)
p <- 0.01

allResInf1 <- allRes[allRes$weightFisher < p,]

if (nrow(allResInf1) == 0) {
  stop(
    "No GO term has been found above the p value. Stopping now.
    Please try to increase the required p-value (option p)."
  )
}
aRevigorer = "napus_aRevigorer.txt"
utils::write.table(
  allResInf1[, c("GO.ID", "weightFisher")],
  file = aRevigorer,
  quote = F,
  row.names = F,
  col.names = F
)


# Oleracea


background <- gzfile('Oleracea_vs_Swissprot_best_score_only.csv_go_terms.csv.gz', 'rt')
foreground <-
  "NewPAV_Table_with_Oleracea_Rapa_Only_genes_OleraceaIndsOnly_OleraceaGenesOnly_VariableGenes.txt"

annAT <- readMappings(background, sep = " ", IDsep = ";")
allgenes <- unique(unlist(annAT))

# give file with your genes of interest, one gen_eid per line
mygenes <- scan(foreground , what = "")
geneList <- factor(as.integer(allgenes %in% mygenes))
names(geneList) <- allgenes

GOdata <- new (
  "topGOdata",
  ontology = GOterm,
  
  allGenes = geneList,
  #nodeSize = top, # REVIGO DOES NOT USE THIS?
  annot = annFUN.GO2genes,
  GO2genes = annAT
)

test.stat <-
  methods::new(
    'weightCount',
    testStatistic = topGO::GOFisherTest,
    name = 'Fisher test',
    sigratio = 'ratio'
  )
resultWeight <- topGO::getSigGroups(GOdata, test.stat)
topGO::geneData(resultWeight)
allRes <- topGO::GenTable(
  GOdata,
  weightFisher = resultWeight,
  orderBy = 'weightFisher',
  ranksOf = 'weightFisher',
  topNodes = length(topGO::score(resultWeight))
)

utils::write.table(
  allRes,
  file = 'oleracea_weightFisher.csv',
  quote = F,
  row.names = F,
  col.names = T
)
p <- 0.01

allResInf1 <- allRes[allRes$weightFisher < p,]

if (nrow(allResInf1) == 0) {
  stop(
    "No GO term has been found above the p value. Stopping now.
    Please try to increase the required p-value (option p)."
  )
}
aRevigorer = "oleracea_aRevigorer.txt"
utils::write.table(
  allResInf1[, c("GO.ID", "weightFisher")],
  file = aRevigorer,
  quote = F,
  row.names = F,
  col.names = F
)
