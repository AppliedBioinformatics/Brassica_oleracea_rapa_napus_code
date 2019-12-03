The files in this folder collect how I ran the GO-enrichment analysis and visualisation.


I split up all input protein files for the three annotations using seqkit split -s 2000, and download swissprot

    wget --continue ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_swissprot.fasta.gz

swissprot from the same location 
UniProtKB/Swiss-Prot Release 2019_04 of 08-May-2019

Then run blastp:

    diamond blastp --query /scratch/pawsey0149/pbayer/BLAST_NR_Brassica_PanPan/Oleracea_EVM.without_RNASeq_FINAL.No_contaminants_no_small.gff3_prot.fa.split/Oleracea_EVM.without_RNASeq_FINAL.No_contaminants_no_small.gff3_prot.part_064.fa --db /scratch/pawsey0149/pbayer/BLAST_NR_Brassica_PanPan/uniprot_both.fasta --outfmt 6 --evalue 1e-20 --out /scratch/pawsey0149/pbayer/BLAST_NR_Brassica_PanPan/Oleracea_EVM.without_RNASeq_FINAL.No_contaminants_no_small.gff3_prot.fa.split/Oleracea_EVM.without_RNASeq_FINAL.No_contaminants_no_small.gff3_prot.part_064.fa_vs_uniprot.csv

Then downloaded GOA from ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz last modified 9/05/2019 (May)

Kept only the top blast hit for all three species:

    awk '{ if(!x[$1]++) {print $0; bitscore=($14-1)} else { if($14>bitscore) print $0} }' Napus_all.csv > Napus_best_hit_only.csv

And ran this script to get the GO term by gene:
<pre>
from glob import glob
filename_dict = {} # key: filename, value: {uniprot_gene: [napus_gene1, napus_gene2...]}
filename_dict_genes = {} # key: filename, value: {go_term1:[napus_gene1, napus_gene2...]}

for filename in glob('*only.csv'):
    this_dict = {}
    for line in open(filename):
        ll = line.split()
        napus_gene, uniprot_hit = ll[0], ll[1]
        # tr|M4DY10|M4DY10_BRARP to M4DY10
        uniprot_hit = uniprot_hit.split('|')[1]
        if uniprot_hit in this_dict:
            this_dict[uniprot_hit].append(napus_gene)
        else:
            this_dict[uniprot_hit] = [napus_gene]
    filename_dict[filename] = this_dict
    filename_dict_genes[filename] = {}


# final aim:
# GO:0003700 gene_A;gene_B;gene_C
for line in open('goa_uniprot_all.gaf'):
    if line.startswith('!'):
        continue
    ll = line.split('\t')
    uniprot = ll[1]
    go = ll[4]

    for f in filename_dict:
        if uniprot in filename_dict[f]:
            # f = Napus_best_hit_only.csv
            # uniprot = A041231
            # filename_dict[f][uniprot] = ['evm.model.chrA05.1728', 'evm.model.chrA05.1766', 'evm.model.unplaced_contigs.475']
            hits = filename_dict[f][uniprot]
            if go in filename_dict_genes[f]:
                filename_dict_genes[f][go] += hits
            else:
                filename_dict_genes[f][go] = hits

for f in filename_dict_genes:
    out = open(f + '_go_terms.csv', 'w')
    for go in filename_dict_genes[f]:
        out.write('%s %s\n'%(go, ';'.join(filename_dict_genes[f][go])))
</pre>


I gzipped those output files, and then used topGO for enrichment:

<pre>
library('topGO')
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
  ontology = GOterm,pwd

  allGenes = geneList,
  #nodeSize = top, # REVIGO DOES NOT USE THIS?
  annot = annFUN.GO2genes,
  GO2genes = annAT
)GOdata <- new (
  "topGOdata",
  ontology = GOterm,

  allGenes = geneList,
  #nodeSize = top, # REVIGO DOES NOT USE THIS?
  annot = annFUN.GO2genes,
  GO2genes = annAT
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
</pre>


And uploaded that napus_aRevigorer.txt to REVIGO, it looks like this:

<pre>
GO:0006491 0.00013
GO:0006517 0.00013
GO:0009682 0.00016
GO:0071332 0.00023
GO:0002239 0.00027
GO:0051013 0.00028
GO:0009409 0.00029
GO:0006979 0.00030
....
</pre>

Then, I installed CirGO: https://github.com/IrinaVKuznetsova/CirGO

and downloaded the results from REVIGO under the TreeMap tab, saved as csv, and ran CirGO:

     python2 CirGO.py -inputFile ../../napus_REVIGO.csv -outputFile Visual_BP.svg -fontSize 6.5 -numCat 40 -legend "Name & Proportion of Biological process (inner ring)"

I changed some code within CirGOVisualGO.py:

around line 390, I added some code so that longer GO-terms have line breaks inserted:

<pre>
390         texts[nr].set_rotation_mode("anchor")
391         texts[nr].set_x(x)
392         texts[nr].set_y(y)
393         temp = []
394         import textwrap
395         for o in outer_lab:
396             newo = []
397             for i in o:
398                 newo.append('\n'.join(textwrap.wrap(i, 40)))
399             temp.append(newo)
400         outer_lab = temp
401
402         texts[nr].set_text(np.array(outer_lab).flatten()[nr])
403                                                                   
</pre>

and I commented out the legend plotting part:
<pre>
469     #plt.legend(p_inner[ ::-1],
470     #           ["%s, %1.1f%%" % (l, s) for l, s in zip(np.array(legend_labeling), pie_chart_percent )[ ::-1]],
471     #           ncol = 1,
472     #           loc = "upper left",
473     #           bbox_to_anchor = (1.25, 1.25),
474     #           frameon = False,
475     #           prop = {"size": font_size },
476     #           title = legend_name,                                                                                                                                                                           477     #           bbox_transform = plt.gcf().transFigure)
478
</pre>

and added the legend text into the figure legend in the Word document.

to get a png:

     rsvg-convert -d 300 -p 300 -o output.png Visual_BP.svg


done! 
