I aligned all reads with all pangenomes using bowtie2 --end-to-end --sensitive

Then, I had to merge all read alignments into one library per BioSample, NOT per SRA accession -we can have X libraries per individual after all.

if you have only the RUN IDs, this R-script will use the SRA API to get the BioSample and the Species:

library(rentrez)
library(assertthat)
library(readr)
search_ind <- function(term){
  # get the IDs for a run ID
  # ERR457868 searched, returns 1011219
  results <- entrez_search(db="sra", term=term)$ids
  assert_that(length(results) == 1)
  results
}

search_links <- function(term) {
  # get the Links in the SRA for a specific run's ID
  # then pull out only the biosample from there
  results <- entrez_link(dbfrom='sra', id=term, db='biosample')$links$sra_biosample
  # there should be only one biosample per ID
  assert_that(length(results) == 1)
  results
}

search_summary <- function(term) {
  # Get the summary for a BioSample, then return only the accession
  summary <- entrez_summary(db='biosample', id=term)
  accession <- summary$accession
  accession
}

get_species <- function(term) {
  summary <- entrez_summary(db='biosample', id=term)
  organism <- summary$organism
  organism
}
names <- read.table('./Names_only.txt', head=F)
#ERR457868 ERR475358 ERR475359 ERR475360 ERR475361 ERR479604 ....
ids <- sapply(names$V1, search_ind, USE.NAMES=F)
#[1] "1011219" "1533032" "984971"  "984969"  "984970"  "1533033"....
links <- sapply(ids, search_links, USE.NAMES=F)
# "3087115" "3769628" "3031276" "3031277" "3031278" "3769630"
accession_ids <- sapply(links, search_summary, USE.NAMES=F)
# [1] "SAMEA2399445" "SAMEA2445339" "SAMEA2729910" "SAMEA2445340" "SAMEA2445341" "SAMEA2467095"
names$V2 <- accession_ids
organisms <- sapply(links, get_species, USE.NAMES=F)
names$V3 <- organisms
write_csv(names, 'IDs_to_Samples.csv')

Here's the code to make the samtools merge commands:

The table Reads_to_Samples.csv looks like this:

ERR479609       SAMEA2467100
SRR5007238      SAMN05756140
SRR5007229      SAMN05756140
SRR3032032      SAMN04349657
ERR479616       SAMEA2467107
ERR479611       SAMEA2467102
ERR479612       SAMEA2467103
ERR479634       SAMEA2467125
ERR479613       SAMEA2467104
ERR479635       SAMEA2467126
ERR479629       SAMEA2467120
SRR3032030      SAMN04349655
ERR479608       SAMEA2467099
ERR479606       SAMEA2467097
ERR457868       SAMEA2399445
ERR479636       SAMEA2467127
SRR3420372      SAMN04857229
ERR457816       SAMEA2399444
ERR457887       SAMEA2399444
SRR5007251      SAMN05756141
SRR5007230      SAMN05756140
ERR479628       SAMEA2467119
...

from collections import defaultdict
from glob import glob
cultivar_by_read = defaultdict(list) # key: cultivar / biosample, value:
for line in open('Reads_to_Samples.csv'):
    ll = line.split()
    library, sample=  ll
    cultivar_by_read[sample].append(library)

refs = ['Darmor_NRGene_pseudomolecules_v1_chloroplast_mitochondrium_both_rounds', 'Brassica_oleracea.v2.1.dna.toplevel_mito_chloro_both', 'BrapaV2.1PacBio.Chr_mito_chloro_with_first_round_and_second_round']

moveout = open('Move_us_away.txt','w')

for s in cultivar_by_read:
    libs = cultivar_by_read[s]
    if len(libs) == 1: 
        # leave unmerged
        continue

    skip = False
    all_files = [ list(), list(), list() ]
    for l in libs:
        for counter, r in enumerate(refs):
            this_file = glob('*%s*%s*sorted.bam'%(l, r))
            if len(this_file) != 1:
                skip = True
                print('*%s*%s*sorted.bam'%(l,r))

            this_file = this_file[0]
            all_files[counter].append(this_file)

    # Napus_SRR5007238_1.fastq_Darmor_NRGene_pseudomolecules_v1_chloroplast_mitochondrium_both_rounds_sorted.bam to Napus
    thistype = all_files[0][0].split('_')[0]
    outnames = []
    for r in refs:
        thisname = '%s_MERGED_%s_%s_sorted.bam'%(thistype, s, r)
        outnames.append(thisname)

    for counter, f in enumerate(all_files):
        out = outnames[counter]
        f = ['/scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/PAV_Alignments/' + x for x in f]
        cmd = '"samtools merge -f /scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/PAV_Alignments/%s %s"'%(out, ' '.join(f))
        for x in f:
            moveout.write(x + '\n')
            moveout.write(x + '.bai\n')
        print(cmd)

And deleted the files in Move_us_away.txt.

Then with he merged libraries I ran mosdepth.

I ran mosdepth using mostly standard parameters so I get a per-base bed file.

I made a gff3 file with only the exons (grep -P '\texon\t' newPan.gff3) and ran bedtools intersect with my mosdepth output and the exon-only gff3:

   bedtools intersect -a /scratch/pawsey0149/pbayer/DOwnloadPanPanData/fastq/10081_Harsh1_S1_L005_R1_001.fastq_newPan_sorted_MERGED.bam_mosdepth.per-base.bed.gz -b /scratch/pawsey0149/pbayer/DOwnloadPanPanData/fastq/newPan_exon_only.gff3 -wao | grep -v '\-1' | awk '{if (\$4 >=2) print}' > /scratch/pawsey0149/pbayer/DOwnloadPanPanData/fastq/10081_Harsh1_S1_L005_R1_001.fastq_newPan_sorted_MERGED.bam_mosdepth.per-base.bed.gz_vs_exons.bed

Then I wrote a script which takes a output bed file, and writes a single row into the whole table of PAV. I first used bedops' gff2bed to get a bed file of my genes.

   gff2bed < newPan.gff3 > newPan.bed

so I don't get confused with my start and end positions, as bed counts from 0 and gff3 counts from 1.

Script to call PAV - this will write two lines per input mosdepth output file:

import sys
from collections import defaultdict
gene_dict = defaultdict(int) # key: gene name, value: total length of gene based on exons

with open('/scratch/pawsey0149/pbayer/DOwnloadPanPanData/fastq/newPan.bed') as fh:
    '''chrA01  2829    3118    evm.model.chrA01.1.exon6        .       -       EVM     exon    .       ID=evm.model.chrA01.1.exon6;Parent=evm.model.chrA01.1
    chrA01  3500    3700    cds.evm.model.chrA01.1  .       -       EVM     CDS     1       ID=cds.evm.model.chrA01.1;Parent=evm.model.chrA01.1
    chrA01  3500    3700    evm.model.chrA01.1.exon5        .       -       EVM     exon    .       ID=evm.model.chrA01.1.exon5;Parent=evm.model.chrA01.1
    chrA01  3776    3909    cds.evm.model.chrA01.1  .       -       EVM     CDS     2       ID=cds.evm.model.chrA01.1;Parent=evm.model.chrA01.1
    chrA01  3776    3909    evm.model.chrA01.1.exon4        .       -       EVM     exon    .       ID=evm.model.chrA01.1.exon4;Parent=evm.model.chrA01.1
    chrA01  4021    4385    cds.evm.model.chrA01.1  .       -       EVM     CDS     0       ID=cds.evm.model.chrA01.1;Parent=evm.model.chrA01.1'''

    for line in fh:
        ll = line.split()
        if len(ll) < 3: continue
        if ll[7] != 'exon': continue
        start, end = map(int, ll[1:3])
        length = end - start
        name = ll[-1].split(';')[-1]
        assert 'Parent' in name
        name = name.replace('Parent=','')
        gene_dict[name] += length

f = sys.argv[1]
overlap_dict = defaultdict(int)
with open(f) as fh:
    '''chrA01  2526    2527    927     chrA01  EVM     exon    2527    2599    .       -       .       ID=evm.model.chrA01.1.exon7;Parent=evm.model.chrA01.1   1
    chrA01  2530    2532    907     chrA01  EVM     exon    2527    2599    .       -       .       ID=evm.model.chrA01.1.exon7;Parent=evm.model.chrA01.1   2
    '''
    for line in fh:
        ll = line.split()
        cov = int(ll[3])
        if cov < 2: continue
        start_overlap, end_overlap = map(int, ll[1:3])
        length_overlap = end_overlap - start_overlap
        name = ll[-2].split(';')[-1]
        assert 'Parent' in name
        name = name.replace('Parent=','')
        overlap_dict[name] += length_overlap

genes = sorted(gene_dict)
header = ['Individual'] + genes
print('\t'.join(header))
ll = [f]
for g in genes:
    overlap_percentage = overlap_dict[g]/float(gene_dict[g])*100
    if overlap_percentage <= 5:
        # LOST
        ll.append('0')
    else:
        ll.append('1')

print('\t'.join(ll))

Usage:


   python makePAV.py Brassica_rapa_SAMN04498195_MERGED.bam_mosdepth.per-base.bed.gz_vs_exons.bed > Brassica_rapa_SAMN04498195_MERGED.bam_mosdepth.per-base.bed.gz_vs_exons_PAV.tsv

Then, run that script for all bedtools output, then

   head -1 Brassica_napus_SAMN04218413_MERGED.bam_mosdepth.per-base.bed.gz_vs_exons_PAV.tsv > Header
   cat *PAV.tsv | grep -v Ind >> Header 
   mv Header Final_PAV_Table.tsv

