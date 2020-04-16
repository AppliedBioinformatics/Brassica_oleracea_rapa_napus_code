First, I used BWA to estimate the insert sizes of all libraries.

I took first 100000 reads:

   for l in *paired.fastq; do head -100000 $l > ${l}_head; done

   for l in `pwd`/*_1.fastq; do echo \"/scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/bwa-0.7.16a/bwa aln -n 15 -k 3 -t 4 /scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/Napus/Darmor_NRGene_pseudomolecules_v1_chloroplast_mitochondrium.fasta $l \> ${l}.sai\" ;done

   for l in `pwd`/*_1_paired.fastq_head.sai; do second=${l/_1_paired/_2_paired}; echo \" /scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/bwa-0.7.16a/bwa sampe /scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/Oleracea/Brassica_oleracea.v2.1.dna.toplevel_mito_chloro.fa $l $second ${l%%.sai} ${second%%.sai} \> /dev/null 2\> ${l}_stats\"; done

Files ending in _stats have the insert sizes.

I removed all libraries with large (>1000bp) insert sizes, these are mistaken mate-paired libraries.

Alignment command for bowtie2 before masurca:

   for l in `pwd`/*_1_paired.fastq; do echo \"bowtie2 -I 0 -X 1000 -x `pwd`/../Darmor_NRGene_pseudomolecules_v1_chloroplast_mitochondrium -1 $l -2 ${l/1_paired.fastq/2_paired.fastq} --end-to-end --sensitive --threads 4  2\> ${l}_bowtie2stats \| samtools view -Sb - \| samtools sort - ${l}_sorted\"; done

For all of the first samples' list above, extracted forward and reverse unmapped reads:

   "samtools view -h -b -f 68 /scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/Rapa/fastq/SRR1592657_1_paired.fastq_sorted.bam > /scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/Rapa/fastq/SRR1592657_1_paired.fastq_sorted_unmapped_first.bam"
   "samtools view -h -b -f 132 /scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/Rapa/fastq/SRR1592657_1_paired.fastq_sorted.bam > /scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/Rapa/fastq/SRR1592657_1_paired.fastq_sorted_unmapped_second.bam"

etc.

Then merged those:

   "samtools merge /scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/Rapa/fastq/SRR3926923_1_paired.fastq_sorted_unmapped_merged.bam /scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/Rapa/fastq/SRR3926923_1_paired.fastq_sorted_unmapped_first.bam /scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/Rapa/fastq/SRR3926923_1_paired.fastq_sorted_unmapped_second.bam"

and sorted by name:

   samtools sort -n /scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/Napus/fastq/ERR479650_1_paired.fastq_sorted_unmapped_merged.bam /scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/Napus/fastq/ERR479650_1_paired.fastq_sorted_unmapped_merged_sortedName

and converted to fastq:

   "bamtools convert -in /scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/Rapa/fastq/SRR3926923_1_paired.fastq_sorted_unmapped_merged_sortedName.bam -out /scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/Rapa/fastq/SRR3926923_1_paired.fastq_sorted_unmapped_merged_sortedName.fastq -format fastq"

and split and synchronised the fastq using this script:

<pre>
import sys
from Bio import SeqIO

out_paired_1 = sys.argv[1] + "_R1.fastq"
out_paired_2 = sys.argv[1] + "_R2.fastq"
out_unpaired = sys.argv[1] + "_unpaired.fastq"

out_paired_1 = open(out_paired_1, "w")
out_paired_2 = open(out_paired_2, "w")
out_unpaired = open(out_unpaired, "w")

read_dict = {} # { key: read name, value: {"R1":SeqRecord, "R2":SeqRecord} }

for l in SeqIO.parse(sys.argv[1], "fastq"):
    short = l.id.split('/')[0]
    if short not in read_dict:
        read_dict[short] = {"R1":False, "R2":False}

    if "/1" in l.id:
        read_dict[short]["R1"] = l
    else:
        read_dict[short]["R2"] = l

    # do we have a proper pair?
    if read_dict[short]["R2"] and read_dict[short]["R1"]:
        out_paired_1.write(read_dict[short]["R1"].format("fastq"))
        out_paired_2.write(read_dict[short]["R2"].format("fastq"))
        # we will never see this again, save some memory
        del read_dict[short]

        # now write out the remaining singles, it's virtually guaranteed that
        # these won't have another pair
        for r in read_dict:
            r1 = read_dict[r]["R1"]
            r2 = read_dict[r]["R2"]
            if r1:
                out_unpaired.write(r1.format("fastq"))
            if r2:
                out_unpaired.write(r2.format("fastq"))
        del read_dict
        read_dict = {}

</pre>

Usage:

   "python /scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/splitFiles.py /scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/Rapa/fastq/SRR3926922_1_paired.fastq_sorted_unmapped_merged_sortedName.fastq"

and used those files in masurca.

Here's my script to take those files and write three masurca control files:


<pre>
insert_sizes = {}

for line in open('Napus_Rapa_Oleracea_Insertsizes_all.csv'):
    ll = line.split()
    name, insert, std = ll[0], ll[1], ll[3]
    insert_sizes[name] = (insert, std)

template = open('test.conf')
template_head='''# example configuration file 

# DATA is specified as type {PE,JUMP,OTHER,PACBIO} and 5 fields:
# 1)two_letter_prefix 2)mean 3)stdev 4)fastq(.gz)_fwd_reads
# 5)fastq(.gz)_rev_reads. The PE reads are always assumed to be
# innies, i.e. --->.<---, and JUMP are assumed to be outties
# <---.--->. If there are any jump libraries that are innies, such as
# longjump, specify them as JUMP and specify NEGATIVE mean. Reverse reads
# are optional for PE libraries and mandatory for JUMP libraries. Any
# OTHER sequence data (454, Sanger, Ion torrent, etc) must be first
# converted into Celera Assembler compatible .frg files (see
# http://wgs-assembler.sourceforge.com)
DATA
'''
template_bottom = '''
END

PARAMETERS
#this is k-mer size for deBruijn graph values between 25 and 127 are supported, auto will compute the optimal size based on the read data and GC content
GRAPH_KMER_SIZE = auto
#set this to 1 for all Illumina-only assemblies
#set this to 1 if you have less than 20x long reads (454, Sanger, Pacbio) and less than 50x CLONE coverage by Illumina, Sanger or 454 mate pairs
#otherwise keep at 0
USE_LINKING_MATES = 1
#this parameter is useful if you have too many Illumina jumping library mates. Typically set it to 60 for bacteria and 300 for the other organisms 
LIMIT_JUMP_COVERAGE = 300
#these are the additional parameters to Celera Assembler.  do not worry about performance, number or processors or batch sizes -- these are computed automatically. 
#set cgwErrorRate=0.25 for bacteria and 0.1<=cgwErrorRate<=0.15 for other organisms.
CA_PARAMETERS =  cgwErrorRate=0.15
#minimum count k-mers used in error correction 1 means all k-mers are used.  one can increase to 2 if Illumina coverage >100
KMER_COUNT_THRESHOLD = 1
#whether to attempt to close gaps in scaffolds with Illumina data
CLOSE_GAPS=1
#auto-detected number of cpus to use
NUM_THREADS = 16
#this is mandatory jellyfish hash size -- a safe value is estimated_genome_size*estimated_coverage
JF_SIZE = 358916659
#set this to 1 to use SOAPdenovo contigging/scaffolding module.  Assembly will be worse but will run faster. Useful for very large (>5Gbp) genomes
SOAP_ASSEMBLY=0
END
'''


import glob

from itertools import product
from string import ascii_lowercase
keywords = [''.join(i) for i in product(ascii_lowercase, repeat = 2)]

ol_out = open('masurca_olera.conf','w')
ra_out = open('masurca_rapa.conf','w')
na_out = open('masurca_napus.conf','w')

ol_out.write(template_head)
ra_out.write(template_head)
na_out.write(template_head)


counter = 0
for r1 in glob.glob('/scratch/pawsey0149/pbayer/Brassica_rapa_oleracea_data/*/fastq/*sortedName.fastq_R1.fastq'):
    r2 = r1.replace('R1.fastq','R2.fastq')
    unp = r1.replace('R1.fastq', 'unpaired.fastq')
    sra = r1.split('/')[-1].split('_')[0]
    isize, std = insert_sizes[sra]
    type = r1.split('/')[-3]
    keyword = keywords[counter]
    counter += 1
    # PE= pe 180 20  /FULL_PATH/frag_1.fastq  /FULL_PATH/frag_2.fastq
    line1 = 'PE= %s %s %s %s %s\n'%(keyword, isize, std, r1, r2)
    keyword = keywords[counter]
    line2 = 'PE= %s %s %s %s\n'%(keyword, isize, std, unp) 
    counter += 1
    if type == 'Oleracea':
        ol_out.write(line1)
        ol_out.write(line2)
    elif type == 'Rapa':
        ra_out.write(line1)
        ra_out.write(line2)
    elif type =='Napus':
        na_out.write(line1)
        na_out.write(line2)
    else:
        print type
ol_out.write(template_bottom)
ra_out.write(template_bottom)
na_out.write(template_bottom)
</pre>

This will write the control files for the three pangenomes.

Be aware that masurca always makes a file called 'assemble.sh' in the current working directory, so running three jobs at a time will result in problems, which is why I made subdirectories for all of these files.
