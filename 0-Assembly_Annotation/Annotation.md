
I first used RepeatModeler to call repeats de novo, using only contigs bigger than 10kb.

    BuildDatabase -name Brassica_oleracea.v2.1.dna.toplevel_mito_chloro_both_rounds_10kb_filtered Brassica_oleracea.v2.1.dna.toplevel_mito_chloro_both_rounds_10kb_filtered.fasta
    BuildDatabase -name Darmor_NRGene_pseudomolecules_v1_chloroplast_mitochondrium_both_rounds_10kb_filtered Darmor_NRGene_pseudomolecules_v1_chloroplast_mitochondrium_both_rounds_10kb_filtered.fasta

    RepeatModeler -pa 12 -database Brassica_oleracea.v2.1.dna.toplevel_mito_chloro_both_rounds_10kb_filtered > OLERACEA_OUT 2> OLERACEA_ERR 
    RepeatModeler -pa 12 -database Darmor_NRGene_pseudomolecules_v1_chloroplast_mitochondrium_both_rounds_10kb_filtered > NAPUS_OUT 2> NAPUS_ERR 

Then I used hisat2 to align all RNASeq reads (see below commands) and merged the bam files using picard MergeSamFiles (important: merge sequence dictionaries!)

    java -jar picard.jar MergeSamFiles OUTPUT=merged.bam MERGE_SEQUENCE_DICTIONARIES=true INPUT=ERR431516_1.fastq_sorted.bam  .....

Then I used RepeatMasker with the predicted repeats as -lib argument, IMPORTANT: RepeatMasker needs to be run with xsmall to make softmasked references!

    RepeatMasker -lib repeatmodeler_output -pa 12 -xsmall reference.fasta

and ran BRAKER the following way:

    braker.pl --genome=genome.fa --species=speciesname --bam=accepted hits.bam --softmasked

Then I used EvidenceModeler to merge the final output of genemark and augustus. I removed all gene models without RNASeq support.

<pre>
perl EVidenceModeler-1.1.1/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl augustus.gtf > augustus_evm.gtf

genemark_gtf2gff3 genemark.gtf > genemark.gff3
perl EVidenceModeler-1.1.1/EvmUtils/gff3_gene_prediction_file_validator.pl genemark.gff3

cat augustus_evm.gtf genemark.gff3 > augustus_evm_and_genemark.gff3

grep -v '^#' stringtie_all_merged.gtf | ./EVidenceModeler-1.1.1/EvmUtils/misc/cufflinks_gtf_to_alignment_gff3.pl - >  oleracea_stringtie.gff3

# | in sequence names cause errors
sed -i.bak 's/|/_/g' augustus_evm_and_genemark.gff3
sed -i.bak 's/|/_/g' Brassica_oleracea.v2.1.dna.toplevel_mito_chloro_both_rounds.fasta.masked
sed -i.bak 's/|/_/g' oleracea_stringtie.gff3 &
sed -i.bak 's/|/_/g' Brassica_oleracea.v2.1.dna.toplevel_mito_chloro_both_rounds.fasta.gff3

$EVM_HOME/EvmUtils/partition_EVM_inputs.pl --repeats Brassica_oleracea.v2.1.dna.toplevel_mito_chloro_both_rounds.fasta.gff3 --gene_predictions augustus_evm_and_genemark.gff3 --transcript_alignments oleracea_stringtie.gff3 --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out --genome Brassica_oleracea.v2.1.dna.toplevel_mito_chloro_both_rounds.fasta.masked

$EVM_HOME/EvmUtils/write_EVM_commands.pl --repeats Brassica_oleracea.v2.1.dna.toplevel_mito_chloro_both_rounds.fasta.gff3 --gene_predictions augustus_evm_and_genemark.gff3 --transcript_alignments oleracea_stringtie.gff3 --genome Brassica_oleracea.v2.1.dna.toplevel_mito_chloro_both_rounds.fasta.masked --weights `pwd`/weights.txt --output_file_name evm.out --partitions partitions_list.out > commands.list

$EVM_HOME/EvmUtils/execute_EVM_commands.pl commands.list | tee run.log

$EVM_HOME/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out

python removeNonRNASeq.py
</pre>

The script is here:
<pre>
import subprocess

run_command = 'find . -maxdepth 2 -name evm.out'
process = subprocess.Popen(run_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
for c in iter(lambda: process.stdout.readline(), ''):
    out = open(c.rstrip().replace('.out','_nonRNARemoved.out'), 'w')
    this_gene = []
    for line in open(c.rstrip()):
        if line.startswith('!!'):
            out.write(line)
            continue
        if line.startswith('#'):
            if this_gene:
                # check whether this_gene has RNASeq in it
                keep = False
                for l in this_gene:
                    if 'Cufflinks' in l:
                        keep = True
                        break

                if not keep: continue

                for l in this_gene:
                    out.write(l)

            this_gene = [line]
        else:
            this_gene.append(line)
    # the last gene still has to be printed
    if this_gene:
        keep = False
        for l in this_gene:
            if 'Cufflinks' in l:
                keep = True
                break

        if not keep: continue

        for l in this_gene:
            out.write(l)
</pre>

This will iterate over all evm.out files and delete genes without RNASeq support and write new evm.out files ending in 'nonRNARemoved'. Check whether they have content where needed and nothing got borked, then rename them to evm.out.

<pre>
for l in $(find . -maxdepth 2 -name evm_nonRNARemoved.out); do mv $l ${l%%_nonRNARemoved.out}.out; done

$EVM_HOME/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome Brassica_oleracea.v2.1.dna.toplevel_mito_chloro_both_rounds.fasta.masked

find . -regex ".*evm.out.gff3" -exec cat {} \; > EVM.all_without_RNASeq.gff3
</pre>

That single gff3 file wll have the final gene models. Use Stringtie's cuffread to make protein and cds fastas from that.

