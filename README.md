
# END-VoC WP3 - D3.2 SARS-CoV-2 NGS Bioinformatics pipeline 




**######################### Gather the raw read dataset ########################**

 849 SARS-CoV-2 sequencing projects from 2019/11/01 to 2020/05/12 were downloaded from the NCBI SRA raw read database
SRA can be searched in various ways, including (i) on a web browser at https://www.ncbi.nlm.nih.gov/sra using a query such as "("2019/11/01"[Publication Date] : "2020/05/12"[Publication Date]) AND txid2697049[Organism:noexp] NOT 00000000000[Mbases]" or using the e-utility dedicated command line tool with commands such as "elink -target sra -db taxonomy -id 2697049 | efetch -mode xml > 00.sra.xml "

In the same manner, downloading the actual fastq files of the dataset can be done either manualy or by using appropriate command line tools (Aspera, fasterq-dump ....
Detailed scripts are not provided here because the reader should have its own dataset to analyse.

the list of datasets studied as well as information about the samples is provided in Table_1.

**######################### Produce the consensus sequences of each sample ########################**

The dataset comprises sequencing data produced with an Illumina amplicon strategy and with a Nanopore amplicon strategy 

**Analysis of the Illumina based sequencing projects**
Details about how to specify your input files paths in ILLUMINA_AMPLICON.samplesheet.csv are provided at https://nf-co.re/viralrecon/2.5/usage

    bash nextflow run nf-core/viralrecon -r 2.4.1  \
        --input ILLUMINA_AMPLICON.samplesheet.csv \
        --outdir ./output_ILLUMINA_AMPLICON \
        --platform illumina \
        --protocol amplicon \
        --genome 'MN908947.3' \
        --primer_set artic \
        --primer_set_version 1 \
        --skip_kraken2 \
        --skip_assembly \
        -profile singularity \
        --consensus_caller ivar \
        --max_memory '30.GB' \
        --max_cpus 30


**Analysis of the Nanopore based sequencing projects**
Details about how to organize your input file folder architecture and specify your input files paths in NANOPORE_AMPLICON.samplesheet.cs are provided at https://nf-co.re/viralrecon/2.5/usage

bash nextflow run nf-core/viralrecon -r 2.4.1  \
--input NANOPORE_AMPLICON.samplesheet.csv \
--outdir ./output_NANOPORE_AMPLICON/ \
--platform nanopore \
--genome 'MN908947.3' \
--primer_set_version 1 \
--fastq_dir ./OXFORD_NANOPORE_AMPLICON/ \
--artic_minion_caller medaka \
--artic_minion_medaka_model r941_min_high_g360 \
--artic_minion_aligner minimap2 \
-profile singularity \
--skip_nanoplot \
--max_memory '30.GB' \
--max_cpus 30

**######################### Produce an alignment for subsequent analysis ########################**

Locate your consensus sequences.
In the above case, consensus sequences of the Illumina data are in 

> ./output_ILLUMINA_AMPLICON/variants/ivar/consensus/ivar/*.consensus.fa

and those of Nanopore data are in 

> ./output_NANOPORE_AMPLICON/*.consensus.fasta

all sequences are aligned against the Wuhan-Hu-1 reference sequence (EPI_ISL_402125)

    cat EPI_ISL_402125.fasta > total.aln
    cat ./output_ILLUMINA_AMPLICON/variants/ivar/consensus/ivar/*.consensus.fa >> total.aln
    cat ./output_NANOPORE_AMPLICON/*.consensus.fasta >> total.aln

The actual alignment is performed using Mafft, embedded into the augur bioinformatic toolkit

    augur align --sequences total.aln --reference-name 'EPI_ISL_402125' --fill-gaps --output total_aligned.aln --nthreads 8



Page generated using https://stackedit.io/app


