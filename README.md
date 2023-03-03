



# END-VoC WP3 - D3.2
# SARS-CoV-2 NGS Bioinformatics pipeline v1


# **## Gather the raw read dataset ##**

 849 SARS-CoV-2 sequencing projects generated from the 2019/11/01 to the 2020/05/12 were downloaded from the NCBI SRA raw read database.
SRA can be searched in various ways, including (i) on a web browser at https://www.ncbi.nlm.nih.gov/sra using a query such as `("2019/11/01"[Publication Date] : 2020/05/12"[Publication Date]) AND txid2697049[Organism:noexp] NOT 00000000000[Mbases]` or using the e-utility dedicated command line tool with commands such as `elink -target sra -db taxonomy -id 2697049 | efetch -mode xml > 00.sra.xml`

In the same manner, downloading the actual fastq files of the dataset can be done either manualy or by using appropriate command line tools (Aspera, fasterq-dump ...)
Detailed scripts are not provided here because the reader should have his own dataset to analyse.

The list of genome accessions studied as well as their metadata information are provided in Table_1.

**########## 2. Produce the consensus sequences of each sample ##########**

The dataset comprises sequencing data produced with both Illumina and Nanopore amplicon strategies. For its adaptability, we chose to analyse it using the nf-core/viralrecon workflow. Extensive documentation about it, including its installation procedure, usage and numerous possible outputs can be found at https://nf-co.re/viralrecon

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

**########## 3. Produce an alignment for subsequent analysis ##########**

Locate your consensus sequences.
In the above case, consensus sequences of the Illumina samples locate in 

> ./output_ILLUMINA_AMPLICON/variants/ivar/consensus/ivar/*.consensus.fa

and those of Nanopore samples in 

> ./output_NANOPORE_AMPLICON/*.consensus.fasta

All sequences are aligned against the Wuhan-Hu-1 reference sequence (EPI_ISL_402125, available on NCBI here: https://www.ncbi.nlm.nih.gov/search/all/?term=wuhan-hu-1), the first genome assembly published from the pandemic using the multi-sequence alignment tool MAFFT implemented via the AUGUR pipeline (https://github.com/nextstrain/augur).

    cat ./data/EPI_ISL_402125.fasta > ./data/00.total.aln
    cat ./output_ILLUMINA_AMPLICON/variants/ivar/consensus/ivar/*.consensus.fa >> ./data/00.total.aln
    cat ./output_NANOPORE_AMPLICON/*.consensus.fasta >> ./data/00.total.aln


    augur align --sequences ./data/00.total.aln --reference-name 'EPI_ISL_402125' --fill-gaps --output ./data/01.total_aligned.aln --nthreads 8

The beginning and ends of alignments can often be noisy. Hence we mask the first 55 and last 100 positions of the alignment. Sites flagged as possible sequencing errors also need to be masked. Masking consists in replacing the position with 'N'. The up to date list of sites to mask that we use is available at  https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf

You can read more about potential sequencing errors in SARS-CoV-2 genomes here: https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473 These are important to be aware of as sequencing errors will appear homoplasic so may signpost to selection or recombination artefactually. 

    perl -F'\t' -ane 'if( $F[6] eq "mask"){print $F[1] . "\n"}' ./data/problematic_sites_sarsCov2.vcf > ./data/problematic_sites_sarsCov2_inline.vcf
    python3 ./scripts/mask-alignment_v2.py --alignment ./data/01.total_aligned.aln --mask-from-beginning 55 --mask-from-end 100 --mask-sites ./data/problematic_sites_sarsCov2_inline.vcf --output ./data/02.total.mask.aln

We now want to get rid of sequences that are of too low quality. We exclude any sequences that display more than 1500 Ns

    perl ./scripts/Count_N_in_seq.pl -i ./data/02.total.mask.aln -o ./data/03.N_nb.tsv
    perl -F'\t' -ane 'if($F[-1] < 1500){print}' ./data/03.N_nb.tsv > ./data/03.N_nb_passed.tsv
    perl ./scripts/Extract_specific_sequences_from_fastq_or_fasta_v2.pl -c ./data/03.N_nb_passed.tsv -i ./data/02.total.mask.aln -f 0 -o ./data/04.total.N.aln -s "\t"


We also want to exclude sequences that contain too many SNPs: with a mean rate of evolution of 2 snp per month and our sampling spanning less than 5 months, we don't except more than 20 SNP


    perl ./scripts/Pairwise_SNP_distance_from_fasta_v3_two_files.pl -i ./data/EPI_ISL_402125.fasta -i2 ./data/04.total.N.aln -o ./data/05.total.snp_count.csv
    perl -F',' -ane 'if($F[-1] > 20){print }' ./data/05.total.snp_count.csv > ./data/06.total.high.snp_count.csv
    perl ./scripts/Extract_specific_sequences_from_fastq_or_fasta_v2.pl -c ./data/06.total.high.snp_count.csv -i ./data/04.total.N.aln -f 1 -r -o ./data/07.total_high_qual.aln -s "," 

At this stage the dataset comprises 596 high quality sequences (from 849 !).

**########## 4. Build a Maximum likelihood tree ##########**

Building a maximum likelihood phylogenetic tree on the masked alignment using the tree builder RaxML (https://cme.h-its.org/exelixis/web/software/raxml/) again implemented via the AUGUR pipeline.

`augur tree --alignment ./data/07.total_high_qual.aln --method raxml --nthreads 8 --output ./data/08.raxml.tree`


**########## 5. Phylogenetics reconstruction using R ##########**

This section is a modified version of Lucy van Dorp's SARS-CoV-2 Phylogenomics Practical

### Introduction 

In this practical, you will explore some of the R packages that can be used to analyse genome sequencing data. The intention is to give a grounding in the basics of reading, plotting and making inference from phylogenetic trees using genomes from SARS-CoV-2, the agent of COVID-19, as an example.


By the end of the practical, you should be able to: 

* Read a whole genome alignment into R
* Obtain metrics from this alignment 
* Plot a phylogenetic tree
* Estimate a temporal regression
* Perform a simple tip-dating analysis

The raw data used to apply these steps are available in the `./data` folder.

### Reading in a whole genome alignment

First, we will load in some R packages. These contain useful functions which supplement the basic R functions. If you do not have these, you will have to install them with e.g. `install.packages('adegenet')`.

```{r libraries, message=FALSE}
library(adegenet) # Population genetics R package containing tools for extracting SNP counts and assessing diversity over the alignment
library(adephylo) # Phylogenetics R package containing tools for querying phylogenetic trees
library(ape) # The most widely used package for reading, writing, customising and plotting phylogenetic trees
library(treedater) # Method for estimating mutation rates and divergence estimates over phylogenetic trees
library(lubridate) # Tool for manipulating dates from calender format to decimal formats.

setwd("/home/folder_of_your_choice/")
```

We will then read in the fasta format alignment:

```{r align}
sequences <- read.dna('./data/07.total_high_qual.aln',format='fasta')
```

Some explanation of the above code: 

* `read.dna` is a R function (i.e. comes with the R package Ape) which reads in fasta format alignment.
* `format='fasta'` is an argument which tells `read.dna` that this is an aligned fasta file.

The resulting alignment is stored as a `DNAbin` object. Printing this object gives us some information about the alignement.
  
```{r inspect_dataset_1, message=FALSE}
sequences
```
  
  Now lets look at the samples in the alignment (these are the rows of the alignment, one aligned sequence per row). We can use the function `head()` to view the first ten rows:
```{r inspect_dataset_2}
head(rownames(sequences))
```

We can check how many positions there are in the alignment by querying its dimensions using the function `dim()`. Dim reports dimensions, the first number being the number of rows and the second the number of columns. In this case the number of columns corresponds to the number of nucleotides. 



* We have profile aligned to the SARS-CoV-2 reference genome Wuhan-Hu-1. This genome is 29903 bases in length so the alignment should be exactly this long too.
</details>

```{r inspect_dataset_3}
dim(sequences)
```

SARS-CoV-2 is actually large for a virus. For example the Influenza A virus genome is only around 13,588 bases.
    
We can also ask how many positions in the alignment vary - the single nucleotide polymorphisms (SNPs). We can do this using a function from adegenet `seg.sites()` which reports all variant positions. We can assess the length of this string using `length()`:
      
```{r inspect_dataset_4}
nb_SNP <- length(seg.sites(sequences))
nb_SNP
```

This tells us the number of positions in the 29903 base pair alignment vary in at least one of the included sequences. In fact in this alignment the mean SNP difference between any two assemblies is **only 8 mutations**. Remember even today SARS-CoV-2 is a young and genetically homogenous virus.

Storing the alignment like this also allows us to query any particular position. 

A very commonly discussed mutation in the context of the early SARS-CoV-2 pandemic is at position 23403 corresponding to a protein change D614G in the SARS-CoV-2 spike protein. Lets see how many of the samples in this alignment vary at this position. We can do this by querying the position and computing the base frequencies using the `base.freq()` function. Lets just look at the first 10 sequences using `head()`:
  
```{r inspect_dataset_5}
base.freq(sequences[,23403],freq=T)
```

Some explanation of the above code: 

* `sequences[,23403`] queries this position in the DNAbin object (our alignment).
* `freq=T` is an argument to `base.freq()` which requests the absolute counts of A, C, G and T in the alignment to be returned (rather than their proportions)

* 140 of the genomes have an `a` at this position (the ancestral type). 
* 452 genomes have a `g` at the position, which correponds to the nonsynonymous change to 614G. * This corresponds to ~30% of the included genomes carrying the mutation.
</details>

Now the frequency of D614G is essentially fixed. This is likely due to a combination of founder introductions of viruses carrying this mutation and the functional impact of this mutation in modulating transmissibility.

### Reading in metadata
Sequence data is only useful for epidemiological inference if we have **associated metadata on where and when genomes were collected**. Unfortunately we often don't have information on patient status, though this would be useful. 

We can read the metadata table into R to see where the samples come from and use this to understand better the emergence and spread of the virus.

A common data format is **t**ab **s**eperated **v**alue (.tsv) format because it can both be opened in Excel but also easily read by different programming languages in R. We can read in tab separated files using the `read.delim()` function in R: 

```{r read_meta}
meta <- read.delim("Table_1_metadata.tsv",header=T,sep='\t',as.is=T)
```

Some explanation of the above code: 

* `read.delim` is a base R function which reads tsv files.
* `header='T'` specifies that we have a header in the .tsv file.
* `sep='\t'` specifies that entries in the file are seperated by a tab (encoded by `\t` in R). Note if we had a **c**omma **s**eperated **v**alue file (.csv) we could specify that `sep=','`.
* `as.is=T` ensures factors are read in as characters.

We can have a look at the first ten rows of the file using `head()`:

```{r view_meta}
head(meta)
```

We can query any of the variables in the metadata by pulling off different variables of the data frame using the dollar symbol `$`. 

Lets take a look at the Country represented. We can do this using the `$` sign and the `table()` command as follows:  

```{r explore_meta_1}
table(meta$Country)
```

At this early stage of the pandemic most samples were from Asia and Europe. However this rapidly shifted. European sequencing has typically dominated sequencing efforts. This is largely due to the COG-UK consortium (https://www.cogconsortium.uk).



The "country" field of the metadata table needs to be modified as it also contain more precise information (eg. USA: Connecticut ). Although this is surprinsingly not the case here, metadata tables typically contain missplells and errors in the country names and those need to be searched for and corrected manually. Remember `R` can only assess the variables it reads so does not consider categories with typos or different spellings (eg. capitalised and not capitalised) as the same level (e.g. "USA" vs. "usa" vs. "U.S.A" vs. "United States of America" ...).

```
meta$Country <- sapply(strsplit(meta$Country,":"), `[`, 1)
```

The metadata table doesn't contain any continent information so we will add it
```
countries_continents <- read.delim("./data/countries_continents.csv",header=T,sep=',',as.is=T)
meta$Continent <- countries_continents$Continent[match(meta$Country,countries_continents$Country)]

table(meta$Continent)
```

We now attribute a color to each continent and add this information to the metadata table
```
continent_colors <- c('pink','aquamarine','cornflowerblue','darkmagenta','brown1','darkgoldenrod')
names(continent_colors) <- c('Africa','Asia','Europe','North America','South America','Oceania')

meta$Colour <- continent_colors[meta$Continent]
```

### Reading in a phylogenetic tree

We will next use the Ape package from R to read in a maximum likelihood (ML) phylogenetic tree over the 596 SARS-CoV-2 genome assemblies. Trees are often stored as `Newick` file types, which are strings of text which denote the tree structure (via brackets) and the branch lengths measured in substitutions per site.
  
We can do this using the function `read.tree()`:
  
```{r read_tree}
mltree <- read.tree("./data/08.raxml.tree")
```
  
This is now stored as a tree object. Lets see what it looks like:
  
```{r query_tree}
mltree
```  

R efficiently stores this tree as a `phylo format` object. You can see this tree has 596 tips and each tip corresponds to a genome in our dataset. We can look at the order of tips by extracting variables from this object. For example to print the first ten tip labels:
  
```{r query_tip_tree}
head(mltree$tip.label)
```   

### Excluding low-quality data

Gathering dataset from distinct laboratories on public databases often yields datasets very heterogeneous in quality. It's therefore important to check the quality of the sequences. In this case they have been prefiltered before starting the R analysis in order to be able to create the phylogenetic tree, but this is still worth a look.


```
# We start by creating an object that will receive the identifiers and the number of N of each sequence
quality_of_seqs <- data.frame(ncol=2,nrow=nrow(sequences))

# and we run a loop on each sequence to compute its proportion of N
for(i in seq(nrow(sequences))) {
  # Get prop. of Ns
  prop.Ns <- base.freq(sequences[i,], all = TRUE)[["n"]]
 quality_of_seqs[i,] <- c(rownames(sequences)[i],prop.Ns)
}
hist(as.numeric(quality_of_seqs[,2]),breaks=200,xlab="Proportion of Ns in each sequence")
```

### Rooting the phylogenetic tree

It is important to root phylogenetic trees so as to orientate the topology. This is particularly relevant for phylogenetic dating which relies on **root-to-tip distances** to compute temporal signal. Typically we would use a phylogenetic `out-group`. This is a genome sequence which is known to be distant from all samples in the tree but not too distant so as to not share large parts of the alignment.
  
However, for SARS-CoV-2 there remains no appropriate out-group genome.

* the closest available genomes from bats remain **decades divergent** from SARS-CoV-2. 
* they are too distant to be used as out-groups and the genomes are fairly low quality.

Instead we use the reference (and first) sequence of SARS-CoV-2 **Wuhan-Hu-1**. 
  
We can root the tree using the `root` function available in the Ape package:
  
```{r root_tree}
mltree.root <- root(mltree,"EPI_ISL_402125",resolve.root = T)
```

### Plotting a phylogenetic tree
Lets have a first go at plotting the tree. As we installed Ape we can use the `plot` command which will already interpret that we have a `phylo object`. To speed up the plot time we chose not to display tip-labels and use `show.tip.label=FALSE`:
  
```{r plot_tree_1}
plot(mltree.root,show.tip.label=FALSE)
```

Trees are not very useful without annotation. We previously added a colour for each continent to the metadata table. By matching the tip labels of the tree to the metadata using the `match()` function we can create a colour vector for all of the tips in the tree and provide this to the plot function.

```{r plot_tree_2}
col.vec <- meta$Colour[match(mltree.root$tip.label,meta$SRA_Run)]
plot(mltree.root,show.tip.label=FALSE)
tiplabels(pch=20,col=col.vec,frame='none')
```

Even at this early stage of the pandemic there is not much in the way of geographic structure. What does this mean? 
* the early pandemic first wave demonstrates that many samples from world-wide regions can be found across the global phylogenetic tree.
* you can think of the SARS-CoV-2 population as panmictic - any sample can represent any geographic region
* this is consistent with many introductions of the virus to many geographic regions at multiple times meaning there is no one patient zero.
* though note this situation changed by the time of the second COVID-19 wave when restrictions on travel supporting the emergence of more local SARS-CoV-2 clusters.

![schematics](https://github.com/END-VOC/WP3_D3.2_SARS-CoV-2_NGS_bioinformatics_pipeline/blob/main/tree_1.png)

### Dating a phylogenetic tree

Many uncertainties surround the age of SARS-CoV-2. Using genomic datasets we can estimate the **time to the most recent common ancestor (tMRCA)** of the sequences we include in our alignment.

First we need to obtain the sampling dates of all of the samples in our alignment. We can do this by querying the metadata dataframe for the variable collection date (`meta$Date`), which we store as a variable `date.vec`:
  
```{r date_1}
date.vec <- meta$Date[match(mltree.root$tip.label,meta$SRA_Run)]

```

We already stored the collection dates as `date.vec` but most tools can not read dates in this format. We therefore convert calender dates to decimal dates. There is a useful package for doing this called `lubridate` in R. Note that sometimes date entires are entered wrong - often months and days are flipped - it is important to check this carefully.

```{r date_3}
date.dec.vec <- decimal_date(as.Date(date.vec,'%Y-%m-%d'))
names(date.dec.vec) <- mltree.root$tip.label
```

There are many possible tools for **phylogenetic tip-dating**. For this dataset we will use an approach called **TreeDater** which fits a strict or relaxed molecular clock to a phylogenetic tree and estimates evolutionary rates and times of common ancestry. It is convenient as it is readily implementable as an R package.

You can read more in the TreeDater paper here: https://academic.oup.com/ve/article/3/2/vex025/4100592

For bespoke tools like TreeDater learning how to apply them requires careful reading of associated documentation. TreeDater has an associated manual available here: https://cran.r-project.org/web/packages/treedater/index.html

We can see from this manual that we require as input our rooted phylogenetic tree `mltree.root`, the associated dates in decimal years `date.dec.vec` and the number of SNPs informing the phylogenetic tree (calculated earlier and stored in the "nb_SNP" variable). 

For simiplicity we will assume a **strict clock model**. This is a reasonable assumption for the earliest SARS-CoV-2 samples but other more complex models exist and there is now some evidence of the rate of evolution changing subtly over the course of the pandemic.

We can then run TreeDater. We store the output as a variable termed `dtr`. Note this may take five minutes or so depending on your laptop:
```{r date_4}
dtr <- dater(mltree.root,date.dec.vec,nb_SNP,clock='strict')
```

Before interpreting the results lets check we have meaningful temporal signal in the alignment by performing a **simple regression** of root-to-tip distances and sampling time. TreeDater handily provides a function for this which can be applied to the `dtr` object:

```{r date_5}
rootToTipRegressionPlot(dtr)
```

Our regression has a p-value of <1e-22. 

The p-value of the temporal regression is highly significant. What does this mean? 
* the regression analysis supports measurably evolution over the time-span of our dataset
* temporal signal is often not this strong due to eg. noisy alignments, erroneous SNP calls, and (more likely) the effect of recombination and/or structural variants
* note there is now good evidence for recombination in SARS-CoV-2 which is in line with the behaviour of other coronaviruses.


We can therefore interpret the results of the analysis by querying elements of the `dtr` object. The one we are most interested in is the **estimated age** of all of these genomes:

```{r date_6}
dtr$timeOfMRCA
```

Our analysis suggests a decimal date of 2019.997. This translates to a calender estimate of the **30th December 2019**.

This is a bit more recent than we might expect. Some of the earliest cases reported in China date to mid November 2019. 

 Why might our inferred date be too recent?

* we are only considering a small number of genomes from a large viral resevoir early in the pandemic. Incorporating more diversity may alter the estimated rates.
* this analysis employs a strict clock model (constant rate). Allowing for relaxation of these parameters may shift the date, with the likely outcome being to older times.
* this is only a point estimate and estimates such as these should always be presented with appropriate uncertainty.
* phylgenetic dating analyses capture this uncertainty either through boot-strapping or by sampling sets of posterior trees (the latter is the approach used in the BEAST framework which you may have come across).

![schematics](https://github.com/END-VOC/WP3_D3.2_SARS-CoV-2_NGS_bioinformatics_pipeline/blob/main/root_to_tip.png)

### Plotting a dated phylogenetic tree

Finally we can plot the phylogenetic tree with branch lengths now **scaled by unit time** rather than by the number of substitutions. This aids interpretability, for example when we have identical sequences that were sampled weeks a part. Actually this is very common for SARS-CoV-2 which **mutates more slowly than it transmits**.

Plotting a dated tree can be achieved using the function `axisPhylo()` and specifying the results from the TreeDater `dtr` object analysis. We can colour the tip points by continent just as before:

```{r plot_date}
plot(dtr,no.mar=F,cex=0.2,show.tip.label=F)
axisPhylo(side=1,root.time=dtr$timeOfMRCA,backward=F)
tiplabels(pch=20,col=col.vec)
legend("topleft",inset=c(0.1,0),c('Africa','Asia','Europe','NorthAmerica','SouthAmerica','Oceania'),col=c('pink','aquamarine','cornflowerblue','darkmagenta','brown1','darkgoldenrod'),pch=20,bty='n',cex=1.2)
```

![schematics](https://github.com/END-VOC/WP3_D3.2_SARS-CoV-2_NGS_bioinformatics_pipeline/blob/main/tree_2.png)

You can now see that the tree is scaled by unit time. By looking at the x-axis we are able to estimate both the tMRCA but also the **approximate age of different nodes** on the phylogenetic tree. This can be useful for estimating the age of introduced lineages as well as the age of the appearance of particular sets of mutations. Analyses such as these suggest D614G first emerged in mid January 2019.

To see an up to date dated phylogenetic tree on a subsampled dataset of SARS-CoV-2 you can check the NextStrain web interface https://nextstrain.org/ncov/global?c=region.


### Conclusions

Hopefully this has served as an introduction to the flexibility of R for analysing phylogenetic data. This analysis demonstrates the feasibility of estimating the age of SARS-CoV-2 using a small dataset encompassing the genomes uploaded during the first three months of the COVID-19 genomics sequencing efforts.

Phylogenetic datasets of SARS-CoV-2 continue to grow and one challenge is scaling the kinds of methods presented here to hundreds of thousands of viral sequences.

To sum up, we have been able to:

* Read in and evalute the diversity of an alignment
* Read in, root and plot a phylogenetic tree
* Perform a simple phylogenetic tip-dating analysis

### Extensions

The commands used should all be quite portable and flexible. If you're interested, I encourage you to try them for yourself with variations. 

There are also many further tools that can be used in R for phylogenetics:

* **ggtree** - an impressive suite of tools for visulisation and plotting of phylogenetic trees allowing for additional customisations that are not possible in Ape.

* **BactDating** - another R package for phylogenetic tip-dating.

* **phangorn** - another widely used R package for phylogenetic analyses and visulations.

There are also many online tutorials. A nice phylogenetics tutorial in R which touches on some of the concepts mentioned in this practical is available here:
  * **Introduction to phylogenies in R** (http://www.phytools.org/Cordoba2017/ex/2/Intro-to-phylogenies.html)

And if you're excited by phylogenetic analyses of SARS-CoV-2 some of our recent papers (also mentioned in the lecture) are below:

* Emergence of genomic diversity and recurrent mutations in SARS-CoV-2 (https://www.sciencedirect.com/science/article/abs/pii/S1567134820301829)

* No evidence for increased transmissibility from recurrent mutations in SARS-CoV-2 (https://www.nature.com/articles/s41467-020-19818-2)

* The genomic epidemiology of SARS-CoV-2 in Palestine (https://www.biorxiv.org/content/10.1101/2020.10.26.355677v3)

* Recurrent mutations in SARS-CoV-2 genomes isolated from mink point to rapid host-adaptation (https://www.biorxiv.org/content/10.1101/2020.11.16.384743v1)

* Pre-existing T cell-mediated cross-reactivity to SARS-CoV-2 cannot solely be explained by prior exposure to endemic human coronaviruses (https://www.biorxiv.org/content/10.1101/2020.12.08.415703v1)

