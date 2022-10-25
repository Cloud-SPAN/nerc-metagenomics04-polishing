---
title: "Polishing an assembly"
teaching: 30
exercises: 10
questions:
- "Why do assemblies need to be polished?"
- "What are the different purposes of polishing with short and long reads?"
- "What software can we used to do long and short read polishing?"

objectives:
- "Understand why polishing metagenomes is important."  
- "Understand the different programs used to do short and long read polishing."
keypoints:
- "Short reads have a higher base accuracy than long reads and can be used to remove errors in assemblies generated with long reads."
- "Long reads have a lower accuracy but help generate a more contiguous (less fragmented) assembly, so are used to get the structure of the metagenome, but may have small misassemblies or single nucleotide polymorphisms (SNPs)"
- "Medaka is used to polish an assembly with long reads."
- "Pilon is used to polish an assembly with short reads."
---



In the [previous episode](https://cloud-span.github.io/metagenomics01-qc-assembly/03-assembly/index.html) we generated a draft assembly using Flye from our long read Nanopore data. 

Long reads can span regions which would would difficult to assemble with short reads such as regions with large repeats. Despite this, some long reads will be misassembled. In addition, the base accuracy of long reads is lower than that short reads and some bases will be incorrectly assigned. Consequently it is common to correct these errors known as "polishing" an assembly. We will use two polishing strategies:
1. Polishing with long reads using [Medaka](https://github.com/nanoporetech/medaka). This maps the raw long reads to the assembly to identify contigs that have been incorrectly joined.
2. Polishing with short reads using [Pilon](https://github.com/broadinstitute/pilon) which uses the more accurate short read data to correct incorrectly called bases in the assembly.
More detail on the advantages and disadvantages of short and long read sequencing is covered in our [Statistically useful experimental design](https://cloud-span.github.io/experimental_design00-overview/) workshop in [Platform choice](https://cloud-span.github.io/experimental_design01-principles/01-platform/index.html).

<img align="left" width="525" height="307" src="{{ page.root }}/fig/04_polishing_diagram_v1.png" alt="Diagram showing overlap of reads for polishing" /> &nbsp; &nbsp; &nbsp;

**Polishing with short reads**
In the diagram, the long read assembly is shown at the top. The four short reads shown below the assembly have been aligned to it. The third position in the assembly is `A` but the three short reads that cover this region contain a `T`. The assembly probably contains a miscalled base but it can be corrected - or polished - at this position with the higher accuracy short read information. Typically short read polishing would be carried out three times with the base accuracy increasing each time. However, to reduce the compute requirements and the time required to finish the assembly, we will be performing it just once.

<br clear="left"/>

## Why bother polishing?
How important polishing is to your analysis will depend on what you need it for. Usually we generate metagenome assemblies so we can compare the sequences to a database and find out what taxa they belong to.

You might NOT need to polish your assembly if:
- you only need to taxa to the genus level (meaning single incorrect bases are not important)

You DO need to polish your assembly if:
- you want to identify taxa to the species level (if possible). This is a common requirement since one of the main advantages of whole genome sequencing over amplicon sequencing is that you can assign annotations to the species level.  We will cover [Taxonomic annotations](https://cloud-span.github.io/metagenomics03-taxonomic-anno/) later in the course.
- you want to generate protein predictions or identify protein structure domains to determine the functionality of metagenomes. This is discussed in more detail in Watson and Warr (2019): [Errors in long-read assemblies can critically affect protein prediction](https://www.nature.com/articles/s41587-018-0004-z).   

## Polishing an assembly with long reads

First we will polish the draft Flye assembly using the filtered raw long reads. As with the assembly, we need to use polishing software that is especially written for long read raw reads.

[Medaka](https://github.com/nanoporetech/medaka) is a command line tool built by Oxford Nanopore Technologies which will polish an assembly by generating a consensus from raw Nanopore sequences using a recurrent neural network.

We will be using one Medaka command, `medaka_consensus`. This pipeline will first align the raw reads to the draft assembly, then process this alignment to generate a pileup. The pileup is presented to a recurrent neural network in order to produce a consensus sequence.

Medaka is installed on the AWS instance. Look at the help page for `medaka_consensus`:

~~~
medaka_consensus -h
~~~
{: .bash}

> ## `medaka_consensus` Help
>
> ~~~
> medaka 1.7.0
>
> Assembly polishing via neural networks. Medaka is optimized
> to work with the Flye assembler.
>
> medaka_consensus [-h] -i <fastx> -d <fasta>
>
>     -h  show this help text.
>     -i  fastx input basecalls (required).
>     -d  fasta input assembly (required).
>     -o  output folder (default: medaka).
>     -g  don't fill gaps in consensus with draft sequence.
>     -r  use gap-filling character instead of draft sequence (default: None)
>     -m  medaka model, (default: r941_min_hac_g507).
>         Choices: r103_fast_g507 r103_hac_g507 r103_min_high_g345 r103_min_high_g360 r103_prom_high_g360 r103_sup_g507 r1041_e82_400bps_fast_g615 r1041_e82_400bps_hac_g615 r1041_e82_400bps_sup_g615 r104_e81_fast_g5015 r104_e81_hac_g5015 r104_e81_sup_g5015 r104_e81_sup_g610 r10_min_high_g303 r10_min_high_g340 r941_e81_fast_g514 r941_e81_hac_g514 r941_e81_sup_g514 r941_min_fast_g303 r941_min_fast_g507 r941_min_hac_g507 r941_min_high_g303 r941_min_high_g330 r941_min_high_g340_rle r941_min_high_g344 r941_min_high_g351 r941_min_high_g360 r941_min_sup_g507 r941_prom_fast_g303 r941_prom_fast_g507 r941_prom_hac_g507 r941_prom_high_g303 r941_prom_high_g330 r941_prom_high_g344 r941_prom_high_g360 r941_prom_high_g4011 r941_prom_sup_g507 r941_sup_plant_g610
>         Alternatively a .tar.gz/.hdf file from 'medaka train'.
>     -f  Force overwrite of outputs (default will reuse existing outputs).
>     -x  Force recreation of alignment index.
>     -t  number of threads with which to create features (default: 1).
>     -b  batchsize, controls memory use (default: 100).
> ~~~
> {: .output}
{: .solution}

* The usage is `medaka_consensus [-h] -i <fastx> -d <fasta>` indicating that the `-i` and `-d` flags are mandatory.
  - `-i` indicates the input basecalls _i.e._, the Nanopore raw-reads, what we are polishing with.
  - `-d` indicates the assembly we are polishing
* Other flags are optional
  - `-m` allows to select an apropriate recurrent neural network model. The [documentation](https://github.com/nanoporetech/medaka#models) describes the models which are named to indicate i) the pore type, ii) the sequencing device (MinION or PromethION), iii) the basecaller variant, and iv) the basecaller version, with the format: `{pore}_{device}_{caller variant}_{caller version}`. Medaka doesn't offer an exact model for our dataset. We will use the closest available model: `r941_prom_fast_g303`. It is also possible specify a bespoke model.  
  - `-o` allows us to specify the output directory  
  - `-t` allows us to specify the number of threads so we can speed the process up

The `medaka_consensus` polishing will take about 30 mins so we will run it in the background and redirect the output to a file.
Make sure you are in the `analysis` folder and run `medaka_consensus` on `assembly.fasta`:
~~~
cd analysis/
medaka_consensus -i ../data/nano_fastq/ERR3152367_sub5_filtered.fastq -d assembly/assembly.fasta -m r941_prom_fast_g303 -o medaka -t 8 &> medaka.out &
~~~
{: .bash}
We have added `&> medaka.out &` to redirect the output to `medaka.out `and run the command in the background.  

We can check the command is running using `jobs`:
~~~
jobs
~~~
{: .bash}

If it is successfully running you should see an output like:
~~~
[1]+  Running                 medaka_consensus -i ../data/nano_fastq/ERR3152367_sub5_filtered.fastq -d assembly/assembly.fasta -m r941_prom_fast_g303 -o medaka -t 8 &> medaka.out &
~~~
{: .output}

We can also look in the output file (`medaka.out`) to check the progress of the command.
~~~
less medaka.out
~~~
{: .bash}
If the Medaka command has been run correctly you will see something like this at the start of the output:
~~~
Checking program versions
This is medaka 1.7.0
Program    Version    Required   Pass     
bcftools   1.15.1     1.11       True     
bgzip      1.15.1     1.11       True     
minimap2   2.24       2.11       True     
samtools   1.15.1     1.11       True     
tabix      1.15.1     1.11       True     
Aligning basecalls to draft
Constructing minimap index.
[M::mm_idx_gen::0.515*0.99] collected minimizers
[M::mm_idx_gen::0.648*1.40] sorted minimizers
[M::main::0.877*1.29] loaded/built the index for 146 target sequence(s)
[M::mm_idx_stat] kmer size: 15; skip: 10; is_hpc: 0; #seq: 146
[M::mm_idx_stat::0.910*1.28] distinct minimizers: 2598367 (94.68% are singletons); average occurrences: 1.076; average spacing: 5.350; total length: 14953273
~~~
{: .output}

> ## Help!
> Medaka may give you a warning along the lines of:
> ~~~
> 2022-10-25 09:07:35.970532: W tensorflow/stream_executor/platform/default/dso_loader.cc:64] Could not load dynamic library 'libcudart.so.11.0'; dlerror: libcudart.so.11.0: cannot open shared object file: No such file or directory
> 2022-10-25 09:07:35.970583: I tensorflow/stream_executor/cuda/cudart_stub.cc:29] Ignore above cudart dlerror if you do not have a GPU set up on your machine.
> 2022-10-25 09:07:39.935310: W tensorflow/stream_executor/platform/default/dso_loader.cc:64] Could not load dynamic library 'libcudart.so.11.0'; dlerror: libcudart.so.11.0: cannot open shared object file: No such file or directory
> 2022-10-25 09:07:39.935346: I tensorflow/stream_executor/cuda/cudart_stub.cc:29] Ignore above cudart dlerror if you do not have a GPU set up on your machine.
> ~~~
> {: .output}
> Don't worry, you can safely ignore this warning wherever it pops up. It is telling us that it couldn't load a library required for parallel computing using GPUs. We are not using a GPU setup and so this warning is irrelevant. 
{: .callout} 

Medaka first looks for the other programs that it needs (known as dependencies) and their versions. These dependencies are installed on the AWS instance. Once it confirms they are present, it begins by aligning the raw reads (basecalls) to the assembly using minimap.

<kbd>q</kbd> will quit from `less`.

Once Medaka has completed the end of the file will contain something like:
~~~
less medaka.out
~~~
{: .bash}

<kbd>G</kbd>  will take you to the end   
~~~
[20:51:16 - DataIndx] Loaded 1/1 (100.00%) sample files.
[20:51:16 - DataIndx] Loaded 1/1 (100.00%) sample files.
[20:51:16 - DataIndx] Loaded 1/1 (100.00%) sample files.
[20:51:16 - DataIndx] Loaded 1/1 (100.00%) sample files.
[20:51:16 - DataIndx] Loaded 1/1 (100.00%) sample files.
[20:51:16 - DataIndx] Loaded 1/1 (100.00%) sample files.
[20:51:16 - DataIndx] Loaded 1/1 (100.00%) sample files.
[20:51:17 - DataIndx] Loaded 1/1 (100.00%) sample files.
[20:51:17 - DataIndx] Loaded 1/1 (100.00%) sample files.
[20:51:17 - DataIndx] Loaded 1/1 (100.00%) sample files.
Polished assembly written to medaka/consensus.fasta, have a nice day.
~~~
{: .output}

Once Medaka has completed we can navigate into the output directory and look at the files it has generated.

~~~
cd medaka
ls
~~~
{: .bash}

~~~
calls_to_draft.bam  calls_to_draft.bam.bai  consensus.fasta  consensus.fasta.gaps_in_draft_coords.bed  consensus_probs.hdf
~~~
{: .output}

Medaka has created multiple files:

* `calls_to_draft.bam` - a BAM file containing the alignment of the raw reads (basecalls) to the draft assembly
* `calls_to_draft.bam.bai` - an index file of the above BAM file
* `consensus.fasta` - the consensus sequence, or polished assembly in our case in FASTA format
* `consensus.fasta.gaps_in_draft_coords.bed` - a BED file containing information about the location of any gaps in the consensus sequence which can be used when visualising the assembly
* `consensus_probs.hdf` - a file that contains the output of the neural network calculations and is not an output for end-users, so we don't need to worry about this file

In our case we're interested in the polished assembly, so we want the `consensus.fasta` file.

> ## BAM and SAM Files
> A [SAM file](https://genome.sph.umich.edu/wiki/SAM), is a tab-delimited text file that contains information for each individual read and its alignment to the genome. The paper by [Heng Li et al.](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) provides the full specification.
>
> The compressed binary version of SAM is called a BAM file. We use this version to reduce size and to allow for indexing, which enables efficient random access of the data contained within the file.  
>
> The file begins with a header, which is optional. The header describes the source of data, reference sequence, method of alignment etc. - these will change depending on the aligner being used. Following the header is the alignment section. Each line that follows corresponds to alignment information for a single read. There are 11 mandatory fields for essential mapping information and a variable number of other fields for aligner specific information.
>
> See [Genomics - Variant Calling](https://cloud-span.github.io/04genomics/01-variant_calling/index.html) for a deeper dive.  
{: .callout}

## Polishing with short reads

We will be using the program [Pilon](https://github.com/broadinstitute/pilon) to further polish the draft assembly using the raw short reads. Pilon will improve a draft assembly by filling gaps, fixing misassemblies and correcting bases. You can read more about how it works in the paper [Pilon: An Integrated Tool for Comprehensive Microbial Variant Detection and Genome Assembly Improvement](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0112963).

Bioinformatics programs are not built equally. Some programs, like Flye or Medaka, require very few input files as they will generate any that they need within the pipeline. Some programs however, require a lot of user input to generate the input files that are needed.  

Pilon is in the latter group of bioinformatics software, so we will need to do some pre-processing using other programs to create some of the inputs needed.

### Generating the Pilon input files

We will first use the program [BWA](https://github.com/lh3/bwa) to generate an alignment of the raw short reads against the draft genome in consensus.fasta. The steps will be:
1. indexing the polished assembly, consensus.fasta with `bwa index`. Indexing allows the aligner to quickly find potential alignment sites for query sequences in a genome, which saves time during alignment.
2. creating a directory for the outputs of Pilon
3. aligning the short reads (the illumina data) to the assembly, consensus.fasta with `bwa mem`
4. converting the short read alignment alignment to a BAM file `samtools view`
5. sorting the short read alignment with `samtools sort`
6. indexing the short read alignment with `samtools index`

Make sure are in the `analysis` folder and index the consensus assembly:
~~~
cd ~/cs_course/analysis/
bwa index medaka/consensus.fasta
~~~
{: .bash}
This should only take a few seconds to complete so we don't need to run the job in the background.
Once the indexing is complete you should see an output like:
~~~
[bwa_index] Pack FASTA... 0.51 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 5.86 seconds elapse.
[bwa_index] Update BWT... 0.10 sec
[bwa_index] Pack forward-only FASTA... 0.10 sec
[bwa_index] Construct SA from BWT and Occ... 1.81 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index consensus.fasta
[main] Real time: 8.704 sec; CPU: 8.395 sec
~~~
{: .output}

This will generate five additional files in the `medaka` directory with the file extensions `.pac`, `.bwt`, `.ann`, `.amb` and `.sa`. These files are used by BWA in step 3.

Next we will make an output directory and move into it:

~~~
mkdir pilon
cd pilon
~~~
{: .bash}

We will now do steps 3, 4 and 5 in one go by chaining them together with pipes.

> ## Chaining together commands with a pipe
> It is possible to chain together commands in unix using a process known as "piping". This allows the output from one command to be directly passed as input to the next without the need for intermediate files. This is useful when the intermediate file is not needed and keeps your workspace tidy (and unfull). The pipe character is `|` and obtainined with <kbd>⇧ Shift</kbd> + <kbd>\</kbd> on most keyboards.
>
> You can use multiple pipes in one command but data will only go from the left to the right:
>
> `command1 | command2 | command3 | .... |`
{: .callout}

We will be using two pipes to join three separate steps. First we will align the raw reads to the draft assembly, then convert the output to BAM format, before finally sorting this alignment to generate a sorted BAM file. Chaining the steps together together will only generate one final output file, avoiding the need to generate large intermediary files we don't need again between the other two steps.

3. we will align the short reads (the illumina data) to the assembly, consensus.fasta with `bwa mem`:
   `bwa mem -t 8 consensus.fasta ../data/illumina_fastq/ERR2935805.fastq`
4. convert the short read alignment alignment to a BAM file `samtools view`:
   `samtools view - -Sb`
5. sort the short read alignment with `samtools sort`:
   `samtools sort - -@4 -o short_read_alignment.bam`  

This will take around 60 minutes so we will use `&> alignment.out &` to redirect the process to a file and to run the command in the background. We will also wrap our whole command in brackets so we run all three steps in the background.

Add the pipes between these commands and run:
~~~
(bwa mem -t 8 ../medaka/consensus.fasta ../../data/illumina_fastq/ERR2935805.fastq | samtools view - -Sb | samtools sort - -@4 -o short_read_alignment.bam) &> alignment.out &
~~~
{: .bash}

Once the command is running, you can check the process of this job by looking at the `alignment.out` file
~~~
less alignment.out
~~~
{: .bash}

~~~
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 8758 sequences (40002050 bp)...
[M::process] read 8872 sequences (40001127 bp)...
[M::mem_process_seqs] Processed 8758 reads in 134.856 CPU sec, 33.785 real sec
[M::process] read 8774 sequences (40005021 bp)...
[M::mem_process_seqs] Processed 8872 reads in 138.347 CPU sec, 34.631 real sec
[M::process] read 8736 sequences (40009331 bp)...
[M::mem_process_seqs] Processed 8774 reads in 130.086 CPU sec, 32.545 real sec
[M::process] read 8848 sequences (40002710 bp)...
[M::mem_process_seqs] Processed 8736 reads in 132.912 CPU sec, 33.277 real sec
[M::process] read 8884 sequences (40009423 bp)...
[M::mem_process_seqs] Processed 8848 reads in 134.169 CPU sec, 33.883 real sec
[M::process] read 8902 sequences (40003755 bp)...
[M::mem_process_seqs] Processed 8884 reads in 129.038 CPU sec, 32.410 real sec
[M::process] read 8760 sequences (40000601 bp)...
~~~
{: .output}
Once completed, the end of the `alignment.out` file should contain something like:
~~~
[M::mem_process_seqs] Processed 8862 reads in 117.475 CPU sec, 29.369 real sec
[M::process] read 4795 sequences (23206538 bp)...
[M::mem_process_seqs] Processed 8610 reads in 125.241 CPU sec, 31.397 real sec
[M::mem_process_seqs] Processed 4795 reads in 87.225 CPU sec, 23.866 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -t 8 consensus.fasta ../ERR3152367_sub5_filtered.fastq
[main] Real time: 2454.659 sec; CPU: 9557.109 sec
[bam_sort_core] merging from 4 files and 4 in-memory blocks...
~~~
{: .output}
We have now generated the `short_read_alignment.bam` file - this is a binary file (meaning it's not human readable) so we won't be checking its contents.

Now carry out step 6, index the alignment:

~~~
samtools index short_read_alignment.bam
~~~
{: .bash}

> ## Something to think about
> Why didn't we include this command in the sequence of pipes in the previous step? The answer is that we will need access to the BAM file produced for further analysis.
> If we included this step as part of a pipe the intermediate BAM file would not be saved.
{: .callout}

This command will take around one or two minutes so we don't need to run it in the background.

Once your prompt has returned you should also have a file named `short_read_alignment.bam.bai` which is the index.

### Running Pilon
Now we have generated the necessary input files we can finally run Pilon.

First navigate up(`..`)  to the `analysis` directory
Pilon is installed on the AWS instance and we can view the help documentation using:

~~~
cd ..
pilon --help
~~~
{: .bash}

> ## Pilon help Documentation
> ~~~
> Pilon version 1.24 Thu Jan 28 13:00:45 2021 -0500
>
>    Usage: pilon --genome genome.fasta [--frags frags.bam] [--jumps jumps.bam] [--unpaired unpaired.bam]
>                 [...other options...]
>           pilon --help for option details
>
>
>         INPUTS:
>           --genome genome.fasta
>              The input genome we are trying to improve, which must be the reference used
>              for the bam alignments.  At least one of --frags or --jumps must also be given.
>           --frags frags.bam
>              A bam file consisting of fragment paired-end alignments, aligned to the --genome
>              argument using bwa or bowtie2.  This argument may be specified more than once.
>           --jumps jumps.bam
>              A bam file consisting of jump (mate pair) paired-end alignments, aligned to the
>              --genome argument using bwa or bowtie2.  This argument may be specified more than once.
>           --unpaired unpaired.bam
>              A bam file consisting of unpaired alignments, aligned to the --genome argument
>              using bwa or bowtie2.  This argument may be specified more than once.
>           --bam any.bam
>              A bam file of unknown type; Pilon will scan it and attempt to classify it as one
>              of the above bam types.
>           --nanopore ont.bam
>              A bam file containing Oxford Nanopore read alignments. Experimental.
>           --pacbio pb.bam
>              A bam file containing Pacific Biosciences read alignments. Experimental.
>         OUTPUTS:
>           --output prefix
>              Prefix for output files
>           --outdir directory
>              Use this directory for all output files.
>           --changes
>              If specified, a file listing changes in the <output>.fasta will be generated.
>           --vcf
>              If specified, a vcf file will be generated
>           --vcfqe
>               If specified, the VCF will contain a QE (quality-weighted evidence) field rather
>               than the default QP (quality-weighted percentage of evidence) field.
>           --tracks
>               This options will cause many track files (*.bed, *.wig) suitable for viewing in
>               a genome browser to be written.
>         CONTROL:
>           --variant
>              Sets up heuristics for variant calling, as opposed to assembly improvement;
>              equivalent to "--vcf --fix all,breaks".
>           --chunksize
>              Input FASTA elements larger than this will be processed in smaller pieces not to
>              exceed this size (default 10000000).
>           --diploid
>              Sample is from diploid organism; will eventually affect calling of heterozygous SNPs
>           --fix fixlist
>              A comma-separated list of categories of issues to try to fix:
>                "snps": try to fix individual base errors;
>                "indels": try to fix small indels;
>                "gaps": try to fill gaps;
>                "local": try to detect and fix local misassemblies;
>                "all": all of the above (default);
>                "bases": shorthand for "snps" and "indels" (for back compatibility);
>                "none": none of the above; new fasta file will not be written.
>              The following are experimental fix types:
>                "amb": fix ambiguous bases in fasta output (to most likely alternative);
>                "breaks": allow local reassembly to open new gaps (with "local");
>                "circles": try to close circular elements when used with long corrected reads;
>                "novel": assemble novel sequence from unaligned non-jump reads.
>           --dumpreads
>              Dump reads for local re-assemblies.
>           --duplicates
>              Use reads marked as duplicates in the input BAMs (ignored by default).
>           --iupac
>              Output IUPAC ambiguous base codes in the output FASTA file when appropriate.
>           --nonpf
>              Use reads which failed sequencer quality filtering (ignored by default).
>           --targets targetlist
>              Only process the specified target(s).  Targets are comma-separated, and each target
>              is a fasta element name optionally followed by a base range.
>              Example: "scaffold00001,scaffold00002:10000-20000" would result in processing all of
>              scaffold00001 and coordinates 10000-20000 of scaffold00002.
>              If "targetlist" is the name of a file, each line will be treated as a target
>              specification.
>           --verbose
>              More verbose output.
>           --debug
>              Debugging output (implies verbose).
>           --version
>              Print version string and exit.
>         HEURISTICS:
>           --defaultqual qual
>              Assumes bases are of this quality if quals are no present in input BAMs (default 10).
>           --flank nbases
>              Controls how much of the well-aligned reads will be used; this many bases at each
>              end of the good reads will be ignored (default 10).
>           --gapmargin
>              Closed gaps must be within this number of bases of true size to be closed (100000)
>           --K
>              Kmer size used by internal assembler (default 47).
>           --mindepth depth
>              Variants (snps and indels) will only be called if there is coverage of good pairs
>              at this depth or more; if this value is >= 1, it is an absolute depth, if it is a
>              fraction < 1, then minimum depth is computed by multiplying this value by the mean
>              coverage for the region, with a minumum value of 5 (default 0.1: min depth to call
>              is 10% of mean coverage or 5, whichever is greater).
>           --mingap
>              Minimum size for unclosed gaps (default 10)
>           --minmq
>              Minimum alignment mapping quality for a read to count in pileups (default 0)
>           --minqual
>              Minimum base quality to consider for pileups (default 0)
>           --nostrays
>              Skip making a pass through the input BAM files to identify stray pairs, that is,
>              those pairs in which both reads are aligned but not marked valid because they have
>              inconsistent orientation or separation. Identifying stray pairs can help fill gaps
>              and assemble larger insertions, especially of repeat content.  However, doing so
>              sometimes consumes considerable memory.
> ~~~
> {: .output}
{: .solution}

You can read more about the possible outputs Pilon can produce in the [Wiki](https://github.com/broadinstitute/pilon/wiki/Output-File-Descriptions).

We can see there are many different options for pilon. We will be using the defaults for our assembly.
* `--genome` - this will be the output assembly from Medaka
* `--unpaired` - the short reads we used to create the BAM alignment were unpaired, so we need to specify this using this flag
* `--outdir` - this will generate a directory for all the output

Check you are in the `analysis` folder and run

~~~
cd ~/cs_course/analysis/
pilon --genome medaka/consensus.fasta --unpaired pilon/short_read_alignment.bam --outdir pilon &> pilon.out &

~~~
{: .bash}

We can again keep track of the analysis by looking at the `pilon.out` file with `less`.

The top of the file:  
~~~
Pilon version 1.24 Thu Jan 28 13:00:45 2021 -0500
Genome: ../medaka/consensus.fasta
Fixing snps, indels, gaps, local
Input genome size: 14961385
Scanning BAMs
short_read_alignment.bam: 97413971 reads, 0 filtered, 97027227 mapped, 0 proper, 0 stray, Unpaired 100% 147+/-55, max 313
Processing contig_47:1-12303
unpaired short_read_alignment.bam: coverage 3
Total Reads: 532, Coverage: 3, minDepth: 5
Confirmed 3506 of 12303 bases (28.50%)
Corrected 25 snps; 0 ambiguous bases; corrected 18 small insertions totaling 22 bases, 28 small deletions totaling 32 bases
# Attempting to fix local continuity breaks
# fix break: contig_47:10552-11713 0 -0 +0 NoSolution
contig_47:1-12303 log:
Finished processing contig_47:1-12303
~~~
{: .output}

When Pilon finishes (aroud 20 mins) the end of the file will contain something like:
~~~
Writing updated contig_86_pilon to pilon.fasta
Writing updated contig_90_pilon to pilon.fasta
Writing updated contig_30_pilon to pilon.fasta
Writing updated contig_136_pilon to pilon.fasta
Writing updated contig_107_pilon to pilon.fasta
Mean unpaired coverage: 535
Mean total coverage: 535
~~~
{: .output}

Navigate into the `pilon` directory and have a look at the output files Pilon has produced.
~~~
cd pilon
ls
~~~
{: .bash}
~~~
alignment.out  pilon.fasta short_read_alignment.bam  short_read_alignment.bam.bai
~~~
{: .output}

We can see pilon has produced a fasta file `pilon.fasta`, which is the newly polished assembly.
This file is now our assembly.

In the next episode we will assess the quality of this assembly and compare its quality to that of the unpolished assemly.

> ## Recommended reading:
> While you're waiting for the polishing to finish, here's some things you might want to read about:
> * Comparison of combined assembly and polishing method [Trycycler: consensus long-read assemblies for bacterial genomes](https://link.springer.com/article/10.1186/s13059-021-02483-z)
> * Polishing strategy for ONT and Pacbio Hifi reads [Polishing high-quality genome assemblies](https://www.nature.com/articles/s41592-022-01515-1)
> * Comparison of polishing of ONT data with alignment free tool Jasper compared to POLCA,NextPolish and ntEdit [JASPER: a fast genome polishing tool that improves accuracy and creates population-specific reference genomes](https://www.biorxiv.org/content/10.1101/2022.06.14.496115v1.full)
> * Comparison of short read polishers including pilon to the polisher Polypolish [Polypolish: Short-read polishing of long-read bacterial genome assemblies](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009802#)
> * Pilon short read polisher paper [Pilon: An Integrated Tool for Comprehensive Microbial Variant Detection and Genome Assembly Improvement](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0112963)
> * Accuracy of polishers including medaka for nanopore data [Nanpore consensus quality](https://github.com/rrwick/August-2019-consensus-accuracy-update#racon)
> * Comparison of nanopore polishing tools [Comparative evaluation of Nanopore polishing tools for microbial genome assembly and polishing strategies for downstream analysis](https://www.nature.com/articles/s41598-021-00178-w)
{: .callout}  
