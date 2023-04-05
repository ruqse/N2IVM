

# RNA-Seq Nextflow Pipeline for Ivermectin-Exposed C. elegans (N2)

## Creating Decoy Sequence headers and Gentrome for C. elegans RNA-Seq Analysis



```
# Create decoy file
# The C. elegans genome file: c_elegans.PRJNA13758.WS283.genomic.fa.gz
# The C. elegans transcriptome file: c_elegans.PRJNA13758.WS283.all_transcript.fa.gz

grep "^>" <(gunzip -c c_elegans.PRJNA13758.WS283.genomic.fa.gz) | cut -d " " -f 1 > decoys.txt

sed -i.bak -e 's/>//g' decoys.txt
```


Concatenate the C. elegans transcriptome and genome files to create the gentrome file:

```
cat c_elegans.PRJNA13758.WS283.all_transcript.fa.gz c_elegans.PRJNA13758.WS283.genomic.fa.gz > gentrome.fa.gz

```

## Pipeline processes raw sequencing data, performs quality control, and quantifies transcript abundance


*Input parameters*

params.reads: Raw sequencing reads in FASTQ format. 

params.decoys: File containing decoy sequences headers.

params.gentrome: Concatenated genome and transcriptome fasta file.

params.adapters: Adapter sequences for trimming. 

params.sortmerna_ref: SortMeRNA reference database. 

*Output directories*

params.multiqc: Directory for the MultiQC report.

params.outdir: Base output directory. 

params.trimmooutdir: Subdirectory for Trimmomatic results. 

params.salmonoutdir: Subdirectory for Salmon results.

params.fastqcdir: Subdirectory for FastQC results. 

params.sortmernaoutdir: Subdirectory for SortMeRNA results. 

Processes
1. SORTMERNA: Removes residual ribosomal RNA using SortMeRNA with a curated reference database.

2. TRIMMOMATIC: Trims poor-quality read ends and removes adapter sequences using Trimmomatic.

3. SALMON_INDEX: Creates a binary index of the transcriptome file using Salmon.

4. SALMON_QUANTIFICATION: Performs quantification of transcript abundance using Salmon.

5. FASTQC: Assesses the quality of raw reads, rRNA-removed reads, and trimmed reads using FastQC.

6. MULTIQC: Aggregates results from Salmon and FastQC in a single report using MultiQC.

**Workflow**

The pipeline processes input files in parallel and aggregates results using the specified output directories. The final output includes the MultiQC report, which summarizes quality control and quantification results.

Once you have completed quantifying the salmon, you can move forward with downstream analysis using differential expression tools such as DESeq2 or edgeR. With the help of the tximport package, you can easily import transcript-level quantifications from salmon and choose to aggregate them to the gene level for gene-level differential expression analysis.


