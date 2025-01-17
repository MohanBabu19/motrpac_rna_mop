# MoTrPAC RNA-seq pipeline

**NOTE:** (31 Jan 19) This MOP is currently outdated. Refer to the most recent version of the MoTrPAC RNA-seq MOP, which implements the pipeline available here: https://github.com/yongchao/motrpac_rnaseq


**Contact:** Nicole Gay (nicolerg@stanford.edu), Brunilda Balliu (bballiu@stanford.edu), Laure Fresard (lfresard@stanford.edu)

Sourced from MoTrPAC RNA-seq MOP: https://docs.google.com/document/d/1oz8jZAY9Rq4uqenp-0RMkQBhMtjsMLKlxfjlmRExSQ0/edit?ts=5b04a52e#  
**Authors:**  

* Brunilda Balliu, Laure Fresard, Stephen Montgomery Lab @ Stanford  
* Surendra Dasari @ Mayo  
* Tuuli Lappalainen @ NYGC  
* Elena Zaslavsky, Yongchao Ge @ Mtn Sinai  

This repository contains code to run all components of the MoTrPAC RNA-seq pipeline. The following software is required:

* STAR 2.6.1b https://github.com/alexdobin/STAR  
* Python 2.7 https://www.python.org/  
* cutadapt 1.18 https://pypi.python.org/pypi/cutadapt  
* Picard tools 2.18.16 https://broadinstitute.github.io/picard/ and https://github.com/broadinstitute/picard  
* Samtools 1.3.1 http://www.htslib.org/   
* RSEM 1.3.1 https://github.com/deweylab/RSEM/ (download link)   
* MultiQC 1.6 http://multiqc.info/  
* bowtie2 2.3.4.3 http://bowtie-bio.sourceforge.net/bowtie2/index.shtml  
* FastQC 0.11.8 https://www.bioinformatics.babraham.ac.uk/projects/fastqc/. FastQC is a java-based tool. Please ensure that the correct version of java is installed.   

## Outline  

##### A. Pre-alignment sample processing and QC

* A.1 Pre-alignment QC with MultiQC and FASTQC
* A1.1 Run FASTQC on each pair of FASTQ files
* A.1.2 MultiQC
* A.2 Adapter trimming  

##### B. Genome reference-specific analyses 

* B.1 Build genome index with STAR
* B.2 Alignment
* B.3 Quantification
* B.3.1 Prepare RSEM reference
* B.3.2 Run quantification 
* B.4 Post-alignment QC
* B.4.1 RNA-seq metrics with Picard CollectRnaSeqMetrics
* B.4.2 Count alignments by segment 
* B.4.3 Calculate % contamininants (rRNA and globin)
* B.4.3.1 Build bowtie indices for rRNA and globin FASTA files 
* B.4.3.2 Align raw FASTQ files to rRNA and globin bowtie2 indices  

##### C. Compile important metrics

* C.1 Summarize alignment metrics from STAR with MultiQC 
* C.2 Mark PCR duplicates
* C.3 Compile important metrics into a single report 

##### D. Flag problematic samples 
##### E. Post-Quantification QC 

# A. Pre-alignment sample processing and QC

Section A is **not** genome reference-specific. The same steps can be run on both human and rat samples simultaneously.

## A.1 Pre-alignment QC with MultiQC and FASTQC

### A1.1 Run FASTQC on each pair of FASTQ files

This step will generate a single zipped folder and html file for each FASTQ file. The folder contains metrics about read quality. 

Parameters:  

* `FASTQ_DIR`: Input directory that contains FASTQ files. The following code assumes that there are gzipped, paired R1 and R2 files for each sample, and it excludes "Undetermined" FASTQ files generated by `bcl2fastq`. 
* `OUTPUT_DIR`: Output folder for FASTQC reports
* `ADAPTER_FILE`: A tab-delimited file that contains expected adapter sequences. See `data/adapter.tsv` for a usable example.  

The following code assumes that the paired FASTQ files are named with the format `ANY_SAMPLE_NAME_R1_001.fastq.gz` and `ANY_SAMPLE_NAME_R2_001.fastq.gz`.  

```bash
for prefix in `ls ${FASTQ_DIR} | grep "fastq.gz" |  grep -v "_I1_" | grep -v "Undetermined" | sed "s/_R[0-9]_001\..*//" | uniq`; do 
	r1=`ls ${FASTQ_DIR} | grep ${prefix} | grep "R1_001.fastq.gz"`
	r2=`ls ${FASTQ_DIR} | grep ${prefix} | grep "R2_001.fastq.gz"`
	fastqc -o ${OUTPUT_DIR} -t 12 -f fastq -a $ADAPTER_FILE ${FASTQ_DIR}/${r1} ${FASTQ_DIR}/${r2}
done
```

### A.1.2 MultiQC

MultiQC compiles all log files from analysis and processing tools and combines them into a single report. Run MultiQC on all directories containing outputs from FastQC:

Parameters:  

* `OUTPUT_DIR`: Output folder for MultiQC report
* `FASTQC_DIR_N`: Directory/directories containing FASTQC reports from the previous step. Include as many directories as needed (minimum 1).

Options:  

* `-d`: Prepend the directory to file names (useful for keeping track of samples in multiple runs)  
* `-f`: Overwrite existing MultiQC reports in `OUTPUT_DIR`  

```bash
multiqc \
	-d \
	-f \
	-n $REPORT_NAME \
	-o $OUTPUT_DIR \
	$FASTQC_DIR_1 \
	$FASTQC_DIR_2 
```

## A.2 Adapter trimming

For each paired-end FASTQ file (`${SAMPLE}\_R1_001.fastq.gz` and `${SAMPLE}\_R2_001.fastq.gz`), remove adapters (set using INDEXED_ADAPTER_PREFIX and UNIVERSAL_ADAPTER) using cutadapt v1.16. Eliminate reads that are too short after removing adapters and save trimmed FASTQ files (`${SAMPLE}\_R1_001.trimmed.fastq.gz` and `${SAMPLE}\_R2_001.trimmed.fastq.gz`). Save reads that were trimmed because they were too short in `${prefix}_R[1-2]_001.tooshort.fastq.gz`.  

Parameters:  

* `FASTQ_DIR`: Input directory that contains FASTQ files. The following code assumes that it contains gzipped, paired R1 and R2 files for each sample.  
* `OUTPUT_DIR`: Output directory for trimmed FASTQ files  
* `m`: Remove reads that are this short after removing adapters, e.g. m=50 bp for 150 bp reads and m=20 bp for 75 bp reads

```bash
# for NUGEN kits (Illumina TruSeq adapters):
INDEXED_ADAPTER_PREFIX=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
UNIVERSAL_ADAPTER=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

for prefix in `ls ${FASTQ_DIR} | grep "fastq.gz" |  grep -v "_I1_" | sed "s/_R[0-9]_001\..*//" | uniq`; do 
	r1=`ls ${indir} | grep ${prefix} | grep "R1_001.fastq.gz"`
	r2=`ls ${indir} | grep ${prefix} | grep "R2_001.fastq.gz"`

	cutadapt \
		-a $INDEXED_ADAPTER_PREFIX \
		-A $UNIVERSAL_ADAPTER \
		-o ${OUTPUT_DIR}/${prefix}_R1_001.trimmed.fastq.gz -p ${OUTPUT_DIR}/${prefix}_R2_001.trimmed.fastq.gz \
		-m 20 \
		--too-short-output ${OUTPUT_DIR}/${prefix}_R1_001.tooshort.fastq.gz \
		--too-short-paired-output ${OUTPUT_DIR}/${prefix}_R2_001.tooshort.fastq.gz \
 		${FASTQ_DIR}/${r1} ${FASTQ_DIR}/${r2} 
done
```

# B. Genome reference-specific analyses

Section B requires different inputs for rat and human samples.  

## B.1 Build genome index with STAR

For rat (Rnor_6.0, Release 94), use the following files:  

* Genome build (`GENOME_FASTA`): ftp://ftp.ensembl.org/pub/release-94/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz     
* GTF (`GTF_FILE`): ftp://ftp.ensembl.org/pub/release-94/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.94.gtf.gz  

For human (GRCh38.p12, Release 29), use the following files:  

* Genome build (`GENOME_FASTA`): ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.p12.genome.fa.gz  
* GTF (`GTF_FILE`): ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz  

Parameters:  

* `NUM_THREADS`: Number of threads to run on
* `INDEX_DIRECTORY`: Output folder for the genome index
* `OVERHANG`: splice-junction-data-base-overhang (sjdbOverhang) parameter for STAR. Should have a value of `read length – 1`, e.g. for read length 75 it should be set to 74.  

Note that `GENOME_FASTA` and `GTF_FILE` cannot be in `.gz` format for compatibility with STAR (i.e. gunzip the files before running the following command).

```bash
NUM_THREADS=10
OVERHANG=74 # for 2x75 NextSeq and MiSeq runs

STAR \
	--runThreadN $NUM_THREADS \
	--runMode genomeGenerate \
	--genomeDir $INDEX_DIRECTORY \
	--genomeFastaFiles $GENOME_FASTA \
	--sjdbGTFfile $GTF_FILE \
	--sjdbOverhang $OVERHANG
```

## B.2 Alignment 

Align each trimmed FASTQ file to the reference using STAR.  

Parameters:  

* `GTF_FILE`: Same file as in step B.1
* `INDEX_DIRECTORY`: Path to genome index from step B.1  
* `TRIMMED_FASTQ_DIR`: Folder containing trimmed FASTQ files (step A.2) 
* `OUTPUT_DIR`: Directory for output of STAR alignment, including sorted BAM files

```bash
for prefix in `ls ${TRIMMED_FASTQ_DIR} | grep "fastq.gz" | sed "s/_R[0-9]_001\..*//" | uniq`; do 

	r1=`ls ${TRIMMED_FASTQ_DIR} | grep ${prefix} | grep "R1_001.trimmed.fastq.gz"`
	r2=`ls ${TRIMMED_FASTQ_DIR} | grep ${prefix} | grep "R2_001.trimmed.fastq.gz"`

	STAR \
		--genomeDir $INDEX_DIRECTORY \
		--readFilesIn ${TRIMMED_FASTQ_DIR}/${r1} ${TRIMMED_FASTQ_DIR}/${r2} \
		--outFileNamePrefix ${OUTPUT_DIR}/${prefix}_ \
		--readFilesCommand zcat \
		--outSAMattributes NH HI AS NM MD nM\
		--outFilterType BySJout \
		--runThreadN 8 \
		--outSAMtype BAM SortedByCoordinate \
		--quantMode TranscriptomeSAM\
		--genomeLoad LoadAndKeep \
		--limitBAMsortRAM 15000000000 
done
```

## B.3 Quantification

Use RSEM v. 1.3.1 for quantification to obtain gene-level and transcript-level quantification

### B.3.1 Prepare RSEM reference

First, prepare the RSEM reference. This only needs to be done once.  

Parameters:  

* `OUTPUT_PREFIX`: Output prefix for RSEM reference, including output directory. i.e. `OUTPUT_PREFIX=${OUTPUT_DIR}/${PREFIX}`  
* `GENOME_FASTA` and `GTF_FILE` are the same as in step B.1  

```bash 
$PATH_TO_RSEM/rsem-prepare-reference \
	--gtf $GTF_FILE \
	$GENOME_FASTA \
	$OUTPUT_PREFIX
```

### B.3.2 Run quantification 

Now run quantification with the sorted `*Aligned.toTranscriptome.sorted.bam` files.  

Parameters:  

* `ALIGNED_IN`: Input directory of aligned BAM files  
* `RSEM_REFERENCE`: Same as `OUTPUT_PREFIX` in the previous step (B.3.1)  
* `OUTPUT_DIR`: Output directory for quantification stats  

```bash 
for sample in `ls ${ALIGNED_IN} | grep "Aligned.toTranscriptome.out.bam" | sed "s/_Aligned.*//" | sort | uniq`; do

	bam_file=${ALIGNED_IN}/${sample}_Aligned.toTranscriptome.out.bam

	$PATH_TO_RSEM/rsem-calculate-expression \
		-p 12 \
		--bam \
		--paired-end \
		--no-bam-output \
		--forward-prob 0.5 \
		--seed 12345 \
		${bam_file} \
		${RSEM_REFERENCE} \
		${OUTPUT_DIR}/${sample}

	done
done
```

## B.4 Post-alignment QC

### B.4.1 RNA-seq metrics with Picard CollectRnaSeqMetrics

Compute post alignment QC metrics, including % mapped to coding, intron, intergenic, UTR, and % correct strand, and 5’ to 3’ bias.  

refFlat files are required for this step, which can be generated from GTF files. Human (gencode.v29.annotation.refFlat) and rat (Rattus_norvegicus.Rnor_6.0.94.refFlat) refFlat files are provided in the `data` directory of this repostitory, generated from the GTF files specified in step B.1. See the README file in the `data` directory for more information on how these files were generated.

Parameters:  

* `REF_FLAT`: Gene annotations in refFlat format  
* `BAM_DIR`: Path to BAM files output by STAR in step B.2 
* `OUTPUT_DIR`: Output directory for RNA-seq metrics   

```bash
for sample in `ls $BAM_DIR | grep "Aligned.sortedByCoord.out.bam" | sed "s/_Aligned.*//" | sort | uniq`; do

	bam=${BAM_DIR}/${sample}_Aligned.sortedByCoord.out.bam

	java -Xmx10g -jar picard.jar CollectRnaSeqMetrics \
		I=$bam  \
		O=$OUTPUT_DIR/${sample}_output.RNA_Metrics \
		REF_FLAT=$REF_FLAT \
		STRAND=FIRST_READ_TRANSCRIPTION_STRAND 
done
```

### B.4.2 Count alignments by segment 

Calculate the number of reads mapped to each genomic segment (chromosomes and contigs).  

This section outputs a single file in `${INDIR}` called `${sample}.idxstats.log` with the following columns:  

* `SAMPLE`: Prefix of ${PREFIX}\_Aligned.sortedByCoord.out.bam
* `TOTAL_READS`: Total mapped reads
* `CHROM_READS`: Number of reads mapped to chromosomes (autosomal and sex)
* `MITO_READS`: Number of reads mapped to the mitochondrial chromosome
* `CONTIG_READS`: Number of reads mapped to contigs 

Parameters:  

* `INDIR`: Directory containing ${SAMPLE}\_Aligned.sortedByChrom.out.bam files from step B.2  

```bash
INDIR=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/aligned_bam
mkdir -p ${INDIR}/idxstats
for bam in `ls ${INDIR} | grep -E "Aligned.sortedByCoord.out.bam$"`; do 
	# index 
	samtools index ${INDIR}/$bam
	# get summary
	echo "seq_name	seq_length	n_mapped_reads	n_unmapped_reads" > ${INDIR}/idxstats/${sample}.idxstats.log 
	samtools idxstats ${INDIR}/$bam >> ${INDIR}/idxstats/${sample}.idxstats.log
done
first=1
# summarize in one file 
for log in `ls ${INDIR}/idxstats | grep "idxstats.log"`; do 
	if [ "$first" == 1 ]; then 
		cut -f 1,2 ${INDIR}/idxstats/${log} > ${INDIR}/idxstats/tmp1
		first=0
	fi
	sample_prefix=`echo ${log} | sed "s/\.idxstats\.log//"`
	cut -f 3 ${INDIR}/idxstats/${log} | sed "1 s/^.*$/$sample_prefix/" > ${INDIR}/idxstats/tmp2
	paste ${INDIR}/idxstats/tmp1 ${INDIR}/idxstats/tmp2 > ${INDIR}/idxstats/tmp3
	rm ${INDIR}/idxstats/tmp1 ${INDIR}/idxstats/tmp2
	mv ${INDIR}/idxstats/tmp3 ${INDIR}/idxstats/tmp1
done
mv ${INDIR}/idxstats/tmp1 ${INDIR}/idxstats/idxstats.merged.log 
```

### B.4.3 Calculate % contamininants (rRNA and globin)

Map the original FASTQ files against the species-specific ribosomal RNA and hemoglobin sequences to get the number of reads that map to rRNA and globin. 

See the `*.fa` files in the `data` directory of this repository for rat and human rRNA and globin FASTA files. The README in the `data` directory details how these files were generated.  

#### B.4.3.1 Build bowtie indices for rRNA and globin FASTA files  

Run the following command to build bowtie2 indices for each FASTA file, where `INDEX_PREFIX` includes the path to the output directory plus a prefix for the index (e.g. `/PATH/TO/OUTDIR/rn6_globin`):  

```bash
bowtie2-build $FASTA_FILE $INDEX_PREFIX 
``` 
#### B.4.3.2 Align raw FASTQ files to rRNA and globin bowtie2 indices  

Once the indices are built for both globin and rRNA, run the following code to map each BAM file to the contaminant FASTA files. This step will generate a SAM file and an alignment stat file for each sample. The *merge results* section of the code produces `merged_rRNA.log` and `merged_globin.log` summary files to provide % reads mapped to each sample.  

The following codeassumes that the paired FASTQ files are named with the format `ANY_SAMPLE_NAME_R1_001.fastq.gz` and `ANY_SAMPLE_NAME_R2_001.fastq.gz`. 

Parameters:  

* `RRNA_REF`: Path to rRNA bowtie2 index (same as `$INDEX_PREFIX` above)  
* `GLOBIN_REF`: Path to globin bowtie2 index (same as `$INDEX_PREFIX` above)  
* `OUTDIR`: Output directory for SAM file and alignment stats  
* `INDIR`: Directory with raw FASTQ files  

```bash
for prefix in `ls ${INDIR} | grep "fastq.gz" | grep -v "Undetermined" | grep -v "_I1_" | sed "s/_R[0-9]_001\..*//" | uniq`; do
	r1=${INDIR}/${prefix}_R1_001.fastq.gz
	r2=${INDIR}/${prefix}_R2_001.fastq.gz

	bowtie2 --local -p 20 -s --sam-nohead -x $RRNA_REF \
		-1 $r1 -2 $r2 >${OUTDIR}/${prefix}_rRNA.sam  2>${OUTDIR}/${prefix}_rRNA.log

	bowtie2 --local -p 20 -s --sam-nohead -x $GLOBIN_REF \
		-1 $r1 -2 $r2 >${OUTDIR}/${prefix}_globin.sam  2>${OUTDIR}/${prefix}_globin.log
done

# merge results

for file in `ls $OUTDIR | grep "rRNA.log"`; do 
	prefix=`echo $file | sed "s/_rRNA.*//"`
	alignment=`tail -1 $file`; echo "$prefix	$alignment" >> $OUTDIR/merged_rRNA.log
done

for file in `ls $OUTDIR | grep "globin.log"`; do 
	prefix=`echo $file | sed "s/_globin.*//"`
	alignment=`tail -1 $file`; echo "$prefix	$alignment" >> $OUTDIR/merged_globin.log
done
```

# C. Compile important metrics

The metrics in Section C can be determined the same way for both human and rat samples from outputs from the preceding steps.

## C.1 Summarize alignment metrics from STAR with MultiQC

Compile the metrics from the Log.final.out files into a single MultiQC report.  

These metrics include:  

* Number of total reads  
* Average read length  
* Average mapping length  
* Number (or percentage) of uniquely mapped reads  
* Multi-mapped reads( too many loci, multiple loci)  
* Unmapped reads (too many mismatches, too short, other)  
* Chimeric reads  

Parameters:  

* `OUTPUT_DIR`: Output folder for MultiQC report
* `FASTQC_DIR_N`: Directory/directories containing STAR log files from alignment. Include as many directories as needed (minimum 1).

Options:  

* `-d`: Prepend the directory to file names (useful for keeping track of samples in multiple runs)  
* `-f`: Overwrite existing MultiQC reports in `OUTPUT_DIR`  

```bash
multiqc \
	-d \
	-f \
	-n $REPORT_NAME \
	-o $OUTPUT_DIR \
	$FASTQC_DIR_1 \
	$FASTQC_DIR_2 
```
## C.2 Mark PCR duplicates

For each bam file, mark PCR duplicates using Picard tools. 

Parameters:  

* `BAM_DIR`: Directory containing ${SAMPLE}\_Aligned.sortedByChrom.out.bam files from step B.2  
* `OUTDIR`: Output directory for metrics file  

```bash 
for sample in `ls $BAM_DIR | grep "Aligned.sortedByCoord.out.bam" | sed "s/_Aligned.*//" | sort | uniq`; do

	bam=${BAM_DIR}/${sample}_Aligned.sortedByCoord.out.bam

	java -Xmx10g -jar picard.jar MarkDuplicates \
		INPUT=$bam \
		OUTPUT=${BAM_DIR}/${sample}.markedDup \
		CREATE_INDEX=true \
		VALIDATION_STRINGENCY=SILENT \
		METRICS_FILE=${OUTDIR}/${sample}.markedDup.out \
		REMOVE_DUPLICATES=false 
done 
```

## C.3 Compile important metrics into a single report 

### C.3.1 FASTQ metrics (raw and trimmed)

* PCT_GC: FASTQC metric from raw FASTQ files (take the average of the two values from the read pair)
* READ_CNT_RAW: number of reads in raw FASTQs for each sample (i.e. the line count of SAMPLE_R1_001.fastq divided by 2)
* READ_CNT_TRIMMED: number of reads in trimmed FASTQs for each sample (i.e. the line count of SAMPLE_R1_001.trimmed.fastq divided by 2)

### C.3.2 Contamination metrics (step B.4.4.2)

* PCT_RRNA: from merged_rRNA.log
* PCT_GLOBIN: from merged_globin.log

### C.3.3 Post-alignment RNA-seq metrics from Picard's CollectRnaSeqMetrics output files (step B.4.1)

* PCT_EXONIC: (PCT_CODING_BASES + PCT_UTR_BASES)\*100 
* PCT_INTRONIC: PCT_INTRONIC_BASES 
* PCT_INTERGENIC: PCT_INTERGENIC_BASES
* PCT_CORRECT_STRAND: PCT_CORRECT_STRAND_READS
* 5_3_BIAS: MEDIAN_5PRIME_TO_3PRIME_BIAS

### C.3.4 Alignment metrics from STAR's ${SAMPLE}\_Log.final.out files (step B.2)

* READ_CNT_PAIRED_INPUT: `Number of input reads` 
* READ_CNT_MAPPED_UNIQ: `Uniquely mapped reads number`
* PCT_MAPPED_UNIQ: `Uniquely mapped reads %`
* READ_CNT_MAPPED_MULTI: `Number of reads mapped to multiple loci`
* PCT_MAPPED_MULTI: `% of reads mapped to multiple loci`
* READ_COUNT_TOO_MANY_MULTI: `Number of reads mapped to too many loci`
* PCT_TOO_MANY_MULTI: `% of reads mapped to too many loci`

### C.3.5 Alignment metrics from RSEM ${SAMPLE}.stat/${SAMPLE}.cnt files (step B.3.2)

* N_ALIGNABLE: line 1, col 2 (space-delimited)
* N_UNIQUE: line 2, col 1 (space-delimited)
* N_MULTI: line 2, col 2 (space-delimited)
* N_UNCERTAIN: line 2, col 3 (space-delimited)
* N_TOTAL_ALIGNMENTS: line 3, col 1 (space-delimited)

### C.3.6 Alignment metrics from merged.regional.read.counts.txt (step B.4.3)

* TOTAL_MAPPED: TOTAL_READS
* MAPPED_CHROM: CHROM_READS
* MAPPED_MT: MITO_READS
* MAPPED_CONTIG: CONTIG_READS 

### C.3.7 PCR duplicate metrics from ${SAMPLE}.markedDup.out (step C.2)

* TBD 

# D. Flag problematic samples

*These metrics are very liberal to remove samples of very low quality. Additional samples might be removed downstream using more stringent criteria.*

Flag samples with 

1. Low RIN scores, e.g. RIN < 6
2. Low number of sequenced reads after trimming adaptors, e.g. samples with < 20M reads
3. Abnormal GC content, e.g. GC% > 80% or < 20%
4. High percentage of rRNA, e.g. rRNA > 20%
5. Low number of mapped reads, e.g. < 60%
6. Low number of % exonic read, % exonic read < 50%
7. Mapped read count < 50% of average mapped read count per sample  

# E. Post-Quantification QC 

*These steps should be updated when new samples become available because they are depend on all samples.*

1. Confirm sample concordance with other omics genotypes from same individual using VerifyBAMID.
2. Obtain filtered and normalized expression data and filter outliers:  

  * Remove genes with 0 counts in every sample. 
  * Filter lowly expressed genes (for QC purposes only, definition of expressed gene might change for each analyses): genes with mean counts < 5 and zero counts in more than 20% of samples.
  * Perform library size correction and variance stabilization using DESeq2 or other methods.
  * Compute D-statistic for each sample based on filtered and normalized expression data, i.e. average Spearman’s correlation with other samples from the same tissue and time point. Flag samples with D<.80.
  * Perform PCA based on filtered and normalized expression data across and within tissues. Mark outlier samples (outside +/- 3 sd).  
