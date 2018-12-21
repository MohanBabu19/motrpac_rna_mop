#!/bin/bash 

# QC on raw FASTQ files #############################################################################################

# FASTQC

ADAPTER_FILE=adapter.tsv

for run in RUN1 RUN2 RUN3; do

	indir=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_${run}
	outdir=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_${run}/fastqc

	if [ ! -e "$outdir" ]; then
		mkdir ${outdir}
	fi

	for prefix in `ls ${indir} | grep "fastq.gz" |  grep -v "_I1_" | grep -v "Undetermined" | sed "s/_R[0-9]_001\..*//" | uniq`; do 
		r1=`ls ${indir} | grep ${prefix} | grep "R1_001.fastq.gz"`
		r2=`ls ${indir} | grep ${prefix} | grep "R2_001.fastq.gz"`
		taskset -c 0-11 fastqc -o ${outdir} -t 12 -f fastq -a ${ADAPTER_FILE} ${indir}/${r1} ${indir}/${r2} 
	done
done

wait

# MultiQC

taskset -c 0-11 multiqc \
	-d \
	-f \
	-n motrpac.rna.raw.fastq.multiqc \
	-o /mnt/lab_data/montgomery/nicolerg/motrpac/rna/multiqc \
	/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN1/fastqc \
	/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN2/fastqc \
	/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN3/fastqc

# adapter trimming #############################################################################################

INDEXED_ADAPTER_PREFIX=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
UNIVERSAL_ADAPTER=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

for run in {1..3}; do 

	indir=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN${run}
	outdir=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN${run}/cutadapt

	if [ ! -e "$outdir" ]; then
		mkdir ${outdir}
	fi

	for prefix in `ls ${indir} | grep "fastq.gz" |  grep -v "_I1_" | sed "s/_R[0-9]_001\..*//" | uniq`; do 
		r1=`ls ${indir} | grep ${prefix} | grep "R1_001.fastq.gz"`
		r2=`ls ${indir} | grep ${prefix} | grep "R2_001.fastq.gz"`

		taskset -c 0-23 cutadapt \
			-a $INDEXED_ADAPTER_PREFIX \
			-A $UNIVERSAL_ADAPTER \
			-o ${outdir}/${prefix}_R1_001.trimmed.fastq.gz -p ${outdir}/${prefix}_R2_001.trimmed.fastq.gz \
			-m 20 \
			--too-short-output ${outdir}/${prefix}_R1_001.tooshort.fastq.gz \
			--too-short-paired-output ${outdir}/${prefix}_R2_001.tooshort.fastq.gz \
			${indir}/${r1} ${indir}/${r2} 
	done
done
