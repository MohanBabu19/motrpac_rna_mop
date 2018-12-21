#!/bin/bash

export PICARD='/mnt/lab_data/montgomery/nicolerg/tools/picard.jar'
GENOME_FASTA=/mnt/lab_data/montgomery/nicolerg/motrpac/hsapiens_ref/GRCh38.p12.genome.fa
GTF_FILE=/mnt/lab_data/montgomery/nicolerg/motrpac/hsapiens_ref/gencode.v29.annotation.gtf
INDEX_DIRECTORY=/mnt/lab_data/montgomery/nicolerg/motrpac/hsapiens_ref/STAR

hs_fastq_prefix=( /mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN1/PAX2_thaw_1x_S1 \
	/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN1/PAX2_thaw_2x_S3 \
	/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN1/PAX3_thaw_1x_S2 \
	/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN1/PAX3_thaw_2x_S4 \
	/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN2/RD165_S7 \
	/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN2/RD166_S8 \
	/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN2/RD167_S9 )

hs_sample=( PAX2_thaw_1x_S1 \
	PAX2_thaw_2x_S3 \
	PAX3_thaw_1x_S2 \
	PAX3_thaw_2x_S4 \
	RD165_S7 \
	RD166_S8 \
	RD167_S9 )

# build genome index with STAR ##############################################################################

if [ ! -e "$INDEX_DIRECTORY" ]; then
	mkdir $INDEX_DIRECTORY
fi

/usr/local/bin/STAR-2.5.2b \
	--runThreadN 10 \
	--runMode genomeGenerate \
	--genomeDir $INDEX_DIRECTORY \
	--genomeFastaFiles $GENOME_FASTA \
	--sjdbGTFfile $GTF_FILE \
	--sjdbOverhang 74

# alignment ##################################################################################################

for run in {1..3}; do 

	indir=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN${run}/cutadapt
	outdir=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN${run}/aligned

	if [ ! -e "$outdir" ]; then
		mkdir ${outdir}
	fi

	for prefix in `ls ${indir} | grep "fastq.gz" | sed "s/_R[0-9]_001\..*//" | uniq`; do 

		r1=`ls ${indir} | grep ${prefix} | grep "R1_001.trimmed.fastq.gz"`
		r2=`ls ${indir} | grep ${prefix} | grep "R2_001.trimmed.fastq.gz"`

		out=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN${run}/aligned/${prefix}

		if [ ! -e "$out" ]; then
			mkdir ${out}
		fi

		taskset -c 24-35 /usr/local/bin/STAR-2.5.2b \
			--genomeDir $INDEX_DIRECTORY \
			--readFilesIn ${indir}/${r1} ${indir}/${r2} \
			--outFileNamePrefix ${out}/${prefix}_ \
			--readFilesCommand zcat \
			--outSAMattributes NH HI AS NM MD nM\
			--outFilterType BySJout \
			--runThreadN 8 \
			--outSAMtype BAM SortedByCoordinate \
			--quantMode TranscriptomeSAM\
			--genomeLoad LoadAndKeep \
			--limitBAMsortRAM 15000000000 

	done
done

# quantification ##################################################################################################

# first, build prepare the RSEM reference:

OUTPUT=/mnt/lab_data/montgomery/nicolerg/motrpac/hsapiens_ref/RSEM_reference/hsapiens

taskset -c 12-23 /mnt/lab_data/montgomery/nicolerg/tools/RSEM-1.3.1/rsem-prepare-reference \
	--gtf $GTF_FILE \
	$GENOME_FASTA \
	$OUTPUT

outdir=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/quant

# now quantification: 

for prefix in "${hs_fastq_prefix[@]}"; do
	pre=`echo $prefix | sed "s/.*\///"`
	base=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/aligned_bam
	bam_file=${base}/${pre}_Aligned.toTranscriptome.out.bam

	taskset -c 12-23 /mnt/lab_data/montgomery/nicolerg/tools/RSEM-1.3.1/rsem-calculate-expression \
		-p 12 \
		--bam \
		--paired-end \
		--no-bam-output \
		--forward-prob 0.5 \
		--seed 12345 \
		$bam_file \
		/mnt/lab_data/montgomery/nicolerg/motrpac/hsapiens_ref/RSEM_reference/hsapiens \
		${outdir}/${pre} &

	running=`ps -ef | grep "nicolerg" | grep "RSEM-1.3.1" | wc -l`
	while [ $running -gt 4 ]; do
		sleep 30
		running=`ps -ef | grep "nicolerg" | grep "RSEM-1.3.1" | wc -l`
	done

done

# collect RNA-seq metrics ########################################################################################

refdir=/mnt/lab_data/montgomery/nicolerg/motrpac/hsapiens_ref
flat=$refdir/gencode.v29.annotation.refFlat
bam_base=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/aligned_bam

for file in "${hs_fastq_prefix[@]}"; do
	prefix=`echo $file | sed "s/.*\///"`
	bam_file=${bam_base}/${prefix}_Aligned.sortedByCoord.out.bam
	rna_metrics=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/metrics/${prefix}_output.RNA_Metrics

	java -Xmx10g -jar $PICARD CollectRnaSeqMetrics \
		I=$bam_file \
		O=$rna_metrics \
		REF_FLAT=$flat \
		STRAND=FIRST_READ_TRANSCRIPTION_STRAND 
done

# count alignments by segment #######################################################################################

for sample in "${hs_sample[@]}"; do 

	file=${sample}_Aligned.sortedByCoord.out.bam
	# get table of reads per chromosome
	samtools view ${indir}/${file} | cut -f3 | sort | uniq -c > ${indir}/tmp.${sample}.read.cnts
	# get total
	total=`sed "s/^[ \t]*//" ${indir}/tmp.${sample}.read.cnts | cut -d' ' -f1 | awk '{total += $0} END{print total}'`
	# get num chromosomal
	chrom=`grep "chr[1-9XY]" ${indir}/tmp.${sample}.read.cnts | sed "s/^[ \t]*//" | cut -d' ' -f1 | awk '{total += $0} END{print total}'`
	# get num mitochondrial
	mito=`grep "chrM" ${indir}/tmp.${sample}.read.cnts | sed "s/^[ \t]*//" | cut -d' ' -f1 | awk '{total += $0} END{print total}'`
	# get contigs
	contig=`grep -v "chr" ${indir}/tmp.${sample}.read.cnts | sed "s/^[ \t]*//" | cut -d' ' -f1 | awk '{total += $0} END{print total}'`

	echo "$sample	$total	$chrom	$mito	$contig" >> ${indir}/merged.regional.read.counts.txt

done

# calculate % contaminants #######################################################################################

RRNA_REF=/mnt/lab_data/montgomery/nicolerg/motrpac/hsapiens_ref/hsapiens_rRNA
GLOBIN_REF=/mnt/lab_data/montgomery/nicolerg/motrpac/hsapiens_ref/hsapiens_hemoglobin
OUTDIR=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/QC

# first, build bowtie2 indices

refdir=/mnt/lab_data/montgomery/nicolerg/motrpac/hsapiens_ref
bowtie2-build $refdir/gencode.v29.annotation.rRNA.fa $refdir/hsapiens_rRNA
bowtie2-build $refdir/gencode.v29.annotation.hemoglobin.fa $refdir/hsapiens_hemoglobin

# now map raw reads

for run in {1..3}; do
	indir=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN${run}
	for prefix in `ls ${indir} | grep "fastq.gz" | grep -v "Undetermined" | grep -v "_I1_" | grep -v "PAX" | sed "s/_R[0-9]_001\..*//" | uniq`; do
		
		r1=${indir}/${prefix}_R1_001.fastq.gz
		r2=${indir}/${prefix}_R2_001.fastq.gz

		bowtie2 --local -p 20 -s --sam-nohead -x $RRNA_REF \
			-1 $r1 -2 $r2 >${OUTDIR}/${prefix}_rRNA.sam  2>${OUTDIR}/${prefix}_rRNA.log

		bowtie2 --local -p 20 -s --sam-nohead -x $GLOBIN_REF \
			-1 $r1 -2 $r2 >${OUTDIR}/${prefix}_globin.sam  2>${OUTDIR}/${prefix}_globin.log
	done
done

# merge results

for file in `ls $OUTDIR | grep "rRNA.log"`; do 
	prefix=`echo $file | sed "s/_rRNA.*//"`
	alignment=`tail -1 $file`; echo "$prefix	$alignment" >> $OUTDIR/merged_rRNA.log
done

for file in `ls $OUTDIR | grep "globin.log"`; do 
	prefix=`echo $file | sed "s/_globin.*//"`
	alignment=`tail -1 $file`; echo "$prefix	$alignment" >> merged_globin.log
done
