#!/bin/bash

export PICARD='/mnt/lab_data/montgomery/nicolerg/tools/picard.jar'
GENOME_FASTA=/mnt/lab_data/montgomery/nicolerg/motrpac/reference/ensembl/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa
GTF_FILE=/mnt/lab_data/montgomery/nicolerg/motrpac/reference/ensembl/Rattus_norvegicus.Rnor_6.0.94.gtf
INDEX_DIRECTORY=/mnt/lab_data/montgomery/nicolerg/motrpac/reference/STAR

# build genome index with STAR ##############################################################################

# gunzip /mnt/lab_data/montgomery/nicolerg/motrpac/reference/ensembl/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
# gunzip /mnt/lab_data/montgomery/nicolerg/motrpac/reference/ensembl/Rattus_norvegicus.Rnor_6.0.94.gtf.gz

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
			--quantMode TranscriptomeSAM \
			--genomeLoad LoadAndKeep \
			--limitBAMsortRAM 15000000000 

	done
done

# quantification ##################################################################################################

# first, build prepare the RSEM reference:

OUTPUT_DIR=/mnt/lab_data/montgomery/nicolerg/motrpac/reference/RSEM_reference/rn6

if [ ! -d "$OUTPUT_DIR" ]; then
	mkdir ${OUTPUT_DIR}
fi

taskset -c 12-23 /mnt/lab_data/montgomery/nicolerg/tools/RSEM-1.3.1/rsem-prepare-reference \
	--gtf $GTF_FILE \
	$GENOME_FASTA \
	$OUTPUT_DIR

# now quantification: 

for run in {1..3}; do 
	basedir=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN${run}/aligned
	outdir=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/quant
	for dir in `ls ${basedir}`; do 
		if [ ! -d "$basedir/$dir" ]; then
			continue
		fi
		bam_file=`ls ${basedir}/${dir} | grep "Aligned.toTranscriptome.out.bam"`
		bam=${basedir}/${dir}/${bam_file}

		taskset -c 12-23 /mnt/lab_data/montgomery/nicolerg/tools/RSEM-1.3.1/rsem-calculate-expression \
			-p 12 \
			--bam \
			--paired-end \
			--no-bam-output \
			--forward-prob 0.5 \
			--seed 12345 \
			$bam \
			/mnt/lab_data/montgomery/nicolerg/motrpac/reference/RSEM_reference/rn6 \
			${outdir}/$dir 
	done
done

# collect RNA-seq metrics ########################################################################################

flat=/mnt/lab_data/montgomery/nicolerg/motrpac/reference/ensembl/Rattus_norvegicus.Rnor_6.0.94.refFlat
bam_base=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/aligned_bam

for sample in `ls $bam_base | grep "sortedByCoord" | sed "s/_Aligned.*//" | sort | uniq`; do
	bam_file=${bam_base}/${sample}_Aligned.sortedByCoord.out.bam
	rna_metrics=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/metrics/${sample}_output.RNA_Metrics

	java -Xmx10g -jar $PICARD CollectRnaSeqMetrics \
		I=$bam_file \
		O=$rna_metrics \
		REF_FLAT=$flat \
		STRAND=FIRST_READ_TRANSCRIPTION_STRAND 

done

# count alignments by segment #######################################################################################

for sample in `ls $bam_base | grep "sortedByCoord" | sed "s/_Aligned.*//" | sort | uniq`; do 

	file=${sample}_Aligned.sortedByCoord.out.bam

	# get table of reads per chromosome
	samtools view ${indir}/${file} | cut -f3 | sort | uniq -c > ${indir}/tmp.${sample}.read.cnts
	# get total
	total=`sed "s/^[ \t]*//" ${indir}/tmp.${sample}.read.cnts | cut -d' ' -f1 | awk '{total += $0} END{print total}'`
	# get num chromosomal
	head -20 ${indir}/tmp.${sample}.read.cnts > ${indir}/tmp.${sample}.chrom
	tail -2 ${indir}/tmp.${sample}.read.cnts >> ${indir}/tmp.${sample}.chrom
	chrom=`sed "s/^[ \t]*//" ${indir}/tmp.${sample}.chrom | cut -d' ' -f1 | awk '{total += $0} END{print total}'`
	# get num mitochondrial
	mito=`grep -E " MT$" ${indir}/tmp.${sample}.read.cnts | sed "s/^[ \t]*//" | cut -d' ' -f1 | awk '{total += $0} END{print total}'`
	# get contigs
	contig=`tail -n +21 ${indir}/tmp.${sample}.read.cnts | head -n -3 | sed "s/^[ \t]*//" | cut -d' ' -f1 | awk '{total += $0} END{print total}'`
	echo "$sample	$total	$chrom	$mito	$contig" >> ${indir}/merged.regional.read.counts.txt

done

# calculate % contaminants #######################################################################################

RRNA_REF=/mnt/lab_data/montgomery/nicolerg/motrpac/rat_ref/ensembl/rn6_rRNA
GLOBIN_REF=/mnt/lab_data/montgomery/nicolerg/motrpac/reference/ensembl/rn6_hemoglobin
OUTDIR=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/QC

# first, build bowtie2 indices

refdir=/mnt/lab_data/montgomery/nicolerg/motrpac/rat_ref
bowtie2-build $refdir/Rattus_norvegicus.Rnor_6.0.94.rRNA.fa $refdir/rn6_rRNA
bowtie2-build $refdir/Rattus_norvegicus.Rnor_6.0.94.hemoglobin.fa $refdir/rn6_hemoglobin

# now map raw reads

for run in {1..3}; do
	indir=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN${run}
	for prefix in `ls ${indir} | grep "fastq.gz" | grep -v "Undetermined" | grep -v "_I1_" | sed "s/_R[0-9]_001\..*//" | uniq`; do

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
