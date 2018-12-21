#!/bin/bash 

# collect alignment metrics from STAR (MultiQC) ########################################################################

# Collect the metrics from the Log.final.out from star alignment
# This includes the number of total reads,
# 	avg read length, 
# 	av mapping length, 
# 	number (or percentage) of uniquely mapped reads, 
# 	multi-mapped reads( too many loci, multiple loci), 
# 	unmapped reads (too many mismatches, too short, other) 
# 	and chimeric reads

# # make master file with outputs from STAR alignment (OPTIONAL)

# for run in {1..3}; do
# 	first=1
# 	basedir=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN${run}/aligned
# 	for dir in `ls ${basedir}`; do 
# 		if [ ! -d "$basedir/$dir" ]; then
# 			continue
# 		fi
# 		if [ "$first" -eq "1" ]; then
# 			echo "STAT" > ${basedir}/STAR.merged.rows
# 			cut -f1 ${basedir}/${dir}/${dir}_Log.final.out | sed "s/^[ \t]*//" | sed "s/ |//" >> ${basedir}/STAR.merged.rows
# 			first=0
# 		fi
# 		echo "$dir" > ${basedir}/${dir}.tmp
# 		cut -f2 ${basedir}/${dir}/${dir}_Log.final.out | sed "s/^[ \t]*//" >> ${basedir}/${dir}.tmp
# 	done
# 	paste ${basedir}/STAR.merged.rows ${basedir}/*.tmp > ${basedir}/STAR.merged.log.out 
# 	rm ${basedir}/STAR.merged.rows ${basedir}/*.tmp
# 	sed -i "/UNIQUE READS/d" ${basedir}/STAR.merged.log.out 
# 	sed -i "/MULTI-MAPPING READS/d" ${basedir}/STAR.merged.log.out 
# 	sed -i "/UNMAPPED READS/d" ${basedir}/STAR.merged.log.out 
# 	sed -i "/CHIMERIC READS/d" ${basedir}/STAR.merged.log.out 
# 	sed -i "/^\s*$/d" ${basedir}/STAR.merged.log.out 
# done

# run it through MultiQC

taskset -c 0-11 multiqc \
	-d \
	-f \
	--ignore STAR.merged* \
	-n motrpac.rna.postalignment.report \
	-o /mnt/lab_data/montgomery/nicolerg/motrpac/rna/multiqc \
	/mnt/lab_data/montgomery/nicolerg/motrpac/rna/aligned 

# mark PCR duplicates #####################################################################################

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

# compile metrics into a single report ##########################################################################

## the following code just prints values to stdout, not in a particularly organized way
## this is mainly to show how to systematically access values of interest

# get % exonic reads from RNA metrics (should be >50%)
for file in `ls /mnt/lab_data/montgomery/nicolerg/motrpac/rna/metrics`; do 
	echo $file
	frac_coding=`head -8 /mnt/lab_data/montgomery/nicolerg/motrpac/rna/metrics/$file | tail -1 | cut -f17`
	frac_utr=`head -8 /mnt/lab_data/montgomery/nicolerg/motrpac/rna/metrics/$file | tail -1 | cut -f18`

	total_frac=`awk -v frac_coding="$frac_coding" -v frac_utr="$frac_utr" 'BEGIN{print frac_utr + frac_coding}'`

	awk -v total_frac="$total_frac" 'BEGIN{print total_frac * 100}'
done

# get number of reads in raw FASTQ files

for run in {1..3}; do
	indir=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN${run}
	for file in `ls ${indir} | grep "R1_001.fastq.gz"`; do 
		sample=`echo $file | sed "s/_R1_001.fastq.gz//"`
		count=`zcat ${indir}/${file} | wc -l`
		y=2 # line count/4 == reads in one file; multiply by 2 for total number of reads in PE files
		echo $sample
		echo $((count / y))
	done
done

# get number of reads in trimmed FASTQ files

for run in {1..3}; do
	indir=/mnt/lab_data/montgomery/nicolerg/motrpac/rna/FASTQ_RUN${run}/cutadapt
	for file in `ls ${indir} | grep "R1_001.trimmed.fastq.gz"`; do 
		sample=`echo $file | sed "s/_R1_001.trimmed.fastq.gz//"`
		count=`zcat ${indir}/${file} | wc -l`
		y=2 # line count/4 == reads in one file; multiply by 2 for total number of reads in PE files
		echo $sample
		echo $((count / y))
	done
done
