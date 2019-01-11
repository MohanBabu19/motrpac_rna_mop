# Ancillary files 

The files in this directory are intermediate files necessary to run all steps of the RNA pipeline. They are provided for the convenience of the user. This document provides additional details about how these files were generated.

### Table of contents:  

* adapter.tsv

refFlat files:

  * Rattus_norvegicus.Rnor_6.0.94.refFlat
  * gencode.v29.annotation.refFlat

globin files:  

  * rat_hemoglobin_genes.txt
  * human_hemoglobin_genes.txt
  * Rattus_norvegicus.Rnor_6.0.94.hemoglobin.fa
  * gencode.v29.annotation.hemoglobin.fa

ribosomal RNA files:  

  * rn_rRNA.fa
  * hg_rRNA.fa

### adapter.tsv

This file contains several adapter sequences used for adapter trimming. Not all of them are relevant for our data, but including additional adapters does not compromise adapter trimming. The most necessary sequences for adapter trimming for the NuGEN kits are `TruSeq1` and `TruSeq Universal`. 

### refFlat files 

* `Rattus_norvegicus.Rnor_6.0.94.refFlat`: refFlat file for rat (Rnor_6.0, Release 94)  
* `gencode.v29.annotation.refFlat`: refFlat for for human (GRCh38.p12, Release 29)  

refFlat files are a file format similar to the UCSC genePred format (read more here: https://genome.ucsc.edu/FAQ/FAQformat.html#format9 and here: https://genome.ucsc.edu/goldenpath/gbdDescriptions.html). They can be generated from GTF files using the UCSC gtfToGenePred tool (https://bioconda.github.io/recipes/ucsc-gtftogenepred/README.html) plus a little extra reformatting. The following code was used to generate refFlat files using the GTF files specified in step B.1 of the RNA MOP: 

Parameters:  

* `GTF_FILE`: Same as in steps 1 and 4.1. 
* `OUTPUT_DIR`: Output directory for refFlat file 
* `REF_FLAT`: Path for output refFlat file 

```bash 
# first, convert GTF to genePred

gtfToGenePred -genePredExt $GTF_FILE genePred.tmp

# now convert genePred to refFlat

awk 'BEGIN{FS=OFS="\t"};{print $12, $1, $2,$3,$4,$5,$6,$7,$8,$9,$10}' genePred.tmp > $REF_FLAT
rm genePred.tmp
```

### globin files 

* `rat_hemoglobin_genes.txt`: Table of Ensembl rat hemoglobin genes used to generate a FASTA file of rat globin sequences
* `human_hemoglobin_genes.txt`: Table of Ensembl human hemoglobin genes used to generate a FASTA file of human globin sequences 
* `Rattus_norvegicus.Rnor_6.0.94.hemoglobin.fa`: FASTA file of rat globin sequences
* `gencode.v29.annotation.hemoglobin.fa`: FASTA file of human globin sequences 

Only the alpha and beta globin subunit sequences are included because they are sufficient to detect globin contamination. The other subunits account for a much smaller fraction of globin-derived reads in non-globin depleted blood samples. 

There are some differences in annotation between Ensembl and NCBI rat globin genes. The gene name corresponding to ENSRNOG00000029886 is actually Hba-a2, but this annotation overlaps with NCBI Hba1 annotation, and Ensembl does not otherwise have an Hba1 annotation, so ENSRNOG00000029886 is given the Hba1 designation in `rat_hemoglobin_genes.txt`. 

`rat_hemoglobin_genes.txt` and `human_hemoglobin_genes.txt` are not necessary to run the pipeline. However, the FASTA files are. The following code was used to generate the FASTA files:

```bash
# for rat ############################################

outfile=Rattus_norvegicus.Rnor_6.0.94.hemoglobin.fa
for gene in ENSRNOG00000058105 ENSRNOG00000029886 ENSRNOG00000047321; do
	wget -q --header='Content-type:text/x-fasta' "https://rest.ensembl.org/sequence/id/${gene}?type=genomic;species=rat"  -O - >> $outfile
done

# for human ##########################################

outfile=gencode.v29.annotation.hemoglobin.fa
for gene in ENSG00000206172 ENSG00000188536 ENSG00000244734; do
	wget -q --header='Content-type:text/x-fasta' "https://rest.ensembl.org/sequence/id/${gene}?type=genomic;species=human"  -O - >> $outfile
done
```

### ribosomal RNA files

* `rn_rRNA.fa`: FASTA file of 28S, 45S, 18S, 5S ribosomal RNA rat sequences 
* `hg_rRNA.fa`: FASTA file of 12S, 16S, 45S, 18S, 5.8S, 28S, 5S human rRNA sequences 

These FASTA files are NCBI sequences provided by Yongchao Ge.
