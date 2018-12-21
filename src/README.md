These scripts provide EXAMPLES of how to run all parts of the RNA MOP. These scripts will NOT run without file paths being changed. Other user-specific changes may need to be made as well. 

`1_rna_pipeline.sh` includes the blocks that correspond to Section A in the MOP. These blocks can be run on all raw FASTQ files, regardless of species/sample.  

`2_human_rna_pipeline.sh` and `2_rat_rna_pipeline` are reference-specific and correspond to Section B in the MOP.  

`3_rna_pipeline.sh` roughly corresponds to Section C in the MOP. These blocks can be run on all outputs of previous steps, regardless of species, to collect statistics for all samples. 
