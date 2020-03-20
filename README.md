# HNSCC-WES-SNV-Workflow
Workflow for identifying single-nucleotide-variants in Head and Neck Squamous Cell Carcinoma Whole Exome Sequencing samples

Version Notes
- These analyses were carried out on the OHSU cluster computing system (Exacloud) using CentOS 7.5.1804 unless otherwise noted

Workflow
Step 1) Alignment of sequencingn reads to the Hg19 
This submit script uses the program BWA and the function mem to align sequencing reads to a the hg19 human
reference genome.  This one aligns sequencing reads from each lane individually.  They will later be combined
into a single .bam file for each sample.

The "-M" flag tells the program to Mark shorter split hits
as secondary (for Picard compatibility)

The "-t" flag tells the program how many threads to use.

The output file is specified in "--output" above

'''
fasta_dir=<path to fasta file>
alignment_run=<sample ID>
ID="DNA170613MK.L00"${SLURM_ARRAY_TASK_ID}
PL="Illumina"
LB=$alignment_run"T"
SM=$alignment_run"T"
input_F=DNA170613MK_DNA_Tissue_"$alignment_run"_S24_L00"${SLURM_ARRAY_TASK_ID}"_R1_001.fastq
input_R=DNA170613MK_DNA_Tissue_"$alignment_run"_S24_L00"${SLURM_ARRAY_TASK_ID}"_R2_001.fastq

bwa-0.7.15/bwa mem -M -R "@RG\tID:$ID\tPL:$PL\tLB:$LB\tSM:$SM" -t 32 $fasta_dir/hg19.fa $input_F $input_R
'''
