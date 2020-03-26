# HNSCC-WES-SNV-Workflow
Workflow for identifying single-nucleotide-variants (SNV) in Head and Neck Squamous Cell Carcinoma (HNSCC) Whole Exome Sequencing (WES) samples

# Version Notes
- These analyses were carried out on the OHSU cluster computing system (Exacloud) using CentOS 7.5.1804 unless otherwise noted
- Exacloud uses the job scheduler, Slurm, for job submissions.  See separate files for Slurm submit scripts. 
- Alignment of sequencing reads was accomplished using the Burrows-Wheeler Aligner.  The version used was bwa-0.7.15
- GATK version 3.6 (Picard included)
- All python scripts were run on python version 2.7.13

# Workflow
**Step 1) Alignment of sequencingn reads to the hg19** 
The program BWA is used to align sequencing reads to a the hg19 human reference genome.  Sequencing reads from each lane are aligned individually.  They will later be combined into a single .bam file for each sample.  Separate alignments were carried out for the tumor and matched normal samples.  

The "-M" flag tells the program to mark shorter split hits as secondary (for Picard compatibility)

The "-t" flag tells the program how many threads to use.

```
fasta_dir=<path to fasta file>
alignment_run=<sample ID>
ID = <Read group identifier>
PL="Illumina"
LB=$alignment_run"T"
SM=$alignment_run"T"
input_F=<fastq file for foward reads>
input_R=<fastq file for reverse reads>

bwa-0.7.15/bwa mem -M -R "@RG\tID:$ID\tPL:$PL\tLB:$LB\tSM:$SM" -t 32 $fasta_dir/hg19.fa $input_F $input_R
```
Notes on read group fields:
- ID: This is the read group identifier.  It designates which read group each read belongs to.  With Illumina data, ID's are a combination of the flowcell identifier and the lane identifier.  
- PL: This is the platform identifier.  It tells which sequencing technology was used.  In our case, we used Illumina for all our sequencing.
- LB: This is the DNA preparation library identifier.  The LB field can be used to identify molecular duplicates in case the same DNA library was sequenced on multiple lanes.
- SM: This is the sample identifier.  This is the unique name given the the individual sample.  GATK tools treat all read groups with the same SM value as containing sequencing data for the same sample, and this is also used for the sample column in the VCF file.

**Step 2) Merge BAM files**
Picard Tools -- utilizing Java -- is used to take the individual alignemnts from each of the WES lanes and merge them into one complete .bam file.  Separate merges were carried out for the tumor and matched normal samples.

```
#for n lanes

ALIGNMENT_RUN=<Sample ID>
file_1=aligned_lane_1.bam"
file_2=aligned_lane_2.bam"
#            .
#            .
#            .
file_<n>=aligned_lane_<n>.bam"

java -Xmx8G -jar picard.jar MergeSamFiles \
    I=$file_1 \
    I=$file_2 \
#           .
#           .
#           .
    I=$file_<n> \
    O=aligned.bam
```

**Step 3) Mark Duplicates**
Picard Tools -- utilizing Java -- is used to mark any duplicate reads from a sequence alignment file.  Separate MarkDuplicates runs were carried out for the tumor and matched normal samples.

The CREATE_INDEX=true command will create an index for the outputed .bam file.  This is needed for downstream GATK tools.

The VALIDATION_STRINGENCY=SILENT is set per GATK pipeline

A Temp Directory (TMP_DIR) can be designated if needed

```
ALIGNMENT_RUN=<Sample ID>

java -Xmx8G -jar $PICARD_DIR/picard.jar MarkDuplicates \
I=aligned.bam \
O=rg_added_aligned_MarkedDup.bam \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
M=Markdup_metrics.txt \
TMP_DIR=<path to appropriate location for a temp directory>
```

**Step 4) Base Quality Score Recalibration**
The GATK tool, "BaseRecalibrator" is used to take a .bam file, the GATK reference for hg19 as well as a file containing the known sites of variation in hg19 and produces a recal_data.table as the output.  This consists of several tables.
- The list of arguments
- The quantized qualities talbe
- The recalibration table by read group
- The recalibration talbe by quality score
- The recalibration table for all the optional covariates
Separate Base Quality Score Recalibration (BQSR) runs were carried out for the tumor and matched normal samples.

```
ALIGNMENT_RUN=<Sample ID>

REF="hg19.fa"
INPUT_FILE="rg_added_aligned_MarkedDup.bam"
KNOWN_SITES_1="sorted_dbsnp_138.hg19.vcf"
KNOWN_SITES_2="Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
INTERVALS="nexterarapidcapture_exome_targetedregions_v1.2.bed"

java -jar GenomeAnalysisTK.jar -T BaseRecalibrator \
    -R $REF \
    -I $INPUT_FILE \
    -knownSites $KNOWN_SITES_1 \
    -knownSites $KNOWN_SITES_2 \
    -L $INTERVALS \
    --interval_padding 50 \
    -o recal_data.table
```

**Step 5) Apply the BQSR to the variants file**
The GATK tool, "PrintReads" was used to take the report from the GATK BaseRecalibrator and actually recalibrate the quality scores in the called bases in the alignment .bam file.  Separate PrintReads runs were carried out for the tumor and matched normal samples.

```
ALIGNMENT_RUN=<Sample ID>

REF="hg19.fa"
INPUT_FILE="rg_added_aligned_MarkedDup.bam"
TABLE="recal_data.table"

java -jar GenomeAnalysisTK.jar -T PrintReads -R $REF -I $INPUT_FILE -BQSR $TABLE -o recal_reads.bam
```

**Step 6) Tumor-Only Mutect2 run to create prePON file**
The GATK tool, "MuTect2" was used to prepare the normal samples to be combined into a panel of normals.  Each normal sample is run separately.

```
ALIGNMENT_RUN=<Sample ID>

REF="hg19.fa"
INPUT_FILE="recal_reads.bam"
DBSNP="sorted_dbsnp_138.hg19.vcf"
COSMIC="sorted_CosmicCodingMuts.vcf"
INTERVALS="nexterarapidcapture_exome_targetedregions_v1.2.bed"

java -jar GenomeAnalysisTK.jar -T MuTect2 \
    -R $REF \
    -I:tumor $INPUT_FILE \
    --dbsnp $DBSNP \
    --cosmic $COSMIC \
    -L $INTERVALS \
    --interval_padding 50 \
    --artifact_detection_mode \
    -o prePON.vcf
```

**Step 7) Combine prePON files into a Panel of Normals**
The GATK tool, "CombineVariants" was used to create a PONs

```
#For n samples

REF="hg19.fa"

java -jar GenomeAnalysisTK.jar -T CombineVariants \
    -R $REF \
    -V <path to directory for Sample 1>/prePON.vcf \
    -V <path to directory for Sample 2>/prePON.vcf \
    #                       .
    #                       .
    #                       .
    -V <path to directory for Sample n>/prePON.vcf \
    -minN 2 \
    --setKey \"null\" \
    --filteredAreUncalled \
    --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
    -o MuTect2_PON.vcf
```

**Step 8) Estimate level of contamination in tumor samples**
The GATK tool, "ContEst" was used to determine the percent contamination of an input bam.  This method uses the nomral bam for genotypeing on-the-fly

```
ALIGNMENT_RUN=<Sample ID>

REF="hg19.fa"
TUMOR=<path to directory for tumor sample>/recal_reads.bam"
NORMAL=<path to directory for matched normal sample>/recal_reads.bam"
POPFILE="hapmap_3.3_hg19_pop_stratified_af.vcf"

java -jar GenomeAnalysisTK.jar -T ContEst \
    -R $REF \
    -I:eval $TUMOR \
    -I:genotype $NORMAL \
    --popfile $POPFILE \
    -o sample.conest
```

**Step 9) Call variants in tumor sample relative to the matched normal sample**
The GATK tool, "MuTect2" was used to call mutations in the tumor samples relative to the normals.  Each sample was run separately.  MuTect2 has the ability to use COSMIC data in conjunction with dbSNP to adjust the threshold for evidence of a variant in the normal.  If a variant is present in dbSNP, but not in COSMIC then more evidence is required from the normal sample to prove the variant is not present in germline.  rsIDs from the dbSNP file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate.  dbSNP overlap is only used to require more evidence of absence in the normal if the variant in question has been seen before in germline.

Note: The contamination results from Step 8 are reported as percentages.  If the sample.conest reports a contaminatio of 0.3, then it should be entered into MuTect2 as 0.003.

```
ALIGNMENT_RUN=<Sample ID>

REF="hg19.fa"
TUMOR=<path to directory for tumor sample>/recal_reads.bam"
NORMAL=<path to directory for normal sample>/recal_reads.bam"
DBSNP="sorted_dbsnp_138.hg19.vcf"
COSMIC="sorted_CosmicCodingMuts.vcf"
INTERVALS="nexterarapidcapture_exome_targetedregions_v1.2.bed"
CONTAMINATION=<results from sample.conest>
PANOFNORMS="MuTect2_PON.vcf"

java -jar GenomeAnalysisTK.jar -T MuTect2 \
    -R $REF \
    -I:tumor $TUMOR \
    -I:normal $NORMAL \
    --dbsnp $DBSNP \
    --cosmic $COSMIC \
    -L $INTERVALS \
    --interval_padding 50 \
    --contamination_fraction_to_filter $CONTAMINATION \
    -PON $PANOFNORMS \
    -o variants.vcf
```

**Step 10) Filter variant file by variants that only have the "PASS" designation in their filtering field**
The python program, "filter_VCF_by_PASS.py", processes the .vcf file from Mutect2 output and filters it so only the variants with a "PASS" in their filter feild will be carried forward to the new file.

```
python filter_VCF_by_PASS.py
```

**Step 11) Format the variants for processing via the Provean website**
The pyhton program, "preProvean.py", takes the variants file and formats each variant so they can be evaluated by Provean

```
python preProvean.py
```

**Step 12) Add allele frequency to the Provean output**
The python program, "add_af_to_output.py", takes the Provean results and adds in the af from the .vcf file to the output

```
python add_af_to_output.py
```

**Step 13) Add COSMIC and ExAC information**
The Python program, "add_COSMIC_ExAC_to_output.py", takes the Provean results (with af added) and adds in any COSMIC and/or ExAC annotations if they exist

```
python add_COSMIC_ExAC_to_output.py
```

# Refernces
1) Read Groups.  https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
2) Li H. and Durbin R. Fast and accurate short read alignment with Burrows-Wheeler Transform.  Bioinformatics 2009;25:1754-60.
