### Allele specific expression study with Phaser

#### Steps to phase and impute using 1000 Genome and Impute2 to study allelic imbalance  (BETA version)


Download,
```
wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz

wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3_chrX.tgz

tar -xvzf *.tgz

export thousandGenome=$PATH:/path/to/1000GP_Phase3/

```

#### This is an example on a selected region of Chromosome 1

*For simplicity and demo, testing only on one chromosome.* 

#### Step 1: Convert to plink map and ped format from bed bim format

```
plink --bfile MSchip_RNAseq --recode --tab --out MSchip_RNAseq
```

#### Step 2: Select based on chromosome

```
plink --bfile MSchip_RNAseq --recode --out MSchip_RNAseq_recode_chr1 --chr 1
```

#### Step 3: Convert to GEN file (genotype file) 

```
gtool -P --ped MSchip_RNAseq_recode_chr1.ped --map MSchip_RNAseq_recode_chr1.map --og MSchip_RNAseq_recode_chr1.gen
```

#### Step 4: Prephase on the selected region 

```
impute2 -prephase_g -m $thousandGenome/genetic_map_chr1_combined_b37.txt -g MSchip_RNAseq_recode_chr1.gen -int 1 249250621  -Ne 20000 -o MSchip_RNAseq_chr1.prephasing.impute2 -allow_large_regions
```

*File "hg19_chrLen" for chromosome size can be found in this repository*

##### Strand information for different chips (eg.Illumina) can be found here,
```
http://www.well.ox.ac.uk/~wrayner/strand/index.html#Illumina
```

#### Step 5: Imputation into pre-phased haplotypes

```
impute2 -use_prephased_g -m $thousandGenome/genetic_map_chr1_combined_b37.txt -h $thousandGenome/1000GP_Phase3_chr1.hap.gz -l $thousandGenome/1000GP_Phase3_chr1.legend.gz  -known_haps_g MSchip_RNAseq_chr1.prephasing.impute2_haps  -int 1 10000000 -Ne 20000  -o MSchip_RNAseq_chr1.phased.chunk1.impute2 -phase -allow_large_regions -strand_g MSchip.strand
```



[[Note in case of the following error(seldom get into this error)
*"ERROR: Individual 47 (1-indexed) has invalid alleles '00.333' at position 215314870 "*
If SNPs not determined by prephasing, it confuses this step of imputation  and result in following error.One way to fix using the following grep to remove the confused phases and continue to previous step of prephasing (step 4)]] 
###### grep -v '0.333' MSchip_RNAseq_chr1.prephasing.impute2_haps > MSchip_RNAseq_chr1.prephasing.impute2_haps_corrected

The following script can concatenate all the pieces(eg. 24 pieces) of a single chromosome(for example chromosome). This can reduce the memory usage,
###### cat MSchip_RNAseq_chr1.phased.chunk{1..24}.impute2 |sort -k3n |sed 's/---/chr1/g' >MSchip_RNAseq_chr1.phased.ALLchunk.merged.impute2

###### cat MSchip_RNAseq_chr1.phased.chunk{1..24}.impute2_haps |sort -k3n |sed 's/---/chr1/g' >MSchip_RNAseq_chr1.phased.ALLchunk.merged.impute2_haps.haps

#### Step 6: Converting haps into VCF format
```
shapeit -convert --input-haps MSchip_RNAseq_chr1.phased.ALL_chunks.impute2_haps --output-vcf MSchip_RNAseq_chr1.phased.ALL_chunks.impute2_haps.vcf.gz --thread 8
tabix -f -p vcf MSchip_RNAseq_chr1.phased.ALL_chunks.impute2_haps.vcf.gz
```
*It requires a sample file for example, "MSchip_RNAseq_chr1.phased.ALL_chunks.impute2_haps.sample"*

Sample format can be found here http://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#hapsample

#### Step 7: Subseting vcf file based on Sample (for Phaser ASE algorithm)
*for example, sample "N000271" on chromosome 1*
```
vcf-subset -c N000271 MSchip_RNAseq_chr1.phased.ALL_chunks.impute2_haps.vcf.gz|bgzip -c >N000271_MSchip_RNAseq_chr1.phased.ALL_chunks.impute2_haps.vcf.gz

tabix -f -p vcf N000271_MSchip_RNAseq_chr1.phased.ALL_chunks.impute2_haps.vcf.gz
```

#### Step 8: ASE test using phASER tool

Reference to Phaser and package https://github.com/secastel/phaser/tree/master/phaser

'chr' label corrected version for hla bed file in the repository or download,
###### wget ftp://ftp.nygenome.org/sec/phaser/annot_files_v100/hg19_hla.bed.gz
###### wget ftp://ftp.nygenome.org/sec/phaser/annot_files_v100/hg19_haplo_count_blacklist.bed.gz

###### samtools index N000271_Aligned_sorted_sort.bam 

```
python /phaser/phaser.py --vcf N000271_MSchip_RNAseq_chr1.phased.ALL_chunks.impute2_haps.vcf.gz --bam  N000271_Aligned_sorted_sort.bam --paired_end 1 --mapq 255 --baseq 10 --sample N000271 --blacklist hg19_hla.chr.bed --haplo_count_blacklist hg19_haplo_count_blacklist.chr.bed --threads 4 --o N000271_MSrepASE
```
#### Step 9: Haplotype expression quantifications using phASER Gene AE tool

'chr' label corrected version in the repository or download
###### wget ftp://ftp.nygenome.org/sec/phaser/hg19_ensembl.bed.gz

```
python /phaser_gene_ae/phaser_gene_ae.py --haplotypic_counts N000271_MSrepASE.haplotypic_counts.txt --features hg19_ensembl.chr.bed --o N000271_Phaser_gene_ae.txt
```
R script "imbalance_FDR.R" file can be downloaded from this Repository. This test uses a binomial test to determine if genes have significant allelic imbalance (deviation from 50/50 expression).
```
Rscript --verbose imbalance_FDR.R N000271_Phaser_gene_ae.txt N000271_Phaser_gene_ae_FDR.txt 
```

###### cat *_Phaser_gene_ae_FDR.txt >Phaser_result_allsamples.txt
