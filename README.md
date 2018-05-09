## Phase and Impute using 1000 Genome and Impute2


Download,
```
wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz

wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3_chrX.tgz

tar -xvzf *.tgz

export thousandGenome=$PATH:/path/to/1000GP_Phase3/

```

### This is an example on a selected region of Chromosome 1

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

*Repository file "hg19_chrLen" for chromosome size*

*if SNPs not determined by prephasing, it confuses imputation (next) step and result in following error in step 5*
*"ERROR: Individual 47 (1-indexed) has invalid alleles '00.333' at position 215314870"*
###### grep -v '0.333' MSchip_RNAseq_chr1.prephasing.impute2_haps > MSchip_RNAseq_chr1.prephasing.impute2_haps_corrected
*Most cases this error doesnt occur,so ignore this correction with grep!*


##### Strand information for different chips (eg.Illumina) can be found here,
```
http://www.well.ox.ac.uk/~wrayner/strand/index.html#Illumina
```

#### Step 5: Imputation into pre-phased haplotypes

```
impute2 -use_prephased_g -m $thousandGenome/genetic_map_chr1_combined_b37.txt -h $thousandGenome/1000GP_Phase3_chr1.hap.gz -l $thousandGenome/1000GP_Phase3_chr1.legend.gz  -known_haps_g MSchip_RNAseq_chr1.prephasing.impute2_haps  -int 1 10000000 -Ne 20000  -o MSchip_RNAseq_chr1.phased.chunk1.impute2 -phase -allow_large_regions -strand_g MSchip.strand
```
The following script can concatenate all the pieces(eg. 24 pieces) of a single chromosome(for example chromosome). This can reduce the memory usage,
###### ls MSchip_RNAseq_chr1.phased.chunk{1..24}.impute2 |cat |sort -k3n |sed 's/---/1/g' >MSchip_RNAseq_chr1.phased.ALLchunk.merged.impute2

###### ls MSchip_RNAseq_chr1.phased.chunk{1..24}.impute2_haps |cat |sort -k3n |sed 's/---/1/g'>MSchip_RNAseq_chr1.phased.ALLchunk.merged.impute2_haps.haps

#### Step 6: Converting haps into VCF format
```
shapeit -convert --input-haps MSchip_RNAseq_chr1.phased.ALL_chunks.impute2_haps --output-vcf MSchip_RNAseq_chr1.phased.ALL_chunks.impute2_haps.vcf --thread 8
```
*It requires a sample file for example, "MSchip_RNAseq_chr1.phased.ALL_chunks.impute2_haps.sample"*

Sample format can be found here http://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#hapsample
