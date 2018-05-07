## Phase and Impute using 1000 Genome and Impute2


Download,
```
wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz

wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3_chrX.tgz

tar -xvzf *.tgz

export 1000Genome=$PATH:/path/to/1000Genome

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

#### Step 3: Convert to gen file 

```
gtool -P --ped MSchip_RNAseq_recode_chr1.ped --map MSchip_RNAseq_recode_chr1.map --og MSchip_RNAseq_recode_chr1.gen
```

#### Step 4: Prephase on the selected region 

```
impute2 -prephase_g -m $1000Genome/genetic_map_chr1_combined_b37.txt -g MSchip_RNAseq_recode_chr1.gen -int 652566 249218992  -Ne 81 -o MSchip_RNAseq_chr1.prephasing.impute2 -allow_large_regions
```

*250000 bps on both sides of the gene level….since I am focusing on the exome/RNAseq regions*


#### Step 5: Imputation into pre-phased haplotypes

```
impute2 -use_prephased_g -m $1000Genome/genetic_map_chr1_combined_b37.txt -h $1000Genome/1000GP_Phase3_chr1.hap.gz -l $1000Genome/1000GP_Phase3_chr1.legend.gz  -known_haps_g MSchip_RNAseq_chr1.prephasing.impute2_haps  -int 652566 249218992 -Ne 81  -o MSchip_RNAseq_chr1.phased.chunk1.impute2 -phase -allow_large_regions

```
