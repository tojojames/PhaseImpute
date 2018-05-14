### source: https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/

rm(list = ls()); gc();
args = commandArgs(trailingOnly=TRUE)
# load haplotypic counts
phaser_test_case_gene_ae = read.delim(args[1])
#phaser_test_case_gene_ae = read.delim("phaser_test_case_gene_ae.txt")
# select only genes with sufficient coverage
cov8 = subset(phaser_test_case_gene_ae, phaser_test_case_gene_ae$totalCount >= 8)
# perform a binomial test for deviation from 0.5
cov8$binom_p = apply(cov8[,c("aCount","bCount")], 1, function(x) binom.test(x[1],x[1]+x[2],p=0.5)$p.value)
# perform multiple testing correction with FDR
cov8$binom_q = p.adjust(cov8$binom_p, method = "fdr")
# plot haplotype A versus haplotype, highlight sites with significant imbalance (FDR < 5%)
#plot(cov8$bCount, cov8$aCount, log="xy", col=(cov8$binom_q<0.05)+1, ylab="Haplotype B Count", xlab="Haplotype A Count")
#abline(0,1,col="grey")
#legend("topleft",c("No significant imbalance","Significant imbalance"),pch=c(15,15),col=c(1,2))
write.table(cov8,file=args[2],sep = "\t",row.names =F)
