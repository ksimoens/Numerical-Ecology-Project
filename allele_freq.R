library(poppr)

# Load data with all the SNPs and individuals
load('Input/lobster_1278ind_79snps_40pop.RData')

# Combine different temporal samples
popNames(data_filt) <- gsub("Sar13", "Sar", popNames(data_filt))
popNames(data_filt) <- gsub("Sar17", "Sar", popNames(data_filt))
popNames(data_filt) <- gsub("Idr16", "Idr", popNames(data_filt))
popNames(data_filt) <- gsub("Idr17", "Idr", popNames(data_filt))
print(summary(data_filt$pop))

# Calculate allele frequencies per population from genotype
al.freq.tot <- data.frame(rraf(data_filt, by_pop=TRUE, correction = FALSE))
# Only retain one allele (q = 1 - p)
al.freq.tot <- al.freq.tot[, seq(1, dim(al.freq.tot)[2], 2)]

print(al.freq.tot)

write.csv(al.freq.tot,"Output/allele_freq_total.csv")

# Outlier SNPs
outliers <- c("11291","15581","32358","42395","53935","58053",
              "65064","65576")
snpdata.sel <- data_filt[loc=outliers]

al.freq.sel <- data.frame(rraf(snpdata.sel, by_pop=TRUE, correction = FALSE))
al.freq.sel <- al.freq.sel[, seq(1, dim(al.freq.sel)[2], 2)]

print(al.freq.sel)

write.csv(al.freq.sel, "Output/allele_freq_sel.csv")

# Neutral SNPs = not outliers
neutral <- rownames(as.matrix(data_filt@all.names))
neutral <- neutral[!(neutral %in% outliers)]
snpdata.neut = data_filt[loc=neutral]

al.freq.neut <- data.frame(rraf(snpdata.neut, by_pop=TRUE, correction = FALSE))
al.freq.neut <- al.freq.neut[, seq(1, dim(al.freq.neut)[2], 2)]

print(al.freq.neut)

write.csv(al.freq.neut, "Output/allele_freq_neut.csv")
