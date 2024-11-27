library(vcfR)
library(tidyr)
library(ggplot2)
library(dplyr)

#here we read in the filtered VCF, which we have filtered for biallelic sites
data <- read.vcfR("combined_final_tagged.vcf.gz")

#now we extract the Allele Depth, which we will use to filter heterozygotes
ad <- extract.gt(data,"AD")

#we need to extract the reference depth, which is before the comma, and then convert to numeric
ref <- as.data.frame(gsub(",.*","",ad))
ref <- apply(ref,2,as.numeric)
#head(ref)

#we need to extract the alternate depth, which is after the comma, and then convert to numeric
alt <- as.data.frame(gsub(".*,","",ad))
alt <- apply(alt,2,as.numeric)
#head(alt)

#now we convert to a ratio, between 0 and 1
ratio <- ref/(ref+alt)
head(ratio)

#we calculate the number of heterozygous sites per position
num_mixed <- rowSums(ratio > 0.25 & ratio < 0.75,na.rm=T)

print("Does lengths of data@gt and ratio match, they should NOT")
print(length(data@gt) == length(ratio))

#we want to replace values in data@gt based on values in ratio, but the problem is that data@gt has an extra column of "FORMAT"
FORMAT <- data@gt[,"FORMAT"]
data@gt <- data@gt[,-1]

print("Does lengths of data@gt and ratio match, they should currently")
print(length(data@gt) == length(ratio))

#now replace all genotypes with NA if ratio is bad
#data@gt[ratio > 0.25 & ratio < 0.75] <- "0:0,0:0:0:0,0"
#and we want to exclude those with more than two "heterozygotes", since that is not possibe
data@fix <- data@fix[num_mixed < 2,]
data@gt  <- data@gt[num_mixed < 2,]
FORMAT <- FORMAT[num_mixed <2]

#now we need to add FORMAT column back
data@gt <- cbind(FORMAT,data@gt)

print("Does lengths of data@gt and ratio match, they should not anymore")
print(length(data@gt) == length(ratio))

#this file should now be smaller
head(data@gt)

write.vcf(data,file="combined_final_tagged_nohets.vcf.gz")
