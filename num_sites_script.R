
library(tidyr)
library(dplyr)
file <- read.table("number.sites.analysis.tsv",strip.white=T,header=F)
skip_last <- head(file,-1)
new_file <- separate(skip_last,V2,c("V2","V3"),sep="/")
wide_file <- pivot_wider(new_file,values_from="V1",names_from="V3")
write.table(wide_file,"number.sites.long.tsv",sep="\t",row.names=F)
