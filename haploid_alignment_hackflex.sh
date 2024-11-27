#!/bin/bash
#Tells which version of bash to use, put it on the start of each script.

REFERENCE=/mnt/scratch/auxie001/outgroup/ref_genome/Af293_combined.fna
BAMS=/mnt/LTR_userdata/auxie001/groeniii/ftps.baseclear.com/141184/raw_sequences
#then we need to index the genome file 
#bwa-mem2 index $REFERENCE


#then we need to map the reads to the reference genome 

#For loop for the analysis of all 24 samples at the same time 

#for sample in 102C11 11A2 11A5 11A6 16C32 29C8 38C31 39B10 40A20 48A4 48A6 55C1 55C4 58C25 60C3 82C21 86C22 88C20 98C11 98C12 9C22
#do echo $sample
#/mnt/LTR_userdata/auxie001/programs/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 8 -R '@RG\tID:'$sample'\tSM:'$sample'\tLB:hack_flex' $REFERENCE $BAMS/$sample*R1*.gz $BAMS/$sample*R2*.gz | samtools view -b | samtools fixmate -@ 6 -m - - | samtools sort -@ 8 -m 3G - | samtools markdup -@ 6 - $sample.sorted.bam
#samtools index $sample.sorted.bam
#/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk HaplotypeCaller -I $sample.sorted.bam -R $REFERENCE  -O $sample.g.vcf.gz -ERC GVCF -ploidy 1
#done


#for sample in S1  S10  S11  S12  S13  S14  S15  S16  S17  S18  S19  S2  S20  S21  S22  S23  S24  S3  S4  S5  S6  S7  S8  S9
#do echo $sample
#/mnt/LTR_userdata/auxie001/programs/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 8 -R '@RG\tID:'$sample'\tSM:'$sample'\tLB:hack_flex' $REFERENCE $sample/*_1.fq.gz $sample/*_2.fq.gz | samtools view -b | samtools fixmate -@ 6 -m - - | samtools sort -@ 8 -m 3G - | samtools markdup -@ 6 - $sample.sorted.bam
#samtools index $sample.sorted.bam
#/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk HaplotypeCaller -I $sample.sorted.bam -R $REFERENCE  -O $sample.g.vcf.gz -ERC GVCF -ploidy 1
#done

#The GVCF files need to be combined into one output that lists where in the genome there are variants

#find *.g.vcf.gz > haploid_variant_files.list

#/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk CombineGVCFs -R $REFERENCE -O haploid_test_hack_flex.g.vcf.gz --variant haploid_variant_files.list
#Joined genotyping - the combined variants will be spit out by this
#/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk GenotypeGVCFs -R $REFERENCE -O haploid_test_hack_flex.vcf.gz -V haploid_test_hack_flex.g.vcf.gz

#The variants will be filtered based on their likelyhood  (see criteria in script); filter out low quality variants
#/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk SplitVcfs -I haploid_test_hack_flex.vcf.gz --SNP_OUTPUT combined_SNPs.vcf.gz --INDEL_OUTPUT combined_INDELs.vcf.gz --STRICT false
#echo "starting filtering SNPs"
#bcftools filter -O u -s "QD" -m + -e "QD < 2.0" combined_SNPs.vcf.gz | bcftools filter -O u -s "MQ" -m + -e "MQ < 50.0" |bcftools filter -O u -s "FS" -m + -e "FS > 60.0" | bcftools filter -O u -s "MQRank" -m + -e "MQRankSum < -12.5" | bcftools filter -O u -s "ReadPos" -m + -e "ReadPosRankSum < -8.0" | bcftools view -O z -f "PASS" > combined_SNPs_filtered.vcf.gz
#echo "finished filtering SNPs"
#bcftools index combined_SNPs_filtered.vcf.gz
#echo "starting filtering INDELs"
#bcftools filter -O u -s "QD" -m + -e "QD < 2.0" combined_INDELs.vcf.gz | bcftools filter -O u -s "FS" -m + -e "FS > 200.0" | bcftools filter -O u -s "ReadPos" -m + -e "ReadPosRankSum < -20.0" | bcftools view -O z -f "PASS" > combined_INDELs_filtered.vcf.gz
#echo "finished filtering INDELs"
#bcftools index combined_INDELs_filtered.vcf.gz
#bcftools concat -O z -a combined_INDELs_filtered.vcf.gz combined_SNPs_filtered.vcf.gz > combined_filtered.vcf.gz
#/mnt/LTR_userdata/auxie001/programs/gatk-4.2.6.1/gatk IndexFeatureFile -I combined_filtered.vcf.gz

#rm combined_INDELs* combined_SNPs*

#bcftools stats -s - combined_filtered.vcf.gz > combined_filtered_stats.txt

#Sample	ng/ul	 	Sample	ng/ul
#1. 11A2 	13. 40A20	2.1
#2. 11A6 	14. 48A4	3.3
#3. 48A6 	15. 98C11	3.2
#4. 55C1 	16. 11A5	3.4
#5. 55C4	17. 102C11	4.1
#6. 60C3 	18. 29C8	2.04
#7. 82C21 	19. 60C3	2.5
#8. 86C22 	20. 88C20	1.9
#9. AfIR964 	21. 16C32	3.8
#10. AfIR974 	22. 39B10	3.7
#11. 9C22 	23. 58C25	6.0
#12. 38C31 	24. 98C12	1.5

#Note that after looking at the data, sample 23 and sample 24 were mixed up.


#bcftools view -m 2 -M 2 -s ^S9,S10 combined_filtered.vcf.gz | bcftools filter -e "AN<38" | bcftools filter -e "AF>0.98" -o combined_final_tagged.vcf.gz

#bcftools stats -s - combined_final_tagged.vcf.gz > combined_final_tagged_stats.txt

#now process in R, currently on Ben's laptop

#Rscript bo_filtering.R

#zcat combined_final_tagged_nohets.vcf.gz | sed "s|\./\.|\.:0,0:0:0:0|g" | gzip > combined_final_tagged_nohets_tagged.vcf.gz

#sep 5; the R script is filtering too strictly, so try re-running without removing heterozygyotes by switching from combined_final_tagged_nohets_tagged.vcf.gz to combined_final_tagged.vcf.gz

#now extract pairwise
for i in S1 S2 S3 S4 S5 S6 S7 S8 S11 S12 S13 S14 S15 S16 S17 S18 S19 S20 S21 S22 S23 S24
do bcftools view -s $i -c 1 -O z -o temp_hack_$i.vcf.gz combined_final_tagged_nohets.vcf.gz
done

for i in 102C11 11A2 11A5 11A6 16C32 29C8 38C31 39B10 40A20 48A4 48A6 55C1 55C4 58C25 60C3 82C21 86C22 88C20 98C11 98C12 9C22
do bcftools view -s $i -c 1 -O z -o temp_base_$i.vcf.gz combined_final_tagged_nohets.vcf.gz
done

for i in temp*.gz; do bcftools index $i; done

bcftools isec temp_hack_S1.vcf.gz  temp_base_11A2.vcf.gz   -p compare_S1
bcftools isec temp_hack_S2.vcf.gz  temp_base_11A6.vcf.gz   -p compare_S2
bcftools isec temp_hack_S3.vcf.gz  temp_base_48A6.vcf.gz   -p compare_S3
bcftools isec temp_hack_S4.vcf.gz  temp_base_55C1.vcf.gz   -p compare_S4
bcftools isec temp_hack_S5.vcf.gz  temp_base_55C4.vcf.gz   -p compare_S5
bcftools isec temp_hack_S6.vcf.gz  temp_base_60C3.vcf.gz   -p compare_S6
bcftools isec temp_hack_S7.vcf.gz  temp_base_82C21.vcf.gz  -p compare_S7
bcftools isec temp_hack_S8.vcf.gz  temp_base_86C22.vcf.gz  -p compare_S8
bcftools isec temp_hack_S11.vcf.gz temp_base_9C22.vcf.gz   -p compare_S11
bcftools isec temp_hack_S12.vcf.gz temp_base_38C31.vcf.gz  -p compare_S12
bcftools isec temp_hack_S13.vcf.gz temp_base_40A20.vcf.gz  -p compare_S13
bcftools isec temp_hack_S14.vcf.gz temp_base_48A4.vcf.gz   -p compare_S14
bcftools isec temp_hack_S15.vcf.gz temp_base_98C11.vcf.gz  -p compare_S15
bcftools isec temp_hack_S16.vcf.gz temp_base_11A5.vcf.gz   -p compare_S16
bcftools isec temp_hack_S17.vcf.gz temp_base_102C11.vcf.gz -p compare_S17
bcftools isec temp_hack_S18.vcf.gz temp_base_29C8.vcf.gz   -p compare_S18
bcftools isec temp_hack_S19.vcf.gz temp_base_60C3.vcf.gz   -p compare_S19
bcftools isec temp_hack_S20.vcf.gz temp_base_88C20.vcf.gz  -p compare_S20
bcftools isec temp_hack_S21.vcf.gz temp_base_16C32.vcf.gz  -p compare_S21
bcftools isec temp_hack_S22.vcf.gz temp_base_39B10.vcf.gz  -p compare_S22
bcftools isec temp_hack_S23.vcf.gz temp_base_98C12.vcf.gz  -p compare_S23
bcftools isec temp_hack_S24.vcf.gz temp_base_58C25.vcf.gz  -p compare_S24

#since bcftools isec prints out the header, we need to strip this to easily count the variant positions
for i in compare_*/*.vcf; do sed -i "/^#/d" $i; done

#rm temp*.gz temp*.csi

wc -l compare*/* > number.sites.analysis.tsv

Rscript num_sites_script.R
sed -i '1 i\#column V2 is the sample being compared, 0000.vcf is # variants unique to hackflex, 0001.vcf is # variants unique to baseclear, and 0002 and 0003 are common to both' number.sites.long.tsv

#make plink file for plotting pca
/mnt/LTR_userdata/auxie001/programs/plink --pca --allow-extra-chr --vcf combined_final_tagged_nohets.vcf.gz
