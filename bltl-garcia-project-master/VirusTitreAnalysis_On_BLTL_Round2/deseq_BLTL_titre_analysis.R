library("DESeq2")

# read in raw counts and metadata from Round1
inputmatrix_round1= as.matrix(read.csv("/Users/julieragsdale/Desktop/covid19/deseq2/covid19_counts_tximport.csv", row.names = "X"))
head(inputmatrix_round1)

mode(inputmatrix_round1) <- "integer"
head(inputmatrix_round1)

# read in raw counts and metadata from Round2
inputmatrix_round2= as.matrix(read.csv("/Users/julieragsdale/Desktop/covid19/BLT-L/200729_UNC41-A00434_0134_AHNMGMDMXX_counts_tximport.csv", row.names = "X"))
head(inputmatrix_round2)

# transfer matrix from float to integer
mode(inputmatrix_round2) <- "integer"
head(inputmatrix_round2)

#read in all Metadata for BLTL mice

metadt = read.csv("/Users/julieragsdale/Desktop/covid19/BLT-L/BTL_mice_titre_metadata.csv",row.names= "sampleName")
head(metadt)

# select data for baseline contrast
#baseinput = inputmatrix[, c("EBP_COVID_26_F283_CAACACCT_TGGTAGCT_S11_L001","EBP_COVID_27_F283_CGGCTAAT_AGAACGAG_S12_L001",
#                            "EBP_COVID_31_F289_ATATGCGC_CTGATCGT_S16_L001","EBP_COVID_25_F276_AAGCACTG_GTTGACCT_S10_L001",
#                            "EBP_COVID_28_F285_GATAGGCT_TGTGACTG_S13_L001","EBP_COVID_29_F285_ACTCCATC_CCTTGATC_S14_L001")]  

#baseinput = inputmatrix[, c("EBP_COVID_26_F283_CAACACCT_TGGTAGCT_S11_L001","EBP_COVID_27_F283_CGGCTAAT_AGAACGAG_S12_L001",
#                            "EBP_COVID_24_F168_GCTATCCT_CACCTGTT_S9_L001","EBP_COVID_30_F288_ATCACACG_TACGCTAC_S15_L001")]  

baseinput_round1 = inputmatrix_round1[, c("COVID_07.E744","COVID_08.E855","COVID_09.E861","COVID_18.E862","COVID_19.E862","COVID_20.E863",
                            "COVID_21.E863")] 

baseinput_round2 = inputmatrix_round2[, c("EBP_COVID_30_F288_ATCACACG_TACGCTAC_S15_L001","EBP_COVID_24_F168_GCTATCCT_CACCTGTT_S9_L001",
                                          "EBP_COVID_20_E863_AGTGGATC_GTCATCGA_S5_L001","EBP_COVID_21_E863_CAATGTGG_GTCGAAGA_S6_L001",
                                          "EBP_COVID_18_E862_ATGCCTGT_CGCAATCT_S3_L001","EBP_COVID_19_E862_GATCGAGT_AGTTCGTC_S4_L001",
                                          "EBP_COVID_26_F283_CAACACCT_TGGTAGCT_S11_L001","EBP_COVID_27_F283_CGGCTAAT_AGAACGAG_S12_L001")]  

head(baseinput_round1)
head(baseinput_round2)

library(dplyr)

#Merge counts files from Rounds 1 and 2 into object baseinput_all

baseinput_all<- merge(baseinput_round1,baseinput_round2,by=0) 
head(baseinput_all)

#Rename the rows of object baseinput_all to Ensembl IDs and make into baseinputall2 object

baseinput_all2 <- baseinput_all[,-1]
rownames(baseinput_all2) <- baseinput_all[,1]

head(baseinput_all2)


#Pull out the relevant BLTL samples' metadata and columns round, time_point_days, and estimated_titre
#Estimated_titre is raw covid reads divided by raw human reads counts 

head(metadt)

basemeta = metadt[c("COVID_07.E744","COVID_08.E855","COVID_09.E861","COVID_18.E862","COVID_19.E862","COVID_20.E863",
                    "COVID_21.E863","EBP_COVID_30_F288_ATCACACG_TACGCTAC_S15_L001","EBP_COVID_24_F168_GCTATCCT_CACCTGTT_S9_L001",
                    "EBP_COVID_20_E863_AGTGGATC_GTCATCGA_S5_L001","EBP_COVID_21_E863_CAATGTGG_GTCGAAGA_S6_L001",
                    "EBP_COVID_18_E862_ATGCCTGT_CGCAATCT_S3_L001","EBP_COVID_19_E862_GATCGAGT_AGTTCGTC_S4_L001",
                    "EBP_COVID_26_F283_CAACACCT_TGGTAGCT_S11_L001","EBP_COVID_27_F283_CGGCTAAT_AGAACGAG_S12_L001"), 
                 c("round","Time_point_days","estimated_titre")]
head(basemeta)

#Play around with R data types to satisfy DeSeq2

basemeta$round <- as.factor(basemeta$round)
basemeta$Time_point_days <- as.factor(basemeta$Time_point_days)
basemeta$estimated_titre <- as.double(basemeta$estimated_titre)

head(basemeta)

#New basemeta 2 object with two levels of virus: lots of virus and no virus 

basemeta2 <- cbind(basemeta, Virus_Presence="virus")
head(basemeta2)
basemeta2$Virus_Presence <- as.character(basemeta2$Virus_Presence)
basemeta2['COVID_09.E861','Virus_Presence'] <- "lots_of_virus"
basemeta2['EBP_COVID_30_F288_ATCACACG_TACGCTAC_S15_L001','Virus_Presence'] <- "lots_of_virus"
basemeta2['EBP_COVID_24_F168_GCTATCCT_CACCTGTT_S9_L001','Virus_Presence'] <- "lots_of_virus"

basemeta2['EBP_COVID_27_F283_CGGCTAAT_AGAACGAG_S12_L001','Virus_Presence'] <- "no_virus"
basemeta2['EBP_COVID_26_F283_CAACACCT_TGGTAGCT_S11_L001','Virus_Presence'] <- "no_virus"
basemeta2['EBP_COVID_19_E862_GATCGAGT_AGTTCGTC_S4_L001','Virus_Presence'] <- "no_virus"
basemeta2['EBP_COVID_18_E862_ATGCCTGT_CGCAATCT_S3_L001','Virus_Presence'] <- "no_virus"
basemeta2['EBP_COVID_21_E863_CAATGTGG_GTCGAAGA_S6_L001','Virus_Presence'] <- "no_virus"
basemeta2['EBP_COVID_20_E863_AGTGGATC_GTCATCGA_S5_L001','Virus_Presence'] <- "no_virus"
basemeta2$Virus_Presence <- as.factor(basemeta2$Virus_Presence)


#Trying a log transformation of estimated_titre to normalize the distribution, investigate distribution with bar plot to see if it's multimodal,
#so as to inform the bining of virus titre categories 

basemeta2$Viral_Distr <- log2(basemeta2$estimated_titre)

barplot(basemeta2$estimated_titre)
barplot(abs(basemeta2$Viral_Distr[-15]),width =0.5)

hist(abs(basemeta2$Viral_Distr[-15]),breaks=50)


#Trying another basemeta object (basemeta3) that has Three Virus Titre Categories: High, Zero, and Low 

basemeta3 <- basemeta2

basemeta3$Virus_Presence <- as.character(basemeta3$Virus_Presence)
basemeta3['COVID_09.E861','Virus_Presence'] <- "High"
basemeta3['EBP_COVID_30_F288_ATCACACG_TACGCTAC_S15_L001','Virus_Presence'] <- "High"
basemeta3['EBP_COVID_24_F168_GCTATCCT_CACCTGTT_S9_L001','Virus_Presence'] <- "High"

basemeta3[4:7,'Virus_Presence'] <- "Zero"
basemeta3[10:15,'Virus_Presence'] <- "Zero"
basemeta3[1:2,'Virus_Presence'] <- "Low"

#Trying another basemeta object (basemeta4) that has Four Virus Titre Categories: High, Zero, and Low, and Moderate 


basemeta4 <- basemeta3
basemeta4$Virus_Presence <- as.character(basemeta4$Virus_Presence)
basemeta4[8:9,'Virus_Presence'] <- "Moderate"
basemeta4$Virus_Presence <- as.factor(basemeta4$Virus_Presence)

basemeta4$combined_time_virus <- with(basemeta4,paste0(Time_point_days,Virus_Presence))

#Try another way to normalize distribution of estimated virus titre 

basemeta4$normalized_virus_titre <- normalize(basemeta4$estimated_titre, method = "standardize")


# load data into DESeq object
# Design: design = ~ Model + infection means that deseq2 will test the effect of the infection(the last factor), controlling the effect of Model so that the algorithm returns the fold change result only from the effect of time. design = ~infection means the algorithm will return the fold change that result from infection without correcting for fold change that result from Model.

#Several Different Stats Models that were tried:


#WITH the basemeta object with two virus levels: High and Low

#deseqdt = DESeqDataSetFromMatrix(countData = baseinput_all2, colData = basemeta2, design = ~ Time_point_days + round + Virus_Presence)

#WITH basemeta object with four levels of virus titre: High, Low, Zero, and Moderate

#deseqdt = DESeqDataSetFromMatrix(countData = baseinput_all2, colData = basemeta4, design = ~ Time_point_days + round + Virus_Presence )

#deseqdt = DESeqDataSetFromMatrix(countData = baseinput_all2, colData = basemeta2, design = ~ Time_point_days + round + Viral_Distr)


#THIS WAS THE FINAL RESULTING STATS MODEL !!!! that was used to avoid linear variables, time point + virus level were combined

#deseqdt = DESeqDataSetFromMatrix(countData = baseinput_all2, colData = basemeta4, design = ~ combined_time_virus + round )

#deseqdt = DESeqDataSetFromMatrix(countData = baseinput_all2, colData = basemeta4, design = ~ Time_point_days + round + normalized_virus_titre)


dds = DESeq(deseqdt)

# Testing Different Contrasts:

#res2=results(dds, contrast=c("round","stage1treatmentX", "stage2treatmentX"))
#res3=results(dds, contrast=c("round","1", "2"))

#head(results(dds, contrast=c("Time_point_days","6", "Naive_Control")))

#res2=results(dds, contrast=c("combined_time_virus","Naive_ControlZero","6Low"))
#res2=results(dds, contrast=c("combined_time_virus","6Low","6Moderate"))
#res2=results(dds, contrast=c("combined_time_virus","6Moderate","6High"))

res2=results(dds, contrast=c("combined_time_virus","Naive_ControlZero","6High"))
res2Cutoff = res2[which(res2$padj< 0.05),]
head(res2Cutoff)
res2Cutoff = resCutoff[which(abs(res2Cutoff$log2FoldChange)>0.58),]
head(res2Cutoff)


# Create volcano plot

library("EnhancedVolcano")

EnhancedVolcano(res2, lab = rownames(res2), x="log2FoldChange", y= "pvalue", title = "BLT-L Day 6 Low Virus vs. Moderate Virus Titre", subtitle="fold change cutoff 1.5", xlim=c(-10, 10), ylim=c(0,45), colAlpha=0.7, transcriptPointSize=2, FCcutoff=1.5)



#write.table(resCutoffOrdered, "/Users/julieragsdale/Desktop/covid19/BLT-L/BLTL_Day6_vs_NC_CombinedTimeVirusFactor_resCutOffOrdered_092120.txt", sep="\t")

write.table(res2Cutoff,"/Users/julieragsdale/Desktop/covid19/BLT-L/Deseq2_VirusTitre_Results/BLTL_Day6vsNC_contrast_NCvsHighVirusTitre_resCutOff_0922.txt",sep="\t")
