library("DESeq2")

# read in raw counts and metadata (pay attention to metadata's format and internal structure on column names and first column)
inputmatrix= as.matrix(read.csv("/Users/julieragsdale/Desktop/covid19/BLT-L/200729_UNC41-A00434_0134_AHNMGMDMXX_counts_tximport.csv", row.names = "X"))
head(inputmatrix)

# transfer matrix from float to integer
mode(inputmatrix) <- "integer"
head(inputmatrix)

metadt = read.csv("/Users/julieragsdale/Desktop/covid19/BLT-L/BLT-L_metadata_newdataonly.csv",row.names= "SampleID")
head(metadt)

# select data for baseline contrast
baseinput = inputmatrix[, c("EBP_COVID_26_F283_CAACACCT_TGGTAGCT_S11_L001","EBP_COVID_27_F283_CGGCTAAT_AGAACGAG_S12_L001",
                            "EBP_COVID_31_F289_ATATGCGC_CTGATCGT_S16_L001","EBP_COVID_25_F276_AAGCACTG_GTTGACCT_S10_L001",
                            "EBP_COVID_28_F285_GATAGGCT_TGTGACTG_S13_L001","EBP_COVID_29_F285_ACTCCATC_CCTTGATC_S14_L001")]  
head(baseinput)
basemeta = metadt[c("EBP_COVID_26_F283_CAACACCT_TGGTAGCT_S11_L001","EBP_COVID_27_F283_CGGCTAAT_AGAACGAG_S12_L001",
                            "EBP_COVID_31_F289_ATATGCGC_CTGATCGT_S16_L001","EBP_COVID_25_F276_AAGCACTG_GTTGACCT_S10_L001",
                            "EBP_COVID_28_F285_GATAGGCT_TGTGACTG_S13_L001","EBP_COVID_29_F285_ACTCCATC_CCTTGATC_S14_L001"), 
                  c("Model","Time_point_days")]
head(basemeta)

basemeta$Time_point_days <- as.factor(basemeta$Time_point_days)

# load data into DESeq object
# Design: design = ~ Model + infection means that deseq2 will test the effect of the infection(the last factor), controlling the effect of Model so that the algorithm returns the fold change result only from the effect of time. design = ~infection means the algorithm will return the fold change that result from infection without correcting for fold change that result from Model.
deseqdt = DESeqDataSetFromMatrix(countData = baseinput, colData = basemeta, design = ~ Time_point_days)
dds = DESeq(deseqdt)


# ouput in human readable format
res = results(dds)

# --> rename gene symbols here (see tutorial)

# create volcano plot
library("EnhancedVolcano")

EnhancedVolcano(res, lab = rownames(res), x="log2FoldChange", y= "pvalue", title = "BLT-L Day 2 vs. Naive Control", subtitle="fold change cutoff 1.5", xlim=c(-10, 10), ylim=c(0,45), colAlpha=0.7, transcriptPointSize=2, FCcutoff=1.5)

# filter and select top differentially expressed genes
# The lfc.cutoff is set to 0.58; remember that we are working with log2 fold changes so this translates to an actual fold change of 1.5 which is pretty reasonable.
# https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html

resCutoff = res[which(res$padj< 0.05),]
head(resCutoff)
resCutoff = resCutoff[which(abs(resCutoff$log2FoldChange)>0.58),]
head(resCutoff)

# sort by log2foldchange value
resCutoffOrdered = resCutoff[order(resCutoff$log2FoldChange),]
head(resCutoffOrdered)

#write.table(res,"/Users/julieragsdale/Desktop/covid19/BLT-L/BLTL_Day2_vs_NC_unfiltered_res_2.txt",sep="\t")
#write.table(resCutoffOrdered, "/Users/julieragsdale/Desktop/covid19/BLTL_Day2_vs_NC_resCutOffOrdered.txt", sep="\t")

