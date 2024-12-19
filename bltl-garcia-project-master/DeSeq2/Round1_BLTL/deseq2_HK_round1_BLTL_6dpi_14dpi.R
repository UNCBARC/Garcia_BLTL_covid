library("DESeq2")

# read in raw counts and metadata (pay attention to metadata's format and internal structure on column names and first column)
inputmatrix= as.matrix(read.csv("/Users/julieragsdale/Desktop/covid19/deseq2/covid19_counts_tximport.csv", row.names = "X"))
head(inputmatrix)

# transfer matrix from float to integer
mode(inputmatrix) <- "integer"
head(inputmatrix)

metadt = read.csv("/Users/julieragsdale/Desktop/covid19/deseq2/colmeta.csv", row.names = "X")
head(metadt)

# select data for baseline contrast
#baseinput = inputmatrix[, c("COVID_07.E744", "COVID_08.E855", "COVID_09.E861", "COVID_18.E862","COVID_19.E862","COVID_20.E863","COVID_21.E863")] #for 2dpi BLT
#baseinput = inputmatrix[, c("COVID_04.E130","COVID_05.E434","COVID_06.E829","COVID_22.E128","COVID_23.E838","COVID_16.E840","COVID_17.E844")] #6dpi vs. NC LoM 
baseinput = inputmatrix[, c("COVID_07.E744","COVID_08.E855","COVID_09.E861","COVID_18.E862","COVID_19.E862","COVID_20.E863","COVID_21.E863")] #Round 1 BLTL 6dpi vs NC

baseinput = inputmatrix[, c("COVID_13.E856","COVID_14.E858","COVID_15.E860","COVID_18.E862","COVID_19.E862","COVID_20.E863","COVID_21.E863")] #Round 1 BLTL 14dpi vs. NC

head(baseinput)
basemeta = metadt[c("COVID_13.E856","COVID_14.E858","COVID_15.E860","COVID_18.E862","COVID_19.E862","COVID_20.E863","COVID_21.E863"), c("Model","days","infection")]
head(basemeta)


# load data into DESeq object
# Design: design = ~ Model + infection means that deseq2 will test the effect of the infection(the last factor), controlling the effect of Model so that the algorithm returns the fold change result only from the effect of time. design = ~infection means the algorithm will return the fold change that result from infection without correcting for fold change that result from Model.
deseqdt = DESeqDataSetFromMatrix(countData = baseinput, colData = basemeta, design = ~infection)
dds = DESeq(deseqdt)


# ouput in human readable format
res = results(dds)

# --> rename gene symbols here (see tutorial)



# create volcano plot
library("EnhancedVolcano")

EnhancedVolcano(res, lab = rownames(res), x="log2FoldChange", y= "pvalue", title = "Round 1 BLTL Day 14", subtitle="fold change cutoff 1.5", xlim=c(-10, 10), ylim=c(0,45), colAlpha=0.7, transcriptPointSize=2, FCcutoff=1.5)

# filter and select top differentially expressed genes
# The lfc.cutoff is set to 0.58; remember that we are working with log2 fold changes so this translates to an actual fold change of 1.5 which is pretty reasonable.
# https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html

write.table(res,"/Users/julieragsdale/Desktop/covid19/BLT-L/BLTL_Round1_14dpi_v_NaiveControl_Deseq2_ALL_results.txt",sep="\t")

resCutoff = res[which(res$padj< 0.05),]
head(resCutoff)
resCutoff = resCutoff[which(abs(resCutoff$log2FoldChange)>0.58),]
head(resCutoff)

# sort by log2foldchange value
resCutoffOrdered = resCutoff[order(resCutoff$log2FoldChange),]
head(resCutoffOrdered)

write.table(resCutoffOrdered, "/Users/julieragsdale/Desktop/covid19/BLT-L/BLTL_Round1_14dpi_v_NaiveControl_Deseq2_res_cutofforder.txt", sep="\t")

