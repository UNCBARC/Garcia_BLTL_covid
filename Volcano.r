#module add r/4.1.0
library(DESeq2)
library("EnhancedVolcano")
library(org.Hs.eg.db)

metadt = read.csv("victor_deseq_meta.csv", row.names = "sampleName")
inputmatrix= as.matrix(read.csv("batch_corr_victor_2020_counts_tximport.csv", row.names = "geneid"))

mode(inputmatrix) <- "integer"

basemeta = metadt[metadt$Model == "BLTL", ]
basemeta$type = c("E","E","E","E","E","E","E","E","E","E","E","E","E","E","F","F","F","F","F","F","F","F")
metadt = basemeta

############################################################
#Day2:
basemeta = metadt[(metadt$days == "Control" | metadt$days == "2") & metadt$FCs == "FC2", ]
baseinput = inputmatrix[, rownames(basemeta)]

# match gene name to ensemble id

ens <- rownames(baseinput)

symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(baseinput), names(symbols))]
rownames(baseinput) <- symbols
keep <- !is.na(rownames(baseinput))
baseinput <- baseinput[keep,]

deseqdt = DESeqDataSetFromMatrix(countData = baseinput, colData = basemeta, design = ~type + days)

dds = DESeq(deseqdt)
resultsNames(dds)
res = results(dds, contrast=c("days", "2", "Control"))
#res <- lfcShrink(dds, contrast=c("days", "2", "Control"), type="ashr")

sig = res[which(res$padj < 0.05 & abs(res$log2FoldChange) >= 0.585), ]
#sig = res[which(res$padj < 0.05 & abs(res$log2FoldChange) >= 1), ]

#write.table(res, "BLTL_BatchCorr_day2vsControl_FC2_designEFdays.csv", quote=F, sep=",", col.names=NA)
#write.table(sig, "BLTL_BatchCorr_day2vsControl_FC2_designEFdays_sig.csv", quote=F, sep=",", col.names=NA)

# get volcano plot

# set the base colour as 'black'
keyvals <- rep('grey', nrow(res))

# set the base name/label as 'NS'
names(keyvals) <- rep('NS', nrow(res))

# modify keyvals for variables with fold change > 1.5
keyvals[which(res$log2FoldChange > 0.585 & res$padj < 0.05)] <- '#E50000'
names(keyvals)[which(res$log2FoldChange > 0.585 & res$padj < 0.05)] <- 'Up'

# modify keyvals for variables with fold change < -1.5
keyvals[which(res$log2FoldChange < -0.585 & res$padj < 0.05)] <- 'blue'
names(keyvals)[which(res$log2FoldChange < -0.585 & res$padj < 0.05)] <- 'Down'

# Selected genes to represent
#s_genes <- c('IGHA2', 'CLDN22')
#selectLab = s_genes,
# , ylim = c(0, 10)

EnhancedVolcano(res, #xlim = c(-10, 10), ylim = c(0, 20), 
                colCustom = keyvals,
                lab = rownames(res),  
                #lab = NA, 
                pointSize = 4.5,
                x="log2FoldChange", y= "padj", title = "BLTL - Day2 vs Control", subtitle = "Fold change cutoff 1.5, padj value cutoff 0.05", 
                pCutoff = 0.05, FCcutoff =0.585, cutoffLineType = 'longdash', cutoffLineCol = 'black') 



####################################################################
# Day6
basemeta = metadt[(metadt$days == "Control" | metadt$days == "6")  & metadt$Model == "BLTL", ]
baseinput = inputmatrix[, rownames(basemeta)]

# match gene name to ensemble id

ens <- rownames(baseinput)

symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(baseinput), names(symbols))]
rownames(baseinput) <- symbols
keep <- !is.na(rownames(baseinput))
baseinput <- baseinput[keep,]

deseqdt = DESeqDataSetFromMatrix(countData = baseinput, colData = basemeta, design = ~FCs + type + days)

dds = DESeq(deseqdt)
res = results(dds, contrast=c("days", "6", "Control"))
#res = lfcShrink(dds=dds, contrast=c("days", "6", "Control"), type="ashr")

sig = res[which(res$padj < 0.05 & abs(res$log2FoldChange) >= 0.585), ]

#write.table(res, "BLTL_BatchCorr_day6vsControl_FC1FC2_designFCsEFdays.csv", quote=F, sep=",", col.names=NA)
#write.table(sig, "BLTL_BatchCorr_day6vsControl_FC1FC2_designFCsEFdays_sig.csv", quote=F, sep=",", col.names=NA)

# get volcano plot

# set the base colour as 'black'
keyvals <- rep('grey', nrow(res))

# set the base name/label as 'NS'
names(keyvals) <- rep('NS', nrow(res))

# modify keyvals for variables with fold change > 1.5
keyvals[which(res$log2FoldChange > 0.585 & res$padj < 0.05)] <- '#E50000'
names(keyvals)[which(res$log2FoldChange > 0.585 & res$padj < 0.05)] <- 'Up'

# modify keyvals for variables with fold change < -1.5
keyvals[which(res$log2FoldChange < -0.585 & res$padj < 0.05)] <- 'blue'
names(keyvals)[which(res$log2FoldChange < -0.585 & res$padj < 0.05)] <- 'Down'

# Selected genes to represent
#s_genes <- c('IGHA2', 'CLDN22')
#selectLab = s_genes,
#, ylim = c(0, 10)

EnhancedVolcano(res, #xlim = c(-10, 10), ylim = c(0, 20), 
                colCustom = keyvals,
                lab = rownames(res),  
                #lab = NA, 
                pointSize = 4.5,
                x="log2FoldChange", y= "padj", title = "BLTL - Day6 vs Control", subtitle = "Fold change cutoff 1.5, p value cutoff 0.05", 
                pCutoff = 0.05, FCcutoff =0.585, cutoffLineType = 'longdash', cutoffLineCol = 'black')




###############################################################
# Day14
basemeta = metadt[(metadt$days == "Control" | metadt$days == "14") & metadt$FCs == "FC1", ]
baseinput = inputmatrix[, rownames(basemeta)]

# match gene name to ensemble id

ens <- rownames(baseinput)

symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(baseinput), names(symbols))]
rownames(baseinput) <- symbols
keep <- !is.na(rownames(baseinput))
baseinput <- baseinput[keep,]

deseqdt = DESeqDataSetFromMatrix(countData = baseinput, colData = basemeta, design = ~days)

dds = DESeq(deseqdt)
res = results(dds, contrast=c("days", "14", "Control"))
#res = lfcShrink(dds=dds, contrast=c("days", "14", "Control"), type="ashr")
summary(res)

sig = res[which(res$padj < 0.05 & abs(res$log2FoldChange) >= 0.585), ]

#write.table(res, "BLTL_BatchCorr_day14vsControl_FC1_designdays.csv", quote=F, sep=",", col.names=NA)
#write.table(sig, "BLTL_BatchCorr_day14vsControl_FC1_designdays_sig.csv", quote=F, sep=",", col.names=NA)

# get volcano plot


# set the base colour as 'black'
keyvals <- rep('grey', nrow(res))

# set the base name/label as 'NS'
names(keyvals) <- rep('NS', nrow(res))

# modify keyvals for variables with fold change > 1.5
keyvals[which(res$log2FoldChange > 0.585 & res$padj < 0.05)] <- '#E50000'
names(keyvals)[which(res$log2FoldChange > 0.585 & res$padj < 0.05)] <- 'Up'

# modify keyvals for variables with fold change < -1.5
keyvals[which(res$log2FoldChange < -0.585 & res$padj < 0.05)] <- 'blue'
names(keyvals)[which(res$log2FoldChange < -0.585 & res$padj < 0.05)] <- 'Down'

# Selected genes to represent
#s_genes <- c('IGHA2', 'CLDN22')
#selectLab = s_genes,
#, ylim = c(0, 10)

EnhancedVolcano(res, #xlim = c(-10, 10), 
                #ylim = c(0, 20), 
                colCustom = keyvals,
                lab = rownames(res), 
                #lab = NA, 
                pointSize = 4.5,
                x="log2FoldChange", y= "padj", title = "BLTL - Day14 vs Control", subtitle = "Fold change cutoff 1.5, p value cutoff 0.05", 
                pCutoff = 0.05, FCcutoff =0.585, cutoffLineType = 'longdash', cutoffLineCol = 'black')

