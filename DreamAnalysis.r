library(Biobase)
library(variancePartition)  # BiocManager::install("variancePartition")
library(openxlsx)  # BiocManager::install("openxlsx")
library('edgeR')

set.seed(12345)

# data input
meta <- read.csv("v3_GeoMx_metadata_61.csv", row.names = "workname")
count <- read.csv("v3_GeoMx_count_61.csv")
q3 <- read.csv("Q3_v3_select61_rmAll1.csv")

# sort the dataframe to make sure count and q3 have the same order of gene list
count <- count[order(count$TargetName),]
q3 <- q3[order(q3$TargetName),]

# check
rownames(q3) <- NULL
head(q3[1:5,1:2])
dim(q3)

rownames(count) <- NULL
head(count[1:5,1:2])
dim(count)

#rownames(meta) <- meta$work_name
rownames(count) <- count$TargetName
rownames(q3) <- q3$TargetName

#count <- count[, -c(1:5)]
count <- count[, -c(1)]
count <- count[,rownames(meta)]

count <- as.matrix(count)
mode(count) <- "integer"


# remove first column
q3 <- q3[,-c(1)]
q3<-q3[, rownames(meta)]

dim(count)
dim(q3)

mean<-rowMeans(q3)
neg.mean<-rowMeans(q3[rownames(q3)=="NegProbe-WTX",])

# filtering counts lower than negative probe
count<-count[mean>neg.mean,]
q3<-q3[mean>neg.mean,]

# naive samples

# naive: E835, E734, E828
# day2: E836, E133
# day6: E731, E130, E434


naive <- c("E835", "E734", "E828")
day2 <- c("E836", "E133")
day6 <- c("E731", "E130", "E434")


#basemeta <- meta[meta$mouse%in% naive,]
#basecount <- count[,rownames(basemeta)]
basemeta <- meta
basecount <- count


# Do model fit and DE analysis - full model

form <- ~ 0 + day + (1|mouse) + (1|Scan_name) + (1|cell) + (1|infection)

# Normalize data by voom with dream.

vm <- voomWithDreamWeights(DGEList(basecount), form, basemeta)

# Using contrasts to compare coefficients
L = makeContrastsDream(form, basemeta, 
                       contrasts = c("dayday2 - daynaive", 
                                     "dayday6 - daynaive"))

# "dayinfected - dayuninfected"

fitmm <- dream(vm, form, basemeta, L)
fitmm = eBayes(fitmm)


# check columns
head(coef(fitmm))
colnames((fitmm))

# Examine design matrix
head(fitmm$design, 3)

#colnames(coef(fitmm))[1] <- "Intercept" # avoid warning: Renaming (Intercept) to Intercept
#colnames(fitmm$design)[1] <- "Intercept"

#day2 - naive 

#top.table <- topTable(fitmm, sort.by = "P", n = Inf)

topTable(fitmm, coef="dayday2 - daynaive", number=5)

top.table <- topTable(fitmm, coef="dayday2 - daynaive", number=Inf, adjust.method="BH")

top.table$Gene <- rownames(top.table)

head(top.table, 5)

# length(which(top.table$adj.P.Val < 0.05)) # the number of DE genes are there

sig = top.table[which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) >= 0.585), ]
dim(sig)

write.table(top.table, file='Diana_day2_vs_naive_countNeg_voom_dream.csv', sep=',', col.names=T, row.names=F, quote=F)
write.table(sig, file='Diana_day2_vs_naive_countNeg_voom_dream_sig.csv', sep=',', col.names=T, row.names=F, quote=F)



#day6 - naive

topTable(fitmm, coef="dayday6 - daynaive", number=5)

top.table <- topTable(fitmm, coef="dayday6 - daynaive", number=Inf, adjust.method="BH")

# length(which(top.table$adj.P.Val < 0.05)) # the number of DE genes are there

top.table$Gene <- rownames(top.table)

head(top.table, 5)

sig = top.table[which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) >= 0.585), ]
dim(sig)

write.table(top.table, file='Diana_day6_vs_naive_countNeg_voom_dream.csv', sep=',', col.names=T, row.names=F, quote=F)
write.table(sig, file='Diana_day6_vs_naive_countNeg_voom_dream_sig.csv', sep=',', col.names=T, row.names=F, quote=F)






# Do model fit and DE analysis - full model

form <- ~ 0 + infection + (1|mouse) + (1|Scan_name) + (1|cell) + (1|day)

# Normalize data by voom with dream.

vm <- voomWithDreamWeights(DGEList(basecount), form, basemeta)

# Using contrasts to compare coefficients
L = makeContrastsDream(form, basemeta, 
                       contrasts = c("infectioninfected - infectionuninfected"))

# "dayinfected - dayuninfected"

fitmm <- dream(vm, form, basemeta, L)
fitmm = eBayes(fitmm)

#infected - uninfected

topTable(fitmm, coef="infectioninfected - infectionuninfected", number=5)

top.table <- topTable(fitmm, coef="infectioninfected - infectionuninfected", number=Inf, adjust.method="BH")

# length(which(top.table$adj.P.Val < 0.05)) # the number of DE genes are there

top.table$Gene <- rownames(top.table)

head(top.table, 5)

sig = top.table[which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) >= 0.585), ]
dim(sig)

write.table(top.table, file='Diana_infected_vs_uninfected_countNeg_voom_dream.csv', sep=',', col.names=T, row.names=F, quote=F)
write.table(sig, file='Diana_infected_vs_uninfected_countNeg_voom_dream_sig.csv', sep=',', col.names=T, row.names=F, quote=F)




### compare all samples:
###
###
###
###

# Do model fit and DE analysis - full model

form <- ~ 0 + cell + (1|mouse) + (1|Scan_name)

# Normalize data by voom with dream.

vm <- voomWithDreamWeights(DGEList(count), form, meta)

# Using contrasts to compare coefficients
L = makeContrastsDream(form, meta, 
                       contrasts = c("cellFoam - cellInter", 
                                     "cellFoam - cellControl",
                                     "cellInter - cellControl"))

fitmm <- dream(vm, form, meta, L)
fitmm = eBayes(fitmm)


# check columns
head(coef(fitmm))
colnames((fitmm))

# Examine design matrix
head(fitmm$design, 3)

#colnames(coef(fitmm))[1] <- "Intercept" # avoid warning: Renaming (Intercept) to Intercept
#colnames(fitmm$design)[1] <- "Intercept"

#Foamy vs Inter 

#top.table <- topTable(fitmm, sort.by = "P", n = Inf)

topTable(fitmm, coef="cellFoam - cellInter", number=5)

top.table <- topTable(fitmm, coef="cellFoam - cellInter", number=Inf, adjust.method="BH")

top.table$Gene <- rownames(top.table)

head(top.table, 5)

# length(which(top.table$adj.P.Val < 0.05)) # the number of DE genes are there

write.table(top.table, file='Garcia_v3_AllSamples_Foam_vs_Inter_countNeg_voom_dream.csv', sep=',', col.names=T, row.names=F, quote=F)

#Inter vs control

topTable(fitmm, coef="cellInter - cellControl", number=5)

top.table <- topTable(fitmm, coef="cellInter - cellControl", number=Inf, adjust.method="BH")

# length(which(top.table$adj.P.Val < 0.05)) # the number of DE genes are there

top.table$Gene <- rownames(top.table)

head(top.table, 5)

write.table(top.table, file='Garcia_v3_AllSamples_Inter_vs_Control_countNeg_voom_dream.csv', sep=',', col.names=T, row.names=F, quote=F)

#Foamy vs control

topTable(fitmm, coef="cellFoam - cellControl", number=5)

top.table <- topTable(fitmm, coef="cellFoam - cellControl", number=Inf, adjust.method="BH")

# length(which(top.table$adj.P.Val < 0.05)) # the number of DE genes are there

top.table$Gene <- rownames(top.table)

head(top.table, 5)

write.table(top.table, file='Garcia_v3_AllSamples_Foam_vs_Control_countNeg_voom_dream.csv', sep=',', col.names=T, row.names=F, quote=F)
