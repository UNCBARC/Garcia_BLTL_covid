library(tidyverse) 
library(fgsea)
# to install fgsea:
# library(devtools)
# install_github("ctlab/fgsea")
# FGSEA tutorial: https://stephenturner.github.io/deseq-to-fgsea/

# downloaded from http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
pathways.hallmark <- gmtPathways("c5.all.v2023.2.Hs.symbols.gmt") 


# Day2 vs naive:
res <- read_csv("Diana_day2_vs_naive_countNeg_voom_dream.csv")
#res <- read_csv("Diana_day2_vs_naive_countNeg_voom_dream_sig.csv")

res2 <- res %>% 
  dplyr::select(Gene, t) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(Gene) %>% 
  summarize(stat=mean(t))

ranks <- deframe(res2)

#fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=100000) #fgseaSimple
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks) #fgseaMultilevel

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy$leadingEdge <- sapply(fgseaResTidy$leadingEdge, paste, collapse=",") 

write.table(fgseaResTidy, "Diana_day2_vs_naive_countNeg_voom_dream_fgseaMultilevel.txt", sep="\t", col.names=NA, quote=F)





# Day6 vs naive:

res <- read_csv("Diana_day6_vs_naive_countNeg_voom_dream.csv")

res2 <- res %>% 
  dplyr::select(Gene, t) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(Gene) %>% 
  summarize(stat=mean(t))

ranks <- deframe(res2)

#fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=100000) #fgseaSimple
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks) #fgseaMultilevel

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy$leadingEdge <- sapply(fgseaResTidy$leadingEdge, paste, collapse=",") 

write.table(fgseaResTidy, "Diana_day6_vs_naive_countNeg_voom_dream_fgseaMultilevel.txt", sep="\t", col.names=NA, quote=F)





# infected vs uninfected:
res <- read_csv("Diana_infected_vs_uninfected_countNeg_voom_dream.csv")

res2 <- res %>% 
  dplyr::select(Gene, t) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(Gene) %>% 
  summarize(stat=mean(t))

ranks <- deframe(res2)

#fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=100000) #fgseaSimple
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks) #fgseaMultilevel

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy$leadingEdge <- sapply(fgseaResTidy$leadingEdge, paste, collapse=",") 

write.table(fgseaResTidy, "Diana_infected_vs_uninfected_countNeg_voom_dream_fgseaMultilevel.txt", sep="\t", col.names=NA, quote=F)



