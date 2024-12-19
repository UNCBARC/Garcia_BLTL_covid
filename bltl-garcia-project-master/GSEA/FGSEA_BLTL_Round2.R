library(tidyverse) 

res <- read_csv("/Users/julieragsdale/Desktop/covid19/diff_expr_gene_list_baselineContrastp0.05lfc0.58_with_hgnc.csv") 

#res <- read.csv("/Users/julieragsdale/Desktop/covid19/BLT-L/BLTL_Day2_vs_NC_resCutOffOrdered_w_GeneSymbols.txt",sep="\t")

res2 <- res %>% 
  dplyr::select(GeneSymbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(GeneSymbol) %>% 
  summarize(stat=mean(stat))

head(res2)

library(fgsea)

ranks <- deframe(res2)
head(ranks, 20)

pathways.hallmark <- gmtPathways("/Users/julieragsdale/Desktop/covid19/gsea/c5.all.v7.1.symbols.gmt")

pathways.hallmark %>% 
  head() %>% 
  lapply(head)

inflammatory.pathways.hallmark <- pathways.hallmark$GO_INFLAMMATORY_RESPONSE

#fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=100000)
#fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

#fgseaResTidy %>% 
#  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
#  arrange(padj) %>% 
#  DT::datatable()

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

fgseaResTidy$leadingEdge <- sapply(fgseaResTidy$leadingEdge, paste, collapse=",") 

#fgseaResTidy2<- p05_fgseaResTidy[order(p05_fgseaResTidy$pval),]

p.unadjusted <- fgseaResTidy$pval
p_adjusted <- p.adjust(p.unadjusted, method = "BH")

write.table(fgseaResTidy,"/Users/julieragsdale/Desktop/covid19/LoM_Day2_all_FGSEA_results_nperm_100thl.txt",sep="\t")

inflammatory_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("INFLAMMATORY",pathway))

#PADJ cutoff Inflammatory GSEA Plot 
ggplot(inflammatory_fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Inflammatory pathways NES from GSEA") + 
  theme_minimal()
  
viral_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("VIR",pathway))

#P-VALUE cutoff Inflammatory/Viral GSEA Plot 
ggplot(inflammatory_fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Viral GOs NES from GSEA") + 
  scale_fill_manual(values=c("#0000FF","#FF0000")) +
  theme_minimal()

#custom_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("GO_CELLULAR_RESPONSE_TO_INTERFERON",pathway))

#custom1_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("GO_RESPONSE_TO_TYPE_I_INTERFERON",pathway))

#custom2_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("GO_CYTOKINE_SECRETION",pathway))[1,]

#custom3_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("GO_ACTIVATION_OF_INNATE_IMMUNE_RESPONSE",pathway))

#custom4_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("GO_INNATE_IMMUNE_RESPONSE",pathway))[1,]

#custom5_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("GO_ACUTE_INFLAMMATORY_RESPONSE",pathway))[1,]

#custom6_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("GO_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE",pathway))

#custom7_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("GO_CYTOKINE_ACTIVITY",pathway))



custom2_1_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("GO_INNATE_IMMUNE_RESPONSE",pathway))[1,]

custom2_2_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("CYTOKINE_MEDIATED_SIGNALING_PATHWAY",pathway))

custom2_3_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("GO_REGULATION_OF_RESPONSE_TO_STRESS",pathway))[1,]

custom2_4_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("INFLAMMATORY_RESPONSE",pathway))[1,]

custom2_5_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("CYTOKINE_PRODUCTION",pathway))[2,]

custom2_6_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("RESPONSE_TO_VIRUS",pathway))[1,]

custom2_7_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("RESPONSE_TO_TYPE_I_INTERFERON",pathway))

custom2_8_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("GO_NIK_NF_KAPPAB_SIGNALING",pathway))

custom2_9_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("GO_REGULATION_OF_CELL_DEATH",pathway))

custom2_10_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("ACUTE_INFLAMMATORY_RESPONSE",pathway))[1,]

custom2_11_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("GO_COAGULATION",pathway))

custom2_12_fgseaResTidy<- dplyr::filter(fgseaResTidy, grepl("GO_COMPLEMENT_ACTIVATION",pathway))[1,]




custom2_fgseaResTidy <-rbind(custom2_1_fgseaResTidy,custom2_2_fgseaResTidy,custom2_3_fgseaResTidy,custom2_4_fgseaResTidy,
                             custom2_5_fgseaResTidy,custom2_6_fgseaResTidy,custom2_7_fgseaResTidy,custom2_8_fgseaResTidy,
                             custom2_9_fgseaResTidy,custom2_10_fgseaResTidy,custom2_11_fgseaResTidy,custom2_12_fgseaResTidy)



ggplot(custom2_fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="NES from GSEA") + 
  scale_fill_manual(values=c("#0000FF","#FF0000")) +
  theme_minimal()

ggplot(custom2_fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="NES from GSEA") + 
  scale_fill_manual(values=c("#FF8106","#FF0000")) +
  theme_minimal()

#p05_fgseaResTidy <-fgseaResTidy[apply(fgseaResTidy < 0.05, 1, all), ]

p05_fgseaResTidy <- filter(fgseaResTidy, padj < 0.05)

#top_p05_fgseaResTidy <- head(p05_fgseaResTidy,100)

ggplot(top_p05_fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="All p<0.05 GO pathways NES from GSEA") 
  theme_minimal()

p05_fgseaResTidy$leadingEdge <- sapply(p05_fgseaResTidy$leadingEdge, paste, collapse=",") 

#p05_fgseaResTidy2<- p05_fgseaResTidy[order(p05_fgseaResTidy$pval),]
p05_fgseaResTidy<- p05_fgseaResTidy[order(p05_fgseaResTidy$padj),]

#write.table(fgseaResTidy,"/Users/julieragsdale/Desktop/covid19/BLT-L/BLT_Round2_Day2_all_FGSEA_results_1.txt",sep="\t")

write.table(p05_fgseaResTidy,file="/Users/julieragsdale/Desktop/covid19/LoM_Ro_Day2_padj05_FGSEA_results_nperm_100th.txt",sep="\t",col.names = TRUE)

