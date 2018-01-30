# install.packages("feather")
library(feather)
library(tidyverse)

apobec_hnscc_df <- read_feather("~/Documents/manuscripts/APOBEC_HNSCC/nfkb.feather")

head(apobec_hnscc_df)


apobec_hnscc_df %>%
  group_by(categ) %>%
  summarize(mean(S2),mean(S13))

# APOBEC_expression <- read.table(file = "expression/hnsc_tcga_pub_hnsc_tcga_pub_rna_seq_v2_mrna_median_Zscores_we_want.txt",skip=1,header = T,stringsAsFactors = F,sep=" ")
# downloaded from http://www.cbioportal.org/index.do?session_id=5a0f48fa498e5df2e2986845&show_samples=false&
APOBEC_expression <- read.table(file = "expression/hnsc_tcga_hnsc_tcga_rna_seq_v2_mrna_median_Zscores_tab.txt",header = T,stringsAsFactors = F,sep="\t")


load("HNSC_MAF.RData")
HNSC.MAF <- HNSC.MAF[which(HNSC.MAF$Variant_Type=="SNP"),]

head(HNSC.MAF)
length(unique(HNSC.MAF$Unique_patient_identifier))
unique(names(APOBEC_expression))
unique(HNSC.MAF$Unique_patient_identifier)

names(APOBEC_expression)

# Need to only keep names that correspond to APOBEC

APOBEC_expression <- APOBEC_expression[startsWith(x = names(APOBEC_expression),prefix = "TCGA")]
APOBEC_expression <- APOBEC_expression[endsWith(x = names(APOBEC_expression),suffix = ".01")]

# Need to replace all . with - and delete last 3

names(APOBEC_expression) <- gsub("\\.", "-", names(APOBEC_expression))

first.12 <- function(string_to_12){
  return(paste(unlist(strsplit(string_to_12,split = ""))[1:12],collapse = ""))
}

names(APOBEC_expression) <- unlist(lapply(names(APOBEC_expression),first.12))


APOBEC3B.expression.and.hpv <- as.data.frame(matrix(data = NA,nrow = length(unique(HNSC.MAF$Unique_patient_identifier)),ncol=3))
rownames(APOBEC3B.expression.and.hpv) <- unique(HNSC.MAF$Unique_patient_identifier)
colnames(APOBEC3B.expression.and.hpv) <- c("HPV_status","APOBEC3B_expression_Z_score","SNV_count")

for(i in 1:nrow(APOBEC3B.expression.and.hpv)){
  APOBEC3B.expression.and.hpv[i,"SNV_count"] <- length(which(HNSC.MAF$Unique_patient_identifier==rownames(APOBEC3B.expression.and.hpv)[i]))
  APOBEC3B.expression.and.hpv[i,"HPV_status"] <- apobec_hnscc_df$categ[which(apobec_hnscc_df$patient_id==rownames(APOBEC3B.expression.and.hpv)[i])]
  if(length(which(names(APOBEC_expression)==rownames(APOBEC3B.expression.and.hpv)[i]))>0){
    APOBEC3B.expression.and.hpv[i,"APOBEC3B_expression_Z_score"] <- as.numeric(APOBEC_expression[which(names(APOBEC_expression)==rownames(APOBEC3B.expression.and.hpv)[i])])
  }
}


APOBEC.data <- subset(APOBEC3B.expression.and.hpv, HPV_status!="ambiguous")

library(ggplot2)

APOBEC.plot <- ggplot(data = APOBEC.data, aes(x=APOBEC3B_expression_Z_score,y=SNV_count)) + 
  geom_point(aes(x=APOBEC3B_expression_Z_score,y=SNV_count)) + 
  geom_smooth(method=lm) + 
  theme_bw() +
  facet_grid(HPV_status ~ .) 

APOBEC.plot
ggsave(filename = "Figures/SNV_vs_APOBEC3B_expression.png")

APOBEC.plot.subset <- ggplot(data = subset(APOBEC.data,SNV_count<2000), aes(x=APOBEC3B_expression_Z_score,y=SNV_count)) + 
  geom_point(aes(x=APOBEC3B_expression_Z_score,y=SNV_count)) + 
  geom_smooth(method=lm) + 
  theme_bw() +
  facet_grid(HPV_status ~ .) 

APOBEC.plot.subset
ggsave(filename = "Figures/subset_SNV_vs_APOBEC3B_expression.png",plot = APOBEC.plot.subset)


summary(lm(APOBEC.data$SNV_count ~ APOBEC.data$APOBEC3B_expression_Z_score))

summary(lm(APOBEC.data$SNV_count[which(APOBEC.data$HPV_status=="HPV+")] ~ APOBEC.data$APOBEC3B_expression_Z_score[which(APOBEC.data$HPV_status=="HPV+")]))

summary(lm(APOBEC.data$SNV_count[which(APOBEC.data$HPV_status=="HPV-")] ~ APOBEC.data$APOBEC3B_expression_Z_score[which(APOBEC.data$HPV_status=="HPV-")]))


summary(lm(APOBEC.data$SNV_count[which(APOBEC.data$HPV_status=="HPV-" & APOBEC.data$SNV_count<2000)] ~ APOBEC.data$APOBEC3B_expression_Z_score[which(APOBEC.data$HPV_status=="HPV-"  & APOBEC.data$SNV_count<2000)]))





APOBEC.plot <- ggplot(data = APOBEC.data, aes(y=APOBEC3B_expression_Z_score,x=SNV_count)) + 
  geom_point(aes(y=APOBEC3B_expression_Z_score,x=SNV_count)) + 
  geom_smooth(method=lm) + 
  theme_bw() +
  facet_grid(HPV_status ~ .) 

APOBEC.plot

APOBEC.plot.subset <- ggplot(data = subset(APOBEC.data,SNV_count<2000), aes(x=APOBEC3B_expression_Z_score,x=SNV_count)) + 
  geom_point(aes(y=APOBEC3B_expression_Z_score,x=SNV_count)) + 
  geom_smooth(method=lm) + 
  theme_bw() +
  facet_grid(HPV_status ~ .) 

APOBEC.plot.subset


summary(lm(APOBEC.data$APOBEC3B_expression_Z_score ~APOBEC.data$SNV_count ))

summary(lm( APOBEC.data$APOBEC3B_expression_Z_score[which(APOBEC.data$HPV_status=="HPV+")]~APOBEC.data$SNV_count[which(APOBEC.data$HPV_status=="HPV+")] ))

summary(lm(APOBEC.data$APOBEC3B_expression_Z_score[which(APOBEC.data$HPV_status=="HPV-")] ~ APOBEC.data$SNV_count[which(APOBEC.data$HPV_status=="HPV-")] ))


summary(lm( APOBEC.data$APOBEC3B_expression_Z_score[which(APOBEC.data$HPV_status=="HPV-"  & APOBEC.data$SNV_count<2000)]~APOBEC.data$SNV_count[which(APOBEC.data$HPV_status=="HPV-" & APOBEC.data$SNV_count<2000)] ))





# SNV count vs. APOBEC  vs. HPV+ 

APOBEC.data$APOBEC_pos_neg <- NA
for(i in 1:nrow(APOBEC.data)){
  APOBEC.data$APOBEC_pos_neg[i] <- apobec_hnscc_df$has_apobec[which(apobec_hnscc_df$patient_id==rownames(APOBEC.data)[i])]
}

require(ggplot2)
ggplot(data=APOBEC.data,aes(y=log10(SNV_count), x=APOBEC_pos_neg))+ geom_boxplot() + geom_jitter(alpha=0.4,width = 0.1)  + facet_grid(HPV_status ~.) + coord_flip() + labs(x="APOBEC signal",y="log10(SNV count)")

APOBEC.data %>% group_by(HPV_status, APOBEC_pos_neg) %>% count()

hpv.pos.apo.pos <- APOBEC.data$SNV_count[which(APOBEC.data$HPV_status=="HPV+" & APOBEC.data$APOBEC_pos_neg==T)]
hpv.pos.apo.neg <- APOBEC.data$SNV_count[which(APOBEC.data$HPV_status=="HPV+" & APOBEC.data$APOBEC_pos_neg==F)]

wilcox.test(hpv.pos.apo.neg,hpv.pos.apo.pos,paired = F,conf.int = T)
wilcox.test(hpv.pos.apo.pos,hpv.pos.apo.neg,paired = F,conf.int = T)
median(hpv.pos.apo.neg);median(hpv.pos.apo.pos)

t.test(hpv.pos.apo.neg,hpv.pos.apo.pos)

hpv.neg.apo.pos <- APOBEC.data$SNV_count[which(APOBEC.data$HPV_status=="HPV-" & APOBEC.data$APOBEC_pos_neg==T)]
hpv.neg.apo.neg <- APOBEC.data$SNV_count[which(APOBEC.data$HPV_status=="HPV-" & APOBEC.data$APOBEC_pos_neg==F)]

wilcox.test(hpv.neg.apo.neg,hpv.neg.apo.pos,paired = F,conf.int = T)
wilcox.test(hpv.neg.apo.pos,hpv.neg.apo.neg,paired = F,conf.int = T)

t.test(hpv.neg.apo.neg,hpv.neg.apo.pos)


median(hpv.neg.apo.neg);median(hpv.neg.apo.pos)

# wilcox.test(log10(hpv.pos.apo.neg),log10(hpv.pos.apo.pos))

# wilcox.test(log10(hpv.neg.apo.neg),log10(hpv.neg.apo.pos))
# wilcox.test(hpv.neg.apo.pos,hpv.neg.apo.neg)

wilcox.test(hpv.neg.apo.pos,hpv.pos.apo.pos)
wilcox.test(hpv.pos.apo.neg,hpv.neg.apo.neg)

#Testing 
# mean(apobec_hnscc_df$S2[which(apobec_hnscc_df$categ=="HPV+")])


hpv.neg <- APOBEC.data$SNV_count[which(APOBEC.data$HPV_status=="HPV-")]

wilcox.test(hpv.pos.apo.pos,hpv.neg,paired = F,conf.int = T)
# ggplot(data = APOBEC.data,aes(y=log10(SNV_count),x=APOBEC_pos_neg))


wilcox.test(hpv.pos.apo.neg,hpv.neg,paired = F,conf.int = T)






# Need to catalogue all mutations as APOBEC mutations or not. 

load("~/Documents/Selection_analysis/HNSC/selection_output/HNSC_selection_output.RData")

genes <- unique(selection.output$complete_mutation_data$Gene)

nrow(selection.output$complete_mutation_data)

selection.output$complete_mutation_data$TCW_TKW <- NA

# weight.and.trinuc.df <- as.data.frame(matrix(nrow=nrow(selection.output$complete_mutation_data),ncol=5,data=NA))
# colnames(weight.and.trinuc.df) <- c("gene","mutation","tumor","APOBEC_weight","TCW_TKW")

if(length(which(is.na(selection.output$complete_mutation_data$Nucleotide_trinuc_context)))>0){
  selection.output$complete_mutation_data <- selection.output$complete_mutation_data[-which(is.na(selection.output$complete_mutation_data$Nucleotide_trinuc_context)),] 
}


for(i in 1:nrow(selection.output$complete_mutation_data)){
  if(((selection.output$complete_mutation_data$Nucleotide_trinuc_context[i]=="TCA" & selection.output$complete_mutation_data$Alternative_Nucleotide[i]=="T") | 
      (selection.output$complete_mutation_data$Nucleotide_trinuc_context[i]=="TGA" & selection.output$complete_mutation_data$Alternative_Nucleotide[i]=="A")) | 
     ((selection.output$complete_mutation_data$Nucleotide_trinuc_context[i]=="TCT" & selection.output$complete_mutation_data$Alternative_Nucleotide[i]=="T") | 
      (selection.output$complete_mutation_data$Nucleotide_trinuc_context[i]=="AGA" & selection.output$complete_mutation_data$Alternative_Nucleotide[i]=="A")) |
     ((selection.output$complete_mutation_data$Nucleotide_trinuc_context[i]=="TCA" & selection.output$complete_mutation_data$Alternative_Nucleotide[i]=="G") | 
      (selection.output$complete_mutation_data$Nucleotide_trinuc_context[i]=="TGA" & selection.output$complete_mutation_data$Alternative_Nucleotide[i]=="C")) |
     ((selection.output$complete_mutation_data$Nucleotide_trinuc_context[i]=="TCT" & selection.output$complete_mutation_data$Alternative_Nucleotide[i]=="G") | 
      (selection.output$complete_mutation_data$Nucleotide_trinuc_context[i]=="AGA" & selection.output$complete_mutation_data$Alternative_Nucleotide[i]=="C"))){
    selection.output$complete_mutation_data$TCW_TKW[i] <- 1
  }else{
    selection.output$complete_mutation_data$TCW_TKW[i] <- 0
  }
}



APOBEC.data$percent_TCW_TKW <- NA

for(i in 1:nrow(APOBEC.data)){
  APOBEC.data$percent_TCW_TKW[i] <- sum(selection.output$complete_mutation_data$TCW_TKW[which(selection.output$complete_mutation_data$Unique_patient_identifier==rownames(APOBEC.data)[i])])/length(which(selection.output$complete_mutation_data$Unique_patient_identifier==rownames(APOBEC.data)[i]))
}

head(APOBEC.data)


library(ggplot2)
apobec.signal.vs.percent <- ggplot(data=APOBEC.data,aes(y=percent_TCW_TKW, x=APOBEC_pos_neg))+ geom_boxplot() + geom_jitter(alpha=0.4,width = 0.1)  + facet_grid(HPV_status ~.) + coord_flip() + labs(x="APOBEC signal",y="Percent of SNV that are TCW --> TKW") + theme_bw()
ggsave(filename = "~/Documents/manuscripts/APOBEC_HNSCC/Figures/APOBEC_signal_vs_percent.png",plot = apobec.signal.vs.percent)


tcw.percent.apo.pos.hpv.pos <- APOBEC.data$percent_TCW_TKW[which(APOBEC.data$HPV_status=="HPV+" & APOBEC.data$APOBEC_pos_neg==T)]
tcw.percent.apo.pos.hpv.neg <- APOBEC.data$percent_TCW_TKW[which(APOBEC.data$HPV_status=="HPV-" & APOBEC.data$APOBEC_pos_neg==T)]

median(tcw.percent.apo.pos.hpv.pos)
median(tcw.percent.apo.pos.hpv.neg)

wilcox.test(tcw.percent.apo.pos.hpv.pos, tcw.percent.apo.pos.hpv.neg, paired = F, conf.int = T)



tcw.percent.apo.neg.hpv.pos <- APOBEC.data$percent_TCW_TKW[which(APOBEC.data$HPV_status=="HPV+" & APOBEC.data$APOBEC_pos_neg==F)]
tcw.percent.apo.neg.hpv.neg <- APOBEC.data$percent_TCW_TKW[which(APOBEC.data$HPV_status=="HPV-" & APOBEC.data$APOBEC_pos_neg==F)]

wilcox.test(tcw.percent.apo.neg.hpv.pos,tcw.percent.apo.neg.hpv.neg,conf.int = T,paired = F)


# from https://www.nature.com/articles/ncomms3513#supplementary-information 
# Tang et al 2013 The landscape of viral expression and host gene fusion and adaptation in human cancer
# Where Hinderson 2014 got their viral data from 
tang.2013 <- read.csv(file = "input_data/ncomms3513-s2.txt",skip = 3,sep="\t",header = T,stringsAsFactors = F)

head(tang.2013)

tang.2013 <- subset(tang.2013, Cancer=="HNSC")

tang.2013$Sample.barcode[1]

names.split <- strsplit(tang.2013$Sample.barcode,split="-")

tang.2013$sample.short <- NA
for(i in 1:nrow(tang.2013)){
 tang.2013$sample.short[i] <- paste(names.split[[i]][1:3],collapse = "-")
}

both.datasets <- union(tang.2013$sample.short,apobec_hnscc_df$patient_id)

HPV.comparison <- as.data.frame(matrix(data = NA,nrow=length(both.datasets),ncol=3))
colnames(HPV.comparison) <- c("tumor","tang2013","ours")

HPV.comparison$tumor <- both.datasets

for(i in 1:nrow(HPV.comparison)){
  if(length(which(apobec_hnscc_df$patient_id==HPV.comparison$tumor[i]))==1){
    HPV.comparison$ours[i] <- apobec_hnscc_df$categ[which(apobec_hnscc_df$patient_id==HPV.comparison$tumor[i])]
  }
  if(length(which(tang.2013$sample.short==HPV.comparison$tumor[i]))==1){
   HPV.comparison$tang2013[i] <-  tang.2013$Virus.description[which(tang.2013$sample.short==HPV.comparison$tumor[i])]
   if(HPV.comparison$tang2013[i]=="N/A"){ HPV.comparison$tang2013[i] <- "not_detected"}
  }
  
}


HPV.comparison







