

# Used RSEM data from cbioportal (http://www.cbioportal.org/)
# Steps: 
# click the DOWNLOAD DATA tab
# specify "Head and NEck Squamous Cell Carcinoma (TCGA, PanCancer Atlas)  (523 samples as of Sept 25, 2018)
# specify "mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)
# Enter all APOBEC3* genes
# Click download

APOBEC_RSEM_data <- read.delim(file = "expression/cbioportal/cBioPortal_data.txt",header = T,stringsAsFactors = F)
APOBEC_RSEM_data <- APOBEC_RSEM_data[,-1]

APOBEC_RSEM_data <- APOBEC_RSEM_data[,-525]

APOBEC_RSEM_data_melt <- melt(APOBEC_RSEM_data,id.vars = "COMMON")

APOBEC_RSEM_data_melt$variable <- gsub(pattern = "\\.",replacement = "-",x = APOBEC_RSEM_data_melt$variable)

APOBEC_RSEM_data_melt$variable <- gsub(pattern = "-01",replacement = "",x = APOBEC_RSEM_data_melt$variable)




APOBEC_RSEM_data_melt$HPV_classification <- NA

for(i in 1:nrow(APOBEC_RSEM_data_melt)){
  if(APOBEC_RSEM_data_melt$variable[i] %in% HPV.classification$Tumor_name){
    APOBEC_RSEM_data_melt$HPV_classification[i] <- HPV.classification$Classification[which(HPV.classification$Tumor_name==APOBEC_RSEM_data_melt$variable[i])]
  }
}

APOBEC_RSEM_data_melt <- APOBEC_RSEM_data_melt[!is.na(APOBEC_RSEM_data_melt$HPV_classification),]

# head( APOBEC_expression_data[,-1])

# APOBEC_expression_data_melted <- melt(data = APOBEC_expression_data[,-1],id.vars = c("SAMPLE_ID_short","HPV_classification"))

# head(APOBEC_expression_data_melted)

library(ggplot2)
APOBEC_RSEM_data_melt <- APOBEC_RSEM_data_melt[-which(is.na(APOBEC_RSEM_data_melt$value )),]

plot_text_size <- 12

apobec_plots <- ggplot(data = APOBEC_RSEM_data_melt, aes(x=HPV_classification,y=value)) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha=0.2,width = 0.2)  + facet_wrap(~COMMON,scales = "free") + theme_classic()  + labs(x="HPV classification", y="mRNA expression") + theme(strip.text = element_text(size=plot_text_size), axis.title = element_text(size=plot_text_size),axis.text = element_text(size=plot_text_size))

ggsave(apobec_plots,filename = "Figures/apobec_expression.png",width = 6.5,height = 6.5)

ggsave(filename = "Figures/FigS1_A3_expression.eps",plot = apobec_plots,width = 6.5,height = 6.5,dpi = 300,device=cairo_ps, fallback_resolution = 300)

t.test(data=subset(APOBEC_RSEM_data_melt, COMMON=="APOBEC3A"), value~HPV_classification,alternative="greater")

t.test(data=subset(APOBEC_RSEM_data_melt, COMMON=="APOBEC3B"), value~HPV_classification,alternative="greater")

t.test(data=subset(APOBEC_RSEM_data_melt, COMMON=="APOBEC3C"), value~HPV_classification,alternative="greater")

t.test(data=subset(APOBEC_RSEM_data_melt, COMMON=="APOBEC3D"), value~HPV_classification,alternative="greater")

t.test(data=subset(APOBEC_RSEM_data_melt, COMMON=="APOBEC3F"), value~HPV_classification,alternative="greater")

t.test(data=subset(APOBEC_RSEM_data_melt, COMMON=="APOBEC3G"), value~HPV_classification,alternative="greater")

t.test(data=subset(APOBEC_RSEM_data_melt, COMMON=="APOBEC3H"), value~HPV_classification,alternative="greater")

unique_tumors <- data.frame(tumors=unique(APOBEC_RSEM_data_melt$variable), HPV_status=NA,stringsAsFactors = F)

for(i in 1:nrow(unique_tumors)){
 unique_tumors$HPV_status[i] <- as.character(APOBEC_expression_data_melted$HPV_classification)[which(as.character(APOBEC_expression_data_melted$SAMPLE_ID_short)==as.character(unique_tumors$tumors)[i])[1] ]
}

head(unique_tumors)


length(which(unique_tumors$HPV_status=="HPV−"))
length(which(unique_tumors$HPV_status=="HPV+"))


# wilcox.test(data=subset(APOBEC_RSEM_data_melted, variable=="APOBEC3A"), value~HPV_classification,alternative="greater")

