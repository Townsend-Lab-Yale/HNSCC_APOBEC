TCGA_HPVstatus <- read.table(file = "~/Documents/AcademicPapers/TCGA_HNSC_nature_supp/1.2_HPVstatus.txt",header = T,stringsAsFactors = F,sep = "\t")

head(TCGA_HPVstatus)


load("output_data/HNSC_selection_with_APOBEC.RData")

head(HNSC.selection.output)

TCGA_positives <- TCGA_HPVstatus$Barcode[which(TCGA_HPVstatus$Final_HPV_Status=="Positive")] 

Us_positives <- HNSC.selection.output$Unique_patient_identifier[which(HNSC.selection.output$HPV_call=="HPV+")]

TCGA_positives %in% Us_positives

TCGA_positives[!(TCGA_positives %in% Us_positives)]

TCGA_HPVstatus[which(TCGA_HPVstatus$Barcode=="TCGA-CR-7404"),]

HNSC.selection.output$Unique_patient_identifier[which(HNSC.selection.output$Unique_patient_identifier %in% TCGA_positives & HNSC.selection.output$Gene=="TP53")]


HNSC.selection.output$Unique_patient_identifier[which(HNSC.selection.output$Unique_patient_identifier %in% Us_positives & HNSC.selection.output$Gene=="TP53")]


BROAD.clinical <- read.csv("input_data/BROAD/gdac.broadinstitute.org_HNSC.Merge_Clinical.Level_1.2016012800.0.0/HNSC.clin.merged.txt",header = F,stringsAsFactors = F,sep="\t",row.names = 1)

head(BROAD.clinical)
BROAD_patients <- toupper(BROAD.clinical["patient.bcr_patient_barcode",])

toupper(BROAD.clinical["patient.hpv_test_results.hpv_test_result.hpv_status",])

BROAD_HPV_status <- matrix(data = c(toupper(BROAD.clinical["patient.bcr_patient_barcode",]),toupper(BROAD.clinical["patient.hpv_test_results.hpv_test_result.hpv_status",]),toupper(BROAD.clinical["patient.hpv_test_results.hpv_test_result.hpv_call_1",]),toupper(BROAD.clinical["patient.hpv_test_results.hpv_test_result.hpv_call_2",]),toupper(BROAD.clinical["patient.hpv_test_results.hpv_test_result.hpv_call_3",])),nrow = length(toupper(BROAD.clinical["patient.bcr_patient_barcode",])),ncol = 5,byrow = F)

BROAD_HPV_status <- BROAD_HPV_status[BROAD_HPV_status[,1] %in% unique(HNSC.MAF$Unique_patient_identifier),] # just the tumors with sequencing data

head(BROAD_HPV_status)
nrow(BROAD_HPV_status)

BROAD.positive.tumors <- BROAD_HPV_status[,1][which(BROAD_HPV_status[,2]=="POSITIVE")]

rownames(BROAD_HPV_status) <- BROAD_HPV_status[,1]

BROAD_HPV_status[TCGA_positives,]

BROAD_HPV_status

nrow(BROAD_HPV_status[which(!is.na(BROAD_HPV_status[,3])),])
BROAD_HPV_status[which(BROAD_HPV_status[,2]=="NEGATIVE"),]
BROAD_HPV_status[which(BROAD_HPV_status[,2]=="POSITIVE"),]

length(which(BROAD_HPV_status[,2]=="POSITIVE"))


# testing for TP53
BROAD_tp53 <- data.frame(tumor= BROAD.positive.tumors, TP53_mutation=NA)

for(i in 1:nrow(BROAD_tp53)){
  if(length(which(HNSC.MAF$Unique_patient_identifier==BROAD_tp53$tumor[i]))>0){
    BROAD_tp53$TP53_mutation[i] <- length(which(HNSC.MAF$Hugo_Symbol[which(HNSC.MAF$Unique_patient_identifier==BROAD_tp53$tumor[i])]=="TP53"))
  }
}

nrow(BROAD_tp53)

length(which(BROAD_tp53$TP53_mutation>0))

unique(apobec_hnscc_df$categ)
length(unique(apobec_hnscc_df$patient_id))

# virusscan.calls <- matrix(data = c(unique(apobec_hnscc_df$patient_id),rep(NA,508)),nrow = 508,ncol=2,byrow = F)
# for(i in 1:nrow(virusscan.calls)){
#   
# }


virusscan.calls <- apobec_hnscc_df[,c("patient_id","categ")]

virusscan.calls[which(is.na(virusscan.calls)),]

length(which(virusscan.calls[,"categ"]=="HPV+"))



virusscan.calls.standalone <- read.csv(file = "~/Documents/slack_downloads/virusscan_hnsc_hpv_reads.csv",header = T,stringsAsFactors = F)

nrow(virusscan.calls.standalone)

length(which(virusscan.calls.standalone$patient_id %in% apobec_hnscc_df$patient_id))

missing <- as.data.frame(apobec_hnscc_df[which(!(apobec_hnscc_df$patient_id %in% virusscan.calls.standalone$patient_id)),])

unique(missing$categ)

which(virusscan.calls.standalone$patient_id=="TCGA-CQ-7064")
virusscan.calls.standalone$one_over_1000 <- NA

write.table(x = missing,file = "output_data/missing_from_virusscan_but_in_categorization_file.txt",sep = "\t",quote = F,row.names = F)




head(as.data.frame(virusscan.calls),510)
