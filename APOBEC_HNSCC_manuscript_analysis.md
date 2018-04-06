Analysis for APOBEC-induced mutations and their cancer effect size in head and neck squamous cell carcinoma
================
Vincent L. Cannataro

-   [processing MAF file and splitting into HPV pos and neg](#processing-maf-file-and-splitting-into-hpv-pos-and-neg)
-   [Trinucleotide heatmaps](#trinucleotide-heatmaps)
-   [Gene-level mutation rates](#gene-level-mutation-rates)
-   [Calculating trinucleotide context with weights for each tumor](#calculating-trinucleotide-context-with-weights-for-each-tumor)
-   [Loading selection output and merging with APOBEC and HPV status](#loading-selection-output-and-merging-with-apobec-and-hpv-status)
-   [Selection and tornado plots](#selection-and-tornado-plots)
    -   [tornado plots](#tornado-plots)
-   [heatmap and dendrogram of the signatures](#heatmap-and-dendrogram-of-the-signatures)
-   [HPV vs SNP count](#hpv-vs-snp-count)
-   [Logistic regression of APOBEC and HPV](#logistic-regression-of-apobec-and-hpv)
    -   [Mutations that have the highest prevalence and selection intensity in HNSCC](#mutations-that-have-the-highest-prevalence-and-selection-intensity-in-hnscc)

processing MAF file and splitting into HPV pos and neg
======================================================

``` r
# import data from the VirusScan manuscript --- Cao, S., Wendl, M. C., Wyczalkowski, M. A., Wylie, K., Ye, K., Jayasinghe, R., … Ding, L. (2016). Divergent viral presentation among human tumors and adjacent normal tissues. Scientific Reports, 6(May), 28294. https://doi.org/10.1038/srep28294 
virusscan.data <- read.table(file = "input_data/virusscan/vscan_counts.tsv",sep = "\t",header = T,stringsAsFactors = F)
```

``` r
# load in the MAF file from the NCI 
HNSC.MAF <- read.csv("input_data/NCI/gdc_download_20180201_160847/1aa33f25-3893-4f37-a6a4-361c9785d07e/TCGA.HNSC.mutect.1aa33f25-3893-4f37-a6a4-361c9785d07e.DR-10.0.somatic.maf",skip=5,header = T,sep = "\t",stringsAsFactors = F)

# source custom functions
source("R/flip_function.R")
source("R/unique_tumor_addition.R")
source("R/hg39_to_hg19_converter.R")
```

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, cbind, colMeans,
    ##     colnames, colSums, do.call, duplicated, eval, evalq, Filter,
    ##     Find, get, grep, grepl, intersect, is.unsorted, lapply,
    ##     lengths, Map, mapply, match, mget, order, paste, pmax,
    ##     pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce,
    ##     rowMeans, rownames, rowSums, sapply, setdiff, sort, table,
    ##     tapply, union, unique, unsplit, which, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## Loading required package: GenomeInfoDb

``` r
source("R/DNP_remover.R")
source("R/tumor_allele_adder.R")

# convert the MAF to hg19 coordinates 
HNSC.MAF <- hg38.to.hg19.converter(chain = "input_data/hg38Tohg19.chain",hg38_maf = HNSC.MAF)
```

    ## Loading in specified MAF...

    ## Number of rows in the MAF that failed to convert:  3

``` r
NCI.MAF <- read.csv("input_data/NCI/gdc_download_20180201_160847/1aa33f25-3893-4f37-a6a4-361c9785d07e/TCGA.HNSC.mutect.1aa33f25-3893-4f37-a6a4-361c9785d07e.DR-10.0.somatic.maf",skip=5,header = T,sep = "\t",stringsAsFactors = F)
NCI.MAF <- unique.tumor.addition.function(MAF.file = NCI.MAF,non.TCGA.characters.to.keep = 'all')
```

    ## Summary statistics of the number of mutations per unique tumor:

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    1.00   88.75  139.00  201.40  207.25 4050.00

``` r
message("Number of NCI tumors in the VirusScan dataset:")
```

    ## Number of NCI tumors in the VirusScan dataset:

``` r
length(which(unique(NCI.MAF$Unique_patient_identifier) %in% virusscan.data$patient_id))
```

    ## [1] 477

``` r
# Processed data from Hedberg et al., 2016 (doi 10.1172/JCI82066) using https://github.com/Townsend-Lab-Yale/HNSCC_APOBEC/blob/master/process_yale_data.ipynb
Hedberg.2016 <- read.csv(file = "output_data/yale_filtered.maf",header = T,sep = "\t",stringsAsFactors = F)
as.character(unique(Hedberg.2016$Tumor_Sample_Barcode)) # Tumor names from Hedberg et al., 2016. 
```

    ##  [1] "PY-10T" "PY-11T" "PY-12T" "PY-13T" "PY-14T" "PY-15T" "PY-16T"
    ##  [8] "PY-17T" "PY-19T" "PY-1T"  "PY-20T" "PY-21T" "PY-22T" "PY-23T"
    ## [15] "PY-24T" "PY-25T" "PY-3T"  "PY-4T"  "PY-5T"  "PY-6T"  "PY-7T" 
    ## [22] "PY-8T"

``` r
Hedberg.2016 <- Hedberg.2016[-which(Hedberg.2016$Tumor_Sample_Barcode %in% c("PY-19T","PY-1T","PY-14T","PY-13T","PY-7T")),] # filtering out tumors that are also within the TCGA data set


source("R/merging_NCI_and_local_MAF_files.R")

HNSC.MAF <- merging_TCGA_and_local_MAFdata_function(NCI_data = HNSC.MAF,Local_data = Hedberg.2016)
```

    ## These are the important headers that need to be contained in both files:

    ##  [1] "Hugo_Symbol"            "Chromosome"            
    ##  [3] "Tumor_Seq_Allele2"      "Variant_Classification"
    ##  [5] "Variant_Type"           "trv_type"              
    ##  [7] "transcript_error"       "Reference_Allele"      
    ##  [9] "Start_Position"         "strand"                
    ## [11] "Tumor_Sample_Barcode"   "t_ref_count"           
    ## [13] "t_alt_count"

    ## Important headers not in NCI_data:

    ## [1] "trv_type"         "transcript_error"

    ## Important headers not in Local_data:

    ## [1] "trv_type"         "transcript_error" "strand"

    ## Making sure all the essential column headers are the same so they can be properly merged...

    ## [1] "Local_data is missing column name header strand"
    ## [1] "Could not find strand header in  Local_data . You need to manually find the appropriate header and change it to strand"

    ## Still a problem and need to be fixed manually:

    ## Important headers not in NCI_data:

    ## [1] "trv_type"         "transcript_error"

    ## Important headers not in Local_data:

    ## [1] "trv_type"         "transcript_error" "strand"

    ## Merging the data frames along their common headers...

    ## Merging Completed

``` r
# tail(HNSC.MAF)

HNSC.MAF <- unique.tumor.addition.function(MAF.file = HNSC.MAF,non.TCGA.characters.to.keep = 'all')
```

    ## Summary statistics of the number of mutations per unique tumor:

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##     1.0    86.0   137.0   197.9   206.0  4049.0

``` r
# remove potential DNP 
HNSC.MAF <- DNP.remover(MAF = HNSC.MAF)
```

    ## Removing possible DNP

    ## Total count of potential DNP removed:  2128

    ## DNP removal complete

    ## Deleting any mutations detected in TCGA recurrent tumors

``` r
# add a column that is just the variant allele 
HNSC.MAF <- tumor.allele.adder(MAF = HNSC.MAF)
# HNSC.MAF <- unique.tumor.addition.function(MAF.file = HNSC.MAF,non.TCGA.characters.to.keep = 'all')

HNSC.MAF$HPV_call <- NA
HNSC.MAF$Virusscan_counts <- NA

for(i in 1:length(unique(HNSC.MAF$Unique_patient_identifier))){
  if(length(which(virusscan.data$patient_id==unique(HNSC.MAF$Unique_patient_identifier)[i]))>0){
    HNSC.MAF$HPV_call[HNSC.MAF$Unique_patient_identifier==unique(HNSC.MAF$Unique_patient_identifier)[i]] <- 
      #apobec_hnscc_df$categ[which(apobec_hnscc_df$patient_id==unique(HNSC.MAF$Unique_patient_identifier)[i])]
      ifelse(virusscan.data$VScan_counts[which(virusscan.data$patient_id==unique(HNSC.MAF$Unique_patient_identifier)[i])]>100,
             "HPV+",
             "HPV−")
    HNSC.MAF$Virusscan_counts[HNSC.MAF$Unique_patient_identifier==unique(HNSC.MAF$Unique_patient_identifier)[i]] <- 
      virusscan.data$VScan_counts[which(virusscan.data$patient_id==unique(HNSC.MAF$Unique_patient_identifier)[i])]
  }
}

# assigning HPV status to PY tumors
HNSC.MAF$HPV_call[which(startsWith(x = HNSC.MAF$Unique_patient_identifier,prefix = "PY"))] <- "HPV−"
HNSC.MAF$HPV_call[which(HNSC.MAF$Unique_patient_identifier=="PY-16T")] <- "HPV+"

# assigning HPV status to tumors with our VirusScan results
# HNSC.MAF$HPV_call[which(HNSC.MAF$Unique_patient_identifier=="TCGA-QK-A6IF")] <- "HPV+"
our_VirusScan_results <- read.csv(file = "input_data/VirusScan_RPHM_on_24_additional_TCGA_patients.txt",header = T,sep = "\t",stringsAsFactors = F)

for(i in 1:nrow(our_VirusScan_results)){
 if(length(which(HNSC.MAF$Unique_patient_identifier==our_VirusScan_results$patient_id[i]))>0){
   HNSC.MAF$HPV_call[which(HNSC.MAF$Unique_patient_identifier==our_VirusScan_results$patient_id[i])] <- 
     ifelse(our_VirusScan_results$rphm[i]>100,"HPV+","HPV−")
   
 }
  
}


# results of two tumors that did not have RNA-seq data (we could not run Virusscan) but did have consistent results from 
# p16 and ISH clinical data
HNSC.MAF$HPV_call[which(HNSC.MAF$Unique_patient_identifier=="TCGA-CN-A63Y")] <- "HPV+"
HNSC.MAF$HPV_call[which(HNSC.MAF$Unique_patient_identifier=="TCGA-IQ-A61L")] <- "HPV−"

HNSC.MAF.hpvpos <- HNSC.MAF[which(HNSC.MAF$HPV_call=="HPV+"),]
MAF_for_analysis <- HNSC.MAF.hpvpos
length(unique(HNSC.MAF.hpvpos$Unique_patient_identifier))
```

    ## [1] 69

``` r
save(MAF_for_analysis, file="output_data/HNSC_HPVpos_MAF.RData")

HNSC.MAF.hpvneg <- HNSC.MAF[which(HNSC.MAF$HPV_call=="HPV−"),]
MAF_for_analysis <- HNSC.MAF.hpvneg
length(unique(HNSC.MAF.hpvneg$Unique_patient_identifier))
```

    ## [1] 451

``` r
save(MAF_for_analysis, file="output_data/HNSC_HPVneg_MAF.RData")


MAF_for_analysis <- HNSC.MAF
length(unique(HNSC.MAF$Unique_patient_identifier))
```

    ## [1] 525

``` r
save(MAF_for_analysis, file="output_data/HNSC_MAF.RData")

message("HPV positive tumor with TP53 mutation:") 
```

    ## HPV positive tumor with TP53 mutation:

``` r
HNSC.MAF.hpvpos$Unique_patient_identifier[which(HNSC.MAF.hpvpos$Hugo_Symbol=="TP53")]
```

    ## [1] "TCGA-CR-7368"

``` r
HNSC.MAF.hpvpos[which(HNSC.MAF.hpvpos$Hugo_Symbol=="TP53"),]
```

    ##       Chromosome Start_Position End_Position Hugo_Symbol Entrez_Gene_Id
    ## 57340         17        7577538      7577538        TP53           7157
    ##       Center                    NCBI_Build Variant_Classification
    ## 57340     BI Converted_from_GRCh38_to_hg19      Missense_Mutation
    ##       Variant_Type Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2
    ## 57340          SNP                C                 C                 T
    ##         dbSNP_RS      dbSNP_Val_Status         Tumor_Sample_Barcode
    ## 57340 rs11540652 byCluster;byFrequency TCGA-CR-7368-01A-11D-2129-08
    ##        Matched_Norm_Sample_Barcode Match_Norm_Seq_Allele1
    ## 57340 TCGA-CR-7368-10A-01D-2129-08                   <NA>
    ##       Match_Norm_Seq_Allele2 Tumor_Validation_Allele1
    ## 57340                   <NA>                       NA
    ##       Tumor_Validation_Allele2 Match_Norm_Validation_Allele1
    ## 57340                       NA                            NA
    ##       Match_Norm_Validation_Allele2 Verification_Status Validation_Status
    ## 57340                            NA                  NA                NA
    ##       Mutation_Status Sequencing_Phase Sequence_Source Validation_Method
    ## 57340         Somatic               NA              NA                NA
    ##       Score BAM_File           Sequencer
    ## 57340    NA       NA Illumina HiSeq 2000
    ##                          Tumor_Sample_UUID
    ## 57340 4b194ab3-d213-4a7a-be46-909b4f0c7291
    ##                   Matched_Norm_Sample_UUID    HGVSc       HGVSp
    ## 57340 2f1b122e-205d-44be-a6e1-1fc7f2271d99 c.743G>A p.Arg248Gln
    ##       HGVSp_Short   Transcript_ID Exon_Number t_depth t_ref_count
    ## 57340     p.R248Q ENST00000269305        7/11      77          50
    ##       t_alt_count n_depth n_ref_count n_alt_count
    ## 57340          27      98        <NA>        <NA>
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 all_effects
    ## 57340 TP53,missense_variant,p.R248Q,ENST00000617185,NM_001126114.2,c.743G>A,MODERATE,,deleterious(0.02),probably_damaging(0.999),-1;TP53,missense_variant,p.R248Q,ENST00000420246,,c.743G>A,MODERATE,,deleterious(0.02),probably_damaging(0.999),-1;TP53,missense_variant,p.R209Q,ENST00000622645,NM_001276696.1,c.626G>A,MODERATE,,deleterious(0.05),probably_damaging(0.999),-1;TP53,missense_variant,p.R209Q,ENST00000610292,NM_001126118.1,c.626G>A,MODERATE,,deleterious(0.02),probably_damaging(1),-1;TP53,missense_variant,p.R209Q,ENST00000610538,NM_001276695.1,c.626G>A,MODERATE,,deleterious(0.02),probably_damaging(0.996),-1;TP53,missense_variant,p.R248Q,ENST00000455263,NM_001126113.2,c.743G>A,MODERATE,,deleterious(0.02),probably_damaging(0.996),-1;TP53,missense_variant,p.R248Q,ENST00000269305,NM_000546.5&NM_001126112.2,c.743G>A,MODERATE,YES,deleterious(0),probably_damaging(1),-1;TP53,missense_variant,p.R209Q,ENST00000620739,NM_001276761.1&NM_001276760.1,c.626G>A,MODERATE,,deleterious(0.02),probably_damaging(1),-1;TP53,missense_variant,p.R209Q,ENST00000619485,,c.626G>A,MODERATE,,deleterious(0.02),probably_damaging(1),-1;TP53,missense_variant,p.R248Q,ENST00000445888,,c.743G>A,MODERATE,,deleterious(0),probably_damaging(1),-1;TP53,missense_variant,p.R116Q,ENST00000510385,NM_001126116.1,c.347G>A,MODERATE,,tolerated(0.06),probably_damaging(0.999),-1;TP53,missense_variant,p.R89Q,ENST00000618944,NM_001276698.1,c.266G>A,MODERATE,,tolerated(0.07),probably_damaging(0.999),-1;TP53,missense_variant,p.R89Q,ENST00000610623,NM_001276699.1,c.266G>A,MODERATE,,tolerated(0.07),probably_damaging(0.996),-1;TP53,missense_variant,p.R116Q,ENST00000504290,NM_001126117.1,c.347G>A,MODERATE,,tolerated(0.06),probably_damaging(0.996),-1;TP53,missense_variant,p.R89Q,ENST00000619186,NM_001276697.1,c.266G>A,MODERATE,,deleterious(0.03),probably_damaging(1),-1;TP53,missense_variant,p.R116Q,ENST00000504937,NM_001126115.1,c.347G>A,MODERATE,,deleterious(0.03),probably_damaging(1),-1;TP53,missense_variant,p.R248Q,ENST00000359597,,c.743G>A,MODERATE,,deleterious(0.02),probably_damaging(0.999),-1;TP53,missense_variant,p.R237Q,ENST00000615910,,c.710G>A,MODERATE,,deleterious(0),probably_damaging(0.992),-1;TP53,missense_variant,p.R248Q,ENST00000413465,,c.743G>A,MODERATE,,deleterious(0.02),probably_damaging(0.999),-1;TP53,missense_variant,p.R116Q,ENST00000509690,,c.347G>A,MODERATE,,tolerated(0.07),probably_damaging(1),-1;TP53,missense_variant,p.R155Q,ENST00000514944,,c.464G>A,MODERATE,,deleterious(0.04),probably_damaging(0.999),-1;TP53,downstream_gene_variant,,ENST00000508793,,,MODIFIER,,,,-1;TP53,downstream_gene_variant,,ENST00000604348,,,MODIFIER,,,,-1;TP53,downstream_gene_variant,,ENST00000503591,,,MODIFIER,,,,-1;TP53,upstream_gene_variant,,ENST00000576024,,,MODIFIER,,,,-1;TP53,downstream_gene_variant,,ENST00000574684,,,MODIFIER,,,,-1;TP53,downstream_gene_variant,,ENST00000505014,,,MODIFIER,,,,-1
    ##       Allele            Gene         Feature Feature_type      Consequence
    ## 57340      T ENSG00000141510 ENST00000269305   Transcript missense_variant
    ##       cDNA_position CDS_position Protein_position Amino_acids  Codons
    ## 57340      933/2579     743/1182          248/393         R/Q cGg/cAg
    ##               Existing_variation ALLELE_NUM DISTANCE SYMBOL SYMBOL_SOURCE
    ## 57340 rs11540652;TP53_g.13380G>A          1       NA   TP53          HGNC
    ##          HGNC_ID        BIOTYPE CANONICAL        CCDS            ENSP
    ## 57340 HGNC:11998 protein_coding       YES CCDS11118.1 ENSP00000269305
    ##       SWISSPROT TREMBL UNIPARC                     RefSeq           SIFT
    ## 57340    P04637 K7PPA8         NM_000546.5;NM_001126112.2 deleterious(0)
    ##                   PolyPhen EXON INTRON
    ## 57340 probably_damaging(1) 7/11       
    ##                                                                      DOMAINS
    ## 57340 Pfam_domain:PF00870;Prints_domain:PR00386;Superfamily_domains:SSF49417
    ##       GMAF AFR_MAF AMR_MAF ASN_MAF EAS_MAF EUR_MAF SAS_MAF AA_MAF EA_MAF
    ## 57340   NA      NA      NA      NA      NA      NA      NA     NA     NA
    ##         CLIN_SIG SOMATIC
    ## 57340 pathogenic        
    ##                                                      PUBMED MOTIF_NAME
    ## 57340 25032700;20377871;18798306;15450681;25105660;21264207         NA
    ##       MOTIF_POS HIGH_INF_POS MOTIF_SCORE_CHANGE   IMPACT PICK
    ## 57340        NA           NA                 NA MODERATE    1
    ##       VARIANT_CLASS TSL HGVS_OFFSET PHENO MINIMISED   ExAC_AF ExAC_AF_Adj
    ## 57340           SNV   1          NA   1;0         1 5.765e-05   5.768e-05
    ##       ExAC_AF_AFR ExAC_AF_AMR ExAC_AF_EAS ExAC_AF_FIN ExAC_AF_NFE
    ## 57340           0           0   0.0002312           0   7.493e-05
    ##       ExAC_AF_OTH ExAC_AF_SAS GENE_PHENO           FILTER
    ## 57340           0           0         NA panel_of_normals
    ##       Unique_patient_identifier Tumor_allele HPV_call Virusscan_counts
    ## 57340              TCGA-CR-7368            T     HPV+            31647

The number of tumors in the whole dataset: 525

The number of tumors in the HPV+ dataset: 69

The number of tumors in the HPV− dataset: 451

Tumors that are HPV+ and also have a mutation in TP53: TCGA-CR-7368

Tumors that are from TCGA and HPV+: 68

Tumors that are from TCGA and HPV+: 435

Trinucleotide heatmaps
======================

The SNV selection intensity pipeline was run on the HPV data. The pipeline may be found here: <https://github.com/Townsend-Lab-Yale/SNV_selection_intensity> and the associated manuscript is here: <https://doi.org/10.1101/229724>

``` r
library(ggplot2)
load("input_data/selection_from_cluster/HNSC_HPVpos/trinuc_output/trinuc_data_HNSC_HPVpos.RData")
HPV.pos.trinuc.mutation_data <- trinuc.mutation_data
HPV.pos.trinuc.heatmap <- ggplot(data=HPV.pos.trinuc.mutation_data, aes(Downstream, Upstream)) +
  geom_tile(aes(fill = proportion*100), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue", name="Percent")
HPV.pos.trinuc.heatmap <- HPV.pos.trinuc.heatmap + facet_grid(.~section_labels, labeller = label_parsed) 
HPV.pos.trinuc.heatmap <- HPV.pos.trinuc.heatmap +  geom_text(aes(label = round(proportion, 4)*100),size=2)
# p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
# panel.background = element_blank(), axis.line = element_line(colour = "black"))
HPV.pos.trinuc.heatmap <- HPV.pos.trinuc.heatmap + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                      panel.grid.minor = element_blank(),
                                                                      axis.ticks = element_blank(),
                                                                      strip.text=element_text(size=15),
                                                                      axis.title.x = element_text(size=15),
                                                                      axis.title.y = element_text(size=15),
                                                                      axis.text.x = element_text(size=12),
                                                                      axis.text.y=element_text(size=12),plot.title = element_text(hjust = 0.5),legend.position = "left") 
# ggtitle(paste("Trinucleotide profile for ",tumor.name,sep=""))
# p


HPV.pos.trinuc.heatmap
```

![](APOBEC_HNSCC_manuscript_analysis_files/figure-markdown_github/Figure%20trinuc%20context%20HPVpos-1.png)

``` r
ggsave(paste("Figures/","HPVpos","_trinuc_heatmap.png",sep=""),height = 1.5,width = 7,plot = HPV.pos.trinuc.heatmap,dpi=300)
```

``` r
load("input_data/selection_from_cluster/HNSC_HPVneg/trinuc_output/trinuc_data_HNSC_HPVneg.RData")
HPV.neg.trinuc.mutation_data <- trinuc.mutation_data
HPV.neg.trinuc.heatmap <- ggplot(data=HPV.neg.trinuc.mutation_data, aes(Downstream, Upstream)) +
  geom_tile(aes(fill = proportion*100), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue", name="Percent")
HPV.neg.trinuc.heatmap <- HPV.neg.trinuc.heatmap + facet_grid(.~section_labels, labeller = label_parsed) 
HPV.neg.trinuc.heatmap <- HPV.neg.trinuc.heatmap +  geom_text(aes(label = round(proportion, 4)*100),size=2)
# p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
# panel.background = element_blank(), axis.line = element_line(colour = "black"))
HPV.neg.trinuc.heatmap <- HPV.neg.trinuc.heatmap + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                      panel.grid.minor = element_blank(),
                                                                      axis.ticks = element_blank(),
                                                                      strip.text=element_text(size=15),
                                                                      axis.title.x = element_text(size=15),
                                                                      axis.title.y = element_text(size=15),
                                                                      axis.text.x = element_text(size=12),
                                                                      axis.text.y=element_text(size=12),plot.title = element_text(hjust = 0.5),legend.position = "left") 
# ggtitle(paste("Trinucleotide profile for ",tumor.name,sep=""))
# p

HPV.neg.trinuc.heatmap
```

![](APOBEC_HNSCC_manuscript_analysis_files/figure-markdown_github/Figure%20trinuc%20context%20HPVneg-1.png)

``` r
ggsave(filename = paste("Figures/","HPVneg","_trinuc_heatmap.png",sep=""),height = 1.5,width = 7,plot = HPV.neg.trinuc.heatmap,dpi = 300)
```

Gene-level mutation rates
=========================

``` r
HPV.pos.mut.rates <- read.csv(file = "input_data/selection_from_cluster/HNSC_HPVpos/mutsig_output/MAF_HNSC_HPVpos.txt.gene_rates.txt",header = T,sep = "\t",stringsAsFactors = F)
# head(HPV.pos.mut.rates)

HPV.neg.mut.rates <- read.csv(file = "input_data/selection_from_cluster/HNSC_HPVneg/mutsig_output/MAF_HNSC_HPVneg.txt.gene_rates.txt",header = T,sep = "\t",stringsAsFactors = F)
# head(HPV.neg.mut.rates)

# HPV.neg.mut.rates


all.equal(HPV.pos.mut.rates$gene,HPV.neg.mut.rates$gene)
```

    ## [1] TRUE

``` r
# HPV.neg.mut.rates$r_x_X[which(HPV.neg.mut.rates$gene=="PIK3CA")]
# HPV.pos.mut.rates$r_x_X[which(HPV.pos.mut.rates$gene=="PIK3CA")]

mutation_rates <- data.frame(gene=HPV.pos.mut.rates$gene,positive_mut_rates=HPV.pos.mut.rates$r_x_X,negative_mut_rates=HPV.neg.mut.rates$r_x_X)
library(ggrepel)
```

    ## Warning: package 'ggrepel' was built under R version 3.4.2

``` r
# mutation_rates[which(mutation_rates$gene=="PIK3CA"),]

source("R/fancy_scientific_code.R")




mutation_rates_full_scatter <- ggplot(data = mutation_rates, aes(x=positive_mut_rates,y=negative_mut_rates)) +
  geom_point(alpha=0.2,col="black",size=0.5) + 
  geom_smooth(method='lm',color="red") + 
  # geom_text_repel(data=subset(mutation_rates, positive_mut_rates > 5e-5 | negative_mut_rates > 3e-5 ), aes(label=gene),size=3) + 
  geom_point(data = subset(mutation_rates, gene %in% c("FBXW7","PIK3CA","MAPK1")),aes(x=positive_mut_rates,y=negative_mut_rates),size=3,col="blue") + 
  geom_text_repel(data = subset(mutation_rates, gene %in% c("FBXW7","PIK3CA","MAPK1")),aes(x=positive_mut_rates,y=negative_mut_rates,label=gene),size=5,col="blue",fontface="bold") +
  labs(x="Mutation rate in HPV+ tumors",y="Mutation rate in HPV− tumors") + 
  coord_equal(ratio=1) + 
  theme_bw() +
  geom_abline(slope=1, intercept=0) + 
  scale_x_continuous(labels=fancy_scientific) + 
  scale_y_continuous(labels=fancy_scientific) + theme(plot.margin = margin(r=.2,unit = "in"))


mutation_rates_reduced <- ggplot(data = mutation_rates, aes(x=positive_mut_rates,y=negative_mut_rates)) + geom_point(alpha=0.2,col="black",size=0.5) + geom_point(data = subset(mutation_rates, gene %in% c("FBXW7","PIK3CA","MAPK1")),aes(x=positive_mut_rates,y=negative_mut_rates),size=3,col="blue") + geom_text_repel(data = subset(mutation_rates, gene %in% c("FBXW7","PIK3CA","MAPK1")),aes(x=positive_mut_rates,y=negative_mut_rates,label=gene),size=6,col="blue",fontface="bold") +
  labs(x="Mutation rate in HPV+ tumors",y="Mutation rate in HPV− tumors") + coord_equal(ratio=1,xlim=c(0,0.5e-5),ylim=c(0,0.5e-5)) + theme_bw() + geom_smooth(method='lm',formula=y~x,color="red") + geom_abline(slope=1, intercept=0) + scale_x_continuous(labels=fancy_scientific) + scale_y_continuous(labels=fancy_scientific) 

# mutation_rates_full_scatter <- ggplot(data = mutation_rates, aes(y=positive_mut_rates,x=negative_mut_rates)) +
#   geom_point(alpha=0.2,col="black",size=0.5) + 
#     geom_smooth(method='lm',color="red") + 
#   # geom_text_repel(data=subset(mutation_rates, positive_mut_rates > 5e-5 | negative_mut_rates > 3e-5 ), aes(label=gene),size=3) + 
#   geom_point(data = subset(mutation_rates, gene %in% c("FBXW7","PIK3CA","MAPK1")),aes(y=positive_mut_rates,x=negative_mut_rates),size=3,col="blue") + 
#   geom_text_repel(data = subset(mutation_rates, gene %in% c("FBXW7","PIK3CA","MAPK1")),aes(y=positive_mut_rates,x=negative_mut_rates,label=gene),size=5,col="blue",fontface="bold",box.padding = 5) +
#   labs(y=bquote("Mutation rate in "~HPV^{"-"}~ "\n tumors"), x=bquote("Mutation rate in "~HPV^{"+"}~ "tumors"))+ 
#   coord_equal(ratio=1) + 
#   theme_bw() + 
#   geom_abline(slope=1, intercept=0) + 
#   scale_x_continuous(labels=fancy_scientific) + 
#   scale_y_continuous(labels=fancy_scientific) + theme(plot.margin = margin(0, 1, 0, 0,unit = "in"))

mutation_rates_full_scatter
```

![](APOBEC_HNSCC_manuscript_analysis_files/figure-markdown_github/mutation%20rate%20of%20HPV%20positive%20vs%20negative-1.png)

``` r
summary(lm(data = mutation_rates,formula = positive_mut_rates~negative_mut_rates))
```

    ## 
    ## Call:
    ## lm(formula = positive_mut_rates ~ negative_mut_rates, data = mutation_rates)
    ## 
    ## Residuals:
    ##        Min         1Q     Median         3Q        Max 
    ## -7.416e-06 -6.240e-07 -1.450e-07  3.970e-07  1.136e-04 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        1.097e-06  2.740e-08   40.04   <2e-16 ***
    ## negative_mut_rates 1.682e-01  1.074e-02   15.66   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2.364e-06 on 18860 degrees of freedom
    ## Multiple R-squared:  0.01283,    Adjusted R-squared:  0.01278 
    ## F-statistic: 245.1 on 1 and 18860 DF,  p-value: < 2.2e-16

``` r
ggsave(plot = mutation_rates_full_scatter,filename = "Figures/mutation_rates_full_scatter.png",height = 3,width = 8,dpi=300)

ggsave(plot = mutation_rates_reduced,filename = "Figures/mutation_rates_reduced.png",height = 2.5,width = 2.5)


# ggplot(data = mutation_rates, aes(x = positive_mut_rates, y=negative_mut_rates)) + geom_point() + geom_smooth() + coord_equal()
# ggplot(data = mutation_rates, aes(y = positive_mut_rates, x=negative_mut_rates)) + geom_point() + geom_smooth() + coord_equal()
```

Calculating trinucleotide context with weights for each tumor
=============================================================

``` r
library(reshape2)
```

    ## Warning: package 'reshape2' was built under R version 3.4.3

``` r
load("output_data/HNSC_MAF.RData")
HNSC.MAF <- MAF_for_analysis
source("R/trinuc_signatures_w_weights.R")
trinuc.contexts <- trinuc.profile.function_withweights(input.MAF = HNSC.MAF,
                                                       signature.choice = "signatures.cosmic", 
                                                       minimum.mutations.per.tumor=50,save.figs=F)
```

    ## Loading required package: deconstructSigs

    ## Removing all recurrent mutations...

    ## Finding the number of mutations per tumor

    ## Number of tumors over specified minimum mutation number of 50: 461

    ## Cleaning input to only contain tumors above the minimum...

    ## Calculating trinucleotide mutation counts...

    ## Warning in mut.to.sigs.input(mut.ref = input.MAF, sample.id = "Unique_patient_identifier", : Check ref bases -- not all match context:
    ##   TCGA-CV-7418:chr1:148346596:C:G

    ## Calculating individual tumor mutational signatures...

    ## No id variables; using all as measure variables

    ## Statistical summary of the proportion of the mutational signature in each tumor sample that is 'unknown'

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ## 0.00000 0.05560 0.09802 0.10846 0.15214 0.42909

``` r
length(trinuc.contexts[[2]]) #number of tumors with >50 mutations
```

    ## [1] 461

``` r
save(trinuc.contexts,file="output_data/trinuc_contexts_HNSC_MAF.RData")
```

Loading selection output and merging with APOBEC and HPV status
===============================================================

``` r
library(tidyverse)
```

    ## Warning: package 'tidyverse' was built under R version 3.4.2

    ## ── Attaching packages ──────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ tibble  1.4.2     ✔ purrr   0.2.4
    ## ✔ tidyr   0.8.0     ✔ dplyr   0.7.4
    ## ✔ readr   1.1.1     ✔ stringr 1.2.0
    ## ✔ tibble  1.4.2     ✔ forcats 0.2.0

    ## Warning: package 'tibble' was built under R version 3.4.3

    ## Warning: package 'tidyr' was built under R version 3.4.3

    ## Warning: package 'purrr' was built under R version 3.4.2

    ## Warning: package 'dplyr' was built under R version 3.4.2

    ## ── Conflicts ─────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::collapse()   masks IRanges::collapse()
    ## ✖ dplyr::combine()    masks BiocGenerics::combine()
    ## ✖ dplyr::desc()       masks IRanges::desc()
    ## ✖ tidyr::expand()     masks S4Vectors::expand()
    ## ✖ dplyr::filter()     masks stats::filter()
    ## ✖ dplyr::first()      masks S4Vectors::first()
    ## ✖ dplyr::lag()        masks stats::lag()
    ## ✖ ggplot2::Position() masks BiocGenerics::Position(), base::Position()
    ## ✖ purrr::reduce()     masks GenomicRanges::reduce(), IRanges::reduce()
    ## ✖ dplyr::rename()     masks S4Vectors::rename()
    ## ✖ dplyr::slice()      masks IRanges::slice()

``` r
load("input_data/selection_from_cluster/HNSC/selection_output/HNSC_all_selection_output.RData")
HNSC.selection.output <- as.tibble(selection.output$complete_mutation_data) %>%
  select(Gene, starts_with("Nucleo"), 
         Chromosome, starts_with("Reference"),
         starts_with("Alternative"), Tumor_origin, 
         Unique_patient_identifier, starts_with("Amino"), Codon_position, synonymous.mu, trinucs, Gamma_epistasis)  
# HNSC.selection.output

# Assigning HPV status
load("output_data/HNSC_MAF.RData")
HNSC.MAF <- MAF_for_analysis
HNSC.selection.output$HPV_call <- NA
for(i in 1:length(unique(HNSC.selection.output$Unique_patient_identifier))){
  HNSC.selection.output$HPV_call[which(HNSC.selection.output$Unique_patient_identifier == unique(HNSC.selection.output$Unique_patient_identifier)[i])] <- HNSC.MAF$HPV_call[which(HNSC.MAF$Unique_patient_identifier == unique(HNSC.selection.output$Unique_patient_identifier)[i])[1]]
}


load("output_data/trinuc_contexts_HNSC_MAF.RData")

weights.df <- trinuc.contexts$signature.weights[[1]]$weights
for(i in 2:length(trinuc.contexts$signature.weights)){
  weights.df <- rbind(weights.df,trinuc.contexts$signature.weights[[i]]$weights)
}

HNSC.selection.output$Sig_2 <- NA
HNSC.selection.output$Sig_13 <- NA


for(i in 1:length(unique(HNSC.selection.output$Unique_patient_identifier))){
  if(length(which(rownames(weights.df) == unique(HNSC.selection.output$Unique_patient_identifier)[i]))>0){
    HNSC.selection.output$Sig_2[which(HNSC.selection.output$Unique_patient_identifier == unique(HNSC.selection.output$Unique_patient_identifier)[i])] <- weights.df$Signature.2[which(rownames(weights.df) == unique(HNSC.selection.output$Unique_patient_identifier)[i])]
    HNSC.selection.output$Sig_13[which(HNSC.selection.output$Unique_patient_identifier == unique(HNSC.selection.output$Unique_patient_identifier)[i])] <- weights.df$Signature.13[which(rownames(weights.df) == unique(HNSC.selection.output$Unique_patient_identifier)[i])]
  }
}


HNSC.selection.output <- mutate(.data = HNSC.selection.output, APOBEC_weight = Sig_2 + Sig_13)

# head(HNSC.selection.output)


# Assigning TCW --> TKW trinucleotide context 



HNSC.selection.output$TCW_TKW <- NA
HNSC.selection.output$TCN_TKN <- NA

for(j in 1:nrow(HNSC.selection.output)){
  if(is.na(HNSC.selection.output$Nucleotide_trinuc_context[j])){
    #need to make a call of the same nucleotide position and amino acid alternative, find which is NOT NA, and make this the new "j" 
    matches <- which(HNSC.selection.output$Nucleotide_chromosome_position == HNSC.selection.output$Nucleotide_chromosome_position[j] &
                       HNSC.selection.output$Alternative_Nucleotide == HNSC.selection.output$Alternative_Nucleotide[j])
    
    new.j <- matches[which(!is.na(HNSC.selection.output$Nucleotide_trinuc_context[matches]))]
    if(((HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="TCA" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="T") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="TGA" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="A")) | 
       ((HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="TCT" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="T") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="AGA" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="A")) |
       ((HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="TCA" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="G") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="TGA" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="C")) |
       ((HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="TCT" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="G") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="AGA" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="C"))){
      HNSC.selection.output$TCW_TKW[j] <- 1
    }else{
      HNSC.selection.output$TCW_TKW[j] <- 0
    }
    
  }else{
    if(((HNSC.selection.output$Nucleotide_trinuc_context[j]=="TCA" & HNSC.selection.output$Alternative_Nucleotide[j]=="T") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[j]=="TGA" & HNSC.selection.output$Alternative_Nucleotide[j]=="A")) | 
       ((HNSC.selection.output$Nucleotide_trinuc_context[j]=="TCT" & HNSC.selection.output$Alternative_Nucleotide[j]=="T") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[j]=="AGA" & HNSC.selection.output$Alternative_Nucleotide[j]=="A")) |
       ((HNSC.selection.output$Nucleotide_trinuc_context[j]=="TCA" & HNSC.selection.output$Alternative_Nucleotide[j]=="G") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[j]=="TGA" & HNSC.selection.output$Alternative_Nucleotide[j]=="C")) |
       ((HNSC.selection.output$Nucleotide_trinuc_context[j]=="TCT" & HNSC.selection.output$Alternative_Nucleotide[j]=="G") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[j]=="AGA" & HNSC.selection.output$Alternative_Nucleotide[j]=="C"))){
      HNSC.selection.output$TCW_TKW[j] <- 1
    }else{
      HNSC.selection.output$TCW_TKW[j] <- 0
    }
  }
  
  # Now, mutations that could be TCN --> TKN
  if(is.na(HNSC.selection.output$Nucleotide_trinuc_context[j])){
    #need to make a call of the same nucleotide position and amino acid alternative, find which is NOT NA, and make this the new "j" 
    matches <- which(HNSC.selection.output$Nucleotide_chromosome_position == HNSC.selection.output$Nucleotide_chromosome_position[j] &
                       HNSC.selection.output$Alternative_Nucleotide == HNSC.selection.output$Alternative_Nucleotide[j])
    
    new.j <- matches[which(!is.na(HNSC.selection.output$Nucleotide_trinuc_context[matches]))]
    if(((HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="TCA" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="T") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="TGA" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="A")) | 
       ((HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="TCT" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="T") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="AGA" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="A")) |
       ((HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="TCA" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="G") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="TGA" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="C")) |
       ((HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="TCT" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="G") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="AGA" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="C")) |
       
       ((HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="TCC" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="T") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="GGA" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="A")) | 
       ((HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="TCG" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="T") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="CGA" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="A")) |
       ((HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="TCC" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="G") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="GGA" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="C")) |
       ((HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="TCG" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="G") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[new.j]=="CGA" & HNSC.selection.output$Alternative_Nucleotide[new.j]=="C"))
    ){
      HNSC.selection.output$TCN_TKN[j] <- 1
    }else{
      HNSC.selection.output$TCN_TKN[j] <- 0
    }
    
  }else{
    if(((HNSC.selection.output$Nucleotide_trinuc_context[j]=="TCA" & HNSC.selection.output$Alternative_Nucleotide[j]=="T") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[j]=="TGA" & HNSC.selection.output$Alternative_Nucleotide[j]=="A")) | 
       ((HNSC.selection.output$Nucleotide_trinuc_context[j]=="TCT" & HNSC.selection.output$Alternative_Nucleotide[j]=="T") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[j]=="AGA" & HNSC.selection.output$Alternative_Nucleotide[j]=="A")) |
       ((HNSC.selection.output$Nucleotide_trinuc_context[j]=="TCA" & HNSC.selection.output$Alternative_Nucleotide[j]=="G") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[j]=="TGA" & HNSC.selection.output$Alternative_Nucleotide[j]=="C")) |
       ((HNSC.selection.output$Nucleotide_trinuc_context[j]=="TCT" & HNSC.selection.output$Alternative_Nucleotide[j]=="G") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[j]=="AGA" & HNSC.selection.output$Alternative_Nucleotide[j]=="C"))|
       
       ((HNSC.selection.output$Nucleotide_trinuc_context[j]=="TCC" & HNSC.selection.output$Alternative_Nucleotide[j]=="T") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[j]=="GGA" & HNSC.selection.output$Alternative_Nucleotide[j]=="A")) | 
       ((HNSC.selection.output$Nucleotide_trinuc_context[j]=="TCG" & HNSC.selection.output$Alternative_Nucleotide[j]=="T") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[j]=="CGA" & HNSC.selection.output$Alternative_Nucleotide[j]=="A")) |
       ((HNSC.selection.output$Nucleotide_trinuc_context[j]=="TCC" & HNSC.selection.output$Alternative_Nucleotide[j]=="G") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[j]=="GGA" & HNSC.selection.output$Alternative_Nucleotide[j]=="C")) |
       ((HNSC.selection.output$Nucleotide_trinuc_context[j]=="TCG" & HNSC.selection.output$Alternative_Nucleotide[j]=="G") | 
        (HNSC.selection.output$Nucleotide_trinuc_context[j]=="CGA" & HNSC.selection.output$Alternative_Nucleotide[j]=="C"))){
      HNSC.selection.output$TCN_TKN[j] <- 1
    }else{
      HNSC.selection.output$TCN_TKN[j] <- 0
    }
  }
  
  
}


HNSC.selection.output$Name <- NA
for(i in 1:nrow(HNSC.selection.output)){
  HNSC.selection.output$Name[i] <- paste(HNSC.selection.output$Gene[i]," ",ifelse(!is.na(HNSC.selection.output$Amino_acid_reference[i]),paste(HNSC.selection.output$Amino_acid_reference[i],HNSC.selection.output$Amino_acid_position[i],HNSC.selection.output$Amino_acid_alternative[i]," ", HNSC.selection.output$Chromosome[i],"_",HNSC.selection.output$Nucleotide_chromosome_position[i],HNSC.selection.output$Alternative_Nucleotide[i],sep=""),paste(HNSC.selection.output$Reference_Nucleotide[i],HNSC.selection.output$Nucleotide_chromosome_position[i],HNSC.selection.output$Alternative_Nucleotide[i],"NCSNV")),sep="")
}

HNSC.selection.output$Name_short <- NA
for(i in 1:nrow(HNSC.selection.output)){
  HNSC.selection.output$Name_short[i] <- paste(HNSC.selection.output$Gene[i]," ",ifelse(!is.na(HNSC.selection.output$Amino_acid_reference[i]),paste(HNSC.selection.output$Amino_acid_reference[i],HNSC.selection.output$Amino_acid_position[i],HNSC.selection.output$Amino_acid_alternative[i],sep=""),paste(HNSC.selection.output$Reference_Nucleotide[i],HNSC.selection.output$Nucleotide_chromosome_position[i],HNSC.selection.output$Alternative_Nucleotide[i],"NCSNV")),sep="")
}



# HNSC.selection.output
save(HNSC.selection.output,file = "output_data/HNSC_selection_with_APOBEC.RData")

HNSC.selection.output.recur <- subset(HNSC.selection.output, Nucleotide_change_tally>1)
save(HNSC.selection.output.recur, file = "output_data/HNSC_selection_with_APOBEC_recur.RData")
```

``` r
trinuc.w.HPV <- weights.df

trinuc.w.HPV$HPV <- NA

for(i in 1:nrow(trinuc.w.HPV)){
  trinuc.w.HPV$HPV[i] <- HNSC.selection.output$HPV_call[which(HNSC.selection.output$Unique_patient_identifier==rownames(weights.df)[i])[1]]
}

# Among tumors that had enough substitutions that we could calculate mutational signatures ...
length(which(trinuc.w.HPV$HPV=="HPV+")) # ...how many tumors are HPV+ 
```

    ## [1] 47

``` r
length(which(trinuc.w.HPV$HPV=="HPV−")) # ...how many tumors are HPV- 
```

    ## [1] 409

``` r
length(which(is.na(trinuc.w.HPV$HPV)))
```

    ## [1] 5

``` r
save(trinuc.w.HPV,file = "output_data/signature_weights_w_HPV.RData")


length(unique(HNSC.selection.output$Unique_patient_identifier[which(HNSC.selection.output$HPV_call=="HPV+" & HNSC.selection.output$APOBEC_weight > 0)])) # Out of all tumors, how many were HPV+ and had an APOBEC signature 
```

    ## [1] 39

``` r
# length(which((trinuc.w.HPV$`Signature.2` > 0 | trinuc.w.HPV$`Signature.13`>0) & trinuc.w.HPV$HPV=="HPV+"))



length(which((trinuc.w.HPV$`Signature.2` > 0 | trinuc.w.HPV$`Signature.13`>0) & trinuc.w.HPV$HPV=="HPV−"))/length(which(trinuc.w.HPV$HPV=="HPV−")) # proportion of HPV- tumors with enough substitutions to measure signatures that have APOBEC signature
```

    ## [1] 0.6748166

``` r
length(which((trinuc.w.HPV$`Signature.2` > 0 | trinuc.w.HPV$`Signature.13`>0) & trinuc.w.HPV$HPV=="HPV+"))/length(which(trinuc.w.HPV$HPV=="HPV+")) # proportion of HPV+ tumors with enough substitutions to measure signatures that have APOBEC signature
```

    ## [1] 0.8297872

``` r
# mean weights
mean(trinuc.w.HPV$`Signature.2`[which(trinuc.w.HPV$HPV=="HPV+")])
```

    ## [1] 0.2111875

``` r
mean(trinuc.w.HPV$`Signature.13`[which(trinuc.w.HPV$HPV=="HPV+")])
```

    ## [1] 0.1601071

``` r
mean(trinuc.w.HPV$`Signature.2`[which(trinuc.w.HPV$HPV=="HPV−")])
```

    ## [1] 0.07625592

``` r
mean(trinuc.w.HPV$`Signature.13`[which(trinuc.w.HPV$HPV=="HPV−")])
```

    ## [1] 0.1105457

Number of tumors with greater than 50 mutations and a trinucleotide signature: 461

``` r
recur.hpv.pos <- subset(HNSC.selection.output.recur, HPV_call=="HPV+")
recur.hpv.pos <- subset(recur.hpv.pos, Name_short %in% names(table(recur.hpv.pos$Name_short))[which(table(recur.hpv.pos$Name_short)>1)]) #just the recurrent mutations 
hpv.pos.names <- unique(recur.hpv.pos$Name_short)

recur.hpv.neg <- subset(HNSC.selection.output.recur, HPV_call=="HPV−")
recur.hpv.neg <- subset(recur.hpv.neg, Name_short %in% names(table(recur.hpv.neg$Name_short))[which(table(recur.hpv.neg$Name_short)>1)]) #just the recurrent mutations 
hpv.neg.names <- unique(recur.hpv.neg$Name_short)

intersect(hpv.pos.names,hpv.neg.names)
```

    ## [1] "FBXW7 R505G"  "PIK3CA E545K" "PIK3CA E542K" "MAPK1 E322K"

``` r
both.names <- union(hpv.neg.names,hpv.pos.names)

#make a dataframe with nrows == union of names, and columnes the number in each type of tumor and then the prevalence in each type

prevalence.df <- as.data.frame(matrix(data = NA, nrow = length(both.names), ncol=5))
colnames(prevalence.df) <- c("Name","tally_HPVpos","tally_HPVneg","prev_HPVpos","prev_HPVneg")

prevalence.df$Name <- both.names
prevalence.df$tally_HPVpos <- 0
prevalence.df$tally_HPVneg <- 0
for(i in 1:nrow(prevalence.df)){
  if(length(which(hpv.pos.names == prevalence.df$Name[i]))>0){
    prevalence.df$tally_HPVpos[i] <- length(which(recur.hpv.pos$Name_short == prevalence.df$Name[i]))
  }
  if(length(which(hpv.neg.names == prevalence.df$Name[i]))>0){
    prevalence.df$tally_HPVneg[i] <- length(which(recur.hpv.neg$Name_short == prevalence.df$Name[i]))
  }
}

prevalence.df$prev_HPVpos <- prevalence.df$tally_HPVpos/length(unique(HNSC.selection.output$Tumor_origin[which(HNSC.selection.output$HPV_call=="HPV+")]))
prevalence.df$prev_HPVneg <- prevalence.df$tally_HPVneg/length(unique(HNSC.selection.output$Tumor_origin[which(HNSC.selection.output$HPV_call=="HPV−")])) 

library(ggrepel)
```

``` r
library(ggrepel)
common.text.size <- 6
library(cowplot)
```

    ## Warning: package 'cowplot' was built under R version 3.4.3

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

``` r
# getting them to match up! 
prev_plot <- ggplot(data = prevalence.df) + 
  geom_point(aes(y = prev_HPVneg, x = prev_HPVpos),alpha=0.5,size=1.5,color="red") + 
  geom_text_repel(data = subset(prevalence.df, (prev_HPVneg>0.025 | prev_HPVpos>0.04) | (prev_HPVneg>0.0 & prev_HPVpos>0)),aes(y = prev_HPVneg, x = prev_HPVpos,label=Name),box.padding =0.6,size=common.text.size*(5/14),color="black",segment.alpha = 0.3) + 
  theme_classic() + 
  expand_limits(x = 0, y = c(0,max(prevalence.df$prev_HPVneg)+.005)) + 
  geom_segment(aes(x = 0, xend =max(prevalence.df$prev_HPVneg)+.005, y = 0, yend =max(prevalence.df$prev_HPVneg)+.005),linetype="dashed") +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  labs(y=bquote("Prevalence in "~HPV^{"−"}~ "\n tumors"), x=bquote("Prevalence in "~HPV^{"+"}~ "tumors")) + theme(axis.text=element_text(size=common.text.size), axis.title=element_text(size=common.text.size,face="bold"),plot.margin = margin(r=.2,t=.0,l=.1,b=0.1,unit = "in")) + 
  # coord_equal() +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(margin = margin(t = -.2, r = 0, b = 0, l = 0)))  #+ geom_label(aes(x = prev_HPVneg, y = prev_HPVpos, label=Name))
# prev_plot

prev_plot_gridoff <- ggplot_gtable(ggplot_build(prev_plot))
prev_plot_gridoff$layout$clip[prev_plot_gridoff$layout$name == "panel"] <- "off"
library(grid)
library(gridExtra)
```

    ## Warning: package 'gridExtra' was built under R version 3.4.1

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

``` r
prev_plot_grob <- arrangeGrob(prev_plot_gridoff)
ggsave(filename = "Figures/prevalence_plot.png",plot = prev_plot_grob,height = 1.2,width = 3.25,dpi = 600)


mutation_rates_full_scatter <- ggplot(data = mutation_rates, aes(x=positive_mut_rates,y=negative_mut_rates)) +
  geom_point(alpha=0.2,col="red",size=0.5) + 
  geom_smooth(method='lm',color="red",size=.5) + 
  # geom_text_repel(data=subset(mutation_rates, positive_mut_rates > 5e-5 | negative_mut_rates > 3e-5 ), aes(label=gene),size=3) + 
  geom_point(data = subset(mutation_rates, gene %in% c("FBXW7","PIK3CA","MAPK1")),aes(x=positive_mut_rates,y=negative_mut_rates),size=1.5,col="black") + 
  geom_text_repel(data = subset(mutation_rates, gene %in% c("FBXW7","PIK3CA","MAPK1")),aes(x=positive_mut_rates,y=negative_mut_rates,label=gene),size=common.text.size*(5/14),segment.alpha = 1,col="black",box.padding = 0.6) +
  labs(y=bquote("Mutation rate in "~HPV^{"−"}~ "\n tumors"), x=bquote("Mutation rate in "~HPV^{"+"}~ "tumors")) + 
  # coord_equal() + 
  # theme_classic() + 
  expand_limits(x = 0, y = c(0,max(mutation_rates$negative_mut_rates)+1e-6)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_segment(aes(x = 0, xend =max(mutation_rates$negative_mut_rates)+1e-6, y = 0, yend =max(mutation_rates$negative_mut_rates)+1e-6),linetype="dashed") + 
  # geom_abline(slope=1, intercept=0,linetype = "dashed") + 
  scale_x_continuous(labels=fancy_scientific,expand = c(0, 0)) + 
  scale_y_continuous(labels=fancy_scientific,expand = c(0, 0)) + theme(axis.text=element_text(size=common.text.size),
                                                                       axis.title=element_text(size=common.text.size,face="bold"),
                                                                       plot.margin = margin(r=.2,t=.0,b=.0,l=.1,unit = "in")) + theme(axis.title.x = element_text(margin = margin(t = -.2, r = 0, b = 0, l = 0)))  

mutation_rates_full_scatter
```

![](APOBEC_HNSCC_manuscript_analysis_files/figure-markdown_github/prevalence%20and%20mutation%20rate%20plots-1.png)

``` r
mutation_rates_gridoff <- ggplot_gtable(ggplot_build(mutation_rates_full_scatter))
mutation_rates_gridoff$layout$clip[mutation_rates_gridoff$layout$name == "panel"] <- "off"
# library(grid)
# library(gridExtra)

mutation_rates_plot_grob <- arrangeGrob(mutation_rates_gridoff)


ggsave(plot = mutation_rates_plot_grob,filename = "Figures/mutation_rates_full_scatter.png",height = 1.2,width = 3.25,dpi=600)


library(cowplot)

combined.plot.prev.and.muts <- plot_grid(prev_plot,mutation_rates_full_scatter,nrow = 2,align='v')

ggsave(plot = combined.plot.prev.and.muts,filename = "Figures/selection_and_mutrates_combined.png",height = 1.2*2,width = 3.25,dpi=600)

HPV.neg.trinuc.heatmap <- ggplot(data=HPV.neg.trinuc.mutation_data, aes(Downstream, Upstream)) +
  geom_tile(aes(fill = proportion*100), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue", name="Percent") + facet_grid(.~section_labels, labeller = label_parsed)  +  #geom_text(aes(label = round(proportion, 4)*100),size=.8) +
  theme_bw() + theme(panel.grid.major = element_blank(),plot.margin =  margin(r=0,t=5,b=.0,l=0),
                     panel.grid.minor = element_blank(),
                     axis.ticks = element_blank(),
                     strip.text=element_text(size=common.text.size),
                     axis.title.x = element_text(size=common.text.size),
                     axis.title.y = element_text(size=common.text.size),
                     axis.text.x = element_text(size=common.text.size,margin=margin(t=0)),
                     legend.text = element_text(size=common.text.size,margin=margin(l=-2)),
                     legend.title = element_text(size=common.text.size),
                     axis.text.y=element_text(size=common.text.size,margin=margin(r=0)),
                     plot.title = element_text(hjust = 0.5),
                     legend.position = "right") +   
  guides(fill = guide_colorbar(barwidth = .5, barheight = 2)) +
  theme(strip.text.x = element_text(margin = margin(0.4,0,0.4,0))) + 
  theme(legend.margin=margin(t = 0,r=0,b=0,l=-7),legend.box.margin=margin(0,0,0,0)) + 
  theme(panel.spacing = unit(.1, "lines")) + 
  theme(axis.title.x = element_text(margin = margin(t = -.2, r = 0, b = 0, l = 0)))  

HPV.pos.trinuc.heatmap <- ggplot(data=HPV.pos.trinuc.mutation_data, aes(Downstream, Upstream)) +
  geom_tile(aes(fill = proportion*100), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue", name="Percent",labels=c(13,10,7,4,1),breaks=c(13,10,7,4,1)) + facet_grid(.~section_labels, labeller = label_parsed)  +  #geom_text(aes(label = round(proportion, 4)*100),size=.8) +
  theme_bw() + theme(panel.grid.major = element_blank(),plot.margin =  margin(r=0,t=5,b=.0,l=0),
                     panel.grid.minor = element_blank(),
                     axis.ticks = element_blank(),
                     strip.text=element_text(size=common.text.size),
                     axis.title.x = element_text(size=common.text.size),
                     axis.title.y = element_text(size=common.text.size),
                     axis.text.x = element_text(size=common.text.size,margin=margin(t=0)),
                     legend.text = element_text(size=common.text.size,margin=margin(l=-2)),
                     legend.title = element_text(size=common.text.size),
                     axis.text.y=element_text(size=common.text.size,margin=margin(r=0)),
                     plot.title = element_text(hjust = 0.5),
                     legend.position = "right") +   
  guides(fill = guide_colorbar(barwidth = .5, barheight = 2)) +
  theme(strip.text.x = element_text(margin = margin(0.4,0,0.4,0))) +
  theme(legend.margin=margin(t = 0,r=0,b=0,l=-7),legend.box.margin=margin(0,0,0,0)) + 
  theme(panel.spacing = unit(.1, "lines")) +
  theme(axis.title.x = element_text(margin = margin(t = -.2, r = 0, b = 0, l = 0)))  
# HPV.pos.trinuc.heatmap



combined.plot.prev.and.muts.heat <- plot_grid(prev_plot_gridoff,mutation_rates_gridoff,HPV.neg.trinuc.heatmap,HPV.pos.trinuc.heatmap,nrow = 4,align='v',axis='l',rel_heights = c(1,1,.4,.4),labels = "AUTO",label_size = 6)

ggsave(plot = combined.plot.prev.and.muts.heat,filename = "Figures/selection_and_mutrates_heat_combined.png",height = 1*4,width = 3.25,dpi=600)
```

Selection and tornado plots
===========================

``` r
load("input_data/selection_from_cluster/HNSC_HPVpos/selection_output/HNSC_HPVpos_selection_output.RData")
HPV.pos.selectionoutput <- selection.output
HPV.pos.selection_minrecur <- subset(HPV.pos.selectionoutput$all_mutations, freq>1)


load("input_data/selection_from_cluster/HNSC_HPVneg/selection_output/HNSC_HPVneg_selection_output.RData")
HPV.neg.selectionoutput <- selection.output
HPV.neg.selection_minrecur <- subset(HPV.neg.selectionoutput$all_mutations, freq>1)


HPV.pos.selection_minrecur$Name_short <- NA
for(i in 1:nrow(HPV.pos.selection_minrecur)){
  HPV.pos.selection_minrecur$Name_short[i] <- paste(HPV.pos.selection_minrecur$Gene[i]," ",ifelse(!is.na(HPV.pos.selection_minrecur$AA_Ref[i]),paste(HPV.pos.selection_minrecur$AA_Ref[i],HPV.pos.selection_minrecur$AA_Pos[i],HPV.pos.selection_minrecur$AA_Change[i],sep=""),paste(HPV.pos.selection_minrecur$Nuc_Ref[i],HPV.pos.selection_minrecur$Nucleotide_position[i],HPV.pos.selection_minrecur$Nuc_Change[i],"NCSNV")),sep="")
}


HPV.neg.selection_minrecur$Name_short <- NA
for(i in 1:nrow(HPV.neg.selection_minrecur)){
  HPV.neg.selection_minrecur$Name_short[i] <- paste(HPV.neg.selection_minrecur$Gene[i]," ",ifelse(!is.na(HPV.neg.selection_minrecur$AA_Ref[i]),paste(HPV.neg.selection_minrecur$AA_Ref[i],HPV.neg.selection_minrecur$AA_Pos[i],HPV.neg.selection_minrecur$AA_Change[i],sep=""),paste(HPV.neg.selection_minrecur$Nuc_Ref[i],HPV.neg.selection_minrecur$Nucleotide_position[i],HPV.neg.selection_minrecur$Nuc_Change[i],"NCSNV")),sep="")
}


gamma.df <- cbind(prevalence.df,NA,NA)
colnames(gamma.df)[6:7] <- c("gamma_HPVpos","gamma_HPVneg")

for(i in 1:nrow(gamma.df)){
  if(length(which(HPV.pos.selection_minrecur$Name_short == gamma.df$Name[i]))>0){
    gamma.df$gamma_HPVpos[i] <- HPV.pos.selection_minrecur$gamma_epistasis[which(HPV.pos.selection_minrecur$Name_short == gamma.df$Name[i])][1]
  }
  if(length(which(HPV.neg.selection_minrecur$Name_short == gamma.df$Name[i]))>0){
    gamma.df$gamma_HPVneg[i] <- HPV.neg.selection_minrecur$gamma_epistasis[which(HPV.neg.selection_minrecur$Name_short == gamma.df$Name[i])][1]
  }
}

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  l <- gsub("0e\\+00","0",l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}



library(ggrepel)


gamma_plot <- ggplot(data = gamma.df) + geom_point(aes(x = gamma_HPVneg, y = gamma_HPVpos),alpha=0.5,size=2) + geom_text_repel(aes(x = gamma_HPVneg, y = gamma_HPVpos,label=Name),box.padding = 0.8,size=common.text.size*(5/14),color="darkred",fontface="bold") + theme_bw() + labs(x=bquote("Selection intensity in "~HPV^{"−"}~ "tumors"), y=bquote("Selection intensity in "~HPV^{"+"}~ "tumors")) + theme(axis.text=element_text(size=common.text.size), axis.title=element_text(size=common.text.size,face="bold"))  + scale_y_log10(labels = fancy_scientific) + scale_x_log10(labels = fancy_scientific) + coord_equal() #+ coord_equal(xlim = c(0,),ylim = c(0,))

gamma_plot
```

    ## Warning: Removed 293 rows containing missing values (geom_point).

    ## Warning: Removed 293 rows containing missing values (geom_text_repel).

![](APOBEC_HNSCC_manuscript_analysis_files/figure-markdown_github/gamma%20gamma%20plot-1.png)

``` r
ggsave(filename = "Figures/gamma_gamma_plot.png",plot = gamma_plot,height = 2,width = 2,dpi=600)
```

    ## Warning: Removed 293 rows containing missing values (geom_point).

    ## Warning: Removed 293 rows containing missing values (geom_text_repel).

tornado plots
-------------

``` r
### Figures for manuscript

# Selection plots with the same colors 
# setwd("~/Documents/Selection_analysis/HNSC/NCI/")
# source("~/GitHub/site_specific_selection_intensity/functions_script.R")

# load in the HPV neg, load in the HPV pos, find all unique genes and generate color palette 
# This will be sorted by selection, not frequency. 
# Neg 
load("input_data/selection_from_cluster/HNSC_HPVneg/selection_output/HNSC_HPVneg_selection_output.RData")

selection.subset.neg <- selection.output$all_mutations[which(selection.output$all_mutations$freq>1),]
selection.subset.neg <- selection.subset.neg[which(!is.na(selection.subset.neg$Gene)),]
selection.subset.neg <- selection.subset.neg[order(-selection.subset.neg$gamma_epistasis),]

if(nrow(selection.subset.neg)>25){
  selection.subset.neg.ordered <- selection.subset.neg[1:25,] #
}else{
  selection.subset.neg.ordered <- selection.subset.neg 
}
selection.subset.neg.ordered <- selection.subset.neg.ordered[order(selection.subset.neg.ordered$gamma_epistasis),]
selection.subset.neg.ordered$Name <- NA

for(i in 1:nrow(selection.subset.neg.ordered)){
  selection.subset.neg.ordered$Name[i] <- paste(selection.subset.neg.ordered$Gene[i]," ",ifelse(!is.na(selection.subset.neg.ordered$AA_Ref[i]),paste(selection.subset.neg.ordered$AA_Ref[i],selection.subset.neg.ordered$AA_Pos[i],selection.subset.neg.ordered$AA_Change[i],sep=""),"NCSNV"))
}



# selection.subset.neg.ordered$Name <- paste(selection.subset.neg.ordered$Gene," ",selection.subset.neg.ordered$AA_Ref,selection.subset.neg.ordered$AA_Pos,selection.subset.neg.ordered$AA_Change,sep="")

length(unique(selection.subset.neg.ordered$Name))
```

    ## [1] 25

``` r
if(length(which(selection.subset.neg.ordered$Name=="TP53   NCSNV"))>1){
  selection.subset.neg.ordered$Name[which(selection.subset.neg.ordered$Name=="TP53   NCSNV")] <- paste("TP53    NCSNV",letters[length((which(selection.subset.neg.ordered$Name=="TP53   NCSNV"))):1],sep="")
}
selection.subset.neg.ordered$Name <- factor(selection.subset.neg.ordered$Name, levels=selection.subset.neg.ordered$Name)
# selection.subset.neg.ordered
```

``` r
# Pos

load("input_data/selection_from_cluster/HNSC_HPVpos/selection_output/HNSC_HPVpos_selection_output.RData")
selection.subset.pos <- selection.output$all_mutations[which(selection.output$all_mutations$freq>1),]
selection.subset.pos <- selection.subset.pos[which(!is.na(selection.subset.pos$Gene)),]
selection.subset.pos <- selection.subset.pos[order(-selection.subset.pos$gamma_epistasis),]

if(nrow(selection.subset.pos)>25){
  selection.subset.pos.ordered <- selection.subset.pos[1:25,] #
}else{
  selection.subset.pos.ordered <- selection.subset.pos 
}
selection.subset.pos.ordered <- selection.subset.pos.ordered[order(selection.subset.pos.ordered$gamma_epistasis),]
selection.subset.pos.ordered$Name <- NA

for(i in 1:nrow(selection.subset.pos.ordered)){
  selection.subset.pos.ordered$Name[i] <- paste(selection.subset.pos.ordered$Gene[i]," ",ifelse(!is.na(selection.subset.pos.ordered$AA_Ref[i]),paste(selection.subset.pos.ordered$AA_Ref[i],selection.subset.pos.ordered$AA_Pos[i],selection.subset.pos.ordered$AA_Change[i],sep=""),"NCSNV"))
}



# selection.subset.pos.ordered$Name <- paste(selection.subset.pos.ordered$Gene," ",selection.subset.pos.ordered$AA_Ref,selection.subset.pos.ordered$AA_Pos,selection.subset.pos.ordered$AA_Change,sep="")
selection.subset.pos.ordered$Name <- factor(selection.subset.pos.ordered$Name, levels=selection.subset.pos.ordered$Name)
# selection.subset.pos.ordered
```

``` r
##Find colors used for unique genes in combined build. 
source("R/fancy_scientific_code.R")
selection.subset.combined.ordered <- rbind(selection.subset.neg.ordered,selection.subset.pos.ordered)

unique.genes <- unique(selection.subset.combined.ordered$Gene)
#Make a dataframe of just the unique genes
selection.subset.combined.ordered.unique <- as.data.frame(matrix(nrow=0,ncol=ncol(selection.subset.combined.ordered)))
colnames(selection.subset.combined.ordered.unique) <- colnames(selection.subset.combined.ordered)
for(i in 1:length(unique.genes)){
  selection.subset.combined.ordered.unique <- rbind(selection.subset.combined.ordered.unique,selection.subset.combined.ordered[which(selection.subset.combined.ordered$Gene==unique.genes[i])[1],])
}

#from http://stackoverflow.com/questions/18265941/two-horizontal-bar-charts-with-shared-axis-in-ggplot2-similar-to-population-pyr
library('grid')
library('gridExtra')


g.mid <- ggplot(selection.subset.combined.ordered.unique,aes(x=1,y=Name)) +
  geom_text(aes(label=Name),size=7) +
  # geom_segment(aes(x=0.94,xend=0.96,yend=Name)) +
  # geom_segment(aes(x=1.04,xend=1.065,yend=Name)) +
  ggtitle("") +
  ylab(NULL) + 
  scale_x_continuous(expand=c(0,0),limits=c(0.94,1.065)) + ggtitle("HNSCC HPV+") +
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(1,-1,1,-1), "mm"))
# g.mid

g1 <- ggplot(data=selection.subset.combined.ordered.unique,aes(x=Name,y=mu, fill=Gene)) +
  geom_bar(stat="identity") + ggtitle("Mutation rate") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        plot.margin = unit(c(1,-1,1,10), "mm")) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major =element_line(color="lightgrey"),panel.grid.minor =element_line(color="lightgrey")) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) +
  geom_text(aes(label=round(mu*1e5,2)), vjust=-0.5, hjust=0.5, position=position_dodge(width=0.9),size=4,angle=90) +
  geom_text(aes(label=freq,y=-2e-7), position=position_dodge(width=0.9),size=6,angle=0,color="black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = 'none',axis.text.x = element_text(size=12)) +
  scale_y_reverse(labels=fancy_scientific) + coord_flip() 
# g1

g2 <- ggplot(data=selection.subset.combined.ordered.unique, aes(x=Name,y=gamma_epistasis,fill=Gene)) +
  geom_bar(stat="identity") + ggtitle("Selection intensity") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = unit(c(1,2,1,-1), "mm")) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major =element_line(color="lightgrey"),panel.grid.minor =element_line(color="lightgrey")) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size=12)) +
  coord_flip() + theme(legend.position = c(0.85, .5),legend.text = element_text(size=10)) + theme(legend.position="none")
# g2

gg1 <- ggplot_gtable(ggplot_build(g1))
gg2 <- ggplot_gtable(ggplot_build(g2))
gg.mid <- ggplot_gtable(ggplot_build(g.mid))

# grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(4.25/10,2/10,3.5/10))


gg.combined <- arrangeGrob(gg1,gg.mid,gg2,ncol=3,widths=c(5/10,3.5/10,5/10))

# ggplot_build(g1)$data
# ggplot_build(g.mid)$data

color.data <- ggplot_build(g1)$data 

gene.colors <- selection.subset.combined.ordered.unique$Gene
names(gene.colors) <- color.data[[1]]$fill


color.df <- as.data.frame(matrix(data = NA, nrow=length(unique(selection.subset.combined.ordered.unique$Gene)),ncol=2))
colnames(color.df) <- c("gene","color")
color.df$gene <- unique(selection.subset.combined.ordered.unique$Gene)
for(i in 1:nrow(selection.subset.combined.ordered.unique)){
  if(is.na(color.df$color[which(color.df$gene==selection.subset.combined.ordered.unique$Gene[i])])){
    color.df$color[which(color.df$gene==selection.subset.combined.ordered.unique$Gene[i])] <- 
      color.data[[1]]$fill[i]
  }
}

# color.df

color.vec <- color.df$color
names(color.vec) <- color.df$gene
```

``` r
### Now, for the separate plots
# HPV negative 
library(ggplot2)
library(cowplot)


neg.colors <- color.vec[names(color.vec) %in% selection.subset.neg.ordered$Gene]

g.mid <- ggplot(selection.subset.neg.ordered,aes(x=1,y=Name)) +
  geom_text(aes(label=Name),size=common.text.size*(5/14)) +
  # geom_segment(aes(x=0.94,xend=0.96,yend=Name)) +
  # geom_segment(aes(x=1.04,xend=1.065,yend=Name)) +
  ggtitle("") +
  ylab(NULL) + 
  # scale_x_continuous(expand=c(0,0),limits=c(0.94,1.065)) + ggtitle("HNSCC HPV-") +
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),axis.line=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(1,-1,1,-1), "mm")) 
# g.mid

g1 <- ggplot(data=selection.subset.neg.ordered,aes(x=Name,y=mu, fill=Gene)) +
  geom_bar(stat="identity") + ggtitle("Mutation rate") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), axis.line=element_blank(),
        axis.ticks.y = element_blank(), panel.border = element_blank(),
        plot.margin = unit(c(1,-1,1,3.7), "mm")) +
  theme(panel.background = element_blank()) + scale_fill_manual("Legend", values = neg.colors) +
  theme(panel.grid.major =element_line(color="lightgrey"),panel.grid.minor =element_line(color="lightgrey")) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) +
  # geom_text(aes(label=round(mu*1e5,2)), vjust=-0.5, hjust=0.5, position=position_dodge(width=0.9),size=4,angle=90) +
  geom_text(aes(label=freq,y=-8e-8), position=position_dodge(width=0.9),size=common.text.size*(5/14),angle=0,color="black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = 'none',axis.text.x = element_text(size=common.text.size,angle = 45, hjust = 1)) +
  scale_y_reverse(labels=fancy_scientific) + coord_flip() + 
  theme(axis.text.x = element_text(size = common.text.size))+ theme(plot.title = element_text(size=common.text.size))
# g1

g2 <- ggplot(data=selection.subset.neg.ordered, aes(x=Name,y=gamma_epistasis,fill=Gene)) +
  geom_bar(stat="identity") + ggtitle("Selection intensity") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line=element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),panel.border = element_blank(),
        plot.margin = unit(c(1,-1,1,-1), "mm")) +
  theme(panel.background = element_blank()) + scale_fill_manual("Legend", values = neg.colors) +
  theme(panel.grid.major =element_line(color="lightgrey"),panel.grid.minor =element_line(color="lightgrey")) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size=common.text.size,angle = 45, hjust = 1)) +
  coord_flip() + theme(legend.position = c(0.85, .5),legend.text = element_text(size=common.text.size)) + theme(legend.position="none") + scale_y_continuous(labels=fancy_scientific,breaks=c(0,1e5,2e5,3e5),minor_breaks = seq(0,3e5,5e4)) + theme(axis.text.x = element_text(size =common.text.size)) + theme(plot.title = element_text(size=common.text.size))
# g2

# gg1 <- ggplot_gtable(ggplot_build(g1))
# gg2 <- ggplot_gtable(ggplot_build(g2))
# gg.mid <- ggplot_gtable(ggplot_build(g.mid))
# library(grid);library(gridExtra)
# grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(3/10,3.5/10,3/10))

library(cowplot)

combined_neg_gamma_plot <- plot_grid(g1,g.mid,g2,align = 'h',axis="t",nrow=1,rel_widths = c(1,.7,1))
combined_neg_gamma_plot
```

![](APOBEC_HNSCC_manuscript_analysis_files/figure-markdown_github/hpv%20neg%20tornado%20plot-1.png)

``` r
ggsave(filename = "Figures/combined_neg_gamma.png",plot = combined_neg_gamma_plot,width = 3.25,dpi = 600,height = 2.8)


# gg.combined <- arrangeGrob(gg1,gg.mid,gg2,ncol=3,widths=c(3/10,2/10,3/10))
# ggsave(gg.combined, filename = paste("Figures/same_color_gamma_epistasis_plot_","HNSC_HPV-.png",sep=""),units = "in",height=8,width = 12,dpi = 400)
```

``` r
# HPV positive



pos.colors <- color.vec[names(color.vec) %in% selection.subset.pos.ordered$Gene]

g.mid.pos <- ggplot(selection.subset.pos.ordered,aes(x=1,y=Name)) +
  geom_text(aes(label=Name),size=common.text.size*(5/14)) +
  # geom_segment(aes(x=0.94,xend=0.96,yend=Name)) +
  # geom_segment(aes(x=1.04,xend=1.065,yend=Name)) +
  ggtitle("") +
  ylab(NULL) + 
  # scale_x_continuous(expand=c(0,0),limits=c(0.94,1.065)) + ggtitle("HNSCC HPV-") +
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),axis.line=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(1,-1,1,-1), "mm")) 
# g.mid

g1.pos <- ggplot(data=selection.subset.pos.ordered,aes(x=Name,y=mu, fill=Gene)) +
  geom_bar(stat="identity") + ggtitle("Mutation rate") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), axis.line=element_blank(),
        axis.ticks.y = element_blank(), panel.border = element_blank(),
        plot.margin = unit(c(1,-1,1,1), "mm")) +
  theme(panel.background = element_blank()) + scale_fill_manual("Legend", values = pos.colors) +
  theme(panel.grid.major =element_line(color="lightgrey"),panel.grid.minor =element_line(color="lightgrey")) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) +
  # geom_text(aes(label=round(mu*1e5,2)), vjust=-0.5, hjust=0.5, position=position_dodge(width=0.9),size=4,angle=90) +
  geom_text(aes(label=freq,y=-3.5e-7), position=position_dodge(width=0.9),size=common.text.size*(5/14),angle=0,color="black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = 'none',axis.text.x = element_text(size=common.text.size,angle = 45, hjust = 1)) +
  scale_y_reverse(labels=fancy_scientific) + coord_flip() + 
  theme(axis.text.x = element_text(size = common.text.size))+ theme(plot.title = element_text(size=common.text.size))
# g1

g2.pos <- ggplot(data=selection.subset.pos.ordered, aes(x=Name,y=gamma_epistasis,fill=Gene)) +
  geom_bar(stat="identity") + ggtitle("Selection intensity") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line=element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),panel.border = element_blank(),
        plot.margin = unit(c(1,r=1,1,-1), "mm")) +
  theme(panel.background = element_blank()) + scale_fill_manual("Legend", values = pos.colors) +
  theme(panel.grid.major =element_line(color="lightgrey"),panel.grid.minor =element_line(color="lightgrey")) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size=common.text.size,angle = 45, hjust = 1)) +
  coord_flip() + theme(legend.position = c(0.85, .5),legend.text = element_text(size=common.text.size)) + theme(legend.position="none") + scale_y_continuous(labels=fancy_scientific,breaks=c(0,2.5e5,5e5,7.5e5,1e6,1.2e6)) + theme(axis.text.x = element_text(size =common.text.size)) + theme(plot.title = element_text(size=common.text.size))
# g2



library(cowplot)

combined_pos_gamma_plot <- plot_grid(g1.pos,g.mid.pos,g2.pos,align = 'h',axis="t",nrow=1,rel_widths = c(1,.7,1))
combined_pos_gamma_plot
```

![](APOBEC_HNSCC_manuscript_analysis_files/figure-markdown_github/hpv%20pos%20tornado%20plot-1.png)

``` r
ggsave(filename = "Figures/combined_pos_gamma.png",plot = combined_pos_gamma_plot,width = 3.25,dpi = 600,height = 2.8)

combined_both_gamma <- plot_grid(g1,g.mid,g2,g1.pos,g.mid.pos,g2.pos,align = 'h',axis="t",nrow=1,rel_widths = c(1,.6,1,1,.6,1))
ggsave(filename = "Figures/combined_both_gamma.png",plot = combined_both_gamma,width = 3.25*2,dpi = 600,height = 2.8)

#
```

``` r
library(ggplot2);library(ggrepel)
gamma_plot <- ggplot(data = gamma.df) + 
  geom_point(aes(x = gamma_HPVneg, y = gamma_HPVpos),alpha=0.5,size=1.5,col="red") + 
  geom_text_repel(aes(x = gamma_HPVneg, y = gamma_HPVpos,label=Name),
                  box.padding =.7,
                  size=common.text.size*(5/14),segment.alpha=0.2,
                  color="black") + 
  labs(x=bquote("Selection intensity in "~HPV^{"−"}~ "tumors"), y=bquote("Selection intensity in "~HPV^{"+"}~ "tumors")) +
  theme(axis.text=element_text(size=common.text.size), 
        axis.title=element_text(size=common.text.size,face="bold"))  + 
  expand_limits(x = c(1,(max(gamma.df$gamma_HPVpos,na.rm = T)+2e7)), y = c(1,(max(gamma.df$gamma_HPVpos,na.rm = T)+2e7))) + 
  # scale_x_continuous(expand = c(0, 0)) + 
  # scale_y_continuous(expand = c(0, 0)) + 
  scale_y_log10(labels = c(expression(10^0),expression(10^2),expression(10^4),expression(10^6)),breaks=c(1,1e2,1e4,1e6),expand = c(0,0)) +
  scale_x_log10(labels = c(expression(10^0),expression(10^2),expression(10^4),expression(10^6)),breaks=c(1,1e2,1e4,1e6),expand = c(0,0)) + coord_equal() #+
# coord_equal() #+ coord_equal(xlim = c(0,),ylim = c(0,))
gamma_plot
```

    ## Warning: Removed 293 rows containing missing values (geom_point).

    ## Warning: Removed 293 rows containing missing values (geom_text_repel).

![](APOBEC_HNSCC_manuscript_analysis_files/figure-markdown_github/adding%20gamma%20gamma%20plot%20to%20gamma%20tornado%20plots-1.png)

``` r
library(viridis)
```

    ## Loading required package: viridisLite

``` r
library(cowplot)
# plot_grid(combined_both_gamma,NULL, nrow=1,ncol=2,align = 'h',rel_widths = c(1,.1))
combined_both_gamma_forcombined <- plot_grid(g1,g.mid,g2,g1.pos,g.mid.pos,g2.pos,NULL,align = 'h',axis="t",nrow=1,rel_widths = c(1,.7,1,1,.7,1,1),ncol=7)
Fig2 <- ggdraw() + 
  draw_plot(combined_both_gamma_forcombined) +
  geom_rect(data=data.frame(x=0.715,y=0.23),aes(xmin=x,xmax=x+.28,ymin=y,ymax=y+.52),color="white",fill="white") + 
  draw_plot(gamma_plot,  x=0.57, y=0.19,width = .55,height = .55) + 
  draw_plot_label(c("A", "B","C"), c(0, 0.43,.715), c(.985, .985,(.23+.52)), size = common.text.size)
```

    ## Warning: Removed 293 rows containing missing values (geom_point).

    ## Warning: Removed 293 rows containing missing values (geom_text_repel).

``` r
Fig2
```

![](APOBEC_HNSCC_manuscript_analysis_files/figure-markdown_github/adding%20gamma%20gamma%20plot%20to%20gamma%20tornado%20plots-2.png)

``` r
# Fig2
ggsave(filename = "Figures/Fig4_tornadoplots.png",dpi = 600,width = 3.25*2,height =3,plot = Fig2)
```

heatmap and dendrogram of the signatures
========================================

``` r
# library(plyr)
library(reshape2)
# library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggdendro)
library(gridExtra)
library(dendextend)
```

    ## Warning: package 'dendextend' was built under R version 3.4.2

    ## 
    ## ---------------------
    ## Welcome to dendextend version 1.6.0
    ## Type citation('dendextend') for how to cite the package.
    ## 
    ## Type browseVignettes(package = 'dendextend') for the package vignette.
    ## The github page is: https://github.com/talgalili/dendextend/
    ## 
    ## Suggestions and bug-reports can be submitted at: https://github.com/talgalili/dendextend/issues
    ## Or contact: <tal.galili@gmail.com>
    ## 
    ##  To suppress this message use:  suppressPackageStartupMessages(library(dendextend))
    ## ---------------------

    ## 
    ## Attaching package: 'dendextend'

    ## The following object is masked from 'package:ggdendro':
    ## 
    ##     theme_dendro

    ## The following object is masked from 'package:stats':
    ## 
    ##     cutree

``` r
# useful: https://stackoverflow.com/questions/42047896/joining-a-dendrogram-and-a-heatmap
# and https://plot.ly/ggplot2/ggdendro-dendrograms/


load("output_data/trinuc_contexts_HNSC_MAF.RData")
weights.df <- trinuc.contexts$signature.weights[[1]]$weights
for(i in 2:length(trinuc.contexts$signature.weights)){
  weights.df <- rbind(weights.df,trinuc.contexts$signature.weights[[i]]$weights)
}


# prevalence of weights
weight.prevalence <- rep(NA,30)

for(i in 1:length(weight.prevalence)){
  weight.prevalence[i] <- length(which(weights.df[,i]>0))
}

message("Order of signatures in decreasing prevalence")
```

    ## Order of signatures in decreasing prevalence

``` r
order(weight.prevalence,decreasing = T)
```

    ##  [1]  1 13  2  7  6  3 16  4 24  5 15 10 18 19 20 30 22 11 29 25 12 23  8
    ## [24] 17 21 26 14 27  9 28

``` r
# Dendrogram data
# not scaling because all data should be on the same scale

colnames(weights.df) <- 1:30

dend <- as.dendrogram(hclust(dist(weights.df)))
dend_data <- dendro_data(dend)

dend_t <- as.dendrogram(hclust(dist(t(weights.df))))
dend_data_t <- dendro_data(dend_t)

# Setup segment data

segment_data_tumors <- with(
  segment(dend_data), 
  data.frame(x = y, y = x, xend = yend, yend = xend))

segment_data_signatures <- with(
  segment(dend_data_t), 
  data.frame(x = y, y = x, xend = yend, yend = xend))

# position data
tumor_pos_table <- with(
  dend_data$labels,
  data.frame(ycenter = x, tumor=as.character(label), height=1))

signature_pos_table <- with(
  dend_data_t$labels,
  data.frame(xcenter = x, signature=as.character(label), width=1))

tumor_names <- rownames(weights.df)
signature_names <- colnames(weights.df)




heatmap_data <- weights.df %>%
  add_column(id = rownames(weights.df)) %>%
  reshape2::melt(value.name = "Signature weight") %>%
  dplyr::rename(tumor = id, signature = variable) %>%
  left_join(signature_pos_table) %>%
  left_join(tumor_pos_table)
```

    ## Using id as id variables

    ## Joining, by = "signature"

    ## Warning: Column `signature` joining factors with different levels, coercing
    ## to character vector

    ## Joining, by = "tumor"

    ## Warning: Column `tumor` joining character vector and factor, coercing into
    ## character vector

``` r
# Limits for the vertical axes
tumor_axis_limits <- with(
  tumor_pos_table, 
  c(min(ycenter - 0.5 * height), max(ycenter + 0.5 * height))
) + 
  0.1 * c(-1, 1) # extra spacing: 0.1
signature_axis_limits <- with(
  signature_pos_table, 
  c(min(xcenter - 0.5 * width), max(xcenter + 0.5 * width))
) + 
  0.1 * c(-1, 1) # extra spacing: 0.1

# ADD HPV status
segment_data_tumors$color <- "black"
# length(which(segment_data_tumors$xend==0));nrow(tumor_pos_table)

# segment_data_tumors$y[which(segment_data_tumors$xend==0)]
segment_data_tumors$tumor <- NA
segment_data_tumors$tumor[which(segment_data_tumors$xend==0)] <- as.character(tumor_pos_table$tumor)
segment_data_tumors$HPV_status <- NA



for(i in 1:nrow(segment_data_tumors[which(segment_data_tumors$xend==0),])){
  segment_data_tumors$HPV_status[which(segment_data_tumors$xend==0)][i] <- HNSC.selection.output$HPV_call[which(HNSC.selection.output$Unique_patient_identifier==segment_data_tumors$tumor[which(segment_data_tumors$xend==0)][i])[1]]
}

# which(is.na(segment_data_tumors$color) & !is.na(segment_data_tumors$tumor)) #all accounted for

# head(segment_data_tumors,10)

for(i in 1:nrow(segment_data_tumors)){
  if(!is.na(segment_data_tumors$HPV_status[i])){
    
    if(segment_data_tumors$HPV_status[i]=="ambiguous"){
      segment_data_tumors$color[i] <- "gray" 
    }
    if(segment_data_tumors$HPV_status[i]=="HPV+"){
      segment_data_tumors$color[i] <- "blue" 
    }
    if(segment_data_tumors$HPV_status[i]=="HPV−"){
      segment_data_tumors$color[i] <- "orange" 
    }
  }
  
}
```

``` r
dendro_size <- .12
# Heatmap plot
plt_hmap <- ggplot(heatmap_data, 
                   aes(x = xcenter, y = ycenter, fill = `Signature weight`, 
                       height = height, width = width)) + 
  geom_tile() +
  scale_fill_gradient2("Signature weight", high = "steelblue", low = "white") +
  scale_x_continuous(breaks = signature_pos_table$xcenter, 
                     labels = signature_pos_table$signature, 
                     expand = c(0, 0),position = "top") + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  scale_y_continuous(breaks = tumor_pos_table[, "ycenter"], 
                     labels = rep("", nrow(tumor_pos_table)),
                     limits = tumor_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "Signature", y = "") +
  theme_bw() +
  theme(axis.text.x.top = element_text(size = common.text.size, vjust=0.5,hjust=0.5 ,  angle = 90),axis.title.y = element_blank(),
        # margin: top, right, bottom, and left
        plot.margin = unit(c(-.051, -.02, -.051, -.051), "cm"),
        # plot.margin = unit(c(-1, 0.2, 0.2, 2), "cm"),
        panel.grid.minor = element_blank(),legend.position = "none",axis.title = element_blank()) 

# try and make a heatmap that is just HPV status... 

segment_data_subset <- subset(segment_data_tumors,color!="black")
segment_data_subset$x <- 1
plt_dendr_subset <- ggplot(segment_data_subset) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend,color=color),size=dendro_size) + 
  scale_x_reverse(expand = c(0, 0.01)) + 
  scale_y_continuous(breaks = tumor_pos_table$y_center, 
                     labels = tumor_pos_table$tumor, 
                     limits = tumor_axis_limits, 
                     expand = c(0, 0)) + scale_colour_identity() +
  labs(x = "", y = "Tumors", colour = "", size = "") + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(),axis.title.x  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),panel.border = element_blank(),panel.grid.major = element_blank(),
        plot.margin = unit(c(-.1, -.1, -.1, -.1), "cm"))



# Dendrogram plot
plt_dendr <- ggplot(segment_data_tumors) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend,color=color),size=dendro_size) + 
  scale_x_reverse(expand = c(0, 0.01)) + 
  scale_y_continuous(breaks = tumor_pos_table$y_center, 
                     labels = tumor_pos_table$tumor, 
                     limits = tumor_axis_limits, 
                     expand = c(0, 0)) + scale_colour_identity() +
  labs(x = "", y = "Tumors", colour = "", size = "") + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(),axis.title.x  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),panel.border = element_blank(),panel.grid.major = element_blank(),plot.margin = unit(c(-.1, -.1, -.1, -.1), "cm"))

plt_dendr_signature <- ggplot(segment_data_signatures) + 
  geom_segment(aes(x = xend, y = yend, xend = x, yend = y)) + 
  # scale_x_reverse(expand = c(0, 0.5)) + 
  scale_y_continuous(breaks = signature_pos_table$y_center, 
                     labels = signature_pos_table$signature, 
                     limits = signature_axis_limits, 
                     expand = c(0, 0),position="top") +
  labs(x = "Distance", y = "Signatures", colour = "", size = "") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),axis.text = element_blank(),axis.ticks.y = element_blank(),axis.title = element_blank(),panel.border = element_blank(),panel.grid.major = element_blank(),plot.margin = unit(c(-.1, -.1, -.1, -.1), "cm")) + coord_flip()

library(cowplot)
heatmap_and_dendro <- plot_grid(NULL,plt_dendr_signature,plt_dendr, plt_hmap, align = 'hv', rel_widths = c(0.2, 1), rel_heights = c(.4,2))

heatmap_and_dendro_wbar <- plot_grid(NULL,NULL,plt_dendr_signature,plt_dendr,plt_dendr_subset, plt_hmap, align = 'hv', rel_widths = c(0.2,.05, 1), rel_heights = c(.4,2))
# heatmap_and_dendro



ggsave(filename = "Figures/heatmap_and_dendro.png",plot = heatmap_and_dendro,height = 5,width = 7)
# ggsave(filename = "Figures/heatmap_and_dendro.pdf",plot = heatmap_and_dendro,height = 9,width = 6)
# 
# 

dendro_and_heatmap_labs <- ggdraw() + 
  draw_plot(plot_grid(NULL,heatmap_and_dendro,rel_widths = c(.05,1),nrow=1)) +
  geom_text(x=.72,y=.95,label="Signatures",size=common.text.size*(5/14)) + 
  geom_text(x=.03,y=.42,label="Tumors",size=common.text.size*(5/14),angle=90)

dendro_and_heatmap_labs
```

![](APOBEC_HNSCC_manuscript_analysis_files/figure-markdown_github/generate%20heatmap-1.png)

``` r
ggsave(filename = "Figures/heatmap_and_dendro_w_labs.png",plot = dendro_and_heatmap_labs,height = 4,width = 6.5)

# plt_dendr_axis <- axis_canvas(plt_hmap,axis = 'y') + ggplot(segment_data_tumors) + 
#     geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
#     scale_x_reverse(expand = c(0, 0.01)) + 
#     scale_y_continuous(breaks = tumor_pos_table$y_center, 
#                        labels = tumor_pos_table$gene, 
#                        limits = tumor_axis_limits, 
#                        expand = c(0, 0)) + 
#     labs(x = "Distance", y = "", colour = "", size = "") +
#     theme_bw() + 
#     theme(panel.grid.minor = element_blank(),axis.title.x  = element_blank(),plot.margin = unit(c(0, 0.0, 0.0, -0.7), "cm"))
```

``` r
library(cowplot)
# combined.plot.prev.and.muts.heat

# dendro_and_heatmap_labs

combined_prev_mut_heat_dendro <- plot_grid(combined.plot.prev.and.muts.heat,NULL,nrow = 2,ncol=1,rel_heights = c(1,.6),labels = c("","E"),label_size = common.text.size,rel_widths = c(1,5))

# save_plot(combined_prev_mut_heat_dendro,filename = "Figures/prev_mut_heat_dendro.png",base_height =  8,base_width  = 3.25,dpi=600)

combined_prev_mut_heat_dendro2 <- ggdraw(combined_prev_mut_heat_dendro) + draw_plot(heatmap_and_dendro,y = 0,x = .05,height = .37,width = .94) +
  geom_text(x=.55,y=.35,label="Signatures",size=common.text.size*(5/14)) + 
  geom_text(x=.05,y=.17,label="Tumors",size=common.text.size*(5/14),angle=90)

combined_prev_mut_heat_dendro2
```

![](APOBEC_HNSCC_manuscript_analysis_files/figure-markdown_github/combining%20prevalence%20mutation%20and%20heatmaps-1.png)

``` r
save_plot(combined_prev_mut_heat_dendro2,filename = "Figures/Fig1_prev_mut_heat_dendro.png",base_height =  8,base_width  = 3.25,dpi=600)
```

HPV vs SNP count
================

``` r
library(ggplot2)
load("output_data/HNSC_selection_with_APOBEC.RData")

HNSC.selection.HPVknown <- HNSC.selection.output[which(HNSC.selection.output$HPV_call != "ambiguous"),] #also removes NA


SNV.APOBEC.df <- as.data.frame(matrix(data = NA, nrow=length(unique(HNSC.selection.HPVknown$Unique_patient_identifier)),ncol=5))
colnames(SNV.APOBEC.df) <- c("Tumor","SNV_count","APOBEC","APOBEC_weight","HPV_status")

SNV.APOBEC.df$Tumor <- unique(HNSC.selection.HPVknown$Unique_patient_identifier)

for(i in 1:nrow(SNV.APOBEC.df)){
  SNV.APOBEC.df$SNV_count[i] <- length(which(HNSC.selection.HPVknown$Unique_patient_identifier==SNV.APOBEC.df$Tumor[i]))
  if(!is.na(HNSC.selection.HPVknown$APOBEC_weight[which(HNSC.selection.HPVknown$Unique_patient_identifier==SNV.APOBEC.df$Tumor[i])[1]])){
    SNV.APOBEC.df$APOBEC[i] <- if(HNSC.selection.HPVknown$APOBEC_weight[which(HNSC.selection.HPVknown$Unique_patient_identifier==SNV.APOBEC.df$Tumor[i])[1]]>0){"Yes"
    }else{if(HNSC.selection.HPVknown$APOBEC_weight[which(HNSC.selection.HPVknown$Unique_patient_identifier==SNV.APOBEC.df$Tumor[i])[1]]==0){"No"}}}
  SNV.APOBEC.df$HPV_status[i] <- HNSC.selection.HPVknown$HPV_call[which(HNSC.selection.HPVknown$Unique_patient_identifier==SNV.APOBEC.df$Tumor[i])[1]]
  SNV.APOBEC.df$APOBEC_weight[i] <- HNSC.selection.HPVknown$APOBEC_weight[which(HNSC.selection.HPVknown$Unique_patient_identifier==SNV.APOBEC.df$Tumor[i])[1]]
}

# View(SNV.APOBEC.df)

write.table(x = SNV.APOBEC.df,file = "output_data/SNV_APOBEC_HPV_status.txt",sep="\t",quote = F,row.names = F)

source("R/fancy_scientific_code.R")

HPV.vs.APOBEC <- ggplot(data = SNV.APOBEC.df) + geom_boxplot(aes(y=SNV_count,x=APOBEC),color="dark red") + geom_jitter(aes(y=SNV_count,x=APOBEC),width= 0.2,alpha=0.5) + facet_grid(.~HPV_status) + theme_bw() + scale_y_log10(labels=fancy_scientific) + labs(x="APOBEC signature detected",y="SNV count")


ggsave(filename = "Figures/HPV_status_and_APOBEC_vs_SNV.png",plot = HPV.vs.APOBEC)
```

    ## Saving 7 x 5 in image

``` r
wilcox.test(data=subset(SNV.APOBEC.df,HPV_status=="HPV−"),SNV_count~APOBEC,conf.int=T) # SNV count vs. APOBEC signature in HPV - tumors
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  SNV_count by APOBEC
    ## W = 17816, p-value = 0.6309
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -17.00006  11.99994
    ## sample estimates:
    ## difference in location 
    ##              -3.999963

``` r
wilcox.test(data=subset(SNV.APOBEC.df,HPV_status=="HPV+"),SNV_count~APOBEC, conf.int = T,alternative="two.sided") # SNV count vs. APOBEC signature in HPV - tumors
```

    ## Warning in wilcox.test.default(x = c(74L, 75L, 141L, 119L, 137L, 73L,
    ## 64L, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(x = c(74L, 75L, 141L, 119L, 137L, 73L,
    ## 64L, : cannot compute exact confidence intervals with ties

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  SNV_count by APOBEC
    ## W = 76.5, p-value = 0.02531
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -140.999963   -4.000086
    ## sample estimates:
    ## difference in location 
    ##              -53.99995

``` r
hpv.pos.snv.apobec <- subset(SNV.APOBEC.df,HPV_status=="HPV+") 
hpv.neg.snv.apobec <- subset(SNV.APOBEC.df,HPV_status=="HPV−")

mean(hpv.pos.snv.apobec$SNV_count[which(hpv.pos.snv.apobec$APOBEC=="Yes")])-mean(hpv.pos.snv.apobec$SNV_count[which(hpv.pos.snv.apobec$APOBEC=="No")]) #SNV count difference in mean
```

    ## [1] 105.359

``` r
median(hpv.pos.snv.apobec$SNV_count[which(hpv.pos.snv.apobec$APOBEC=="Yes")])-median(hpv.pos.snv.apobec$SNV_count[which(hpv.pos.snv.apobec$APOBEC=="No")]) #SNV count difference in median
```

    ## [1] 67.5

``` r
# fold difference of median 
median(hpv.pos.snv.apobec$SNV_count[which(hpv.pos.snv.apobec$APOBEC=="Yes")])/median(hpv.pos.snv.apobec$SNV_count[which(hpv.pos.snv.apobec$APOBEC=="No")])
```

    ## [1] 1.90604

``` r
mean(hpv.neg.snv.apobec$SNV_count[which(hpv.neg.snv.apobec$APOBEC=="Yes")])-mean(hpv.neg.snv.apobec$SNV_count[which(hpv.neg.snv.apobec$APOBEC=="No")])
```

    ## [1] -79.35159

``` r
median(hpv.neg.snv.apobec$SNV_count[which(hpv.neg.snv.apobec$APOBEC=="Yes")])-median(hpv.neg.snv.apobec$SNV_count[which(hpv.neg.snv.apobec$APOBEC=="No")])
```

    ## [1] 9

``` r
summary(hpv.pos.snv.apobec$SNV_count[which(hpv.pos.snv.apobec$APOBEC=="Yes")])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    50.0    87.0   142.0   197.4   245.0   694.0

``` r
summary(hpv.pos.snv.apobec$SNV_count[which(hpv.pos.snv.apobec$APOBEC=="No")])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   53.00   70.75   74.50   92.00  123.50  141.00

``` r
summary(hpv.neg.snv.apobec$SNV_count[which(hpv.neg.snv.apobec$APOBEC=="Yes")])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    46.0   100.0   135.0   170.2   195.5   950.0

``` r
summary(hpv.neg.snv.apobec$SNV_count[which(hpv.neg.snv.apobec$APOBEC=="No")])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    46.0    91.0   126.0   249.5   233.0  3490.0

``` r
# 
# #HPV vs. APOBEC signature plot for the manuscript
# 
# 
# HPV.vs.APOBEC_ms <- ggplot(data = subset(SNV.APOBEC.df, !is.na(SNV.APOBEC.df$APOBEC))) + geom_boxplot(aes(y=SNV_count,x=APOBEC),color="dark red") + geom_jitter(aes(y=SNV_count,x=APOBEC),width= 0.2,alpha=0.5) + facet_grid(.~HPV_status) + theme_bw() + scale_y_log10(labels=fancy_scientific) + labs(x="APOBEC signature detected",y="SNV count")  + theme(strip.text.x = element_text(size=18)) +theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
# 
# HPV.vs.APOBEC_ms
# 
# ggsave(filename = "Figures/SNV_count_vs_HPV_and_APOBEC.png",plot = HPV.vs.APOBEC_ms)
```

Logistic regression of APOBEC and HPV
=====================================

``` r
load("output_data/HNSC_selection_with_APOBEC.RData")

SNV.APOBEC.df <- read.table("output_data/SNV_APOBEC_HPV_status.txt",header = T)

SNV.APOBEC.df.knownAPOBEC <- SNV.APOBEC.df[which(!is.na(SNV.APOBEC.df$APOBEC)),] 
# length(which(SNV.APOBEC.df.knownAPOBEC$HPV_status=="HPV+" & SNV.APOBEC.df.knownAPOBEC$APOBEC_weight==0))
SNV.APOBEC.df.knownAPOBEC$HPV_status_value <- ifelse(SNV.APOBEC.df.knownAPOBEC$HPV_status=="HPV+",1,0)

common.text.size.large <- 10

HPV_status_vs_APOBEC_weight <- ggplot() + 
  geom_point(data = SNV.APOBEC.df.knownAPOBEC,aes(y=HPV_status_value,x=APOBEC_weight),alpha=0.5,size=3,shape="|",color="blue") + 
  geom_abline(slope=0,intercept = 1) +
  stat_smooth(data = SNV.APOBEC.df.knownAPOBEC,aes(y=HPV_status_value,x=APOBEC_weight),method="glm", method.args=list(family="binomial"), se=T, fullrange=TRUE) + 
  geom_abline(slope=0,intercept = 1) + geom_segment(aes(x=1,xend=1,y=0,yend=1)) + 
  labs(y="Probability HPV positive",x="APOBEC signature weight") + theme_classic() +
  theme(text = element_text(size=common.text.size.large),plot.margin = margin(r=10,t=10,l=1)) + 
  expand_limits(x = c(0,1), y = 0) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))

HPV_vs_APOBEC_gridoff <- ggplot_gtable(ggplot_build(HPV_status_vs_APOBEC_weight))
HPV_vs_APOBEC_gridoff$layout$clip[HPV_vs_APOBEC_gridoff$layout$name == "panel"] <- "off"
library(grid)
library(gridExtra)

g <- arrangeGrob(HPV_vs_APOBEC_gridoff)
ggsave(plot = g, filename ="Figures/HPV_status_APOBEC_weight.png",width = 3.25,height = 2,dpi = 600)

# library(ggplot2)
HPV_status_vs_APOBEC_weight
```

![](APOBEC_HNSCC_manuscript_analysis_files/figure-markdown_github/logistic%20regression%20HPV%20APOBEC-1.png)

``` r
ggsave(plot = HPV_status_vs_APOBEC_weight, filename ="Figures/Fig2_HPVvsAPOBEC.png",width = 3.25,height = 2,dpi = 600)

log_regression <- glm(HPV_status_value~APOBEC_weight,data = SNV.APOBEC.df.knownAPOBEC,family = binomial)

summary(log_regression)
```

    ## 
    ## Call:
    ## glm(formula = HPV_status_value ~ APOBEC_weight, family = binomial, 
    ##     data = SNV.APOBEC.df.knownAPOBEC)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.1109  -0.4700  -0.3489  -0.2994   2.5009  
    ## 
    ## Coefficients:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)    -3.0824     0.2713 -11.361  < 2e-16 ***
    ## APOBEC_weight   3.3794     0.6494   5.204 1.95e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 302.58  on 455  degrees of freedom
    ## Residual deviance: 275.44  on 454  degrees of freedom
    ## AIC: 279.44
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
# nrow(SNV.APOBEC.df.knownAPOBEC)
# which(is.na(SNV.APOBEC.df.knownAPOBEC$APOBEC_weight))
# which(is.na(SNV.APOBEC.df.knownAPOBEC$APOBEC))

# length(unique(SNV.APOBEC.df.knownAPOBEC$Tumor))==nrow(SNV.APOBEC.df.knownAPOBEC)
message("HPV+ tumors with APOBEC signal")
```

    ## HPV+ tumors with APOBEC signal

``` r
length(which(SNV.APOBEC.df.knownAPOBEC$APOBEC=="Yes" & SNV.APOBEC.df.knownAPOBEC$HPV_status=="HPV+"))/length(which(SNV.APOBEC.df.knownAPOBEC$HPV_status=="HPV+"))
```

    ## [1] 0.8297872

``` r
message("HPV− tumors with APOBEC signal")
```

    ## HPV− tumors with APOBEC signal

``` r
length(which(SNV.APOBEC.df.knownAPOBEC$APOBEC=="Yes" & SNV.APOBEC.df.knownAPOBEC$HPV_status=="HPV−"))/length(which(SNV.APOBEC.df.knownAPOBEC$HPV_status=="HPV−"))
```

    ## [1] 0.6748166

``` r
# Among samples with an APOBEC signal, total mutation load was... 

# SNV.APOBEC.yes <- subset(SNV.APOBEC.df,APOBEC=="Yes")

wilcox.test(data = subset(SNV.APOBEC.df,APOBEC=="Yes"),SNV_count~HPV_status,conf.int=T)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  SNV_count by HPV_status
    ## W = 5414.5, p-value = 0.9521
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -23.99999  29.99998
    ## sample estimates:
    ## difference in location 
    ##              0.9999928

``` r
wilcox.test(data = subset(SNV.APOBEC.df,APOBEC=="No"),SNV_count~HPV_status,conf.int=T)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  SNV_count by HPV_status
    ## W = 278.5, p-value = 0.02414
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -108.000020   -3.000002
    ## sample estimates:
    ## difference in location 
    ##              -41.06115

``` r
# SNV.APOBEC.df$APOBEC <- factor(SNV.APOBEC.df$APOBEC,levels = c("No APOBEC Signature","APOBEC Signature","NA"))
SNV.APOBEC.df.comparingAPOBEC <- SNV.APOBEC.df

SNV.APOBEC.df.comparingAPOBEC$APOBEC <- as.character(SNV.APOBEC.df.comparingAPOBEC$APOBEC)
SNV.APOBEC.df.comparingAPOBEC$APOBEC[which(SNV.APOBEC.df.comparingAPOBEC$APOBEC=="Yes")] <- "APOBEC present"
SNV.APOBEC.df.comparingAPOBEC$APOBEC[which(SNV.APOBEC.df.comparingAPOBEC$APOBEC=="No")] <- "APOBEC absent"

SNV_count_comparing_APOBEC_plot <- ggplot(data = subset(SNV.APOBEC.df.comparingAPOBEC, !is.na(SNV.APOBEC.df.comparingAPOBEC$APOBEC))) + geom_boxplot(aes(y=SNV_count,x=HPV_status),color="dark red") + geom_jitter(aes(y=SNV_count,x=HPV_status),width= 0.2,alpha=0.5) + facet_grid(.~APOBEC) + theme_bw() + scale_y_log10(labels=fancy_scientific) + labs(x="HPV status",y="SNV count")  + theme(strip.text.x = element_text(size=15)) +theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))

# SNV_count_comparing_APOBEC_plot

ggsave(filename = "Figures/SNV_count_comparing_APOBEC.png",plot = SNV_count_comparing_APOBEC_plot)
```

    ## Saving 7 x 5 in image

``` r
library(cowplot)
fancy_scientific_justexp <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  l <- gsub("0e\\+00","0",l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  new.l <- strsplit(x = l,split = "%")[[1]][3]
  # return this as an expression
  parse(text=new.l)
}
common.text.size.large2 <- 12

SNV.APOBEC.df$HPV_status <- factor(SNV.APOBEC.df$HPV_status, labels = c(expression(HPV^{"+"}),expression(HPV^{"−"})))

HPV.vs.APOBEC_ms_forcombine <- ggplot(data = subset(SNV.APOBEC.df, !is.na(SNV.APOBEC.df$APOBEC))) +
  geom_boxplot(aes(y=SNV_count,x=APOBEC),color="black",outlier.shape = NA) + 
  geom_jitter(aes(y=SNV_count,x=APOBEC),width= 0.2,alpha=0.4,color="blue") + 
  facet_grid(.~HPV_status, labeller = "label_parsed") + 
  theme_classic() + 
  expand_limits( y = c(1,(max(SNV.APOBEC.df$SNV_count,na.rm = T)+1e3))) + 
  scale_y_log10(labels=c(expression(10^0),expression(10^01),expression(10^02),expression(10^03)),breaks=c(1,10,100,1000),expand = c(0,0)) + 
  labs(x="APOBEC signature present",y="SNV count")  + 
  theme(strip.text.x = element_text(size=common.text.size.large2)) +
  theme(axis.text=element_text(size=common.text.size.large2), axis.title=element_text(size=common.text.size.large2,face="bold"))
# plot_grid(

# SNV.APOBEC.df.comparingAPOBEC$HPV_status <- factor(SNV.APOBEC.df.comparingAPOBEC$HPV_status, levels = c(expression(HPV^"−"),expression(HPV^"+")))
SNV.APOBEC.df.comparingAPOBEC$APOBEC <- factor(SNV.APOBEC.df.comparingAPOBEC$APOBEC, levels = c("APOBEC present", "APOBEC absent",NA))

SNV_count_comparing_APOBEC_plot_forcombine <- ggplot(data = subset(SNV.APOBEC.df.comparingAPOBEC, !is.na(SNV.APOBEC.df.comparingAPOBEC$APOBEC))) + 
  geom_boxplot(aes(y=SNV_count,x=HPV_status),color="black",outlier.shape = NA) + 
  geom_jitter(aes(y=SNV_count,x=HPV_status),width= 0.2,alpha=0.4,color="blue") + 
  facet_grid(.~APOBEC) + 
  theme_classic() + 
  scale_x_discrete(labels = c(expression(HPV^{"+"}),expression(HPV^{"−"}))) + 
  expand_limits( y = c(1,(max(SNV.APOBEC.df.knownAPOBEC$SNV_count,na.rm = T)+1e3))) + 
  scale_y_log10(labels=c(expression(10^0),expression(10^01),expression(10^02),expression(10^03)),breaks=c(1,10,100,1000),expand = c(0,0)) + 
  labs(x="HPV status",y="SNV count")  + 
  theme(strip.text.x = element_text(size=common.text.size.large2)) + 
  theme(axis.text=element_text(size=common.text.size.large2), axis.title=element_text(size=common.text.size.large2,face="bold"),axis.title.y = element_blank())

combined_HPV_APOBEC_plot <- plot_grid(HPV.vs.APOBEC_ms_forcombine,SNV_count_comparing_APOBEC_plot_forcombine,align='h',ncol = 2,nrow = 1,labels = c("A","B"),label_size = common.text.size.large2)

combined_HPV_APOBEC_plot
```

![](APOBEC_HNSCC_manuscript_analysis_files/figure-markdown_github/combining%20SNV%20HPV%20APOBEC%20plots-1.png)

``` r
ggsave(plot = combined_HPV_APOBEC_plot,filename = "Figures/Fig3_combined_APOBEC_vs_HPV.png",height = 3,width = 3.25*2)


# head(SNV.APOBEC.df)
# SNV.APOBEC.df.HPVpos <- SNV.APOBEC.df[which(SNV.APOBEC.df$HPV_status=='HPV^{\n "+"\n}'),]
library(dplyr)
unique(SNV.APOBEC.df$HPV_status)
```

    ## [1] HPV^{\n    "+"\n} HPV^{\n    "−"\n}
    ## Levels: HPV^{\n    "+"\n} HPV^{\n    "−"\n}

``` r
SNV.APOBEC.df %>% group_by(HPV_status,APOBEC) %>% summarize(median(SNV_count)) 
```

    ## # A tibble: 6 x 3
    ## # Groups:   HPV_status [?]
    ##   HPV_status            APOBEC `median(SNV_count)`
    ##   <fct>                 <fct>                <dbl>
    ## 1 "HPV^{\n    \"+\"\n}" No                    74.5
    ## 2 "HPV^{\n    \"+\"\n}" Yes                  142  
    ## 3 "HPV^{\n    \"+\"\n}" <NA>                  41.0
    ## 4 "HPV^{\n    \"−\"\n}" No                   126  
    ## 5 "HPV^{\n    \"−\"\n}" Yes                  135  
    ## 6 "HPV^{\n    \"−\"\n}" <NA>                  33.0

``` r
SNV.APOBEC.df %>% group_by(HPV_status,APOBEC) %>% summarize(mean(SNV_count)) 
```

    ## # A tibble: 6 x 3
    ## # Groups:   HPV_status [?]
    ##   HPV_status            APOBEC `mean(SNV_count)`
    ##   <fct>                 <fct>              <dbl>
    ## 1 "HPV^{\n    \"+\"\n}" No                  92.0
    ## 2 "HPV^{\n    \"+\"\n}" Yes                197  
    ## 3 "HPV^{\n    \"+\"\n}" <NA>                38.6
    ## 4 "HPV^{\n    \"−\"\n}" No                 250  
    ## 5 "HPV^{\n    \"−\"\n}" Yes                170  
    ## 6 "HPV^{\n    \"−\"\n}" <NA>                31.3

``` r
SNV.APOBEC.df$`Proportion TCW to TKW` <- NA

for(i in 1:nrow(SNV.APOBEC.df)){
  SNV.APOBEC.df$`Proportion TCW to TKW`[i] <- length(which(HNSC.selection.HPVknown$TCW_TKW[which(HNSC.selection.HPVknown$Unique_patient_identifier==SNV.APOBEC.df$Tumor[i])]==1))/(length(which(HNSC.selection.HPVknown$TCW_TKW[which(HNSC.selection.HPVknown$Unique_patient_identifier==SNV.APOBEC.df$Tumor[i])]==1))+length(which(HNSC.selection.HPVknown$TCW_TKW[which(HNSC.selection.HPVknown$Unique_patient_identifier==SNV.APOBEC.df$Tumor[i])]==0)))
}

HPV.vs.APOBEC_proportion <- ggplot(data = SNV.APOBEC.df) + geom_boxplot(aes(y=`Proportion TCW to TKW`,x=APOBEC),color="dark red") + geom_jitter(aes(y=`Proportion TCW to TKW`,x=APOBEC),width= 0.2,alpha=0.5) + facet_grid(.~HPV_status) + theme_bw() + labs(x="APOBEC signature detected")

ggsave(filename = "Figures/HPV_status_and_APOBEC_vs_proportionTCW.png",plot = HPV.vs.APOBEC_proportion)
```

    ## Saving 7 x 5 in image

``` r
# HPV.vs.APOBEC_proportion_ms <- ggplot(data = subset(SNV.APOBEC.df,!is.na(SNV.APOBEC.df$APOBEC))) + geom_boxplot(aes(y=`Proportion TCW to TKW`,x=APOBEC),color="dark red") + geom_jitter(aes(y=`Proportion TCW to TKW`,x=APOBEC),width= 0.2,alpha=0.5) + facet_grid(.~HPV_status) + theme_bw() + labs(x="APOBEC signature detected")

wilcox.test(data=SNV.APOBEC.df[which(SNV.APOBEC.df$APOBEC=="Yes"),],`Proportion TCW to TKW`~HPV_status,conf.int=T)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Proportion TCW to TKW by HPV_status
    ## W = 7338.5, p-value = 0.0002389
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  0.05079741 0.18093822
    ## sample estimates:
    ## difference in location 
    ##              0.1136811

``` r
wilcox.test(data=SNV.APOBEC.df[which(SNV.APOBEC.df$APOBEC=="No"),],`Proportion TCW to TKW`~HPV_status,conf.int=T)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Proportion TCW to TKW by HPV_status
    ## W = 482.5, p-value = 0.6623
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.02411994  0.01835812
    ## sample estimates:
    ## difference in location 
    ##           -0.003107094

``` r
# ggsave(filename = "Figures/HPV_status_and_APOBEC_vs_proportionTCW_forMS.png",plot = HPV.vs.APOBEC_proportion_ms)


SNV.APOBEC.df$total_TCW <- SNV.APOBEC.df$SNV_count * SNV.APOBEC.df$`Proportion TCW to TKW`
SNV.APOBEC.df$total_TCW_x_APOBEC_weight <- SNV.APOBEC.df$total_TCW * SNV.APOBEC.df$APOBEC_weight

SNV.APOBEC.df %>% 
  group_by(HPV_status, APOBEC) %>%
  summarize(sum(total_TCW))
```

    ## # A tibble: 6 x 3
    ## # Groups:   HPV_status [?]
    ##   HPV_status            APOBEC `sum(total_TCW)`
    ##   <fct>                 <fct>             <dbl>
    ## 1 "HPV^{\n    \"+\"\n}" No                 53.0
    ## 2 "HPV^{\n    \"+\"\n}" Yes              3509  
    ## 3 "HPV^{\n    \"+\"\n}" <NA>              120  
    ## 4 "HPV^{\n    \"−\"\n}" No               2514  
    ## 5 "HPV^{\n    \"−\"\n}" Yes             13005  
    ## 6 "HPV^{\n    \"−\"\n}" <NA>              153

``` r
SNV.APOBEC.df %>% 
  group_by(HPV_status, APOBEC) %>%
  summarize(sum(total_TCW_x_APOBEC_weight,na.rm = T))
```

    ## # A tibble: 6 x 3
    ## # Groups:   HPV_status [?]
    ##   HPV_status            APOBEC `sum(total_TCW_x_APOBEC_weight, na.rm = T)`
    ##   <fct>                 <fct>                                        <dbl>
    ## 1 "HPV^{\n    \"+\"\n}" No                                               0
    ## 2 "HPV^{\n    \"+\"\n}" Yes                                           2340
    ## 3 "HPV^{\n    \"+\"\n}" <NA>                                             0
    ## 4 "HPV^{\n    \"−\"\n}" No                                               0
    ## 5 "HPV^{\n    \"−\"\n}" Yes                                           5704
    ## 6 "HPV^{\n    \"−\"\n}" <NA>                                             0

``` r
# nrow(SNV.APOBEC.df)
# unique(SNV.APOBEC.df$APOBEC)
# length(which(SNV.APOBEC.df$APOBEC=="Yes" & SNV.APOBEC.df$HPV_status=="HPV+"))

# probably.wrong.apo.pos.hpv.pos <- as.character(SNV.APOBEC.df$Tumor[which(SNV.APOBEC.df$APOBEC=="Yes" & SNV.APOBEC.df$HPV_status=="HPV+")])



SNV.APOBEC.tib <- as_tibble(SNV.APOBEC.df) %>%
  group_by(APOBEC,HPV_status) %>%
  summarize(n())

# SNV.HPVplus.total <- SNV.APOBEC.tib$`n()`[which(SNV.APOBEC.tib$HPV_status=="HPV+")]

SNV.APOBEC.tib
```

    ## # A tibble: 6 x 3
    ## # Groups:   APOBEC [?]
    ##   APOBEC HPV_status            `n()`
    ##   <fct>  <fct>                 <int>
    ## 1 No     "HPV^{\n    \"+\"\n}"     8
    ## 2 No     "HPV^{\n    \"−\"\n}"   133
    ## 3 Yes    "HPV^{\n    \"+\"\n}"    39
    ## 4 Yes    "HPV^{\n    \"−\"\n}"   276
    ## 5 <NA>   "HPV^{\n    \"+\"\n}"    22
    ## 6 <NA>   "HPV^{\n    \"−\"\n}"    41

``` r
# 7/(37+7)
# 37/(37+7)
# 
# 256/(256+124)
# 124/(256+124)

# fisher.test(matrix(data = c(2,2,250,250),nrow=2))
fisher.test(matrix(data = SNV.APOBEC.tib$`n()`[1:4],nrow=2))
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  matrix(data = SNV.APOBEC.tib$`n()`[1:4], nrow = 2)
    ## p-value = 0.03043
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  0.1673275 0.9591682
    ## sample estimates:
    ## odds ratio 
    ##  0.4263654

Mutations that have the highest prevalence and selection intensity in HNSCC
---------------------------------------------------------------------------

``` r
load("output_data/HNSC_selection_with_APOBEC_recur.RData")
just.ALL <- HNSC.selection.output.recur
ALL.names <- unique(just.ALL$Name)

ALL.prev.df <- as.data.frame(matrix(data = NA, nrow=length(ALL.names),ncol=5))
colnames(ALL.prev.df) <- c("Mutation","Frequency", "Trinuc context", "Alternative Nucleotide","Selection intensity")

ALL.prev.df$Mutation <- ALL.names

for(i in 1:nrow(ALL.prev.df)){
  ALL.prev.df$Frequency[i] <- length(which(just.ALL$Name == ALL.prev.df$Mutation[i]))
  ALL.prev.df$`Trinuc context`[i] <- just.ALL$Nucleotide_trinuc_context[which(just.ALL$Name == ALL.prev.df$Mutation[i])[1]]
  ALL.prev.df$`Alternative Nucleotide`[i] <- just.ALL$Alternative_Nucleotide[which(just.ALL$Name == ALL.prev.df$Mutation[i])[1]]
  ALL.prev.df$`Selection intensity`[i] <- round(just.ALL$Gamma_epistasis[which(just.ALL$Name == ALL.prev.df$Mutation[i])[1]],2)
}

ALL.prev.df <- ALL.prev.df[order(ALL.prev.df$Frequency,decreasing = T),]
```

    ## 3 substitutions occur in only tumors
    ## that did not have enough mutations to measure trinucleotide context
    ## and will be removed from the analysis

    ## 3 substitutions occur in only tumors
    ## that did not have enough mutations to measure trinucleotide context
    ## and will be removed from the analysis

``` r
# Only one amino acid substitution caused by two different nucleotide substitutions, but neither are TCW --> TKW
# table(ALL.prev.original$short_name)[order(table(ALL.prev.original$short_name),decreasing = T)]
# ALL.prev.original[which(ALL.prev.original$short_name=="CDKN2A W110*"),]

head(HNSC.selection.output)
```

    ## # A tibble: 6 x 31
    ##   Gene   Nucleotide_Gene_P… Nucleotide_chromosome_po… Nucleotide_change_t…
    ##   <chr>               <dbl>                     <int>                <dbl>
    ## 1 RPA2                   NA                  28240999                 1.00
    ## 2 RPA2                  292                  28233480                 1.00
    ## 3 SLC6A9                480                  44475695                 1.00
    ## 4 SLC6A9                130                  44477352                 1.00
    ## 5 SLC6A9                621                  44474213                 1.00
    ## 6 SLC6A9               1181                  44467985                 1.00
    ## # ... with 27 more variables: Nucleotide_mutation_rate <dbl>,
    ## #   Nucleotide_trinuc_context <chr>, Chromosome <chr>,
    ## #   Reference_Nucleotide <chr>, Reference_Count <chr>,
    ## #   Alternative_Nucleotide <chr>, Alternative_Count <chr>,
    ## #   Tumor_origin <chr>, Unique_patient_identifier <chr>,
    ## #   Amino_acid_position <int>, Amino_acid_codon <chr>,
    ## #   Amino_acid_reference <chr>, Amino_acid_alternative <chr>,
    ## #   Amino_acid_mutation_rate <dbl>, Amino_acid_change_tally <int>,
    ## #   Codon_position <dbl>, synonymous.mu <dbl>, trinucs <chr>,
    ## #   Gamma_epistasis <dbl>, HPV_call <chr>, Sig_2 <dbl>, Sig_13 <dbl>,
    ## #   APOBEC_weight <dbl>, TCW_TKW <dbl>, TCN_TKN <dbl>, Name <chr>,
    ## #   Name_short <chr>

``` r
unique(ALL.prev.original$`Trinuc context`)
```

    ##  [1] "TGA" "CCG" "CGC" "CAT" "TCG" "CGT" "TAT" "CGG" "GGT" "CGA" "GCG"
    ## [12] "TCC" "TAC" "AGA" "GGG" "TGG" "AAG" "TCA" "TGT" "GGC" "GGA" "ACG"
    ## [23] "GTA" "CCC" "TGC" "GCA" "ACC" "CCA" "TTG" "GAG" "AGG" "ATG" "GCT"
    ## [34] "TCT" "GCC" "TAG" "ACA" "CCT" "AAT" "CAG" "AAC" "AGC" "ACT" "AGT"
    ## [45] "GAA"

``` r
ALL.prev.df$`APOBEC type TCW` <- "False"

for( i in 1:nrow(ALL.prev.df)){
  if(length(which(HNSC.selection.output$Name_short==ALL.prev.df$`Short name`[i]))==0){print(paste("NO MATCH!",i))}
  if(HNSC.selection.output$TCW_TKW[which(HNSC.selection.output$Name_short==ALL.prev.df$`Short name`[i])[1]]==1){
    ALL.prev.df$`APOBEC type TCW`[i] <- "True"
  }
}



wilcox.test(data=ALL.prev.df, `Selection intensity`~`APOBEC type TCW`)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Selection intensity by APOBEC type TCW
    ## W = 15126, p-value = 2.069e-07
    ## alternative hypothesis: true location shift is not equal to 0

``` r
common.text.size.large3 <- 20

distribution.of.APOBEC.TCW <- ggplot(data = ALL.prev.df) + 
  geom_boxplot(aes(y=`Selection intensity`,x=`APOBEC type TCW`),width=0.5,outlier.shape = NA) + 
  geom_jitter(aes(y=`Selection intensity`,x=`APOBEC type TCW`),alpha=0.4,width  = 0.1,size=1.5,col="red") + 
  theme_classic() + 
  labs(y="Selection intensity") + 
  expand_limits( y = c(1,(max(ALL.prev.df$`Selection intensity`,na.rm = T)+1e5))) + 
  scale_y_log10(labels=c(expression(10^0),expression(10^01),expression(10^02),expression(10^03),expression(10^04),expression(10^05)),breaks=c(1,10,100,1000,1e4,1e5),expand = c(0,0)) + 
  # scale_y_log10(labels=c(expression(10^03),expression(10^05)),breaks=c(1000,100000)) + 
  theme(axis.text=element_text(size=common.text.size.large3), axis.title=element_text(size=common.text.size.large3,face="bold")) + 
  labs(x=expression(TCW %->% TKW))

distribution.of.APOBEC.TCW
```

![](APOBEC_HNSCC_manuscript_analysis_files/figure-markdown_github/Selection%20intensity%20vs%20TCW%20to%20TKW-1.png)

``` r
ggsave(filename = "Figures/Fig5_selection_vs_TCW.png",plot = distribution.of.APOBEC.TCW,width = 3.25,height = 3.25,dpi = 600)

message("Number of recurrent mutations:");nrow(ALL.prev.df)
```

    ## Number of recurrent mutations:

    ## [1] 368

``` r
length(which(ALL.prev.df$`APOBEC type TCW`=="True"))
```

    ## [1] 74

``` r
# ALL.prev.df[which(ALL.prev.df$`APOBEC type TCW`=="True"),]
```

``` r
# head(signature.breakdown.ALL)
# head(ALL.prev.df)

ALL.prev.df.split <- as.data.frame(matrix(data = NA,nrow = 2*nrow(ALL.prev.df),ncol=5))
colnames(ALL.prev.df.split) <- c("Name","Selection from mutational context","APOBEC context", "TCW to TKW", "TCN to TKN")
ALL.prev.df.split$Name <- c(ALL.prev.df$`Short name`,ALL.prev.df$`Short name`)
ALL.prev.df.split$`APOBEC context` <-  c(rep("APOBEC processes",nrow(ALL.prev.df.split)/2),rep("non-APOBEC processes",nrow(ALL.prev.df.split)/2))


for(i in 1:length(unique(ALL.prev.df.split$Name))){
  these.pos <- which(ALL.prev.df.split$Name ==  unique(ALL.prev.df.split$Name)[i])
  # ALL.prev.df.split$`TCW to TKW`[these.pos] <- ALL.prev.df$`APOBEC type TCW`[i]
  # ALL.prev.df.split$`TCN to TKN`[these.pos] <- ALL.prev.df$`APOBEC type TCN`[i]
  ALL.prev.df.split$`Selection from mutational context`[these.pos] <- c(ALL.prev.df$`Selection from APOBEC`[i],ALL.prev.df$`Selection from NonAPOBEC`[i])
}

head(ALL.prev.df.split)
```

    ##           Name Selection from mutational context   APOBEC context
    ## 1 PIK3CA E545K                            720.63 APOBEC processes
    ## 2   TP53 R283P                            640.30 APOBEC processes
    ## 3   TP53 R196P                            385.67 APOBEC processes
    ## 4 PIK3CA E542K                            340.02 APOBEC processes
    ## 5  EPHA2 S419*                            217.82 APOBEC processes
    ## 6 EP300 D1399N                            172.51 APOBEC processes
    ##   TCW to TKW TCN to TKN
    ## 1         NA         NA
    ## 2         NA         NA
    ## 3         NA         NA
    ## 4         NA         NA
    ## 5         NA         NA
    ## 6         NA         NA

``` r
# ALL.prev.df.split$short_name <- NA
# for(i in 1:nrow(ALL.prev.df.split)){
#   ALL.prev.df.split$short_name[i] <- paste(strsplit(x = ALL.prev.df.split$short_name[i], split=" ")[[1]][1:2],collapse = " ")
# }

# this quick fix messes up the non-coding mutations, fix the one we want labeled 

ALL.prev.df.split$Name[which(ALL.prev.df.split$Name=="CCDC50 A 191098615 C NCSNV")] <- "CCDC50 S.S."

# ggplot(data = ALL.prev.df.split) +geom_violin(aes(x=`APOBEC context`, y = `Selection from mutational context`))  + geom_jitter(aes(x=`APOBEC context`, y = `Selection from mutational context`))

selection.from.contexts <- ggplot(data = ALL.prev.df.split) +geom_violin(aes(x=`APOBEC context`, y = log10(`Selection from mutational context`)))  + geom_jitter(aes(x=`APOBEC context`, y = log10(`Selection from mutational context`))) + labs(x="Mutational signature attributed to mutation")
source("R/fancy_scientific_code.R")

# selection.from.contexts <- ggplot(data = ALL.prev.df.split) +geom_violin(aes(x=`APOBEC context`, y = `Selection from mutational context`))  + geom_jitter(aes(x=`APOBEC context`, y = `Selection from mutational context`)) + scale_y_log10(labels=fancy_scientific) + labs(x="Mutational signatures attributed to mutation", y="Proportion of selection intensity attributed to signature and \nweighted by mutation prevalence among tumors") + theme_bw()

ALL.prev.df.split$`APOBEC context` <- factor(ALL.prev.df.split$`APOBEC context`, levels = c("APOBEC processes","non-APOBEC processes"),labels = c("APOBEC","non-APOBEC"))

library(ggrepel)

selection.from.contexts.labels <- ggplot(data = ALL.prev.df.split) + 
  geom_point(aes(x=`APOBEC context`, y = `Selection from mutational context`),alpha=0.2,col="red") + 
  geom_text_repel(data = subset(ALL.prev.df.split,`APOBEC context`=="APOBEC" & `Selection from mutational context` > 100),aes(x=`APOBEC context`, y = `Selection from mutational context`,label=Name),color="black",segment.alpha = 0.4,size=common.text.size*(5/14)) + #scale_y_log10(labels=fancy_scientific) + 
  geom_text_repel(data = subset(ALL.prev.df.split,`APOBEC context`=="non-APOBEC" & `Selection from mutational context` > 650),aes(x=`APOBEC context`, y = `Selection from mutational context`,label=Name),box.padding = .31,segment.alpha=0.4,color="black",size=common.text.size*(5/14)) + 
  labs(y="Net realized selection intensity",x="Mutation process") + theme_classic() + 
  theme(axis.text=element_text(size=common.text.size.large), axis.title=element_text(size=common.text.size.large2,face="bold")) +
  expand_limits(x = 0, y = c(0,max(ALL.prev.df.split$`Selection from mutational context`)+100)) + 
  # scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + scale_x_discrete(expand=c(0,1)) # expand axis so symmetrical 


selection.from.contexts.labels
```

![](APOBEC_HNSCC_manuscript_analysis_files/figure-markdown_github/Mutation%20type%20vs%20gamma%20from%20APOBEC%20or%20not-1.png)

``` r
# ggsave(filename = "Figures/selection_from_process.png",plot = selection.from.contexts)
# View(ALL.prev.df.split)
ggsave(filename = "Figures/Fig6_selection_from_process_labels.png",plot = selection.from.contexts.labels,width = 3.25,height = 3.25,dpi = 600)

wilcox.test(`Selection from mutational context`~`APOBEC context`,data = ALL.prev.df.split)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Selection from mutational context by APOBEC context
    ## W = 32770, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0

``` r
# HNSC.MAF[which(HNSC.MAF$Start_Position==191098615),]
```

``` r
save.image(file = "output_data/completed_analysis_APOBEC_HNSCC.RData")
sessionInfo(package = NULL)
```

    ## R version 3.4.0 (2017-04-21)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: OS X El Capitan 10.11.6
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ##  [1] grid      parallel  stats4    stats     graphics  grDevices utils    
    ##  [8] datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] dendextend_1.6.0      ggdendro_0.1-20       viridis_0.5.0        
    ##  [4] viridisLite_0.3.0     gridExtra_2.3         cowplot_0.9.2        
    ##  [7] bindrcpp_0.2          forcats_0.2.0         stringr_1.2.0        
    ## [10] dplyr_0.7.4           purrr_0.2.4           readr_1.1.1          
    ## [13] tidyr_0.8.0           tibble_1.4.2          tidyverse_1.2.1      
    ## [16] deconstructSigs_1.8.0 reshape2_1.4.3        ggrepel_0.7.0        
    ## [19] ggplot2_2.2.1         rtracklayer_1.36.6    GenomicRanges_1.28.6 
    ## [22] GenomeInfoDb_1.12.3   IRanges_2.10.5        S4Vectors_0.14.7     
    ## [25] BiocGenerics_0.22.1  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-131                      bitops_1.0-6                     
    ##  [3] matrixStats_0.53.0                lubridate_1.7.2                  
    ##  [5] httr_1.3.1                        prabclus_2.2-6                   
    ##  [7] rprojroot_1.3-2                   tools_3.4.0                      
    ##  [9] backports_1.1.2                   utf8_1.1.3                       
    ## [11] R6_2.2.2                          lazyeval_0.2.1                   
    ## [13] colorspace_1.3-2                  trimcluster_0.1-2                
    ## [15] nnet_7.3-12                       mnormt_1.5-5                     
    ## [17] compiler_3.4.0                    cli_1.0.0                        
    ## [19] rvest_0.3.2                       Biobase_2.36.2                   
    ## [21] xml2_1.2.0                        DelayedArray_0.2.7               
    ## [23] labeling_0.3                      diptest_0.75-7                   
    ## [25] scales_0.5.0                      DEoptimR_1.0-8                   
    ## [27] robustbase_0.92-8                 mvtnorm_1.0-7                    
    ## [29] psych_1.7.8                       digest_0.6.15                    
    ## [31] Rsamtools_1.28.0                  foreign_0.8-69                   
    ## [33] rmarkdown_1.8                     BSgenome.Hsapiens.UCSC.hg19_1.4.0
    ## [35] XVector_0.16.0                    pkgconfig_2.0.1                  
    ## [37] htmltools_0.3.6                   BSgenome_1.44.2                  
    ## [39] rlang_0.1.6                       readxl_1.0.0                     
    ## [41] rstudioapi_0.7                    bindr_0.1                        
    ## [43] jsonlite_1.5                      mclust_5.4                       
    ## [45] BiocParallel_1.10.1               RCurl_1.95-4.10                  
    ## [47] magrittr_1.5                      modeltools_0.2-21                
    ## [49] GenomeInfoDbData_0.99.0           Matrix_1.2-12                    
    ## [51] Rcpp_0.12.15                      munsell_0.4.3                    
    ## [53] stringi_1.1.6                     whisker_0.3-2                    
    ## [55] yaml_2.1.16                       MASS_7.3-48                      
    ## [57] SummarizedExperiment_1.6.5        zlibbioc_1.22.0                  
    ## [59] flexmix_2.3-14                    plyr_1.8.4                       
    ## [61] crayon_1.3.4                      lattice_0.20-35                  
    ## [63] Biostrings_2.44.2                 haven_1.1.1                      
    ## [65] hms_0.4.1                         knitr_1.19                       
    ## [67] pillar_1.1.0                      fpc_2.1-11                       
    ## [69] XML_3.98-1.9                      glue_1.2.0                       
    ## [71] evaluate_0.10.1                   modelr_0.1.1                     
    ## [73] cellranger_1.1.0                  gtable_0.2.0                     
    ## [75] kernlab_0.9-25                    assertthat_0.2.0                 
    ## [77] broom_0.4.3                       class_7.3-14                     
    ## [79] GenomicAlignments_1.12.2          cluster_2.0.6
