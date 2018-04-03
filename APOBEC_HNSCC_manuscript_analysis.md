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

# assigning HPV status
HNSC.MAF$HPV_call[which(startsWith(x = HNSC.MAF$Unique_patient_identifier,prefix = "PY"))] <- "HPV−"
HNSC.MAF$HPV_call[which(HNSC.MAF$Unique_patient_identifier=="PY-16T")] <- "HPV+"


HNSC.MAF.hpvpos <- HNSC.MAF[which(HNSC.MAF$HPV_call=="HPV+"),]
MAF_for_analysis <- HNSC.MAF.hpvpos
length(unique(HNSC.MAF.hpvpos$Unique_patient_identifier))
```

    ## [1] 67

``` r
save(MAF_for_analysis, file="output_data/HNSC_HPVpos_MAF.RData")

HNSC.MAF.hpvneg <- HNSC.MAF[which(HNSC.MAF$HPV_call=="HPV−"),]
MAF_for_analysis <- HNSC.MAF.hpvneg
length(unique(HNSC.MAF.hpvneg$Unique_patient_identifier))
```

    ## [1] 427

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

Trinucleotide heatmaps
======================

The SNV selection intensity pipeline was run on the HPV data. The pipeline may be found here: <https://github.com/Townsend-Lab-Yale/SNV_selection_intensity> and the associated manuscript is here: <https://doi.org/10.1101/229724>

``` r
library(ggplot2)
load("output_data/selection_from_cluster/HNSC_HPVpos/trinuc_output/trinuc_data_HNSC_HPVpos.RData")
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
load("output_data/selection_from_cluster/HNSC_HPVneg/trinuc_output/trinuc_data_HNSC_HPVneg.RData")
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
HPV.pos.mut.rates <- read.csv(file = "output_data/selection_from_cluster/HNSC_HPVpos/mutsig_output/MAF_HNSC_HPVpos.txt.gene_rates.txt",header = T,sep = "\t",stringsAsFactors = F)
# head(HPV.pos.mut.rates)

HPV.neg.mut.rates <- read.csv(file = "output_data/selection_from_cluster/HNSC_HPVneg/mutsig_output/MAF_HNSC_HPVneg.txt.gene_rates.txt",header = T,sep = "\t",stringsAsFactors = F)
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
    ## -7.747e-06 -6.350e-07 -1.470e-07  4.020e-07  1.170e-04 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        1.128e-06  2.756e-08   40.94   <2e-16 ***
    ## negative_mut_rates 1.667e-01  1.109e-02   15.04   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2.418e-06 on 18860 degrees of freedom
    ## Multiple R-squared:  0.01185,    Adjusted R-squared:  0.01179 
    ## F-statistic: 226.1 on 1 and 18860 DF,  p-value: < 2.2e-16

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

    ## ── Attaching packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ tibble  1.4.2     ✔ purrr   0.2.4
    ## ✔ tidyr   0.8.0     ✔ dplyr   0.7.4
    ## ✔ readr   1.1.1     ✔ stringr 1.2.0
    ## ✔ tibble  1.4.2     ✔ forcats 0.2.0

    ## Warning: package 'tibble' was built under R version 3.4.3

    ## Warning: package 'tidyr' was built under R version 3.4.3

    ## Warning: package 'purrr' was built under R version 3.4.2

    ## Warning: package 'dplyr' was built under R version 3.4.2

    ## ── Conflicts ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
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
load("output_data/selection_from_cluster/HNSC/selection_output/HNSC_all_selection_output.RData")
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

    ## [1] 46

``` r
length(which(trinuc.w.HPV$HPV=="HPV−")) # ...how many tumors are HPV- 
```

    ## [1] 386

``` r
length(which(is.na(trinuc.w.HPV$HPV)))
```

    ## [1] 29

``` r
save(trinuc.w.HPV,file = "output_data/signature_weights_w_HPV.RData")


length(unique(HNSC.selection.output$Unique_patient_identifier[which(HNSC.selection.output$HPV_call=="HPV+" & HNSC.selection.output$APOBEC_weight > 0)])) # Out of all tumors, how many were HPV+ and had an APOBEC signature 
```

    ## [1] 38

``` r
# length(which((trinuc.w.HPV$`Signature.2` > 0 | trinuc.w.HPV$`Signature.13`>0) & trinuc.w.HPV$HPV=="HPV+"))



length(which((trinuc.w.HPV$`Signature.2` > 0 | trinuc.w.HPV$`Signature.13`>0) & trinuc.w.HPV$HPV=="HPV−"))/length(which(trinuc.w.HPV$HPV=="HPV−")) # proportion of HPV- tumors with enough substitutions to measure signatures that have APOBEC signature
```

    ## [1] 0.6735751

``` r
length(which((trinuc.w.HPV$`Signature.2` > 0 | trinuc.w.HPV$`Signature.13`>0) & trinuc.w.HPV$HPV=="HPV+"))/length(which(trinuc.w.HPV$HPV=="HPV+")) # proportion of HPV+ tumors with enough substitutions to measure signatures that have APOBEC signature
```

    ## [1] 0.826087

``` r
# mean weights
mean(trinuc.w.HPV$`Signature.2`[which(trinuc.w.HPV$HPV=="HPV+")])
```

    ## [1] 0.2081386

``` r
mean(trinuc.w.HPV$`Signature.13`[which(trinuc.w.HPV$HPV=="HPV+")])
```

    ## [1] 0.1579057

``` r
mean(trinuc.w.HPV$`Signature.2`[which(trinuc.w.HPV$HPV=="HPV−")])
```

    ## [1] 0.07588167

``` r
mean(trinuc.w.HPV$`Signature.13`[which(trinuc.w.HPV$HPV=="HPV−")])
```

    ## [1] 0.1089066

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
  geom_text_repel(data = subset(mutation_rates, gene %in% c("FBXW7","PIK3CA","MAPK1")),aes(x=positive_mut_rates,y=negative_mut_rates,label=gene),size=common.text.size*(5/14),segment.alpha = 0.3,col="black",box.padding = 0.6) +
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
load("output_data/selection_from_cluster/HNSC_HPVpos/selection_output/HNSC_HPVpos_selection_output.RData")
HPV.pos.selectionoutput <- selection.output
HPV.pos.selection_minrecur <- subset(HPV.pos.selectionoutput$all_mutations, freq>1)


load("output_data/selection_from_cluster/HNSC_HPVneg/selection_output/HNSC_HPVneg_selection_output.RData")
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

    ## Warning: Removed 258 rows containing missing values (geom_point).

    ## Warning: Removed 258 rows containing missing values (geom_text_repel).

![](APOBEC_HNSCC_manuscript_analysis_files/figure-markdown_github/gamma%20gamma%20plot-1.png)

``` r
ggsave(filename = "Figures/gamma_gamma_plot.png",plot = gamma_plot,height = 2,width = 2,dpi=600)
```

    ## Warning: Removed 258 rows containing missing values (geom_point).

    ## Warning: Removed 258 rows containing missing values (geom_text_repel).

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
load("output_data/selection_from_cluster/HNSC_HPVneg/selection_output/HNSC_HPVneg_selection_output.RData")

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

load("output_data/selection_from_cluster/HNSC_HPVpos/selection_output/HNSC_HPVpos_selection_output.RData")
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
        plot.margin = unit(c(1,-1,1,1.5), "mm")) +
  theme(panel.background = element_blank()) + scale_fill_manual("Legend", values = neg.colors) +
  theme(panel.grid.major =element_line(color="lightgrey"),panel.grid.minor =element_line(color="lightgrey")) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) +
  # geom_text(aes(label=round(mu*1e5,2)), vjust=-0.5, hjust=0.5, position=position_dodge(width=0.9),size=4,angle=90) +
  geom_text(aes(label=freq,y=-4e-8), position=position_dodge(width=0.9),size=common.text.size*(5/14),angle=0,color="black") +
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
                  box.padding =.36,
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

    ## Warning: Removed 258 rows containing missing values (geom_point).

    ## Warning: Removed 258 rows containing missing values (geom_text_repel).

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

    ## Warning: Removed 258 rows containing missing values (geom_point).

    ## Warning: Removed 258 rows containing missing values (geom_text_repel).

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
    ## W = 15940, p-value = 0.6689
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -17.00005  12.00001
    ## sample estimates:
    ## difference in location 
    ##              -3.000012

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
    ## W = 76.5, p-value = 0.02972
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -149.00000   -3.00002
    ## sample estimates:
    ## difference in location 
    ##              -49.81802

``` r
hpv.pos.snv.apobec <- subset(SNV.APOBEC.df,HPV_status=="HPV+") 
hpv.neg.snv.apobec <- subset(SNV.APOBEC.df,HPV_status=="HPV−")

mean(hpv.pos.snv.apobec$SNV_count[which(hpv.pos.snv.apobec$APOBEC=="Yes")])-mean(hpv.pos.snv.apobec$SNV_count[which(hpv.pos.snv.apobec$APOBEC=="No")]) #SNV count difference in mean
```

    ## [1] 105.2895

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

    ## [1] -67.69072

``` r
median(hpv.neg.snv.apobec$SNV_count[which(hpv.neg.snv.apobec$APOBEC=="Yes")])-median(hpv.neg.snv.apobec$SNV_count[which(hpv.neg.snv.apobec$APOBEC=="No")])
```

    ## [1] 7

``` r
summary(hpv.pos.snv.apobec$SNV_count[which(hpv.pos.snv.apobec$APOBEC=="Yes")])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    50.0    86.0   142.0   197.3   252.5   694.0

``` r
summary(hpv.pos.snv.apobec$SNV_count[which(hpv.pos.snv.apobec$APOBEC=="No")])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   53.00   70.75   74.50   92.00  123.50  141.00

``` r
summary(hpv.neg.snv.apobec$SNV_count[which(hpv.neg.snv.apobec$APOBEC=="Yes")])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    46.0   100.0   133.5   168.5   197.0   950.0

``` r
summary(hpv.neg.snv.apobec$SNV_count[which(hpv.neg.snv.apobec$APOBEC=="No")])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   46.00   91.25  126.50  236.20  232.25 3490.00

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
    ## -1.1257  -0.4833  -0.3552  -0.3074   2.4802  
    ## 
    ## Coefficients:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)    -3.0284     0.2729 -11.098  < 2e-16 ***
    ## APOBEC_weight   3.3583     0.6617   5.075 3.87e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 292.98  on 431  degrees of freedom
    ## Residual deviance: 267.16  on 430  degrees of freedom
    ## AIC: 271.16
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

    ## [1] 0.826087

``` r
message("HPV− tumors with APOBEC signal")
```

    ## HPV− tumors with APOBEC signal

``` r
length(which(SNV.APOBEC.df.knownAPOBEC$APOBEC=="Yes" & SNV.APOBEC.df.knownAPOBEC$HPV_status=="HPV−"))/length(which(SNV.APOBEC.df.knownAPOBEC$HPV_status=="HPV−"))
```

    ## [1] 0.6735751

``` r
# Among samples with an APOBEC signal, total mutation load was... 

# SNV.APOBEC.yes <- subset(SNV.APOBEC.df,APOBEC=="Yes")

wilcox.test(data = subset(SNV.APOBEC.df,APOBEC=="Yes"),SNV_count~HPV_status,conf.int=T)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  SNV_count by HPV_status
    ## W = 4927.5, p-value = 0.9807
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -25.00004  28.00002
    ## sample estimates:
    ## difference in location 
    ##             -0.9999884

``` r
wilcox.test(data = subset(SNV.APOBEC.df,APOBEC=="No"),SNV_count~HPV_status,conf.int=T)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  SNV_count by HPV_status
    ## W = 266.5, p-value = 0.02603
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -108.000089   -2.999956
    ## sample estimates:
    ## difference in location 
    ##                -42.524

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
    ## 3 "HPV^{\n    \"+\"\n}" <NA>                  42.0
    ## 4 "HPV^{\n    \"−\"\n}" No                   126  
    ## 5 "HPV^{\n    \"−\"\n}" Yes                  134  
    ## 6 "HPV^{\n    \"−\"\n}" <NA>                  34.0

``` r
SNV.APOBEC.df %>% group_by(HPV_status,APOBEC) %>% summarize(mean(SNV_count)) 
```

    ## # A tibble: 6 x 3
    ## # Groups:   HPV_status [?]
    ##   HPV_status            APOBEC `mean(SNV_count)`
    ##   <fct>                 <fct>              <dbl>
    ## 1 "HPV^{\n    \"+\"\n}" No                  92.0
    ## 2 "HPV^{\n    \"+\"\n}" Yes                197  
    ## 3 "HPV^{\n    \"+\"\n}" <NA>                38.9
    ## 4 "HPV^{\n    \"−\"\n}" No                 236  
    ## 5 "HPV^{\n    \"−\"\n}" Yes                169  
    ## 6 "HPV^{\n    \"−\"\n}" <NA>                31.7

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
    ## W = 6708, p-value = 0.0003675
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  0.04709324 0.17647540
    ## sample estimates:
    ## difference in location 
    ##              0.1076937

``` r
wilcox.test(data=SNV.APOBEC.df[which(SNV.APOBEC.df$APOBEC=="No"),],`Proportion TCW to TKW`~HPV_status,conf.int=T)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Proportion TCW to TKW by HPV_status
    ## W = 456.5, p-value = 0.6589
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.02406108  0.01768338
    ## sample estimates:
    ## difference in location 
    ##            -0.00342477

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
    ## 2 "HPV^{\n    \"+\"\n}" Yes              3408  
    ## 3 "HPV^{\n    \"+\"\n}" <NA>              116  
    ## 4 "HPV^{\n    \"−\"\n}" No               2089  
    ## 5 "HPV^{\n    \"−\"\n}" Yes             11852  
    ## 6 "HPV^{\n    \"−\"\n}" <NA>              150

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
    ## 2 "HPV^{\n    \"+\"\n}" Yes                                           2278
    ## 3 "HPV^{\n    \"+\"\n}" <NA>                                             0
    ## 4 "HPV^{\n    \"−\"\n}" No                                               0
    ## 5 "HPV^{\n    \"−\"\n}" Yes                                           5048
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
    ## 2 No     "HPV^{\n    \"−\"\n}"   126
    ## 3 Yes    "HPV^{\n    \"+\"\n}"    38
    ## 4 Yes    "HPV^{\n    \"−\"\n}"   260
    ## 5 <NA>   "HPV^{\n    \"+\"\n}"    21
    ## 6 <NA>   "HPV^{\n    \"−\"\n}"    40

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
    ## p-value = 0.04208
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  0.1701955 0.9831748
    ## sample estimates:
    ## odds ratio 
    ##  0.4351451

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

    ## Warning: package 'DT' was built under R version 3.4.3

<!--html_preserve-->

<script type="application/json" data-for="htmlwidget-55446c8e80f8d4729eb2">{"x":{"filter":"none","data":[["1","98","83","3","35","29","167","266","346","14","19","25","42","183","321","23","13","89","24","22","95","36","77","33","60","78","240","320","104","350","330","136","46","7","75","37","110","281","354","43","27","245","88","48","93","55","10","38","152","52","2","4","96","209","216","76","65","309","333","232","288","80","206","200","73","238","284","109","151","192","349","101","369","115","113","20","279","248","116","287","97","147","342","268","189","223","295","137","345","298","340","311","62","179","105","50","54","41","64","117","17","244","230","215","28","31","174","341","148","15","74","253","220","224","111","99","318","365","194","178","286","122","205","70","359","11","12","6","103","335","368","339","265","249","364","228","310","184","322","218","366","204","69","134","306","353","59","18","348","150","139","213","324","16","292","85","243","79","123","133","273","356","314","131","176","187","328","145","159","127","186","307","146","177","86","84","236","172","231","274","299","302","94","102","181","325","168","323","81","315","327","34","316","158","210","21","195","362","264","347","237","107","108","68","188","58","261","132","296","361","87","119","191","246","278","72","169","290","221","162","199","226","138","155","329","128","135","212","352","263","173","283","143","285","300","259","234","182","241","257","130","197","203","154","153","317","277","225","308","334","367","262","100","106","156","51","71","114","112","120","165","180","255","360","157","160","260","331","344","207","227","282","304","118","267","301","336","351","363","129","144","247","121","164","250","280","332","5","8","9","26","30","32","39","40","44","45","47","49","53","56","57","61","63","66","67","82","90","91","92","124","125","126","140","141","142","149","161","163","166","170","171","175","185","190","193","198","201","202","208","211","217","219","222","229","233","235","239","242","251","252","254","256","258","269","270","271","272","275","276","289","291","293","294","303","305","312","313","319","326","337","338","343","355","357","358","370","371"],["PIK3CA E545K","TP53 R283P","TP53 R196P","PIK3CA E542K","EPHA2 S419*","EP300 D1399N","EPHA6 D952H","HIST1H4C R68P","VAMP1 G 6579803 C NCSNV","TP53 E298*","FBXW7 R505G","TP53 Y163C","TP53 E294*","FAT1 E1994*","EP300 E1514K","TP53 P278S","HRAS G13V","TP53 E56*","TP53 H193L","MAPK1 E322K","TP53 S106R","TGFBR2 E544K","TP53 E285K","NFE2L2 E79K","RHOA E40Q","TP53 Q144*","HLA-B Q250E","EP300 C1164Y","CDKN2A E69*","CCDC50 A 191098615 C NCSNV","LAMB4 R1585P","NCOR1 E2306K","CDKN2A E88*","TP53 R213*","TP53 E68*","CASP8 Q524*","KIAA1109 I225M","HIST2H2BE F71L","ATP6AP2 E119Q","TP53 Q192*","CDKN2A G 21968242 A NCSNV","ELF4 S415L","TP53 G266E","NFE2L2 D29H","TP53 H193P","FAT1 S3373*","TP53 Y220C","TP53 V173M","ZNF43 G431A","MYC S161L","CDKN2A R80*","TP53 R175H","TP53 Q331*","CAPN2 Q413E","EPHA3 T802R","TP53 R158L","MAP2 E865*","TSHZ3 E841*","ZZZ3 R5Q","ANKHD1-EIF4EBP3 L376V","SMTN E963K","TP53 Q136*","IL17A E83Q","DEPDC5 C 32179890 T NCSNV","TP53 R249S","CNNM3 E422K","C3orf24 F33L","MUC16 E1409Q","ANKRD11 S769*","TGFBR2 S320*","PARP2 M432I","CDKN2A D108Y","FLOT1 G 30698343 C NCSNV","HUWE1 E3771K","AP3D1 Q1077Q","RAC1 A178V","ESRRA D219N","HCRTR2 D100Y","HUWE1 E4177K","ALPK3 Q860Q","TP53 C242F","FAM21C D460N","ARHGAP9 G 57868406 A NCSNV","NRAP R1205W","PIK3CA E726K","AP1G2 D243N","CCDC135 E313K","HSD17B1 C 40704915 T NCSNV","RRP12 E1204K","MRAP E28E","FABP3 C 31838696 T NCSNV","KRT14 G 39740666 A NCSNV","SLC44A2 F459F","FAT1 S1200*","RNASE3 R128T","MYH9 E530K","NOTCH1 E455K","TP53 G245V","PADI3 P481P","NFE2L2 E79Q","TP53 G245S","SLC27A1 F620F","CASP8 R292W","GIGYF2 S221*","CDKN2A W110*","TP53 H179R","WTAP E346K","RTKN R556C","FAM21C E1034K","TP53 R273C","TP53 C176Y","MYH7 E149K","HRAS Q61L","B2M E97K","DNAH5 Q1797E","ZNF253 C 20003622 G NCSNV","ARMC5 C 31477389 G NCSNV","ANXA6 R231Q","ITPR3 I2503I","FAT1 R885*","LRRN3 L190V","MGAM S703F","ALPK1 Q551*","APC R283Q","TMCO6 C 140021368 T NCSNV","TP53 R196*","TP53 R248Q","CDKN2A R58*","CDKN2A G 21971208 A NCSNV","PSG4 E263K","DDX23 L710L","LDLR M708I","UMODL1 L581L","PRKDC F1244F","ZNF710 F388F","CASP8 T 202137501 C NCSNV","SLC8A2 L195L","YSK4 Q798H","CDT1 C 88870221 T NCSNV","WDR70 E8K","RAP1B G12E","PSG8 S395F","ITGAV F51F","NCOR1 P1081L","RHBDF1 S285L","VCX2 A 8138014 T NCSNV","CASP8 R494*","TP53 R248W","OR4F15 R183W","PLCG2 L1038L","PGBD1 L412L","DCLK1 R23Q","DUSP6 E81K","TP53 R306*","MED29 L2L","TP53 R280G","HOXB5 E170K","TP53 P151S","IFFO1 L292L","NRXN3 G 80158584 A NCSNV","CDH7 Q225K","SMPDL3B F275F","ZNF234 G 44662400 C NCSNV","TMEM215 S168L","PLEKHG4 E744K","PIK3CA M1043V","NCOA1 L1361L","MYH9 E457K","PIK3CG R544*","FBXW7 R543G","PIK3CA R88Q","IFT140 E664K","ZNF510 S321*","NLRP5 E751K","TP53 Y234C","TP53 V143M","IKZF1 Q156Q","TNN D1091N","BTBD8 F178F","FCRL3 V288V","ECE1 A243A","TDRD6 E1997*","TP53 H179Y","CDKN2A G 21970900 T NCSNV","FAT1 Q1694*","C7orf42 C286R","TSPYL5 S193L","MBOAT7 R424W","TP53 C238F","MAGEB3 V75A","C9orf135 T222I","HIST1H3C K37M","TBX21 E494K","AKAP9 Q3402K","PAMR1 I580I","HRAS G12S","LRRC4C E637Q","SP140L D419N","IFNGR1 R123*","CNKSR1 L343L","KIAA0100 V1227V","FLG A2503V","PHLDB2 V1139V","DPP10 T 116257182 A NCSNV","PIK3CA Q546R","CASP8 R472*","SLMAP A 57743311 G NCSNV","DENND5B G988E","DNAH6 R442H","CEP68 C 65296555 T NCSNV","TP53 Q331Q","PRDM9 R77Q","NRD1 E181E","PCDHGA12 A329A","OR2D3 R70C","KCNA4 E280K","OR4C15 M351L","DHX15 R399H","HRAS G12D","PDE3A L275L","ZSCAN22 A321A","MCF2L R928Q","FUNDC2 F70F","NSD1 R788*","GNPDA2 G 44724233 A NCSNV","HEATR7B2 V1189V","NCOR1 R1561Q","WSCD2 P43P","KCNK9 L122L","MAP3K7 Y300C","PCDHA10 A481T","ITGB1 D158N","PKHD1L1 R1654Q","CCDC110 R111C","PNKD G 219137576 A NCSNV","BCAN R579Q","CTNNA2 A67T","FAT1 R3400*","HLA-B R7*","ATP10B R311Q","C6orf118 T328T","MYH13 K1852N","GPR158 L9L","KIAA1407 R673W","USH2A N3137N","MFSD10 R138Q","HLA-A R7*","HTR1A L166L","SMARCA4 P913L","NXPH1 T25T","CORO2B A156T","EPHB3 R420H","GRM8 G50D","MYO15A R3294W","THSD7A C728F","PCDHA6 S608S","COL22A1 P1601L","ASXL1 R693*","CSMD3 R3127Q","CDH10 R128R","ZNF208 S437S","FAT1 R937*","SPRR2B P69P","PTOV1 L99L","THSD7A R1046C","PIK3CG A156V","HHIPL2 F67F","VPS39 P567P","SFI1 R821Q","LY96 D100D","FCAR V233M","MYO3A C 26500926 T NCSNV","NADSYN1 D265D","NRXN1 R1256H","TAF1L R1397C","C9 D312N","RPA1 R31H","MTTP G722G","PAQR8 F255F","HCN1 V381V","KSR2 T555M","CDH9 T737M","PCDHA13 D481D","SALL1 P1109P","MMP16 R154H","C12orf60 T120T","GFRA3 A208A","PIK3CA H1047R","TP53 R282W","TP53 R273H","TP53 V157F","TP53 R110L","TP53 Y236C","TP53 R337L","TP53 R337C","TP53 T125T","TP53 G 7578176 T NCSNV","CDKN2A E120*","KIR3DL2 V405I","MAGEB1 P274P","ADAMTS18 H78H","HRAS G13R","FGFR3 S249C","PRB3 R39*","GABRA2 V74M","NLRP8 V428I","TP53 L201*","TP53 V272M","TP53 T 7577017 G NCSNV","TP53 G 7579311 T NCSNV","PFAS Y1202Y","MAEL A303V","TNIK R949H","CREBBP R1446C","SLCO1C1 S362R","ADCY1 R346H","CEP164 N79N","GUCY2C S606S","TMTC2 T409R","PREX1 T1402M","NOTCH1 A465T","NOTCH1 R353C","ZFHX4 G754R","PLXND1 G913R","BZRAP1 R1376K","TGFBR2 R553C","MYOM1 R63Q","SGCE R137H","CNTNAP2 C745C","KDM6A R519*","OR4C6 A214A","PDZD2 R953*","LINGO2 P410T","HRAS G12C","CASP8 R127*","TMEM132D A56V","EDIL3 Q173K","HIST1H3I K37M","IVD G408W","ATP7B V831I","OR11H12 Q30K","MYH7 V411I","PCDHA4 N451N","PCDH17 N579N","SDR9C7 A104A","KCNT2 R541*","FBXL13 R308W","TCF21 A124A","MASP1 T215A","SLC24A4 E137K","FAM190B Y766C","PPP4R4 T89M","CHCHD4 T92M","SLC34A1 R507H","TRHR V272L","DUOXA1 R332C","OR2T33 D153D","MIA C 41281377 T NCSNV","SLC5A7 G336C","HPN A 35557003 T NCSNV","SYT6 R164H","OR8I2 A67S","FXYD5 R145W","H3F3B K37M","OBP2A T107M","POLI V117L","GABBR2 W149*","CRYGS R36H"],[0.9,0.52,0.88,0.75,1,0.38,0.76,0.96,0.76,0.03,0.03,0.07,0.33,0.53,0.51,0.29,0.09,0.65,0.02,0.36,0.38,0.96,0.87,0.82,0.98,0.36,0.99,0.09,0.68,0.02,0.02,0.97,0.29,0.03,0.37,0.89,1,0.92,0.99,0.4,0.2,0.27,0.44,0.97,0.02,0.98,0.01,0.07,0.71,0.79,0.01,0.01,0.64,1,0.06,0.03,0.11,0.17,0.15,0.99,0.62,0.15,1,0.88,0.06,0.48,0.98,0.98,0.96,0.99,0.9,0.05,0.88,0.95,0.91,0.02,0.95,0.56,0.88,0.89,0.06,0.73,0.13,0.63,0.81,0.34,0.43,0.97,0.56,0.9,0.32,0.36,0.3,0.99,0.74,0.23,0.36,0.02,0.29,0.93,0.02,0.5,0.52,0.5,0.02,0.01,0.71,0.33,0.28,0.02,0.07,0.53,0.07,0.96,0.98,0.92,0.5,0.41,0.54,0.5,0.97,0.33,0.8,0.24,0.79,0.01,0.01,0.01,0.28,0.88,0.98,0.72,0.41,0.33,0.21,0.06,0.94,0.99,0.28,0.32,0.32,0.4,0.34,0.07,0.17,0.03,0.12,0.01,0.42,0.65,0.99,0.41,0.39,0.01,0.99,0.01,0.26,0.02,0.89,0.49,0.09,0.28,0.45,0.21,0.19,0.03,0.74,0.17,0.45,0.02,0.12,0.16,0.36,0.47,0.01,0.02,0.64,0.35,0.47,0.89,0.13,0.04,0.02,0.03,0.26,0.02,0.15,0.16,0.01,0.02,0.04,0.01,0.14,0.02,0.29,0.01,0.49,0.07,0.06,0.19,0.35,0.04,0.03,0.01,0.01,0.03,0.01,0.03,0.01,0.02,0.01,0.21,0.01,0.08,0.09,0.24,0.03,0.05,0.03,0.61,0.19,0.03,0.06,0.02,0.02,0.53,0.01,0.12,0.32,0.01,0.04,0.04,0.06,0.02,0.03,0.04,0.05,0.03,0.01,0.02,0.05,0.03,0.06,0.02,0.01,0.01,0.02,0.44,0.01,0.06,0.02,0.01,0.01,0.02,0.01,0.03,0.01,0.01,0.03,0.02,0.02,0.01,0.19,0.03,0.03,0.02,0.01,0.01,0.01,0.02,0.02,0.01,0.02,0.02,0.01,0.01,0.01,0.01,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.1,0.48,0.12,0.25,0,0.62,0.24,0.04,0.24,0.97,0.97,0.93,0.67,0.47,0.49,0.71,0.91,0.35,0.98,0.64,0.62,0.04,0.13,0.18,0.02,0.64,0.01,0.91,0.32,0.98,0.98,0.03,0.71,0.97,0.63,0.11,0,0.08,0.01,0.6,0.8,0.73,0.56,0.03,0.98,0.02,0.99,0.93,0.29,0.21,0.99,0.99,0.36,0,0.94,0.97,0.89,0.83,0.85,0.01,0.38,0.85,0,0.12,0.94,0.52,0.02,0.02,0.04,0.01,0.1,0.95,0.12,0.05,0.09,0.98,0.05,0.44,0.12,0.11,0.94,0.27,0.87,0.37,0.19,0.66,0.57,0.03,0.44,0.1,0.68,0.64,0.7,0.01,0.26,0.77,0.64,0.98,0.71,0.07,0.98,0.5,0.48,0.5,0.98,0.99,0.29,0.67,0.72,0.98,0.93,0.47,0.93,0.04,0.02,0.08,0.5,0.59,0.46,0.5,0.03,0.67,0.2,0.76,0.21,0.99,0.99,0.99,0.72,0.12,0.02,0.28,0.59,0.67,0.79,0.94,0.06,0.01,0.72,0.68,0.68,0.6,0.66,0.93,0.83,0.97,0.88,0.99,0.58,0.35,0.01,0.59,0.61,0.99,0.01,0.99,0.74,0.98,0.11,0.51,0.91,0.72,0.55,0.79,0.81,0.97,0.26,0.83,0.55,0.98,0.88,0.84,0.64,0.53,0.99,0.98,0.36,0.65,0.53,0.11,0.87,0.96,0.98,0.97,0.74,0.98,0.85,0.84,0.99,0.98,0.96,0.99,0.86,0.98,0.71,0.99,0.51,0.93,0.94,0.81,0.65,0.96,0.97,0.99,0.99,0.97,0.99,0.97,0.99,0.98,0.99,0.79,0.99,0.92,0.91,0.76,0.97,0.95,0.97,0.39,0.81,0.97,0.94,0.98,0.98,0.47,0.99,0.88,0.68,0.99,0.96,0.96,0.94,0.98,0.97,0.96,0.95,0.97,0.99,0.98,0.95,0.97,0.94,0.98,0.99,0.99,0.98,0.56,0.99,0.94,0.98,0.99,0.99,0.98,0.99,0.97,0.99,0.99,0.97,0.98,0.98,0.99,0.81,0.97,0.97,0.98,0.99,0.99,0.99,0.98,0.98,0.99,0.98,0.98,0.99,0.99,0.99,0.99,0.98,0.99,0.99,0.99,0.99,0.99,0.99,0.98,0.99,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.99,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[24,2,2,18,4,5,2,2,2,6,6,5,3,2,2,5,7,2,5,6,2,4,2,4,3,2,2,2,2,2,2,2,3,10,2,4,2,2,2,3,5,2,2,3,2,3,7,3,2,3,19,13,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,3,3,3,3,2,6,2,2,2,9,4,2,2,2,6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,7,7,11,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,6,2,2,2,2,2,6,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,4,2,2,2,6,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,13,9,7,5,4,4,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],[17515.33,323226.81,115043.21,13223.12,28589.29,47668.12,53738.03,30784.96,26410.66,219242.74,166055.56,69814.95,22640.27,21038.79,21432.92,14971.91,31888.02,12325.21,153732.99,6153.15,16080.05,2552.98,5522.17,2865.77,3194.62,12137.46,4294.43,46787.85,6091.72,206418.36,206283.91,4164.17,9127.34,26447.83,10498.82,2175.09,3605.73,3205.09,2847.6,4635.1,5541.62,10108.6,6026.17,1630.6,117702.99,1596.37,60496.55,19691.55,2898.88,1709.38,21258.05,30734.49,3096.51,1963.68,32552.33,64994.97,17374.77,11052.37,12130.12,1833.55,2864.03,11685.29,1679.69,1903.81,26801.24,3315.13,1617.23,1595.08,1600.66,1541.44,1678.99,29678.94,1681.61,1525.65,1552.54,23088.04,1437.62,2410.54,1525.65,1495.62,21743.94,1782.63,9817.14,1946.21,1495.75,3557.52,2799.15,1227.22,2123.09,1316.97,3699.86,3248.44,2570.96,1065.46,1419.72,3027.06,1934.13,34595.5,2376.21,1088.16,16843.45,1991.28,1892.55,1956.02,10822.43,44655.18,1251.09,2657.09,3095.7,14388.51,12137.46,1592.61,11934.87,836.29,806.26,849.98,1560.08,1904.33,1429.41,1535.28,749.75,2188.13,888.57,2878.54,875.44,19610.37,19610.37,12415.61,2427.09,771.94,673.83,858.66,1495.23,1840.78,2842.35,9889.85,612.87,572.12,2028.67,1736.89,1724.25,1372.94,1585.33,7613.17,3081.51,17342.04,2835.87,16843.45,1167.38,737.77,480.87,1122.1,1132.67,14388.51,436.07,42892.58,1618.03,20403.18,435.14,794.26,4288.86,1330.64,818.47,1740.47,1899.76,11775.96,468.72,2020.04,764.31,15662.01,2597.51,1952.45,857.82,619.89,28100.4,13155.07,405.59,733.26,548.04,281.05,1902.58,6197.34,11685.29,7714.5,884.08,11362.32,1423.9,1327.92,21119.89,10016.62,4924.95,9530.63,1372.57,9441.81,614.14,5648.71,346.82,2272.45,2638.59,818.98,427.89,3682.48,4633.13,13524.42,13758.91,2985.56,13049.69,4274.27,12470.13,6256.21,12280.1,573.03,11679.36,1451.58,1223.52,446.99,3535.53,2025.31,3316.8,158.53,505.79,3134.44,1491.82,4351.44,4033.01,145.07,7613.17,644.21,237.9,6597.42,1431.04,1424.19,909.37,2807.35,1851.01,1331.62,975.91,1616.32,4542.18,2136.96,774.28,1343.63,660.85,1902.67,3294.25,3398.95,1379.12,58.15,2656.89,420.35,1326.32,2350.61,2119.12,1065.91,1989.79,402.54,1737.84,1823.87,504.15,740.02,812.98,1616.32,87.31,554.16,442.55,687.38,1280.53,1411.44,1364.79,479.5,534.04,1013.58,496.02,369.99,700.09,873.6,755.5,898.77,337.32,565.48,626.5,414.09,193.52,337.27,341.4,139.93,373.43,69531.79,25110.48,12381.96,213699.71,171313.73,55967.6,97289.62,7239.03,12797.78,28671.43,9127.34,1129.34,258.31,455.16,90309.54,49563.7,976.75,819.33,726.7,83207.15,20403.18,117702.99,19154.13,1175.59,613.95,808.04,1917.38,1316.24,1036.12,1369,1369.85,44243.24,1466.55,1160.91,1358.98,602.19,3787.61,2499.03,2338.38,2161.61,2092.99,452.09,1857.84,103.4,2007.29,3045.71,27938.83,1992.45,586.71,2796.98,4317.45,12032.61,1124.58,2297.22,1058.65,632.08,322.51,939.08,611.5,1376.57,128.15,11122.7,1130.04,61282.9,1149.48,1180.4,4328.89,4169.39,1660.33,116.07,1516.82,2861.91,9038.72,1223.59,3995.85,1612.55,22104.92,814.01,9363.17,1647.51,1395.01],[720.63,640.3,385.67,340.02,217.82,172.51,155.58,112.58,76.47,75.17,56.93,46.54,42.69,42.48,41.64,41.35,38.27,30.52,29.28,25.32,23.28,18.67,18.3,17.9,17.89,16.65,16.2,16.04,15.78,15.73,15.72,15.39,15.13,15.11,14.8,14.75,13.74,11.23,10.74,10.59,10.56,10.4,10.1,9.04,8.97,8.94,8.07,7.88,7.84,7.72,7.69,7.61,7.55,7.48,7.44,7.43,7.28,7.16,6.93,6.92,6.76,6.68,6.4,6.38,6.13,6.06,6.04,5.95,5.85,5.81,5.76,5.65,5.64,5.52,5.38,5.28,5.2,5.14,5.11,5.07,4.97,4.96,4.86,4.67,4.62,4.61,4.59,4.53,4.53,4.52,4.51,4.46,4.41,4.02,4,3.98,3.98,3.95,3.94,3.86,3.85,3.79,3.75,3.73,3.71,3.4,3.38,3.34,3.3,3.29,3.24,3.22,3.18,3.06,3.01,2.98,2.97,2.97,2.94,2.92,2.77,2.75,2.71,2.63,2.63,2.61,2.61,2.6,2.59,2.59,2.52,2.36,2.34,2.31,2.27,2.26,2.19,2.16,2.16,2.12,2.1,2.09,2.05,2.03,2,1.98,1.94,1.92,1.87,1.83,1.81,1.75,1.68,1.64,1.64,1.63,1.6,1.55,1.48,1.48,1.47,1.42,1.4,1.39,1.38,1.35,1.32,1.31,1.31,1.19,1.19,1.19,1.18,1.11,1.07,1,0.99,0.98,0.98,0.95,0.94,0.94,0.89,0.88,0.88,0.87,0.81,0.81,0.8,0.76,0.75,0.73,0.73,0.72,0.68,0.65,0.65,0.61,0.6,0.59,0.57,0.56,0.53,0.52,0.52,0.51,0.5,0.49,0.48,0.48,0.47,0.46,0.44,0.44,0.42,0.41,0.4,0.39,0.38,0.37,0.37,0.36,0.34,0.33,0.31,0.29,0.29,0.29,0.29,0.25,0.22,0.22,0.21,0.21,0.21,0.2,0.19,0.18,0.17,0.16,0.15,0.15,0.15,0.14,0.13,0.13,0.11,0.1,0.1,0.1,0.1,0.09,0.08,0.08,0.08,0.07,0.07,0.07,0.06,0.06,0.06,0.06,0.06,0.06,0.05,0.05,0.05,0.05,0.05,0.04,0.04,0.04,0.04,0.03,0.03,0.03,0.03,0.03,0.03,0.02,0.02,0.02,0.01,0.01,0.01,0.01,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[80.07,591.04,52.59,113.34,0,281.47,49.13,4.69,24.15,2430.46,1840.84,618.36,86.68,37.67,40.01,101.24,386.91,16.43,1434.84,45.01,37.98,0.78,2.73,3.93,0.37,29.59,0.16,162.2,7.43,770.63,770.13,0.48,37.03,488.66,25.2,1.82,0,0.98,0.11,15.89,42.22,28.11,12.86,0.28,439.42,0.18,798.55,104.65,3.2,2.05,761.65,753.43,4.25,0,116.57,240.17,58.91,34.95,39.28,0.07,4.15,37.84,0,0.87,95.97,6.57,0.12,0.12,0.24,0.06,0.64,107.41,0.77,0.29,0.53,258.59,0.27,4.04,0.7,0.63,77.86,1.83,32.54,2.74,1.08,8.94,6.08,0.14,3.56,0.5,9.58,7.92,10.28,0.04,1.41,13.32,7.07,193.73,9.64,0.29,188.65,3.79,3.46,3.73,181.82,336.83,1.38,6.78,8.49,161.15,43,2.85,42.28,0.13,0.06,0.26,2.97,4.28,2.5,2.92,0.09,5.58,0.68,8.33,0.7,258.86,258.86,257.54,6.66,0.35,0.05,0.92,3.36,4.7,8.55,35.42,0.14,0.02,5.56,4.5,4.47,3.14,3.99,26.97,9.74,64.08,14.26,190.57,2.58,0.98,0.02,2.52,2.63,162.8,0.02,161.77,4.56,76.17,0.18,1.54,14.87,3.65,1.71,5.24,5.86,43.51,0.46,6.39,1.6,58.47,8.71,6.25,2.09,1.25,105.98,49.11,0.56,1.82,1.11,0.12,6.31,22.66,43.63,28.51,2.49,42.42,4.61,4.25,79.65,37.4,18.01,71.89,4.5,35.25,1.66,63.91,0.67,8.05,9.45,2.53,1.06,13.47,17.12,51.01,51.89,16.55,49.22,15.79,47.03,23.36,46.31,1.72,44.05,5.09,4.24,1.29,13.06,7.33,12.26,0.24,1.56,11.58,5.34,16.25,15.06,0.26,28.71,2.16,0.62,24.88,5.23,5.21,3.26,10.48,6.84,4.87,3.53,5.97,17.13,7.98,2.8,4.97,2.37,7.1,12.42,12.82,5.15,0.12,10.02,1.51,4.95,8.87,7.99,3.98,7.5,2.23,6.55,6.88,1.86,2.76,3.04,6.1,0.27,2.05,1.64,2.57,4.83,5.32,5.15,1.79,1.99,3.82,1.85,1.38,2.64,3.29,2.85,3.39,1.26,2.13,2.36,1.56,0.73,1.27,1.29,0.52,1.41,1721.74,430.47,165.09,2035.24,1305.25,426.42,555.94,41.37,73.13,163.84,52.16,6.45,1.48,2.6,516.05,283.22,5.58,3.12,2.77,316.98,77.73,448.39,72.97,4.48,2.34,3.08,7.3,5.01,3.95,5.22,5.22,168.55,5.59,4.42,5.18,2.29,14.43,9.52,8.91,8.23,7.97,1.72,7.08,0.39,7.65,11.6,106.43,7.59,2.24,10.66,16.45,45.84,4.28,8.75,4.03,2.41,1.23,3.58,2.33,5.24,0.49,42.37,4.3,233.46,4.38,4.5,16.49,15.88,6.33,0.44,5.78,10.9,34.43,4.66,15.22,6.14,84.21,3.1,35.67,6.28,5.31]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Short name<\/th>\n      <th>Median APOBEC contribution<\/th>\n      <th>Median NonAPOBEC contribution<\/th>\n      <th>Frequency<\/th>\n      <th>Selection intensity<\/th>\n      <th>Selection from APOBEC<\/th>\n      <th>Selection from NonAPOBEC<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
<!--/html_preserve-->
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
    ##  [1] DT_0.4                dendextend_1.6.0      ggdendro_0.1-20      
    ##  [4] viridis_0.5.0         viridisLite_0.3.0     gridExtra_2.3        
    ##  [7] cowplot_0.9.2         bindrcpp_0.2          forcats_0.2.0        
    ## [10] stringr_1.2.0         dplyr_0.7.4           purrr_0.2.4          
    ## [13] readr_1.1.1           tidyr_0.8.0           tibble_1.4.2         
    ## [16] tidyverse_1.2.1       deconstructSigs_1.8.0 reshape2_1.4.3       
    ## [19] ggrepel_0.7.0         ggplot2_2.2.1         rtracklayer_1.36.6   
    ## [22] GenomicRanges_1.28.6  GenomeInfoDb_1.12.3   IRanges_2.10.5       
    ## [25] S4Vectors_0.14.7      BiocGenerics_0.22.1  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] colorspace_1.3-2                  class_7.3-14                     
    ##  [3] modeltools_0.2-21                 mclust_5.4                       
    ##  [5] rprojroot_1.3-2                   XVector_0.16.0                   
    ##  [7] rstudioapi_0.7                    flexmix_2.3-14                   
    ##  [9] mvtnorm_1.0-7                     lubridate_1.7.2                  
    ## [11] xml2_1.2.0                        mnormt_1.5-5                     
    ## [13] robustbase_0.92-8                 knitr_1.19                       
    ## [15] jsonlite_1.5                      Rsamtools_1.28.0                 
    ## [17] broom_0.4.3                       cluster_2.0.6                    
    ## [19] kernlab_0.9-25                    shiny_1.0.5                      
    ## [21] compiler_3.4.0                    httr_1.3.1                       
    ## [23] backports_1.1.2                   assertthat_0.2.0                 
    ## [25] Matrix_1.2-12                     lazyeval_0.2.1                   
    ## [27] cli_1.0.0                         htmltools_0.3.6                  
    ## [29] tools_3.4.0                       gtable_0.2.0                     
    ## [31] glue_1.2.0                        GenomeInfoDbData_0.99.0          
    ## [33] Rcpp_0.12.15                      Biobase_2.36.2                   
    ## [35] cellranger_1.1.0                  trimcluster_0.1-2                
    ## [37] Biostrings_2.44.2                 nlme_3.1-131                     
    ## [39] fpc_2.1-11                        crosstalk_1.0.0                  
    ## [41] psych_1.7.8                       rvest_0.3.2                      
    ## [43] mime_0.5                          XML_3.98-1.9                     
    ## [45] DEoptimR_1.0-8                    zlibbioc_1.22.0                  
    ## [47] MASS_7.3-48                       scales_0.5.0                     
    ## [49] BSgenome_1.44.2                   hms_0.4.1                        
    ## [51] SummarizedExperiment_1.6.5        yaml_2.1.16                      
    ## [53] stringi_1.1.6                     BiocParallel_1.10.1              
    ## [55] rlang_0.1.6                       pkgconfig_2.0.1                  
    ## [57] prabclus_2.2-6                    matrixStats_0.53.0               
    ## [59] bitops_1.0-6                      evaluate_0.10.1                  
    ## [61] lattice_0.20-35                   bindr_0.1                        
    ## [63] GenomicAlignments_1.12.2          htmlwidgets_1.0                  
    ## [65] labeling_0.3                      BSgenome.Hsapiens.UCSC.hg19_1.4.0
    ## [67] plyr_1.8.4                        magrittr_1.5                     
    ## [69] R6_2.2.2                          DelayedArray_0.2.7               
    ## [71] pillar_1.1.0                      haven_1.1.1                      
    ## [73] whisker_0.3-2                     foreign_0.8-69                   
    ## [75] RCurl_1.95-4.10                   nnet_7.3-12                      
    ## [77] modelr_0.1.1                      crayon_1.3.4                     
    ## [79] utf8_1.1.3                        rmarkdown_1.8                    
    ## [81] readxl_1.0.0                      digest_0.6.15                    
    ## [83] diptest_0.75-7                    xtable_1.8-2                     
    ## [85] httpuv_1.3.5                      munsell_0.4.3
