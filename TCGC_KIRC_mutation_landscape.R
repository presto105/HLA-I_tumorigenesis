library("maftools")
library("dplyr")
library("ggplot2")
library("ggpubr")
library("RColorBrewer")

####Input_data_setting####
df_FILE="TCGA_ann.tsv"
raw_df <- read.csv(df_FILE, header=T,sep='\t')
names(raw_df)[names(raw_df)=="age_at_initial_pathologic_diagnosis"] <- "age_at_diagnosis"


########## Data preprocessing ##########
filtered_df <- raw_df %>% 
  mutate(event = 1,
         Homo_at_least_1_locus = as.factor(ifelse(Homozygous_gene_count>=1,1,
                                                  ifelse(Homozygous_gene_count==0,0,NA))),
         Homo_at_least_2_locus = as.factor(ifelse(Homozygous_gene_count>=2,1,
                                                  ifelse(Homozygous_gene_count==0,0,NA))),
         Homo_at_3_locus = as.factor(ifelse(Homozygous_gene_count==3,1,
                                            ifelse(Homozygous_gene_count==0,0,NA))),
         Germline_variant = as.factor(ifelse(P_total_GV_count==0,0,1)),
         age_at_diagnosis = as.numeric(as.character(age_at_diagnosis)),
         race2 = ifelse(race=="WHITE",1,0),
         gender = as.factor(as.character(gender)),
         tumor_stage = as.factor(as.character(tumor_stage)),
         prior_malignancy = as.factor(as.character(prior_malignancy)),
         tumor = as.factor(as.character(project_id)),
         stage = ifelse(grepl("IV",tumor_stage),4,
                        ifelse(grepl("III",tumor_stage),3,
                               ifelse(grepl("I/II",tumor_stage),NA,
                                      ifelse(grepl("II",tumor_stage),2,
                                             ifelse(grepl("I",tumor_stage),1,
                                                    ifelse(grepl(0,tumor_stage),1,NA)))))),
         VHL = relevel(as.factor(ifelse(is.na(VHL_mut_type), NA,
                                        ifelse(VHL_mut_type=='None','None','Mutation'))),ref='None'),
         PBRM1 = relevel(as.factor(ifelse(is.na(PBRM1_mut_type),NA,
                                          ifelse(PBRM1_mut_type=='None','None','Mutation'))),ref='None'),
         chr3p_VHL_binned = relevel(as.factor(ifelse(is.na(chr3p_VHL_binned), NA,
                                                     ifelse(chr3p_VHL_binned=='-1','-1','0'))),ref='0'),
         chr3p_PBRM1_binned = relevel(as.factor(ifelse(is.na(chr3p_PBRM1_binned), NA,
                                                       ifelse(chr3p_PBRM1_binned=='-1','-1','0'))),ref='0'),
         dNdS=ifelse(non_silent_TMB==0|silent_TMB==0, NA, non_silent_TMB/silent_TMB)
  ) %>% 
  filter(stage=='1'&prior_malignancy=="no"&!is.na(age_at_diagnosis)&!is.na(Homozygous_gene_count)) %>% 
  filter(Germline_variant==0&(ILC_Infection_status=="X"|is.na(ILC_Infection_status))) %>% 
  filter(tumor=='KIRC')

#################### 3p loss ##########################
c3p_loss_filtered_df <- filtered_df %>%
  mutate(VHL_PBRM1_class = ifelse(VHL=='Mutation'&PBRM1=='Mutation'&chr3p_PBRM1_binned==-1,'Both_biallelic',
                                  ifelse(VHL=='Mutation'&PBRM1!='Mutation'&chr3p_PBRM1_binned==-1,'VHL_biallelic',
                                         ifelse(VHL!='Mutation'&PBRM1=='Mutation'&chr3p_PBRM1_binned==-1,'PBRM1_biallelic',
                                                ifelse((VHL!='Mutation'&PBRM1!='Mutation'&chr3p_PBRM1_binned==-1),'3p_Haplo_loss',
                                                       ifelse(VHL!='Mutation'&PBRM1!='Mutation'&chr3p_PBRM1_binned!=-1, 'None',NA))))),
         VHL_biallelic = factor(ifelse(VHL=='Mutation'&chr3p_PBRM1_binned==-1,'Mutation','None'),levels=c('None','Mutation')),
         PBRM1_biallelic = factor(ifelse(PBRM1=='Mutation'&chr3p_PBRM1_binned==-1,'Mutation','None'),levels=c('None','Mutation'))) %>% 
  filter(!is.na(VHL)&!is.na(chr3p_PBRM1_binned)&!(chr3p_VHL_binned==-1&chr3p_PBRM1_binned!=-1))




#### Filtered ####
maf_data <- paste('Prefiltered_KIRC_mc3.v0.2.8.CONTROLLED',".maf",sep='')
gene_data <- paste("KIRC.markSMGs.tsv",sep='')

#### get gene list ####  
maf = read.maf(maf = maf_data, clinicalData = gene_data)

filtered_maf = subsetMaf(maf=maf, tsb=c3p_loss_filtered_df$submitter_id)
hotspot <- colnames(filtered_maf@clinical.data)[2:length(colnames(filtered_maf@clinical.data))]
hotspot <- c("X3p_loss", "Homozygous_cnt","diagnosis_age")

#### Clinical data_annotation color setting ####
col = RColorBrewer::brewer.pal(n = 9, name = 'Paired')
GF_cols <- c("#cccccc", "#56876D", "#56876D", "#2272A9", "#E2C044",
             "#040000", "#8E7DBE", "#EF8354", "#A8763E", "#A8763E")
names(GF_cols) = c("None","Frame_Shift_Del","Frame_Shift_Ins", "Missense_Mutation","Nonsense_Mutation",
                   "Multi_Hit", "Nonstop_Mutation","Splice_Site","In_Frame_Ins","In_Frame_Del")
anno_cols = list()
for (geneName in hotspot){
  anno_cols[[geneName]] = GF_cols
}

homo_cnt_cols <- c("#D2AB99", "#8DB38B", "#56876D", "#04724D")
names(homo_cnt_cols) <- c("0","1","2","3")
anno_cols[["Homozygous_cnt"]]  = homo_cnt_cols


min_age <- as.numeric(min(filtered_maf@clinical.data$diagnosis_age))
max_age <- as.numeric(max(filtered_maf@clinical.data$diagnosis_age))
interval <- max_age-min_age+1

rbPal <- colorRampPalette(c('red','blue'))
grad_color <- rbPal(interval)
names(grad_color) <- as.character(seq(min_age,max_age,1))
anno_cols[["diagnosis_age"]]  = grad_color

loh <- as.character(unique(filtered_maf@clinical.data$X3p_loss))
col <- c("#C93E4C","#B3B3B3","#B3B3B3")
names(col) <- loh
anno_cols[["X3p_loss"]] = col

#Shows sample summry.
getSampleSummary(filtered_maf)

pdf("KIRC_mutation_landscape.pdf", width=18, height=12)
  oncoplot(maf = filtered_maf, top = 2, removeNonMutated = F, drawRowBar = T, drawColBar = F,
           clinicalFeatures = hotspot, annotationColor = anno_cols, colors = GF_cols,
           sortByAnnotation = T, logColBar =F, writeMatrix = TRUE)
dev.off()