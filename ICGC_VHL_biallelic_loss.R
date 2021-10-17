rm(list=ls())

library(forestplot)
library(moonBook)
require(rms)
library(survival)
library(survminer)
library(dplyr)
library(Hmisc)
library(timereg)
library(openxlsx)
library(ggplot2)
library(coin)
library(tidyr)
library(KMsurv)

########## Set working directory ########## 
ICGC_HLA_genotype_FILE="ICGC_ann.tsv"

########## Data preprocessing ##########
annotated_df<- read.csv(ICGC_HLA_genotype_FILE, header=T,sep='\t')

filtered_df <- annotated_df %>% 
  mutate(event = 1,
         Homo_at_least_1_locus = as.factor(ifelse(Homozygous_gene_count>=1,1,
                                                  ifelse(Homozygous_gene_count==0,0,NA))),
         Homo_at_least_2_locus = as.factor(ifelse(Homozygous_gene_count>=2,1,
                                                  ifelse(Homozygous_gene_count==0,0,NA))),
         Homo_at_3_locus = as.factor(ifelse(Homozygous_gene_count==3,1,
                                            ifelse(Homozygous_gene_count==0,0,NA))),
         Homo_only_A = as.factor(ifelse(A_zygosity=='Homozygous'&B_zygosity!='Homozygous'&C_zygosity!='Homozygous',1,0)),
         Homo_only_B = as.factor(ifelse(A_zygosity!='Homozygous'&B_zygosity=='Homozygous'&C_zygosity!='Homozygous',1,0)),
         Homo_only_C = as.factor(ifelse(A_zygosity!='Homozygous'&B_zygosity!='Homozygous'&C_zygosity=='Homozygous',1,0)),
         age_at_diagnosis = as.numeric(as.character(age_at_diagnosis)),
         gender = as.factor(as.character(gender))
         ) %>% 
  filter(!is.na(Homozygous_gene_count)&!is.na(age_at_diagnosis))

########## Function ##########
cumulative_onset_plot <- function(data, allele, dir, filename, annotation=NULL){
  age_index_min <<-  floor(min(as.numeric(as.character(data$age_at_diagnosis)),na.rm=T))
  age_index_max <<- ceiling(max(as.numeric(as.character(data$age_at_diagnosis)),na.rm=T))
  interval <<- age_index_max-age_index_min
  
  form <- formula(paste0("Surv(age_at_diagnosis,event) ~ ",allele))
  fit <- do.call(survfit, args=list(formula=form, data=data))
  
  if(dir.exists(dir)==F){
    dir.create(dir)
  }
  
  surv_p <<- ggsurvplot(fit, conf.int=FALSE, fun="event", xlim=c(age_index_min,age_index_max+1), break.time.by=10,
                        pval=F, size = 1.2,#pval.coord=c(age_index_max-interval/4.2, 0.01), pval.size=7, 
                        fontsize=6, xlab="Age", ylab="Cumulative onset",  
                        font.title=15, legend.title=allele, legend=c(0.80,0.15), font.legend=c(15),
                        font.tickslab=c(15, "black"), font.x = c(15, "black"), font.y = c(15, "black"),
                        palette=c("#56B4E9","#E69F00"), cumevents=T) # count table
  if(!is.null(annotation)) surv_p$plot <- surv_p$plot + annotation
  
  pdf(paste0(dir,filename,".pdf"), width=12, height=10, onefile=F)
  print(surv_p)
  dev.off()
  
  tiff(paste0(dir,filename,".tiff"), width=12, height=10, units="in", res=100)
  print(surv_p)
  dev.off()
}

######### Accelerated failure time model (AFT model) #########
expected_time_weibull <- function(beta, log_scale){
  return(exp(beta)*gamma(1+exp(log_scale)))}

dist_select <- function(data, var){
  
  AFT_data <- data
  
  form <- formula(paste("Surv(age_at_diagnosis,event)~",paste(var,collapse="+")))
  
  AIC_res <- c()
  dist_class <- c("weibull", "exponential", "lognormal", "loglogistic")
  for(dist in dist_class){
    AIC_res <- c(AIC_res, AIC(survreg(formula=form, data=AFT_data, dist=dist)))
  }
  AIC_df <- data.frame(AIC_res, row.names=dist_class)
  min_dist <- dist_class[which(AIC_res==min(AIC_res))]
  
  return(list(min_dist=min_dist,AIC_df=AIC_df))
}

AFT_model <- function(data, tumor, var, ad_var, dist){
  
  AFT_data <- data
  
  
  form <- formula(paste("Surv(age_at_diagnosis,event)~",paste(c(var,ad_var),collapse="+")))
  model <- summary(survreg(formula=form, data=AFT_data, dist=dist))
  
  unadjust <- as.data.frame(model$table)
  unadjust[] <- lapply(unadjust,sprintf,fmt="%.3f")
  unadjust$p <- ifelse(unadjust$p<0.001,"<0.001",unadjust$p)
  
  form <- formula(paste("Surv(age_at_diagnosis,event)~",paste(c(var,ad_var),collapse="+")))
  model <- summary(survreg(formula=form, data=AFT_data, dist=dist))
  adjust <- as.data.frame(model$table)
  adjust[] <- lapply(adjust,sprintf,fmt="%.3f")
  adjust$p <- ifelse(adjust$p<0.001,"<0.001",adjust$p)
  
  return(list(unadjust=unadjust,adjust=adjust))
}

FR_merged_cumulative_onset_plot <- function(data, dist, dir, filename){

  AFT_data <- data
  if (length(unique(AFT_data$tumor_type))==1) {
    form <- formula(paste("Surv(age_at_diagnosis,event)~Homo_at_least_1_locus+gender"))
  } else {
    form <- formula(paste("Surv(age_at_diagnosis,event)~Homo_at_least_1_locus+gender+tumor_type"))  
  }
  
  form_LR <- formula(paste("Surv(age_at_diagnosis,event)~Homo_at_least_1_locus"))
  
  model <- summary(survreg(formula=form, data=AFT_data, dist=dist))$table
  
  n <- nrow(AFT_data)
  homo_n <- nrow(AFT_data %>% filter(Homo_at_least_1_locus=='1'))
  hetero_n <- nrow(AFT_data %>% filter(Homo_at_least_1_locus=='0'))
  
  time_ratio <- exp(model[2,1])
  lower <- exp(model[2,1]-qnorm(0.975)*model[2,2])
  upper <- exp(model[2,1]+qnorm(0.975)*model[2,2])
  hetero_age <- round(expected_time_weibull(model[1,1], model[3,1]),2)
  homo_age <- round(expected_time_weibull(model[1,1]+model[2,1], model[3,1]),2)
  diff_age <- hetero_age-homo_age
  
  p_value <- model[2,4]
  
  time_ratio <- sprintf("%.3f", time_ratio)
  lower <- sprintf("%.3f", lower)
  upper <- sprintf("%.3f", upper)
  
  time_ratio_with_CI <- paste0(time_ratio," (",lower,",",upper,")")
  p_value <- ifelse(p_value<0.001," < 0.001",sprintf("%.3f",p_value))
  
  LR <- survdiff(form_LR, AFT_data)
  p <- pchisq(LR$chisq,1,lower.tail=F)
  p <- ifelse(p_value<0.001," < 0.001",paste0("= ",sprintf("%.3f",p)))
  
  cumulative_onset_plot(data=data, allele="Homo_at_least_1_locus", 
                        dir=dir, filename=filename,
                        annotation=annotate("text", x = age_index_max-interval/8, y = 0.25, label = 
                                              paste0("TR (95% CI) = ",time_ratio_with_CI,"\n p ",p), parse = F, size = 6))
  res <- data.frame(homo_n, hetero_n, time_ratio_with_CI, p_value, hetero_age, homo_age, diff_age)
  return(res)
}

Boxplot <- function(data, start, end, dir, filename){
  box_p <- data %>% 
    ggplot(aes(x=Homo_at_least_1_locus,y=age_at_diagnosis,fill=Homo_at_least_1_locus)) +
    geom_boxplot(width=0.8,outlier.shape = NA) +
    geom_point(position =position_jitterdodge(jitter.width=0.2)) +
    stat_compare_means(size=6, label.x = 1.4, method = "wilcox.test",label = "p.format") +
    coord_cartesian(ylim = c(start, end)) +
    scale_fill_manual(values=c("#56B4E9","#E69F00")) +
    theme_classic() + theme(legend.position="bottom",text = element_text(size=15))
  
  pdf(paste0(dir, filename,".pdf"), width=4, height=8, onefile=F)
  print(box_p)
  dev.off()
  
  tiff(paste0(dir, filename,".tiff"), width=4, height=8, units="in", res=100)
  print(box_p)
  dev.off()
}


## homo median age
filtered_df %>%
  filter(VHL_biallelic_mut=='O') %>% 
  group_by(Homo_at_least_1_locus) %>% summarise(median = median(age_at_diagnosis, na.rm = TRUE))

######### Total stage #############
Homo_cnt_over10 <- filtered_df %>% group_by(tumor_type, Homo_at_least_1_locus) %>% 
  summarise(n=n()) %>% filter(Homo_at_least_1_locus==1&n>=10)
# tumor_types <- sort(as.character(unique(Homo_cnt_over10$tumor_type)))

ST_df <- filtered_df %>%
  filter(!is.na(age_at_diagnosis)&!is.na(Homozygous_gene_count))

y_start <- min(ST_df$age_at_diagnosis)
y_end <- max(ST_df$age_at_diagnosis)

AFT_regression <- data.frame()
############ group 1 ###############
ST_G1 <- ST_df %>% filter(tumor_type=='RECA-EU')
dist <- dist_select(ST_G1, "Homo_at_least_1_locus")[[1]]
Boxplot(data=ST_G1, start=y_start, end=y_end, dir="Boxplot/", filename="ICGC_ST_G1_Box")
ST_G1_tmp <- FR_merged_cumulative_onset_plot(data=ST_G1, dist=dist, dir="plot/", filename="ICGC_ST_G1_KM")
ST_G1_tmp['Group'] <- 'stage_Total_G1'
AFT_regression <- rbind(AFT_regression, ST_G1_tmp)

############ group 2 ###############
ST_G2 <- ST_df %>% filter(tumor_type!='RECA-EU')
dist <- dist_select(ST_G2, "Homo_at_least_1_locus")[[1]]
Boxplot(data=ST_G2, start=y_start, end=y_end, dir="Boxplot/", filename="ICGC_ST_G2_Box")
ST_G2_tmp <- FR_merged_cumulative_onset_plot(data=ST_G2, dist=dist, dir="plot/", filename="ICGC_ST_G2_KM")
ST_G2_tmp['Group'] <- 'stage_Total_G2'
AFT_regression <- rbind(AFT_regression, ST_G2_tmp)

############ group 3 ###############
ST_G3 <- ST_df %>% filter(VHL_biallelic_mut=='O')
dist <- dist_select(ST_G3, "Homo_at_least_1_locus")[[1]]
Boxplot(data=ST_G3, start=y_start, end=y_end, dir="Boxplot/", filename="ICGC_ST_G3_Box")
ST_G3_tmp <- FR_merged_cumulative_onset_plot(data=ST_G3, dist=dist, dir="plot/", filename="ICGC_ST_G3_KM")
ST_G3_tmp['Group'] <- 'stage_Total_G3'
AFT_regression <- rbind(AFT_regression, ST_G3_tmp)

ST_G3 %>%
  group_by(Homo_at_least_1_locus) %>% summarise(median = median(age_at_diagnosis, na.rm = TRUE))

############ group 4 ###############
ST_G4 <- ST_df %>% filter(is.na(VHL_biallelic_mut))
dist <- dist_select(ST_G4, "Homo_at_least_1_locus")[[1]]
Boxplot(data=ST_G4, start=y_start, end=y_end, dir="Boxplot/", filename="ICGC_ST_G4_Box")
ST_G4_tmp <- FR_merged_cumulative_onset_plot(data=ST_G4, dist=dist, dir="plot/", filename="ICGC_ST_G4_KM")
ST_G4_tmp['Group'] <- 'stage_Total_G4'
AFT_regression <- rbind(AFT_regression, ST_G4_tmp)

############ group 5 ###############
ST_G5 <- ST_df %>%
  filter(tumor_type=='RECA-EU') %>% 
  filter(is.na(VHL_biallelic_mut))
dist <- dist_select(ST_G5, "Homo_at_least_1_locus")[[1]]
Boxplot(data=ST_G5, start=y_start, end=y_end, dir="Boxplot/", filename="ICGC_ST_G5_Box")
ST_G5_tmp <- FR_merged_cumulative_onset_plot(data=ST_G5, dist=dist, dir="plot/", filename="ICGC_ST_G5_KM")
ST_G5_tmp['Group'] <- 'stage_Total_G5'
AFT_regression <- rbind(AFT_regression, ST_G5_tmp)


############ group 6 ###############
ST_G6 <- ST_df %>%
  filter(tumor_type!='RECA-EU') %>% 
  filter(VHL_biallelic_mut=='O')
dist <- dist_select(ST_G6, "Homo_at_least_1_locus")[[1]]
Boxplot(data=ST_G6, start=y_start, end=y_end, dir="Boxplot/", filename="ICGC_ST_G6_Box")
ST_G6_tmp <- FR_merged_cumulative_onset_plot(data=ST_G6, dist=dist, dir="plot/", filename="ICGC_ST_G6_KM")
ST_G6_tmp['Group'] <- 'stage_Total_G6'
AFT_regression <- rbind(AFT_regression, ST_G6_tmp)

write.xlsx(AFT_regression, "Stage_total_group_AFT_result.xlsx")


######### Stage 1 #############
sup_stage1 <- c(0,'I','IA','IB','Ic','T1aN0M0','UICC:1a','UICC:1b','BCLC:0','BCLC:A')
stage1 <- c('T1/N0/M0', 'T1N0M0', '1', 'IA/IB', 'IB',
            'IA', 'I', 'T1N0', '1B', '1A', 'T1cN0M0',
            'T1cN0(sn)', 'T0N0M0', 'I,N', '0', 'IAE',
            'Stage 1','T1bN0M0', 'T1aN0M0','T1NXM0')

S1_df <- filtered_df %>%
  filter(stage_sup %in% sup_stage1)

tmp_df <- filtered_df %>%
  filter(!stage_sup %in% sup_stage1) %>% 
  filter(tumor_stage %in% stage1)
  
S1_df <- rbind(S1_df,tmp_df) %>%
  filter(!is.na(age_at_diagnosis)&!is.na(Homozygous_gene_count))

y_start <- min(S1_df$age_at_diagnosis)
y_end <- max(S1_df$age_at_diagnosis)

AFT_regression <- data.frame()
############ group 1 ###############
S1_G1 <- S1_df %>% filter(tumor_type=='RECA-EU')
dist <- dist_select(S1_G1, "Homo_at_least_1_locus")[[1]]
Boxplot(data=S1_G1, start=y_start, end=y_end, dir="Boxplot/", filename="ICGC_S1_G1_Box")
S1_G1_tmp <- FR_merged_cumulative_onset_plot(data=S1_G1, dist=dist, dir="plot/", filename="ICGC_S1_G1_KM")
S1_G1_tmp['Group'] <- 'stage1_G1'
AFT_regression <- rbind(AFT_regression, S1_G1_tmp)

############ group 2 ###############
S1_G2 <- S1_df %>% filter(tumor_type!='RECA-EU')
dist <- dist_select(S1_G2, "Homo_at_least_1_locus")[[1]]
Boxplot(data=S1_G2, start=y_start, end=y_end, dir="Boxplot/", filename="ICGC_S1_G2_Box")
S1_G2_tmp <- FR_merged_cumulative_onset_plot(data=S1_G2, dist=dist, dir="plot/", filename="ICGC_S1_G2_KM")
S1_G2_tmp['Group'] <- 'stage1_G2'
AFT_regression <- rbind(AFT_regression, S1_G2_tmp)

############ group 3 ###############
S1_G3 <- S1_df %>% filter(VHL_biallelic_mut=='O')
dist <- dist_select(S1_G3, "Homo_at_least_1_locus")[[1]]
Boxplot(data=S1_G3, start=y_start, end=y_end, dir="Boxplot/", filename="ICGC_S1_G3_Box")
S1_G3_tmp <- FR_merged_cumulative_onset_plot(data=S1_G3, dist=dist, dir="plot/", filename="ICGC_S1_G3_KM")
S1_G3_tmp['Group'] <- 'stage1_G3'
AFT_regression <- rbind(AFT_regression, S1_G3_tmp)

S1_G3 %>%
  group_by(Homo_at_least_1_locus) %>% summarise(median = median(age_at_diagnosis, na.rm = TRUE))

############ group 4 ###############
S1_G4 <- S1_df %>% filter(is.na(VHL_biallelic_mut))
dist <- dist_select(S1_G4, "Homo_at_least_1_locus")[[1]]
Boxplot(data=S1_G4, start=y_start, end=y_end, dir="Boxplot/", filename="ICGC_S1_G4_Box")
S1_G4_tmp <- FR_merged_cumulative_onset_plot(data=S1_G4, dist=dist, dir="plot/", filename="ICGC_S1_G4_KM")
S1_G4_tmp['Group'] <- 'stage1_G4'
AFT_regression <- rbind(AFT_regression, S1_G4_tmp)

############ group 5 ###############
S1_G5 <- S1_df %>%
  filter(tumor_type=='RECA-EU') %>% 
  filter(is.na(VHL_biallelic_mut))
dist <- dist_select(S1_G5, "Homo_at_least_1_locus")[[1]]
Boxplot(data=S1_G5, start=y_start, end=y_end, dir="Boxplot/", filename="ICGC_S1_G5_Box")
S1_G5_tmp <- FR_merged_cumulative_onset_plot(data=S1_G5, dist=dist, dir="plot/", filename="ICGC_S1_G5_KM")
S1_G5_tmp['Group'] <- 'stage1_G5'
AFT_regression <- rbind(AFT_regression, S1_G5_tmp)

S1_G5 %>%
  group_by(Homo_at_least_1_locus) %>% summarise(median = median(age_at_diagnosis, na.rm = TRUE))

############ group 6 ###############
S1_G6 <- S1_df %>%
  filter(tumor_type!='RECA-EU') %>% 
  filter(VHL_biallelic_mut=='O')
dist <- dist_select(S1_G6, "Homo_at_least_1_locus")[[1]]
Boxplot(data=S1_G6, start=y_start, end=y_end, dir="Boxplot/", filename="ICGC_S1_G6_Box")
S1_G6_tmp <- FR_merged_cumulative_onset_plot(data=S1_G6, dist=dist, dir="plot/", filename="ICGC_S1_G6_KM")
S1_G6_tmp['Group'] <- 'stage1_G6'
AFT_regression <- rbind(AFT_regression, S1_G6_tmp)

write.xlsx(AFT_regression, "Stage1_group_AFT_result.xlsx")
