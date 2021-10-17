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
annotated_df <- read.csv(ICGC_HLA_genotype_FILE, header=T,sep='\t')

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
         gender = as.factor(as.character(gender)),
         mutation_cnt = Coding_SNV_cnt+Coding_indel_cnt,
         ) %>% 
  filter(!is.na(Homozygous_gene_count))

########## Function ##########
cumulative_onset_plot <- function(data, tumor, allele, dir, filename, annotation=NULL){
  
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
                        fontsize=6, xlab="Age", ylab="Cumulative onset", title=paste(tumor), 
                        font.title=15, legend.title=allele, legend=c(0.80,0.15), font.legend=c(15),
                        font.tickslab=c(15, "black"), font.x = c(15, "black"), font.y = c(15, "black"),
                        palette=c("#56B4E9","#E69F00"), cumevents=T) # count table
  if(!is.null(annotation)) surv_p$plot <- surv_p$plot + annotation
  
  pdf(paste0(dir,paste("ICGC",filename,sep="_"),".pdf"), width=12, height=10, onefile=F)
  print(surv_p)
  dev.off()
  
  tiff(paste0(dir,paste("ICGC",filename,sep="_"),".tiff"), width=12, height=10, units="in", res=100)
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

FR_merged_cumulative_onset_plot <- function(data, tumor, dist, filename){
  
  AFT_data <- data
  form <- formula(paste("Surv(age_at_diagnosis,event)~Homo_at_least_1_locus+gender"))
  form_LR <- formula(paste("Surv(age_at_diagnosis,event)~Homo_at_least_1_locus"))
  
  model <- summary(survreg(formula=form, data=AFT_data, dist=dist))$table
  
  n <- nrow(AFT_data)
  
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
  p <- ifelse(p<0.001," < 0.001",paste0("= ",sprintf("%.3f",p)))
  cumulative_onset_plot(data=data, tumor=tumor, allele="Homo_at_least_1_locus", 
                        dir="plot/", filename=filename,
                        annotation=annotate("text", x = age_index_max-interval/8, y = 0.25, label = 
                                              paste0("TR (95% CI) = ",time_ratio_with_CI,"\n p ",p), parse = F, size = 6))
}


###################### analysis ##########################
######################## total ###########################
t_type <- 'RECA-EU'
filtered_df <- filtered_df %>% filter(tumor_type == t_type & (tumor_stage=="T1N0M0"|tumor_stage=="T1NXM0"))

##########################################################
filtered_df %>% filter(VHL_biallelic_mut=='O')
table(filtered_df$Homo_at_least_1_locus)

dist <- dist_select(filtered_df, "Homo_at_least_1_locus")[[1]]
FR_merged_cumulative_onset_plot(data=filtered_df, tumor=t_type, dist=dist, filename="Total_KIRC_KM_plot")

y_start <- min(filtered_df$age_at_diagnosis)
y_end <- max(filtered_df$age_at_diagnosis)

###########################################################
#################### VHL biallelic mut ####################
VHL_filtered_df <- filtered_df %>% filter(VHL_biallelic_mut=='O')
table(VHL_filtered_df$Homo_at_least_1_locus)

dist <- dist_select(VHL_filtered_df, "Homo_at_least_1_locus")[[1]]
FR_merged_cumulative_onset_plot(data=VHL_filtered_df, tumor=t_type, dist=dist, filename="VHL_biallelic_loss_Matced_KIRC_KM_plot")

p <- VHL_filtered_df %>%
  ggplot(aes(x=Homo_at_least_1_locus,y=age_at_diagnosis,fill=Homo_at_least_1_locus)) +
  geom_boxplot(width=0.8,outlier.shape = NA) +
  geom_point(position =position_jitterdodge(jitter.width=0.2)) +
  stat_compare_means(size=4, label.x = 1.4, method = "wilcox.test",label = "p.format") +
  coord_cartesian(ylim = c(y_start, y_end)) +
  scale_fill_manual(values=c("#56B4E9","#E69F00")) +
  theme_classic() + theme(legend.position="bottom",text = element_text(size=10))
print(p)

VHL_filtered_df %>% group_by(Homo_at_least_1_locus) %>% summarise(median=median(age_at_diagnosis, na.rm=T))
print(p)

VHL_filtered_df %>% filter(VHL_biallelic_mut=='O') %>% 
  group_by(Homo_at_least_1_locus) %>% summarise(n=n())


###########################################################
#################### PBRM1 biallelic mut ####################
PBRM1_filtered_df <- filtered_df %>% filter(PBRM1_biallelic_mut=='O')
table(PBRM1_filtered_df$Homo_at_least_1_locus)

dist <- dist_select(PBRM1_filtered_df, "Homo_at_least_1_locus")[[1]]
FR_merged_cumulative_onset_plot(data=PBRM1_filtered_df, tumor=t_type, dist=dist, filename="PBRM1_biallelic_loss_Matced_KIRC_KM_plot")

p <- PBRM1_filtered_df %>%
  ggplot(aes(x=Homo_at_least_1_locus,y=age_at_diagnosis,fill=Homo_at_least_1_locus)) +
  geom_boxplot(width=0.8,outlier.shape = NA) +
  geom_point(position =position_jitterdodge(jitter.width=0.2)) +
  stat_compare_means(size=4, label.x = 1.4, method = "wilcox.test",label = "p.format") +
  coord_cartesian(ylim = c(y_start, y_end)) +
  scale_fill_manual(values=c("#56B4E9","#E69F00")) +
  theme_classic() + theme(legend.position="bottom",text = element_text(size=10))
print(p)

PBRM1_filtered_df %>% filter(PBRM1_biallelic_mut=='O') %>% 
  group_by(Homo_at_least_1_locus) %>% summarise(n=n())


###################################################
#################### Wild type ####################
WT_filtered_df <- filtered_df %>% filter(is.na(VHL_biallelic_mut)&is.na(PBRM1_biallelic_mut))
table(WT_filtered_df$Homo_at_least_1_locus)

dist <- dist_select(WT_filtered_df, "Homo_at_least_1_locus")[[1]]
FR_merged_cumulative_onset_plot(data=WT_filtered_df, tumor=t_type, dist=dist, filename="WT_KIRC_KM_plot")

p <- WT_filtered_df %>%
  ggplot(aes(x=Homo_at_least_1_locus,y=age_at_diagnosis,fill=Homo_at_least_1_locus)) +
  geom_boxplot(width=0.8,outlier.shape = NA) +
  geom_point(position =position_jitterdodge(jitter.width=0.2)) +
  stat_compare_means(size=4, label.x = 1.4, method = "wilcox.test",label = "p.format") +
  coord_cartesian(ylim = c(y_start, y_end)) +
  scale_fill_manual(values=c("#56B4E9","#E69F00")) +
  theme_classic() + theme(legend.position="bottom",text = element_text(size=10))
print(p)

WT_filtered_df %>% 
  group_by(Homo_at_least_1_locus) %>% summarise(n=n())

###############################################################
###############################################################
none_VHL_filtered_df <- filtered_df %>% filter(VHL_biallelic_mut!='O'|is.na(VHL_biallelic_mut))
table(none_VHL_filtered_df$Homo_at_least_1_locus)

dist <- dist_select(none_VHL_filtered_df, "Homo_at_least_1_locus")[[1]]
FR_merged_cumulative_onset_plot(data=none_VHL_filtered_df, tumor=t_type, dist=dist, filename="None_VHL_biallelic_loss_total_KIRC_KM_plot")

VHL_loss <- none_VHL_filtered_df %>%
  ggplot(aes(x=Homo_at_least_1_locus,y=mutation_cnt,fill=Homo_at_least_1_locus)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position =position_jitterdodge()) +
  stat_compare_means(method = "wilcox.test",label = "p.format") +
  # ylim(0, 200) +
  theme_classic()


#############################################################
#################### PBRM1 biallelic mut ####################
PBRM1_filtered_df <- filtered_df %>% filter(PBRM1_biallelic_mut=='O')
table(PBRM1_filtered_df$Homo_at_least_1_locus)

dist <- dist_select(PBRM1_filtered_df, "Homo_at_least_1_locus")[[1]]
FR_merged_cumulative_onset_plot(data=PBRM1_filtered_df, tumor=t_type, dist=dist, filename="PBRM1_biallelic_loss_Matced_KIRC_KM_plot")

none_PBRM1_filtered_df <- filtered_df %>% filter(PBRM1_biallelic_mut!='O'|is.na(PBRM1_biallelic_mut))
table(none_PBRM1_filtered_df$Homo_at_least_1_locus)

dist <- dist_select(none_PBRM1_filtered_df, "Homo_at_least_1_locus")[[1]]
FR_merged_cumulative_onset_plot(data=none_PBRM1_filtered_df, tumor=t_type, dist=dist, filename="None_PBRM1_biallelic_loss_total_KIRC_KM_plot")


PBRM1_loss <- none_VHL_filtered_df %>%
  ggplot(aes(x=Homo_at_least_1_locus,y=mutation_cnt,fill=Homo_at_least_1_locus)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position =position_jitterdodge()) +
  stat_compare_means(method = "wilcox.test",label = "p.format") +
  theme_classic()

###########################################################
#################### both biallelic mut ####################
VHL_PBRM1_filtered_df <- filtered_df %>% filter(VHL_biallelic_mut=='O'&PBRM1_biallelic_mut=='O')
AFT_model_summary(data=VHL_PBRM1_filtered_df, stage=1, tumor=t_type, 
                  var=c("Homo_at_least_1_locus"), ad_var=c("gender"))
cumulative_onset_plot(data=VHL_PBRM1_filtered_df, tumor=t_type, allele="Homo_at_least_1_locus")
