require(rms)
library(survival)
library(survminer)
library(dplyr)
library(tibble)
library(Hmisc)
library(timereg)
library(openxlsx)
library(ggplot2)
library(coin)
library(tidyr)
library(KMsurv)
library(metafor)

########## Data preprocessing ##########
df_FILE="TCGA_ann.tsv"
raw_df <- read.csv(df_FILE, header=T,sep='\t')

names(raw_df)[names(raw_df)=="age_at_initial_pathologic_diagnosis"] <- "age_at_diagnosis"
filtered_df <- raw_df %>% 
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
         Germline_variant = as.factor(ifelse(P_total_GV_count==0,0,1)),
         age_at_diagnosis = as.numeric(as.character(age_at_diagnosis)),
         ethnicity = as.factor(as.character(ethnicity)),
         race2 = as.factor(ifelse(race=="WHITE",1,0)),
         gender = as.factor(as.character(gender)),
         tumor_stage = as.factor(as.character(tumor_stage)),
         prior_malignancy = as.factor(as.character(prior_malignancy)),
         tumor = as.factor(as.character(project_id)),
         staging_classification = ifelse(grepl("ajcc_pathological_staging",staging_method),"Pathological",
                                         ifelse(grepl("NA",staging_method),NA,
                                                ifelse(grepl("stag",staging_method), "Clinical",NA
                                                ))),
         stage = ifelse(grepl("IV",tumor_stage),4,
                        ifelse(grepl("III",tumor_stage),3,
                               ifelse(grepl("I/II",tumor_stage),NA,
                                      ifelse(grepl("II",tumor_stage),2,
                                             ifelse(grepl("I",tumor_stage),1,
                                                    ifelse(grepl(0,tumor_stage),1,NA))))))) %>% 
  filter(stage=='1'&prior_malignancy=="no"&!is.na(age_at_diagnosis)&!is.na(Homozygous_gene_count)) %>% 
  filter(Germline_variant==0&(ILC_Infection_status=="X"|is.na(ILC_Infection_status)))

######### Survival Curve ##########
cumulative_onset_plot <- function(data, stage, tumor, allele, dir, annotation=NULL){

  if(grepl('total',tumor)){
    data <- data
  }else{
    data <- data[which(data$stage==stage&data$tumor==tumor),]
  }
  
  
  age_index_min <-  floor(min(as.numeric(as.character(data$age_at_diagnosis)),na.rm=T))
  age_index_max <- ceiling(max(as.numeric(as.character(data$age_at_diagnosis)),na.rm=T))
  interval <- age_index_max-age_index_min
  
  form <- formula(paste0("Surv(age_at_diagnosis,event) ~ ",allele))
  fit <- do.call(survfit, args=list(formula=form, data=data))
  
  if(dir.exists(dir)==F){
    dir.create(dir)
  }
  
  pdf(paste0(dir,paste("Stage",stage,tumor,allele,"type.pdf",sep="_")), width=12, height=10, onefile=F)
  p <- ggsurvplot(fit, conf.int=FALSE, fun="event", xlim=c(age_index_min,age_index_max+1),
                    break.time.by=10, pval=F, size = 1.2, fontsize=6,
                    xlab="Age", ylab="Cumulative onset", title=paste(tumor,"stage",stage),
                    legend.title=allele, font.title=15, legend=c(0.80,0.15), font.legend=c(15),
                    font.tickslab=c(15, "black"), font.x = c(15, "black"), font.y = c(15, "black"),
                    palette=c("#56B4E9","#E69F00"), cumevents=T) # count table
  if(!is.null(annotation)) p$plot <- p$plot + annotation
    print(p)
  dev.off()
}

expected_time_weibull <- function(beta, log_scale){
  return(exp(beta)*gamma(1+exp(log_scale)))}

freq <- filtered_df %>% group_by(tumor,Homo_at_least_1_locus) %>% filter(stage == 1) %>% 
  summarise(n=n()) 
n_table <- table(filtered_df$tumor)[table(filtered_df$tumor)!= 0]

Homo_cnt_over10 <- filtered_df %>% group_by(tumor, Homo_at_least_1_locus) %>% 
  summarise(n=n()) %>% filter(Homo_at_least_1_locus==1&n>=10)

Homo_cnt_less10 <- filtered_df %>% group_by(tumor, Homo_at_least_1_locus) %>% 
  summarise(n=n()) %>% filter(Homo_at_least_1_locus==1&n<10)
Homo_cnt_less10_type <- c('BLCA', unique(as.character(Homo_cnt_less10$tumor)))

nrow(filtered_df %>% filter(tumor %in% Homo_cnt_less10_type))

tumor_types <- sort(as.character(unique(Homo_cnt_over10$tumor)))
n <- as.numeric(unname(n_table[tumor_types]))

for(i in 1:length(tumor_types)){
  tumor <- tumor_types[i]
  AFT_data <- filtered_df[which(filtered_df$stage==1&filtered_df$tumor==tumor_types[i]),]
  form <- formula(paste("Surv(age_at_diagnosis,event)~Homo_at_least_1_locus+gender+race2"))
  
  m1 <- AIC(survreg(formula=form, data=AFT_data, dist="weibull"))
  m2 <- AIC(survreg(formula=form, data=AFT_data, dist="exponential"))
  m3 <- AIC(survreg(formula=form, data=AFT_data, dist="lognormal"))
  m4 <- AIC(survreg(formula=form, data=AFT_data, dist="loglogistic"))
  
  tmp <- data.frame(tumor, weibull=m1, exponential=m2, lognormal=m3, loglogistic=m4)
  if(i==1){
    res <- tmp
  }else{
    res <- rbind(res, tmp)
  }
}
which_dist <- apply(res[,2:5],1,which.min)

dist_list <- c("weibull","exponential","lognormal","loglogistic")


for(i in 1:length(tumor_types)){

  order <- 'Each_tumor'
  AFT_data <- filtered_df[which(filtered_df$stage==1&filtered_df$tumor==tumor_types[i]),]
  form <- formula(paste("Surv(age_at_diagnosis,event)~Homo_at_least_1_locus+gender+race2"))
  form_LR <- formula(paste("Surv(age_at_diagnosis,event)~Homo_at_least_1_locus"))
  
  model <- summary(survreg(formula=form, data=AFT_data, dist=dist_list[which_dist[i]]))$table

  homo_n <- nrow(AFT_data %>% filter(Homo_at_least_1_locus=='1'))
  hetero_n <- nrow(AFT_data %>% filter(Homo_at_least_1_locus=='0'))
  
  tumor <- tumor_types[i]
  time_ratio <- exp(model[2,1])
  lower <- exp(model[2,1]-qnorm(0.975)*model[2,2])
  upper <- exp(model[2,1]+qnorm(0.975)*model[2,2])
  hetero_age <- round(expected_time_weibull(model[1,1], model[length(rownames(model)),1]),2)
  homo_age <- round(expected_time_weibull(model[1,1]+model[2,1], model[length(rownames(model)),1]),2)
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

  age_index_min <<-  floor(min(as.numeric(as.character(AFT_data$age_at_diagnosis)),na.rm=T))
  age_index_max <<- ceiling(max(as.numeric(as.character(AFT_data$age_at_diagnosis)),na.rm=T))
  interval <<- age_index_max-age_index_min

  cumulative_onset_plot(data=filtered_df, stage=1, tumor=tumor_types[i], allele="Homo_at_least_1_locus",
                         dir="plot/", #annotation=NA)
                         annotation=annotate("text", x = age_index_max-interval/8, y = 0.25, label =
                                               paste0("TR (95% CI) = ",time_ratio_with_CI,"\n p ",p), parse = F, size = 5))

  tmp <- data.frame(tumor, homo_n, hetero_n, time_ratio_with_CI, p_value, hetero_age, homo_age, diff_age, order)
  if(i==1){
    uni_tumor_res <- tmp
  }else{
    uni_tumor_res <- rbind(uni_tumor_res, tmp)
  }
}

uni_tumor_res %>% arrange(desc(time_ratio_with_CI),diff_age)

uni_tumor_res <- uni_tumor_res %>% arrange(desc(time_ratio_with_CI))
write.xlsx(uni_tumor_res,"Univariable_homo_at_1_least_tumor_AFT_200929.xlsx")

######### Accelerated failure time model (AFT model) for multiple tumors #########
AFT_model <- function(data, stage, tumor, var, ad_var, tumor_level, dist){

  AFT_data <- data[which(data$stage==stage&data$tumor%in%tumor),]
  AFT_data$tumor <- factor(as.character(AFT_data$tumor), levels=tumor_level)
  
  form <- formula(paste("Surv(age_at_diagnosis,event)~",paste(c("tumor",var,ad_var),collapse="+")))
  model <- summary(survreg(formula=form, data=AFT_data, dist=dist))$table

  time_ratio <- exp(model[,1])
  lower <- exp(model[,1]-qnorm(0.975)*model[,2])
  upper <- exp(model[,1]+qnorm(0.975)*model[,2])
  
  p_value <- model[,4]
  
  time_ratio <- sprintf("%.3f", time_ratio)
  lower <- sprintf("%.3f", lower)
  upper <- sprintf("%.3f", upper)
  
  time_ratio_with_CI <- paste0(time_ratio," (",lower,",",upper,")")
  p_value <- ifelse(p_value<0.001," < 0.001",sprintf("%.3f",p_value))
  
  adjust <- data.frame(time_ratio_with_CI, p_value)

  return(list(adjust, BIC(survreg(formula=form, data=AFT_data, dist=dist)),AFT_data))
}

########### age difference ############
Age_diff_CI_cal <- function(data, input_set){
  
  AFT_data <- data[which(data$stage==1),]
  AFT_data$tumor <- as.factor(as.character(AFT_data$tumor))
  
  form <- formula(paste("Surv(age_at_diagnosis,event)~Homo_at_least_1_locus+gender+race2+tumor"))
  form_LR <- formula(paste("Surv(age_at_diagnosis,event)~Homo_at_least_1_locus"))
  
  model <- summary(survreg(formula=form, data=AFT_data, dist=dist_list[which_dist[i]]))$table
  
  homo_n <- nrow(AFT_data %>% filter(Homo_at_least_1_locus=='1'))
  hetero_n <- nrow(AFT_data %>% filter(Homo_at_least_1_locus=='0'))
  
  for(i in 2:4){
    
    tumor <- rownames(model)[i]
    time_ratio <- exp(model[i,1])
    lower <- exp(model[i,1]-qnorm(0.975)*model[i,2])
    upper <- exp(model[i,1]+qnorm(0.975)*model[i,2])
    hetero_age <- round(expected_time_weibull(model[1,1], model[length(rownames(model)),1]),2)
    homo_age <- round(expected_time_weibull(model[1,1]+model[i,1], model[length(rownames(model)),1]),2)
    diff_age <- hetero_age-homo_age
    p_value <- model[i,4]
    
    time_ratio <- sprintf("%.3f", time_ratio)
    lower <- sprintf("%.3f", lower)
    upper <- sprintf("%.3f", upper)
    
    time_ratio_with_CI <- paste0(time_ratio," (",lower,",",upper,")")
    p_value <- ifelse(p_value<0.001," < 0.001",sprintf("%.3f",p_value))
    
    tmp <- data.frame(tumor, homo_n, hetero_n, time_ratio_with_CI, p_value, hetero_age, homo_age, diff_age, order)
    
    if(i==2){
      LR <- survdiff(form_LR, AFT_data)
      p <- pchisq(LR$chisq,1,lower.tail=F)
      p <- ifelse(p<0.001," < 0.001",paste0("= ",sprintf("%.3f",p)))
      
      age_index_min <<-  floor(min(as.numeric(as.character(AFT_data$age_at_diagnosis)),na.rm=T))
      age_index_max <<- ceiling(max(as.numeric(as.character(AFT_data$age_at_diagnosis)),na.rm=T))
      interval <<- age_index_max-age_index_min
      
      cumulative_onset_plot(data=AFT_data, stage=1, tumor=input_set, allele="Homo_at_least_1_locus",
                            dir="plot/", #annotation=NA)
                            annotation=annotate("text", x = age_index_max-interval/8, y = 0.25, label =
                                                  paste0("FR (95% CI) = ",time_ratio_with_CI,"\n p ",p), parse = F, size = 5))
      multi_res <- tmp
    }else{
      multi_res <- rbind(multi_res, tmp)
    }
  }
  return(multi_res)
}

###############################################
positive_b_t_types <- as.character(unique(uni_tumor_res$tumor))
representative_res <- AFT_model(data=filtered_df, stage=1, 
                                 tumor=positive_b_t_types, 
                                 tumor_level=positive_b_t_types, # 첫번째 tumor가 reference tumor
                                 var=c("Homo_at_least_1_locus"), 
                                 ad_var=c("gender","race2"), dist="loglogistic")



Pan_cancer_df <- filtered_df %>% filter(tumor %in% positive_b_t_types)
Pan_cancer_res <- Age_diff_CI_cal(Pan_cancer_df, 'total')
Pan_cancer_res[1,1] <- "Pan_cancer_HLA"
Pan_cancer_res[1,which(colnames(Pan_cancer_res)=="order")] <- "Pan_cancer_HLA"
Pan_cancer_res[2:3,which(colnames(Pan_cancer_res)=="order")] <- "Pan_cancer_adjusting_variable"

Pan_cancer_res[2,2] <- nrow(Pan_cancer_df %>% filter(gender=='MALE'))
Pan_cancer_res[2,3] <- nrow(Pan_cancer_df %>% filter(gender=='FEMALE'))
Pan_cancer_res[3,2] <- nrow(Pan_cancer_df %>% filter(race2=='1'))
Pan_cancer_res[3,3] <- nrow(Pan_cancer_df %>% filter(race2=='0'))


########### Sample count ################
Sample_count_table <- t(as.data.frame(representative_res[[3]] %>%
  group_by(tumor, Homo_at_least_1_locus) %>%
  summarise(n=n()) %>% 
  arrange(desc(n)) %>% 
  pivot_wider(names_from = tumor, values_from = n)))

Sample_count_table[1,which(Sample_count_table[1,]=="0")] <- 'Heterozygous'
Sample_count_table[1,which(Sample_count_table[1,]=="1")] <- 'Homozygous'

colnames(Sample_count_table) <- Sample_count_table[1,]
Sample_count_table <- as.data.frame(Sample_count_table[-1,])

Sample_count_table <- Sample_count_table %>% 
  rownames_to_column('Tumor_type') %>% 
  mutate(
  Heterozygous = as.numeric(Heterozygous),
  Homozygous = as.numeric(Homozygous),
  Total = as.numeric(Heterozygous)+as.numeric(Homozygous),
  )

write.xlsx(Sample_count_table,"Final_analysis_dataset_count_table.xlsx", rowNames = FALSE)
write.xlsx(representative_res[[1]],"All_Representative_result_Adjusted_multivariate_homo_at_least_1.xlsx", rowNames = TRUE)

################ without ccRCC ##############################
negative_b_t_types <- c('KIRC')
positive_b_t_types <- as.character(unique(uni_tumor_res$tumor)[!unique(uni_tumor_res$tumor) %in% negative_b_t_types])

representative_res <- AFT_model(data=filtered_df, stage=1, 
           tumor=positive_b_t_types, 
           tumor_level=positive_b_t_types, # 첫번째 tumor가 reference tumor
           var=c("Homo_at_least_1_locus"), 
           ad_var=c("gender","race2"), dist="loglogistic")
representative_res[[1]]

Pan_cancer_wo_ccrcc_types_df <- filtered_df %>% filter(tumor %in% positive_b_t_types)
Pan_cancer_wo_ccrcc_res <- Age_diff_CI_cal(data=Pan_cancer_wo_ccrcc_types_df, input_set='total_wo_ccRCC')[1,]
Pan_cancer_wo_ccrcc_res[1,1] <- "Pan_cancer_wo_ccRCC_HLA" 
Pan_cancer_wo_ccrcc_res[1,which(colnames(Pan_cancer_res)=="order")] <- "Pan_cancer_wo_ccRCC_HLA" 

write.xlsx(representative_res[[1]],"Representative_result_Adjusted_multivariate_homo_at_least_1.xlsx", rowNames = TRUE)

merged_res <- rbind(Pan_cancer_res[1,],Pan_cancer_wo_ccrcc_res,uni_tumor_res,Pan_cancer_res[2:3,])
rownames(merged_res) <- 1:length(rownames(merged_res))

########## forest plot ########## 
merged_res$FR <- as.numeric(substr(merged_res$time_ratio_with_CI,1,5))
merged_res$lower <- as.numeric(substr(merged_res$time_ratio_with_CI,8,12))
merged_res$upper <- as.numeric(substr(merged_res$time_ratio_with_CI,14,18))
merged_res$homo_age <- as.character(sprintf("%.2f",round(as.numeric(merged_res$homo_age),2)))
merged_res$hetero_age <- as.character(sprintf("%.2f",round(as.numeric(merged_res$hetero_age),2)))
merged_res$diff_age <- as.character(round(as.numeric(merged_res$diff_age),2))
merged_res$order <- as.factor(merged_res$order)
merged_res <- rev(merged_res %>% arrange(merged_res$order))

pdf("AFT_forest_plot.pdf", width=16, height=10, onefile=F)
  forest(x=merged_res$FR, ci.lb=merged_res$lower, ci.ub=merged_res$upper, slab=merged_res$tumor, 
         # rows=c(1:4,6:16),
         ilab=cbind(merged_res$homo_n,merged_res$hetero_n,merged_res$homo_age, merged_res$hetero_age,merged_res$diff_age),
         xlim=c(0,2), ylim=c(1, 17),
         ilab.xpos=c(0.5,0.7,1.3,1.45,1.6), xlab='',
         digits=c(3,3), cex=1.5, efac=1.5, refline=1, top=1.5)
  text(c(0.1,0.5,0.7,1.3,1.45,1.6,1.86), 17.1, c('Tumor types','Homo (N)','Hetero (N)',
                                                 'Homo age','Hetero age','Diffage',
                                                 'Failure Rate [95% CI]'), cex=1.3, font=2)
dev.off()


data <- filtered_df

tumor_types <- uni_tumor_res$tumor
data <- data %>% filter(tumor %in% tumor_types)
data$tumor <- factor(as.character(data$tumor), levels=rev(tumor_types))

temp <- data %>% filter(!is.na(get("Homo_at_least_1_locus"))) %>% group_by(tumor, get("Homo_at_least_1_locus")) %>%
  summarise(median=median(age_at_diagnosis, na.rm=T),
            Q1=quantile(age_at_diagnosis, c(0.25), na.rm=T),
            Q3=quantile(age_at_diagnosis, c(0.75), na.rm=T),
            n=n())

temp$med_IQR <- paste0(sprintf("%.1f",temp$median),"(",sprintf("%.1f",temp$Q1),"-",sprintf("%.1f",temp$Q3),")")

for(i in 1:length(unique(temp$`get("Homo_at_least_1_locus")`))){
  if(i==1){
    res <- temp %>% filter(`get("Homo_at_least_1_locus")`==unique(temp$`get("Homo_at_least_1_locus")`)[i]) %>% select(tumor, n, med_IQR)
    names(res)[-1] <- paste0(c("n","med_IQR"),"_",unique(temp$`get(type)`)[i])
  }else{
    tmp <- temp %>% filter(`get("Homo_at_least_1_locus")`==unique(temp$`get("Homo_at_least_1_locus")`)[i]) %>% select(tumor, n, med_IQR)
    names(tmp)[-1] <- paste0(c("n","med_IQR"),"_",unique(temp$`get("Homo_at_least_1_locus")`)[i])
    res <- merge(res, tmp, by="tumor", all=T)
  }
}

res$p_value <- NA
for(i in 1:nrow(res)){
  p_value <- tryCatch(
    pvalue(kruskal_test(age_at_diagnosis~get("Homo_at_least_1_locus"), data[data$tumor==res$tumor[i],])),
    error=function(e){NA})
  p_value <- ifelse(p_value<0.001,"<0.001",sprintf("%.3f",p_value))
  res$p_value[i] <- p_value
}

write.xlsx(res, paste0("Homo_at_least_1_locus.xlsx"), row.names=F)

######### Box plot #########
p <- Pan_cancer_df %>% filter(!is.na(get("Homo_at_least_1_locus"))) %>%
  ggplot(aes(Homo_at_least_1_locus, age_at_diagnosis, fill=get("Homo_at_least_1_locus"))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#56B4E9","#E69F00")) +
  scale_y_continuous(breaks=seq(25, 100, 10), limit = c(25, 90))+
  coord_flip() +
  labs(x = "Tumor type", y = "Onset age") +
  guides(fill=guide_legend(title="Homo_at_least_1_locus"))+
  theme_classic() + theme(legend.position="bottom",text = element_text(size=30))

pdf("Pancancer_boxplot.pdf", width=16, height=10)
  print(p)
dev.off()

p <- Pan_cancer_wo_ccrcc_types_df %>% filter(!is.na(get("Homo_at_least_1_locus"))) %>%
  ggplot(aes(Homo_at_least_1_locus, age_at_diagnosis, fill=get("Homo_at_least_1_locus"))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#56B4E9","#E69F00")) +
  scale_y_continuous(breaks=seq(25, 100, 10), limit = c(25, 90))+
  coord_flip() +
  labs(x = "Tumor type", y = "Onset age") +
  guides(fill=guide_legend(title="Homo_at_least_1_locus"))+
  theme_classic() + theme(legend.position="bottom",text = element_text(size=30))

pdf("Pan_wo_ccRCC_boxplot.pdf", width=16, height=10)
  print(p)
dev.off()

p <- data %>% filter(!is.na(get("Homo_at_least_1_locus"))) %>%
  ggplot(aes(tumor, age_at_diagnosis, fill=get("Homo_at_least_1_locus"))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#56B4E9","#E69F00")) +
  scale_y_continuous(breaks=seq(25, 100, 10), limit = c(25, 90))+
  coord_flip() +
  labs(x = "Tumor type", y = "Onset age") +
  guides(fill=guide_legend(title="Homo_at_least_1_locus"))+
  theme_classic() + theme(legend.position="bottom",text = element_text(size=30))

pdf(paste0("Tumor_types_boxplot.pdf"), width=16, height=10)
  print(p)
dev.off()