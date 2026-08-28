################################################################################
################################################################################
###################  Study 1 Code ##############################################
################################################################################



## Libraries -------------------------------------------------------------------
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cli)
library(ztable)
library(forcats)
library(poLCA)
library(purrr)

## Upload Data from ABCD 5.1 ---------------------------------------------------
data <- read.csv("data/5.1/mh_p_cbcl.csv", header=TRUE)
subset_data <- subset(data, eventname == 'baseline_year_1_arm_1' & cbcl_q02_p < 1 & cbcl_q105_p < 1)
subset_data <- subset_data[, c('src_subject_id','cbcl_scr_syn_anxdep_t','cbcl_scr_syn_withdep_t','cbcl_scr_syn_somatic_t','cbcl_scr_syn_social_t','cbcl_scr_syn_thought_t','cbcl_scr_syn_attention_t','cbcl_scr_syn_rulebreak_t','cbcl_scr_syn_aggressive_t')]
subset_data <- na.omit(subset_data)

# Modify CBCL data to binary scores based on T Score cutoff
subset_data$cbcl_anxdep_binary <- ifelse(subset_data$cbcl_scr_syn_anxdep_t >= 65, 2, 1)
subset_data$cbcl_withdep_binary <- ifelse(subset_data$cbcl_scr_syn_withdep_t >= 65, 2, 1)
subset_data$cbcl_somatic_binary <- ifelse(subset_data$cbcl_scr_syn_somatic_t >= 65, 2, 1)
subset_data$cbcl_social_binary <- ifelse(subset_data$cbcl_scr_syn_social_t >= 65, 2, 1)
subset_data$cbcl_thought_binary <- ifelse(subset_data$cbcl_scr_syn_thought_t >= 65, 2, 1)
subset_data$cbcl_attention_binary <- ifelse(subset_data$cbcl_scr_syn_attention_t >= 65, 2, 1)
subset_data$cbcl_rulebreak_binary <- ifelse(subset_data$cbcl_scr_syn_rulebreak_t >= 65, 2, 1)
subset_data$cbcl_aggressive_binary <- ifelse(subset_data$cbcl_scr_syn_aggressive_t >= 65, 2, 1)



## LCA -------------------------------------------------------------------------
f <- with(subset_data, cbind(cbcl_anxdep_binary,cbcl_withdep_binary,cbcl_somatic_binary,cbcl_social_binary,cbcl_thought_binary,cbcl_attention_binary,cbcl_rulebreak_binary,cbcl_aggressive_binary)~1)

set.seed(01012)
lc2<-poLCA(f, data=subset_data, nclass=2, na.rm = FALSE, nrep=30, maxiter=100000) 
lc3<-poLCA(f, data=subset_data, nclass=3, na.rm = FALSE, nrep=30, maxiter=100000)
lc4<-poLCA(f, data=subset_data, nclass=4, na.rm = FALSE, nrep=30, maxiter=100000)
lc5<-poLCA(f, data=subset_data, nclass=5, na.rm = FALSE, nrep=30, maxiter=100000) 
lc6<-poLCA(f, data=subset_data, nclass=6, na.rm = FALSE, nrep=30, maxiter=100000)


results <- data.frame(Model=c("Model 1"),                                                 
                      AIC = lc2$aic,                                                    
                      log_likelihood = lc2$llik,                                        
                      df = lc2$resid.df,                                                
                      BIC = lc2$bic,                                                   
                      ABIC = (-2*lc2$llik) + ((log((lc2$N + 2/24)) * lc2$npar)),        
                      CAIC = (-2*lc2$llik) + lc2$npar * (1 + log(lc2$N)),               
                      likelihood_ratio = lc2$Gsq)  

results$Model <- as.integer(results$Model)

results[1,1]<-c("2 class")
results[2,1]<-c("3 class")
results[3,1]<-c("4 class")
results[4,1]<-c("5 class")
results[5,1]<-c("6 class")

#Add AIC values in column 2

results[2,2]<-lc3$aic
results[3,2]<-lc4$aic
results[4,2]<-lc5$aic
results[5,2]<-lc6$aic

#Add log-likelihood values to rows in column 3 

results[2,3]<-lc3$llik
results[3,3]<-lc4$llik
results[4,3]<-lc5$llik
results[5,3]<-lc6$llik

#Add df values to rows in column 4

results[2,4]<-lc3$resid.df
results[3,4]<-lc4$resid.df
results[4,4]<-lc5$resid.df
results[5,4]<-lc6$resid.df

#Add BIC values to rows in column 5 

results[2,5]<-lc3$bic
results[3,5]<-lc4$bic
results[4,5]<-lc5$bic
results[5,5]<-lc6$bic

#Add ABIC values to rows in column 6 - calculation for ABIC included 

results[2,6]<-(-2*lc3$llik) + ((log((lc3$N + 2)/24)) * lc3$npar)
results[3,6]<-(-2*lc4$llik) + ((log((lc4$N + 2)/24)) * lc4$npar)
results[4,6]<-(-2*lc5$llik) + ((log((lc5$N + 2)/24)) * lc5$npar)
results[5,6]<-(-2*lc6$llik) + ((log((lc6$N + 2)/24)) * lc6$npar)

#Add CAIC values to rows in column 7 - calculation for CAIC included 

results[2,7]<- (-2*lc3$llik) + lc3$npar * (1 + log(lc3$N))
results[3,7]<- (-2*lc4$llik) + lc4$npar * (1 + log(lc4$N))
results[4,7]<- (-2*lc5$llik) + lc5$npar * (1 + log(lc5$N))
results[5,7]<- (-2*lc6$llik) + lc6$npar * (1 + log(lc6$N))

#Add likelihood ratio values to rows in column 8 

results[2,8]<-lc3$Gsq
results[3,8]<-lc4$Gsq
results[4,8]<-lc5$Gsq
results[5,8]<-lc6$Gsq

#Create entropy as a function 

entropy<-function (p) sum(-p*log(p))

#Add column 9 to dataframe for entropy 

results[1,9] <- c("-")

#Calculate entropy for Model 1
error_prior<-entropy(lc2$P)
error_post<-mean(apply(lc2$posterior,1, entropy),na.rm = TRUE)
results[1,9]<-round(((error_prior-error_post) / error_prior),3)


#Calculate entropy for Model 2 

error_prior<-entropy(lc3$P)
error_post<-mean(apply(lc3$posterior,1, entropy),na.rm = TRUE)
results[2,9]<-round(((error_prior-error_post) / error_prior),3)

#Calculate entropy for Model 3

error_prior<-entropy(lc4$P) 
error_post<-mean(apply(lc4$posterior,1, entropy),na.rm = TRUE)
results[3,9]<-round(((error_prior-error_post) / error_prior),3)

#Calculate entropy for Model 4 

error_prior<-entropy(lc5$P) 
error_post<-mean(apply(lc5$posterior,1, entropy),na.rm = TRUE)
results[4,9]<-round(((error_prior-error_post) / error_prior),3)

#Calculate entropy for Model 5 

error_prior<-entropy(lc6$P)
error_post<-mean(apply(lc6$posterior,1, entropy),na.rm = TRUE)
results[5,9]<-round(((error_prior-error_post) / error_prior),3)

#Combine results to a dataframe 

colnames(results)<-c("Model","AIC", "log-likelihood","resid. df","BIC","aBIC","cAIC","likelihood-ratio","Entropy")

#Create table for Model_results 


ztable::ztable(results)

#Assign Model column as a factor under forcats 

results$Model <- as_factor(results$Model) 

#Convert to long format

results2 <- tidyr::gather(results,Kriterium,Guete,5:8)
results2 <- results2 %>%
  mutate(Kriterium = ifelse(Kriterium == 'likelihood-ratio','like.\nratio',Kriterium))

#Plot aBIC, BIC, cAIC and likelihood ratio data for each class 

fit.plot <- ggplot(results2) + geom_point(aes(x=Model, y=Guete), size=3) + geom_line(aes(Model, Guete, group = 1)) + theme_bw() + labs(x = "", y="", title = "") + facet_grid(Kriterium ~. ,scales = "free") + theme_bw(base_size = 16, base_family = "") + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour="grey", linewidth=0.5), legend.title = element_text(size = 16, face = 'bold'), axis.text = element_text(size = 16), axis.title = element_text(size = 16), legend.text = element_text(size=16), axis.line = element_line(colour = "black")) 
fit.plot


## Average posterior probabilities
get_app_polca <- function(model, model_name) {
  
  # Posterior probability matrix
  probs <- as.matrix(model$posterior)
  
  # Assigned classes
  assigned_class <- model$predclass
  
  # Remove missing assignments if present
  valid_rows <- !is.na(assigned_class)
  
  probs <- probs[valid_rows, , drop = FALSE]
  assigned_class <- assigned_class[valid_rows]
  
  # Extract assigned-class probability
  assigned_prob <- probs[
    cbind(seq_len(nrow(probs)), assigned_class)
  ]
  
  # Summarise APP by class
  out <- data.frame(
    class = assigned_class,
    assigned_prob = assigned_prob
  ) %>%
    group_by(class) %>%
    summarise(
      APP = mean(assigned_prob, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(model = model_name)
  
  return(out)
}

models <- list(
  lc2 = lc2,
  lc3 = lc3,
  lc4 = lc4,
  lc5 = lc5,
  lc6 = lc6
)

pprob_results <- imap_dfr(models, get_app_polca)





#Extract starting values from best model (4 class)

probs.start<-lc4$probs.start

#Re-run the best model (class 4) with graphs=true

set.seed(123)

lc<-poLCA(f, subset_data, nclass=4, probs.start=probs.start,graphs=TRUE, na.rm=TRUE, nrep=30, maxiter=10000)

#Assign values to poLCA posterior states

poLCA.posterior.states <- data.frame(lc$predclass,lc$posterior)

#Add state numbers to columns for posterior states

names(poLCA.posterior.states) <- c('state','s1','s2', 's3', 's4')

#Assign state column as factor/categorical

poLCA.posterior.states$state <- as.factor(poLCA.posterior.states$state)

#Change order of factors for purpose of the graph
#1 is universal, 2 is low, 3 is internal, 4 is external

poLCA.posterior.states$state <- factor(poLCA.posterior.states$state, levels = c("2","3",'4',"1"))

#Calculate average posterior probabilities for each class 

poLCA.posterior.states.s2 <- poLCA.posterior.states[poLCA.posterior.states$state == '2',]

round(colMeans(lc$posterior)*100,2)

#Plot data for posterior states
colnames(subset_data)[colnames(subset_data) %in% c('cbcl_anxdep_binary','cbcl_withdep_binary','cbcl_somatic_binary','cbcl_social_binary','cbcl_thought_binary','cbcl_attention_binary','cbcl_rulebreak_binary','cbcl_aggressive_binary')] <-
  c('Anxious','Withdrawn','Somatic','Social','Thought','Attention','Delinquent','Aggressive')
plot.data <- cbind(subset_data, poLCA.posterior.states$state, poLCA.posterior.states$s1, poLCA.posterior.states$s2, poLCA.posterior.states$s3,poLCA.posterior.states$s4) %>% gather(key="measure", value="value", 'Anxious','Withdrawn','Somatic','Social','Thought','Attention','Delinquent','Aggressive')
plot.data <- plot.data[, c('src_subject_id','poLCA.posterior.states$state','poLCA.posterior.states$s1','poLCA.posterior.states$s2','poLCA.posterior.states$s3','poLCA.posterior.states$s4','measure','value')]
#Change colnames

colnames(plot.data) <- c("Key", "class", "s1", "s2", "s3","s4", "measure", "value")

#Assign individual to class based on maximum probability 

max_lca_probs <- cbind(plot.data, poLCA.posterior.states$state, poLCA.posterior.states$s1, poLCA.posterior.states$s2, poLCA.posterior.states$s3,poLCA.posterior.states$s4) %>% 
  mutate(id=1:nrow(.)) %>% 
  gather("state_prob", "value", starts_with("S", ignore.case=FALSE)) %>% 
  group_by(id) %>% 
  summarize(max_prob=max(value, na.rm = TRUE))

#Probability of cases assigned to classes 

summary(max_lca_probs)

## LCA plots -------------------------------------------------------------------

#Change 1 scores to 0 and 2 scores to 1 for probability counts under 1

plot.data[,8][plot.data[,8] ==1] <- 0
plot.data[,8][plot.data[,8] ==2] <- 1

#Create summary plot with class, measure and calculation of mean probability of class responding 'yes' to measure

summary.plot.data <- plot.data %>% group_by(class, measure) %>% summarize(probability=mean(value))
desired_order <- c('Anxious','Withdrawn','Somatic','Social','Thought','Attention','Delinquent','Aggressive')  # Specify the desired order

# Convert 'measure' to a factor with the desired order of levels
summary.plot.data$measure <- factor(summary.plot.data$measure, levels = desired_order)


#Check frequency scores using 'count' - how many participants said 'yes' to a measure in each class

summary.plot.data2 <- plot.data %>% group_by(class, measure) %>% count(value)

#Create plot

ggplot(summary.plot.data, aes(y=probability, x=measure, group=class, color=class)) + geom_point() + geom_line(size = 1) +  ggtitle("LCA plot") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Change graph 

ggplot(summary.plot.data, aes(y=probability, x=measure, group=class, color=class)) + geom_point() + geom_line(linewidth = 0.6) +  ggtitle("Figure 1") + theme(axis.text.x = element_text(size = 9, angle = 90, hjust = 1)) + scale_color_manual(labels = c('Low Symptom (86.20%)', 'Internalizing (8.86%)', 'Externalizing (2.79%)','Universal (2.15%)'), values = c("#212237","#ec5552", "#60dc60","#0000FF")) + scale_y_continuous(limits = c(0, 1))

#Change to APA format

ggplot(summary.plot.data, aes(x=measure, y=probability, group=class)) + geom_point(aes(shape=class), size = 3) + geom_line() + scale_shape_manual(labels=c('Low Symptom', 'Internalizing', 'Externalizing','Universal'), values = c(15, 16, 17, 18)) + scale_x_discrete(limits = c('Anxious','Withdrawn','Somatic','Social','Thought','Attention','Delinquent','Aggressive')) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + labs(x = "CBCL Subscale", y = "Probability", title = "Figure 1") + theme_classic() + theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 12), axis.text = element_text(size = 10), plot.margin = unit(c(1,1,0.5,1), "cm"), axis.text.x = element_text(angle = 45, hjust = 1))

## Build final data set for SPSS analysis --------------------------------------
#Plot participant key and class assignment
final.data <- plot.data[, c('Key','class')]
final.data <- final.data[!duplicated(final.data$Key),]

#Open demographic data
data.demo <- read.csv('data/5.1/abcd_p_demo.csv', header=TRUE)
data.demo <- subset(data.demo, eventname == 'baseline_year_1_arm_1')
colnames(data.demo)[colnames(data.demo) %in% c('demo_race_a_p___10', 'demo_race_a_p___11', 'demo_race_a_p___12','demo_race_a_p___13',
                                               'demo_race_a_p___14', 'demo_race_a_p___15', 'demo_race_a_p___16', 'demo_race_a_p___17',
                                               'demo_race_a_p___18', 'demo_race_a_p___19', 'demo_race_a_p___20', 'demo_race_a_p___21',
                                               'demo_race_a_p___22', 'demo_race_a_p___23', 'demo_race_a_p___24', 'demo_race_a_p___25')] <-
  c('White', 'Black', 'Native', 'Alaskan','Hawaiian','Guamanian', 'Samoan', 'Pacific Islander', 'Indian', 'Chinese', 'Filipino', 'Japanese', 'Korean', 'Vietnamese',
    'Other Asian', 'Other Race')
colnames(data.demo)[colnames(data.demo)=='demo_ethn_v2'] <- 'Latino'
colnames(data.demo)[colnames(data.demo)=='src_subject_id'] <- 'Key'
data.demo <- data.demo %>%
  mutate('demo_sex_v2' = ifelse(demo_sex_v2 == 1, 0, 1), 'Latino' = ifelse(Latino == 1, 1, 0))
colnames(data.demo)[colnames(data.demo)=='demo_sex_v2'] <- 'Female'
data.demo <- data.demo[, c('Key','Female','White', 'Black', 'Native', 'Alaskan','Hawaiian','Guamanian', 'Samoan', 'Pacific Islander', 'Indian', 'Chinese', 'Filipino', 'Japanese', 'Korean', 'Vietnamese',
                           'Other Asian', 'Other Race', 'Latino')]
data.demo <- data.demo %>%
  mutate(Latino = ifelse(Latino == 1, 1, 0),
         White_NL = ifelse(Latino == 1, 0, White),
         Black_NL = ifelse(Latino == 1, 0, Black),
         Native_NL = ifelse(Latino == 1, 0, Native+Alaskan),
         PacificIslander_NL = ifelse(Latino == 1, 0, `Pacific Islander`+Hawaiian+Guamanian+Samoan),
         Asian_NL = ifelse(Latino == 1, 0, Indian+Chinese+Filipino+Japanese+Korean+Vietnamese+`Other Asian`),
         OtherRace_NL = ifelse(Latino == 1, 0, `Other Race`))
data.demo <- data.demo[, c('Key','Female','White_NL', 'Black_NL', 'Native_NL', 'PacificIslander_NL', 'Asian_NL', 'OtherRace_NL', 'Latino')]

#Merge class data with demo data
final.data <- merge(final.data, data.demo, by ='Key', all.x = TRUE)



#Open location data
data.loca <- read.csv('abcd_y_lt.csv', header=TRUE)
data.loca <- subset(data.loca, eventname == 'baseline_year_1_arm_1')
colnames(data.loca)[colnames(data.loca) %in% c('src_subject_id', 'site_id_l', 'interview_age')] <- c('Key', 'Location', 'Age')
data.loca$Location <- as.integer(gsub("site", "", data.loca$Location))
data.loca$RML <- ifelse(data.loca$Location == 1 | data.loca$Location == 2 | data.loca$Location == 6 | data.loca$Location == 7 | data.loca$Location == 8
                        | data.loca$Location == 9 | data.loca$Location == 10 | data.loca$Location == 12 | data.loca$Location == 13 | data.loca$Location == 14
                        | data.loca$Location == 17 | data.loca$Location == 19 | data.loca$Location == 20 | data.loca$Location == 21, 1, ifelse(data.loca$Location == 3 | data.loca$Location == 4 | data.loca$Location == 11 | data.loca$Location == 15 | data.loca$Location == 16, '', 0))

data.loca$MML <- ifelse(data.loca$Location == 3 | data.loca$Location == 4 | data.loca$Location == 11 | data.loca$Location == 15 | data.loca$Location == 16, 1, 0)
data.loca <- data.loca [,c('Key', 'Location', 'Age', 'RML','MML')]

#Merge location data with final data
final.data <- merge(final.data, data.loca, by = 'Key', all.x = TRUE)

#Split LCA into multiple columns
#1 is universal, 2 is low, 3 is internal, 4 is external

final.data <- final.data %>%
  mutate('Universal' = ifelse(class == 1, 1, 0),   # C1 is 1 if 'class' is 2, else 0
         'Internalizing' = ifelse(class == 3, 1, 0),   # C2 is 1 if 'class' is 3, else 0
         'Externalizing' = ifelse(class == 4, 1, 0))

#Same process for drug data

#Alcohol
data.alc <- read.csv('su_y_alc_exp.csv', header=TRUE)
data.alc <- subset(data.alc, eventname == '1_year_follow_up_y_arm_1')
data.alc <- data.alc[,c('src_subject_id',"aeq_positive_expectancies_ss", "aeq_negative_expectancies_ss")]
colnames(data.alc)[colnames(data.alc) %in% c('src_subject_id', "aeq_positive_expectancies_ss", "aeq_negative_expectancies_ss")] <-
  c('Key', 'Alc_P1', 'Alc_N1')
final.data <- merge(final.data, data.alc, by = 'Key', all.x=TRUE)

data.alc <- read.csv('su_y_alc_exp.csv', header=TRUE)
data.alc <- subset(data.alc, eventname == '2_year_follow_up_y_arm_1')
data.alc <- data.alc[,c('src_subject_id',"aeq_positive_expectancies_ss", "aeq_negative_expectancies_ss")]
colnames(data.alc)[colnames(data.alc) %in% c('src_subject_id', "aeq_positive_expectancies_ss", "aeq_negative_expectancies_ss")] <-
  c('Key', 'Alc_P2', 'Alc_N2')
final.data <- merge(final.data, data.alc, by = 'Key', all.x=TRUE)

data.alc <- read.csv('su_y_alc_exp.csv', header=TRUE)
data.alc <- subset(data.alc, eventname == '3_year_follow_up_y_arm_1')
data.alc <- data.alc[,c('src_subject_id',"aeq_positive_expectancies_ss", "aeq_negative_expectancies_ss")]
colnames(data.alc)[colnames(data.alc) %in% c('src_subject_id', "aeq_positive_expectancies_ss", "aeq_negative_expectancies_ss")] <-
  c('Key', 'Alc_P3', 'Alc_N3')
final.data <- merge(final.data, data.alc, by = 'Key', all.x=TRUE)




#Cannabis
data.can <- read.csv('su_y_can_exp.csv', header=TRUE)
data.can <- subset(data.can, eventname == '1_year_follow_up_y_arm_1')
data.can <- data.can[,c('src_subject_id',"meeq_positive_expectancies_ss", "meeq_negative_expectancies_ss")]
colnames(data.can)[colnames(data.can) %in% c('src_subject_id', "meeq_positive_expectancies_ss", "meeq_negative_expectancies_ss")] <-
  c('Key', 'Mar_P1', 'Mar_N1')
final.data <- merge(final.data, data.can, by = 'Key', all.x=TRUE)

data.can <- read.csv('su_y_can_exp.csv', header=TRUE)
data.can <- subset(data.can, eventname == '2_year_follow_up_y_arm_1')
data.can <- data.can[,c('src_subject_id',"meeq_positive_expectancies_ss", "meeq_negative_expectancies_ss")]
colnames(data.can)[colnames(data.can) %in% c('src_subject_id', "meeq_positive_expectancies_ss", "meeq_negative_expectancies_ss")] <-
  c('Key', 'Mar_P2', 'Mar_N2')
final.data <- merge(final.data, data.can, by = 'Key', all.x=TRUE)

data.can <- read.csv('su_y_can_exp.csv', header=TRUE)
data.can <- subset(data.can, eventname == '3_year_follow_up_y_arm_1')
data.can <- data.can[,c('src_subject_id',"meeq_positive_expectancies_ss", "meeq_negative_expectancies_ss")]
colnames(data.can)[colnames(data.can) %in% c('src_subject_id', "meeq_positive_expectancies_ss", "meeq_negative_expectancies_ss")] <-
  c('Key', 'Mar_P3', 'Mar_N3')
final.data <- merge(final.data, data.can, by = 'Key', all.x=TRUE)

#Nicotine
data.nic <- read.csv('su_y_nic_exp.csv', header=TRUE)
data.nic <- subset(data.nic, eventname == '1_year_follow_up_y_arm_1')
data.nic <- data.nic %>% 
  mutate('Nic_P1' = ascq_section_q01+ascq_section_q02+ascq_section_q03+ascq_section_q04,
         'Nic_N1' = ascq_section_q05+ascq_section_q06)
colnames(data.nic)[colnames(data.nic)=='src_subject_id'] <- 'Key'
data.nic <- data.nic[,c('Key', "Nic_P1", "Nic_N1")]
final.data <- merge(final.data, data.nic, by='Key', all.x=TRUE)

data.nic <- read.csv('su_y_nic_exp.csv', header=TRUE)
data.nic <- subset(data.nic, eventname == '2_year_follow_up_y_arm_1')
data.nic <- data.nic %>% 
  mutate('Nic_P2' = ascq_section_q01+ascq_section_q02+ascq_section_q03+ascq_section_q04,
         'Nic_N2' = ascq_section_q05+ascq_section_q06)
colnames(data.nic)[colnames(data.nic)=='src_subject_id'] <- 'Key'
data.nic <- data.nic[,c('Key', "Nic_P2", "Nic_N2")]
final.data <- merge(final.data, data.nic, by='Key', all.x=TRUE)

data.nic <- read.csv('su_y_nic_exp.csv', header=TRUE)
data.nic <- subset(data.nic, eventname == '3_year_follow_up_y_arm_1')
data.nic <- data.nic %>% 
  mutate('Nic_P3' = ascq_section_q01+ascq_section_q02+ascq_section_q03+ascq_section_q04,
         'Nic_N3' = ascq_section_q05+ascq_section_q06)
colnames(data.nic)[colnames(data.nic)=='src_subject_id'] <- 'Key'
data.nic <- data.nic[,c('Key', "Nic_P3", "Nic_N3")]
final.data <- merge(final.data, data.nic, by='Key', all.x=TRUE)

#Vape
data.vape <- read.csv('su_y_vap_exp.csv', header=TRUE)
data.vape <- subset(data.vape, eventname == '3_year_follow_up_y_arm_1')
data.vape <- data.vape %>% 
  mutate('Vape_P3' = su_vape_expect_01+su_vape_expect_02+su_vape_expect_03+su_vape_expect_04,
         'Vape_N3' = su_vape_expect_05+su_vape_expect_06+su_vape_expect_07+su_vape_expect_08)
colnames(data.vape)[colnames(data.vape)=='src_subject_id'] <- 'Key'
data.vape <- data.vape[,c('Key', "Vape_P3", "Vape_N3")]
final.data <- merge(final.data, data.vape, by='Key', all.x=TRUE)

final.data <- final.data %>%
  mutate('Alc_P1.A' = Alc_P1 / 3,
         'Alc_P2.A' = Alc_P2 / 3,
         'Alc_P3.A' = Alc_P3 / 3,
         'Alc_N1.A' = Alc_N1 / 3,
         'Alc_N2.A' = Alc_N2 / 3,
         'Alc_N3.A' = Alc_N3 / 3,
         'Nic_P1.A' = Nic_P1 / 4,
         'Nic_P2.A' = Nic_P2 / 4,
         'Nic_P3.A' = Nic_P3 / 4,
         'Nic_N1.A' = Nic_N1 / 2,
         'Nic_N2.A' = Nic_N2 / 2,
         'Nic_N3.A' = Nic_N3 / 2,
         'Mar_P1.A' = Mar_P1 / 3,
         'Mar_P2.A' = Mar_P2 / 3,
         'Mar_P3.A' = Mar_P3 / 3,
         'Mar_N1.A' = Mar_N1 / 3,
         'Mar_N2.A' = Mar_N2 / 3,
         'Mar_N3.A' = Mar_N3 / 3)


final.data <- merge(final.data, df, by='Key')
colnames(final.data)[colnames(final.data) %in% c('Alc_P1.y', 'Alc_P2.y', 'Alc_P3.y',
                                                 'Alc_N1.y', 'Alc_N2.y', 'Alc_N3.y',
                                                 'Nic_P1.y', 'Nic_P2.y', 'Nic_P3.y',
                                                 'Nic_N1.y', 'Nic_N2.y', 'Nic_N3.y',
                                                 'Mar_P1.y', 'Mar_P2.y', 'Mar_P3.y',
                                                 'Mar_N1.y', 'Mar_N2.y', 'Mar_N3.y')] <-
  c('Alc_P1.A', 'Alc_P2.A', 'Alc_P3.A',
    'Alc_N1.A', 'Alc_N2.A', 'Alc_N3.A',
    'Nic_P1.A', 'Nic_P2.A', 'Nic_P3.A',
    'Nic_N1.A', 'Nic_N2.A', 'Nic_N3.A',
    'Mar_P1.A', 'Mar_P2.A', 'Mar_P3.A',
    'Mar_N1.A', 'Mar_N2.A', 'Mar_N3.A')



final.data <- final.data %>%
  mutate('Alc_P.C21' = Alc_P2.x - Alc_P1.x,
         'Alc_P.C32' = Alc_P3.x - Alc_P2.x,
         'Alc_P.C31' = Alc_P3.x - Alc_P1.x,
         'Alc_N.C21' = Alc_N2.x - Alc_N1.x,
         'Alc_N.C32' = Alc_N3.x - Alc_N2.x,
         'Alc_N.C31' = Alc_N3.x - Alc_N1.x,
         'Nic_P.C21' = Nic_P2.x - Nic_P1.x,
         'Nic_P.C32' = Nic_P3.x - Nic_P2.x,
         'Nic_P.C31' = Nic_P3.x - Nic_P1.x,
         'Nic_N.C21' = Nic_N2.x - Nic_N1.x,
         'Nic_N.C32' = Nic_N3.x - Nic_N2.x,
         'Nic_N.C31' = Nic_N3.x - Nic_N1.x,
         'Mar_P.C21' = Mar_P2.x - Mar_P1.x,
         'Mar_P.C32' = Mar_P3.x - Mar_P2.x,
         'Mar_P.C31' = Mar_P3.x - Mar_P1.x,
         'Mar_N.C21' = Mar_N2.x - Mar_N1.x,
         'Mar_N.C32' = Mar_N3.x - Mar_N2.x,
         'Mar_N.C31' = Mar_N3.x - Mar_N1.x)

#Save final.data
write.csv(final.data, file = "output/study1/final_data.csv", row.names = TRUE)

## NOTE: Remaining statistical analysis performed on SPSS and not included here-----