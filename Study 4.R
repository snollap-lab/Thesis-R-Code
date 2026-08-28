################################################################################
################################################################################
################### Study 4 Code ###############################################
################################################################################


# Load Packages ----------------------------------------------------------------

library(clipr)
library(dplyr)
library(stringr)
library(tidyr)
library(broom)
library(tidyverse)
library(psych) 
library(DescTools) 
library(Hmisc) 
library(rcompanion) 
library(statmod) 
library(skimr) 
library(caret)
library(rsample)
library(mice) 
library(doParallel) 
library(doRNG)
library(recipes)
library(miselect) 
library(ggplot2) 
library(pROC)
library(diffdf)
library(data.table)
library(moments)
library(lmtest)
library(gtsummary)
library(arsenal) 
library(rstatix)
library(fastDummies)
library(randomForest)
library(cvAUC)
library(svMisc)
library(purrr)
library(blastula)


blastula::create_smtp_creds_key(
  id = "email_smtp",
  user = " email here", #  email here, for automated updates
  host = "smtp.gmail.com",
  port = 587,
  use_ssl = TRUE,
  overwrite = TRUE
)

# Prepare substance use data for Mplus  ----------------------------------------
dict <- read.csv('data/6.0/datadictionary.csv', header=TRUE)

su <- read.csv('data/6.0/su_y_su.csv', header=TRUE)
sui <- su %>% 
  select(participant_id, session_id, starts_with('su_y_sui')) %>%
  filter(endsWith(session_id, "A")) %>%
  rename(age0 = su_y_sui_age)

use_6mo <- su %>% 
  select(participant_id, session_id, starts_with('su_y_mysu')) %>%
  filter(endsWith(session_id, "M")) %>%
  subset(session_id != 'ses-05M') %>%
  select(participant_id, session_id, contains('use_') & !contains('branch')
         & !contains('001__01')) %>%
  select(-contains('_1mo_'), -contains('_002'), 
         -contains('_1wk_'), -contains('_first_')) %>%
  mutate( 
    across(.cols = starts_with('su_y_mysu') & ends_with('001'), 
           .fns = ~replace(.x, .x == 777, NA)))

use_yr <- sui %>%
  select(participant_id, session_id, contains('use_') & !contains('branch')
         & !contains('001__01')) %>%
  subset(session_id != 'ses-06A')


dict_names <- c(names(use_yr %>% select(-participant_id, -session_id, -contains('_l'),
                                        -contains('_l_'))), names(
                                          use_6mo %>% select(-participant_id, -session_id, -contains('_l'),
                                                             -contains('_l_'))))

use_dict <- data.frame(
  name = dict_names
) %>%
  left_join(dict[,c('name', 'description')], by = 'name')

#### Summarise use across drug categories
use_yr <- use_yr %>%
  mutate(
    alc_use = pmax(
      su_y_sui__use__alc_001,
      su_y_sui__use__alc__sip_001,
      su_y_sui__use__alc_001__l,
      su_y_sui__use__alc__sip_001__l,
      na.rm = TRUE),
    mj_use = pmax(
      su_y_sui__use__mj__blunt_001,        
      su_y_sui__use__mj__blunt_001__l,    
      su_y_sui__use__mj__conc_001,      
      su_y_sui__use__mj__conc__smoke_001__l,
      su_y_sui__use__mj__conc__vape_001__l,
      su_y_sui__use__mj__delta_001__l,  
      su_y_sui__use__mj__drink_001,     
      su_y_sui__use__mj__drink_001__l, 
      su_y_sui__use__mj__edbl_001,    
      su_y_sui__use__mj__edbl_001__l,
      su_y_sui__use__mj__puff_001,       
      su_y_sui__use__mj__puff_001__l,     
      su_y_sui__use__mj__smoke_001,       
      su_y_sui__use__mj__smoke_001__l,   
      su_y_sui__use__mj__synth_001,   
      su_y_sui__use__mj__synth_001__l,    
      su_y_sui__use__mj__tinc_001,      
      su_y_sui__use__mj__tinc_001__l,   
      su_y_sui__use__mj__vape_001__l,       
      na.rm = TRUE),
    nic_use = pmax(
      su_y_sui__use__nic__chew_001,
      su_y_sui__use__nic__chew_001__l,
      su_y_sui__use__nic__cig_001,
      su_y_sui__use__nic__cig_001__l,
      su_y_sui__use__nic__cigar_001,
      su_y_sui__use__nic__cigar_001__l,
      su_y_sui__use__nic__hookah_001,
      su_y_sui__use__nic__hookah_001__l,
      su_y_sui__use__nic__pipe_001,
      su_y_sui__use__nic__pipe_001__l,
      su_y_sui__use__nic__puff_001,
      su_y_sui__use__nic__puff_001__l,
      su_y_sui__use__nic__rplc_001,
      su_y_sui__use__nic__rplc_001__l,
      su_y_sui__use__nic__rplc_001__v01__l,
      su_y_sui__use__nic__vape_001,
      su_y_sui__use__nic__vape_001__l,
      na.rm = TRUE),
    oth_use = pmax(
      su_y_sui__use__coc_001,
      su_y_sui__use__coc_001__l,
      su_y_sui__use__cath_001,
      su_y_sui__use__cath_001__l,
      su_y_sui__use__cbd_001__l,
      su_y_sui__use__dxm_001,
      su_y_sui__use__dxm_001__l,
      su_y_sui__use__ghb_001,
      su_y_sui__use__ghb_001__l,
      su_y_sui__use__hall_001,
      su_y_sui__use__hall_001__l,
      su_y_sui__use__inh_001,
      su_y_sui__use__inh_001__l,
      su_y_sui__use__inh__sniff_001,
      su_y_sui__use__ket_001,
      su_y_sui__use__ket_001__l,
      su_y_sui__use__mdma_001,
      su_y_sui__use__mdma_001__l,
      su_y_sui__use__meth_001,
      su_y_sui__use__meth_001__l,
      su_y_sui__use__opi_001,
      su_y_sui__use__opi_001__l,
      su_y_sui__use__othdrg_001,
      su_y_sui__use__othdrg_001__l,
      su_y_sui__use__othdrg_002__01__l___1,
      su_y_sui__use__othdrg_002__01__l___2,
      su_y_sui__use__othdrg_002__01__l___3,
      su_y_sui__use__othdrg_002__01__l___4,
      su_y_sui__use__othdrg_002__01__l___5,
      su_y_sui__use__othdrg_002__01__l___6,
      su_y_sui__use__othdrg_002__01__l___7,
      su_y_sui__use__othdrg_002__01__l___8,
      su_y_sui__use__othdrg_002__01__l___9,
      su_y_sui__use__othdrg_002__l,
      su_y_sui__use__roid_001 ,
      su_y_sui__use__roid_001__l ,
      su_y_sui__use__rxopi_001 ,
      su_y_sui__use__rxopi_001__l ,
      su_y_sui__use__rxsed_001 ,
      su_y_sui__use__rxsed_001__l ,
      su_y_sui__use__rxstim_001 ,
      su_y_sui__use__rxstim_001__l ,
      su_y_sui__use__salv_001 ,
      su_y_sui__use__salv_001__l ,
      su_y_sui__use__shroom_001 ,
      su_y_sui__use__shroom_001__l ,
      su_y_sui__use__vape__flav_001__l ,
      su_y_sui__use__vape__oth_001__l ,
      na.rm = TRUE),
    qc_use = pmax(
      su_y_sui__use__qc_001 ,
      su_y_sui__use__qc_001__l ,
      na.rm = TRUE)
  )



use_loop <- c('alc_use', 'nic_use', 'mj_use', 'oth_use')

for (i in use_loop) {
  cat("\n---", i, "---\n")
  temp <- table(use_yr$session_id, use_yr[[i]])
  print(temp)
}


use_6mo <- use_6mo %>%
  mutate(
    alc_use = pmax(
      su_y_mysu__use__alc__sip__6mo_001,
      su_y_mysu__use__alc__6mo_001,
      su_y_mysu__use__alc_001,
      su_y_mysu__use__alc__sip_001,
      na.rm = TRUE),
    mj_use = pmax(
      su_y_mysu__use__mj__6mo_001,
      su_y_mysu__use__mj__conc__pipe__6mo_001,
      su_y_mysu__use__mj__conc__vape__6mo_001,
      su_y_mysu__use__mj__drink__tinc__6mo_001,
      su_y_mysu__use__mj__edbl__6mo_001,
      su_y_mysu__use__mj__smoke__6mo_001,
      su_y_mysu__use__mj__synth__6mo_001,
      su_y_mysu__use__mj__vape__6mo_001,       
      na.rm = TRUE),
    nic_use = pmax(
      su_y_mysu__use__nic__6mo_001,
      su_y_mysu__use__nic__chew__6mo_001,
      su_y_mysu__use__nic__cig__6mo_001,
      su_y_mysu__use__nic__cigar__6mo_001,
      su_y_mysu__use__nic__oth_001,
      su_y_mysu__use__nic__puff_001,
      su_y_mysu__use__nic__puff__6mo_001,
      su_y_mysu__use__nic__vape__6mo_001,
      na.rm = TRUE),
    oth_use = pmax(
      su_y_mysu__use__cath__6mo_001,
      su_y_mysu__use__cbd_001,
      su_y_mysu__use__coc__6mo_001,
      su_y_mysu__use__dxm__6mo_001,
      su_y_mysu__use__ghb__6mo_001,
      su_y_mysu__use__hall__6mo_001,
      su_y_mysu__use__illdrg_001,
      su_y_mysu__use__illdrg__6mo_001,
      su_y_mysu__use__inh__6mo_001,
      su_y_mysu__use__ket__6mo_001,
      su_y_mysu__use__mdma__6mo_001,
      su_y_mysu__use__meth__6mo_001,su_y_mysu__use__opi__6mo_001,
      su_y_mysu__use__qc__6mo_001,
      su_y_mysu__use__roid__6mo_001,
      su_y_mysu__use__rxopi__6mo_001,
      su_y_mysu__use__rxsed__6mo_001,
      su_y_mysu__use__rxstim__6mo_001,
      su_y_mysu__use__salv__6mo_001,
      su_y_mysu__use__shroom__6mo_001,
      su_y_mysu__use__vape__flav__6mo_001,
      na.rm = TRUE)
  )

for (i in use_loop) {
  cat("\n---", i, "---\n")
  temp <- table(use_6mo$session_id, use_6mo[[i]])
  print(temp)
}

use_yr <- use_yr %>%
  select(participant_id, session_id, ends_with('_use')) %>%
  select(-qc_use) %>%
  mutate(
    session_id = 
      case_when(
        session_id == 'ses-00A' ~ 0,
        session_id == 'ses-01A' ~ 2,
        session_id == 'ses-02A' ~ 4,
        session_id == 'ses-03A' ~ 6,
        session_id == 'ses-04A' ~ 8,
        session_id == 'ses-05A' ~ 10
      )) %>%
  rename(time_mplus = session_id)

use_6mo <- use_6mo %>%
  select(participant_id, session_id, ends_with('_use')) %>%
  mutate(
    session_id = 
      case_when(
        session_id == 'ses-00M' ~ 1,
        session_id == 'ses-01M' ~ 3,
        session_id == 'ses-02M' ~ 5,
        session_id == 'ses-03M' ~ 7,
        session_id == 'ses-04M' ~ 9
      )) %>%
  rename(time_mplus = session_id)

data.long <- rbind(use_yr, use_6mo) %>% arrange(participant_id, time_mplus)

participants <- unique(data.long$participant_id)
times <- c(0,1,2,3,4,5,6,7,8,9,10)
skeleton <- expand.grid(participant_id = participants, time_mplus = times)

data.long <- skeleton %>%
  left_join(data.long, by = c("participant_id", "time_mplus")) %>%
  arrange(participant_id, time_mplus)

rm(skeleton, use_6mo, use_yr, i, dict_names, participants, temp, times)

####### Basic Exclusions -------------------------------------------------------
# Identify baseline users
baseline_use <- data.long %>%
  filter(time_mplus == 0) %>%
  select(-time_mplus) %>%
  mutate(
    alc_baseline = alc_use,
    mj_baseline = mj_use,
    nic_baseline = nic_use,
    oth_baseline = oth_use
  ) %>%
  select(-ends_with('use'))



#Identify participants with >5 missing timepoints
use_vars <- c("alc_use", "nic_use", "mj_use", "oth_use")

exclude_ids <- data.long %>%
  mutate(all_missing = if_all(all_of(use_vars), ~ is.na(.x))) %>%  # mark rows where *all* uses are NA
  group_by(participant_id) %>%
  summarise(n_all_missing = sum(all_missing), .groups = "drop") %>%
  filter(n_all_missing > 5) %>%
  pull(participant_id)

#940 participants with >5 missing timepoints, 11868 - 940 = 10928

# exclude those participants entirely
data.long <- data.long %>%
  filter(!participant_id %in% exclude_ids)

id_temp <- data.long %>%
  pull(participant_id)


## MRI Info + Exclusions -------------------------------------------------------
mri_info <- read.csv('data/6.0/mr_y_qc__mot.csv', header = TRUE) %>%
  filter(session_id == 'ses-00A') %>%
  filter(participant_id %in% id_temp) %>%
  select(participant_id, contains('mot_mean'), contains('indicator')) %>%
  select(-contains('unwarp'), -contains('disp'))


temp <- names(mri_info %>% select(contains('indicator')))
for(i in 1:length(temp)) {
  col <- temp[i]
  col2 <- mri_info[,col]
  cat('---- ', col, ' ----\n')
  print(table(col2, useNA = 'ifany'))
  cat('----------------\n')
}

### Largest % of missing MRI data from tfmri

mri_info <- mri_info %>%
  select(-contains('tfmri'))

#Make exclusion summary variables
mri_info <- mri_info %>%
  rowwise() %>%
  mutate(
    exclude_mot = ifelse(max(mr_y_qc__mot__rsfmri__mot_mean , na.rm = TRUE) > 0.5, 1, 0),
    exclude_qc  = ifelse(min(c_across(contains('indicator')), na.rm = TRUE) < 1,   1, 0)
  ) %>%
  ungroup() %>%
  select(participant_id, exclude_mot, exclude_qc)

data.long <- data.long %>%
  left_join(mri_info, by = 'participant_id')

table(data.long$exclude_mot, useNA = 'ifany')
table(data.long$exclude_qc, useNA = 'ifany')
table(data.long$exclude_qc, data.long$exclude_mot, useNA = 'ifany')

#7456 participants with no MRI exclusions

data.long <- data.long %>%
  filter(exclude_qc == 0 & exclude_mot == 0) %>%
  select(-exclude_qc, -exclude_mot)

rm(mri_info, col, col2, i, id_temp, temp)

## MPLUS ID's ------------------------------------------------------------------

for (i in use_loop) {
  cat("\n---", i, "---\n")
  temp <- table(data.long$time_mplus, data.long[[i]], useNA = 'ifany')
  print(temp)
}

#id's for mplus
id_lookup <- data.frame(
  mplus_id = as.numeric(factor(unique(data.long$participant_id))),
  participant_id = unique(data.long$participant_id)
)
write.csv(id_lookup, "output/study4/mPlus/id_lookup.csv", row.names = FALSE)



#Pivot Wider for Mplus

#Add age for covariation
age_bl <- sui %>%
  filter(session_id == 'ses-00A') %>%
  select(participant_id, age0) %>%
  mutate(age0 = as.numeric(scale(age0)))

# Pivot to wide
data.wide <- data.long %>%
  pivot_wider(
    id_cols = participant_id,
    names_from = time_mplus,
    values_from = c(alc_use, mj_use, nic_use, oth_use),
    names_sep = "_"
  ) %>%
  arrange(participant_id) %>%
  left_join(age_bl, by = 'participant_id') %>%
  select(participant_id, age0, everything())

names(data.wide) <- gsub("_use_", "", names(data.wide))



#Replace NA with distinguishable code for MPlus
data.wide[is.na(data.wide)] <- -9999

data.mplus <- data.wide %>%
  select(-alc0, -mj0, -nic0, -oth0) %>%
  left_join(id_lookup, by = 'participant_id') %>%
  select(mplus_id, everything()) %>%
  select(-participant_id)

write.table(
  data.mplus,
  file = "output/study4/mPlus/poly_mplus_wBL.dat",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)





rm(baseline_use, exclude_ids, i, temp, age_bl, data.mplus, sui, use_dict, use_loop, use_vars)


#### Results -------------------------------------------------------------------

#ADD aBIC, DO LRT FOR NG = 3:5
##Information Criteria
mplus_ic <- data.frame(
  ng = 1:10,
  AIC = NA,
  BIC = NA,
  aBIC = NA,
  logLike = NA,
  entropy = NA
)

# Combine all ng vectors into a list, info manually retrieved from mplus output
ng_list <- list(
  c(118816.961, 118886.129, 118854.351, -59398.481, NA),  # ng1
  c(84893.908, 85018.405, 84961.205, -42428.954, 0.880),  # ng2
  c(80820.450, 81014.112, 80925.134, -40382.225, 0.878),  # ng3
  c(79076.220, 79339.047, 79218.291, -39500.110, 0.827),  # ng4
  c(78092.239, 78424.231, 78271.697, -38998.119, 0.799),  # ng5
  c(77486.450, 77887.607, 77703.295, -38685.225, 0.809),  # ng6
  c(77136.650, 77606.972, 77390.883, -38500.325, 0.799),  # ng7
  c(76889.399, 77428.887, 77181.020, -38366.700, 0.805),  # ng8
  c(76792.801, 77401.453, 77121.808, -38308.400, 0.835),  # ng9
  c(76656.283, 77334.101, 77022.678, -38230.141, 0.687)   # ng10
)

# Loop through each class and fill mplus_ic
for (i in 1:10) {
  mplus_ic[mplus_ic$ng == i, c("AIC","BIC",'aBIC',"logLike", 'entropy')] <- ng_list[[i]]
}

mplus_poly_ic_long <- mplus_ic %>%
  pivot_longer(
    cols = c(AIC, BIC, aBIC, logLike, entropy),
    names_to = "Criterion",
    values_to = "Value"
  )

mplus_ic <- mplus_ic %>%
  arrange(ng) %>%
  mutate(
    deltaAIC = ifelse( ng == 1, NA, AIC - lag(AIC)),
    deltaBIC = ifelse( ng == 1, NA, BIC - lag(BIC)),
    deltaABIC = ifelse(ng == 1, NA, aBIC - lag(aBIC)),
    deltalogLike = ifelse( ng == 1, NA, logLike - lag(logLike)),
    deltaEntropy = ifelse(ng == 1, NA, entropy - lag(entropy))
  ) %>%
  select(ng, AIC, deltaAIC, BIC, deltaBIC, aBIC, deltaABIC, logLike, deltalogLike,
         entropy, deltaEntropy)



# Define the file pattern and empty list for results
ng_files <- sprintf("output/study4/mPlus/poly_wBL/poly_ng%d_mplus_wMRI.dat", 1:10)
assignment_counts <- list()

assignment_counts <- lapply(ng_files, function(file) {
  message("Reading: ", file)
  dat <- read.table(file, header = TRUE, fill = TRUE)
  lst <- ncol(dat)
  last_col <- dat[[lst]]
  table(last_col)
})

names(assignment_counts) <- paste0("ng", 1:10)
assignment_counts

mplus_ic$class_dist <- assignment_counts

combos <- lapply(1:10, function(num) {
  choose(num, 2)
})

mplus_ic$unique_tests <- combos

mplus_ic <- mplus_ic %>%
  select(ng, unique_tests, class_dist, everything())

mplus_ic$class_dist <- sapply(mplus_ic$class_dist, 
                              function(x) paste(x, collapse = "; "))

mplus_ic$unique_tests <- as.character(mplus_ic$unique_tests)


write.csv(mplus_ic, 'output/study4/mplus_information_wBL.csv', row.names = FALSE)


#### Add best model -----------------------------------------------------------

four.model <- read.table("mPlus/poly_wBL/poly_ng4_mplus_wMRI.dat", header = TRUE)
comb_ng4_assignments <- data.frame(
  mplus_id = as.numeric(four.model$X1),
  ng4_class = factor(four.model[[ncol(four.model)]])
)

id_lookup <- id_lookup %>%
  left_join(comb_ng4_assignments, by = 'mplus_id')


data.long <- data.long %>%
  left_join(id_lookup[,c('participant_id', 'ng4_class')], by = 'participant_id')



data.wide <- data.wide %>%
  left_join(id_lookup[,c('participant_id', 'ng4_class')], by = 'participant_id') %>%
  select(participant_id, everything())

data.wide[data.wide == -9999] <- NA

rm(four.model, comb_ng4_assignments, id_lookup)

saveRDS(data.long, 'output/study4/data_long.rds')


# Wrangle explanatory data --------------------------------------------------------------
dict <- read.csv('data/6.0/datadictionary.csv', header=TRUE)

##Independent Variables: Self and Peer Involvement with Substance Use #####
su <- read.csv('data/6.0/su_y_su.csv', header=TRUE)

#   substance-related predictors
su_gating_pred <- su %>% 
  select(participant_id, session_id, starts_with('su_y_sui')) %>% 
  subset(session_id == 'ses-00A') %>% 
  select(
    participant_id, contains('hrd')) 

colnames(su_gating_pred) <- paste(
  colnames(su_gating_pred), 'base', sep = '.')

su_gating_pred <- su_gating_pred %>% 
  rename(participant_id = participant_id.base)

################################################################################
#                           Peer Group Deviance                                #
################################################################################

# import questionnaire
su_y_peerdevia <- su %>% 
  select(participant_id, session_id, starts_with('su_y_pgd')) %>% 
  subset(session_id == 'ses-00A') %>% 
  select(participant_id,
         su_y_pgd_001, su_y_pgd_002, 
         su_y_pgd_003, su_y_pgd_004, 
         su_y_pgd_005, su_y_pgd_006,
         su_y_pgd_007, su_y_pgd_008, 
         su_y_pgd_009) %>%
  left_join(su_gating_pred, by = 'participant_id')

# re-code due to gating criteria
su_y_peerdevia <- su_y_peerdevia %>% 
  mutate(
    peer_dev_1 = ifelse(
      su_y_sui__hrd__nic_001.base == 0 & is.na(su_y_pgd_001), 0, su_y_pgd_001),
    peer_dev_2 = ifelse(
      su_y_sui__hrd__alc_001.base == 0 & is.na(su_y_pgd_002), 0, su_y_pgd_002),
    peer_dev_3 = ifelse(
      su_y_sui__hrd__alc_001.base == 0 & is.na(su_y_pgd_003), 0, su_y_pgd_003),
    peer_multi = ifelse(
      su_y_sui__hrd__alc_001.base == 1 | su_y_sui__hrd__nic_001.base == 1 | su_y_sui__hrd__mj_001.base == 1 | 
        su_y_sui__hrd__mj__synth_001.base == 1 | su_y_sui__hrd__coc_001.base == 1 | 
        su_y_sui__hrd__inh_001.base == 1, 1, 0),
    peer_dev_4 = ifelse(
      peer_multi == 0 & is.na(su_y_pgd_004), 0, su_y_pgd_004),
    peer_dev_5 = ifelse(
      su_y_sui__hrd__mj_001.base == 0 & is.na(su_y_pgd_005), 0, su_y_pgd_005),
    peer_dev_6 = ifelse(
      su_y_sui__hrd__inh_001.base == 0 & is.na(su_y_pgd_006), 0, su_y_pgd_006),
    peer_dev_7 = ifelse(
      su_y_sui__hrd_001.base == 0 & is.na(su_y_pgd_007), 0, su_y_pgd_007),
    peer_dev_8 = ifelse(
      su_y_sui__hrd_001.base == 0 & is.na(su_y_pgd_008), 0, su_y_pgd_008),
    peer_dev_9 = ifelse(
      su_y_sui__hrd__nic_001.base == 0 & is.na(su_y_pgd_009), 0, su_y_pgd_009))

# create summary scores
su_y_peerdevia <- su_y_peerdevia %>% 
  mutate(
    peer_alc = peer_dev_2 + peer_dev_3,
    peer_tob = peer_dev_1 + peer_dev_9,
    peer_cb = peer_dev_5, 
    peer_other = peer_dev_6 + peer_dev_7, 
    peer_prob = peer_dev_4 + peer_dev_8) 

# subset relevant variables 
su_y_peerdevia <- su_y_peerdevia %>% 
  select(participant_id, peer_alc, peer_tob, peer_cb, peer_other, peer_prob)

# merge w/full dataset
dv_names <- c('alc1', 'alc2', 'alc3',
              'alc4', 'alc5', 'alc6', 
              'alc7', 'alc8', 'alc9', 'alc10',
              'nic1', 'nic2', 'nic3',
              'nic4', 'nic5', 'nic6', 
              'nic7', 'nic8', 'nic9', 'nic10',
              'mj1', 'mj2', 'mj3',
              'mj4', 'mj5', 'mj6', 
              'mj7', 'mj8', 'mj9', 'mj10',
              'oth1', 'oth2', 'oth3',
              'oth4', 'oth5', 'oth6', 
              'oth7', 'oth8', 'oth9', 'oth10')

data.2 <- data.wide %>% 
  left_join(su_y_peerdevia, by = 'participant_id') %>%
  select(-all_of(dv_names), -age0) %>%
  mutate(alc0 = ifelse(is.na(alc0), 0, alc0),
         mj0 = ifelse(is.na(mj0), 0, mj0),
         nic0 = ifelse(is.na(nic0), 0, nic0),
         oth0 = ifelse(is.na(oth0), 0, oth0)) %>%
  mutate(alc0 = factor(alc0, levels = c(0,1), labels = c('No', 'Yes')),
         mj0 = factor(mj0, levels = c(0,1), labels = c('No', 'Yes')),
         nic0 = factor(nic0, levels = c(0,1), labels = c('No', 'Yes')),
         oth0 = factor(oth0, levels = c(0,1), labels = c('No', 'Yes')))

rm(su_y_peerdevia)

################################################################################
#                             Intent to Use
################################################################################

# import questionnaire
su_y_path_intuse <- su %>% 
  select(participant_id, session_id, starts_with('su_y_itu')) %>% 
  subset(session_id == 'ses-00A') %>% 
  select(participant_id,
         su_y_itu__alc_001, su_y_itu__alc_002, su_y_itu__alc_003, su_y_itu__mj_001,
         su_y_itu__mj_002, su_y_itu__mj_003, su_y_itu__nic_001, su_y_itu__nic_002,
         su_y_itu__nic_003) %>%  
  left_join(su_gating_pred, by = 'participant_id')


# re-code due to gating criteria
su_y_path_intuse <- su_y_path_intuse %>% 
  mutate(
    path_1 = ifelse(
      (su_y_sui__hrd__nic_001.base == 0 & is.na(su_y_itu__nic_001)), 4, su_y_itu__nic_001),
    path_2 = ifelse(
      (su_y_sui__hrd__alc_001.base == 0 & is.na(su_y_itu__alc_001)), 4, su_y_itu__alc_001),
    path_3 = ifelse(
      (su_y_sui__hrd__mj_001.base == 0 & is.na(su_y_itu__mj_001)), 4, su_y_itu__mj_001), 
    path_4 = ifelse(
      (su_y_sui__hrd__nic_001.base == 0 & is.na(su_y_itu__nic_002)), 4, su_y_itu__nic_002),
    path_5 = ifelse(
      (su_y_sui__hrd__alc_001.base == 0 & is.na(su_y_itu__alc_002)), 4, su_y_itu__alc_002),
    path_6 = ifelse(
      (su_y_sui__hrd__mj_001.base == 0 & is.na(su_y_itu__mj_002)), 4, su_y_itu__mj_002),
    path_7 = ifelse(
      (su_y_sui__hrd__nic_001.base == 0 & is.na(su_y_itu__nic_003)), 4, su_y_itu__nic_003),
    path_8 = ifelse(
      (su_y_sui__hrd__alc_001.base == 0 & is.na(su_y_itu__alc_003)), 4, su_y_itu__alc_003),
    path_9 = ifelse(
      (su_y_sui__hrd__mj_001.base == 0 & is.na(su_y_itu__mj_003)), 4, su_y_itu__mj_003)) 


# re-code Don't Know and Refuse to Answer
su_y_path_intuse <- su_y_path_intuse %>% 
  mutate(
    across(.cols = path_1:path_9,  
           .fns = ~replace(.x, .x %in% c(5, 6, 999, 777), NA))) %>% 
  # re-code to be on the same scale as peer deviance questions 
  # greater scores, greater intent to use
  mutate(
    path_1 = 5 - path_1,
    path_2 = 5 - path_2,
    path_3 = 5 - path_3,
    path_4 = 5 - path_4,
    path_5 = 5 - path_5,
    path_6 = 5 - path_6,
    path_7 = 5 - path_7,
    path_8 = 5 - path_8,
    path_9 = 5 - path_9)

# create summary scores
su_y_path_intuse <- su_y_path_intuse %>% 
  mutate(
    path_alc = path_2 + path_5 + path_8, 
    path_tob = path_1 + path_4 + path_7, 
    path_cb = path_3 + path_6 + path_9)

# subset relevant variables 
su_y_path_intuse <- su_y_path_intuse %>% 
  select(participant_id, path_alc, path_tob, path_cb)

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(su_y_path_intuse, by = 'participant_id')

rm(su_y_path_intuse)

## Independent Variables: Parenting Behaviors  ----------------------------------

################################################################################
#             Community Risk and Protective Factors Survey
################################################################################
fce <- read.csv('data/6.0/fc_p_aclt.csv')

# import questionnaire
su_p_crpf <- su %>% 
  select(participant_id, session_id, starts_with('su_p_crpf'))  %>% 
  subset(session_id == 'ses-00A') %>% 
  select(participant_id,
         su_p_crpf__alc_001, su_p_crpf__illdrg_001, su_p_crpf__mj_001,
         su_p_crpf__nic__cig_001, su_p_crpf__nic__vape_001)
# re- code values of (4) 'Don't Know' to NA
su_p_crpf <- su_p_crpf %>% 
  mutate(
    across(
      .cols = starts_with('su_p_crpf'),
      .fns = ~ ifelse(.x == 4, NA, .x))) %>% 
  mutate(
    crpf = su_p_crpf__alc_001 + su_p_crpf__illdrg_001 + su_p_crpf__mj_001 + 
      su_p_crpf__nic__cig_001 + su_p_crpf__nic__vape_001) %>% 
  select(participant_id, crpf)

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(su_p_crpf, by = 'participant_id') 

rm(su_p_crpf)

################################################################################
#                           Parent Rules
################################################################################

su_p_pr <- su %>% 
  select(participant_id, session_id, starts_with('su_p_rule'))  %>% 
  subset(session_id == 'ses-00A') %>% 
  select(participant_id,su_p_rule__alc_001, su_p_rule__mj_001, su_p_rule__nic__cig_001,) %>% 
  mutate(par_rules = su_p_rule__alc_001 + su_p_rule__mj_001 + su_p_rule__nic__cig_001) 


# retain total score
su_p_pr <- su_p_pr %>% 
  select(participant_id, par_rules)

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(su_p_pr, by = 'participant_id')

rm(su_p_pr, su_gating, su_gating_pred, su_y_mypi, su_y_sui)

################################################################################
#                     Parental Monitoring Questionnaire                        #
################################################################################

ce_y_pm <- fce %>% 
  select(participant_id, session_id, starts_with('fc_y_pm')) %>% 
  subset (session_id == 'ses-00A') %>% 
  select(participant_id, fc_y_pm_mean) 

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(ce_y_pm, by = 'participant_id') 

rm(ce_y_pm)

################################################################################
#                       Family Environment Scale                               #
################################################################################

ce_y_fes <- fce %>% 
  select(participant_id, session_id, starts_with('fc_y_fes')) %>% 
  subset (session_id == 'ses-00A') %>% 
  select(participant_id, fc_y_fes__confl_mean) 

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(ce_y_fes, by = 'participant_id') 

rm(ce_y_fes)

################################################################################
#                            CRPBI                                             #
################################################################################

ce_y_crpbi <- fce %>% 
  select(participant_id, session_id, starts_with('fc_y_crpbi')) %>%
  subset (session_id == 'ses-00A') %>% 
  rename(
    crpbi_y_ss_parent = fc_y_crpbi__cg1_mean
  ) %>%
  select(participant_id, crpbi_y_ss_parent) 

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(ce_y_crpbi, by = 'participant_id') 

rm(ce_y_crpbi)


## Independent Variables: Demographics  -----------------------------------------

################################################################################
#                         ABCD General: Baseline Age
################################################################################
dem <- read.csv('data/6.0/ab_g_dyn.csv')

abcd_y_lt <- dem %>% 
  select(participant_id, session_id, starts_with('ab_g_dyn'))   %>% 
  subset (session_id == 'ses-00A') %>% 
  rename(age_baseline = ab_g_dyn__visit_age) %>% 
  select(participant_id, age_baseline) 
abcd_y_lt$age_baseline <- round(abcd_y_lt$age_baseline, digits = 2)

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(abcd_y_lt, by = 'participant_id')

rm(abcd_y_lt)

################################################################################
#                 ABCD Parent Demographics Questionnaire
################################################################################

abcd_p_demo <- dem %>% 
  select(participant_id, session_id, starts_with('ab_p_demo')) %>% 
  subset (session_id == 'ses-00A') %>%
  left_join(dem %>% select(participant_id, session_id, starts_with('ab_g_stc')) %>% 
              subset (session_id == 'ses-00A'), by = 'participant_id' ) %>%
  select(-contains('session_id')) %>%
  select(participant_id, ab_p_demo_age, ab_g_stc__cohort_sex,
         ab_p_demo__income__hhold_001,
         ab_p_demo__race_001___10:ab_p_demo__race_001___25, ab_p_demo__race_001___0, 
         ab_p_demo__race_001___777, ab_p_demo__race_001___999, 
         ab_p_demo__ethn_001, ab_p_demo__relig_001, ab_p_demo__edu__slf_001) %>%
  select(-contains('ntvlang'), -contains('prtnr')) %>% 
  mutate(
    ab_p_demo__income__hhold_001 = case_when(
      ab_p_demo__income__hhold_001 <= 4 ~ 1,
      ab_p_demo__income__hhold_001 >= 5 & ab_p_demo__income__hhold_001 < 7 ~ 2,
      ab_p_demo__income__hhold_001 >= 7 & ab_p_demo__income__hhold_001 < 9 ~ 3,
      ab_p_demo__income__hhold_001 >= 9  & ab_p_demo__income__hhold_001 < 777 ~ 4,
      ab_p_demo__income__hhold_001 == 999 | ab_p_demo__income__hhold_001 == 777 ~ NA
    )
  ) %>%
  mutate(
    sex = factor(ab_g_stc__cohort_sex, 
                 levels = c(1,2), 
                 labels = c('Male', 'Female')),
    eth_hisp = factor(ab_p_demo__ethn_001, 
                      levels = c(0,1), 
                      labels = c('non_Hispanic', 'Hispanic')),
    income = factor(ab_p_demo__income__hhold_001, 
                    levels = c(1, 2, 3, 4),
                    labels = c('inc_1', 'inc_2', 'inc_3', 'inc_4')),
    religion = factor(ab_p_demo__relig_001,
                      levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 
                                 14, 15, 16, 17),
                      labels = c('rp_1', 'rp_2', 'rp_3', 'rp_4', 'rp_5',
                                 'rp_6', 'rp_7', 'rp_8', 'rp_9', 'rp_10',
                                 'rp_11', 'rp_12', 'rp_13', 'rp_14', 
                                 'rp_15', 'rp_16', 'rp_17')))


# recode values of (777) Refuse to Answer and (999) as NA
abcd_p_demo <- abcd_p_demo %>%  
  mutate(
    across(.cols = -participant_id,  
           .fns = ~replace(.x, .x %in% c(999, 777), NA)))

# create summary education variable 
abcd_p_demo <- abcd_p_demo %>% 
  mutate(
    p_edu = case_when(
      ab_p_demo__edu__slf_001 < 13 ~ 'Less_than_HS_Degree_GED_Equivalent',
      ab_p_demo__edu__slf_001 == 13 ~ 'HS_Graduate_GED_Equivalent',
      ab_p_demo__edu__slf_001 == 14 ~ 'HS_Graduate_GED_Equivalent',
      ab_p_demo__edu__slf_001 == 15 ~ 'Some_College_or_Associates_Degree',
      ab_p_demo__edu__slf_001 == 16 ~ 'Some_College_or_Associates_Degree',
      ab_p_demo__edu__slf_001 == 17 ~ 'Some_College_or_Associates_Degree',
      ab_p_demo__edu__slf_001 == 18 ~ 'Bachelors_Degree',
      ab_p_demo__edu__slf_001 == 19 ~ 'Masters_Degree',
      ab_p_demo__edu__slf_001 >= 20 ~ 'Professional_School_or_Doctoral_Degree')) %>% 
  mutate(p_edu = as.factor(p_edu))
table(abcd_p_demo$p_edu, useNA = 'ifany')

# create summary race variables           
abcd_p_demo <- abcd_p_demo %>%   
  mutate(
    White = ifelse(ab_p_demo__race_001___10 == 1, 1, 0),
    Black = ifelse(ab_p_demo__race_001___11 == 1, 1, 0),
    Asian = ifelse(ab_p_demo__race_001___18 == 1 | ab_p_demo__race_001___19 == 1 | 
                     ab_p_demo__race_001___20 == 1 | ab_p_demo__race_001___21 == 1 | 
                     ab_p_demo__race_001___22 == 1 | ab_p_demo__race_001___23 == 1 |
                     ab_p_demo__race_001___24 == 1, 1, 0),
    AIAN = ifelse( #AIAN: American Indian Alaskan Native
      ab_p_demo__race_001___12 == 1 | ab_p_demo__race_001___13 == 1, 1, 0),
    NHPI = ifelse( #NHPI: Native Hawaiian and Other Pacific
      ab_p_demo__race_001___14 == 1 | ab_p_demo__race_001___15 == 1 | 
        ab_p_demo__race_001___16 == 1 | ab_p_demo__race_001___17, 1, 0),
    Other = ifelse(ab_p_demo__race_001___25 == 1 | ab_p_demo__race_001___0 == 1, 1, 0)) %>% 
  mutate(
    MultiRacial = ifelse(
      (White + Black + Asian + AIAN + NHPI + Other) > 1, 1, 0)) %>% 
  mutate(
    race_4l = case_when(
      White == 1 & (Black == 0 & Asian == 0 & AIAN == 0 & NHPI == 0 & Other == 0
                    & MultiRacial == 0) ~ 'White',
      Black == 1 & (White == 0 & Asian == 0 & AIAN == 0 & NHPI == 0 & Other == 0
                    & MultiRacial == 0) ~ 'Black',
      Asian == 1 & (White == 0 & Black == 0 & AIAN == 0 & NHPI == 0 & Other == 0
                    & MultiRacial == 0) ~ 'Asian',
      AIAN == 1 | NHPI == 1 | Other == 1 | MultiRacial == 1 ~ 
        'Other_MultiRacial'),
    race_4l = factor(race_4l)) 

# subset relevant variables
abcd_p_demo <- abcd_p_demo %>% 
  select(participant_id, 
         sex, race_4l, eth_hisp, income, religion, p_edu)
str(abcd_p_demo)

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(abcd_p_demo, by = 'participant_id')

rm(abcd_p_demo)

data.2 <- data.2 %>% 
  rename(sex_2l = sex) %>%  
  mutate(sex_2l = as.factor(sex_2l))


## Independent Variables: Mental Health  ----------------------------------------
mh <- read.csv('data/6.0/mh_p_cbcl.csv')


################################################################################
#                    Child Behavior Checklist Scores (CBCL)                    #
################################################################################

mh_p_cbcl <- mh %>% 
  select(participant_id, session_id, starts_with('mh_p_cbcl')) %>% 
  subset (session_id == 'ses-00A') %>% 
  select(participant_id, starts_with('mh_p_cbcl__synd') & ends_with('_tscore'))

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(mh_p_cbcl, by = 'participant_id')

rm(mh_p_cbcl)

################################################################################
#                           Family History                                     #
################################################################################

mh_p_fhx <- mh %>% 
  select(participant_id, session_id, starts_with('mh_p_famhx')) %>% 
  subset (session_id == 'ses-00A') %>% 
  select(participant_id,
         mh_p_famhx__alc__fath_001___0, mh_p_famhx__alc__pat__gfath_001___0,
         mh_p_famhx__alc__pat__gmoth_001___0, mh_p_famhx__alc__moth_001___0,
         mh_p_famhx__alc__mat__gfath_001___0, mh_p_famhx__alc__mat__gmoth_001___0,
         mh_p_famhx__drg__fath_001___0, mh_p_famhx__drg__pat__gfath_001___0,
         mh_p_famhx__drg__pat__gmoth_001___0, mh_p_famhx__drg__moth_001___0,
         mh_p_famhx__drg__mat__gfath_001___0, mh_p_famhx__drg__mat__gmoth_001___0,)



mh_p_fhx <- mh_p_fhx %>% 
  mutate(
    across(.cols = starts_with('mh_p_famhx__alc__fath') | 
             starts_with('mh_p_famhx__drg__fath') |
             starts_with('mh_p_famhx__alc__moth') |
             starts_with('mh_p_famhx__drg__moth'), 
           .fns = ~replace(.x, .x == 1, .50))) %>% 
  mutate(
    across(.cols = starts_with('mh_p_famhx__alc__pat') |
             starts_with('mh_p_famhx__drg__pat') |
             starts_with('mh_p_famhx__alc__mat') | starts_with('mh_p_famhx__drg__mat'), 
           .fns = ~replace(.x, .x == 1, .25))) 

mh_p_fhx$count_na <- rowSums(is.na(select(mh_p_fhx, starts_with('mh_p_famhx'))))

# create density score
mh_p_fhx$mh_density <- rowSums(
  (select(mh_p_fhx, starts_with('mh_p_famhx'))), na.rm = TRUE)
table(mh_p_fhx$mh_density, useNA = 'ifany')

# if missing all 12 variables used in density score, code as NA 
mh_p_fhx <- mh_p_fhx %>% 
  mutate(mh_density = ifelse(count_na == 12, NA, mh_density))
table(mh_p_fhx$mh_density, useNA = 'ifany') # (NA) = 127

# merge w/full dataset 
mh_p_fhx <- mh_p_fhx %>%
  select(participant_id, mh_density)

data.2 <- data.2 %>% 
  left_join(mh_p_fhx, by = 'participant_id')

rm(mh_p_fhx)

################################################################################
#                       Prodromal Psychosis Symptoms                           #
################################################################################

mh_y_pps <- mh %>% 
  select(participant_id, session_id, starts_with('mh_y_pps')) %>% 
  subset (session_id == 'ses-00A') %>% 
  select(participant_id, mh_y_pps__severity_score)

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(mh_y_pps, by = 'participant_id')

rm(mh_y_pps)

################################################################################
#                               UPPS-P                                         #
################################################################################

mh_y_upps <- mh %>% 
  select(participant_id, session_id, starts_with('mh_y_upps')) %>% 
  subset (session_id == 'ses-00A') %>% 
  select(participant_id, mh_y_upps__nurg_sum,
         mh_y_upps__pers_sum, mh_y_upps__plan_sum,
         mh_y_upps__sens_sum, mh_y_upps__purg_sum)

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(mh_y_upps, by = 'participant_id')

rm(mh_y_upps)

################################################################################
#                               BIS/BAS                                        #
################################################################################

mh_y_bisbas <- mh %>% 
  select(participant_id, session_id, starts_with('mh_y_bisbas')) %>% 
  subset (session_id == 'ses-00A') %>% 
  select(participant_id, mh_y_bisbas__bis_sum, mh_y_bisbas__bas__dr_sum, 
         mh_y_bisbas__bas__fs_sum, mh_y_bisbas__bas__rr_sum)

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(mh_y_bisbas, by = 'participant_id')

rm(mh_y_bisbas)

## Independent Variables: Physical Health  --------------------------------------

ph <- read.csv('data/6.0/ph_p_pds.csv')

################################################################################
#                 Sports Activities Involvement Questionnaire                  #
################################################################################

ph_p_saiq <- ph %>% 
  select(participant_id, session_id, starts_with('ph_p_saiq')) %>% 
  subset (session_id == 'ses-00A') %>% 
  select(participant_id, starts_with('ph_p_saiq__activs') & ends_with('001')) 


# recode values of '999' to NA
ph_p_saiq <- ph_p_saiq %>% 
  mutate(
    across(.cols = starts_with('ph_p_saiq__activs'), 
           .fns = ~replace(.x, .x == 999, NA)))

# count number missing recreational activity vars (n=29 activities)
ph_p_saiq$count_na <- rowSums(is.na(select(ph_p_saiq, starts_with('ph_p_saiq__activs'))))

# create continuous and binary recreational activity variables
ph_p_saiq <- ph_p_saiq %>% 
  mutate(
    rec_con = rowSums((select(ph_p_saiq, starts_with('ph_p_saiq__activs'))), 
                      na.rm = TRUE),
    rec_con = ifelse(count_na == 29, NA, rec_con),
    rec_bin = ifelse(rec_con == 0, 0, 1),
    rec_bin = factor(rec_bin,
                     levels = c(0, 1),
                     labels = c('No', 'Yes'))) %>% 
  select(participant_id, rec_con, rec_bin)

table(ph_p_saiq$rec_con, useNA = 'ifany') # NA = 1
table(ph_p_saiq$rec_bin, useNA = 'ifany') # NA = 1

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(ph_p_saiq, by = 'participant_id')

rm(ph_p_saiq)

################################################################################
#                             Screen Time Survey (STQ)                         #
################################################################################

nt_p_stq <- read.csv('data/6.0/nt_p_yst.csv')  %>% 
  subset (session_id == 'ses-00A') %>% 
  select(participant_id, nt_p_yst__wkdy__hr_001, nt_p_yst__wknd__hr_001) 

# recode if >= 18 to NA
nt_p_stq <- nt_p_stq %>% 
  mutate(
    nt_p_yst__wkdy__hr_001 = ifelse(
      nt_p_yst__wkdy__hr_001 >= 18, NA, nt_p_yst__wkdy__hr_001),
    nt_p_yst__wknd__hr_001 = ifelse(
      nt_p_yst__wknd__hr_001 >= 18, NA, nt_p_yst__wknd__hr_001))

nt_p_stq <-  nt_p_stq %>% 
  mutate(screentime = nt_p_yst__wkdy__hr_001 + nt_p_yst__wknd__hr_001) %>% 
  select(participant_id, screentime)

table(nt_p_stq$screentime, useNA = 'ifany')

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(nt_p_stq, by = 'participant_id')

rm(nt_p_stq)

################################################################################
#         Youth Risk Behavior Survey (YRB) - Exercise Physical Activity        #
################################################################################

ph_y_yrb <- ph %>% 
  select(participant_id, session_id, starts_with('ph_y_pa')) %>% 
  subset (session_id == 'ses-00A') %>%
  select(participant_id, ph_y_pa_001, ph_y_pa_002, 
         ph_y_pa_003)%>% 
  mutate(act1 = ph_y_pa_001,
         # recode physical_activity2_y to be on the same scale as question 1
         act2 = ph_y_pa_002 - 1, 
         act5 = ph_y_pa_003) %>% 
  select(participant_id, act1, act2, act5)

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(ph_y_yrb, by = 'participant_id')

rm(ph_y_yrb)

################################################################################
#                 Developmental History Questionnaire                          #
################################################################################

ph_p_dhx <- ph %>% 
  select(participant_id, session_id, starts_with('ph_p_dhx')) %>% 
  subset (session_id == 'ses-00A') %>% 
  select(participant_id, ends_with('a')) %>%
  select(-contains('rx'))

# re-code values of 999 to NA for non-caffeine variables
ph_p_dhx <- ph_p_dhx %>% 
  mutate(
    across(.cols = starts_with('ph_p_dhx'),  
           .fns = ~replace(.x, .x %in% c(999), NA))) 


# create binary substance exposure variable 
ph_p_dhx <- ph_p_dhx %>% 
  mutate(exp_sub = rowSums((select(
    ph_p_dhx, starts_with('ph_p_dhx') )),
    na.rm = TRUE),
    count_na = rowSums((is.na(select(
      ph_p_dhx, starts_with('ph_p_dhx') ))))) %>%
  mutate(
    exp_sub = ifelse(exp_sub == 0, 0, 1),
    exp_sub = ifelse(count_na == 10, NA, exp_sub),   
    exp_sub = factor(exp_sub,
                     levels = c(0, 1),
                     labels = c('No', 'Yes')))

table(ph_p_dhx$exp_sub, useNA = 'ifany')# NA = 170

# create caffeine use variable
ph_p_dhx <- ph_p_dhx %>% 
  left_join(ph[ph$session_id == 'ses-00A',] %>% 
              select(participant_id, ph_p_dhx__caff_001), by = 'participant_id') %>%
  mutate(
    exp_caf = ifelse(
      (ph_p_dhx__caff_001 == -1 | ph_p_dhx__caff_001 == 999), NA, 
      ph_p_dhx__caff_001)) %>% 
  mutate(
    exp_caf_rec = case_when(
      exp_caf == 0 ~ '0', # no caffeine
      exp_caf == 3 ~ '1', # less than once a week
      exp_caf == 2 ~ '2', # not daily but more than 1x/week
      exp_caf == 1 ~ '3')) %>%  # daily
  mutate(exp_caf_rec = as.numeric(exp_caf_rec))  

# merge w/full dataset 
ph_p_dhx <- ph_p_dhx %>% 
  select(participant_id, exp_sub, exp_caf_rec)

data.2 <- data.2 %>% 
  left_join(ph_p_dhx, by = 'participant_id')

rm(ph_p_dhx)

################################################################################
#                            Sleep Disorder Scale                              #
################################################################################

ph_p_sds <- ph %>% 
  select(participant_id, session_id, starts_with('ph_p_sds')) %>% 
  subset (session_id == 'ses-00A') %>%
  mutate(
    across(
      # select all columns that match "__{dimension}_"
      .cols = matches("__da_|__dims_|__does_|__hyphy_|__sbd_|__swtd_"),
      .fns = as.numeric # ensure numeric just in case
    )
  ) %>%
  rowwise() %>%
  mutate(
    ph_p_sds__da_sum     = sum(c_across(starts_with("ph_p_sds__da_")), na.rm = TRUE),
    ph_p_sds__dims_sum   = sum(c_across(starts_with("ph_p_sds__dims_")), na.rm = TRUE),
    ph_p_sds__does_sum   = sum(c_across(starts_with("ph_p_sds__does_")), na.rm = TRUE),
    ph_p_sds__hyphy_sum  = sum(c_across(starts_with("ph_p_sds__hyphy_")), na.rm = TRUE),
    ph_p_sds__sbd_sum    = sum(c_across(starts_with("ph_p_sds__sbd_")), na.rm = TRUE),
    ph_p_sds__swtd_sum   = sum(c_across(starts_with("ph_p_sds__swtd_")), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(participant_id, ends_with('sum'))

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(ph_p_sds, by = 'participant_id')

rm(ph_p_sds)

################################################################################
#                         Puberty Development Scale                            #
################################################################################

ph_p_pds <- ph %>% 
  select(participant_id, session_id, starts_with('ph_p_pds')) %>% 
  subset (session_id == 'ses-00A') %>% 
  select(participant_id, ph_p_pds__f_categ, ph_p_pds__m_categ)

# create z-scores for each sex
df_sex <- data.2 %>% 
  select(participant_id, sex_2l) 

ph_p_pds <- ph_p_pds %>% 
  left_join(df_sex, by = 'participant_id') %>% 
  mutate(pds_female_z = 
           ifelse(sex_2l == 'Female', scale(ph_p_pds__f_categ), NA),
         pds_male_z = 
           ifelse(sex_2l == 'Male', scale(ph_p_pds__m_categ), NA)) 

# combine z-scores for each sex into a single variable
ph_p_pds$pds <- rowSums(ph_p_pds[5:6], na.rm = TRUE)

# if sex is missing PDS, recode to NA
ph_p_pds <- ph_p_pds %>% 
  mutate(
    pds = ifelse(
      (sex_2l == 'Female' & is.na(ph_p_pds__f_categ) |
         sex_2l == 'Male' & is.na(ph_p_pds__m_categ) ), NA, pds))


# merge w/full dataset 
ph_p_pds <- ph_p_pds %>% 
  select(participant_id, pds)

data.2 <- data.2 %>% 
  left_join(ph_p_pds, by = 'participant_id')

rm(ph_p_pds, df_sex)

################################################################################
#                  Ohio State Traumatic Brain Injury                           #
################################################################################

ph_p_otbi <- ph %>% 
  select(participant_id, session_id, starts_with('ph_p_otbi')) %>% 
  subset (session_id == 'ses-00A') %>% 
  select(participant_id, ph_p_otbi_001, ph_p_otbi_002, ph_p_otbi_003, ph_p_otbi_004) 

ph_p_otbi <- ph_p_otbi %>% 
  mutate(
    tbi_injury = ifelse(
      (ph_p_otbi_001 == 1 | ph_p_otbi_002 == 1 | ph_p_otbi_003 == 1 | ph_p_otbi_004 == 1), 'Yes', 'No'),
    tbi_injury = factor(tbi_injury)) %>% 
  select(participant_id, tbi_injury) 

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(ph_p_otbi, by = 'participant_id')

rm(ph_p_otbi)

## Independent Variables: Culture & Environment  --------------------------------

################################################################################
#                        Acculturation Survey  (ACC)                           #
################################################################################

ce_y_acc <- fce %>%
  select(participant_id, session_id, starts_with('fc_y_aclt')) %>%
  subset (session_id == 'ses-00A')

# recode (999) Don't Know and (777) Refuse to answer as NA
ce_y_acc <- ce_y_acc %>% 
  select(participant_id, fc_y_aclt_002) %>% 
  mutate(
    across(.cols = fc_y_aclt_002,  
           .fns = ~replace(.x, .x %in% c(999, 777), NA))) 

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(ce_y_acc, by = 'participant_id')

rm(ce_y_acc)

################################################################################
#                   Neighborhood Safety/Crime Survey NSC                       #
################################################################################

ce_y_nsc <- fce %>%
  select(participant_id, session_id, starts_with('fc_y_nsc')) %>% 
  subset (session_id == 'ses-00A') %>% 
  select(participant_id, fc_y_nsc__ns_003) 

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(ce_y_nsc, by = 'participant_id')

rm(ce_y_nsc)

################################################################################
#             Diagnostic Interview for DSM-5 Background Items                  #
################################################################################

mh_p_ksads_bg <- read.csv('data/6.0/mh_p_kbi.csv') %>% 
  subset (session_id == 'ses-00A') %>% 
  select(participant_id, mh_p_kbi__school_006, mh_p_kbi__school_003,
         mh_p_kbi__school_007___1:mh_p_kbi__school_007___10) 

# recode variables as needed
mh_p_ksads_bg <- mh_p_ksads_bg %>% 
  mutate(
    across(.cols = mh_p_kbi__school_006,  
           .fns = ~replace(.x, .x %in% c(6, -1), NA)), 
    across(.cols = mh_p_kbi__school_003,  
           .fns = ~replace(.x, .x %in% c(777, 3), NA))) %>% 
  mutate(mh_p_kbi__school_006 = 
           factor(mh_p_kbi__school_006,
                  levels = c(1, 2, 3, 4, 5),
                  labels = c('Grade_A', 'Grade_B', 'Grade_C', 
                             'Grade_D', 'Grade_F'))) %>% 
  mutate(
    mh_p_kbi__school_003 = case_when(
      mh_p_kbi__school_003 == 1 ~ 'Yes',
      mh_p_kbi__school_003 == 0 ~ 'No'), 
    mh_p_kbi__school_003 = as.factor(mh_p_kbi__school_003))

# check if multiple special education categories are selected 
mh_p_ksads_bg$se_multi <- rowSums(
  (select(mh_p_ksads_bg, starts_with('mh_p_kbi__school_007'))), na.rm = TRUE)

table(mh_p_ksads_bg$se_multi, useNA = 'ifany')

# create special education services groups
mh_p_ksads_bg <- mh_p_ksads_bg %>% 
  mutate(
    se_services = case_when(
      se_multi == 0 ~ as.character(NA),
      se_multi == 1 & 
        (mh_p_kbi__school_007___1 == 1 | mh_p_kbi__school_007___2 == 1 |
           mh_p_kbi__school_007___3 == 1 | mh_p_kbi__school_007___4 == 1 |
           mh_p_kbi__school_007___5 == 1 | mh_p_kbi__school_007___6 == 1 |
           mh_p_kbi__school_007___7 == 1 ) ~ 'Emotion_or_Learning_Support',
      se_multi == 1 & mh_p_kbi__school_007___8 == 1 ~ 'Gifted', 
      se_multi == 1 & mh_p_kbi__school_007___9 == 1 ~ 'Other', 
      se_multi == 1 & mh_p_kbi__school_007___10 == 1 ~ 'None', 
      se_multi >= 2 & (mh_p_kbi__school_007___1 == 1 | 
                         mh_p_kbi__school_007___2 == 1 |
                         mh_p_kbi__school_007___3 == 1 | mh_p_kbi__school_007___4 == 1 |
                         mh_p_kbi__school_007___5 == 1 | mh_p_kbi__school_007___6 == 1 |
                         mh_p_kbi__school_007___7 == 1) & (mh_p_kbi__school_007___8 == 0 & 
                                                             mh_p_kbi__school_007___9 == 0 & mh_p_kbi__school_007___10 == 0) 
      ~ 'Emotion_or_Learning_Support',
      se_multi >= 2 & (mh_p_kbi__school_007___1 == 1 | 
                         mh_p_kbi__school_007___2 == 1 |
                         mh_p_kbi__school_007___3 == 1 | mh_p_kbi__school_007___4 == 1 |
                         mh_p_kbi__school_007___5 == 1 | mh_p_kbi__school_007___6 == 1 |
                         mh_p_kbi__school_007___7 == 1) & 
        (mh_p_kbi__school_007___8 == 1 & mh_p_kbi__school_007___9 == 0) & 
        (mh_p_kbi__school_007___10 == 0) 
      ~ 'Combined Services_1', 
      se_multi >= 2 & (mh_p_kbi__school_007___1 == 1 | 
                         mh_p_kbi__school_007___2 == 1 |
                         mh_p_kbi__school_007___3 == 1 | mh_p_kbi__school_007___4 == 1 |
                         mh_p_kbi__school_007___5 == 1 | mh_p_kbi__school_007___6 == 1 |
                         mh_p_kbi__school_007___7 == 1) & 
        (mh_p_kbi__school_007___8 == 0 & mh_p_kbi__school_007___9 == 1) & 
        (mh_p_kbi__school_007___10 == 0) 
      ~ 'Combined Services_2', 
      se_multi >= 2 & (mh_p_kbi__school_007___1 == 0 & 
                         mh_p_kbi__school_007___2 == 0 &
                         mh_p_kbi__school_007___3 == 0 & mh_p_kbi__school_007___4 == 0 &
                         mh_p_kbi__school_007___5 == 0 & mh_p_kbi__school_007___6 == 0 &
                         mh_p_kbi__school_007___7 == 0) & 
        (mh_p_kbi__school_007___8 == 1 & mh_p_kbi__school_007___9 == 1) & 
        (mh_p_kbi__school_007___10 == 0) 
      ~ 'Combined Services_3', 
      se_multi >= 2 & (mh_p_kbi__school_007___1 == 1 | 
                         mh_p_kbi__school_007___2 == 1 |
                         mh_p_kbi__school_007___3 == 1 | mh_p_kbi__school_007___4 == 1 |
                         mh_p_kbi__school_007___5 == 1 | mh_p_kbi__school_007___6 == 1 |
                         mh_p_kbi__school_007___7 == 1) & 
        (mh_p_kbi__school_007___8 == 1 & mh_p_kbi__school_007___9 == 1) & 
        (mh_p_kbi__school_007___10 == 0)
      ~ 'Combined Services_4',
      se_multi >= 2 & mh_p_kbi__school_007___10 == 1 ~ 'None'))  

table(mh_p_ksads_bg$se_services, useNA = 'ifany')
# Combined services_1 = Emotion or learning support + Gifted 
# Combined services_2 = Emotion or learning support + Other
# Combined services_3 = Gifted + Other
# Combined services_4 = Emotion or learning support + Gifted + Other  

# breakdown of those who indicated multiple groups
table(mh_p_ksads_bg$se_multi, mh_p_ksads_bg$se_services, useNA = 'ifany')

# combine combined services groups into larger category
mh_p_ksads_bg <- mh_p_ksads_bg %>% 
  mutate(
    se_services = if_else(
      se_services == 'Combined Services_1' | 
        se_services == 'Combined Services_2' |
        se_services == 'Combined Services_3' |
        se_services == 'Combined Services_4', 'Combined_Services', se_services))
table(mh_p_ksads_bg$se_services, useNA = 'ifany') # (NA) = 101



# merge w/full dataset 
mh_p_ksads_bg <- mh_p_ksads_bg %>%  
  rename(det_susp = mh_p_kbi__school_003) %>% 
  select(participant_id, mh_p_kbi__school_006, det_susp, se_services) %>% 
  mutate(se_services = factor(se_services))

data.2 <- data.2 %>% 
  left_join(mh_p_ksads_bg, by = 'participant_id') 

# create summary score for failing grades (D and F)
data.2 <- data.2 %>% 
  mutate(mh_p_kbi__school_006 = case_when(
    mh_p_kbi__school_006 == 'Grade_A' ~ 'Grade_A',
    mh_p_kbi__school_006 == 'Grade_B' ~ 'Grade_B',
    mh_p_kbi__school_006 == 'Grade_C' ~ 'Grade_C',
    mh_p_kbi__school_006 == 'Grade_D' ~ 'Grade_Fail',
    mh_p_kbi__school_006 == 'Grade_F' ~ 'Grade_Fail')) %>% 
  mutate(mh_p_kbi__school_006 = as.factor(mh_p_kbi__school_006))
class(data.2$mh_p_kbi__school_006) # factor

rm(mh_p_ksads_bg)

################################################################################
#                                  SRPF                                        #
################################################################################

ce_y_srpf <- fce %>%
  select(participant_id, session_id, starts_with('fc_y_srpf')) %>% 
  subset (session_id == 'ses-00A') %>%
  select(participant_id, fc_y_srpf__env_mean, fc_y_srpf__involv_mean, fc_y_srpf__dis_mean) 

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(ce_y_srpf, by = 'participant_id') 

rm(ce_y_srpf)

## Independent Variables: Biospecimens  -----------------------------------------

################################################################################
#                   Hormone Saliva Salimetric Scores                           #
################################################################################

ph_y_sal_horm <- read.csv('data/6.0/ph_y_phs.csv') %>% 
  subset (session_id == 'ses-00A') %>% 
  select(participant_id, ph_y_phs__ert_mean, ph_y_phs__dhea_mean)

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(ph_y_sal_horm, by = 'participant_id') 

rm(ph_y_sal_horm)

## Independent Variables: Neurocognitive Factors  -------------------------------
nc <- read.csv('data/6.0/nc_y_cct.csv')

################################################################################
#                                   RAVLT                                      #
################################################################################

nc_y_ravlt <- nc %>%
  select(participant_id, session_id, starts_with('nc_y_ravlt')) %>%
  subset (session_id == 'ses-00A') %>% 
  select(participant_id, nc_y_ravlt__trial1__crct_count, 
         nc_y_ravlt__trial5__crct_count,  nc_y_ravlt__trial6__sd__crct_count,
         nc_y_ravlt__trial7__ld__crct_count) 

nc_y_ravlt <- nc_y_ravlt %>% 
  mutate(pea_ravlt_learn = 
           nc_y_ravlt__trial5__crct_count - nc_y_ravlt__trial1__crct_count) %>%  
  select(participant_id, pea_ravlt_learn, nc_y_ravlt__trial6__sd__crct_count, 
         nc_y_ravlt__trial7__ld__crct_count)

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(nc_y_ravlt, by = 'participant_id') 

rm(nc_y_ravlt)

################################################################################
#                                  WISC-V                                      #
################################################################################

nc_y_wisc <- nc %>%
  select(participant_id, session_id, starts_with('nc_y_wisc')) %>% 
  subset (session_id == 'ses-00A') %>% 
  select(participant_id, nc_y_wisc__scaled_score)

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(nc_y_wisc, by = 'participant_id') 

rm(nc_y_wisc)

################################################################################
#                         Cash Choice Task                                     #
################################################################################

nc_y_cct <- nc %>%
  select(participant_id, session_id, starts_with('nc_y_cct')) %>% 
  subset (session_id == 'ses-00A') %>% 
  select(participant_id, nc_y_cct_001)

# recode values of (3) Don't Know to NA
nc_y_cct <- nc_y_cct %>% 
  mutate(
    across(.cols = nc_y_cct_001,  
           .fns = ~replace(.x, .x == 3, NA))) %>% 
  mutate(
    cct = factor(nc_y_cct_001, levels = c(1,2), 
                 labels = c('Immediate', 'Delayed'))) %>% 
  select(participant_id, cct)

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(nc_y_cct, by = 'participant_id')

rm(nc_y_cct)

################################################################################
#                           Little Man Task                                    #
################################################################################

nc_y_lmt <- nc %>%
  select(participant_id, session_id, starts_with('nc_y_lmt')) %>% 
  subset (session_id == 'ses-00A') %>% 
  select(participant_id, nc_y_lmt__crct_count, nc_y_lmt__wrong_count, 
         nc_y_lmt__crct_acc, nc_y_lmt__crct_rt, nc_y_lmt_effncy)

nc_y_lmt <- nc_y_lmt %>% 
  mutate(
    lmt_acc = (nc_y_lmt__crct_count)/(nc_y_lmt__crct_count + nc_y_lmt__wrong_count),
    lmt_acc = ifelse(
      (nc_y_lmt__crct_count == 0 & nc_y_lmt__wrong_count == 0), NA, lmt_acc)) %>% 
  select(participant_id, lmt_acc, nc_y_lmt__crct_rt, nc_y_lmt_effncy)

# merge w/full dataset  
data.2 <- data.2 %>% 
  left_join(nc_y_lmt, by = 'participant_id')

rm(nc_y_lmt)

################################################################################
#                                NIH ToolBox                                   #
################################################################################

nc_y_nihtb <- nc %>%
  select(participant_id, session_id, starts_with('nc_y_nihtb')) %>% 
  subset (session_id == 'ses-00A') %>%
  select(participant_id, nc_y_nihtb__picvcb__agecor_score, 
         nc_y_nihtb__flnkr__agecor_score, nc_y_nihtb__lswmt__agecor_score,
         nc_y_nihtb__crdst__agecorr_score, nc_y_nihtb__pttcp__agecor_score, 
         nc_y_nihtb__picsq__agecor_score, nc_y_nihtb__readr__agecor_score)

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(nc_y_nihtb, by = 'participant_id')

rm(nc_y_nihtb)


## MRI Data ---------------------------------------------------------------------
mri <- read.csv('data/6.0/mr_y_dti__fs__fa__aseg.csv')

################################################################################
#                                DTI                                  #
################################################################################
dti <- mri %>%
  select(participant_id, session_id, starts_with('mr_y_dti')) %>%
  subset(session_id == 'ses-00A') %>%
  select(-session_id, -contains('_mean'))   # Keep only weighted means

data.2 <- data.2 %>% 
  left_join(dti, by = 'participant_id')

rm(dti)

################################################################################
#                                SMRI                                  #
################################################################################
smri <- read.csv('data/6.0/mr_y_smri__sulc__dsk.csv') %>%
  select(participant_id, session_id, starts_with('mr_y_smri')) %>%
  subset(session_id == 'ses-00A') %>%
  select(-session_id)

data.2 <- data.2 %>% 
  left_join(smri, by = 'participant_id')

rm(smri)

################################################################################
#                                rsfMRI                                  #
################################################################################
rsfmri <- mri %>%
  select(participant_id, session_id, starts_with('mr_y_rsfmri')) %>%
  subset(session_id == 'ses-00A') %>%
  select(-session_id)

cort <- rsfmri %>%
  select(-contains('aseg'))

subcort <- rsfmri %>%
  select(participant_id, contains('aseg'))

# Get all rsfmri variable names
rsfmri_vars <- names(cort) %>% 
  str_subset("^mr_y_rsfmri__corr__gpnet__")

# Extract network pairs
pair_df <- tibble(var = rsfmri_vars) %>%
  mutate(
    net1 = str_match(var, "__gpnet__(.*?)__")[,2],
    net2 = str_match(var, "__gpnet__.*?__(.*?)_mean")[,2]
  ) %>%
  # Sort the two names alphabetically within each pair
  rowwise() %>%
  mutate(pair_id = paste(sort(c(net1, net2)), collapse = "__")) %>%
  ungroup() %>%
  # Keep only one per unique pair
  distinct(pair_id, .keep_all = TRUE)

# These are  unique variable names
unique_rsfmri_vars <- pair_df$var

# Optional: subset the dataframe to just those
cort <- cort %>% select(participant_id, all_of(unique_rsfmri_vars))

data.2 <- data.2 %>% 
  left_join(cort, by = 'participant_id') %>%
  left_join(subcort, by = 'participant_id')

rm(rsfmri, cort, subcort)

# Sample Characteristics: Numeric Variables ------------------------------------

# drop objects no longer needed
data.2 <- data.2 %>% filter(!is.na(ng4_class))

rm(dem, fce, mh, nc, ph, su, sui, mri)

################################################################################
#                      create table & domain names 
################################################################################

# create table and domain names
temp <- data.2 %>% 
  select_if(is.factor)
factors <- names(temp)
rm(temp)

dv_names <- c('alc1', 'alc2', 'alc3',
              'alc4', 'alc5', 'alc6', 
              'alc7', 'alc8', 'alc9', 'alc10',
              'nic1', 'nic2', 'nic3',
              'nic4', 'nic5', 'nic6', 
              'nic7', 'nic8', 'nic9', 'nic10',
              'mj1', 'mj2', 'mj3',
              'mj4', 'mj5', 'mj6', 
              'mj7', 'mj8', 'mj9', 'mj10',
              'oth1', 'oth2', 'oth3',
              'oth4', 'oth5', 'oth6', 
              'oth7', 'oth8', 'oth9', 'oth10')

table_names_numeric <- data.2 %>% 
  select(-participant_id, -all_of(factors)) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'variable') %>% 
  select(variable) %>%
  mutate(
    table_name = 
      case_when(
        # Self and Peer Involvement with Substance Use 
        variable == 'peer_alc' ~ 'Peer Substance Use - Alcohol',
        variable == 'peer_tob' ~ 'Peer Substance Use - Nicotine',
        variable == 'peer_cb' ~ 'Peer Substance Use - Cannabis',
        variable == 'peer_other' ~ 'Peer Substance Use - Other Substance Use',
        variable == 'peer_prob' ~ 'Peer Substance Use - Problems with Substance Use',
        variable == 'path_alc' ~ 'Intent to Use - Alcohol',
        variable == 'path_tob' ~ 'Intent to Use - Nicotine',
        variable == 'path_cb' ~ 'Intent to Use - Cannabis',
        
        # Parenting Behaviors
        variable == 'crpf' ~ 'CRPF',
        variable == 'par_rules' ~ 'Parent Rules',
        variable == 'fc_y_pm_mean' ~ 'PMQ - Parental Monitoring',
        variable == 'fc_y_fes__confl_mean' ~ 'FES - Conflict Subscale',
        variable == 'crpbi_y_ss_parent' ~ 'CRPBI - Acceptance Subscale',
        
        # Demographics
        variable == 'age_baseline' ~ 'Age at Baseline',
        
        # Mental Health
        variable == 'mh_p_cbcl__synd__anxdep_tscore' ~ 'CBCL - Anxiety / Depression',
        variable == 'mh_p_cbcl__synd__wthdep_tscore' ~ 'CBCL - Withdraw / Depression',
        variable == 'mh_p_cbcl__synd__som_tscore' ~ 'CBCL - Somatic Symptoms',
        variable == 'mh_p_cbcl__synd__soc_tscore' ~ 'CBCL - Social Problems',
        variable == 'mh_p_cbcl__synd__tho_tscore' ~ 'CBCL - Thought Problems',
        variable == 'mh_p_cbcl__synd__attn_tscore' ~ 'CBCL - Attention Problems',
        variable == 'mh_p_cbcl__synd__rule_tscore' ~ 'CBCL - Rule-breaking Behavior',
        variable == 'mh_p_cbcl__synd__aggr_tscore' ~ 'CBCL - Aggressive Behavior',
        variable == 'mh_p_cbcl__synd__int_tscore' ~ 'CBCL - Internalizing Disorders',
        variable == 'mh_p_cbcl__synd__ext_tscore' ~ 'CBCL - Externalizing Disorders',
        variable == 'mh_density' ~ 'Family History of Substance Use',
        variable == 'mh_y_pps__severity_score' ~ 'PPS - Positive Symptoms Endorsed',
        variable == 'mh_y_upps__nurg_sum' ~ 'UPPS-P - Negative Urgency',
        variable == 'mh_y_upps__pers_sum' ~ 'UPPS-P - Lack of Perseverance',
        variable == 'mh_y_upps__plan_sum' ~ 'UPPS-P - Lack of Planning',
        variable == 'mh_y_upps__sens_sum' ~ 'UPPS-P - Sensation Seeking',
        variable == 'mh_y_upps__purg_sum' ~ 'UPPS-P - Positive Urgency',
        variable == 'mh_y_bisbas__bis_sum' ~ 'BIS/BAS - BIS Sum',
        variable == 'mh_y_bisbas__bas__dr_sum' ~ 'BIS/BAS - Drive',
        variable == 'mh_y_bisbas__bas__fs_sum' ~ 'BIS/BAS - Fun Seeking',
        variable == 'mh_y_bisbas__bas__rr_sum' ~ 'BIS/BAS - Reward Responsiveness',
        
        # Physical Health
        variable == 'rec_con' ~ 'Recreational Activities (continuous)',
        variable == 'screentime' ~ 'Screentime',
        variable == 'act1' ~ 'Physical Activity - Days of 60-min.',
        variable == 'act2' ~ 'Physical Activity - Days of Strength Building',
        variable == 'act5' ~ 'Physical Activity - Days of P.E. Course',
        variable == 'exp_caf_rec' ~ 'Prenatal Exposure - Caffeine',
        variable == 'ph_p_sds__dims_sum' ~ 'SDS - Disorders of Initiating and Maintaining Sleep',
        variable == 'ph_p_sds__sbd_sum' ~ 'SDS - Sleep Breathing disorders',
        variable == 'ph_p_sds__da_sum' ~ 'SDS - Disorder of Arousal',
        variable == 'ph_p_sds__swtd_sum' ~ 'SDS - Sleep-Wake transition Disorders',
        variable == 'ph_p_sds__does_sum' ~ 'SDS - Disorders of Excessive Somnolence',
        variable == 'ph_p_sds__hyphy_sum' ~ 'SDS - Sleep Hyperhydrosis',
        variable == 'pds' ~ 'PDS - Puberty Development',
        
        # Culture & Environment
        variable == 'fc_y_aclt_002' ~ 'Bilingual Status',
        variable == 'fc_y_nsc__ns_003' ~ 'Neighborhood Safety',
        variable == 'fc_y_srpf__env_mean' ~ 'SRPF - School Environment Subscale',
        variable == 'fc_y_srpf__involv_mean' ~ 'SRPF - School Involvement Subscale',
        variable == 'fc_y_srpf__dis_mean' ~ 'SRPF - School Disengagement Subscale',
        
        # Biospecimens
        variable == 'ph_y_phs__ert_mean' ~ 'Hormones - Testosterone',
        variable == 'ph_y_phs__dhea_mean' ~ 'Hormones - DHEA',
        
        # Neurocog
        variable == 'pea_ravlt_learn' ~ 'RAVLT - Learning Score',
        variable == 'nc_y_ravlt__trial6__sd__crct_count' ~ 'RAVLT - Short Delay',
        variable == 'nc_y_ravlt__trial7__ld__crct_count' ~ 'RAVLT - Long Delay',
        variable == 'nc_y_wisc__scaled_score' ~ 'WISC-V - Matrix Reasoning',
        variable == 'lmt_acc' ~ 'LMT - Accuracy',
        variable == 'nc_y_lmt__crct_rt' ~ 'LMT - Reaction Time for Correct Trials',
        variable == 'nc_y_lmt_effncy' ~ 'LMT - Efficiency',
        variable == 'nc_y_nihtb__picvcb__agecor_score' ~ 'NIH Toolbox - Picture Vocabulary Test',
        variable == 'nc_y_nihtb__flnkr__agecor_score' ~ 'NIH Toolbox - Flanker Inhibitory Control and Attention Test',
        variable == 'nc_y_nihtb__lswmt__agecor_score' ~ 'NIH Toolbox - List Sorting Working Memory Test',
        variable == 'nc_y_nihtb__crdst__agecorr_score' ~ 'NIH Toolbox - Dimensional Change Card Sort Test',
        variable == 'nc_y_nihtb__pttcp__agecor_score' ~ 'NIH Toolbox - Pattern Comparison Processing Speed Test',
        variable == 'nc_y_nihtb__picsq__agecor_score' ~ 'NIH Toolbox - Picture Sequence Memory Test',
        variable == 'nc_y_nihtb__readr__agecor_score' ~ 'NIH Toolbox - Oral Reading Recognition Test')) 

domain_names_numeric <- data.2 %>% 
  select(-participant_id, -all_of(factors)) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "variable") %>% 
  select(variable) %>% 
  mutate(
    domain_name = case_when(
      
      # Self and Peer Involvement with Substance Use
      substr(variable, 1, 4) == 'peer' ~ 'Self and Peer Involvement with Substance Use',
      substr(variable, 1, 4) == 'path' ~ 'Self and Peer Involvement with Substance Use',
      
      # Parenting Behaviors
      variable == 'crpf' ~ 'Parenting Behaviors',
      variable == 'par_rules' ~ 'Parenting Behaviors',
      variable == 'fc_y_pm_mean' ~ 'Parenting Behaviors',
      variable == 'fc_y_fes__confl_mean' ~ 'Parenting Behaviors',
      variable == 'crpbi_y_ss_parent' ~ 'Parenting Behaviors',
      
      # Demographics
      variable == 'age_baseline' ~ 'Demographics',
      
      # Mental Health
      substr(variable, 1, 9) == 'mh_p_cbcl' ~ 'Mental Health',
      variable == 'mh_density' ~ 'Mental Health',
      variable == 'mh_y_pps__severity_score' ~ 'Mental Health',
      substr(variable, 1, 9) == 'mh_y_upps' ~ 'Mental Health',
      substr(variable, 1, 8) == 'mh_y_bis' ~ 'Mental Health',
      
      # Physical Health
      variable == 'rec_con' ~ 'Physical Health',
      variable == 'screentime' ~ 'Physical Health',
      variable == 'act1' ~ 'Physical Health',
      variable == 'act2' ~ 'Physical Health',
      variable == 'act5' ~ 'Physical Health',
      variable == 'exp_caf_rec' ~ 'Physical Health',
      substr(variable, 1, 8) == 'ph_p_sds' ~ 'Physical Health',
      variable == 'pds' ~ 'Physical Health',
      
      # Culture & Environment
      variable == 'fc_y_aclt_002' ~ 'Culture & Environment',
      variable == 'fc_y_nsc__ns_003' ~ 'Culture & Environment',
      substr(variable, 1, 9) == 'fc_y_srpf' ~ 'Culture & Environment',
      
      # Biospecimens
      substr(variable, 1, 8) == 'ph_y_phs' ~ 'Hormones',
      
      # Neurocog
      substr(variable, 1, 5) == 'nc_y_' ~ 'Neurocognitive Factors',
      variable == 'lmt_acc' | variable == 'pea_ravlt_learn' ~ 'Neurocognitive Factors'))

table_names_numeric <- table_names_numeric %>% 
  full_join(domain_names_numeric, by = 'variable') 
rm(domain_names_numeric)

# select factor variables
factor_vars <- data.2 %>% select(-ng4_class, -ng3_class) %>%
  mutate(across(where(is.factor), ~ fct_na_value_to_level(., "Missing"))) %>% select_if(is.factor)

# create factor dictionary by expanding levels
table_names_factors <- factor_vars %>%
  summarise(across(everything(), ~list(levels(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "levels") %>%
  unnest(levels) %>%
  mutate(
    variable_level = paste0(variable, '_', levels)) %>%
  select(-variable, -levels) %>%
  rename(variable = variable_level) %>%
  mutate(
    table_name = case_when(
      variable == 'alc0_Yes' ~ 'Baseline Use - Alcohol',
      variable == 'alc0_No' ~ 'No Baseline Use - Alcohol',
      variable == 'mj0_Yes' ~ 'Baseline Use - Marijuana',
      variable == 'mj0_No' ~ 'No Baseline Use - Marijuana',
      variable == 'nic0_Yes' ~ 'Baseline Use - Nicotine',
      variable == 'nic0_No' ~ 'No Baseline Use - Nicotine',
      variable == 'oth0_Yes' ~ 'Baseline Use - Other',
      variable == 'oth0_No' ~ 'No Baseline Use - Other',
      variable == 'sex_2l_Male' ~ 'Sex - Male',
      variable == 'sex_2l_Female' ~ 'Sex - Female',
      variable == 'sex_2l_Missing' ~ 'Sex - No Response',
      variable == 'race_4l_Asian' ~ 'Race - Asian',
      variable == 'race_4l_Black' ~ 'Race - Black',
      variable == 'race_4l_Other_MultiRacial' ~ 'Race - Other/Multi-Racial',
      variable == 'race_4l_White' ~ 'Race - White',
      variable == 'race_4l_Missing' ~ 'Race - No Response',
      variable == 'eth_hisp_non_Hispanic' ~ 'Ethnicity - Non-Hispanic',
      variable == 'eth_hisp_Hispanic' ~ 'Ethnicity - Hispanic',
      variable == 'eth_hisp_Missing' ~ 'Ethnicity - No Response',
      variable == 'income_inc_1' ~ 'Primary Caregiver Income - Up to $24,999',
      variable == 'income_inc_2' ~ 'Primary Caregiver Income - Between $25,000 and $49,999',
      variable == 'income_inc_3' ~ 'Primary Caregiver Income - Between $50,000 and $99,999',
      variable == 'income_inc_4' ~ 'Primary Caregiver Income - Greater than $100,000',
      variable == 'income_Missing' ~ 'Primary Caregiver Income - No Response',
      variable == 'religion_rp_1' ~ 'Religious Affiliation - Mainline Protestant',
      variable == 'religion_rp_2' ~ 'Religious Affiliation - Evangelical Protestant',
      variable == 'religion_rp_3' ~ 'Religious Affiliation - Historically Black Church',
      variable == 'religion_rp_4' ~ 'Religious Affiliation - Roman Catholic',
      variable == 'religion_rp_5' ~ 'Religious Affiliation - Jewish (Judaism)',
      variable == 'religion_rp_6' ~ 'Religious Affiliation - LDS Church (Mormon)',
      variable == 'religion_rp_7' ~ "Religious Affiliation - Jehovah's Witness",
      variable == 'religion_rp_8' ~ 'Religious Affiliation - Muslim (Islam)',
      variable == 'religion_rp_9' ~ 'Religious Affiliation - Buddhist',
      variable == 'religion_rp_10' ~ 'Religious Affiliation - Hindu',
      variable == 'religion_rp_11' ~ 'Religious Affiliation - Orthodox Christian',
      variable == 'religion_rp_12' ~ 'Religious Affiliation - Unitarian (Universalist)',
      variable == 'religion_rp_13' ~ 'Religious Affiliation - Other Christian',
      variable == 'religion_rp_14' ~ 'Religious Affiliation - Atheist',
      variable == 'religion_rp_15' ~ 'Religious Affiliation - Agnostic',
      variable == 'religion_rp_16' ~ 'Religious Affiliation - Other',
      variable == 'religion_rp_17' ~ 'Religious Affiliation - Nothing in Particular',
      variable == 'religion_Missing' ~ 'Religious Affiliation - No Response',
      variable == 'p_edu_Bachelors_Degree' ~ "Primary Caregiver's Maximum Education - Bachelor's Degree",
      variable == 'p_edu_HS_Graduate_GED_Equivalent' ~ "Primary Caregiver's Maximum Education - High School or GED",
      variable == 'p_edu_Less_than_HS_Degree_GED_Equivalent' ~ "Primary Caregiver's Maximum Education - Less than High School or GED",
      variable == 'p_edu_Masters_Degree' ~ "Primary Caregiver's Maximum Education - Master's Degree",
      variable == 'p_edu_Professional_School_or_Doctoral_Degree' ~ "Primary Caregiver's Maximum Education - Doctoral Degree",
      variable == 'p_edu_Some_College_or_Associates_Degree' ~ "Primary Caregiver's Maximum Education - Some College or Associate's Degree",
      variable == 'p_edu_Missing' ~ "Primary Caregiver's Maximum Education - No Response",
      variable == 'rec_bin_No' ~ 'Recreational Activities (binary) - No',
      variable == 'rec_bin_Yes' ~ 'Recreational Activities (binary) - Yes',
      variable == 'rec_bin_Missing' ~ 'Recreational Activities (binary) - No Response',
      variable == 'exp_sub_No' ~ 'Prenatal Exposure - Substance Use - No',
      variable == 'exp_sub_Yes' ~ 'Prenatal Exposure - Substance Use - Yes',
      variable == 'exp_sub_Missing' ~ 'Prenatal Exposure - Substance Use - No Response',
      variable == 'tbi_injury_No' ~ 'TBI Injury - No',
      variable == 'tbi_injury_Yes' ~ 'TBI Injury - Yes',
      variable == 'tbi_injury_Missing' ~ 'TBI Injury - No Response',
      variable == 'mh_p_kbi__school_006_Grade_A' ~ 'Grade - A',
      variable == 'mh_p_kbi__school_006_Grade_B' ~ 'Grade - B',
      variable == 'mh_p_kbi__school_006_Grade_C' ~ 'Grade - C',
      variable == 'mh_p_kbi__school_006_Grade_Fail' ~ 'Grade - Fail',
      variable == 'mh_p_kbi__school_006_Missing' ~ 'Grade - No Response',
      variable == 'det_susp_No' ~ 'Detention / Suspension - No',
      variable == 'det_susp_Yes' ~ 'Detention / Suspension - Yes',
      variable == 'det_susp_Missing' ~ 'Detention / Suspension - No Response',
      variable == 'se_services_Combined_Services' ~ 'Special Education Services - Combined Services',
      variable == 'se_services_Emotion_or_Learning_Support' ~ 'Special Education Services - Emotion or Learning Support',
      variable == 'se_services_Gifted' ~ 'Special Education Services - Gifted',
      variable == 'se_services_None' ~ 'Special Education Services - None',
      variable == 'se_services_Other' ~ 'Special Education Services - Other',
      variable == 'se_services_Missing' ~ 'Special Education Services - No Response',
      variable == 'cct_Immediate' ~ 'Cash Choice Task - Immediate',
      variable == 'cct_Delayed' ~ 'Cash Choice Task - Delayed',
      variable == 'cct_Missing' ~ 'Cash Choice Task - No Response'
    ),
    domain_name = case_when(
      variable == 'alc0_Yes' ~ 'Self and Peer Involvement with Substance Use',
      variable == 'alc0_No' ~ 'Self and Peer Involvement with Substance Use',
      variable == 'mj0_Yes' ~ 'Self and Peer Involvement with Substance Use',
      variable == 'mj0_No' ~ 'Self and Peer Involvement with Substance Use',
      variable == 'nic0_Yes' ~ 'Self and Peer Involvement with Substance Use',
      variable == 'nic0_No' ~ 'Self and Peer Involvement with Substance Use',
      variable == 'oth0_Yes' ~ 'Self and Peer Involvement with Substance Use',
      variable == 'oth0_No' ~ 'Self and Peer Involvement with Substance Use',
      variable == 'sex_2l_Male' ~ 'Demographics',
      variable == 'sex_2l_Female' ~ 'Demographics',
      variable == 'sex_2l_Missing' ~ 'Demographics',
      variable == 'race_4l_Asian' ~ 'Demographics',
      variable == 'race_4l_Black' ~ 'Demographics',
      variable == 'race_4l_Other_MultiRacial' ~ 'Demographics',
      variable == 'race_4l_White' ~ 'Demographics',
      variable == 'race_4l_Missing' ~ 'Demographics',
      variable == 'eth_hisp_non_Hispanic' ~ 'Demographics',
      variable == 'eth_hisp_Hispanic' ~ 'Demographics',
      variable == 'eth_hisp_Missing' ~ 'Demographics',
      variable == 'income_inc_1' ~ 'Demographics',
      variable == 'income_inc_2' ~ 'Demographics',
      variable == 'income_inc_3' ~ 'Demographics',
      variable == 'income_inc_4' ~ 'Demographics',
      variable == 'income_Missing' ~ 'Demographics',
      variable == 'religion_rp_1' ~ 'Demographics',
      variable == 'religion_rp_2' ~ 'Demographics',
      variable == 'religion_rp_3' ~ 'Demographics',
      variable == 'religion_rp_4' ~ 'Demographics',
      variable == 'religion_rp_5' ~ 'Demographics',
      variable == 'religion_rp_6' ~ 'Demographics',
      variable == 'religion_rp_7' ~ 'Demographics',
      variable == 'religion_rp_8' ~ 'Demographics',
      variable == 'religion_rp_9' ~ 'Demographics',
      variable == 'religion_rp_10' ~ 'Demographics',
      variable == 'religion_rp_11' ~ 'Demographics',
      variable == 'religion_rp_12' ~ 'Demographics',
      variable == 'religion_rp_13' ~ 'Demographics',
      variable == 'religion_rp_14' ~ 'Demographics',
      variable == 'religion_rp_15' ~ 'Demographics',
      variable == 'religion_rp_16' ~ 'Demographics',
      variable == 'religion_rp_17' ~ 'Demographics',
      variable == 'religion_Missing' ~ 'Demographics',
      variable == 'p_edu_Bachelors_Degree' ~ "Demographics",
      variable == 'p_edu_HS_Graduate_GED_Equivalent' ~ "Demographics",
      variable == 'p_edu_Less_than_HS_Degree_GED_Equivalent' ~ "Demographics",
      variable == 'p_edu_Masters_Degree' ~ "Demographics",
      variable == 'p_edu_Professional_School_or_Doctoral_Degree' ~ "Demographics",
      variable == 'p_edu_Some_College_or_Associates_Degree' ~ "Demographics",
      variable == 'p_edu_Missing' ~ "Demographics",
      variable == 'rec_bin_No' ~ 'Physical Health',
      variable == 'rec_bin_Yes' ~ 'Physical Health',
      variable == 'rec_bin_Missing' ~ 'Physical Health',
      variable == 'exp_sub_No' ~ 'Physical Health',
      variable == 'exp_sub_Yes' ~ 'Physical Health',
      variable == 'exp_sub_Missing' ~ 'Physical Health',
      variable == 'tbi_injury_No' ~ 'Physical Health',
      variable == 'tbi_injury_Yes' ~ 'Physical Health',
      variable == 'tbi_injury_Missing' ~ 'Physical Health',
      variable == 'mh_p_kbi__school_006_Grade_A' ~ 'Culture & Environment',
      variable == 'mh_p_kbi__school_006_Grade_B' ~ 'Culture & Environment',
      variable == 'mh_p_kbi__school_006_Grade_C' ~ 'Culture & Environment',
      variable == 'mh_p_kbi__school_006_Grade_Fail' ~ 'Culture & Environment',
      variable == 'mh_p_kbi__school_006_Missing' ~ 'Culture & Environment',
      variable == 'det_susp_No' ~ 'Culture & Environment',
      variable == 'det_susp_Yes' ~ 'Culture & Environment',
      variable == 'det_susp_Missing' ~ 'Culture & Environment',
      variable == 'se_services_Combined_Services' ~ 'Culture & Environment',
      variable == 'se_services_Emotion_or_Learning_Support' ~ 'Culture & Environment',
      variable == 'se_services_Gifted' ~ 'Culture & Environment',
      variable == 'se_services_None' ~ 'Culture & Environment',
      variable == 'se_services_Other' ~ 'Culture & Environment',
      variable == 'se_services_Missing' ~ 'Culture & Environment',
      variable == 'cct_Immediate' ~ 'Neurocognitive Factors',
      variable == 'cct_Delayed' ~ 'Neurocognitive Factors',
      variable == 'cct_Missing' ~ 'Neurocognitive Factors')
  ) %>%
  select(variable, table_name, domain_name)

table_names_all <- bind_rows(table_names_numeric, table_names_factors)

mri_names <- table_names_all %>%
  filter(is.na(table_name)) %>%
  select(variable) %>%
  left_join(dict %>%
              select(name, description, domain) %>%
              rename(variable = name), by = 'variable') %>%
  rename(domain_name = domain, table_name = description) %>%
  select(variable, table_name, domain_name)

table_names_all <- table_names_all %>%
  filter(!is.na(table_name))

table_names_all <- bind_rows(table_names_all, mri_names)

table_names_all$table_name <- gsub(" \\(full shell DTI\\) in AtlasTrack fiber tract", "",
                                   table_names_all$table_name, fixed = TRUE)


data.3 <- data.2 %>%
  droplevels()

saveRDS(data.3, 'output/study4/data_management/wBL_fulldata.rds')
saveRDS(table_names_all, 'output/study4/data_management/table_names_wBL.rds')


################################################################################
#                    descriptives split by DV
################################################################################

# Split by ng4_class
demos_by_class <- data.2 %>%
  filter (!is.na(ng4_class)) %>%
  select(-participant_id) %>%
  group_split(ng4_class)

# Create summary stats per class
obs_num_by_class <- map_dfr(
  demos_by_class,
  ~ get_summary_stats(.x, type = "common") %>%
    select(variable, n, mean, sd, min, max) %>%
    left_join(table_names_numeric, by = "variable") %>%
    relocate(domain_name, table_name, variable) %>%
    mutate(class = unique(.x$ng4_class)),
  .id = "class_id"
)

# Spread wide so each class’s summary stats appear as columns
obs_num_by_class_wide <- obs_num_by_class %>%
  select(-class_id) %>%
  pivot_wider(
    names_from = class,
    values_from = c(n, mean, sd, min, max),
    names_glue = "{.value}_{class}"
  )

# One-way ANOVA (numeric vars only)
anova_results <- data.2 %>%
  filter(!is.na(ng4_class)) %>%
  select(-participant_id) %>%
  select(where(is.numeric)) %>%
  map_df(~ broom::tidy(aov(. ~ data.2$ng4_class)), .id = "variable") %>%
  select(variable, statistic, p.value) %>%
  rename(F.statistic = statistic) %>%
  mutate(across(c(F.statistic, p.value), \(x) round(x, 3))) %>%   # <- updated form
  left_join(table_names_numeric, by = "variable") %>%
  relocate(domain_name, table_name, variable, F.statistic, p.value) %>%
  mutate(p.value = ifelse(p.value == 0.000, "<0.001", p.value)) %>%
  filter(!is.na(F.statistic) & !is.na(p.value))


# Merge summaries + ANOVA results
obs_num_ng4 <- obs_num_by_class_wide %>%
  full_join(anova_results, by = c("domain_name", "table_name", "variable"))

# add headers
temp_su <- obs_num_ng4 %>% 
  subset(domain_name == 'Self and Peer Involvement with Substance Use') %>% 
  add_row(table_name = 'Self and Peer Involvement with Substance Use', 
          .before = 1) 

temp_pb <- obs_num_ng4 %>% 
  subset(domain_name == 'Parenting Behaviors') %>% 
  add_row(table_name = 'Parenting Behaviors', .before = 1) 

temp_mh <- obs_num_ng4 %>% 
  subset(domain_name == 'Mental Health') %>% 
  add_row(table_name = 'Mental Health', .before = 1) 

temp_ph <- obs_num_ng4 %>% 
  subset(domain_name == 'Physical Health') %>% 
  add_row(table_name = 'Physical Health', .before = 1) 

temp_ce <- obs_num_ng4 %>% 
  subset(domain_name == 'Culture & Environment') %>% 
  add_row(table_name = 'Culture & Environment', .before = 1) 

temp_hor <- obs_num_ng4 %>% 
  subset(domain_name == 'Hormones') %>% 
  add_row(table_name = 'Hormones', .before = 1) 

temp_nc <- obs_num_ng4 %>% 
  subset(domain_name == 'Neurocognitive Factors') %>% 
  add_row(table_name = 'Neurocognitive Factors', .before = 1) 

temp_mri <- obs_num_ng4 %>% 
  subset(domain_name == 'Imaging') %>% 
  add_row(table_name = 'Imaging', .before = 1) 

obs_num_DV <- rbind(
  temp_su, temp_pb, temp_mh, temp_ph, temp_ce, temp_hor, temp_nc, temp_mri) %>% 
  select(-domain_name)

rm(
  temp_su, temp_pb, temp_mh, temp_ph, temp_ce, temp_hor, temp_nc, temp_mri)

# export
write.csv(obs_num_DV, 'output/study4/data_management/obs_num_DV_wBL.csv') 
rm(obs_num_DV, obs_num_ng4, obs_num_by_class, obs_num_by_class_wide, anova_results)

################################################################################
#                    descriptives for whole sample
################################################################################

obs_num <- get_summary_stats(data.2, type = "common") %>% 
  select(variable, n, mean, sd, min, max) %>% 
  rename(
    n_obs = n,
    mean_obs = mean,
    sd_obs = sd,
    min_obs = min,
    max_obs = max) %>% 
  slice(-1) %>% 
  left_join(table_names_numeric, by = 'variable') %>% 
  relocate(domain_name, table_name, variable)

# add headers
temp_demo <- obs_num %>% 
  subset(domain_name == 'Demographics') %>% 
  add_row(table_name = 'Demographics', .before = 1) 

temp_su <- obs_num %>% 
  subset(domain_name == 'Self and Peer Involvement with Substance Use') %>% 
  add_row(table_name = 'Self and Peer Involvement with Substance Use', 
          .before = 1) 

temp_pb <- obs_num %>% 
  subset(domain_name == 'Parenting Behaviors') %>% 
  add_row(table_name = 'Parenting Behaviors', .before = 1) 

temp_mh <- obs_num %>% 
  subset(domain_name == 'Mental Health') %>% 
  add_row(table_name = 'Mental Health', .before = 1) 

temp_ph <- obs_num %>% 
  subset(domain_name == 'Physical Health') %>% 
  add_row(table_name = 'Physical Health', .before = 1) 

temp_ce <- obs_num %>% 
  subset(domain_name == 'Culture & Environment') %>% 
  add_row(table_name = 'Culture & Environment', .before = 1) 

temp_hor <- obs_num %>% 
  subset(domain_name == 'Hormones') %>% 
  add_row(table_name = 'Hormones', .before = 1) 

temp_nc <- obs_num %>% 
  subset(domain_name == 'Neurocognitive Factors') %>% 
  add_row(table_name = 'Neurocognitive Factors', .before = 1) 

temp_mri <- obs_num %>% 
  subset(domain_name == 'Imaging') %>% 
  add_row(table_name = 'Imaging', .before = 1) 

obs_num <- rbind(
  temp_demo, temp_su, temp_pb, temp_mh, temp_ph, temp_ce, temp_hor, temp_nc, temp_mri) %>% 
  select(-domain_name)

rm(
  temp_demo, temp_su, temp_pb, temp_mh, temp_ph, temp_ce, temp_hor, temp_nc, temp_mri)

write.csv(obs_num, 'output/study4/data_management/obs_num_wBL.csv') 

# Sample Characteristics: Categorical Variables --------------------------------
# Function to summarize categorical variable frequencies and proportions
factor_sum <- function(x) {
  n <- table(x, useNA = "no")
  proportion <- prop.table(n)
  percentage <- proportion * 100
  as.data.frame(cbind(n, proportion = round(proportion, 4), percentage = round(percentage, 2)))
}

# Split dataset by ng4_class
demos_by_class <- data.2 %>%
  filter(!is.na(ng4_class)) %>%
  select(all_of(factors)) %>%
  group_split(ng4_class)

# Compute summary tables for each class and each categorical variable
obs_factor_ng4 <- map_dfr(
  demos_by_class,
  function(df_class) {
    class_id <- unique(df_class$ng4_class)
    map_dfr(
      factors[factors != "ng4_class"],
      function(varname) {
        tab <- factor_sum(df_class[[varname]]) %>%
          rownames_to_column("category") %>%
          mutate(variable = varname, class = class_id) %>%
          relocate(class, variable, category)
      }
    )
  }
)

# Pivot wider so each class has separate n/proportion/percentage columns
obs_factor_ng4_wide <- obs_factor_ng4 %>%
  pivot_wider(
    names_from = class,
    values_from = c(n, proportion, percentage),
    names_glue = "{.value}_{class}"
  )

# Missing Data in Observed Sample ----------------------------------------------

# names of variables w/missing data
NA_obs_numeric <- data.2 %>% 
  select_if(~ any(is.na(.))) %>% 
  select_if(~is.numeric(.))

NA_obs_numeric <- NA_obs_numeric %>%
  Desc(., plotit = FALSE) 
capture.output(
  NA_obs_numeric, 
  file = 'output/study4/data_management/supplement_wBL/NA_obs_numeric.txt')

NA_obs_factor <- data.2 %>% 
  select_if(~ any(is.na(.))) %>% 
  select_if(~is.factor(.))

NA_obs_factor <- NA_obs_factor %>% 
  Desc(., plotit = FALSE) 
capture.output(
  NA_obs_factor, 
  file = 'output/study4/data_management/supplement_wBL/NA_obs_factor.txt')

rm(NA_obs_numeric, NA_obs_factor)

# Fully automated imputation and elastic net pipeline --------------------------
system("powercfg -change standby-timeout-ac 0")  # Disables sleep


tryCatch({
  # ==============================================================================
  # MICE PREPROCESSING FOR 4-LEVEL OUTCOME (PAIRWISE COMPARISONS)
  # ==============================================================================
  # This script automates MICE imputation for all pairwise comparisons of a
  # 4-level outcome variable (6 total comparisons: 4 choose 2)
  # ==============================================================================
  
  # Load Required Packages -------------------------------------------------------
  library(dplyr)
  library(tidyr)
  library(mice)
  library(rsample)
  library(recipes)
  library(doParallel)
  library(doRNG)
  library(rstatix)
  library(fastDummies)
  library(blastula)
  library(keyring)
  
  
  
  # Configure Parallelization ----------------------------------------------------
  nCores <- min(parallel::detectCores()) - 1
  options(mc.cores = nCores, cores = nCores)
  doParallel::registerDoParallel(cores = nCores)
  cat("### Using", foreach::getDoParWorkers(), "cores\n")
  
  # Helper Functions -------------------------------------------------------------
  factor_sum <- function(x) {
    n <- table(x)
    proportion <- round(prop.table(n), 4)
    percentage <- (round(prop.table(n), 4)) * 100
    OUT <- cbind(n, proportion, percentage)
  }
  
  # Generate All Pairwise Comparisons --------------------------------------------
  ng4_levels <- 1:4
  pairs <- combn(ng4_levels, 2, simplify = FALSE)
  
  cat("\n### Pairwise Comparisons to Process:\n")
  for (i in seq_along(pairs)) {
    cat(sprintf("Comparison %d: Level %d vs Level %d\n", 
                i, pairs[[i]][1], pairs[[i]][2]))
  }
  
  # Master Function: Process One Split for One Comparison -----------------------
  process_split_comparison <- function(data, split_seed, split_num, 
                                       comparison_num, level1, level2,
                                       min_class_size = 30) {
    
    cat(sprintf("\n========================================\n"))
    cat(sprintf("Processing Split %d, Comparison %d (%d vs %d)\n", 
                split_num, comparison_num, level1, level2))
    cat(sprintf("========================================\n"))
    
    # Subset data for this comparison
    data_subset <- data %>%
      filter(ng4_class %in% c(level1, level2)) %>%
      mutate(
        # Create binary outcome: level1 = 1, level2 = 0
        # Keep as numeric throughout to avoid factor conversion issues
        ng4_class_binary = ifelse(ng4_class == level2, 1, 0),
        ng4_class = ng4_class_binary
      ) %>%
      select(-ng4_class_binary)
    
    # Check total sample size
    n_level2 <- sum(data_subset$ng4_class == 1)
    n_level1 <- sum(data_subset$ng4_class == 0)
    
    cat(sprintf("Total sample sizes - Level %d: %d, Level %d: %d\n",
                level1, n_level1, level2, n_level2))
    cat(sprintf("Binary coding: Level %d = 0, Level %d = 1\n", level1, level2))
    
    # Check if either class is too small
    if (n_level1 < min_class_size || n_level2 < min_class_size) {
      warning(sprintf("⚠ Insufficient sample size for comparison %d vs %d (n1=%d, n2=%d). Skipping.\n",
                      level1, level2, n_level1, n_level2))
      return(NULL)
    }
    
    # Calculate minimum expected training size
    min_train_expected <- min(n_level1, n_level2) * 0.75
    
    if (min_train_expected < 20) {
      warning(sprintf("⚠ Expected training class size too small (%d). Adjusting split proportion.\n",
                      floor(min_train_expected)))
      prop_train <- max(0.85, 25 / min(n_level1, n_level2))
      prop_train <- min(prop_train, 0.95)
    } else {
      prop_train <- 0.75
    }
    
    cat(sprintf("Using train proportion: %.2f\n", prop_train))
    
    # Create train/test split with adjusted proportion
    # Create a temporary factor for stratification only
    data_subset_temp <- data_subset %>%
      mutate(ng4_class_factor = factor(ng4_class))
    
    set.seed(split_seed)
    split_obj <- initial_split(data_subset_temp, prop = prop_train, strata = ng4_class_factor)
    
    trn <- rsample::training(split_obj) %>% select(-ng4_class_factor)
    tst <- rsample::testing(split_obj) %>% select(-ng4_class_factor)
    
    # Verify minimum class sizes in training set
    trn_level1 <- sum(trn$ng4_class == 0)
    trn_level2 <- sum(trn$ng4_class == 1)
    
    if (trn_level1 < 10 || trn_level2 < 10) {
      warning(sprintf("⚠ Training set has insufficient observations (n1=%d, n2=%d). Skipping.\n",
                      trn_level1, trn_level2))
      return(NULL)
    }
    
    cat(sprintf("Training sizes - Level %d (coded as 1): %d, Level %d (coded as 0): %d\n",
                level1, trn_level1, level2, trn_level2))
    
    # Validation checks
    validate_split(trn, tst, level1, level2)
    
    # Run MICE on training and test
    trn_imputed <- run_mice_pipeline(trn, split_num, comparison_num, "trn")
    tst_imputed <- run_mice_pipeline(tst, split_num, comparison_num, "tst")
    
    # Return results
    list(
      trn_obs = trn,
      tst_obs = tst,
      trn_imp = trn_imputed$imp_df,
      tst_imp = tst_imputed$imp_df,
      trn_mice = trn_imputed$mice_obj,
      tst_mice = tst_imputed$mice_obj,
      sample_sizes = list(
        total = c(n_level1, n_level2),
        train = c(trn_level1, trn_level2),
        test = c(sum(tst$ng4_class == 0), sum(tst$ng4_class == 1))
      )
    )
  }
  
  # Validation Function ----------------------------------------------------------
  validate_split <- function(trn, tst, level1, level2) {
    
    # ID check
    id_overlap <- any(trn$participant_id %in% tst$participant_id)
    if (!id_overlap) {
      cat("✓ No ID overlap between train and test\n")
    } else {
      stop("✗ ID overlap detected between train and test")
    }
    
    # Outcome distribution check
    trn_prop <- mean(trn$ng4_class == level1)
    tst_prop <- mean(tst$ng4_class == level1)
    
    if (trn_prop >= 0.4 && trn_prop <= 0.6 && tst_prop >= 0.4 && tst_prop <= 0.6) {
      cat(sprintf("✓ Outcome proportions reasonable - Train: %.2f, Test: %.2f\n",
                  trn_prop, tst_prop))
    } else {
      warning(sprintf("⚠ Outcome proportions imbalanced - Train: %.2f, Test: %.2f\n",
                      trn_prop, tst_prop))
    }
    
    # Factor variable check
    factor_vars <- c('race_4l', 'eth_hisp', 'income', 'religion', 'p_edu', 
                     'sex_2l', 'rec_bin', 'exp_sub', 'tbi_injury', 
                     'mh_p_kbi__school_006', 'det_susp', 'se_services', 'cct')
    
    all_represented <- TRUE
    for (var in factor_vars) {
      if (var %in% names(trn)) {
        trn_levels <- length(unique(trn[[var]]))
        tst_levels <- length(unique(tst[[var]]))
        if (trn_levels != tst_levels) {
          warning(sprintf("⚠ %s has different levels: Train=%d, Test=%d\n",
                          var, trn_levels, tst_levels))
          all_represented <- FALSE
        }
      }
    }
    
    if (all_represented) {
      cat("✓ All factor levels represented in both sets\n")
    }
  }
  
  # MICE Pipeline Function -------------------------------------------------------
  run_mice_pipeline <- function(data, split_num, comparison_num, dataset_type) {
    
    cat(sprintf("\n--- MICE for %s dataset ---\n", toupper(dataset_type)))
    
    # ============================================================================
    # FIX: Drop unused factor levels to prevent dimension mismatches
    # ============================================================================
    factor_cols <- sapply(data, is.factor)
    if (any(factor_cols)) {
      data[factor_cols] <- lapply(data[factor_cols], droplevels)
      cat("✓ Dropped unused factor levels\n")
    }
    
    # ============================================================================
    # FIX: Check for and remove constant/near-constant variables
    # ============================================================================
    n_unique <- sapply(data, function(x) {
      if (is.factor(x)) {
        length(levels(x))
      } else {
        length(unique(na.omit(x)))
      }
    })
    
    constant_vars <- names(n_unique)[n_unique <= 1]
    
    if (length(constant_vars) > 0) {
      cat(sprintf("⚠ Removing %d constant/single-level variables: %s\n", 
                  length(constant_vars), 
                  paste(constant_vars, collapse = ", ")))
      data <- data[, !names(data) %in% constant_vars, drop = FALSE]
    }
    
    # Extract variable names
    all_vars <- names(data)
    miss_vars <- names(data)[colSums(is.na(data)) > 0]
    
    cat(sprintf("Variables with missingness: %d\n", length(miss_vars)))
    
    # Build predictor matrix
    predictor_matrix <- matrix(0, ncol = length(all_vars), nrow = length(all_vars))
    rownames(predictor_matrix) <- all_vars
    colnames(predictor_matrix) <- all_vars
    
    # Specify imputer variables
    imputer_vars <- c('age_baseline', 'sex_2l', 'race_4l', 'eth_hisp', 
                      'income', 'religion', 'p_edu')
    imputer_vars <- intersect(unique(imputer_vars), all_vars)
    
    # Specify imputed variables
    imputed_only_vars <- miss_vars
    imputed_vars <- intersect(unique(c(imputed_only_vars, imputer_vars)), miss_vars)
    
    # ============================================================================
    # FIX: Check for collinearity and reduce predictor matrix if needed
    # ============================================================================
    # Compute correlation matrix for numeric variables
    numeric_vars <- names(data)[sapply(data, is.numeric)]
    numeric_vars <- setdiff(numeric_vars, "ng4_class")  # Exclude outcome
    
    if (length(numeric_vars) > 1) {
      cor_matrix <- cor(data[, numeric_vars], use = "pairwise.complete.obs")
      
      # Find highly correlated pairs (|r| > 0.95)
      high_cor <- which(abs(cor_matrix) > 0.95 & abs(cor_matrix) < 1, arr.ind = TRUE)
      
      if (nrow(high_cor) > 0) {
        high_cor_pairs <- unique(t(apply(high_cor, 1, sort)))
        high_cor_pairs <- high_cor_pairs[!duplicated(high_cor_pairs), , drop = FALSE]
        
        cat(sprintf("⚠ Found %d highly correlated variable pairs (|r| > 0.95)\n", 
                    nrow(high_cor_pairs)))
        
        # Remove one variable from each highly correlated pair
        vars_to_remove <- c()
        for (i in 1:nrow(high_cor_pairs)) {
          var1 <- rownames(cor_matrix)[high_cor_pairs[i, 1]]
          var2 <- rownames(cor_matrix)[high_cor_pairs[i, 2]]
          
          # Keep the variable with less missingness
          miss1 <- sum(is.na(data[[var1]]))
          miss2 <- sum(is.na(data[[var2]]))
          
          if (miss1 > miss2) {
            vars_to_remove <- c(vars_to_remove, var1)
          } else {
            vars_to_remove <- c(vars_to_remove, var2)
          }
        }
        
        vars_to_remove <- unique(vars_to_remove)
        
        if (length(vars_to_remove) > 0) {
          cat(sprintf("  Removing %d collinear variables: %s\n", 
                      length(vars_to_remove),
                      paste(vars_to_remove, collapse = ", ")))
          data <- data[, !names(data) %in% vars_to_remove, drop = FALSE]
          
          # Update variable lists
          all_vars <- names(data)
          miss_vars <- names(data)[colSums(is.na(data)) > 0]
          imputer_vars <- intersect(imputer_vars, all_vars)
          imputed_vars <- intersect(imputed_vars, all_vars)
          
          # Rebuild predictor matrix with new dimensions
          predictor_matrix <- matrix(0, ncol = length(all_vars), nrow = length(all_vars))
          rownames(predictor_matrix) <- all_vars
          colnames(predictor_matrix) <- all_vars
        }
      }
    }
    
    # Construct predictor matrix
    imputer_matrix <- predictor_matrix
    imputer_matrix[, imputer_vars] <- 1
    
    imputed_matrix <- predictor_matrix
    imputed_matrix[imputed_vars, ] <- 1
    
    predictor_matrix <- imputer_matrix * imputed_matrix
    diag(predictor_matrix) <- 0
    
    # ============================================================================
    # FIX: Limit number of predictors per variable to avoid over-specification
    # ============================================================================
    max_predictors <- floor(nrow(data) / 3)  # Rule of thumb: n/3 predictors max
    
    for (var in rownames(predictor_matrix)) {
      n_predictors <- sum(predictor_matrix[var, ])
      
      if (n_predictors > max_predictors) {
        cat(sprintf("⚠ Variable '%s' has %d predictors (max: %d), reducing...\n", 
                    var, n_predictors, max_predictors))
        
        # Keep only the most important predictors (imputer_vars + highest correlation)
        predictors_for_var <- names(which(predictor_matrix[var, ] == 1))
        
        # Always keep imputer_vars
        keep_vars <- intersect(predictors_for_var, imputer_vars)
        
        # Add other predictors based on correlation
        other_predictors <- setdiff(predictors_for_var, imputer_vars)
        
        if (length(other_predictors) > 0 && length(keep_vars) < max_predictors) {
          if (is.numeric(data[[var]])) {
            # Calculate correlations with numeric predictors
            numeric_predictors <- intersect(other_predictors, numeric_vars)
            
            if (length(numeric_predictors) > 0) {
              cors <- abs(cor(data[[var]], data[, numeric_predictors], 
                              use = "pairwise.complete.obs"))
              top_predictors <- names(sort(cors[1,], decreasing = TRUE)[1:(max_predictors - length(keep_vars))])
              keep_vars <- c(keep_vars, top_predictors)
            }
          }
          
          # Fill remaining slots with other predictors if needed
          if (length(keep_vars) < max_predictors) {
            remaining_slots <- max_predictors - length(keep_vars)
            keep_vars <- c(keep_vars, other_predictors[1:min(remaining_slots, length(other_predictors))])
          }
        }
        
        # Update predictor matrix
        predictor_matrix[var, ] <- 0
        predictor_matrix[var, keep_vars] <- 1
      }
    }
    
    # ============================================================================
    # FIX: Wrap dry run in tryCatch to catch dimension errors early
    # ============================================================================
    dry_mice <- tryCatch({
      mice(data = data, m = 1, predictorMatrix = predictor_matrix, maxit = 0)
    }, error = function(e) {
      cat("✗ Error in dry MICE run:\n")
      cat(sprintf("  %s\n", e$message))
      
      # Try with default predictor matrix as fallback
      cat("  Attempting with quickpred predictor matrix...\n")
      pred_quick <- quickpred(data, mincor = 0.1, minpuc = 0.1)
      mice(data = data, m = 1, predictorMatrix = pred_quick, maxit = 0)
    })
    
    predictor_matrix <- dry_mice$predictorMatrix
    dry_mice$method[setdiff(all_vars, imputed_vars)] <- ""
    
    # ============================================================================
    # FIX: Add verbose output for first imputation to catch issues
    # ============================================================================
    M <- 5
    set.seed(808915)
    
    # Run first imputation with verbose to diagnose issues
    cat("Running first imputation (verbose)...\n")
    mice_first <- tryCatch({
      mice::mice(data = data, m = 1, print = FALSE,
                 predictorMatrix = predictor_matrix,
                 method = dry_mice$method, maxit = 20)
    }, error = function(e) {
      cat(sprintf("✗ Error in first MICE imputation: %s\n", e$message))
      
      # Print diagnostic information
      cat("\nDiagnostic Information:\n")
      cat(sprintf("  Dataset dimensions: %d rows x %d cols\n", nrow(data), ncol(data)))
      cat("\n  Factor variables and their levels:\n")
      for (var in names(data)[sapply(data, is.factor)]) {
        cat(sprintf("    %s: %d levels - %s\n", 
                    var, 
                    nlevels(data[[var]]),
                    paste(levels(data[[var]]), collapse = ", ")))
      }
      stop(e)
    })
    
    # Continue with remaining imputations in parallel
    if (M > 1) {
      cat(sprintf("Running remaining %d imputations in parallel...\n", M - 1))
      mice_rest <- foreach(i = seq_len(M - 1), .combine = ibind) %dorng% {
        mice::mice(data = data, m = 1, print = FALSE,
                   predictorMatrix = predictor_matrix,
                   method = dry_mice$method, maxit = 20)
      }
      
      mice_out <- ibind(mice_first, mice_rest)
    } else {
      mice_out <- mice_first
    }
    
    # Check for logged events
    if (is.null(mice_out$loggedEvents)) {
      cat("✓ No logged events during imputation\n")
    } else {
      warning("⚠ Logged events detected:\n")
      print(mice_out$loggedEvents)
    }
    
    # Create long format dataframe
    imp_df <- mice::complete(mice_out, "long", include = TRUE)
    
    cat(sprintf("✓ MICE completed: %d imputations\n", M))
    
    list(
      mice_obj = mice_out,
      imp_df = imp_df
    )
  }
  # Center and Scale Function ----------------------------------------------------
  center_scale_data <- function(imp_df) {
    
    cat("\n--- Center and Scale ---\n")
    
    # IMPORTANT: ng4_class should already be numeric (0/1) at this point
    # Verify this
    cat(sprintf("ng4_class unique values before C&S: %s\n", 
                paste(unique(imp_df$ng4_class[imp_df$.imp != 0]), collapse=", ")))
    
    # Store ng4_class and remove observed data
    imp_df_cs <- imp_df %>% subset(.imp != 0)
    
    # Store observed data separately
    imp_df_imp0 <- imp_df %>% subset(.imp == 0) %>% select(-participant_id)
    
    # Store .imp, .id, pds, and ng4_class (don't normalize these)
    imp_df_imp15 <- imp_df %>% 
      subset(.imp != 0) %>% 
      select(.imp, .id, pds, ng4_class)
    
    # Build recipe - exclude ng4_class from normalization
    recipe_prep <- imp_df_cs %>%
      recipe(ng4_class ~ .) %>%
      step_rm(ends_with('.id') | ends_with('.imp') | participant_id | pds) %>%
      update_role(ng4_class, new_role = "outcome") %>%
      step_normalize(all_numeric_predictors()) %>%
      prep()
    
    # Apply recipe
    imp_df_cs_2 <- recipe_prep %>% bake(new_data = imp_df_cs)
    
    # Remove ng4_class from centered/scaled data (we'll add it back from imp_df_imp15)
    imp_df_cs_2 <- imp_df_cs_2 %>% select(-ng4_class)
    
    # Combine: add back .imp, .id, pds, and ng4_class
    imp_df_cs_3 <- cbind(imp_df_imp15, imp_df_cs_2)
    
    # Add back observed data
    imp_df_cs_final <- rbind(imp_df_imp0, imp_df_cs_3)
    
    # Validate ng4_class is still 0/1
    unique_outcome <- unique(imp_df_cs_final$ng4_class[imp_df_cs_final$.imp != 0])
    cat(sprintf("ng4_class unique values after C&S: %s\n", 
                paste(sort(unique_outcome), collapse=", ")))
    
    if (!all(unique_outcome %in% c(0, 1))) {
      stop("ERROR: ng4_class is not binary (0/1) after center and scale!")
    }
    
    # Validate centering and scaling of predictors
    cs_check <- imp_df_cs_final %>%
      subset(.imp != 0) %>%
      select_if(is.numeric) %>%
      select(-.imp, -.id, -ng4_class) %>%
      get_summary_stats(type = "common") %>%
      select(variable, mean, sd)
    
    mean_ok <- all(cs_check$mean >= -0.10 & cs_check$mean <= 0.10)
    sd_ok <- all(cs_check$sd >= 0.95 & cs_check$sd <= 1.05)
    
    if (mean_ok && sd_ok) {
      cat("✓ Center and scale validated\n")
    } else {
      warning("⚠ Center and scale validation failed\n")
    }
    
    imp_df_cs_final
  }
  
  # Create Dummy Variables -------------------------------------------------------
  create_dummy_variables <- function(imp_df_cs_final) {
    
    cat("\n--- Creating Dummy Variables ---\n")
    
    # Extract factor variables
    factor_vars <- imp_df_cs_final %>% select_if(is.factor)
    
    # Create dummy codes for ALL categories
    factor_dummy <- fastDummies::dummy_cols(factor_vars, remove_first_dummy = FALSE)
    
    # Define expected dummy variable names (excluding reference groups)
    expected_dummies <- c(
      #Baseline Use
      'alc0_Yes', 'mj0_Yes', 'nic0_Yes', 'oth0_Yes',
      # Demographics
      "race_4l_Asian", "race_4l_Black", "race_4l_Other_MultiRacial",
      "eth_hisp_Hispanic",
      "income_inc_1", "income_inc_2", "income_inc_3",
      "religion_rp_1", "religion_rp_2", "religion_rp_3", "religion_rp_4",
      "religion_rp_5", "religion_rp_6", "religion_rp_7", "religion_rp_8",
      "religion_rp_9", "religion_rp_10", "religion_rp_11", "religion_rp_12",
      "religion_rp_13", "religion_rp_14", "religion_rp_15", "religion_rp_16",
      "p_edu_Less_than_HS_Degree_GED_Equivalent",
      "p_edu_HS_Graduate_GED_Equivalent",
      "p_edu_Some_College_or_Associates_Degree",
      "p_edu_Masters_Degree",
      "p_edu_Professional_School_or_Doctoral_Degree",
      "sex_2l_Female",
      # Physical Health
      "rec_bin_Yes", "exp_sub_Yes", "tbi_injury_Yes",
      # Culture & Environment
      "mh_p_kbi__school_006_Grade_B", "mh_p_kbi__school_006_Grade_C",
      "mh_p_kbi__school_006_Grade_Fail",
      "det_susp_Yes",
      "se_services_Emotion_or_Learning_Support", "se_services_Gifted",
      "se_services_Other", "se_services_Combined_Services",
      # Neurocognitive
      "cct_Immediate"
    )
    
    # Select only the expected dummies that exist, create zeros for missing ones
    dummy_selected <- data.frame(matrix(0, nrow = nrow(factor_dummy), ncol = length(expected_dummies)))
    colnames(dummy_selected) <- expected_dummies
    
    # Fill in values for dummies that exist
    for (col in expected_dummies) {
      if (col %in% colnames(factor_dummy)) {
        dummy_selected[[col]] <- factor_dummy[[col]]
      } else {
        cat(sprintf("  ⚠ Creating zero column for missing dummy: %s\n", col))
      }
    }
    
    # Remove original factor variables
    imp_df_final <- imp_df_cs_final %>%
      select(-race_4l, -eth_hisp, -income, -religion, -p_edu, -sex_2l,
             -rec_bin, -exp_sub, -tbi_injury,
             -mh_p_kbi__school_006, -det_susp, -se_services, -cct)
    
    # Add dummy codes
    imp_df_final <- cbind(imp_df_final, dummy_selected)
    
    cat(sprintf("✓ Total variables after dummy coding: %d\n", ncol(imp_df_final)))
    cat(sprintf("✓ Total dummy variables: %d\n", length(expected_dummies)))
    
    imp_df_final
  }
  
  
  # Create MICE Object -----------------------------------------------------------
  create_mice_object <- function(imp_df_final) {
    
    obs_df <- imp_df_final %>% subset(.imp == 0) %>% select(-.imp, -.id)
    
    where_matrix <- matrix(TRUE, nrow = nrow(obs_df), ncol = ncol(obs_df))
    colnames(where_matrix) <- colnames(obs_df)
    
    mice_obj <- as.mids(imp_df_final, where = where_matrix, .imp = ".imp", .id = ".id")
    
    mice_obj
  }
  
  # Master Processing Loop -------------------------------------------------------
  # This processes all splits for all comparisons
  
  # Define split seeds
  split_seeds <- c(11897, 95039, 48378, 21773, 48339, 
                   88593, 87388, 12993, 52701, 34885)
  
  # Initialize storage
  results <- list()
  data.3 <- readRDS('output/study4/data_management/wBL_fulldata.rds')
  table_names_all <- readRDS('output/study4/data_management/table_names_wBL.rds')
  
  for (i in seq_along(pairs)) {
    
  }
  
  # Loop through each comparison
  for (comp_idx in seq_along(pairs)) {
    job::job({
      p <- pairs[[comp_idx]]
      level_order <- p[order(-table(data.3[data.3$ng4_class %in% p,]$ng4_class)[p])]
      level1 <- level_order[1]
      level2 <- level_order[2]
      
      cat(sprintf("\n\n"))
      cat(sprintf("================================================================================\n"))
      cat(sprintf("COMPARISON %d: LEVEL %d vs LEVEL %d\n", comp_idx, level1, level2))
      cat(sprintf("================================================================================\n"))
      
      # Loop through each split
      for (split_idx in 1:10) {
        
        tryCatch({
          
          # Process split
          split_results <- process_split_comparison(
            data = data.3 ,  # original dataset
            split_seed = split_seeds[split_idx],
            split_num = split_idx,
            comparison_num = comp_idx,
            level1 = level1,
            level2 = level2
          )
          
          # Center and scale
          split_results$trn_imp_cs <- center_scale_data(split_results$trn_imp)
          split_results$tst_imp_cs <- center_scale_data(split_results$tst_imp)
          
          # Create dummy variables
          split_results$trn_final <- create_dummy_variables(split_results$trn_imp_cs)
          split_results$tst_final <- create_dummy_variables(split_results$tst_imp_cs)
          
          # Create final MICE objects
          split_results$trn_mice_final <- create_mice_object(split_results$trn_final)
          split_results$tst_mice_final <- create_mice_object(split_results$tst_final)
          
          # Store results with naming convention
          comp_name <- sprintf("comp%d_split%d", comp_idx, split_idx)
          results[[comp_name]] <- split_results
          
          # Assign to global environment with meaningful names
          assign(sprintf("trn_%d_%dvs%d_obs", split_idx, level1, level2),
                 split_results$trn_obs, envir = .GlobalEnv)
          assign(sprintf("tst_%d_%dvs%d_obs", split_idx, level1, level2),
                 split_results$tst_obs, envir = .GlobalEnv)
          assign(sprintf("trn_%d_%dvs%d_final", split_idx, level1, level2),
                 split_results$trn_final, envir = .GlobalEnv)
          assign(sprintf("tst_%d_%dvs%d_final", split_idx, level1, level2),
                 split_results$tst_final, envir = .GlobalEnv)
          assign(sprintf("trn_%d_%dvs%d_mice", split_idx, level1, level2),
                 split_results$trn_mice_final, envir = .GlobalEnv)
          assign(sprintf("tst_%d_%dvs%d_mice", split_idx, level1, level2),
                 split_results$tst_mice_final, envir = .GlobalEnv)
          
          cat(sprintf("\n✓✓✓ Successfully completed Split %d, Comparison %d ✓✓✓\n",
                      split_idx, comp_idx))
          
        }, error = function(e) {
          cat(sprintf("\n✗✗✗ ERROR in Split %d, Comparison %d: %s ✗✗✗\n",
                      split_idx, comp_idx, e$message))
        })}
    })
  }
  
  # Save Results -----------------------------------------------------------------
  result_list <- c(
    paste0("trn_", 1:10, "_1vs2_mice"),
    paste0("trn_", 1:10, "_1vs3_mice"),
    paste0("trn_", 1:10, "_2vs3_mice"),
    paste0("trn_", 1:10, "_1vs4_mice"),
    paste0("trn_", 1:10, "_2vs4_mice"),
    paste0("trn_", 1:10, "_3vs4_mice")
  )
  while (any(!sapply(result_list, exists)) ||
         any(sapply(result_list, function(x) exists(x) && length(get(x)) == 0))) {
    message("Imputation incomplete. Waiting 30 minutes...")
    Sys.sleep(30 * 60)  # sleep for 30 minutes
    # refresh list here if needed, e.g.
    # my_list <- get_new_data()
  }
  
  cat("\n\n=== Saving Results ===\n")
  save(list = ls(pattern = "^(trn_|tst_)"), file = 'output/study4/MICE_prep/imputed_data_wBL.RData')
  cat("✓ Saved to mice_pairs.RData\n")
  
  
  
  
  email <- compose_email(
    body = "Imputation complete."
  )
  
  smtp_send(
    email,
    from = " email here",
    to   = "destination email here",
    subject = "Imputation step done",
    credentials = creds_key("gmail_smtp")
  )
  
  # Summary ----------------------------------------------------------------------
  cat("\n\n=== PROCESSING COMPLETE ===\n")
  cat(sprintf("Total comparisons processed: %d\n", length(pairs)))
  cat(sprintf("Total splits per comparison: 10\n"))
  cat(sprintf("Total datasets created: %d\n", length(results)))
  cat("\nObject naming convention:\n")
  cat("  - trn_<split>_<level1>vs<level2>_obs    : Training observed data\n")
  cat("  - tst_<split>_<level1>vs<level2>_obs    : Test observed data\n")
  cat("  - trn_<split>_<level1>vs<level2>_final  : Training imputed/processed\n")
  cat("  - tst_<split>_<level1>vs<level2>_final  : Test imputed/processed\n")
  cat("  - trn_<split>_<level1>vs<level2>_mice   : Training MICE object\n")
  cat("  - tst_<split>_<level1>vs<level2>_mice   : Test MICE object\n")
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(glmnet)
  library(rsample)
  library(purrr)
  library(pROC)
  library(tidyverse)
  library(job)
  library(ggsci)
  
  # #if necessary:
  #load("~/output/study4/MICE_prep/imputed_data_wMRI.RData")
  
  # Function: Run Elastic Net on MICE Data ---------------------------------------
  run_elastic_net_mice <- function(mice_obj, obs_df, alpha = 0.5, nfolds = 10) {
    
    cat("\n--- Running Elastic Net with MICE ---\n")
    
    # Extract number of imputations
    m <- mice_obj$m
    
    # Check sample size for CV folds
    n <- nrow(obs_df)
    min_class_size <- min(table(obs_df$ng4_class))
    
    # Adjust nfolds if necessary
    if (min_class_size < nfolds * 2) {
      nfolds_adjusted <- max(3, floor(min_class_size / 2))
      cat(sprintf("⚠ Adjusting nfolds from %d to %d due to small class size (%d)\n",
                  nfolds, nfolds_adjusted, min_class_size))
      nfolds <- nfolds_adjusted
    }
    
    # Further check: minimum 2 observations per fold per class
    if (floor(min_class_size / nfolds) < 2) {
      nfolds <- max(2, floor(min_class_size / 2))
      cat(sprintf("⚠ Further adjusting nfolds to %d to ensure min 2 obs/fold\n", nfolds))
    }
    
    # Storage for models and predictions
    models <- vector("list", m)
    lambdas <- numeric(m)
    max_aucs <- numeric(m)
    kept_columns <- vector("list", m)  # NEW: Track which columns are kept
    
    # Run elastic net on each imputation
    for (i in 1:m) {
      
      tryCatch({
        # Complete data for this imputation
        complete_data <- mice::complete(mice_obj, action = i)
        
        # Prepare X and Y
        y <- complete_data$ng4_class
        X <- complete_data %>% 
          select(-ng4_class, -age_baseline) %>%
          as.matrix()
        
        #Check for zero variance columns
        col_vars <- apply(X, 2, var)
        zero_var_cols <- which(col_vars == 0 | is.na(col_vars))
        
        if (length(zero_var_cols) > 0) {
          cat(sprintf("  Removing %d zero-variance columns\n", length(zero_var_cols)))
          X <- X[, -zero_var_cols, drop = FALSE]
        }
        
        # Store which columns were kept
        kept_columns[[i]] <- colnames(X)
        
        # Run cross-validated elastic net
        set.seed(123 + i)
        cv_model <- cv.glmnet(
          x = X,
          y = y,
          family = "binomial",
          alpha = alpha,
          nfolds = nfolds,
          type.measure = "auc",
          parallel = FALSE,  # Disable parallel to avoid nested parallelization issues
          standardize = FALSE  # Already standardized
        )
        
        models[[i]] <- cv_model
        lambdas[i] <- cv_model$lambda.min
        max_aucs[i] <- max(cv_model$cvm)
        
        cat(sprintf("  Imputation %d: Lambda = %.6f, Max AUC = %.4f\n",
                    i, cv_model$lambda.min, max(cv_model$cvm)))
        
      }, error = function(e) {
        cat(sprintf("  ✗ Error in imputation %d: %s\n", i, e$message))
        models[[i]] <- NULL
        lambdas[i] <- NA
        max_aucs[i] <- NA
      })
    }
    
    # Remove failed models
    valid_idx <- !sapply(models, is.null)
    
    if (sum(valid_idx) == 0) {
      stop("All imputations failed during elastic net fitting")
    }
    
    if (sum(valid_idx) < m) {
      cat(sprintf("⚠ Warning: %d out of %d imputations failed\n", m - sum(valid_idx), m))
      models <- models[valid_idx]
      lambdas <- lambdas[valid_idx]
      max_aucs <- max_aucs[valid_idx]
    }
    
    # Use median lambda across successful imputations
    lambda_final <- median(lambdas, na.rm = TRUE)
    cat(sprintf("\n✓ Final lambda (median): %.6f\n", lambda_final))
    cat(sprintf("✓ Mean training AUC: %.4f (SD: %.4f)\n", 
                mean(max_aucs, na.rm = TRUE), sd(max_aucs, na.rm = TRUE)))
    
    list(
      models = models,
      lambda_final = lambda_final,
      lambdas = lambdas,
      max_aucs = max_aucs,
      nfolds_used = nfolds,
      kept_columns = kept_columns  # NEW: Return column info
    )
  }
  
  # Function: Get Pooled Coefficients (Rubin's Rules) ---------------------------
  pool_coefficients <- function(enet_results, mice_obj, lambda = NULL) {
    
    if (is.null(lambda)) {
      lambda <- enet_results$lambda_final
    }
    
    m <- length(enet_results$models)
    
    # Extract coefficients from each imputation
    coef_list <- lapply(enet_results$models, function(model) {
      as.matrix(coef(model, s = lambda))
    })
    
    # Ensure all have same variables (elastic net may select different variables)
    all_vars <- unique(unlist(lapply(coef_list, rownames)))
    
    # Create matrix with all coefficients
    coef_matrix <- matrix(0, nrow = length(all_vars), ncol = m)
    rownames(coef_matrix) <- all_vars
    
    for (i in 1:m) {
      vars_i <- rownames(coef_list[[i]])
      coef_matrix[vars_i, i] <- coef_list[[i]][, 1]
    }
    
    # Pool using Rubin's rules (simple average for elastic net)
    pooled_coef <- rowMeans(coef_matrix)
    
    # Calculate stability (% of imputations where coefficient is non-zero)
    stability <- rowMeans(coef_matrix != 0) * 100
    
    # Create results dataframe
    results_df <- data.frame(
      variable = names(pooled_coef),
      coefficient = pooled_coef,
      stability_pct = stability,
      stringsAsFactors = FALSE
    ) %>%
      filter(coefficient != 0) %>%
      arrange(desc(abs(coefficient)))
    
    list(
      pooled_coef = pooled_coef,
      coef_matrix = coef_matrix,
      results_df = results_df
    )
  }
  
  # Function: Make Predictions on Test Data -------------------------------------
  predict_elastic_net_mice <- function(enet_results, pooled_coef, 
                                       test_mice_obj, test_obs_df) {
    
    m <- test_mice_obj$m
    lambda <- enet_results$lambda_final
    
    # Storage for predictions from each imputation
    pred_list <- vector("list", m)
    
    for (i in 1:m) {
      
      # Complete test data
      complete_test <- mice::complete(test_mice_obj, action = i)
      
      # Prepare X
      X_test <- complete_test %>%
        select(-ng4_class) %>%
        as.matrix()
      
      # Get columns that were kept during training
      kept_cols <- enet_results$kept_columns[[i]]
      
      # Find missing and extra columns
      missing_cols <- setdiff(kept_cols, colnames(X_test))
      extra_cols <- setdiff(colnames(X_test), kept_cols)
      
      # Handle missing columns
      if (length(missing_cols) > 0) {
        cat(sprintf("  ⚠ Test imputation %d missing %d columns, filling with zeros: %s\n",
                    i, length(missing_cols), 
                    paste(head(missing_cols, 3), collapse = ", ")))
        
        # Create matrix with missing columns (filled with 0)
        missing_matrix <- matrix(0, 
                                 nrow = nrow(X_test), 
                                 ncol = length(missing_cols))
        colnames(missing_matrix) <- missing_cols
        
        # Combine with existing test data
        X_test <- cbind(X_test, missing_matrix)
      }
      
      # Remove extra columns and reorder to match training
      X_test <- X_test[, kept_cols, drop = FALSE]
      
      # Verify dimensions match
      if (ncol(X_test) != length(kept_cols)) {
        stop(sprintf("Dimension mismatch after alignment: expected %d cols, got %d",
                     length(kept_cols), ncol(X_test)))
      }
      
      # Make predictions
      pred_list[[i]] <- predict(
        enet_results$models[[i]],
        newx = X_test,
        s = lambda,
        type = "response"
      )[, 1]
    }
    
    # Average predictions across imputations
    pred_matrix <- do.call(cbind, pred_list)
    pooled_predictions <- rowMeans(pred_matrix)
    
    # Get true outcomes from observed data
    y_true <- test_obs_df$ng4_class
    
    list(
      predictions = pooled_predictions,
      pred_matrix = pred_matrix,
      y_true = y_true
    )
  }
  
  # Function: Calculate Performance Metrics --------------------------------------
  calculate_metrics <- function(predictions, y_true) {
    
    # ROC curve and AUC
    roc_obj <- pROC::roc(y_true, predictions, quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    
    # Find optimal threshold using Youden's index
    coords_obj <- pROC::coords(roc_obj, "best", best.method = "youden")
    threshold <- coords_obj$threshold
    
    # Predictions at optimal threshold
    pred_class <- ifelse(predictions >= threshold, 1, 0)
    
    # Confusion matrix
    cm <- table(Predicted = pred_class, Actual = y_true)
    
    # Calculate metrics
    accuracy <- sum(diag(cm)) / sum(cm)
    sensitivity <- cm[2, 2] / sum(cm[, 2])
    specificity <- cm[1, 1] / sum(cm[, 1])
    ppv <- cm[2, 2] / sum(cm[2, ])
    npv <- cm[1, 1] / sum(cm[1, ])
    
    list(
      auc = auc_val,
      threshold = threshold,
      accuracy = accuracy,
      sensitivity = sensitivity,
      specificity = specificity,
      ppv = ppv,
      npv = npv,
      confusion_matrix = cm,
      roc_obj = roc_obj
    )
  }
  
  # Master Function: Complete Elastic Net Pipeline ------------------------------
  run_complete_pipeline <- function(split_num, comp_idx, level1, level2, 
                                    alpha = 0.5, nfolds = 10) {
    
    
    cat("\n")
    cat("========================================================================\n")
    cat(sprintf("ELASTIC NET: Split %d, Comparison %d (%d vs %d)\n", 
                split_num, comp_idx, level1, level2))
    cat("========================================================================\n")
    
    # Get object names
    trn_mice_name <- sprintf("trn_%d_%dvs%d_mice", split_num, level1, level2)
    tst_mice_name <- sprintf("tst_%d_%dvs%d_mice", split_num, level1, level2)
    trn_obs_name <- sprintf("trn_%d_%dvs%d_obs", split_num, level1, level2)
    tst_obs_name <- sprintf("tst_%d_%dvs%d_obs", split_num, level1, level2)
    
    # Get objects from global environment
    trn_mice <- get(trn_mice_name, envir = .GlobalEnv)
    tst_mice <- get(tst_mice_name, envir = .GlobalEnv)
    trn_obs <- get(trn_obs_name, envir = .GlobalEnv)
    tst_obs <- get(tst_obs_name, envir = .GlobalEnv)
    
    # 1. Train elastic net
    enet_results <- run_elastic_net_mice(trn_mice, trn_obs, alpha, nfolds)
    
    # 2. Pool coefficients
    pooled <- pool_coefficients(enet_results, trn_mice)
    
    cat("\n--- Selected Variables ---\n")
    print(head(pooled$results_df, 20))
    
    # 3. Predict on test set
    predictions <- predict_elastic_net_mice(enet_results, pooled$pooled_coef,
                                            tst_mice, tst_obs)
    
    # 4. Calculate metrics
    metrics <- calculate_metrics(predictions$predictions, predictions$y_true)
    
    cat("\n--- Test Set Performance ---\n")
    cat(sprintf("AUC: %.4f\n", metrics$auc))
    cat(sprintf("Accuracy: %.4f\n", metrics$accuracy))
    cat(sprintf("Sensitivity: %.4f\n", metrics$sensitivity))
    cat(sprintf("Specificity: %.4f\n", metrics$specificity))
    cat(sprintf("PPV: %.4f\n", metrics$ppv))
    cat(sprintf("NPV: %.4f\n", metrics$npv))
    
    # Return all results
    list(
      split = split_num,
      comparison = comp_idx,
      levels = c(level1, level2),
      enet_results = enet_results,
      pooled_coef = pooled,
      predictions = predictions,
      metrics = metrics
    )
  }
  
  # Run All Elastic Net Models ---------------------------------------------------
  cat("\n\n")
  cat("================================================================================\n")
  cat("RUNNING ELASTIC NET ON ALL COMPARISONS\n")
  cat("================================================================================\n")
  
  # Storage for all results
  enet_all_results <- list()
  
  # Loop through comparisons and splits
  for (comp_idx in seq_along(pairs)) {
    
    p <- pairs[[comp_idx]]
    level_order <- p[order(-table(data.3[data.3$ng4_class %in% p,]$ng4_class)[p])]
    level1 <- level_order[1]
    level2 <- level_order[2]
    
    job({
      for (split_idx in 1:10) {
        
        tryCatch({
          
          # Run complete pipeline
          result <- run_complete_pipeline(
            split_num = split_idx,
            comp_idx = comp_idx,
            level1 = level1,
            level2 = level2,
            alpha = 0.5,  # Elastic net (0 = ridge, 1 = lasso)
            nfolds = 10
          )
          
          # Store result
          result_name <- sprintf("comp%d_split%d", comp_idx, split_idx)
          enet_all_results[[result_name]] <- result
          
          # Assign to global environment
          assign(sprintf("enet_%d_%dvs%d", split_idx, level1, level2),
                 result, envir = .GlobalEnv)
          
        }, error = function(e) {
          cat(sprintf("\n✗ ERROR in Split %d, Comparison %d: %s\n",
                      split_idx, comp_idx, e$message))
        })
      }
      
      assign(sprintf('enet_results_comp%d', comp_idx), enet_all_results,
             envir = .GlobalEnv)}, title = sprintf("GLMNET_Comp%d", comp_idx),
      import = c(ls(pattern = paste0("^(trn_|tst_)\\d+_", level1, "vs", level2, "_(mice|obs)$")),
                 'comp_idx', 'level1', 'level2', 'calculate_metrics', 'pool_coefficients',
                 'predict_elastic_net_mice', 'run_complete_pipeline', 'run_elastic_net_mice',
                 'enet_all_results', 'pairs'))
  }
  
  # Save all elastic net results
  result_list <- c('enet_results_comp1', 'enet_results_comp2', 'enet_results_comp3',
                   'enet_results_comp4', 'enet_results_comp5', 'enet_results_comp6')
  while (any(!sapply(result_list, exists)) ||
         any(sapply(result_list, function(x) exists(x) && length(get(x)) == 0))) {
    message("Regressions incomplete. Waiting 30 minutes...")
    Sys.sleep(30 * 60)  # sleep for 30 minutes
    # refresh  list here if needed, e.g.
    # my_list <- get_new_data() 
  }
  enet_all_results = c(enet_results_comp1, enet_results_comp2, enet_results_comp3,
                       enet_results_comp4, enet_results_comp5, enet_results_comp6)
  save(list = c(ls(pattern = "^(enet_)"),'table_names_all', 'data.3'), file = "output/study4/en_results/claude_elastic_net_results_wBL.RData")
  cat("\n✓ All elastic net results saved to elastic_net_results.RData\n")
  
  email <- compose_email(
    body = "eNet complete."
  )
  
  smtp_send(
    email,
    from = " email here",
    to   = "destination email here",
    subject = "Elastic Net step 1 done",
    credentials = creds_key("gmail_smtp")
  )
  
  # Summary Statistics -----------------------------------------------------------
  cat("\n\n=== ELASTIC NET SUMMARY ===\n")
  
  # Extract AUCs for each comparison
  summary_df <- data.frame()
  
  for (comp_idx in seq_along(pairs)) {
    
    p <- pairs[[comp_idx]]
    level_order <- p[order(-table(data.3[data.3$ng4_class %in% p,]$ng4_class)[p])]
    level1 <- level_order[1]
    level2 <- level_order[2]
    
    aucs <- numeric(10)
    
    for (split_idx in 1:10) {
      result_name <- sprintf("comp%d_split%d", comp_idx, split_idx)
      if (result_name %in% names(enet_all_results)) {
        aucs[split_idx] <- enet_all_results[[result_name]]$metrics$auc
      }
    }
    
    summary_df <- rbind(summary_df, data.frame(
      comparison = sprintf("%d vs %d", level1, level2),
      mean_auc = mean(aucs),
      sd_auc = sd(aucs),
      min_auc = min(aucs),
      max_auc = max(aucs)
    ))
  }
  
  print(summary_df)
  
  # ==============================================================================
  # EXTRACT AND SUMMARIZE ELASTIC NET RESULTS FROM ENVIRONMENT
  # ==============================================================================
  # Extracts results from enet_X_XvsX objects in the global environment
  # and creates comprehensive summaries and visualizations
  # ==============================================================================
  
  # Define pairwise comparisons
  pairs <- combn(1:4, 2, simplify = FALSE)
  
  # Initialize storage
  all_results <- list()
  
  cat("=== Extracting Elastic Net Results from Environment ===\n\n")
  
  # Loop through all comparisons and splits
  for (comp_idx in seq_along(pairs)) {
    
    p <- pairs[[comp_idx]]
    level_order <- p[order(-table(data.3[data.3$ng4_class %in% p,]$ng4_class)[p])]
    level1 <- level_order[1]
    level2 <- level_order[2]
    
    cat(sprintf("\nComparison %d: %d vs %d\n", comp_idx, level1, level2))
    
    # Storage for this comparison's splits
    comp_results <- list()
    
    for (split_idx in 1:10) {
      
      # Get object name
      enet_name <- sprintf("enet_%d_%dvs%d", split_idx, level1, level2)
      
      # Check if object exists
      if (exists(enet_name, envir = .GlobalEnv)) {
        
        enet_obj <- get(enet_name, envir = .GlobalEnv)
        
        # Extract key information
        comp_results[[paste0("split_", split_idx)]] <- list(
          split = split_idx,
          comparison = sprintf("%d vs %d", level1, level2),
          levels = c(level1, level2),
          auc = enet_obj$metrics$auc,
          accuracy = enet_obj$metrics$accuracy,
          sensitivity = enet_obj$metrics$sensitivity,
          specificity = enet_obj$metrics$specificity,
          ppv = enet_obj$metrics$ppv,
          npv = enet_obj$metrics$npv,
          threshold = enet_obj$metrics$threshold,
          lambda_min = enet_obj$enet_results$lambda_final,
          n_predictors = nrow(enet_obj$pooled_coef$results_df),
          coefficients = enet_obj$pooled_coef$results_df,
          roc_object = enet_obj$metrics$roc_obj,
          confusion_matrix = enet_obj$metrics$confusion_matrix
        )
        
        cat(sprintf("  Split %d: AUC = %.3f, Accuracy = %.3f, %d predictors\n",
                    split_idx, 
                    enet_obj$metrics$auc,
                    enet_obj$metrics$accuracy,
                    nrow(enet_obj$pooled_coef$results_df)))
        
      } else {
        cat(sprintf("  Split %d: NOT FOUND\n", split_idx))
      }
    }
    
    # Store results for this comparison
    if (length(comp_results) > 0) {
      all_results[[sprintf("comp_%d_%dvs%d", comp_idx, level1, level2)]] <- comp_results
    }
  }
  
  # Check if we have any results
  if (length(all_results) == 0) {
    stop("No elastic net results found in environment!")
  }
  
  cat("\n\n=== Creating Summary Tables ===\n")
  
  # ============================================================================
  # 1. PERFORMANCE SUMMARY BY SPLIT
  # ============================================================================
  
  performance_by_split <- list()
  
  for (comp_name in names(all_results)) {
    comp_results <- all_results[[comp_name]]
    
    for (split_name in names(comp_results)) {
      split_result <- comp_results[[split_name]]
      
      performance_by_split[[paste(comp_name, split_name, sep = "_")]] <- data.frame(
        comparison = split_result$comparison,
        split = split_result$split,
        auc = split_result$auc,
        accuracy = split_result$accuracy,
        sensitivity = split_result$sensitivity,
        specificity = split_result$specificity,
        ppv = split_result$ppv,
        npv = split_result$npv,
        lambda_min = split_result$lambda_min,
        n_predictors = split_result$n_predictors,
        stringsAsFactors = FALSE
      ) %>%
        mutate(comparison = case_when(
          comparison == '1 vs 2' ~ 'Low Use (Reference) v Low Risk - Alcohol Use',
          comparison == '1 vs 3' ~ 'Low Use (Reference) v High Risk - Alcohol Use',
          comparison == '2 vs 3' ~ 'Low Risk - Alcohol Use (Reference) v High Risk - Alcohol Use',
          comparison == '1 vs 4' ~ 'Low Use (Reference) v Polysubstance Use',
          comparison == '2 vs 4' ~ 'Low Risk - Alcohol Use (Reference) v Polysubstance Use',
          comparison == '3 vs 4' ~ 'High Risk - Alcohol Use (Reference) v Polysubstance Use'
        ))
    }
  }
  
  performance_df <- bind_rows(performance_by_split) %>%
    arrange(comparison, split)
  
  cat("\n--- Performance by Split ---\n")
  print(head(performance_df, 20))
  
  write.csv(performance_df,
            "output/study4/en_results/performance_by_split_wMRI.csv",
            row.names = FALSE)
  
  # ============================================================================
  # 2. AGGREGATED PERFORMANCE SUMMARY (MEAN ACROSS SPLITS)
  # ============================================================================
  
  performance_summary <- performance_df %>%
    group_by(comparison) %>%
    summarise(
      n_splits = n(),
      mean_auc = mean(auc, na.rm = TRUE),
      sd_auc = sd(auc, na.rm = TRUE),
      min_auc = min(auc, na.rm = TRUE),
      max_auc = max(auc, na.rm = TRUE),
      mean_accuracy = mean(accuracy, na.rm = TRUE),
      sd_accuracy = sd(accuracy, na.rm = TRUE),
      mean_sensitivity = mean(sensitivity, na.rm = TRUE),
      sd_sensitivity = sd(sensitivity, na.rm = TRUE),
      mean_specificity = mean(specificity, na.rm = TRUE),
      sd_specificity = sd(specificity, na.rm = TRUE),
      mean_n_predictors = mean(n_predictors, na.rm = TRUE),
      sd_n_predictors = sd(n_predictors, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_auc))
  
  cat("\n--- Aggregated Performance Summary ---\n")
  print(performance_summary)
  
  write.csv(performance_summary,
            "output/study4/en_results/performance_summary_aggregated_wMRI.csv",
            row.names = FALSE)
  
  # ============================================================================
  # 3. TOP PREDICTORS BY COMPARISON (ACROSS ALL SPLITS)
  # ============================================================================
  
  cat("\n--- Extracting Coefficients ---\n")
  
  all_coefficients <- list()
  
  for (comp_name in names(all_results)) {
    comp_results <- all_results[[comp_name]]
    
    for (split_name in names(comp_results)) {
      split_result <- comp_results[[split_name]]
      
      if (nrow(split_result$coefficients) > 0) {
        coefs <- split_result$coefficients %>%
          mutate(
            comparison = split_result$comparison,
            split = split_result$split
          ) %>%
          select(comparison, split, variable, coefficient, stability_pct)
        
        all_coefficients[[paste(comp_name, split_name, sep = "_")]] <- coefs
      }
    }
  }
  
  all_coefs_df <- bind_rows(all_coefficients)
  
  cat(sprintf("Total coefficient entries: %d\n", nrow(all_coefs_df)))
  
  # Calculate variable importance: how often selected and average coefficient
  variable_importance <- all_coefs_df %>%
    mutate(comparison = case_when(
      comparison == '1 vs 2' ~ 'Low Use (Reference) v Low Risk - Alcohol Use',
      comparison == '1 vs 3' ~ 'Low Use (Reference) v High Risk - Alcohol Use',
      comparison == '2 vs 3' ~ 'Low Risk - Alcohol Use (Reference) v High Risk - Alcohol Use',
      comparison == '1 vs 4' ~ 'Low Use (Reference) v Polysubstance Use',
      comparison == '2 vs 4' ~ 'Low Risk - Alcohol Use (Reference) v Polysubstance Use',
      comparison == '3 vs 4' ~ 'High Risk - Alcohol Use (Reference) v Polysubstance Use'
    )) %>%
    group_by(comparison, variable) %>%
    summarise(
      n_splits_selected = n(),
      mean_coefficient = mean(coefficient, na.rm = TRUE),
      sd_coefficient = sd(coefficient, na.rm = TRUE),
      mean_stability = mean(stability_pct, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(abs_mean_coef = abs(mean_coefficient)) %>%
    arrange(comparison, desc(n_splits_selected), desc(abs_mean_coef))
  
  write.csv(variable_importance,
            "output/study4/en_results/variable_importance_by_comparison_wMRI.csv",
            row.names = FALSE)
  
  # Top predictors (selected in >=7 splits with high coefficients)
  top_predictors <- variable_importance %>%
    filter(n_splits_selected >= 7) %>%
    filter(variable != '(Intercept)') %>%
    group_by(comparison) %>%
    slice_max(order_by = abs_mean_coef, n = 20) %>%
    ungroup() %>%
    left_join(table_names_all[,c('variable', 'table_name', 'domain_name')], by = 'variable') %>%
    select(comparison, domain_name, table_name, variable, everything()) %>%
    arrange(comparison, rev(abs_mean_coef))
  
  cat("\n--- Top Predictors (selected in >=7 splits) ---\n")
  print(head(top_predictors, 30))
  
  write.csv(top_predictors,
            "output/study4/en_results/top_predictors_consistent_wMRI.csv",
            row.names = FALSE)
  
  
  # ============================================================================
  # 4. COEFFICIENT STABILITY HEATMAP DATA
  # ============================================================================
  
  # Create wide format: variables x comparisons, showing mean coefficient
  coef_wide <- variable_importance %>%
    filter(n_splits_selected >= 5) %>%  # Only variables selected in >=5 splits
    select(comparison, variable, mean_coefficient) %>%
    pivot_wider(
      names_from = comparison,
      values_from = mean_coefficient,
      values_fill = 0
    )
  
  write.csv(coef_wide,
            "output/study4/en_results/coefficients_wide_stable_wMRI.csv",
            row.names = FALSE)
  
  cat(sprintf("\nStable variables (selected in >=5 splits): %d\n", nrow(coef_wide)))
  
  # ============================================================================
  # 5. ROC CURVES VISUALIZATION
  # ============================================================================
  
  cat("\n--- Creating ROC Curves ---\n")
  
  roc_data_list <- list()
  
  for (comp_name in names(all_results)) {
    comp_results <- all_results[[comp_name]]
    
    for (split_name in names(comp_results)) {
      split_result <- comp_results[[split_name]]
      
      if (!is.null(split_result$roc_object)) {
        roc_obj <- split_result$roc_object
        
        roc_data_list[[paste(comp_name, split_name, sep = "_")]] <- data.frame(
          comparison = split_result$comparison,
          split = split_result$split,
          sensitivity = roc_obj$sensitivities,
          specificity = roc_obj$specificities,
          auc = split_result$auc,
          stringsAsFactors = FALSE
        ) %>%
          mutate(comparison = case_when(
            comparison == '1 vs 2' ~ 'Low Use (Reference) v Low Risk - Alcohol Use',
            comparison == '1 vs 3' ~ 'Low Use (Reference) v High Risk - Alcohol Use',
            comparison == '2 vs 3' ~ 'Low Risk - Alcohol Use (Reference) v High Risk - Alcohol Use',
            comparison == '1 vs 4' ~ 'Low Use (Reference) v Polysubstance Use',
            comparison == '2 vs 4' ~ 'Low Risk - Alcohol Use (Reference) v Polysubstance Use',
            comparison == '3 vs 4' ~ 'High Risk - Alcohol Use (Reference) v Polysubstance Use'
          ))
      }
    }
  }
  
  roc_data <- bind_rows(roc_data_list) %>%
    mutate(
      comparison_split = paste0(comparison, " (Split ", split, ")"),
      fpr = 1 - specificity
    )
  
  # Plot: One facet per comparison, showing all splits
  p_roc_facets <- ggplot(roc_data, 
                         aes(x = fpr, y = sensitivity, 
                             color = factor(split), group = comparison_split)) +
    geom_line(alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    facet_wrap(~ comparison, ncol = 2) +
    labs(
      title = "ROC Curves by Comparison (All Splits)",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = "Split"
    ) +
    theme_minimal() +
    theme(legend.position = 'none')
  
  ggsave("output/study4/en_results/roc_curves_by_comparison_all_splits_wMRI.png",
         p_roc_facets, width = 12, height = 10, dpi = 300)
  
  # Plot: Mean ROC curve per comparison
  mean_roc_data <- roc_data %>%
    group_by(comparison) %>%
    summarise(
      mean_auc = mean(auc),
      .groups = "drop"
    ) %>%
    distinct()
  
  # For mean curves, we need to interpolate at common FPR points
  fpr_seq <- seq(0, 1, by = 0.01)
  
  # Interpolate each split's ROC curve to common FPR grid
  mean_roc_curves <- roc_data %>%
    group_by(comparison, split) %>%
    arrange(fpr) %>%
    reframe(
      fpr_grid = fpr_seq,
      interp_sens = approx(fpr, sensitivity, xout = fpr_seq, rule = 2)$y
    ) %>%
    group_by(comparison, fpr_grid) %>%
    summarise(
      mean_sensitivity = mean(interp_sens, na.rm = TRUE),
      sd_sensitivity = sd(interp_sens, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(mean_roc_data, by = "comparison") %>%
    mutate(comparison_auc = sprintf("%s\n(Mean AUC = %.3f)", comparison, mean_auc))
  
  p_roc_mean <- ggplot(mean_roc_curves, 
                       aes(x = fpr_grid, y = mean_sensitivity, color = comparison_auc)) +
    geom_line(linewidth = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    labs(
      title = "Mean ROC Curves Across Splits",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = "Comparison"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 9)
    )
  
  ggsave("output/study4/en_results/roc_curves_mean_by_comparison_wMRI.png",
         p_roc_mean, width = 10, height = 8, dpi = 300)
  
  cat("ROC curves saved\n")
  
  # ============================================================================
  # 6. PERFORMANCE BOXPLOTS
  # ============================================================================
  
  cat("\n--- Creating Performance Boxplots ---\n")
  
  p_auc <- ggplot(performance_df, aes(x = comparison, y = auc, fill = comparison)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(
      title = "AUC Distribution Across Splits",
      x = "Comparison",
      y = "AUC"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    coord_cartesian(ylim = c(0.5, 1))
  
  ggsave("output/study4/en_results/auc_boxplot_by_comparison.png",
         p_auc, width = 10, height = 6, dpi = 300)
  
  p_accuracy <- ggplot(performance_df, aes(x = comparison, y = accuracy, fill = comparison)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(
      title = "Accuracy Distribution Across Splits",
      x = "Comparison",
      y = "Accuracy"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    coord_cartesian(ylim = c(0.5, 1))
  
  ggsave("output/study4/en_results/accuracy_boxplot_by_comparison_wMRI.png",
         p_accuracy, width = 10, height = 6, dpi = 300)
  
  cat("Performance boxplots saved\n")
  
  # ============================================================================
  # 7. VARIABLE SELECTION FREQUENCY
  # ============================================================================
  
  cat("\n--- Analyzing Variable Selection Frequency ---\n")
  
  var_selection_freq <- all_coefs_df %>%
    group_by(comparison, variable) %>%
    summarise(
      times_selected = n(),
      .groups = "drop"
    ) %>%
    arrange(comparison, desc(times_selected))
  
  # Variables selected in all 10 splits
  always_selected <- var_selection_freq %>%
    filter(times_selected == 10)
  
  cat(sprintf("Variables selected in all 10 splits: %d\n", nrow(always_selected)))
  
  write.csv(always_selected,
            "output/study4/en_results/variables_always_selected_wMRI.csv",
            row.names = FALSE)
  
  # Histogram of selection frequency
  p_freq <- ggplot(var_selection_freq %>%
                     mutate(comparison = case_when(
                       comparison == '1 vs 2' ~ 'Low Use (Reference) v Low Risk - Alcohol Use',
                       comparison == '1 vs 3' ~ 'Low Use (Reference) v High Risk - Alcohol Use',
                       comparison == '2 vs 3' ~ 'Low Risk - Alcohol Use (Reference) v High Risk - Alcohol Use',
                       comparison == '1 vs 4' ~ 'Low Use (Reference) v Polysubstance Use',
                       comparison == '2 vs 4' ~ 'Low Risk - Alcohol Use (Reference) v Polysubstance Use',
                       comparison == '3 vs 4' ~ 'High Risk - Alcohol Use (Reference) v Polysubstance Use'
                     )), 
                   aes(x = times_selected, fill = comparison)) +
    geom_histogram(binwidth = 1, color = "black", alpha = 0.7) +
    facet_wrap(~ comparison, ncol = 2) +
    labs(
      title = "Variable Selection Frequency Across Splits",
      x = "Number of Splits Variable Was Selected",
      y = "Count of Variables"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  ggsave("output/study4/en_results/variable_selection_frequency_wMRI.png",
         p_freq, width = 10, height = 8, dpi = 300)
  
  cat("Variable selection frequency saved\n")
  
  # ============================================================================
  # 8. CONFUSION MATRIX SUMMARY
  # ============================================================================
  
  cat("\n--- Summarizing Confusion Matrices ---\n")
  
  confusion_summary <- list()
  
  for (comp_name in names(all_results)) {
    comp_results <- all_results[[comp_name]]
    
    for (split_name in names(comp_results)) {
      split_result <- comp_results[[split_name]]
      
      cm <- split_result$confusion_matrix
      
      confusion_summary[[paste(comp_name, split_name, sep = "_")]] <- data.frame(
        comparison = split_result$comparison,
        split = split_result$split,
        true_neg = cm[1, 1],
        false_pos = cm[2, 1],
        false_neg = cm[1, 2],
        true_pos = cm[2, 2],
        stringsAsFactors = FALSE
      ) %>%
        mutate(comparison = case_when(
          comparison == '1 vs 2' ~ 'Low Use (Reference) v Low Risk - Alcohol Use',
          comparison == '1 vs 3' ~ 'Low Use (Reference) v High Risk - Alcohol Use',
          comparison == '2 vs 3' ~ 'Low Risk - Alcohol Use (Reference) v High Risk - Alcohol Use',
          comparison == '1 vs 4' ~ 'Low Use (Reference) v Polysubstance Use',
          comparison == '2 vs 4' ~ 'Low Risk - Alcohol Use (Reference) v Polysubstance Use',
          comparison == '3 vs 4' ~ 'High Risk - Alcohol Use (Reference) v Polysubstance Use'
        ))
    }
  }
  
  confusion_df <- bind_rows(confusion_summary)
  
  write.csv(confusion_df,
            "output/study4/en_results/confusion_matrices_all_splits_wMRI.csv",
            row.names = FALSE)
  
  # ============================================================================
  # FINAL SUMMARY
  # ============================================================================
  
  cat("\n\n")
  cat("================================================================================\n")
  cat("=== ELASTIC NET ANALYSIS COMPLETE ===\n")
  cat("================================================================================\n\n")
  
  cat("Files saved to output/study4/en_results_claude/:\n\n")
  cat("Performance Metrics:\n")
  cat("  - performance_by_split.csv\n")
  cat("  - performance_summary_aggregated.csv\n")
  cat("  - confusion_matrices_all_splits.csv\n\n")
  
  cat("Variable/Coefficient Information:\n")
  cat("  - variable_importance_by_comparison.csv\n")
  cat("  - top_predictors_consistent.csv\n")
  cat("  - coefficients_wide_stable.csv\n")
  cat("  - variables_always_selected.csv\n\n")
  
  cat("Visualizations:\n")
  cat("  - roc_curves_by_comparison_all_splits.png\n")
  cat("  - roc_curves_mean_by_comparison.png\n")
  cat("  - auc_boxplot_by_comparison.png\n")
  cat("  - accuracy_boxplot_by_comparison.png\n")
  cat("  - variable_selection_frequency.png\n\n")
  
  cat("Summary Statistics:\n")
  cat(sprintf("  Total comparisons analyzed: %d\n", length(all_results)))
  cat(sprintf("  Total splits analyzed: %d\n", nrow(performance_df)))
  cat(sprintf("  Mean AUC (across all): %.3f (SD: %.3f)\n", 
              mean(performance_df$auc, na.rm = TRUE),
              sd(performance_df$auc, na.rm = TRUE)))
  cat(sprintf("  Mean Accuracy (across all): %.3f (SD: %.3f)\n",
              mean(performance_df$accuracy, na.rm = TRUE),
              sd(performance_df$accuracy, na.rm = TRUE)))
  cat(sprintf("  Total unique variables selected: %d\n",
              length(unique(all_coefs_df$variable))))
  cat(sprintf("  Variables selected in all 10 splits: %d\n",
              nrow(always_selected)))
  
  cat("\n")
  cat("Best performing comparison (by mean AUC):\n")
  best_comp <- performance_summary %>% slice_max(mean_auc, n = 1)
  cat(sprintf("  %s: Mean AUC = %.3f (SD: %.3f)\n",
              best_comp$comparison, best_comp$mean_auc, best_comp$sd_auc))
  
  email <- compose_email(
    body = "eNet complete."
  )
  
  smtp_send(
    email,
    from = " email here",
    to   = "destination email here",
    subject = "Four Class Model Complete",
    credentials = creds_key("gmail_smtp")
  )
  
  
  
  rm(list = ls())
  
  
  
  
}, error = function(e) {
  # Send email on error
  smtp_send(
    compose_email(
      body = paste("Error occurred:\n\n", e$message, "\n\nTraceback:\n",
                   paste(capture.output(traceback()), collapse = "\n"))
    ),
    from = "your email here",
    to   = "destination email here",
    subject = "Error Occurred",
    credentials = creds_key("gmail_smtp")
  )
  
  stop(e)  # Re-throw the error if you want it to still fail
})

rm(list = ls())
system("powercfg -change standby-timeout-ac 15")  # Resets sleep timeout
system('shutdown -l')

