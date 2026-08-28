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

