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
