################################################################################
################################################################################
################### Study 3 Code ###############################################
################################################################################


# Load Packages -----------------------------------------------------------
library(mediation)
library(dplyr)
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
library(recipes)
library(ggplot2) 
library(pROC)
library(diffdf)
library(data.table)
library(moments)
library(lmtest)
library(gtsummary)
library(arsenal) 
library(rstatix)
library(cvAUC)
library(svMisc)
library(reshape2)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("sva")
BiocManager::install("neuroconductor")
library(sva)
#library(summarytools)
library(patchwork)
library(lavaan)
library(pscl)
library(broom) 
library(fastDummies)
library(flextable)
library(officer)
library(showtext)
library("ggpubr")
library(stringr)
library("plotrix")
library(gridExtra)
library(magick)
library(showtext)
library(ggeffects)
library(openxlsx)
library(grid)
library(purrr)
library(ggtext)
library(rtools)

showtext_auto()
# Times New Roman is usually installed as "Times New Roman"
font_add(
  "Times New Roman",
  regular    = "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
  bold       = "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf",
  italic     = "/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf",
  bolditalic = "/System/Library/Fonts/Supplemental/Times New Roman Bold Italic.ttf"
)
font_add("Times New Roman", regular = "C:/WINDOWS/Fonts/times.ttf")

showtext_opts(dpi = 320)

# Create final dataset based on inclusion criteria & define outcome ------------

# import substance use data (baseline and annual follow-up)
su_y_sui <- read_csv('data/5.1/su_y_sui.csv',
                     col_types = list(tlfb_blunt_use = col_double(),
                                      tlfb_mdma_use = col_double(),
                                      tlfb_ket_use = col_double())) %>% 
  type_convert()

# subset baseline substance use variables 
su_base <- su_y_sui %>% 
  subset(eventname == 'baseline_year_1_arm_1') %>% 
  select(
    src_subject_id,
    tlfb_age_calc_inmonths,
    # primary substance use variables 
    isip_1b_yn,
    tlfb_tob_puff, tlfb_chew_use, tlfb_cigar_use, tlfb_hookah_use,
    tlfb_pipes_use, tlfb_nicotine_use,
    tlfb_mj_puff, tlfb_blunt_use, tlfb_mj_conc_use, tlfb_mj_drink_use,
    tlfb_tincture_use, 
    tlfb_mj_synth_use, tlfb_coc_use, tlfb_bsalts_use, tlfb_meth_use, 
    tlfb_mdma_use, tlfb_ket_use, tlfb_ghb_use, tlfb_opi_use, tlfb_hall_use, 
    tlfb_shrooms_use, tlfb_salvia_use, tlfb_steroids_use, tlfb_sniff_use, 
    tlfb_amp_use, tlfb_tranq_use, tlfb_vicodin_use, tlfb_cough_use, 
    tlfb_other_use,
    # gating criteria variables
    tlfb_alc, tlfb_alc_sip, tlfb_tob, tlfb_mj, tlfb_mj_synth, 
    tlfb_list_yes_no, tlfb_list___1:tlfb_list___12, tlfb_inhalant, 
    tlfb_rx_misuse) 

# add suffix to indicate timepoint for baseline variables 
colnames(su_base) <- paste(
  colnames(su_base), 'base', sep = '.')
su_base <- su_base %>% 
  rename(src_subject_id = src_subject_id.base)

# subset annual follow-up substance use variables 
su_yr <- su_y_sui %>% 
  arrange(src_subject_id) %>% 
  subset(
    eventname == '1_year_follow_up_y_arm_1' |
      eventname == '2_year_follow_up_y_arm_1' | 
      eventname == '3_year_follow_up_y_arm_1') %>%  
  mutate(
    eventname = 
      case_when(
        eventname == '1_year_follow_up_y_arm_1' ~ '1_year',
        eventname == '2_year_follow_up_y_arm_1' ~ '2_year',
        eventname == '3_year_follow_up_y_arm_1' ~ '3_year')) %>% 
  select(
    src_subject_id, eventname,
    # primary substance use variables 
    isip_1b_yn_l, 
    tlfb_tob_puff_l, tlfb_chew_use_l, tlfb_hookah_use_l, tlfb_nicotine_use_l, 
    tlfb_mj_puff_l, tlfb_blunt_use_l, tlfb_mj_conc_use_l, tlfb_mj_drink_use_l, 
    tlfb_tincture_use_l, 
    tlfb_mj_synth_use_l, tlfb_coc_use_l, tlfb_bsalts_use_l, 
    tlfb_meth_use_l, tlfb_mdma_use_l, tlfb_ket_use_l, tlfb_ghb_use_l, 
    tlfb_opi_use_l, tlfb_lsd_use_l, tlfb_shrooms_use_l, tlfb_salvia_use_l, 
    tlfb_steroids_use_l, tlfb_inhalant_use_l, 
    tlfb_amp_use_l, tlfb_tranq_use_l, 
    tlfb_vicodin_use_l, tlfb_cough_use_l, tlfb_other_use_l,
    # gating criteria
    tlfb_alc_l, tlfb_alc_sip_l, tlfb_tob_l, 
    tlfb_mj_l, tlfb_mj_synth_l, tlfb_list_yes_no_l, tlfb_list_l___1, 
    tlfb_list_l___2, tlfb_list_l___3, tlfb_list_l___4, tlfb_list_l___5, 
    tlfb_list_l___6, tlfb_list_l___7, tlfb_list_l___8, tlfb_list_l___9, 
    tlfb_list_l___10, tlfb_list_l___11, tlfb_list_l___12, tlfb_inhalant_l, 
    tlfb_rx_misuse_l)

# import mid-year substance use data 
su_y_mypi <- read_delim('data/5.1//su_y_mypi.csv') %>% 
  type_convert()

# subset mid-year substance use variables 
su_mid <- su_y_mypi %>% 
  subset(
    eventname == '6_month_follow_up_arm_1' | 
      eventname == '18_month_follow_up_arm_1' |
      eventname == '30_month_follow_up_arm_1') %>% 
  mutate(
    eventname = 
      case_when(
        eventname == '6_month_follow_up_arm_1' ~ '6_month',
        eventname == '18_month_follow_up_arm_1' ~ '18_month',
        eventname == '30_month_follow_up_arm_1' ~ '30_month')) %>% 
  select(
    src_subject_id, eventname,
    # primary substance use variables 
    mypi_alc_sip_1b, 
    mypi_ecig, mypi_cigar_used, mypi_flavoring, mypi_chew_pst_used, 
    mypi_mj_used, mypi_mj_edible, mypi_mj_oils, mypi_mj_tinc_used, mypi_mj_vape,
    mypi_mj_oils_vaped, 
    mypi_mj_synth_used, mypi_coke_used, mypi_meth_used, mypi_ghb_used, 
    mypi_heroin_used, mypi_sniff_used, mypi_pills_used, mypi_pills_dep_used, 
    mypi_pr_used, mypi_cold_used, mypi_high_other_used,
    # gating criteria
    mypi_alc, mypi_alc_full_drink, mypi_alc_sip,
    mypi_tob, mypi_tob_used, mypi_chew, mypi_mj, mypi_mj_30, mypi_sniff, mypi_pills, 
    mypi_high_other)

# transpose annual and mid-year substance use data from long to wide 
su_yr_wide <- su_yr %>%
  pivot_wider(id_cols = 'src_subject_id',
              names_from = 'eventname', 
              names_sep = ".",
              values_from = c(
                # primary substance use variables
                isip_1b_yn_l, 
                tlfb_tob_puff_l, tlfb_chew_use_l, tlfb_hookah_use_l, 
                tlfb_nicotine_use_l, 
                tlfb_mj_puff_l, tlfb_blunt_use_l, tlfb_mj_conc_use_l, 
                tlfb_mj_drink_use_l, tlfb_tincture_use_l, 
                tlfb_mj_synth_use_l, tlfb_coc_use_l, tlfb_bsalts_use_l, 
                tlfb_meth_use_l, tlfb_mdma_use_l, tlfb_ket_use_l, 
                tlfb_ghb_use_l, tlfb_opi_use_l, tlfb_lsd_use_l, 
                tlfb_shrooms_use_l, tlfb_salvia_use_l, tlfb_steroids_use_l, 
                tlfb_inhalant_use_l, tlfb_amp_use_l, tlfb_tranq_use_l, 
                tlfb_vicodin_use_l, tlfb_cough_use_l, tlfb_other_use_l,
                # gating criteria variables 
                tlfb_alc_l, tlfb_alc_sip_l, tlfb_tob_l, tlfb_mj_l, 
                tlfb_mj_synth_l, tlfb_list_yes_no_l,  tlfb_inhalant_l, 
                tlfb_rx_misuse_l))

su_mid_wide <- su_mid %>%
  pivot_wider(id_cols = 'src_subject_id',
              names_from = 'eventname', 
              names_sep = ".",
              values_from = c(
                # primary substance use variables
                mypi_alc_sip_1b, 
                mypi_ecig, mypi_cigar_used, mypi_flavoring, mypi_chew_pst_used, 
                mypi_mj_used, mypi_mj_edible, mypi_mj_oils, mypi_mj_tinc_used, 
                mypi_mj_vape, mypi_mj_oils_vaped, 
                mypi_mj_synth_used, mypi_coke_used, mypi_meth_used, 
                mypi_ghb_used, mypi_heroin_used, mypi_sniff_used, 
                mypi_pills_used, mypi_pills_dep_used, mypi_pr_used, 
                mypi_cold_used, mypi_high_other_used,
                # gating criteria variables
                mypi_alc, mypi_alc_full_drink, mypi_alc_sip, mypi_tob, 
                mypi_tob_used, mypi_chew, mypi_mj, mypi_mj_30, mypi_sniff, 
                mypi_pills, mypi_high_other))

# merge all substance use data 
data.1 <- su_base %>%
  full_join(su_yr_wide, by = 'src_subject_id') %>% 
  full_join(su_mid_wide, by = 'src_subject_id') 

names(data.1)

# drop objects no longer needed
rm(su_mid, su_mid_wide, su_yr, su_yr_wide)


###### Inclusion Criteria #1 ###################################################

# (1) criteria: 
# - remove those who have insufficient substance use data (i.e., missing gating 
#   criteria heard of variables) at baseline and follow-up timepoints

data.1 <- data.1 %>% 
  mutate(
    exclude_1 = ifelse(
      # baseline
      is.na(tlfb_alc.base) & is.na(tlfb_tob.base) & is.na(tlfb_mj.base) &
        is.na(tlfb_mj_synth.base) & is.na(tlfb_list_yes_no.base) &
        is.na(tlfb_rx_misuse.base) & is.na(tlfb_inhalant.base) &
        
        # 6-month
        is.na(mypi_alc.6_month) & is.na(mypi_tob.6_month) & 
        is.na(mypi_chew.6_month) & is.na(mypi_mj.6_month) & 
        is.na(mypi_sniff.6_month) & is.na(mypi_pills.6_month) &
        is.na(mypi_high_other.6_month) &
        
        # 18-month
        is.na(mypi_alc.18_month) & is.na(mypi_tob.18_month) & 
        is.na(mypi_chew.18_month) & is.na(mypi_mj.18_month) & 
        is.na(mypi_sniff.18_month) & is.na(mypi_pills.18_month) & 
        is.na(mypi_high_other.18_month) &
        
        # 30-month 
        is.na(mypi_alc.30_month) & is.na(mypi_tob.30_month) & 
        is.na(mypi_chew.30_month) & is.na(mypi_mj.30_month) & 
        is.na(mypi_sniff.30_month) & is.na(mypi_pills.30_month) &
        is.na(mypi_high_other.30_month) &
        
        # 1-year
        is.na(tlfb_alc_l.1_year) & is.na(tlfb_tob_l.1_year) &
        is.na(tlfb_mj_l.1_year) & is.na(tlfb_mj_synth_l.1_year) & 
        is.na(tlfb_list_yes_no_l.1_year) & is.na(tlfb_inhalant_l.1_year) & 
        is.na(tlfb_rx_misuse_l.1_year) & 
        
        # 2-year
        is.na(tlfb_alc_l.2_year) & is.na(tlfb_tob_l.2_year) &
        is.na(tlfb_mj_l.2_year) & is.na(tlfb_mj_synth_l.2_year) & 
        is.na(tlfb_list_yes_no_l.2_year) & is.na(tlfb_inhalant_l.2_year) & 
        is.na(tlfb_rx_misuse_l.2_year) & 
        
        # 3-year  
        is.na(tlfb_alc_l.3_year) & is.na(tlfb_tob_l.3_year) &
        is.na(tlfb_mj_l.3_year) & is.na(tlfb_mj_synth_l.3_year) & 
        is.na(tlfb_list_yes_no_l.3_year) & is.na(tlfb_inhalant_l.3_year) & 
        is.na(tlfb_rx_misuse_l.3_year),
      'Drop', 'Keep'))

table(data.1$exclude_1, useNA = 'ifany') 
# n = 1 missing gating at baseline and follow-up

data.exclude <- data.1 %>%
  subset(exclude_1 != 'Keep') %>%
  select(src_subject_id, exclude_1) %>%
  mutate(exclude_1 = ifelse(exclude_1 == 'Drop', TRUE, FALSE))

data.1 <- data.1 %>% 
  subset(exclude_1 == 'Keep')

###### expected sample size -----#
# n = 11867 (11868 - 1)
###### expected sample size -----#

# (2) criteria: 
# - remove those who are missing baseline but have follow-up gating, on the 
#   basis that we cannot establish initial use at baseline if they are missing
#   gating criteria

data.1 <- data.1 %>% 
  mutate(
    exclude_2 = ifelse(
      is.na(tlfb_alc.base) & is.na(tlfb_tob.base) & is.na(tlfb_mj.base) &
        is.na(tlfb_mj_synth.base) & is.na(tlfb_list_yes_no.base) &
        is.na(tlfb_rx_misuse.base) & is.na(tlfb_inhalant.base), 'Drop', 'Keep'))

table(data.1$exclude_2, useNA = 'ifany') 
# n = 8 to drop

###### double check ---#
# confirm no baseline substance use data present

# indicator of missing baseline data
data.1 <- data.1 %>% 
  mutate(
    base_NA = ifelse(
      is.na(isip_1b_yn.base) & 
        is.na(tlfb_tob_puff.base) & is.na(tlfb_chew_use.base) & 
        is.na(tlfb_cigar_use.base) & is.na(tlfb_hookah_use.base) & 
        is.na(tlfb_pipes_use.base) & is.na(tlfb_nicotine_use.base) & 
        is.na(tlfb_mj_puff.base) & is.na(tlfb_blunt_use.base) & 
        is.na(tlfb_mj_conc_use.base) & is.na(tlfb_tincture_use.base) & 
        is.na(tlfb_mj_synth_use.base) & is.na(tlfb_coc_use.base) & 
        is.na(tlfb_bsalts_use.base) & is.na(tlfb_meth_use.base) & 
        is.na(tlfb_mdma_use.base) & is.na(tlfb_ket_use.base) & 
        is.na(tlfb_ghb_use.base) & is.na(tlfb_opi_use.base) & 
        is.na(tlfb_hall_use.base) & is.na(tlfb_shrooms_use.base) & 
        is.na(tlfb_salvia_use.base) & is.na(tlfb_steroids_use.base) & 
        is.na(tlfb_sniff_use.base) & is.na(tlfb_amp_use.base) &
        is.na(tlfb_tranq_use.base) & is.na(tlfb_vicodin_use.base) & 
        is.na(tlfb_cough_use.base) & is.na(tlfb_other_use.base), 
      'Missing All Baseline', 'Not Missing All Baseline'))

data.1 <- data.1 %>% 
  mutate(
    exclude_2_rev = case_when(
      exclude_2 == 'Keep' ~ 'Keep',
      exclude_2 == 'Drop' & base_NA == 'Missing All Baseline' ~ 
        'Gating Missing & SU Data Missing ',
      exclude_2 == 'Drop' & base_NA == 'Not Missing All Baseline' ~ 
        'Gating Missing & SU Data Present'))

table(data.1$exclude_2_rev, useNA = 'ifany') 
# n = 2 to retain due to SU data present
###### double check ---#


data.exclude <- data.exclude %>%
  full_join(data.1 %>%  
              subset(exclude_2_rev != 'Keep' & 
                       exclude_2_rev != 'Gating Missing & SU Data Present') %>% 
              mutate(exclude_2_rev = TRUE) %>% select(src_subject_id, exclude_2_rev) %>%
              rename(exclude_2 = exclude_2_rev),
            by = 'src_subject_id')


data.1 <- data.1 %>% 
  subset(exclude_2_rev == 'Keep' | 
           exclude_2_rev == 'Gating Missing & SU Data Present')

###### expected sample size -----#
# n = 11861 (11867 - 6)
###### expected sample size -----#

# create summary variables for substances w/multiple sub-categories
data.1 <- data.1 %>% 
  mutate(
    tlfb_nicotine.base = ifelse(
      (tlfb_tob_puff.base == 1 | tlfb_chew_use.base == 1 | 
         tlfb_cigar_use.base == 1 | tlfb_hookah_use.base == 1 | 
         tlfb_pipes_use.base == 1 | tlfb_nicotine_use.base == 1), 1, 0),
    tlfb_cannabis.base = ifelse(
      (tlfb_mj_puff.base == 1 | tlfb_blunt_use.base == 1 | 
         tlfb_mj_conc_use.base == 1 | tlfb_mj_drink_use.base == 1 | 
         tlfb_tincture_use.base == 1), 1, 0),
    tlfb_rx.base = ifelse(
      (tlfb_amp_use.base == 1 | tlfb_tranq_use.base == 1 | 
         tlfb_vicodin_use.base == 1 | tlfb_cough_use.base == 1), 1, 0),
    tlfb_hall.base = 
      ifelse(
        (tlfb_hall_use.base == 1 | tlfb_shrooms_use.base == 1 | 
           tlfb_salvia_use.base == 1), 1, 0)) 

# check for any cases where indicated sip but missing religious context
table_alc.base <- table(
  data.1$tlfb_alc_sip.base, data.1$isip_1b_yn.base, useNA = 'ifany')
table_alc.1y <- table(
  data.1$tlfb_alc_sip_l.1_year, data.1$isip_1b_yn_l.1_year, useNA = 'ifany')
table_alc.2y <- table(
  data.1$tlfb_alc_sip_l.2_year, data.1$isip_1b_yn_l.2_year, useNA = 'ifany')
table_alc.3y <- table(
  data.1$tlfb_alc_sip_l.3_year, data.1$isip_1b_yn_l.3_year, useNA = 'ifany')

table_alc.6m <- table(
  data.1$mypi_alc_sip.6_month, data.1$mypi_alc_sip_1b.6_month, useNA = 'ifany')
table_alc.18m <- table(
  data.1$mypi_alc_sip.18_month, data.1$mypi_alc_sip_1b.18_month, useNA = 'ifany')
table_alc.30m <- table(
  data.1$mypi_alc_sip.30_month, data.1$mypi_alc_sip_1b.30_month, useNA = 'ifany')

names(dimnames(table_alc.base)) <- c('tlfb_alc_sip', 'isip_1b_yn')
names(dimnames(table_alc.1y)) <- c('tlfb_alc_sip_l', 'isip_1b_yn_l')
names(dimnames(table_alc.2y)) <- c('tlfb_alc_sip_l', 'isip_1b_yn_l')
names(dimnames(table_alc.3y)) <- c('tlfb_alc_sip_l', 'isip_1b_yn_l')
names(dimnames(table_alc.6m)) <- c('mypi_alc_sip', 'mypi_alc_sip_1b')
names(dimnames(table_alc.18m)) <- c('mypi_alc_sip', 'mypi_alc_sip_1b')
names(dimnames(table_alc.30m)) <- c('mypi_alc_sip', 'mypi_alc_sip_1b')

table_alc.base 
table_alc.1y
table_alc.2y
table_alc.3y
table_alc.6m
table_alc.18m
table_alc.30m
# no cases where indicated sip but missing religious context

rm(table_alc.base, table_alc.1y, table_alc.2y, table_alc.3y,
   table_alc.6m, table_alc.18m, table_alc.30m)

# subset participants w/no substance use at baseline
data.1 <- data.1 %>% 
  mutate(
    .cols = if_any(
      isip_1b_yn.base | tlfb_tob_puff.base |
        tlfb_chew_use.base | tlfb_cigar_use.base | tlfb_hookah_use.base |
        tlfb_pipes_use.base | tlfb_nicotine_use.base |
        tlfb_mj_puff.base | tlfb_blunt_use.base | tlfb_mj_conc_use.base | 
        tlfb_mj_drink_use.base | tlfb_tincture_use.base |
        tlfb_mj_synth_use.base | tlfb_coc_use.base | tlfb_bsalts_use.base | 
        tlfb_meth_use.base | tlfb_mdma_use.base | tlfb_ket_use.base | 
        tlfb_ghb_use.base | tlfb_opi_use.base | tlfb_hall_use.base | 
        tlfb_shrooms_use.base | tlfb_salvia_use.base | tlfb_steroids_use.base |
        tlfb_sniff_use.base | tlfb_amp_use.base | tlfb_tranq_use.base | 
        tlfb_vicodin_use.base | tlfb_cough_use.base | tlfb_other_use.base,
      .fns = ~ifelse(.x == 1, 1, 0)),
    .names = NULL) %>% 
  rename(su_base_use = .cols)

table(data.1$su_base_use, useNA = 'ifany') # TRUE = 2110

data.exclude <- data.exclude %>%
  full_join(data.1 %>% select(src_subject_id, su_base_use) %>%
              subset(su_base_use == TRUE), by = 'src_subject_id') %>%
  rename(exclude_3 = su_base_use)

data.1 <- data.1 %>%
  subset(is.na(su_base_use)) %>%
  mutate(su_base_use = case_when(
    is.na(su_base_use) ~ FALSE,
    su_base_use == TRUE ~ TRUE
  )) %>%
  rename(DV_BL = su_base_use)

###### expected sample size -----#
# n = 11861 - 2110 = 9751
###### expected sample size -----#


# recode 'refuse to answer' responses at mid-year assessment to NA
data.1 <- data.1 %>%  
  mutate( 
    across(.cols = starts_with('mypi_') & ends_with('.6_month'), 
           .fns = ~replace(.x, .x == 777, NA))) %>% 
  mutate( 
    across(.cols = starts_with('mypi_') & ends_with('.18_month'), 
           .fns = ~replace(.x, .x == 777, NA))) %>% 
  mutate( 
    across(.cols = starts_with('mypi_') & ends_with('.30_month'), 
           .fns = ~replace(.x, .x == 777, NA))) 


###### Inclusion Criteria #2 ###################################################



# check number missing age (& substance use variables of interest at 3-year) 
# from age variable in longitudinal TLFB questionnaire (su_y_sui)

criteria_3 <- su_y_sui %>% 
  subset(eventname == '3_year_follow_up_y_arm_1') %>% 
  mutate(age_3yr = tlfb_age_l) %>% 
  select(src_subject_id, age_3yr)

data.1 <- data.1 %>% 
  left_join(criteria_3, by = 'src_subject_id')

table(data.1$age_3yr, useNA = 'ifany') 
# n = 1330 missing age at 3-year

temp <- data.1 %>% 
  subset(is.na(age_3yr)) %>% 
  select(
    isip_1b_yn_l.3_year, 
    tlfb_tob_puff_l.3_year, tlfb_chew_use_l.3_year, tlfb_hookah_use_l.3_year, 
    tlfb_nicotine_use_l.3_year, 
    tlfb_mj_puff_l.3_year, tlfb_blunt_use_l.3_year, tlfb_mj_conc_use_l.3_year, 
    tlfb_mj_drink_use_l.3_year, tlfb_tincture_use_l.3_year, 
    tlfb_mj_synth_use_l.3_year, tlfb_coc_use_l.3_year, tlfb_bsalts_use_l.3_year, 
    tlfb_meth_use_l.3_year, tlfb_mdma_use_l.3_year, tlfb_ket_use_l.3_year, 
    tlfb_ghb_use_l.3_year, tlfb_opi_use_l.3_year, tlfb_lsd_use_l.3_year, 
    tlfb_shrooms_use_l.3_year, tlfb_salvia_use_l.3_year, 
    tlfb_steroids_use_l.3_year, tlfb_inhalant_use_l.3_year, 
    tlfb_amp_use_l.3_year, tlfb_tranq_use_l.3_year, tlfb_vicodin_use_l.3_year, 
    tlfb_cough_use_l.3_year, tlfb_other_use_l.3_year)

na_count <- temp %>% 
  summarise(
    across(everything(), ~ sum(is.na(.)))) %>% 
  t() 

# - of the n = 1330 missing age at the 3-year follow-up (majority missing most 
#   substance use variables, n's = 1329 or 1330)

rm(temp, na_count, criteria_3)

data.exclude <- data.exclude %>%
  full_join(data.1 %>% select(src_subject_id, age_3yr) %>%
              subset(is.na(age_3yr)) %>%
              rename(exclude_4 = age_3yr) %>%
              mutate(exclude_4 = TRUE),
            by = 'src_subject_id')

# subset participants >= age 12
data.1 <- data.1 %>%
  subset(!is.na(age_3yr)) 
# n = 1330 exclude: missing age (& substance use variables of interest at 3-year)

###### expected sample size -----#
# n = 841 (9751 - 1330)
###### expected sample size -----#

# recode 'refuse to answer' responses at mid-year assessment to NA
data.1 <- data.1 %>%  
  mutate( 
    across(.cols = starts_with('mypi_') & ends_with('.6_month'), 
           .fns = ~replace(.x, .x == 777, NA))) %>% 
  mutate( 
    across(.cols = starts_with('mypi_') & ends_with('.18_month'), 
           .fns = ~replace(.x, .x == 777, NA))) %>% 
  mutate( 
    across(.cols = starts_with('mypi_') & ends_with('.30_month'), 
           .fns = ~replace(.x, .x == 777, NA))) 


###### Inclusion Criteria #4a                               #


# criteria: 
# - substance use data present for at least 1 follow-up timepoint, exclude if 
#   missing substance use variables across all timepoints 

# subset gating criteria heard of use variables 
su_gating <- data.1 %>% 
  select(
    src_subject_id,
    
    # base
    tlfb_alc.base, tlfb_tob.base, tlfb_mj.base, tlfb_mj_synth.base, 
    tlfb_list_yes_no.base, tlfb_inhalant.base, tlfb_rx_misuse.base,  
    
    # annual
    tlfb_alc_l.1_year, tlfb_alc_l.2_year, tlfb_alc_l.3_year, 
    tlfb_tob_l.1_year, tlfb_tob_l.2_year, tlfb_tob_l.3_year, 
    tlfb_mj_l.1_year, tlfb_mj_l.2_year, tlfb_mj_l.3_year, 
    tlfb_mj_synth_l.1_year, tlfb_mj_synth_l.2_year, tlfb_mj_synth_l.3_year, 
    tlfb_list_yes_no_l.1_year, tlfb_list_yes_no_l.2_year, 
    tlfb_list_yes_no_l.3_year , 
    tlfb_inhalant_l.1_year, tlfb_inhalant_l.2_year, tlfb_inhalant_l.3_year, 
    tlfb_rx_misuse_l.1_year, tlfb_rx_misuse_l.2_year, tlfb_rx_misuse_l.3_year, 
    
    # mid-year
    mypi_alc.6_month, mypi_alc.18_month, mypi_alc.30_month, 
    mypi_tob.6_month, mypi_tob.18_month, mypi_tob.30_month,  
    mypi_chew.6_month, mypi_chew.18_month, mypi_chew.30_month, 
    mypi_mj.6_month, mypi_mj.18_month, mypi_mj.30_month,
    mypi_sniff.6_month, mypi_sniff.18_month, mypi_sniff.30_month, 
    mypi_pills.6_month, mypi_pills.18_month, mypi_pills.30_month, 
    mypi_high_other.6_month, mypi_high_other.18_month, mypi_high_other.30_month)

# subset primary substance use variables 
data.1 <- data.1 %>% 
  select(
    src_subject_id,
    tlfb_age_calc_inmonths.base,
    
    # baseline (n = 30)
    isip_1b_yn.base,
    tlfb_tob_puff.base, tlfb_chew_use.base, tlfb_cigar_use.base, 
    tlfb_hookah_use.base, tlfb_pipes_use.base, tlfb_nicotine_use.base,
    tlfb_mj_puff.base, tlfb_blunt_use.base, tlfb_mj_conc_use.base, 
    tlfb_mj_drink_use.base, tlfb_tincture_use.base, 
    tlfb_mj_synth_use.base, tlfb_coc_use.base, tlfb_bsalts_use.base, 
    tlfb_meth_use.base, tlfb_mdma_use.base, tlfb_ket_use.base, 
    tlfb_ghb_use.base, tlfb_opi_use.base, tlfb_hall_use.base, 
    tlfb_shrooms_use.base, tlfb_salvia_use.base, tlfb_steroids_use.base, 
    tlfb_sniff_use.base, tlfb_amp_use.base, tlfb_tranq_use.base, 
    tlfb_vicodin_use.base, tlfb_cough_use.base, tlfb_other_use.base,
    
    # annual follow-up (n = 28)
    isip_1b_yn_l.1_year, isip_1b_yn_l.2_year, isip_1b_yn_l.3_year, 
    tlfb_tob_puff_l.1_year, tlfb_tob_puff_l.2_year, tlfb_tob_puff_l.3_year, 
    tlfb_chew_use_l.1_year, tlfb_chew_use_l.2_year, tlfb_chew_use_l.3_year, 
    tlfb_hookah_use_l.1_year, tlfb_hookah_use_l.2_year, tlfb_hookah_use_l.3_year, 
    tlfb_nicotine_use_l.1_year, tlfb_nicotine_use_l.2_year, tlfb_nicotine_use_l.3_year, 
    tlfb_mj_puff_l.1_year, tlfb_mj_puff_l.2_year, tlfb_mj_puff_l.3_year, 
    tlfb_blunt_use_l.1_year, tlfb_blunt_use_l.2_year, tlfb_blunt_use_l.3_year, 
    tlfb_mj_conc_use_l.1_year, tlfb_mj_conc_use_l.2_year, tlfb_mj_conc_use_l.3_year, 
    tlfb_mj_drink_use_l.1_year, tlfb_mj_drink_use_l.2_year, tlfb_mj_drink_use_l.3_year, 
    tlfb_tincture_use_l.1_year, tlfb_tincture_use_l.2_year, tlfb_tincture_use_l.3_year, 
    tlfb_mj_synth_use_l.1_year, tlfb_mj_synth_use_l.2_year, tlfb_mj_synth_use_l.3_year, 
    tlfb_coc_use_l.1_year, tlfb_coc_use_l.2_year, tlfb_coc_use_l.3_year, 
    tlfb_bsalts_use_l.1_year, tlfb_bsalts_use_l.2_year, tlfb_bsalts_use_l.3_year, 
    tlfb_meth_use_l.1_year, tlfb_meth_use_l.2_year, tlfb_meth_use_l.3_year, 
    tlfb_mdma_use_l.1_year, tlfb_mdma_use_l.2_year, tlfb_mdma_use_l.3_year,  
    tlfb_ket_use_l.1_year, tlfb_ket_use_l.2_year, tlfb_ket_use_l.3_year, 
    tlfb_ghb_use_l.1_year, tlfb_ghb_use_l.2_year, tlfb_ghb_use_l.3_year, 
    tlfb_opi_use_l.1_year, tlfb_opi_use_l.2_year, tlfb_opi_use_l.3_year, 
    tlfb_lsd_use_l.1_year, tlfb_lsd_use_l.2_year, tlfb_lsd_use_l.3_year, 
    tlfb_shrooms_use_l.1_year, tlfb_shrooms_use_l.2_year, tlfb_shrooms_use_l.3_year, 
    tlfb_salvia_use_l.1_year, tlfb_salvia_use_l.2_year, tlfb_salvia_use_l.3_year,
    tlfb_steroids_use_l.1_year, tlfb_steroids_use_l.2_year, tlfb_steroids_use_l.3_year,
    tlfb_inhalant_use_l.1_year, tlfb_inhalant_use_l.2_year, tlfb_inhalant_use_l.3_year,
    tlfb_amp_use_l.1_year, tlfb_amp_use_l.2_year, tlfb_amp_use_l.3_year,
    tlfb_tranq_use_l.1_year, tlfb_tranq_use_l.2_year, tlfb_tranq_use_l.3_year,
    tlfb_vicodin_use_l.1_year, tlfb_vicodin_use_l.2_year, tlfb_vicodin_use_l.3_year,
    tlfb_cough_use_l.1_year, tlfb_cough_use_l.2_year, tlfb_cough_use_l.3_year,
    tlfb_other_use_l.1_year, tlfb_other_use_l.2_year, tlfb_other_use_l.3_year,
    
    # mid-year follow-up (n = 21)
    mypi_alc_sip_1b.6_month, mypi_alc_sip_1b.18_month, mypi_alc_sip_1b.30_month,
    mypi_ecig.6_month, mypi_ecig.18_month, mypi_ecig.30_month, 
    mypi_cigar_used.6_month, mypi_cigar_used.18_month, mypi_cigar_used.30_month,
    mypi_flavoring.6_month, mypi_flavoring.18_month, mypi_flavoring.30_month,
    mypi_chew_pst_used.6_month, mypi_chew_pst_used.18_month, mypi_chew_pst_used.30_month,
    mypi_mj_used.6_month, mypi_mj_used.18_month, mypi_mj_used.30_month,
    mypi_mj_edible.6_month, mypi_mj_edible.18_month, mypi_mj_edible.30_month,
    mypi_mj_oils.6_month, mypi_mj_oils.18_month, mypi_mj_oils.30_month, 
    mypi_mj_tinc_used.6_month, mypi_mj_tinc_used.18_month, mypi_mj_tinc_used.30_month, 
    mypi_mj_vape.6_month, mypi_mj_vape.18_month, mypi_mj_vape.30_month,
    mypi_mj_oils_vaped.6_month, mypi_mj_oils_vaped.18_month, mypi_mj_oils_vaped.30_month,
    mypi_mj_synth_used.6_month, mypi_mj_synth_used.18_month, mypi_mj_synth_used.30_month, 
    mypi_coke_used.6_month, mypi_coke_used.18_month, mypi_coke_used.30_month, 
    mypi_meth_used.6_month, mypi_meth_used.18_month, mypi_meth_used.30_month, 
    mypi_ghb_used.6_month, mypi_ghb_used.18_month, mypi_ghb_used.30_month,
    mypi_heroin_used.6_month, mypi_heroin_used.18_month, mypi_heroin_used.30_month, 
    mypi_sniff_used.6_month, mypi_sniff_used.18_month, mypi_sniff_used.30_month, 
    mypi_pills_used.6_month, mypi_pills_used.18_month, mypi_pills_used.30_month, 
    mypi_pills_dep_used.6_month, mypi_pills_dep_used.18_month, mypi_pills_dep_used.30_month, 
    mypi_pr_used.6_month, mypi_pr_used.18_month, mypi_pr_used.30_month, 
    mypi_cold_used.6_month, mypi_cold_used.18_month, mypi_cold_used.30_month, 
    mypi_high_other_used.6_month, mypi_high_other_used.18_month, mypi_high_other_used.30_month,
    
    DV_BL)

names(data.1)

# check for missing substance use at each timepoint
time_6m <- data.1 %>% 
  select(src_subject_id, ends_with('.6_month')) %>% 
  mutate(tlfb_6m_NA = rowSums(is.na(.))) %>% 
  select(src_subject_id, tlfb_6m_NA)
time_18m <- data.1 %>% 
  select(src_subject_id, ends_with('.18_month')) %>% 
  mutate(tlfb_18m_NA = rowSums(is.na(.))) %>% 
  select(src_subject_id, tlfb_18m_NA)
time_30m <- data.1 %>% 
  select(src_subject_id, ends_with('.30_month')) %>% 
  mutate(tlfb_30m_NA = rowSums(is.na(.))) %>% 
  select(src_subject_id, tlfb_30m_NA)
time_1y <- data.1 %>% 
  select(src_subject_id, ends_with('.1_year')) %>% 
  mutate(tlfb_1y_NA = rowSums(is.na(.))) %>% 
  select(src_subject_id, tlfb_1y_NA)
time_2y <- data.1 %>% 
  select(src_subject_id, ends_with('.2_year')) %>% 
  mutate(tlfb_2y_NA = rowSums(is.na(.))) %>% 
  select(src_subject_id, tlfb_2y_NA)
time_3y <- data.1 %>% 
  select(src_subject_id, ends_with('.3_year')) %>% 
  mutate(tlfb_3y_NA = rowSums(is.na(.))) %>% 
  select(src_subject_id, tlfb_3y_NA)

data.1 <- data.1 %>% 
  left_join(time_6m, by = 'src_subject_id') %>% 
  left_join(time_18m, by = 'src_subject_id') %>% 
  left_join(time_30m, by = 'src_subject_id') %>% 
  left_join(time_1y, by = 'src_subject_id') %>% 
  left_join(time_2y, by = 'src_subject_id') %>% 
  left_join(time_3y, by = 'src_subject_id') 
rm(time_6m, time_18m, time_30m, time_1y, time_2y, time_3y)

# create indicator of missing substance use data across all timepoints    
data.1 <- data.1 %>%   
  mutate(
    tlfb_NA = ifelse ((
      tlfb_6m_NA == 22 & tlfb_18m_NA == 22 & tlfb_30m_NA == 22 & 
        tlfb_1y_NA == 28 & tlfb_2y_NA == 28 & tlfb_3y_NA == 28), 1, 0)) 
table(data.1$tlfb_NA, useNA = 'ifany')
# - no IDs missing all substance use variables across timepoints of interest


###### Create Dependent Variable ###############################################


# calculate substance use at each timepoint 
data.1 <- data.1 %>%  
  mutate(
    .cols = if_any(ends_with('6_month'),
                   .fns = ~ifelse (.x == 1, TRUE, FALSE))) %>% 
  rename(DV_6m = .cols) %>% 
  mutate(
    .cols = if_any(ends_with('18_month'),
                   .fns = ~ifelse (.x == 1, TRUE, FALSE))) %>% 
  rename(DV_18m = .cols) %>% 
  mutate(
    .cols = if_any(ends_with('30_month'),
                   .fns = ~ifelse (.x == 1, TRUE, FALSE))) %>% 
  rename(DV_30m = .cols) %>% 
  mutate(
    .cols = if_any(ends_with('1_year'),
                   .fns = ~ifelse (.x == 1, TRUE, FALSE))) %>% 
  rename(DV_1y = .cols) %>% 
  mutate(
    .cols = if_any(ends_with('2_year'),
                   .fns = ~ifelse (.x == 1, TRUE, FALSE))) %>% 
  rename(DV_2y = .cols) %>% 
  mutate(
    .cols = if_any(ends_with('3_year'),
                   .fns = ~ifelse (.x == 1, TRUE, FALSE))) %>% 
  rename(DV_3y = .cols) %>% 
  replace_na(
    list(DV_6m = FALSE, DV_18m = FALSE, DV_30m = FALSE,
         DV_1y = FALSE, DV_2y = FALSE, DV_3y = FALSE))

# create single dependent variable for any substance use across all timepoints 
data.1 <- data.1 %>% 
  mutate(
    DV = ifelse(DV_6m == TRUE | DV_1y == TRUE | 
                  DV_18m == TRUE | DV_2y == TRUE |
                  DV_30m == TRUE | DV_3y == TRUE, 1, 0))
table(data.1$DV, useNA = 'ifany')

#                 Breakdown of SU at Each Timepoint                            #


# create summary variables for substances w/multiple sub-categories
# - nicotine (annual): 4 variables
# - nicotine (mid-year): 4 variables
# - cannabis (annual): 5 variables
# - cannabis (mid-year): 6 variables
# - Rx (annual): 4 variables
# - Rx (mid-year): 4 variables

data.1 <- data.1 %>% 
  mutate(
    nicotine.1_year = 
      ifelse(
        (tlfb_tob_puff_l.1_year == 1 | tlfb_chew_use_l.1_year == 1 | 
           tlfb_hookah_use_l.1_year == 1 | tlfb_nicotine_use_l.1_year == 1), 
        1, 0),
    nicotine.2_year = 
      ifelse(
        (tlfb_tob_puff_l.2_year == 1 | tlfb_chew_use_l.2_year == 1 | 
           tlfb_hookah_use_l.2_year == 1 | tlfb_nicotine_use_l.2_year == 1), 
        1, 0),
    nicotine.3_year = 
      ifelse(
        (tlfb_tob_puff_l.3_year == 1 | tlfb_chew_use_l.3_year == 1 | 
           tlfb_hookah_use_l.3_year == 1 | tlfb_nicotine_use_l.3_year == 1), 
        1, 0),
    
    cannabis.1_year = 
      ifelse(
        (tlfb_mj_puff_l.1_year == 1 | tlfb_blunt_use_l.1_year == 1 | 
           tlfb_mj_conc_use_l.1_year == 1 | tlfb_mj_drink_use_l.1_year == 1 | 
           tlfb_tincture_use_l.1_year == 1), 
        1, 0),
    cannabis.2_year = 
      ifelse(
        (tlfb_mj_puff_l.2_year == 1 | tlfb_blunt_use_l.2_year == 1 | 
           tlfb_mj_conc_use_l.2_year == 1 | tlfb_mj_drink_use_l.2_year == 1 | 
           tlfb_tincture_use_l.2_year == 1), 
        1, 0),
    cannabis.3_year = 
      ifelse(
        (tlfb_mj_puff_l.3_year == 1 | tlfb_blunt_use_l.3_year == 1 | 
           tlfb_mj_conc_use_l.3_year == 1 | tlfb_mj_drink_use_l.3_year == 1 | 
           tlfb_tincture_use_l.3_year == 1), 
        1, 0),
    
    hall.1_year = 
      ifelse(
        (tlfb_lsd_use_l.1_year == 1 | tlfb_shrooms_use_l.1_year == 1 | 
           tlfb_salvia_use_l.1_year == 1), 
        1, 0),
    
    hall.2_year = 
      ifelse(
        (tlfb_lsd_use_l.2_year == 1 | tlfb_shrooms_use_l.2_year == 1 | 
           tlfb_salvia_use_l.2_year == 1), 
        1, 0),
    
    hall.3_year = 
      ifelse(
        (tlfb_lsd_use_l.3_year == 1 | tlfb_shrooms_use_l.3_year == 1 | 
           tlfb_salvia_use_l.3_year == 1), 
        1, 0),
    
    Rx.1_year = 
      ifelse(
        (tlfb_amp_use_l.1_year == 1 | tlfb_tranq_use_l.1_year == 1 | 
           tlfb_vicodin_use_l.1_year == 1 | tlfb_cough_use_l.1_year == 1), 
        1, 0),
    Rx.2_year = 
      ifelse(
        (tlfb_amp_use_l.2_year == 1 | tlfb_tranq_use_l.2_year == 1 | 
           tlfb_vicodin_use_l.2_year == 1 | tlfb_cough_use_l.2_year == 1), 
        1, 0),
    Rx.3_year = 
      ifelse(
        (tlfb_amp_use_l.3_year == 1 | tlfb_tranq_use_l.3_year == 1 | 
           tlfb_vicodin_use_l.3_year == 1 | tlfb_cough_use_l.3_year == 1), 
        1, 0),
    
    nicotine.6_month = 
      ifelse(
        (mypi_ecig.6_month == 1 | mypi_cigar_used.6_month == 1 | 
           mypi_flavoring.6_month == 1 | mypi_chew_pst_used.6_month == 1), 
        1, 0),
    nicotine.18_month = 
      ifelse(
        (mypi_ecig.18_month == 1 | mypi_cigar_used.18_month == 1 | 
           mypi_flavoring.18_month == 1 | mypi_chew_pst_used.18_month == 1), 
        1, 0),
    nicotine.30_month = 
      ifelse(
        (mypi_ecig.30_month == 1 | mypi_cigar_used.30_month == 1 | 
           mypi_flavoring.30_month == 1 | mypi_chew_pst_used.30_month == 1), 
        1, 0),
    
    cannabis.6_month = 
      ifelse(
        (mypi_mj_used.6_month == 1 | mypi_mj_edible.6_month == 1 | 
           mypi_mj_oils.6_month == 1 | mypi_mj_tinc_used.6_month == 1 | 
           mypi_mj_vape.6_month == 1 | mypi_mj_oils_vaped.6_month == 1), 
        1, 0),
    cannabis.18_month = 
      ifelse(
        (mypi_mj_used.18_month == 1 | mypi_mj_edible.18_month == 1 | 
           mypi_mj_oils.18_month == 1 | mypi_mj_tinc_used.18_month == 1 | 
           mypi_mj_vape.18_month == 1 | mypi_mj_oils_vaped.18_month == 1), 
        1, 0),
    cannabis.30_month = 
      ifelse(
        (mypi_mj_used.30_month == 1 | mypi_mj_edible.30_month == 1 | 
           mypi_mj_oils.30_month == 1 | mypi_mj_tinc_used.30_month == 1 | 
           mypi_mj_vape.30_month == 1 | mypi_mj_oils_vaped.30_month == 1), 
        1, 0),
    
    Rx.6_month = 
      ifelse(
        (mypi_pills_used.6_month == 1 | mypi_pills_dep_used.6_month == 1 | 
           mypi_pr_used.6_month == 1 | mypi_cold_used.6_month == 1), 
        1, 0),
    Rx.18_month = 
      ifelse(
        (mypi_pills_used.18_month == 1 | mypi_pills_dep_used.18_month == 1 | 
           mypi_pr_used.18_month == 1 | mypi_cold_used.18_month == 1), 
        1, 0),
    Rx.30_month = 
      ifelse(
        (mypi_pills_used.30_month == 1 | mypi_pills_dep_used.30_month == 1 | 
           mypi_pr_used.30_month == 1 | mypi_cold_used.30_month == 1), 
        1, 0))



data.1 <- data.1 %>%
  rename(age = tlfb_age_calc_inmonths.base) %>%
  select(src_subject_id, age, starts_with('DV'))

#############fMRI Exclusions----------------------------------------------------
###### Inclusion Criteria #5  ##################################################


# criteria:  
# - exclude if does not meet ABCD quality control inclusion criteria
# - exclude if missing all neuroimaging variables of primary interest

# import resting state f-MRI (n = 272 variables)
net_to_net_fc <- read_csv('data/5.1/mri_y_rsfmr_cor_gp_gp.csv') %>%
  subset(eventname == 'baseline_year_1_arm_1') %>%
  select(src_subject_id, rsfmri_c_ngd_dt_ngd_dt,
         rsfmri_c_ngd_dt_ngd_fo, rsfmri_c_ngd_dt_ngd_sa,
         rsfmri_c_ngd_fo_ngd_fo, rsfmri_c_ngd_fo_ngd_sa, rsfmri_c_ngd_sa_ngd_sa)


# import Quality Control - Recommended Image Inclusion
mri_qc_incl <- read_csv('data/5.1/mri_y_qc_incl.csv') %>% 
  subset(eventname == 'baseline_year_1_arm_1') %>%
  select(src_subject_id, imgincl_rsfmri_include)

mri_qc_fd <- read.csv('data/5.1/mri_y_qc_raw_rsfmr.csv') %>% 
  subset(eventname == 'baseline_year_1_arm_1') %>%
  select(src_subject_id, iqc_rsfmri_all_mean_motion) %>%
  mutate(mean_motion_qc = case_when(
    iqc_rsfmri_all_mean_motion > 0.5 | is.na(iqc_rsfmri_all_mean_motion) ~ 0,
    iqc_rsfmri_all_mean_motion <= 0.5 ~ 1
  )) %>%
  select(src_subject_id, mean_motion_qc, iqc_rsfmri_all_mean_motion)


data.1 <- data.1 %>% 
  left_join(net_to_net_fc, by = 'src_subject_id') %>%
  left_join(mri_qc_incl, by = 'src_subject_id')  %>%
  left_join(mri_qc_fd, by = 'src_subject_id')

data.exclude <- data.exclude %>%
  full_join(data.1 %>% 
              subset(imgincl_rsfmri_include != 1 | mean_motion_qc != 1 |
                       is.na(imgincl_rsfmri_include) | is.na(mean_motion_qc)) %>%
              select(src_subject_id, mean_motion_qc) %>%
              rename(exclude_5 = mean_motion_qc) %>%
              mutate(exclude_5 = TRUE), by = 'src_subject_id')

#Subset data by QC
data.1 <- data.1 %>% 
  subset(imgincl_rsfmri_include == 1 & mean_motion_qc == 1) %>%
  select(-imgincl_rsfmri_include, -age, -mean_motion_qc) %>%
  rename (mean_fd = iqc_rsfmri_all_mean_motion)

#Remove NA rsFC
data.1 <- data.1 %>%
  mutate(exclude_na = ifelse(src_subject_id %in% data.1[!complete.cases(data.1), ]$src_subject_id,
                             TRUE, FALSE))

data.exclude <- data.exclude %>%
  full_join(data.1[!complete.cases(data.1), ] %>% select(src_subject_id, exclude_na),
            by = 'src_subject_id') %>%
  rename(exclude_6 = exclude_na)

data.1 <- na.omit(data.1)

###### expected sample size -----#
# n = 9369
###### expected sample size -----#

rm( net_to_net_fc, net_to_subcort_fc, mri_qc_incl, mri_qc_fd)

data.ex2 <- data.exclude %>%
  pivot_longer(cols = starts_with('exclude'),
               names_to = 'exclusion',
               names_prefix = 'exclude_',
  ) %>%
  filter(value) %>%  # keep only the TRUE one
  mutate(exclusion = readr::parse_number(exclusion)) %>% # get just the number
  select(src_subject_id, exclusion)

########## Initiation over time ------------------------------------------------
data.1 <- data.1 %>%
  mutate(
    DV_1y = ifelse(DV_6m == TRUE, TRUE, DV_1y),
    DV_18m = ifelse(DV_1y == TRUE, TRUE, DV_18m),
    DV_2y = ifelse(DV_18m == TRUE, TRUE, DV_2y),
    DV_30m = ifelse(DV_2y == TRUE, TRUE, DV_30m),
    DV_3y = ifelse(DV_30m == TRUE, TRUE, DV_3y)
  )

data.1 <- data.1 %>%
  select(src_subject_id, DV, starts_with('rsfmri'), mean_fd)

rm(su_gating, su_y_mypi, su_y_sui)

# NeuroCombat for batch effects -----------------------------------------------
# Control for family, site and scanner effects 

## Site and Family effects

site <- read.csv('data/5.1/abcd_y_lt.csv') %>%
  subset(eventname == 'baseline_year_1_arm_1') %>%
  select(src_subject_id, site_id_l, interview_age)

data.ex2 <- data.ex2 %>%
  left_join(site %>% select(src_subject_id, interview_age), by = 'src_subject_id')

## Scanner effects

scanner <- read.csv('data/5.1/mri_y_adm_info.csv') %>%
  subset(eventname == 'baseline_year_1_arm_1') %>%
  select(src_subject_id, mri_info_manufacturer, mri_info_deviceserialnumber)

##Sex effects
sex <- read.csv('data/5.1/abcd_p_demo.csv') %>% 
  subset(eventname == 'baseline_year_1_arm_1') %>%
  select(src_subject_id,demo_sex_v2) %>% 
  mutate(
    demo_sex_v2 = case_when(demo_sex_v2 == 1 | demo_sex_v2 == 3 ~ 1,
                            demo_sex_v2 == 2 | demo_sex_v2 == 4 ~ 2),
    demo_sex_v2 = factor(demo_sex_v2, 
                         levels = c(1,2), 
                         labels = c('Male', 'Female'))) %>%
  rename(sex = demo_sex_v2)

controls <- site %>%
  left_join(scanner, by = 'src_subject_id') %>%
  left_join(sex, by = 'src_subject_id')


## BASELINE RESIDUALS

# 1. Prepare data
data.1 <- na.omit(data.1) %>%
  left_join(controls, by = 'src_subject_id') %>%
  rename(site = site_id_l,
         brand = mri_info_manufacturer,
         device = mri_info_deviceserialnumber) %>%
  mutate(batch = interaction(site, brand, device, drop = TRUE)) %>%
  select(-site, -brand, -device)

# 2. Extract features and set rownames to ID
features <- data.1 %>%
  select(src_subject_id, starts_with('rsfmri'))

feature_matrix <- features %>%
  select(-src_subject_id) %>%
  as.matrix()

rownames(feature_matrix) <- features$src_subject_id  # critical step

# 3. Model matrix
mod <- model.matrix(~1, data = data.1)

# 4. Run ComBat
combat_output <- ComBat(
  dat = t(feature_matrix),
  batch = data.1$batch,
  mod = mod,
  par.prior = TRUE,
  prior.plots = FALSE
)

# 5. Transpose back and keep rownames as IDs
features_resid <- as.data.frame(t(combat_output))
features_resid$src_subject_id <- rownames(features_resid)

# 6. Merge back explicitly by ID
data.2 <- data.1 %>%
  select(-starts_with('rsfmri'), -batch) %>%
  left_join(features_resid, by = "src_subject_id")




# Visualize Variance ----------------------------------------------------------
raw.plot <- ggplot(data.1 %>%
                     select(-batch) %>%
                     pivot_longer(
                       cols = starts_with('rsfmri'),
                       names_to = 'Measure',
                       values_to = 'rsfmri_value'
                     ) %>%
                     mutate(
                       Measure = factor(Measure,
                                        levels = c(
                                          'rsfmri_c_ngd_dt_ngd_dt',
                                          'rsfmri_c_ngd_fo_ngd_fo',
                                          'rsfmri_c_ngd_sa_ngd_sa',
                                          'rsfmri_c_ngd_dt_ngd_sa',
                                          'rsfmri_c_ngd_dt_ngd_fo',
                                          'rsfmri_c_ngd_fo_ngd_sa'
                                        ), labels = c(
                                          'DMN:DMN',
                                          'FPN:FPN',
                                          'SN:SN',
                                          'DMN:SN',
                                          'DMN:FPN',
                                          'FPN:SN'
                                        )),
                       DV = factor(DV,
                                   levels = c(0,1),
                                   labels = c('Uninitiated', 'Initiated'))
                     ),
                   aes(x = DV, y = rsfmri_value, color = factor(DV))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~ Measure, nrow = 1, scales = "fixed") +  # 6 columns, 1 row
  labs(
    title = "Raw rsfMRI Values by Initiation Status Across Measures",
    x = "Initiation Status",
    y = "rsfMRI Value"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
raw.plot


ggsave('output/study3/study 3 baseline raw variance plot.png', raw.plot, 'png',
       path = NULL,
       scale = 1,
       width = 11,
       height = 8.5,
       units = 'in',
       dpi = 320,
       limitsize = FALSE,
       bg = 'white',
       create.dir = FALSE)

data.long <- data.2 %>%
  pivot_longer(
    cols = starts_with('rsfmri'),
    names_to = 'Measure',
    values_to = 'rsfmri_value'
  )

resid.plot <- ggplot(data.long %>%
                       mutate(
                         Measure = factor(Measure,
                                          levels = c(
                                            'rsfmri_c_ngd_dt_ngd_dt',
                                            'rsfmri_c_ngd_fo_ngd_fo',
                                            'rsfmri_c_ngd_sa_ngd_sa',
                                            'rsfmri_c_ngd_dt_ngd_sa',
                                            'rsfmri_c_ngd_dt_ngd_fo',
                                            'rsfmri_c_ngd_fo_ngd_sa'
                                          ), labels = c(
                                            'DMN:DMN',
                                            'FPN:FPN',
                                            'SN:SN',
                                            'DMN:SN',
                                            'DMN:FPN',
                                            'FPN:SN'
                                          )),
                         DV = factor(DV,
                                     levels = c(0,1),
                                     labels = c('Uninitiated', 'Initiated'))
                       ),
                     aes(x = DV, y = rsfmri_value, color = factor(DV))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~ Measure, nrow = 1, scales = "fixed") +  # 6 columns, 1 row
  labs(
    title = "Residual rsfMRI Values by Initiation Status Across Measures",
    x = "Initiation Status",
    y = "rsfMRI Value"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
resid.plot

ggsave('output/study3/study 3 baseline resid variance plot.png', resid.plot, 'png',
       path = NULL,
       scale = 1,
       width = 11,
       height = 8.5,
       units = 'in',
       dpi = 320,
       limitsize = FALSE,
       bg = 'white',
       create.dir = FALSE)

##### Demographics ----------------------------------------------------------


#                 ABCD Parent Demographics Questionnaire


abcd_p_demo <- read.csv('data/5.1/abcd_p_demo.csv', header = TRUE) %>% 
  subset (eventname == 'baseline_year_1_arm_1') %>% 
  select(src_subject_id, demo_prnt_age_v2, demo_sex_v2, demo_gender_id_v2, 
         demo_comb_income_v2,
         demo_race_a_p___10:demo_race_a_p___25, demo_race_a_p___0, 
         demo_race_a_p___77, demo_race_a_p___99, 
         demo_ethn_v2, demo_relig_v2, demo_prnt_ed_v2) %>% 
  mutate(
    demo_comb_income_v2 = case_when(
      demo_comb_income_v2 <= 4 ~ 1,
      demo_comb_income_v2 >= 5 & demo_comb_income_v2 < 7 ~ 2,
      demo_comb_income_v2 >= 7 & demo_comb_income_v2 < 9 ~ 3,
      demo_comb_income_v2 >= 9  & demo_comb_income_v2 < 777 ~ 4,
      demo_comb_income_v2 == 999 | demo_comb_income_v2 == 777 ~ NA
    )
  ) %>% 
  mutate(
    sex = factor(demo_sex_v2, 
                 levels = c(1,2,3,4), 
                 labels = c('Male', 'Female', 
                            'Intersex-Male', 'Intersex-Female')),
    eth_hisp = factor(demo_ethn_v2, 
                      levels = c(1,2), 
                      labels = c('Hispanic', 'non_Hispanic')),
    gender = factor(demo_gender_id_v2, 
                    levels = c(1,2,3,4,5,6), 
                    labels = c('Male', 'Female', 'Trans male', 'Trans female',
                               'Gender queer', 'Different')),
    income = factor(demo_comb_income_v2, 
                    levels = c(1, 2, 3, 4),
                    labels = c('inc_1', 'inc_2', 'inc_3', 'inc_4')),
    religion = factor(demo_relig_v2,
                      levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 
                                 14, 15, 16, 17),
                      labels = c('rp_1', 'rp_2', 'rp_3', 'rp_4', 'rp_5',
                                 'rp_6', 'rp_7', 'rp_8', 'rp_9', 'rp_10',
                                 'rp_11', 'rp_12', 'rp_13', 'rp_14', 
                                 'rp_15', 'rp_16', 'rp_17')))

# recode values of (777) Refuse to Answer and (999) as NA
abcd_p_demo <- abcd_p_demo %>%  
  mutate(
    across(.cols = ends_with('_v2'),  
           .fns = ~replace(.x, .x %in% c(999, 777), NA)))

# create summary education variable 
abcd_p_demo <- abcd_p_demo %>% 
  mutate(
    p_edu = case_when(
      demo_prnt_ed_v2 < 13 ~ 'Less_than_HS_Degree_GED_Equivalent',
      demo_prnt_ed_v2 == 13 ~ 'HS_Graduate_GED_Equivalent',
      demo_prnt_ed_v2 == 14 ~ 'HS_Graduate_GED_Equivalent',
      demo_prnt_ed_v2 == 15 ~ 'Some_College_or_Associates_Degree',
      demo_prnt_ed_v2 == 16 ~ 'Some_College_or_Associates_Degree',
      demo_prnt_ed_v2 == 17 ~ 'Some_College_or_Associates_Degree',
      demo_prnt_ed_v2 == 18 ~ 'Bachelors_Degree',
      demo_prnt_ed_v2 == 19 ~ 'Masters_Degree',
      demo_prnt_ed_v2 >= 20 ~ 'Professional_School_or_Doctoral_Degree')) %>% 
  mutate(p_edu = as.factor(p_edu))
table(abcd_p_demo$p_edu, useNA = 'ifany')

# create summary race variables           
abcd_p_demo <- abcd_p_demo %>%   
  mutate(
    White = ifelse(demo_race_a_p___10 == 1, 1, 0),
    Black = ifelse(demo_race_a_p___11 == 1, 1, 0),
    Asian = ifelse(demo_race_a_p___18 == 1 | demo_race_a_p___19 == 1 | 
                     demo_race_a_p___20 == 1 | demo_race_a_p___21 == 1 | 
                     demo_race_a_p___22 == 1 | demo_race_a_p___23 == 1 |
                     demo_race_a_p___24 == 1, 1, 0),
    AIAN = ifelse( #AIAN: American Indian Alaskan Native
      demo_race_a_p___12 == 1 | demo_race_a_p___13 == 1, 1, 0),
    NHPI = ifelse( #NHPI: Native Hawaiian and Other Pacific
      demo_race_a_p___14 == 1 | demo_race_a_p___15 == 1 | 
        demo_race_a_p___16 == 1 | demo_race_a_p___17, 1, 0),
    Other = ifelse(demo_race_a_p___25 == 1 | demo_race_a_p___0 == 1, 1, 0)) %>% 
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
  select(src_subject_id, 
         sex, gender, race_4l, eth_hisp, income, religion, p_edu) %>% 
  mutate(
    sex_2l = case_when(
      sex == 'Intersex-Male' ~ 'Male',
      sex == 'Male' ~ 'Male', 
      sex == 'Female' ~ 'Female')) %>% 
  select(-sex) %>% 
  mutate(sex_2l = as.factor(sex_2l))
str(abcd_p_demo)

# merge w/full dataset 
data.2 <- data.2 %>% 
  left_join(abcd_p_demo, by = 'src_subject_id')

data.ex2 <- data.ex2 %>%
  left_join(abcd_p_demo, by = 'src_subject_id') %>%
  select(-gender, -religion)

rm(abcd_p_demo)

data.2 <- data.2 %>% 
  select(-gender, -religion)

#Scale residuals and factor DV
data.2 <- data.2 %>%
  rename(
    DMN_DMN = rsfmri_c_ngd_dt_ngd_dt,
    FPN_FPN = rsfmri_c_ngd_fo_ngd_fo,
    SN_SN = rsfmri_c_ngd_sa_ngd_sa,
    DMN_SN = rsfmri_c_ngd_dt_ngd_sa,
    DMN_FPN = rsfmri_c_ngd_dt_ngd_fo,
    FPN_SN = rsfmri_c_ngd_fo_ngd_sa
  ) %>%
  mutate(
    DV = factor(DV,
                levels = c(0,1),
                labels = c('Uninitiated', 'Initiated'))
  )



# Exclusion Analysis -------------------------------------------------------

data.ex2 <- data.ex2 %>%
  mutate(exclusion = TRUE)

data.temp <- data.2 %>%
  mutate(exclusion = FALSE) %>%
  select(src_subject_id, exclusion, 
         race_4l, eth_hisp, income, p_edu, sex_2l, interview_age)

data.ex2 <- rbind(data.ex2, data.temp)

# Step 1: Original categorical summary (unchanged)
vars_cat <- c("race_4l", "eth_hisp", "income", "p_edu", "sex_2l")

summary_tables_cat <- lapply(vars_cat, function(var) {
  data.ex2 %>%
    group_by(exclusion, .data[[var]]) %>%
    summarise(n = n(), .groups = 'drop') %>%
    mutate(variable = var,
           level = as.character(.data[[var]])) %>%
    select(variable, level, exclusion, n)
}) %>%
  bind_rows()

total_by_exclusion <- summary_tables_cat %>%
  group_by(exclusion, variable) %>%
  summarise(total = sum(n), .groups = 'drop')

summary_df_cat <- summary_tables_cat %>%
  left_join(total_by_exclusion, by = c("exclusion", "variable")) %>%
  mutate(percent = round(100 * n / total, 1)) %>%
  select(variable, level, exclusion, n, percent) %>%
  mutate(exclusion = case_when(exclusion == TRUE ~ 1,
                               exclusion == FALSE ~ 0))

summary_df_cat_wide <- summary_df_cat %>%
  pivot_wider(
    names_from = exclusion,
    values_from = c(n, percent),
    names_glue = "exclusion_{exclusion}_{.value}"
  ) %>%
  rename(
    `Freq (Included)` = exclusion_0_n,
    `% within included` = exclusion_0_percent,
    `Freq (Excluded)` = exclusion_1_n,
    `% within excluded` = exclusion_1_percent
  ) %>%
  replace_na(list(
    `Freq (Included)` = 0, `% within included` = 0,
    `Freq (Excluded)` = 0, `% within excluded` = 0
  ))

# Step 2: Run chi-square/Fisher tests for categorical vars
run_test <- function(varname) {
  tab <- table(data.ex2[[varname]], data.ex2$exclusion)
  if (any(tab < 5) && all(dim(tab) == c(2,2))) {
    test <- fisher.test(tab)
    return(paste0("Fisher p = ", signif(test$p.value, 3)))
  } else if (all(dim(tab) == c(2,2))) {
    test <- chisq.test(tab, correct = TRUE)
    return(paste0("Chi-sq (Yates) p = ", signif(test$p.value, 3)))
  } else {
    test <- chisq.test(tab)
    return(paste0("Chi-sq p = ", signif(test$p.value, 3)))
  }
}

test_results_cat <- tibble(
  variable = vars_cat,
  `Freq. Test (p-value)` = map_chr(vars_cat, run_test)
)

cat_table <- summary_df_cat_wide %>%
  left_join(test_results_cat, by = "variable")

rm(data.temp)
# Step 3: Summary for continuous var: interview_age_1
data.ex2 <- data.ex2 %>%
  mutate(exclusion = case_when(exclusion == TRUE ~ 1,
                               exclusion == FALSE ~ 0),
         exclusion = factor(exclusion,
                            levels = c(0,1)))

age_summary <- read.csv('data/5.1/abcd_y_lt.csv') %>%
  filter(eventname == 'baseline_year_1_arm_1') %>%
  select(src_subject_id, interview_age) %>%
  right_join(data.ex2 %>%
               select(src_subject_id, exclusion), by = 'src_subject_id') %>%
  group_by(exclusion) %>%
  summarise(
    mean_age = mean(interview_age, na.rm = TRUE),
    sd_age = sd(interview_age, na.rm = TRUE),
    min_age = min(interview_age, na.rm = TRUE),
    max_age = max(interview_age, na.rm = TRUE),
  ) %>%
  pivot_wider(names_from = exclusion, values_from = c(mean_age, sd_age, min_age, max_age), names_prefix = "exclusion_") %>%
  mutate(
    variable = "interview_age",
    `Freq (Included)` = sprintf("%.1f (%.1f)", mean_age_exclusion_0, sd_age_exclusion_0),
    `Freq (Excluded)` = sprintf("%.1f (%.1f)", mean_age_exclusion_1, sd_age_exclusion_1),
    `% within included` = sprintf("Min: %.1f\nMax: %.1f", min_age_exclusion_0, max_age_exclusion_0),
    `% within excluded` = sprintf("Min: %.1f\nMax: %.1f", min_age_exclusion_1, max_age_exclusion_1),
    `Freq. Test (p-value)` = sprintf("t-test: p=%.4f", 
                                     t.test(data.ex2$interview_age ~ data.ex2$exclusion)$p.value)
  ) %>%
  select(variable, `Freq (Included)`, `Freq (Excluded)`, `% within included`, `% within excluded`, `Freq. Test (p-value)`)


# Coerce relevant columns to character in both tables before binding
cat_table_fixed <- cat_table %>%
  mutate(across(c(`Freq (Included)`, `Freq (Excluded)`,
                  `% within included`, `% within excluded`,
                  `Freq. Test (p-value)`), as.character))

age_summary_fixed <- age_summary %>%
  mutate(across(c(`Freq (Included)`, `Freq (Excluded)`,
                  `% within included`, `% within excluded`,
                  `Freq. Test (p-value)`), as.character))

# Now bind rows safely
final_table <- bind_rows(cat_table_fixed, age_summary_fixed)

data.ex2 <- data.ex2 %>%
  mutate(exclusion = factor(exclusion,
                            levels = c(0,1),
                            labels = c('Included', 'Excluded')))

write.csv(final_table,'output/study3/Study 3 Exclusion Analysis.csv')

rm(scanner, sex, site, su_base, summary_df_cat, summary_df_cat_wide, summary_tables_cat,
   test_results_cat, total_by_exclusion, data.exclude, feature_matrix, combat_output, cat_table,
   cat_table_fixed, age_summary, age_summary_fixed, vars_cat, run_test, features, features_resid,
   controls, mod)

######## Mean Differences ------------------------------------------------------

data.long <- data.2 %>%
  pivot_longer(
    cols = c("DMN_DMN", "FPN_FPN", "SN_SN",
             "DMN_FPN", "DMN_SN", "FPN_SN"),
    names_to = 'Measure',
    values_to = 'rsfmri_value'
  )

run_test <- function(data, covariates, dv = "DV") {
  bind_rows(
    lapply(covariates, function(var) {
      x <- data[[var]]
      
      # Categorical
      if (is.factor(x) || is.character(x)) {
        tab <- table(x, data[[dv]])
        prop_initiated <- prop.table(tab, margin = 2)[, "Initiated"]
        prop_uninitiated <- prop.table(tab, margin = 2)[, "Uninitiated"]
        total_initiated <- tab[, "Initiated"]
        total_uninitiated <- tab[, "Uninitiated"]
        
        tab_total <- table(x)
        prop_total <- prop.table(tab_total)
        total_all <- as.numeric(tab_total)
        
        if (all(dim(tab) == c(2, 2))) {
          test <- fisher.test(tab)
          test_type <- "Fisher"
          statistic <- NA_real_
          df <- NA_real_
          p.value <- as.numeric(test$p.value)
          residuals_str <- NA_character_  # Fisher's test doesn't have residuals
        } else {
          test <- suppressWarnings(chisq.test(tab))
          test_type <- "Chi-sq"
          statistic <- as.numeric(test$statistic)
          df <- as.numeric(test$parameter)
          p.value <- as.numeric(test$p.value)
          
          # Extract standardized residuals
          std_residuals <- test$stdres
          # Format residuals as a string (e.g., "group1: 2.3, group2: -1.5")
          residuals_str <- paste(
            paste0(rownames(std_residuals), ": ", 
                   round(std_residuals[, "Initiated"], 2)),
            collapse = "; "
          )
        }
        
        tibble(
          GroupingVariable = var,
          group = names(prop_initiated),
          prop_initiated = as.character(round(prop_initiated, 3)),
          prop_uninitiated = as.character(round(prop_uninitiated, 3)),
          total_initiated = as.character(total_initiated),
          total_uninitiated = as.character(total_uninitiated),
          prop_total = as.character(round(prop_total, 3)),
          total_all = as.character(total_all),
          test_type = test_type,
          statistic = statistic,
          df = df,
          p.value = p.value,
          std_residuals = residuals_str
        )
        
        # Continuous
      } else if (is.numeric(x)) {
        init_vals <- data[data[[dv]] == "Initiated", var, drop = TRUE]
        uninit_vals <- data[data[[dv]] == "Uninitiated", var, drop = TRUE]
        all_vals <- data[[var]]
        
        test <- t.test(init_vals, uninit_vals)
        
        tibble(
          GroupingVariable = var,
          group = NA_character_,
          prop_initiated = paste(range(init_vals, na.rm = TRUE), collapse = ","),
          prop_uninitiated = paste(range(uninit_vals, na.rm = TRUE), collapse = ","),
          total_initiated = paste0(as.character(round(mean(init_vals, na.rm = TRUE), 3)),' (',as.character(round(sd(init_vals, na.rm = TRUE), 2)), ')'),
          total_uninitiated = paste0(as.character(round(mean(uninit_vals, na.rm = TRUE), 3)),' (',as.character(round(sd(uninit_vals, na.rm = TRUE), 2)), ')'),
          prop_total = paste(range(all_vals, na.rm = TRUE), collapse = ","),
          total_all = paste0(as.character(round(mean(all_vals, na.rm = TRUE), 3)),
                             ' (', as.character(round(sd(all_vals, na.rm = TRUE), 2)), ')'),
          test_type = "t-test",
          statistic = as.numeric(test$statistic),
          df = as.numeric(test$parameter),
          p.value = test$p.value,
          std_residuals = NA_character_  # t-tests don't have residuals
        )
      } else {
        NULL
      }
    })
  )
}

group_vars <- c("race_4l", "eth_hisp", "sex_2l", "income", "p_edu", 'interview_age',
                'DMN_DMN', 'FPN_FPN', 'SN_SN', 'DMN_FPN', 'DMN_SN', 'FPN_SN')
mean_results <- bind_rows(lapply(group_vars, function(v) run_test(data.2, v)))

write.csv(mean_results,'output/study3/Study 3 Baseline Mean Differences V3.csv')

data.2 <- data.2 %>%
  mutate(interview_age = as.numeric(scale(interview_age)),
         mean_fd = as.numeric(scale(mean_fd)))
