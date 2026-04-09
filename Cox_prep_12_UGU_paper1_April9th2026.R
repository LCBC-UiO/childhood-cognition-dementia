#Revised script for dementia prediction paper
#by KBW April9th 2026
library (tidyverse)
library (tidyverse)
library(ggplot2)
library(survival)
library(ragg)
library(haven)
library(tidyr)

#first run scripts get NPR_inpatient4 and getCDR..5
#as of Jan24th, run NPR_inpatient4 and getCDR..5 (fixing all DODSDAT)
#the scripts creates dementia diagnoses
source("/safe/data/KBW/get_NPR_inpatient_dementia_diagnoses4.R")
source("/safe/data/KBW/get_CDR_dementia_diagnoses5_jan24th.R")

library (tidyverse)
#my signs: §|~ %>% $ {} [ ]
#NPR_dem_only already contains only
  #select(ID, NPR_dementia, NPR_dem_date)

class(NPR_dem_only$NPR_dem_date)
head(NPR_dem_only$NPR_dem_date)

NPR_dem_only <-NPR_dem_only %>%
  mutate(NPR_dem_date = as.Date(NPR_dem_date))

CDR_dem <- CDR%>%
  select(ID, CDR_dementia, CDR_dem_date, CDR_death_date)

#take earliest dementia_date across sources
dem_all <-NPR_dem_only%>%
  full_join(CDR_dem, by = "ID")%>%
             mutate(
               NPR_dementia = coalesce (NPR_dementia, 0L),
               #had mistake here CDR_dmenetia
               CDR_dementia =coalesce(CDR_dementia, 0L),
               dementia_any =if_else(
                 NPR_dementia ==1L | CDR_dementia == 1L,
                 1L, 0L
               ),
               dementia_date = case_when(
                 NPR_dementia ==1L & CDR_dementia == 1L ~ pmin (NPR_dem_date, CDR_dem_date),
                 NPR_dementia ==1L ~ NPR_dem_date,
                 CDR_dementia ==1L ~ CDR_dem_date,
                 TRUE ~ as.Date (NA_character_)
               )
             )
               
table (dem_all$dementia_any)

end_of_follow_up <-as.Date("2025-11-03") #the date Andreas got the registry data

dem_all <-dem_all%>%
  mutate(
    censor_date = case_when(
      !is.na(CDR_death_date) ~ CDR_death_date,
      TRUE ~ end_of_follow_up))
      
#check how many unique cases, n =415
  
all_dementia <- dem_all%>%
  select ("ID", "dementia_any")%>%
    filter(dementia_any ==1)%>%
      distinct(ID, .keep_all = TRUE) 

#read in ugu1948
library (haven)
ugu1948 <-read_sav ("/safe/data/UGU-raw-SPSS/ugu1948_English_v2.4_2024_11_20_LIFE.sav")
#View(ugu1948)
#11945 obs of 386 vars
ugu1948 <-ugu1948%>%
  mutate (ID = as.integer (ugukod) )

############take away registry opt-out IDs########## 
#as informed of March 12th 2026
#IDs with opt-out-of-registry

#the script registry_opt_out.R creates optout_ids
source("/safe/data/KBW/registry_opt_out.R")

#exclude those who have opted out from ugu1948
ugu1948 <-ugu1948%>%
  filter(!ID %in% optout_ids)


ugu1948 <-ugu1948%>%
  filter(!ID %in% optout_ids)

#full data
cox_ugu1948 <- ugu1948%>% 
left_join(dem_all, by = "ID")

#replace missing dementia_any with 0
cox_ugu1948 <- cox_ugu1948%>%
  mutate(
    dementia_any = if_else (is.na(dementia_any), 0L, dementia_any)
    
  )
#ensure everyone has a censor_date (death or end_of_follow_up, whichever comes first)
cox_ugu1948 <-cox_ugu1948%>%
  mutate(
    censor_date = case_when(
      !is.na(CDR_death_date) ~ CDR_death_date,
      TRUE ~ end_of_follow_up
    )
  )

##create birth_date using RSYEAR and RSMONTH
#use the 15th of the months as standard mind-month approximation
cox_ugu1948 <- cox_ugu1948%>%
  mutate(
    birth_month = sprintf("%02d",RSMONTH),
    birth_date = as.Date(paste0(RSYEAR, "-", birth_month, "-15"))
  )

#define start date (all entered in 1961)
#it says tests took place in the period May8th-27th 1961
cox_ugu1948 <- cox_ugu1948%>%
  mutate(start_date = as.Date("1961-05-08"))
#compute start_age
cox_ugu1948 <-cox_ugu1948%>%
  mutate(
    start_age = as.numeric((start_date-birth_date)/365.25)
  )

#######cause of reused personnummere#####
####check that dementia_date is not after CDR_death_date####
#It says FALSE, so we+re clean
any(cox_ugu1948$dementia_date > cox_ugu1948$CDR_death_date, na.rm = TRUE)



#compute end date and end_age(dementia or censoring)
cox_ugu1948 <-cox_ugu1948%>%
  mutate(
    event_date = if_else(dementia_any ==1L, dementia_date, censor_date),
    end_age = as.numeric((event_date -birth_date)/365.25)
  )
##########do some stuff on ugu1948 again to prep########
# I will scale the tests again
cox_ugu1948 <- cox_ugu1948%>%
  mutate(
    z_TS6IITP = as.numeric(scale (TS6IITP)),
    z_TS6ISTP = as.numeric(scale (TS6ISTP)),
    z_TS6IVOTP = as.numeric(scale (TS6IVOTP)),
  )
#and I make sex into 0 and 1, keeping NA where necessary, so here, 
#if my raw var RSSEX uses as in codebook 1 = male, 2 = female,
#then I recode, to keep direction consistent, 0 = male (reference), 1 = female
cox_ugu1948 <-cox_ugu1948%>%
  mutate(
    sex01 =ifelse(RSSEX ==2,1,
                  ifelse(RSSEX ==1,0, NA))
  )
#recode 9s to NA in mother father edu
cox_ugu1948 <- cox_ugu1948 %>%
  mutate(
    RFSUN4 = na_if(RFSUN4, 9),
    RMSUN4   = na_if(RMSUN4, 9)
  )
#make mean parental edu, and use either val if one missing
cox_ugu1948 <- cox_ugu1948 %>%
  mutate(
    parentedu_mean = rowMeans(select(., RFSUN4, RMSUN4), na.rm = TRUE)
  )
#standardize parental edu var
cox_ugu1948 <- cox_ugu1948%>%
  mutate(
    z_parentedu_mean = as.numeric(scale (parentedu_mean)),
  )

#some might argue we should correct for 3 tests, 
#but they are so related
#so, we do another analysis to alleviate this concern
#Now, instead of testing test1: beta1 =0, test 2; beta2=, test 3: beta3 =0,
#which gives a multiple testing problem, we test one gathered hypothesis
#H0=betaIITP =BetaISTP =BetaIVOTP =0
#against
#H1: at least one of them is not null
#this is an omnibus-likelihood ratio-test 
#df =3
#controlling for all other covariates
#avoids multiple testing (it is one hypothesis, not three)
#gives evidence for cognitive ability as a block
#is statistically decent, when subtests are part of one same model

# I filter out make sure we have same dataset on all vars
model_data48 <- cox_ugu1948 %>%
  select(ID, dementia_any, sex01, z_parentedu_mean,
         z_TS6IITP, z_TS6ISTP, z_TS6IVOTP,
         start_age, end_age,) %>%
  filter(!is.na(start_age)) %>%
  filter(!is.na(end_age)) %>%
  filter(!is.na(dementia_any)) %>%
  filter(!is.na(sex01)) %>%
  filter(!is.na(z_parentedu_mean)) %>%
  filter(!is.na(z_TS6IITP))%>%
  filter(!is.na(z_TS6ISTP))%>%
  filter(!is.na(z_TS6IVOTP))

#my signs: §|~ %>% $ {} [ ]


table (model_data48$sex01)


cox_full <-coxph(formula =
                   Surv(start_age, end_age, dementia_any) ~
                   z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                   z_parentedu_mean +sex01,
                 data = model_data48
)
summary(cox_full)

#now I am making a reduced model w/o cognition:

cox_reduced <-coxph(formula =
                      Surv(start_age, end_age, dementia_any) ~
                      z_parentedu_mean +sex01,
                    data = model_data48
)
summary(cox_reduced)

#then
anova(cox_reduced, cox_full, test = "LRT")


########now do the cox on ugu1948######

library(survival)

#with opt-out excluded, n= 10539, number of events= 287; was N = 10547 w/opt-outs

cox_all <-coxph(formula =
  Surv(start_age, end_age, dementia_any) ~
    z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
    z_parentedu_mean +sex01,
  data = cox_ugu1948
)

summary(cox_all)



########Make HR plot#######

#now make a forest plot of teh standardized effects, figure 1:
#extract ORs and CIs
HRs <-exp(coef(cox_all))
CIs <-exp(confint.default(cox_all))

#make clean dataframe for plotting
#remember I recoded RSSEX, whic was  1= Male, and 2 = female, to 0 = male and 1 is female
plot_df <-data.frame(
  Predictor = c("Inductive reasoning", "Spatial ability", "Verbal ability",
                "Parental education", "Sex (male vs. female)"),
  HR =HRs,
  CI_low = CIs [,1],
  CI_high = CIs[,2]
)

library (ggplot2)


#inductive on top
plot_df$Predictor <- factor(plot_df$Predictor,
                            levels = c("Sex (male vs. female)",
                                       "Parental education",
                                       "Verbal ability",
                                       "Spatial ability", 
                                       "Inductive reasoning"
                            ))

#I think I want one color for the predictors and one for the covariates, to cheer up the plot also
plot_df$Type <- factor(
  c("Cognitive", "Cognitive", "Cognitive",
    "Covariate", "Covariate"),
  levels = c("Cognitive","Covariate"))


#specify where to put plot
outdir <- "/safe/data/KBW/"

library(ragg)

agg_tiff(
  file.path(outdir, "Figure1_HR_plot_600_dpi.tif"),
  width = 7,
  height =5,
  units = "in",
  res = 600,
  compression ="lzw"
)

print(
  ggplot(plot_df, aes(x= Predictor, y = HR, color = Type))+
  geom_point(size =3)+
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.1, linewidth =0.8) +
  geom_hline(yintercept = 1, linetype = "dashed")+
  coord_flip()+
  theme_classic(base_size =18)+
  scale_color_manual(values = c("Covariate" = "#d62728", #red
                                "Cognitive" = "#1f77b4")) + #blue
  labs(title = "Childhood cognitive scores 
       and dementia risk",
       y = "Hazard ratio (per 1SD)",
       x = "",
       color = ""
  ))

dev.off()

#add sex interaction
#this model is what we do in the paper
cox_all_sex <-coxph(formula =
                      Surv(start_age, end_age, dementia_any) ~
                      z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                      z_parentedu_mean +sex01 +
                      z_TS6IITP:sex01,
                    data = cox_ugu1948
)

summary(cox_all_sex)


#########Jan 6th 2026 - try to add their own education######

#Add their own edu in 1990
FOB_1990_data <-read.csv ("/safe/data/RAW_DATA/SCB/FOB/FoB_1990.csv")

#data from 1990 = edu_42 for 1948 and edu_37 for 1953 
select_FOB_1990 <- FOB_1990_data%>%
  mutate (edu_FOB_1990 = UtbNiva)%>%
  select (ID, edu_FOB_1990)

#can check SUN2, SUN3, SUN5, seems also edu-vars (look in codebook)
#using UtbNiva
#distribitions for utbildniva are 1-7 and then 9, check, codebook says:
#"UtbNiva": {
#  "Question_label": "UtbildningsnivÃ¥",
#  "Item_label": "",
#  "Response_alternatives": 
#"0 = Alla fÃ¶dda efter 1974 och fÃ¶re 1926\
#n1 = FÃ¶rgymnasial utbildning kortare Ã¤n 9 Ã¥r\
#n2 = FÃ¶rgymnasial utbildning 9 Ã¥r\
#n3 = Gymnasial utbildning hÃ¶gst 2 Ã¥r\
#n4 = Gymnasial utbildning lÃ¤ngre Ã¤n 2 Ã¥r men max 3 Ã¥r\
#n5 = Eftergymnasial utbildning kortare Ã¤n 3 Ã¥r (inkl. 4-Ã¥rigt gymn.)
#n6 = Eftergymnasial utbildning 3 Ã¥r  eller lÃ¤ngre (exkl. forskarutb)
#n7 = Forskarutbildning\n9 = ospecificerad nivÃ¥",
#> table (cox_ugu1948$edu_FOB_1990)
#1    2    3    4    5    6    7    9 
#1921 1979 2426 1404 1385 1692   76  509

#recode 9 to NA in FOB_1990
select_FOB_1990$edu_FOB_1990[select_FOB_1990$edu_FOB_1990 ==9] <- NA


########precautions for reused personnummer######

#to take out the one (Jan28th) that has swopped personnumer, first in excel file in vault
#and any other such
#I will simply check that none have dem_date of FOB1990 vals after death_date

cox_ugu1948 <- cox_ugu1948%>% 
  left_join(select_FOB_1990, by = "ID")

table (cox_ugu1948$edu_FOB_1990)

#the following lists all that have edu 1990 but died before in DORS/CDR,
#and it is just this one person
cox_ugu1948[
  !is.na(cox_ugu1948$edu_FOB_1990) &
  !is.na(cox_ugu1948$CDR_death_date) & 
  cox_ugu1948$CDR_death_date < as.Date("1990-01-01"),
  c("ID", "CDR_death_date", "edu_FOB_1990")
]

#now I will replace such impossible cases with NA in principle (here it is just one) 
#post this code in channel, as others may want to use same type of standard in UGU 
cox_ugu1948$edu_FOB_1990[
  !is.na(cox_ugu1948$edu_FOB_1990) &
  !is.na(cox_ugu1948$CDR_death_date) & 
    cox_ugu1948$CDR_death_date < as.Date("1990-01-01")
] <-NA

#and sanity check for no dementia date after death date in CDR, luckily not the case
#but remember to check same for uGU1953 later
any(cox_ugu1948$dementia_date > cox_ugu1948$CDR_death_date, na.rm = TRUE)

#standardize own edu
cox_ugu1948 <- cox_ugu1948%>%
  mutate(
    z_edu_FOB_1990 = as.numeric(scale (edu_FOB_1990)),
  )

########get descriptoves for those having own edu######

# I filter out make sure we have same dataset on all vars
model_data48_own_edu <- cox_ugu1948 %>%
  select(ID, dementia_any, sex01, z_parentedu_mean, 
         z_edu_FOB_1990,edu_FOB_1990,
         z_TS6IITP, z_TS6ISTP, z_TS6IVOTP,
         start_age, end_age,) %>%
  filter(!is.na(start_age)) %>%
  filter(!is.na(end_age)) %>%
  filter(!is.na(dementia_any)) %>%
  filter(!is.na(sex01)) %>%
  filter(!is.na(z_parentedu_mean)) %>%
  filter(!is.na(z_edu_FOB_1990)) %>%
  filter(!is.na(z_TS6IITP))%>%
  filter(!is.na(z_TS6ISTP))%>%
  filter(!is.na(z_TS6IVOTP))

table (model_data48_own_edu$edu_FOB_1990)

#Table of midlife education proportions by sex
prop.table(table(model_data48_own_edu$edu_FOB_1990, model_data48_own_edu$sex01), 2)

#add own edu in cox
cox_all_edu <-coxph(formula =
                      Surv(start_age, end_age, dementia_any) ~
                      z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                      z_parentedu_mean + z_edu_FOB_1990 +sex01,
                    data = model_data48_own_edu
)

summary(cox_all_edu)




#######Make higher res figure 2 with own edu Jan 6th for paper########
#make a forest plot of the standardized effects:
#extract ORs and CIs
HRs <-exp(coef(cox_all_edu))
CIs <-exp(confint.default(cox_all_edu))
#make clean dataframe for plotting
#recoded RSSEX, which was  1= Male, and 2 = female, to 0 = male and 1 is female
plot_df <-data.frame(
  Predictor = c("Inductive reasoning", "Spatial ability", "Verbal ability",
                "Parental education", "Midlife education", "Sex (male vs. female)"),
  HR =HRs,
  CI_low = CIs [,1],
  CI_high = CIs[,2]
)

library (ggplot2)


#inductive top
plot_df$Predictor <- factor(plot_df$Predictor,
                            levels = c("Sex (male vs. female)",
                                       "Midlife education",
                                       "Parental education",
                                       "Verbal ability",
                                       "Spatial ability", 
                                       "Inductive reasoning"
                            ))

#I think I want one color for the predictors and one for the covariates, to cheer up the plot also
plot_df$Type <- factor(
  c("Cognitive", "Cognitive", "Cognitive",
    "Covariate", "Covariate", "Covariate"),
  levels = c("Cognitive","Covariate"))
#specify where to put plot
outdir <- "/safe/data/KBW/"

library(ragg)

agg_tiff(
  file.path(outdir, "Figure2_w_own_1990_edu_HR_plot_600_dpi.tif"),
  width = 7,
  height =5,
  units = "in",
  res = 600,
  compression ="lzw"
)

print(
  ggplot(plot_df, aes(x= Predictor, y = HR, color = Type))+
    geom_point(size =3)+
    geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.1, linewidth =0.8) +
    geom_hline(yintercept = 1, linetype = "dashed")+
    coord_flip()+
    theme_classic(base_size =18)+
    scale_color_manual(values = c("Covariate" = "#d62728", #red
                                  "Cognitive" = "#1f77b4")) + #blue
    labs(title = "Childhood cognitive scores, 
education and dementia risk",
         y = "Hazard ratio (per 1SD)",
         x = "",
         color = ""
    ))

dev.off()



cox_all_sex <-coxph(formula =
                      Surv(start_age, end_age, dementia_any) ~
                      z_TS6IITP + sex01 +
                      z_TS6IITP:sex01,
                    data = cox_ugu1948
)
summary(cox_all_sex)
#makes no sense to do it just on the interaction term
cox_sex_int <-coxph(formula =
                      Surv(start_age, end_age, dementia_any) ~
                      z_TS6IITP:sex01,
                    data = cox_ugu1948
)
summary(cox_sex_int)

#######Now load ugu1953 and merge ######

#LOAD UGU1953
ugu1953 <-read_sav ("/safe/data/UGU-raw-SPSS/ugu1953_English_v2.3_2024-11-20.sav")
#View(ugu1953)
#9929 obs of 270 vars
ugu1953 <-ugu1953%>%
  mutate (ID = as.integer (ugukod) )

#exclude those who have opted out from ugu1948
ugu1953 <-ugu1953%>%
  filter(!ID %in% optout_ids)


#and I make sex into 0 and 1, keeping NA where necessary
#here, rssex ==1 remains 1, where as 2 becomes 0
#ugu1953 <-ugu1953%>%
#  mutate(
#    sex01 =ifelse(RSSEX ==1,1,
#                  ifelse(RSSEX ==2,0, NA))
#  )

#and I make sex into 0 and 1, keeping NA where necessary, so here, 
#if my raw var RSSEX uses as in codebook 1 = male, 2 = female,
#then I recode, to keep direction consistent, 0 = male (reference), 1 = female
ugu1953 <-ugu1953%>%
  mutate(
    sex01 =ifelse(RSSEX ==2,1,
                  ifelse(RSSEX ==1,0, NA))
  )

#full data
cox_ugu1953 <- ugu1953%>% 
  left_join(dem_all, by = "ID")

#replace missing dementia_any with 0
cox_ugu1953 <- cox_ugu1953%>%
  mutate(
    dementia_any = if_else (is.na(dementia_any), 0L, dementia_any)
    
  )
#ensure everyone has a censor_date (death or end_of_follow_up, whichever comes first)
cox_ugu1953 <-cox_ugu1953%>%
  mutate(
    censor_date = case_when(
      !is.na(CDR_death_date) ~ CDR_death_date,
      TRUE ~ end_of_follow_up
    )
  )

##create birth_date using RSYEAR and RSMONTH
#use the 15th of the months as standard mind-month approximation
cox_ugu1953 <- cox_ugu1953%>%
  mutate(
    birth_month = sprintf("%02d",RSMONTH),
    birth_date = as.Date(paste0(RSYEAR, "-", birth_month, "-15"))
  )

#define start date (all entered in 1966)
#it says tests took place in the period May8th-27th 1966????
cox_ugu1953 <- cox_ugu1953%>%
  mutate(start_date = as.Date("1966-05-09"))
#compute start_age
cox_ugu1953 <-cox_ugu1953%>%
  mutate(
    start_age = as.numeric((start_date-birth_date)/365.25)
  )

#compute end date and end_age(dementia or censoring)
cox_ugu1953 <-cox_ugu1953%>%
  mutate(
    event_date = if_else(dementia_any ==1L, dementia_date, censor_date),
    end_age = as.numeric((event_date -birth_date)/365.25)
  )

# I will scale the tests within cohort
cox_ugu1953 <- cox_ugu1953%>%
  mutate(
    z_TS6IITP = as.numeric(scale (TS6IITP)),
    z_TS6ISTP = as.numeric(scale (TS6ISTP)),
    z_TS6IVOTP = as.numeric(scale (TS6IVOTP)),
  )

#######cause of reused personnummere#####
####check that dementia_date is not after CDR_death_date####
#It says FALSE, so we+re clean
any(cox_ugu1953$dementia_date > cox_ugu1953$CDR_death_date, na.rm = TRUE)

#add missing columns from 1948 into 1953 as NA, without adding new data
#taht did not work, damn!
#vars_1948 <-names (ugu1948)
#ugu1953 <-ugu1953%>%
#  mutate (across(setdiff(vars_1948, names(.)),~ NA ))

#Just take the overlapping columns then
common_cols <-intersect(names(cox_ugu1948), names(cox_ugu1953))
ugu1948_small <-cox_ugu1948 %>% select(all_of(common_cols))
ugu1953_small <-cox_ugu1953 %>%select(all_of(common_cols))

ugu_merged <-bind_rows(ugu1948_small, ugu1953_small)
#note, there are a series of conflicting labels when merging, use below to see
warnings() #seems OK to ignore for now, but read up on this!

#maybe best not to restandardize TS6IITP across cohorts
#Cox_ugu4853z <- ugu_merged%>%
#  mutate(
#    z_TS6IITP = as.numeric(scale (TS6IITP)),
#    z_TS6ISTP = as.numeric(scale (TS6ISTP)),
#    z_TS6IVOTP = as.numeric(scale (TS6IVOTP)),
#  )

Cox_ugu4853x <- ugu_merged


cox_all_merged <-coxph(formula =
                  Surv(start_age, end_age, dementia_any) ~
                  z_TS6IITP + z_TS6ISTP + z_TS6IVOTP,
                data = Cox_ugu4853x
)

summary(cox_all_merged)

cox_all_merged <-coxph(formula =
                         Surv(start_age, end_age, dementia_any) ~
                         z_TS6IITP + z_TS6ISTP + z_TS6IVOTP + sex01,
                       data = Cox_ugu4853x
)

summary(cox_all_merged)

cox_all_merged <-coxph(formula =
                         Surv(start_age, end_age, dementia_any) ~
                         z_TS6IITP + sex01,
                       data = Cox_ugu4853x
)

summary(cox_all_merged)




##########test if non-linear in merged cohort#####
#there is no significant smooth effect of inductive reasoning on dementia risk
model_data4853 <- Cox_ugu4853x %>%
  select(ID, dementia_any, TS6IITP, sex01) %>%
  filter(!is.na(dementia_any)) %>%
  filter(!is.na(sex01)) %>%
  filter(!is.na(TS6IITP))

# now n= 19919, number of events= 367 (n was 19942)
#is significant

model_data4853z <- Cox_ugu4853x %>%
  select(ID, dementia_any, z_TS6IITP, sex01) %>%
  filter(!is.na(dementia_any)) %>%
  filter(!is.na(sex01)) %>%
  filter(!is.na(z_TS6IITP))

library (mgcv)

#is there a non-linear effect (smooth) across both cohorts?
#note taht here we do no use standardized scores
model_gam_large <-gam(dementia_any ~s(TS6IITP),
                      data = model_data4853, 
                      family = binomial)
summary(model_gam_large)

#if we do use standardized scores, the smooth term is not significant
model_gam_large_z <-gam(dementia_any ~s(z_TS6IITP),
                      data = model_data4853z, 
                      family = binomial)
summary(model_gam_large_z)
#simple fig, looks like hell, but clear, quite linear this
#surprises me actually, may be Swedes are more linear than Norwegians, I wonder why...
plot (model_gam_large, shade = TRUE, rug = TRUE)
#this is obvious, but I have to do a formal test for the ms, I think
#compare to logistic linear model
model_lin_large <-glm(dementia_any ~ TS6IITP,
                      data = model_data4853, 
                      family = binomial)
anova(model_lin_large, model_gam_large, test = "Chisq")

#make nicer plot for paper
tiff("/safe/data/KBW/Log_odds.tiff", width=6, height=5, units = "in", res = 600)
par(mar = c(5,5,4,2), cex =1.4, cex.lab = 1.6, cex.axis = 1.4)

plot (model_gam_large, 
      shade = TRUE,
      shade.col = "gray90",
      rug = FALSE,
      seWithMean= TRUE,
      xlab = "Childhood inductive reasoning",
      ylab = "Smooth effect on log-odds of dementia",
      main = "",
      bty = "n",
      ylim = c(-0.7, 0.7),
      scheme =1)

dev.off ()

getwd()

#####Make Figure 3 600dpi#####
outdir <- "/safe/data/KBW/"

library (ragg)
agg_tiff(
  file.path(outdir, "Figure3_smooth_600_dpi.tif"),
  width = 7,
  height =5,
  units = "in",
  res = 600,
  compression ="lzw"
)

par(
  mar = c(5,8,4,2), 
  cex =1.4, 
  cex.lab = 1.6, 
  cex.axis = 1.4)

plot (model_gam_large, 
      shade = TRUE,
      shade.col = "gray90",
      rug = FALSE,
      seWithMean= TRUE,
      xlab = "Childhood inductive reasoning",
      ylab = "Log-odds of dementia",
      main = "",
      bty = "n",
      ylim = c(-0.7, 0.7),
      scheme =1)

dev.off ()


######remake figures with higher resolution#######

#for 1948

model_data1948 <- cox_ugu1948 %>%
  select(ID, dementia_any, TS6IITP, TS6ISTP, TS6IVOTP, RSSEX, RFSUN4, RMSUN4, parentedu_mean) %>%
  filter(!is.na(dementia_any)) %>%
  filter(!is.na(RSSEX)) %>%
  filter(!is.na(parentedu_mean)) %>%
  filter(!is.na(TS6IITP) & !is.na(TS6ISTP) & !is.na(TS6IVOTP))

#how many females and males in this model (RSSEX: 1 = Male, 2 = female)
table (model_data1948$RSSEX)

#check parent edu in this model

#frequenccies of parental education
table (model_data1948$RFSUN4)
table (model_data1948$RMSUN4)


summary (model_data1948$parentedu_mean)
sd (model_data1948$parentedu_mean)

#show the subtest correlations
cor(model_data1948$TS6IITP, model_data1948$TS6ISTP, use = "pairwise.complete.obs")
cor(model_data1948$TS6IITP, model_data1948$TS6IVOTP, use = "pairwise.complete.obs")
cor(model_data1948$TS6IVOTP, model_data1948$TS6ISTP, use = "pairwise.complete.obs")

#plot the distributions of scores from the subtests
library (ggplot2)
library (tidyr)

model_data1948_long <-model_data1948%>%
  select (TS6IVOTP, TS6ISTP, TS6IITP)%>%
  pivot_longer(everything(), names_to = "Test", values_to = "Score")%>%
  mutate(Test = recode(Test, 
                       "TS6IVOTP" = "Verbal",
                       "TS6ISTP" = "Spatial",
                       "TS6IITP" = "Inductive"))

#specify where to put plot
outdir <- "/safe/data/KBW/"

tiff(
  file.path(outdir, "subtest_distributions_age13.tif"),
  width = 7,
  height =5,
  units = "in",
  res = 300,
  compression ="lzw"
  )

print(
  ggplot(model_data1948_long, aes(x= Score, fill = Test, color =Test))+
  geom_density(alpha = 0.3, size =1)+
  labs(x = "Test score",
       y = "Density",
       title = "Distribution of subtest scores, age 13")+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette ="Dark2")+
  theme_minimal(base_size =14)+
  theme(
    panel.grid = element_blank(),
    panel.border =element_blank(),
    axis.line = element_line(),
    plot.title =element_text(size =16, face = "bold")
    
  ))

dev.off()

#or use ragg (for better text rendering?)
library(ragg)

agg_tiff(
  file.path(outdir, "subtest_distributions_age13_600_dpi.tif"),
  width = 7,
  height =5,
  units = "in",
  res = 600,
  compression ="lzw"
)

print(
  ggplot(model_data1948_long, aes(x= Score, fill = Test, color =Test))+
    geom_density(alpha = 0.3, size =1)+
    labs(x = "Test score",
         y = "Density",
         title = "Distribution of subtest scores, age 13")+
    scale_fill_brewer(palette = "Dark2")+
    scale_color_brewer(palette ="Dark2")+
    theme_minimal(base_size =14)+
    theme(
      panel.grid = element_blank(),
      panel.border =element_blank(),
      axis.line = element_line(),
      plot.title =element_text(size =16, face = "bold")
      
    ))

dev.off()


##now do som additional numbers for Lancet Helathy Longevity version Dec 21st

#check who has dementia from which source across all registry data, total cases seems 403
#check overlap for groups, 

#the blow shows that aomng teh 415 cases identified in all of UGU cross cohorts,
#inldung thsoe with other missing daat (so irrespectove of cognitive tests)
#329 were in Ugu48, of which 276 in NPR, and 53 in CDR only
#86 were in UGU53, of which 8 in CDR only,
#but remember, when doing teh merged analysis, we only have 376 dementia cases
"So for uG1948"
# 0-0, not in either, 
#0-1: 53 have no dementia in NPR, have in CDR only
# 1-0: 58: have dementia in NPR and not CDR, 
# 1-1:97 overlap, have in both;
# 1-NA: 121: have dementia in NPR, no CDR data (likely still alive at end of follow-u)
#NA-NA: 8795 -no info in either registry, people never Hospitalized for dementia/nit died/no dementia recorded at death
#(NPR only:228)   NA-1 (CDR only: 61) NA-NA  (none of them:4366)
table (paste (cox_ugu1948$NPR_dementia,
              cox_ugu1948$CDR_dementia, 
              sep="-"))

# for ugu53, corresponding numbers are:
#0-1:1532, 0-1:8, 1-0: 15, 1-1:29, 1-NA:34, NA-NA 8309
table (paste (cox_ugu1953$NPR_dementia,
              cox_ugu1953$CDR_dementia, 
              sep="-"))

#in full sample, regardless of cognitives tests
#of 11945 ugu48, 329 dementia cases and 11616 no
#of 9927 in ugu53, 86 dementia cases, and 9841 no
table (cox_ugu1948$dementia_any)
table (cox_ugu1953$dementia_any)

# check w/o missings:
model_data1948 <- cox_ugu1948 %>%
  select(ID, dementia_any, TS6IITP, TS6ISTP, TS6IVOTP, RSSEX, RFSUN4, RMSUN4, 
         NPR_dementia, CDR_dementia, parentedu_mean) %>%
  filter(!is.na(dementia_any)) %>%
  filter(!is.na(RSSEX)) %>%
  filter(!is.na(parentedu_mean)) %>%
  filter(!is.na(TS6IITP) & !is.na(TS6ISTP) & !is.na(TS6IVOTP))

table (model_data1948$dementia_any)

table (paste (model_data1948$NPR_dementia,
              model_data1948$CDR_dementia, 
              sep="-"))

# check w/o missing TS6IITP and parent edu:
model_data1948b <- cox_ugu1948 %>%
  select(ID, dementia_any, TS6IITP, TS6ISTP, TS6IVOTP, RSSEX, RFSUN4, RMSUN4, 
         NPR_dementia, CDR_dementia, parentedu_mean) %>%
  filter(!is.na(dementia_any)) %>%
  filter(!is.na(RSSEX)) %>%
  filter(!is.na(TS6IITP))

table (model_data1948b$dementia_any)

table (paste (model_data1948b$NPR_dementia,
              model_data1948b$CDR_dementia, 
              sep="-"))

#check for ugu1953
model_data1953 <- cox_ugu1953 %>%
  select(ID, dementia_any, TS6IITP, RSSEX,  
         NPR_dementia, CDR_dementia) %>%
  filter(!is.na(dementia_any)) %>%
  filter(!is.na(RSSEX)) %>%
  filter(!is.na(TS6IITP))

table (model_data1953$dementia_any)

table (paste (model_data1953$NPR_dementia,
              model_data1953$CDR_dementia, 
              sep="-"))
#check which version of R, is 4.3.3 (2024-02-29)
R.version.string

##########after Lancet desk rejection, do some checks for possible add ons ######

#first, we can considier adding parental edu for ugu1953
#will not be a real value, and no reason for rejection, bu still
#fix parental_edu also with different variable sin ugu1953
#recode 0s = no response to NA in father (QS68) and mother (QS69) edu
ugu1953 <- ugu1953 %>%
  mutate(
    QS68 = na_if(QS68, 0),
    QS69   = na_if(QS69, 0)
  )
#make mean parental edu, and use either val if one missing
#note that here, primary school = 1, junior secondary school =2,
#graduated or equivalent = 3, academic degree = 4, check: can it be directly comparable to 1948?
ugu1953 <- ugu1953 %>%
  mutate(
    parentedu_mean = rowMeans(select(., QS68, QS69), na.rm = TRUE)
  )

#standardize parental edu var
ugu1953 <- ugu1953%>%
  mutate(
    z_parentedu_mean = as.numeric(scale (parentedu_mean)),
  )

##########Add analyses on other somatic morbidities for JAMA ######

# I have run KBW_Charlson_R... and saved teh data, now I load that

charlson <-read.csv ("/safe/data/KBW/charlsondata.csv")
names(charlson)

#check
table (charlson$CCIunw, useNA="ifany")
table (charlson$CCIw, useNA="ifany")

#compute CCI-nodemntia
charlson$CCI_nodementia <- charlson$CCIunw - charlson$Dementia
charlson$CCI_event_nodementia <-ifelse (charlson$CCI_nodementia >0,1,0)

table (charlson$CCI_event_nodementia)

#merge, keeping all 1948 cases
cox_ugu1948 <-merge(
  cox_ugu1948,
  charlson,
  by = "ID",
  all.x =TRUE
)



#I need to get end_age for charlson non-dementia events

#first identify all non-dementia date columns

#my signs: §|~ %>% $ {} [ ] \\

#list all date vars
date_vars <-names(cox_ugu1948)[grepl("^date\\.", names(cox_ugu1948))]
date_vars

#remove dementia
date_vars_nodementia <- date_vars [date_vars != "date.Dementia"]
#compute earliest non-dementia date
cox_ugu1948$first_CCI_date <-
  do.call(pmin, c(cox_ugu1948 [date_vars_nodementia], na.rm = TRUE))

#create event indicator
cox_ugu1948$CCI_event_timebased <-
  ifelse(is.na(cox_ugu1948$first_CCI_date),0,1)

#compute CCI_date and end_age(CCI or censoring)
cox_ugu1948 <-cox_ugu1948%>%
  mutate(
    first_CCI_date_date = as.Date(as.character(first_CCI_date), "%Y%m%d"),
      event_date_CCI = if_else(
        CCI_event_timebased ==1L, 
        first_CCI_date_date,
       censor_date),
    #bug_fix April 9th
    #end_age_CCI = as.numeric((event_date -birth_date)/365.25)
    end_age_CCI = as.numeric((event_date_CCI -birth_date)/365.25)
  )


#######cause of reused personnummere#####
####check that charlson_date is not after CDR_death_date####
#If it says FALSE, we+re clean, but it says TRUE, for 6
#check if they ar ein the models

cox_ugu1948Y <- cox_ugu1948 %>%
  #select(ID, dementia_any, sex01, z_parentedu_mean,
  #       z_TS6IITP, z_TS6ISTP, z_TS6IVOTP,
  #       start_age, end_age,) %>%
  filter(!is.na(start_age)) %>%
  filter(!is.na(end_age)) %>%
  filter(!is.na(dementia_any)) %>%
  filter(!is.na(sex01)) %>%
  filter(!is.na(z_parentedu_mean)) %>%
  filter(!is.na(z_TS6IITP))%>%
  filter(!is.na(z_TS6ISTP))%>%
  filter(!is.na(z_TS6IVOTP))%>%
filter(!is.na(CCI_event_timebased))

any(cox_ugu1948Y$first_CCI_date_date > cox_ugu1948Y$CDR_death_date, na.rm = TRUE)
sum(cox_ugu1948Y$first_CCI_date_date > cox_ugu1948Y$CDR_death_date, na.rm = TRUE)

#from this, I get 4 cases with discrepancy of 1-2 days between death date and CCI date
#but also two cases with several years, those two need to go out
cox_ugu1948YX <- cox_ugu1948Y%>%
  filter (first_CCI_date_date > CDR_death_date)%>%
  select(ID, first_CCI_date_date, CDR_death_date)

#my signs: §|~ %>% $ {} [ ] \\

#I write some lines to exclude all persons that have 3 days or more discrepancy
#where CCI dates in date.* format is more than 2 days from death (utdatum)
#this here should select all Charlson date variables(date.*)
#if event date > death date + 3 days, I set to NA
#but thsi si scary, and does not work, I think it drops all who did not die
#cox_ugu1948[ , grep("date\\.", names(cox_ugu1948))] <-
#  lapply(cox_ugu1948[ , grep("date\\.", names(cox_ugu1948))],
#         function(x)  {x [x > cox_ugu1948$CDR_death_date + 3] <- NA; x })

#I will ratehr do this in tidyverse

#then I redo my stuff so far

#list all date vars
date_vars <-names(cox_ugu1948)[grepl("^date\\.", names(cox_ugu1948))]
date_vars

#identify problem IDs
problem_ids <- cox_ugu1948 %>%
  filter(!is.na(CDR_death_date)) %>%
  filter(if_any(all_of(date_vars),
                ~as.Date(as.character(.x), "%Y%m%d") > CDR_death_date +3 ))
#count
nrow(problem_ids)

#my signs: §|~ %>% $ {} [ ] \\
#set impossible cci dates, more than 3 days to na (across diagnoses)
cox_ugu1948 <- cox_ugu1948%>%
  mutate(across(all_of(date_vars),
                ~ {
                  #convert numeric YYYYMMDD to Date for comparison
                  d <-as.Date(as.character(.x), "%Y%m%d")
                  #identify impossible post-mortem dates
                  bad <-!is.na(CDR_death_date) &
                    !is.na(d) &
                    d > CDR_death_date +3
                  #set those to NA in original numeric variable
                  .x[bad] <- NA
                  #return cleaned numeric vector
                  .x}))

#remove dementia
date_vars_nodementia <- date_vars [date_vars != "date.Dementia"]
#compute earliest non-dementia date
cox_ugu1948$first_CCI_date <-
  do.call(pmin, c(cox_ugu1948 [date_vars_nodementia], na.rm = TRUE))

#create event indicator
cox_ugu1948$CCI_event_timebased <-
  ifelse(is.na(cox_ugu1948$first_CCI_date),0,1)

#compute CCI_date and end_age(CCI or censoring)
cox_ugu1948 <-cox_ugu1948%>%
  mutate(
    first_CCI_date_date = as.Date(as.character(first_CCI_date), "%Y%m%d"),
    event_date_CCI = if_else(
      CCI_event_timebased ==1L, 
      first_CCI_date_date,
      censor_date),
    #bug fix April9th
    #end_age_CCI = as.numeric((event_date -birth_date)/365.25)
    end_age_CCI = as.numeric((event_date_CCI -birth_date)/365.25)
  )



#redo cox

cox_no_dementia <-coxph(formula =
                  Surv(start_age, end_age_CCI, CCI_event_timebased) ~
                  z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                  z_parentedu_mean +sex01,
                data = cox_ugu1948)

summary(cox_no_dementia)

cox_no_dementia_edu <-coxph(formula =
                          Surv(start_age, end_age_CCI, CCI_event_timebased) ~
                          z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                          z_parentedu_mean +z_edu_FOB_1990 +sex01,
                        data = cox_ugu1948)

summary(cox_no_dementia_edu)

#compare again
cox_all <-coxph(formula =
                  Surv(start_age, end_age, dementia_any) ~
                  z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                  z_parentedu_mean +sex01,
                data = cox_ugu1948
)

summary(cox_all)


#there is a difference in magnitude of risk, so I will formally test if it is statistically significant
#stack the data
#d1 <-transform(cox_ugu1948, event = dementia_any, end = end_age, outcome =1)
#d0 <-transform(cox_ugu1948, event = CCI_event_timebased, end =end_age_CCI, outcome = 0) 
#stacked <-rbind(d1, d0)

#run interaction model, there is not a statistically significant interaction
# I will not include this in teh paper, it is just to see
#cox_compare <-coxph(formula =
#                      Surv(start_age, end, event) ~
#                      z_TS6IITP * outcome +
#                      z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
#                      z_parentedu_mean +sex01 +
#                      strata(outcome),
#                    data = stacked)

#summary(cox_compare)

#cox_compare <-coxph(formula =
#                      Surv(start_age, end, event) ~
#                      z_TS6IITP * outcome +
#                      z_parentedu_mean +sex01 +
#                      strata(outcome),
#                    data = stacked)

#summary(cox_compare)

######analyze cardiovascular events specifically######

#check how many events I have - this includes people not havig cogtests and all vars
table (!is.na(cox_ugu1948$date.Cerebrovascular_disease))
table (!is.na(cox_ugu1948$date.Myocardial_infarction))
table (!is.na(cox_ugu1948$date.Congestive_heart_failure))

#make events and event ages 
#check with CCI defintion of dementia versus ours, 
#but remember ours include quite a few cases from CDR also (n = 43) 
#will vary from that used (only 242 cases, p = .0575, HR = .86)
#They include:
#ICD 9: 290, 294B, 331A, 331B, 331C, 331X 
#(think 294B is persistent mental disorders due to brain damage, and 331x is other cerebral degenerations)
#ICD10:F00, F01, F02, F03, F051, G30, G311, G319
#I think F051 is delirum superimposed on dementia
#and G311 is sSenile degenration of brain, not elsewhere classified

#We include (n = 287 cases):
#ICD9:290 and 331A, 331B, 331C
#ICD10:F00, F01, F02, F03,G30, G310, G318, G319
#our G310 is frontotemporal dementia
#and our G318 is other specified degenerative diseases

#cox_ugu1948$CCI_dementia_event <-ifelse(!is.na(cox_ugu1948$date.Dementia), 1, 0)

#cox_ugu1948$CCI_dementia_age <- as.numeric(
#  as.Date(as.character(cox_ugu1948$date.Dementia), "%Y%m%d") -cox_ugu1948$birth_date)/365.25

#cox_ugu1948$CCI_dementia_endage <- ifelse(cox_ugu1948$CCI_dementia_event==1,
#                                             cox_ugu1948$CCI_dementia_age,
#                                             cox_ugu1948$end_age)

#cox_CCI_dementia <-coxph(formula =
#                              Surv(start_age, CCI_dementia_endage, CCI_dementia_event) ~
#                              z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
#                              z_parentedu_mean +sex01,
#                            data = cox_ugu1948)

#summary(cox_CCI_dementia)

#cerebrovascular_disease
cox_ugu1948$cerebrovascular_event <-ifelse(!is.na(cox_ugu1948$date.Cerebrovascular_disease), 1, 0)

cox_ugu1948$cerebrovascular_age <- as.numeric(
  as.Date(as.character(cox_ugu1948$date.Cerebrovascular_disease), "%Y%m%d") -cox_ugu1948$birth_date)/365.25

cox_ugu1948$cerebrovascular_endage <- ifelse(cox_ugu1948$cerebrovascular_event==1,
                                             cox_ugu1948$cerebrovascular_age,
                                             cox_ugu1948$end_age)
  
cox_cerebrovascular <-coxph(formula =
                          Surv(start_age, cerebrovascular_endage, cerebrovascular_event) ~
                          z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                          z_parentedu_mean +sex01,
                        data = cox_ugu1948)

summary(cox_cerebrovascular)


#myocardial infarction
cox_ugu1948$myocardial_event <-ifelse(!is.na(cox_ugu1948$date.Myocardial_infarction), 1, 0)

cox_ugu1948$myocardial_age <- as.numeric(
  as.Date(as.character(cox_ugu1948$date.Myocardial_infarction), "%Y%m%d") -cox_ugu1948$birth_date)/365.25

cox_ugu1948$myocardial_endage <- ifelse(cox_ugu1948$myocardial_event==1,
                                             cox_ugu1948$myocardial_age,
                                             cox_ugu1948$end_age)

cox_myocardial <-coxph(formula =
                              Surv(start_age, myocardial_endage, myocardial_event) ~
                              z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                              z_parentedu_mean +sex01,
                            data = cox_ugu1948)

summary(cox_myocardial)


#congestive heart failure
cox_ugu1948$heartfailure_event <-ifelse(!is.na(cox_ugu1948$date.Congestive_heart_failure), 1, 0)

cox_ugu1948$heartfailure_age <- as.numeric(
  as.Date(as.character(cox_ugu1948$date.Congestive_heart_failure), "%Y%m%d") -cox_ugu1948$birth_date)/365.25

cox_ugu1948$heartfailure_endage <- ifelse(cox_ugu1948$heartfailure_event==1,
                                        cox_ugu1948$heartfailure_age,
                                        cox_ugu1948$end_age)

cox_heartfailure <-coxph(formula =
                         Surv(start_age, heartfailure_endage, heartfailure_event) ~
                         z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                         z_parentedu_mean +sex01,
                       data = cox_ugu1948)

summary(cox_heartfailure)


#my signs: §|~ %>% $ {} [ ] \\
#any cardiovascular event
#define composite

cox_ugu1948$cvd_any_event <-ifelse(
    cox_ugu1948$cerebrovascular_event==1 |
    cox_ugu1948$myocardial_event==1 |
    cox_ugu1948$heartfailure_event==1,1,0)

#define earliest of the three
cox_ugu1948$cvd_any_age <-pmin(
  cox_ugu1948$cerebrovascular_age,
    cox_ugu1948$myocardial_age,
    cox_ugu1948$heartfailure_age,
    na.rm = TRUE)

#define end_age for composite
cox_ugu1948$cvd_any_endage <-ifelse(
  cox_ugu1948$cvd_any_event==1,
    cox_ugu1948$cvd_any_age,
    cox_ugu1948$end_age)

#model child cog - composite CVD
cox_any_CVD <-coxph(formula =
            Surv(start_age, cvd_any_endage, cvd_any_event) ~
            z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
            z_parentedu_mean +sex01,
            data = cox_ugu1948)

summary(cox_any_CVD)

cox_any_CVD_edu <-coxph(formula =
                      Surv(start_age, cvd_any_endage, cvd_any_event) ~
                      z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                      z_parentedu_mean + z_edu_FOB_1990 +sex01,
                    data = cox_ugu1948)

summary(cox_any_CVD_edu)


######run diabetes too######

date_vars
#diabetes with and without chronic complications combined

#compute earliest diabetes date
cox_ugu1948$diabetes_date <- pmin(
  cox_ugu1948$date.Diabetes_without_chronic_complication,
  cox_ugu1948$date.Diabetes_with_chronic_complication,
  na.rm = TRUE)

#my signs: §|~ %>% $ {} [ ] \\

cox_ugu1948$diabetes_date[is.infinite(cox_ugu1948$diabetes_date)] <-NA

cox_ugu1948$diabetes_event <-ifelse(!is.na(cox_ugu1948$diabetes_date), 1, 0)

cox_ugu1948$diabetes_age <- as.numeric(
  as.Date(as.character(cox_ugu1948$diabetes_date), "%Y%m%d") -cox_ugu1948$birth_date)/365.25

cox_ugu1948$diabetes_endage <- ifelse(cox_ugu1948$diabetes_event==1,
                                          cox_ugu1948$diabetes_age,
                                          cox_ugu1948$end_age)


#######Build time-varying covariate dataset using tmerge######
# April 2026: Fix time-ordering so morbidity covariates only "switch on"
# at the age they actually occur, properly handling dementia as competing event.
# Models where CVD/diabetes/CCI is the OUTCOME are unchanged above.
# Models below use tmerge long format: each person contributes intervals
# with covariate = 0 before their event and = 1 after.

# Fix Inf in cvd_any_age (pmin returns Inf when all components are NA with na.rm=TRUE)
cox_ugu1948$cvd_any_age[is.infinite(cox_ugu1948$cvd_any_age)] <- NA

# Compute CCI event age as a standalone variable for tdc()
cox_ugu1948 <- cox_ugu1948 %>%
  mutate(
    CCI_event_age = if_else(
      CCI_event_timebased == 1L,
      as.numeric((first_CCI_date_date - birth_date) / 365.25),
      NA_real_)
  )

# Build long-format dataset: one row per person-interval, timescale = age
cox_ugu1948_tv <- tmerge(
  data1  = cox_ugu1948,
  data2  = cox_ugu1948,
  id     = ID,
  tstart = start_age,
  tstop  = end_age,
  dementia_tv = event(end_age, dementia_any)
)

# Add each morbidity as a time-dependent covariate (tdc)
# tdc(t): covariate switches from 0 to 1 at age t; NA t = never switches
cox_ugu1948_tv <- tmerge(cox_ugu1948_tv, cox_ugu1948, id = ID,
  CCI_tv             = tdc(CCI_event_age))
cox_ugu1948_tv <- tmerge(cox_ugu1948_tv, cox_ugu1948, id = ID,
  cvd_any_tv         = tdc(cvd_any_age))
cox_ugu1948_tv <- tmerge(cox_ugu1948_tv, cox_ugu1948, id = ID,
  diabetes_tv        = tdc(diabetes_age))
cox_ugu1948_tv <- tmerge(cox_ugu1948_tv, cox_ugu1948, id = ID,
  cerebrovascular_tv = tdc(cerebrovascular_age))
cox_ugu1948_tv <- tmerge(cox_ugu1948_tv, cox_ugu1948, id = ID,
  myocardial_tv      = tdc(myocardial_age))
cox_ugu1948_tv <- tmerge(cox_ugu1948_tv, cox_ugu1948, id = ID,
  heartfailure_tv    = tdc(heartfailure_age))

# Diagnostics
nrow(cox_ugu1948_tv)           # should exceed nrow(cox_ugu1948)
attr(cox_ugu1948_tv, "tcount") # split-row summary per covariate

#######Dementia models with time-varying morbidity covariates######
# All use Surv(tstart, tstop, dementia_tv) on cox_ugu1948_tv

# Somatic morbidity (CCI, time-varying) predicting dementia
cox_dementia_any_somatic <- coxph(formula =
                               Surv(tstart, tstop, dementia_tv) ~
                               z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                               z_parentedu_mean + sex01 +
                               CCI_tv,
                             data = cox_ugu1948_tv)

summary(cox_dementia_any_somatic)

cox_dementia_any_somatic_edu <- coxph(formula =
                                   Surv(tstart, tstop, dementia_tv) ~
                                   z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                                   z_parentedu_mean + z_edu_FOB_1990 + sex01 +
                                   CCI_tv,
                                 data = cox_ugu1948_tv)

summary(cox_dementia_any_somatic_edu)

# Cerebrovascular disease (time-varying) predicting dementia
cox_dementia_CVD <- coxph(formula =
                  Surv(tstart, tstop, dementia_tv) ~
                  z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                  z_parentedu_mean + sex01 +
                  cerebrovascular_tv,
                data = cox_ugu1948_tv)

summary(cox_dementia_CVD)

# Myocardial infarction (time-varying) predicting dementia
cox_dementia_myocardial <- coxph(formula =
                           Surv(tstart, tstop, dementia_tv) ~
                           z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                           z_parentedu_mean + sex01 +
                           myocardial_tv,
                         data = cox_ugu1948_tv)

summary(cox_dementia_myocardial)

# Congestive heart failure (time-varying) predicting dementia
cox_dementia_heartfailure <- coxph(formula =
                                  Surv(tstart, tstop, dementia_tv) ~
                                  z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                                  z_parentedu_mean + sex01 +
                                  heartfailure_tv,
                                data = cox_ugu1948_tv)

summary(cox_dementia_heartfailure)

# Composite CVD (time-varying) predicting dementia
cox_dementia_any_CVD <- coxph(formula =
            Surv(tstart, tstop, dementia_tv) ~
            z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
            z_parentedu_mean + sex01 +
            cvd_any_tv,
            data = cox_ugu1948_tv)

summary(cox_dementia_any_CVD)

cox_dementia_any_CVD_edu <- coxph(formula =
                               Surv(tstart, tstop, dementia_tv) ~
                               z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                               z_parentedu_mean + z_edu_FOB_1990 + sex01 +
                               cvd_any_tv,
                             data = cox_ugu1948_tv)

summary(cox_dementia_any_CVD_edu)

# Diabetes (time-varying) predicting dementia
cox_dementia_diabetes <- coxph(formula =
                                    Surv(tstart, tstop, dementia_tv) ~
                                    z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                                    z_parentedu_mean + sex01 +
                                    diabetes_tv,
                                  data = cox_ugu1948_tv)

summary(cox_dementia_diabetes)

cox_dementia_diabetes_edu <- coxph(formula =
                                Surv(tstart, tstop, dementia_tv) ~
                                z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                                z_parentedu_mean + z_edu_FOB_1990 + sex01 +
                                diabetes_tv,
                              data = cox_ugu1948_tv)

summary(cox_dementia_diabetes_edu)

cox_diabetes <-coxph(formula =
                           Surv(start_age, diabetes_endage, diabetes_event) ~
                           z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                           z_parentedu_mean +sex01,
                         data = cox_ugu1948)

summary(cox_diabetes)

cox_diabetes_edu <-coxph(formula =
                       Surv(start_age, diabetes_endage, diabetes_event) ~
                       z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                       z_parentedu_mean + z_edu_FOB_1990 +sex01,
                     data = cox_ugu1948)

summary(cox_diabetes_edu)


table(cox_ugu1948$Diabetes_without_chronic_complication)
table(cox_ugu1948$Diabetes_with_chronic_complication)



#######Make higher res new figure 3 with own edu Jan 6th for paper########
#make a forest plot of the standardized effects:
#extract ORs and CIs
HRs <-exp(coef(cox_dementia_any_CVD_edu))
CIs <-exp(confint.default(cox_dementia_any_CVD_edu))
#make clean dataframe for plotting
#recoded RSSEX, which was  1= Male, and 2 = female, to 0 = male and 1 is female
plot_df <-data.frame(
  Predictor = c("Inductive reasoning", "Spatial ability", "Verbal ability",
                "Parental education", "Midlife education", "Cardiovascular morbidity",
                "Sex (male vs. female)"),
  HR =HRs,
  CI_low = CIs [,1],
  CI_high = CIs[,2]
)

library (ggplot2)



#inductive top
plot_df$Predictor <- factor(plot_df$Predictor,
                            levels = c("Sex (male vs. female)",
                                       "Cardiovascular morbidity",
                                       "Midlife education",
                                       "Parental education",
                                       "Verbal ability",
                                       "Spatial ability", 
                                       "Inductive reasoning"
                            ))

#I want one color for the predictors and one for the covariates, 
plot_df$Type <- factor(
  c("Cognitive", "Cognitive", "Cognitive",
    "Covariate", "Covariate", "Covariate", "Covariate"),
  levels = c("Cognitive","Covariate"))
#specify where to put plot
outdir <- "/safe/data/KBW/"

library(ragg)

agg_tiff(
  file.path(outdir, "Figure3_w_CVD_1990_edu_HR_plot_600_dpi.tif"),
  width = 12,
  height =5,
  units = "in",
  res = 600,
  compression ="lzw"
)

print(
  ggplot(plot_df, aes(x= Predictor, y = HR, color = Type))+
    geom_point(size =3)+
    geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.1, linewidth =0.8) +
    geom_hline(yintercept = 1, linetype = "dashed")+
    coord_flip()+
    #scale_y_log10(breaks = c(0.5, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 2, 2.5, 3))+
    scale_y_continuous(limits = c(0.5, 3),
                       breaks = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0))+
    theme_classic(base_size =18)+
    scale_color_manual(values = c("Covariate" = "#d62728", #red
                                  "Cognitive" = "#1f77b4")) + #blue
    labs(title = "Childhood cognitive scores, education, cardiovascular \nmorbidity and dementia risk",
         y = "Hazard ratio (per 1SD)",
         x = "",
         color = "",
        theme (plot.title = element_text(hjust = 0.5),
               plot.totle.positon = "plot")
    ))


#\n
#reset graphics device (warnings)
dev.off()

#warnings()

###make a table:

#what I hve already run:
# I filter out make sure we have same dataset on all vars
model_data48_table <- cox_ugu1948 %>%
  #select(ID, dementia_any, sex01, z_parentedu_mean,
        # z_TS6IITP, z_TS6ISTP, z_TS6IVOTP,
        # start_age, end_age,) %>%
  filter(!is.na(start_age)) %>%
  filter(!is.na(end_age)) %>%
  filter(!is.na(dementia_any)) %>%
  filter(!is.na(sex01)) %>%
  filter(!is.na(z_parentedu_mean)) %>%
  filter(!is.na(z_TS6IITP))%>%
  filter(!is.na(z_TS6ISTP))%>%
  filter(!is.na(z_TS6IVOTP))
#my signs: §|~ %>% $ {} [ ] \\
summary (model_data48_table$end_age)


#make death and dementia ages

model_data48_table <- model_data48_table%>%
  mutate(
    dementia_age = if_else(
      dementia_any ==1L,
      as.numeric((dementia_date-birth_date)/365.25),
                 NA_real_),
      death_age = if_else(
        !is.na(CDR_death_date),
        as.numeric((CDR_death_date-birth_date)/365.25),
                   NA_real_),
        followup_years = end_age -start_age)
#check

summary(model_data48_table$dementia_age, na.rm =TRUE)
summary(model_data48_table$death_age, na.rm =TRUE)



sum(!is.na(model_data48_table$dementia_age)) == sum(model_data48_table$dementia_any ==1)

#my signs: §|~ %>% $ {} [ ] \\
#table_template
make_jama_table <-function(data, group_var, event_age_var) {
  data%>%
    group_by({{group_var}})%>%
    summarise(
      n=n(),
      female = paste0(
        sum(sex01==1, na.rm =TRUE),
        "(",
        round(mean(sex01==1, na.rm =TRUE)*100,1),
        "%)"),
      inductive =paste0(
        round(mean(TS6IITP, na.rm =TRUE),2),
        "(",
        round(sd(TS6IITP, na.rm=TRUE),2),
        ")"),
      spatial =paste0(
        round(mean(TS6ISTP, na.rm =TRUE),2),
        "(",
        round(sd(TS6ISTP, na.rm=TRUE),2),
        ")"),
      verbal =paste0(
        round(mean(TS6IVOTP, na.rm =TRUE),2),
        "(",
        round(sd(TS6IVOTP, na.rm=TRUE),2),
        ")"),
      parentedu =paste0(
        round(mean(parentedu_mean, na.rm =TRUE),2),
        "(",
        round(sd(parentedu_mean, na.rm=TRUE),2),
        ")"),
      midlifeedu =paste0(
        round(mean(edu_FOB_1990, na.rm =TRUE),2),
        "(",
        round(sd(edu_FOB_1990, na.rm=TRUE),2),
        ")"),
      event_age =paste0(
        round(median({{event_age_var}}, na.rm =TRUE),2),
        "(",
        round(IQR({{event_age_var}}, na.rm=TRUE),2),
        ")"),
      death_age =paste0(
        round(median(death_age, na.rm =TRUE),2),
        "(",
        round(IQR(death_age, na.rm=TRUE),2),
        ")"),
      followup =paste0(
        round(median(followup_years, na.rm =TRUE),2),
        "(",
        round(IQR(followup_years, na.rm=TRUE),2),
        ")"
      ),
  .groups="drop")} 

#run for dementia
make_jama_table(model_data48_table, dementia_any, dementia_age)
#make var for missin midlife edu
model_data48_table <- model_data48_table%>%
  mutate(midlife_missing = if_else(is.na(edu_FOB_1990), 1,0))
#run for edu_missing
make_jama_table(model_data48_table, midlife_missing, dementia_age)
#run for major somatic morbidity
make_jama_table(model_data48_table, CCI_event_timebased, dementia_age)
#run for major somatic morbidity
make_jama_table(model_data48_table, cvd_any_event, dementia_age)
#run for diabetes
make_jama_table(model_data48_table, diabetes_event, dementia_age)

#do reliability corrections on z_scores
model_data48_table <-model_data48_table %>%
  mutate(z_TS6IITP_corr =z_TS6IITP / sqrt(0.92),
         z_TS6ISTP_corr =z_TS6ISTP / sqrt(0.87),
         z_TS6IVOTP_corr =z_TS6IVOTP / sqrt(0.92))

cox_all_corr <-coxph(formula =
                  Surv(start_age, end_age, dementia_any) ~
                  z_TS6IITP_corr + z_TS6ISTP_corr + z_TS6IVOTP_corr +
                  z_parentedu_mean +sex01,
                data = model_data48_table)

summary(cox_all_corr)

#or restandardize after correction so HRs are interpretable per 1 SD
model_data48_table <-model_data48_table %>%
  mutate(z_TS6IITP_corr =as.numeric(scale(z_TS6IITP / sqrt(0.92))),
         z_TS6ISTP_corr =as.numeric(scale(z_TS6ISTP / sqrt(0.87))),
         z_TS6IVOTP_corr =as.numeric(scale(z_TS6IVOTP / sqrt(0.92))))

cox_all_corr <-coxph(formula =
                       Surv(start_age, end_age, dementia_any) ~
                       z_TS6IITP_corr + z_TS6ISTP_corr + z_TS6IVOTP_corr +
                       z_parentedu_mean +sex01,
                     data = model_data48_table)

summary(cox_all_corr)

#address possible issue of competing risk by death

#my signs: §|~ %>% $ {} [ ] \\

library(cmprsk)
#create time variable (end_age - start_age), 
#not really needed, as I defined censor date as detah or endo f follow-up,
#and then event_date as dementia, death, or admin end/end-of follow-up,
#where en_age is dementia aor censoring (and censoring is end of followup or death)
#model_data48_table$ftime <-model_data48_table$end_age -model_data48_table$start_age

#create competing risk status variable for thsi dataset 
#(1 if dementia before death, 2 if death occurs before dementia, 0 otherwise)

model_data48_table <-model_data48_table%>%
  mutate(
    status_cr = case_when(
      dementia_any ==1L ~ 1,           #dementia 
    !is.na(CDR_death_date)~ 2,         #death before dementia
  TRUE ~ 0                             #alive at end of follow-up
  ),  
  
  ftime = end_age-start_age)

#run Fine-Gray Model
fg_model <- crr(
  ftime = model_data48_table$ftime,
  fstatus = model_data48_table$status_cr,
  cov1 = as.matrix(model_data48_table[, c("z_TS6IITP",
                                        "z_parentedu_mean",
                                        "sex01")]))

summary (fg_model)

#warninsg pertian to different labels across ugu1948 and ugu1953
warnings ()
