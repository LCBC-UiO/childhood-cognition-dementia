# ----- Buka descriptives: deaths and age at dementia diagnosis ---------------
# Drop in after the cvd_any_event / diabetes_event definitions in
# Cox_prep_8_UGU_paper1_Feb_2026.R. Uses cox_ugu1948 with variables you
# already have: dementia_any, dementia_date, CDR_death_date, birth_date,
# CCI_event_timebased, cvd_any_event, diabetes_event.

library(dplyr)

# 1. Death indicator and ages -------------------------------------------------
cox_ugu1948 <- cox_ugu1948 %>%
  mutate(
    died          = as.integer(!is.na(CDR_death_date)),
    age_at_death  = as.numeric(CDR_death_date - birth_date) / 365.25,
    age_at_dem_dx = as.numeric(dementia_date  - birth_date) / 365.25
  )

# 2. Helper: per-exposure summary --------------------------------------------
desc_by <- function(df, group_var) {
  df %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      n            = n(),
      n_dead       = sum(died),
      pct_dead     = round(100 * mean(died), 1),
      age_death_md = round(median(age_at_death,  na.rm = TRUE), 2),
      age_death_iqr= round(IQR   (age_at_death,  na.rm = TRUE), 2),
      n_dem        = sum(dementia_any),
      age_dem_md   = round(median(age_at_dem_dx, na.rm = TRUE), 2),
      age_dem_iqr  = round(IQR   (age_at_dem_dx, na.rm = TRUE), 2),
      .groups = "drop"
    ) %>%
    rename(group = 1) %>%
    mutate(exposure = group_var, .before = 1)
}

# 3. Run for the four exposure splits -----------------------------------------
desc_tbl <- bind_rows(
  desc_by(cox_ugu1948, "dementia_any"),
  desc_by(cox_ugu1948, "CCI_event_timebased"),  # major somatic morbidity
  desc_by(cox_ugu1948, "cvd_any_event"),         # composite CVD
  desc_by(cox_ugu1948, "diabetes_event")         # diabetes
)

print(desc_tbl, n = Inf)
# write.csv(desc_tbl, "Buka_descriptives_byExposure.csv", row.names = FALSE)
