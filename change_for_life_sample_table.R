library(tidyverse)

# ── Sample characteristics table ──────────────────────────────────────────────
# Produces Table 1: characteristics by survival status (Survived vs Died)

make_sample_table <- function(data) {
  
  data %>%
    mutate(status = if_else(death == 0, "Survived", "Died")) %>%
    group_by(status) %>%
    summarise(
      n                = n(),
      cohort_1948_n    = sum(cohort == 1948),
      cohort_1948_pct  = mean(cohort == 1948) * 100,
      cohort_1953_n    = sum(cohort == 1953),
      cohort_1953_pct  = mean(cohort == 1953) * 100,
      child_g_mean     = mean(child_g_z),
      child_g_sd       = sd(child_g_z),
      army_g_mean      = mean(army_g_z),
      army_g_sd        = sd(army_g_z),
      resid_mean       = mean(army_g_resid),
      resid_sd         = sd(army_g_resid),
      pedu_mean        = mean(parentedu_mean_z),
      pedu_sd          = sd(parentedu_mean_z),
      age_death_mean   = mean(stop_age[death == 1]),
      age_death_sd     = sd(stop_age[death == 1]),
      followup_mean    = mean(stop_age - start_age),
      followup_sd      = sd(stop_age - start_age),
      .groups = "drop"
    )
}

# Cause-specific death counts
make_cause_counts <- function(data) {
  data %>%
    mutate(status = if_else(death == 0, "Survived", "Died")) %>%
    group_by(status) %>%
    summarise(
      n              = n(),
      cancer_n       = sum(ev_cancer == 1),
      cancer_pct     = mean(ev_cancer == 1) * 100,
      cvd_n          = sum(ev_cvd == 1),
      cvd_pct        = mean(ev_cvd == 1) * 100,
      resp_n         = sum(ev_resp == 1),
      resp_pct       = mean(ev_resp == 1) * 100,
      external_n     = sum(ev_external == 1),
      external_pct   = mean(ev_external == 1) * 100,
      accidental_n   = sum(ev2_accidental == 1),
      accidental_pct = mean(ev2_accidental == 1) * 100,
      suicide_n      = sum(ev_suicide == 1),
      suicide_pct    = mean(ev_suicide == 1) * 100,
      other_n        = sum(ev_other == 1),
      other_pct      = mean(ev_other == 1) * 100,
      .groups = "drop"
    )
}

# Total column
make_totals <- function(data) {
  data %>%
    summarise(
      status           = "Total",
      n                = n(),
      cohort_1948_n    = sum(cohort == 1948),
      cohort_1948_pct  = mean(cohort == 1948) * 100,
      cohort_1953_n    = sum(cohort == 1953),
      cohort_1953_pct  = mean(cohort == 1953) * 100,
      child_g_mean     = mean(child_g_z),
      child_g_sd       = sd(child_g_z),
      army_g_mean      = mean(army_g_z),
      army_g_sd        = sd(army_g_z),
      resid_mean       = mean(army_g_resid),
      resid_sd         = sd(army_g_resid),
      pedu_mean        = mean(parentedu_mean_z),
      pedu_sd          = sd(parentedu_mean_z),
      age_death_mean   = mean(stop_age[death == 1]),
      age_death_sd     = sd(stop_age[death == 1]),
      followup_mean    = mean(stop_age - start_age),
      followup_sd      = sd(stop_age - start_age)
    )
}

# Run
tab      <- make_sample_table(model_army_all)
causes   <- make_cause_counts(model_army_all)
totals   <- make_totals(model_army_all)

# Print readable summary
cat("\n=== SAMPLE CHARACTERISTICS ===\n")
print(tab)
cat("\n=== CAUSE-SPECIFIC DEATHS ===\n")
print(causes)
cat("\n=== TOTALS ===\n")
print(totals)

# Format for table (paste into ms)
format_mean_sd <- function(m, s) sprintf("%.2f (%.2f)", m, s)
format_n_pct   <- function(n, p) sprintf("%d (%.1f%%)", n, p)

cat("\n=== FORMATTED TABLE 1 ===\n")
cat(sprintf("%-40s %-25s %-25s %-25s\n",
    "Characteristic",
    paste0("Survived (n=", tab$n[tab$status=="Survived"], ")"),
    paste0("Died (n=", tab$n[tab$status=="Died"], ")"),
    paste0("Total (n=", totals$n, ")")))
cat(strrep("-", 115), "\n")

rows <- list(
  c("Cohort 1948, n (%)",
    format_n_pct(tab$cohort_1948_n[tab$status=="Survived"], tab$cohort_1948_pct[tab$status=="Survived"]),
    format_n_pct(tab$cohort_1948_n[tab$status=="Died"],     tab$cohort_1948_pct[tab$status=="Died"]),
    format_n_pct(totals$cohort_1948_n, totals$cohort_1948_pct)),
  c("Cohort 1953, n (%)",
    format_n_pct(tab$cohort_1953_n[tab$status=="Survived"], tab$cohort_1953_pct[tab$status=="Survived"]),
    format_n_pct(tab$cohort_1953_n[tab$status=="Died"],     tab$cohort_1953_pct[tab$status=="Died"]),
    format_n_pct(totals$cohort_1953_n, totals$cohort_1953_pct)),
  c("Child g, mean (SD)",
    format_mean_sd(tab$child_g_mean[tab$status=="Survived"], tab$child_g_sd[tab$status=="Survived"]),
    format_mean_sd(tab$child_g_mean[tab$status=="Died"],     tab$child_g_sd[tab$status=="Died"]),
    format_mean_sd(totals$child_g_mean, totals$child_g_sd)),
  c("Army g, mean (SD)",
    format_mean_sd(tab$army_g_mean[tab$status=="Survived"], tab$army_g_sd[tab$status=="Survived"]),
    format_mean_sd(tab$army_g_mean[tab$status=="Died"],     tab$army_g_sd[tab$status=="Died"]),
    format_mean_sd(totals$army_g_mean, totals$army_g_sd)),
  c("Adolescent change (residual), mean (SD)",
    format_mean_sd(tab$resid_mean[tab$status=="Survived"], tab$resid_sd[tab$status=="Survived"]),
    format_mean_sd(tab$resid_mean[tab$status=="Died"],     tab$resid_sd[tab$status=="Died"]),
    format_mean_sd(totals$resid_mean, totals$resid_sd)),
  c("Parental education, mean (SD)",
    format_mean_sd(tab$pedu_mean[tab$status=="Survived"], tab$pedu_sd[tab$status=="Survived"]),
    format_mean_sd(tab$pedu_mean[tab$status=="Died"],     tab$pedu_sd[tab$status=="Died"]),
    format_mean_sd(totals$pedu_mean, totals$pedu_sd)),
  c("Age at death, mean (SD)",
    "—",
    format_mean_sd(tab$age_death_mean[tab$status=="Died"], tab$age_death_sd[tab$status=="Died"]),
    format_mean_sd(totals$age_death_mean, totals$age_death_sd)),
  c("Follow-up time, mean (SD)",
    format_mean_sd(tab$followup_mean[tab$status=="Survived"], tab$followup_sd[tab$status=="Survived"]),
    format_mean_sd(tab$followup_mean[tab$status=="Died"],     tab$followup_sd[tab$status=="Died"]),
    format_mean_sd(totals$followup_mean, totals$followup_sd))
)

for (r in rows) {
  cat(sprintf("%-40s %-25s %-25s %-25s\n", r[1], r[2], r[3], r[4]))
}

cat("\nDeaths by cause (among those who died):\n")
cause_rows <- list(
  c("  Cancer, n (%)",
    "—",
    format_n_pct(causes$cancer_n[causes$status=="Died"],     causes$cancer_pct[causes$status=="Died"]),
    format_n_pct(causes$cancer_n[causes$status=="Died"],     causes$cancer_pct[causes$status=="Died"])),
  c("  Cardiovascular, n (%)",
    "—",
    format_n_pct(causes$cvd_n[causes$status=="Died"],        causes$cvd_pct[causes$status=="Died"]),
    format_n_pct(causes$cvd_n[causes$status=="Died"],        causes$cvd_pct[causes$status=="Died"])),
  c("  Respiratory, n (%)",
    "—",
    format_n_pct(causes$resp_n[causes$status=="Died"],       causes$resp_pct[causes$status=="Died"]),
    format_n_pct(causes$resp_n[causes$status=="Died"],       causes$resp_pct[causes$status=="Died"])),
  c("  External, n (%)",
    "—",
    format_n_pct(causes$external_n[causes$status=="Died"],   causes$external_pct[causes$status=="Died"]),
    format_n_pct(causes$external_n[causes$status=="Died"],   causes$external_pct[causes$status=="Died"])),
  c("    Accidental, n (%)",
    "—",
    format_n_pct(causes$accidental_n[causes$status=="Died"], causes$accidental_pct[causes$status=="Died"]),
    format_n_pct(causes$accidental_n[causes$status=="Died"], causes$accidental_pct[causes$status=="Died"])),
  c("    Suicide, n (%)",
    "—",
    format_n_pct(causes$suicide_n[causes$status=="Died"],    causes$suicide_pct[causes$status=="Died"]),
    format_n_pct(causes$suicide_n[causes$status=="Died"],    causes$suicide_pct[causes$status=="Died"])),
  c("  Other, n (%)",
    "—",
    format_n_pct(causes$other_n[causes$status=="Died"],      causes$other_pct[causes$status=="Died"]),
    format_n_pct(causes$other_n[causes$status=="Died"],      causes$other_pct[causes$status=="Died"]))
)

for (r in cause_rows) {
  cat(sprintf("%-40s %-25s %-25s %-25s\n", r[1], r[2], r[3], r[4]))
}
