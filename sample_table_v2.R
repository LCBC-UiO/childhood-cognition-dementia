library(tidyverse)

# ── Helper functions ───────────────────────────────────────────────────────────
fmt_mean_sd <- function(x) sprintf("%.2f (%.2f)", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
fmt_n_pct   <- function(x, n_total) sprintf("%d (%.1f%%)", sum(x, na.rm=TRUE), sum(x, na.rm=TRUE)/n_total*100)

# ── Split data ─────────────────────────────────────────────────────────────────
survived <- model_army_all %>% filter(death == 0)
died     <- model_army_all %>% filter(death == 1)
total    <- model_army_all

n_surv <- nrow(survived)
n_died <- nrow(died)
n_tot  <- nrow(total)

cat(sprintf("\n%-45s %-22s %-22s %-22s\n",
    "Characteristic",
    paste0("Survived (n=", n_surv, ")"),
    paste0("Died (n=", n_died, ")"),
    paste0("Total (n=", n_tot, ")")))
cat(strrep("-", 111), "\n")

# ── Cohort ─────────────────────────────────────────────────────────────────────
cat(sprintf("%-45s %-22s %-22s %-22s\n", "Cohort 1948, n (%)",
    fmt_n_pct(survived$cohort == 1948, n_surv),
    fmt_n_pct(died$cohort == 1948,     n_died),
    fmt_n_pct(total$cohort == 1948,    n_tot)))

cat(sprintf("%-45s %-22s %-22s %-22s\n", "Cohort 1953, n (%)",
    fmt_n_pct(survived$cohort == 1953, n_surv),
    fmt_n_pct(died$cohort == 1953,     n_died),
    fmt_n_pct(total$cohort == 1953,    n_tot)))

# ── Cognitive measures ─────────────────────────────────────────────────────────
cat(sprintf("%-45s %-22s %-22s %-22s\n", "Childhood cognitive level, mean (SD)",
    fmt_mean_sd(survived$child_g_z),
    fmt_mean_sd(died$child_g_z),
    fmt_mean_sd(total$child_g_z)))

cat(sprintf("%-45s %-22s %-22s %-22s\n", "Army g (age 18), mean (SD)",
    fmt_mean_sd(survived$army_g_z),
    fmt_mean_sd(died$army_g_z),
    fmt_mean_sd(total$army_g_z)))

cat(sprintf("%-45s %-22s %-22s %-22s\n", "Adolescent cognitive change, mean (SD)",
    fmt_mean_sd(survived$army_g_resid),
    fmt_mean_sd(died$army_g_resid),
    fmt_mean_sd(total$army_g_resid)))

# ── Parental education ─────────────────────────────────────────────────────────
cat(sprintf("%-45s %-22s %-22s %-22s\n", "Parental education (z-score), mean (SD)",
    fmt_mean_sd(survived$parentedu_mean_z),
    fmt_mean_sd(died$parentedu_mean_z),
    fmt_mean_sd(total$parentedu_mean_z)))

# ── Age / follow-up ────────────────────────────────────────────────────────────
cat(sprintf("%-45s %-22s %-22s %-22s\n", "Age at death, mean (SD)",
    "—",
    fmt_mean_sd(died$stop_age),
    fmt_mean_sd(total$stop_age[total$death == 1])))

cat(sprintf("%-45s %-22s %-22s %-22s\n", "Follow-up time (years), mean (SD)",
    fmt_mean_sd(survived$stop_age - survived$start_age),
    fmt_mean_sd(died$stop_age - died$start_age),
    fmt_mean_sd(total$stop_age - total$start_age)))

# ── Deaths by cause (% of full sample) ────────────────────────────────────────
cat(strrep("-", 111), "\n")
cat(sprintf("%-45s %-22s %-22s %-22s\n", "Deaths by cause, n (% of full sample)", "", "", ""))

causes <- list(
  "Cancer"        = "ev_cancer",
  "Cardiovascular"= "ev_cvd",
  "Respiratory"   = "ev_resp",
  "External"      = "ev_external",
  "  Accidental"  = "ev2_accidental",
  "  Suicide"     = "ev_suicide",
  "Other"         = "ev_other"
)

for (label in names(causes)) {
  var <- causes[[label]]
  if (var %in% names(total)) {
    cat(sprintf("%-45s %-22s %-22s %-22s\n", label,
        "—",
        fmt_n_pct(died[[var]],  n_tot),
        fmt_n_pct(total[[var]], n_tot)))
  } else {
    cat(sprintf("  [Variable %s not found in data]\n", var))
  }
}

cat(strrep("-", 111), "\n")
cat("Note: Cognitive measures and parental education are standardized (z-scores) within cohort.\n")
cat("Adolescent cognitive change = residual of army g after regressing on childhood g.\n")
cat("Death cause percentages are of the full analytic sample (n =", n_tot, ").\n")
