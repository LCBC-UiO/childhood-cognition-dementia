# =============================================================================
# Cox_prep_13_diagnostics_April2026.R
# Supplementary diagnostic analyses for Walhovd et al. dementia paper
# Run AFTER Cox_prep_12_UGU_paper1_April9th2026.R has completed successfully.
# All objects created in that script (cox_ugu1948, model_data48,
# cox_ugu1948_tv, cox_all, cox_full, etc.) are assumed to be in the environment.
#
# New analyses added here (not in Cox_prep_12):
#   1. VIF / multicollinearity diagnostics for cognitive predictors
#   2. Individual subtest models (each domain entered alone — addresses R1)
#   3. Schoenfeld residual test of proportional hazards assumption
#   4. Birth month sensitivity analysis
#   5. Emigration sensitivity analysis (template — requires emigration data)
#
# All sections print output to the console. Save outputs to /safe/data/KBW/
# by redirecting with sink() or capture.output() as needed.
# =============================================================================

# Ensure all required packages are loaded
library(survival)
library(car)      # for vif() on coxph objects
library(tidyverse)


# =============================================================================
# SECTION 1: VIF / MULTICOLLINEARITY DIAGNOSTICS
# =============================================================================
# Why: Reviewer 1 raised the "perils of partialling" concern — entering
# correlated cognitive predictors jointly inflates standard errors and may
# make unique-variance estimates unreliable. VIF (Variance Inflation Factor)
# quantifies how much each predictor's variance is inflated by collinearity
# with the others. Rule of thumb: VIF < 5 is acceptable; VIF > 10 is severe.
# With inter-correlations of .43–.56 we expect VIFs around 1.3–1.6.
# =============================================================================

cat("\n=== SECTION 1: VIF / MULTICOLLINEARITY DIAGNOSTICS ===\n")

# Use the main Cox model (cox_full or cox_all — both have the same cognitive
# predictors; cox_full was built on model_data48 which is the cleaner
# filtered dataset, so use that for the VIF calculation).
# car::vif() handles coxph objects via a GVIF (generalised VIF) calculation.

cat("\n--- VIF for cognitive predictors in joint dementia model (car::vif) ---\n")
cat("Using cox_full (3 cognitive tests + parental education + sex)\n\n")

vif_result <- vif(cox_full)
print(vif_result)

# INTERPRETATION GUIDE (print to output):
cat("\nInterpretation guide:
  VIF = 1.0: no collinearity with other predictors
  VIF 1–2  : low collinearity, standard errors modestly inflated
  VIF 2–5  : moderate collinearity, worth noting
  VIF > 5  : high collinearity, potentially problematic
  VIF > 10 : severe collinearity\n")

# --- Manual approach (backup if car is not available in vault) ---
# Calculate VIF from the predictor correlation matrix.
# For a predictor X_j, VIF_j = 1 / (1 - R^2_j), where R^2_j is the
# R-squared from regressing X_j on all other predictors.

cat("\n--- Manual VIF calculation from predictor correlations (backup method) ---\n")

# Use the same model_data48 as in the main Cox models
pred_matrix <- model_data48 %>%
  select(z_TS6IITP, z_TS6ISTP, z_TS6IVOTP, z_parentedu_mean, sex01) %>%
  na.omit()

# Inductive reasoning (z_TS6IITP) regressed on all other predictors
r2_inductive <- summary(lm(z_TS6IITP ~ z_TS6ISTP + z_TS6IVOTP +
                              z_parentedu_mean + sex01,
                            data = pred_matrix))$r.squared
vif_inductive <- 1 / (1 - r2_inductive)

# Spatial ability (z_TS6ISTP)
r2_spatial <- summary(lm(z_TS6ISTP ~ z_TS6IITP + z_TS6IVOTP +
                            z_parentedu_mean + sex01,
                          data = pred_matrix))$r.squared
vif_spatial <- 1 / (1 - r2_spatial)

# Verbal ability (z_TS6IVOTP)
r2_verbal <- summary(lm(z_TS6IVOTP ~ z_TS6IITP + z_TS6ISTP +
                           z_parentedu_mean + sex01,
                         data = pred_matrix))$r.squared
vif_verbal <- 1 / (1 - r2_verbal)

# Parental education
r2_parentedu <- summary(lm(z_parentedu_mean ~ z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
                              sex01,
                            data = pred_matrix))$r.squared
vif_parentedu <- 1 / (1 - r2_parentedu)

vif_manual <- data.frame(
  Predictor = c("Inductive reasoning (z_TS6IITP)",
                "Spatial ability (z_TS6ISTP)",
                "Verbal ability (z_TS6IVOTP)",
                "Parental education (z_parentedu_mean)"),
  R2_with_others = round(c(r2_inductive, r2_spatial, r2_verbal, r2_parentedu), 3),
  VIF = round(c(vif_inductive, vif_spatial, vif_verbal, vif_parentedu), 3)
)

print(vif_manual)

# Also print the condition number of the predictor correlation matrix.
# Condition numbers > 30 indicate problematic collinearity.
cond_num <- kappa(cor(pred_matrix))
cat(sprintf("\nCondition number of predictor correlation matrix: %.2f\n", cond_num))
cat("(Values < 30 are generally acceptable)\n")


# =============================================================================
# SECTION 2: INDIVIDUAL SUBTEST MODELS (each domain entered separately)
# =============================================================================
# Why: When predictors are correlated, the joint model partials out shared
# variance. An alternative is to enter each domain alone — this shows the
# total (not unique) association of each domain with dementia.
# These go in the Supplement as sensitivity analyses addressing R1.
# Expected finding: inductive reasoning should also predict dementia when
# entered alone; verbal and spatial should be weaker or null.
# =============================================================================

cat("\n=== SECTION 2: INDIVIDUAL SUBTEST MODELS (each domain alone) ===\n")
cat("These are Supplement eTables addressing the 'perils of partialling' concern.\n")
cat("Model: Surv(start_age, end_age, dementia_any) ~ one cognitive test + parental edu + sex\n\n")

# --- 2a. Inductive reasoning alone ---
cat("--- 2a: Inductive reasoning alone ---\n")
cox_inductive_only <- coxph(
  Surv(start_age, end_age, dementia_any) ~
    z_TS6IITP + z_parentedu_mean + sex01,
  data = model_data48
)
summary(cox_inductive_only)

# --- 2b. Spatial ability alone ---
cat("\n--- 2b: Spatial ability alone ---\n")
cox_spatial_only <- coxph(
  Surv(start_age, end_age, dementia_any) ~
    z_TS6ISTP + z_parentedu_mean + sex01,
  data = model_data48
)
summary(cox_spatial_only)

# --- 2c. Verbal ability alone ---
cat("\n--- 2c: Verbal ability alone ---\n")
cox_verbal_only <- coxph(
  Surv(start_age, end_age, dementia_any) ~
    z_TS6IVOTP + z_parentedu_mean + sex01,
  data = model_data48
)
summary(cox_verbal_only)

# --- 2d. Comparison summary table ---
cat("\n--- 2d: Summary — HRs for inductive reasoning across model specifications ---\n")
cat("This table is useful for the Supplement to show robustness.\n\n")

# Extract HR and CI for inductive reasoning (z_TS6IITP) from each model
extract_hr <- function(model, predictor, label) {
  coefs <- summary(model)$coefficients
  conf  <- exp(confint.default(model))
  if (predictor %in% rownames(coefs)) {
    hr  <- round(exp(coefs[predictor, "coef"]), 3)
    lo  <- round(conf[predictor, 1], 3)
    hi  <- round(conf[predictor, 2], 3)
    p   <- round(coefs[predictor, "Pr(>|z|)"], 4)
    cat(sprintf("%-50s HR = %.3f (%.3f–%.3f), p = %.4f\n", label, hr, lo, hi, p))
  } else {
    cat(sprintf("%-50s [predictor not in model]\n", label))
  }
}

extract_hr(cox_inductive_only, "z_TS6IITP",
           "Inductive alone (no partialling)")
extract_hr(cox_full,           "z_TS6IITP",
           "Inductive joint (all 3 tests, partialled)")
extract_hr(cox_all_edu,        "z_TS6IITP",
           "Inductive joint + midlife education")

# Same for spatial and verbal (to show they truly add nothing):
cat("\n")
extract_hr(cox_spatial_only, "z_TS6ISTP",
           "Spatial alone")
extract_hr(cox_full,         "z_TS6ISTP",
           "Spatial joint (partialled)")
cat("\n")
extract_hr(cox_verbal_only,  "z_TS6IVOTP",
           "Verbal alone")
extract_hr(cox_full,         "z_TS6IVOTP",
           "Verbal joint (partialled)")


# =============================================================================
# SECTION 3: PROPORTIONAL HAZARDS ASSUMPTION CHECK (Schoenfeld residuals)
# =============================================================================
# Why: The Cox model assumes that the hazard ratio is constant over time
# (proportional hazards). If this is violated, the HR estimate averages
# over a time-varying effect and the p-values may be misleading.
# The cox.zph() test regresses scaled Schoenfeld residuals on time.
# A significant p-value indicates a violation for that predictor.
# This check is good practice and should be reported in the Supplement.
# =============================================================================

cat("\n=== SECTION 3: PROPORTIONAL HAZARDS ASSUMPTION (Schoenfeld residuals) ===\n")
cat("H0: hazard ratio is constant over time (proportional hazards holds)\n")
cat("Significant p < .05 suggests a time-varying hazard ratio for that predictor.\n\n")

ph_test <- cox.zph(cox_full)
print(ph_test)

cat("\n")
# Plot the smooth of Schoenfeld residuals over time for visual inspection.
# Save to file for potential inclusion in supplement.
tiff("/safe/data/KBW/PH_assumption_Schoenfeld_residuals.tiff",
     width = 10, height = 8, units = "in", res = 300)
par(mfrow = c(2, 3))  # layout for up to 5 predictors + global test
plot(ph_test)
dev.off()
cat("Schoenfeld residual plots saved to /safe/data/KBW/PH_assumption_Schoenfeld_residuals.tiff\n")

# Note: if any cognitive predictor shows a significant violation,
# consider adding a time*predictor interaction term, or stratifying by
# a time-based variable. For this paper, a brief Methods note is sufficient
# if the violation is minor.


# =============================================================================
# SECTION 4: BIRTH MONTH SENSITIVITY ANALYSIS
# =============================================================================
# Why: Within the 1948 cohort, all children were tested in May 1961,
# but birth months varied. Younger children (born later in 1948) were
# slightly younger at testing and may score marginally lower.
# This is the "relative age effect" concern raised by Martin Fischer.
# If birth month is unrelated to dementia after adjusting for test score,
# it is not a confounder and need not be included in main models.
# If it does attenuate the cognition-dementia association, it is worth noting.
#
# RSMONTH is already in the data (used to construct birth_date).
# We simply add it as an additional covariate.
# =============================================================================

cat("\n=== SECTION 4: BIRTH MONTH SENSITIVITY ANALYSIS ===\n")
cat("Testing whether birth month (relative age within cohort) confounds results.\n")
cat("RSMONTH: 1 = January, 12 = December (all born in 1948, tested May 1961).\n")
cat("Younger children at test (RSMONTH = 4–12, born closer to test date) may score lower.\n\n")

# Check distribution of birth months
cat("Birth month distribution in analytic sample:\n")
print(table(model_data48$RSMONTH, useNA = "ifany"))

# Main dementia model + birth month as covariate
# We include RSMONTH as a continuous covariate (months 1–12)
# standardise it so the effect is on the same scale as other predictors

# Add RSMONTH to model_data48 if not already present
model_data48_bm <- cox_ugu1948 %>%
  select(ID, dementia_any, sex01, z_parentedu_mean,
         z_TS6IITP, z_TS6ISTP, z_TS6IVOTP,
         start_age, end_age, RSMONTH) %>%
  filter(!is.na(start_age)) %>%
  filter(!is.na(end_age)) %>%
  filter(!is.na(dementia_any)) %>%
  filter(!is.na(sex01)) %>%
  filter(!is.na(z_parentedu_mean)) %>%
  filter(!is.na(z_TS6IITP)) %>%
  filter(!is.na(z_TS6ISTP)) %>%
  filter(!is.na(z_TS6IVOTP)) %>%
  filter(!is.na(RSMONTH)) %>%
  mutate(z_RSMONTH = as.numeric(scale(RSMONTH)))

cat(sprintf("\nN in birth month model: %d (events: %d)\n",
            nrow(model_data48_bm),
            sum(model_data48_bm$dementia_any)))

cat("\n--- Main model without birth month (reference) ---\n")
cox_bm_ref <- coxph(
  Surv(start_age, end_age, dementia_any) ~
    z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
    z_parentedu_mean + sex01,
  data = model_data48_bm   # same N as birth month model for fair comparison
)
summary(cox_bm_ref)

cat("\n--- Model with birth month added ---\n")
cox_bm_sensitivity <- coxph(
  Surv(start_age, end_age, dementia_any) ~
    z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
    z_parentedu_mean + sex01 + z_RSMONTH,
  data = model_data48_bm
)
summary(cox_bm_sensitivity)

# Key thing to check: does the HR for z_TS6IITP change materially?
cat("\n--- Change in inductive reasoning HR after birth month adjustment ---\n")
hr_without_bm <- exp(coef(cox_bm_ref)["z_TS6IITP"])
hr_with_bm    <- exp(coef(cox_bm_sensitivity)["z_TS6IITP"])
cat(sprintf("HR z_TS6IITP without birth month: %.3f\n", hr_without_bm))
cat(sprintf("HR z_TS6IITP with birth month:    %.3f\n", hr_with_bm))
cat(sprintf("Attenuation: %.1f%%\n",
            100 * (hr_without_bm - hr_with_bm) / (hr_without_bm - 1)))
# Interpretation: if attenuation is < 5–10%, birth month is not a meaningful confounder.


# =============================================================================
# SECTION 5: EMIGRATION SENSITIVITY ANALYSIS (template)
# =============================================================================
# Why: Individuals who emigrated from Sweden after 1961 are administratively
# censored in NPR-IPR and CDR when they leave. If emigration is correlated
# with cognitive ability (e.g., higher-ability individuals more likely to
# emigrate), this creates informative censoring — people who emigrate are
# not at risk of being registered with dementia in Sweden, but may develop
# dementia abroad. This would bias associations toward the null if
# higher-ability individuals emigrate more.
#
# Standard approach: obtain emigration dates from Statistics Sweden's
# emigration register (EMIG), censor at emigration date instead of
# end of follow-up, and rerun primary models. Compare HRs.
#
# NOTE: This template assumes the emigration data is available as a file
# called emig_data.csv in /safe/data/RAW_DATA/SCB/EMIG/ (adjust path as needed).
# Variable names assumed: ID, emigration_date (YYYY-MM-DD format).
# =============================================================================

cat("\n=== SECTION 5: EMIGRATION SENSITIVITY ANALYSIS (template) ===\n")
cat("*** Requires emigration data from Statistics Sweden. ***\n")
cat("*** Adjust file path and variable names before running. ***\n\n")

# --- Step 1: Load emigration data ---
# emig_raw <- read.csv("/safe/data/RAW_DATA/SCB/EMIG/emig_data.csv")
#
# # Keep only first emigration per person (some may re-migrate)
# emig_first <- emig_raw %>%
#   mutate(emigration_date = as.Date(emigration_date)) %>%
#   group_by(ID) %>%
#   slice_min(emigration_date, with_ties = FALSE) %>%
#   ungroup() %>%
#   select(ID, emigration_date)
#
# cat(sprintf("Emigration events found in cohort: %d\n", nrow(emig_first)))

# --- Step 2: Add emigration date to 1948 cohort ---
# cox_ugu1948_emig <- cox_ugu1948 %>%
#   left_join(emig_first, by = "ID")
#
# # Redefine censor date: earliest of (death, emigration, end_of_follow_up)
# cox_ugu1948_emig <- cox_ugu1948_emig %>%
#   mutate(
#     censor_date_emig = pmin(
#       CDR_death_date,
#       emigration_date,
#       end_of_follow_up,
#       na.rm = TRUE
#     ),
#     event_date_emig = if_else(dementia_any == 1L, dementia_date, censor_date_emig),
#     end_age_emig    = as.numeric((event_date_emig - birth_date) / 365.25)
#   )
#
# # Exclude anyone whose dementia date is AFTER the emigration date
# # (they would have been registered abroad, so we cannot count them)
# cox_ugu1948_emig <- cox_ugu1948_emig %>%
#   mutate(
#     dementia_any_emig = if_else(
#       dementia_any == 1L & !is.na(emigration_date) &
#         dementia_date > emigration_date,
#       0L,          # treat as censored at emigration
#       dementia_any
#     )
#   )
#
# cat(sprintf("N after emigration censoring: %d\n",
#             sum(!is.na(cox_ugu1948_emig$end_age_emig))))
# cat(sprintf("Dementia events after emigration censoring: %d\n",
#             sum(cox_ugu1948_emig$dementia_any_emig == 1L, na.rm = TRUE)))

# --- Step 3: Rebuild model_data with emigration censoring ---
# model_data_emig <- cox_ugu1948_emig %>%
#   select(ID, dementia_any_emig, sex01, z_parentedu_mean,
#          z_TS6IITP, z_TS6ISTP, z_TS6IVOTP,
#          start_age, end_age_emig) %>%
#   rename(dementia_any = dementia_any_emig,
#          end_age      = end_age_emig) %>%
#   filter(!is.na(start_age)) %>%
#   filter(!is.na(end_age)) %>%
#   filter(!is.na(dementia_any)) %>%
#   filter(!is.na(sex01)) %>%
#   filter(!is.na(z_parentedu_mean)) %>%
#   filter(!is.na(z_TS6IITP)) %>%
#   filter(!is.na(z_TS6ISTP)) %>%
#   filter(!is.na(z_TS6IVOTP))

# --- Step 4: Rerun primary dementia model with emigration censoring ---
# cox_emig_sensitivity <- coxph(
#   Surv(start_age, end_age, dementia_any) ~
#     z_TS6IITP + z_TS6ISTP + z_TS6IVOTP +
#     z_parentedu_mean + sex01,
#   data = model_data_emig
# )
# summary(cox_emig_sensitivity)
#
# # Compare HR for inductive reasoning
# hr_no_emig  <- exp(coef(cox_full)["z_TS6IITP"])
# hr_emig     <- exp(coef(cox_emig_sensitivity)["z_TS6IITP"])
# cat(sprintf("\nHR z_TS6IITP primary analysis:              %.3f\n", hr_no_emig))
# cat(sprintf("HR z_TS6IITP censoring at emigration:       %.3f\n", hr_emig))
# cat("If similar, emigration is not materially biasing results.\n")


# =============================================================================
# SECTION 6: ADDITIONAL SENSITIVITY — INDUCTIVE REASONING IN TMERGE MODELS
# =============================================================================
# Check that the inductive reasoning HR is stable across all tmerge models.
# This is already implicit from the main results but a summary table is
# useful for the Supplement.
# =============================================================================

cat("\n=== SECTION 6: INDUCTIVE REASONING HR ACROSS ALL TMERGE MODELS ===\n")
cat("Confirms robustness of inductive reasoning effect after time-varying adjustment.\n\n")

tmerge_models <- list(
  "Primary (no morbidity covariate)"      = cox_full,
  "CCI as time-varying covariate"         = cox_dementia_any_somatic,
  "Composite CVD as time-varying"         = cox_dementia_any_CVD,
  "Cerebrovascular as time-varying"       = cox_dementia_CVD,
  "Myocardial infarction as time-varying" = cox_dementia_myocardial,
  "Heart failure as time-varying"         = cox_dementia_heartfailure,
  "Diabetes as time-varying"              = cox_dementia_diabetes,
  "CVD + midlife education"               = cox_dementia_any_CVD_edu,
  "CCI + midlife education"               = cox_dementia_any_somatic_edu
)

cat(sprintf("%-47s %6s  %16s  %7s\n",
            "Model", "HR", "95% CI", "p"))
cat(strrep("-", 80), "\n")

for (label in names(tmerge_models)) {
  m <- tmerge_models[[label]]
  coefs <- tryCatch(summary(m)$coefficients, error = function(e) NULL)
  conf  <- tryCatch(exp(confint.default(m)), error = function(e) NULL)

  if (!is.null(coefs) && "z_TS6IITP" %in% rownames(coefs)) {
    hr  <- round(exp(coefs["z_TS6IITP", "coef"]), 3)
    lo  <- round(conf["z_TS6IITP", 1], 3)
    hi  <- round(conf["z_TS6IITP", 2], 3)
    p   <- round(coefs["z_TS6IITP", "Pr(>|z|)"], 4)
    cat(sprintf("%-47s %6.3f  [%5.3f – %5.3f]  %7.4f\n", label, hr, lo, hi, p))
  } else {
    cat(sprintf("%-47s  [not available]\n", label))
  }
}

cat("\nNote: 'Primary' model uses Surv(start_age, end_age, dementia_any).\n")
cat("All other models use Surv(tstart, tstop, dementia_tv) from cox_ugu1948_tv.\n")
cat("Small differences in HR across models indicate stable inductive reasoning association.\n")


# =============================================================================
# DONE
# =============================================================================
cat("\n=== DIAGNOSTICS COMPLETE ===\n")
cat("Objects created in this script:\n")
cat("  vif_result           — car::vif() output on cox_full\n")
cat("  vif_manual           — manual VIF data frame\n")
cat("  cox_inductive_only   — inductive reasoning alone (dementia)\n")
cat("  cox_spatial_only     — spatial ability alone (dementia)\n")
cat("  cox_verbal_only      — verbal ability alone (dementia)\n")
cat("  ph_test              — cox.zph() proportional hazards test\n")
cat("  cox_bm_ref           — primary model on birth-month N (for comparison)\n")
cat("  cox_bm_sensitivity   — primary model + RSMONTH as covariate\n")
cat("  (emigration section is commented template — activate when data available)\n\n")
cat("Saved files:\n")
cat("  /safe/data/KBW/PH_assumption_Schoenfeld_residuals.tiff\n")
