### Measurement-error calibrated HR figure

library(tidyverse)
library(survival)
library(ggplot2)

outdir <- "/safe/data/KBW/"

# event variables (already made)
model_army_all <- model_army_all %>%
  mutate(
    ev_allcause = death,
    ev_external = as.integer(death == 1 & cause_group == "External"),
    ev_cancer = as.integer(death == 1 & cause_group == "Cancer"),
    ev_cvd = as.integer(death == 1 & cause_group == "Cardiovascular"),
    ev_resp = as.integer(death == 1 & cause_group == "Respiratory"),
    ev_other = as.integer(death == 1 & cause_group == "Other"),
    ev2_accidental = as.integer(death == 1 & cause_group2 == "Accidental"),
    ev2_suicide = as.integer(death == 1 & cause_group2 == "Suicide")
  )
# Matrix correction from Ole's Breit-anchored model
K <- matrix(
  c(1.085, 0,
    -0.626, 1.882),
  nrow = 2,
  byrow = TRUE
)

fit_and_correct <- function(data, event_var, outcome_name) {
  
  f <- as.formula(
    paste0(
      "Surv(start_age, stop_age, ", event_var, ") ~ ",
      "child_g_z + army_g_resid + z_parentedu_mean + strata(cohort)"
    )
  )
  
  fit <- coxph(f, data = data)
  
  b <- coef(fit)[c("child_g_z", "army_g_resid")]
  V <- vcov(fit)[c("child_g_z", "army_g_resid"),
                 c("child_g_z", "army_g_resid")]
  
  se <- sqrt(diag(V))
  
  raw <- tibble(
    outcome = outcome_name,
    estimate_type = "Observed-score",
    predictor = c("Childhood cognitive level",
                  "Adolescent cognitive change"),
    beta = as.numeric(b),
    se = as.numeric(se),
    HR = exp(as.numeric(b)),
    CI_low = exp(as.numeric(b) - 1.96 * as.numeric(se)),
    CI_high = exp(as.numeric(b) + 1.96 * as.numeric(se)),
    events = sum(data[[event_var]], na.rm = TRUE)
  )
  
  b_corr <- as.vector(K %*% b)
  V_corr <- K %*% V %*% t(K)
  se_corr <- sqrt(diag(V_corr))
  
  corrected <- tibble(
    outcome = outcome_name,
    estimate_type = "Matrix-corrected",
    predictor = c("Childhood cognitive level",
                  "Adolescent cognitive change"),
    beta = b_corr,
    se = se_corr,
    HR = exp(b_corr),
    CI_low = exp(b_corr - 1.96 * se_corr),
    CI_high = exp(b_corr + 1.96 * se_corr),
    events = sum(data[[event_var]], na.rm = TRUE)
  )
  
  bind_rows(raw, corrected)
}

me_results <- bind_rows(
  fit_and_correct(model_army_all, "ev_allcause", "All-cause"),
  fit_and_correct(model_army_all, "ev_cancer", "Cancer"),
  fit_and_correct(model_army_all, "ev_cvd", "Cardiovascular"),
  fit_and_correct(model_army_all, "ev_resp", "Respiratory"),
  fit_and_correct(model_army_all, "ev_other", "Other"),
  fit_and_correct(model_army_all, "ev_external", "External"),
  fit_and_correct(model_army_all, "ev2_accidental", "Accidental"),
  fit_and_correct(model_army_all, "ev2_suicide", "Suicide")
)

me_table <- me_results %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f–%.2f)", HR, CI_low, CI_high)
  ) %>%
  select(outcome, events, estimate_type, predictor, HR_CI)

write.table(
  me_table,
  file = file.path(outdir, "measurement_error_calibrated_HRs.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

print(me_table)

plotdat <- me_results %>%
  mutate(
    outcome_lab = paste0(outcome, " (n = ", events, ")"),
    outcome_lab = factor(
      outcome_lab,
      levels = rev(c(
        paste0("All-cause (n = ", sum(model_army_all$ev_allcause, na.rm = TRUE), ")"),
        paste0("Cancer (n = ", sum(model_army_all$ev_cancer, na.rm = TRUE), ")"),
        paste0("Cardiovascular (n = ", sum(model_army_all$ev_cvd, na.rm = TRUE), ")"),
        paste0("Respiratory (n = ", sum(model_army_all$ev_resp, na.rm = TRUE), ")"),
        paste0("Other (n = ", sum(model_army_all$ev_other, na.rm = TRUE), ")"),
        paste0("External (n = ", sum(model_army_all$ev_external, na.rm = TRUE), ")"),
        paste0("Accidental (n = ", sum(model_army_all$ev2_accidental, na.rm = TRUE), ")"),
        paste0("Suicide (n = ", sum(model_army_all$ev2_suicide, na.rm = TRUE), ")")
      ))
    ),
    estimate_type = factor(
      estimate_type,
      levels = c("Observed-score", "Matrix-corrected")
    ),
    predictor = factor(
      predictor,
      levels = c("Childhood cognitive level",
                 "Adolescent cognitive change")
    )
  )

p_me <- ggplot(
  plotdat,
  aes(
    x = HR,
    y = outcome_lab,
    xmin = CI_low,
    xmax = CI_high,
    color = predictor
  )
) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey25") +
  geom_pointrange(
    position = position_dodge(width = 0.45),
    linewidth = 0.65,
    size = 0.55
  ) +
  facet_wrap(~ estimate_type, nrow = 1) +
  scale_x_log10(
    limits = c(0.35, 1.35),
    breaks = c(0.5, 0.7, 1.0),
    labels = c("0.5", "0.7", "1.0")
  ) +
  scale_color_manual(
    values = c(
      "Childhood cognitive level" = "#00BFC4",
      "Adolescent cognitive change" = "#F8766D"
    )
  ) +
  labs(
    title = "Observed-score and measurement-error calibrated mortality associations",
    x = "Hazard ratio (log scale)",
    y = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0),
    panel.spacing = unit(1.2, "lines")
  )

print(p_me)

ggsave(
  filename = file.path(outdir, "measurement_error_calibrated_HRs.png"),
  plot = p_me,
  width = 9,
  height = 5.5,
  dpi = 300
)