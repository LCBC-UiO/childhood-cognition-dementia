############################################################
## Intrinsic / external / accidental mortality
## + formal heterogeneity tests
##
## Uses existing variables in model_army_all:
## start_age
## stop_age
## cohort
## child_g_z
## army_g_resid
## z_parentedu_mean
## ev_intrinsic
## ev_external
## ev2_accidental
############################################################

library(survival)

dat <- model_army_all
dat$.id_stack <- seq_len(nrow(dat))

############################################################
## 1. Standard cause-specific Cox models
############################################################

m_intrinsic <- coxph(
  Surv(start_age, stop_age, ev_intrinsic) ~
    child_g_z + army_g_resid + z_parentedu_mean + strata(cohort),
  data = dat
)

m_external <- coxph(
  Surv(start_age, stop_age, ev_external) ~
    child_g_z + army_g_resid + z_parentedu_mean + strata(cohort),
  data = dat
)

m_accidental <- coxph(
  Surv(start_age, stop_age, ev2_accidental) ~
    child_g_z + army_g_resid + z_parentedu_mean + strata(cohort),
  data = dat
)

cat("\n\n### INTRINSIC ###\n")
print(summary(m_intrinsic))

cat("\n\n### EXTERNAL ###\n")
print(summary(m_external))

cat("\n\n### ACCIDENTAL ###\n")
print(summary(m_accidental))


############################################################
## 2. Word-ready HR table
############################################################

get_hr <- function(model, outcome) {
  s <- summary(model)
  
  x <- data.frame(
    Outcome = outcome,
    Predictor = rownames(s$coef),
    HR = s$conf.int[, "exp(coef)"],
    L95 = s$conf.int[, "lower .95"],
    U95 = s$conf.int[, "upper .95"],
    p = s$coef[, "Pr(>|z|)"],
    row.names = NULL
  )
  
  x <- x[x$Predictor %in% c("child_g_z", "army_g_resid"), ]
  
  x$Predictor <- ifelse(
    x$Predictor == "child_g_z",
    "Childhood cognitive level",
    "Adolescent cognitive change"
  )
  
  x$HR_95_CI <- sprintf("%.2f (%.2f–%.2f)", x$HR, x$L95, x$U95)
  
  x$p_value <- ifelse(
    x$p < .001,
    "<0.001",
    sprintf("%.3f", x$p)
  )
  
  x[, c("Outcome", "Predictor", "HR_95_CI", "p_value")]
}

tab_main <- rbind(
  get_hr(m_intrinsic, "Intrinsic"),
  get_hr(m_external, "External"),
  get_hr(m_accidental, "Accidental")
)

cat("\n\n### TABLE FOR MANUSCRIPT ###\n")
print(tab_main, row.names = FALSE)

cat("\n\n### TAB-SEPARATED VERSION ###\n")
write.table(
  tab_main,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


############################################################
## 3. Function for stacked heterogeneity tests
############################################################

make_stack <- function(dat, event_a, event_b, name_a, name_b) {
  
  d_a <- dat
  d_a$cause_group <- name_a
  d_a$event_stack <- d_a[[event_a]]
  
  d_b <- dat
  d_b$cause_group <- name_b
  d_b$event_stack <- d_b[[event_b]]
  
  out <- rbind(d_a, d_b)
  out$cause_group <- factor(out$cause_group, levels = c(name_a, name_b))
  
  out
}

test_interactions <- function(model) {
  
  cn <- names(coef(model))
  
  child_int <- cn[grepl("child_g_z:cause_group|cause_group.*:child_g_z", cn)]
  resid_int <- cn[grepl("army_g_resid:cause_group|cause_group.*:army_g_resid", cn)]
  
  run_test <- function(terms) {
    b <- coef(model)[terms]
    V <- vcov(model)[terms, terms, drop = FALSE]
    chi <- as.numeric(t(b) %*% solve(V) %*% b)
    df <- length(terms)
    p <- pchisq(chi, df = df, lower.tail = FALSE)
    
    data.frame(
      Chi_square = chi,
      df = df,
      p = p
    )
  }
  
  cat("\nChildhood cognitive level difference across causes:\n")
  print(run_test(child_int))
  
  cat("\nAdolescent cognitive change difference across causes:\n")
  print(run_test(resid_int))
  
  cat("\nJoint test of both cognitive differences across causes:\n")
  print(run_test(c(child_int, resid_int)))
}


############################################################
## 4. Intrinsic vs external heterogeneity test
############################################################

stack_intr_ext <- make_stack(
  dat = dat,
  event_a = "ev_intrinsic",
  event_b = "ev_external",
  name_a = "Intrinsic",
  name_b = "External"
)

m_intr_ext_het <- coxph(
  Surv(start_age, stop_age, event_stack) ~
    child_g_z * cause_group +
    army_g_resid * cause_group +
    z_parentedu_mean +
    strata(cohort, cause_group) +
    cluster(.id_stack),
  data = stack_intr_ext
)

cat("\n\n### HETEROGENEITY: INTRINSIC VS EXTERNAL ###\n")
print(summary(m_intr_ext_het))
test_interactions(m_intr_ext_het)


############################################################
## 5. Intrinsic vs accidental heterogeneity test
############################################################

stack_intr_acc <- make_stack(
  dat = dat,
  event_a = "ev_intrinsic",
  event_b = "ev2_accidental",
  name_a = "Intrinsic",
  name_b = "Accidental"
)

m_intr_acc_het <- coxph(
  Surv(start_age, stop_age, event_stack) ~
    child_g_z * cause_group +
    army_g_resid * cause_group +
    z_parentedu_mean +
    strata(cohort, cause_group) +
    cluster(.id_stack),
  data = stack_intr_acc
)

cat("\n\n### HETEROGENEITY: INTRINSIC VS ACCIDENTAL ###\n")
print(summary(m_intr_acc_het))
test_interactions(m_intr_acc_het)


############################################################
## 6. Interpretation text helper
############################################################

cat("\n\n### HOW TO REPORT ###\n")

cat("\nIf heterogeneity tests are NOT significant:\n")
cat("Associations showed qualitatively different tendencies across causes of death, ")
cat("but formal heterogeneity tests did not provide strong evidence that the coefficients ")
cat("differed statistically across cause groups.\n")

cat("\nIf intrinsic vs external heterogeneity IS significant:\n")
cat("Formal heterogeneity tests supported differences between intrinsic and external mortality, ")
cat("with childhood cognitive level relatively more strongly associated with intrinsic mortality ")
cat("and adolescent cognitive change relatively more strongly associated with external mortality.\n")

cat("\n\nDone.\n")