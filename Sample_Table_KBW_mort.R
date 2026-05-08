############################################################
### TABLE 1: SAMPLE CHARACTERISTICS BY ALIVE / DEAD STATUS
### mean (SD), no age interactions
############################################################

library(dplyr)
library(tidyr)

fmt_mean_sd <- function(x, digits = 2) {
  if (all(is.na(x))) return("—")
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  paste0(round(m, digits), " (", round(s, digits), ")")
}

fmt_n_pct <- function(x, level) {
  n <- sum(x == level, na.rm = TRUE)
  d <- sum(!is.na(x))
  paste0(n, " (", round(100 * n / d, 1), "%)")
}

make_row_cont <- function(dat, var, label, digits = 2) {
  dat %>%
    group_by(Mortality_status) %>%
    summarise(value = fmt_mean_sd(.data[[var]], digits), .groups = "drop") %>%
    pivot_wider(names_from = Mortality_status, values_from = value) %>%
    mutate(Characteristic = label, .before = 1)
}

make_row_cat <- function(dat, var, level, label) {
  dat %>%
    group_by(Mortality_status) %>%
    summarise(value = fmt_n_pct(.data[[var]], level), .groups = "drop") %>%
    pivot_wider(names_from = Mortality_status, values_from = value) %>%
    mutate(Characteristic = label, .before = 1)
}

tabdat <- model_army_all %>%
  mutate(
    Mortality_status = if_else(death == 1, "Dead", "Alive"),
    cohort = as.character(cohort),
    follow_up_years = stop_age - start_age,
    age_at_death = if_else(death == 1, stop_age, NA_real_)
  )

n_row <- tabdat %>%
  count(Mortality_status) %>%
  mutate(value = as.character(n)) %>%
  select(Mortality_status, value) %>%
  pivot_wider(names_from = Mortality_status, values_from = value) %>%
  mutate(Characteristic = "N", .before = 1)

table1_dead_alive <- bind_rows(
  n_row,
  make_row_cat(tabdat, "cohort", "1948", "1948 cohort, n (%)"),
  make_row_cat(tabdat, "cohort", "1953", "1953 cohort, n (%)"),
  make_row_cont(tabdat, "child_g_z", "Childhood cognitive level, mean (SD)"),
  make_row_cont(tabdat, "army_g_z", "Late-adolescent cognitive ability, mean (SD)"),
  make_row_cont(tabdat, "army_g_resid", "Adolescent cognitive change, mean (SD)"),
  make_row_cont(tabdat, "z_parentedu_mean", "Parental education, mean (SD)"),
  make_row_cont(tabdat, "age_at_death", "Age at death, mean (SD), y"),
  make_row_cont(tabdat, "follow_up_years", "Follow-up time, mean (SD), y")
) %>%
  select(Characteristic, Alive, Dead)

cat("\n\n================ TABLE 1 DEAD VS ALIVE ================\n\n")
print(as.data.frame(table1_dead_alive), row.names = FALSE)
cat("\n================ END TABLE 1 DEAD VS ALIVE ================\n\n")

write.table(
  table1_dead_alive,
  file = file.path(outdir, "Table1_dead_alive_sample_characteristics.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)