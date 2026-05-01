# ----- Overlap diagram: dementia, somatic morbidity, any CVD, diabetes -------
# Drop in after the cvd_any_event / diabetes_event definitions in
# Cox_prep_8_UGU_paper1_Feb_2026.R. Uses cox_ugu1948 with variables you
# already have: dementia_any, CCI_event_timebased, cvd_any_event, diabetes_event.

library(dplyr)

# 1. Build a 0/1 indicator matrix of the four cases ---------------------------
overlap_df <- cox_ugu1948 %>%
  transmute(
    Dementia          = as.integer(dementia_any           == 1),
    `Somatic morbidity` = as.integer(CCI_event_timebased    == 1),
    `Any CVD`         = as.integer(cvd_any_event          == 1),
    Diabetes          = as.integer(diabetes_event         == 1)
  ) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0L, .)))

# 2. Cross-tab of all 16 possible combinations (counts) -----------------------
# This is the underlying data that drives any Venn/UpSet plot.
overlap_counts <- overlap_df %>%
  count(Dementia, `Somatic morbidity`, `Any CVD`, Diabetes,
        name = "n", sort = TRUE)

print(overlap_counts, n = Inf)
# write.csv(overlap_counts, "overlap_counts_4sets.csv", row.names = FALSE)


# 3a. UpSet plot (recommended for 4+ sets, readable when sizes differ a lot) --
# install.packages("UpSetR")
library(UpSetR)

upset_input <- as.data.frame(overlap_df)

pdf("Figure_Sx_overlap_upset.pdf", width = 8, height = 5)
upset(
  upset_input,
  sets = c("Somatic morbidity", "Any CVD", "Diabetes", "Dementia"),
  keep.order = TRUE,
  order.by = "freq",
  mainbar.y.label = "Number of individuals",
  sets.x.label    = "Cases per condition",
  text.scale = c(1.3, 1.2, 1.2, 1.0, 1.2, 1.1)
)
dev.off()


# 3b. Alternative: proportional Euler diagram (eulerr) ------------------------
# install.packages("eulerr")
library(eulerr)

# Build a named vector of region sizes from overlap_counts.
# Region names use "&" as separator (eulerr's convention).
region <- function(d, s, c, db) {
  parts <- c()
  if (d  == 1) parts <- c(parts, "Dementia")
  if (s  == 1) parts <- c(parts, "Somatic morbidity")
  if (c  == 1) parts <- c(parts, "Any CVD")
  if (db == 1) parts <- c(parts, "Diabetes")
  if (length(parts) == 0) return(NA_character_)  # drop the all-zero region
  paste(parts, collapse = "&")
}

eul_vec <- overlap_counts %>%
  mutate(label = mapply(region,
                        Dementia, `Somatic morbidity`, `Any CVD`, Diabetes)) %>%
  filter(!is.na(label)) %>%
  { setNames(.$n, .$label) }

fit <- euler(eul_vec, shape = "ellipse")

pdf("Figure_Sx_overlap_euler.pdf", width = 7, height = 5)
plot(fit,
     quantities = list(fontsize = 9),
     labels     = list(fontsize = 11),
     fills      = list(alpha = 0.5),
     edges      = TRUE,
     main       = "Overlap of cases during follow-up (1948 cohort, n=10,539)")
dev.off()


# 4. (Optional) 3-set version focused on the headline outcomes ----------------
# Often more readable in print. Drops "Somatic morbidity" since CVD and Diabetes
# are largely subsumed within the Charlson index anyway.
eul3 <- overlap_df %>%
  count(Dementia, `Any CVD`, Diabetes, name = "n") %>%
  mutate(label = mapply(function(d,c,db) {
    parts <- c()
    if (d  == 1) parts <- c(parts, "Dementia")
    if (c  == 1) parts <- c(parts, "Any CVD")
    if (db == 1) parts <- c(parts, "Diabetes")
    if (length(parts) == 0) return(NA_character_)
    paste(parts, collapse = "&")
  }, Dementia, `Any CVD`, Diabetes)) %>%
  filter(!is.na(label)) %>%
  { setNames(.$n, .$label) }

pdf("Figure_Sx_overlap_euler_3sets.pdf", width = 6, height = 5)
plot(euler(eul3, shape = "ellipse"),
     quantities = TRUE,
     fills      = list(alpha = 0.55),
     main       = "Overlap of dementia, any CVD, and diabetes cases\n(1948 cohort, n=10,539)")
dev.off()
