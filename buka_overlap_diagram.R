# ----- Overlap diagram: dementia, any CVD, diabetes --------------------------
# Drop in after the cvd_any_event / diabetes_event definitions in
# Cox_prep_8_UGU_paper1_Feb_2026.R. Uses cox_ugu1948 with variables you
# already have: dementia_any, cvd_any_event, diabetes_event.

library(dplyr)

# 0. Where to put the PDFs ----------------------------------------------------
# Set this to whatever you like, e.g. the same folder as your R script.
# The default just uses the current working directory (run getwd() to see it).
out_dir <- getwd()
cat("Figures will be written to:", out_dir, "\n")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Build a 0/1 indicator matrix of the three case types ---------------------
overlap_df <- cox_ugu1948 %>%
  transmute(
    Dementia  = as.integer(dementia_any   == 1),
    `Any CVD` = as.integer(cvd_any_event  == 1),
    Diabetes  = as.integer(diabetes_event == 1)
  ) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0L, .)))

# 2. Cross-tab of all 8 possible combinations (counts) ------------------------
# (Printed as a plain data.frame to avoid the tibble 'na.print' issue.)
overlap_counts <- overlap_df %>%
  count(Dementia, `Any CVD`, Diabetes, name = "n", sort = TRUE) %>%
  as.data.frame()

print(overlap_counts)
# write.csv(overlap_counts, "overlap_counts_3sets.csv", row.names = FALSE)


# 3a. UpSet plot --------------------------------------------------------------
# install.packages("UpSetR")
library(UpSetR)

upset_input <- as.data.frame(overlap_df)

upset_pdf <- file.path(out_dir, "Figure_Sx_overlap_upset.pdf")
pdf(upset_pdf, width = 7, height = 4.5)
upset(
  upset_input,
  sets = c("Any CVD", "Diabetes", "Dementia"),
  keep.order = TRUE,
  order.by = "freq",
  mainbar.y.label = "Number of individuals",
  sets.x.label    = "Cases per condition",
  text.scale = c(1.3, 1.2, 1.2, 1.0, 1.2, 1.1)
)
dev.off()
cat("Wrote:", upset_pdf, "\n")


# 3b. Proportional Euler diagram (alternative) --------------------------------
# install.packages("eulerr")
library(eulerr)

eul_vec <- overlap_counts %>%
  mutate(label = mapply(function(d, c, db) {
    parts <- c()
    if (d  == 1) parts <- c(parts, "Dementia")
    if (c  == 1) parts <- c(parts, "Any CVD")
    if (db == 1) parts <- c(parts, "Diabetes")
    if (length(parts) == 0) return(NA_character_)  # drop the all-zero region
    paste(parts, collapse = "&")
  }, Dementia, `Any CVD`, Diabetes)) %>%
  filter(!is.na(label)) %>%
  { setNames(.$n, .$label) }

fit <- euler(eul_vec, shape = "ellipse")

euler_pdf <- file.path(out_dir, "Figure_Sx_overlap_euler.pdf")
pdf(euler_pdf, width = 6, height = 5)
plot(fit,
     quantities = list(fontsize = 10),
     labels     = list(fontsize = 11),
     fills      = list(alpha = 0.55),
     edges      = TRUE,
     main       = "Overlap of dementia, any CVD, and diabetes cases\n(1948 cohort, n=10,539)")
dev.off()
cat("Wrote:", euler_pdf, "\n")
