# =============================================================================
# flowchart_UGU1948_April2026.R
# STROBE-style participant flow diagram for Walhovd et al.
#
# Run AFTER Cox_prep_12 (or Cox_prep_13) so that cox_ugu1948 and
# model_data48 are already in the environment.
# The box counts are computed from the live data so they update
# automatically when N changes (e.g. after further exclusions).
#
# Output: TIFF at 600 dpi suitable for journal submission
#   /safe/data/KBW/Figure_flowchart_600dpi.tif
# =============================================================================

library(ggplot2)
library(ragg)

# =============================================================================
# STEP 1: COMPUTE COUNTS FROM DATA
# =============================================================================
# These are derived from objects already in the environment.
# If you need to hard-code them instead, just replace each with a number.

# Box 1: All individuals in 1948 cohort with registry data
# (before any predictor missingness or opt-out exclusions)
n_registry <- nrow(ugu1948) + length(optout_ids)   # add back the 8 excluded opt-outs

# Box 2: After excluding registry opt-outs
n_after_optout <- nrow(ugu1948)
n_excluded_optout <- length(optout_ids)

# Box 3: Analytic sample — those with complete cognitive, sex, parental education data
n_analytic <- nrow(model_data48)
n_excluded_missing <- n_after_optout - n_analytic

# Box 4: Dementia cases in analytic sample
n_dementia <- sum(model_data48$dementia_any == 1L)
n_no_dementia <- n_analytic - n_dementia

# Box 5: Subsample with midlife (1990) education
# model_data48_own_edu is already filtered for complete midlife education
n_midlife_edu <- nrow(model_data48_own_edu)
n_missing_midlife <- n_analytic - n_midlife_edu
n_dementia_edu_subsample <- sum(model_data48_own_edu$dementia_any == 1L)

cat("=== Flow diagram counts ===\n")
cat(sprintf("All with registry data (incl. opt-outs): %d\n", n_registry))
cat(sprintf("Excluded (registry opt-outs, March 2026): %d\n", n_excluded_optout))
cat(sprintf("After opt-out exclusion: %d\n", n_after_optout))
cat(sprintf("Excluded (missing cognition/sex/parental edu): %d\n", n_excluded_missing))
cat(sprintf("Primary analytic sample: %d\n", n_analytic))
cat(sprintf("  Dementia cases: %d (%.1f%%)\n", n_dementia, 100*n_dementia/n_analytic))
cat(sprintf("  No dementia:    %d\n", n_no_dementia))
cat(sprintf("Midlife education subsample: %d\n", n_midlife_edu))
cat(sprintf("  Missing midlife education: %d\n", n_missing_midlife))
cat(sprintf("  Dementia cases in subsample: %d\n", n_dementia_edu_subsample))

# =============================================================================
# STEP 2: BUILD THE DIAGRAM USING ggplot2 RECTANGLES AND TEXT
# =============================================================================
# Layout uses a normalised coordinate system (x: 0–10, y: 0–10).
# Main flow goes down the centre (x = 5).
# Exclusion boxes sit to the right (x = 7.8).

# Helper: draw a box with wrapped text
# We build all elements as data frames and layer them with geom_rect + geom_text.

# --- Box definitions ---
# Each box: xmin, xmax, ymin, ymax, label (may contain "\n")

boxes <- data.frame(
  id    = c("B1", "B2", "B3", "B4", "B5", "E1", "E2"),
  xmin  = c(2.5,  2.5,  2.5,  2.5,  2.5,  6.0,  6.0),
  xmax  = c(7.5,  7.5,  7.5,  7.5,  7.5,  9.8,  9.8),
  ymin  = c(8.8,  6.8,  4.8,  2.8,  0.8,  7.2,  5.2),
  ymax  = c(9.8,  7.8,  5.8,  3.8,  1.8,  7.8,  5.8),
  label = c(
    # B1 — top: all with registry data
    sprintf("1948 cohort with registry linkage\nn = %s", format(n_registry, big.mark = ",")),
    # B2 — after opt-out exclusion
    sprintf("After excluding registry opt-outs\nn = %s", format(n_after_optout, big.mark = ",")),
    # B3 — primary analytic sample
    sprintf("Primary analytic sample\n(complete cognitive, sex, parental education data)\nn = %s  |  dementia cases: %d (%.1f%%)",
            format(n_analytic, big.mark = ","), n_dementia, 100*n_dementia/n_analytic),
    # B4 — follow-up summary
    sprintf("Follow-up to November 2025 (>60 years)\nDementia: %d  |  Censored: %d",
            n_dementia, n_no_dementia),
    # B5 — midlife education subsample
    sprintf("Midlife education subsample\n(1990 Census data available)\nn = %s  |  dementia cases: %d",
            format(n_midlife_edu, big.mark = ","), n_dementia_edu_subsample),
    # E1 — exclusion box 1
    sprintf("Excluded: %d\n(registry opt-outs, March 2026)", n_excluded_optout),
    # E2 — exclusion box 2
    sprintf("Excluded: %s\n(missing cognitive tests,\nparental education, or sex)",
            format(n_excluded_missing, big.mark = ","))
  ),
  is_exclusion = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
  stringsAsFactors = FALSE
)

# Midpoint of each box for text placement
boxes$xmid <- (boxes$xmin + boxes$xmax) / 2
boxes$ymid <- (boxes$ymin + boxes$ymax) / 2

# --- Arrow definitions ---
# Vertical arrows connecting main boxes (B1->B2->B3->B4->B5)
# Horizontal arrows from main flow to exclusion boxes (E1, E2)

arrows_v <- data.frame(
  x    = 5,
  xend = 5,
  y    = c(8.8, 6.8, 4.8, 2.8),   # start = bottom of each box
  yend = c(7.8, 5.8, 3.8, 1.8)    # end   = top of next box
)

# Horizontal arrows branching right to exclusion boxes
# Branch point is midway between the two vertical boxes
arrows_h <- data.frame(
  x    = 5,
  xend = c(6.0, 6.0),
  y    = c(7.3, 5.3),
  yend = c(7.3, 5.3)
)

# =============================================================================
# STEP 3: RENDER AND SAVE
# =============================================================================

outdir <- "/safe/data/KBW/"

agg_tiff(
  file.path(outdir, "Figure_flowchart_600dpi.tif"),
  width  = 9,
  height = 11,
  units  = "in",
  res    = 600,
  compression = "lzw"
)

p <- ggplot() +

  # --- Main boxes ---
  geom_rect(
    data = boxes[!boxes$is_exclusion, ],
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill  = "white",
    color = "black",
    linewidth = 0.6
  ) +

  # --- Exclusion boxes (lighter border) ---
  geom_rect(
    data = boxes[boxes$is_exclusion, ],
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill  = "#f5f5f5",
    color = "grey50",
    linewidth = 0.5,
    linetype  = "dashed"
  ) +

  # --- Box labels ---
  geom_text(
    data = boxes,
    aes(x = xmid, y = ymid, label = label),
    size     = 3.2,
    lineheight = 1.3,
    hjust    = 0.5,
    vjust    = 0.5,
    fontface = ifelse(boxes$is_exclusion, "plain", "plain")
  ) +

  # --- Vertical arrows (main flow) ---
  geom_segment(
    data = arrows_v,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow     = arrow(length = unit(0.2, "cm"), type = "closed"),
    linewidth = 0.6,
    color     = "black"
  ) +

  # --- Horizontal arrows to exclusion boxes ---
  geom_segment(
    data = arrows_h,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow     = arrow(length = unit(0.2, "cm"), type = "closed"),
    linewidth = 0.5,
    color     = "grey40",
    linetype  = "solid"
  ) +

  # --- Bottom annotation ---
  annotate("text",
           x = 5, y = 0.35,
           label = "Note: 1953 replication cohort (n = 9,403; 80 dementia cases) analysed separately.",
           size  = 2.8,
           color = "grey40",
           hjust = 0.5) +

  # --- Clean theme ---
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 10.5)) +
  theme_void() +
  theme(
    plot.margin = margin(10, 10, 10, 10)
  )

print(p)
dev.off()

cat(sprintf("\nFlow diagram saved to %s\n",
            file.path(outdir, "Figure_flowchart_600dpi.tif")))


# =============================================================================
# OPTIONAL: lower-res preview version (for checking on screen)
# =============================================================================

tiff(
  file.path(outdir, "Figure_flowchart_preview_150dpi.tif"),
  width  = 9,
  height = 11,
  units  = "in",
  res    = 150
)
print(p)
dev.off()
cat("Preview version (150 dpi) also saved.\n")
