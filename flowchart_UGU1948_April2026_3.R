# =============================================================================
# flowchart_UGU1948_April2026.R
# Full STROBE-style participant flow diagram — Walhovd et al.
#
# Run AFTER Cox_prep_12 so that ugu1948, model_data48, and
# model_data48_own_edu are in the environment.
#
# Hard-coded numbers from UGU documentation:
#   https://www4.gu.se/compeat/FUR/UGU/Variabellistor/kohort_1948.html
#   n_target = 12,166 (sampling frame)
#   n_test_booklets = 10,670 (returned aptitude test booklets)
#   n_registry = 11,945 (rows in ugu1948_English_v2.4 — the analysis file)
#
# All other counts are computed dynamically from the live data.
#
# Output:
#   /safe/data/KBW/Figure_flowchart_600dpi.tif  (600 dpi, submission quality)
#   /safe/data/KBW/Figure_flowchart_preview.tif (150 dpi, screen preview)
# =============================================================================

library(ggplot2)
library(ragg)

# =============================================================================
# STEP 1: COMPUTE ALL COUNTS
# =============================================================================

# --- Hard-coded (fixed historical facts from UGU documentation) ---
n_target        <- 12166  # sampling frame: born 5th/15th/25th each month 1948
n_test_booklets <- 10670  # aptitude test booklets returned (cognitive data used in analyses)
n_registry      <- 11945  # rows in ugu1948_English_v2.4 (the starting analysis file)
# Note: n_registry > n_test_booklets because the analysis file is based on
# school administrative data (collected for more children than returned test booklets).
# Children without a returned test booklet appear in the file but are later
# excluded due to missing cognitive scores.

# --- Computed dynamically from analysis objects ---
# Opt-out exclusion: use difference rather than length(optout_ids), because
# optout_ids covers both the 1948 and 1953 cohorts combined.
n_after_optout    <- nrow(ugu1948)
n_excluded_optout <- n_registry - n_after_optout

# Analytic sample
n_analytic         <- nrow(model_data48)
n_excluded_missing <- n_after_optout - n_analytic

# Dementia outcomes
n_dementia    <- sum(model_data48$dementia_any == 1L)
n_no_dementia <- n_analytic - n_dementia

# Midlife education subsample
n_midlife_edu  <- nrow(model_data48_own_edu)
n_miss_midlife <- n_analytic - n_midlife_edu
n_dem_midlife  <- sum(model_data48_own_edu$dementia_any == 1L)

cat("=== Flow diagram counts — verify before finalising manuscript ===\n")
cat(sprintf("Sampling frame (selection days, 1948):    %d\n", n_target))
cat(sprintf("  Test booklets returned:                 %d\n", n_test_booklets))
cat(sprintf("UGU 1948 analysis file (n_registry):      %d\n", n_registry))
cat(sprintf("  Opt-outs excluded:                      %d\n", n_excluded_optout))
cat(sprintf("After opt-out exclusion:                  %d\n", n_after_optout))
cat(sprintf("  Missing predictors excluded:            %d\n", n_excluded_missing))
cat(sprintf("Primary analytic sample:                  %d\n", n_analytic))
cat(sprintf("  Dementia: %d (%.1f%%)  Censored: %d\n",
            n_dementia, 100*n_dementia/n_analytic, n_no_dementia))
cat(sprintf("Midlife education subsample:              %d\n", n_midlife_edu))
cat(sprintf("  Dementia cases in subsample:            %d\n", n_dem_midlife))


# =============================================================================
# STEP 2: DEFINE LAYOUT
# =============================================================================
# Coordinate system: x 0–13, y 0–16 (higher y = higher on page).
# Main flow: x 0.4–8.0, centred at x_c = 4.2
# Exclusion boxes: x 8.8–12.8, centred at x_e = 10.8

x_c  <- 4.2;  xl  <- 0.4;  xr  <- 8.0
x_e  <- 10.8; xe_l <- 8.8; xe_r <- 12.8

# Y positions [ymin, ymax] — boxes are tall to give text room to breathe
yB0 <- c(13.8, 15.5)   # sampling frame
yB1 <- c(11.4, 13.0)   # UGU 1948 analysis file (n_registry)
yB2 <- c( 9.0, 10.6)   # after opt-out exclusion
yB3 <- c( 5.5,  8.2)   # primary analytic sample (tallest)
yB4 <- c( 3.0,  4.6)   # midlife education subsample

# Branch points for horizontal arrows = midpoint of each vertical gap
br1 <- mean(c(yB0[1], yB1[2]))   # gap between B0 and B1
br2 <- mean(c(yB1[1], yB2[2]))   # gap between B1 and B2
br3 <- mean(c(yB2[1], yB3[2]))   # gap between B2 and B3

# Exclusion box half-height
eh <- 0.85

# Title / body text offsets within each box
t_off  <- 0.43   # title above midpoint
b_off  <- 0.32   # body below midpoint
t_B3   <- 0.75   # B3 is taller — larger offsets
b_B3   <- 0.60


# =============================================================================
# STEP 3: BUILD DATA FRAMES
# =============================================================================

# Main box rectangles
main_boxes <- data.frame(
  id   = c("B0","B1","B2","B3","B4"),
  xmin = xl, xmax = xr,
  ymin = c(yB0[1], yB1[1], yB2[1], yB3[1], yB4[1]),
  ymax = c(yB0[2], yB1[2], yB2[2], yB3[2], yB4[2]),
  fill = c("#FAFBFC","#FAFBFC","#FAFBFC","#EBF5FB","#FAFBFC"),
  stringsAsFactors = FALSE
)
main_boxes$ymid <- (main_boxes$ymin + main_boxes$ymax) / 2

# Exclusion box rectangles (3 exclusion steps)
excl_boxes <- data.frame(
  id   = c("E1","E2","E3"),
  xmin = xe_l, xmax = xe_r,
  ymin = c(br1, br2, br3) - eh,
  ymax = c(br1, br2, br3) + eh,
  stringsAsFactors = FALSE
)
excl_boxes$ymid <- (excl_boxes$ymin + excl_boxes$ymax) / 2

# Title labels (bold)
titles <- data.frame(
  x = x_c,
  y = c(
    mean(yB0) + t_off,
    mean(yB1) + t_off,
    mean(yB2) + t_off,
    mean(yB3) + t_B3,
    mean(yB4) + t_off
  ),
  label = c(
    "UGU 1948 sampling frame",
    "UGU 1948 registry-linked dataset",
    "After excluding registry opt-outs",
    "Primary analytic sample",
    "Midlife education subsample"
  ),
  stringsAsFactors = FALSE
)

# Body labels (regular)
bodies <- data.frame(
  x = x_c,
  y = c(
    mean(yB0) - b_off,
    mean(yB1) - b_off,
    mean(yB2) - b_off,
    mean(yB3) - b_B3,
    mean(yB4) - b_off
  ),
  label = c(
    # B0 — sampling frame
    sprintf("Individuals born 5th, 15th, 25th of each month in 1948 in Sweden\nn = %s  |  aptitude test booklets returned: %s",
            format(n_target, big.mark=","), format(n_test_booklets, big.mark=",")),
    # B1 — registry-linked file
    sprintf("Linked to NPR-IPR and Cause-of-Death Register\nFollow-up through November 2025 (>60 years)\nn = %s",
            format(n_registry, big.mark=",")),
    # B2 — after opt-outs
    sprintf("n = %s", format(n_after_optout, big.mark=",")),
    # B3 — analytic sample
    sprintf("Complete aptitude tests, sex, and parental education\nn = %s  |  follow-up >60 years\nDementia cases: %d (%.1f%%)  |  Censored: %s",
            format(n_analytic, big.mark=","),
            n_dementia, 100*n_dementia/n_analytic,
            format(n_no_dementia, big.mark=",")),
    # B4 — midlife education
    sprintf("1990 Census education data available\nn = %s  |  dementia cases: %d",
            format(n_midlife_edu, big.mark=","), n_dem_midlife)
  ),
  stringsAsFactors = FALSE
)

# Exclusion box labels
excl_text <- data.frame(
  x = x_e,
  y = excl_boxes$ymid,
  label = c(
    # E1 — gap between sampling frame and analysis file
    sprintf("Not in registry-linked dataset\nn = %s\n(no registry linkage or administrative data)",
            format(n_target - n_registry, big.mark=",")),
    # E2 — opt-outs
    sprintf("Excluded: %d\n(opted out of registry linkage,\nMarch 2026)",
            n_excluded_optout),
    # E3 — missing predictors
    sprintf("Excluded: %s\n(missing aptitude tests,\nparental education, or sex)",
            format(n_excluded_missing, big.mark=","))
  ),
  stringsAsFactors = FALSE
)

# Arrows
arrows_v <- data.frame(
  x=x_c, xend=x_c,
  y    = c(yB0[1], yB1[1], yB2[1], yB3[1]),
  yend = c(yB1[2], yB2[2], yB3[2], yB4[2])
)

arrows_h <- data.frame(
  x=x_c, xend=xe_l,
  y    = c(br1, br2, br3),
  yend = c(br1, br2, br3)
)


# =============================================================================
# STEP 4: RENDER
# =============================================================================

p <- ggplot() +

  # Main box fills
  geom_rect(
    data = main_boxes,
    aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=I(fill)),
    color="#1A252F", linewidth=0.75
  ) +

  # Extra border on analytic sample box (B3) to make it stand out
  geom_rect(
    data = main_boxes[main_boxes$id == "B3", ],
    aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
    fill=NA, color="#1A5276", linewidth=1.3
  ) +

  # Exclusion box fills
  geom_rect(
    data = excl_boxes,
    aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
    fill="#F4F6F7", color="#85929E", linewidth=0.6, linetype="dashed"
  ) +

  # Title text (bold)
  geom_text(
    data = titles,
    aes(x=x, y=y, label=label),
    size=4.7, fontface="bold", hjust=0.5, vjust=0.5,
    color="#1A252F", lineheight=1.3
  ) +

  # Body text (regular)
  geom_text(
    data = bodies,
    aes(x=x, y=y, label=label),
    size=4.1, fontface="plain", hjust=0.5, vjust=0.5,
    color="#2C3E50", lineheight=1.4
  ) +

  # Exclusion box text
  geom_text(
    data = excl_text,
    aes(x=x, y=y, label=label),
    size=3.8, hjust=0.5, vjust=0.5,
    color="#4D5656", lineheight=1.3
  ) +

  # Vertical arrows
  geom_segment(
    data=arrows_v,
    aes(x=x, xend=xend, y=y, yend=yend),
    arrow=arrow(length=unit(0.28,"cm"), type="closed"),
    linewidth=0.8, color="#1A252F"
  ) +

  # Horizontal arrows to exclusion boxes
  geom_segment(
    data=arrows_h,
    aes(x=x, xend=xend, y=y, yend=yend),
    arrow=arrow(length=unit(0.22,"cm"), type="open"),
    linewidth=0.6, color="#7F8C8D"
  ) +

  # Footnote
  annotate("text",
    x=x_c, y=2.2,
    label="Note: UGU 1953 replication cohort analysed separately (n \u2248 9,403; dementia cases: 80).",
    size=3.8, color="#7F8C8D", hjust=0.5, fontface="italic"
  ) +

  coord_cartesian(xlim=c(0,13), ylim=c(1.6, 16.3)) +
  theme_void() +
  theme(
    plot.margin     = margin(14, 14, 14, 14),
    plot.background = element_rect(fill="white", color=NA)
  )

# Save 600 dpi
outdir <- "/safe/data/KBW/"
agg_tiff(file.path(outdir,"Figure_flowchart_600dpi.tif"),
         width=11, height=13, units="in", res=600, compression="lzw")
print(p)
dev.off()
cat(sprintf("Saved: %s\n", file.path(outdir,"Figure_flowchart_600dpi.tif")))

# Save 150 dpi preview
tiff(file.path(outdir,"Figure_flowchart_preview.tif"),
     width=11, height=13, units="in", res=150)
print(p)
dev.off()
cat(sprintf("Preview: %s\n", file.path(outdir,"Figure_flowchart_preview.tif")))
