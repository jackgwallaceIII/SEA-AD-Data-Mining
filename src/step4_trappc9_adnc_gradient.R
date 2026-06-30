# =============================================================================
# Step 4: TRAPPC9 Expression Across ADNC Severity Gradient
# =============================================================================
# Run this AFTER step2_TRAPPC9_seurat_analysis.R has completed and `seu`
# is loaded in your R environment with TRAPPC9 normalized expression.
#
# Tests whether TRAPPC9 expression changes progressively across the full
# ADNC severity gradient: Not AD -> Low -> Intermediate -> High
# =============================================================================

library(ggplot2)
library(ggpubr)
library(dplyr)

# Make sure ADNC is an ordered factor in the right order
seu$ADNC <- factor(seu$ADNC,
                   levels = c("Not AD", "Low", "Intermediate", "High"),
                   ordered = TRUE)

adnc_df <- data.frame(
  TRAPPC9  = GetAssayData(seu, layer = "data")["TRAPPC9", ],
  ADNC     = seu$ADNC,
  CellType = seu$Subclass
)

# =============================================================================
# 1. Pooled across cell types — boxplot across all 4 ADNC groups
# =============================================================================

p_adnc_pooled <- ggplot(adnc_df, aes(x = ADNC, y = TRAPPC9)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = NA,
               color = "black", linewidth = 1.1, staplewidth = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.06, size = 0.4, color = "grey60") +
  stat_compare_means(method = "kruskal.wallis", label = "p.format") +
  theme_classic(base_size = 12) +
  labs(y = "TRAPPC9 expression (log-normalized)", x = "ADNC severity",
       title = "TRAPPC9 expression across ADNC severity (pooled)")

ggsave(file.path(OUTPUT_DIR, "TRAPPC9_ADNC_gradient_pooled.pdf"),
       p_adnc_pooled, width = 6, height = 5, device = cairo_pdf)

p_adnc_pooled


# =============================================================================
# 2. Split by cell type — does the trend differ between Astrocyte/Immune?
# =============================================================================

p_adnc_ct <- ggplot(adnc_df, aes(x = ADNC, y = TRAPPC9)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = NA,
               color = "black", linewidth = 1.1, staplewidth = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.06, size = 0.4, color = "grey60") +
  stat_compare_means(method = "kruskal.wallis", label = "p.format") +
  facet_wrap(~ CellType) +
  theme_classic(base_size = 12) +
  labs(y = "TRAPPC9 expression (log-normalized)", x = "ADNC severity",
       title = "TRAPPC9 expression across ADNC severity by cell type")

ggsave(file.path(OUTPUT_DIR, "TRAPPC9_ADNC_gradient_by_celltype.pdf"),
       p_adnc_ct, width = 9, height = 5, device = cairo_pdf)

p_adnc_ct


# =============================================================================
# 3. Summary stats table — mean/median expression per ADNC group
# =============================================================================

adnc_summary <- adnc_df %>%
  group_by(CellType, ADNC) %>%
  summarise(
    n_cells        = n(),
    n_expressing   = sum(TRAPPC9 > 0),
    pct_expressing = round(n_expressing / n_cells * 100, 2),
    mean_expr      = round(mean(TRAPPC9), 4),
    median_expr    = round(median(TRAPPC9), 4),
    .groups = "drop"
  ) %>%
  arrange(CellType, ADNC)

print(adnc_summary)

write.csv(adnc_summary,
          file.path(OUTPUT_DIR, "TRAPPC9_ADNC_gradient_summary.csv"),
          row.names = FALSE)


# =============================================================================
# 4. Trend test — is there a monotonic trend across severity? (Jonckheere-Terpstra)
# =============================================================================
# Tests for an ordered trend rather than just "groups differ somewhere"
# More appropriate than Kruskal-Wallis when ADNC has a natural order

if (!requireNamespace("clinfun", quietly = TRUE)) {
  install.packages("clinfun")
}
library(clinfun)

message("\nJonckheere-Terpstra trend test (pooled):")
jt_pooled <- jonckheere.test(adnc_df$TRAPPC9, as.numeric(adnc_df$ADNC))
print(jt_pooled)

for (ct in unique(adnc_df$CellType)) {
  message("\nJonckheere-Terpstra trend test — ", ct, ":")
  sub <- adnc_df %>% filter(CellType == ct)
  jt_ct <- jonckheere.test(sub$TRAPPC9, as.numeric(sub$ADNC))
  print(jt_ct)
}