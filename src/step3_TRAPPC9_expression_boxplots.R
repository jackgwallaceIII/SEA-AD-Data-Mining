# =============================================================================
# Step 3: TRAPPC9 Control vs AD Boxplot
# =============================================================================
# Run this AFTER step2_TRAPPC9_seurat_analysis.R has completed and `seu`
# is loaded in your R environment with TRAPPC9 normalized expression.
# =============================================================================

library(ggplot2)
library(ggpubr)

# Subset to just Not AD and High for the simple two-group comparison
seu_box <- subset(seu, ADNC %in% c("Not AD", "High"))
seu_box$Disease_Group2 <- ifelse(seu_box$ADNC == "Not AD", "control", "AD")

box_df <- data.frame(
  TRAPPC9 = GetAssayData(seu_box, layer = "data")["TRAPPC9", ],
  Group   = factor(seu_box$Disease_Group2, levels = c("control", "AD"))
)

p_box <- ggplot(box_df, aes(x = Group, y = TRAPPC9)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, fill = NA,
               color = "black", linewidth = 1.2,
               staplewidth = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.08, size = 0.4, color = "grey60") +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_classic(base_size = 12) +
  labs(y = "TRAPPC9 expression (log-normalized)", x = "")

ggsave(file.path(OUTPUT_DIR, "TRAPPC9_boxplot_control_vs_AD.pdf"),
       p_box, width = 5, height = 5, device = cairo_pdf)

p_box

# ---- Optional: split by cell type (Immune vs Astrocyte) -------------------
box_df_ct <- data.frame(
  TRAPPC9 = GetAssayData(seu_box, layer = "data")["TRAPPC9", ],
  Group   = factor(seu_box$Disease_Group2, levels = c("control", "AD")),
  CellType = seu_box$Subclass
)

p_box_ct <- ggplot(box_df_ct, aes(x = Group, y = TRAPPC9)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, fill = NA,
               color = "black", linewidth = 1.2,
               staplewidth = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.08, size = 0.4, color = "grey60") +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  facet_wrap(~ CellType) +
  theme_classic(base_size = 12) +
  labs(y = "TRAPPC9 expression (log-normalized)", x = "")

ggsave(file.path(OUTPUT_DIR, "TRAPPC9_boxplot_control_vs_AD_by_celltype.pdf"),
       p_box_ct, width = 8, height = 5, device = cairo_pdf)

p_box_ct