# =============================================================================
# Violin plot with overlaid points and p-values by cell type and ADNC
# =============================================================================
# Run after seu is loaded with normalized TRAPPC9 data and Subclass/ADNC
# metadata has been decoded.
# =============================================================================

library(ggplot2)
library(ggpubr)

seu$ADNC <- factor(seu$ADNC,
                   levels = c("Not AD", "Low", "Intermediate", "High"),
                   ordered = TRUE)

violin_df <- data.frame(
  TRAPPC9  = GetAssayData(seu, layer = "data")["TRAPPC9", ],
  ADNC     = seu$ADNC,
  CellType = seu$Subclass
)

ADNC_COLORS <- c(
  "Not AD"       = "#2166AC",
  "Low"          = "#74ADD1",
  "Intermediate" = "#F4A582",
  "High"         = "#D6604D"
)

p_violin <- ggplot(violin_df, aes(x = ADNC, y = TRAPPC9, fill = ADNC)) +
  geom_violin(scale = "width", trim = TRUE, color = "black", linewidth = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.05, size = 0.3, color = "black") +
  stat_compare_means(method = "kruskal.test", label = "p.format",
                     label.y = 3.9) +
  facet_wrap(~ CellType) +
  scale_fill_manual(values = ADNC_COLORS) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(y = "TRAPPC9 expression (log-normalized)", x = "",
       title = "TRAPPC9 expression by cell type and ADNC — SEA-AD Hippocampus")

ggsave(file.path(OUTPUT_DIR, "TRAPPC9_violin_points_pvalue.pdf"),
       p_violin, width = 10, height = 6, device = cairo_pdf)

p_violin