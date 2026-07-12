# =============================================================================
# SEA-AD TRAPPC9 Analysis — Step 3: Cross-Region Bar Plot
# =============================================================================
# Prerequisites: step2 must have completed for HIP, MTG, and PFC.
# Output: results/TRAPPC9_cross_region_barplot.pdf
# =============================================================================

library(ggplot2)
library(dplyr)
library(readr)

REGIONS <- c("HIP", "MTG", "PFC")
GENE    <- "TRAPPC9"

df <- bind_rows(lapply(REGIONS, function(r) {
  f <- file.path("results", r, paste0(GENE, "_expression_by_celltype_ADNC.csv"))
  if (!file.exists(f)) { message("Missing: ", f); return(NULL) }
  read_csv(f, show_col_types = FALSE) %>% mutate(Region = r)
}))

plot_df <- df %>%
  filter(Subclass %in% c("Astrocyte", "Immune"),
         ADNC     %in% c("Not AD", "High")) %>%
  mutate(
    ADNC   = factor(ADNC, levels = c("Not AD", "High"),
                    labels = c("Control", "AD (High)")),
    Region = factor(Region, levels = REGIONS)
  )

p <- ggplot(plot_df, aes(x = Region, y = pct_expressing, fill = ADNC)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = paste0(round(pct_expressing, 1), "%"),
                y = pct_expressing + 1.8),
            position = position_dodge(width = 0.8),
            size = 3.4, fontface = "bold", color = "grey20") +
  facet_wrap(~ Subclass) +
  scale_fill_manual(values = c("Control" = "#2166AC", "AD (High)" = "#D6604D")) +
  scale_y_continuous(limits = c(0, 92), expand = c(0, 0)) +
  theme_classic(base_size = 13) +
  theme(
    legend.position   = "right",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 11),
    strip.text        = element_text(face = "bold", size = 13),
    strip.background  = element_rect(fill = "grey93", color = NA),
    axis.text         = element_text(color = "black", size = 11),
    axis.title        = element_text(size = 12),
    plot.title        = element_text(face = "bold", size = 14),
    plot.subtitle     = element_text(color = "grey45", size = 11),
    panel.spacing     = unit(1.2, "lines")
  ) +
  labs(
    title    = paste("% cells expressing", GENE, "by brain region"),
    subtitle = "SEA-AD snRNA-seq — Control vs AD (High ADNC)",
    x        = "Brain Region",
    y        = paste0("% ", GENE, "+ cells")
  )

out <- file.path("results", paste0(GENE, "_cross_region_barplot.pdf"))
ggsave(out, p, width = 9, height = 5)
message("Saved: ", out)
