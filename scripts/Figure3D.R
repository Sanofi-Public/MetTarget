#!/usr/bin/env Rscript

options(echo = TRUE)
Sys.time()

library(openxlsx)
library(dplyr)
library(ggplot2)

df <- as.data.frame(read.xlsx("../data/data.Figure3D.xlsx", colNames = TRUE))
df[1:3, ]

#" ----------
#"
#" # Plot

UC_top100_thr <- 0.724218777
CD_top100_thr <- 0.799496205

has_label <- !is.na(df$label) & nzchar(df$label)
df[, "type"] <- ifelse(has_label & df$label != "SPHK1", "#22538D", "#C53B32")

label_df <- df %>%
  filter(has_label) %>%
  group_by(UC = round(UC, 3)) %>%
  mutate(
    label_nudge_x = 0.02,
    label_nudge_y = (row_number() - (n() + 1) / 2) * 0.025
  ) %>%
  ungroup()

p <- ggplot(df, aes(x = UC, y = CD)) +
  geom_point(aes(color = type), size = 0.8, alpha = 0.7, shape = 1, show.legend = FALSE) +
  scale_color_manual(values = c("#22538D" = "#22538D", "#C53B32" = "#C53B32")) +
  geom_label(
    data = label_df,
    aes(x = UC, y = CD, label = label, color = type,
        nudge_x = label_nudge_x, nudge_y = label_nudge_y),
    inherit.aes = FALSE,
    size = 1.8,
    family = "Sans",
    label.size = 0.3,
    label.padding = unit(0.5, "lines"),
    fill = alpha("white", 0.35),
    alpha = 0.5,
    show.legend = FALSE
  ) +
  xlim(c(0.2, 1.01)) +
  ylim(c(0.2, 1.01)) +
  geom_vline(aes(xintercept = UC_top100_thr), linewidth = 0.2, linetype = "dashed", color = "#717172") +
  geom_hline(aes(yintercept = CD_top100_thr), linewidth = 0.2, linetype = "dashed", color = "#717172") +
  ylab("MetTarget score (CD)") +
  xlab("MetTarget score (UC)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6, family = "Sans"),
        axis.text.y = element_text(size = 6, family = "Sans"),
        axis.title.y = element_text(size = 7, family = "Sans"),
        axis.title.x = element_text(size = 7, family = "Sans"),
        axis.line.x.bottom = element_line(linewidth = 0.2),
        axis.line.y.left = element_line(linewidth = 0.2),
        plot.title = element_text(size = 7, family = "Sans", hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "top")

png(filename = "Figure3D.png",
    width = 1200, height = 1200, units = "px", res = 300)
print(p)
dev.off()

#" ----------
#"
#" # Done.
