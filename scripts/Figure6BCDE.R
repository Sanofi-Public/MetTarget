#!/usr/bin/env Rscript

options(echo = TRUE)
Sys.time()

library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)
library(reshape2)

theme_pubr_border <- function() {
  theme_classic() +
    theme(
      axis.line = element_line(linewidth = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
}

pf543 <- read.xlsx("../data/Supplemental.Data.2.xlsx", sheet = 1, rows = c(1:5))
colnames(pf543)[1] <- "Compound"

fty720 <- read.xlsx("../data/Supplemental.Data.2.xlsx", sheet = 1, rows = c(7:11))
colnames(fty720)[1] <- "Compound"

df <- bind_rows(pf543, fty720)

df_lf <- melt(df, id.vars = "Compound", value.name = "Log Ratio")
df_lf <- df_lf %>%
  separate_wider_delim(
    cols = variable,
    delim = ":",
    names = c("System", "Biomarker")
  )

df_lf$System <- factor(
  df_lf$System,
  levels = c("3C", "4H", "LPS", "SAg", "BT", "BF4T",
             "BE3C", "CASM3C", "HDF3CGF", "KF3CT", "MyoF", "lMphg")
)

plot_figure6 <- function(data, colors, out_file, width = 24, height = 8) {
  p <- ggplot(data, aes(x = Biomarker, y = `Log Ratio`, color = Compound, group = Compound)) +
    geom_line(linewidth = 1) +
    facet_grid(~System, scales = "free_x", space = "free_x") +
    theme_pubr_border() +
    ylim(c(-1.5, 1.0)) +
    scale_color_manual(values = colors) +
    labs(color = "Profiles") +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  svg(out_file, width = width, height = height)
  print(p)
  dev.off()
}

plot_figure6(
  df_lf[grepl("SAR_Compound", df_lf$Compound), , drop = FALSE],
  c("red3", "orange2", "gold", "green4"),
  "Figure6B.svg"
)

plot_figure6(
  df_lf[grepl("SAR_Compound 3, 10000 nM|FTY720, 1100 nM", df_lf$Compound), , drop = FALSE],
  c("black", "red3"),
  "Figure6C.svg"
)

plot_figure6(
  df_lf[grepl("SAR_Compound", df_lf$Compound) & df_lf$System %in% c("BT", "HDF3CGF"), , drop = FALSE],
  c("red3", "orange2", "gold", "green4"),
  "Figure6D.svg",
  width = 12
)

plot_figure6(
  df_lf[grepl("SAR_Compound 3, 10000 nM|FTY720, 1100 nM", df_lf$Compound) &
          df_lf$System %in% c("BT", "HDF3CGF"), , drop = FALSE],
  c("black", "red3"),
  "Figure6E.svg",
  width = 12
)

sessionInfo()

# Done.
