#!/usr/bin/env Rscript

options(echo = TRUE)
Sys.time()

plot_title <- "Ulcerative colitis"

genes <- "SPHK1,SPHK2,S1PR1"
genes <- strsplit(genes, ",")[[1]]

options(scipen = 100, digits = 2)

#" ----------
#"
#" # Load R packages

library(openxlsx)
library(dplyr)
library(ggplot2)

source_colors <- c(
  "GSE16879" = "#3B4992",
  "GSE73661" = "#EE0000",
  "GSE23597 10 mg/kg" = "#008B45",
  "GSE23597 5 mg/kg" = "#631879",
  "pooled colon" = "#008280"
)

#" ----------
#" 
#" # Input data

data_df <- as.data.frame(read.xlsx("../data/data.Figure4F.xlsx", sheet = 1, colNames = TRUE))

data_df[, "source"] <- gsub("Meta-analysis", "pooled colon", data_df[, "source"])

data_df[, "source"] <- factor(data_df[, "source"], 
                              levels = c("GSE16879", "GSE73661", "GSE23597 10 mg/kg", "GSE23597 5 mg/kg", "pooled colon"))

data_df[, "gene"] <- factor(data_df[, "gene"], levels = genes)

data_df <- data_df %>% mutate(label = 
                                case_when(source == "pooled colon" ~
                                            paste0("P = ", ifelse(gene %in% c("SPHK1", "SPHK2", "S1PR1"), format(P, scientific = T), round(P, 3)), 
                                                   "\nadj. P = ", ifelse(gene %in% c("SPHK1", "SPHK2", "S1PR1"), format(FDR, scientific = T), round(FDR, 3)),
                                                   "\n95% CI = [", format(lower, digits = 2), 
                                                   ", ", format(upper, digits = 2), "]"), 
                                          P < 0.001 ~ paste0("P = ", format(P, scientific = T)), 
                                          TRUE ~ paste0("P = ", round(P, 3))))

x_max <- max(1, max(data_df[, "log2FC"]))
data_df <- data_df %>% mutate(
  label_nudge_x = if_else(log2FC > x_max - 0.05, -0.18, 0.12)
)

p <- ggplot(data=data_df, aes(y = source, x = log2FC, color = source)) +
  geom_point(aes(size = -log10(P)), alpha = 0.8, show.legend = TRUE) +
  geom_vline(xintercept=0, color="black", linetype="dotted", alpha=.5) +
  geom_label(
    aes(label = label, nudge_x = label_nudge_x),
    size = 2,
    alpha = 0.5,
    label.padding = unit(1.2, "lines")
  ) +
  scale_color_manual(values = source_colors) +
  scale_size_continuous(breaks = c(2, 4, 6), labels = c("2", "4", "6"), range = c(1, 9)) +
  scale_y_discrete(breaks = c("GSE16879", "GSE73661", "GSE23597 10 mg/kg", "GSE23597 5 mg/kg", "pooled colon"),
                   labels = c("colon\nGSE16879\n(n=23)", "colon\nGSE73661\n(n=23)", "colon\nGSE23597\n10 mg/kg\n(n=14)", "colon\nGSE23597\n5 mg/kg\n(n=13)", "pooled colon\n(n=73)")) +
  ylab("") +
  xlab("\nLog2 Fold Change") +
  xlim(c(-1, x_max)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=7, family="Sans"),
        axis.text.y = element_text(size=7, family="Sans"),
        axis.title.y = element_text(size=7, family="Sans"), 
        axis.title.x = element_text(size=7, family="Sans"), 
        axis.line.x.bottom = element_line(size=0.2),
        axis.line.y.left = element_line(size=0.2),
        plot.title = element_text(size=8, family="Sans", hjust = 0.5), 
        strip.text = element_text(size=7, family="Sans"),
        panel.grid.major.y = element_blank(),
        panel.spacing = unit(1, "lines"), 
        legend.text = element_text(size=7, family="Sans"),
        legend.title = element_text(size=7, family="Sans"),
        legend.position = "bottom") +
  ggtitle(plot_title) +
  guides(color = "none") + 
  facet_wrap(~ gene, ncol = 3)

if (nrow(data_df) > 0) {
  png(filename="Figure4F.png", width=2400, height=1500, units="px", res=300)
  print(p)
  dev.off()
}

sessionInfo()

#" ----------
#"
#" # Done.
