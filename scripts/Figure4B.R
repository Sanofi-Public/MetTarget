#!/usr/bin/env Rscript

options(echo = TRUE)
Sys.time()

library(ggplot2)
library(dplyr)
library(openxlsx)

#" ----------
#" 
#" # set parameters

select_pathway <- c('Metabolism of lipids')

mycolors <- c("#E06B3C")

plot_figure4b <- function(data_df, prefix, ridge_color) {
  data_df[, "label"] <- sapply(data_df[, "label"], function(x) {
    ifelse(is.na(x), NA, paste0("italic(", x, ")"))
  })

  dens <- stats::density(data_df$score, bw = 0.02, from = 0, to = 1, n = 512)
  max_y <- max(dens$y)
  label_df <- data_df[!is.na(data_df$label), , drop = FALSE]
  label_df$y <- stats::approx(dens$x, dens$y, xout = label_df$score)$y + 0.3 * max_y

  ggplot(data_df, aes(x = score)) +
    geom_density(
      fill = ridge_color,
      color = ridge_color,
      alpha = 0.5,
      bw = 0.02
    ) +
    geom_label(
      data = label_df,
      aes(x = score, y = y, label = label),
      inherit.aes = FALSE,
      parse = TRUE,
      size = 2,
      family = "Sans",
      label.size = 0.3,
      label.padding = unit(0.8, "lines"),
      fill = alpha("white", 0.7),
      color = "#424345"
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ) +
    ylab("Density") +
    xlab(paste0("MetTarget score (", prefix, ")")) +
    theme_bw() +
    ggtitle(paste0(
      ifelse(prefix == "CD", "Crohn's disease", "Ulcerative colitis"),
      ": ",
      select_pathway
    )) +
    theme(
      axis.text.x = element_text(size = 7, family = "Sans"),
      axis.text.y = element_text(size = 7, family = "Sans"),
      legend.text = element_text(size = 7, family = "Sans"),
      axis.title.y = element_text(size = 7, family = "Sans"),
      axis.title.x = element_text(size = 8, family = "Sans"),
      axis.line.x.bottom = element_line(size = 0.2),
      axis.line.y.left = element_line(size = 0.2),
      plot.title = element_text(size = 7, family = "Sans", hjust = 0.5),
      legend.position = "none"
    )
}

#" ----------
#"
#" # Plot CD

prefix <- "CD"
data_df <- as.data.frame(read.xlsx("../data/data.Figure4B.xlsx", sheet = prefix, colNames = TRUE))
p <- plot_figure4b(data_df, prefix, mycolors[1])

png(filename = paste0("Figure4B.", prefix, ".png"),
    width = 1400, height = 400, units = "px", res = 300)
print(p)
dev.off()

#" ----------
#"
#" # Plot UC

prefix <- "UC"
data_df <- as.data.frame(read.xlsx("../data/data.Figure4B.xlsx", sheet = prefix, colNames = TRUE))
p <- plot_figure4b(data_df, prefix, mycolors[1])

png(filename = paste0("Figure4B.", prefix, ".png"),
    width = 1400, height = 400, units = "px", res = 300)
print(p)
dev.off()

closeAllConnections()

sessionInfo()

#" ----------
#"
#" # Done.
