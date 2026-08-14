#!/usr/bin/env Rscript

options(echo = TRUE)
Sys.time()

library(ggplot2)

cols <- c("#D44F43", "#70A1C9")

format_pvalue_labels <- function(stat_df) {
  stat_df$label <- ""
  sig <- stat_df$FDR < 0.05
  if (any(sig)) {
    e <- floor(log10(stat_df$FDR[sig]))
    m <- round(stat_df$FDR[sig] / 10^e, digits = 1)
    stat_df$label[sig] <- paste0(
      "adj. ~ italic(P) == ",
      sprintf("%g %%*%% 10^{%+d}", m, e)
    )
  }
  stat_df
}

plot_figure1c <- function(bar_data, stat_df, title, x_limit, ylab_text, out_file) {
  stat_df <- format_pvalue_labels(stat_df)

  p <- ggplot(bar_data, aes(x = percent, y = cell)) +
    geom_line(color = "#E59740", linewidth = 0.5) +
    geom_point(size = 4, alpha = 0.5, aes(fill = metabolic), color = "black", shape = 21) +
    geom_text(
      data = subset(stat_df, label != ""),
      aes(x = x.position, y = cell, label = label),
      parse = TRUE,
      size = 2.2,
      vjust = 0.5,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = x_limit) +
    scale_y_discrete(limits = rev) +
    labs(fill = "", x = ylab_text, y = "") +
    ggtitle(title) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 7, vjust = 0.5, hjust = 1, family = "Sans"),
      axis.text.y = element_text(size = 8, family = "Sans"),
      legend.text = element_text(size = 8, family = "Sans"),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 8, family = "Sans"),
      legend.position = "bottom"
    )

  svg(out_file, width = 4, height = 6)
  print(p)
  dev.off()
}

bar_data <- readRDS("../data/metabvsnonmetab.percent.rds")
stat_data <- readRDS("../data/metabvsnonmetab.fisherexacttest.rds")

plot_figure1c(
  bar_data[["CD"]], stat_data[["CD"]],
  "CD ileum", c(0, 0.085), "% of DEGs", "Figure1C.CD.svg"
)

plot_figure1c(
  bar_data[["UC"]], stat_data[["UC"]],
  "UC colon", c(0, 0.14), "% of DEGs", "Figure1C.UC.svg"
)

sessionInfo()

# Done.
