#!/usr/bin/env Rscript

options(echo = TRUE)
Sys.time()

library(ggplot2)
library(openxlsx)

#' ----------
#'
#' # set parameters

select_pathway <- c("Biological oxidations",
                    "Metabolism of lipids",
                    "Metabolism of vitamins and cofactors",
                    "Metabolism of porphyrins",
                    "Metabolism of amino acids and derivatives",
                    "Metabolism of carbohydrates",
                    "Metabolism of nucleotides",
                    "TCA cycle and respiratory electron transport",
                    "Transport of small molecules",
                    "Integration of energy metabolism",
                    "  ",
                    "metabolic genes",
                    "non-metabolic genes",
                    " ",
                    "")

N <- length(select_pathway) - 3

mycolors <- c("#C53B32",
              "#E06B3C",
              "#ECB14A",
              "#CCB871",
              "#72BA76",
              "#78ADAA",
              "#5CBEE6",
              "#3D3D77",
              "#992E75",
              "#F5D6A7",
              "#DA4554",
              "#3A76B7")

category_colors <- setNames(rep(mycolors, length.out = length(select_pathway)), select_pathway)

prepare_figure3a_data <- function(sheet) {
  data_df <- as.data.frame(read.xlsx("../data/data.Figure3A.xlsx", sheet = sheet, colNames = TRUE))
  data_df[which(is.na(data_df[, "category"])), "category"] <- ""
  data_df <- rbind.data.frame(
    data_df,
    data.frame("score" = NA, "category" = " "),
    data.frame("score" = NA, "category" = "  ")
  )
  data_df[, "category"] <- factor(data_df[, "category"], levels = select_pathway)
  data_df
}

prepare_ridge_df <- function(data_df, categories, bw = 0.02, rel_min_height = 0.01, scale = 0.75) {
  ridge_list <- lapply(categories, function(cat) {
    scores <- data_df$score[data_df$category == cat & !is.na(data_df$score)]
    if (length(scores) < 2) {
      return(NULL)
    }
    dens <- stats::density(scores, bw = bw, from = 0, to = 1, n = 512)
    hmax <- max(dens$y, na.rm = TRUE)
    if (!is.finite(hmax) || hmax <= 0) {
      return(NULL)
    }
    height <- dens$y / hmax * scale
    height[height / scale < rel_min_height] <- 0
    data.frame(
      category = cat,
      x = dens$x,
      height = height,
      stringsAsFactors = FALSE
    )
  })
  ridge_df <- do.call(rbind, ridge_list)
  ridge_df$category <- factor(ridge_df$category, levels = categories)
  ridge_df$y <- as.numeric(ridge_df$category)
  ridge_df
}

plot_figure3a <- function(data_df, prefix) {
  ridge_df <- prepare_ridge_df(data_df, select_pathway)

  ggplot(ridge_df, aes(x = x, group = category)) +
    geom_ribbon(
      aes(ymin = y, ymax = y + height, fill = category, colour = category),
      alpha = 0.5,
      linewidth = 0.2
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ) +
    scale_y_continuous(breaks = seq_along(select_pathway), labels = select_pathway) +
    scale_fill_manual(values = category_colors, guide = "none") +
    scale_colour_manual(values = category_colors, guide = "none") +
    ylab("") +
    xlab(paste0("MetTarget score (", prefix, ")")) +
    theme_bw() +
    ggtitle(ifelse(prefix == "CD", "Crohn's disease", "Ulcerative colitis")) +
    theme(
      axis.text.x = element_text(size = 7, family = "Sans"),
      axis.text.y = element_text(size = 7, family = "Sans"),
      legend.text = element_text(size = 7, family = "Sans"),
      axis.title.y = element_text(size = 7, family = "Sans"),
      axis.title.x = element_text(size = 8, family = "Sans"),
      axis.line.x.bottom = element_line(linewidth = 0.2),
      axis.line.y.left = element_line(linewidth = 0.2),
      plot.title = element_text(size = 7, family = "Sans", hjust = 0.5),
      legend.position = "none"
    )
}

#' ----------
#'
#' # Plot CD

prefix <- "CD"
data_df <- prepare_figure3a_data(prefix)
p <- plot_figure3a(data_df, prefix)

png(filename = paste0("Figure3A.", prefix, ".png"),
    width = 1400, height = 1600, units = "px", res = 300)
print(p)
dev.off()

#' ----------
#'
#' # Plot UC

prefix <- "UC"
data_df <- prepare_figure3a_data(prefix)
p <- plot_figure3a(data_df, prefix)

png(filename = paste0("Figure3A.", prefix, ".png"),
    width = 1400, height = 1600, units = "px", res = 300)
print(p)
dev.off()

closeAllConnections()

sessionInfo()

#' ----------
#'
#' # Done.
