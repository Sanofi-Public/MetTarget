#!/usr/bin/env Rscript

options(echo = TRUE)
Sys.time()

#" ----------
#"
#" # Load R packages

library(ggplot2)
library(openxlsx)

#' ----------
#'
#' # Input data

pathway_df <- as.data.frame(read.xlsx("../data/data.Figure1D.xlsx", sheet = 1, colNames = TRUE))

for (col in c("padj", "imbalance")) pathway_df[, col] <- as.numeric(pathway_df[, col])
pathway_df[, "disease"] <- factor(pathway_df[, "disease"], levels = c("CD ileum", "CD colon", "UC colon"))

stopifnot(max(pathway_df[, "imbalance"]) <= 1)
stopifnot(min(pathway_df[, "imbalance"]) >= -1)

pathway_df[order(pathway_df[, 'pathway'], decreasing = T), 'order'] <- 1:nrow(pathway_df)

p <- ggplot(pathway_df, 
            aes(x = -log10(padj), y = reorder(pathway, order))) + 
  geom_bar(stat='identity', aes(fill = imbalance), width = 0.5, size = 0.2, col = '#424345', linewidth = 0.2) +
  scale_fill_gradient2(low = '#548CD8', high = '#EA3323', mid = 'white', 
                       breaks = c(-1, 0, 1), labels = c('all down', 'mixed up/down', 'all up')) +
  geom_vline(xintercept = -log10(0.05), linetype = 'dashed', size = 0.3, col = '#424345') + 
  facet_wrap(~disease) +
  xlim(0, 3.5) +
  ylab('') +
  xlab('-log10(FDR)') +
  labs(fill = 'leading-edge\nimbalance\n') +
  theme_bw() +
  theme(axis.text.x = element_text(size=5, family='Sans'),
        axis.text.y = element_text(size=5, family='Sans'),
        axis.title.y = element_text(size=6, family='Sans'), 
        axis.title.x = element_text(size=6, family='Sans'), 
        plot.title = element_text(size=6, family='Sans'),
        panel.ontop = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill=NA, colour=NA),
        strip.text = element_text(size=7, family='Sans'),
        legend.title = element_text(size=6, family='Sans'),
        legend.text = element_text(size=6, family='Sans'), 
        legend.key.height = unit(0.2, 'cm'), 
        legend.key.width = unit(0.3, 'cm'),
        legend.position = 'right')

png(filename=file.path(paste0('Figure1D.png')),
    width=2000, height=1400, units='px', res=300)
print(p)
dev.off()

sessionInfo()

#' ----------
#'
#' # Done.
