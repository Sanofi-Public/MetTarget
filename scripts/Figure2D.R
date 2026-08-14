#!/usr/bin/env Rscript

options(echo = TRUE)
Sys.time()

#' ----------
#'
#' # Load R packages

library(tidyr, "1.3.1")
library(ggplot2, "4.0.0")

#' ----------
#'
#' # Input data

# UC

df <- c(32.909, 16.821, 11.139, 10.068, 2.962, 
        21.253, 10.159, 6.305, 7.231, 1.443, 
        50.958, 27.854, 19.681, 14.020, 6.083)

dim(df) <- c(5, 3)

df

colnames(df) <- c('OR', 'Lower', 'Upper')
rownames(df) <- c('MetTarget top 200 (unseen)', 'Priority index top 200', 'Open Targets top 200', 'DisGeNET', 'GWAS')
df <- cbind.data.frame('method' = rownames(df), df)

df_long <- gather(df, condition, measurement, OR:Upper, factor_key=TRUE)
df_long

UC_df_long <- df

# CD

df <- c(26.7, 14.238, 4.297, 12.521, 0.910, 
        17.020, 8.413, 1.876, 9.023, 0.289, 
        41.981, 24.095, 9.844, 17.376, 2.864)

dim(df) <- c(5, 3)

df

colnames(df) <- c('OR', 'Lower', 'Upper')
rownames(df) <- c('MetTarget top 200 (unseen)', 'Priority index top 200', 'Open Targets top 200', 'DisGeNET', 'GWAS')
df <- cbind.data.frame('method' = rownames(df), df)

library(tidyr)

df_long <- gather(df, condition, measurement, OR:Upper, factor_key=TRUE)
df_long

CD_df_long <- df

# combine

comb_df_long <- rbind.data.frame(data.frame('disease' = 'CD', CD_df_long),
                                 data.frame('disease' = 'UC', UC_df_long))

comb_df_long[, 'group'] <- paste0(comb_df_long[, 'method'], ' ', comb_df_long[, 'disease'])
# comb_df_long[, 'rank'] <- c(6, 5, 4, 2, 1, 0)
comb_df_long[, 'rank'] <- c(10, 9, 8, 7, 6, 4, 3, 2, 1, 0)
comb_df_long[, 'CI'] <- paste0('(', round(comb_df_long[, 'Lower'], digits = 1), ', ', round(comb_df_long[, 'Upper'], digits = 1), ')   ')

comb_df_long <- rbind.data.frame(comb_df_long, 
                                 data.frame('rank' = 11, 'OR' = NA, 'Lower' = 0.2, 'Upper' = NA, 
                                            'disease' = NA, 'method' = NA, 'group' = NA, 'CI' = NA))

p <- ggplot(comb_df_long, aes(x = OR, y = rank, color = method, label = method)) + 
  geom_rect(aes(xmin = 0, xmax = 1500, ymin = -0.5, ymax = 4.5), fill = "pink", alpha = 0.03, color = "white") +
  geom_rect(aes(xmin = 0, xmax = 1500, ymin = 5.5, ymax = 10.5), fill = "#A6CEE4", alpha = 0.03, color = "white") +
  geom_vline(aes(xintercept = 0), size = 0.1, linetype = "dashed", color = '#717172') + 
  geom_vline(aes(xintercept = 1), size = 0.1, linetype = "dashed", color = '#717172') + 
  geom_vline(aes(xintercept = 10), size = 0.1, linetype = "dashed", color = '#717172') + 
  geom_vline(aes(xintercept = 100), size = 0.1, linetype = "dashed", color = '#717172') +
  geom_vline(aes(xintercept = 1000), size = 0.1, linetype = "dashed", color = '#717172') +
  geom_errorbarh(aes(xmax = Upper, xmin = Lower), height = 0.5, width = 0.05, size = 0.2, alpha = 1) +
  geom_point(size = 1, shape=21, fill="white") +
  xlim(0, 1500) +
  coord_fixed(clip = 'off', ratio = 0.15) + 
  scale_x_continuous(trans='log10') +
  scale_color_manual(values = c('#C63A29', '#CB7631', '#EBAE44', '#4185B0', '#4EACB9'), 
                     breaks = c('MetTarget top 200 (unseen)', 'Priority index top 200', 'Open Targets top 200', 'DisGeNET', 'GWAS')) +
  annotate("text", x = 1700, 
           y=comb_df_long$rank, 
           label = round(comb_df_long$OR, 1), 
           size = 5/.pt, family='Sans', hjust = 0) +
  annotate("text", x = 5050, 
           y=comb_df_long$rank, 
           label = comb_df_long$CI, 
           size = 5/.pt, family='Sans', hjust = 0) +
  scale_y_continuous(breaks = c(2, 8), labels=c('UC', 'CD')) +
  theme_classic() +
  xlab('Odds Ratio (95% CI)') +
  ylab('') +
  labs(col='') +
  theme(axis.text.x = element_text(size=6, family='Sans'),
        axis.text.y = element_text(size=6, family='Sans'),
        axis.title.y = element_text(size=6, family='Sans'), 
        axis.title.x = element_text(size=6, family='Sans'), 
        plot.title = element_text(size=6, family='Sans'),
        panel.ontop = element_blank(),
        panel.background = element_blank(),
        strip.background =element_rect(fill=NA, colour=NA),
        strip.text = element_text(size=6, family='Sans'),
        legend.title = element_text(size=6, family='Sans'),
        legend.text = element_text(size=5, family='Sans'), 
        legend.key.height = unit(0.05, 'cm'), 
        legend.key.width = unit(0.5, 'cm'),
        legend.position = 'top',
        legend.box = 'vertical') +
  guides(color = guide_legend(nrow = 5))

png(filename=paste0('Figure2D.png'),
    width=1200, height=700, units='px', res=300)
print(p)
dev.off()

#' ----------
#'
#' # Done.
