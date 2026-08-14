#!/usr/bin/env Rscript

options(echo = TRUE)
Sys.time()

#" ----------
#"
#" # Load R packages

library(ggplot2)
library(openxlsx)

#" ----------
#"
#" # Input data

tmp <- as.data.frame(read.xlsx("../data/data.Figure1B.xlsx", sheet = 1, colNames = TRUE))

selected_cols <- c("disease", "db", "X.squared", "adj.P.Value")

data_df <- NULL
for (i in 1:nrow(tmp)) {
  data_df <- rbind.data.frame(data_df, 
                             data.frame("class" = "metabolic", "percentage" = tmp[i, "pct_metabolic"], tmp[i, selected_cols]),
                             data.frame("class" = "non-metabolic", "percentage" = tmp[i, "pct_non-metabolic"], tmp[i, selected_cols]))
}

for (col in c("percentage", "X.squared", "adj.P.Value")) data_df[, col] <- as.numeric(data_df[, col])
data_df[, "disease"] <- factor(data_df[, "disease"], levels = c("CD ileum", "CD colon", "UC colon"))

ii <- which(data_df[, "class"] == "metabolic")

# label 1
data_df[ii, "label"] <- paste0("chi^2 == ", data_df[ii, "X.squared"])

# label 2
e <- floor(log10(data_df[ii, "adj.P.Value"]))
m <- round(data_df[ii, "adj.P.Value"] / 10^e, digits = 1)
data_df[ii, "label2"] <- paste0(
  "adj. ~ italic(P) == ",
  sprintf("%g %%*%% 10^{%+d}", m, e)
)

p <- ggplot(data_df, aes(x = db, y = percentage)) + 
  geom_line(color="#E59740", size=0.5, ) + 
  geom_point(size=3, alpha=1, aes(fill=class), color="black", shape=21) +
  scale_fill_manual(values = c("#D35043", "#61A3CD")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0.14, 0.42)) +
  geom_text(aes(label = label), parse = TRUE, size = 1.5, hjust = -0.05, nudge_x = -0.3, nudge_y = 0.04) +
  geom_text(aes(label = label2), parse = TRUE, size = 1.5, hjust = -0.05, nudge_x = -0.5, nudge_y = 0.025) +
  facet_grid(cols = vars(disease)) +
  ylab("% of deregulated pathway genes") +
  labs(fill="") +
  theme_bw() +
  theme(axis.text.x = element_text(size=6, family="Sans"),
        axis.text.y = element_text(size=6, family="Sans"),
        legend.text = element_text(size=7, family="Sans"), 
        axis.title.y = element_text(size=7, family="Sans"), 
        axis.title.x = element_blank(),
        legend.position = "bottom")

png(filename=file.path(paste0("Figure1B.png")),
    width=1200, height=900, units="px", res=300)
print(p)
dev.off()

sessionInfo()

#" ----------
#"
#" # Done.
