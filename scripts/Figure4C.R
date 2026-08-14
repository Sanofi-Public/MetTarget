#!/usr/bin/env Rscript

options(echo = TRUE)
Sys.time()

library(dplyr)
library(openxlsx)
library(ggplot2)

genes <- c("SPHK2", "SPHK1")

width <- length(genes) * 600

#" ----------
#" 
#" # Input data

UC_long_df <- as.data.frame(read.xlsx("../data/data.Figure4C.xlsx", sheet = "UC", colNames = TRUE))
CD_long_df <- as.data.frame(read.xlsx("../data/data.Figure4C.xlsx", sheet = "CD", colNames = TRUE))

for (col in c("value", "norm_value")) {
  UC_long_df[, col] <- as.numeric(UC_long_df[, col])
  CD_long_df[, col] <- as.numeric(CD_long_df[, col])
}

#" ----------
#" 
#" # Plot UC

set.seed(2)

prefix <- "UC"
long_df <- UC_long_df

long_df[, "gene"] <- factor(long_df[, "gene"], levels = genes)

ranked_evidence <- unique(long_df[order(long_df[, "value"], decreasing = F), "evidence"])
ranked_evidence[1:5]

long_df[, "evidence_id"] <- match(long_df[, "evidence"], ranked_evidence)

shared_features <- union(UC_long_df[, "evidence"], CD_long_df[, "evidence"])

temp_cmap <- sample(grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)], 
                    length(shared_features))

cmap <- temp_cmap[match(rev(ranked_evidence), shared_features)]

long_df <- long_df[order(long_df[, "value"], decreasing = T), ]

long_df[order(long_df[, "value"], decreasing = T), ]

num_col <- ifelse(prefix == "CD", 4, 3)
lg_height <- ifelse(prefix == "CD", 34, 5)

p <- ggplot(long_df[order(long_df[, "norm_value"], decreasing = T), ], 
            aes(fill=factor(evidence, levels = rev(ranked_evidence)), 
                y=gene, x=norm_value)) + 
  geom_bar(position="stack", stat="identity", width = 0.5, color = "black", size = 0.05, linewidth = 0.1) + 
  scale_fill_manual(values = cmap) +
  ylab("") +
  xlim(c(0, 0.8)) +
  xlab(paste0("feature scores (", prefix, ")")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=7, family="Sans"),
        axis.text.y = element_text(size=7, family="Sans"),
        legend.text = element_text(size=4, family="Sans"), 
        legend.title = element_text(size=5, family="Sans"), 
        axis.title.y = element_text(size=7, family="Sans"), 
        axis.title.x = element_text(size=7, family="Sans"), 
        axis.line.x.bottom = element_line(size=0.2),
        axis.line.y.left = element_line(size=0.2),
        plot.title = element_text(size=7, family="Sans", hjust = 0.5), 
        strip.text = element_text(size=7, family="Sans"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.key.height = unit(2, "mm"), 
        legend.key.width = unit(1, "mm"), 
        legend.key.spacing.y = unit(0, "pt"),
        legend.box.margin = margin(0, 0, lg_height, 0), 
        legend.position = "bottom") +
  guides(fill=guide_legend(ncol=num_col, byrow=FALSE, title="")) 

png(filename=paste0("Figure4C.", prefix, ".png"),
    width=1500, height=600, units="px", res = 300)
print(p)
dev.off()

#" ----------
#" 
#" # Plot CD

set.seed(2)

prefix <- "CD"
long_df <- CD_long_df

long_df[, "gene"] <- factor(long_df[, "gene"], levels = genes)

ranked_evidence <- unique(long_df[order(long_df[, "value"], decreasing = F), "evidence"])
ranked_evidence[1:5]

long_df[, "evidence_id"] <- match(long_df[, "evidence"], ranked_evidence)

shared_features <- union(UC_long_df[, "evidence"], CD_long_df[, "evidence"])

temp_cmap <- sample(grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)], 
                    length(shared_features))

cmap <- temp_cmap[match(rev(ranked_evidence), shared_features)]

long_df <- long_df[order(long_df[, "value"], decreasing = T), ]

long_df[order(long_df[, "value"], decreasing = T), ]

num_col <- ifelse(prefix == "CD", 4, 3)
lg_height <- ifelse(prefix == "CD", 34, 5)

p <- ggplot(long_df[order(long_df[, "norm_value"], decreasing = T), ], 
            aes(fill=factor(evidence, levels = rev(ranked_evidence)), 
                y=gene, x=norm_value)) + 
  geom_bar(position="stack", stat="identity", width = 0.5, color = "black", size = 0.05, linewidth = 0.1) + 
  scale_fill_manual(values = cmap) +
  ylab("") +
  xlim(c(0, 0.8)) +
  xlab(paste0("feature scores (", prefix, ")")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=7, family="Sans"),
        axis.text.y = element_text(size=7, family="Sans"),
        legend.text = element_text(size=4, family="Sans"), 
        legend.title = element_text(size=5, family="Sans"), 
        axis.title.y = element_text(size=7, family="Sans"), 
        axis.title.x = element_text(size=7, family="Sans"), 
        axis.line.x.bottom = element_line(size=0.2),
        axis.line.y.left = element_line(size=0.2),
        plot.title = element_text(size=7, family="Sans", hjust = 0.5), 
        strip.text = element_text(size=7, family="Sans"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.key.height = unit(2, "mm"), 
        legend.key.width = unit(1, "mm"), 
        legend.key.spacing.y = unit(0, "pt"),
        legend.box.margin = margin(0, 0, lg_height, 0), 
        legend.position = "bottom") +
  guides(fill=guide_legend(ncol=num_col, byrow=FALSE, title="")) 

png(filename=paste0("Figure4C.", prefix, ".png"),
    width=1500, height=600, units="px", res = 300)
print(p)
dev.off()

sessionInfo()

#" ----------
#"
#" # Done.
