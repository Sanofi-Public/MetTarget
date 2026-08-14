rm(list = ls())
param_path = "."
data_path = "../data"

library(ggplot2)

micro.df <- readRDS(file.path(data_path,"cellMicro.metab.DEG.rds"))
CD.micro.df <- micro.df[["CD"]]
UC.micro.df <- micro.df[["UC"]]

colors <- CD.micro.df$colors
names(colors) <- CD.micro.df$cell

CD.bar <- ggplot(CD.micro.df, aes(x = reorder(cell,count), y = count, fill = cell)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  theme_minimal() + 
  ggtitle("Crohn's disease") + 
  scale_fill_manual(values = colors) + 
  labs(y = "Number of DEGs \nin top 100 metabolic targets",x="") + 
  theme(plot.title = element_text(size = 12, family = 'Sans'),
        axis.title.x = element_text(size = 12, family = 'Sans'),
        axis.title.y = element_text(size = 7, family = 'Sans'),
        text=element_text(size=12, family = 'Sans'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5) )

colors <- UC.micro.df$colors
names(colors) <- UC.micro.df$cell
UC.bar <- ggplot(UC.micro.df, aes(x = reorder(cell,count), y = count, fill = cell)) +
  geom_bar(stat = 'identity') +
  ggtitle("Ulcerative colitis") + 
  coord_flip() +
  theme_minimal() + 
  scale_fill_manual(values = colors) + 
  labs(y = "Number of DEGs \nin top 100 metabolic targets",x="") + 
  theme(plot.title = element_text(size = 12, family = 'Sans'),
        axis.title.x = element_text(size = 12, family = 'Sans'),
        axis.title.y = element_text(size = 7, family = 'Sans'),
        text=element_text(size=12, family = 'Sans'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 



svg(file.path(param_path,paste0(Sys.Date(),'_CD_cellMicro_meta_DEGS.svg')), width=6.6,height=5.4)
print(CD.bar)
dev.off()

svg(file.path(param_path,paste0(Sys.Date(),'_UC_cellMicro_meta_DEGS.svg')), width=6.6,height=5.4)
print(UC.bar)
dev.off()
