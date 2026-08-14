rm(list = ls())
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
param_path = "."
data_path = "../data"
#plot CD GSEA results 
plotmat_cd = readRDS(file.path(data_path,"CD.singlecell.GSEA.fig1E.rds"))
filtered_plotmat_cd <- plotmat_cd[sort(rownames(plotmat_cd)),]

hmCellFun <- function(j,i,x,y,w,h,col){
  if (!is.na(filtered_plotmat_cd[i,j])){
    if (abs(filtered_plotmat_cd[i,j]) > -log10(0.05)){
      grid.text('*',x,y,
                gp = grid::gpar(fontsize = 10),
                just = c(0.5, 0.8))
    }
  }
}

cell_cat <- c(rep('Nonimmune',3), rep('Immune', ncol(plotmat_cd)-4), 'Nonimmune')
names(cell_cat) <- colnames(plotmat_cd)

cat_labels <- factor(c('Immune','Nonimmune'), levels = c('Immune','Nonimmune'))
ha_cd <- HeatmapAnnotation(
  category = anno_block(gp = gpar(fill = c("#EFC02D","#7BD4EE")),
                        which = 'column',
                        height = unit(rep(6,3), 'mm'),
                        labels = cat_labels,
                        labels_gp = gpar(cex = 0.7))
)
ht_opt("TITLE_PADDING" = unit(0, "mm"))
breakpoints <- c(-3, -2, -1, 0, 1, 2, 3)
colors <- rev(brewer.pal(7,'RdBu'))
color_mapping <- colorRamp2(breakpoints, colors)

pushViewport(viewport(gp = gpar(fontfamily = "Sans")))
ht_cd <- Heatmap(filtered_plotmat_cd,
              col = color_mapping,
              cluster_columns = FALSE,
              column_names_side = 'bottom',
              row_names_side = 'left',
              column_split = cell_cat,
              column_gap = unit(1,'mm'),
              rect_gp = gpar(col = '#424345', lwd = 0.5),
              border =  '#424345',
              row_gap = unit(1,'mm'),
              column_names_gp = grid::gpar(fontsize = 10),
              row_names_gp = grid::gpar(fontsize = 10),
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              na_col = '#424345',
              cell_fun = hmCellFun,
              heatmap_legend_param = list(direction = 'horizontal',
                                          title = 'signed -log10(p)',
                                          title_position = 'topcenter',
                                          title_gp = gpar(fontsize = 9),
                                          labels_gp = gpar(fontsize = 9)),
              top_annotation = ha_cd,
              column_title_gp = gpar(fontsize = 10, fontface = "bold"),
              column_title = "CD ileum"
)




#plot UC GSEA results
plotmat_uc = readRDS(file.path(data_path,"UC.singlecell.GSEA.fig1E.rds"))
filtered_plotmat_uc <- plotmat_uc[sort(rownames(plotmat_uc)),]

hmCellFun <- function(j,i,x,y,w,h,col){
  if (!is.na(filtered_plotmat_uc[i,j])){
    if (abs(filtered_plotmat_uc[i,j]) > -log10(0.05)){
      grid.text('*',x,y,
                gp = grid::gpar(fontsize = 10),
                just = c(0.5, 0.8))
    }
  }
}

cell_cat <- c(rep('Nonimmune',3), rep('Immune', ncol(plotmat_uc)-3))
names(cell_cat) <- colnames(plotmat_uc)

cat_labels <- factor(c('Immune','Nonimmune'), levels = c('Immune','Nonimmune'))
ha_uc <- HeatmapAnnotation(
  category = anno_block(gp = gpar(fill = c("#EFC02D","#7BD4EE")),
                        which = 'column',
                        height = unit(rep(6,3), 'mm'),
                        labels = cat_labels,
                        labels_gp = gpar(cex = 0.7))
)


pushViewport(viewport(gp = gpar(fontfamily = "Sans")))

if(!all(rownames(filtered_plotmat_cd)==rownames(filtered_plotmat_uc))){
  warnings("rownames of CD and UC matrix are not identical!")
}

ht_uc <- Heatmap(filtered_plotmat_uc,
              col = color_mapping,
              cluster_columns = FALSE,
              column_names_side = 'bottom',
              row_names_side = 'left',
              column_split = cell_cat,
              column_gap = unit(1,'mm'),
              rect_gp = gpar(col = '#424345', lwd = 0.5),
              border =  '#424345',
              row_gap = unit(1,'mm'),
              column_names_gp = grid::gpar(fontsize = 10),
              row_names_gp = grid::gpar(fontsize = 10),
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              na_col = '#424345',
              cell_fun = hmCellFun,
              show_row_names = F,
              #heatmap_legend_param = list(direction = 'vertical',
              #                            title = '-log10(p)',
              #                           title_position = 'topcenter',
              #                            title_gp = gpar(fontsize = 9),
              #                            labels_gp = gpar(fontsize = 9)),
              show_heatmap_legend = F,
              top_annotation = ha_uc,
              column_title_gp = gpar(fontsize = 10, fontface = "bold"),
              column_title = "UC colon"
)

svg(file.path(param_path,paste0(Sys.Date(),'Fig1E.svg')),  width=14,height=9.3)
draw(ht_cd+ht_uc,heatmap_legend_side = "bottom",padding = unit(c(2, 40, 2, 2), "mm"))
dev.off()
