library(Seurat)
library(Census)
library(cowplot)

peng_sub = readRDS('peng_sub.RDS')

predlist = census_main(peng_sub, organ = 'Pancreas', predict_cancer = T)

peng_sub$Census = predlist$pred$celltype

plot_grid(DimPlot(peng_sub, group.by = 'Census', label = T) + NoLegend(),
          DimPlot(peng_sub, group.by = 'celltype', label = T) + NoLegend())
