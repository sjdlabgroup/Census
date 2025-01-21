Census: automated and hierarchical cell-type annotation for scRNA-seq
================
Bassel Ghaddar
10/20/2022

## Contents

1.  Introduction
2.  Annotating normal tissues
3.  Annotating cancer tissues
4.  Training and using a custom model
5.  Advanced options
6.  Reference

Note: a demo script (demo.r) and dataset (peng_sub.RDS) are also included in this directory.

Also note: Census requires Seurat version 4 to run. Please load Seurat v4 prior to loading the Census library in R. 

## 1. Introduction

Census is a fast and fully automated hierarchical cell-type identification method in R for single-cell RNA-seq (scRNA-seq). Briefly, Census implements a collection of hierarchically organized gradient-boosted decision tree models that successively classify individual cells according to a predefined cell-type hierarchy (Fig. 1a). New datasets are annotated using the pretrained models followed by a custom-developed label-stabilizing algorithm (Fig. 1b). This involves using clusters and prediction contours in UMAP space to propagate higher confidence binary prediction outcomes. The Census model starts at the root node of the cell-type hierarchy and successively classifies cells according to their node identities until terminal classifications are reached. The main Census model is trained on the Tabula Sapiens and Cancer Cell Line Encyclopedia and can identify 175 cell-types from 24 organs. It can also identify cancer cells and their likely cell of origin. Custom models can also be easily trained using other references. The main functionality is described below.

To install Census, please run:

``` r
remotes::install_github('sjdlabgroup/Census', ref = 'main')
```

![](Figure%201.png) <font size="2"> **Figure 1.** A schematic representation of the Census training (a) and predicting (b) algorithms

<font size="3">

To view the main Census model's cell-type hierarchy:

``` r
library(Census)
ggsave(plot_cell_hierarchy(Census:::census_ts_model$hierarchy_mat, 
                           axis.limits = c(15,-15),
                           edge.thickness = 0.5,
                           label.line.padding = 0.1,
                           text.size = 1.5, 
                           label.size = 0.75
), 
filename = 'plot.pdf', width = 3, height = 8)                 
```

<img src="Figure%202.png" height="300" />

<font size="2"> **Figure 2.** Census main model cell-type hierarchy. See file Figure 2.pdf for a high-def version.

<font size="3">

## 2. Annotating normal tissues

Fully-automated cell-type identifcation is done with the `census_main` function, which uses the pretrained Tabula Sapiens model. The input data is a Seurat object that contains a UMAP dimensionality reduction. Users should also specify the tissue organ(s), which can be or more of: Bladder, Blood, Bone\_Marrow, Eye, Fat, Heart, Kidney, Large\_Intestine, Liver, Lung, Lymph\_Node, Mammary, Muscle, Pancreas, Protate, Salivary\_Gland, Skin, Small\_Intestine, Spleen, Thymus, Tongue, Trachea, Uterus, Vasculature.

This example uses Census to predict cell-types in lung cell atlas data from Travaglini et al, Nature 2020. First the data needs to be processed using the standard Seurat pipeline to produce a UMAP dimensionality reduction. Next, run the `census_main` function for cell-type annotation.

``` r
library(Census)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

# Seurat
obj = obj %>% CreateSeuratObject() %>% NormalizeData() %>% ScaleData() %>%
  FindVariableFeatures() %>% RunPCA() %>% RunUMAP(dims=1:50)

# Census
res = census_main(obj, organ = 'Lung')

# Plot Census annotations
ggplot(res$pred, aes(umap1, umap2, color = celltype)) + 
  geom_point(size = 0.1, shape = 16) + 
  geom_label_repel(data = res$pred %>% 
                     group_by(celltype) %>% 
                     summarize(u1 = mean(umap1), u2 = mean(umap2)), 
                   force = 100, 
                   aes(u1,u2, label=celltype), 
                   size = 2, color = 'black', 
                   label.padding = 0.1, 
                   max.overlaps = 20) +
  scale_color_manual(values=colorRampPalette(brewer.pal(10, "Spectral"))(length(unique(res$pred$celltype)))) +
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(color = 'black', size = 7))
```

![](Figure%203.png)

<font size="2"> **Figure 3.** Census predictions on Travaglini et al., Nature 2020 lung cell atlas. Left, original author annotations; Right, Census predictions.

<font size="3">

## 2. Annotating cancer tissues

Cancer tissues can also be annotated using the `census_main` function. The process is the same as for normal tissues, but afterwards an organ specific model predicts whether classified epithelial cells are normal or malignant. The normal cell-type is also retained as the predicted cell of origin. Make sure to set the parameter `predict_cancer=T`. If more than one organ is specified in the `organ` parameter, the first one is used for the cancer model. Currently available cancer models include: Kidney, Large\_Intestine, Liver, Lung, Mammary, Pancreas.

This example uses Census to predict cell-types in a pancreatic cancer dataset data from Peng et al., Nature Cell Research 2019.

``` r
library(Census)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

# Seurat
obj = obj %>% CreateSeuratObject() %>% NormalizeData() %>% ScaleData() %>%
  FindVariableFeatures() %>% RunPCA() %>% RunUMAP(dims=1:50)

# Census
res = census_main(obj, organ = 'Pancreas', predict_cancer = T)

# Plot Census annotations
ggplot(res$pred, aes(umap1, umap2, color = celltype)) + 
  geom_point(size = 0.1, shape = 16) + 
  geom_label_repel(data = res$pred %>% 
                     group_by(celltype) %>% 
                     summarize(u1 = mean(umap1), u2 = mean(umap2)), 
                   force = 100, 
                   aes(u1,u2, label=celltype), 
                   size = 2, color = 'black', 
                   label.padding = 0.1, 
                   max.overlaps = 20) +
  scale_color_manual(values=colorRampPalette(brewer.pal(10, "Spectral"))(length(unique(res$pred$celltype)))) +
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(color = 'black', size = 7))

# Plot predicted cell-type of origin for cancer cells
ggplot(res$pred, aes(umap1, umap2, color = cancer_origin)) + 
  geom_point(size = 0.1, shape = 16) + 
  geom_label_repel(data = res$pred %>% 
                     group_by(celltype) %>% 
                     summarize(u1 = mean(umap1), u2 = mean(umap2)), 
                   force = 100, 
                   aes(u1,u2, label=cancer_origin), 
                   size = 2, color = 'black', 
                   label.padding = 0.1, 
                   max.overlaps = 20) +
  scale_color_manual(values=colorRampPalette(brewer.pal(10, "Spectral"))(length(unique(res$pred$celltype)))) +
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(color = 'black', size = 7))
```

![](Figure%204.png)

<font size="2"> **Figure 4.** Census predictions on the Peng et al., Nature Cell Research 2019 pancreatic cancer dataset. Left, original author annotations; Middle, Census predictions; Right, predicted cancer cell of origin

<font size="3">

## 4. Training and using a custom model

While the core Census model is trained on the Tabula Sapiens, with one line of code users can seamlessly train their own models with other references for customized applications. In this example we use lung cell atlas data from Travaglini et al., Nature 2020 to train a Census model, and we use that model to annotate lung data from Madissoon et al., Genome Biology 2020. Please see function help or Section 5 for advanced parameters on model training parameters.

The only inputs for model training are a Seurat object containing training data (UMAP not required) and a character vector of training cell-type labels. The `census_train` function computes the cell-type hierarchy and carries out node marker selection and hierarchical model training. Models, markers, the cell-type hierarchy, and model AUCs are all returned. You can view the resulting cell-type hierarchy using `plot_cell_hierarchy`.

``` r
library(Census)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

# train is a Seurat object containing training counts data from Travaglini et al., Nature 2020
# train_celltypes is a character vector of cell-types 
# see Seurat::FindMarkers for optional parameters influence marker selection and algorithm speed

# train model
model = census_train(train, celltypes = train_celltypes)

# view cell-type hierarchy
plot_cell_hierarchy(model$hierarchy_mat)

# predict new data
# test is a counts matrix of lung data from Madissoon et al., Genome Biology 2020
test = test %>% CreateSeuratObject() %>% NormalizeData() %>% ScaleData() %>% 
  FindVariableFeatures() %>% RunPCA() %>% RunUMAP(dims=1:50)

res = census_predict(test, model)

# Plot Census annotations
ggplot(res$pred, aes(umap1, umap2, color = celltype)) + 
  geom_point(size = 0.1, shape = 16) + 
  geom_label_repel(data = res$pred %>% 
                     group_by(celltype) %>% 
                     summarize(u1 = mean(umap1), u2 = mean(umap2)), 
                   force = 100, 
                   aes(u1,u2, label=celltype), 
                   size = 2, color = 'black', 
                   label.padding = 0.1, 
                   max.overlaps = 20) +
  scale_color_manual(values=colorRampPalette(brewer.pal(10, "Spectral"))(length(unique(res$pred$celltype)))) +
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(color = 'black', size = 7))
```

## 5. Advanced options

For custom model training, the `Seurat::FindMarkers` function is used to find differentially expressed genes at each node of the cell-type hierarchy. Please see its documentation for relevant parameters that can influence speed and output. We've found that adjusting the max.cells.per.ident and logfc.threshold to be useful. Model training also utilizes the `xgboost` package. Please see its documentation for adjustable parameters.

During prediction of new data, the outcome is influenced by the number of internally calculated UMAP clusters (which are included in the function output). To adjust clustering users may alter the `resolution` parameter from the `Seurat::FindClusters` function. Higher values will lead to more clusters. Also during annotation of new data when using `census_main`, the parameter "test" can be either "all" or "tabula\_sapiens". "all" will allow prediciton of all theoretically possible cell-types in the given organ. "tabula\_sapiens" will limit predictions to cell-types that were directly observed in the given organ in the Tabula Sapiens. Lastly, if `contour_data = T` during prediction with either `census_main` or `census_predict`, all prediction contours are retained and outputed from each node and round of prediction.

## 6. Reference

For any questions, please contact Bassel Ghaddar: bassel dot ghaddar at gmail dot com

Please also see the following preprint for more information!

Ghaddar, B. and De, S. (2022). Census: accurate, automated, deep, fast, and hierarchical scRNA-seq cell-type annotation. BioRxiv.
