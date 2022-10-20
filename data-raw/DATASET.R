## code to prepare `DATASET` dataset goes here

library(usethis)
library(stringr)
library(dplyr)

use_data_raw()

# usethis::use_data(DATASET, overwrite = TRUE)

## make cell_organ_df, a data frame that describes which cell types can occur in which organs

# 1 = cell-types that were detected in a given organ in Tabula Sapiens
# 2 = cell-types that could occur in a given organ

ts.meta = readRDS('./data-raw/ts.meta.RDS')
cell_organ_df = table(ts.meta$cell_ontology_class %>%
                        str_remove_all('\\,') %>%
                        str_replace('t follicular helper t cells', 'cd4-positive follicular helper t cell') %>%
                        str_replace('cd4-positive helper t cells', 'cd4-positive helper t cell') %>%
                        str_replace('alpha-beta t cell', 't cell') %>%
                        str_replace('leucocyte', 'leukocyte') %>%
                        str_replace('pancreatic stellate cell', 'stellate cell') %>%
                        str_replace('lymphatic endothelial cells', 'endothelial cell of lymphatic vessel') %>%
                        str_replace('artery endothelial cell', 'endothelial cell of artery'),
                      ts.meta$organ_tissue) %>% as.data.frame.matrix()
cell_organ_df[cell_organ_df > 0] = 2

cell_organ_df['adipocyte', cell_organ_df['adipocyte', ] != 2] = 1
cell_organ_df['b cell', cell_organ_df['b cell', ] != 2] = 1
cell_organ_df['basophil', cell_organ_df['basophil', ] != 2] = 1
cell_organ_df['capillary aerocyte', cell_organ_df['capillary aerocyte', ] != 2] = 1
cell_organ_df['capillary endothelial cell', cell_organ_df['capillary endothelial cell', ] != 2] = 1
j = str_which(rownames(cell_organ_df), '^cd'); for(i in j){cell_organ_df[i, cell_organ_df[i, ] != 2] = 1 }
j = str_which(rownames(cell_organ_df), 'monocyte'); for(i in j){cell_organ_df[i, cell_organ_df[i, ] != 2] = 1 }
j = str_which(rownames(cell_organ_df), 'dendritic'); for(i in j){cell_organ_df[i, cell_organ_df[i, ] != 2] = 1 }
cell_organ_df['endothelial cell', cell_organ_df['endothelial cell', ] != 2] = 1
cell_organ_df['endothelial cell of artery', cell_organ_df['endothelial cell of artery', ] != 2] = 1
cell_organ_df['endothelial cell of lymphatic vessel', cell_organ_df['endothelial cell of lymphatic vessel', ] != 2] = 1
cell_organ_df['endothelial cell of vascular tree', cell_organ_df['endothelial cell of vascular tree', ] != 2] = 1
cell_organ_df['erythrocyte', cell_organ_df['erythrocyte', ] != 2] = 1
cell_organ_df['fibroblast', cell_organ_df['fibroblast', ] != 2] = 1
cell_organ_df['immature natural killer cell', cell_organ_df['immature natural killer cell', ] != 2] = 1
cell_organ_df['immune cell', cell_organ_df['immune cell', ] != 2] = 1
cell_organ_df['innate lymphoid cell', cell_organ_df['innate lymphoid cell', ] != 2] = 1
cell_organ_df['keratinocyte', 'Skin'] = 1
cell_organ_df['leukocyte', cell_organ_df['leukocyte', ] != 2] = 1
cell_organ_df['macrophage', cell_organ_df['macrophage', ] != 2] = 1
cell_organ_df['mast cell', cell_organ_df['mast cell', ] != 2] = 1
cell_organ_df['mature nk t cell', cell_organ_df['mature nk t cell', ] != 2] = 1
cell_organ_df['memory b cell', cell_organ_df['memory b cell', ] != 2] = 1
j = str_which(rownames(cell_organ_df), 'myeloid'); for(i in j){cell_organ_df[i, cell_organ_df[i, ] != 2] = 1 }
cell_organ_df['myofibroblast cell', cell_organ_df['myofibroblast cell', ] != 2] = 1
cell_organ_df['neutrophil', cell_organ_df['neutrophil', ] != 2] = 1
cell_organ_df['stellate cell', 'Liver'] = 1
cell_organ_df['pericyte cell', cell_organ_df['pericyte cell', ] != 2] = 1
cell_organ_df['plasma cell', cell_organ_df['plasma cell', ] != 2] = 1
cell_organ_df['plasmablast', cell_organ_df['plasmablast', ] != 2] = 1
cell_organ_df['platelet', cell_organ_df['platelet', ] != 2] = 1
cell_organ_df['regulatory t cell', cell_organ_df['regulatory t cell', ] != 2] = 1
cell_organ_df['schwann cell', cell_organ_df['schwann cell', ] != 2] = 1
cell_organ_df['smooth muscle cell', cell_organ_df['smooth muscle cell', ] != 2] = 1
cell_organ_df['t cell', cell_organ_df['t cell', ] != 2] = 1
cell_organ_df['t follicular helper cell', cell_organ_df['t follicular helper cell', ] != 2] = 1
cell_organ_df['type i nk t cell', cell_organ_df['type i nk t cell', ] != 2] = 1
cell_organ_df['vascular associated smooth muscle cell', cell_organ_df['vascular associated smooth muscle cell', ] != 2] = 1
cell_organ_df['vein endothelial cell', cell_organ_df['vein endothelial cell', ] != 2] = 1
cell_organ_df['vein or capillary endothelial cell', cell_organ_df['vein or capillary endothelial cell', ] != 2] = 1


## Main Tabula Sapiens model

ts.model = readRDS('./data-raw/ts.model7300.pctrank.RDS')
x = ts.model$hierarchy_mat
x = apply(x, 2, function(y) str_replace(y, 'pancreatic stellate cell', 'stellate cell'))
colnames(x) = str_replace(colnames(x), 'pancreatic stellate cell', 'stellate cell')
x = data.frame(x, stringsAsFactors = F, check.names = F)
ts.model$hierarchy_mat = x

census_ts_model = ts.model

## Cancer models
cancer_models = readRDS('./data-raw/cancer_models.RDS')

census_cancer_models = cancer_models

use_data(cell_organ_df, census_ts_model, census_cancer_models, internal = T, overwrite = T)



