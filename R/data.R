#' @title Possible cell-types per organ
#'
#' @description A data frame describing the organs in which each cell-type may occur. A value of '1' indicates the cell-type was directly detected in the given organ in Tabula Sapiens. A value of '2' indicates that the cell-type may be detected in other organs but was not detected in the Tabula Sapiens. This generally applies to most immune, stromal, and endothelial cell-types.
#'
#' @format A cell-type x organ data frame
"cell_organ_df"

#' @title Census Tabula Sapiens model
#'
#' @description Hierarchical gradient-boosted tree model trained on the Tabula Sapiens
#'
#' @format A list containing standard model outputs: models (xgboost tree model for each node), markers (marker genes used in each model), auc (model performance on training data), hierarchy_mat (cell-type hierarchy used by model)
"census_ts_model"

#' @title Census cancer models
#'
#' @description Collection of models to detect cancer cells in the following organs: breast, colon, kidney, liver, lung, pancreas
#'
#' @format A list of lists for each model. The model list format is the same as ts.model
"census_cancer_models"
