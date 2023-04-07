#' Create cell hierarchy
#'
#' @import Seurat
#' @import tidyverse
#' @import dplyr
#' @import stringr
#' @import Matrix
#' @rawNamespace import(xgboost, except = slice)
#' @importFrom ggdendro "dendro_data"
#' @importFrom dendextend "partition_leaves"
#' @import ggraph
#' @import ggplot2
#' @importFrom igraph "graph_from_data_frame"
#' @importFrom igraph "V"
#' @importFrom MASS "kde2d"
#' @importFrom sp "point.in.polygon"
#' @import foreach
#' @importFrom matrixStats "rowMaxs"
#'
#' @param counts Gene by cell counts matrix
#' @param celltypes Character vector of cell-type for each barcode in counts
#' @param counts.bulk Optional pseudo-bulk cell-type by gene counts matrix. May be useful if counts is very large and hard to transpose
#'
#' @return A data frame containing the cell hiearchy. Each column describes the node lineage for a given cell-type.
#' @export
#'
cell_hierarchy = function(counts, celltypes, counts.bulk = NULL){

  # counts is a gene by cell matrix
  # celltypes is a character vector of cell-types for each barcode in counts

  counts = Matrix::t(counts)

  if(is.null(counts.bulk)){
    bulk = list()
    for(i in unique(celltypes)){
      bulk[[i]] = data.frame(celltype = i, gene = colnames(counts), n = counts[celltypes == i, ] %>% Matrix::colSums())
    }
    counts = bulk %>%
      dplyr::bind_rows() %>%
      tidyr::pivot_wider(id_cols = celltype, names_from = gene, values_from = n, values_fill = list(n=0)) %>%
      tibble::column_to_rownames('celltype')

    counts.bulk = log2(counts/rowSums(counts)*10000+1)
  }

  h = hclust(dist(counts.bulk), method = 'ward.D2')
  dat = ggdendro::dendro_data(as.dendrogram(h), type = 'rectangle')

  subtrees = dendextend::partition_leaves(as.dendrogram(h))
  leaves = subtrees[[1]]
  pathRoutes = function(leaf){which(sapply(subtrees, function(x) leaf %in% x))}
  paths = lapply(leaves, function(x) pathRoutes(x) %>% data.frame() %>% t() %>% data.frame()) %>%
    dplyr::bind_rows() %>%
    data.frame()

  for(i in 1:nrow(paths)){
    paths[i, is.na(paths[i,])] = max(paths[i,], na.rm = T)
  }

  nr = nrow(paths); nc = ncol(paths)

  nodes = paths %>% as.matrix() %>% as.vector()

  paths = paths %>% as.matrix() %>%
    as.vector() %>%
    factor(levels = unique(nodes)) %>%
    as.numeric() %>%
    matrix(nrow = nr, ncol = nc)

  for(i in 1:nrow(paths)){
    paths[i, ] = replace(paths[i,], paths[i,] == paths[i,ncol(paths)], leaves[i])
  }

  paths = paths %>% t() %>% data.frame(stringsAsFactors = F)
  colnames(paths) = paths[nrow(paths), ]

  paths
}

#' Plot the cell hierarchy
#'
#' @param hierarchy_mat Cell-type hierarchy dataframe produced by the cell_hierarchy function
#' @param circular Plots a circular dendorgram if True
#' @param color Colors to use for cell lables. Can be a single color or multiple colors. If multiple colors are given they will be matched the major hierarchy branches.
#' @param edge.thickness Thickness of dendrogram lines
#' @param text.size Cell-type label text size
#' @param label.size Node label text size
#' @param label.line.padding Controls size of the node label outline
#' @param axis.limits Controls horizontal size of dendrogram. Takes arguements in the form "c(x,y)" where x and y can be positive or negative numbers
#'
#' @return A ggraph object plotting the cell-type hiearchy
#' @export
#'
plot_cell_hierarchy = function(hierarchy_mat, circular = F,
                               color = c('#000004FF', '#420A68FF', '#DD513AFF', '#FCA50AFF'),
                               edge.thickness = 0.5, text.size = 2, label.size = 2,
                               label.line.padding = 0.25, axis.limits = NULL){

  edges = list()
  for(j in 1:ncol(hierarchy_mat)){
    xx = hierarchy_mat[,j] %>% unique() %>% as.character()
    d = c()
    for(i in 1:(length(xx)-1)){
      d = rbind(d, data.frame(from = xx[i], to = xx[i+1], check.names = F))
    }
    edges[[j]] = d
  }
  edges = bind_rows(edges)

  graph = igraph::graph_from_data_frame(edges)

  paths = apply(hierarchy_mat, 2, function(x) paste(unique(x), collapse = ',')) %>% unlist()

  if(length(color) == 1){igraph::V(graph)$color = color}else{
    n = apply(hierarchy_mat, 1, function(x) length(unique(x)))
    n = which.min(abs(n - length(color)))[1]
    color = color[1:length(unique(unlist(hierarchy_mat[n,])))]
    color = rep(color, table(unlist(hierarchy_mat[n, ])))
    igraph::V(graph)$color = color[match(igraph::V(graph)$name, colnames(hierarchy_mat))]
  }

  if(circular == T){
    ggraph::ggraph(graph, layout = 'dendrogram', circular = circular) +
      ggraph::geom_edge_elbow(width = edge.thickness) +
      ggraph::geom_node_text(ggplot2::aes(filter = leaf,
                                          label = name,
                                          angle = -((-node_angle(x, y)+90)%%180)+90,
                                          color = color),
                             size = text.size, hjust = 'outward') +
      ggraph::geom_node_label(ggplot2::aes(filter = leaf == F, label = name), size = label.size,
                              label.padding = ggplot2::unit(label.line.padding, 'lines')) +
      ggplot2::scale_color_identity() +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position="none",
                     plot.margin=ggplot2::unit(c(0,0,0,0),"cm")) +
      ggplot2::expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))
  }
  else{
    ggraph::ggraph(graph, layout = 'dendrogram', circular = circular) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_reverse(limits = axis.limits) +
      ggraph::geom_edge_elbow(width = edge.thickness) +
      ggraph::geom_node_text(ggplot2::aes(filter = leaf,
                                          label = name,
                                          color = color),
                             size = text.size,
                             hjust = 'left') +
      ggraph::geom_node_label(ggplot2::aes(filter = leaf == F, label = name), size = label.size,
                              label.padding = ggplot2::unit(label.line.padding, 'lines')) +
      ggplot2::scale_color_identity() +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position="none",
                     plot.margin=ggplot2::unit(c(0,0,0,0),"cm")) +
      ggplot2::expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))
  }
}

#' Internal function for training a classifier on a specific node in the cell hierarchy
#'
#' @param obj Seurat object containing training data
#' @param node Character giving node identity to train
#' @param hierarchy_mat Cell-type hierarchy dataframe
#' @param celltypes Character vector of cell-type identities for all barcodes in the training data
#' @param metadata Optional metadata to use for proportional weighting of training data. Default weighting of training data is to weight inversely proportional to class size
#' @param markers Optional cell-type markers to be used for model training. Must be a dataframe with a column named "gene". If absent, Seurat's FindMarkers is used.
#' @param cell_node_match Optional custom matching of cell-type labels to cell-type hierarchy labels. This should seldom be used.
#' @param sparsity Fraction of values to replace with NA/missing values. Can be a vector or NULL
#' @param ptile Number of ntile groups to create for gene expression values; e.g. if 100, will percentile rank genes per cell, if 4 will quartile rank genes
#' @param topn Optional number of marker genes to use per model. Will selec the topn genes by fold change ranking. If NULL will use all statistically significant genes
#' @param nrounds Number of boosting rounds for model classification
#' @param eval_metric Evaluation metric for xgbosot classification algorithm
#' @param colsample Fraction of features to sample per tree. See xgboost documentation for more details for this and other possible parameters
#' @param ... Arguments passed to other methods
#'
#' @return Returns a list containing the xgboost model, a dataframe containing the marker genes data, and the classification AUC
#' @export
#'
train_node = function(obj,
                      node,
                      hierarchy_mat,
                      celltypes,
                      metadata = NULL,
                      markers = NULL,
                      cell_node_match = NULL,
                      sparsity = c(0.9),
                      ptile = 100,
                      topn = NULL,
                      nrounds = 20,
                      eval_metric = 'auc',
                      colsample = 0.2,
                      ...){

  obj$celltype = celltypes

  idx = which(hierarchy_mat == node, arr.ind = T)
  new_ids = hierarchy_mat[idx[1,1] + 1, idx[,2]] %>% as.character() %>% unique()
  if(is.null(cell_node_match)){cell_node_match = match(obj$celltype, colnames(hierarchy_mat))}
  Seurat::Idents(obj) = hierarchy_mat[idx[1,1] + 1, cell_node_match]
  obj = subset(obj, idents = new_ids)

  if(is.null(markers)){markers = Seurat::FindAllMarkers(obj, only.pos = T, ...)}
  if(is.integer(topn)){markers = markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(topn, avg_logFC)}

  x = obj@assays$RNA@counts[markers$gene, ] %>% as.matrix()
  new_x = x
  y = Seurat::Idents(obj) %>% factor(levels = Seurat::Idents(obj) %>% unique() %>% as.character() %>% sort()) %>% as.numeric() - 1

  # add sparsity
  if(!is.null(sparsity)){
    for(i in 1:length(sparsity)){
      f = replicate(ncol(new_x), rbinom(nrow(new_x), 1, 1-sparsity[i]))
      x = cbind(x, new_x*f)
    }
  }

  x[x == 0] = NA
  x = apply(x, 2, dplyr::ntile, ptile) %>% t()
  colnames(x) = unique(markers$gene)

  x = x/matrixStats::rowMaxs(x, na.rm = T)

  # set weights
  w = table(y)/sum(table(y))
  w = w[match(y, names(w))]

  if(!is.null(metadata)){
    metadata = subset(metadata, rownames(metadata) %in% colnames(obj))
    for(i in 1:ncol(metadata)){
      f = table(metadata[,i])/sum(table(metadata[,i]))
      f = f[match(metadata[,i], names(f))]
      w = w*f
    }
  }
  w = 1/w

  # make w and y match the length of x
  y = rep(y, nrow(x)/ncol(obj))
  w = rep(w, nrow(x)/ncol(obj))

  xg = xgboost::xgboost(data = x, label = y, weight = w, nrounds = nrounds, eval_metric = eval_metric, verbosity = 0,
                        colsample_bylevel=colsample, colsample_bynode=colsample, colsample_bytree=colsample, ...)

  list(model = xg, markers = markers, auc = xg$evaluation_log[nrow(xg$evaluation_log), 2] %>% as.numeric())

}

#' Train a hierarchical Census model given reference counts data and cell-type labels
#'
#' @param obj Seurat object containing training data
#' @param celltypes Character vector of cell-type identities for all barcodes in the training data
#' @param hierarchy_mat Cell-type hierarchy dataframe
#' @param metadata Optional metadata to use for proportional weighting of training data. Default weighting of training data is to weight inversely proportional to class size
#' @param markers.list Optional list containing cell-type markers to be used for model training. Each list element must contain a dataframe with a column named "gene". If absent, Seurat's FindMarkers is used. The list must be named according to the cell-type hierarchy node names.
#' @param cell_node_match Optional custom matching of cell-type labels to cell-type hierarchy labels. This should seldom be used.
#' @param sparsity Fraction of values to replace with NA/missing values. Can be a vector or NULL
#' @param ptile Number of ntile groups to create for gene expression values; e.g. if 100, will percentile rank genes per cell, if 4 will quartile rank genes
#' @param topn Optional number of marker genes to use per model. Will selec the topn genes by fold change ranking. If NULL will use all statistically significant genes
#' @param nrounds Number of boosting rounds for model classification
#' @param eval_metric Evaluation metric for xgbosot classification algorithm
#' @param colsample Fraction of features to sample per tree. See xgboost documentation for more details for this and other possible parameters
#' @param ... Arguments passed to other methods. The most relevant methods include train_node, Seurat::FindAllMarkers, and xgboost::xgboost
#'
#' @return Returns a list containing the hierarchical models, a list of dataframes containing the marker genes data for each node model, the classification AUC for each node model, and the cell-type hierarchy used
#' @export
#'
census_train = function(obj,
                        celltypes,
                        hierarchy_mat = NULL,
                        metadata = NULL,
                        markers.list = NULL,
                        sparsity = c(0.9),
                        ptile = 100,
                        nrounds = 20,
                        eval_metric = 'auc',
                        colsample = 0.2,
                        verbose = T,
                        ...){
  xg.list = list()
  markers.list = list()
  auc.df = data.frame()

  if(is.null(hierarchy_mat)){
    if(verbose == T){cat('Finding cell-type hierarchy\n')}
    hierarchy_mat = cell_hierarchy(obj@assays$RNA@counts, celltypes)
  }

  cell_node_match = match(celltypes, colnames(hierarchy_mat))

  remaining_ids = '1'
  while(length(remaining_ids) > 0){
    for(i in remaining_ids){
    # loop = foreach(i = remaining_ids, .packages = 'dplyr') %do% {
      if(verbose == T){cat(paste('\rTraining node:', i, '\n'))}
      if(!is.null(markers.list)){markers = markers.list[[i]]}
      node_res = train_node(obj = obj,
                            node = i,
                            hierarchy_mat = hierarchy_mat,
                            celltypes = celltypes,
                            metadata = metadata,
                            markers = markers,
                            cell_node_match = cell_node_match,
                            sparsity = sparsity,
                            # noise = noise,
                            # noise_frac = noise_frac,
                            ptile = ptile,
                            nrounds = nrounds,
                            eval_metric = eval_metric,
                            colsample = colsample,
                            ...)

      xg.list[[i]] = node_res$model
      markers.list[[i]] = node_res$markers
      auc.df = rbind(auc.df, data.frame(node = i, auc = node_res$auc))
    }

    idx = which(hierarchy_mat == i, arr.ind = T)
    ids = hierarchy_mat[idx[1,1]+1,] %>% as.character() %>% unique() %>% sort()
    remaining_ids = which(ids %in% celltypes == F)
    remaining_ids = ids[remaining_ids] %>% sort()
  }
  if(verbose == T){cat('\nDone\n')}
  list(models = xg.list, markers = markers.list, auc = auc.df, hierarchy_mat = hierarchy_mat)
}

#' Label stabilizing algorithm using cluster and prediction contour adjustment. Used in census_predict. First, the average label is propagated within each UMAP SNN cluster. Next, prediction contours are computed on the UMAP plot using the MASS R package. In areas where prediction contours do not overlap, all cells within the contour are given the identity of the contour. In areas where the prediction contours overlap, cells within the overlapping region are given the identity of the most common label in that region. After resolving contour disputes, the most common label is again propagated across each UMAP SNN cluster, and new prediction contours are computed. This process is repeated until either there are no more overlapping prediction contours or until there are no further changes to any cell labels.
#'
#' @param pred_df Initial prediction dataframe from census_predict that contains barcode, umap, and model prediction data.
#' @param contour_level Fraction of points to be included in prediciton contours. Default is 0.999.
#' @param min_prob_diff Probability difference from 0.5 of barcodes to include when calculating prediction contours. The default is to use all data.
#' @param contour_n Number of contour x-y coordinates to estimate.
#'
#' @return Returns a list that contains the final prediction dataframe, the prediction history during the algorithm, and a list containing the contours from each round of the algorithm.
#' @export
#'
contour_adjust = function(pred_df, contour_level = 0.999, min_prob_diff = 0, contour_n = 200){

  getLevel <- function(x,y,prob){
    kk <- MASS::kde2d(x,y)
    dx <- diff(kk$x[1:2])
    dy <- diff(kk$y[1:2])
    sz <- sort(kk$z)
    c1 <- cumsum(sz) * dx * dy
    approx(c1, sz, xout = 1 - prob)$y
  }

  pred_mat = as.matrix(pred_df$pred, ncol = 1)

  # cluster average reset
  clust_pred = data.frame()
  for(i in unique(pred_df$census_clusters)){
    clust_pred = rbind(clust_pred, data.frame(cluster = i, prob = mean(pred_df$prob[pred_df$census_clusters == i])))
  }
  clust_pred = ifelse(clust_pred$prob > 0.5, max(pred_df$pred), min(pred_df$pred))
  names(clust_pred) = unique(pred_df$census_clusters)

  pred_df$pred = clust_pred[match(pred_df$census_clusters, names(clust_pred))]

  # pred_mat = as.matrix(pred_df$pred, ncol = 1)

  pred_mat = cbind(pred_mat, pred_df$pred)
  cl.list = list()
  counter = 0

  overlap_list = list()
  overlap_barcodes = 0

  while(any(duplicated(overlap_list)) == F & length(overlap_barcodes) > 0){
    # print(counter)
    counter = counter + 1
    cl.df = data.frame()

    if(counter > 1 & c(0 %in% overlap_barcodes)){break}

    # get contour for group 1
    pred_df_sub = subset(pred_df, pred == min(pred) & abs(prob - 0.5) > min_prob_diff)
    le = getLevel(pred_df_sub$u1, pred_df_sub$u2, contour_level)
    dens = MASS::kde2d(pred_df_sub$u1,
                       pred_df_sub$u2,
                       lims = c(c(min(pred_df$u1) - 0.25*diff(range(pred_df$u1)), max(pred_df$u1) + 0.25*diff(range(pred_df$u1))),
                                c(min(pred_df$u2) - 0.25*diff(range(pred_df$u2)), max(pred_df$u2) + 0.25*diff(range(pred_df$u2)))),
                       n = contour_n)
    cl = contourLines(dens, levels = le)
    mat0 = data.frame()
    if(length(cl) > 0){
      for(i in 1:length(cl)){
        cl.df = rbind(cl.df, data.frame(n = counter, level = contour_level, type = min(pred_df$pred), x = cl[[i]]$x, y = cl[[i]]$y))
        inner = sp::point.in.polygon(pred_df$u1, pred_df$u2, cl[[i]]$x, cl[[i]]$y)
        idx = which(inner == 1)
        if(length(idx) == 0){next}
        mat0 = rbind(mat0, data.frame(barcode = pred_df$barcode[idx],
                                      pred = pred_df$pred[idx],
                                      cont.type = min(pred_df$pred),
                                      contour = i,
                                      stringsAsFactors = F))
      }
    }


    # get contour for group 2
    pred_df_sub = subset(pred_df, pred == max(pred_df$pred) & abs(prob-0.5) > min_prob_diff)
    le = getLevel(pred_df_sub$u1, pred_df_sub$u2, contour_level)
    dens = MASS::kde2d(pred_df_sub$u1,
                       pred_df_sub$u2,
                       lims = c(c(min(pred_df$u1) - 0.25*diff(range(pred_df$u1)), max(pred_df$u1) + 0.25*diff(range(pred_df$u1))),
                                c(min(pred_df$u2) - 0.25*diff(range(pred_df$u2)), max(pred_df$u2) + 0.25*diff(range(pred_df$u2)))),
                       n = contour_n)
    cl2 = contourLines(dens, levels = le)
    mat = data.frame()
    if(length(cl2) > 0){
      for(i in 1:length(cl2)){
        cl.df = rbind(cl.df, data.frame(n = counter, level = contour_level, type = max(pred_df$pred), x = cl2[[i]]$x, y = cl2[[i]]$y))
        inner = sp::point.in.polygon(pred_df$u1, pred_df$u2, cl2[[i]]$x, cl2[[i]]$y)
        idx = which(inner == 1)
        if(length(idx) == 0){next}
        mat = rbind(mat, data.frame(barcode = pred_df$barcode[idx],
                                    pred = pred_df$pred[idx],
                                    cont.type = max(pred_df$pred),
                                    contour = i,
                                    stringsAsFactors = F))
      }
      if(nrow(mat) == 0 | nrow(mat0) == 0){
        if(nrow(mat) == 0 & nrow(mat0) == 0){
          break
        }
        mat = ifelse(nrow(mat) > 0, mat, mat0)
      } else {mat = full_join(mat, mat0, by = 'barcode') %>% tibble()}
    } else{mat = mat0}

    new_pred = dplyr::select(pred_df, barcode, pred)

    # adjust predictions in non-overlapping contours
    switch_bc = mat$barcode[which(mat$pred.x != mat$cont.type.x & is.na(mat$pred.y))]
    if(length(switch_bc) > 0){
      new_pred$pred[new_pred$barcode %in% switch_bc] = max(pred_df$pred)
    }
    switch_bc = mat$barcode[which(mat$pred.y != mat$cont.type.y & is.na(mat$pred.x))]
    if(length(switch_bc) > 0){
      new_pred$pred[new_pred$barcode %in% switch_bc] = min(pred_df$pred)
    }

    # adjust predictions in overlapping contours
    if('pred.x' %in% colnames(mat)){
      overlap_mat = subset(mat, !is.na(pred.x) & !is.na(pred.y))

      for(i in unique(overlap_mat$contour.x)){
        xx = subset(overlap_mat, contour.x == i)
        for(j in unique(xx$contour.y)){
          xy = subset(xx, contour.y == j)
          t = table(xy$pred.x)
          new_pred$pred[new_pred$barcode %in% xy$barcode] = names(t)[which.max(t)] %>% as.character()
        }
      }

      for(i in unique(overlap_mat$contour.y)){
        xx = subset(overlap_mat, contour.y == i)
        for(j in unique(xx$contour.x)){
          xy = subset(xx, contour.x == j)
          t = table(xy$pred.x)
          new_pred$pred[new_pred$barcode %in% xy$barcode] = names(t)[which.max(t)] %>% as.character()
        }
      }

      overlap_list[[counter]] = overlap_mat$barcode
      overlap_barcodes = overlap_mat$barcode
      pred_df$pred = new_pred$pred
    }

    cl.df$type = as.character(cl.df$type)
    cl.list[[counter]] = cl.df

    # cluster adjust
    for(i in unique(pred_df$census_clusters)){
      t = table(pred_df$pred[pred_df$census_clusters == i])
      pred_df$pred[pred_df$census_clusters == i] = names(t)[which.max(t)] %>% as.character()
    }

    pred_mat = cbind(pred_mat, pred_df$pred)

  }

  rownames(pred_mat) = pred_df$barcode
  list(final_pred_df = pred_df, pred_mat = pred_mat, cl.list = cl.list)
}

#' Predict cell identities for a specific node in the hierarchical model. Used internally census_predict
#'
#' @param obj Seurat object containing training data
#' @param model Model list that contains the hierarchcial models, markers, cell-type hierarchy that is the output of census_train.
#' @param node Character identity of the node to predict
#' @param get_prob If False, the main Tabula Sapiens model will bypass models if one node is not an allowed node. If True, the model will still perform classification at that node and will report probabilities, and unallowed nodes will be adjusted.
#' @param allowed_nodes For the Tabula Sapiens model, if any nodes are not allowed for a specific organ their prediction will accordingly be adjusted.
#'
#' @return Returns the output from the label stabilizing algorithm: a list that contains the final prediction dataframe, the prediction history during the algorithm, and a list containing the contours from each round of the algorithm.
#' @export
#'
predict_node = function(obj, model, node, get_prob = T, allowed_nodes = NULL){

  idx = which(model$hierarchy_mat == node, arr.ind = T)
  new_ids = model$hierarchy_mat[idx[1,1] + 1, idx[,2]] %>% as.character() %>% unique()

  if(!is.null(allowed_nodes) & any(new_ids %in% allowed_nodes == F)){
    obj = obj[, Idents(obj) == node]
    pred_res = data.frame(barcode = colnames(obj),
                          pred = new_ids[new_ids %in% allowed_nodes],
                          prob = NA,
                          census_clusters = obj$census_clusters,
                          u1 = obj@reductions$umap@cell.embeddings[,1],
                          u2 = obj@reductions$umap@cell.embeddings[,2])
  }

  if(any(new_ids %in% allowed_nodes == F) & get_prob == F){
    return(pred_res)
  } else {
    g = intersect(model$markers[[node]]$gene, rownames(obj))
    temp_seurat = obj[g, Idents(obj) == node]

    x2 = temp_seurat@assays$RNA@counts %>% as.matrix()

    # add missing genes
    g = setdiff(model$markers[[node]]$gene, rownames(obj))
    if(length(g) > 0){
      x2 = rbind(x2, matrix(0, nrow = length(g), ncol = ncol(x2)))
      rownames(x2) = c(rownames(temp_seurat), g)
    }

    x2 = x2[model$markers[[node]]$gene, ]

    x2[x2 == 0] = NA

    x2 = apply(x2, 2, dplyr::ntile, 100) %>% t()
    colnames(x2) = model$markers[[node]]$gene

    x2 = x2/matrixStats::rowMaxs(x2, na.rm = T)

    pred.prob = predict(model$models[[node]], x2)
    pred.class = ifelse(pred.prob > 0.5, max(new_ids), min(new_ids))

    pred_df = data.frame(barcode = colnames(temp_seurat),
                         pred = pred.class,
                         prob = pred.prob,
                         census_clusters = temp_seurat$census_clusters,
                         u1 = temp_seurat@reductions$umap@cell.embeddings[,1],
                         u2 = temp_seurat@reductions$umap@cell.embeddings[,2],
                         stringsAsFactors = F)

    # if(!is.null(allowed_nodes) & any(pred_df$pred %in% allowed_nodes == F)){
    #   pred_df$pred = new_ids[new_ids %in% allowed_nodes]
    # }

    pred_res = contour_adjust(pred_df)

    return(pred_res)
  }
}

#' Automated cell-type annotation using the main Census model trained on the Tabula Sapiens and Cancer Cell Line Encyclopedia
#'
#' @param obj Seurat object containing training data
#' @param organ A vector specifying the organ(s) of the tissue. Can be more than one organ, and if no organ is chosen then all organs and cell-types will be possible prediction results. Allowed organs include: Bladder, Blood, Bone_Marrow, Eye, Fat, Heart, Kidney, Large_Intestine, Liver, Lung, Lymph_Node, Mammary, Muscle, Pancreas, Protate, Salivary_Gland, Skin, Small_Intestine, Spleen, Thymus, Tongue, Trachea, Uterus, Vasculature
#' @param test A character specifying either "all" or "tabula_sapiens". "all" will allow prediciton of all theoretically possible cell-types in the given organ. "tabula_sapiens" will limit predictions to cell-types that were directly observed in the given organ in the Tabula Sapiens.
#' @param predict_cancer If False, then only normal cell types will be predicted. If cancer cells are expected, set to True. In this case, all cells initially undergo cell-type prediction using the Tabula Sapiens model. Then for those that are classified as epithelial, the organ specific cancer model is applied to predict cancer cells. If multiple organs are specified, the first organ will be used for the cancer model. Cancer models are available for the following organs: Kidney, Large_Intestine, Liver, Lung, Mammary, Pancreas.
#' @param controu_data If True, contour data will be retained, if False it will not be output
#' @param verbose If True, some outputs will be printed to the console.
#'
#' @param ... Arguments passed to other methods. The most useful may be the "resolution" parameter for clustering using Seruat::FindClusters.
#'
#' @return Returns a list that contains the prediction results. The list cnontains the final prediction dataframe, class history of all barcodes, probability history, otional prediction contour data, and a record of adjusted nodes
#' @export
#'
census_main = function(obj, organ = NULL, test = 'all', predict_cancer = F, contour_data = F, verbose = T, ...){

  if('umap' %in% names(obj@reductions) == F){
    if(verbose == T){cat('Running Seurat \n')}
    obj = obj %>% NormalizeData(verbose = F) %>% ScaleData(verbose = F) %>%
      FindVariableFeatures(verbose = F) %>% RunPCA(verbose = F) %>%
      RunUMAP(dims=1:50, verbose = F) %>%
      FindNeighbors(reduction='umap', dims=1:2, verbose = F) %>%
      FindClusters(verbose = F)
  }

  # get umap clusters
  if('census_clusters' %in% colnames(obj@meta.data) == F){
    if(verbose == T){cat('Finding umap clusters \n')}
    obj = obj %>% FindNeighbors(reduction='umap', dims=1:2, verbose = F) %>%
      FindClusters(verbose = F, ...)
    obj$census_clusters = obj$seurat_clusters
  }


  beginning.time = Sys.time()
  if(verbose == T){cat('Started cell-type annotation \n')}


  # subset object for speed
  g = lapply(census_ts_model$markers, function(y) y$gene) %>% unlist() %>% unique()
  obj = obj[intersect(rownames(obj), g), ]

  pred_res = list()
  class_hist = data.frame(barcode = colnames(obj))
  prob_hist = data.frame(barcode = colnames(obj))
  final_pred = data.frame(barcode = colnames(obj), cluster = obj$census_clusters, celltype = NA)

  obj$census_celltype = '1'
  Idents(obj) = obj$census_celltype
  remaining_ids = '1'

  if(!is.null(organ)){
    if(test == 'all'){
      possible_celltypes = rownames(cell_organ_df)[as.matrix(which(cell_organ_df[, organ] > 0, arr.ind = T))[,1]]
      allowed_nodes = census_ts_model$hierarchy_mat[, possible_celltypes] %>% unlist() %>% unname() %>% unique()
    }
    if(test == 'tabula_sapiens'){
      possible_celltypes = rownames(cell_organ_df)[as.matrix(which(cell_organ_df[, organ] > 1, arr.ind = T))[,1]]
      allowed_nodes = census_ts_model$hierarchy_mat[, possible_celltypes] %>% unlist() %>% unname() %>% unique()
    }
    if(test %in% c('all', 'tabula_sapiens') == F){stop("test must be either 'all' or 'tabula_sapiens'")}
  } else {allowed_nodes = NULL}

  # str = '\rNodes predicted: 1'
  while(length(remaining_ids) > 0){
    for(i in remaining_ids){
    # loop = foreach(i = remaining_ids, .packages = 'dplyr') %do% {
      # if(verbose == T & i > 1){cat(paste(str, '   ')); str = paste0(str, ', ', i)}
      if(verbose == T){cat(paste('\rPredicting node:', i, '   '))}
      pred_res[[i]] = predict_node(obj = obj, model = census_ts_model, node = i, allowed_nodes = allowed_nodes)
      obj$census_celltype[pred_res[[i]]$final_pred_df$barcode] = pred_res[[i]]$final_pred_df$pred
      final_pred$celltype[final_pred$barcode %in% pred_res[[i]]$final_pred_df$barcode] = pred_res[[i]]$final_pred_df$pred
      class_hist = left_join(class_hist, data.frame(obj$census_celltype) %>% tibble::rownames_to_column('barcode'), by = 'barcode')
      prob_hist = left_join(prob_hist, dplyr::select(pred_res[[i]]$final_pred_df, barcode, prob), by = 'barcode')
      Idents(obj) = obj$census_celltype
    }
    remaining_ids = obj$census_celltype[obj$census_celltype %in% colnames(census_ts_model$hierarchy_mat) == F] %>% unique() %>% gtools::mixedsort()
    if(!is.null(allowed_nodes)){remaining_ids = intersect(remaining_ids, allowed_nodes)}
  }

  # orginial_pred_res = pred_res
  if(!is.null(allowed_nodes) & any(final_pred$celltype %in% allowed_nodes == F)){
    obj2 = obj[, which(obj$census_celltype %in% allowed_nodes == F)]
  }

  adjusted_nodes = data.frame(stringsAsFactors = F)
  while(!is.null(allowed_nodes) & any(final_pred$celltype %in% allowed_nodes == F)){
    # switch ids
    switch_ids = obj2$census_celltype[which(obj2$census_celltype %in% allowed_nodes == F)] %>% unique() %>% as.character()
    for(i in switch_ids){
      idx = which(census_ts_model$hierarchy_mat == i, arr.ind = T)
      old_id = census_ts_model$hierarchy_mat[idx[1,1] - 1, idx[,2]] %>% as.character() %>% unique()
      idx = which(census_ts_model$hierarchy_mat == old_id, arr.ind = T)
      new_ids = census_ts_model$hierarchy_mat[idx[1,1] + 1, idx[,2]] %>% as.character() %>% unique()
      new_ids = new_ids[new_ids %in% allowed_nodes]

      adjusted_nodes = rbind(adjusted_nodes, data.frame(node = old_id,
                                                        cluster = pred_res[[old_id]]$final_pred_df$census_clusters %>% unique(),
                                                        old = i,
                                                        new = new_ids,
                                                        stringsAsFactors = F))
      pred_res[[old_id]]$final_pred_df$pred = new_ids

      class_hist = left_join(class_hist, data.frame(obj2$census_celltype, stringsAsFactors = F) %>%
                               tibble::rownames_to_column('barcode'),
                             by = 'barcode')

      final_pred$celltype = str_replace(final_pred$celltype, paste0('^', i, '$'), new_ids)
      obj2$census_celltype = str_replace(obj2$census_celltype, paste0('^', i, '$'), new_ids)
    }

    # subset and re-evaluate
    remaining_ids = obj2$census_celltype[obj2$census_celltype %in% colnames(census_ts_model$hierarchy_mat) == F] %>% unique() %>% gtools::mixedsort()
    remaining_ids = intersect(remaining_ids, allowed_nodes)

    if(length(remaining_ids) == 0){break}

    Idents(obj2) = obj2$census_celltype

    while(length(remaining_ids) > 0){
      for(i in remaining_ids){
      # loop = foreach(i = remaining_ids, .packages = 'dplyr') %do% {
        # if(verbose == T){cat(paste(str, '   ')); str = paste0(str, ', ', i)}
        if(verbose == T){cat(paste('\rPredicting node:', i, '   '))}
        pred_res[[i]] = predict_node(obj = obj2, model = census_ts_model, node = i, allowed_nodes = allowed_nodes)
        obj2$census_celltype[pred_res[[i]]$final_pred_df$barcode] = pred_res[[i]]$final_pred_df$pred
        final_pred$celltype[final_pred$barcode %in% pred_res[[i]]$final_pred_df$barcode] = pred_res[[i]]$final_pred_df$pred
        Idents(obj2) = obj2$census_celltype
        class_hist = left_join(class_hist, data.frame(obj2$census_celltype) %>% tibble::rownames_to_column('barcode'), by = 'barcode')
        prob_hist = left_join(prob_hist, dplyr::select(pred_res[[i]]$final_pred_df, barcode, prob), by = 'barcode')
      }
      remaining_ids = obj2$census_celltype[obj2$census_celltype %in% colnames(census_ts_model$hierarchy_mat) == F] %>% unique() %>% gtools::mixedsort()
      if(!is.null(allowed_nodes)){remaining_ids = intersect(remaining_ids, allowed_nodes)}
    }

  }

  if(predict_cancer == T){
    # get epithelial cells
    if(verbose == T){cat(paste('\rPrediciting node:', 'cancer   '))}
    i = which(pred_res[['3']]$final_pred_df$pred == '7')
    if(length(i) > 0){
      obj = obj[, which(colnames(obj) %in% pred_res[['3']]$final_pred_df$barcode[i])]
      obj$census_celltype = '0'
      Idents(obj) = obj$census_celltype
      cancer_pred = predict_node(obj = obj, model = census_cancer_models[[organ[1]]], node = '0')
      final_pred$cancer_origin = final_pred$celltype
      if(length(which(cancer_pred$final_pred_df$pred == 'cancer cell')) > 0){
        cancer_barcodes = cancer_pred$final_pred_df$barcode[which(cancer_pred$final_pred_df$pred == 'cancer cell')]
        final_pred$celltype[final_pred$barcode %in% cancer_barcodes] = 'cancer cell'
        final_pred$cancer_origin[final_pred$celltype != 'cancer cell'] = NA
      }
    }
  }

  if(verbose == T){cat(paste('\nDone; total prediction time =', format(Sys.time() - beginning.time, digits=3), '\n'))}

  prob_hist = tibble::column_to_rownames(prob_hist, 'barcode')
  colnames(prob_hist) = 1:ncol(prob_hist)

  class_hist = tibble::column_to_rownames(class_hist, 'barcode')
  colnames(class_hist) = 1:ncol(class_hist)

  final_pred = data.frame(final_pred,
                          umap1 = pred_res[['1']]$final_pred_df$u1,
                          umap2 = pred_res[['1']]$final_pred_df$u2)

  if(contour_data == F){pred_res = NULL}

  list(pred = final_pred, class_hist = class_hist, prob_hist = prob_hist, pred_data = pred_res, adjusted_nodes = adjusted_nodes)
}


#' Use a Census model to annotate new datasets
#'
#' @param obj Seurat object containing training data
#' @param model Hierarchical models list that is the output of census_train. Contains the models, markers, and cell-type hierarchy
#' @param contour_data If True, contour data are outputed
#' @param verbose if True, messages are printed to the console during prediction.
#' @param ... Arguments passed to other methods
#'
#' @return A ggraph object plotting the cell hiearchy
#' @export
#'
census_predict = function(obj, model, contour_data = F, verbose = T, ...){

  if('umap' %in% names(obj@reductions) == F){
    if(verbose == T){cat('Running Seurat \n')}
    obj = obj %>% NormalizeData(verbose = F) %>% ScaleData(verbose = F) %>%
      FindVariableFeatures(verbose = F) %>% RunPCA(verbose = F) %>%
      RunUMAP(dims=1:50, verbose = F) %>%
      FindNeighbors(reduction='umap', dims=1:2,graph.name='census', verbose = F) %>%
      FindClusters(graph.name='census', verbose = F)
  }

  # get umap clusters
  if('census_clusters' %in% colnames(obj@meta.data) == F){
    if(verbose == T){cat('Finding umap clusters \n')}
    obj = obj %>% FindNeighbors(reduction='umap', dims=1:2, verbose = F) %>%
      FindClusters(verbose = F, ...)
    obj$census_clusters = obj$seurat_clusters
  }


  beginning.time = Sys.time()
  if(verbose == T){cat('Started cell-type annotation \n')}


  # subset object for speed
  g = lapply(model$markers, function(y) y$gene) %>% unlist() %>% unique()
  obj = obj[intersect(rownames(obj), g), ]

  pred_res = list()
  class_hist = data.frame(barcode = colnames(obj))
  prob_hist = data.frame(barcode = colnames(obj))
  final_pred = data.frame(barcode = colnames(obj), cluster = obj$census_clusters, celltype = NA)

  obj$census_celltype = '1'
  Idents(obj) = obj$census_celltype
  remaining_ids = '1'
  allowed_nodes = NULL

  # str = '\rNodes predicted: 1'
  while(length(remaining_ids) > 0){
    for(i in remaining_ids){
      # loop = foreach(i = remaining_ids, .packages = 'dplyr') %do% {
      # if(verbose == T & i > 1){cat(paste(str, '   ')); str = paste0(str, ', ', i)}
      if(verbose == T){cat(paste('\rPredicting node:', i, '   '))}
      pred_res[[i]] = predict_node(obj = obj, model = model, node = i, allowed_nodes = allowed_nodes)
      obj$census_celltype[pred_res[[i]]$final_pred_df$barcode] = pred_res[[i]]$final_pred_df$pred
      final_pred$celltype[final_pred$barcode %in% pred_res[[i]]$final_pred_df$barcode] = pred_res[[i]]$final_pred_df$pred
      class_hist = left_join(class_hist, data.frame(obj$census_celltype) %>% tibble::rownames_to_column('barcode'), by = 'barcode')
      prob_hist = left_join(prob_hist, dplyr::select(pred_res[[i]]$final_pred_df, barcode, prob), by = 'barcode')
      Idents(obj) = obj$census_celltype
    }
    remaining_ids = obj$census_celltype[obj$census_celltype %in% colnames(model$hierarchy_mat) == F] %>% unique() %>% gtools::mixedsort()
  }

  if(verbose == T){cat(paste('\nDone; total prediction time =', format(Sys.time() - beginning.time, digits=3), '\n'))}

  prob_hist = tibble::column_to_rownames(prob_hist, 'barcode')
  colnames(prob_hist) = 1:ncol(prob_hist)

  class_hist = tibble::column_to_rownames(class_hist, 'barcode')
  colnames(class_hist) = 1:ncol(class_hist)

  final_pred = data.frame(final_pred,
                          umap1 = pred_res[['1']]$final_pred_df$u1,
                          umap2 = pred_res[['1']]$final_pred_df$u2)

  if(contour_data == F){pred_res = NULL}

  list(pred = final_pred, class_hist = class_hist, prob_hist = prob_hist, pred_data = pred_res)
}



