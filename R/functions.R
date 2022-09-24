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
#' @importFrom igraph "graph_from_data_frame"
#' @importFrom MASS "kde2d"
#' @importFrom sp "point.in.polygon"
#' @import foreach
#' @importFrom matrixStats "rowMaxs"
#'
#' @param ... Arguments passed to other methods
#'
#' @return A data frame containing the cell hiearchy. Each column describes the node lineage for a given cell-type.
#' @export
#'
cell_hierarchy = function(counts, celltypes, counts.bulk = NULL){

  # counts is a cell by gene matrix
  # celltypes is a character vector of cell-types for each barcode in counts

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
#' @param ... Arguments passed to other methods
#'
#' @return A ggraph object plotting the cell hiearchy
#' @export
#'
plot_cell_hierarchy = function(hierarchy.mat, circular = F, n.color.per.group = NULL,
                               color = c('firebrick1', 'firebrick4', 'dodgerblue1', 'dodgerblue4',
                                         'seagreen3','darkgreen', 'slateblue1', 'slateblue4'),
                               edge.thickness = 0.5, text.size = 2, label.size = 2,
                               label.line.padding = 0.25, axis.limits = NULL){

  edges = list()
  for(j in 1:ncol(hierarchy.mat)){
    xx = hierarchy.mat[,j] %>% unique() %>% as.character()
    d = c()
    for(i in 1:(length(xx)-1)){
      d = rbind(d, data.frame(from = xx[i], to = xx[i+1], check.names = F))
    }
    edges[[j]] = d
  }
  edges = bind_rows(edges)

  graph = graph_from_data_frame(edges)

  paths = apply(hierarchy.mat, 2, function(x) paste(unique(x), collapse = ',')) %>% unlist()

  if(length(color) == 1){V(graph)$color = color}else{
    color = matrix(color, nrow = 2) %>% t()

    major_group = apply(hierarchy.mat, 2, function(x) unique(x)[nrow(color)/2+1])
    minor_group = colnames(hierarchy.mat)
    t = table(major_group, minor_group)

    graph.color = c()
    for(i in 1:nrow(t)){
      idx = which(t[i,] > 0)
      if(is.null(n.color.per.group)){n.color.per.group = length(idx)}
      c = colorRampPalette(color[i,])(n.color.per.group)
      if(length(c) < length(idx)){c = rep(c, length.out = length(idx))}

      ii = minor_group[minor_group %in% names(idx)]
      ii = ii[order(colnames(hierarchy.mat)[colnames(hierarchy.mat) %in% ii])]
      names(c) = ii

      graph.color = c(graph.color, c)
    }
    graph.color = graph.color[match(minor_group, names(graph.color))]

    names(graph.color) = minor_group
    graph.color = graph.color[match(V(graph)$name, names(graph.color))]
    V(graph)$color = graph.color
  }

  if(circular == T){
    ggraph(graph, layout = 'dendrogram', circular = circular) +
      geom_edge_elbow(width = edge.thickness) +
      geom_node_text(aes(filter = leaf,
                         label = name,
                         angle = -((-node_angle(x, y)+90)%%180)+90,
                         color = color),
                     size = text.size, hjust = 'outward') +
      geom_node_label(aes(filter = leaf == F, label = name), size = label.size,
                      label.padding = unit(label.line.padding, 'lines')) +
      scale_color_identity() +
      theme_void() +
      theme(legend.position="none",
            plot.margin=unit(c(0,0,0,0),"cm")) +
      expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))
  }
  else{
    ggraph(graph, layout = 'dendrogram', circular = circular) +
      coord_flip() +
      scale_y_reverse(limits = axis.limits) +
      geom_edge_elbow(width = edge.thickness) +
      geom_node_text(aes(filter = leaf,
                         label = name,
                         color = color),
                     size = text.size,
                     hjust = 'left') +
      geom_node_label(aes(filter = leaf == F, label = name), size = label.size,
                      label.padding = unit(label.line.padding, 'lines')) +
      scale_color_identity() +
      theme_void() +
      theme(legend.position="none",
            plot.margin=unit(c(0,0,0,0),"cm")) +
      expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))
  }
}

#' Train a classifier on a specific node in the cell hierarchy
#'
#' @param ... Arguments passed to other methods
#'
#' @return A ggraph object plotting the cell hiearchy
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
                      noise = c(-0.5, 0.5),
                      noise_frac = 0.5,
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

  # add noise
  if(!is.null(noise)){
    f = replicate(ncol(new_x), runif(nrow(new_x), noise[1], noise[2]))
    r = replicate(ncol(new_x), rbinom(nrow(new_x), 1, noise_frac))
    f = f*r; f[f == 0] = 1
    x = cbind(x, new_x*f)
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

#' Train a hierarchical Census model
#'
#' @param ... Arguments passed to other methods
#'
#' @return A ggraph object plotting the cell hiearchy
#' @export
#'
census_train = function(obj,
                        hierarchy_mat,
                        celltypes,
                        metadata = NULL,
                        markers.list = NULL,
                        sparsity = c(0.9),
                        noise = c(-0.5, 0.5),
                        noise_frac = 0.5,
                        ptile = 100,
                        nrounds = 20,
                        eval_metric = 'auc',
                        colsample = 0.2,
                        print_node = T,
                        ...){
  xg.list = list()
  markers.list = list()
  auc.df = data.frame()

  cell_node_match = match(celltypes, colnames(hierarchy_mat))

  remaining_ids = '1'
  while(length(remaining_ids) > 0){
    for(i in remaining_ids){
    # loop = foreach(i = remaining_ids, .packages = 'dplyr') %do% {
      if(print_node == T){print(i)}
      if(!is.null(markers.list)){markers = markers.list[[i]]}
      node_res = train_node(obj = obj,
                            node = i,
                            hierarchy_mat = hierarchy_mat,
                            celltypes = celltypes,
                            metadata = metadata,
                            markers = markers,
                            cell_node_match = cell_node_match,
                            sparsity = sparsity,
                            noise = noise,
                            noise_frac = noise_frac,
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

  list(models = xg.list, markers = markers.list, auc = auc.df, hierarchy_mat = hierarchy_mat)
}

#' Label stabilizing algorithm using cluster and prediction contour adjustment
#'
#' @param ... Arguments passed to other methods
#'
#' @return A ggraph object plotting the cell hiearchy
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

#' Predict cell identities for a specific node in the hierarchical model
#'
#' @param ... Arguments passed to other methods
#'
#' @return A ggraph object plotting the cell hiearchy
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

#' Automated cell-type annotation using the main Census model
#'
#' @param ... Arguments passed to other methods
#'
#' @return A ggraph object plotting the cell hiearchy
#' @export
#'
census_predict = function(obj, organ = NULL, test = 'all', predict_cancer = F, contour_data = F, verbose = T, ...){

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
    obj = obj %>% FindNeighbors(reduction='umap', dims=1:2, graph.name='census', verbose = F) %>%
      FindClusters(graph.name='census', verbose = F, ...)
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
    # if(verbose == T){str = paste0(str, ', ', 'cancer\n   '); cat(str)}
    if(verbose == T){cat(paste('\rPrediciting node:', 'cancer   '))}
    i = which(pred_res[['3']]$final_pred_df$pred == '7')
    obj = obj[, which(colnames(obj) %in% pred_res[['3']]$final_pred_df$barcode[i])]
    obj$census_celltype = '0'
    Idents(obj) = obj$census_celltype
    cancer_pred = predict_node(obj = obj, model = census_cancer_models[[organ]], node = '0')
    final_pred$cancer_origin = final_pred$celltype
    if(length(which(cancer_pred$final_pred_df$pred == 'cancer cell')) > 0){
      cancer_barcodes = cancer_pred$final_pred_df$barcode[which(cancer_pred$final_pred_df$pred == 'cancer cell')]
      final_pred$celltype[final_pred$barcode %in% cancer_barcodes] = 'cancer cell'
      final_pred$cancer_origin[final_pred$celltype != 'cancer cell'] = NA
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


