# TODO Implement for non-celda_C result lists
#' Calculate the silhouette width for a given celda model.
#' 
#' @param celda.mod A celda_C, celda_G, or celda_CG object 
#' @param dis.matrix A dissimilarity matrix, generated from the counts matrix used to generate the celda.mod parameter
#' @param param Whether to calculate based off the z (K) or y (L) clustering for a celda_CG object. Defaults to "z".
#' @return An object of class silhouette
#' @export
calculate_silhouette_width = function(celda.mod, dis.matrix, param="z") {
  if (!(param %in% c("z", "y"))) stop("Invalid parameter provided for 'param', must be 'z' or 'y'")
  
  clusters = switch(class(celda.mod),
                    "celda_C" = celda.mod$z,
                    "celda_G" = celda.mod$y,
                    "celda_CG"= celda.mod[[param]])
  return(cluster::silhouette(x=clusters, dist=dis.matrix))
}


#' Calculate a dissimilarity matrix between columns of a counts matrix.
#' 
#' In order to make distances between cells more comparable despite an invidiual cell's depth of sequencing,
#' this function requires that the counts matrix be normalized within cells (to CPM), 
#' and subsequently across cells. This is performed automatically (and toggled by the _normalize_ parameter) 
#' on the provided counts matrix, unless otherwise specified.
#' 
#' @param counts A counts matrix used to run celda
#' @param normalize Whether to normalize the provided matrix automatically. Defaults to TRUE.
#' @param dist.metric Which dissimilarity metric to use. Options are "euclidian", "cosine". Defaults to "euclidian"
#' @export
calculate_dissimilarity_matrix = function(counts, normalize=T, dist.metric="euclidian") {
  normalized.counts = t(scale(t(normalizeCounts(counts)), center=T, scale=T))
    
  if (dist.metric == "euclidian") {
    dis.matrix = as.matrix(t(dist(t(normalized.counts), method="euclidian")))
  }
  else if (dist.metric == "cosine") {
    dis.matrix = cosine_dissimilarity_matrix(normalized.counts)
  }
  else stop("Invalid dissimilarity metric provided.")
  return(dis.matrix)
}


# Calculate a cosine disimilarity matrix for all columns of a counts matrix
cosine_dissimilarity_matrix = function(counts) {
  # Normalize columns of counts to between 0 and 1:
  normalized.counts = apply(counts, 1, 
                            function(col) {
                              max.val = max(col)
                              min.val = min(col)
                              denom = max.val - min.val
                              normalized.col = sapply(col, function(i) {
                                return((i - min.val)/denom)
                              })
                              return(normalized.col)
                            })
  
  column.index.combinations = expand.grid(i=1:ncol(counts),
                                          j=1:ncol(counts))
  dissim.mat = matrix(apply(column.index.combinations, 1,
                      function(idx) {
                        x = counts[, idx[1]]
                        y = counts[, idx[2]]
                        return(cosine_dissimilarity(x, y))
                      }),
                      ncol(counts), ncol(counts))
  return(dissim.mat)
}


# Calculate the cosine (dis)similarity for two vectors
cosine_dissimilarity = function(x, y) {
  dissimilarity = x %*% y / sqrt((x %*% x) * (y %*% y))
  return(dissimilarity)
}


# TODO Implement for non-celda_C result lists
#' Visualize the average cluster silhouette width for all models in a celda_list.
#' 
#' Returns a ggplot of boxplots, breaking out silhouette widths by each value of K/L
#' (or combination thereof) present in the celda_list, across all chains.
#' 
#' @param celda.res A celda_list returned from celda()
#' @param counts The counts matrix used to generate celda.res
#' @param normalize.counts Whether to normalize the counts matrix, defaults to TRUE
#' @param dist.metric The distance metric to use in calculating the dissimilarity matrix. Choices are "euclidian", "cosine". Defaults to "euclidian".
#' @param title Plot title to use when visualizing silhouette width distributions
#' @export
visualize_silhouette_widths = function(celda.res, counts, normalize.counts=T,
                                       dist.metric="euclidian",
                                       title="Average Silhouette Width By Cluster Size") {
  
  if (class(celda.res) != "celda_list") stop("Must provide a celda_list object for celda.res parameter")
  if (!(compare_count_matrix(counts, celda.res$count.checksum))) stop("MD5 checksum of counts matrix doesn't match celda.res checksum.")
  
  dissimilarity.matrix = switch(celda.res$content.type,
                                "celda_C" = calculate_dissimilarity_matrix(counts, dist.metric),
                                "celda_G" = calculate_dissimilarity_matrix(t(counts), dist.metric),
                                "celda_CG" = list("cells"=calculate_dissimilarity_matrix(counts, dist.metric),
                                                  "genes"=calculate_dissimilarity_matrix(t(counts), dist.metric)))
  
  # For models with a single clustering parameter:
  if (celda.res$content.type != "celda_CG") {
    silhouette.widths = lapply(celda.res$res.list, calculate_silhouette_width,
                               dissimilarity.matrix)
    plot.df = celda.res$run.params
    avg.sil.width = lapply(silhouette.widths, function(sil) { return(summary(sil)$avg.width) })
    plot.df$avg.sil.width = as.numeric(avg.sil.width)
    
    cluster.label = ifelse(celda.res$content.type == "celda_C", "K", "L")
    return(ggplot2::ggplot(plot.df, ggplot2::aes(x=factor(K), y=avg.sil.width)) + 
            ggplot2::geom_boxplot(outlier.color=NA, fill=NA) +
            ggplot2::geom_point(position=ggplot2::position_jitter(width=0.1, height=0)) +
            ggplot2::xlab(cluster.label) + ggplot2::ylab("Average Silhouette Width") + 
            ggplot2::ggtitle(title) + ggplot2::theme_bw())
  }
  
  # For models with multiple levels of clustering (e.g. celda_CG):
  else {
    k.sil.widths = lapply(celda.res$res.list, calculate_silhouette_width,
                          dissimilarity.matrix$cells, "z")
    l.sil.widths = lapply(celda.res$res.list, calculate_silhouette_width,
                          dissimilarity.matrix$genes, "y")
    
    # TODO Determine best summary metric... right now get average width for K & L, then average
    plot.df = celda.res$run.params
    avg.k.sil.width = lapply(k.sil.widths,  function(sil) { return(summary(sil)$avg.width) })
    avg.l.sil.width = lapply(l.sil.widths,  function(sil) { return(summary(sil)$avg.width) })
    plot.df$mean.overall.sil.width = (unlist(avg.k.sil.width) + unlist(avg.l.sil.width)) / 2
    plot.df$key = paste(plot.df$K, plot.df$L, sep=",")
    
    # Order models by L, then K, in the boxplots:
    sorted.plot.df = dplyr::arrange(plot.df, L, K)  # Sort by L, then K within L
    plot.df$key = factor(plot.df$key, levels=unique(sorted.plot.df$key))
    
    return(ggplot2::ggplot(plot.df, ggplot2::aes(x=key, y=mean.overall.sil.width)) + 
            ggplot2::geom_boxplot(outlier.color=NA, fill=NA) +
            ggplot2::geom_point(position=ggplot2::position_jitter(width=0.1, height=0)) +
            ggplot2::xlab("K,L Combination") + 
            ggplot2::ylab("Average of Average K / L Silhouette Widths") + 
            ggplot2::theme(axis.text.x=element_text(angle=90, hjust=1)) +
            ggplot2::theme(axis.text.y=element_text(hjust=1)) +
            ggplot2::ggtitle(title) + ggplot2::theme_bw())
  } 
}



