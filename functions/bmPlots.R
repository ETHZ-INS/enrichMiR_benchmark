#' TPvFP
#'
#' @param df benchmark dataframe
#' @param labels vector of method labels to be displayed in the plot; default NULL
#' @param perm binary; add when enrichMiR was run with miRNA target annotation permutations; default FALSE
#' @param label.all binary; add when you want all methods to be labelled; default FALSE
#' @param full.scale binary; FALSE: automatic y-axis scale adjustment; TRUE: from 0 to 1
#'
#' @return a ggplot of TP rates vs. total number of FPs (at FDR<.05)
#'
TPvFP <- function(df, labels=NULL, perm=FALSE, label.all=FALSE, full.scale=TRUE, features=NULL){
  
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
  })
  
  if(label.all){
    labels <- unique(df$method) 
  }
  # aggregate to get a single scores per method (& permutation)
  if(perm){
    df.agg <- aggregate(df[,c("FP.atFDR05","TP.atFDR05")], by=df[,c("prop","method",features)], FUN=mean)
    label.data <- subset(df.agg, method %in% labels) %>% filter(prop=="original")
  } else {
    if(is.null(features)){
      df.agg <- aggregate(df[,c("FP.atFDR05","TP.atFDR05")], by=list(method=df$method), FUN=mean)
      label.data <- subset(df.agg, method %in% labels)
      prop <- NULL
    } else {
      df.agg <- aggregate(df[,c("FP.atFDR05","TP.atFDR05")], by=df[,c("method",features)], FUN=mean)
      label.data <- subset(df.agg, method %in% labels)
      prop <- NULL
    }
  }
  
  if(full.scale){
    scale <- scale_y_continuous(limits = c(0,1))
  } else scale <- NULL
  
  if(!is.null(features)){
    wrap <- facet_wrap(as.formula(paste0("~", paste(features, collapse="+"))))
  } else wrap <- NULL
  
  # create plot
  set.seed(123)
  ggplot(df.agg, aes(x=FP.atFDR05, y=TP.atFDR05, color=method, shape=prop, group=method)) + 
    geom_line() + 
    geom_point(size=2.5) + 
    ggrepel::geom_label_repel(data=label.data, 
                     aes(label=method), 
                     force=20, 
                     size=2.5,
                     ylim=c(.05,.95),
                     segment.color="black",
                     segment.alpha=.2, 
                     arrow=arrow(length = unit(0,"npc"), type = "open", ends = "last"),
                     box.padding = unit(.45, "lines"),
                     point.padding = unit(.45, 'lines'),
                     nudge_y=-.15, nudge_x=1
    ) +
    geom_point() +
    scale_x_sqrt() +
    cowplot::theme_cowplot() +
    theme(panel.border=element_rect(colour = "grey35", fill = NA),
          panel.grid.major=element_line(colour = "grey85"), 
          axis.line=element_blank(),
          strip.background=element_rect(colour = "grey35"),
          axis.ticks=element_line(colour = "grey35"),
          text = element_text(size=10), axis.text = element_text(size=9)
    ) +
    scale + wrap
}

#--------------------------------------------------------------------------------

#' HM
#'
#' @param df benchmark dataframe
#' @param metric string; score metric to display (has to be a column name in df)
#' @param form function for df reshaping
#' @param fn function; arithmetic operation for the reshaping of numeric df columns 
#' @param cluster.cols binary; if heatmap columns should be clustered
#' @param exclude string vector; define which rows should be dropped
#' @param ... additional options for the ComplexHeatmap::Heatmap function
#'
#' @return a Heatmap class object
#'
HM <- function(df, metric="detPPV", form=method~prop, fn=mean, cluster.cols=FALSE, exclude=NULL, ...){
  
  # reshape data
  m <- reshape2::dcast(df, form=form, value.var=metric, na.rm=TRUE, fun.aggregate=fn)
  
  # generate matrix & annotations
  rownames(m) <- m[,1]
  m[,1] <- NULL
  # get rid of rows with all NaNs
  m <- na.omit(m)
  
  # exclusion of defined rows
  m <- subset(m, !rownames(m) %in% exclude)
  
  # heatmap
  hm <- ComplexHeatmap::Heatmap(m, cluster_columns=cluster.cols, name=metric,
                                rect_gp=grid::gpar(col="white",lwd =.5), ...)
  return(hm)
}


