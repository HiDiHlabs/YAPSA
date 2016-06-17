#' Plot the exposures of a cohort
#'
#' The exposures \code{H}, determined by NMF or by \code{\link{LCD}}, are
#' displayed as a stacked barplot by calling
#' \itemize{
#'  \item \code{\link[ggplot2]{geom_bar}} and optionally
#'  \item \code{\link[ggplot2]{geom_text}}.
#' }
#' The x-axis displays the PIDs (patient identifier or sample), the y-axis
#' the counts attributed to the different signatures with their respective
#' colours per PID. Is called by \code{\link{plot_relative_exposures}}.
#'
#' @param in_exposures_df
#'  Numerical data frame encoding the exposures \code{H}, i.e. which
#'  signature contributes how much to which PID (patient identifier or sample).
#' @param in_signatures_ind_df
#'  A data frame containing meta information about the signatures
#' @param in_subgroups_df
#'  A data frame indicating which PID (patient or sample identifyier) belongs
#'  to which subgroup
#' @param in_sum_ind
#'  Index vector influencing the order in which the PIDs are going to be
#'  displayed
#' @param in_subgroups.field
#'  String indicating the column name in \code{in_subgroups_df} to take the
#'  subgroup information from.
#' @param in_title
#'  Title for the plot to be created.
#' @param in_labels
#'  Flag, if \code{TRUE} the PIDs are displayed on the x-axis
#' @param in_show_subgroups
#'  Flag, if \code{TRUE} then PIDs are grouped by subgroups
#' @param legend_height
#'  How many signatures should be displayed in one column together at most.
#'
#' @return The generated barplot - a ggplot2 plot
#'
#' @examples
#'  data(lymphoma_cohort_LCD_results)
#'  plot_exposures(lymphoma_Nature2013_COSMIC_cutoff_exposures_df,
#'                 chosen_signatures_indices_df,
#'                 COSMIC_subgroups_df)
#'
#' @seealso \code{\link{LCD}}
#' @seealso \code{\link[ggplot2]{geom_bar}}
#' @seealso \code{\link[ggplot2]{geom_text}}
#' @seealso \code{\link{plot_relative_exposures}}
#' 
#' @import ggplot2
#' @import reshape2
#' @export
#' 
plot_exposures <- function(in_exposures_df,
                           in_signatures_ind_df,
                           in_subgroups_df=NULL,
                           in_sum_ind=NULL,
                           in_subgroups.field="subgroup",
                           in_title="",
                           in_labels=TRUE,
                           in_show_subgroups=TRUE,
                           legend_height=10) {
  .e <- environment()
  temp_exposures_df <- in_exposures_df
  number_of_sigs <- dim(in_exposures_df)[1]
  number_of_PIDs <- dim(in_exposures_df)[2]
  ## prepare for output
  temp_exposures_df$sig_index <- 
    in_signatures_ind_df$index[match(rownames(temp_exposures_df),
                                     in_signatures_ind_df$sig)]
  if(!is.null(in_subgroups_df)){
    subgroups_ind <- which(names(in_subgroups_df)==in_subgroups.field)
    my_sum_ind <- order(in_subgroups_df[,subgroups_ind],
                        in_subgroups_df$compl_sum) 
    my_PIDs_vector <- in_subgroups_df$PID
  } else {
    total_counts <- colSums(in_exposures_df)
    my_sum_ind <- order(total_counts)
    my_PIDs_vector <- names(in_exposures_df)
  }
  if(!is.null(in_sum_ind)) {
    my_sum_ind <- in_sum_ind
  }
  order_sum_ind <- order(my_sum_ind)
  names(order_sum_ind) <- names(in_exposures_df)
  ## visualize as stacked barplot with ggplot2 for comparison with matlab NMF output
  exposures_df_melt <- melt(temp_exposures_df,id.vars="sig_index")
  names(exposures_df_melt)[2] <- "PID"
  exposures_df_melt <- exposures_df_melt[order(exposures_df_melt$PID,
                                               exposures_df_melt$sig_index),]
  exposures_df_melt$PID_index <- order_sum_ind[exposures_df_melt$PID]
  exposures_df_melt$PID_index <- as.factor(exposures_df_melt$PID_index)
  exposures_df_melt$sig_index <- as.factor(exposures_df_melt$sig_index)
  number_of_legend_cols <- 
    ifelse((dim(in_signatures_ind_df)[1] %% legend_height)==0,0,1) +
    (dim(in_signatures_ind_df)[1] %/% legend_height)
  p0 <- ggplot(environment = .e) + 
    geom_bar(data=exposures_df_melt, 
             aes_string(x="PID_index",y="value",fill="sig_index"),
             stat='identity') +
    scale_fill_manual(labels=in_signatures_ind_df$sig,
                      values=in_signatures_ind_df$colour) + 
    scale_x_discrete("PID",breaks=seq(1,number_of_PIDs),
                     labels=my_PIDs_vector[my_sum_ind])
  p <- p0 + theme(axis.text.x=element_text(angle=90,size=7,vjust=0.5)) +
    guides(fill=guide_legend(ncol=number_of_legend_cols))
  q <- p0 + theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    labs(x="", y="")
  if(!is.null(in_subgroups_df) & in_show_subgroups){
    pid_sums <- colSums(in_exposures_df)
    pid_ind <- match(my_PIDs_vector,names(pid_sums))
    pid_sums <- pid_sums[pid_ind]
    upper_y_lim <- 1.1*max(pid_sums)
    in_subgroups_df$pid_sums <- pid_sums
    in_subgroups_df <- in_subgroups_df[my_sum_ind,]
    in_subgroups_df$my_x <- seq(1,dim(in_subgroups_df)[1])
    p1 <- p +
      geom_text(data=in_subgroups_df,
                aes_string(label=in_subgroups.field,x="my_x",y="pid_sums"),
                angle=90,size=3,hjust=0) + ylim(0,upper_y_lim) +
      labs(title=in_title, y="")    
  } else {
    p1 <- p +
      labs(y="")
  }
  if(!in_labels) {p1 <- q + theme(legend.position = "none")}
  return(p1)
}



plot_exposures_old <- function(in_exposures_df,
                           in_signatures_ind_df,
                           in_subgroups_df,
                           in_sum_ind=NULL,
                           in_subgroups.field="subgroup",
                           in_title="",
                           in_labels=TRUE,
                           in_show_subgroups=TRUE) {
  .e <- environment()
  subgroups_ind <- which(names(in_subgroups_df)==in_subgroups.field)
  if(!is.null(in_sum_ind)) {
    my_sum_ind <- in_sum_ind
  } else {
    my_sum_ind <- order(in_subgroups_df[,subgroups_ind],
                        in_subgroups_df$compl_sum)
  }
  order_sum_ind <- order(my_sum_ind)
  legend_height <- 10
  ## prepare for output
  temp_exposures_df <- in_exposures_df
  temp_exposures_df$sig_index <- 
    in_signatures_ind_df$index[match(rownames(temp_exposures_df),
                                     in_signatures_ind_df$sig)]
  number_of_sigs <- dim(in_exposures_df)[1]
  ## visualize as stacked barplot with ggplot2 for comparison with matlab NMF output
  exposures_df_melt <- melt(temp_exposures_df,id.vars="sig_index")
  names(exposures_df_melt)[2] <- "PID"
  exposures_df_melt <- exposures_df_melt[order(exposures_df_melt$PID,
                                               exposures_df_melt$sig_index),]
#   exposures_df_melt$PID_index <- 
#     in_subgroups_df$index[match(exposures_df_melt$PID,in_subgroups_df$PID)]
  exposures_df_melt$PID_index <- order_sum_ind[match(exposures_df_melt$PID,
                                                     in_subgroups_df$PID)]
  exposures_df_melt$PID_index <- as.factor(exposures_df_melt$PID_index)
  exposures_df_melt$sig_index <- as.factor(exposures_df_melt$sig_index)
  number_of_legend_cols <- 
    ifelse((dim(in_signatures_ind_df)[1] %% legend_height)==0,0,1) +
    (dim(in_signatures_ind_df)[1] %/% legend_height)
  p <- ggplot(environment = .e) + 
    geom_bar(data=exposures_df_melt, 
             aes_string(x="PID_index",y="value",fill="sig_index"),
             stat='identity') +
    scale_fill_manual(labels=in_signatures_ind_df$sig,
                      values=in_signatures_ind_df$colour) + 
    scale_x_discrete("PID",breaks=seq(1,dim(in_subgroups_df)[1]),
                     labels=in_subgroups_df$PID[my_sum_ind]) +
    theme(axis.text.x=element_text(angle=90,size=7,vjust=0.5)) +
    guides(fill=guide_legend(ncol=number_of_legend_cols))
  q <- ggplot(environment = .e) + 
    geom_bar(data=exposures_df_melt,
             aes_string(x="PID_index",y="value",fill="sig_index"),
             stat='identity') +
    scale_fill_manual(labels=in_signatures_ind_df$sig,
                      values=in_signatures_ind_df$colour) + 
    scale_x_discrete(breaks=seq(1,dim(in_subgroups_df)[1]),label=NULL) +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    labs(x="", y="")
  pid_sums <- colSums(in_exposures_df)
  pid_ind <- match(in_subgroups_df$PID,names(pid_sums))
  pid_sums <- pid_sums[pid_ind]
  upper_y_lim <- 1.1*max(pid_sums)
  in_subgroups_df$pid_sums <- pid_sums
  in_subgroups_df <- in_subgroups_df[my_sum_ind,]
  in_subgroups_df$my_x <- seq(1,dim(in_subgroups_df)[1])
  if(in_show_subgroups) {
    p1 <- p +
      geom_text(data=in_subgroups_df,
                aes_string(label=in_subgroups.field,x="my_x",y="pid_sums"),
                angle=90,size=3,hjust=0) + ylim(0,upper_y_lim) +
      labs(title=in_title, y="")
  } else {
    p1 <- p +
      labs(y="")
  }
  if(!in_labels) {p1 <- q + theme(legend.position = "none")}
  return(p1)
}


#' Plot the relative exposures of a cohort
#'
#' Plot the relative or normalized exposures of a cohort. This function first
#' normalizes its input and then sends the normalized data to
#' \code{\link{plot_relative_exposures}}.
#'
#' @param in_exposures_df
#'  Numerical data frame encoding the exposures \code{H}, i.e. which
#'  signature contributes how much to which PID (patient identifier or sample).
#' @param in_signatures_ind_df
#'  A data frame containing meta information about the signatures
#' @param in_subgroups_df
#'  A data frame indicating which PID (patient or sample identifyier) belongs
#'  to which subgroup
#' @param in_sum_ind
#'  Index vector influencing the order in which the PIDs are going to be
#'  displayed
#' @param in_subgroups.field
#'  String indicating the column name in \code{in_subgroups_df} to take the
#'  subgroup information from.
#' @param in_title
#'  Title for the plot to be created.
#' @param in_labels
#'  Flag, if \code{TRUE} the PIDs are displayed on the x-axis
#' @param in_show_subgroups
#'  Flag, if \code{TRUE} then PIDs are grouped by subgroups
#'
#' @return The generated barplot - a ggplot2 plot
#'
#' @examples
#'  data(lymphoma_cohort_LCD_results)
#'  plot_relative_exposures(lymphoma_Nature2013_COSMIC_cutoff_exposures_df,
#'                          chosen_signatures_indices_df,
#'                          COSMIC_subgroups_df)
#'
#' @seealso \code{\link{plot_exposures}}
#' 
#' @export
#' 
plot_relative_exposures <- function(in_exposures_df,in_signatures_ind_df,in_subgroups_df,in_sum_ind=NULL,
                                    in_subgroups.field="subgroup",in_title="",in_labels=TRUE,
                                    in_show_subgroups=TRUE) {
  t_temp_temp_exposures_df <- as.data.frame(t(in_exposures_df))
  t_rowSums <- rowSums(t_temp_temp_exposures_df)
  t_rowSums[which(t_rowSums==0)] <- 1e06
  t_temp_exposures_df <- t_temp_temp_exposures_df/t_rowSums
  this_exposures_df <- as.data.frame(t(t_temp_exposures_df))
  return(plot_exposures(this_exposures_df,in_signatures_ind_df,in_subgroups_df,in_sum_ind,in_subgroups.field,in_title,in_labels,in_show_subgroups))
}


#' Plot the exposures of a cohort as a ComplexHeatmap
#'
#' The exposures \code{H}, determined by NMF or by \code{\link{LCD}}, are
#' displayed as a stacked barplot by calling \code{\link[ComplexHeatmap]{Heatmap}}.
#' The x-axis displays the PIDs (patient identifier or sample), the y-axis
#' the counts attributed to the different signatures with their respective
#' colours per PID. It is analogous to \code{\link{plot_exposures}}.
#'
#' @param mat
#'  Numerical data frame encoding the exposures \code{H}, i.e. which
#'  signature contributes how much to which PID (patient identifier or sample).
#' @param col
#'  A data frame containing meta information about the signatures
#' @param anno
#'  A data frame indicating which PID (patient or sample identifyier) belongs
#'  to which subgroup
#' @param anno_color
#'  Index vector influencing the order in which the PIDs are going to be
#'  displayed
#' @param ylab
#'  String indicating the column name in \code{in_subgroups_df} to take the
#'  subgroup information from.
#' @param title
#'  Title for the plot to be created.
#' @param in_labels
#'  Whether or not to show the names of the samples.
#' @param in_barplot_borders
#'  Whether or not to show border lines in barplot
#' @param in_column_anno_borders
#'  Whether or not to draw separating lines between the fields in the annotation
#'
#' @details
#'  It might be necessary to install the newest version of the development branch
#'  of the packages \pkg{circlize} and \pkg{ComplexHeatmap} by Zuguang Gu:
#'  \code{devtools::install_github("jokergoo/circlize")}
#'  \code{devtools::install_github("jokergoo/ComplexHeatmap")}
#'
#' @return The function doesn't return any value.
#' 
#' @examples
#'  data(lymphoma_cohort_LCD_results)
#'  temp_anno_color <- aggregate(col~subgroup,data=COSMIC_subgroups_df,function(l) return(l[1]))
#'  my_anno_color <- temp_anno_color$col
#'  names(my_anno_color) <- temp_anno_color$subgroup
#'  enhanced_barplot(lymphoma_Nature2013_COSMIC_cutoff_exposures_df[,order(COSMIC_subgroups_df$index)],
#'                   chosen_signatures_indices_df$colour,
#'                   COSMIC_subgroups_df$subgroup[order(COSMIC_subgroups_df$index)],
#'                   my_anno_color)  
#'  
#' @seealso \code{\link[ComplexHeatmap]{HeatmapAnnotation}}
#' @seealso \code{\link[ComplexHeatmap]{Heatmap}}
#' @seealso \code{\link[ComplexHeatmap]{decorate_heatmap_body}}
#' @seealso \code{\link{plot_exposures}}
#' 
#' @import ComplexHeatmap
#' @export
#' 
enhanced_barplot = function(mat, col, anno, anno_color, ylab = NULL,
                            title = "",in_labels=TRUE,
                            in_barplot_borders=TRUE,
                            in_column_anno_borders=FALSE) {
  names(col) = rownames(mat)
  mat_foo = matrix(rep(rownames(mat), ncol(mat)), nrow = nrow(mat))
  
  if(in_labels){
    if(!in_column_anno_borders){
      ha = HeatmapAnnotation(df = data.frame(anno = anno), col = list(anno = anno_color),
                             text = anno_text(colnames(mat), offset = unit(1, "npc"), just = "right", rot = 90),
                             annotation_height = unit.c(unit(5, "mm"), max_text_width(colnames(mat)) ))      
    } else {
      ha = HeatmapAnnotation(df = data.frame(anno = anno), col = list(anno = anno_color),
                             text = anno_text(colnames(mat), offset = unit(1, "npc"), just = "right", rot = 90),
                             annotation_height = unit.c(unit(5, "mm"), max_text_width(colnames(mat)) ),
                             gp = gpar(col="black"))
    }
  } else {
    if(!in_column_anno_borders){
      ha = HeatmapAnnotation(df = data.frame(anno = anno), col = list(anno = anno_color),
                             annotation_height = unit.c(unit(5, "mm"), max_text_width(colnames(mat)) ))    
    } else {
      ha = HeatmapAnnotation(df = data.frame(anno = anno), col = list(anno = anno_color),
                             annotation_height = unit.c(unit(5, "mm"), max_text_width(colnames(mat)) ),
                             gp = gpar(col="black"))          
    }
  }
  ht = Heatmap(mat_foo, col = col, name = "main", cluster_rows = FALSE, cluster_columns = FALSE,
               rect_gp = gpar(type = "none"), bottom_annotation = ha,
               heatmap_legend_param = list(at = rev(names(col)), labels = rev(names(col))))
  
  draw(ht, padding = unit(c(2, 20, 2, 2), "mm"), column_title = title)
  
  decorate_heatmap_body("main", {
    pushViewport(viewport(xscale = c(0, ncol(mat)),
                          yscale = c(0, max(colSums(mat)))))
    if(in_barplot_borders) {
      for(i in seq_len(ncol(mat))) {
        x = i - 0.5
        for(j in seq_len(nrow(mat))) {
          y = sum(mat[1:j, i])
          grid.rect(x, y, width = 0.9, height = mat[j, i], just = c("top"),
                    default.units = "native", gp = gpar(fill = col[j]))
        }
      }
    } else {
      for(i in seq_len(ncol(mat))) {
        x = i - 0.5
        for(j in seq_len(nrow(mat))) {
          y = sum(mat[1:j, i])
          grid.rect(x, y, width = 0.9, height = mat[j, i], just = c("top"),
                    default.units = "native", gp = gpar(fill = col[j],
                                                        col = NA))
        }
      }      
    }
    grid.yaxis(gp=gpar(fontsize=8))
    if(!is.null(ylab)){
      grid.text(ylab, x = unit(0, "npc") - unit(15, "mm"), just = "left", rot = 90)      
    }    
    
    upViewport()
  })
}


#' Wrapper for enhanced_barplot
#' 
#' @param in_exposures_df
#'  Numerical data frame encoding the exposures \code{H}, i.e. which
#'  signature contributes how much to which PID (patient identifier or sample).
#' @param in_signatures_ind_df
#'  A data frame containing meta information about the signatures. If NULL, the
#'  colour information for the signatures is taken from a rainbow palette.
#' @param in_subgroups_df
#'  A data frame indicating which PID (patient or sample identifyier) belongs
#'  to which subgroup. If NULL, it is assumed that all PIDs belong to one common
#'  subgroup. The colour coding for the default subgroup is red.
#' @param in_sum_ind
#'  Index vector influencing the order in which the PIDs are going to be
#'  displayed
#' @param in_subgroups.field
#'  String indicating the column name in \code{in_subgroups_df} to take the
#'  subgroup information from.
#' @param in_title
#'  Title for the plot to be created.
#' @param in_labels
#'  Flag, if \code{TRUE} the PIDs are displayed on the x-axis
#' @param in_show_subgroups
#'  Flag, if \code{TRUE} then PIDs are grouped by subgroups
#' @param ylab
#'  Label of the y-axis on the plot to be generate
#' @param in_barplot_borders
#'  Whether or not to show border lines in barplot
#' @param in_column_anno_borders
#'  Whether or not to draw separating lines between the fields in the annotation
#'
#' @return The generated barplot - a ggplot2 plot
#'
#' @examples
#'  data(lymphoma_cohort_LCD_results)
#'  exposures_barplot(lymphoma_Nature2013_COSMIC_cutoff_exposures_df,
#'                    chosen_signatures_indices_df,
#'                    COSMIC_subgroups_df)  
#'  
#' @seealso \code{\link{enhanced_barplot}}
#' 
#' @export
#' 
exposures_barplot <- function(in_exposures_df,in_signatures_ind_df=NULL,
                              in_subgroups_df=NULL,in_sum_ind=NULL,
                              in_subgroups.field="subgroup",in_title="",
                              in_labels=TRUE,in_show_subgroups=TRUE,ylab=NULL,
                              in_barplot_borders=TRUE,
                              in_column_anno_borders=FALSE) {
  if(!is.null(in_subgroups_df)){
    temp_anno_color <- aggregate(col~subgroup,data=in_subgroups_df,
                                 function(l) return(l[1]))
    my_anno_color <- temp_anno_color$col
    names(my_anno_color) <- temp_anno_color$subgroup
    order_ind <- order(in_subgroups_df$index)
    if(!is.null(in_sum_ind)) order_ind <- in_sum_ind
    subgroup_vector <- in_subgroups_df$subgroup[order_ind]
  } else {
    total_counts <- colSums(in_exposures_df)
    order_ind <- rev(order(total_counts))
    if(!is.null(in_sum_ind)) order_ind <- in_sum_ind
    my_anno_color <- c("#FF0000FF")
    names(my_anno_color) <- "subgroup"
    subgroup_vector <- factor(rep("subgroup",dim(in_exposures_df)[2]))
  }
  if(!is.null(in_signatures_ind_df)){
    signature_colour_vector <- in_signatures_ind_df$colour
  } else {
    signature_colour_vector <- rainbow(dim(in_exposures_df)[1])
  }
  title  <- in_title
  enhanced_barplot(in_exposures_df[,order_ind],
                   signature_colour_vector,
                   subgroup_vector,
                   my_anno_color,in_labels=in_labels,
                   in_barplot_borders=in_barplot_borders,
                   in_column_anno_borders=in_column_anno_borders)
}



#' Plot the exposures of a cohort with different layers of annotation
#'
#' The exposures \code{H}, determined by NMF or by \code{\link{LCD}}, are
#' displayed as a stacked barplot by calling \code{\link[ComplexHeatmap]{Heatmap}}.
#' The x-axis displays the PIDs (patient identifier or sample), the y-axis
#' the counts attributed to the different signatures with their respective
#' colours per PID. It is analogous to \code{\link{plot_exposures}}. As many
#' layers of information as desired can be added via an annotation data frame.
#' The annotation data is handled in a way similar to
#' \code{\link{annotation_heatmap_exposures}}. This function calls: 
#' \itemize{
#'  \item \code{\link[ComplexHeatmap]{rowAnnotation}},
#'  \item \code{\link[ComplexHeatmap]{HeatmapAnnotation}} and
#'  \item \code{\link[ComplexHeatmap]{Heatmap}}
#' }
#'
#' @details
#'  It might be necessary to install the newest version of the development branch
#'  of the packages \pkg{circlize} and \pkg{ComplexHeatmap} by Zuguang Gu:
#'  \code{devtools::install_github("jokergoo/circlize")}
#'  \code{devtools::install_github("jokergoo/ComplexHeatmap")}
#'
#' @param in_exposures_df
#'  Numerical data frame encoding the exposures \code{H}, i.e. which
#'  signature contributes how much to which PID (patient identifier or sample).
#' @param in_signatures_ind_df
#'  A data frame containing meta information about the signatures
#' @param in_subgroups_df
#'  A data frame indicating which PID (patient or sample identifyier) belongs
#'  to which subgroup
#' @param in_annotation_df
#'  A data frame indicating which PID (patient or sample identifyier) belongs
#'  to which subgroup for all layers of annotation
#' @param in_annotation_col
#'  A list indicating colour attributions for all layers of annotation
#' @param in_signatures_ind_df
#'  A data frame containing meta information about the signatures, especially the
#'  asserted colour
#' @param ylab
#'  String indicating the column name in \code{in_subgroups_df} to take the
#'  subgroup information from.
#' @param title
#'  Title for the plot to be created.
#' @param in_labels
#'  Whether or not to show the names of the samples.
#' @param in_barplot_borders
#'  Whether or not to show border lines in barplot
#' @param in_column_anno_borders
#'  Whether or not to draw separating lines between the fields in the annotation
#' @param in_annotation_legend_side
#'  Where to put the legends of the annotation df, default is right.
#'
#' @details
#'  It might be necessary to install the newest version of the development branch
#'  of the packages \pkg{circlize} and \pkg{ComplexHeatmap} by Zuguang Gu:
#'  \code{devtools::install_github("jokergoo/circlize")}
#'  \code{devtools::install_github("jokergoo/ComplexHeatmap")}
#'
#' @return The function doesn't return any value.
#' 
#' @examples
#'  NULL  
#'  
#' @seealso \code{\link[ComplexHeatmap]{HeatmapAnnotation}}
#' @seealso \code{\link[ComplexHeatmap]{Heatmap}}
#' @seealso \code{\link[ComplexHeatmap]{decorate_heatmap_body}}
#' @seealso \code{\link{annotation_heatmap_exposures}}
#' @seealso \code{\link{plot_exposures}}
#' 
#' @import ComplexHeatmap
#' @import circlize
#' @export
#' 
annotation_exposures_barplot <- function(in_exposures_df,
                                         in_signatures_ind_df,
                                         in_subgroups_df,
                                         in_annotation_df,
                                         in_annotation_col,
                                         ylab = NULL,
                                         title = "",in_labels=FALSE,
                                         in_barplot_borders=TRUE,
                                         in_column_anno_borders=FALSE,
                                         in_annotation_legend_side="right") {
  order_ind <- order(in_subgroups_df$index)
  mat <- in_exposures_df[,order_ind]
  this_annotation_df <- in_annotation_df[order_ind,]
  col <- in_signatures_ind_df$colour
  
  names(col) = rownames(mat)
  mat_foo = matrix(rep(rownames(mat), ncol(mat)), nrow = nrow(mat))
  
  if(in_labels){
    if(!in_column_anno_borders){
      ha = HeatmapAnnotation(df = this_annotation_df,
                             col = in_annotation_col,
                             text = anno_text(colnames(mat), offset = unit(1, "npc"), just = "right", rot = 90))      
    } else {
      ha = HeatmapAnnotation(df = this_annotation_df,
                             col = in_annotation_col,
                             text = anno_text(colnames(mat), offset = unit(1, "npc"), just = "right", rot = 90),
                             gp = gpar(col="black"))
    }
  } else {
    if(!in_column_anno_borders){
      ha = HeatmapAnnotation(df = this_annotation_df,
                             col = in_annotation_col)    
    } else {
      ha = HeatmapAnnotation(df = this_annotation_df,
                             col = in_annotation_col,
                             gp = gpar(col="black"))          
    }
  }
  ht = Heatmap(mat_foo, col = col, name = "main", cluster_rows = FALSE, cluster_columns = FALSE,
               rect_gp = gpar(type = "none"), bottom_annotation = ha,
               heatmap_legend_param = list(at = rev(names(col)), labels = rev(names(col))))
  
  draw(ht, padding = unit(c(2, 20, 2, 2), "mm"), column_title = title,
       annotation_legend_side=in_annotation_legend_side)
  
  decorate_heatmap_body("main", {
    pushViewport(viewport(xscale = c(0, ncol(mat)),
                          yscale = c(0, max(colSums(mat)))))
    if(in_barplot_borders) {
      for(i in seq_len(ncol(mat))) {
        x = i - 0.5
        for(j in seq_len(nrow(mat))) {
          y = sum(mat[1:j, i])
          grid.rect(x, y, width = 0.9, height = mat[j, i], just = c("top"),
                    default.units = "native", gp = gpar(fill = col[j]))
        }
      }
    } else {
      for(i in seq_len(ncol(mat))) {
        x = i - 0.5
        for(j in seq_len(nrow(mat))) {
          y = sum(mat[1:j, i])
          grid.rect(x, y, width = 0.9, height = mat[j, i], just = c("top"),
                    default.units = "native", gp = gpar(fill = col[j],
                                                        col = NA))
        }
      }      
    }
    grid.yaxis(gp=gpar(fontsize=8))
    if(!is.null(ylab)){
      grid.text(ylab, x = unit(0, "npc") - unit(15, "mm"), just = "left", rot = 90)      
    }    
    
    upViewport()
  })
}


#' Make an annotation plot for subgroup information
#' 
#' For cohort-wide plots like those produced by \code{\link{plot_exposures}}
#' or \code{\link{plot_relative_exposures}}, annotation plots for per-PID
#' information can be put on top. This function produces an annotation plot
#' for subgroup inforamtion.
#' 
#' @param in_subgroup_df
#'  Data frame carrying the subgroup information
#' @param mycol
#'  Vector containing the colour information for the subgroups
#' @param in_subgroup.field
#'  Indicates the name of the column in which the subgroup information
#'  is encoded
#' @param in_index.field
#'  Indicates the name of the column in which the index information
#'  is encoded
#' @param in_size
#'  Numerical value indicating the size unit (important when uniting different
#'  plots)
#' 
#' @return The plot as a ggplot2.
#' 
#' @examples
#'  data(lymphoma_cohort_LCD_results)
#'  temp_anno_color <- aggregate(col~subgroup,data=COSMIC_subgroups_df,function(l) return(l[1]))
#'  plot_colour_legend_subgroups(COSMIC_subgroups_df,temp_anno_color$col)
#'
#' @seealso \code{\link{plot_exposures}}
#' @seealso \code{\link{plot_relative_exposures}}
#' 
#' @import ggplot2
#' @export
#'
plot_colour_legend_subgroups <- function(in_subgroup_df,mycol,
                                         in_subgroup.field="subgroup",
                                         in_index.field="index",
                                         in_size=0.9) {
  .e <- environment()
  p2 <- ggplot(in_subgroup_df, aes_string(x=in_index.field, y=1, fill=in_subgroup.field,
                                          width = in_size, height = in_size), environment = .e) + 
    geom_tile() + 
    theme_bw(base_size = 14) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_blank()) +
    scale_fill_manual(values = mycol, 
                      guide = guide_legend(title = '', keywidth = 0.7, keyheight = 0.7, nrow = 1)) + 
    theme(legend.position = "top") +
    xlab('') +
    ylab('') + 
    coord_fixed(ratio=1)
  return(p2)
}



#' Plot results of the Stratification of a Mutational Catalogue
#'
#' Plot a big composite figure with 3 columns: in the left column the per-PID
#' absolute exposures will be shown, in the middle column the per_PID
#' relative or normalized exposures will be shown, in the right column the
#' cohort-wide exposures are shown (averaged over PIDs).
#'
#' @param number_of_strata
#'  Number of strata as deduced from \code{link{SMC}}
#' @param output_path
#'  Path to file where the results are going to be stored. If NULL, the results will
#'  be plotted to the running environment.
#' @param decomposition_method
#'  String for the filename of the generated barplot.
#' @param number_of_sigs
#'  Number of signatures
#' @param name_list
#'  Names of the contructed strata.
#' @param exposures_strata_list
#'  The list of \code{s} strata specific exposures Hi, all are numerical data 
#'  frames with \code{l} rows and \code{m} columns, \code{l} being the number
#'  of signatures and \code{m} being the number of samples
#' @param this_signatures_ind_df
#'  A data frame containing meta information about the signatures
#' @param this_subgroups_df
#'  A data frame indicating which PID (patient or sample identifyier) belongs
#'  to which subgroup
#' @param in_strata_order_ind
#'  Index vector defining reordering of the strata
#' @param exposures_both_rel_df_list
#'  A list of \code{s} strata specific cohortwide (i.e. averaged over cohort)
#'  normalized exposures
#' @param cohort_method_flag
#'  Either or several of \code{c("all_PIDs","cohort","norm_PIDs")}, representing
#'  alternative ways to average over the cohort.
#'  
#' @examples
#' NULL
#'  
#' @return The function doesn't return any value.
#' 
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @export
#' 
plot_SMC_old <- function(number_of_strata,output_path,decomposition_method,number_of_sigs,name_list,
                         exposures_strata_list,this_signatures_ind_df,this_subgroups_df,
                         in_strata_order_ind,exposures_both_rel_df_list,cohort_method_flag) {
  plot_list_cohort <- list()
  all_method_colour_vector <- c("grey60","grey40","grey20")
  names(all_method_colour_vector) <- c("all_PIDs","cohort","norm_PIDs")
  if(length(cohort_method_flag)<2) {
    for (i in seq_len(length(exposures_both_rel_df_list))) { 
      temp_df <- exposures_both_rel_df_list[[i]]
      temp_df$sig_index <- this_signatures_ind_df$index[match(temp_df$sig,this_signatures_ind_df$sig)]
      temp_df$sig_index <- as.factor(temp_df$sig_index)
      plot_list_cohort[[i]] <- ggplot() + 
        ggplot2::geom_bar(data=temp_df,aes_string(x="sig_index",y="exposure",fill="sig_index",
                                                  size=0.3),
                          stat='identity',position="dodge",width=.7) + 
        scale_fill_manual(name="sig",labels=this_signatures_ind_df$sig,values=this_signatures_ind_df$colour) +
        geom_errorbar(data=temp_df,aes_string(x="sig_index",ymin="exposure_min",ymax="exposure_max"),width=0.5) +
        labs(x="",y="") +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
              legend.position = "none")
    }    
  } else {
    method_colour_vector <- all_method_colour_vector[cohort_method_flag]
    for (i in seq_len(length(exposures_both_rel_df_list))) { 
      temp_df <- exposures_both_rel_df_list[[i]]
      temp_df$sig_index <- this_signatures_ind_df$index[match(temp_df$sig,this_signatures_ind_df$sig)]
      temp_df$sig_index <- as.factor(temp_df$sig_index)
      plot_list_cohort[[i]] <- ggplot() + 
        ggplot2::geom_bar(data=temp_df,aes_string(x="sig_index",y="exposure",colour="method",
                                                  fill="sig_index",size=0.3),
                          stat='identity',position="dodge",width=.7) + 
        scale_fill_manual(name="sig",labels=this_signatures_ind_df$sig,values=this_signatures_ind_df$colour) +
        scale_colour_manual(values=method_colour_vector) +
        geom_errorbar(data=temp_df,aes_string(x="sig_index",ymin="exposure_min",ymax="exposure_max"),width=0.5) +
        labs(x="",y="") +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              legend.position = "none")
    }    
  }
  all_cohort_plot <- ggplotGrob(plot_list_cohort[[1]] + theme(legend.position="right") + 
                                  guides(size=FALSE))$grobs
  common_legend <- all_cohort_plot[[which(sapply(all_cohort_plot, function(x) x$name) == "guide-box")]]
  
  plot_list_abs <- list()
  plot_list_rel <- list()
  plot_list_abs[[1]] <- plot_exposures(exposures_strata_list$exposures_all_df,
                                       this_signatures_ind_df,this_subgroups_df,
                                       in_title="new",in_labels=FALSE)
  plot_list_rel[[1]] <- plot_relative_exposures(exposures_strata_list$exposures_all_df,
                                                this_signatures_ind_df,this_subgroups_df,
                                                in_title="new",in_labels=FALSE)
  for (i in seq_len(number_of_strata)) {
    plot_list_abs[[i+1]] <- plot_exposures(exposures_strata_list$sub_exposures_list[[i]],
                                           this_signatures_ind_df,this_subgroups_df,
                                           in_title="new",in_labels=FALSE)
    plot_list_rel[[i+1]] <- plot_relative_exposures(exposures_strata_list$sub_exposures_list[[i]],
                                                    this_signatures_ind_df,this_subgroups_df,
                                                    in_title="new",in_labels=FALSE)
  }
  strata_order_ind <- c(1,in_strata_order_ind+1)
  plot_list_abs <- plot_list_abs[strata_order_ind]
  plot_list_rel <- plot_list_rel[strata_order_ind]
  plot_list_cohort <- plot_list_cohort[strata_order_ind]
  plot_list_all <- list()
  for (i in seq_len(number_of_strata+1)) {
    plot_list_all[[2*(i-1)+1]] <- plot_list_abs[[i]]
    plot_list_all[[2*(i-1)+2]] <- plot_list_rel[[i]]
  }
  
  temp_name_list <- c("all",name_list)
  name_df <- data.frame(t(data.frame(t(temp_name_list[strata_order_ind]))))
  names(name_df) <- "names"
  name_df$names <- as.character(name_df$names)
  name_df$x_pos <- 0
  name_df$y_pos <- 0
  name_df$new_names <- gsub("non_","non ",name_df$names)
  name_df$new_names <- gsub("_","\n",name_df$new_names)
  
  horizontal_label_width <- 1
  horizontal_big_element_width <- 10
  horizontal_small_element_width <- 4
  horizontal_legend_width <- 2
  number_of_horizontal_units <- horizontal_label_width + 2*horizontal_big_element_width + horizontal_small_element_width + horizontal_legend_width
  vertical_legend_height <- 1
  vertical_element_height <- 5
  number_of_vertical_units <- vertical_legend_height + (number_of_strata+1)*vertical_element_height
  horizontal_figure_factor <- 60
  vertical_figure_factor <- 40
  
  if(!is.null(output_path)){
    png(output_path,
        width=number_of_horizontal_units*horizontal_figure_factor,
        height=number_of_vertical_units*vertical_figure_factor)    
  }
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(number_of_vertical_units, number_of_horizontal_units)))
  vertical_temp_stop <- vertical_legend_height
  for (i in seq_len(number_of_strata+1)) {
    horizontal_temp_start <- horizontal_label_width + 1
    horizontal_temp_stop <- (horizontal_temp_start - 1) + horizontal_big_element_width
    vertical_temp_start <- vertical_temp_stop + 1
    vertical_temp_stop <- (vertical_temp_start - 1) + vertical_element_height    
    print(plot_list_abs[[i]], vp = vplayout(vertical_temp_start:vertical_temp_stop, horizontal_temp_start:horizontal_temp_stop))
    p <- ggplot() + geom_text(data=name_df[i,],aes_string(label="new_names",x="x_pos",y="y_pos"),angle=90,size=3) +
      theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_blank()) +
      labs(x="", y="")  
    print(p, vp = vplayout(vertical_temp_start:vertical_temp_stop, 1:horizontal_label_width))
    horizontal_temp_start <- horizontal_temp_stop + 1
    horizontal_temp_stop <- (horizontal_temp_start - 1) + horizontal_big_element_width
    print(plot_list_rel[[i]], vp = vplayout(vertical_temp_start:vertical_temp_stop, horizontal_temp_start:horizontal_temp_stop))
    horizontal_temp_start <- horizontal_temp_stop + 1
    horizontal_temp_stop <- (horizontal_temp_start - 1) + horizontal_small_element_width
    print(plot_list_cohort[[i]], vp = vplayout(vertical_temp_start:vertical_temp_stop, horizontal_temp_start:horizontal_temp_stop))
  }
  legend_vp = vplayout(1:number_of_vertical_units, (horizontal_temp_stop + 1):(horizontal_temp_stop + horizontal_legend_width))
  pushViewport(legend_vp)
  grid.draw(common_legend)
  if(!is.null(output_path)){
    dev.off()
  }
  return()
}


plot_SMC_PID_facet <- function(in_abs_df_list,in_rel_df_list,
                               in_signatures_ind_df,in_subgroups_df,
                               in_name_list,in_sum_ind=NULL,
                               in_subgroups.field="subgroup",
                               in_strata_order_ind=seq_len(length(in_name_list)),
                               in_label_orientation="turn"){
  .e <- environment()
  legend_height <- 10
  subgroups_ind <- which(names(in_subgroups_df)==in_subgroups.field)
  if(!is.null(in_sum_ind)) {
    my_sum_ind <- in_sum_ind
  } else {
    my_sum_ind <- order(in_subgroups_df[,subgroups_ind],in_subgroups_df$compl_sum)
  }
  order_sum_ind <- order(my_sum_ind)
  number_of_sigs <- dim(in_abs_df_list[[1]])[1]
  abs_melt_df_list <- lapply(in_abs_df_list, FUN=function(x) {
    x$sig_index <- in_signatures_ind_df$index[match(rownames(x),in_signatures_ind_df$sig)]
    exposures_df_melt <- melt(x,id.vars="sig_index")
    names(exposures_df_melt)[2] <- "PID"
    exposures_df_melt <- exposures_df_melt[order(exposures_df_melt$PID,exposures_df_melt$sig_index),]
    exposures_df_melt$PID_index <- as.factor(
      order_sum_ind[match(exposures_df_melt$PID,in_subgroups_df$PID)])
    exposures_df_melt$sig_index <- as.factor(exposures_df_melt$sig_index) 
    return(exposures_df_melt)
  })
  for(i in seq_len(length(abs_melt_df_list))){
    abs_melt_df_list[[i]]$stratum <- in_name_list[[i]]
  }
  abs_melt_df_list <- abs_melt_df_list[in_strata_order_ind]
  abs_melt_df <- do.call("rbind",abs_melt_df_list)
  abs_melt_df$stratum <- factor(abs_melt_df$stratum,levels=unique(abs_melt_df$stratum))
  rel_melt_df_list <- lapply(in_rel_df_list, FUN=function(x) {
    x$sig_index <- in_signatures_ind_df$index[match(rownames(x),in_signatures_ind_df$sig)]
    exposures_df_melt <- melt(x,id.vars="sig_index")
    names(exposures_df_melt)[2] <- "PID"
    exposures_df_melt <- exposures_df_melt[order(exposures_df_melt$PID,exposures_df_melt$sig_index),]
    exposures_df_melt$PID_index <- as.factor(
      order_sum_ind[match(exposures_df_melt$PID,in_subgroups_df$PID)])
    exposures_df_melt$sig_index <- as.factor(exposures_df_melt$sig_index) 
    return(exposures_df_melt)
  })
  for(i in seq_len(length(rel_melt_df_list))){
    rel_melt_df_list[[i]]$stratum <- in_name_list[[i]]
  }
  rel_melt_df_list <- rel_melt_df_list[in_strata_order_ind]
  rel_melt_df <- do.call("rbind",rel_melt_df_list)
  rel_melt_df$stratum <- factor(rel_melt_df$stratum,levels=unique(rel_melt_df$stratum))
  number_of_legend_cols <- ifelse((dim(in_signatures_ind_df)[1] %% legend_height)==0,0,1) +
    (dim(in_signatures_ind_df)[1] %/% legend_height)
#   p_abs <- ggplot(abs_melt_df,environment = .e) + 
#     geom_bar(aes_string(x="PID_index",y="value",fill="sig_index"),stat='identity') +
#     facet_grid(stratum ~ ., scales="free_y") +
#     scale_fill_manual(labels=in_signatures_ind_df$sig,values=in_signatures_ind_df$colour) + 
#     scale_x_discrete("PID",breaks=seq(1,dim(in_subgroups_df)[1]),labels=in_subgroups_df$PID[my_sum_ind]) +
#     theme(axis.text.x=element_text(angle=90,size=7,vjust=0.5),
#           panel.background = element_rect(fill=NA, colour="black")) +
#     guides(fill=guide_legend(ncol=number_of_legend_cols))
#   p_rel <- ggplot(rel_melt_df,environment = .e) + 
#     geom_bar(aes_string(x="PID_index",y="value",fill="sig_index"),stat='identity') +
#     facet_grid(stratum ~ ., scales="free_y") +
#     scale_fill_manual(labels=in_signatures_ind_df$sig,values=in_signatures_ind_df$colour) + 
#     scale_x_discrete("PID",breaks=seq(1,dim(in_subgroups_df)[1]),labels=in_subgroups_df$PID[my_sum_ind]) +
#     theme(axis.text.x=element_text(angle=90,size=7,vjust=0.5),
#           panel.background = element_rect(fill=NA, colour="black")) +
#     guides(fill=guide_legend(ncol=number_of_legend_cols))
  abs_melt_df$method <- "abs"
  rel_melt_df$method <- "rel"
  all_melt_df <- rbind(abs_melt_df,rel_melt_df)
  all_melt_df$method <- factor(all_melt_df$method)
  p_all <- ggplot(all_melt_df,environment = .e) + 
    geom_bar(aes_string(x="PID_index",y="value",fill="sig_index"),stat='identity') +
    facet_wrap(~stratum+method, scales="free_y",ncol=2) +
    scale_fill_manual(labels=in_signatures_ind_df$sig,values=in_signatures_ind_df$colour) + 
    scale_x_discrete("PID",breaks=seq(1,dim(in_subgroups_df)[1]),labels=in_subgroups_df$PID[my_sum_ind]) +
    theme(axis.title.y=element_blank(),
          panel.background=element_rect(fill=NA, colour="black"),
          panel.border=element_rect(fill=NA, colour="black"),
          strip.background=element_rect(colour="black")) +
    guides(fill=guide_legend(ncol=number_of_legend_cols))
  if(in_label_orientation=="turn"){
    p_all <- p_all + theme(axis.text.x=element_text(angle=90,size=7,vjust=0.5))
  }
  return(p_all)
}


plot_group_facet <- function(in_df_list,in_name_list,
                             in_signatures_ind_df,
                             in_strata_order_ind=seq_len(length(in_name_list))){
  for(i in seq_len(length(in_df_list))){
    in_df_list[[i]]$stratum <- paste0(in_name_list[[i]],", cohort")
  }
  in_df_list <- in_df_list[in_strata_order_ind]
  all_melt_df <- do.call("rbind",in_df_list)
  all_melt_df$stratum <- factor(all_melt_df$stratum,levels=unique(all_melt_df$stratum))
  p_cohort <- ggplot(all_melt_df) + 
    ggplot2::geom_bar(aes_string(x="sig",y="exposure",fill="sig",
                                              size=0.3),
                      stat='identity',position="dodge",width=.7) + 
    scale_fill_manual(name="sig",labels=in_signatures_ind_df$sig,values=in_signatures_ind_df$colour) +
    geom_errorbar(aes_string(x="sig",ymin="exposure_min",ymax="exposure_max"),width=0.5) +
    facet_wrap(~stratum,ncol=1) +
    theme(axis.title.y=element_blank(),
          panel.background = element_rect(fill=NA, colour="black"),
          panel.border=element_rect(fill=NA, colour="black"),
          strip.background=element_rect(colour="black")) +
    guides(size=FALSE)
  return(p_cohort)
}


#' Plot results of the Stratification of a Mutational Catalogue
#'
#' Plot a big composite figure with 3 columns: in the left column the per-PID
#' absolute exposures will be shown, in the middle column the per_PID
#' relative or normalized exposures will be shown, in the right column the
#' cohort-wide exposures are shown (averaged over PIDs).
#'
#' @param number_of_strata
#'  Number of strata as deduced from \code{link{SMC}}
#' @param output_path
#'  Path to file where the results are going to be stored. If NULL, the results will
#'  be plotted to the running environment.
#' @param decomposition_method
#'  String for the filename of the generated barplot.
#' @param number_of_sigs
#'  Number of signatures
#' @param name_list
#'  Names of the contructed strata.
#' @param exposures_strata_list
#'  The list of \code{s} strata specific exposures Hi, all are numerical data 
#'  frames with \code{l} rows and \code{m} columns, \code{l} being the number
#'  of signatures and \code{m} being the number of samples
#' @param this_signatures_ind_df
#'  A data frame containing meta information about the signatures
#' @param this_subgroups_df
#'  A data frame indicating which PID (patient or sample identifyier) belongs
#'  to which subgroup
#' @param in_strata_order_ind
#'  Index vector defining reordering of the strata
#' @param exposures_both_rel_df_list
#'  A list of \code{s} strata specific cohortwide (i.e. averaged over cohort)
#'  normalized exposures
#' @param cohort_method_flag
#'  Either or several of \code{c("all_PIDs","cohort","norm_PIDs")}, representing
#'  alternative ways to average over the cohort.
#' @param fig_width
#'  Width of the figure to be plotted
#' @param fig_height
#'  Height of the figure to be plotted
#' @param fig_type
#'  png or pdf
#' @param in_label_orientation
#'  Whether or not to turn the labels on the x-axis.
#' @param this_sum_ind
#'  Optional set of indices for reordering the PIDs
#'  
#' @examples
#' NULL
#'  
#' @return The function doesn't return any value.
#' 
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @export
#' 
plot_SMC <- function(number_of_strata,output_path,decomposition_method,number_of_sigs,name_list,
                     exposures_strata_list,this_signatures_ind_df,this_subgroups_df,
                     in_strata_order_ind,exposures_both_rel_df_list,cohort_method_flag,
                     fig_width=1200,fig_height=900,fig_type="png",
                     in_label_orientation="turn",this_sum_ind=NULL) {
  my_strata_order_ind <- c(1,in_strata_order_ind+1)
  in_abs_df_list <- add_as_fist_to_list(exposures_strata_list$sub_exposures_list,
                                        exposures_strata_list$exposures_all_df)
  in_rel_df_list <- add_as_fist_to_list(exposures_strata_list$sub_norm_exposures_list,
                                        exposures_strata_list$norm_exposures_all_df)
  in_name_list <- c("all",name_list)
  p_all <- plot_SMC_PID_facet(in_abs_df_list,in_rel_df_list,
                              this_signatures_ind_df,this_subgroups_df,
                              in_name_list,in_sum_ind=this_sum_ind,
                              in_subgroups.field="subgroup",
                              in_strata_order_ind=my_strata_order_ind,
                              in_label_orientation=in_label_orientation)
  p_all_no_legend <- p_all + guides(fill=FALSE)  
  p_all_short <- p_all +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank())
  p_all_short_no_legend <- p_all_short + guides(fill=FALSE)
  p_cohort <- plot_group_facet(exposures_both_rel_df_list,
                               in_name_list,this_signatures_ind_df,
                               in_strata_order_ind=my_strata_order_ind)
  p_cohort_comb <- p_cohort + theme(strip.text=element_blank())
  p_cohort_short <- p_cohort_comb +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank())    
  if(!is.null(output_path)){
    if(fig_type=="png") png(output_path,width=fig_width,height=fig_height)
    else if(fig_type=="pdf") pdf(output_path,width=fig_width,height=fig_height)
  }
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow =1, ncol = 14)))
  if(in_label_orientation=="turn"){
    print(p_all_short_no_legend, vp = vplayout(1, c(1:10)))
    print(p_cohort_short, vp = vplayout(1, c(11:14)))    
  } else {
    print(p_all_no_legend, vp = vplayout(1, c(1:10)))
    print(p_cohort_comb, vp = vplayout(1, c(11:14)))
  }
  if(!is.null(output_path)){
    dev.off()
  }
  return(list(p_all=p_all,p_cohort=p_cohort))
}




#' Plot averaged signature exposures per subgroup
#' 
#' Plot one averaged signature exposure pattern per subgroup. Uses
#' \code{\link{split_exposures_by_subgroups}}.
#' 
#' @param in_exposures_df
#'  Numerical data frame of the exposures (i.e. contributions of the
#'  different signatures to the number of point mutations per PID)
#' @param in_signatures_ind_df
#'  Data frame carrying additional information on the signatures
#' @param in_subgroups_df
#'  Data frame indicating which PID belongs to which subgroup
#' @param in_subgroups.field
#'  Name indicating which column in \code{in_subgroups_df} contains the
#'  subgroup information
#' @param in_PID.field
#'  Name indicating which column in \code{in_subgroups_df} contains the
#'  PID information 
#'  
#' @return NULL
#' 
#' @seealso \code{\link{split_exposures_by_subgroups}}
#' 
#' @examples
#'  NULL
#'  
#' @export
#' 
stat_plot_subgroups_old <- function(in_exposures_df,in_subgroups_df,
                                    in_signatures_ind_df,
                                    in_subgroups.field="subgroup",
                                    in_PID.field="PID"){
  plots_per_row <- 4
  # split the data with custom function
  exposures_df_list <- split_exposures_by_subgroups(
    in_exposures_df=in_exposures_df,in_subgroups_df=in_subgroups_df,
    in_subgroups.field=in_subgroups.field,
    in_PID.field=in_PID.field)
  # compute statistical quantities to be displayed
  plot_df_list <- lapply(exposures_df_list, FUN=function(x) {
    average_vec <- average_over_present(x,1)
    SEM_vec <- stderrmean_over_present(x,1)
    max_vec <- average_vec + SEM_vec
    min_vec <- average_vec - SEM_vec
    out_df <- data.frame(index=factor(seq_len(length(average_vec))),
                         mean=average_vec,SEM=SEM_vec,
                         max=max_vec,min=min_vec)
    return(out_df)
  })
  # make similar data frame without splitting into subgroups
  average_vec <- average_over_present(in_exposures_df,1)
  SEM_vec <- stderrmean_over_present(in_exposures_df,1)
  max_vec <- average_vec + SEM_vec
  min_vec <- average_vec - SEM_vec
  plot_df_all <- data.frame(index=factor(seq_len(length(average_vec))),
                            mean=average_vec,SEM=SEM_vec,
                            max=max_vec,min=min_vec)
  all_ind <- length(plot_df_list)+1
  plot_df_list[[all_ind]] <- plot_df_all
  names(plot_df_list)[all_ind] <- "all"
  # make the plots
  plot_list_subgroups <- list()
  for(i in seq_len(length(plot_df_list))){
    temp_df <- plot_df_list[[i]]
    plot_list_subgroups[[i]] <- ggplot() + 
      ggplot2::geom_bar(data=temp_df,aes_string(x="index",y="mean",fill="index",
                                                size=0.3),
                        stat='identity',position="dodge",width=.7) + 
      scale_fill_manual(name="sig",labels=in_signatures_ind_df$sig,values=in_signatures_ind_df$colour) +
      geom_errorbar(data=temp_df,aes_string(x="index",ymin="min",ymax="max"),width=0.5) +
      labs(x="",y="",title=names(plot_df_list)[i]) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
            legend.position = "none")    
  }
  # order the plots to fit nicely in big panel
  number_of_subgroups <- length(plot_df_list)
  numer_of_rows <- ceiling(number_of_subgroups/plots_per_row)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow = numer_of_rows, 
                                             ncol = plots_per_row)))
  for(i in seq_len(length(plot_df_list))){
    my_row <- ceiling(i/plots_per_row)
    my_col <- i-((my_row-1)*plots_per_row)
    print(plot_list_subgroups[[i]], vp = vplayout(my_row, my_col))    
  }
}


#' Plot averaged signature exposures per subgroup
#' 
#' Plot one averaged signature exposure pattern per subgroup. Uses
#' \code{\link{split_exposures_by_subgroups}}.
#' 
#' @param in_exposures_df
#'  Numerical data frame of the exposures (i.e. contributions of the
#'  different signatures to the number of point mutations per PID)
#' @param in_signatures_ind_df
#'  Data frame carrying additional information on the signatures
#' @param in_subgroups_df
#'  Data frame indicating which PID belongs to which subgroup
#' @param in_subgroups.field
#'  Name indicating which column in \code{in_subgroups_df} contains the
#'  subgroup information
#' @param in_PID.field
#'  Name indicating which column in \code{in_subgroups_df} contains the
#'  PID information 
#'  
#' @return NULL
#' 
#' @seealso \code{\link{split_exposures_by_subgroups}}
#' 
#' @examples
#'  NULL
#'  
#' @export
#' 
stat_plot_subgroups <- function(in_exposures_df,in_subgroups_df,
                                in_signatures_ind_df,
                                in_subgroups.field="subgroup",
                                in_PID.field="PID",
                                in_colour_vector=NULL){
  .e <- environment()
  plots_per_row <- 4
  # split the data with custom function
  exposures_df_list <- split_exposures_by_subgroups(
    in_exposures_df=in_exposures_df,in_subgroups_df=in_subgroups_df,
    in_subgroups.field=in_subgroups.field,
    in_PID.field=in_PID.field)
  number_of_subgroups <- length(exposures_df_list)
  # compute statistical quantities to be displayed
  my_sig_names <- rownames(exposures_df_list[[1]])
  plot_df_list <- lapply(exposures_df_list, FUN=function(x) {
    average_vec <- average_over_present(x,1)
    SEM_vec <- stderrmean_over_present(x,1)
    max_vec <- average_vec + SEM_vec
    min_vec <- average_vec - SEM_vec
    out_df <- data.frame(index=factor(seq_len(length(average_vec))),
                         mean=average_vec,SEM=SEM_vec,
                         max=max_vec,min=min_vec,
                         sig=my_sig_names)
    return(out_df)
  })
  # make similar data frame without splitting into subgroups
  average_vec <- average_over_present(in_exposures_df,1)
  SEM_vec <- stderrmean_over_present(in_exposures_df,1)
  max_vec <- average_vec + SEM_vec
  min_vec <- average_vec - SEM_vec
  plot_df_all <- data.frame(index=factor(seq_len(length(average_vec))),
                            mean=average_vec,SEM=SEM_vec,
                            max=max_vec,min=min_vec,
                            sig=my_sig_names)
  plot_df_list <- add_as_fist_to_list(plot_df_list,plot_df_all)
  names(plot_df_list)[1] <- "all"
  for(my_stratum in names(plot_df_list)){
    plot_df_list[[my_stratum]]$subgroup <- my_stratum  
  }
  plot_df <- do.call(rbind,plot_df_list)
  plot_df$subgroup <- factor(plot_df$subgroup,
                         levels=unique(plot_df$subgroup))
  plot_df$sig <- factor(plot_df$sig,levels=unique(plot_df$sig))
  # make the plots
  p <- ggplot(plot_df,environment=.e) +
    geom_bar(aes_string(x="sig",y="mean",fill="sig",size=0.3),
             stat='identity',position="dodge",width=.7) +
    scale_fill_manual(name="sig",labels=in_signatures_ind_df$sig,values=in_signatures_ind_df$colour) +
    geom_errorbar(aes_string(x="sig",ymin="min",ymax="max"),width=0.5) +
    facet_wrap(~subgroup)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background = element_rect(fill=NA, colour="black"),
          panel.border=element_rect(fill=NA, colour="black")) +
    guides(size=FALSE)
  p1 <- p + theme(strip.background=element_rect(colour="black"))
  # display as dodged barplots
  my_palette <- c("black",rainbow(number_of_subgroups))
  names(my_palette) <- unique(plot_df$subgroup)
  if(!is.null(in_colour_vector)){
    match_ind <- match(names(in_colour_vector),names(my_palette))
    my_palette[match_ind] <- in_colour_vector
  }
  q <- ggplot(plot_df,environment=.e,aes(sig,y=mean,group=subgroup)) +
    geom_bar(aes(fill=sig,col=subgroup,size=0.3),
             stat='identity',position="dodge",width=.7,size=1) +
    scale_fill_manual(name="sig",labels=in_signatures_ind_df$sig,values=in_signatures_ind_df$colour) +
    scale_colour_manual(values=my_palette) +
    geom_errorbar(aes(ymin=min,ymax=max),width=0.5,position=position_dodge(width=0.7)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background = element_rect(fill=NA, colour="black"),
          panel.border=element_rect(fill=NA, colour="black")) +
    guides(size=FALSE)
  return(list(facet_plot=p1,dodged_plot=q))
}



#' Cluster the PIDs according to their signature exposures
#'
#' The PIDs are clustered according to their signature exposures by calling
#' first creating a distance matrix:
#' \itemize{
#'  \item \code{\link{dist}}, then
#'  \item \code{\link{hclust}} and then
#'  \item \code{\link[dendextend]{labels_colors}} to colour the labels (the
#'    text) of the leaves in the dendrogram.
#' }
#' Typically one colour per subgroup.
#'
#' @param in_exposures_df
#'  Numerical data frame encoding the exposures \code{H}, i.e. which
#'  signature contributes how much to which PID (patient identifier or sample).
#' @param in_subgroups_df
#'  A data frame indicating which PID (patient or sample identifyier) belongs
#'  to which subgroup
#' @param in_method
#'  Method of the clustering to be supplied to \code{\link{dist}}. Can be either
#'  of: \code{euclidean}, \code{maximum}, \code{manhattan}, \code{canberra},
#'  \code{binary} or \code{minkowski}
#' @param in_subgroup_column
#'  Indicates the name of the column in which the subgroup information
#'  is encoded in \code{in_subgroups_df}
#' @param in_palette
#'  Palette with colours or colour codes for the labels (the text) of the leaves
#'  in the dendrogram. Typically one colour per subgroup. If none is specified, a
#'  rainbow palette of the length of the number of subgroups will be used as
#'  default.
#' @param in_cutoff
#'  A numeric value less than 1. Signatures from within \code{W}
#'  with an overall exposure less than \code{in_cutoff} will be
#'  discarded for the clustering.
#' @param in_filename
#'  A path to save the dendrogram. If none is specified, the figure will be plotted
#'  to the running environment.
#' @param in_shift_factor
#'  Graphical parameter to adjust figure to be created
#' @param in_cex
#'  Graphical parameter to adjust figure to be created
#' @param in_title
#'  Title in the figure to be created under \code{in_filename}
#'
#' @return A list with entries
#'  \code{hclust} and
#'  \code{dendrogram}.
#' \itemize{
#'  \item \code{hclust}:
#'    The object created by \code{\link{hclust}}
#'  \item \code{dendrogram}:
#'    The above object wrapped in \code{\link{as.dendrogram}}
#' }
#'
#' @examples
#'  data(lymphoma_cohort_LCD_results)
#'  hclust_exposures(rel_lymphoma_Nature2013_COSMIC_cutoff_exposures_df,
#'                   COSMIC_subgroups_df,
#'                   in_method="manhattan",
#'                   in_subgroup_column="subgroup")
#'
#' @seealso \code{\link{hclust}}
#' @seealso \code{\link{dist}}
#' @seealso \code{\link[dendextend]{labels_colors}}
#' 
#' @import dendextend
#' @export
#' 
hclust_exposures <- function(in_exposures_df,in_subgroups_df,
                             in_method="manhattan",
                             in_subgroup_column="subgroup",
                             in_palette=NULL,in_cutoff=0,in_filename=NULL,
                             in_shift_factor=0.3,in_cex=0.2,in_title="") {
  ## 1. choose only signatures with exposures above cutoff
  exposures_sum_df <- data.frame(sum=apply(in_exposures_df,1,sum))
  exposures_sum_df$sum_norm <- exposures_sum_df$sum/sum(exposures_sum_df$sum)
  sig_choice_ind <- which(exposures_sum_df$sum_norm >= in_cutoff)
  reduced_exposures_df <- in_exposures_df[sig_choice_ind,]
  ## 2. adapt for hierarchic clustering
  my_exposures_df <- as.data.frame(t(reduced_exposures_df))
  ## 3. run clustering
  my_exposures_hierarchy <- hclust(dist(my_exposures_df,method=in_method))
  my_exposures_hc <- as.dendrogram(my_exposures_hierarchy)
  ## 4. colour the labels of the leaves according to subgroups
  subgroup_column_index <- which(tolower(names(in_subgroups_df))==tolower(in_subgroup_column))
  number_of_subgroups <- length(unique(in_subgroups_df[,subgroup_column_index]))
  if(is.null(in_palette)){
    in_palette=rainbow(number_of_subgroups)
  }
  colorCodes <- in_palette[seq_len(number_of_subgroups)]
  labels_colors(my_exposures_hc) <- colorCodes[factor(in_subgroups_df[,subgroup_column_index])][order.dendrogram(my_exposures_hc)]
  ## 5. prepare for plotting and plot
  max_height <- max(my_exposures_hierarchy$height)
  if (!is.null(in_filename)) {
    png(in_filename,width=1000,height=400)
  }
  plot(my_exposures_hc,cex=in_cex,ylim=c((-1)*in_shift_factor*max_height,max_height),main=in_title)
  if (!is.null(in_filename)) {
    dev.off()    
  }
  return(list(hclust=my_exposures_hierarchy,dendrogram=my_exposures_hc))
}


#' Create heatmap to cluster the PIDs according to their signature exposures
#'
#' The PIDs are clustered according to their signature exposures. The procedure
#' is analogous to \code{\link{hclust_exposures}}, but the graphical output is
#' more detailed due to the heatmap. This function calls: 
#' \itemize{
#'  \item \code{\link[gplots]{heatmap.2}} for all statistics.
#' }
#'
#' @param in_exposures_df
#'  Numerical data frame encoding the exposures \code{H}, i.e. which
#'  signature contributes how much to which PID (patient identifier or sample).
#' @param in_subgroups_df
#'  A data frame indicating which PID (patient or sample identifyier) belongs
#'  to which subgroup
#' @param in_signatures_ind_df
#'  A data frame containing meta information about the signatures, especially the
#'  asserted colour
#' @param in_method
#'  Method of the clustering to be supplied to \code{\link{dist}}. Can be either
#'  of: \code{euclidean}, \code{maximum}, \code{manhattan}, \code{canberra},
#'  \code{binary} or \code{minkowski}
#' @param in_subgroup_column
#'  Indicates the name of the column in which the subgroup information
#'  is encoded in \code{in_subgroups_df}
#' @param in_subgroup_colour_column
#'  Indicates the name of the column in which the colour information for
#'  subgroups is encoded in \code{in_subgroups_df}. If NULL, a rainbow palette
#'  is used instead.
#' @param in_palette
#'  Palette with colours or colour codes for the labels (the text) of the leaves
#'  in the dendrogram. Typically one colour per subgroup. If none is specified, a
#'  rainbow palette of the length of the number of subgroups will be used as
#'  default.
#' @param in_cutoff
#'  A numeric value less than 1. Signatures from within \code{W}
#'  with an overall exposure less than \code{in_cutoff} will be
#'  discarded for the clustering.
#' @param in_filename
#'  A path to save the heatmap. If none is specified, the figure will be plotted
#'  to the running environment.
#' @param in_title
#'  Title in the figure to be created under \code{in_filename}
#'
#' @return The function doesn't return any value.
#' 
#' @examples
#'  data(lymphoma_cohort_LCD_results)
#'  heatmap_exposures(rel_lymphoma_Nature2013_COSMIC_cutoff_exposures_df,
#'                    COSMIC_subgroups_df,
#'                    chosen_signatures_indices_df,
#'                    in_subgroup_colour_column="col",
#'                    in_method="manhattan",
#'                    in_subgroup_column="subgroup")
#'
#' @seealso \code{\link{hclust}}
#' @seealso \code{\link{dist}}
#' @seealso \code{\link[dendextend]{labels_colors}}
#' 
#' @importFrom gplots heatmap.2
#' @export
#' 
heatmap_exposures <- function(in_exposures_df,in_subgroups_df,in_signatures_ind_df,
                              in_method="manhattan",in_subgroup_column="subgroup",
                              in_subgroup_colour_column=NULL,
                              in_palette=NULL,in_cutoff=0,in_filename=NULL,in_title="") {
  ## 1. choose only signatures with exposures above cutoff
  exposures_sum_df <- data.frame(sum=apply(in_exposures_df,1,sum))
  exposures_sum_df$sum_norm <- exposures_sum_df$sum/sum(exposures_sum_df$sum)
  sig_choice_ind <- which(exposures_sum_df$sum_norm >= in_cutoff)
  reduced_exposures_df <- in_exposures_df[sig_choice_ind,]
  ## 2. reformat for heatmap.2 function
  reduced_exposures_matrix <- as.matrix(reduced_exposures_df)
  ## 3. custom colour code for the heatmap
  heatmap_palette <- c("lightyellow","yellow","orange","red","red","red","darkred","darkred","darkred")
  myColorRange <- colorRampPalette(heatmap_palette)(n=100)
  ## 4. colour the labels of the leaves according to subgroups
  subgroup_column_index <- which(tolower(names(in_subgroups_df))==tolower(in_subgroup_column))
  number_of_subgroups <- length(unique(in_subgroups_df[,subgroup_column_index]))
  if(is.null(in_subgroup_colour_column)){
    if(is.null(in_palette)){
      in_palette=rainbow(number_of_subgroups)
    }
    colorCodes <- in_palette[seq_len(number_of_subgroups)] 
    subgroups_colour_vector <- colorCodes[factor(in_subgroups_df[,subgroup_column_index])]
  } else {
    subgroup_colour_ind <- which(names(in_subgroups_df)==in_subgroup_colour_column)
    colorCodes <- unique(in_subgroups_df[,subgroup_colour_ind])
    subgroups_colour_vector<- in_subgroups_df[,subgroup_colour_ind]
  }
  colorNames <- sapply(colorCodes,function(l) {
    my_ind <- min(which(subgroups_colour_vector==l))
    return(in_subgroups_df[my_ind,subgroup_column_index])
                                              })
  sig_colour_vector <- in_signatures_ind_df$colour
  if (!is.null(in_filename)) {
    png(in_filename,width=800,height=800)
  }
  heatmap.2(reduced_exposures_matrix, distfun=function(i)
    dist(i,method=in_method), hclustfun=function(i)
      hclust(i,method="complete"),
    trace="none",main=in_title,col=myColorRange,margins=c(12,8),srtCol=45,
    ColSideColors=subgroups_colour_vector,
    RowSideColors=sig_colour_vector)
  legend("topright",
         legend = colorNames,
         col = colorCodes,
         lty= 1,
         lwd = 10)
  if (!is.null(in_filename)) {
    dev.off()    
  }
  return()
}


#' Heatmap to cluster the PIDs on their signature exposures (with ComplexHeatmap)
#'
#' The PIDs are clustered according to their signature exposures. The procedure
#' is analogous to \code{\link{heatmap_exposures}}, but unsing the package
#' \pkg{ComplexHeatmap} by Zuguang Gu instead. This function calls: 
#' \itemize{
#'  \item \code{\link[ComplexHeatmap]{rowAnnotation}},
#'  \item \code{\link[ComplexHeatmap]{HeatmapAnnotation}} and
#'  \item \code{\link[ComplexHeatmap]{Heatmap}}
#' }
#'
#' @param in_exposures_df
#'  Numerical data frame encoding the exposures \code{H}, i.e. which
#'  signature contributes how much to which PID (patient identifier or sample).
#' @param in_subgroups_df
#'  A data frame indicating which PID (patient or sample identifyier) belongs
#'  to which subgroup
#' @param in_signatures_ind_df
#'  A data frame containing meta information about the signatures, especially the
#'  asserted colour
#' @param in_data_type
#'  Title in the figure
#' @param in_method
#'  Method of the clustering to be supplied to \code{\link{dist}}. Can be either
#'  of: \code{euclidean}, \code{maximum}, \code{manhattan}, \code{canberra},
#'  \code{binary} or \code{minkowski}
#' @param in_subgroup_column
#'  Indicates the name of the column in which the subgroup information
#'  is encoded in \code{in_subgroups_df}
#' @param in_subgroup_colour_column
#'  Indicates the name of the column in which the colour information for
#'  subgroups is encoded in \code{in_subgroups_df}. If NULL, a rainbow palette
#'  is used instead.
#' @param in_palette
#'  Palette with colours for the heatmap. Default is
#'  \code{colorRamp2(c(0, 0.2, 0.4, 0.6), c('white', 'yellow', 'orange', 'red'))}
#' @param in_cutoff
#'  A numeric value less than 1. Signatures from within \code{W}
#'  with an overall exposure less than \code{in_cutoff} will be
#'  discarded for the clustering.
#' @param in_filename
#'  A path to save the heatmap. If none is specified, the figure will be plotted
#'  to the running environment.
#' @param in_column_anno_borders
#'  Whether or not to draw separating lines between the fields in the annotation
#' @param in_row_anno_borders
#'  Whether or not to draw separating lines between the fields in the annotation
#'
#' @details
#'  It might be necessary to install the newest version of the development branch
#'  of the packages \pkg{circlize} and \pkg{ComplexHeatmap} by Zuguang Gu:
#'  \code{devtools::install_github("jokergoo/circlize")}
#'  \code{devtools::install_github("jokergoo/ComplexHeatmap")}
#'
#' @return The function doesn't return any value.
#' 
#' @examples
#'  data(lymphoma_cohort_LCD_results)
#'  complex_heatmap_exposures(rel_lymphoma_Nature2013_COSMIC_cutoff_exposures_df,
#'                            COSMIC_subgroups_df,
#'                            chosen_signatures_indices_df,
#'                            in_data_type="norm exposures",
#'                            in_subgroup_colour_column="col",
#'                            in_method="manhattan",
#'                            in_subgroup_column="subgroup")
#'
#' @seealso \code{\link[ComplexHeatmap]{Heatmap}}
#' @seealso \code{\link{heatmap_exposures}}
#'
#' @import ComplexHeatmap
#' @import circlize
#' @export
#' 
complex_heatmap_exposures <- function(in_exposures_df,in_subgroups_df,in_signatures_ind_df,
                                      in_data_type="norm exposures",
                                      in_method="manhattan",in_subgroup_column="subgroup",
                                      in_subgroup_colour_column=NULL,
                                      in_palette=colorRamp2(c(0,0.2,0.4,0.6),c('white','yellow','orange','red')),
                                      in_cutoff=0,in_filename=NULL,
                                      in_column_anno_borders=FALSE,
                                      in_row_anno_borders=FALSE){
  if(!in_row_anno_borders){
    row_anno = rowAnnotation(df = data.frame(signature = in_signatures_ind_df$index),
                           col = list(signature = structure(in_signatures_ind_df$colour, 
                                                      names = as.character(in_signatures_ind_df$index))),
                           show_legend = FALSE)
  } else {  
    row_anno = rowAnnotation(df = data.frame(signature = in_signatures_ind_df$index),
                             col = list(signature = structure(in_signatures_ind_df$colour, 
                                                              names = as.character(in_signatures_ind_df$index))),
                             show_legend = FALSE,
                             gp = gpar(col="black"))
  }
  if(!in_column_anno_borders){
    column_anno = HeatmapAnnotation(df = data.frame(subtype = in_subgroups_df[[in_subgroup_column]]), 
                                    col = list(subtype = structure(unique(in_subgroups_df[[in_subgroup_colour_column]]),
                                                                   names = unique(in_subgroups_df[[in_subgroup_column]]))))
  } else {
    column_anno = HeatmapAnnotation(df = data.frame(subtype = in_subgroups_df[[in_subgroup_column]]), 
                                    col = list(subtype = structure(unique(in_subgroups_df[[in_subgroup_colour_column]]),
                                                                   names = unique(in_subgroups_df[[in_subgroup_column]]))),
                                    gp = gpar(col="black"))   
  }
  ht_list = row_anno + Heatmap(in_exposures_df, col = in_palette,
                               top_annotation = column_anno, 
                               clustering_distance_rows = in_method,
                               clustering_distance_columns = in_method,
                               heatmap_legend_param = list(title = in_data_type))
  draw(ht_list, row_dend_side = 'left')
}


#' Heatmap to cluster the PIDs on their signature exposures (with ComplexHeatmap)
#'
#' The PIDs are clustered according to their signature exposures. The procedure
#' is analogous to \code{\link{complex_heatmap_exposures}}, but enabling more
#' than one annotation row for the PIDs. This function calls: 
#' \itemize{
#'  \item \code{\link[ComplexHeatmap]{rowAnnotation}},
#'  \item \code{\link[ComplexHeatmap]{HeatmapAnnotation}} and
#'  \item \code{\link[ComplexHeatmap]{Heatmap}}
#' }
#'
#' @param in_exposures_df
#'  Numerical data frame encoding the exposures \code{H}, i.e. which
#'  signature contributes how much to which PID (patient identifier or sample).
#' @param in_annotation_df
#'  A data frame indicating which PID (patient or sample identifyier) belongs
#'  to which subgroup for all layers of annotation
#' @param in_annotation_col
#'  A list indicating colour attributions for all layers of annotation
#' @param in_signatures_ind_df
#'  A data frame containing meta information about the signatures, especially the
#'  asserted colour
#' @param in_data_type
#'  Title in the figure
#' @param in_method
#'  Method of the clustering to be supplied to \code{\link{dist}}. Can be either
#'  of: \code{euclidean}, \code{maximum}, \code{manhattan}, \code{canberra},
#'  \code{binary} or \code{minkowski}
#' @param in_palette
#'  Palette with colours or colour codes for the heatmap. Default is
#'  \code{colorRamp2(c(0, 0.2, 0.4, 0.6), c('white', 'yellow', 'orange', 'red'))}
#' @param in_cutoff
#'  A numeric value less than 1. Signatures from within \code{W}
#'  with an overall exposure less than \code{in_cutoff} will be
#'  discarded for the clustering.
#' @param in_filename
#'  A path to save the heatmap. If none is specified, the figure will be plotted
#'  to the running environment.
#' @param in_column_anno_borders
#'  Whether or not to draw separating lines between the fields in the annotation
#' @param in_row_anno_borders
#'  Whether or not to draw separating lines between the fields in the annotation
#' @param in_show_PIDs
#'  Whether or not to show the PIDs on the x-axis
#' @param in_annotation_legend_side
#'  Where to put the legends of the annotation df, default is right.
#'
#' @details
#'  It might be necessary to install the newest version of the development branch
#'  of the packages \pkg{circlize} and \pkg{ComplexHeatmap} by Zuguang Gu:
#'  \code{devtools::install_github("jokergoo/circlize")}
#'  \code{devtools::install_github("jokergoo/ComplexHeatmap")}
#'
#' @return The function doesn't return any value.
#' 
#' @examples
#'  NULL
#'
#' @seealso \code{\link[ComplexHeatmap]{Heatmap}}
#' @seealso \code{\link{heatmap_exposures}}
#' @seealso \code{\link{complex_heatmap_exposures}}
#'
#' @import ComplexHeatmap
#' @import circlize
#' @export
#' 
annotation_heatmap_exposures <- function(in_exposures_df,
                                         in_annotation_df,
                                         in_annotation_col,
                                         in_signatures_ind_df,
                                         in_data_type="norm exposures",
                                         in_method="manhattan",
                                         in_palette=colorRamp2(c(0, 0.2, 0.4, 0.6), 
                                                               c('white', 'yellow', 'orange', 'red')),
                                         in_cutoff=0,in_filename=NULL,
                                         in_column_anno_borders=FALSE,
                                         in_row_anno_borders=FALSE,
                                         in_show_PIDs=TRUE,
                                         in_annotation_legend_side="right"){
  if(!in_row_anno_borders){
    row_anno = rowAnnotation(df = data.frame(signature = in_signatures_ind_df$index),
                             col = list(signature = structure(in_signatures_ind_df$colour, 
                                                              names = as.character(in_signatures_ind_df$index))),
                             show_legend = FALSE)        
  } else {
    row_anno = rowAnnotation(df = data.frame(signature = in_signatures_ind_df$index),
                             col = list(signature = structure(in_signatures_ind_df$colour, 
                                                              names = as.character(in_signatures_ind_df$index))),
                             show_legend = FALSE,
                             gp = gpar(col="black"))    
  }
  if(!in_column_anno_borders){
    column_anno = HeatmapAnnotation(df = in_annotation_df, 
                                    col = in_annotation_col)
  } else {
    column_anno = HeatmapAnnotation(df = in_annotation_df, 
                                    col = in_annotation_col,
                                    gp = gpar(col="black"))    
  }
  if(in_show_PIDs){
    ht_list = row_anno + Heatmap(in_exposures_df, col = in_palette,
                                 top_annotation = column_anno, 
                                 clustering_distance_rows = in_method,
                                 clustering_distance_columns = in_method,
                                 heatmap_legend_param = list(title = in_data_type))
  } else {
    ht_list = row_anno + Heatmap(in_exposures_df, col = in_palette,
                                 top_annotation = column_anno, 
                                 clustering_distance_rows = in_method,
                                 clustering_distance_columns = in_method,
                                 heatmap_legend_param = list(title = in_data_type),
                                 show_column_names = FALSE)
  }
  draw(ht_list, row_dend_side = 'left',annotation_legend_side=in_annotation_legend_side)
}


#' Create a rainfall plot in a trellis structure
#'
#' A trellis is a plot structure which allows space optimized multi-panel multi
#' track plots. This function uses the package \pkg{gtrellis} developed by
#' Zuguang Gu, also available at 
#' \url{http://www.bioconductor.org/packages/release/bioc/html/gtrellis.html}.
#' The graphics in the tracks within a gtrellis plot are mostly drawn with
#' functions from the package \pkg{grid}. Note that for technical reasons,
#' the column indicating the chromosome MUST have the name \emph{chr}
#' and be the first column in the data frame supplied to the gtrellis
#' functions. Therefore reformatting is performed in this function before
#' calling gtrellis functions.
#'
#' @param in_rainfall_dat
#'  Data frame which has to contain at least columns for chromosome,
#'  position, intermutational distance and colour information
#' @param in_point_size
#'  size of the points in the rainfall plot to be created has to be
#'  provided with appropriate units, e.g. in_point_size=unit(0.5,"mm")
#' @param in_rect_list
#'  Optional argument, if present, will lead to highlighting of specified
#'  regions by coloured but transparent rectangles
#' @param in_title
#'  Title in the figure to be created.
#' @param in_CHROM.field
#'  String indicating which column of \code{in_rainfall_dat} carries the
#'  chromosome information
#' @param in_POS.field
#'  String indicating which column of \code{in_rainfall_dat} carries the
#'  position information
#' @param in_dist.field
#'  String indicating which column of \code{in_rainfall_dat} carries the
#'  intermutational distance information
#' @param in_col.field
#'  String indicating which column of \code{in_rainfall_dat} carries the
#'  colour information encoding the nucleotide exchange
#'
#' @return The function doesn't return any value.
#' 
#' @examples
#' NULL
#' \dontrun{
#'  data(lymphoma_test)
#'  choice_PID <- "4121361"
#'  PID_df <- subset(lymphoma_test_df,PID==choice_PID)
#'  trellis_rainfall_plot(PID_df,in_point_size=unit(0.5,"mm"))
#' }
#' @seealso \code{\link[gtrellis]{gtrellis_layout}}
#' @seealso \code{\link[gtrellis]{add_track}}
#' @seealso \code{\link[grid]{grid.points}}
#' 
#' @import gtrellis
#' @export
#' 
trellis_rainfall_plot <- function(in_rainfall_dat,in_point_size=unit(1,"mm"),
                                  in_rect_list=NULL,in_title="",
                                  in_CHROM.field="CHROM",in_POS.field="POS",
                                  in_dist.field="dist",in_col.field="col") {
  my_alpha <- 0.1
  CHROM_ind <- which(names(in_rainfall_dat)==in_CHROM.field)
  POS_ind <- which(names(in_rainfall_dat)==in_POS.field)
  dist_ind <- which(names(in_rainfall_dat)==in_dist.field)
  col_ind <- which(names(in_rainfall_dat)==in_col.field)
  names(in_rainfall_dat)[CHROM_ind] <- "chr"
  rainfall_dat <- in_rainfall_dat[,c(CHROM_ind,POS_ind,dist_ind,col_ind)]
  gtrellis_layout(n_track = 1, ncol = 5, byrow = FALSE,
                  title=in_title,
                  track_axis = TRUE,
                  track_height = unit.c(unit(1, "null")), 
                  track_ylim = c(0, 8),
                  track_ylab = c("intermut dist"),
                  add_name_track = TRUE, add_ideogram_track = TRUE)
  ## draw partially transparent rectangular regions if desired
  if(!is.null(in_rect_list)){
    for(i in seq_len(length(in_rect_list))){
      temp_list <- in_rect_list[[i]]
      temp_color <- temp_list$col
      temp_df <- temp_list$df
      temp_alpha <- temp_list$alpha
      add_track(temp_df, track = 2, panel.fun = function(gr) {
        grid.rect(gr[[2]], unit(0, "npc"), width=gr[[4]], height=unit(1, "npc"),
                  default.units = "native", hjust=0, vjust=0, gp = gpar(col=temp_color,fill =temp_color,alpha=temp_alpha))        
      })
    }
  }
  # track for rainfall plots
  add_track(rainfall_dat, track = 2, panel.fun = function(gr) {
    x = gr[[2]]
    y = log10(gr[[3]])
    grid.points(x, y, pch = 16, size = in_point_size, gp = gpar(col = gr[[4]]))
  }) 
}


#' Plot the spectra of nucleotide exchanges
#'
#' Plots the spectra of nucleotide exchanges in their triplet contexts. If
#' several columns are present in the input data frame, the spectra are plotted
#' for every column separately.
#'
#' @param in_catalogue_df
#'  Numerical data frame encoding the exchange spectra to be displayed, either
#'  a mutational catalogue \code{V} or a signatures matrix \code{W}.
#' @param in_colour_vector
#'  Specifies the colours of the 6 nucleotide exchanges if non-null.
#' @param in_show_triplets
#'  Whether or not to show the triplets on the x-axis
#' @param in_show_axis_title
#'  Whether or not to show the name of the y-axis
#'
#' @return The generated barplot - a ggplot2 plot
#'
#' @examples
#'  NULL
#'
#' @seealso \code{\link[ggplot2]{geom_bar}}
#' @seealso \code{\link[ggplot2]{facet_grid}}
#' 
#' @import ggplot2
#' @import reshape2
#' @export
#' 
plotExchangeSpectra <- function(in_catalogue_df,
                                in_colour_vector=NULL,
                                in_show_triplets=FALSE,
                                in_show_axis_title=FALSE){
  .e <- environment()
  in_catalogue_df$triplet_exchange <- rownames(in_catalogue_df)
  in_catalogue_df$nuc_exchange <- gsub(" .+$","",in_catalogue_df$triplet_exchange)
  in_catalogue_df$triplet <- gsub("^.+ ","",in_catalogue_df$triplet_exchange)
  catalogue_df_melt <- melt(in_catalogue_df,id.vars=c("nuc_exchange","triplet",
                                                      "triplet_exchange"))
  my_palette <- c("lightblue","black","red","grey","green","pink")
  names(my_palette) <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  if(!is.null(in_colour_vector)) my_palette <- in_colour_vector
  p <- ggplot(catalogue_df_melt,environment=.e) +
    geom_bar(aes(x=triplet,y=value,fill=nuc_exchange),
             stat='identity') +
    scale_fill_manual(name="exchange",values=my_palette) +
    facet_grid(variable~nuc_exchange,scales="free_x")
  p1 <- p +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  if(in_show_triplets) {
    p1 <- p +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_text(angle=90,vjust=0.5),
            axis.ticks.x=element_blank())
  }
  if(!in_show_axis_title) p1 <- p1 + theme(axis.title.y=element_blank())
  return(p1)
}