#' Group strata from different stratification axes
#'
#' For a comparison of the strata from different orthogonal stratification
#' axes, i.e. othogonal SMCs, the strata have to be grouped and reformatted.
#' This function does this task for the comparison by cosine similarity of
#' signature exposures. Output of this function is the basis for applying
#' \code{\link{plot_strata}} and \code{\link{make_comparison_matrix}}. It
#' is called by the wrapper functions \code{\link{compare_SMCs}},
#' \code{\link{run_plot_strata_general}} or 
#' \code{\link{run_comparison_general}}.
#'
#' @param in_stratification_lists_list
#'  List of lists with entries from different (orthogonal) stratification 
#'  axes or SMCs
#' @param in_remove_signature_ind
#'  Omit one of the signatures in \code{in_signatures_ind_df} for the 
#'  comparison if non-NULL. The parameter specifies the index of the
#'  signature to be removed.
#' @param in_additional_stratum
#'  Include an additionally supplied stratum in comparison in non-NULL.
#'  
#' @return A list with entries
#'  \code{strata_df},
#'  \code{number_of_SMCs},
#'  \code{number_of_strata}.
#' \itemize{
#'  \item \code{strata_df}:
#'   Pasted numerical data frame of all strata (these are going to be
#'   compared e.g. by \code{\link{make_comparison_matrix}}).
#'  \item \code{number_of_SMCs}:
#'    Number of orthogonal stratifications in
#'    \code{in_stratification_lists_list} and additional ones.
#'  \item \code{number_of_strata}:
#'    Cumulative number of strata (sum over the numbers of strata of
#'    the different stratifications in \code{in_stratification_lists_list})
#'    and additional ones.
#' }
#' 
#' @examples
#'  NULL
#'  
#' @seealso \code{\link{plot_strata}}
#' @seealso \code{\link{make_comparison_matrix}}
#' @seealso \code{\link{compare_SMCs}}
#' @seealso \code{\link{run_plot_strata_general}}
#' @seealso \code{\link{run_comparison_general}}
#' 
#' @export
#' 
make_strata_df <- function(in_stratification_lists_list,
                           in_remove_signature_ind=NULL,
                           in_additional_stratum=NULL) {
  my_number_of_signatures <- 
    dim(in_stratification_lists_list[[1]]$cohort[[1]])[1]
  my_SMCs <- names(in_stratification_lists_list)
  my_strata <- c()
  for(SMC in my_SMCs) {
    temp_strata_names <- 
      paste0(SMC,"_",unlist(in_stratification_lists_list[[SMC]]$name_list))
    temp_all <- paste0(SMC,"_all")
    my_strata <- c(my_strata,temp_all,temp_strata_names)
  }
  out_strata_df <- repeat_df(0,my_number_of_signatures,length(my_strata))
  stratum_counter <- 0
  SMC_counter <- 0
  my_number_of_strata <- c()
  for(SMC in my_SMCs) {
    SMC_counter <- SMC_counter + 1
    my_number_of_strata[SMC_counter] <- 
      length(in_stratification_lists_list[[SMC]]$cohort)
    for(stratum in seq_len(my_number_of_strata[SMC_counter])) {
      stratum_counter <- stratum_counter + 1
      temp_df <- in_stratification_lists_list[[SMC]]$cohort[[stratum]]
      out_strata_df[,stratum_counter] <- temp_df$exposure
      names(out_strata_df)[stratum_counter] <- my_strata[stratum_counter]
    }
  }
  rownames(out_strata_df) <- paste0("S",rownames(out_strata_df))
  my_number_of_SMCs <- length(my_SMCs)
  if(!(is.null(in_remove_signature_ind))) {
    out_strata_df <- out_strata_df[(-1)*in_remove_signature_ind,]
  }
  if(!(is.null(in_additional_stratum))) {
    out_strata_df$compare <- in_additional_stratum
    my_number_of_SMCs <- my_number_of_SMCs+1
    my_number_of_strata <- c(my_number_of_strata,1)
  }
  return(list(strata_df=out_strata_df,
              number_of_SMCs=my_number_of_SMCs,
              number_of_strata=my_number_of_strata))
}


#' Group strata from different stratification axes
#'
#' For a comparison of the strata from different orthogonal stratification
#' axes, i.e. othogonal SMCs, the strata have to be grouped and reformatted.
#' This function does this task for the comparison by cosine similarity of
#' mutational catalogues. Output of this function is the basis for applying
#' \code{\link{make_comparison_matrix}}. It is called by the wrapper function
#' \code{\link{run_comparison_catalogues}}.
#'
#' @param in_stratification_lists_list
#'  List of lists with entries from different (orthogonal) stratification 
#'  axes or SMCs
#' @param in_additional_stratum
#'  Include an additionally supplied stratum in comparison in non-NULL.
#'  
#' @return A list with entries
#'  \code{strata_df},
#'  \code{number_of_SMCs},
#'  \code{number_of_strata}.
#' \itemize{
#'  \item \code{strata_df}:
#'   Pasted numerical data frame of all strata (these are going to be
#'   compared e.g. by \code{\link{make_comparison_matrix}}).
#'  \item \code{number_of_SMCs}:
#'    Number of orthogonal stratifications in
#'    \code{in_stratification_lists_list} and additional ones.
#'  \item \code{number_of_strata}:
#'    Cumulative number of strata (sum over the numbers of strata of
#'    the different stratifications in \code{in_stratification_lists_list})
#'    and additional ones.
#' }
#' 
#' @examples
#'  NULL
#'  
#' @seealso \code{\link{plot_strata}}
#' @seealso \code{\link{make_comparison_matrix}}
#' @seealso \code{\link{run_comparison_catalogues}}
#' 
#' @export
#' 
make_catalogue_strata_df <- function(in_stratification_lists_list,
                                     in_additional_stratum=NULL) {
  my_number_of_PIDs <- 
    dim(in_stratification_lists_list[[1]]$catalogues_list[[1]])[2]
  my_number_of_features <- 
    dim(in_stratification_lists_list[[1]]$catalogues_list[[1]])[1]
  my_SMCs <- names(in_stratification_lists_list)
  my_strata <- c()
  ## normalize the catalogues
  normalized_catalogues_lists_list <- list()
  for(SMC in my_SMCs) {
    temp_strata_names <- 
      paste0(SMC,"_",unlist(in_stratification_lists_list[[SMC]]$name_list))
    temp_all <- paste0(SMC,"_all")
    my_strata <- c(my_strata,temp_all,temp_strata_names)
    my_catalogues_list <- in_stratification_lists_list[[SMC]]$catalogues_list
    sum_catalogue <- sum_over_list_of_df(my_catalogues_list)
    temp_length <- length(my_catalogues_list)
    my_catalogues_list[[temp_length+1]] <- sum_catalogue
    permut_ind <- c(temp_length+1,seq_len(temp_length))
    this_catalogues_list <- my_catalogues_list[permut_ind]
    normalized_catalogues_lists_list[[SMC]] <- 
      base::lapply(this_catalogues_list,
                   function(l) normalize_df_per_dim(l,2))
  }
  out_strata_df <- repeat_df(0,my_number_of_features,length(my_strata))
  stratum_counter <- 0
  SMC_counter <- 0
  my_number_of_strata <- c()
  for(SMC in my_SMCs) {
    SMC_counter <- SMC_counter + 1
    my_number_of_strata[SMC_counter] <- 
      length(in_stratification_lists_list[[SMC]]$cohort)
    for(stratum in seq_len(my_number_of_strata[SMC_counter])) {
      stratum_counter <- stratum_counter + 1
      temp_df <- normalized_catalogues_lists_list[[SMC]][[stratum]]
      out_strata_df[,stratum_counter] <- average_over_present(temp_df,1)
      names(out_strata_df)[stratum_counter] <- my_strata[stratum_counter]
    }
  }
  rownames(out_strata_df) <- 
    rownames(in_stratification_lists_list[[1]]$catalogues_list[[1]])
  my_number_of_SMCs <- length(my_SMCs)
  if(!(is.null(in_additional_stratum))) {
    out_strata_df$compare <- in_additional_stratum
    my_number_of_SMCs <- my_number_of_SMCs+1
    my_number_of_strata <- c(my_number_of_strata,1)
  }
  return(list(strata_df=out_strata_df,
              number_of_SMCs=my_number_of_SMCs,
              number_of_strata=my_number_of_strata))
}


#' Plot all strata from different stratification axes together
#'
#' Plot the cohort wide signature exposures of all strata from 
#' different stratification axes together. Naturally called by
#' \code{\link{compare_SMCs}}.
#'
#' @param in_strata_list
#'  Data structure created by \code{make_strata_df} or
#'  \code{make_catalogue_strata_df} in which the strata from different
#'  orthogonal stratification axes are reorganized in a consistent
#'  structure.
#' @param in_signatures_ind_df
#'  A data frame containing meta information about the signatures
#' @param output_path
#'  Path to directory where the results, especially the figure produced,
#'  are going to be stored.
#' @param in_attribute
#'  Additional string for the file name where the figure output
#'  is going to be stored.
#'
#' @return The function doesn't return any value.
#' 
#' @examples
#'  NULL
#'  
#' @seealso \code{\link{compare_SMCs}}
#'
#' @import ggplot2
#' @export
#' 
plot_strata <- function(in_strata_list,
                        in_signatures_ind_df,
                        output_path=NULL,
                        in_attribute="") {
  in_strata_df <- in_strata_list$strata_df
  number_of_SMCs <- in_strata_list$number_of_SMCs
  number_of_strata <- in_strata_list$number_of_strata
  plot_list <- list()
  stratum_counter <- 0
  for(stratum in names(in_strata_df)) {
    stratum_counter <- stratum_counter + 1
    temp_df <- data.frame(sig=rownames(in_strata_df),
                          exposure=in_strata_df[,stratum])
    plot_list[[stratum_counter]] <- ggplot() + 
      ggplot2::geom_bar(data=temp_df,
                        aes_string(x="sig",y="exposure",fill="sig",size=0.3),
                        stat='identity',position="dodge",width=.7) + 
      scale_fill_manual(values=in_signatures_ind_df$colour) + ## colour_vector needs to be supplied!!
      labs(x="",y="",title=stratum) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            legend.position = "none")
  }
  horizontal_element_width <- 10
  number_of_horizontal_units <- number_of_SMCs*horizontal_element_width
  vertical_element_height <- 10
  number_of_vertical_units <- (max(number_of_strata))*vertical_element_height
  horizontal_figure_factor <- 60
  vertical_figure_factor <- 40
  
  if(!is.null(output_path)){
    fileName <- file.path(output_path,paste0(in_attribute,"_all_strata.png"))
    png(fileName,
      width=number_of_horizontal_units*horizontal_figure_factor,
      height=number_of_vertical_units*vertical_figure_factor)
  }
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(number_of_vertical_units, 
                                             number_of_horizontal_units)))
  this_offsets <- c(0,cumsum(number_of_strata)[1:length(number_of_strata)-1])
  vertical_temp_stop <- 0
  for (i in seq_len(sum(number_of_strata))) {
    this_SMC <- findInterval(i,cumsum(number_of_strata)+1)+1
    this_stratum <- i - this_offsets[this_SMC]
    horizontal_temp_start <- (this_SMC-1)*horizontal_element_width+1
    horizontal_temp_stop <- 
      (horizontal_temp_start - 1) + horizontal_element_width
    vertical_temp_start <- vertical_temp_stop + 1
    vertical_temp_start <- (this_stratum-1)*vertical_element_height+1
    vertical_temp_stop <- (vertical_temp_start - 1) + vertical_element_height    
    print(plot_list[[i]], 
          vp = vplayout(vertical_temp_start:vertical_temp_stop, 
                        horizontal_temp_start:horizontal_temp_stop))
  }
  if(!is.null(output_path)){
    dev.off()
  }
  return()
}


#' Compute a similarity matrix for different strata
#'
#' Compute and plot a similarity matrix for different strata from 
#' different stratification axes together. First, \code{\link{compare_sets}} is
#' called on \code{in_strata_df} with itself, yielding a distance matrix (a
#' numerical data frame) \code{dist_df} of the strata. The corresponding
#' similarity matrix \code{1-dif_df} is then passed to
#' \code{\link[corrplot]{corrplot}}.
#' 
#'
#' @param in_strata_df
#'  Numerical data frame of all strata to be compared.
#' @param output_path
#'  Path to directory where the results, especially the figure produced by
#'  \code{\link[corrplot]{corrplot}} is going to be stored.
#' @param in_nrect
#'  Number of clusters in the clustering procedure provided by
#'  \code{\link[corrplot]{corrplot}}
#' @param in_attribute
#'  Additional string for the file name where the figure produced by
#'  \code{\link[corrplot]{corrplot}} is going to be stored.
#' @param in_palette
#'  Colour palette for the matrix
#'  
#' @return The comparison matrix of cosine similarities.
#' 
#' @examples
#' data(sigs)
#' make_comparison_matrix(
#'  AlexCosmicValid_sig_df,in_nrect=9,
#'  in_palette=colorRampPalette(c("blue","green","red"))(n=100))
#' 
#' @seealso \code{\link{compare_SMCs}}
#' 
#' @importFrom grDevices colorRampPalette png dev.off
#' @importFrom corrplot corrplot
#' @export
#' 
make_comparison_matrix <- function(in_strata_df,output_path=NULL,in_nrect=5,
                                   in_attribute="",in_palette=NULL) {
  if(is.null(in_palette)){
    myColorRange <- colorRampPalette(c("blue","green","green","green","green",
                                       "green","green","green","green","green",
                                       "red"))(n=100)    
  } else {
    myColorRange <- in_palette
  }
  compare_list <- compare_sets(in_strata_df,in_strata_df)
  dist_df <- compare_list$distance
  if(!is.null(output_path)){
    fileName <- file.path(output_path,
                          paste0(in_attribute,"_comparison_matrix.png"))
    png(fileName,width=500,height=500)
  }
  corrplot(as.matrix(1-dist_df),method="color",order="hclust",
           addrect=in_nrect,col=myColorRange)
  if(!is.null(output_path)){
    dev.off()
  }
  return(as.matrix(1-dist_df))
}


#' Compare all strata from different stratifications
#'
#' Compare all strata from different orthogonal stratification axes, i.e.
#' othogonal SMCs by cosine similarity of signature exposures. First calls
#' \itemize{
#'  \item \code{make_strata_df}, then
#'  \item \code{\link{plot_strata}} and finally
#'  \item \code{\link{make_comparison_matrix}}
#' }
#'
#' @param in_stratification_lists_list
#'  List of lists with entries from different (orthogonal) stratification 
#'  axes or SMCs
#' @param in_signatures_ind_df
#'  A data frame containing meta information about the signatures
#' @param output_path
#'  Path to directory where the results, especially the figure produced by
#'  \code{\link[corrplot]{corrplot}} is going to be stored.
#' @param in_nrect
#'  Number of clusters in the clustering procedure provided by
#'  \code{\link[corrplot]{corrplot}}
#' @param in_attribute
#'  Additional string for the file name where the figure produced by
#'  \code{\link[corrplot]{corrplot}} is going to be stored.
#'  
#' @return The comparison matrix of cosine similarities.
#' 
#' @examples
#'  NULL
#'  
#' @seealso \code{\link{plot_strata}}
#' @seealso \code{\link{make_comparison_matrix}}
#' 
#' @export
#' 
compare_SMCs <- function(in_stratification_lists_list,
                         in_signatures_ind_df,
                         output_path,
                         in_nrect=5,
                         in_attribute="") {
  strata_list <- make_strata_df(in_stratification_lists_list)
  strata_df <- strata_list$strata_df
  plot_strata(strata_list,in_signatures_ind_df,output_path,in_attribute)
  reduced_strata_df <- strata_df
  names(reduced_strata_df)[1] <- "all"
  remove_ind <- grep("_all$",names(reduced_strata_df))
  reduced_strata_df <- reduced_strata_df[,-remove_ind]
  my_matrix <- make_comparison_matrix(reduced_strata_df,output_path,
                                      in_nrect,in_attribute)
  return(my_matrix)
}

#' Wrapper function for \code{plot_strata}
#'
#' First calls
#' \itemize{
#'  \item \code{make_strata_df}, then
#'  \item \code{\link{plot_strata}}
#' }
#'
#' @param in_stratification_lists_list
#'  List of lists with entries from different (orthogonal) stratification 
#'  axes or SMCs
#' @param in_signatures_ind_df
#'  A data frame containing meta information about the signatures
#' @param output_path
#'  Path to directory where the results, especially the figure produced by
#'  \code{\link{plot_strata}} is going to be stored.
#' @param in_attribute
#'  Additional string for the file name where the figure produced by
#'  \code{\link{plot_strata}} is going to be stored.
#' @param in_remove_signature_ind
#'  Omit one of the signatures in \code{in_signatures_ind_df} for the 
#'  comparison if non-NULL. The parameter specifies the index of the
#'  signature to be removed.
#' @param in_additional_stratum
#'  Include an additionally supplied stratum in comparison in non-NULL.
#'  
#' @return The function doesn't return any value.
#' 
#' @examples
#'  NULL
#'  
#' @seealso \code{\link{plot_strata}}
#' 
#' @export
#' 
run_plot_strata_general <- function(in_stratification_lists_list,
                                    in_signatures_ind_df,
                                    output_path=NULL,
                                    in_attribute="",
                                    in_remove_signature_ind=NULL,
                                    in_additional_stratum=NULL) {
  strata_list <- make_strata_df(in_stratification_lists_list,
                                in_remove_signature_ind,
                                in_additional_stratum)
  plot_strata(strata_list$strata_df,in_signatures_ind_df,
              output_path,in_attribute)
}


#' Compare all strata from different stratifications
#'
#' Compare all strata from different orthogonal stratification axes, i.e.
#' othogonal SMCs by cosine similarity of signature exposures. Function 
#' similar to \code{\link{compare_SMCs}}, but without calling
#' \code{\link{plot_strata}}. First calls
#' \itemize{
#'  \item \code{make_strata_df}, then
#'  \item \code{\link{make_comparison_matrix}}
#' }
#'
#' @param in_stratification_lists_list
#'  List of lists with entries from different (orthogonal) stratification 
#'  axes or SMCs
#' @param output_path
#'  Path to directory where the results, especially the figure produced by
#'  \code{\link[corrplot]{corrplot}} is going to be stored.
#' @param in_nrect
#'  Number of clusters in the clustering procedure provided by
#'  \code{\link[corrplot]{corrplot}}
#' @param in_attribute
#'  Additional string for the file name where the figure produced by
#'  \code{\link[corrplot]{corrplot}} is going to be stored.
#' @param in_remove_signature_ind
#'  Omit one of the signatures in \code{in_signatures_ind_df} for the 
#'  comparison if non-NULL. The parameter specifies the index of the
#'  signature to be removed.
#' @param in_additional_stratum
#'  Include an additionally supplied stratum in comparison in non-NULL.
#'  
#' @return The comparison matrix of cosine similarities.
#' 
#' @examples
#'  NULL
#'  
#' @seealso \code{\link{make_comparison_matrix}}
#' @seealso \code{\link{compare_SMCs}}
#' @seealso \code{\link{run_comparison_catalogues}}
#' 
#' @export
#' 
run_comparison_general <- function(in_stratification_lists_list,
                                   output_path=NULL,
                                   in_nrect=5,
                                   in_attribute="",
                                   in_remove_signature_ind=NULL,
                                   in_additional_stratum=NULL) {
  strata_list <- make_strata_df(in_stratification_lists_list,
                                in_remove_signature_ind,
                                in_additional_stratum)
  strata_df <- strata_list$strata_df
  reduced_strata_df <- strata_df
  names(reduced_strata_df)[1] <- "all"
  remove_ind <- grep("_all$",names(reduced_strata_df))
  if(length(remove_ind)>0) {
    reduced_strata_df <- reduced_strata_df[,-remove_ind]
  }
  my_matrix <- make_comparison_matrix(reduced_strata_df,
                                      output_path,
                                      in_nrect,
                                      in_attribute)
  return(my_matrix)
}

#' Compare all strata from different stratifications
#'
#' Compare all strata from different orthogonal stratification axes, i.e.
#' othogonal SMCs by cosine similarity of mutational catalogues. Function 
#' similar to \code{\link{run_comparison_general}}. First calls
#' \itemize{
#'  \item \code{make_catalogue_strata_df}, then
#'  \item \code{\link{make_comparison_matrix}}
#' }
#'
#' @param in_stratification_lists_list
#'  List of lists with entries from different (orthogonal) stratification 
#'  axes or SMCs
#' @param output_path
#'  Path to directory where the results, especially the figure produced by
#'  \code{\link[corrplot]{corrplot}} is going to be stored.
#' @param in_nrect
#'  Number of clusters in the clustering procedure provided by
#'  \code{\link[corrplot]{corrplot}}
#' @param in_attribute
#'  Additional string for the file name where the figure produced by
#'  
#' @return The comparison matrix of cosine similarities.
#' 
#' @examples
#'  NULL
#'  
#' @seealso \code{\link{make_comparison_matrix}}
#' @seealso \code{\link{run_comparison_general}}
#' 
#' @export
#' 
run_comparison_catalogues <- function(in_stratification_lists_list,
                                      output_path=NULL,
                                      in_nrect=5,
                                      in_attribute="") {
  catalogue_strata_list <- 
    make_catalogue_strata_df(in_stratification_lists_list)
  catalogue_strata_df <- catalogue_strata_list$strata_df
  reduced_strata_df <- catalogue_strata_df
  names(reduced_strata_df)[1] <- "all"
  remove_ind <- grep("_all$",names(reduced_strata_df))
  if(length(remove_ind)>0) {
    reduced_strata_df <- reduced_strata_df[,-remove_ind]
  }
  my_matrix <- make_comparison_matrix(reduced_strata_df,
                                      output_path,
                                      in_nrect,
                                      in_attribute)
  return(my_matrix)
}
