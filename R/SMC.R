#' Wrapper function for the Stratification of a Mutational Catalogue
#'
#'\code{run_SMC} takes as input a big dataframe constructed from a vcf-like 
#'file of a whole cohort. This wrapper function calls custom functions to
#'construct a mutational catalogue and stratify it according to categories
#'indicated by a special column in the input dataframe:
#'\itemize{
#'  \item \code{\link{create_mutation_catalogue_from_df}}
#'  \item \code{\link{stratify_and_create_mutational_catalogue}}
#'  \item \code{\link{adjust_number_of_columns_in_list_of_catalogues}}
#'}
#'This stratification
#'yields a collection of stratified mutational catalogues, these are
#'reformatted and sent to the custom function \code{\link{SMC}} and thus
#'indirectly to \code{\link{LCD_SMC}} to perform a signature analysis
#'of the stratified mutational catalogues. The result is then handed over
#'to \code{\link{plot_SMC}} for visualization.
#'
#' @param my_table
#'  A big dataframe constructed from a vcf-like file of a whole cohort. The 
#'  first columns are those of a standard vcf file, followed by an arbitrary
#'  number of custom or user defined columns. One of these must carry a PID
#'  (patient or sample identifyier) and one must be the category used for
#'  stratification.
#' @param this_signatures_df
#'  A numeric data frame \code{W} in with \code{n} rows and \code{l} columns,
#'  \code{n} being the number of features and \code{l} being the number of
#'  signatures
#' @param this_signatures_ind_df
#'  A data frame containing meta information about the signatures
#' @param this_subgroups_df
#'  A data frame indicating which PID (patient or sample identifyier) belongs
#'  to which subgroup
#' @param column_name
#'  Name of the column in \code{my_table} which is going to be used for
#'  stratification
#' @param refGenome
#'  FaFile of the reference genome to extract the motif context of the
#'  variants in \code{my_table}
#' @param cohort_method_flag
#'  Either or several of \code{c("all_PIDs","cohort","norm_PIDs")}, representing
#'  alternative ways to average over the cohort.
#' @param in_strata_order_ind
#'  Index vector defining reordering of the strata
#' @param wordLength
#'  Integer number defining the length of the features or motifs, e.g. 3 for
#'  tripletts or 5 for pentamers
#' @param verbose_flag
#'  Verbose if \code{verbose_flag=1}
#' @param target_dir
#'  Path to directory where the results of the stratification procedure are
#'  going to be stored if non-NULL.
#' @param strata_dir
#'  Path to directory where the mutational catalogues of the different strata
#'  are going to be stored if non-NULL
#' @param output_path
#'  Path to directory where the results, especially the figures produced by
#'  \code{\link{plot_SMC}} are going to be stored.
#' @param in_all_exposures_df
#'  Optional argument, if specified, \code{H}, i.e. the overall exposures
#'  without stratification, is set to equal \code{in_all_exposures_df}. This 
#'  is equivalent to forcing the LCD_SMC procedure to use e.g. the exposures 
#'  of a previously performed NMF decomposition.
#' @param in_rownames
#'  Optional parameter to specify rownames of the mutational catalogue \code{V}
#'  i.e. the names of the features.
#' @param in_norms
#'  If specified, vector of the correction factors for every motif due to
#'  differing trinucleotide content. If null, no correction is applied.
#' @param in_label_orientation
#'  Whether or not to turn the labels on the x-axis.
#' @param this_sum_ind
#'  Optional set of indices for reordering the PIDs
#'   
#' @return A list with entries
#'  \code{exposures_list},
#'  \code{catalogues_list},
#'  \code{cohort} and
#'  \code{name_list}.
#' \itemize{
#'  \item \code{exposures_list}:
#'    The list of \code{s} strata specific exposures Hi, all are numerical data 
#'    frames with \code{l} rows and \code{m} columns, \code{l} being the number
#'    of signatures and \code{m} being the number of samples
#'  \item \code{catalogues_list}:
#'    A list of \code{s} strata specific cohortwide (i.e. averaged over cohort)
#'    normalized exposures
#'  \item \code{cohort}:
#'    \code{subgroups_df} adjusted for plotting
#'  \item \code{name_list}:
#'    Names of the contructed strata.
#' }
#' 
#' @examples
#' NULL
#' \dontrun{
#'  data(sigs)
#'  data(lymphoma_test)
#'  data(lymphoma_cohort_LCD_results)
#'  strata_list <- cut_breaks_as_intervals(lymphoma_test_df$random_norm,
#'                                         in_outlier_cutoffs=c(-4,4),
#'                                         in_cutoff_ranges_list=list(c(-2.5,-1.5),c(0.5,1.5)),
#'                                         in_labels=c("small","intermediate","big"))
#'  lymphoma_test_df$random_cat <- strata_list$category_vector
#'  choice_ind <- (names(lymphoma_Nature2013_COSMIC_cutoff_exposures_df) 
#'                 %in% unique(lymphoma_test_df$PID))
#'  lymphoma_test_exposures_df <- lymphoma_Nature2013_COSMIC_cutoff_exposures_df[,choice_ind]
#'  temp_subgroups_df <- make_subgroups_df(lymphoma_test_exposures_df,lymphoma_test_df)
#'  lymphoma_test_signatures_df <- Alex_COSMIC_signatures_df[,chosen_signatures_indices_df$index]
#'  mut_density_list <- run_SMC(lymphoma_test_df,
#'                              lymphoma_test_signatures_df,
#'                              chosen_signatures_indices_df,
#'                              temp_subgroups_df,
#'                              column_name="random_cat",
#'                              refGenome=BSgenome.Hsapiens.UCSC.hg19,
#'                              cohort_method_flag="norm_PIDs",
#'                              in_rownames = rownames(Alex_COSMIC_signatures_df))
#' }
#' 
#' @seealso \code{\link{create_mutation_catalogue_from_df}}
#' @seealso \code{\link{stratify_and_create_mutational_catalogue}}
#' @seealso \code{\link{adjust_number_of_columns_in_list_of_catalogues}}
#' @seealso \code{\link{normalizeMotifs_otherRownames}}
#' @seealso \code{\link{SMC}}
#' @seealso \code{\link{plot_SMC}}
#' @seealso \code{\link{LCD_SMC}}
#' 
#' @export
#' 
run_SMC <- function(my_table,
                    this_signatures_df,
                    this_signatures_ind_df,
                    this_subgroups_df,
                    column_name,
                    refGenome,
                    cohort_method_flag="all_PIDs",
                    in_strata_order_ind=seq_len(length(unique(my_table[,column_name]))),
                    wordLength=3,verbose_flag=1,
                    target_dir=NULL,strata_dir=NULL,output_path=NULL,
                    in_all_exposures_df=NULL,
                    in_rownames = c(),
                    in_norms=NULL,
                    in_label_orientation="turn",
                    this_sum_ind=NULL) {
  
  refGenome_Seqinfo <- seqinfo(refGenome)
  if (verbose_flag==1) {cat("\nYAPSA:::run_SMC::calling create_mutation_catalogue_from_df...\n");}
  merged_results_list <- create_mutation_catalogue_from_df(my_table,
                                                          refGenome_Seqinfo,
                                                          this_seqnames.field="CHROM",
                                                          this_start.field="POS",
                                                          this_end.field="POS",
                                                          this_PID.field="PID",
                                                          this_subgroup.field="subgroup",
                                                          this_rownames = in_rownames,
                                                          refGenome,wordLength,verbose_flag)
  temp_rownames <- rownames(merged_results_list$matrix)
  if (verbose_flag==1) {cat("YAPSA:::run_SMC::calling stratify_and_create_mutational_catalogue...\n");}
  all_list <- stratify_and_create_mutational_catalogue(my_table,column_name,target_dir,strata_dir,
                                                       refGenome_Seqinfo,
                                                       our_seqnames.field="CHROM",
                                                       our_start.field="POS",
                                                       our_end.field="POS",
                                                       our_PID.field="PID",
                                                       our_subgroup.field="subgroup",
                                                       refGenome,wordLength,verbose_flag,
                                                       temp_rownames)
  
  name_list <- all_list$name_list
  df_list <- adjust_number_of_columns_in_list_of_catalogues(all_list,merged_results_list)
  mutation_catalogue_all_df <- as.data.frame(merged_results_list$matrix)
  ## adjust for trinucleotide content if necessary
  if(!is.null(in_norms)) {
    if (verbose_flag==1) {cat("\nYAPSA:::run_SMC::adapting to different kmer distribution by calling normalizeMotifs_otherRownames...\n");}
    df_list <- lapply(df_list, function(l) normalizeMotifs_otherRownames(l,in_norms))
    mutation_catalogue_all_df <- normalizeMotifs_otherRownames(mutation_catalogue_all_df,in_norms)
  }
  number_of_strata <- length(df_list)
  number_of_sigs <- dim(this_signatures_df)[2]
  
  if (verbose_flag==1) {cat("\nYAPSA:::run_SMC::calling SMC...\n");}
  SMC_list <- SMC(df_list,this_signatures_df,in_all_exposures_df,number_of_strata,number_of_sigs,
                  name_list,this_subgroups_df,mutation_catalogue_all_df,cohort_method_flag,
                  in_verbose=verbose_flag)
  exposures_strata_list <- SMC_list$exposures_strata_list
  exposures_both_rel_df_list <- SMC_list$exposures_both_rel_df_list
  this_subgroups_df <- SMC_list$this_subgroups_df
  decomposition_method <- SMC_list$decomposition_method
  
  ## plot
  if (verbose_flag==1) {cat("YAPSA:::run_SMC::calling plot_SMC...\n");}
  my_plot_list <- plot_SMC(number_of_strata,output_path,decomposition_method,
                           number_of_sigs,name_list,exposures_strata_list,
                           this_signatures_ind_df,this_subgroups_df,
                           in_strata_order_ind,exposures_both_rel_df_list,
                           cohort_method_flag,
                           in_label_orientation=in_label_orientation,
                           this_sum_ind=this_sum_ind)
  return(list(exposures_list=exposures_strata_list,catalogues_list=df_list,
              cohort=exposures_both_rel_df_list,name_list=name_list))
}


#' Stratification of a Mutational Catalogue
#'
#'\code{SMC} takes a given collection of stratified mutational catalogues \code{Vi},
#'sends them to perform a mutational signatures decomposition by Linear
#'Combination Decomposition (LCD) with the functions \code{\link{LCD_SMC}}
#'with known signatures \code{W}. It subsequently performs some useful
#'statistics and preparation for plotting with the function \code{\link{plot_SMC}}.
#'\code{SMC} is naturally called by \code{\link{run_SMC}}.
#'
#' @param df_list
#'  A list of \code{s} stratified mutational catalogues \code{Vi} \(numeric data frames\)
#'  with \code{n} rows and \code{m} columns each, \code{n} being 
#'  the number of features and \code{m} being the number of samples.
#'  This list is naturally provided in \code{\link{run_SMC}}.
#' @param this_signatures_df
#'  A numeric data frame \code{W} in with \code{n} rows and \code{l} columns,
#'  \code{n} being the number of features and \code{l} being the number of
#'  signatures
#' @param in_all_exposures_df
#'  The overall exposures \code{H} without stratification, a numeric data frame with 
#'  \code{l} rows and \code{m} columns, \code{l} being
#'  the number of signatures and \code{m} being the number
#'  of samples
#' @param number_of_strata
#'  The length of the list \code{df_list}
#' @param number_of_sigs
#'  The number of signatures used in the current decomposition.
#' @param name_list
#'  A list of names of the different strata
#'@param this_subgroups_df
#'  A data frame indicating which PID (patient or sample identifyier) belongs
#'  to which subgroup
#' @param mutation_catalogue_all_df
#'  The overall mutational catalogue \code{V} without stratification.
#' @param cohort_method_flag
#'  Either or several of \code{c("all_PIDs","cohort","norm_PIDs")}, representing
#'  alternative ways to average over the cohort.
#' @param in_verbose
#'  Verbose if \code{in_verbose=1}
#'  
#' @return A list with entries
#'  \code{exposures_strata_list},
#'  \code{exposures_both_rel_df_list},
#'  \code{this_subgroups_df},
#'  \code{subgroup_ind} and
#'  \code{decomposition_method}.
#' \itemize{
#'  \item \code{exposures_strata_list}:
#'    The list of \code{s} strata specific exposures Hi, all are numerical data 
#'    frames with \code{l} rows and \code{m} columns, \code{l} being the number
#'    of signatures and \code{m} being the number of samples
#'  \item \code{exposures_both_rel_df_list}:
#'    A list of \code{s} strata specific cohortwide (i.e. averaged over cohort)
#'    normalized exposures
#'  \item \code{this_subgroups_df}:
#'    \code{subgroups_df} adjusted for plotting
#'  \item \code{subgroup_ind}:
#'    Index of the subgroups chosen and relevant for plotting.
#'  \item \code{decomposition_method}:
#'    String telling whether LCD or NMF was used, relevant only for handing over
#'    to \code{\link{plot_SMC}}.
#' }
#' 
#' @examples
#' NULL
#' 
#' @seealso \code{\link{run_SMC}}
#' @seealso \code{\link{plot_SMC}}
#' @seealso \code{\link{LCD_SMC}}
#' 
#' @import reshape2
#' @export
#' 
SMC <- function(df_list,this_signatures_df,in_all_exposures_df,number_of_strata,number_of_sigs,
                name_list,this_subgroups_df,mutation_catalogue_all_df,cohort_method_flag,in_verbose=1) {
  ## this is the LCD_SMC procedure
  if (in_verbose==1) {cat("YAPSA:::SMC::calling LCD_SMC on per-PID data...\n");}
  exposures_strata_list <- LCD_SMC(df_list,this_signatures_df,in_all_exposures_df)
  exposures_all_df <- exposures_strata_list$exposures_all_df
  
  sum_over_exposures_all_per_pid <- apply(exposures_all_df,2,sum)
  sum_df <- data.frame(sum=sum_over_exposures_all_per_pid,PID=colnames(exposures_all_df))
  
  ## make cohort-wide statistics
  ## 1. compute exposures in cohort as sum over PIDs
  mean_over_exposures_all_per_sig <- apply(exposures_all_df,1,mean)
  stderrmean_over_exposures_all_per_sig <- apply(exposures_all_df,1,stderrmean)
  exposures_all_PIDs_df <- data.frame(all_abs=mean_over_exposures_all_per_sig,
                                      all_stderrmean=stderrmean_over_exposures_all_per_sig)
  temp_denominator <- sum(exposures_all_PIDs_df$all_abs)
  exposures_all_PIDs_df$all_rel <- exposures_all_PIDs_df$all_abs/temp_denominator
  exposures_all_PIDs_df$all_rstderrmean <- exposures_all_PIDs_df$all_stderrmean/temp_denominator
  for (i in seq_len(number_of_strata)) {
    my_index <- 4*(i)+1
    temp_exposures_df <- exposures_strata_list$sub_exposures_list[[i]]
    exposures_all_PIDs_df[,my_index] <- apply(temp_exposures_df,1,mean)
    exposures_all_PIDs_df[,my_index+1] <- apply(temp_exposures_df,1,stderrmean)
    temp_denominator <- sum(exposures_all_PIDs_df[,my_index])
    exposures_all_PIDs_df[,my_index+2] <- exposures_all_PIDs_df[,my_index]/temp_denominator
    exposures_all_PIDs_df[,my_index+3] <- exposures_all_PIDs_df[,my_index+1]/temp_denominator
    names(exposures_all_PIDs_df)[my_index] <- paste0(name_list[[i]],"_abs")
    names(exposures_all_PIDs_df)[my_index+1] <- paste0(name_list[[i]],"_stderrmean")
    names(exposures_all_PIDs_df)[my_index+2] <- paste0(name_list[[i]],"_rel")
    names(exposures_all_PIDs_df)[my_index+3] <- paste0(name_list[[i]],"_rstderrmean")
  }
  exposures_all_PIDs_df$sig <- rownames(exposures_all_PIDs_df)
  exposures_all_PIDs_df$method <- "all_PIDs"
  ## 2. compute exposures in cohort by running decomposition on a fused vector
  if (in_verbose==1) {cat("YAPSA:::SMC::compute exposures in cohort by running decomposition on a fused vector...\n");}
  catalogue_in_cohort_all <- apply(mutation_catalogue_all_df,1,sum)
  catalogue_in_cohort_df <- data.frame(all_abs=catalogue_in_cohort_all)
  catalogue_in_cohort_df$all_rel <- catalogue_in_cohort_df$all_abs/sum(catalogue_in_cohort_df$all_abs)
  catalogue_in_cohort_df_list <- list()
  for (i in seq_len(number_of_strata)) {
    my_index <- 2*(i)+1
    temp_catalogue_df <- df_list[[i]]
    catalogue_in_cohort_df[,my_index] <- apply(temp_catalogue_df,1,sum)
    catalogue_in_cohort_df[,my_index+1] <- catalogue_in_cohort_df[,my_index]/sum(catalogue_in_cohort_df[,my_index])
    names(catalogue_in_cohort_df)[my_index] <- paste0(name_list[[i]],"_abs")
    names(catalogue_in_cohort_df)[my_index+1] <- paste0(name_list[[i]],"_rel")
    catalogue_in_cohort_df_list[[i]] <- data.frame(all_abs=catalogue_in_cohort_df[,my_index])
  }
  catalogue_in_cohort_in_exposures <- NULL
  decomposition_method <- "LCD"
  if (!is.null(in_all_exposures_df)) {
    temp_vector <- apply(in_all_exposures_df,1,sum)
    catalogue_in_cohort_in_exposures <- data.frame(all_abs=temp_vector)
    decomposition_method <- "NMF"
  }
  if (in_verbose==1) {cat("YAPSA:::SMC::calling LCD_SMC on cohort-wide data...\n");}
  exposures_in_cohort_strata_list <- LCD_SMC(catalogue_in_cohort_df_list,this_signatures_df,catalogue_in_cohort_in_exposures)
  exposures_in_cohort_df <- data.frame(all_abs=exposures_in_cohort_strata_list$exposures_all_df,
                                       all_stderrmean=rep(0,dim(exposures_in_cohort_strata_list$exposures_all_df)[1]))
  exposures_in_cohort_df$all_rel <- exposures_in_cohort_df$all_abs/sum(exposures_in_cohort_df$all_abs)
  exposures_in_cohort_df$all_rstderrmean <- rep(0,dim(exposures_in_cohort_strata_list$exposures_all_df)[1])
  for (i in seq_len(number_of_strata)) {
    my_index <- 4*(i)+1
    temp_exposures_df <- exposures_in_cohort_strata_list$sub_exposures_list[[i]]
    exposures_in_cohort_df[,my_index] <- temp_exposures_df$all_abs
    exposures_in_cohort_df[,my_index+1] <- 0
    exposures_in_cohort_df[,my_index+2] <- exposures_in_cohort_df[,my_index]/sum(exposures_in_cohort_df[,my_index])
    exposures_in_cohort_df[,my_index+3] <- 0
    names(exposures_in_cohort_df)[my_index] <- paste0(name_list[[i]],"_abs")
    names(exposures_in_cohort_df)[my_index+1] <- paste0(name_list[[i]],"_stderrmean")
    names(exposures_in_cohort_df)[my_index+2] <- paste0(name_list[[i]],"_rel")
    names(exposures_in_cohort_df)[my_index+3] <- paste0(name_list[[i]],"_rstderrmean")
  }
  exposures_in_cohort_df$sig <- rownames(exposures_in_cohort_df)
  exposures_in_cohort_df$method <- "cohort"
  ## 3. average over relative exposures
  ## 3.a) compute relative exposures
  exposures_strata_list$norm_exposures_all_df <- normalize_df_per_dim(exposures_strata_list$exposures_all_df,2)
  exposures_strata_list$sub_norm_exposures_list <- lapply(exposures_strata_list$sub_exposures_list,
                                                         function(l) normalize_df_per_dim(l,2))
  ## 3.b) average and build data structure
  mean_over_norm_exposures_all_per_sig <- apply(exposures_strata_list$norm_exposures_all_df,1,mean)
  stderrmean_over_norm_exposures_all_per_sig <- apply(exposures_strata_list$norm_exposures_all_df,1,stderrmean)
  norm_exposures_all_PIDs_df <- data.frame(all_abs=rep(0,length(mean_over_norm_exposures_all_per_sig)),
                                           all_stderrmean=rep(0,length(mean_over_norm_exposures_all_per_sig)),
                                           all_rel=mean_over_norm_exposures_all_per_sig,
                                           all_rstderrmean=stderrmean_over_norm_exposures_all_per_sig)
  temp_list <- lapply(exposures_strata_list$sub_norm_exposures_list,function(l) apply(l,1,mean))
  for (i in seq_len(number_of_strata)) {
    my_index <- 4*(i)+1
    norm_exposures_all_PIDs_df[,my_index] <- 0
    norm_exposures_all_PIDs_df[,my_index+1] <- 0
    norm_exposures_all_PIDs_df[,my_index+2] <- average_over_present(exposures_strata_list$sub_norm_exposures_list[[i]],1)
    norm_exposures_all_PIDs_df[,my_index+3] <- stderrmean_over_present(exposures_strata_list$sub_norm_exposures_list[[i]],1)
    names(norm_exposures_all_PIDs_df)[my_index] <- paste0(name_list[[i]],"_abs")
    names(norm_exposures_all_PIDs_df)[my_index+1] <- paste0(name_list[[i]],"_stderrmean")
    names(norm_exposures_all_PIDs_df)[my_index+2] <- paste0(name_list[[i]],"_rel")
    names(norm_exposures_all_PIDs_df)[my_index+3] <- paste0(name_list[[i]],"_rstderrmean")
  }
  norm_exposures_all_PIDs_df$sig <- rownames(exposures_all_PIDs_df)
  norm_exposures_all_PIDs_df$method <- "norm_PIDs"
  ## 4 unite the 3 different methods into one dataframe
  exposures_combMethod_df <- rbind(exposures_all_PIDs_df,exposures_in_cohort_df,norm_exposures_all_PIDs_df)
  abs_ind <- grep(".*_abs|^sig$|^method$",names(exposures_combMethod_df))
  rel_ind <- grep(".*_rel|^sig$|^method$",names(exposures_combMethod_df))
  rstderrmean_ind <- grep(".*_rstderrmean|^sig$|^method$",names(exposures_combMethod_df))
  exposures_combMethod_abs_df <- exposures_combMethod_df[,abs_ind]
  exposures_combMethod_rel_df <- exposures_combMethod_df[,rel_ind]
  exposures_combMethod_rstderrmean_df <- exposures_combMethod_df[,rstderrmean_ind]
  
  ## prepare for plotting
  if (in_verbose==1) {cat("YAPSA:::SMC::prepare for plotting...\n");}
  exposures_combMethod_rel_melt_df <- melt(exposures_combMethod_rel_df,id.vars=c("sig","method"),value.name="exposure")
  exposures_combMethod_rstderrmean_melt_df <- melt(exposures_combMethod_rstderrmean_df,id.vars=c("sig","method"),value.name="exposure_rstderrmean")
  exposures_combMethod_rel_melt_df$rstderrmean <- exposures_combMethod_rstderrmean_melt_df$exposure_rstderrmean
  exposures_combMethod_rel_melt_df$exposure_min <- exposures_combMethod_rel_melt_df$exposure - exposures_combMethod_rstderrmean_melt_df$exposure_rstderrmean
  exposures_combMethod_rel_melt_df$exposure_max <- exposures_combMethod_rel_melt_df$exposure + exposures_combMethod_rstderrmean_melt_df$exposure_rstderrmean
  counter <- 0
  exposures_combMethod_rel_df_list <- list()
  exposures_combMethod_rstderrmean_df_list <- list()
  method_choice_ind <- which(exposures_combMethod_rel_melt_df$method %in% cohort_method_flag)
  exposures_combMethod_rel_melt_df <- exposures_combMethod_rel_melt_df[method_choice_ind,]
  for (temp_varname in unique(exposures_combMethod_rel_melt_df$variable)) {
    counter <- counter + 1
    temp_ind <- which(exposures_combMethod_rel_melt_df$variable==temp_varname)
    exposures_combMethod_rel_df_list[[counter]] <- exposures_combMethod_rel_melt_df[temp_ind,]
    exposures_combMethod_rel_df_list[[counter]]$sig <- factor(exposures_combMethod_rel_df_list[[counter]]$sig,
                                                              levels=unique(exposures_combMethod_rel_df_list[[counter]]$sig))
  }  

  ## adapt subgroups data structure
  if(dim(mutation_catalogue_all_df)[2]==dim(this_subgroups_df)[1]) {
    this_subgroups_df$PID <- colnames(mutation_catalogue_all_df)
    this_subgroups_df$sum <- sum_df$sum
    max_total_count <- max(this_subgroups_df$sum)
    this_subgroups_df$compl_sum <- max_total_count - this_subgroups_df$sum
    subgroup_ind <- order(this_subgroups_df$subgroup,this_subgroups_df$compl_sum)
    this_subgroups_df$index <- order(subgroup_ind)
  } else {
    cat("YAPSA:::SMC::Warning: dimension mismatch between mutation_catalogue_df and subgroup_df.\n")
    q(status=5)
  }
  
  return(list(exposures_strata_list=exposures_strata_list,
              exposures_both_rel_df_list=exposures_combMethod_rel_df_list,
              this_subgroups_df=this_subgroups_df,
              subgroup_ind=subgroup_ind,
              decomposition_method=decomposition_method))
}


#' Apply statistical tests to a stratification (SMC)
#'
#' \code{stat_test_SMC} tests for enrichment or depletion in the different strata
#' of a stratification of the mutational catalogue for every signature
#' independently by applying Kruskal Wallis tests. For those signatures where the
#' Kruskal Wallis test gives a significant p-value, pairwise posthoc tests are 
#' carried out by calling \code{\link[PMCMR]{posthoc.kruskal.nemenyi.test}}.
#' Additionally all data is tested for normality by Shapiro Wilk tests, so that
#' the user may apply ANOVA and pairwise posthoc t-test where allowed.
#'
#' @param in_strat_list
#'  A list with entries \code{exposures_list}, \code{catalogues_list},
#'  \code{cohort} and \code{name_list} as in the output of \code{\link{run_SMC}}.
#'  \itemize{
#'    \item \code{exposures_list}:
#'      The list of \code{s} strata specific exposures Hi, all are numerical data 
#'      frames with \code{l} rows and \code{m} columns, \code{l} being the number
#'      of signatures and \code{m} being the number of samples
#'    \item \code{catalogues_list}:
#'      A list of \code{s} strata specific cohortwide (i.e. averaged over cohort)
#'      normalized exposures
#'    \item \code{cohort}:
#'      \code{subgroups_df} adjusted for plotting
#'    \item \code{name_list}:
#'      Names of the contructed strata.
#'  }
#' @param in_flag
#'  If "norm", all tests are performed on normalized exposures, otherwise the
#'  absolute exposures are taken.
#'  
#' @return A list with entries
#'  \code{kruskal_df},
#'  \code{shapiro_df},
#'  \code{kruskal_posthoc_list},
#' \itemize{
#'  \item \code{kruskal_df}:
#'    A data frame containing results (statistic and p values) of the Kruskal Wallis
#'    tests (tests for enrichment or depletion in the different strata for every
#'    signature independently).
#'  \item \code{shapiro_df}:
#'    A data frame containing results (p values) of the Shapiro Wilk tests
#'    (tests for normal distribution in the different strata for every
#'    signature independently).
#'  \item \code{kruskal_posthoc_list}:
#'    A list of results of pairwise posthoc tests carried out for those signatures
#'    where the Kruskal Wallis test yielded a significant p-value (carried out by
#'    \code{\link[PMCMR]{posthoc.kruskal.nemenyi.test}}).
#' }
#' 
#' @seealso \code{\link{run_SMC}}
#' @seealso \code{\link{SMC}}
#' @seealso \code{\link[PMCMR]{posthoc.kruskal.nemenyi.test}}
#' @seealso \code{\link[stats]{kruskal.test}}
#' @seealso \code{\link{shapiro_if_possible}}
#' @seealso \code{\link[stats]{shapiro.test}}
#' 
#' @import PMCMR
#' @export
#' 
stat_test_SMC <- function(in_strat_list,in_flag="norm"){
  my_number_of_strata <- length(in_strat_list$exposures_list$sub_norm_exposures_list)
  my_signatures <- rownames(in_strat_list$exposures_list$sub_norm_exposures_list[[1]])
  my_number_of_signatures <- dim(in_strat_list$exposures_list$sub_norm_exposures_list[[1]])[1]
  my_number_of_PIDs <- dim(in_strat_list$exposures_list$sub_norm_exposures_list[[1]])[2]
  out_stat_df <- repeat_df(NA,my_number_of_signatures,3)
  shapiro_df <- repeat_df(NA,my_number_of_signatures,my_number_of_strata)
  rownames(out_stat_df) <- my_signatures
  rownames(shapiro_df) <- my_signatures
  names(out_stat_df) <- c("Kruskal_statistic","df","Kruskal_p_val")
  posthoc_list <- list()
  for(temp_sig in my_signatures){
    sig_exposures_vector <- c()
    sig_strata_vector <- c()
    for(j in seq_len(my_number_of_strata)){
      if(in_flag=="norm"){
        temp_df <- in_strat_list$exposures_list$sub_norm_exposures_list[[j]]
      } else {
        temp_df <- in_strat_list$exposures_list$sub_exposures_list[[j]]
      }
      PID_choice_ind <- which(colSums(temp_df)>0)
      temp_exposures_vector <- as.numeric(temp_df[temp_sig,PID_choice_ind])
      sig_exposures_vector <- c(sig_exposures_vector,temp_exposures_vector)
      sig_strata_vector <- c(sig_strata_vector,
          rep(in_strat_list$name_list[[j]],length(temp_exposures_vector)))
      ## catch exception for temp_exposures_vector being all equal values (all zero)
      #shapiro_df[temp_sig,j] <- shapiro.test(temp_exposures_vector)$p.value
      shapiro_df[temp_sig,j] <- shapiro_if_possible(temp_exposures_vector)
    }
    sig_strata_vector <- factor(sig_strata_vector)
    sig_kruskal <- kruskal.test(sig_exposures_vector,sig_strata_vector)
    posthoc_list[[temp_sig]] <- posthoc.kruskal.nemenyi.test(x=sig_exposures_vector,g=sig_strata_vector)
    out_stat_df[temp_sig,1] <- sig_kruskal$statistic
    out_stat_df[temp_sig,2] <- sig_kruskal$parameter
    out_stat_df[temp_sig,3] <- sig_kruskal$p.value
  }
  out_stat_df$Kruskal_p_val_BH <- p.adjust(out_stat_df$Kruskal_p_val, method="BH")
  return(list(kruskal_df=out_stat_df,shapiro_df=shapiro_df,kruskal_posthoc_list=posthoc_list))
}


#' Test for differences in average signature exposures between subgroups
#' 
#' Apply Kruskal-Wallis tests to detect differences in the signature exposures
#' between different subgroups. Uses \code{\link{split_exposures_by_subgroups}}.
#' Algorithm analogous to \code{\link{stat_test_SMC}}.
#' 
#' @param in_exposures_df
#'  Numerical data frame of the exposures (i.e. contributions of the
#'  different signatures to the number of point mutations per PID)
#' @param in_subgroups_df
#'  Data frame indicating which PID belongs to which subgroup
#' @param in_subgroups.field
#'  Name indicating which column in \code{in_subgroups_df} contains the
#'  subgroup information
#' @param in_PID.field
#'  Name indicating which column in \code{in_subgroups_df} contains the
#'  PID information 
#'  
#' @return A list with entries
#'  \code{kruskal_df},
#'  \code{kruskal_posthoc_list},
#' \itemize{
#'  \item \code{kruskal_df}:
#'    A data frame containing results (statistic and p values) of the Kruskal Wallis
#'    tests (tests for enrichment or depletion in the different strata for every
#'    signature independently).
#'  \item \code{kruskal_posthoc_list}:
#'    A list of results of pairwise posthoc tests carried out for those signatures
#'    where the Kruskal Wallis test yielded a significant p-value (carried out by
#'    \code{\link[PMCMR]{posthoc.kruskal.nemenyi.test}}).
#' }
#' 
#' @seealso \code{\link{split_exposures_by_subgroups}}
#' @seealso \code{\link{stat_test_SMC}}
#' @seealso \code{\link[PMCMR]{posthoc.kruskal.nemenyi.test}}
#' @seealso \code{\link[stats]{kruskal.test}}
#' 
#' @examples
#'  NULL
#'  
#' @export
#' 
stat_test_subgroups <- function(in_exposures_df,in_subgroups_df,
                                in_subgroups.field="subgroup",
                                in_PID.field="PID"){
  # split the data with custom function
  exposures_df_list <- split_exposures_by_subgroups(
    in_exposures_df=in_exposures_df,in_subgroups_df=in_subgroups_df,
    in_subgroups.field=in_subgroups.field,
    in_PID.field=in_PID.field)
  number_of_subgroups <- length(exposures_df_list)
  my_signatures <- rownames(in_exposures_df)
  number_of_signatures <- length(my_signatures)
  out_stat_df <- repeat_df(NA,number_of_signatures,3)
  shapiro_df <- repeat_df(NA,number_of_signatures,number_of_subgroups)
  rownames(out_stat_df) <- my_signatures
  rownames(shapiro_df) <- my_signatures
  names(out_stat_df) <- c("Kruskal_statistic","df","Kruskal_p_val")
  posthoc_list <- list()
  for(temp_sig in my_signatures){
    sig_exposures_vector <- c()
    sig_subgroups_vector <- c()
    for(j in seq_len(number_of_subgroups)){
      temp_exposures_vector <- as.numeric(exposures_df_list[[j]][temp_sig,])
      sig_exposures_vector <- c(sig_exposures_vector,temp_exposures_vector)
      sig_subgroups_vector <- c(sig_subgroups_vector,
                             rep(names(exposures_df_list)[j],
                                 length(temp_exposures_vector)))
    }
    sig_subgroups_vector <- factor(sig_subgroups_vector)
    sig_kruskal <- kruskal.test(sig_exposures_vector,sig_subgroups_vector)
    posthoc_list[[temp_sig]] <- posthoc.kruskal.nemenyi.test(x=sig_exposures_vector,g=sig_subgroups_vector)
    out_stat_df[temp_sig,1] <- sig_kruskal$statistic
    out_stat_df[temp_sig,2] <- sig_kruskal$parameter
    out_stat_df[temp_sig,3] <- sig_kruskal$p.value
  }
  out_stat_df$Kruskal_p_val_BH <- p.adjust(out_stat_df$Kruskal_p_val, method="BH")
  return(list(kruskal_df=out_stat_df,shapiro_df=shapiro_df,kruskal_posthoc_list=posthoc_list))
}