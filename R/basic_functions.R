#' Compute the cosine distance of two vectors
#'
#' @param a,b Numerical vectors of same length
#' @return The scalar product of the two input vectors divided by the
#'        product of the norms of the two input vectors
#' 
#' @examples
#' ## 1. Orthogonal vectors:
#' cosineDist(c(1,0),c(0,1))
#' ## 2. Non-orthogonal vectors:
#' cosineDist(c(1,0),c(1,1))
#' ## Compare trigonometry:
#' 1-cos(pi/4)
#' 
#' @export
#' 
cosineDist <- function(a,b){
  return(1 - sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) )) 
}


#' Wrapper function for \code{\link[grid]{viewport}}
#'
#' @param x,y
#'  Passed to \code{\link[grid]{viewport}} as \code{layout.pos.row} and
#'  \code{layout.pos.col}
#'  
#' @return The result of the \code{\link[grid]{viewport}} function.
#'
#' @examples
#'  NULL
#'
#' @import grid
#' @export
#' 
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)


#' Compare two sets of signatures by cosine distance
#'
#' Compare two sets of signatures, stored in numerical data frames
#' \code{W1} and \code{W2}, by computing the column-wise cosine distance
#'
#' @param in_df_small,in_df_big Numerical data frames \code{W1} and 
#'  \code{W2}, ideally the bigger one first, both with \code{n} rows
#'  and \code{l1} and \code{l2} columns, \code{n} being the number of
#'  features and \code{l1} and \code{l2} being the respective numbers
#'  of signatures of \code{W1} and \code{W2}
#'  
#' @return A list with entries
#'  \code{distance},
#'  \code{hierarchy_small} and
#'  \code{hierarchy_big}.
#' \itemize{
#'  \item \code{distance}:
#'    A numerical data frame with the cosine distances between the
#'    columns of \code{W1}, indexing the rows, and \code{W2}, indexing
#'    the columns
#'  \item \code{hierarchy_small}:
#'    A data frame carrying the information of ranked similarity between
#'    the signatures in \code{W2} with the signatures in \code{W1}
#'  \item \code{hierarchy_big}:
#'    A data frame carrying the information of ranked similarity between
#'    the signatures in \code{W1} with the signatures in \code{W2}
#' }
#' 
#' @seealso \code{\link{cosineDist}}
#' 
#' @examples
#' sig_1_df <- data.frame(matrix(c(1,0,0,0,0,1,0,0,0,0,1,0),ncol=3))
#' names(sig_1_df) <- paste0("B",seq_len(dim(sig_1_df)[2]))
#' sig_2_df <- data.frame(matrix(c(1,1,0,0,0,0,1,1),ncol=2))
#' compare_sets(sig_1_df,sig_2_df)
#' 
#' @export
#' 
compare_sets <- function(in_df_small,in_df_big) {
  sig_names_small <- colnames(in_df_small)
  sig_names_big <- colnames(in_df_big)
  number_of_sigs_small <- dim(in_df_small)[2]
  number_of_sigs_big <- dim(in_df_big)[2]
  number_of_features_small <- dim(in_df_small)[1]
  number_of_features_big <- dim(in_df_big)[1]
  if(number_of_features_small != number_of_features_big) {
    cat("compare_sets: Error: features mismatch")
    return()
  }
  temp_ind <- match(rownames(in_df_small),rownames(in_df_big))
  df_big <- in_df_big[temp_ind,]
  df_small <- in_df_small
  distance_df <- as.data.frame(matrix(rep(1,number_of_sigs_small*number_of_sigs_big),ncol=number_of_sigs_big))
  rownames(distance_df) <- sig_names_small
  colnames(distance_df) <- sig_names_big
  for (i in 1:number_of_sigs_small) {
    for (j in 1:number_of_sigs_big) {
      distance_df[i,j] <- cosineDist(df_small[,i],df_big[,j])
    }
  }  
  hierarchy_small_df <- t(apply(distance_df,1,function(l) return(colnames(distance_df)[order(l)])))
  hierarchy_big_df <- t(apply(distance_df,2,function(l) return(rownames(distance_df)[order(l)])))
  return(list(distance=distance_df,hierarchy_small=hierarchy_small_df,hierarchy_big=hierarchy_big_df))
}


#' Create a data frame with default values
#'
#' @param in_value
#'  Default entry to be repeated in the data frame
#' @param in_rows,in_cols
#'  Dimensions of the data frame to be created
#' @return The created data frame
#' 
#' @examples
#' ## 1. Initialize with numeric value:
#' repeat_df(1,2,3)
#' ## 2. Initialize with NA value:
#' repeat_df(NA,3,2)
#' ## 3. Initialize with character:
#' repeat_df("a",4,3)
#' 
#' @export
#' 
repeat_df <- function(in_value,in_rows,in_cols) {
  return(as.data.frame(matrix(rep(in_value,in_cols*in_rows),ncol=in_cols)))
}


#' Normalize a numerical data frame over a specified dimension
#'
#' Normalization is carried out by dividing by \code{rowSums} or \code{colSums};
#' for rows with \code{rowSums=0} or columns with \code{colSums=0}, the
#' normalization is left out.
#'
#' @param in_df
#'  Data frame to be normalized
#' @param in_dimension
#'  Dimension along which the normalization will be carried out
#' @return The normalized numerical data frame
#' 
#' @examples
#' test_df <- data.frame(matrix(c(1,2,3,0,5,2,3,4,0,6,0,0,0,0,0,4,5,6,0,7),ncol=4))
#' ## 1. Normalize over rows:
#' normalize_df_per_dim(test_df,1)
#' ## 2. Normalize over columns:
#' normalize_df_per_dim(test_df,2)
#' 
#' @export
#' 
normalize_df_per_dim <- function(in_df,in_dimension) {
  out_df <- repeat_df(in_value = 0,in_rows = dim(in_df)[1],in_cols = dim(in_df)[2])
  if(in_dimension==1) {
    choice_ind <- which(rowSums(in_df)>0)
    out_df[choice_ind,] <- in_df[choice_ind,]/rowSums(in_df)[choice_ind]      
  } else if (in_dimension==2) {
    t_df <- t(in_df)
    choice_ind <- which(rowSums(t_df)>0)
    temp_df <- repeat_df(in_value = 0,in_rows = dim(in_df)[2],in_cols = dim(in_df)[1])
    temp_df[choice_ind,] <- t_df[choice_ind,]/rowSums(t_df)[choice_ind]
    out_df <- as.data.frame(t(temp_df))
  } else {
    return(NULL)
  }
  colnames(out_df) <- colnames(in_df)
  rownames(out_df) <- rownames(in_df)
  return(out_df)
}


#' Average a data frame over a specified dimension
#'
#' If averaging over columns, zero rows (i.e. those with \code{rowSums=0})
#' are left out, if averaging over rows, zero columns (i.e. those with
#' \code{colSums=0}) are left out.
#'
#' @param in_df
#'  Data frame to be normalized
#' @param in_dimension
#'  Dimension along which the averaging will be carried out
#' @return A vector of the means
#' 
#' @examples
#' test_df <- data.frame(matrix(c(1,2,3,0,5,2,3,4,0,6,0,0,0,0,0,4,5,6,0,7),ncol=4))
#' ## 1. Average over non-zero rows:
#' average_over_present(test_df,1)
#' ## 2. Average over non-zero columns:
#' average_over_present(test_df,2)
#' 
#' @export
#' 
average_over_present <- function(in_df,in_dimension) {
  out_vector <- c()
  if(in_dimension==2) {
    row_ind <- which(rowSums(in_df)>0)
    reduced_df <- in_df[row_ind,]
    if (length(row_ind)==0) {
      return(NULL)
    } else if (length(row_ind)==1) {
      out_vector <- reduced_df
    } else {
      out_vector <- apply(reduced_df,2,mean)      
    }
  } else if (in_dimension==1) {
    col_ind <- which(colSums(in_df)>0)
    reduced_df <- in_df[,col_ind]
    if (length(col_ind)==0) {
      return(NULL)
    } else if (length(col_ind)==1) {
      out_vector <- reduced_df
    } else {
      out_vector <- apply(reduced_df,1,mean)      
    }
  } else {
    return(NULL)
  }
  return(out_vector)
}


#' Standard deviation of a data frame over a specified dimension
#'
#' If computing the standard deviation over columns, zero rows
#' (i.e. those with \code{rowSums=0}) are left out, if computing
#' the standard deviation over rows, zero columns (i.e. those with
#' \code{colSums=0}) are left out.
#'
#' @param in_df
#'  Data frame the standard deviation is computed of
#' @param in_dimension
#'  Dimension along which the computation will be carried out
#' @return A vector of the standard deviations
#' 
#' @examples
#' test_df <- data.frame(matrix(c(1,2,3,0,5,2,3,4,0,6,0,0,0,0,0,4,5,6,0,7),ncol=4))
#' ## 1. Compute standard deviation over non-zero rows:
#' sd_over_present(test_df,1)
#' ## 2. Compute standard deviation over non-zero columns:
#' sd_over_present(test_df,2)
#' 
#' @export
#' 
sd_over_present <- function(in_df,in_dimension) {
  out_vector <- c()
  if(in_dimension==2) {
    row_ind <- which(rowSums(in_df)>0)
    reduced_df <- in_df[row_ind,]
    if (length(row_ind)==0) {
      return(NULL)
    } else if (length(row_ind)==1) {
      out_vector <- rep(0,dim(in_df)[1])
    } else {
      out_vector <- apply(reduced_df,2,sd)      
    }
  } else if (in_dimension==1) {
    col_ind <- which(colSums(in_df)>0)
    reduced_df <- in_df[,col_ind]
    if (length(col_ind)==0) {
      return(NULL)
    } else if (length(col_ind)==1) {
      out_vector <- rep(0,dim(in_df)[2])
    } else {
      out_vector <- apply(reduced_df,1,sd)      
    }
  } else {
    return(NULL)
  }
  return(out_vector)
}


#' Elementwise sum over a list of (numerical) data frames
#'
#' @param in_df_list
#'  List of (numerical) data frames
#' @return A numerical data frame with the same dimensions as the entries of
#'  \code{in_df_list} with elementwise sums
#' 
#' @examples
#' A <- data.frame(matrix(c(1,1,1,2,2,2),ncol=2))
#' B <- data.frame(matrix(c(3,3,3,4,4,4),ncol=2))
#' df_list <- list(A=A,B=B)
#' sum_over_list_of_df(df_list)
#' 
#' @export
#' 
sum_over_list_of_df <- function(in_df_list) {
  sum_df <- repeat_df(0,dim(in_df_list[[1]])[1],dim(in_df_list[[1]])[2])
  for (i in seq_len(length(in_df_list))) {
    sum_df <- sum_df + in_df_list[[i]]
  }
  colnames(sum_df) <- colnames(in_df_list[[1]])
  rownames(sum_df) <- rownames(in_df_list[[1]])
  return(sum_df)
}


#' Compute the standard error of the mean
#'
#' This function returns the standard deviation of an input numerical vector 
#' divided by the square root of the length of the input vector
#' 
#' @param x
#'  A numerical vector
#'  
#' @return
#'  Standard deviation of an input numerical vector divided by the square
#'  root of the length of the input vector
#'
#' @examples
#' A <- c(1,2,3)
#' sd(A)
#' stderrmean(A)
#'  
#' @export
#'  
stderrmean <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))


#' Standard error of the mean of a data frame over a specified dimension
#'
#' If computing the standard error of the mean over columns, zero rows
#' (i.e. those with \code{rowSums=0}) are left out, if computing the
#' standard error of the mean over rows, zero columns (i.e. those with
#' \code{colSums=0}) are left out. Uses the function \code{\link{stderrmean}}
#'
#' @param in_df
#'  Data frame the standard error of the mean is computed of
#' @param in_dimension
#'  Dimension along which the computation will be carried out
#' @return A vector of the standard errors of the mean
#' 
#' @seealso \code{\link{stderrmean}}
#' 
#' @examples
#' test_df <- data.frame(matrix(c(1,2,3,0,5,2,3,4,0,6,0,0,0,0,0,4,5,6,0,7),ncol=4))
#' ## 1. Compute standard deviation over non-zero rows:
#' stderrmean_over_present(test_df,1)
#' ## 2. Compute standard deviation over non-zero columns:
#' stderrmean_over_present(test_df,2)
#' 
#' @export
#' 
stderrmean_over_present <- function(in_df,in_dimension) {
  out_vector <- c()
  if(in_dimension==2) {
    row_ind <- which(rowSums(in_df)>0)
    reduced_df <- in_df[row_ind,]
    if (length(row_ind)==0) {
      return(NULL)
    } else if (length(row_ind)==1) {
      out_vector <- rep(0,dim(in_df)[1])
    } else {
      out_vector <- apply(reduced_df,2,stderrmean)      
    }
  } else if (in_dimension==1) {
    col_ind <- which(colSums(in_df)>0)
    reduced_df <- in_df[,col_ind]
    if (length(col_ind)==0) {
      return(NULL)
    } else if (length(col_ind)==1) {
      out_vector <- rep(0,dim(in_df)[2])
    } else {
      out_vector <- apply(reduced_df,1,stderrmean)      
    }
  } else {
    return(NULL)
  }
  return(out_vector)
}


#' Translate chromosome names to the hg19 naming convention
#'
#' In hg19 naming convention, chromosome names start with the prefix \emph{chr} and
#' the gonosomes are called \emph{X} and \emph{Y}. If data analysis is performed e.g.
#' with \code{\link[BSgenome.Hsapiens.UCSC.hg19]{BSgenome.Hsapiens.UCSC.hg19}}, this
#' naming convention is needed. The inverse transform is done with
#' \code{\link{translate_to_1kG}}.
#'
#' @param in_df
#'  Data frame which carries one column with chromosome information to be reformatted
#' @param in_CHROM.field
#'  String indicating which column of \code{in_df} carries the chromosome information
#'  
#' @return A data frame identical to \code{in_df}, but with the names in the
#'  chromosome column replaced
#' 
#' @seealso \code{\link{translate_to_1kG}}
#' 
#' @examples
#' test_df <- data.frame(CHROM=c(1,2,23,24),POS=c(100,120000000,300000,25000),dummy=c("a","b","c","d"))
#' hg19_df <- translate_to_hg19(test_df, in_CHROM.field = "CHROM")
#' hg19_df
#' 
#' @export
#' 
translate_to_hg19 <- function(in_df,in_CHROM.field="CHROM") {
  out_df <- in_df
  chrom_ind <- which(names(out_df)==in_CHROM.field)
  #names(out_df)[chrom_ind] <- "chr"
  out_df[,chrom_ind] <- as.character(out_df[,chrom_ind])
  out_df[,chrom_ind] <- gsub("23","X",out_df[,chrom_ind])
  out_df[,chrom_ind] <- gsub("24","Y",out_df[,chrom_ind])
  logical_chr_vector <- grepl("^chr.*",out_df[,chrom_ind])
  if(!(all(logical_chr_vector))) {
    replace_ind <- which(!logical_chr_vector)
    out_df[replace_ind,chrom_ind] <- paste0("chr",out_df[replace_ind,chrom_ind])
  }
  return(out_df)
}


#' Translate chromosome names to the 1kG naming convention
#'
#' In 1kG, i.e. 1000 genomes naming convention, chromosome names have no prefix
#' \emph{chr} and the gonosomes are called \emph{23} for \emph{X} and \emph{24}
#' for \emph{Y}. If data analysis is performed e.g. with \code{hs37d5.fa}, this
#' naming convention is needed. The inverse transform is done with
#' \code{\link{translate_to_hg19}}.
#'
#' @param in_df
#'  Data frame which carries one column with chromosome information to be reformatted
#' @param in_CHROM.field
#'  String indicating which column of \code{in_df} carries the chromosome information
#'  
#' @return A data frame identical to \code{in_df}, but with the names in the
#'  chromosome column replaced
#' 
#' @seealso \code{\link{translate_to_1kG}}
#' 
#' @examples
#' test_df <- data.frame(CHROM=c(1,2,23,24),POS=c(100,120000000,300000,25000),dummy=c("a","b","c","d"))
#' hg19_df <- translate_to_hg19(test_df, in_CHROM.field = "CHROM")
#' onekG_df <- translate_to_1kG(hg19_df, in_CHROM.field = "CHROM")
#' onekG_df
#' 
#' @export
#' 
translate_to_1kG <- function(in_df,in_CHROM.field="chr") {
  out_df <- in_df
  chrom_ind <- which(names(out_df)==in_CHROM.field)
  #names(out_df)[chrom_ind] <- "CHROM"
  out_df[,chrom_ind] <- as.character(out_df[,chrom_ind])
  out_df[,chrom_ind] <- gsub("X","23",out_df[,chrom_ind])
  out_df[,chrom_ind] <- gsub("Y","24",out_df[,chrom_ind])
  logical_chr_vector <- grepl("^[0-9]+$",out_df[,chrom_ind])
  if(!(all(logical_chr_vector))) {
    replace_ind <- which(!logical_chr_vector)
    out_df[replace_ind,chrom_ind] <- gsub("chr","",out_df[replace_ind,chrom_ind])
  }
  return(out_df)
}


#' Attribute the nucleotide exchange for an SNV
#'
#' SNVs are grouped into 6 different categories (12/2 as reverse complements are summed
#' over). This function defines the attribution.
#'
#' @param in_dat
#'  Data frame which carries one column for the reference base and one column for the
#'  variant base
#' @param in_REF.field
#'  String indicating which column of \code{in_dat} carries the reference base
#' @param in_ALT.field
#'  String indicating which column of \code{in_dat} carries the variant base
#'  
#' @return A character vector with as many rows as there are in \code{in_dat} which
#'  can be annotated (i.e. appended) to the input data frame.
#' 
#' @examples
#' test_df <- data.frame(CHROM=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5),
#'                        POS=c(1,2,3,4,5,6,1,2,3,4,5,6,7,8),
#'                        REF=c("C","C","C","T","T","T","A","A","A","G","G","G","N","A"),
#'                        ALT=c("A","G","T","A","C","G","C","G","T","A","C","T","A","N"))
#' test_df$change <- attribute_nucleotide_exchanges(test_df,in_REF.field = "REF",in_ALT.field = "ALT")
#' test_df
#' 
#' @export
#' 
attribute_nucleotide_exchanges <- function(in_dat,in_REF.field="REF",in_ALT.field="ALT") {
  name_list <- names(in_dat)
  ## exception handling for input fields
  if(tolower(in_REF.field) %in% tolower(name_list)) {
    cat("attribute_nucleotide_exchanges::in_REF.field found. Retrieving REF information.\n")
    REF_ind <- min(which(tolower(name_list)==tolower(in_REF.field)))
  } else {
    cat("attribute_nucleotide_exchanges::error: in_REF.field not found. Return NULL.\n")
    return(NULL)
  }
  if(tolower(in_ALT.field) %in% tolower(name_list)) {
    cat("attribute_nucleotide_exchanges::in_ALT.field found. Retrieving ALT information.\n")
    ALT_ind <- min(which(tolower(name_list)==tolower(in_ALT.field)))
  } else {
    cat("attribute_nucleotide_exchanges::error: in_ALT.field not found. Return NULL.\n")
    return(NULL)
  }
  ## adapt nomenclature for nucleotide exchanges
  complement <- c('A' = 'T', 'C' = 'G','G' = 'C', 'T' = 'A', 'N' = 'N')
  change_vec <- rep("dummy",dim(in_dat)[1])
  sel <- which(in_dat[,REF_ind] == 'C' | in_dat[,REF_ind] == 'T')
  change_vec[sel] <- paste0(as.character(in_dat[sel,REF_ind]),as.character(in_dat[sel,ALT_ind]))
  sel <- which(in_dat[,REF_ind] == 'A' | in_dat[,REF_ind] == 'G')
  change_vec[sel] <- paste0(complement[as.character(in_dat[sel,REF_ind])],complement[as.character(in_dat[sel,ALT_ind])])  
  change_vec <- factor(change_vec, levels = c("CA", "CG", "CT", "TA", "TC", "TG"))
  return(change_vec)
}


#' Annotate the intermutation distance of variants per PID
#'
#' The function annotates the intermutational distance to a PID wide data frame by
#' applying \code{\link[circlize]{rainfallTransform}} to every chromosome-specific
#' subfraction of the PID wide data.
#'
#' @param in_df
#'  Data frame which carries (at least) one column for the chromosome and one column
#'  for the position.
#' @param in_CHROM.field
#'  String indicating which column of \code{in_df} carries the chromosome information
#' @param in_POS.field
#'  String indicating which column of \code{in_df} carries the position information
#' @param in_mode
#'  String passed to \code{\link[circlize]{rainfallTransform}} indicating which
#'  method to choose for the computation of the intermutational distance.
#'  
#' @return A data frame identical to \code{in_df}, but with the intermutation distance
#' annotated as an additional column on the right named \code{dist}.
#'  
#' @seealso \code{\link{annotate_intermut_dist_cohort}}
#' @seealso \code{\link[circlize]{rainfallTransform}}
#' 
#' @examples
#' test_df <- data.frame(CHROM=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5),
#'                        POS=c(1,2,4,4,6,9,1,4,8,10,20,40,100,200),
#'                        REF=c("C","C","C","T","T","T","A","A","A","G","G","G","N","A"),
#'                        ALT=c("A","G","T","A","C","G","C","G","T","A","C","T","A","N"))
#' min_dist_df <- annotate_intermut_dist_PID(test_df,in_CHROM.field="CHROM",in_POS.field="POS",
#'                                            in_mode="min")
#' max_dist_df <- annotate_intermut_dist_PID(test_df,in_CHROM.field="CHROM",in_POS.field="POS",
#'                                            in_mode="max")
#' min_dist_df
#' max_dist_df
#' 
#' @importFrom circlize rainfallTransform
#' @export
#' 
annotate_intermut_dist_PID <- function(in_df,in_CHROM.field="CHROM",in_POS.field="POS",in_mode="min") {
  CHROM_ind <- which(names(in_df)==in_CHROM.field)
  POS_ind <- which(names(in_df)==in_POS.field)
  inflate_df <- data.frame()
  for(my_chrom in unique(in_df[,CHROM_ind])) {
    my_df <- in_df[which(in_df[,CHROM_ind]==my_chrom),]
    region_df <- data.frame(start=my_df[,POS_ind],end=my_df[,POS_ind])
    if(dim(my_df)[1]<2) {
      temp_df <- region_df
      temp_df$dist <- c(100000000)
    } else {
      temp_df <- rainfallTransform(region_df, mode = in_mode)      
    }
    temp_df$CHROM <- my_chrom
    inflate_df <- rbind(inflate_df,temp_df)
  }
#   df_list <- lapply(unique(in_df[,CHROM_ind]),function(l) {
#     my_df <- in_df[which(in_df[,CHROM_ind]==l),]
#     region_df <- data.frame(start=my_df[,POS_ind],end=my_df[,POS_ind])
#     if(dim(my_df)[1]<2) {
#       temp_df <- region_df
#       temp_df$dist <- c(100000000)
#     } else {
#       temp_df <- rainfallTransform(region_df, mode = in_mode)      
#     }
#     temp_df$CHROM <- l
#     return(temp_df)
#   })
#   inflate_df <- do.call("rbind",df_list)
  out_df <- in_df
  out_df$dist <- inflate_df[,dim(inflate_df)[2]-1]
  return(out_df)
}


#' Annotate the intermutation distance of variants cohort-wide
#'
#' The function annotates the intermutational distance to a cohort wide data frame by
#' applying \code{\link{annotate_intermut_dist_PID}} to every PID-specific subfraction
#' of the cohort wide data. Note that \code{\link{annotate_intermut_dist_PID}} calls
#' \code{\link[circlize]{rainfallTransform}}. If the PID information is missing,
#' \code{\link{annotate_intermut_dist_PID}} is called directly for the whole input.
#'
#' @param in_df
#'  Data frame which carries (at least) one column for the chromosome and one column
#'  for the position. Optionally, a column to specify the PID can be provided.
#' @param in_CHROM.field
#'  String indicating which column of \code{in_df} carries the chromosome information
#' @param in_POS.field
#'  String indicating which column of \code{in_df} carries the position information
#' @param in_PID.field
#'  String indicating which column of \code{in_df} carries the PID information
#' @param in_mode
#'  String passed through \code{\link{annotate_intermut_dist_PID}} to
#'  \code{\link[circlize]{rainfallTransform}} indicating which method to choose for
#'  the computation of the intermutational distance.
#'  
#' @return A data frame identical to \code{in_df}, but with the intermutation distance
#' annotated as an additional column on the right named \code{dist}.
#'  
#' @seealso \code{\link{annotate_intermut_dist_PID}}
#' @seealso \code{\link[circlize]{rainfallTransform}}
#' 
#' @examples
#' test_df <- data.frame(CHROM=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5),
#'                        POS=c(1,2,4,4,6,9,1,4,8,10,20,40,100,200),
#'                        REF=c("C","C","C","T","T","T","A","A","A","G","G","G","N","A"),
#'                        ALT=c("A","G","T","A","C","G","C","G","T","A","C","T","A","N"),
#'                        PID=c(1,1,1,2,2,2,1,1,2,2,2,1,1,2))
#' test_df <- test_df[order(test_df$PID,test_df$CHROM,test_df$POS),]
#' min_dist_df <- annotate_intermut_dist_cohort(test_df,in_CHROM.field="CHROM",in_POS.field="POS",
#'                                              in_PID.field="PID",in_mode="min")
#' max_dist_df <- annotate_intermut_dist_cohort(test_df,in_CHROM.field="CHROM",in_POS.field="POS",
#'                                              in_PID.field="PID",in_mode="max")
#' min_dist_df
#' max_dist_df
#' 
#' @export
#' 
annotate_intermut_dist_cohort <- function(in_df,in_CHROM.field="CHROM",in_POS.field="POS",
                                          in_PID.field=NULL,in_mode="min") {
  if(!is.null(in_PID.field)) {
    PID_ind <- which(names(in_df)==in_PID.field)
#     inflate_df <- data.frame()
#     for(my_PID in unique(in_df[,PID_ind])) {
#       temp_df <- annotate_intermut_dist_PID(in_df[which(in_df[,PID_ind]==my_PID),],
#                                             in_CHROM.field=in_CHROM.field,
#                                             in_POS.field=in_POS.field,in_mode=in_mode)
#       temp_df$PID <- my_PID
#       inflate_df <- rbind(inflate_df,temp_df)
#     }
    df_list <- lapply(unique(in_df[,PID_ind]),function(l) {
      temp_df <- annotate_intermut_dist_PID(in_df[which(in_df[,PID_ind]==l),],
                                            in_CHROM.field=in_CHROM.field,
                                            in_POS.field=in_POS.field,in_mode=in_mode)
      temp_df$PID <- l
      return(temp_df)
    })
    inflate_df <- do.call("rbind",df_list)
    out_df <- in_df
    out_df$dist <- inflate_df$dist
  } else {
    out_df <- annotate_intermut_dist_PID(in_df,in_CHROM.field=in_CHROM.field,
                                         in_POS.field=in_POS.field,in_mode=in_mode)
  }
  return(out_df)
}


#' Make a custom data structure for subgroups
#'
#' Creates a data frame carrying the subgroup information and the order in which the PIDs
#' have to be displayed. Calls \code{\link{aggregate}} on \code{in_vcf_like_df}.
#'
#' @param in_exposures_df
#'  Data frame with the signature exposures
#' @param in_vcf_like_df
#'  vcf-like data frame with point mutation calls
#' @param in_palette
#'  Palette for colour attribution to the subgroups if nun-NULL
#' @param in_subgroup.field
#'  String indicating which column of \code{in_vcf_like_df} carries the subgroup information
#' @param in_PID.field
#'  String indicating which column of \code{in_vcf_like_df} and of \code{in_exposures_df}
#'  carries the PID information
#'  
#' @return subgroups_df:
#'  A data frame carrying the subgroup and rank information.
#' 
#' @examples
#'  data(lymphoma_test)
#'  data(lymphoma_cohort_LCD_results)
#'  choice_ind <- (names(lymphoma_Nature2013_COSMIC_cutoff_exposures_df) 
#'                 %in% unique(lymphoma_test_df$PID))
#'  lymphoma_test_exposures_df <- lymphoma_Nature2013_COSMIC_cutoff_exposures_df[,choice_ind]
#'  make_subgroups_df(lymphoma_test_exposures_df,lymphoma_test_df)
#' 
#' @seealso \code{\link{aggregate}}
#' 
#' @export
#' 
make_subgroups_df <- function(in_exposures_df,
                              in_vcf_like_df,
                              in_palette=NULL,
                              in_subgroup.field="SUBGROUP",
                              in_PID.field="PID"){
  ## 1. rename the colnames in in_vcf_like_df if necessary to be able to run the
  ## aggregate command later
  subgroup_ind <- which(names(in_vcf_like_df)==in_subgroup.field)
  names(in_vcf_like_df)[subgroup_ind] <- "SUBGROUP"
  PID_ind <- which(names(in_vcf_like_df)==in_PID.field)
  names(in_vcf_like_df)[PID_ind] <- "PID"
  ## 2. start work
  this_sum_df <- data.frame(sum=apply(in_exposures_df,2,sum))
  out_subgroups_df <- aggregate(SUBGROUP~PID,data=in_vcf_like_df,
                                   function(l) return(l[1]))
  out_subgroups_df <- merge(out_subgroups_df,this_sum_df,by.x=in_PID.field,by.y=0)
  max_total_count <- max(out_subgroups_df$sum)
  out_subgroups_df$compl_sum <- max_total_count - out_subgroups_df$sum
  temp_ind <- order(out_subgroups_df$SUBGROUP,out_subgroups_df$compl_sum)
  out_subgroups_df$index <- order(temp_ind)
  names(out_subgroups_df)[2] <- "subgroup"
  if(is.null(in_palette)){
    in_palette <- rainbow(length(unique(out_subgroups_df$subgroup)))
  }
  out_subgroups_df$col <- in_palette[factor(out_subgroups_df$subgroup)]
  out_subgroups_df$subgroup <- as.character(out_subgroups_df$subgroup)
  return(out_subgroups_df)
}


#' Wrapper for Shapiro test but allow for all identical values
#' 
#' @param in_vector
#'  Numerical vector the Shapiro-Wilk test is computed on
#' 
#' @return p-value of the Shapiro-Wilk test, zero if all entries in the input
#'  vector \code{in_vector} are identical.
#' 
#' @seealso \code{\link[stats]{shapiro.test}}
#' 
#' @examples
#'  shapiro_if_possible(runif(100,min=2,max=4))
#'  shapiro_if_possible(rnorm(100,mean=5,sd=3))
#'  shapiro_if_possible(rep(4.3,100))
#'  shapiro_if_possible(c("Hello","World"))
#' 
#' @export
#' 
shapiro_if_possible <- function(in_vector){
  if(is.numeric(in_vector)){
    if(length(unique(in_vector)) > 1){
      return(shapiro.test(in_vector)$p.value)
    } else {
      return(0.0)
    }
  } else {
    cat("YAPSA:::shapiro_if_possible::error: Non-numeric input")
    return(NULL)
  }
}


#' Split an exposures data frame by subgroups
#' 
#' If a cohort consists of different subgroups, this function enables
#' to split the data frame storing the signature exposures into a list
#' of data frames with signature exposures, one per subgroup. This
#' functionality is needed for \code{\link{stat_test_subgroups}} and
#' \code{\link{stat_plot_subgroups}}
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
#' @return List of data frames with the subgroup specific signature
#'  exposures.
#' 
#' @seealso \code{\link{stat_test_subgroups}}
#' @seealso \code{\link{stat_plot_subgroups}}
#' 
#' @examples
#'  NULL
#' 
#' @export
#' 
split_exposures_by_subgroups <- function(in_exposures_df,in_subgroups_df,
                                         in_subgroups.field="subgroup",
                                         in_PID.field="PID"){
  subgroup_index <- which(names(in_subgroups_df)==in_subgroups.field)
  PID_index <- which(names(in_subgroups_df)==in_PID.field)
  # split in_exposures_df
  subgroups_list <- split(in_subgroups_df,in_subgroups_df[,subgroup_index])
  exposures_list <- lapply(subgroups_list, FUN=function(x) {
    PID_vector <- x[,PID_index]
    choice_ind <- which(names(in_exposures_df) %in% PID_vector)
    if(length(choice_ind)>1){
      out_df <- in_exposures_df[,choice_ind]
    } else {
      out_df <- data.frame(temp=in_exposures_df[,choice_ind])
      names(out_df)[1] <- PID_vector[1]
    }
    return(out_df)
  })
  return(exposures_list)
}



#' Add an element as first entry to a list
#' 
#' Works for all types of lists and inputs
#' 
#' @param in_list
#'  List to which an element is to be added
#' @param in_element
#'  Element to be added
#' 
#' @return List with input element as first entry.
#' 
#' @examples
#'  NULL
#' 
#' @export
#' 
add_as_fist_to_list <- function(in_list,in_element){
  out_list <- in_list
  number_of_elements <- length(in_list)
  out_list[[number_of_elements+1]] <- in_element
  my_ind <- c(number_of_elements+1,seq_len(number_of_elements))
  return(out_list[my_ind])
}


#' Merge exposure data frames
#' 
#' Merges with the special feature of preserving the signatures and signature 
#' order.
#' 
#' @param in_exposures_list
#'  List of data frames (carrying information on exposures).
#' @param in_signatures_df
#'  Data frame \code{W} in which the columns represent the signatures.
#' 
#' @return A data frame with the merged exposures.
#' 
#' @examples
#'  NULL
#' 
#' @export
#' 
merge_exposures <- function(in_exposures_list,
                            in_signatures_df){
  exposures_list <- lapply(in_exposures_list,FUN=function(current_df) {
    current_df$sig <- rownames(current_df)
    return(current_df)
  })
  exposures_df <- Reduce(function(...) merge(..., by="sig",all=TRUE,sort=FALSE),
                         exposures_list)
  row_order <- order(match(exposures_df$sig,names(in_signatures_df)))
  exposures_df <- exposures_df[row_order,]
  rownames(exposures_df) <- exposures_df$sig
  exposures_df$sig <- NULL
  exposures_df[is.na(exposures_df)] <- 0
  return(exposures_df)
}


#' Generically melts exposure data frames
#' 
#' Melt an exposure data frame with signatures as ID variables.
#' 
#' @param in_df
#'  Numeric data frame with exposures.
#' 
#' @return A data frame with the molten exposures.
#' 
#' @examples
#'  NULL
#' 
#' @export
#' 
melt_exposures <- function(in_df){
  in_df$sig <- rownames(in_df)
  in_df_melt <- melt(in_df,id.vars = "sig")
  names(in_df_melt) <- gsub("variable","PID",names(in_df_melt))
  return(in_df_melt)
}


#' Compares alternative exposures
#' 
#' Compares exposures computed by two alternative approaches for the same cohort
#' 
#' @param in_exposures1_df
#'  Numeric data frame with exposures, ideally the smaller exposure data is
#'  supplied first.
#' @param in_exposures2_df
#'  Numeric data frame with exposures, ideally the bigger exposure data is
#'  supplied second.
#' @param deselect_flag
#'  Wehther signatures absent in both exposure data frames should be removed.
#' 
#' @return A list with entries \code{merge_df}, \code{all_cor.coeff},
#'  \code{all_p.value}, \code{cor.coeff_vector}, \code{p.value_vector},
#'  \code{all_cor.test}, and \code{cor.test_list}.
#' \itemize{
#'  \item \code{merge_df}:
#'    Merged molten input exposure data frames
#'  \item \code{all_cor.coeff}:
#'    Pearson correlation coefficient for all data 
#'    points, i.e. taken all signatures together
#'  \item \code{all_p.value}:
#'    P-value of the Pearson test for all data 
#'    points, i.e. taken all signatures together
#'  \item \code{cor.coeff_vector}:
#'    A vector of Pearson correlation coefficients 
#'    evaluated for every signature independently
#'  \item \code{p.value_vector}:
#'    A vector of p-values of the Pearson tests 
#'    evaluated for every signature independently
#'  \item \code{all_cor.test}:
#'    A data structure as returned by \code{\link{cor.test}} for all data 
#'    points, i.e. taken all signatures together
#'  \item \code{cor.test_list}:
#'    A list of data structures as returned by \code{\link{cor.test}}, but 
#'    evaluated for every signature independently
#' }
#' 
#' @examples
#'  NULL
#' 
#' @export
#' 
compare_exposures <- function(in_exposures1_df,
                              in_exposures2_df,
                              deselect_flag=TRUE){
  signatures1 <- rownames(in_exposures1_df)
  signatures2 <- rownames(in_exposures2_df)
  common_signatures <- intersect(signatures1,signatures2)
  if(length(common_signatures) == 0) return(NULL)
  common_signatures <- union(signatures2,signatures1)
  exposures1_df_melt <- melt_exposures(in_exposures1_df)
  names(exposures1_df_melt) <- gsub("value","exposures1",
                                    names(exposures1_df_melt))
  exposures2_df_melt <- melt_exposures(in_exposures2_df)
  names(exposures2_df_melt) <- gsub("value","exposures2",
                                    names(exposures2_df_melt))
  merge_df <- merge(exposures1_df_melt,exposures2_df_melt,by=c("sig","PID"),
                    sort=FALSE,all=TRUE)
  merge_df$exposures1[which(!is.finite(merge_df$exposures1))] <- 0
  merge_df$exposures2[which(!is.finite(merge_df$exposures2))] <- 0
  if(deselect_flag){
    deselect_ind <- which(merge_df$exposures1==0 & merge_df$exposures2==0)
    merge_df <- merge_df[-deselect_ind,]    
  }
  all_cor.test <- cor.test(merge_df$exposures1,merge_df$exposures2)
  cor.test_list <- lapply(
    split(merge_df,f = merge_df$sig),
    FUN=function(current_df){
      levels1 <- length(unique(current_df$exposures1))
      levels2 <- length(unique(current_df$exposures2))
      if(levels1 > 1 & levels2 > 1) {
        out.test <- cor.test(current_df$exposures1,current_df$exposures2)
      } else {
        out.test <- NULL       
      } 
    })
  cor.coeff_vector <- unlist(lapply(
    cor.test_list,
    function(current_test){
      if(is.null(current_test)) return(0)
      return(as.numeric(current_test$estimate))
    }))
  p.value_vector <- unlist(lapply(
    cor.test_list,
    function(current_test){
      if(is.null(current_test)) return(1)
      return(as.numeric(current_test$p.value))
    }))
  match_ind <- match(common_signatures,names(cor.test_list))
  cor.test_list <- cor.test_list[match_ind]
  cor.coeff_vector <- cor.coeff_vector[match_ind]
  p.value_vector <- p.value_vector[match_ind]
  return(list(merge_df=merge_df,
              all_cor.coeff=all_cor.test$estimate,
              all_p.value=all_cor.test$p.value,
              cor.coeff_vector=cor.coeff_vector,
              p.value_vector=p.value_vector,
              all_cor.test=all_cor.test,
              cor.test_list=cor.test_list))
}


#' Information on initially published mutational signatures
#' 
#' Data frame of the signatures published initially by Alexandrov et al. (Nature 
#' 2013). There are 27 signatures which constitute the columns, 22 of which were
#' validated by an orhtogonal sequencing technology. These 22 are in the first 22
#' columns of the data frame. The column names are \emph{A} pasted to the number
#' of the signature, e.g. \emph{A5}. The nonvalidated signatures have an
#' additional letter in their naming convention: either \emph{AR1} - \emph{AR3}
#' or \emph{AU1} - \emph{AU2}. The rownames are the features, i.e. an encoding of
#' the nucleotide exchanges in their trinucleotide context, e.g. \emph{C>A ACA}.
#' In total there are 96 different features and therefore 96 rows when dealing
#' with a trinucleotide context. This data was downloaded from
#' \url{ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/signatures.txt} and
#' reformatted.
#' 
#' @seealso \code{\link{AlexInitialArtif_sigInd_df}}
#' @seealso \code{\link{AlexCosmicValid_sig_df}}
#' @seealso \code{\link{AlexInitialValid_sig_df}}
#' 
#' @docType data
#' @name AlexInitialArtif_sig_df
#' @usage data(sigs)
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' @references Alexandrov et al. (Nature 2013)
#' 
NULL


#' Meta-information on initially published mutational signatures
#' 
#' A data frame with as many rows as there are signatures in
#' \code{AlexInitialValid_sig_df} and several columns: 
#' \itemize{
#'  \item \code{sig}: signature name
#'  \item \code{index}: corresponding to the row index of the signature
#'  \item \code{colour}: colour for visualization in stacked barplots
#'  \item \code{process}: asserted biological process
#'  \item \code{cat.coarse}: categorization of the signatures according
#'   to the asserted biological processes at low level of detail
#'  \item \code{cat.medium}: categorization of the signatures according
#'   to the asserted biological processes at intermediate level of detail
#'  \item \code{cat.high}: categorization of the signatures according
#'   to the asserted biological processes at high level of detail
#'  \item \code{cat.putative}: categorization of the signatures according
#'   to the asserted biological processes based on clustering and inference
#'  }
#' 
#' @seealso \code{\link{AlexInitialArtif_sig_df}}
#' @seealso \code{\link{AlexInitialValid_sig_df}}
#' 
#' @docType data
#' @name AlexInitialArtif_sigInd_df
#' @usage data(sigs)
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' 
NULL


#' Sub-data frame of the initially published mutational signatures
#' 
#' Data frame of only the validated signatures published initially by Alexandrov
#' et al. (Nature 2013), corresponding to the first 22 columns of
#' \code{AlexInitialArtif_sig_df}
#' 
#' @seealso \code{\link{AlexInitialArtif_sig_df}}
#' 
#' @docType data
#' @name AlexInitialValid_sig_df
#' @usage data(sigs)
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' @references Alexandrov et al. (Nature 2013)
#' 
NULL


#' Information on updated mutational signatures
#' 
#' Data frame of the updated signatures list maintained by Ludmil Alexandrov at
#' \url{http://cancer.sanger.ac.uk/cosmic/signatures}. The data was downloaded from
#' \url{http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt}
#' and reformatted. The column names are \emph{AC} pasted to the number of the
#' signature, e.g. \emph{AC5}. The naming convention for the rows is as described
#' for \code{\link{AlexInitialArtif_sig_df}}.
#' 
#' @seealso \code{\link{AlexInitialArtif_sig_df}}
#' 
#' @docType data
#' @name AlexCosmicValid_sig_df
#' @usage data(sigs)
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' @references Alexandrov et al. (Nature 2013)
#' 
NULL


#' Meta-information on updated mutational signatures
#' 
#' A data frame with as many rows as there are signatures in
#' \code{AlexCosmicValid_sig_df} and several columns: 
#' \itemize{
#'  \item \code{sig}: signature name
#'  \item \code{index}: corresponding to the row index of the signature
#'  \item \code{colour}: colour for visualization in stacked barplots
#'  \item \code{process}: asserted biological process
#'  \item \code{cat.coarse}: categorization of the signatures according
#'   to the asserted biological processes at low level of detail
#'  \item \code{cat.medium}: categorization of the signatures according
#'   to the asserted biological processes at intermediate level of detail
#'  \item \code{cat.fine}: categorization of the signatures according
#'   to the asserted biological processes at high level of detail
#'  \item \code{cat.putative}: categorization of the signatures according
#'   to the asserted biological processes based on clustering and inference
#'  }
#' 
#' @seealso \code{\link{AlexCosmicValid_sig_df}}
#' 
#' @docType data
#' @name AlexCosmicValid_sigInd_df
#' @usage data(sigs)
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' 
NULL


#' Subgroup information for some samples in the vignette
#' 
#' A data frame carrying subgroup information for a subcohort of samples used in the vignette.
#' Data in the vignette is downloaded from \url{ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Lymphoma B-cell/Lymphoma B-cell_clean_somatic_mutations_for_signature_analysis.txt}.
#' In the file available under that link somatic point mutation calls from several samples are 
#' listed in a vcf-like format. One column encodes the sample the variant was found in. In the
#' vignette we want to restrict the analysis to only a fraction of these involved samples. The
#' data frame \code{lymphoma_PID_df} carries the sample identifiers (PID) as rownames and the
#' attributed subgroup in a column called \code{subgroup}.
#' 
#' @docType data
#' @name lymphoma_PID_df
#' @usage data(lymphoma_PID)
#' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/23945592}
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' 
NULL


#' Test data for complex functions
#' 
#' A data frame carrying point mutation calls. It represents a subset of the data stored in
#' \url{ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Lymphoma B-cell/Lymphoma B-cell_clean_somatic_mutations_for_signature_analysis.txt}.
#' In the file available under that link somatic point mutation calls from several samples are 
#' listed in a vcf-like format. One column encodes the sample the variant was found in. The
#' data frame \code{lymphoma_test_df} has only the variants occuring in the sample identifiers
#' (PIDs) 4112512, 4194218 and 4121361.
#' 
#' @docType data
#' @name lymphoma_test_df
#' @usage data(lymphoma_test)
#' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/23945592}
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' 
#' @examples
#' data(lymphoma_test)
#' head(lymphoma_test_df)
#' dim(lymphoma_test_df)
#' table(lymphoma_test_df$PID)
#' 
NULL


#' Test exposures for plot functions
#' 
#' Data frame with exposures for testing the plot functions. Data taken from
#' \url{ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Lymphoma B-cell/Lymphoma B-cell_clean_somatic_mutations_for_signature_analysis.txt}.
#' 
#' @docType data
#' @name lymphoma_Nature2013_COSMIC_cutoff_exposures_df
#' @usage data(lymphoma_cohort_LCD_results)
#' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/23945592}
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' 
#' @seealso \code{\link{rel_lymphoma_Nature2013_COSMIC_cutoff_exposures_df}}
#' @seealso \code{\link{COSMIC_subgroups_df}}
#' @seealso \code{\link{chosen_AlexInitialArtif_sigInd_df}}
#' 
NULL


#' Test normalized exposures for plot functions
#' 
#' Data frame with normalized or relative exposures for testing the plot functions.
#' Data taken from
#' \url{ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Lymphoma B-cell/Lymphoma B-cell_clean_somatic_mutations_for_signature_analysis.txt}.
#' 
#' @docType data
#' @name rel_lymphoma_Nature2013_COSMIC_cutoff_exposures_df
#' @usage data(lymphoma_cohort_LCD_results)
#' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/23945592}
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' 
#' @seealso \code{\link{lymphoma_Nature2013_COSMIC_cutoff_exposures_df}}
#' @seealso \code{\link{COSMIC_subgroups_df}}
#' @seealso \code{\link{chosen_AlexInitialArtif_sigInd_df}}
#' 
NULL


#' Subgroup information for test data for plot functions
#' 
#' Subgroup information for the data stored in
#' \code{\link{lymphoma_Nature2013_COSMIC_cutoff_exposures_df}} and 
#' \code{\link{rel_lymphoma_Nature2013_COSMIC_cutoff_exposures_df}}.
#' 
#' @docType data
#' @name COSMIC_subgroups_df
#' @usage data(lymphoma_cohort_LCD_results)
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' 
#' @seealso \code{\link{lymphoma_Nature2013_COSMIC_cutoff_exposures_df}}
#' @seealso \code{\link{rel_lymphoma_Nature2013_COSMIC_cutoff_exposures_df}}
#' 
NULL


#' Signature information for test data for plot functions
#' 
#' Signature information for the data stored in
#' \code{\link{lymphoma_Nature2013_COSMIC_cutoff_exposures_df}} and 
#' \code{\link{rel_lymphoma_Nature2013_COSMIC_cutoff_exposures_df}}.
#' 
#' @docType data
#' @name chosen_AlexInitialArtif_sigInd_df
#' @usage data(lymphoma_cohort_LCD_results)
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' 
#' @seealso \code{\link{lymphoma_Nature2013_COSMIC_cutoff_exposures_df}}
#' @seealso \code{\link{rel_lymphoma_Nature2013_COSMIC_cutoff_exposures_df}}
#' 
NULL


#' Optimal cutoffs for relative exposures for the COSMIC signatures
#' 
#' Data frame of signature specific cutoff values for the LCD analysis with
#' relative exposures. All values represent optimal cutoffs. The optimal
#' cutoffs were determined for different choices of parameters in the cost
#' function of the optimization. The row index is equivalent to the ratio
#' between costs for false negative attribution and false positive
#' attribution. The columns correspond to the different signatures. To be
#' used with \code{\link{LCD_complex_cutoff}}. Based on the signatures available
#' at http://cancer.sanger.ac.uk/cosmic/signatures .
#' 
#' @seealso \code{\link{cutoffCosmicValid_abs_df}}
#' @seealso \code{\link{cutoffInitialValid_rel_df}}
#' 
#' @docType data
#' @name cutoffCosmicValid_rel_df
#' @usage data(cutoffs)
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' 
NULL


#' Optimal cutoffs for absolute exposures for the COSMIC signatures
#' 
#' Data frame of signature specific cutoff values for the LCD analysis with
#' absolute exposures. All values represent optimal cutoffs. The optimal
#' cutoffs were determined for different choices of parameters in the cost
#' function of the optimization. The row index is equivalent to the ratio
#' between costs for false negative attribution and false positive
#' attribution. The columns correspond to the different signatures. To be
#' used with \code{\link{LCD_complex_cutoff}}. Based on the signatures available
#' at http://cancer.sanger.ac.uk/cosmic/signatures .
#' 
#' @seealso \code{\link{cutoffCosmicValid_rel_df}}
#' @seealso \code{\link{cutoffInitialValid_abs_df}}
#' 
#' @docType data
#' @name cutoffCosmicValid_abs_df
#' @usage data(cutoffs)
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' 
NULL


#' Optimal cutoffs for relative exposures of initial signatures
#' 
#' Data frame of signature specific cutoff values for the LCD analysis with
#' relative exposures. All values represent optimal cutoffs. The optimal
#' cutoffs were determined for different choices of parameters in the cost
#' function of the optimization. The row index is equivalent to the ratio
#' between costs for false negative attribution and false positive
#' attribution. The columns correspond to the different signatures. To be
#' used with \code{\link{LCD_complex_cutoff}}. Based on the signatures published
#' in Alexandrov et al. (Nature 2013).
#' 
#' @seealso \code{\link{cutoffCosmicValid_rel_df}}
#' @seealso \code{\link{cutoffInitialValid_abs_df}}
#' 
#' @docType data
#' @name cutoffInitialValid_rel_df
#' @usage data(cutoffs)
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' 
NULL


#' Optimal cutoffs for absolute exposures of initial signatures
#' 
#' Data frame of signature specific cutoff values for the LCD analysis with
#' absolute exposures. All values represent optimal cutoffs. The optimal
#' cutoffs were determined for different choices of parameters in the cost
#' function of the optimization. The row index is equivalent to the ratio
#' between costs for false negative attribution and false positive
#' attribution. The columns correspond to the different signatures. To be
#' used with \code{\link{LCD_complex_cutoff}}. Based on the signatures published
#' in Alexandrov et al. (Nature 2013).
#' 
#' @seealso \code{\link{cutoffCosmicValid_abs_df}}
#' @seealso \code{\link{cutoffInitialValid_rel_df}}
#' 
#' @docType data
#' @name cutoffInitialValid_abs_df
#' @usage data(cutoffs)
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' 
NULL


#' Optimal adapted cutoffs for absolute exposures of initial signatures
#' 
#' Data frame of signature specific cutoff values for the LCD analysis with
#' absolute exposures. All values represent optimal cutoffs. The optimal
#' cutoffs were determined for different choices of parameters in the cost
#' function of the optimization. The row index is equivalent to the ratio
#' between costs for false negative attribution and false positive
#' attribution. The columns correspond to the different signatures. To be
#' used with \code{\link{LCD_complex_cutoff}}. Based on the signatures published
#' in Alexandrov et al. (Nature 2013).
#' 
#' @seealso \code{\link{cutoffCosmicValid_abs_df}}
#' @seealso \code{\link{cutoffInitialValid_rel_df}}
#' 
#' @docType data
#' @name abs_cutoff_adapt_df
#' @usage data(abs_cutoff_adapt)
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' 
NULL


#' Optimal adapted cutoffs for relative exposures of initial signatures
#' 
#' Data frame of signature specific cutoff values for the LCD analysis with
#' relative exposures. All values represent optimal cutoffs. The optimal
#' cutoffs were determined for different choices of parameters in the cost
#' function of the optimization. The row index is equivalent to the ratio
#' between costs for false negative attribution and false positive
#' attribution. The columns correspond to the different signatures. To be
#' used with \code{\link{LCD_complex_cutoff}}. Based on the signatures published
#' in Alexandrov et al. (Nature 2013).
#' 
#' @seealso \code{\link{cutoffCosmicValid_rel_df}}
#' @seealso \code{\link{cutoffInitialValid_abs_df}}
#' 
#' @docType data
#' @name rel_cutoff_adapt_df
#' @usage data(rel_cutoff_adapt)
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' 
NULL


#' Optimal adapted cutoffs for absolute exposures of COSMIC signatures
#' 
#' Data frame of signature specific cutoff values for the LCD analysis with
#' absolute exposures. All values represent optimal cutoffs. The optimal
#' cutoffs were determined for different choices of parameters in the cost
#' function of the optimization. The row index is equivalent to the ratio
#' between costs for false negative attribution and false positive
#' attribution. The columns correspond to the different signatures. To be
#' used with \code{\link{LCD_complex_cutoff}}. Based on the signatures available
#' at http://cancer.sanger.ac.uk/cosmic/signatures .
#' 
#' @seealso \code{\link{cutoffCosmicValid_abs_df}}
#' @seealso \code{\link{cutoffInitialValid_rel_df}}
#' 
#' @docType data
#' @name abs_cutoff_adapt_COSMIC_df
#' @usage data(abs_cutoff_adapt_COSMIC)
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' 
NULL


#' Optimal adapted cutoffs for relative exposures of COSMIC signatures
#' 
#' Data frame of signature specific cutoff values for the LCD analysis with
#' relative exposures. All values represent optimal cutoffs. The optimal
#' cutoffs were determined for different choices of parameters in the cost
#' function of the optimization. The row index is equivalent to the ratio
#' between costs for false negative attribution and false positive
#' attribution. The columns correspond to the different signatures. To be
#' used with \code{\link{LCD_complex_cutoff}}. Based on the signatures available
#' at http://cancer.sanger.ac.uk/cosmic/signatures .
#' 
#' @seealso \code{\link{cutoffCosmicValid_rel_df}}
#' @seealso \code{\link{cutoffInitialValid_abs_df}}
#' 
#' @docType data
#' @name rel_cutoff_adapt_COSMIC_df
#' @usage data(rel_cutoff_adapt_COSMIC)
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' 
NULL


#' Correction factors for different target capture kits
#' 
#' List of lists with correction factors for different target capture kits.
#' The elements of the overall list are lists, every one carrying information
#' for one target capture kit (and namend after it). The elements of these
#' sublists are 64 dimensional vectors with correction factors for all
#' triplets. They were computed using counts of occurence of the respective
#' triplets in the target capture and in the reference genome and making
#' ratios (either for the counts themselves as in \code{abs_cor} or for the
#' relative occurences in \code{rel_cor}). The information in this data
#' structure may be used as input to
#' \code{\link{normalizeMotifs_otherRownames}}.
#' 
#' @docType data
#' @name targetCapture_cor_factors
#' @usage data(targetCapture_cor_factors)
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' 
NULL


#' Colours codes for displaying SNVs
#' 
#' Vector attributing colours to nucleotide exchanges used when displaying SNV
#' information, e.g. in a rainfall plot.
#' 
#' @docType data
#' @name exchange_colour_vector
#' @usage data(exchange_colour_vector)
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' 
NULL
