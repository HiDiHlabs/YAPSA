# Copyright © 2014-2016  The YAPSA package contributors
# This file is part of the YAPSA package. The YAPSA package is licenced under
# GPL-3

#' Construct a VRanges Object from a data frame
#'
#' In this package, big data frames are generated from cohort wide vcf-like
#' files. This function constructs a VRanges object from such a data frame by
#' using \code{\link{makeGRangesFromDataFrame}} from the package
#' \code{\link{GenomicRanges}}
#'
#' @param in_df 
#'  A big dataframe constructed from a vcf-like file of a whole cohort. The 
#'  first columns are those of a standard vcf file, followed by an arbitrary
#'  number of custom or user defined columns. One of these can carry a PID
#'  (patient or sample identifyier) and one can carry subgroup information.
#' @param in_keep.extra.columns in_seqinfo
#'  Argument passed on to \code{\link{makeGRangesFromDataFrame}}
#' @param in_seqinfo
#'  A seqInfo object, referring to the reference genome used.
#'  Argument passed on to \code{\link{makeGRangesFromDataFrame}}
#' @param in_seqnames.field
#'  Indicates the name of the column in which the chromosome is encoded
#' @param in_start.field
#'  Indicates the name of the column in which the start coordinate is
#'  encoded
#' @param in_end.field
#'  Indicates the name of the column in which the end coordinate is
#'  encoded
#' @param in_PID.field
#'  Indicates the name of the column in which the PID (patient or sample
#'  identifier) is encoded
#' @param in_subgroup.field
#'  Indicates the name of the column in which the subgroup information
#'  is encoded
#' @param in_strand.field
#'  Indicates the name of the column in which the strandedness is encoded
#' @param verbose_flag
#'  Verbose if 1
#'
#' @return
#'  The constructed VRanges object
#' 
#' @examples
#'  data(lymphoma_test)
#'  temp_vr <- makeVRangesFromDataFrame(lymphoma_test_df,
#'                                      in_seqnames.field="CHROM",
#'                                      in_subgroup.field="SUBGROUP",
#'                                      verbose_flag=1)
#' 
#' @seealso \code{\link{makeGRangesFromDataFrame}}
#' 
#' @import GenomicRanges
#' @import VariantAnnotation
#' @importFrom GenomeInfoDb seqlengths seqlengths<-
#' @export
#' 
makeVRangesFromDataFrame <- function(in_df,in_keep.extra.columns=TRUE,
                                     in_seqinfo=NULL,
                                     in_seqnames.field="X.CHROM",
                                     in_start.field="POS",
                                     in_end.field="POS",in_PID.field="PID",
                                     in_subgroup.field="subgroup",
                                     in_strand.field="strand",
                                     verbose_flag=1) {
  out_vr <- NULL
  name_list <- tolower(names(in_df))
  match_list <- c("ref","alt",tolower(in_seqnames.field),
                  tolower(in_start.field),tolower(in_end.field))
  if(length(which(match_list %in% name_list)) == length(match_list)) {
    my_gr <- GenomicRanges::makeGRangesFromDataFrame(in_df,
                                      keep.extra.columns=in_keep.extra.columns,
                                      seqinfo=in_seqinfo,
                                      seqnames.field=in_seqnames.field,
                                      start.field=in_start.field,
                                      end.field=in_end.field,
                                      strand.field=in_strand.field)
    if(!("+" %in% unique(strand(my_gr))) & !("-" %in% unique(strand(my_gr)))){
      strand(my_gr) <- "+"
      if(verbose_flag==1){
        cat("YAPSA:::makeVRangesFromDataFrame::warning:strand ",
            "information missing, set to \"+\".\n")
      }
    }    
    out_vr <- VRanges(seqnames=seqnames(my_gr),ranges=ranges(my_gr),
                      ref=my_gr$REF,alt=my_gr$ALT)
    if(tolower(in_PID.field) %in% name_list) {
      if(verbose_flag==1){
        cat("YAPSA:::makeVRangesFromDataFrame::in_PID.field found. ",
            "Retrieving PID information.\n")
      }
      column_ind <- 
        min(which(tolower(names(mcols(my_gr)))==tolower(in_PID.field)))
      out_vr$PID <- mcols(my_gr)[,column_ind]
      sampleNames(out_vr) <- mcols(my_gr)[,column_ind]
    } else if("pid" %in% name_list) {
      if(verbose_flag==1){
        cat(paste0("YAPSA:::makeVRangesFromDataFrame::warning:in_PID.field ",
                   "not a valid column name,",
                   " but default is valid. Retrieving PID information.\n"))
      }
      out_vr$PID <- mcols(my_gr)[,"PID"]
      sampleNames(out_vr) <- mcols(my_gr)[,column_ind]
    } else {
      if(verbose_flag==1){
        cat("YAPSA:::makeVRangesFromDataFrame::warning:PID information ",
            "missing. Filling up with dummy entries.\n");}
      out_vr$PID <- "dummy_PID"
      sampleNames(out_vr) <- "dummyPID"
    }
    if(tolower(in_subgroup.field) %in% name_list) {
      if(verbose_flag==1){
        cat("YAPSA:::makeVRangesFromDataFrame::in_subgroup.field found. ",
            "Retrieving subgroup information.\n")
      }
      column_ind <- 
        min(which(tolower(names(mcols(my_gr)))==tolower(in_subgroup.field)))
      #out_vr$Type <- mcols(my_gr)[,column_ind]
      mcols(out_vr)[,in_subgroup.field] <- mcols(my_gr)[,column_ind]
    } else if("subgroup" %in% name_list) {
      if(verbose_flag==1){
        cat(paste0("YAPSA:::makeVRangesFromDataFrame::warning:",
                   "in_subgroup.field not a valid column name, ",
                   "but default is valid. Retrieving subgroup information.\n"))
      }
      #out_vr$Type <- my_gr$subgroup
      mcols(out_vr)[,in_subgroup.field] <- my_gr$subgroup
    } else {
      if(verbose_flag==1){cat("YAPSA:::makeVRangesFromDataFrame::warning:",
                              "subgroup information missing. ",
                              "Filling up with dummy entries.\n");}
      #out_vr$Type <- "dummy_subgroup"
      mcols(out_vr)[,in_subgroup.field] <- "dummy_subgroup"
    }
    seqlengths(out_vr) <- seqlengths(my_gr)    
  } else {
    if(verbose_flag==1){
      cat("YAPSA:::makeVRangesFromDataFrame::error:mismatch in column names, ",
          "return NULL.\n")
    }
  }
  return(out_vr)
}


#' Create a Mutational Catalogue from a VRanges Object
#'
#' This function creates a mutational catalogue from a VRanges Object by
#' first calling \code{\link[SomaticSignatures]{mutationContext}} to establish
#' the motif context of the variants in the input VRanges and then calling
#' \code{\link[SomaticSignatures]{motifMatrix}} to build the
#' mutational catalogue \code{V}.
#'
#' @param in_vr 
#'  A VRanges object constructed from a vcf-like file of a whole cohort. The 
#'  first columns are those of a standard vcf file, followed by an arbitrary
#'  number of custom or used defined columns. One of these can carry a PID
#'  (patient or sample identifyier) and one can carry subgroup information.
#' @param in_refGenome
#'  The reference genome handed over to 
#'  \code{\link[SomaticSignatures]{mutationContext}}
#'  and used to extract the motif context of the variants in \code{in_vr}.
#' @param in_wordLength
#'  The size of the motifs to be extracted by 
#'  \code{\link[SomaticSignatures]{mutationContext}}
#' @param in_PID.field
#'  Indicates the name of the column in which the PID (patient or sample
#'  identifier) is encoded
#' @param in_verbose
#'  Verbose if \code{in_verbose=1}
#' @param in_rownames
#'  Optional parameter to specify rownames of the mutational catalogue \code{V}
#'  i.e. the names of the features.
#' @param adapt_rownames
#'  Rownames of the output \code{matrix} will be adapted if 
#'  \code{adapt_rownames=1}
#'
#' @return A list with entries
#'  \code{matrix},
#'  \code{frame},
#' \itemize{
#'  \item \code{matrix}:
#'    The mutational catalogue \code{V}
#'  \item \code{frame}:
#'    Additional and meta information on rownames (features), colnames (PIDs)
#'    and subgroup attribution.
#' }
#' 
#' @examples
#'  library(BSgenome.Hsapiens.UCSC.hg19)
#'  data(lymphoma_test)
#'  data(sigs)
#'  word_length <- 3
#'  temp_vr <- makeVRangesFromDataFrame(
#'    lymphoma_test_df,in_seqnames.field="CHROM",
#'    in_subgroup.field="SUBGROUP",verbose_flag=1)
#'  temp_list <- create_mutation_catalogue_from_VR(
#'    temp_vr,in_refGenome=BSgenome.Hsapiens.UCSC.hg19,
#'    in_wordLength=word_length,in_PID.field="PID",
#'    in_verbose=1)
#'  dim(temp_list$matrix)
#'  head(temp_list$matrix)
#'  test_list <- split(lymphoma_test_df,f=lymphoma_test_df$PID)
#'  other_list <- list()
#'  for(i in seq_len(length(test_list))){
#'    other_list[[i]] <- test_list[[i]][c(1:80),]
#'  }
#'  other_df <- do.call(rbind,other_list)
#'  other_vr <- makeVRangesFromDataFrame(
#'    other_df,in_seqnames.field="CHROM",
#'    in_subgroup.field="SUBGROUP",verbose_flag=1)
#'  other_list <- create_mutation_catalogue_from_VR(
#'    other_vr,in_refGenome=BSgenome.Hsapiens.UCSC.hg19,
#'    in_wordLength=word_length,in_PID.field="PID",
#'    in_verbose=1,in_rownames=rownames(AlexCosmicValid_sig_df)) 
#'  dim(other_list$matrix)
#'  head(other_list$matrix)
#' 
#' @seealso \code{\link[SomaticSignatures]{mutationContext}}
#' @seealso \code{\link[SomaticSignatures]{motifMatrix}}
#' 
#' @importFrom SomaticSignatures mutationContext
#' @importFrom SomaticSignatures motifMatrix
#' @export
#' 
create_mutation_catalogue_from_VR <- function(in_vr,in_refGenome,in_wordLength,
                                              in_PID.field="PID",in_verbose=0,
                                              in_rownames=c(),adapt_rownames=1
                                              ) {
  # adapt levels of in_PID.field
  PID_index <- which(names(mcols(in_vr))==in_PID.field)
  names(mcols(in_vr))[PID_index] <- "PID"
  in_vr$PID <- factor(as.character(in_vr$PID))
  # now start real work of the function
  if ( in_verbose==1 ) {
    cat("YAPSA:::create_mutation_catalogue_from_VR::Extracting mutation ",
        "context for the SNVs from the reference genome...\n");
  }
  my_motifs <- mutationContext(in_vr, in_refGenome, 
                               k=in_wordLength, unify=TRUE)
  if ( in_verbose==1 ) {
    cat("YAPSA:::create_mutation_catalogue_from_VR::Building mutational ",
        "catalogue...\n");
    }
  ## cave: mutationContextMatrix is the old version of the SomaticSignatures 
  ##       package motifMatrix is the new version, but return a list instead 
  ##       of a matrix when the data is sparse
  #my_matrix = mutationContextMatrix(my_motifs, group=in_PID.field,
  #                                  normalize = FALSE)
  my_matrix = motifMatrix(my_motifs, group=in_PID.field,normalize = FALSE)
  ## check if unallowed nucleotides in matrix
  deselect_ind <- grep("[NX]",rownames(my_matrix))
  if (length(deselect_ind) > 0) {
    my_matrix <- my_matrix[-deselect_ind]
  }
  # account for sparsity of the matrix
  if(length(in_rownames) > 0) {
    if(length(in_rownames) < dim(my_matrix)[1]) {
      if(in_verbose==1) {
        cat("YAPSA:::create_mutation_catalogue_from_VR::Warning: supplied ",
            "rownames too short\n");
      }
    } else if (!is.character(in_rownames)) {
      if(in_verbose==1) {
        cat("YAPSA:::create_mutation_catalogue_from_VR::Warning: supplied ",
            "rownames not character format\n");
      }
    } else {
      my_rownames <- transform_rownames_MATLAB_to_R(in_rownames,in_wordLength)
      temp_matrix <- matrix(0,nrow=length(my_rownames),ncol=dim(my_matrix)[2])
      match_ind <- match(rownames(my_matrix),my_rownames)
      if( (length(match_ind) < dim(my_matrix)[1]) | (any(is.na(match_ind)))){
        if(in_verbose==1) {
          cat("YAPSA:::create_mutation_catalogue_from_VR::Warning: my_matrix",
              " contains rownames which are not in supplied rownames.\n");
        }
      } else {
        temp_colnames <- colnames(my_matrix)
        temp_matrix[match_ind,] <- my_matrix
        my_matrix <- temp_matrix
        rownames(my_matrix) <- my_rownames
        colnames(my_matrix) <- temp_colnames
      }
    }
  }
  # string manipulations for export to matlab
  if ( in_verbose==1 ) {
    cat("YAPSA:::create_mutation_catalogue_from_VR::Performing string ",
        "manipulations for output to matlab...\n");
  }
  temp_frame <- data.frame(strsplit(rownames(my_matrix)," "))
  colnames(temp_frame) <- rownames(my_matrix)
  my_frame <- as.data.frame(t(temp_frame[,]))
  colnames(my_frame) <- c("exchange","context")
  my_frame$types <- 
    paste0(substr(my_frame$exchange,1,1),">",substr(my_frame$exchange,2,2))
  ## adapt for word length here
  fraction_length <- floor(in_wordLength/2)
  my_frame$subtypes <- 
    paste0(substr(my_frame$context,1,fraction_length),
           substr(my_frame$exchange,1,1),
           substr(my_frame$context,
                  in_wordLength-fraction_length+1,in_wordLength))  
  my_frame$newnames <- paste0(my_frame$types," ",my_frame$subtypes)
  if(adapt_rownames==1) {
    rownames(my_matrix) <- my_frame$newnames
  }
  return(list(matrix=my_matrix,frame=my_frame))
}


#' Create a Mutational Catalogue from a data frame
#'
#' This function creates a mutational catalogue from a data frame. It is a
#' wrapper function for \code{\link{create_mutation_catalogue_from_VR}}:
#' it first creates a VRanges object from the data frame by
#' \code{\link{makeVRangesFromDataFrame}} and then passes this object on to the
#' above mentioned custom function.
#'
#' @param this_df 
#'  A data frame constructed from a vcf-like file of a whole cohort. The 
#'  first columns are those of a standard vcf file, followed by an arbitrary
#'  number of custom or used defined columns. One of these can carry a PID
#'  (patient or sample identifyier) and one can carry subgroup information.
#' @param this_refGenome_Seqinfo
#'  A seqInfo object, referring to the reference genome used.
#'  Argument passed on to \code{\link{makeGRangesFromDataFrame}} and thus 
#'  indirectly to \code{\link{makeGRangesFromDataFrame}}.
#' @param this_seqnames.field
#'  Indicates the name of the column in which the chromosome is encoded
#' @param this_start.field
#'  Indicates the name of the column in which the start coordinate is
#'  encoded
#' @param this_end.field
#'  Indicates the name of the column in which the end coordinate is
#'  encoded
#' @param this_PID.field
#'  Indicates the name of the column in which the PID (patient or sample
#'  identifier) is encoded
#' @param this_subgroup.field
#'  Indicates the name of the column in which the subgroup information
#'  is encoded
#' @param this_refGenome
#'  The reference genome handed over to 
#'  \code{\link{create_mutation_catalogue_from_VR}}
#'  and indirectly to \code{\link[SomaticSignatures]{mutationContext}}
#'  and used to extract the motif context of the variants in \code{in_vr}.
#' @param this_wordLength
#'  The size of the motifs to be extracted by 
#'  \code{\link[SomaticSignatures]{mutationContext}}
#' @param this_verbose
#'  Verbose if \code{this_verbose=1}
#' @param this_rownames
#'  Optional parameter to specify rownames of the mutational catalogue \code{V}
#'  i.e. the names of the features.
#' @param this_adapt_rownames
#'  Rownames of the output \code{matrix} will be adapted if 
#'  \code{this_adapt_rownames=1}
#'
#' @return A list with entries
#'  \code{matrix} and
#'  \code{frame}
#'  obtained from \code{\link{create_mutation_catalogue_from_VR}}:
#' \itemize{
#'  \item \code{matrix}:
#'    The mutational catalogue \code{V}
#'  \item \code{frame}:
#'    Additional and meta information on rownames (features), colnames (PIDs)
#'    and subgroup attribution.
#' }
#' 
#' @examples
#'  library(BSgenome.Hsapiens.UCSC.hg19)
#'  data(lymphoma_test)
#'  word_length <- 3
#'  temp_list <- create_mutation_catalogue_from_df(
#'    lymphoma_test_df,this_seqnames.field = "CHROM",
#'    this_start.field = "POS",this_end.field = "POS",
#'    this_PID.field = "PID",this_subgroup.field = "SUBGROUP",
#'    this_refGenome = BSgenome.Hsapiens.UCSC.hg19,
#'    this_wordLength = word_length)
#'    dim(temp_list$matrix)
#'    head(temp_list$matrix)
#' 
#' @seealso \code{\link{makeVRangesFromDataFrame}}
#' @seealso \code{\link{create_mutation_catalogue_from_VR}}
#' 
#' @export
#'
create_mutation_catalogue_from_df <- function(this_df,
                                              this_refGenome_Seqinfo=NULL,
                                              this_seqnames.field="X.CHROM",
                                              this_start.field="POS",
                                              this_end.field="POS",
                                              this_PID.field="PID",
                                              this_subgroup.field="subgroup",
                                              this_refGenome,this_wordLength,
                                              this_verbose=1,
                                              this_rownames=c(),
                                              this_adapt_rownames=1) {
  if(is.null(this_refGenome_Seqinfo)){
    this_refGenome_Seqinfo <- seqinfo(this_refGenome)
  }
  my_vr <- makeVRangesFromDataFrame(this_df,
                                    in_keep.extra.columns=TRUE,
                                    in_seqinfo=this_refGenome_Seqinfo,
                                    in_seqnames.field=this_seqnames.field,
                                    in_start.field=this_start.field,
                                    in_end.field=this_end.field,
                                    in_PID.field=this_PID.field,
                                    in_subgroup.field=this_subgroup.field,
                                    verbose_flag=this_verbose)
  my_result_list <- create_mutation_catalogue_from_VR(my_vr,this_refGenome,
                                                      this_wordLength,
                                                      this_PID.field,
                                                      this_verbose,
                                                      this_rownames,
                                                      this_adapt_rownames)
  return(my_result_list)
}


#' Change rownames from one naming convention to another
#'
#' Rownames or names of the features used differ between the different contexts
#' a signature analysis is carried out in. The function
#' \code{transform_rownames_R_to_MATLAB} changes from the convention used in
#' the YAPSA pacakge to the one used by Alexandrov et al. in the MATLAB
#' framework.
#'
#' @param in_rownames 
#'  Character vector of input rownames
#' @param wordLength 
#'  Size of the considered motif context
#'
#' @return A character vector of the translated rownames.
#' 
#' @examples
#' NULL
#' 
#' @export
#'
transform_rownames_R_to_MATLAB <- function(in_rownames,wordLength=3) {
  out_rownames <- in_rownames
  temp_frame <- data.frame(strsplit(in_rownames," "))
  colnames(temp_frame) <- in_rownames
  my_frame <- as.data.frame(t(temp_frame[,]))
  colnames(my_frame) <- c("exchange","context")
  my_frame$types <- paste0(substr(my_frame$exchange,1,1),">",
                           substr(my_frame$exchange,2,2))
  ## adapt for word length here
  fraction_length <- floor(wordLength/2)
  my_frame$subtypes <- paste0(
    substr(my_frame$context,1,fraction_length),
    substr(my_frame$exchange,1,1),
    substr(my_frame$context,wordLength-fraction_length+1,wordLength))  
  out_rownames <- paste0(my_frame$types," ",my_frame$subtypes)
  return(out_rownames)
}


#' Change rownames from one naming convention to another
#'
#' The function \code{transform_rownames_MATLAB_to_R} changes from the 
#' convention used in Alexandrov et al. in the MATLAB framework to the one used
#' by the YAPSA pacakge.
#' 
#' @export
#' @rdname transform_rownames_R_to_MATLAB
#'
transform_rownames_MATLAB_to_R <- function(in_rownames,wordLength=3) {
  out_rownames <- in_rownames
  temp_frame <- data.frame(strsplit(in_rownames," "))
  colnames(temp_frame) <- in_rownames
  my_frame <- as.data.frame(t(temp_frame[,]))
  colnames(my_frame) <- c("exchange","context")
  my_frame$types <- paste0(substr(my_frame$exchange,1,1),
                           substr(my_frame$exchange,3,3))
  ## adapt for word length here
  fraction_length <- floor(wordLength/2)
  my_frame$subtypes <- paste0(
    substr(my_frame$context,1,fraction_length),".",
    substr(my_frame$context,wordLength-fraction_length+1,wordLength))  
  out_rownames <- paste0(my_frame$types," ",my_frame$subtypes)
  return(out_rownames)
}


#' Change rownames from one naming convention to another
#'
#' The function \code{transform_rownames_MATLAB_to_R} changes from the 
#' convention used in stored mutational catalogues by Alexandrov et al. to the 
#' one used by the YAPSA pacakge.
#'
#' @export
#' @rdname transform_rownames_R_to_MATLAB
#'
transform_rownames_nature_to_R <- function(in_rownames,wordLength=3) {
  out_rownames <- in_rownames
  ## adapt for word length here
  fraction_length <- floor(wordLength/2)
  string_length <- nchar(as.character(in_rownames[1]))
  my_frame <- repeat_df(in_value = "a",in_rows = length(in_rownames),
                        in_cols = 2)
  colnames(my_frame) <- c("types","subtypes")
  my_frame$subtypes <- paste0(
    substr(in_rownames,1,fraction_length),
    substr(in_rownames,fraction_length+2,fraction_length+2),
    substr(in_rownames,string_length-fraction_length+1,string_length))  
  my_frame$types <- substr(in_rownames,fraction_length+2,fraction_length+4)
  out_rownames <- paste0(my_frame$types," ",my_frame$subtypes)
  return(out_rownames)
}


#' Change rownames from one naming convention to another
#'
#' The function \code{transform_rownames_YAPSA_to_deconstructSigs} changes from
#' the convention used in the YAPSA package to the one used by the 
#' deconstructSigs package.
#'
#' @export
#' @rdname transform_rownames_R_to_MATLAB
#'
transform_rownames_YAPSA_to_deconstructSigs <- function(in_rownames,
                                                        wordLength=3) {
  out_rownames <- in_rownames
  ## adapt for word length here
  fraction_length <- floor(wordLength/2)
  string_length <- nchar(as.character(in_rownames[1]))
  my_frame <- repeat_df(in_value = "a",in_rows = length(in_rownames),
                        in_cols = 2)
  colnames(my_frame) <- c("types","subtypes")
  my_frame$subtypes <- lapply(in_rownames,function(x) strsplit(x," ")[[1]][2])
  my_frame$types <- lapply(in_rownames,function(x) strsplit(x," ")[[1]][1])
  out_rownames <- paste0(
    substr(my_frame$subtypes,1,fraction_length),"[",
    my_frame$types,"]",
    substr(my_frame$subtypes,fraction_length+2,2*fraction_length+1))
  return(out_rownames)
}


#' Change rownames from one naming convention to another
#'
#' The function \code{transform_rownames_YAPSA_to_deconstructSigs} changes from
#' the convention used in the deconstructSigs package to the one used by the 
#' YAPSA pacakge.
#'
#' @export
#' @rdname transform_rownames_R_to_MATLAB
#'
transform_rownames_deconstructSigs_to_YAPSA <- function(in_rownames,
                                                        wordLength=3) {
  out_rownames <- in_rownames
  ## adapt for word length here
  fraction_length <- floor(wordLength/2)
  string_length <- nchar(as.character(in_rownames[1]))
  my_frame <- repeat_df(in_value = "a",in_rows = length(in_rownames),
                        in_cols = 2)
  colnames(my_frame) <- c("types","subtypes")
  my_frame$subtypes <- paste0(
    substr(in_rownames,1,fraction_length),
    substr(in_rownames,2+fraction_length,2+fraction_length),
    substr(in_rownames,6+fraction_length,5+2*fraction_length))
  my_frame$types <- substr(in_rownames,2+fraction_length,4+fraction_length)
  out_rownames <- paste0(my_frame$types," ",my_frame$subtypes)  
  return(out_rownames)
}


#' Normalize Somatic Motifs with different rownames
#' 
#' This is a wrapper function to 
#' \code{\link[SomaticSignatures]{normalizeMotifs}}. The rownames
#' are first transformed to fit the convention of the 
#' \code{\link[SomaticSignatures]{SomaticSignatures}}
#' package and then passed on to the above mentioned function.
#' 
#' @param in_matrix,in_norms
#'  Arguments to \code{\link[SomaticSignatures]{normalizeMotifs}}
#' @param adjust_counts
#'  Whether to rescale the counts after adaption or not. Default is true.
#'
#' @return
#'  The matrix returned by \code{\link[SomaticSignatures]{normalizeMotifs}}, 
#'  but with rownames transformed back to the convention of the input
#'  
#' @examples
#'  NULL
#'  
#' @importFrom SomaticSignatures normalizeMotifs
#' @export
#'  
normalizeMotifs_otherRownames <- function(in_matrix,
                                          in_norms,
                                          adjust_counts=TRUE) {
  my_matrix <- in_matrix
  rownames(my_matrix) <- transform_rownames_MATLAB_to_R(rownames(in_matrix))
  out_matrix <- SomaticSignatures::normalizeMotifs(my_matrix,in_norms)
  rownames(out_matrix) <- rownames(in_matrix)
  if(adjust_counts){
    total_counts <- colSums(in_matrix)
    intermediate_counts <- colSums(out_matrix)
    mult_factors <- total_counts/intermediate_counts
    out_matrix <- as.data.frame(t(t(out_matrix)*mult_factors))
  }
  return(out_matrix)
}


adjust_number_of_columns_in_list_of_catalogues <- 
  function(in_all_list, in_merged_results_list) {
  out_list <- list()
  number_of_strata <- length(in_all_list$results_lists_list)
  reference_matrix <- in_merged_results_list$matrix
  zero_df <- as.data.frame(reference_matrix-reference_matrix)
  colname_vector <- colnames(reference_matrix)
  for (i in seq_len(number_of_strata)) {
    temp_df <- as.data.frame(in_all_list$results_lists_list[[i]]$matrix)
    match_ind <- match(names(temp_df),colname_vector)
    new_df <- zero_df
    new_df[,match_ind] <- temp_df
    out_list[[i]] <- new_df
  }
  return(out_list)
}


save_mutation_catalogue <- function(in_table,in_matrix,in_frame,
                                    in_subgroup_flag,in_meta_info,in_subDir,
                                    in_verbose) {
  type_df <- aggregate(data=in_table,subgroup~PID, function(l) return(l[1]))
  if (dim(type_df)[1] != length(colnames(in_matrix))) {
    cat("YAPSA:::save_mutation_catalogue::Warning: mismatch between subgroup",
        " information and colnames(in_matrix).\n");
  }
  subgroup_vector <- rep("dummy",dim(in_matrix)[2])
  if (in_subgroup_flag==1) {
    subgroup_vector <- 
      as.character(in_meta_info$Diagnosis[match(colnames(in_matrix),
                                                in_meta_info$PID)])
  }
  if (in_subgroup_flag==2) {
    subgroup_vector <- 
      as.character(type_df$subgroup[match(colnames(in_matrix),type_df$PID)])
  }
  # store to text files
  if ( in_verbose==1 ) {cat("Store to specified files...\n");}
  save(in_table,in_matrix,file=file.path(in_subDir,"MotifMatrix.RData"))
  write.table(in_matrix,sep="\t",file=file.path(in_subDir,"MotifMatrix.txt"),
              col.names=FALSE,row.names=FALSE)
  write(c("PIDs",colnames(in_matrix)),sep="\t",
        file=file.path(in_subDir,"MotifMatrix_colnames.txt"))
  write(c("subgroup",subgroup_vector),sep="\t",
        file=file.path(in_subDir,"MotifMatrix_subgroups.txt"))
  write(c("types",in_frame$types),sep="\t",
        file=file.path(in_subDir,"MotifMatrix_types.txt"))
  write(c("subtypes",in_frame$subtypes),sep="\t",
        file=file.path(in_subDir,"MotifMatrix_subtypes.txt"))
}


stratify_vcf_like_df <- function(in_table,in_column_name,in_verbose=0) {
  out_table_list <- list()
  out_name_list <- list()
  stratification_index <- which(names(in_table)==in_column_name)
  name_list <- unique(in_table[,stratification_index])
  index_range <- sort(name_list)
  counter <- 0
  if(in_verbose == 1) {
    cat("\nYAPSA:::stratify_vcf_like_df::Stratify...\n")
  }
  for (i in index_range) {
    counter <- counter + 1
    temp_ind <- which(in_table[,stratification_index]==i)
    temp_table <- in_table[temp_ind,]
    temp_name <- i
    if(in_verbose == 1) {
      cat("YAPSA:::stratify_vcf_like_df:: temp_name =",temp_name,
          "; number of variants:",dim(temp_table)[1],"; number of PIDs:",
          length(unique(temp_table$PID)),"\n")
    }
    out_table_list[[counter]] <- temp_table
    out_name_list[[counter]] <- temp_name
  }    
  return(list(table_list=out_table_list,name_list=out_name_list))
}


save_stratified_vcf_like_df <- function(in_table_list,in_name_list,
                                        in_target_dir,in_verbose) {
  number_of_strata <- length(in_table_list)
  subDir <- file.path(in_target_dir, "vcf_like")
  strataDir <- file.path(in_target_dir, "strata")
  dir.create(in_target_dir)
  dir.create(subDir)
  dir.create(strataDir)
  if(in_verbose == 1) {
    cat("\nYAPSA:::save_stratified_vcf_like_df::Write to file...\n")
  }
  for (i in seq_len(number_of_strata)) {
    temp_table <- in_table_list[[i]]
    temp_name <- in_name_list[[i]]
    if(in_verbose == 1) {
      cat("YAPSA:::save_stratified_vcf_like_df:: temp_name =",temp_name,
          "; number of variants:",dim(temp_table)[1],"; number of PIDs:",
          length(unique(temp_table$PID)),"\n")
    }
    file_name <- paste0(subDir,"/somatic_SNVs_",temp_name,".tsv")
    names(temp_table) <- sub("X.CHROM","\\#CHROM",names(temp_table))
    write.table(temp_table, file=file_name,sep="\t",
                quote =FALSE, row.names=FALSE)
    temp_dir <- file.path(in_target_dir, "strata", temp_name)
    dir.create(temp_dir)
  }    
  return(1)
}


stratify_and_create_mutational_catalogue <- 
  function(our_table, our_column_name, vcf_target_dir, strata_target_dir,
           our_refGenome_Seqinfo, our_seqnames.field="X.CHROM",
           our_start.field="POS", our_end.field="POS", our_PID.field="PID",
           our_subgroup.field="subgroup", our_refGenome, our_wordLength,
           in_verbose=1, our_rownames=c()) {
  results_lists_list <- list()
  results_list <- stratify_vcf_like_df(our_table,our_column_name,in_verbose)
  table_list <- results_list$table_list
  name_list <- results_list$name_list
  if(!is.null(vcf_target_dir)) {
    status <- save_stratified_vcf_like_df(table_list,name_list,
                                          vcf_target_dir,in_verbose)    
  }
  number_of_strata <- length(table_list)
  if(in_verbose == 1) {
    cat("\nYAPSA:::stratify_and_create_mutational_catalogue::Create ",
        "mutational catalogues...\n")
  }
  for (i in seq_len(number_of_strata)) {
    if(in_verbose == 1) {
      cat("\nYAPSA:::stratify_and_create_mutational_catalogue::Processing ",
          "stratum: ",i,"\n")
    }
    temp_table <- table_list[[i]]
    temp_name <- name_list[[i]]
    temp_results_list <- create_mutation_catalogue_from_df(
      temp_table,
      our_refGenome_Seqinfo,
      this_seqnames.field=our_seqnames.field,
      this_start.field=our_start.field,
      this_end.field=our_end.field,
      this_PID.field=our_PID.field,
      this_subgroup.field=our_subgroup.field,
      our_refGenome,our_wordLength,in_verbose,
      our_rownames)
    if(!is.null(strata_target_dir)) {
      temp_matrix <- temp_results_list$matrix
      temp_frame <- temp_results_list$frame
      MotifMatrix_dir <- file.path(strata_target_dir,temp_name)
      save_mutation_catalogue(temp_table,temp_matrix,temp_frame,2,NULL,
                              MotifMatrix_dir,in_verbose)
    }
    results_lists_list[[i]] <- temp_results_list
  }    
  return(list(table_list=table_list,name_list=name_list,
              results_lists_list=results_lists_list))  
}


#' Wrapper for cut
#' 
#' In this wrapper function for the known \code{\link{cut}} function, the
#' \code{breaks} vector need not be supplied directly, instead, for every
#' break, an interval is supplied and the function optimizes the choice of
#' the breakpoint by chosing a local minimum of the distribution.
#' 
#' @param in_vector
#'  Vector of numerical continuously distributed input
#' @param in_outlier_cutoffs
#'  Interval specifyinf the upper and lower bounds of the range to be
#'  considered
#' @param in_cutoff_ranges_list
#'  List if intervals in which the cutoffs for \code{\link{cut}} have to be
#'  optimized.
#' @param in_labels
#'  Labels assigned to the strata or factors returned
#' @param in_name
#'  String specifying the name of the quantity analyzed (and plotted on the
#'  x-axis of the figure to be created).
#' @param output_path
#'  Path where the figure produced by the density function should be stored
#'  if non-NULL.
#'
#' @return A list with entries
#'  \code{category_vector}, and
#'  \code{density_plot} and
#'  \code{cutoffs}
#' \itemize{
#'  \item \code{category_vector}:
#'    Factor vector of the categories or strata, of the same length as
#'    \code{in_vector}
#'  \item \code{density_plot}:
#'    Density plot produced by the density function and indication of the
#'    chosen cutoffs.
#'  \item \code{cutoffs}:
#'    Vector of the computed optimal cutoffs
#' }
#' 
#' @examples
#'  data(lymphoma_test)
#'  lymphoma_test_df$random_norm <- rnorm(dim(lymphoma_test_df)[1])
#'  temp_list <- cut_breaks_as_intervals(
#'    lymphoma_test_df$random_norm,
#'    in_outlier_cutoffs=c(-4,4),
#'    in_cutoff_ranges_list=list(c(-2.5,-1.5),c(0.5,1.5)),
#'    in_labels=c("small","intermediate","big"))
#'   temp_list$density_plot
#' 
#' @seealso \code{\link{cut}}
#' @seealso \code{\link[stats]{density}}
#' 
#' @export
#' 
cut_breaks_as_intervals <- function(in_vector,
                                    in_outlier_cutoffs=c(0,3000),
                                    in_cutoff_ranges_list=list(c(60,69),
                                                               c(25,32)),
                                    in_labels=c("late","intermediate","early"),
                                    in_name="",
                                    output_path=NULL) {
  ## prepare
  .e <- environment()
  finite_ind <- which(is.finite(in_vector))
  finite_vector <- in_vector[finite_ind]
  choice_ind <- which(finite_vector<=in_outlier_cutoffs[2] &
                        finite_vector>=in_outlier_cutoffs[1])
  my_density <- density(finite_vector[choice_ind])
  my_density_df <- data.frame(x=my_density$x,y=my_density$y)
  ## find appropriate cutoffs in the ranges provided (refine!!)
  number_of_cutoffs <- length(in_cutoff_ranges_list)
  my_cutoff_vector <- rep(0,number_of_cutoffs)
  for(j in seq_len(number_of_cutoffs)) {
    my_cutoff_interval <- in_cutoff_ranges_list[[j]]
    my_range_ind <- which(my_density_df$x<my_cutoff_interval[2] & 
                            my_density_df$x>my_cutoff_interval[1])
    my_min_ind <- which.min(my_density_df$y[my_range_ind])
    my_cutoff_vector[j] <- my_density_df$x[my_range_ind[my_min_ind]]
  }
  ## make density plot
  p <- ggplot(my_density_df,aes_string(x="x",y="y"), environment = .e) + 
    geom_line() + labs(x=in_name, y="density") + 
    geom_vline(xintercept=my_cutoff_vector,colour="red")
  if(!is.null(output_path)) {
    png(output_path,width=400,height=400); p; dev.off()    
  }
  my_cutoff_vector <- c(in_outlier_cutoffs[1],my_cutoff_vector,
                        in_outlier_cutoffs[2])
  ## discretize into categories
  cut_vector <- cut(in_vector,my_cutoff_vector,right=FALSE,labels=in_labels)
  cut_vector <- as.character(cut_vector)
  cut_vector[which(is.na(cut_vector))] <- "undetermined"
  cut_vector <- factor(cut_vector)
  return(list(category_vector=cut_vector,density_plot=p,
              cutoffs=my_cutoff_vector))
}


#' Return those PIDs which have an extreme pattern for signature exposure
#' 
#' For all signatures found in a project, this function returns the sample
#' identifiers (PIDs) with extremely high or extremely low exposures of
#' the respective signatures.
#' 
#' @param in_exposures_df
#'  Data frame with the signature exposures
#' @param in_quantile
#'  Quantile for the amount of extreme PIDs to be selected.
#' 
#' @return A data frame with 4 rows per signature (high PIDs, high
#'  exposures, low PIDs, low exposures); the number of columns depends on
#'  the quantile chosen.
#' 
#' @examples
#'  data(lymphoma_cohort_LCD_results)
#'  get_extreme_PIDs(lymphoma_Nature2013_COSMIC_cutoff_exposures_df,0.05)
#'
#' @import ggplot2
#' @export
#'
get_extreme_PIDs <- function(in_exposures_df,in_quantile=0.03) {
  my_number_of_sigs <- dim(in_exposures_df)[1]
  my_number_of_PIDs <- dim(in_exposures_df)[2]
  my_number_of_extremes <- ceiling(in_quantile*my_number_of_PIDs)
  extreme_PID_df <- repeat_df(0,my_number_of_extremes,my_number_of_sigs*4)
  for(i in seq_len(my_number_of_sigs)) {
    this_sig <- rownames(in_exposures_df)[i]
    my_exposure_vector <- in_exposures_df[i,]
    my_breaks <- as.numeric(quantile(my_exposure_vector,
                                     probs=c(in_quantile,(1-in_quantile))))
    extreme_low_ind <- which(my_exposure_vector<=my_breaks[1])
    extreme_high_ind <- which(my_exposure_vector>=my_breaks[2])
    extreme_PID_df[,4*i-3] <- names(
      my_exposure_vector[extreme_low_ind[seq_len(my_number_of_extremes)]])
    names(extreme_PID_df)[4*i-3] <- paste0(this_sig,"_low_PIDs")
    extreme_PID_df[,4*i-2] <- as.numeric(
      my_exposure_vector[extreme_low_ind[seq_len(my_number_of_extremes)]])
    names(extreme_PID_df)[4*i-2] <- paste0(this_sig,"_low_exp")
    extreme_PID_df[,4*i-1] <- names(
      my_exposure_vector[extreme_high_ind[seq_len(my_number_of_extremes)]])
    names(extreme_PID_df)[4*i-1] <- paste0(this_sig,"_high_PIDs")
    extreme_PID_df[,4*i] <- as.numeric(
      my_exposure_vector[extreme_high_ind[seq_len(my_number_of_extremes)]])
    names(extreme_PID_df)[4*i] <- paste0(this_sig,"_high_exp")
  }
  out_df <- data.frame(t(extreme_PID_df))
  names(out_df) <- paste0("PID_",seq_len(dim(out_df)[2]))
  return(out_df)
}


#' Test if mutated PIDs are enriched in signatures
#' 
#' For all signatures found in a project, this function tests whether PIDs
#' having mutations in a specified list of genes of interest have
#' significantly higher exposures.
#' 
#' @param in_gene_list
#'  List with genes of interest
#' @param in_exposure_df
#'  Data frame with the signature exposures
#' @param in_mut_table
#'  Data frame or table of mutations (derived from vcf-format)
#' @param in_gene.field
#'  Name of the column in which the gene names are to be looked up
#' @param in_p_cutoff
#'  Significance threshold
#' 
#' @return A list with entries
#'  \code{pvals},
#'  \code{exposure_df},
#'  \code{number_of_mutated},
#' \itemize{
#'  \item \code{pvals}:
#'    p-values of the t-tests performed on mutated vs. unmutated PIDs
#'  \item \code{exposure_df}:
#'    Transposed input exposures data frame with additional annotations
#'    for mutation status
#'  \item \code{number_of_mutated}:
#'    Number of PIDs carrying a mutation
#' }
#' 
#' @examples
#'  NULL
#'
#' @export
#'
test_gene_list_in_exposures <- function(in_gene_list,
                                        in_exposure_df,
                                        in_mut_table,
                                        in_gene.field="GENE_short",
                                        in_p_cutoff=0.05) {
  gene_index <- which(names(in_mut_table)==in_gene.field)
  number_of_signatures <- dim(in_exposure_df)[1]
  choice_ind <- which(in_mut_table[,gene_index] %in% in_gene_list)
  #choice_table <- in_mut_table[choice_ind,c(1:5,16:18,20,24,26,29)]
  choice_table <- in_mut_table[choice_ind,]
  this_PIDs <- unique(choice_table$PID)
  temp_ind <- match(this_PIDs,names(in_exposure_df))
  temp_ind_ind <- which(temp_ind>0)
  cat("test_gene_list_in_exposures: number of PIDs in overlap ",
      "between exposures, gene list and mutation table: ",
      length(temp_ind_ind),"\n")
  this_exposure_df <- as.data.frame(t(in_exposure_df))
  this_exposure_df$this_status <- "none"
  this_exposure_df$this_col <- "green"
  this_exposure_df$this_status[temp_ind[temp_ind_ind]] <- "mut"
  this_exposure_df$this_col[temp_ind[temp_ind_ind]] <- "red"
  ## now check for differences in signature counts
  significance_counter <- 0
  out_t_test_vector <- c()
  for (i in seq_len(number_of_signatures)) {
    temp_test_dat <- t.test(
      this_exposure_df[which(this_exposure_df$this_status=="none"),i],
      this_exposure_df[which(this_exposure_df$this_status=="mut"),i],
      alternative="less")
    out_t_test_vector[i] <- temp_test_dat$p.value
    if(out_t_test_vector[i] < in_p_cutoff & !is.na(out_t_test_vector[i])) {
      significance_counter <- significance_counter + 1
    }
  }
  cat("test_gene_list_in_exposures: Number of significancies: ",
      significance_counter,"\n")
  return(list(pvals=out_t_test_vector,
              exposure_df=this_exposure_df,
              number_of_mutated=length(temp_ind_ind)))
}


#' Compare one mutational catalogue to reference mutational catalogues
#' 
#' Compare one mutational catalogue (e.g. of one index patient) to a list of
#' reference mutational catalogues (e.g. from the initial Alexandrov
#' puplication) by cosine similarities
#' 
#' @param in_index_df
#'  Data frame containing the mutational catalogue of interest
#' @param in_comparison_list
#'  List of data frames (ideally named) containing the reference mutational
#'  catalogues
#' 
#' @return A similarity dataframe
#' 
#' @examples
#'  NULL
#'
#' @export
#'
compare_to_catalogues <- function(in_index_df,in_comparison_list){
  sim_vec_list <- lapply(in_comparison_list,FUN=function(x){
    temp_list <- compare_sets(in_index_df,x)
    return(c(as.numeric(1-as.matrix(temp_list$distance))))
  })
  sim_df_list <- list()
  for(i in seq_along(sim_vec_list)) {
    my_entity <- names(sim_vec_list)[i]
    temp_df <- data.frame(entity=my_entity,sim=sim_vec_list[[i]])
    sim_df_list[[i]] <- temp_df
  }
  sim_df <- do.call(rbind,sim_df_list)
  return(sim_df)
}


#' Extract statistical measures for entity comparison
#' 
#' Compare one mutational catalogue (e.g. of one index patient) to a list of
#' reference mutational catalogues (e.g. from the initial Alexandrov
#' puplication) by cosine similarities
#' 
#' @param in_sim_df
#'  A similarity data frame as extracted by \code{\link{compare_to_catalogues}}
#' 
#' @return A dataframe containing statistical measures, prepared for bar plot
#' 
#' @examples
#'  NULL
#'
#' @export
#'
compute_comparison_stat_df <- function(in_sim_df){
  mean_df <- aggregate(sim~entity,data=in_sim_df,FUN=mean)
  median_df <- aggregate(sim~entity,data=in_sim_df,FUN=median)
  sd_df <- aggregate(sim~entity,data=in_sim_df,FUN=sd)
  stat_df <- mean_df
  stat_df$median_sim <- median_df$sim
  stat_df$sd_sim <- sd_df$sim
  names(stat_df)[2] <- "mean_sim"
  return(stat_df)
}


#' Aggregate exposures by category
#' 
#' If a valid category (i.e. it matches to a category specified in 
#' in_sig_ind_df) is supplied, then the exposures are aggregated over this
#' category.
#' 
#' @param in_exposures_df
#'  Input data frame of exposures.
#' @param in_sig_ind_df
#'  Input data frame of meta information on the signatures. Has to match the 
#'  signatures in \code{in_exposures_df}
#' @param in_category
#'  Category to be aggregated over
#' 
#' @return A list with entries:
#' \itemize{
#'  \item \code{exposures}:
#'    The exposures \code{H}, a numeric data frame with 
#'    \code{l} rows and \code{m} columns, \code{l} being
#'    the number of aggregated signatures and \code{m} being the number
#'    of samples
#'  \item \code{norm_exposures}:
#'    The normalized exposures \code{H}, a numeric data frame with 
#'    \code{l} rows and \code{m} columns, \code{l} being
#'    the number of aggregated signatures and \code{m} being the number
#'    of samples
#'  \item \code{out_sig_ind_df}:
#'    Data frame of the type \code{signature_indices_df}, i.e. indicating name,
#'    function and meta-information of the aggregated signatures..
#' }
#' 
#' @seealso \code{\link{LCD_complex_cutoff}}
#' 
#' @examples
#'  NULL
#'
#' @export
#'
aggregate_exposures_by_category <- function(
  in_exposures_df,in_sig_ind_df,in_category){
  if(any(in_category==names(in_sig_ind_df))) {
    ## process exposures
    temp_exposures <- cbind(
      in_exposures_df,in_sig_ind_df[,in_category,drop=FALSE])
    names(temp_exposures) <- gsub(in_category,"cat",names(temp_exposures))
    temp_aggregate_exposures <- aggregate(.~cat,data=temp_exposures,FUN=sum)
    rownames(temp_aggregate_exposures) <- temp_aggregate_exposures$cat
    temp_aggregate_exposures$cat <- NULL
    ## normalize exposures
    norm_aggregate_exposures <- 
      normalize_df_per_dim(temp_aggregate_exposures, 2)
    ## process sig_ind_df
    names(in_sig_ind_df) <- gsub(in_category,"cat",names(in_sig_ind_df))
    out_sig_ind_df <- aggregate(.~cat,data=in_sig_ind_df,FUN=head,1)
    out_sig_ind_df$sig <- out_sig_ind_df$cat
    out_sig_ind_df$cat <- NULL
    out_sig_ind_df$colour <- as.character(out_sig_ind_df$colour)
    return(list(exposures=temp_aggregate_exposures,
                norm_exposures=norm_aggregate_exposures,
                out_sig_ind_df=out_sig_ind_df))
  } else {
    return(NULL)
  }
}
