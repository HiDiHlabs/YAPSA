#' @importFrom GetoptLong qq
#' @importFrom GetoptLong qq.options
#' 
check_perl = function(module = NULL, inc = NULL, perl_bin = "perl") {
  
  op = qq.options("code.pattern")
  qq.options("code.pattern" = "@\\{CODE\\}")
  on.exit(qq.options("code.pattern" = op))
  
  if(is.null(module)) {
    cmd = qq("\"@{perl_bin}\" -v")
  } else if(!is.null(module) && is.null(inc)) {
    cmd = qq("\"@{perl_bin}\" -M@{module} -e \"use @{module}\"")
  } else if(!is.null(module) && !is.null(inc)) {
    cmd = qq("\"@{perl_bin}\" \"-I@{inc}\" -M@{module} -e \"use @{module}\"")
  }
  
  OS = Sys.info()["sysname"]
  if(OS == "Windows") {
    res = system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE, show.output.on.console = FALSE)
  } else {
    res = system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  }
  
  return(ifelse(res, FALSE, TRUE))
}

#' @importFrom GetoptLong qq
#' @importFrom GetoptLong qq.options
#' 
check_bedtools = function(module = NULL, inc = NULL, bedtools_bin = "bedtools") {
  
  op = qq.options("code.pattern")
  qq.options("code.pattern" = "@\\{CODE\\}")
  on.exit(qq.options("code.pattern" = op))
  
  if(is.null(module)) {
    cmd = qq("\"@{bedtools_bin}\" --version")
  }
  
  OS = Sys.info()["sysname"]
  if(OS == "Windows") {
    res = system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE, show.output.on.console = FALSE)
  } else {
    res = system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  }
  
  return(ifelse(res, FALSE, TRUE))
}


## example:
## in_fasta <- "/ibios/co02/reference/Reference_1KG/hs37d5.fa"
## in_word_length <- 3
## project_folder <- "/icgc/dkfzlsdf/analysis/hipo/hipo_028/somaticSignatures/nucleotide_distrib"
## out_csv <- "kmer_frequencies_in_ref.csv"
## 
run_kmer_frequency_pl <- function(in_fasta,in_word_length,project_folder,out_csv) {
  if(!check_perl()) {
    cat("YAPSA:::run_kmer_frequency_pl::Error:unable to find path to perl")
    return(1)
  }
  package_path <- system.file(package='YAPSA')
  perl_command <- paste0("perl")
  script_command <- file.path(package_path,"foreign","kmer_frequencies.pl")
  script_options_command <- paste0("-r ",in_fasta, " -w ",in_word_length)
  output_command <- paste0("> ",file.path(project_folder,out_csv))
  this_command <- paste(perl_command,script_command,script_options_command,output_command,collapse=" ")
  system(this_command,intern=TRUE)
  return(0)
}

## example:
## in_ref_genome <- "/ibios/co02/reference/Reference_1KG/hs37d5.fa"
## in_target_capture_bed <- "/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/targetRegions/Agilent5withUTRs_chr.bed.gz"
## project_folder <- "/icgc/dkfzlsdf/analysis/hipo/hipo_028/somaticSignatures/nucleotide_distrib"
## out_target_capture_fasta <- "hs37d5_Agilent5withUTR_targetCapture.fa"
## 
run_bedtools_getfasta <- function(in_ref_genome,in_target_capture_bed,
                                  project_folder,out_target_capture_fasta) {
  if(!check_bedtools()) {
    cat("YAPSA:::run_bedtools_getfasta::Error:unable to find path to bedtools")
    return(1)
  }
  bedtools_command <- "bedtools getfasta"
  bedtools_options_command <- paste0("-fi ",in_ref_genome, " -bed ",in_target_capture_bed,
                                     " -fo ",file.path(project_folder,out_target_capture_fasta))
  this_command <- paste(bedtools_command,bedtools_options_command,collapse=" ")
  system(this_command)
  return(0)
}

## exmple:
## ${PATH_TO_REFERENCE}=/ibios/co02/reference/Reference_1KG/hs37d5.fa
## ${PATH_TO_PROJECT}=/icgc/dkfzlsdf/analysis/hipo/hipo_028/somaticSignatures/nucleotide_distrib
## ${PATH_TO_TARGET_CAPTURE_BED}=/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/targetRegions/Agilent5withUTRs_chr.bed.gz
##
## step 1: kmers in reference genome:
## 1.1. perl ${PATH_TO_SCRIPT}/kmer_frequencies.pl -r ${PATH_TO_REFERENCE} -w 3 > ${PATH_TO_PROJECT}/kmer_frequencies_in_ref.csv
## step 2: kmers in target capture:
## 2.1. bedtools getfasta -fi ${PATH_TO_REFERENCE} -bed ${PATH_TO_TARGET_CAPTURE_BED} -fo ${PATH_TO_TARGET_CAPTURE_FASTA}
## 2.2. perl ${PATH_TO_SCRIPT}/kmer_frequencies.pl -r ${PATH_TO_TARGET_CAPTURE_FASTA} -w 3 > ${PATH_TO_PROJECT}/kmer_frequencies_in_target_capture.csv
##
## in_ref_genome_fasta <- "/ibios/co02/reference/Reference_1KG/hs37d5.fa"
## in_target_capture_bed <- "/icgc/ngs_share/assemblies/hg19_GRCh37_1000genomes/targetRegions/Agilent5withUTRs_chr.bed.gz"
## in_word_length <- 3
## project_folder <- "/icgc/dkfzlsdf/analysis/hipo/hipo_028/somaticSignatures/nucleotide_distrib"
##

#' Provide normalized correction factors for kmer content
#' 
#' This function is analogous to \code{\link[SomaticSignatures]{normalizeMotifs}}. If an 
#' analysis of mutational signatures is performed on e.g. Whole Exome Sequencing (WES)
#' data, the signatures and exposures have to be adapted to the potentially different kmer
#' (trinucleotide) content of the target capture. The present function takes as arguments
#' paths to the used reference genome and target capture file. It the extracts the
#' sequence of the target capture by calling \code{bedtools getfasta} on the system
#' command prompt. \code{run_kmer_frequency_normalization} then calls a custom made perl
#' script \code{kmer_frequencies.pl} also included in this package to count the occurences
#' of the tripletts in both the whole reference genome and the created target capture
#' sequence. These counts are used for normalization as in
#' \code{\link[SomaticSignatures]{normalizeMotifs}}. Note that
#' \code{\link[SomaticSignatures]{kmerFrequency}} provides a solution to approximate kmer
#' frequencies by random sampling. As opposed to that approach, the function described
#' here deterministically counts all occurences of the kmers in the respective genome.
#' 
#' @param in_ref_genome_fasta
#'  Path to the reference genome fasta file used.
#' @param in_target_capture_bed
#'  Path to a bed file containing the information on the used target capture. May also be
#'  a compressed bed.
#' @param in_word_length
#'  Integer number defining the length of the features or motifs, e.g. 3 for
#'  tripletts or 5 for pentamers
#' @param project_folder
#'  Path where the created files, especially the fasta file with the sequence of the
#'  target capture and the count matrices, can be stored.
#' @param in_verbose
#'  Verbose if \code{in_verbose=1}
#'
#' @return
#'  A numeric vector with correction factors
#'
#' @examples
#'  NULL
#'
#' @seealso \code{\link[SomaticSignatures]{normalizeMotifs}}
#'
#' @export
#' 
run_kmer_frequency_normalization <- function(in_ref_genome_fasta,in_target_capture_bed,in_word_length,
                                             project_folder,in_verbose=1) {
  target_capture_fasta <- "hs37d5_Agilent5withUTR_targetCapture.fa"
  target_capture_kmer_counts_file <- "kmer_frequencies_in_targetCapture.csv"
  reference_genome_kmer_counts_file <- "kmer_frequencies_in_ref.csv"
  if(!file.exists(file.path(project_folder,reference_genome_kmer_counts_file))) {
    if(in_verbose==1) {cat("\nYAPS:::run_kmer_frequency_normalization::reference_genome_kmer_counts_file doesn't exist yet. Create first.\n");}
    run_kmer_frequency_pl(in_ref_genome_fasta,in_word_length,project_folder,reference_genome_kmer_counts_file)
  } else {
    if(in_verbose==1) {cat("\nYAPS:::run_kmer_frequency_normalization::reference_genome_kmer_counts_file already exists.\n");}
  }
  if(!file.exists(file.path(project_folder,target_capture_kmer_counts_file))) {
    if(in_verbose==1) {cat("\nYAPS:::run_kmer_frequency_normalization::target_capture_kmer_counts_file doesn't exist yet. Create first.\n");}
    if(!file.exists(file.path(project_folder,target_capture_fasta))) {
      if(in_verbose==1) {cat("\nYAPS:::run_kmer_frequency_normalization::target_capture_fasta doesn't exist yet. Create first.\n");}
      run_bedtools_getfasta(in_ref_genome_fasta,in_target_capture_bed,
                            project_folder,target_capture_fasta)
    } else {
      if(in_verbose==1) {cat("\nYAPS:::run_kmer_frequency_normalization::target_capture_fasta already exists.\n");}    
    }
    run_kmer_frequency_pl(file.path(project_folder,target_capture_fasta),in_word_length,project_folder,target_capture_kmer_counts_file)
  } else {
    if(in_verbose==1) {cat("\nYAPS:::run_kmer_frequency_normalization::target_capture_kmer_counts_file already exists.\n");}    
  }
  if(in_verbose==1) {cat("\nYAPS:::run_kmer_frequency_normalization::Now read counts into data frames...\n");}
  reference_genome_kmer_frequency_df <- read.csv(file.path(project_folder,reference_genome_kmer_counts_file),header=FALSE,sep="\t")
  names(reference_genome_kmer_frequency_df) <- c("triplet","count")
  reference_genome_total_count <- sum(as.numeric(reference_genome_kmer_frequency_df$count))
  reference_genome_kmer_frequency_df$rel_counts <- as.numeric(reference_genome_kmer_frequency_df$count) / reference_genome_total_count
  target_capture_kmer_frequency_df <- read.csv(file.path(project_folder,target_capture_kmer_counts_file),header=FALSE,sep="\t")
  names(target_capture_kmer_frequency_df) <- c("triplet","count")
  target_capture_total_count <- sum(as.numeric(target_capture_kmer_frequency_df$count))
  target_capture_kmer_frequency_df$rel_counts <- as.numeric(target_capture_kmer_frequency_df$count) / target_capture_total_count
  norms <- reference_genome_kmer_frequency_df$rel_counts / target_capture_kmer_frequency_df$rel_counts
  names(norms) <- reference_genome_kmer_frequency_df$triplet
  return(norms)
}


#' Provide comprehensive correction factors for kmer content
#' 
#' This function is analogous to \code{\link[SomaticSignatures]{normalizeMotifs}}. If an 
#' analysis of mutational signatures is performed on e.g. Whole Exome Sequencing (WES)
#' data, the signatures and exposures have to be adapted to the potentially different kmer
#' (trinucleotide) content of the target capture. The present function takes as arguments
#' paths to the used reference genome and target capture file. It the extracts the
#' sequence of the target capture by calling \code{bedtools getfasta} on the system
#' command prompt. \code{run_kmer_frequency_normalization} then calls a custom made perl
#' script \code{kmer_frequencies.pl} also included in this package to count the occurences
#' of the tripletts in both the whole reference genome and the created target capture
#' sequence. These counts are used for normalization as in
#' \code{\link[SomaticSignatures]{normalizeMotifs}}. Note that
#' \code{\link[SomaticSignatures]{kmerFrequency}} provides a solution to approximate kmer
#' frequencies by random sampling. As opposed to that approach, the function described
#' here deterministically counts all occurences of the kmers in the respective genome.
#' 
#' @param in_ref_genome_fasta
#'  Path to the reference genome fasta file used.
#' @param in_target_capture_bed
#'  Path to a bed file containing the information on the used target capture. May also be
#'  a compressed bed.
#' @param in_word_length
#'  Integer number defining the length of the features or motifs, e.g. 3 for
#'  tripletts or 5 for pentamers
#' @param project_folder
#'  Path where the created files, especially the fasta file with the sequence of the
#'  target capture and the count matrices, can be stored.
#' @param target_capture_fasta
#'  Name of the fasta file of the target capture to be created if not yet existent.
#' @param in_verbose
#'  Verbose if \code{in_verbose=1}
#'
#' @return
#'  A list with 2 entries:
#'  \itemize{
#'    \item \code{rel_cor}:
#'      The correction factors after normalization as in \code{\link{run_kmer_frequency_normalization}}
#'    \item \code{abs_cor}:
#'      The correction factors without normalization.
#'  }
#'
#' @examples
#'  NULL
#'
#' @seealso \code{\link[SomaticSignatures]{normalizeMotifs}}
#'
#' @export
#' 
run_kmer_frequency_correction <- function(in_ref_genome_fasta,in_target_capture_bed,in_word_length,
                                          project_folder,target_capture_fasta="targetCapture.fa",
                                          in_verbose=1) {
  #target_capture_fasta <- "hs37d5_Agilent5withUTR_targetCapture.fa"
  target_capture_kmer_counts_file <- "kmer_frequencies_in_targetCapture.csv"
  reference_genome_kmer_counts_file <- "kmer_frequencies_in_ref.csv"
  dir.create(project_folder,recursive=TRUE)
  if(!file.exists(file.path(project_folder,reference_genome_kmer_counts_file))) {
    if(in_verbose==1) {cat("\nYAPS:::run_kmer_frequency_normalization::reference_genome_kmer_counts_file doesn't exist yet. Create first.\n");}
    run_kmer_frequency_pl(in_ref_genome_fasta,in_word_length,project_folder,reference_genome_kmer_counts_file)
  } else {
    if(in_verbose==1) {cat("\nYAPS:::run_kmer_frequency_normalization::reference_genome_kmer_counts_file already exists.\n");}
  }
  if(!file.exists(file.path(project_folder,target_capture_kmer_counts_file))) {
    if(in_verbose==1) {cat("\nYAPS:::run_kmer_frequency_normalization::target_capture_kmer_counts_file doesn't exist yet. Create first.\n");}
    if(!file.exists(file.path(project_folder,target_capture_fasta))) {
      if(in_verbose==1) {cat("\nYAPS:::run_kmer_frequency_normalization::target_capture_fasta doesn't exist yet. Create first.\n");}
      run_bedtools_getfasta(in_ref_genome_fasta,in_target_capture_bed,
                            project_folder,target_capture_fasta)
    } else {
      if(in_verbose==1) {cat("\nYAPS:::run_kmer_frequency_normalization::target_capture_fasta already exists.\n");}    
    }
    run_kmer_frequency_pl(file.path(project_folder,target_capture_fasta),in_word_length,project_folder,target_capture_kmer_counts_file)
  } else {
    if(in_verbose==1) {cat("\nYAPS:::run_kmer_frequency_normalization::target_capture_kmer_counts_file already exists.\n");}    
  }
  if(in_verbose==1) {cat("\nYAPS:::run_kmer_frequency_normalization::Now read counts into data frames...\n");}
  reference_genome_kmer_frequency_df <- read.csv(file.path(project_folder,reference_genome_kmer_counts_file),header=FALSE,sep="\t")
  names(reference_genome_kmer_frequency_df) <- c("triplet","count")
  reference_genome_total_count <- sum(as.numeric(reference_genome_kmer_frequency_df$count))
  reference_genome_kmer_frequency_df$rel_counts <- as.numeric(reference_genome_kmer_frequency_df$count) / reference_genome_total_count
  target_capture_kmer_frequency_df <- read.csv(file.path(project_folder,target_capture_kmer_counts_file),header=FALSE,sep="\t")
  names(target_capture_kmer_frequency_df) <- c("triplet","count")
  target_capture_total_count <- sum(as.numeric(target_capture_kmer_frequency_df$count))
  target_capture_kmer_frequency_df$rel_counts <- as.numeric(target_capture_kmer_frequency_df$count) / target_capture_total_count
  rel_norms <- reference_genome_kmer_frequency_df$rel_counts / target_capture_kmer_frequency_df$rel_counts
  names(rel_norms) <- reference_genome_kmer_frequency_df$triplet
  abs_norms <- reference_genome_kmer_frequency_df$count / target_capture_kmer_frequency_df$count
  names(abs_norms) <- reference_genome_kmer_frequency_df$triplet
  #out_df <- data.frame(rel_cor=rel_norms,abs_cor=abs_norms)
  #rownames(out_df) <- reference_genome_kmer_frequency_df$triplet
  #return(out_df)
  out_list <- list(rel_cor=rel_norms,abs_cor=abs_norms)
  return(out_list)
}


#' Wrapper function to annotate addition information
#' 
#' Wrapper function to the perl script annotate_vcf.pl which annotates
#' data of a track stored in file_B (may be different formats) to called 
#' variants stored in a vcf-like file_A.
#' 
#' @param in_data_file
#'  Path to the input vcf-like file to be annotated
#' @param in_anno_track_file
#'  Path to the input file containing the annotation track
#' @param in_new_column_name
#'  String indicating the name of the column to be created for annotation.
#' @param out_file
#'  Path where the created files can be stored.
#' @param in_data_file_type
#'  \code{custom} for vcf-like
#' @param in_anno_track_file_type
#'  Type of the file \code{in_anno_track_file} containing the annotation
#'  track.
#' @param in_data_CHROM.field
#'  String indicating which column of \code{in_data_file} contains the
#'  chromosome information.
#' @param in_data_POS.field
#'  String indicating which column of \code{in_data_file} contains the
#'  position information.
#' @param in_data_END.field
#'  String indicating which column of \code{in_data_file} contains the
#'  end information if regions are considered.
#'  
#' @return Return zero if no problems occur.
#' 
#' @examples
#'  NULL
#'  
#' @export
#' 
run_annotate_vcf_pl <- function(in_data_file,
                                in_anno_track_file,
                                in_new_column_name,
                                out_file,
                                in_data_file_type="custom",
                                in_anno_track_file_type="bed",
                                in_data_CHROM.field="CHROM",
                                in_data_POS.field="POS",
                                in_data_END.field="POS") {
  if(!check_perl()) {
    cat("YAPSA:::run_annotate_vcf_pl::Error:unable to find path to perl")
    return(1)
  }
  package_path <- system.file(package='YAPSA')
  perl_command <- paste0("perl")
  script_command <- file.path(package_path,"foreign","annotate_vcf.pl")
  script_options_command <- paste0("-a ",in_data_file,
                                   " --aFileType=",in_data_file_type,
                                   " --aChromColumn ",in_data_CHROM.field,
                                   " --aPosColumn ",in_data_POS.field,
                                   " --aEndColumn ",in_data_END.field,
                                   " -b ",in_anno_track_file,
                                   " --bFileType=",in_anno_track_file_type,
                                   " --columnName ",in_new_column_name)
  output_command <- paste0("> ",file.path(out_file))
  this_command <- paste(perl_command,script_command,script_options_command,output_command,collapse=" ")
  system(this_command,intern=TRUE)
  return(0)
}


#' Build a gene list for a given pathway name
#' 
#' @param in_string
#'  Name or description of the pathway
#' @param in_organism
#'  Name of the taxon to be searched in
#' 
#' @return A character vector of gene names
#' 
#' @examples
#'  NULL
#'  \dontrun{
#'    species <- "hsa"
#'    gene_lists_meta_df <- data.frame(name=c("BER","NHEJ","MMR"),
#'                                     explanation=c("base excision repair",
#'                                                   "non homologous end joining",
#'                                                   "mismatch repair"))
#'    number_of_pathways <- dim(gene_lists_meta_df)[1]
#'    gene_lists_list <- list()
#'    for (i in seq_len(number_of_pathways)) {
#'      temp_list <- build_gene_list_for_pathway(gene_lists_meta_df$explanation[i],species)
#'      gene_lists_list <- c(gene_lists_list,list(temp_list))
#'    }
#'    gene_lists_list
#'  }
#' 
#' @seealso \code{\link[KEGGREST]{keggLink}}
#' @seealso \code{\link[KEGGREST]{keggFind}}
#' @seealso \code{\link{extract_names_from_gene_list}}
#' 
#' @importFrom KEGGREST keggFind
#' @importFrom KEGGREST keggLink
#' @export
#' 
build_gene_list_for_pathway <- function(in_string,in_organism) {
  my_pathway <- keggFind("pathway",in_string)
  my_name <- gsub("map",in_organism,names(my_pathway))
  KEGG_gene_list <- keggLink(in_organism,my_name)
  out_gene_list <- c()
  for (i in seq_len(length(KEGG_gene_list))) {
    out_gene_list[i] <- extract_names_from_gene_list(KEGG_gene_list,i)
  }  
  return(out_gene_list)
}


#' Return gene names from gene lists
#' 
#' @param in_KEGG_gene_list
#'  Gene list to extract names from
#' @param l
#'  Index of the gene to be extracted
#' 
#' @return The gene name.
#' 
#' @examples
#'  NULL
#' 
#' @seealso \code{\link[KEGGREST]{keggGet}}
#' @seealso \code{\link{build_gene_list_for_pathway}}
#' 
#' @importFrom KEGGREST keggGet
#' @export
#' 
extract_names_from_gene_list <- function(in_KEGG_gene_list,l) {
  temp_gene <- keggGet(in_KEGG_gene_list[l])[[1]]
  temp_gene_name <- strsplit(temp_gene$NAME,", ")[[1]][1]  
  return(temp_gene_name)
}
