# Copyright © 2014-2016  The YAPSA package contributors
# This file is part of the YAPSA package. The YAPSA package is licenced under
# GPL-3


#' Data for mutational signatures
#' 
#' The numerical data of the mutational signatures published initially by 
#' Alexandrov et al. (Nature 2013) is stored in data frames with endings 
#' \code{_sig_df}, the associated meta-information is stored in data frames 
#' with endings \code{_sigInd_df}. There are several instances of 
#' \code{_sig_df} and \code{_sigInd_df}, corresponding to results and data 
#' obtained at different times and with different raw data. There always is a 
#' one-to-one correspondence between a \code{_sig_df} and a \code{_sigInd_df}.
#' The data frames of type \code{_sig_df} have as many rows as there are 
#' features, i.e. 96 if analyzing mutational signatures of SNVs in a triplet 
#' context, and as many columns as there are signatures.
#' Data frames of type \code{_sigInd_df} have as many rows as there are 
#' signatures in the corresponding \code{_sig_df} and several columns: 
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
#' @docType data
#' @name sigs
#' @usage data(sigs)
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' @references Alexandrov et al. (Nature 2013)
#' 
NULL


#' Data for initial sigs, including artifacts
#' 
#' \code{AlexInitialArtif_sig_df}: Data frame of the signatures published 
#' initially by Alexandrov et al. 
#' (Nature 2013). There are 27 signatures which constitute the columns, 22 of 
#' which were validated by an orhtogonal sequencing technology. These 22 are in
#' the first 22 columns of the data frame. The column names are \emph{A} pasted
#' to the number of the signature, e.g. \emph{A5}. The nonvalidated signatures
#' have an additional letter in their naming convention: either 
#' \emph{AR1} - \emph{AR3} or \emph{AU1} - \emph{AU2}. The rownames are the 
#' features, i.e. an encoding of the nucleotide exchanges in their 
#' trinucleotide context, e.g. \emph{C>A ACA}. In total there are 96 different 
#' features and therefore 96 rows when dealing with a trinucleotide context. 
#' 
#' @source \code{AlexInitial}: \url{ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/signatures.txt}
#' @name AlexInitialArtif_sig_df
#' @rdname sigs
#' 
NULL


#' Meta-info for initial sigs, including artifacts
#' 
#' \code{AlexInitialArtif_sigInd_df}: Meta-information for 
#' \code{AlexInitialArtif_sig_df} 
#'  
#' @name AlexInitialArtif_sigInd_df
#' @rdname sigs
#' 
NULL


#' Data for initial sigs, only validated
#' 
#' \code{AlexInitialValid_sig_df}: Data frame of only the validated signatures 
#' published initially by Alexandrov et al. (Nature 2013), corresponding to the
#' first 22 columns of \code{AlexInitialArtif_sig_df}
#' 
#' @name AlexInitialValid_sig_df
#' @rdname sigs
#' 
NULL


#' Meta-info for initial sigs, only validated
#' 
#' \code{AlexInitialValid_sigInd_df}: Meta-information for 
#' \code{AlexInitialValid_sig_df} 
#' 
#' @name AlexInitialValid_sigInd_df
#' @rdname sigs
#' 
NULL


#' Data for Cosmic sigs, only validated
#' 
#' \code{AlexCosmicValid_sig_df}: Data frame of the updated signatures list 
#' maintained by Ludmil Alexandrov at
#' \url{http://cancer.sanger.ac.uk/cosmic/signatures}. The column names are 
#' \emph{AC} pasted to the number of the signature, e.g. \emph{AC5}. The naming
#' convention for the rows is as described for 
#' \code{\link{AlexInitialArtif_sig_df}}.
#' 
#' @source \code{AlexCosmic}: \url{http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt}
#' @name AlexCosmicValid_sig_df
#' @rdname sigs
#' 
NULL


#' Meta-info for Cosmic sigs, only validated
#' 
#' \code{AlexCosmicValid_sigInd_df}: Meta-information for 
#' \code{AlexCosmicValid_sig_df} 
#' 
#' @name AlexCosmicValid_sigInd_df
#' @rdname sigs
#' 
NULL


#' Data for Cosmic sigs, including artifacts
#' 
#' \code{AlexCosmicArtif_sig_df}: Data frame of the updated signatures list 
#' maintained by Ludmil Alexandrov at
#' \url{http://cancer.sanger.ac.uk/cosmic/signatures} and complemented by the 
#' artifact signatures from the initial publication, i.e. the last 5 columns of
#' \code{\link{AlexInitialArtif_sig_df}}. The column names are \emph{AC} pasted
#' to the number of the signature, e.g. \emph{AC5}. The naming convention for 
#' the rows is as described for \code{\link{AlexInitialArtif_sig_df}}.
#' 
#' @name AlexCosmicArtif_sig_df
#' @rdname sigs
#' 
NULL


#' Meta-info for Cosmic sigs, including artifacts
#' 
#' \code{AlexCosmicArtif_sigInd_df}: Meta-information for 
#' \code{AlexCosmicArtif_sig_df} 
#' 
#' @name AlexCosmicArtif_sigInd_df
#' @rdname sigs
#' 
NULL


#' Test and example data
#' 
#' Data structures used in examples, tests and the vignette of the YAPSA 
#' package.
#' 
#' @docType data
#' @name exampleYAPSA
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/23945592}
#' 
NULL



#' Subgroup information for some samples in the vignette
#' 
#' \code{lymphoma_PID_df}: A data frame carrying subgroup information for a 
#' subcohort of samples used in the vignette. Data in the vignette is 
#' downloaded from 
#' \url{ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Lymphoma B-cell/Lymphoma B-cell_clean_somatic_mutations_for_signature_analysis.txt}.
#' In the file available under that link somatic point mutation calls from 
#' several samples are listed in a vcf-like format. One column encodes the 
#' sample the variant was found in. In the vignette we want to restrict the 
#' analysis to only a fraction of these involved samples. The data frame 
#' \code{lymphoma_PID_df} carries the sample identifiers (PID) as rownames and 
#' the attributed subgroup in a column called \code{subgroup}.
#' 
#' @name lymphoma_PID_df
#' @usage data(lymphoma_PID)
#' @rdname exampleYAPSA
#' 
NULL


#' Test data for complex functions
#' 
#' \code{lymphoma_test_df}: A data frame carrying point mutation calls. It 
#' represents a subset of the data stored in
#' \url{ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Lymphoma B-cell/Lymphoma B-cell_clean_somatic_mutations_for_signature_analysis.txt}.
#' In the file available under that link somatic point mutation calls from 
#' several samples are listed in a vcf-like format. One column encodes the 
#' sample the variant was found in. The data frame \code{lymphoma_test_df} has 
#' only the variants occuring in the sample identifiers (PIDs) 4112512, 4194218
#' and 4121361.
#' 
#' @name lymphoma_test_df
#' @usage data(lymphoma_test)
#' @rdname exampleYAPSA
#' 
#' @examples
#' data(lymphoma_test)
#' head(lymphoma_test_df)
#' dim(lymphoma_test_df)
#' table(lymphoma_test_df$PID)
#' 
NULL


#' Example data for the vignette
#' 
#' \code{lymphoma_Nature2013_raw_df}: A data frame carrying point mutation 
#' calls. It represents a subset of the data stored in
#' \url{ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Lymphoma B-cell/Lymphoma B-cell_clean_somatic_mutations_for_signature_analysis.txt}.
#' In the file available under that link somatic point mutation calls from 
#' several samples are listed in a vcf-like format. One column encodes the 
#' sample the variant was found in.
#' 
#' @name lymphoma_Nature2013_raw_df
#' @usage data(lymphoma_Nature2013_raw)
#' @rdname exampleYAPSA
#' 
#' @examples
#' data(lymphoma_Nature2013_raw)
#' head(lymphoma_Nature2013_raw_df)
#' dim(lymphoma_Nature2013_raw_df)
#' 
NULL


#' Test exposures for plot functions
#' 
#' \code{lymphoma_Nature2013_COSMIC_cutoff_exposures_df}: Data frame with 
#' exposures for testing the plot functions. Data taken from
#' \url{ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Lymphoma B-cell/Lymphoma B-cell_clean_somatic_mutations_for_signature_analysis.txt}.
#' 
#' @name lymphoma_Nature2013_COSMIC_cutoff_exposures_df
#' @usage data(lymphoma_cohort_LCD_results)
#' @rdname exampleYAPSA
#' 
NULL


#' Test normalized exposures for plot functions
#' 
#' \code{rel_lymphoma_Nature2013_COSMIC_cutoff_exposures_df}: Data frame with 
#' normalized or relative exposures for testing the plot functions. Data taken 
#' from
#' \url{ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Lymphoma B-cell/Lymphoma B-cell_clean_somatic_mutations_for_signature_analysis.txt}.
#' 
#' @name rel_lymphoma_Nature2013_COSMIC_cutoff_exposures_df
#' @usage data(lymphoma_cohort_LCD_results)
#' @rdname exampleYAPSA
#' 
NULL


#' Subgroup information for test data for plot functions
#' 
#' \code{COSMIC_subgroups_df}: Subgroup information for the data stored in
#' \code{\link{lymphoma_Nature2013_COSMIC_cutoff_exposures_df}} and 
#' \code{\link{rel_lymphoma_Nature2013_COSMIC_cutoff_exposures_df}}.
#' 
#' @name COSMIC_subgroups_df
#' @usage data(lymphoma_cohort_LCD_results)
#' @rdname exampleYAPSA
#' 
NULL


#' Sigs info (initial, including artifacts) for test data for plot functions
#' 
#' \code{chosen_AlexInitialArtif_sigInd_df}: Signature information for the data
#' stored in
#' \code{\link{lymphoma_Nature2013_COSMIC_cutoff_exposures_df}} and 
#' \code{\link{rel_lymphoma_Nature2013_COSMIC_cutoff_exposures_df}}.
#' 
#' @name chosen_AlexInitialArtif_sigInd_df
#' @usage data(lymphoma_cohort_LCD_results)
#' @rdname exampleYAPSA
#' 
NULL


#' Sigs info (Cosmic, only validated) for test data for plot functions
#' 
#' \code{chosen_signatures_indices_df}: Signature information for the data
#' stored in
#' \code{\link{lymphoma_Nature2013_COSMIC_cutoff_exposures_df}} and 
#' \code{\link{rel_lymphoma_Nature2013_COSMIC_cutoff_exposures_df}}.
#' 
#' @name chosen_signatures_indices_df
#' @usage data(lymphoma_cohort_LCD_results)
#' @rdname exampleYAPSA
#' 
NULL


#' Cutoffs for a supervised analysis of mutational signatures.
#' 
#' Series of data frames with signature-specific cutoffs. All values represent
#' optimal cutoffs. The optimal cutoffs were determined for different choices 
#' of parameters in the cost function of the optimization. The row index is 
#' equivalent to the ratio between costs for false negative attribution and 
#' false positive attribution. The columns correspond to the different 
#' signatures. To be used with \code{\link{LCD_complex_cutoff}}. 
#' 
#' @docType data
#' @name cutoffs
#' @usage data(cutoffs)
#' @author Daniel Huebschmann \email{huebschmann.daniel@@googlemail.com}
#' 
NULL


#' Opt. cutoffs, rel exposures for the COSMIC sigs, only validated
#' 
#' \code{cutoffCosmicValid_rel_df}: Optimal cutoffs for 
#' \code{\link{AlexCosmicValid_sig_df}}, i.e. COSMIC signatures, only 
#' validated, trained on relative exposures.
#' 
#' @name cutoffCosmicValid_rel_df
#' @rdname cutoffs
#' 
NULL


#' Opt. cutoffs, rel exposures for the COSMIC sigs, including artifacts
#' 
#' \code{cutoffCosmicArtif_rel_df}: Optimal cutoffs for 
#' \code{\link{AlexCosmicArtif_sig_df}}, i.e. COSMIC signatures, including 
#' artifact signatures, trained on relative exposures.
#' 
#' @name cutoffCosmicArtif_rel_df
#' @rdname cutoffs
#' 
NULL


#' Opt. cutoffs, abs exposures for the COSMIC sigs, only validated
#' 
#' \code{cutoffCosmicValid_abs_df}: Optimal cutoffs for 
#' \code{\link{AlexCosmicValid_sig_df}}, i.e. COSMIC signatures, only 
#' validated, trained on absolute exposures.
#' 
#' @name cutoffCosmicValid_abs_df
#' @rdname cutoffs
#' 
NULL


#' Opt. cutoffs, abs exposures for the COSMIC sigs, including artifacts
#' 
#' \code{cutoffCosmicArtif_abs_df}: Optimal cutoffs for 
#' \code{\link{AlexCosmicArtif_sig_df}}, i.e. COSMIC signatures, including 
#' artifact signatures, trained on absolute exposures.
#' 
#' @name cutoffCosmicArtif_abs_df
#' @rdname cutoffs
#' 
NULL


#' Opt. cutoffs, rel exposures for the initial sigs, only validated
#' 
#' \code{cutoffInitialValid_rel_df}: Optimal cutoffs for 
#' \code{\link{AlexInitialValid_sig_df}}, i.e. initially published signatures, 
#' only validated signatures, trained on relative exposures.
#' 
#' @name cutoffInitialValid_rel_df
#' @rdname cutoffs
#' 
NULL


#' Opt. cutoffs, rel exposures for the initial sigs, including artifacts
#' 
#' \code{cutoffInitialArtif_rel_df}: Optimal cutoffs for 
#' \code{\link{AlexInitialArtif_sig_df}}, i.e. initially published signatures, 
#' including artifact signatures, trained on relative exposures.
#' 
#' @name cutoffInitialArtif_rel_df
#' @rdname cutoffs
#' 
NULL


#' Opt. cutoffs, abs exposures for the initial sigs, only validated
#' 
#' \code{cutoffInitialValid_abs_df}: Optimal cutoffs for 
#' \code{\link{AlexInitialValid_sig_df}}, i.e. initially published signatures, 
#' only validated signatures, trained on absolute exposures.
#' 
#' @name cutoffInitialValid_abs_df
#' @rdname cutoffs
#' 
NULL


#' Opt. cutoffs, abs exposures for the initial sigs, including artifacts
#' 
#' \code{cutoffInitialArtif_abs_df}: Optimal cutoffs for 
#' \code{\link{AlexInitialArtif_sig_df}}, i.e. initially published signatures, 
#' including artifact signatures, trained on absolute exposures.
#' 
#' @name cutoffInitialArtif_abs_df
#' @rdname cutoffs
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
#' @return A list of lists of data frames
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
#' @return A named character vector
#' 
NULL
