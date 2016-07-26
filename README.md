# YAPSA
Daniel Huebschmann  
26/08/2015  


# Citation

The package is currently being submitted to 
[Bioconductor](https://www.bioconductor.org/). Please use it once it is 
accepted there and a suitable citation is provided.


# General design

## Introduction

The concept of mutational signatures was introduced in a series of papers by 
Ludmil Alexandrov, Serena Nik-Zainal, Michael Stratton and others (for precise 
citations please refer to the vignette of the package). The general 
approach is as follows:

1. The SNVs are categorized by their nucleotide exchange. In total there are 
`4 x 3 = 12` different nucleotide exchanges, but if summing over reverse 
complements only `12 / 2 = 6` different categories are left. For every SNV 
detected, the motif context around the position of the SNV is extracted. This 
may be a trinucleotide context if taking one base upstream and one base 
downstream of the position of the SNV, but larger motifs may be taken as well 
(e.g. pentamers). Taking into account the motif context increases combinatorial 
complexity: in the case of the trinucleotide context, there are 
`4 x 6 x 4 = 96` different variant categories. These categories are 
called **features** in the following text. The number of features will be 
called `n`.
2. A cohort consists of different samples with the number of samples denoted by 
`m`. For each sample we can count the occurences of each feature, yielding an 
`n`-dimensional vector (`n` being the number of features) per sample. For a 
cohort, we thus get an `n x m` -dimensional matrix, called the 
**mutational catalogue** `V`. It can be understood as a summary indicating 
which sample has how many variants of which category, but omitting the 
information of the genomic coordinates of the variants.
3. The mutational catalogue `V` is quite big and still carries a lot of 
complexity. For many analyses a  reduction of complexity is desirable. One way 
to achieve such a complexity reduction is a matrix decomposition: we would like 
to find two smaller matrices `W` and `H` which if multiplied would span a high 
fraction of the complexity of the big matrix `V` (the mutational catalogue). 
Remember that `V` is an `n x m` -dimensional matrix, `n` being the number 
of features and `m` being the number of samples. `W` in this setting is an 
`n x l` -dimensional matrix and `H` is an `l x m` -dimensional 
matrix. The columns of 
`W` are called the **mutational signatures** and the columns of `H` are called 
**exposures**. `l` denotes the number of mutational signatures. Hence the 
signatures are `n`-dimensional vectors (with `n` being the number of features), 
while the exposures are `l`-dimensional vectors (`l` being the number of 
signatures). Note that as we are dealing with count data, we would like to have 
only positive entries in `W` and `H`. A mathematical method which is able to do 
such a decomposition is the **NMF** (**nonnegative matrix factorization**). It 
basically solves the problem as illustrated in the following figure:

![NMF, Image taken from <https://en.wikipedia.org/wiki/Non-negative_matrix_factorization>](https://upload.wikimedia.org/wikipedia/commons/f/f9/NMF.png)

Of course the whole concept only leads to a reduction in complexity if `l < n`, 
i.e. if the number of signatures is smaller than the number of features, as 
indicated in the above figure. Note that the NMF itself solves the above 
problem for a given number of signatures `l`. Addinional criteria exist to 
evaluate the true number of signatures.

## The YAPSA package

In a context where mutational signatures `W` are already known (because they 
were decribed and published or they are available in a 
database as under <http://cancer.sanger.ac.uk/cosmic/signatures>), we might 
want to just find the exposures `H` for these known signatures in the 
mutational catalogue `V` of a given cohort. 

The **YAPSA**-package (Yet Another Package for Signature Analysis) presented 
here provides the function `LCD` (**l**inear **c**ombination **d**ecomposition) 
to perform this task. The advantage of this method is that there are **no 
constraints on the cohort size**, so `LCD` can be run for as little as one 
sample and thus be used e.g. for signature analysis in personalized oncology. 
In contrast to NMF, `LCD` is very fast and requires very little computational 
resources. The YAPSA package provides additional functions for signature 
analysis, e.g. for stratifying the mutational catalogue to determine signature 
exposures in different strata, part of which is discussed in the vignette of 
the package.


# Install

As long as `YAPSA` is not yet accepted on 
[Bioconductor](https://www.bioconductor.org/) it may be downloaded and installed 
from [github](https://github.com/):


```r
library(devtools)
install_github("huebschm/YAPSA")
```

Of course, `devtools` has to be installed:


```r
install.packages("devtools")
```

`YAPSA` already has preparations to use the newest versions of the pacakges 
`circlize` and `ComplexHeatmap` by Zuguang Gu. These are currently not in the 
release branch of Bioconductor. If you want your system to be ready for the 
next coming update of `YAPSA` you may already now install the newest versions of 
these packages from [github](https://github.com/) as well:


```r
install_github("jokergoo/circlize")
install_github("jokergoo/ComplexHeatmap")
```

If you ran into dependency conflicts before, try rerunning 
`install_github("huebschm/YAPSA")` now.


# Usage

## Example: a cohort of B-cell lymphomas


```r
library(YAPSA)
library(knitr)
opts_chunk$set(echo=TRUE)
opts_chunk$set(fig.show='asis')
```

### Loading example data

Load data in a vcf-like format:


```r
data("lymphoma_Nature2013_raw")
```

Adapt the data structure:


```r
names(lymphoma_Nature2013_raw_df) <- c("PID","TYPE","CHROM","START",
                                       "STOP","REF","ALT","FLAG")
lymphoma_Nature2013_df <- subset(lymphoma_Nature2013_raw_df,TYPE=="subs",
                                 select=c(CHROM,START,REF,ALT,PID))
names(lymphoma_Nature2013_df)[2] <- "POS"
head(lymphoma_Nature2013_df)
```

```
##   CHROM       POS REF ALT      PID
## 1     1 183502381   G   A 07-35482
## 2    18  60985506   T   A 07-35482
## 3    18  60985748   G   T 07-35482
## 4    18  60985799   T   C 07-35482
## 5     2 242077457   A   G 07-35482
## 6     6  13470412   C   T 07-35482
```

Note that there are 48 different samples:


```r
unique(lymphoma_Nature2013_df$PID)
```

```
##  [1] 07-35482       1060           1061           1065          
##  [5] 1093           1096           1102           4101316       
##  [9] 4105105        4108101        4112512        4116738       
## [13] 4119027        4121361        4125240        4133511       
## [17] 4135350        4142267        4158726        4159170       
## [21] 4163639        4175837        4177856        4182393       
## [25] 4189200        4189998        4190495        4193278       
## [29] 4194218        4194891        515            DLBCL-PatientA
## [33] DLBCL-PatientB DLBCL-PatientC DLBCL-PatientD DLBCL-PatientE
## [37] DLBCL-PatientF DLBCL-PatientG DLBCL-PatientH DLBCL-PatientI
## [41] DLBCL-PatientJ DLBCL-PatientK DLBCL-PatientL DLBCL-PatientM
## [45] EB2            FL009          FL-PatientA    G1            
## 48 Levels: 07-35482 1060 1061 1065 1093 1096 1102 4101316 ... G1
```

For convenience later on, we annotate subgroup information to every variant 
(indirectly through the sample it occurs in). For reasons of simplicity, we 
also restrict the analysis to the Whole Genome Sequencing (WGS) datasets:


```r
lymphoma_Nature2013_df$SUBGROUP <- "unknown"
DLBCL_ind <- grep("^DLBCL.*",lymphoma_Nature2013_df$PID)
lymphoma_Nature2013_df$SUBGROUP[DLBCL_ind] <- "DLBCL_other"
MMML_ind <- grep("^41[0-9]+$",lymphoma_Nature2013_df$PID)
lymphoma_Nature2013_df <- lymphoma_Nature2013_df[MMML_ind,]
data(lymphoma_PID)
for(my_PID in rownames(lymphoma_PID_df)) {
  PID_ind <- which(as.character(lymphoma_Nature2013_df$PID)==my_PID)
  lymphoma_Nature2013_df$SUBGROUP[PID_ind] <-
    lymphoma_PID_df$subgroup[which(rownames(lymphoma_PID_df)==my_PID)]
}
lymphoma_Nature2013_df$SUBGROUP <- factor(lymphoma_Nature2013_df$SUBGROUP)
unique(lymphoma_Nature2013_df$SUBGROUP)
```

```
## [1] WGS_D WGS_F WGS_B WGS_I
## Levels: WGS_B WGS_D WGS_F WGS_I
```

As stated [above](#LCD), one of the functions in the YAPSA package (`LCD`) is 
designed to do mutational signatures analysis with known signatures. There are 
(at least) two possible sources for signature data: i) the ones published 
initially by Alexandrov, and ii) an updated and curated 
current set of mutational signatures is maintained by Ludmil Alexandrov at <http://cancer.sanger.ac.uk/cosmic/signatures>.


```r
data(sigs)
current_sig_df <- AlexCosmicValid_sig_df
current_sigInd_df <- AlexCosmicValid_sigInd_df
```

Now we can start using main functions of the YAPSA package: `LCD` and 
`LCD_complex_cutoff`.

### Building a mutational catalogue

Prepare.


```r
library(BSgenome.Hsapiens.UCSC.hg19)
lymphoma_Nature2013_df <- translate_to_hg19(lymphoma_Nature2013_df,"CHROM")
lymphoma_Nature2013_df$change <- 
  attribute_nucleotide_exchanges(lymphoma_Nature2013_df)
lymphoma_Nature2013_df <- 
  lymphoma_Nature2013_df[order(lymphoma_Nature2013_df$PID,
                               lymphoma_Nature2013_df$CHROM,
                               lymphoma_Nature2013_df$POS),]
```

This section uses functions which are to a large extent wrappers for functions 
in the package `SomaticSignatures` by Julian Gehring.


```r
word_length <- 3
lymphomaNature2013_mutCat_list <- 
  create_mutation_catalogue_from_df(
    lymphoma_Nature2013_df,
    this_seqnames.field = "CHROM", this_start.field = "POS",
    this_end.field = "POS", this_PID.field = "PID",
    this_subgroup.field = "SUBGROUP",
    this_refGenome = BSgenome.Hsapiens.UCSC.hg19,
    this_wordLength = word_length)
```

The function `create_mutation_catalogue_from_df` returns a list object with 
several entries. We will use the one called `matrix`.


```r
lymphomaNature2013_mutCat_df <- as.data.frame(
  lymphomaNature2013_mutCat_list$matrix)
```

### LCD analysis with signature-specific cutoffs

When using `LCD_complex_cutoff`, we have to supply a vector of cutoffs with as 
many entries as there are signatures. It may make sense to provide different 
cutoffs for different signatures.


```r
my_cutoff <- 0.06
general_cutoff_vector <- rep(my_cutoff,dim(current_sig_df)[2])
specific_cutoff_vector <- general_cutoff_vector
specific_cutoff_vector[c(1,5)] <- 0
specific_cutoff_vector
```

```
##  [1] 0.00 0.06 0.06 0.06 0.00 0.06 0.06 0.06 0.06 0.06 0.06 0.06 0.06 0.06
## [15] 0.06 0.06 0.06 0.06 0.06 0.06 0.06 0.06 0.06 0.06 0.06 0.06 0.06 0.06
## [29] 0.06 0.06
```

In this example, the cutoff for signatures AC1 and AC5 is thus set to 0, whereas 
the cutoffs for all other signatures remains at 0.06. Running the function
`LCD_complex_cutoff`:


```r
CosmicValid_cutoffSpec_LCDlist <- LCD_complex_cutoff(
  in_mutation_catalogue_df = lymphomaNature2013_mutCat_df,
  in_signatures_df = current_sig_df,
  in_cutoff_vector = specific_cutoff_vector,
  in_sig_ind_df = current_sigInd_df)
```

Some adaptation (extracting and reformatting the information which sample 
belongs to which subgroup):


```r
COSMIC_subgroups_df <- 
  make_subgroups_df(CosmicValid_cutoffSpec_LCDlist$exposures,
                    lymphoma_Nature2013_df)
```

Plotting absolute exposures for visualization:


```r
exposures_barplot(
  in_exposures_df = CosmicValid_cutoffSpec_LCDlist$exposures,
  in_signatures_ind_df = CosmicValid_cutoffSpec_LCDlist$out_sig_ind_df,
  in_subgroups_df = COSMIC_subgroups_df)
```

![Absoute exposures of the COSMIC signatures in the lymphoma mutational catalogues, cutoff of 6% for the LCD (Linear Combination Decomposition).](README_files/figure-html/exposures_barplot_abs_cutoffSpec-1.png)

And relative exposures:


```r
exposures_barplot(
  in_exposures_df = CosmicValid_cutoffSpec_LCDlist$norm_exposures,
  in_signatures_ind_df = CosmicValid_cutoffSpec_LCDlist$out_sig_ind_df,
  in_subgroups_df = COSMIC_subgroups_df)
```

![Relative exposures of the COSMIC signatures in the lymphoma mutational catalogues, cutoff of 6% for the LCD (Linear Combination Decomposition).](README_files/figure-html/exposures_barplot_rel_cutoffSpec-1.png)

Note that the signatures extracted with the signature-specific cutoffs are the 
same in the example displayed here. Depending on the analyzed cohort and the 
choice of cutoffs, the extracted signatures may vary considerably.



## Cluster samples based on their signature exposures

To identify groups of samples which were exposed to similar mutational 
processes, the exposure vectors of the samples can be compared. The YAPSA 
package provides a custom function for this task: `complex_heatmap_exposures`, 
which uses the package *[ComplexHeatmap](http://bioconductor.org/packages/ComplexHeatmap)* by Zuguang Gu. It produces 
output as follows:


```r
complex_heatmap_exposures(CosmicValid_cutoffSpec_LCDlist$norm_exposures,
                          COSMIC_subgroups_df,
                          CosmicValid_cutoffSpec_LCDlist$out_sig_ind_df,
                          in_data_type="norm exposures",
                          in_subgroup_colour_column="col",
                          in_method="manhattan",
                          in_subgroup_column="subgroup")
```

![Clustering of Samples and Signatures based on the relative exposures of the COSMIC signatures in the lymphoma mutational catalogues.](README_files/figure-html/apply_heatmap_exposures-1.png)

If you are interested only in the clustering and not in the heatmap 
information, you could also use `hclust_exposures`:


```r
hclust_list <- 
  hclust_exposures(CosmicValid_cutoffSpec_LCDlist$norm_exposures,
                   COSMIC_subgroups_df,
                   in_method="manhattan",
                   in_subgroup_column="subgroup")
```

The dendrogram produced by either the function `complex_heatmap_exposures` or 
the function `hclust_exposures` can be cut to yield signature exposure specific 
subgroups of the PIDs.


```r
cluster_vector <- cutree(hclust_list$hclust,k=4)
COSMIC_subgroups_df$cluster <- cluster_vector
subgroup_colour_vector <- rainbow(length(unique(COSMIC_subgroups_df$cluster)))
COSMIC_subgroups_df$cluster_col <- 
  subgroup_colour_vector[factor(COSMIC_subgroups_df$cluster)]
complex_heatmap_exposures(CosmicValid_cutoffSpec_LCDlist$norm_exposures,
                          COSMIC_subgroups_df,
                          CosmicValid_cutoffSpec_LCDlist$out_sig_ind_df,
                          in_data_type="norm exposures",
                          in_subgroup_colour_column="cluster_col",
                          in_method="manhattan",
                          in_subgroup_column="cluster")
```

![PIDs labelled by the clusters extracted from the signature analysis.](README_files/figure-html/cluster_PIDs_sig_info-1.png)



## Performing a stratification based on mutation density

This type of analysis is performed using the function `run_SMC` where SMC stands 
for __stratification of the mutational catalogue__. For details on this 
function please consult the vignette.
