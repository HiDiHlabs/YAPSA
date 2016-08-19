test_that("Test cosine_dist with 2-dim vectors", {
  a <- c(1,0)
  b <- c(0,1)
  d <- c(1,1)
  expect_that(YAPSA:::cosineDist(a,b), equals(1))
  expect_that(YAPSA:::cosineDist(a,a), equals(0))
  expect_that(YAPSA:::cosineDist(b,b), equals(0))
  expect_lt(YAPSA:::cosineDist(a,d) - (1-cos(pi/4)),1e-08)
})


test_that("Test makeVRangesFromDataFrame with real data", {
  data(lymphoma_test)
  temp_vr <- makeVRangesFromDataFrame(lymphoma_test_df,
                                      in_seqnames.field="CHROM",
                                      in_subgroup.field="SUBGROUP")
  expect_that(length(temp_vr),equals(dim(lymphoma_test_df)[1]))
  expect_that(as.character(seqnames(temp_vr)),
              is_equivalent_to(lymphoma_test_df$CHROM))
  expect_that(as.data.frame(ranges(temp_vr))[,1],
              is_equivalent_to(lymphoma_test_df$POS))
})


test_that("Test cut_breaks_as_intervals with partially real data",{
  my_precision <- 0.05
  data(lymphoma_test)
  my_outlier_cutoffs <- c(-4,4)
  my_cutoff_ranges_list <- list(c(-2.5,-1.5),c(0.5,1.5))
  lymphoma_test_df$random_norm <- rnorm(dim(lymphoma_test_df)[1])
  temp_list <- 
    cut_breaks_as_intervals(lymphoma_test_df$random_norm,
                            in_outlier_cutoffs=my_outlier_cutoffs,
                            in_cutoff_ranges_list=my_cutoff_ranges_list,
                            in_labels=c("small","intermediate","big"))
  expect_that(temp_list$cutoffs[c(1,4)],equals(my_outlier_cutoffs))
  expect_lt(max(abs(temp_list$cutoffs[c(2,3)]-
                        c(min(my_cutoff_ranges_list[[1]]),
                          max(my_cutoff_ranges_list[[2]])))),
                   my_precision)
})


test_that("Test stratify_vcf_like_df with real data",{
  data(lymphoma_test)
  my_labels <- c("small","intermediate","big")
  strata_list <- 
    cut_breaks_as_intervals(lymphoma_test_df$random_norm,
                            in_outlier_cutoffs=c(-4,4),
                            in_cutoff_ranges_list=list(c(-2.5,-1.5),
                                                       c(0.5,1.5)),
                            in_labels=my_labels)
  lymphoma_test_df$random_cat <- strata_list$category_vector
  stratification_list <- 
    stratify_vcf_like_df(lymphoma_test_df,"random_cat",in_verbose=0)
  expect_that(length(stratification_list$name_list),
              equals(length(unique(lymphoma_test_df$random_cat))))
  expect_that(sort(unlist(stratification_list$name_list)),
              is_equivalent_to(sort(my_labels)))
  expect_that(length(stratification_list$table_list),
              equals(length(unique(lymphoma_test_df$random_cat))))
  sum <- 0
  for (i in seq_len(3)){
    temp_table <- stratification_list$table_list[[i]]
    sum <- sum + dim(temp_table)[1]
    expect_that(as.character(unique(temp_table$random_cat)),
                equals(stratification_list$name_list[[i]]))
  }
  expect_that(sum,equals(dim(lymphoma_test_df)[1]))
})