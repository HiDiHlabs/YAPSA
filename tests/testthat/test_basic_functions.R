test_that("Test cosine_dist with 2-dim vectors", {
  a <- c(1,0)
  b <- c(0,1)
  d <- c(1,1)
  expect_that(YAPSA:::cosineDist(a,b), equals(1))
  expect_that(YAPSA:::cosineDist(a,a), equals(0))
  expect_that(YAPSA:::cosineDist(b,b), equals(0))
  expect_lt(YAPSA:::cosineDist(a,d) - (1-cos(pi/4)),1e-08)
})


test_that("Test compare_sets", {
  ## define two sets of vectors (signatures)
  sig_1_df <- data.frame(matrix(c(1,0,0,0,0,1,0,0,0,0,1,0),ncol=3))
  names(sig_1_df) <- paste0("B",seq_len(dim(sig_1_df)[2]))
  sig_2_df <- data.frame(matrix(c(1,1,0,0,0,0,1,1),ncol=2))
  ## define control parameters
  dummy <- (1-cos(pi/4))
  distance_matrix <- data.frame(matrix(c(dummy,dummy,1,1,1,dummy),ncol=2))
  ## execute function to be tested
  test_list <- YAPSA:::compare_sets(sig_1_df,sig_2_df)
  ## compare
  expect_that(test_list$distance, is_equivalent_to(distance_matrix))
  expect_that(test_list$hierarchy_small[1,], is_equivalent_to(names(sig_2_df)))
  expect_that(test_list$hierarchy_big[1,], is_equivalent_to(names(sig_1_df)))
})


test_that("Test repeat_df", {
  test_df <- data.frame(matrix(rep(1,6),ncol=3))
  expect_that(sum(abs(YAPSA:::repeat_df(1,2,3)-test_df)), equals(0))
  test_df <- data.frame(matrix(rep(NA,6),ncol=2))
  expect_that(YAPSA:::repeat_df(NA,3,2), is_equivalent_to(test_df))
  test_df <- data.frame(matrix(rep("a",12),ncol=3))
  expect_that(YAPSA:::repeat_df("a",4,3), is_equivalent_to(test_df))
})


test_that("Test normalize_df_per_dim", {
  ## define a test data set
  test_df <- data.frame(matrix(c(1,2,3,0,5,2,3,4,0,6,0,0,0,0,0,4,5,6,0,7),ncol=4))
  ## 1. test row-wise
  choice_rows <- which(rowSums(test_df)>0)
  ## 1.a run function
  row_norm_df <- YAPSA:::normalize_df_per_dim(test_df,1)
  ## 1.b compare
  expect_that(rowSums(row_norm_df[choice_rows,]), is_equivalent_to(rep(1,length(choice_rows))))
  expect_that(rowSums(row_norm_df[-choice_rows,]), is_equivalent_to(rep(0,dim(test_df)[1]-length(choice_rows))))
  ## 2. test column-wise
  choice_cols <- which(colSums(test_df)>0)
  ## 2.a run function
  col_norm_df <- YAPSA:::normalize_df_per_dim(test_df,2)
  ## 2.b compare
  expect_that(colSums(col_norm_df[,choice_cols]), is_equivalent_to(rep(1,length(choice_cols))))
  expect_that(sum(col_norm_df[,-choice_cols]), is_equivalent_to(0))
})


test_that("Test average_over_present", {
  ## define a test data set
  test_df <- data.frame(matrix(c(1,2,3,0,5,2,3,4,0,6,0,0,0,0,0,4,5,6,0,7),ncol=4))
  ## 1. test row-wise
  choice_cols <- which(colSums(test_df)>0)
  ## 1.a run function
  row_mean_vec <- YAPSA:::average_over_present(test_df,1)
  ## 1.b compare
  expect_that(apply(test_df[,choice_cols],1,mean), is_equivalent_to(row_mean_vec))
  ## 2. test column-wise
  choice_rows <- which(rowSums(test_df)>0)
  ## 2.a run function
  col_mean_vec <- YAPSA:::average_over_present(test_df,2)
  ## 2.b compare
  expect_that(apply(test_df[choice_rows,],2,mean), is_equivalent_to(col_mean_vec))
})


test_that("Test sd_over_present", {
  ## define a test data set
  test_df <- data.frame(matrix(c(1,2,3,0,5,2,3,4,0,6,0,0,0,0,0,4,5,6,0,7),ncol=4))
  ## 1. test row-wise
  choice_cols <- which(colSums(test_df)>0)
  ## 1.a run function
  row_sd_vec <- YAPSA:::sd_over_present(test_df,1)
  ## 1.b compare
  expect_that(apply(test_df[,choice_cols],1,sd), is_equivalent_to(row_sd_vec))
  ## 2. test column-wise
  choice_rows <- which(rowSums(test_df)>0)
  ## 2.a run function
  col_sd_vec <- YAPSA:::sd_over_present(test_df,2)
  ## 2.b compare
  expect_that(apply(test_df[choice_rows,],2,sd), is_equivalent_to(col_sd_vec))
})


test_that("Test sum_over_list_of_df", {
  ## define a test data set
  A <- data.frame(matrix(c(1,1,1,2,2,2),ncol=2))
  B <- data.frame(matrix(c(3,3,3,4,4,4),ncol=2))
  compare_df <- data.frame(matrix(c(4,4,4,6,6,6),ncol=2))
  ## 1. for a named list
  df_list <- list(A=A,B=B)
  ## compare
  result_df <- sum_over_list_of_df(df_list)
  expect_that(class(result_df), equals("data.frame"))
  expect_that(result_df, is_equivalent_to(compare_df))
  ## 2. for a non-named list
  df_list <- list()
  df_list[[1]] <- A
  df_list[[2]] <- B
  ## compare
  result_df <- sum_over_list_of_df(df_list)
  expect_that(class(result_df), equals("data.frame"))
  expect_that(result_df, is_equivalent_to(compare_df))
})


test_that("Test the stderrmean function with small dummy vector",{
  A <- c(1,2,3)
  expect_lt(abs(stderrmean(A)-(1/(sqrt(3)))),1e-10)
})


test_that("Test stderrmean_over_present", {
  ## define a test data set
  test_df <- data.frame(matrix(c(1,2,3,0,5,2,3,4,0,6,0,0,0,0,0,4,5,6,0,7),ncol=4))
  ## 1. test row-wise
  choice_cols <- which(colSums(test_df)>0)
  ## 1.a run function
  row_stderrmean_vec <- YAPSA:::stderrmean_over_present(test_df,1)
  ## 1.b compare
  expect_that(apply(test_df[,choice_cols],1,stderrmean), is_equivalent_to(row_stderrmean_vec))
  ## 2. test column-wise
  choice_rows <- which(rowSums(test_df)>0)
  ## 2.a run function
  col_stderrmean_vec <- YAPSA:::stderrmean_over_present(test_df,2)
  ## 2.b compare
  expect_that(apply(test_df[choice_rows,],2,stderrmean), is_equivalent_to(col_stderrmean_vec))
})


test_that("Test translate_to_hg19 on very simple synthetic data", {
  ## define a test data set
  test_df <- data.frame(CHROM=c(1,2,23,24),POS=c(100,120000000,300000,25000),dummy=c("a","b","c","d"))
  ## run function
  hg19_df <- translate_to_hg19(test_df, in_CHROM.field = "CHROM")
  ## compare
  expect_that(hg19_df$CHROM,equals(c("chr1","chr2","chrX","chrY")))
  expect_that(hg19_df[,c(2,3)],equals(test_df[,c(2,3)]))
})


test_that("Test translate_to_1kG on very simple synthetic data", {
  ## define a test data set
  test_df <- data.frame(CHROM=c(1,2,23,24),POS=c(100,120000000,300000,25000),dummy=c("a","b","c","d"))
  hg19_df <- translate_to_hg19(test_df, in_CHROM.field = "CHROM")
  ## run function
  onekG_df <- translate_to_1kG(hg19_df, in_CHROM.field = "CHROM")
  ## compare
  expect_that(onekG_df$CHROM,equals(as.character(test_df$CHROM)))
  expect_that(onekG_df[,c(2,3)],equals(test_df[,c(2,3)]))
})


test_that("Test attribute_nucleotide_exchanges on very simple synthetic data", {
  ## define a test data set
  test_df <- data.frame(CHROM=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5),
                        POS=c(1,2,3,4,5,6,1,2,3,4,5,6,7,8),
                        REF=c("C","C","C","T","T","T","A","A","A","G","G","G","N","A"),
                        ALT=c("A","G","T","A","C","G","C","G","T","A","C","T","A","N"))
  ## run function
  test_df$change <- attribute_nucleotide_exchanges(test_df,in_REF.field = "REF",in_ALT.field = "ALT")
  ## compare
  expect_that(test_df$change,
              equals(factor(c("CA","CG","CT","TA","TC","TG","TG","TC","TA","CT","CG","CA",NA,NA),
                            levels=c("CA", "CG", "CT", "TA", "TC", "TG"))))
})


test_that("Test annotate_intermut_dist_PID on very simple synthetic data", {
  ## define a test data set
  test_df <- data.frame(CHROM=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5),
                        POS=c(1,2,4,4,6,9,1,4,8,10,20,40,100,200),
                        REF=c("C","C","C","T","T","T","A","A","A","G","G","G","N","A"),
                        ALT=c("A","G","T","A","C","G","C","G","T","A","C","T","A","N"))
  ## run function
  min_dist_df <- annotate_intermut_dist_PID(test_df,in_CHROM.field="CHROM",in_POS.field="POS",
                                        in_mode="min")
  max_dist_df <- annotate_intermut_dist_PID(test_df,in_CHROM.field="CHROM",in_POS.field="POS",
                                            in_mode="max")
  ## compare
  expect_that(min_dist_df$dist,equals(c(1,1,2,2,2,3,3,3,4,10,10,20,100,100)))
  expect_that(max_dist_df$dist,equals(c(1,2,2,2,3,3,3,4,4,10,20,20,100,100)))
})


test_that("Test annotate_intermut_dist_cohort on very simple synthetic data", {
  ## define a test data set
  test_df <- data.frame(CHROM=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5),
                        POS=c(1,2,4,4,6,9,1,4,8,10,20,40,100,200),
                        REF=c("C","C","C","T","T","T","A","A","A","G","G","G","N","A"),
                        ALT=c("A","G","T","A","C","G","C","G","T","A","C","T","A","N"),
                        PID=c(1,1,1,2,2,2,1,1,2,2,2,1,1,2))
  test_df <- test_df[order(test_df$PID,test_df$CHROM,test_df$POS),]
  ## run function
  min_dist_df <- annotate_intermut_dist_cohort(test_df,in_CHROM.field="CHROM",in_POS.field="POS",
                                               in_PID.field="PID",in_mode="min")
  max_dist_df <- annotate_intermut_dist_cohort(test_df,in_CHROM.field="CHROM",in_POS.field="POS",
                                               in_PID.field="PID",in_mode="max")
  ## compare
  expect_that(min_dist_df$dist,equals(c(1,1,2,3,3,1e+08,1e+08,2,2,3,1e+08,10,10,1e+08)))
  expect_that(max_dist_df$dist,equals(c(1,2,2,3,3,1e+08,1e+08,2,3,3,1e+08,10,10,1e+08)))
})


test_that("Test shapiro_if_possible on very simple synthetic data", {
  significance_threshold <- 0.05
  set.seed(1)
  expect_lt(shapiro_if_possible(runif(100,min=2,max=4)),significance_threshold)
  set.seed(1)
  expect_more_than(shapiro_if_possible(rnorm(100,mean=5,sd=3)),significance_threshold)
  expect_that(shapiro_if_possible(rep(4.3,100)),equals(0))
  expect_null(shapiro_if_possible(c("Hello","World")))
})


test_that("Test make_subgroups_df with real data",{
  data(lymphoma_test)
  data(lymphoma_cohort_LCD_results)
  choice_ind <- (names(lymphoma_Nature2013_COSMIC_cutoff_exposures_df) 
                 %in% unique(lymphoma_test_df$PID))
  lymphoma_test_exposures_df <- lymphoma_Nature2013_COSMIC_cutoff_exposures_df[,choice_ind]
  test_subgroups_df <- make_subgroups_df(lymphoma_test_exposures_df,lymphoma_test_df)
  real_subgroups_df <- aggregate(SUBGROUP~PID,data=lymphoma_test_df,
                                 function(l) return(l[1]))
  real_subgroups_df[,2] <- as.character(real_subgroups_df[,2])
  expect_that(test_subgroups_df[,c(1,2)],is_equivalent_to(real_subgroups_df))
})