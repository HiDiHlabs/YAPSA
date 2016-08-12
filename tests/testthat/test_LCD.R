test_that("Test LCD (linear combination decomposition) with small hand made data frames", {
  ## define raw data
  W_prim <- matrix(c(1,2,3,4,5,6),ncol=2) 
  W_prim_df <- as.data.frame(W_prim)
  W_df <- YAPSA:::normalize_df_per_dim(W_prim_df,2) # corresponds to the signatures
  W <- as.matrix(W_df) 

  ## 1. Simple case: non-negativity already in raw data
  H <- matrix(c(2,5,3,6,1,9,1,2),ncol=4)
  H_df <- as.data.frame(H) # corresponds to the exposures
  V <- W %*% H # matrix multiplication
  V_df <- as.data.frame(V) # corresponds to the mutational catalogue
  ## apply function to be tested
  exposures_df <- YAPSA:::LCD(V_df,W_df)
  ## compare
  expect_less_than(max(abs(exposures_df - H_df)),1e-05)
  
  ## 2. more complicated: raw data already contains negative elements
  ## define indices where sign is going to be swapped
  sign_ind <- c(5,7)
  ## now compute the indices of the other fields in the columns affected
  ## by the sign change
  row_ind <- sign_ind %% dim(H)[1]
  temp_ind <- 2*row_ind -1
  other_ind <- sign_ind + temp_ind
  ## alter the matrix H to yield a new mutational catalogue
  H_compl <- H
  H_compl[sign_ind] <- (-1)*H[sign_ind]
  H_compl_df <- as.data.frame(H_compl) # corresponds to the exposures
  V_compl <- W %*% H_compl # matrix multiplication
  V_compl_df <- as.data.frame(V_compl) # corresponds to the mutational catalogue
  ## apply function to be tested
  exposures_df <- YAPSA:::LCD(V_compl_df,W_df)
  exposures <- as.matrix(exposures_df)
  expect_that(exposures[sign_ind], equals(rep(0,length(sign_ind))))
  expect_that(all(exposures[other_ind]<H_compl[other_ind]),is_true())
})

test_that("Test LCD_cutoff with small hand made data frames", {
  ## define raw data
  W_prim <- matrix(c(1,2,3,4,5,6),ncol=2) 
  W_prim_df <- as.data.frame(W_prim)
  W_df <- YAPSA:::normalize_df_per_dim(W_prim_df,2) # corresponds to the signatures
  W <- as.matrix(W_df) 
  H <- matrix(c(2,5,3,6,1,9,1,2),ncol=4)
  H_df <- as.data.frame(H) # corresponds to the exposures
  V <- W %*% H # matrix multiplication
  V_df <- as.data.frame(V) # corresponds to the mutational catalogue
  ## apply function to be tested
  exposures_small_cutoff_list <- YAPSA:::LCD_cutoff(V_df,W_df,in_cutoff = 0.05)
  exposures_big_cutoff_list <- YAPSA:::LCD_cutoff(V_df,W_df,in_cutoff = 0.4)
  ## compare
  expect_less_than(max(abs(exposures_small_cutoff_list$exposures - H_df)),1e-05)
  expect_that(dim(exposures_big_cutoff_list$exposures)[1], equals(dim(H_df)[1] - 1))
})

test_that("Test LCD_SMC with small hand made data frames", {
  ## define raw data
  W_prim <- matrix(c(1,2,3,4,5,6),ncol=2) 
  W_prim_df <- as.data.frame(W_prim)
  W_df <- YAPSA:::normalize_df_per_dim(W_prim_df,2) # corresponds to the signatures
  W <- as.matrix(W_df) 
  H <- matrix(c(2,5,3,6,1,9,1,2),ncol=4)
  ## define indices where sign is going to be swapped for different strata (perturbation)
  sign_ind_1 <- c(5,7)
  sign_ind_2 <- c(1)
  ## compute the column-wise index complements for the different strata for later use
  col_ind_1 <- sign_ind_1 %/% dim(H)[1]
  other_col_ind_1 <- setdiff(c(0:(dim(H)[2]-1)),col_ind_1)
  other_ind_1 <- rep(other_col_ind_1*dim(H)[1],each=dim(H)[1]) + rep(seq_len(dim(H)[1]),length(other_col_ind_1))
  col_ind_2 <- sign_ind_2 %/% dim(H)[1]
  other_col_ind_2 <- setdiff(c(0:(dim(H)[2]-1)),col_ind_2)
  other_ind_2 <- rep(other_col_ind_2*dim(H)[1],each=dim(H)[1]) + rep(seq_len(dim(H)[1]),length(other_col_ind_2))
  ## alter the matrix H to yield a new mutational catalogue for every stratum (perturbation)
  H_1 <- H
  H_2 <- H
  H_1[sign_ind_1] <- (-1)*H[sign_ind_1]
  H_2[sign_ind_2] <- (-1)*H[sign_ind_2]
  H_1_df <- as.data.frame(H_1)
  H_2_df <- as.data.frame(H_2)
  V_1 <- W %*% H_1 # matrix multiplication
  V_2 <- W %*% H_2 # matrix multiplication
  V <- V_1 + V_2
  V_1_df <- as.data.frame(V_1) # corresponds to the mutational catalogue of stratum 1
  V_2_df <- as.data.frame(V_2) # corresponds to the mutational catalogue of stratum 2
  V_df <- as.data.frame(V) # corresponds to the mutational catalogue of the whole cohort
  V_list <- list()  # make list of data frames
  V_list[[1]] <- V_1_df
  V_list[[2]] <- V_2_df
  ## apply function to be tested
  exposures_strata_list <- YAPSA:::LCD_SMC(V_list,W_df)
  ## compare
  simple_exposures_all_df <- YAPSA:::LCD(V_df,W_df)
  simple_exposures_1_df <- YAPSA:::LCD(V_1_df,W_df)
  simple_exposures_1 <- as.matrix(simple_exposures_1_df)
  simple_exposures_2_df <- YAPSA:::LCD(V_2_df,W_df)
  simple_exposures_2 <- as.matrix(simple_exposures_2_df)
  ## check if the overall decomposition is equal to the result of the normal LCD
  expect_that(exposures_strata_list$exposures_all_df, equals(simple_exposures_all_df))
  ## check if the sum of the exposures of the strata equals the result of the normal LCD
  expect_that(exposures_strata_list$sub_exposures_list[[1]] + exposures_strata_list$sub_exposures_list[[2]],
              equals(simple_exposures_all_df))
  ## check if the positions not affected by the perturbation remain the same
  expect_lt(max(abs(as.matrix(exposures_strata_list$sub_exposures_list[[1]])[other_ind_2]
                           - simple_exposures_1[other_ind_2])),1e-05)
  expect_lt(max(abs(as.matrix(exposures_strata_list$sub_exposures_list[[2]])[other_ind_1]
                           - simple_exposures_2[other_ind_1])),1e-05)
})
