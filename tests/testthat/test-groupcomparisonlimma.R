library(BloodGen3Module)
library(testthat)
library(reshape2)
library(ggplot2)

#Load expression data
library(ExperimentHub)
library(SummarizedExperiment)

dat = ExperimentHub()
res = query(dat , "GSE13015")
GSE13015 = res[["EH5429"]]

test_that("test Groupcomparisonlimma", {

  a  = Group_limma <- Groupcomparisonlimma(GSE13015,
                                           sample_info = NULL,
                                           FC = 1.5,
                                           pval = 0.1 ,
                                           FDR = TRUE,
                                           Group_column = "Group_test",
                                           Test_group = "Test",
                                           Ref_group = "Control")
  b  = Group_limma <- Groupcomparisonlimma(GSE13015,
                                           sample_info = NULL,
                                           FC = 1.5,
                                           pval = 0.1 ,
                                           FDR = TRUE,
                                           Group_column = "Group_test",
                                           Test_group = "Test",
                                           Ref_group = "Control")

  expect_that(a, equals(b))
})
