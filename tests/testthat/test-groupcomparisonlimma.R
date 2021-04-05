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
data_matrix = assay(GSE13015)
sample_ann = data.frame(colData(GSE13015))
head(sample_ann)

test_that("test Groupcomparisonlimma", {

  a  = Group_limma <- Groupcomparisonlimma(data_matrix,
                                           sample_info = sample_ann,
                                           FC = 1.5,
                                           pval = 0.1 ,
                                           FDR = TRUE,
                                           Group_column = "Group_test",
                                           Test_group = "Test",
                                           Ref_group = "Control")
  b  = Group_limma <- Groupcomparisonlimma(data_matrix,
                                           sample_info = sample_ann,
                                           FC = 1.5,
                                           pval = 0.1 ,
                                           FDR = TRUE,
                                           Group_column = "Group_test",
                                           Test_group = "Test",
                                           Ref_group = "Control")

  expect_that(a, equals(b))
})
