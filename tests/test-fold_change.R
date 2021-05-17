library("gtools")
library("testthat")
library("BloodGen3Module")

fold_change = function(df_raw = df_raw,
                       sample_info = sample_info,
                       Group_column = Group_column,
                       Test_group=Test_group,
                       Ref_group=Ref_group){

  FC.group = data.frame(matrix(ncol = 1, nrow = nrow(df_raw)))
  colnames(FC.group) = Test_group
  rownames(FC.group) = rownames(df_raw)

  for (k in 1:nrow(df_raw)) {
    signature = rownames(df_raw)[k]
    test.table <- sample_info
    test.table$scores <- df_raw[k,]
    T2 <- test.table[test.table[, Group_column]==Test_group,]       # "Group_test"; the selected column could be changed to your interested group comparison
    T1 <- test.table[test.table[, Group_column]==Ref_group,]        # "Group_test"; the selected column could be changed to your interested group comparison
    FC.group[signature,] <- foldchange(mean(T2$scores),mean(T1$scores))
  }
  FCgroup <- data.frame(FC.group)
}


library(ExperimentHub)
library(SummarizedExperiment)
dat = ExperimentHub()
res = query(dat , "GSE13015")
GSE13015 = res[["EH5429"]]
data_matrix = assay(GSE13015)
sample_ann = data.frame(colData(GSE13015))


test_that("test fold_change", {
  a  = FCgroup = fold_change(df_raw = data_matrix[c(1:5),],
                                            sample_info = sample_ann,
                                            Group_column = "Group_test",
                                            Test_group="Sepsis",
                                            Ref_group="Control")
  b  = FCgroup = fold_change(df_raw = data_matrix[c(1:5),],
                                            sample_info = sample_ann,
                                            Group_column = "Group_test",
                                            Test_group="Sepsis",
                                            Ref_group="Control")

  expect_that(a, equals(b))
})
