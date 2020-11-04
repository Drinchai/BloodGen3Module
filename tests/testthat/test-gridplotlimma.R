library(BloodGen3Module)
library(testthat)
library(reshape2)
library(ggplot2)

#Load expression data
#Example expression data for package testting

Test_sample =  matrix(data = rexp(1000, rate = 0.01), nrow = 14168, ncol = 20)
control_sample = matrix(data = rexp(1000, rate = 0.1), nrow = 14168, ncol = 10)

data.matrix = data.frame(cbind(Test_sample,control_sample))
data.matrix$Symbol = Module_listGen3$Gene
data.matrix = aggregate(data.matrix,FUN = mean,by=list(data.matrix$Symbol))
rownames(data.matrix) = data.matrix$Group.1
data.matrix$Group.1 = NULL
data.matrix$Symbol = NULL
colnames(data.matrix) = c(paste0(rep("SampleID",30),1:30))

##example sample information
sample_ann = data.frame(SampleID=(colnames(data.matrix)),Group_test = c(rep("Test",20),rep("Control",10)),stringsAsFactors = F)
rownames(sample_ann) = sample_ann$SampleID
rownames(sample_ann) == colnames(data.matrix)
head(sample_ann)

Group_limma <- Groupcomparisonlimma(data.matrix,
                                    sample_info = sample_ann,
                                    FC = 1.5,
                                    pval = 0.1 ,
                                    FDR = TRUE,
                                    Group_column = "Group_test",
                                    Test_group = "Test",
                                    Ref_group = "Control")

test_that("test gridplotlimma", {

  a  = gridplotlimma(Group_limma,
                cutoff = 15,
                Ref_group = "Control",
                filename="Limma_group_comparison")
  b  = gridplotlimma(Group_limma,
                     cutoff = 15,
                     Ref_group = "Control",
                     filename="Limma_group_comparison")

  expect_that(a, equals(b))
})
