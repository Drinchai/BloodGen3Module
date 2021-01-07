#' calculation of Fold-Change
#'
#'A function to calculate fold-change between group comparison; "Test_group" vs "Ref_group"
#'
#' @import   SummarizedExperiment testthat ComplexHeatmap ggplot2 matrixStats gtools reshape2 preprocessCore randomcoloR V8 limma
#' @param x  Matrix of normalized expression data (not Log2 transformed).Genes should be in rows and Sample ID in columns. Row names are required to be valid Gene Symbols
#' @return   A matrix of the fold change comparison between "Test_group" vs ""Ref_group"
#' @examples
#' ## example sample information Example expression
#' ## data for package testting
#'Test_sample = matrix(data = rexp(1000, rate = 0.01),
#'                     nrow = 14168, ncol = 20)
#'control_sample = matrix(data = rexp(1000, rate = 0.1),
#'                        nrow = 14168, ncol = 10)
#'data.matrix = data.frame(cbind(Test_sample, control_sample))
#'data.matrix$Symbol = Module_listGen3$Gene
#'data.matrix = aggregate(data.matrix[,-31], FUN = mean, by = list(data.matrix$Symbol))
#'rownames(data.matrix) = data.matrix$Group.1
#'data.matrix$Group.1 = NULL
#'colnames(data.matrix) = c(paste0(rep("SampleID", 30),
#'                                 1:30))
#'## Example of ample information
#'sample_ann = data.frame(SampleID = (colnames(data.matrix)),
#'                        Group_test = c(rep("Test", 20), rep("Control",
#'                                                            10)), stringsAsFactors = FALSE)
#'rownames(sample_ann) = sample_ann$SampleID
#'
#' FCgroup = fold_change(data.matrix)
#'
#' @author Darawan Rinchai <drinchai@gmail.com>
#' @export

fold_change = function(x){
  df_raw=x
  FC.group = data.frame(matrix(ncol = 1, nrow = nrow(df_raw)))
  colnames(FC.group) = Test_group
  rownames(FC.group) = rownames(df_raw)

  k=1
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






