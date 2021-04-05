#' calculation of Fold-Change
#'
#'A function to calculate fold-change between group comparison; "Test_group" vs "Ref_group"
#'
#' @import   ExperimentHub testthat ComplexHeatmap ggplot2 matrixStats gtools reshape2 preprocessCore randomcoloR V8 limma
#' @param df_raw  Matrix of normalized expression data (not Log2 transformed).Genes should be in rows and Sample ID in columns. Row names are required to be valid Gene Symbols
#' @param sample_info      A dataframe with sample annotation. Sample_info dataframe requires two columns: 1) a column specifying Sample ID (exactly matching Sample ID of data.matrix) and 2) a column specifying group names
#' @return   A matrix of the fold change comparison between "Test_group" vs ""Ref_group"
#' @examples
#'## data could be downloaded from ExperimentHub("GSE13015")
#'library(ExperimentHub)
#'dat = ExperimentHub()
#'res = query(dat , "GSE13015")
#'GSE13015 = res[["EH5429"]]
#'data_matrix = assay(GSE13015)
#'sample_ann = data.frame(colData(GSE13015))
#'
#' FCgroup = fold_change(data.matrix)
#'
#' @author Darawan Rinchai <drinchai@gmail.com>
#' @export

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






