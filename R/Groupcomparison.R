#' Group comparison analysis
#'
#' The Groupcomparison function will perform group comparison analyses and the results are expressed “at the module level” as percent of genes increased or decreased.

#' - Expression matrix and sample annotation files are required to perform this analysis.
#' - The sample annotation file must be loaded using a specific name = "sample_info".
#' - The names of the columns for the conditions used in the analysis must be specified.
#' @import                 ExperimentHub testthat ComplexHeatmap ggplot2 matrixStats gtools reshape2 preprocessCore randomcoloR V8 limma
#' @param data.matrix      Matrix of normalized expression data (not Log2 transformed).Genes should be in rows and Sample ID in columns. Row names are required to be valid Gene Symbols
#' @param sample_info      A dataframe with sample annotation. Sample_info dataframe requires two columns: 1) a column specifying Sample ID (exactly matching Sample ID of data.matrix) and 2) a column specifying group names
#' @param FC               Numeric value specifying the foldchange cut off that will be applied to define increase or decrease of a given transcript compared to the reference group
#' @param pval             Numeric value specifying p-value cut off or False discovery rate	when FDR = TRUE
#' @param FDR              Logical operator to specify whether False discovery rate cut off (using BH-method) should be used
#' @param Group_column		 Character vector identical to the column name from sample_info dataframe that specifies group annotation used for the analysis
#' @param Test_group       Character vector specifying values within the group column (Group_column) that will be used as Test group (samples considered as cases or “intervention” group).
#' @param Ref_group 	     Character vector specifying value within the group column (Group_column) that will be used as Reference group
#' @return                 A matrix of the percentahe of module response in each group comparison
#' @examples
#'## data could be downloaded from ExperimentHub("GSE13015")
#'library(ExperimentHub)
#'library(SummarizedExperiment)
#'dat = ExperimentHub()
#'res = query(dat , "GSE13015")
#'GSE13015 = res[["EH5429"]]
#'data_matrix = assay(GSE13015)
#'sample_ann = data.frame(colData(GSE13015))
#'Group_df = Groupcomparison(data_matrix, sample_info = sample_ann,
#'                           FC = 0, pval = 0.1, FDR = TRUE, Test_group = "Sepsis",
#'                           Group_column = "Group_test", Ref_group = "Control")
#' @author Darawan Rinchai <drinchai@gmail.com>
#' @export
Groupcomparison <- function(data.matrix,
                            sample_info = sample_info,
                            FC = NULL,
                            pval = NULL ,
                            FDR = TRUE,
                            Group_column = NULL,
                            Test_group = "Test_group",
                            Ref_group = "Control"){
  globalVariables("Module_listGen3")
  ### Prepare expression matrix with module list
  df1=Module_listGen3                       # This is module list annotation table
  df2=data.frame(data.matrix)               # expression data (from your own datasets or from step 1)
  df2$Gene = rownames(df2)

  #Annotate gene module to expression matrix
  df.mod = merge(df1,df2,by="Gene",all=FALSE)   # match df1 and df2 by Gene symbol

  rownames(df.mod) = df.mod$Module_gene
  dat.mod.func.Gen3 = df.mod[,c(1:8)]
  dat.mod.Gen3 = df.mod[,-c(1:8)]

  #prepare data for analysis
  ###########
  df_raw = as.matrix(dat.mod.Gen3)          # replace "dat.mod.Gen3" with data_matrix in raw expression data
  mod_func = dat.mod.func.Gen3              # repleace "mod_func" with Gene module annotation table

  #### make sure that expression matrix and sample information are the same order
  df_raw = df_raw[,rownames(sample_info)]
  colnames(df_raw) == rownames(sample_info)

  #############################################
  # Statistic analysis ##
  ############################################
  dat_log2 <- as.matrix(log(df_raw+1,2))      # tranformed data to log 2

  ## prepare entry table
  ########################
  ##### T test
  ########################

  tt_pval = data.frame(matrix(ncol = length(Test_group), nrow = nrow(dat_log2)))
  colnames(tt_pval) = Test_group
  rownames(tt_pval) = rownames(dat_log2)

  # Check if rownames of sample_info and colnames of dat_log2 are in the same order before running loop below
  rownames(sample_info) == colnames(dat_log2)


  for (k in 1:nrow(dat_log2)) {
    signature = rownames(dat_log2)[k]
    test.table <- sample_info
    test.table$scores <- dat_log2[k,]
      T2 <- test.table[test.table[, Group_column] == Test_group,]        # "Group_test"; the selected column could be changed to your interested group comparison
      T1 <- test.table[test.table[, Group_column] == Ref_group,]         # "Group_test"; the selected column could be changed to your interested group comparison
      if(mean(T1$scores) == mean(T2$scores)){
        tt_pval[signature,] = 1
      }else{
        tt_pval[signature,] <- t.test(x =T1$scores,y=T2$scores,paired = FALSE,var.equal = TRUE)$p.value
      }
    }


  pvalue_Group <- data.frame(tt_pval)

  pvalue_Group.FDR <- apply(pvalue_Group,2,function(x) p.adjust(x,method = "fdr"))     ## Apply multiple correction testing
  pvalue_Group.FDR <- as.data.frame(pvalue_Group.FDR)

  if(FDR == "TRUE"){
    Pvalue_cutoff = pvalue_Group.FDR
  }else{
    Pvalue_cutoff = pvalue_Group
  }

  ####################################
  ####calculate fold change ##
  ####################################

  FCgroup = fold_change(df_raw = df_raw,
                        sample_info = sample_info,
                        Group_column = Group_column,
                        Test_group=Test_group,
                        Ref_group=Ref_group)

  #############################################
  # Calculate percentage of response ##
  ############################################

  if (is.null(FC)) {
    FC_cutoff = 0
  }
  else {
    FC_cutoff = as.numeric(FC)
  }
  FC_cutoff = as.numeric(FC)

  if (is.null(pval)) {
    pval = 0.1
  }
  else {
    pval = as.numeric(pval)
  }


  #logical check ##
  mod.up = (FCgroup > FC_cutoff) + (Pvalue_cutoff < pval) == 2             # TRUE Up gene, Both TRUE
  mod.down = (FCgroup < (FC_cutoff*-1)) + (Pvalue_cutoff < pval) == 2      # TRUE down gene, Both TRUE

  ################################################

  ### prepare gene annotation table
  Gene.matrix = mod_func[rownames(mod.up),]

  #####UP GENE#######
  pect_df <- data.frame(Module = row.names(mod.up), mod.up,genes=0)                    # create a new blank table
  pect_df [,] <- NA
  pect_df <- pect_df [-c(2:nrow(pect_df)),]

  for (i in 1:length(unique(Gene.matrix$Module))){                                         # length of module
    module <- unique(as.character(Gene.matrix$Module))[i]                                  # look for only unique module
    sums_up <- colSums(mod.up[Gene.matrix$Module==module,1,drop=FALSE])                    # sum upgene of each column by module
    sums_down <- colSums(mod.down[Gene.matrix$Module==module,1,drop=FALSE])
    sums = sums_up-sums_down
    genes <- nrow(Gene.matrix[Gene.matrix$Module==module,])                                # sum number of gene in each module
    pect_df <- rbind(pect_df,c(module,sums,genes))                                         # paste result into a new fake table
  }

  pect_df <- pect_df [-1,]

  rownames(pect_df) <- pect_df$Module
  pect_df$Module <- NULL
  pect_df.cal <- pect_df
  pect_df <- as.data.frame(lapply(pect_df, as.numeric))                                    # convert data frame to be numberic
  pect_df <- (pect_df/pect_df$genes)*100
  rownames(pect_df) <-rownames(pect_df.cal)
  pect_df <- pect_df[,-ncol(pect_df),drop=FALSE]
  Group_df = pect_df
}
