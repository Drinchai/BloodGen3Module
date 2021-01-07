#' Group comparison analysis using limma
#'
#' The Groupcomparisonlimma function will perform group comparison analyses using "limma" function from "limma R package" and the results are expressed “at the module level” as percent of genes increased or decreased.

#' - Expression matrix and sample annotation files are required to perform this analysis.
#' - The sample annotation file must be loaded using a specific name = "sample_info".
#' - The names of the columns for the conditions used in the analysis must be specified.
#' @import              ComplexHeatmap ggplot2 matrixStats gtools reshape2 preprocessCore randomcoloR V8 limma
#' @param data.matrix   Matrix of normalized expression data (not Log2 transformed). Genes should be in rows and Sample ID in columns.Row names are required to be valid Gene Symbols
#' @param sample_info   A dataframe with sample annotation.
#' @param FC            Numeric value specifying the foldchange cut off that will be applied to define increase or decrease of a given transcript compared to the reference group
#' @param pval          Numeric value specifying p-value cut off or False discovery rate when FDR = TRUE
#' @param FDR           Logical operator  to specify whether False discovery rate cut off (using BH-method) should be used
#' @param Group_column  Character vector identical to the column name from sample_info dataframe that specifies group annotation used for the analysis
#' @param Test_group 		Character vector specifying value within the group column that will be used as Test group
#' @param Ref_group 		Character vector specifying value within the group column that will be used as Reference group
#' @return              A matrix of the percentahe of module response in each group comparison
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
#'Group_limma <- Groupcomparisonlimma(data.matrix, sample_info = sample_ann,
#'FC = 1.5, pval = 0.1, FDR = TRUE, Group_column = "Group_test",
#'Test_group = "Test", Ref_group = "Control")
#' @author Darawan Rinchai <drinchai@gmail.com>
#' @export


Groupcomparisonlimma <- function(data.matrix,
                            sample_info = sample_info,
                            FC = NULL,
                            pval = NULL ,
                            FDR = TRUE,
                            Group_column = NULL,
                            Test_group = "Test_group",
                            Ref_group = "Control"){

  ### Prepare expression matrix with module list
  df1=Module_listGen3                   # This is module list annotation table
  df2=data.frame(data.matrix)                  # expression data (from your own datasets or from step 1)
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
  dat_log2 <- as.matrix(log(df_raw+1,2))      # tranform data to log 2

  ########################
  ##### limma test

  # Analysis
  Expression.matrix = as.data.frame(t(dat_log2))
  Expression.matrix$Group = sample_info[,Group_column][match(rownames(Expression.matrix),rownames(sample_info))]
  levels(factor(Expression.matrix$Group))


  design = model.matrix(~ Expression.matrix$Group)
  Expression.matrix$Group = NULL
  Expression.matrix = as.matrix(Expression.matrix)
  mode(Expression.matrix) = "numeric"

  fit = lmFit(t(Expression.matrix) , design=design)
  fit = eBayes(fit)
  diff.stats = topTable(fit, coef=2,adjust.method = "BH", number = ncol(Expression.matrix))


  pvalue_Group <- data.frame(diff.stats[,"P.Value"])
  rownames(pvalue_Group) = rownames(diff.stats)
  colnames(pvalue_Group) = "pvalue"

  pvalue_Group.FDR <- data.frame(diff.stats[,"adj.P.Val"])
  rownames(pvalue_Group.FDR) = rownames(diff.stats)
  colnames(pvalue_Group.FDR) = "FDR"

  if(FDR == "TRUE"){
    Pvalue_cutoff = pvalue_Group.FDR
  }else{
    Pvalue_cutoff = pvalue_Group
  }

  ####################################
  ####calculate fold change ##
  ####################################

  FCgroup = fold_change(df_raw)

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
  mod.up = (FCgroup > FC_cutoff) + (Pvalue_cutoff < pval) == 2        # TRUE Up gene, Both TRUE

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
    sums_up <- colSums(mod.up[Gene.matrix$Module==module,1,drop=FALSE])                                  # sum upgene of each column by module
    sums_down <- colSums(mod.down[Gene.matrix$Module==module,1,drop=FALSE])
    sums = sums_up-sums_down
    genes <- nrow(Gene.matrix[Gene.matrix$Module==module,])                                # sum number of gene in each module
    pect_df <- rbind(pect_df,c(module,sums,genes))                                 # paste result into a new fake table
  }

  pect_df <- pect_df [-1,]

  rownames(pect_df) <- pect_df$Module
  pect_df$Module <- NULL
  pect_df.cal <- pect_df
  pect_df <- as.data.frame(lapply(pect_df, as.numeric))                                       # convert data frame to be numberic
  pect_df <- (pect_df/pect_df$genes)*100
  rownames(pect_df) <-rownames(pect_df.cal)
  pect_df <- pect_df[,-ncol(pect_df),drop=F]
  Group_df = pect_df
}
