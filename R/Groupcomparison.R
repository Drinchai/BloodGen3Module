#' Group comparison analysis
#'
#' The Groupcomparison function will perform group comparison analyses and the results are expressed “at the module level” as percent of genes increased or decreased.

#' - Expression matrix and sample annotation files are required to perform this analysis.
#' - The sample annotation file must be loaded using a specific name = "sample_info".
#' - The names of the columns for the conditions used in the analysis must be specified.
#' @import                 testthat ComplexHeatmap ggplot2 matrixStats gtools reshape2 preprocessCore randomcoloR V8 limma
#' @param data.matrix      Matrix of normalized expression data (not Log2 transformed).Genes should be in rows and Sample ID in columns. Row names are required to be valid Gene Symbols
#' @param sample_info      A dataframe with sample annotation.
#' @param FC               Numeric value specifying the foldchange cut off that will be applied to define increase or decrease of a given transcript compared to the reference group
#' @param pval             Numeric value specifying p-value cut off or False discovery rate	when FDR = TRUE
#' @param FDR              Logical operator to specify whether False discovery rate cut off (using BH-method) should be used
#' @param Group_column		 Character vector identical to the column name from sample_info dataframe that specifies group annotation used for the analysis
#' @param Ref_group 	     Character vector specifying value within the group column (Group_column) that will be used as Reference group
#' @param Gen3_ann         Dataframe of functional annotation of third generation module
#' @param Module_listGen3  Dataframe of module list annotation
#' @param Color            Character vector color of grid plot

#' @return                 A matrix of the percentahe of module response in each group comparison
#' @examples
#' ## example sample information Example expression
#'## data for package testting
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
#'Group_df = Groupcomparison(data.matrix, sample_info = sample_ann,
#'                           FC = 0, pval = 0.1, FDR = TRUE,
#'                           Group_column = "Group_test", Ref_group = "Control")
#' @author Darawan Rinchai <drinchai@gmail.com>
#' @export
Groupcomparison <- function(data.matrix,
                            sample_info = sample_info,
                            FC = NULL,
                            pval = NULL ,
                            FDR = TRUE,
                            Group_column = NULL,
                            Ref_group = NULL){

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
  group.test = unique(sample_info[, Group_column])

  ########################
  ##### T test
  ########################

  tt_pval = data.frame(matrix(ncol = length(group.test), nrow = nrow(dat_log2)))
  colnames(tt_pval) = group.test
  rownames(tt_pval) = rownames(dat_log2)

  # Check if rownames of sample_info and colnames of dat_log2 are in the same order before running loop below
  rownames(sample_info) == colnames(dat_log2)

  k=1
  for (k in 1:nrow(dat_log2)) {
    signature = rownames(dat_log2)[k]
    test.table <- sample_info
    test.table$scores <- dat_log2[k,]
    i=1
    for (i in 1:length(group.test)) {
      group = group.test[i]
      T2 <- test.table[test.table[, Group_column] == group,]             # "Group_test"; the selected column could be changed to your interested group comparison
      T1 <- test.table[test.table[, Group_column] == Ref_group,]         # "Group_test"; the selected column could be changed to your interested group comparison
      if(mean(T1$scores) == mean(T2$scores)){
        tt_pval[signature,group] = 1
      }else{
        tt_pval[signature,group] <- t.test(x =T1$scores,y=T2$scores,paired = FALSE,var.equal = TRUE)$p.value
      }
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

  FC.group = data.frame(matrix(ncol = length(group.test), nrow = nrow(df_raw)))
  colnames(FC.group) = group.test
  rownames(FC.group) = rownames(df_raw)

  k=1
  for (k in 1:nrow(df_raw)) {
    signature = rownames(df_raw)[k]
    test.table <- sample_info
    test.table$scores <- df_raw[k,]
    for (i in 1:length(group.test)) {
      group = group.test[i]
      T2 <- test.table[test.table[, Group_column]==group,]             # "Group_test"; the selected column could be changed to your interested group comparison
      T1 <- test.table[test.table[, Group_column]== Ref_group,]        # "Group_test"; the selected column could be changed to your interested group comparison
      FC.group[signature,group] <- foldchange(mean(T2$scores),mean(T1$scores))
    }
  }

  FCgroup <- data.frame(FC.group)


  #############################################
  # Calculate percentage of response ##
  ############################################

  if (is.null(FC)) {
    FC_cutoff = 1.5
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
  Group.up = (FCgroup > FC_cutoff) + (Pvalue_cutoff < pval) == 2             # TRUE Up gene, Both TRUE

  Group.down = (FCgroup < (FC_cutoff*-1)) + (Pvalue_cutoff < pval) == 2      # TRUE down gene, Both TRUE


  ################################################
  Gene.matrix <- dat.mod.func.Gen3[rownames(Group.up),]
  Gene.matrix$Module <- as.character(Gene.matrix$Module)

  up.mods.group <- data.frame(matrix(ncol =2+length(group.test), nrow = 1))        # create a new blank table
  colnames(up.mods.group) = c("Module", group.test, "genes")

  i=1
  for (i in 1:length(unique(Gene.matrix$Module))){                                    # length of module
    module <- unique(Gene.matrix$Module)[i]                                           # look for only unique module
    sums <- colSums(Group.up[Gene.matrix$Module==module,])                            # sum upgene of each column by module
    genes <- nrow(dat.mod.func.Gen3[dat.mod.func.Gen3$Module==module,])               # sum number of gene in each module
    up.mods.group <- rbind(up.mods.group,c(module,sums,genes))                        # paste result into a new fake table
  }

  # Calculate percentage of genes that are upregulated in each module
  up.mods.group <-up.mods.group[-1,]
  rownames(up.mods.group) <- up.mods.group$Module
  up.mods.group$Module <- NULL
  up.mods.group.cal <- up.mods.group
  up.mods.group <- as.data.frame(lapply(up.mods.group, as.numeric))                    # convert data frame to be numeric
  up.mods.group <- (up.mods.group/up.mods.group$genes)*100
  rownames(up.mods.group) <-rownames(up.mods.group.cal)
  up.mods.group <- up.mods.group[,-ncol(up.mods.group)]



  #####DOWN GENE#######
  down.mods.group <- data.frame(matrix(ncol =2+length(group.test), nrow = 1))        # create a new blank table
  colnames(down.mods.group) = c("Module", group.test, "genes")

  for (i in 1:length(unique(Gene.matrix$Module))){
    module <- unique(Gene.matrix$Module)[i]
    sums <- colSums (Group.down[Gene.matrix$Module==module,])
    genes <- nrow(dat.mod.func.Gen3[dat.mod.func.Gen3$Module==module,])
    down.mods.group <- rbind(down.mods.group,c(module,sums,genes))
  }
  down.mods.group<-down.mods.group[-1,]

  # Calculate percentage of genes that are downregulated in each module
  rownames(down.mods.group) <- down.mods.group$Module
  down.mods.group$Module <- NULL
  down.mods.group.cal <- down.mods.group
  down.mods.group <- as.data.frame(lapply(down.mods.group, as.numeric))
  down.mods.group <- (down.mods.group/down.mods.group$genes)*100
  rownames(down.mods.group) <- rownames(down.mods.group.cal)
  down.mods.group <- down.mods.group[,-ncol(down.mods.group)]

  ## Prepare data for ploting ##
  res.mods.group <- up.mods.group[,]                                     # prepare a new matrix for new data
  res.mods.group[,1:ncol(res.mods.group)] <- NA                          # Empty matrix

  i=1
  for (i in 1: nrow(up.mods.group)){
    for (j in 1:ncol(up.mods.group)){
      up = up.mods.group[i,j]
      down = down.mods.group[i,j]
      res.mods.group[i,j] = up-down
    }
  }
  Group_df = res.mods.group
}
