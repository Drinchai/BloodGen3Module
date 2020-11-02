#' Group comparison analysis using limma
#'
#'The Groupcomparisonlimma function will perform group comparison analyses using "limma" function from "limma R package" and the results are expressed “at the module level” as percent of genes increased or decreased.

#'- Expression matrix and sample annotation files are required to perform this analysis.
#'- The sample annotation file must be loaded using a specific name = "sample_info".
#'- The names of the columns for the conditions used in the analysis must be specified.
#' @import              ComplexHeatmap ggplot2 matrixStats gtools reshape2 preprocessCore randomcoloR V8 limma
#' @param data.matrix   A table of normalized expression matrix (Important: Expression matrix must be none Log2 transformed as it will be automatic transformed when running this function)
#' @param sample_info   A table of sample information (Important: rownames of sample information must be the same names as in colnames of data.matrix)
#' @param FC            Foldchange cut off to consider th eabundance of a given transcript to be increased or decreased compared to a reference group (Ref_group)
#' @param pval          p-value cut off or False discovery rate when FDR = FALSE
#' @param FDR           False discovery rate cut off using BH method
#' @param Group_column  Name of the columns for the groups used for the analysis
#' @param Test_group 		Characters name of test group or samples that considered as test (Example: Sepsis, Cancer, RSV, Bacteria,.. etc.)
#' @param Ref_group 		Characters name of reference group or samples that considered as control (Example: Control, baseline, Pre-treatment,... etc)
#' @return              A matrix of the percentahe of module response in each group comparison
#' @examples
#' Group_df <- Groupcomparisonlimma(data.matrix,sample_info = sample_info, FC = 0, pval = 0.1, FDR = TRUE, Group_column = "Group_test", Test_group = "Sepsis", Ref_group = "Control")
#' @author
#' Darawan Rinchai <drinchai@gmail.com>
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
  df.mod = merge(df1,df2,by="Gene",all=F)   # match df1 and df2 by Gene symbol

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


  pvalue_Group <- data.frame(diff.stats[,"adj.P.Val"])
  rownames(pvalue_Group) = rownames(diff.stats)
  colnames(pvalue_Group) = "pvalue"

  pvalue_Group.FDR <- data.frame(diff.stats[,"P.Value"])
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

  FC.group = data.frame(matrix(ncol = length(Test_group), nrow = nrow(df_raw)))
  colnames(FC.group) = Test_group
  rownames(FC.group) = rownames(df_raw)

  k=1
  for (k in 1:nrow(df_raw)) {
    signature = rownames(df_raw)[k]
    test.table <- sample_info
    test.table$scores <- df_raw[k,]
    for (i in 1:length(Test_group)) {
      group = Test_group[i]
      T2 <- test.table[test.table[, Group_column]==Test_group,]             
      T1 <- test.table[test.table[, Group_column]== Ref_group,]             
      FC.group[signature,group] <- foldchange(mean(T2$scores),mean(T1$scores))
    }
  }

  FCgroup <- data.frame(FC.group)

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
  Group.up = (FCgroup > FC_cutoff) + (Pvalue_cutoff < pval) == 2        # TRUE Up gene, Both TRUE

  Group.down = (FCgroup < (FC_cutoff*-1)) + (Pvalue_cutoff < pval) == 2      # TRUE down gene, Both TRUE


  ################################################
  Gene.matrix <- dat.mod.func.Gen3[rownames(Group.up),]
  Gene.matrix$Module <- as.character(Gene.matrix$Module)

  up.mods.group <- data.frame(matrix(ncol =2+length(Test_group), nrow = 1))        # create a new blank table
  colnames(up.mods.group) = c("Module", Test_group, "genes")

  i=1
  for (i in 1:length(unique(Gene.matrix$Module))){                                    # length of module
    module <- unique(Gene.matrix$Module)[i]                                           # look for only unique module
    sums <- colSums(Group.up[Gene.matrix$Module==module,1,drop=F])                            # sum upgene of each column by module
    genes <- nrow(dat.mod.func.Gen3[dat.mod.func.Gen3$Module==module,1,drop=F])               # sum number of gene in each module
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
  up.mods.group <- up.mods.group[,-ncol(up.mods.group),drop=F]



  #####DOWN GENE#######
  down.mods.group <- data.frame(matrix(ncol =2+length(Test_group), nrow = 1))        # create a new blank table
  colnames(down.mods.group) = c("Module", Test_group, "genes")

  for (i in 1:length(unique(Gene.matrix$Module))){
    module <- unique(Gene.matrix$Module)[i]
    sums <- colSums (Group.down[Gene.matrix$Module==module,1,drop=F])
    genes <- nrow(dat.mod.func.Gen3[dat.mod.func.Gen3$Module==module,1,drop=F])
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
  down.mods.group <- down.mods.group[,-ncol(down.mods.group),drop=F]

  ## Prepare data for ploting ##
  res.mods.group <- up.mods.group[,1,drop=F]                                     # prepare a new matrix for new data
  res.mods.group[,1:ncol(res.mods.group)] <- NA                                  # Empty matrix

  i=1
  for (i in 1: nrow(up.mods.group)){
    for (j in 1:ncol(up.mods.group)){
      up = up.mods.group[i,j]
      down = down.mods.group[i,j]
      res.mods.group[i,j] = up-down
    }
  }
  Group_limma = res.mods.group
}
