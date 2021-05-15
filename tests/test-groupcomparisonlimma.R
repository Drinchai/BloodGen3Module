library("limma")
library("testthat")
load("./R/sysdata.rda")
Groupcomparisonlimma <- function(data.matrix,
                                 sample_info = NULL,
                                 FC = NULL,
                                 pval = NULL ,
                                 FDR = TRUE,
                                 Group_column = NULL,
                                 Test_group = "Test_group",
                                 Ref_group = "Control",
                                 SummarizedExperiment = TRUE){

  if(is(data.matrix, "SummarizedExperiment")){
    data_matrix = assay(data.matrix)
  }else{
    data_matrix = data.matrix
  }

  #Sample information
  if (is.null(sample_info)) {
    sample_info = data.frame(colData(data.matrix))
  }
  else {
    sample_info = sample_info
  }
  ### Prepare expression matrix with module list
  df1=Module_listGen3                   # This is module list annotation table
  df2=data.frame(data_matrix)                  # expression data (from your own datasets or from step 1)
  df2$Gene = rownames(df2)

  #Annotate gene module to expression matrix
  df.mod = merge(df1,df2,by="Gene",all=FALSE)   # match df1 and df2 by Gene symbol

  rownames(df.mod) = df.mod$Module_gene
  dat.mod.func.Gen3 = df.mod[,c(1:5)]
  dat.mod.Gen3 = df.mod[,-c(1:5)]

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
  pect_df <- pect_df[,-ncol(pect_df),drop=FALSE]
  Group_df = pect_df
  Group_res <- SummarizedExperiment(assays=SimpleList(Percent=as.matrix(Group_df)))

  if (SummarizedExperiment == "TRUE") {
    Group_df = Group_res
  }
  else {
    Group_df = Group_df
  }
}


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
