#' Individual single sample analysis
#'
#'The Individualcomparison function will perform individual sample comparison analysis in reference to a control sample or group of samples, with the results are expressed “at the module level” as percent of genes increased or decreased.
#'
#'- Expression matrix and sample annotation file are required in order to perform this analysis.
#'- The sample annotation file must be loaded using a specific name = "sample.info".
#'- The names of the columns for the conditions used in the analysis must be specified
#'- The default cutoff is set at FC =1.5 and DIFF =10
#'
#' @param data.matrix   Normalized expression data (not Log2 transformed)
#' @param sample_info		 Sample annotation table
#' @param FC               Foldchange cut off to consider the abundance of a given transcript to be increased or decreased compared to a reference group (Ref_group)
#' @param DIFF             Difference cut off to consider the abundance of a given transcript to be increased or decreased compared to a reference group (Ref_group)
#' @param Group_column		 Name of the columns for the groups used for the analysis
#' @param Ref_group 		Reference group or samples that considered as control
#' @return A matrix of the percentahe of module response at individual level
#' @examples
#' Individual_df = Individualcomparison(data.matrix, FC = 1.5, DIFF = 10, Group_column = "Group_test",Ref_group = "Control")
#' @export
#' @author
#' Darawan Rinchai <drinchai@gmail.com>
#
Individualcomparison <- function(data.matrix,
                                 sample_info = sample_info,
                                 FC = NULL,
                                 DIFF = NULL,
                                 Group_column = NULL,
                                 Ref_group = NULL){

  ### Prepare expression matrix with module list
  df1=Module_listGen3                       # This is module list annotation table
  df2=data.frame(data.matrix)               # expression data (from your own datasets or from step 1)
  df2$Gene = rownames(df2)

  #Annotate gene module to expression matrix
  df.mod = merge(df1,df2,by="Gene",all=F)   # match df1 and df2 by Gene symbol

  rownames(df.mod) = df.mod$Module_gene
  dat.mod.func.Gen3 = df.mod[,c(1:8)]
  dat.mod.Gen3 = df.mod[,-c(1:8)]

  ############################
  #prepare data for analysis
  ###########
  df_raw = as.matrix(dat.mod.Gen3)          # replace "dat.mod.Gen3" with data_matrix in raw expression data
  mod_func = dat.mod.func.Gen3              # repleace "mod_func" with Gene module annotation table

  #### make sure that expression matrix and sample information are the same order
  df_raw = df_raw[,rownames(sample_info)]
  colnames(df_raw) == rownames(sample_info)

  # Difference
  Diff.mod.ind.sin <- df_raw[,]
  Diff.mod.ind.sin [,] <- NA

  k=1
  for (k in 1:nrow(df_raw)) {
    signature = rownames(df_raw)[k]
    test.table <- sample_info
    test.table$scores <- df_raw[k,]
    T4 <- test.table
    T3 <- test.table[test.table[, Group_column] == Ref_group,]
    Diff.mod.ind.sin[k,] <- (T4$scores-(mean(T3$scores)))
  }

  Diff.mod.ind.sin <- as.data.frame(Diff.mod.ind.sin)

  ## fold change
  FC.mod.ind.sin <- df_raw[,]
  FC.mod.ind.sin [,] <- NA

  for (k in 1:nrow(df_raw)) {
    signature = rownames(df_raw)[k]
    test.table <- sample_info
    test.table$scores <- df_raw[k,]
    T4 <- test.table
    T3 <- test.table[test.table[, Group_column] == Ref_group,]
    FC.mod.ind.sin[k,] <- foldchange(T4$scores,mean(T3$scores))
  }

  FC.mod.ind.sin <- as.data.frame(FC.mod.ind.sin)


  #############################################
  # Calculate percentage of response ##
  ############################################
  if (is.null(FC)) {
    FC_cutoff = 1.5
  }
  else {
    FC_cutoff = as.numeric(FC)
  }


  if (is.null(DIFF)) {
    DIFF_cutoff = 10
  }
  else {
    DIFF_cutoff = as.numeric(DIFF)
  }
  #logical check ##
  test.up <- (FC.mod.ind.sin > FC_cutoff) + (Diff.mod.ind.sin > DIFF_cutoff) == 2          # TRUE Up gene, Both TRUE
  test.down <- (FC.mod.ind.sin < (FC_cutoff*-1)) + (Diff.mod.ind.sin < (DIFF_cutoff*-1)) == 2      # TRUE down gene, Both TRUE

  ### prepare gene annotation table
  Gene.matrix = dat.mod.func.Gen3[rownames(test.up),]
  #####UP GENE#######
  up.mods.sin <- data.frame(Module = row.names(test.up), test.up,gene=0)                    # create a new blank table
  up.mods.sin [,] <- NA
  up.mods.sin <- up.mods.sin [-c(2:nrow(up.mods.sin)),]

  for (i in 1:length(unique(Gene.matrix$Module))){                                         # length of module
    module <- unique(as.character(Gene.matrix$Module))[i]                                  # look for only unique module
    sums <- colSums(test.up[Gene.matrix$Module==module,])                                  # sum upgene of each column by module
    genes <- nrow(Gene.matrix[Gene.matrix$Module==module,])                                # sum number of gene in each module
    up.mods.sin <- rbind(up.mods.sin,c(module,sums,genes))                                 # paste result into a new fake table
  }

  up.mods.sin <- up.mods.sin [-1,]

  rownames(up.mods.sin) <- up.mods.sin$Module
  up.mods.sin$Module <- NULL
  up.mods.sin.cal <- up.mods.sin
  up.mods.sin <- as.data.frame(lapply(up.mods.sin, as.numeric))                                       # convert data frame to be numberic
  up.mods.sin <- (up.mods.sin/up.mods.sin$gene)*100
  rownames(up.mods.sin) <-rownames(up.mods.sin.cal)
  up.mods.sin <- up.mods.sin[,-ncol(up.mods.sin)]


  #####DOWN GENE#######
  down.mods.sin <- data.frame(Module = row.names(test.down), test.down,gene=0)                         # create a new blank table
  down.mods.sin [,] <- NA                                                                              # create a new blank table
  down.mods.sin <- down.mods.sin [-c(2:nrow(down.mods.sin)),]

  for (i in 1:length(unique(Gene.matrix$Module))){
    module <- unique(Gene.matrix$Module)[i]
    sums <- colSums (test.down[Gene.matrix$Module==module,])
    genes <- nrow(Gene.matrix[Gene.matrix$Module==module,])
    down.mods.sin <- rbind(down.mods.sin,c(module,sums,genes))
  }

  down.mods.sin <- down.mods.sin [-1,]

  rownames(down.mods.sin) <- down.mods.sin$Module
  down.mods.sin$Module <- NULL
  down.mods.sin.cal <- down.mods.sin
  down.mods.sin <- as.data.frame(lapply(down.mods.sin, as.numeric))                                   # convert data frame to be numberic
  down.mods.sin <- (down.mods.sin/down.mods.sin$gene)*100
  rownames(down.mods.sin) <-rownames(down.mods.sin.cal)
  down.mods.sin <- down.mods.sin[,-ncol(down.mods.sin)]


  ## Prepare data for ploting ##
  Sum.mod.sin <- down.mods.sin                                                                       ## prepare a new matrix for new data
  Sum.mod.sin[,] <- NA                                                                               ## Empty matrix

  for (i in 1: nrow(up.mods.sin)){
    for (j in 1:ncol(down.mods.sin)){
      up = up.mods.sin[i,j]
      down = down.mods.sin[i,j]
      Sum.mod.sin[i,j] = up-down
    }
  }
  Individual_df = Sum.mod.sin
}
