#' Individual single sample analysis
#'
#' The Individualcomparison function will perform individual sample comparison analysis in reference to a control sample or group of samples, with the results are expressed “at the module level” as percent of genes increased or decreased.
#'
#' - Expression matrix and sample annotation file are required in order to perform this analysis.
#' - The sample annotation file must be loaded using a specific name = "sample_info".
#' - The names of the columns for the conditions used in the analysis must be specified
#' - The default cutoff is set at FC =1.5 and DIFF =10
#' @import              testthat ComplexHeatmap ggplot2 matrixStats gtools reshape2 preprocessCore randomcoloR V8 limma
#' @param data.matrix   Matrix of normalized expression data (not Log2 transformed).Genes should be in rows and Sample ID in columns. Row names are required to be valid Gene Symbols
#' @param sample_info   A dataframe with sample annotation.
#' @param FC            Numeric value specifying the foldchange cut off that will be applied to define increase or decrease of a given transcript compared to the reference group
#' @param DIFF          Numeric value specifying the difference cut off that will be applied to define increase or decrease of a given transcript compared to the reference group
#' @param Group_column  Character vector identical to the column name from sample_info dataframe that specifies group annotation used for the analysis
#' @param Ref_group 		Character vector specifying value within the group column that will be used as Reference group
#' @return              A matrix of the percentahe of module response at individual level and SummarizedExperiment object
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
#'Individual_df = Individualcomparison(data.matrix, sample_info = sample_ann,
#'                                     FC = 1.5, DIFF = 10, Group_column = "Group_test",
#'                                     Ref_group = "Control")
#' @author Darawan Rinchai <drinchai@gmail.com>
#' @export


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
  df.mod = merge(df1,df2,by="Gene",all=FALSE)   # match df1 and df2 by Gene symbol

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
  mod.up <- (FC.mod.ind.sin > FC_cutoff) + (Diff.mod.ind.sin > DIFF_cutoff) == 2                  # TRUE Up gene, Both TRUE
  mod.down <- (FC.mod.ind.sin < (FC_cutoff*-1)) + (Diff.mod.ind.sin < (DIFF_cutoff*-1)) == 2      # TRUE down gene, Both TRUE

  ### prepare gene annotation table
  Gene.matrix = mod_func[rownames(mod.up),]

  #####UP GENE#######
  pect_df <- data.frame(Module = row.names(mod.up), mod.up,genes=0)                    # create a new blank table
  pect_df [,] <- NA
  pect_df <- pect_df [-c(2:nrow(pect_df)),]

  for (i in 1:length(unique(Gene.matrix$Module))){                                         # length of module
    module <- unique(as.character(Gene.matrix$Module))[i]                                  # look for only unique module
    sums_up <- colSums(mod.up[Gene.matrix$Module==module,])                                  # sum upgene of each column by module
    sums_down <- colSums(mod.down[Gene.matrix$Module==module,])
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
  pect_df <- pect_df[,-ncol(pect_df)]

  Individual_df = pect_df

}
