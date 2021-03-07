#' Fingerprint grid visualization
#'
#' The gridplot function will generate a grid plot as a pdf file. Specific working directory for the analysis need to be specified for saving the file. The result of the plot should be return in the same working directory.
#' The default cut off for visualization is set at 15%, it can be changed to any value between 0-100%.
#' @import            testthat ComplexHeatmap ggplot2 matrixStats gtools reshape2 preprocessCore randomcoloR V8 limma
#' @param Group_df    Dataframe with output generated after running the'Groupcomparison' function
#' @param cutoff 			Numeric value specifying the percentage cut off used for fingerprint visualization (from 0 to 100)
#' @param Ref_group 	Character vector specifying value within the group column that will be used as Reference group
#' @param filename	  Character vector with a name for saving file
#' @return            A pdf file of grid plot
#' @examples
#' #' ## example sample information Example expression
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
#'gridplot(Group_df, cutoff = 15,
#'          Ref_group = "Control",
#'          filename="Group_comparison_cutoff15")
#' @author Darawan Rinchai <drinchai@gmail.com>
#' @export

gridplot = function(Group_df,
                    cutoff = NULL,
                    Ref_group = NULL,
                    filename=NULL){

  ## prepared cluster position
  Group_plot = Group_df
  Group_plot <-Group_plot[rownames(Gen3_ann),1,drop=FALSE]
  rownames(Group_plot)==rownames(Gen3_ann)                         # check if rownames is the same
  rownames(Group_plot) <- Gen3_ann$position
  Group_plot <- as.data.frame(Group_plot)

  ########## An example of DISPLAY DATA > 15 %
  if (is.null(cutoff)) {
    cutoff = 15
  }
  else {
    cutoff = as.numeric(cutoff)
  }
  Group_plot[abs(Group_plot) < cutoff ] <- 0

  ##creat new grid with all filtered cluster##
  mod.group1 <- matrix(nrow=38,ncol=42)
  rownames(mod.group1) = paste0("A",seq_len(38))
  colnames(mod.group1) = paste0("",seq_len(42))

  diseases = colnames(Group_plot)
  N.disease = length(diseases)

  for (i in 1:N.disease){
    disease = diseases[i]
    if(disease == Ref_group){next}
    for (i in 1:nrow(Group_plot)){
      Mx = as.numeric(gsub(x = strsplit(rownames(Group_plot)[i],"\\.")[[1]][[1]],pattern = "A",replacement = ""))
      My = as.numeric(strsplit(rownames(Group_plot)[i],"\\.")[[1]][[2]])
      mod.group1[Mx,My] <- Group_plot[,disease][i]
    }
    mod.group = mod.group1[-c(9:14,19:23),]
    melt_test = melt(mod.group,id.var=c("row.names"))
    colnames(melt_test) = c("Aggregate","Sub_aggregate","%Response")
    pdf(paste0(filename,"_",disease, "vs",Ref_group,"_Gridplot.pdf"), height = 5.5, width = 8.5)
    plot = ggplot(melt_test, aes(Aggregate, as.factor(Sub_aggregate))) +
      geom_tile(color="#E6E6E6" , size = 0.2, fill=color )+
      geom_point(aes(colour=`%Response`),size=4.3)+
      ylab("") +
      xlab("") +
      labs(title= disease)+
      theme(axis.text.x = element_text(angle = -90, hjust = 0))+
      scale_color_gradient2(low = "blue", mid="white", high = "red",limits=c(-100,100), na.value = "#E6E6E6", guide = "colourbar")+
      theme_light() +
      theme(panel.grid.minor = element_line(colour="black", size=0.9))+
      coord_flip() +
      scale_x_discrete(limits = rev(levels(melt_test$Aggregate))) +
      theme(panel.border = element_rect(color = "black",size = 0.5),
            axis.text.x = element_text(colour="black",size=9,angle=0,hjust=0.5,vjust=2,face="plain"),
            axis.text.y = element_text(colour="black",size=9,angle=0,hjust=0.5,vjust=0.5,face="plain"))

    plot(plot)
    dev.off()
  }
}
