#' This is a function for fingerprint grid plot
#'
#' This function will perform grid plot in pdf file.
#' We assumes that you setup the specific working directory for your analysis. The result of the plot should be in your working folder.
#' The result of the plot should be in your working folder.
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
gridplot = function(Group_df, cutoff = NULL){
  Module_ann = Gen3_ann
  colnames(Module_ann) = "Module_col"

  Module_ann  <- cSplit(Module_ann, "Module_col", sep = ",", direction = "wide", fixed = TRUE,
                        drop = TRUE)

  colnames(Module_ann) = c("Module","Cluster","Cluster_location","Function","position")
  Module_ann = as.data.frame(Module_ann)
  rownames(Module_ann) = Module_ann$Module
  Module_ann$Module_color = Module.platte

  ## prepared cluter position
  Group_plot = Group_df
  Group_plot <-Group_plot[rownames(Module_ann),]
  rownames(Group_plot)==rownames(Module_ann)                         # check if rownames is the same
  rownames(Group_plot) <- Module_ann$position
  Group_plot <- as.data.frame(Group_plot)

  ########## An example of DISPLAY DATA > 15 %
  if (is.null(cutoff)) {
    cutoff = 15
  }
  else {
    cutoff = as.numeric(cutoff)
  }
  Group_plot[abs(Group_plot) < cutoff ] <- 0

  #check data
  head(Group_plot)

  # creat new grid with all filtered cluster##
  mod.group1 <- matrix(nrow=38,ncol=42)
  rownames (mod.group1) <- paste0("A",c(1:38))
  colnames (mod.group1) <- paste0("",c(1:42))
  ##

  diseases = colnames(Group_plot)
  N.disease = length(diseases)

  for (i in 1:N.disease){
    disease = diseases[i]
    for (i in 1 : nrow(Group_plot)){
      Mx <- as.numeric(gsub(x = strsplit (rownames(Group_plot)[i],"\\.")[[1]][[1]],pattern = "A",replacement = ""))
      My <- as.numeric(strsplit (rownames(Group_plot)[i],"\\.")[[1]][[2]])
      mod.group1[Mx,My] <- Group_plot[,disease][i]
    }
    mod.group <- mod.group1[-c(9:14,19:23),]
    melt_test <- melt(mod.group,id.var=c("row.names"))
    colnames(melt_test) = c("Aggregate","Sub_aggregate","%Response")
    pdf(paste0("Group_comparison_", disease, "vsControl_Grid.pdf"), height = 5.5, width = 8.5)
    plot = ggplot(melt_test, aes(Aggregate, as.factor(Sub_aggregate))) +
      geom_tile(color="#E6E6E6" , size = 0.2, fill=color )+
      geom_point(aes(colour=`%Response`),size=4.5)+
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
