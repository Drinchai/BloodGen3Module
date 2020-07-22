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
  Module_ann = data.frame(Module_ann)
  Module_ann$Module = as.character(Module_ann$Module)
  Module_ann$Cluster = as.character(Module_ann$Cluster)
  Module_ann$Function = as.character(Module_ann$Function)
  Module_ann$position = as.character(Module_ann$position)

  rownames(Module_ann) = Module_ann$Module
  
  Module.platte= c("#99C68E", "#FFFFFF", "#C226D9", "#00D7D6", "#F85B4D", "#EA7C1E", "#00FFFF", "#8B63FF", "#6BA3D7", "#FF00FF",
                  "#8B63FF", "#FFFFFF", "#EA7C1E", "#FFFFFF", "#FFABA8", "#004181", "#99C68E", "#E5C493", "#FFABA8", "#0077BA",
                  "#8B63FF", "#FFEE07", "#C5A44B", "#B00D29", "#EA7C1E", "#F94A8B", "#B51F83", "#FFFFFF", "#99C68E", "#B00D29",
                  "#E5C493", "#FFFFFF", "#8B63FF", "#000000", "#E5C493", "#FFFFFF", "#FFFFFF", "#E5C493", "#8B63FF", "#FFFFFF",
                  "#B00D29", "#E5C493", "#B497DB", "#E5C493", "#F4CA49", "#00FFFF", "#FFEE07", "#C5A44B", "#99C68E", "#F85B4D",
                  "#FF92D5", "#FFFFFF", "#FFFFFF", "#C5A44B", "#F94A8B", "#FFEE07", "#8B63FF", "#FFABA8", "#EA7C1E", "#FFFFFF",
                  "#99C68E", "#E5C493", "#8B63FF", "#FFABA8", "#E5C493", "#E5C493", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF",
                  "#FFFFFF", "#FFFFFF", "#99C68E", "#B51F83", "#FFABA8", "#FFABA8", "#FFFFFF", "#E5C493", "#FFFFFF", "#B00D29",
                  "#E5C493", "#FFFFFF", "#FFFFFF", "#E5C493", "#B00D29", "#FFFFFF", "#B00D29", "#8B63FF", "#FF92D5", "#FFFFFF",
                  "#B51F83", "#99C68E", "#FFFFFF", "#B497DB", "#E5C493", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#99C68E",
                  "#CFE84D", "#FFFFFF", "#0077BA", "#FFFFFF", "#8B63FF", "#FFFFFF", "#FFFFFF", "#8B63FF", "#B00D29", "#E5C493",
                  "#B00D29", "#F94A8B", "#FFFFFF", "#EA7C1E", "#99C68E", "#E5C493", "#FFFFFF", "#FFABA8", "#8B63FF", "#FFFFFF",
                  "#FFFFFF", "#FFABA8", "#FFFFFF", "#B51F83", "#99C68E", "#004181", "#B00D29", "#E5C493", "#8B63FF", "#E5C493",
                  "#FFFFFF", "#E5C493", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#8B63FF", "#66CCFF", "#E5C493", "#FFFFFF", "#8B63FF",
                  "#8B63FF", "#6BA3D7", "#F4CA49", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#E5C493", "#FFEE07", "#FFABA8", "#99C68E",
                  "#004181", "#E5C493", "#0077BA", "#E5C493", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#4FF300",
                  "#FFFFFF", "#FFFFFF", "#E5C493", "#FFFFFF", "#E5C493", "#FFFFFF", "#000000", "#FFFFFF", "#FFFFFF", "#FF92D5",
                  "#99C68E", "#B51F83", "#000000", "#FFFFFF", "#FFABA8", "#FFFFFF", "#0077BA", "#0077BA", "#FF92D5", "#FFABA8",
                  "#B00D29", "#0077BA", "#FFFFFF", "#E5C493", "#FFFFFF", "#E5C493", "#FFFFFF", "#FFABA8", "#6BA3D7", "#FFFFFF",
                  "#4FF300", "#FFFFFF", "#000000", "#FFFFFF", "#8B63FF", "#FFFFFF", "#EA7C1E", "#FFFFFF", "#C5A44B", "#4F8F00",
                  "#FFFFFF", "#004181", "#FFFFFF", "#FFFFFF", "#004181", "#FFFFFF", "#FFFFFF", "#00FFFF", "#FFFFFF", "#FFFFFF",
                  "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFABA8", "#FFFFFF", "#FFFFFF", "#99C68E", "#EA7C1E", "#FFFFFF", "#FFFFFF",
                  "#FFFFFF", "#FFFFFF", "#F7C599", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#E5C493", "#F4CA49", "#FFFFFF", "#00FFFF",
                  "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFABA8", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF",
                  "#000000", "#0077BA", "#FFABA8", "#EA7C1E", "#FFFFFF", "#C226D9", "#FFFFFF", "#FFFFFF", "#B00D29", "#FFFFFF",
                  "#FFFFFF", "#FFFFFF", "#B00D29", "#99C68E", "#FFFFFF", "#FFFFFF", "#B00D29", "#FFFFFF", "#FFFFFF", "#FFFFFF",
                  "#FFFFFF", "#E5C493", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#E5C493",
                  "#00FFFF", "#FFFFFF", "#8B63FF", "#0077BA", "#FFFFFF", "#FFFFFF", "#004181", "#FFFFFF", "#FFFFFF", "#FFABA8",
                  "#FFFFFF", "#8B63FF", "#FFEE07", "#FFABA8", "#FFFFFF", "#000000", "#FFFFFF", "#E5C493", "#FFFFFF", "#FFFFFF",
                  "#FFFFFF", "#FFFFFF", "#FFFFFF", "#99C68E", "#0077BA", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF",
                  "#4E564D", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#E5C493", "#FFFFFF", "#E5C493", "#FFFFFF", "#FFFFFF", "#FFFFFF",
                  "#004181", "#FFFFFF", "#FFFFFF", "#E5C493", "#8B63FF", "#FFFFFF", "#FFFFFF", "#A7F9D1", "#FFFFFF", "#B00D29",
                  "#FFABA8", "#8B63FF", "#929000", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#8B63FF", "#FFEE07", "#FFFFFF", "#FFFFFF",
                  "#99C68E", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#C7CEF7", "#009051", "#FFFFFF", "#FFFFFF", "#8B63FF", "#942193",
                  "#FFFFFF", "#FFFFFF", "#FFFFFF", "#8B63FF", "#FFFFFF", "#000000", "#FFFFFF", "#4FF300", "#FFEE07", "#8B63FF",
                  "#F4CA49", "#FFFFFF", "#E5C493", "#4F8F00", "#C5A44B", "#FFFFFF", "#F94A8B", "#FFFFFF", "#EA7C1E", "#FFFFFF",
                  "#FFFFFF", "#FFFFFF", "#99C68E", "#FFFFFF", "#FFFFFF", "#FFEE07", "#EA7C1E", "#000000", "#FFFFFF", "#FFABA8",
                  "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#99C68E", "#FFFFFF", "#FFFFFF", "#F94A8B", "#FFFFFF", "#6BA3D7", 
                  "#FFFFFF", "#FFFFFF")
  
  Module_ann$Module_color = Module.platte
  
  #Grid color
  color <- c(rep("white",27),
           rep("white",27),
           rep("white",9),"#E6E6E6","white","#E6E6E6",rep("white",15),
           rep("white",9),"#E6E6E6","white","#E6E6E6",rep("white",15),
           rep("white",9),"#E6E6E6","white","#E6E6E6",rep("white",14),"#E6E6E6",
           rep("white",9),rep("#E6E6E6",3),rep("white",3),"#E6E6E6",rep("white",2),"#E6E6E6",rep("white",3),"#E6E6E6",rep("white",3),"#E6E6E6",
           rep("white",6),"#E6E6E6",rep("white",2),rep("#E6E6E6",3),rep("white",3),rep("#E6E6E6",2),"white","#E6E6E6",rep("white",3),"#E6E6E6","white","#E6E6E6","white","#E6E6E6",
           rep("white",6),"#E6E6E6",rep("white",2),rep("#E6E6E6",3),rep("white",3),rep("#E6E6E6",2),"white","#E6E6E6",rep("white",3),"#E6E6E6","white","#E6E6E6","white","#E6E6E6",
           rep("white",6),rep("#E6E6E6",2),"white",rep("#E6E6E6",3),rep("white",3),rep("#E6E6E6",5),rep("white",2),"#E6E6E6","white","#E6E6E6","white","#E6E6E6",
           rep("white",6),rep("#E6E6E6",2),"white",rep("#E6E6E6",3),rep("white",3),rep("#E6E6E6",5),rep("white",2),"#E6E6E6","white","#E6E6E6","white","#E6E6E6",
           rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",2),"white",rep("#E6E6E6",3),rep("white",3),rep("#E6E6E6",5),rep("white",2),"#E6E6E6","white","#E6E6E6","white","#E6E6E6",
           rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",7),rep("white",2),rep("#E6E6E6",5),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
           rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",8),"white",rep("#E6E6E6",5),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
           rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",8),"white",rep("#E6E6E6",5),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
           rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",8),"white",rep("#E6E6E6",5),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
           rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",14),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
           rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",14),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
           rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",14),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
           rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",14),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
           rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",14),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
           rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",14),rep("white",2),"#E6E6E6","white",rep("#E6E6E6",3),
           rep("white",3),rep("#E6E6E6",2),"white",rep("#E6E6E6",14),rep("white",2),rep("#E6E6E6",5),
           rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",14),rep("white",2),rep("#E6E6E6",5),
           rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",14),rep("white",2),rep("#E6E6E6",5),
           rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",14),rep("white",2),rep("#E6E6E6",5),
           rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",14),rep("white",2),rep("#E6E6E6",5),
           rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",14),rep("white",2),rep("#E6E6E6",5),
           rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",14),rep("white",2),rep("#E6E6E6",5),
           rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",14),rep("white",2),rep("#E6E6E6",5),
           rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",14),rep("white",2),rep("#E6E6E6",5),
           rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",15),"white",rep("#E6E6E6",5),
           rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",15),"white",rep("#E6E6E6",5),
           rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",21),
           rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",21),
           rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",21),
           rep("white",2),rep("#E6E6E6",3),"white",rep("#E6E6E6",21),
           rep("white",2),rep("#E6E6E6",25),
           rep("white",2),rep("#E6E6E6",25),
           "#E6E6E6","white",rep("#E6E6E6",25),
           "#E6E6E6","white",rep("#E6E6E6",25),
           "#E6E6E6","white",rep("#E6E6E6",25),
           "#E6E6E6","white",rep("#E6E6E6",25))
  
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
