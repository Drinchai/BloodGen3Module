#' This is a function for individual fingerprint  plot
#'
#' This function will perform fingerprint heatmap in pdf file.
#' We assumes that you setup the specific working directory for your analysis. The result of the plot should be in your working folder.
#' The result of the plot should be in your working folder.
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
fingerprintplot = function(Individual_df, cutoff = NULL, Aggregate = NULL, height = NULL, width = NULL){
  #Load module annotation
  sample_info = sample.info
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


  Sum.mod.sin = Individual_df
  Sum.mod.sin = Sum.mod.sin[rownames(Module_ann),]
  rownames(Sum.mod.sin) == rownames(Module_ann)

  rownames(Sum.mod.sin) <- paste(Module_ann$Module, Module_ann$Function, sep = ".")

  ############################################################
  ##################### MODULES GEN3 and MODULE WITH FUNCTION DEFINED #######################################
  #modules with function deffined

  Module.list <- unique(Module_ann[,c("Module","Function")])                                                             # creat new dataframe from Module
  Module.list$Modules <- paste(Module.list$Module, Module.list$Function, sep = ".")
  rownames(Module.list) <- Module.list$Modules

  mod.with.function <- Module.list$Modules[which(Module.list$Function!="TBD")]                         # select module that have only function
  Sum.mod.sin.comp.withF <- Sum.mod.sin[rownames(Sum.mod.sin) %in% mod.with.function,]       # selected only modules that have function in this dataset

  ####################################################################################
  ####### DOT Heatmap by complexHeatmap ####

  df_plot = Sum.mod.sin.comp.withF

  ########## An example of DISPLAY DATA > 15 %
  if (is.null(cutoff)) {
    cutoff = 15
  }
  else {
    cutoff = as.numeric(cutoff)
  }

  df_plot[abs(df_plot) < cutoff] <- 0

  df_plot = df_plot[,rownames(sample_info)]

  colnames(df_plot)==rownames(sample_info)

  n.group = length(unique(sample_info$Group_test))

  library(randomcoloR)
  n <- n.group
  palette <- distinctColorPalette(n)

  my.pattle = palette
  names(my.pattle) = unique(sample_info$Group_test)



  col_fun = circlize::colorRamp2(c(-100,0,100), c("blue", "white", "red"))

  ##prepare annotation table
  ####################
  Module_ann$Module_func = paste(Module_ann$Module, Module_ann$Function,sep = ".")


  if (is.null(Aggregate)) {
    anno_table = Module_ann[Module_ann$Module_func%in%rownames(df_plot),]
  }
  else {
    anno_table = Module_ann[grep(Module_ann$Cluster,pattern = Aggregate),]
  }

  rownames(anno_table) == anno_table$Module
  rownames(anno_table) = anno_table$Module_func


  df_plot = df_plot[rownames(anno_table),]


  plate_color = as.character(anno_table$Module_color)
  names(plate_color)=anno_table$Function

  left_ha = rowAnnotation(df = data.frame(Module = anno_table$Function),
                          show_annotation_name = FALSE,simple_anno_size = unit(0.3, "cm"),
                          col = list(Module = plate_color))

  ha_column = HeatmapAnnotation(df = data.frame(Group = sample_info$Group_test),
                                show_annotation_name = FALSE, simple_anno_size = unit(0.3, "cm"),
                                col = list(Group = my.pattle))

  #DOT HEATMAP
  if (is.null(height)) {
    height = 28
  }
  else {
    height = as.numeric(height)
  }

  if (is.null(width)) {
    width = 17
  }
  else {
    width = as.numeric(height)
  }

  pdf(paste0("Gen3_Individual_analysis.pdf"), height = height, width = width)
  ht=Heatmap(df_plot,
             cluster_rows = TRUE,
             cluster_columns = T,
             height = unit(2.1, "mm")*nrow(df_plot),
             width  = unit(2.1, "mm")*ncol(df_plot),
             rect_gp = gpar(type = "none"),
             top_annotation = ha_column,
             left_annotation = left_ha,
             name = "% Response",
             row_names_max_width = unit(10,"in"),
             row_title_gp = gpar(fontsize = 0.1),
             column_names_gp = gpar(fontsize = 4),
             row_names_gp = gpar(fontsize = 5),
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.circle(x = x, y = y, r = unit(0.905, "mm") ,gp = gpar(fill = col_fun(df_plot[i, j]), col = NA))
             }
  )
  draw(ht,heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(2, 20, 2, 2), "mm"))

  dev.off()

}
