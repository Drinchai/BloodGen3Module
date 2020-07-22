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
  Module_ann$Module_color = Module.platte


  Sum.mod.sin = Individual_df
  Sum.mod.sin = Sum.mod.sin[rownames(Module_ann),]
  rownames(Sum.mod.sin) == rownames(Module_ann)

  rownames(Sum.mod.sin) <- paste(Module_ann$Module, Module_ann$Function, sep = ".")

  Sum.mod.sin.comp <- Sum.mod.sin[apply(Sum.mod.sin[,], 1, function(x) !all(x==0)),]                 # Rows sum=0

  ############################################################
  ##################### MODULES GEN3 and MODULE WITH FUNCTION DEFINED #######################################
  #modules with function deffined

  Module.list <- unique(Module_ann[,c("Module","Function")])                                                             # creat new dataframe from Module
  Module.list$Modules <- paste(Module.list$Module, Module.list$Function, sep = ".")
  rownames(Module.list) <- Module.list$Modules

  mod.with.function <- Module.list$Modules[which(Module.list$Function!="TBD")]                         # select module that have only function
  Sum.mod.sin.comp.withF <- Sum.mod.sin.comp[rownames(Sum.mod.sin.comp) %in% mod.with.function,]       # selected only modules that have function in this dataset

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
  rownames(anno_table) == anno_table$Module_func


  df_plot = df_plot[rownames(anno_table),]


  plate_color = as.character(anno_table$Module_color)
  names(plate_color)=anno_table$Function

  left_ha = rowAnnotation(df = data.frame(Module = anno_table$Function),
                          show_annotation_name = FALSE,simple_anno_size = unit(0.3, "cm"),
                          col = list(Module = plate_color))

  ha_column = HeatmapAnnotation(df = data.frame(Group = sample_info$Group_test),
                                show_annotation_name = TRUE,
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


