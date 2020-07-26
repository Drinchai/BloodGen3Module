#' Individual fingerprint visualization
#'
#'The fingerprintplot function will generate fingerprint heatmap plots as a pdf file. The file will be saved in the working directory specified for the analysis.
#'The default cut off for visualization is set at 15%, it can changed to any value between 0-100%.
#'
#' @param Individual_df 		Output table generated after running the 'Individualcomparison' function
#' @param sample_info		 	Sample annotation table
#' @param cutoff 			Sets the percentage cut off used for fingerprint visualization, range of acceptable values from 0 to 100
#' @param rowSplit		 Splits row of heatmap by each aggregate
#' @param Ref_group 		Reference group or samples that considered as control
#' @param Group_column		 Name of the columns for the groups used for the analysis
#' @param show_ref_group	 Plot reference group in the heatmap
#' @param Aggregate		 Selects specific module aggregates for heatmap fingerprint plot
#' @param filename			 Give a name for the file
#' @param height			 Sets height dimension for the heatmap plot
#' @param width		 	Sets width dimension for the heatmap plot
#' @return A heatmap of % of module response in each single sample
#' @examples
#'fingerprintplot(Individual_df, cutoff = 15, rowSplit= TRUE , Group_column= "Group_test", Ref_group =  "Control", Aggregate = c("A28), height = NULL, width = NULL)
#'#' @author
#' Darawan Rinchai <drinchai@gmail.com>
#' @export
fingerprintplot = function(Individual_df,
                           sample_info = sample_info,
                           cutoff = NULL,
                           rowSplit= TRUE ,
                           Ref_group=NULL,
                           show_ref_group = FALSE,
                           Group_column= NULL,
                           Aggregate = NULL,
                           filename = NULL,
                           height = NULL,
                           width = NULL){
  #Load module annotation
  Sum.mod.sin = Individual_df
  Sum.mod.sin = Sum.mod.sin[rownames(Gen3_ann),]
  rownames(Sum.mod.sin) == rownames(Gen3_ann)

  rownames(Sum.mod.sin) <- paste(Gen3_ann$Module, Gen3_ann$Function, sep = ".")

  ############################################################
  ##################### MODULES GEN3 and MODULE WITH FUNCTION DEFINED #######################################
  #modules with function deffined

  Module.list <- unique(Gen3_ann[,c("Module","Function")])                                                             # creat new dataframe from Module
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

  ###remove control sample from plot
  if (show_ref_group == FALSE) {
    sample_info = sample_info[!sample_info[, Group_column]== Ref_group,]
  }


  df_plot = df_plot[,rownames(sample_info)]

  colnames(df_plot)==rownames(sample_info)

  n.group = length(unique(sample_info[, Group_column]))


  n <- n.group
  palette <- distinctColorPalette(n)

  my.pattle = palette
  names(my.pattle) = unique(sample_info[, Group_column])



  col_fun = circlize::colorRamp2(c(-100,0,100), c("blue", "white", "red"))

  ##prepare annotation table
  ####################
  Gen3_ann$Module_func = paste(Gen3_ann$Module, Gen3_ann$Function,sep = ".")


  if (is.null(Aggregate)) {
    anno_table = Gen3_ann[Gen3_ann$Module_func%in%rownames(df_plot),]
  }
  else {
    anno_table = Gen3_ann[grep(Gen3_ann$Cluster,pattern = Aggregate),]
  }

  rownames(anno_table) == anno_table$Module
  rownames(anno_table) = anno_table$Module_func


  df_plot = df_plot[rownames(anno_table),]


  plate_color = as.character(anno_table$Module_color)
  names(plate_color)=anno_table$Function

  left_ha = rowAnnotation(df = data.frame(Module = anno_table$Function),
                          show_annotation_name = FALSE,simple_anno_size = unit(0.3, "cm"),
                          col = list(Module = plate_color))

  ha_column = HeatmapAnnotation(df = data.frame(Group = sample_info[, Group_column]),
                                show_annotation_name = FALSE, simple_anno_size = unit(0.3, "cm"),
                                col = list(Group = my.pattle))
  if (rowSplit == TRUE) {
    rowSplit = anno_table$Cluster
  }
  else {
    rowSplit = NULL
  }

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

  pdf(file = paste0(filename, "_", Aggregate,".pdf"), height = height, width = width)
  ht=Heatmap(df_plot,
             cluster_rows = TRUE,
             cluster_columns = T,
             height = unit(2.1, "mm")*nrow(df_plot),
             width  = unit(2.1, "mm")*ncol(df_plot),
             rect_gp = gpar(type = "none"),
             row_split = rowSplit,
             top_annotation = ha_column,
             left_annotation = left_ha,
             name = "% Response",
             row_names_max_width = unit(10,"in"),
             row_title_gp = gpar(fontsize = 10),
             row_title_rot = 0,
             column_names_gp = gpar(fontsize = 4),
             row_names_gp = gpar(fontsize = 5),
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.circle(x = x, y = y, r = unit(0.905, "mm") ,gp = gpar(fill = col_fun(df_plot[i, j]), col = NA))
             }
  )
  draw(ht,heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(2, 20, 2, 2), "mm"))

  dev.off()

}
