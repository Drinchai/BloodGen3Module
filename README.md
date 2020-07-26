# BloodGen3Module: Modular Repertoire Analysis and Visualization
***
The **BloodGen3Module** package provides functions for R user performing module repertoire analyses and generating fingerprint representations.

Steps involved in module repertoire analysis and visualization include: 

1.	Annotating the gene expression data matrix with module membership information. 
2.	Running statistical tests to determine for each module the proportion of constitutive genes which are differentially expressed.
3.	Expressing results “at the module level” as percent of genes increased or decreased. 
4.	Visualizing results from group comparison as a fingerprint grid and results from individual sample comparisons as a fingerprint heatmap.



## Installation
It is recommended to use the ```install_github``` function from the ```devtools``` package in order to install the R package.

```{r Package installation}
install.packages("devtools")
devtools::install_github("Drinchai/BloodGen3Module")
```

## Usage
```{r setup, warning=FALSE,message=FALSE}
# Load library

library(BloodGen3Module)
library(ComplexHeatmap)
library(ggplot2)
library(reshape2)
library(gtools)
library(matrixStats)
library(preprocessCore)
library(randomcoloR)
```

## Arguments
```{r argument}
data.matrix      Normalized expression data (not Log2 transformed)
sample_info      A table of sample annotation 
FC               Foldchange cut off to consider th eabundance of a given transcript to be increased or decreased compared to a reference group (Ref_group)
DIFF             Difference cut off to consider th eabundance of a given transcript to be increased or decreased compared to a reference group (Ref_group)
pval             p-value cut off or False discovery rate when FDR = FALSE
FDR              False discovery rate cut off
Group_column     Name of the columns for the groups used for the analysis
Ref_group        Reference group or samples that considered as control 
Group_df         Output table generated after running the 'Groupcomparison' function 
Individual_df    Output table generated after running the 'Individualcomparison' function
cutoff           Sets the percentage cut off used for fingerprint visualization, range of acceptable values from 0 to 100
rowSplit         Splits row of heatmap by each aggregate 
show_ref_group	 Plot reference group in the heatmap
Aggregate        Selects specific module aggregates for heatmap fingerprint plot
filename         Give a name for saving file
height           Sets height dimension for the heatmap plot
width            Sets width dimension for the heatmap plot
```



## Input
To perform the modular repertoire analysis, the R package simply requires a sample annotation table and a normalized expression data matrix
For illustrative purposes sample input files can be downloaded here; https://github.com/Drinchai/BloodGen3Module/tree/master/data.

```{r raw data and annotaion preparation}
#Load expression data
data.matrix = data_exp

#Sample annotation file
head(sample_ann)

```

## Group comparison analysis 
The **Groupcomparison** function will perform group comparison analyses and the results are expressed “at the module level” as percent of genes increased or decreased.  
- Expression matrix and sample annotation files are required to perform this analysis. 
- The names of the columns for the conditions used in the analysis must be specified.

```{r group comparison analysis,warning=FALSE}
Group_df <- Groupcomparison(data.matrix,
                            sample_info = sample_ann,
                            FC = 1.5,
                            pval = 0.1 ,
                            FDR = TRUE,
                            Group_column = "Group_test",
                            Ref_group = "Control")
```

## Fingerprint grid visualization 
The **gridplot** function will generate a grid plot as a pdf file. Specific working directory for the analysis need to be specified for saving the file. The result of the plot should be return in the same working directory.

The default cut off for visualization is set at 15%, it can changed to any value between 0-100%. 


```{r grid visulization}

gridplot(Group_df, 
         cutoff = 15, 
         Ref_group = "Control",
         filename= "Group_comparison_")

```
### Grid visualization
![Sepsis vs Control](https://github.com/Drinchai/BloodGen3Module/blob/master/2020%20July26%20Group%20comparison_Fig1.png)

## Individual single sample analysis 
The **Individualcomparison** function will perform individual sample comparison analysis in reference to a control sample or group of samples, with the results are expressed “at the module level” as percent of genes increased or decreased. 

- Expression matrix and sample annotation file are required in order to perform this analysis. 
- The names of the columns for the conditions used in the analysis must be specified
- The default cutoff is set at FC =1.5 and DIFF =10 


```{r individual single sample analysis, warning=FALSE}

Individual_df = Individualcomparison(data.matrix,
                                     sample_info = sample_ann,
                                     FC = 1.5,
                                     DIFF = 10,
                                     Group_column = "Group_test",
                                     Ref_group = "Control")
```

## Individual fingerprint visualization 
The **fingerprintplot** function will generate fingerprint heatmap plots as a pdf file. The file will be saved in the working directory specified for the analysis.

The default cut off for visualization is set at 15%, it can changed to any value between 0-100%.  
 

```{r fingerprint visualization, warning=FALSE}

fingerprintplot(Individual_df,
                sample_info = sample_ann,
                cutoff = 15,
                rowSplit= TRUE ,
                Group_column= "Group_test",
                show_ref_group = FALSE, 
                Ref_group =  "Control",
                Aggregate = "A28",
                filename = "Gen3_Individual_plot",
                height = NULL,
                width = NULL)

```
### Heatmap fingerprint visualization
![Individual_plot](https://github.com/Drinchai/BloodGen3Module/blob/master/2020%20July26%20Individual%20comparison_Fig2.png)


## Publication
A manuscript is currently under consideration for publication, in order to cite the work currently please refer to the bioRxiv preprint:
https://www.biorxiv.org/content/10.1101/2020.07.16.205963v1
