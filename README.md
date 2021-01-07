# BloodGen3Module: Modular Repertoire Analysis and Visualization
***
The **BloodGen3Module** package provides functions for R users to perform module repertoire analyses and generate fingerprint representations.

The steps involved in module repertoire analysis and visualization include:

1.	Annotating the gene expression data matrix with module membership information.
2.	Running statistical tests to determine the proportion of constitutive genes that are differentially expressed for each module.
3.	Expressing results “at the module level” as the percentage of genes that are increased or decreased.
4.	Visualizing results from group comparisons as a fingerprint grid and results from individual sample comparisons as a fingerprint heatmap.


## Installation
It is recommended to use the ```install_github``` function from the ```devtools``` package in order to install the R package.

```{r Package installation}

#Installation from Github

install.packages("devtools")
devtools::install_github("Drinchai/BloodGen3Module")

```

## Usage
```{r setup, warning=FALSE,message=FALSE}
# Load library

library(BloodGen3Module)

```

## Arguments
```{r argument}
data.matrix      Matrix of normalized expression data (that should not be Log2 transformed). Genes should be arranged as rows and Sample ID as columns.Row names are required to be valid Gene Symbols.
sample_info      Dataframe with sample annotation. Sample_info dataframe requires two columns: 1) a column specifying Sample ID (exactly matching Sample ID of data.matrix) and 2) a column specifying group names
FC               Numeric value specifying the fold change cut off that will be applied to define increase or decrease of a given transcript compared to the reference group (Ref_group)
DIFF             Numeric value specifying the absolute difference cut off that will be applied to define increase or decrease of a given transcript compared to the reference group (Ref_group)
pval             Numeric value specifying the p-value cut off or False discovery rate when FDR = TRUE
FDR              Logical operator (TRUE/FALSE) to specify whether False discovery rate cut off (using BH-method) should be used.
Group_column     Character vector identical to the column name from sample_info dataframe that specifies group annotation used for the analysis
Test_group       Character vector specifying values within the group column (Group_column) that will be used as Test group (samples considered as cases or “intervention” group).
Ref_group        Character vector specifying values within the group column (Group_column) that will be used as Reference group (samples considered as control).
Group_df         Dataframe with output generated after running the 'Groupcomparison' function 
Group_limma      Dataframe with output generated after running the 'Groupcomparisonlimma' function
Individual_df    Dataframe with output generated after running the 'Individualcomparison' function
cutoff           Numeric value specifying the percentage cut off used for fingerprint visualization (acceptable values range from 0 to 100).
rowSplit         Logical operator (TRUE/FALSE) to indicate if the rows of the heatmap should be split by aggregate 
show_ref_group	 Logical operator (TRUE/FALSE) to indicate if a reference group should be plotted on the heatmap
Aggregate        Character vector specifying the name of specific module aggregates on the heatmap fingerprint plot
filename         Character vector specifying the name of the saved output file
height           Sets the height of the graphics region in inches. The default values is 28
width            Sets the width of the graphics region in inches. The default values is 17
```


## Input
To perform the modular repertoire analysis, the R package requires a sample annotation table and a normalized expression data matrix. For illustrative purposes, the sample input files can be downloaded from: https://github.com/Drinchai/GSE13015.

```{r raw data and annotaion preparation}
#Load expression data
load(./data_exp.rda)
data.matrix = data_exp

#Sample annotation file
load(./sample_ann.rda)
head(sample_ann)

```

## Group comparison analysis 
The **Groupcomparison** function will perform group comparison analyses. The results are expressed “at the module level” as the percentage of genes that are increased or decreased for a given module.

- Expression matrix and sample annotation files are required to perform this analysis.
- The sample annotation file must be loaded using a specific name = "sample.info".
- The names of the columns for the conditions used in the analysis must be specified.

Using t-test
```{r group comparison analysis,warning=FALSE}
Group_df <- Groupcomparison(data.matrix,
                            sample_info = sample_ann,
                            FC = 1.5,
                            pval = 0.1 ,
                            FDR = TRUE,
                            Group_column = "Group_test",
                            Test_group = "Sepsis",
                            Ref_group = "Control")
```
Using "limma"

```{r group comparison analysis using "limma",warning=FALSE}
Group_limma <- Groupcomparisonlimma(data.matrix,
                                    sample_info = sample_ann,
                                    FC = 1.5,
                                    pval = 0.1 ,
                                    FDR = TRUE,
                                    Group_column = "Group_test",
                                    Test_group = "Sepsis",
                                    Ref_group = "Control")
```


## Fingerprint grid visualization 
The **gridplot** function will generate a grid plot as a PDF file. A specific working directory for the analysis must be specified to save the file. The result of the plot should be returned in the same working directory.


The default cut off for visualization is set at 15%; it can be changed to any value between 0 and 100%.



```{r grid visulization}

gridplot(Group_df, 
         cutoff = 15, 
         Ref_group = "Control",
         filename= "Group_comparison_")

```

### Grid visualization
![Sepsis vs Control](https://github.com/Drinchai/DC_Gen3_Module_analysis/blob/master/2020%20July26%20Group%20comparison_Fig1.png)

## Individual single sample analysis 
The **Individualcomparison** function will perform an individual sample comparison analysis in reference to a control sample or group of samples. The results are expressed “at the module level” as the percentage of genes that are increased or decreased.

- Expression matrix and sample annotation files are required to perform this analysis.
- The sample annotation file must be loaded using a specific name = "sample.info".
- The names of the columns for the conditions used in the analysis must be specified.
- The default cut off is set at fold change (FC) =1.5 and absolute difference (DIFF) =10.



```{r individual single sample analysis, warning=FALSE}

Individual_df = Individualcomparison(data.matrix,
                                     sample_info = sample_ann,
                                     FC = 1.5,
                                     DIFF = 10,
                                     Group_column = "Group_test",
                                     Ref_group = "Control")
```

## Individual fingerprint visualization 
The **fingerprintplot** ffunction will generate fingerprint heatmap plots as a PDF file. The file will be saved in the working directory specified for the analysis.

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
![Individual_plot](https://github.com/Drinchai/DC_Gen3_Module_analysis/blob/master/2020%20July26%20Individual%20comparison_Fig2.png)


## Publication
A manuscript is currently under consideration for publication, in order to cite the work currently please refer to the bioRxiv preprint:
https://www.biorxiv.org/content/10.1101/2020.07.16.205963v1
