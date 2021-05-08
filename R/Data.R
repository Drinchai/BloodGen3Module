#' Functional annotation for visulization
#'
#'A dataframe contain module identifer, functional annotation and color in each specific module
#'
#' @format A data.frame with 5 variables:
#' \code{Modules: Module identifier},
#' \code{Function: Functional annotation},
#' \code{Position: position on fingerprint grid plot},
#' \code{Module_color: specific color of each module for visulization},
#' \code{Cluster: Module cluster membership}
"Gen3_ann"

#' Module identifer and list membership in each module
#'
#' A dataframe contain gene member of 3rd generantion of blood module repertoire construction
#'
#' @format A data.frame with 14168 rows by 5 variables:
#' \code{Module: Module identifier},
#' \code{Gene symbol: gene membership},
#' \code{Module_gene: gene specific module membership},
#' \code{Function: Functional annotation},
#' \code{position: position on fingerprint grid plot}
"Module_listGen3"

#' Color for fingerprint visulization
#'
#' Character vector of color for fingerprint grid plot
#'
#' @format A vector of 1134 character
"color"
