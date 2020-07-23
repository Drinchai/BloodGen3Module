

Gen3_ann <- read.csv("data-raw/Gen3_ann.csv",
                                  stringsAsFactor = FALSE)

rownames(Gen3_ann) = Gen3_ann$Module

#devtools::use_data_raw(Gen3_ann, overwrite = TRUE)
