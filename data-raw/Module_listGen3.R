Module_listGen3 <- read.csv("data-raw/Module_listGen3.csv",
                     stringsAsFactor = FALSE)

rownames(Module_listGen3) = Module_listGen3$Module_gene

usethis::use_data(Module_listGen3, overwrite = TRUE)
