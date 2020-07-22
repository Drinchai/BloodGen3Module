Module_listGen3 <- read.csv("data-raw/Module_listGen3.csv",
                     stringsAsFactor = FALSE)

rownames(Module_listGen3) = Module_listGen3$Module_gene

Module_listGen3 = data.frame(Module_listGen3)

usethis::use_data_raw(Module_listGen3, overwrite = TRUE)
