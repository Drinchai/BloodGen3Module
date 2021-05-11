library(BloodGen3Module)

load("./R/sysdata.rda")
Module_listGen3 = Module_listGen3

Gen3_ann = Gen3_ann

color = color

usethis::use_data(Module_listGen3, Gen3_ann, color, internal = TRUE)
