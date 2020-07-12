library(tidyverse)
library(survey)
library(boot)
library(GGally)
library(broom)

all_scripts <- list.files("R", pattern = ".R", full.names = TRUE)
for(f in all_scripts){
  source(f)
}
rm(f, all_scripts)
