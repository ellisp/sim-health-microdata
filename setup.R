library(tidyverse)
library(survey)
library(boot)
library(GGally)
library(broom)
library(extrafont)
library(frs) # for svg_png()

all_scripts <- list.files("R", pattern = ".R", full.names = TRUE)
for(f in all_scripts){
  source(f)
}
rm(f, all_scripts)

theme_set(theme_minimal(base_family = "Roboto") +
            theme(plot.caption = element_text(colour = "grey50"),
                  strip.text = element_text(size = rel(1), face = "bold"),
                  plot.title = element_text(family = "Sarala")))
  