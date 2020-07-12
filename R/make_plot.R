make_plot<- function(ci){
  ci %>% 
    tidy() %>%
    rename(var = .rownames,
           upper = X97.5..,
           lower = X2.5..) %>%
    mutate(metavar = case_when(
      grepl("continent", var) ~ "Continent compared to 'Eastasia'",
      grepl("weight", var) ~ "Weight compared to 'healthy weight'",
      grepl("race", var) ~ "Species compared to 'Human'",
      grepl("sex", var) ~ "Sex compared to 'Female'",
      grepl("smoking", var) ~ "Smoking status",
      var == "received_treatment" ~ "Treatment",
      var == "comorbidity" ~ "Existing co-morbidity"
    )) %>%
    mutate(var = gsub("^continent", "", var),
           var = gsub("^weight", "", var),
           var = gsub("^race", "", var),
           var = gsub("^sex", "", var),
           var = gsub("^smoking", "", var),
           var = gsub("received_treatment", "Received treatment", var),
           var = gsub("comorbidity", "Existing co-morbidity", var)) %>%
    filter(!var %in% c("(Intercept)", "Missing")) %>%
    mutate(mid = (upper + lower) / 2,
           metavar = fct_reorder(metavar, -mid),
           var = fct_reorder2(var, metavar, mid, .fun = function(x, y){mean(as.numeric(x) - as.numeric(y))})) %>%
    mutate(metavar = fct_relevel(metavar, "Existing co-morbidity", after = Inf),
           var = fct_relevel(var, "Existing co-morbidity")) %>%
    mutate_if(is.numeric, exp) %>%
    ggplot(aes(y = var, yend = var, x = lower, xend = upper, colour = metavar)) +
    geom_vline(xintercept = 1) +
    geom_segment(size = 2) +
    scale_x_log10(breaks = seq(from = 0.75, to = 3, by = 0.25)) +
    scale_colour_brewer(palette = "Set2") +
    labs(x ="95% confidence interval of odds ratio of dying",
         y = "",
         colour = "Type of variable",
         title = "Simulating plausible microdata for clinical studies",
         subtitle = str_wrap("Demonstration of how to create a plausible data for an observational study given marginal totals, 
       models (featuring unobserved variables) for propensity to receive treatment, co-morbidity, and impact on target variable.", 95),
         caption = "http://freerangestats.info")
}