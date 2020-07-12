source("setup.R")

#================Create some marginal totals===============


#---------------Set up variables------------------
continents <- c("Eastasia", "Oceania", "Eurasia")
races <- c("Klingon", "Vulcan", "Human", "Missing")
sexes <- c("Male", "Female", "Other", "Missing")
smoking <- c("Smoking", "Non-smoking", "Missing")
weight <- c("Underweight", "Healthy weight", "Overweight", "Obese", "Missing")

#-------------Race x continent------------------
rc_tot <- expand.grid(race = races, continent = continents) %>%
  mutate(Freq = c(100, 200, 30000, 500,
                  5000, 1000, 2000, 500,
                  2000, 10000,3000, 1000))

c_tot <- rc_tot %>%
  group_by(continent) %>%
  summarise(continent_total = sum(Freq))

r_tot <- rc_tot %>%
  group_by(race) %>%
  summarise(race_total = sum(Freq))

#-------------------Race x sex---------------------
rx_tot <- expand.grid(race = races, sex = sexes) %>%
  left_join(r_tot, by = "race") %>%
  mutate(Freq = case_when(
    sex == "Male" & race == "Klingon" ~ 0.40 * race_total,
    sex == "Female" & race == "Klingon" ~ 0.56 * race_total,
    sex == "Male"                     ~ 0.47 * race_total,
    sex == "Female"                   ~ 0.49 * race_total,
    sex == "Missing"                  ~ 0.03 * race_total,
    sex == "Other"                    ~ 0.01 * race_total
  )) %>%
  select(-race_total)

stopifnot(sum(rc_tot$Freq) == sum(rx_tot$Freq))

x_tot <- rx_tot %>%
  group_by(sex) %>%
  summarise(sex_total = sum(Freq))

#-------------------Sex x smoking----------------

xs_tot <- expand.grid(sex = sexes, smoking = smoking) %>%
  left_join(x_tot, by = "sex") %>%
  mutate(Freq = case_when(
    sex == "Male" & smoking == "Smoking" ~ 0.2 * sex_total,
    sex == "Male" & smoking == "Non-smoking" ~ 0.7 * sex_total,
    sex == "Male" & smoking == "Missing" ~ 0.1 * sex_total,
    sex == "Female" & smoking == "Smoking" ~ 0.15 * sex_total,
    sex == "Female" & smoking == "Non-smoking" ~ 0.72 * sex_total,
    sex == "Female" & smoking == "Missing" ~ 0.13 * sex_total,
    sex %in% c("Other", "Missing") ~ sex_total / 3
  )) %>%
  select(-sex_total)

stopifnot(sum(rc_tot$Freq) == sum(xs_tot$Freq))

s_tot <- xs_tot %>%
  group_by(smoking) %>%
  summarise(smoking_total = sum(Freq))

#---------------------Smoking x weight---------------------

sw_tot <- expand.grid(smoking = smoking, weight = weight) %>%
  left_join(s_tot, by = "smoking") %>%
  mutate(Freq = case_when(
    smoking == "Smoking" & weight == "Obese" ~ 0.2 * smoking_total,
    weight == "Obese" ~ 0.1 * smoking_total,
    weight == "Healthy weight" ~ 0.3 * smoking_total,
    weight == "Overweight" ~ 0.3 * smoking_total,
    weight == "Missing" ~ 0.05 * smoking_total,
    smoking == "Smoking" & weight == "Underweight" ~ 0.15 * smoking_total,
    TRUE ~ 0.25 * smoking_total
  )) %>%
  select(-smoking_total)


stopifnot(sum(rc_tot$Freq) == sum(sw_tot$Freq))

w_tot <- sw_tot %>%
  group_by(weight) %>%
  summarise(weight_total = sum(Freq))

#====================Create microdata of demographic variables that match the marginal totals========

population <- expand.grid(
  continent = continents,
  race = races,
  sex = sexes,
  smoking = smoking,
  weight = weight
) %>%
  as_tibble()

sd <- svydesign(~1, data = population)
sdr <- rake(sd, 
     sample = list(~race + continent,
                   ~race + sex,
                   ~sex + smoking,
                   ~smoking + weight),
     population = list(rc_tot,
                       rx_tot,
                       xs_tot,
                       sw_tot))

population$prob <- weights(sdr) / sum(weights(sdr))

set.seed(123)
the_sample <- sample_n(population, 60000, weight = prob, replace = TRUE) %>%
  select(-prob)


table(the_sample$race)
table(the_sample$continent)
table(the_sample$sex)
table(the_sample$smoking)
table(the_sample$weight)

#=====================Model to get the treatment and the result===========

set.seed(42)
the_sample <- the_sample %>%
  mutate(unobserved1 = rnorm(n()),
         unobserved2 = rnorm(n())) %>%
  mutate(propensity_treatment =
            0.3 * (continent == "Eastasia") +
            0.1 * (race == "Human") +
            0.1 * (smoking == "Non-smoking") +
            0.2 * (weight == "Obese") +
           -0.2 * (weight == "Healthy weight") +
           unobserved1 / 10 +
           -2,
         received_treatment = rbinom(n = n(), size = 1,
                                     prob = inv.logit(propensity_treatment))) %>%
  mutate(propensity_comorbidity = 
           1.1 * (smoking == "Smoking") +
           0.3 * (weight == "Obese") +
           0.1 * (sex == "Male"),
         comorbidity = rbinom(n = n(), size = 1, 
                              prob = inv.logit(propensity_comorbidity))) %>%
  mutate(propensity_death =
            0.3 * (continent == "Oceania") +
           -0.1 * (continent == "Eurasia") +
            0.2 * (race == "Klingon") +
           -0.2 * (race == "Vulcan") +
            0.1 * (smoking == "Missing") +
            0.19 * (sex == "Male") +
            0.1 * (weight == "Obese") +
           -0.2 * (weight == "Healthy weight") +
           -0.05 * (weight == "Under weight") +
            unobserved1 / 12 +
            unobserved2 / 8 +
            0.8  * comorbidity +
            0.2 * received_treatment +
           -3) %>%
  mutate(death = rbinom(n = n(), size = 1,
                        prob = inv.logit(propensity_death))) %>%
  mutate(weight  = fct_relevel(weight, "Healthy weight"),
         smoking = fct_relevel(smoking, "Non-smoking"),
         race    = fct_relevel(race, "Human"),
         sex     = fct_relevel(sex, "Male"))
           

the_sample %>%
  sample_n(size = 2000) %>%
  select(unobserved1, unobserved2, propensity_treatment, propensity_death) %>%
  ggpairs()


the_sample %>%
  select(smoking, weight, received_treatment, death, sex) %>%
  sample_n(size = 2000) %>%
  ggpairs(mapping = aes(fill = sex, colour = sex))

#===============fitting a model============

#---------------Model 1 - with co-morbidity------------

mod1 <- glm(death ~ continent + weight + race + sex + smoking + weight + comorbidity + received_treatment, 
            family = "binomial",
            data = the_sample)

summary(mod1)
anova(mod1, test = "Chi")
ci1 <- confint(mod1)

p1 <- make_plot(ci1)
p1
#---------------Model 2 - without co-morbidity------------

mod2 <- glm(death ~ continent + weight + race + sex + smoking + weight + received_treatment, 
            family = "binomial",
            data = the_sample)

ci2 <- confint(mod2)

p2 <- make_plot(ci2)
p2
