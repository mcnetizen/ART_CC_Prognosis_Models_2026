# ******************************************************************************
dta0 <- readRDS(here::here("list_dta_prognosis_2holdout.rds")) %>% 
  map(~..1 %>% 
        select(id, cohort, kfold, kfold_no,
               fail5, fupys5,
               sexm3, aids2, age10, cd4100, rna, undetect, 
               cd4_vf, cd4_v, rna_v, age, age_f))

# ******************************************************************************
# T0 ----
# ******************************************************************************
dta0t0 <- dta0[[1]]
pipelinet0 <- readRDS(here::here("results_model_t0.rds"))
modelt0 = pipelinet0$model
modelt0@nobs
nrow(dta0t0)

# ..............................................................................
predtime <- seq(0, 5, 5/50)
predtime <- predtime[-1]
dta_plot_survt0 <- vector("list", 50)
for (j in predtime) {
  newtestdata   <- dta0t0 %>% mutate(fupys5 = j)
  x             <- match(j, predtime)
  predlogh      <- predict(modelt0, newtestdata, type = "loghazard")
  preddata      <- predict(modelt0, newtestdata, type = "surv", se.fit = T)
  
  dta_plot_survt0[[x]] <- tibble(kfold    = dta0t0$kfold, 
                                 cd4_vf   = dta0t0$cd4_vf,
                                 aids     = dta0t0$aids2,
                                 age_f    = dta0t0$age_f,
                                 time     = newtestdata$fupys5,
                                 id       = dta0t0$id,
                                 risk     = predlogh,
                                 predsurv = preddata$Estimate, 
                                 lower    = preddata$lower, 
                                 upper    = preddata$upper)
  rm(x, preddata, newtestdata)
}

#...............................................................................
predt0cd <- bind_rows(dta_plot_survt0) %>% 
  select(cd4_vf, time, predsurv, lower, upper) %>% 
  distinct(cd4_vf, .keep_all = T) %>% 
  mutate(time = 0, predsurv = 1, lower = 1, upper = 1) %>% 
  bind_rows(
    bind_rows(dta_plot_survt0)  %>% 
      select(cd4_vf, time, predsurv, lower, upper) %>% 
      group_by(cd4_vf, time) %>% 
      mutate(across(predsurv:upper, ~mean(.))) %>% 
      filter(row_number() == 1) %>% ungroup()) %>% 
  arrange(cd4_vf, time) %>% 
  rename(surv = predsurv, value = cd4_vf) %>% 
  mutate(value = fct_rev(value))

predt0aids <- bind_rows(dta_plot_survt0) %>% 
  select(aids, time, predsurv, lower, upper) %>% 
  distinct(aids, .keep_all = T) %>% 
  mutate(time = 0, predsurv = 1, lower = 1, upper = 1) %>% 
  bind_rows(
    bind_rows(dta_plot_survt0)  %>% 
      select(aids, time, predsurv, lower, upper) %>% 
      group_by(aids, time) %>% 
      mutate(across(predsurv:upper, ~mean(.))) %>% 
      filter(row_number() == 1) %>% ungroup()) %>% 
  arrange(aids, time) %>% 
  rename(surv = predsurv, value = aids) %>% 
  mutate(value = fct_recode(value, Never = "never", Mild = "mild", 
                            "Severe" = "mod_sev"))

predt0age <- bind_rows(dta_plot_survt0) %>% 
  select(age_f, time, predsurv, lower, upper) %>% 
  distinct(age_f, .keep_all = T) %>% 
  mutate(time = 0, predsurv = 1, lower = 1, upper = 1) %>% 
  bind_rows(
    bind_rows(dta_plot_survt0)  %>% 
      select(age_f, time, predsurv, lower, upper) %>% 
      group_by(age_f, time) %>% 
      mutate(across(predsurv:upper, ~mean(.))) %>% 
      filter(row_number() == 1) %>% ungroup()) %>% 
  arrange(age_f, time) %>% 
  rename(surv = predsurv, value = age_f) 

predt0all <- tibble(
  time = 0, predsurv = 1, lower = 1, upper = 1) %>% 
  bind_rows(
    bind_rows(dta_plot_survt0)  %>% 
      select(time, predsurv, lower, upper) %>% 
      group_by(time) %>% 
      mutate(across(predsurv:upper, ~mean(.))) %>% 
      filter(row_number() == 1) %>% ungroup()) %>% 
  arrange(time) %>% 
  rename(surv = predsurv)

predt0 <- bind_rows(ALL = predt0all, CD4 = predt0cd, AIDS = predt0aids, AGE = predt0age,
                    .id = "vbl") %>% 
  mutate(dta = "ART start") %>% 
  select(dta, value, everything())

rm(predt0cd, predt0aids, predt0age, dta_plot_survt0)

#...............................................................................
obst0cd <- survfit(Surv(fupys5, fail5) ~ cd4_vf, dta0t0) %>% 
  summary(., times = seq(0, 5, 5/50)) %>% 
  keep_at(c("time", "surv", "lower", "upper", "strata")) %>% 
  as_tibble() %>% 
  rename(value = strata) %>% 
  mutate(value = str_replace(value, "cd4_vf=", ""), 
         value = factor(value, levels = c('[500,Inf)', '[350,500)', '[200,350)', 
                                          '[100,200)', '[50,100)', '[0,50)')))

obst0aids <- survfit(Surv(fupys5, fail5) ~ aids2, dta0t0) %>% 
  summary(., times = seq(0, 5, 5/50)) %>% 
  keep_at(c("time", "surv", "lower", "upper", "strata")) %>% 
  as_tibble() %>% 
  rename(value = strata) %>% 
  mutate(value = str_replace(value, "aids2=", "")) %>% 
  mutate(value = fct_recode(value, Never = "never", Mild = "mild", 
                            "Severe" = "mod_sev"))

obst0age <- survfit(Surv(fupys5, fail5) ~ age_f, dta0t0) %>% 
  summary(., times = seq(0, 5, 5/50)) %>% 
  keep_at(c("time", "surv", "lower", "upper", "strata")) %>% 
  as_tibble() %>% 
  rename(value = strata) %>% 
  mutate(value = str_replace(value, "age_f=", ""), 
         value = factor(value, levels = c('[16,30)', '[30,40)', '[40,50)', 
                                          '[50,60)', '[60,Inf)')))

obst0all <- survfit(Surv(fupys5, fail5) ~ 1, dta0t0) %>% 
  summary(., times = seq(0, 5, 5/50)) %>% 
  keep_at(c("time", "surv", "lower", "upper")) %>% 
  as_tibble() 

obst0 <- bind_rows(ALL = obst0all, CD4 = obst0cd, AIDS = obst0aids, AGE = obst0age,
                   .id = "vbl") %>% 
  mutate(dta = "ART start") %>% 
  select(dta, vbl, value, everything())

rm(obst0cd, obst0aids, obst0age)

#...............................................................................
t0 <- bind_rows(Predicted = predt0, Observed = obst0, .id = "data")
rm(predt0, obst0)
head(t0)

rm(dta0t0, modelt0, pipelinet0)

# ******************************************************************************
# T1 ----
# ******************************************************************************
dta0t1 <- dta0[[2]]
pipelinet1 <- readRDS(here::here("results_model_t1.rds"))
modelt1 = pipelinet1$model
modelt1@nobs
nrow(dta0t1)

# ..............................................................................
predtime <- seq(0, 5, 5/50)
predtime <- predtime[-1]
dta_plot_survt1 <- vector("list", 50)
for (j in predtime) {
  newtestdata   <- dta0t1 %>% mutate(fupys5 = j)
  x             <- match(j, predtime)
  predlogh      <- predict(modelt1, newtestdata, type = "loghazard")
  preddata      <- predict(modelt1, newtestdata, type = "surv", se.fit = T)
  
  dta_plot_survt1[[x]] <- tibble(kfold    = dta0t1$kfold, 
                                 cd4_vf   = dta0t1$cd4_vf,
                                 aids     = dta0t1$aids2,
                                 age_f    = dta0t1$age_f,
                                 time     = newtestdata$fupys5,
                                 id       = dta0t1$id,
                                 risk     = predlogh,
                                 predsurv = preddata$Estimate, 
                                 lower    = preddata$lower, 
                                 upper    = preddata$upper)
  rm(x, preddata, newtestdata)
}

#...............................................................................
predt1cd <- bind_rows(dta_plot_survt1) %>% 
  select(cd4_vf, time, predsurv, lower, upper) %>% 
  distinct(cd4_vf, .keep_all = T) %>% 
  mutate(time = 0, predsurv = 1, lower = 1, upper = 1) %>% 
  bind_rows(
    bind_rows(dta_plot_survt1)  %>% 
      select(cd4_vf, time, predsurv, lower, upper) %>% 
      group_by(cd4_vf, time) %>% 
      mutate(across(predsurv:upper, ~mean(.))) %>% 
      filter(row_number() == 1) %>% ungroup()) %>% 
  arrange(cd4_vf, time) %>% 
  rename(surv = predsurv, value = cd4_vf)

predt1aids <- bind_rows(dta_plot_survt1) %>% 
  select(aids, time, predsurv, lower, upper) %>% 
  distinct(aids, .keep_all = T) %>% 
  mutate(time = 0, predsurv = 1, lower = 1, upper = 1) %>% 
  bind_rows(
    bind_rows(dta_plot_survt1)  %>% 
      select(aids, time, predsurv, lower, upper) %>% 
      group_by(aids, time) %>% 
      mutate(across(predsurv:upper, ~mean(.))) %>% 
      filter(row_number() == 1) %>% ungroup()) %>% 
  arrange(aids, time) %>% 
  rename(surv = predsurv, value = aids) %>% 
  mutate(value = fct_recode(value, Never = "never", Mild = "mild", 
                            "Severe" = "mod_sev"))

predt1age <- bind_rows(dta_plot_survt1) %>% 
  select(age_f, time, predsurv, lower, upper) %>% 
  distinct(age_f, .keep_all = T) %>% 
  mutate(time = 0, predsurv = 1, lower = 1, upper = 1) %>% 
  bind_rows(
    bind_rows(dta_plot_survt1)  %>% 
      select(age_f, time, predsurv, lower, upper) %>% 
      group_by(age_f, time) %>% 
      mutate(across(predsurv:upper, ~mean(.))) %>% 
      filter(row_number() == 1) %>% ungroup()) %>% 
  arrange(age_f, time) %>% 
  rename(surv = predsurv, value = age_f) 

predt1all <- tibble(
  time = 0, predsurv = 1, lower = 1, upper = 1) %>% 
  bind_rows(
    bind_rows(dta_plot_survt1)  %>% 
      select(time, predsurv, lower, upper) %>% 
      group_by(time) %>% 
      mutate(across(predsurv:upper, ~mean(.))) %>% 
      filter(row_number() == 1) %>% ungroup()) %>% 
  arrange(time) %>% 
  rename(surv = predsurv)

predt1 <- bind_rows(ALL = predt1all, CD4 = predt1cd, AIDS = predt1aids, AGE = predt1age,
                    .id = "vbl") %>% 
  mutate(dta = "ART start6") %>% 
  select(dta, value, everything())

rm(predt1cd, predt1aids, predt1age, dta_plot_survt1)

#...............................................................................
obst1cd <- survfit(Surv(fupys5, fail5) ~ cd4_vf, dta0t1) %>% 
  summary(., times = seq(0, 5, 5/50)) %>% 
  keep_at(c("time", "surv", "lower", "upper", "strata")) %>% 
  as_tibble() %>% 
  rename(value = strata) %>% 
  mutate(value = str_replace(value, "cd4_vf=", ""), 
         value = factor(value, levels = c('[0,50)', '[50,100)', '[100,200)',
                                          '[200,350)', '[350,500)', '[500,Inf)')))

obst1aids <- survfit(Surv(fupys5, fail5) ~ aids2, dta0t1) %>% 
  summary(., times = seq(0, 5, 5/50)) %>% 
  keep_at(c("time", "surv", "lower", "upper", "strata")) %>% 
  as_tibble() %>% 
  rename(value = strata) %>% 
  mutate(value = str_replace(value, "aids2=", "")) %>% 
  mutate(value = fct_recode(value, Never = "never", Mild = "mild", 
                            "Severe" = "mod_sev"))

obst1age <- survfit(Surv(fupys5, fail5) ~ age_f, dta0t1) %>% 
  summary(., times = seq(0, 5, 5/50)) %>% 
  keep_at(c("time", "surv", "lower", "upper", "strata")) %>% 
  as_tibble() %>% 
  rename(value = strata) %>% 
  mutate(value = str_replace(value, "age_f=", ""), 
         value = factor(value, levels = c('[16,30)', '[30,40)', '[40,50)', 
                                          '[50,60)', '[60,Inf)')))

obst1all <- survfit(Surv(fupys5, fail5) ~ 1, dta0t1) %>% 
  summary(., times = seq(0, 5, 5/50)) %>% 
  keep_at(c("time", "surv", "lower", "upper")) %>% 
  as_tibble() 

obst1 <- bind_rows(ALL = obst1all, CD4 = obst1cd, AIDS = obst1aids, AGE = obst1age,
                   .id = "vbl") %>% 
  mutate(dta = "ART start6") %>% 
  select(dta, vbl, value, everything())
head(obst1)

rm(obst1cd, obst1aids, obst1age)

#...............................................................................
t1 <- bind_rows(Predicted = predt1, Observed = obst1, .id = "data")
rm(predt1, obst1)
head(t1)

rm(dta0t1, modelt1, pipelinet1)

# ******************************************************************************
# COMPILE ----
# ******************************************************************************

tplot <- bind_rows(t1= t1, t0 = t0, .id = "dta") %>% 
  mutate(vbl0 = vbl, 
         vbl = str_replace(vbl, "^CD4\\d?$", "CD4 (cells/μL)"), 
         vbl = str_replace(vbl, "^AGE$", "AGE (years)"), 
         vbl = str_replace(vbl, "^ALL$", "OVERALL"), 
         vbl = str_replace(vbl, "^AIDS$", "AIDS EVENT SEVERITY"), 
         dta = str_replace(dta, "t0", "ART start"), 
         dta = str_replace(dta, "t1", "ART start+6m")) 

dta_cal <- bind_rows(t1= t1, t0 = t0, .id = "dta") %>% 
  mutate(data = str_sub(str_to_lower(data), 1, 3)) %>% 
  select(dta, vbl, value, time, surv, data) %>% 
  pivot_wider(id_cols = dta:time, names_from = data, values_from = surv) %>% 
  mutate(across(pre:obs, ~.*100)) %>% 
  mutate(dif = abs(pre - obs), 
         months = time*12)

t_format <- list(t0 = t0, t1 = t1) %>% 
  map(~{..1 %>% 
      mutate(across(lower:upper, ~format(round(., 3), nsmall = 3))) %>% 
      unite(ci, lower, upper, sep = ", ") %>% 
      mutate(data = str_to_lower(str_sub(data, 0, 1)), 
             ci   = str_replace(ci, "^", "\\("),
             ci   = str_replace(ci, "$", "\\)")) %>% 
      pivot_wider(id_cols = c(vbl, value, time), 
                  names_from = data, values_from = surv:ci, 
                  names_vary = "slowest") %>% 
      select(vbl, value, time, ends_with("o"), ends_with("p"))
  }) 

left_join(t_format$t0, t_format$t1, join_by(vbl, value, time)) %>% 
  filter(time %in% seq(0, 5, 0.5)) %>% 
  write_csv(here::here("results_models_calibration.csv"))
  


# ******************************************************************************
# END ----
# ******************************************************************************




