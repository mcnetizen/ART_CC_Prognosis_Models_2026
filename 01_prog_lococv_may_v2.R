# ******************************************************************************
# OUTER CV T0 ----
# ******************************************************************************
dtat0 <- readRDS(here::here("list_dta_prognosis_2moddev.rds"))[[1]]  %>%
  select(id, cohort, kfold, kfold_no,
         fail5, fupys5,
         sex, idu, aidscdc, m_age50, m_cd4, m_rna5)
results1 <- "results_cv_margaret_t0b.rds"
results2 <- "results_cv_margaret_t0b_dtaplot2.rds"

n <- length(unique(dtat0$kfold))
iteration      <- vector("list", n)
iteration_cal <- vector("list", n)
for (i in 1:n){

  dta_out_test  <- subset(dtat0, kfold_no == i)
  dta_out_train <- subset(dtat0, kfold_no != i)

  # *** c-index ----
  mod   <- stpm2(Surv(fupys5, fail5) ~ idu + aidscdc + m_age50 + m_cd4 + m_rna5,
                 link.type = 'PH', smooth.formula = 
                   ~nsx(log(fupys5), df = 1)*idu + nsx(log(fupys5), df = 1)*m_cd4, 
                 dta_out_train)
  pred  <- predict(mod, dta_out_test, type = "loghazard")
  surv  <- with(dta_out_test, Surv(fupys5, fail5))
  w     <- Hmisc::rcorr.cens(1-pred, surv)
  kfold <- unique(dta_out_test$kfold)

  cvout = list(model = summary(mod)@call,
               cstat = tibble(cstat = w[["C Index"]], se = w[['S.D.']]/2, 
                              n = length(dta_out_test$id),
                              kfold_no = unique(dta_out_test$kfold_no), 
                              kfold = kfold))

  iteration[[i]] <-list(cvout = cvout)
  rm(pred, surv, w)
  
  #...............................................................................
  # ** calibration survival plots ----
  predtime <- seq(0, 5, 5/50)
  predtime <- predtime[-1]
  dta_plot_surv <- vector("list", 50)
  for (j in predtime) {
    newtestdata   <- dta_out_test %>% mutate(fupys5 = j)
    x             <- match(j, predtime)
    preddata      <- predict(mod, newtestdata, type = "surv", se.fit = T)
    
    dta_plot_surv[[x]] <- tibble(kfold    = dta_out_test$kfold, 
                                 time     = newtestdata$fupys5,
                                 predsurv = preddata$Estimate, 
                                 lower    = preddata$lower, 
                                 upper    = preddata$upper)
    rm(x, predlogh, preddata, newtestdata)
  }
  
  iteration_cal[[i]] <- list(
    kfold   = unique(dta_out_test$kfold),
    pred    = bind_rows(dta_plot_surv), 
    obsall  = survfit(Surv(fupys5, fail5) ~ 1, dta_out_test))
  
  rm(dta_plot_surv, predtime)
  rm(mod)
  rm(dta_out_test, dta_out_train, dta_plot)
}

cv <- iteration %>% set_names(paste0(rep("i",n), c(1:n)))

#...............................................................................
saveRDS(cv, here::here(results1))
saveRDS(iteration_cal, here::here(results2))

# ******************************************************************************
# OUTER CV T1 ----
# ******************************************************************************
dtat0 <- readRDS(here::here("list_dta_prognosis_2moddev.rds"))[[2]]  %>%
  select(id, cohort, kfold, kfold_no,
         fail5, fupys5,
         sex, idu, aidscdc, m_age50, m_cd4, m_rna5)
results1 <- "results_cv_margaret_t1b.rds"
results2 <- "results_cv_margaret_t1b_dtaplot2.rds"

n <- length(unique(dtat0$kfold))
iteration      <- vector("list", n)
iteration_cal <- vector("list", n)
for (i in 1:n){
  
  dta_out_test  <- subset(dtat0, kfold_no == i)
  dta_out_train <- subset(dtat0, kfold_no != i)
  
  # *** c-index ----
  mod   <- stpm2(Surv(fupys5, fail5) ~ idu + aidscdc + m_age50 + m_cd4 + m_rna5,
                 link.type = 'PH', smooth.formula = 
                   ~nsx(log(fupys5), df = 1)*idu + nsx(log(fupys5), df = 1)*m_cd4, 
                 dta_out_train)
  pred  <- predict(mod, dta_out_test, type = "loghazard")
  surv  <- with(dta_out_test, Surv(fupys5, fail5))
  w     <- Hmisc::rcorr.cens(1-pred, surv)
  kfold <- unique(dta_out_test$kfold)
  
  cvout = list(model = summary(mod)@call,
               cstat = tibble(cstat = w[["C Index"]], se = w[['S.D.']]/2, 
                              n = length(dta_out_test$id),
                              kfold_no = unique(dta_out_test$kfold_no), 
                              kfold = kfold))
  
  iteration[[i]] <-list(cvout = cvout)
  rm(pred, surv, w)
  
  #...............................................................................
  # ** calibration survival plots ----
  predtime <- seq(0, 5, 5/50)
  predtime <- predtime[-1]
  dta_plot_surv <- vector("list", 50)
  for (j in predtime) {
    newtestdata   <- dta_out_test %>% mutate(fupys5 = j)
    x             <- match(j, predtime)
    preddata      <- predict(mod, newtestdata, type = "surv", se.fit = T)
    
    dta_plot_surv[[x]] <- tibble(kfold    = dta_out_test$kfold, 
                                 time     = newtestdata$fupys5,
                                 predsurv = preddata$Estimate, 
                                 lower    = preddata$lower, 
                                 upper    = preddata$upper)
    rm(x, predlogh, preddata, newtestdata)
  }
  
  iteration_cal[[i]] <- list(
    kfold   = unique(dta_out_test$kfold),
    pred    = bind_rows(dta_plot_surv), 
    obsall  = survfit(Surv(fupys5, fail5) ~ 1, dta_out_test))
  
  rm(dta_plot_surv, predtime)
  rm(mod)
  rm(dta_out_test, dta_out_train, dta_plot)
}

cv <- iteration %>% set_names(paste0(rep("i",n), c(1:n)))

#...............................................................................
saveRDS(cv, here::here(results1))
saveRDS(iteration_cal, here::here(results2))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# END ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

