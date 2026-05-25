# ******************************************************************************
# PACAKGES 
# ******************************************************************************

if (!require(here))      library(here)
if (!require(tidyverse)) library(tidyverse)
if (!require(janitor))   library(janitor)
if (!require(survival))  library(survival)
if (!require(survminer)) library(survminer) # survival plots
if (!require(rstpm2))    library(rstpm2)    # stpm2
if (!require(pROC))      library(pROC)      # auc for logit
if (!require(rms))       library(rms)       # rcs
if (!require(Hmisc))     library(Hmisc)     # rcorr.cens 

# ******************************************************************************
# FUNCITONS TO IMPLEMENT LOCOCV 
# ******************************************************************************
# STEP1 ----
# eval_form evaluates the order of the continuous predictors,
# then, used inside select_vblform used to obtain the best functional form for each 
# cstatsumm is used in all steps. It provides the mean summary of the evaluations
fx1_eval_vblform <- function(IFOLD, VBL_H, VBL_O, PRED_H=NULL, PRED_O=NULL){
  dta_train <- subset(dta_out_train, ifold != IFOLD)
  dta_test  <- subset(dta_out_train, ifold == IFOLD)
  
  # the predictors are null to define the ordering of vbls/terms, but
  # these will be used when for the second and third rounds as prior selected
  # terms need to be added
  if (is.null(PRED_H) && is.null(PRED_O)) {
    terms_h <- VBL_H
    terms_o <- VBL_O
  } else {
    terms_h <- paste(PRED_H, VBL_H, sep = " + ")
    terms_o <- paste(PRED_O, VBL_O, sep = " + ")
  }
  
  formula_h <- paste("coxph(Surv(fupys5, fail5) ~ sexm3 + aids2", terms_h, sep = " + ")
  mod_h     <- eval(parse(text = paste(formula_h, "dta_train)", sep = ", ")))
  pred_h    <- predict(mod_h, dta_test, type = "lp")
  
  formula_o <- paste("glm(family = binomial, fail5 ~ sexm3 + aids2", terms_o, sep = " + ")
  mod_o     <- eval(parse(text = paste(formula_o, "dta_train)", sep = ", ")))
  pred_log  <- predict(mod_o, dta_test, type = "link")
  pred_prob <- predict(mod_o, dta_test, type = "response")
  
  surv      <- with(dta_test, Surv(fupys5, fail5))
  
  data.frame(
    option_h = VBL_H,
    option_o = VBL_O,
    itest    = unique(dta_test$ifold),
    itestn   = unique(dta_test$kfold),
    cstat_h  = Hmisc::rcorr.cens(1-pred_h, surv)[["C Index"]],
    cstat_o  = Hmisc::rcorr.cens(1-pred_log, surv)[["C Index"]],
    auc_o    = pROC::auc(dta_test$fail5, pred_prob,
                         levels = c(0, 1), direction = "<")[1], 
    terms_h   = str_replace_all(terms_h, "\\s\\+", ","),
    terms_o   = str_replace_all(terms_o, "\\s\\+", ","))
}

fx_cstatsumm <- function(DATA){
  
  map(c(h = "_h", o = "_o"),
      ~{DATA %>% 
          select(itest, ends_with(.x)) %>% 
          rename_with(~str_remove(., "_\\D{1}")) %>% 
          summarise(meancstat = mean(cstat), no_tests = n(), 
                    .by = c(terms, option)) %>%
          arrange(desc(meancstat))})
}

fx1_select_vblform_t0 <- function(VBL_H, VBL_O, PRED_H=NULL, PRED_O=NULL){
  i       <-  rep(c(1:z), 3)
  xh      <-  c(rep(VBL_H, z), 
                rep(paste0("rcs(", VBL_H, ", 3)"), z), 
                rep(paste0("rcs(", VBL_H, ", 4)"), z)) 
  xo      <-  c(rep(VBL_O, z), 
                rep(paste0("rcs(", VBL_O, ", 3)"), z), 
                rep(paste0("rcs(", VBL_O, ", 4)"), z))
  list    <- list(i, xh, xo)
  
  pmap_df(list, ~fx1_eval_vblform(..1,..2,..3, PRED_H, PRED_O)) %>% 
    fx_cstatsumm()
} 

fx1_select_vblform_t1 <- function(VBL_H, VBL_O, PRED_H=NULL, PRED_O=NULL){
  i       <-  rep(c(1:z), 2)
  
  if (str_detect(VBL_H, "rna") == T) {
    xh      <-  c(rep(VBL_H, z), 
                  rep("undetect + rcs(rna, 3)", z))
  } else {
    xh      <-  c(rep(VBL_H, z), 
                  rep(paste0("rcs(", VBL_H, ", 3)"), z)) 
  }
  
  if (str_detect(VBL_O, "rna") == T) {
    xo      <-  c(rep(VBL_O, z), 
                  rep("undetect+ rcs(rna, 3)", z))
  } else {
    xo      <-  c(rep(VBL_O, z), 
                  rep(paste0("rcs(", VBL_O, ", 3)"), z)) 
  }
  
  list    <- list(i, xh, xo)
  
  pmap_df(list, ~fx1_eval_vblform(..1,..2,..3, PRED_H, PRED_O)) %>% 
    fx_cstatsumm()
}

#...............................................................................
# STEP2&3 ----
# eval_stpmsteps is used in step 2 and 3 
fx2_eval_stpmsteps <- function(IFOLD, TERMS_H, TERMS_O, DF_H, DF_O, X=NULL){
  dta_train <- subset(dta_out_train, ifold != IFOLD)
  dta_test  <- subset(dta_out_train, ifold == IFOLD)
  if (is.null(X)) {X <- "additive"} else {X <- X}
  
  formula_h <- paste("stpm2(Surv(fupys5, fail5)", TERMS_H, sep = " ~ ")
  link_h    <- "link.type = 'PH'"
  base_h    <- paste0("smooth.formula = ~ nsx(log(fupys5), df = ", DF_H, ")")
  
  formula_o <- paste("stpm2(Surv(fupys5, fail5)", TERMS_O, sep = " ~ ")
  link_o    <- "link.type = 'PO'"
  base_o    <- paste0("smooth.formula = ~ nsx(log(fupys5), df = ", DF_O, ")")
  
  mod_h     <- eval(parse(text = paste(formula_h, link_h, base_h, "dta_train)", sep = ", ")))
  pred_h    <- predict(mod_h, dta_test, type = "loghazard")
  
  mod_o     <- eval(parse(text = paste(formula_o, link_o, base_o, "dta_train)", sep = ", ")))
  pred_o    <- predict(mod_o, dta_test, type = "loghazard")
  
  surv      <- with(dta_test, Surv(fupys5, fail5))
  
  data.frame(
    terms_h     = X,
    terms_o     = X,
    option_h    = DF_H,
    option_o    = DF_O,
    itest       = unique(dta_test$ifold),
    itestn      = unique(dta_test$kfold),
    cstat_h     = Hmisc::rcorr.cens(1-pred_h, surv)[["C Index"]],
    cstat_o     = Hmisc::rcorr.cens(1-pred_o, surv)[["C Index"]])
}

fx3_interactions_t0 <- function(PRED){
  p <- str_split_1(PRED, "\\s\\+\\s")
  l <- map(p, ~c(.)) %>% set_names(str_remove_all(p, "rcs\\(|,\\s\\d{1}\\)"))
  # leave the interaction as it will be useful when identifying the best option
  additive     <- PRED
  agexcd4      <- paste(paste(l$sexm3, l$aids2, l$rna, sep = " + "),
                        paste(l$age10, l$cd4100, sep = "*"), sep = " + ")
  aidsxcd4     <- paste(paste(l$sexm3, l$age10, l$rna, sep = " + "),
                        paste(l$aids2,  l$cd4100, sep = "*"), sep = " + ")
  aidsxagexcd4 <- paste(paste(l$sexm3, l$rna, sep = " + "),
                        paste(l$aids2,  l$age10, l$cd4100, sep = "*"), sep = " + ")
  
  list(additive, agexcd4, aidsxcd4, aidsxagexcd4) %>% 
    set_names("additive", "age_cd4", "aids_cd4", "aids_age_cd4")
}

fx3_interactions_t1 <- function(PRED){
  p <- str_split_1(PRED, "\\s\\+\\s")
  l <- map(p, ~c(.)) %>% set_names(str_remove_all(p, "rcs\\(|,\\s\\d{1}\\)"))
  # leave the interaction as it will be useful when identifying the best option
  additive     <- PRED
  agexcd4      <- paste(paste(l$sexm3, l$aids2, l$undetect, l$rna, sep = " + "),
                        paste(l$age10, l$cd4100, sep = "*"), sep = " + ")
  aidsxcd4     <- paste(paste(l$sexm3, l$age10, l$undetect, l$rna, sep = " + "),
                        paste(l$aids2,  l$cd4100, sep = "*"), sep = " + ")
  aidsxagexcd4 <- paste(paste(l$sexm3, l$undetect, l$rna, sep = " + "),
                        paste(l$aids2,  l$age10, l$cd4100, sep = "*"), sep = " + ")
  
  list(additive, agexcd4, aidsxcd4, aidsxagexcd4) %>% 
    set_names("additive", "age_cd4", "aids_cd4", "aids_age_cd4")
}

#...............................................................................
# STEP4 ----
# Note1- the tve interaction is modeled w/ the same functional form used in main effect
# the df for the tve is set to be equal to the baseline complexity
# this will run in the two scales separately because the sequence might stop at different rounds
# c-stat tolerance is fixed at 0.02
fx4_eval_tve <- function(IFOLD, SCALE, DF, TERMS,
                         ROUND, TVE_i, TVE1, TVE2, TVE3, TVE4){
  dta_train <- subset(dta_out_train, ifold != IFOLD)
  dta_test  <- subset(dta_out_train, ifold == IFOLD)
  formula   <- paste("stpm2(Surv(fupys5, fail5)", TERMS, sep = " ~ ")
  link      <- paste0("link.type = '", str_to_upper(SCALE), "'")
  base      <- "smooth.formula = ~ "
  time      <- paste0("nsx(log(fupys5), df = ", DF, ")")
  
  if (ROUND == 1L) {
    
    timetve <- paste(paste0(base, time), TVE_i, sep = "*")
    vbls    <- TVE_i
  } else if (ROUND == 2L) {
    
    timetve1 <- paste(paste0(base, time), TVE1, sep = "*")
    timetve2 <- paste(timetve1, paste(time, TVE_i, sep = "*"), sep = " + ")
    timetve  <- timetve2
    vbls     <- paste(TVE1, TVE_i, sep = ", ")
  } else if (ROUND == 3L) {
    
    timetve1 <- paste(paste0(base, time), TVE1, sep = "*")
    timetve2 <- paste(timetve1, paste(time, TVE2, sep = "*"), sep = " + ")
    timetve3 <- paste(timetve2, paste(time, TVE_i, sep = "*"), sep = " + ")
    timetve  <- timetve3
    vbls     <- paste(TVE1, TVE2, TVE_i, sep = ", ")
  } else if (ROUND == 4L) {
    
    timetve1 <- paste(paste0(base, time), TVE1, sep = "*")
    timetve2 <- paste(timetve1, paste(time, TVE2, sep = "*"), sep = " + ")
    timetve3 <- paste(timetve2, paste(time, TVE3, sep = "*"), sep = " + ")
    timetve4 <- paste(timetve3, paste(time, TVE_i, sep = "*"), sep = " + ")
    timetve  <- timetve4
    vbls     <- paste(TVE1, TVE2, TVE3, TVE_i, sep = ", ")
  } else {
    
    timetve1 <- paste(paste0(base, time), TVE1, sep = "*")
    timetve2 <- paste(timetve1, paste(time, TVE2, sep = "*"), sep = " + ")
    timetve3 <- paste(timetve2, paste(time, TVE3, sep = "*"), sep = " + ")
    timetve4 <- paste(timetve3, paste(time, TVE4, sep = "*"), sep = " + ")
    timetve5 <- paste(timetve4, paste(time, TVE_i, sep = "*"), sep = " + ")
    timetve  <- timetve5
    vbls     <- paste(TVE1, TVE2, TVE3, TVE4, TVE_i, sep = ", ")
  }
  
  mod   <- eval(parse(text = paste(formula, link, timetve, "dta_train)", sep = ", ")))
  pred  <- predict(mod, dta_test, type = "loghazard")
  surv  <- with(dta_test, Surv(fupys5, fail5))
  
  data.frame(
    scale   = SCALE,
    round   = ROUND,
    tve     = vbls,
    itest   = unique(dta_test$ifold),
    itestn  = unique(dta_test$kfold),
    cstat   = Hmisc::rcorr.cens(1-pred, surv)[["C Index"]])
}

fx4_sfs_tve <- function(SCALE, DF, TERMS, CSTAT_PROPORTIONAL){
  tve_options <- unlist(str_split(str_remove(TERMS, "undetect\\s\\+"), pattern = " \\+ |\\*"))
  cstat0      <- CSTAT_PROPORTIONAL
  roundx       <- vector("list", 5)
  tolerance   <- 0.02
  #.............................................................................
  # this function creates the options for eval_tve and summarises the results
  fx4_select_tve <- function(SCALE, DF, TERMS, TVE_OPTIONS, TVE1, TVE2, TVE3, TVE4){
    # tve_options need to be a vector of all the possible options for a given round
    # it also needs to be re-defined in each round
    options   <- list(i     = rep(c(1:z), length(TVE_OPTIONS)),
                      tve_i = unlist(map(TVE_OPTIONS, ~ rep(.x, z))))
    round_i   <- 6 - length(TVE_OPTIONS)
    results   <- pmap(options, 
                      safely(~fx4_eval_tve(..1, SCALE, DF, TERMS, 
                                           round_i, ..2, TVE1, TVE2, TVE3, TVE4)))
    #...........................................................................
    summaries <- results %>% map("result") %>% compact() %>% bind_rows() %>% 
      summarise(meancstat = mean(cstat), no_tests = n(), 
                .by = c(scale, round, tve)) %>%
      arrange(desc(meancstat)) 
    #...........................................................................
    names     <- map2(options[[1]], options[[2]], ~paste(..1,..2, sep = "_"))
    errors    <- results %>% set_names(names) %>% map("error") %>% compact()
    
    list(summaries, errors) %>% set_names(c("summaries", "errors"))
  }
  #.............................................................................
  # ROUND 1
  roundx[[1]]  <- fx4_select_tve(SCALE, DF, TERMS, tve_options)
  cstatchange <- roundx[[1]]$summaries$meancstat[1] - cstat0
  if (cstatchange < tolerance) {
    
    roundx[[1]]$result       <- "proportional model better"
    roundx[[1]]$result$best  <- "none"
    roundx[[1]]$result$cstat <- cstat0 
    return(compact(roundx))
    
  } else {
    # ROUND 2
    tve1        <- roundx[[1]]$summaries$tve[1]
    tve_options <- tve_options[!tve_options %in% tve1]
    roundx[[2]]  <- fx4_select_tve(SCALE, DF, TERMS, tve_options, tve1)
    cstatchange <- roundx[[2]]$summaries$meancstat[1] - roundx[[1]]$summaries$meancstat[1]
    if (cstatchange < tolerance) {
      
      roundx[[2]]$result       <- "round2 shows no improvement"
      roundx[[2]]$result$best  <- tve1
      roundx[[2]]$result$cstat <- roundx[[1]]$summaries$meancstat[1]
      return(compact(roundx))
      
    } else {
      # ROUND 3
      tve2        <- unlist(str_split(roundx[[2]]$summaries$tve[1], ", "))
      tve_options <- tve_options[!tve_options %in% tve2]
      roundx[[3]]  <-  fx4_select_tve(SCALE, DF, TERMS, tve_options, 
                                     tve2[1], tve2[2])
      cstatchange <- roundx[[3]]$summaries$meancstat[1] - roundx[[2]]$summaries$meancstat[1]
      if (cstatchange < tolerance) {
        
        roundx[[3]]$result       <- "round3 shows no improvement"
        roundx[[3]]$result$best  <- paste(tve2)
        roundx[[3]]$result$cstat <- roundx[[2]]$summaries$meancstat[1]
        return(compact(roundx))
        
      } else {
        # ROUND 4
        tve3        <- unlist(str_split(roundx[[3]]$summaries$tve[1], ", "))
        tve_options <- tve_options[!tve_options %in% tve3]
        roundx[[4]]  <-  fx4_select_tve(SCALE, DF, TERMS, tve_options, 
                                       tve3[1], tve3[2], tve3[3])
        cstatchange <- roundx[[4]]$summaries$meancstat[1] - roundx[[3]]$summaries$meancstat[1]
        if (cstatchange < tolerance) {
          
          roundx[[4]]$result       <- "round4 shows no improvement"
          roundx[[4]]$result$best  <- paste(tve3)
          roundx[[4]]$result$cstat <- roundx[[3]]$summaries$meancstat[1]
          return(compact(roundx))
          
        } else {
          # ROUND 5
          tve4        <- unlist(str_split(roundx[[4]]$summaries$tve[1], ", "))
          tve_options <- tve_options[!tve_options %in% tve4]
          roundx[[5]]  <- fx4_select_tve(SCALE, DF, TERMS, tve_options, 
                                        tve4[1], tve4[2], tve4[3], tve4[4])
          cstatchange <- roundx[[5]]$summaries$meancstat[1] - roundx[[4]]$summaries$meancstat[1]
          if (cstatchange < tolerance) {
            
            roundx[[5]]$result       <- "round5 shows no improvement"
            roundx[[5]]$result$best  <- paste(tve4)
            roundx[[5]]$result$cstat <- roundx[[4]]$summaries$meancstat[1]
            return(roundx)
            
          } else {
            roundx[[5]]$result       <- "round5 is best. use all"
            roundx[[5]]$result$best  <- TERMS
            roundx[[5]]$result$cstat <- roundx[[5]]$summaries$meancstat[1]
            return(roundx)
          }
        }
      }
    }
  }
}
