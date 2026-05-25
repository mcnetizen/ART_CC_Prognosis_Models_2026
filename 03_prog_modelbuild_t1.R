# ******************************************************************************
# DATA T1
dta0 <- readRDS(here::here("list_dta_prognosis_2moddev.rds"))[[2]] %>%
  select(id, cohort, kfold, kfold_no,
         fail5, fupys5,
         sexm3, aids2, age10, cd4100, rna, undetect,
         cd4_vf, cd4_v, rna_v, age, age_f)
results1 <- "results_model_t1.rds"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MODEL BUILD ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dta_out_train <- dta0 %>% mutate(ifold = kfold_no) 
dta_out_train %>% count(kfold, kfold_no, cohort, ifold)

# number of inner-kfolds
z <- length(unique(dta_out_train$ifold))
  
#...............................................................................
# step1-part one: sequence of variables for step 1
options <- list(i  = rep(c(1:z), 3), 
                xh = unlist(map(c("cd4100", "age10", "undetect + rna"), ~ rep(.x, z))), 
                xo = unlist(map(c("cd4100", "age10", "undetect + rna"), ~ rep(.x, z))))
step1   <- pmap_df(options, ~fx1_eval_vblform(..1,..2,..3)) %>% fx_cstatsumm()
(order_h   <- step1$h$option)
(order_o   <- step1$h$option)
rm(options)
  
#...............................................................................
# step1- part two: functional form of variables
step1a <- fx1_select_vblform_t1(order_h[1], order_o[1])
(p1    <- c(h = step1a$h$option[1], o = step1a$o$option[1]))
  
step1b <- fx1_select_vblform_t1(order_h[2], order_o[2], p1["h"], p1["o"])
(p2 <- map2(p1, c(step1b$h$option[1], step1b$o$option[1]), ~paste(..1,..2, sep = " + ")))
  
step1c <- fx1_select_vblform_t1(order_h[3], order_o[3], p2["h"], p2["o"])
(p3 <- map2(p2, c(step1c$h$option[1], step1c$h$option[1]), ~paste(..1,..2, sep = " + ")))
  
(pred <- map(p3, ~paste("sexm3 + aids2",..1, sep = " + ")))
rm(p1, p2, p3)
  
#...............................................................................
# step2: baseline hazard
options <- list(i = rep(c(1:z), 3),
                x = c(rep(1, z), rep(2, z), rep(3, z)))
step2   <- pmap_df(options, ~fx2_eval_stpmsteps(..1, pred$h, pred$o,..2,..2),
                   .progress = "step2, model_build") %>% fx_cstatsumm()
rm(options)
  
#...............................................................................
# step3: interactions
(pred2   <- map(pred, ~fx3_interactions_t1(.)))
(pred2   <- map(pred2, ~ as_tibble(.) %>%
                  pivot_longer(everything(), 
                               names_to = "terms", values_to = "terms2")))
options  <- list(i   = rep(c(1:z), 3),
                 # do not re-run the additive
                 ph  = unlist(map(pred2$h$terms2[2:4], ~ rep(.x, z))),
                 po  = unlist(map(pred2$o$terms2[2:4], ~ rep(.x, z))), 
                 x   = unlist(map(pred2$h$terms[2:4],  ~ rep(.x, z))))
results  <- pmap(options, 
                 safely(~fx2_eval_stpmsteps(
                   ..1,..2,..3, step2$h$option[1], step2$o$option[1], ..4)),
                 .progress = "step3, model_build")
step3    <- results %>% map("result") %>% compact() %>% bind_rows() %>% fx_cstatsumm()
# to capture and locate any error
names    <- map2(options[[1]], options[[4]], ~paste(..1,..2, sep = "_"))
(errors3 <- results %>% set_names(names) %>% map("error") %>% compact())
# add the results of the additive from step2 
step3   <- map2(step3, list(step2$h[1, ], step2$o[1, ]),
                ~{bind_rows(.x, .y) %>% arrange(desc(meancstat))})
# add the full terms
step3   <- map2(step3, pred2, ~left_join(.x, .y, join_by(terms)))
step3   <- map(step3, ~ .x %>% rename(df = option, option = terms))
rm(options, results, names)
  
#...............................................................................
# step4: TDE-h
(bestform <- tibble(scale       = c("ph", "po"),
                    basedf      = map_vec(step2,~.x$option[1]), 
                    main_terms  = map_vec(step3,~.x$terms2[1]), 
                    cstat_ph_po = map_vec(step2,~.x$meancstat[1])
                    ))
  
step4 <- pmap(map(bestform, ~c(.)), ~fx4_sfs_tve(..1,..2,..3,..4),
              .progress = "step4, model_build") 
step4 <- step4 %>% set_names(c("h", "o"))
(tve <- (list(h = step4$h %>% map("result") %>% flatten(), 
              o = step4$o %>% map("result") %>% flatten())))
  
#...............................................................................
# scale & model form
(bestform <- tibble(bestform, 
                    tve_terms  = map_vec(map(tve, "best"), ~str_c(., collapse = ", ")),
                    cstat_last = unlist(map(tve, "cstat"))) %>% 
   arrange(desc(cstat_last)))
  
# translate
bestform[1,]
formula   <- paste("stpm2(Surv(fupys5, fail5)", bestform[[1,"main_terms"]], sep = " ~ ")
link      <- paste0("link.type = '", str_to_upper(bestform[[1,"scale"]]), "'")
base      <- "smooth.formula = ~ "
time      <- paste0("nsx(log(fupys5), df = ", bestform[[1,"basedf"]], ")")
if (bestform[[1,"tve_terms"]] == "none") {
  timetve <- paste0(base, time)
} else {
  timetve <- paste(time, str_split(bestform[[1,"tve_terms"]], ",", simplify = T), sep = "*")
  timetve <- paste0(base, paste0(timetve, collapse = " + "))
}
(bestform_mod <- paste(formula, link, timetve, sep = ", "))
rm(formula, link, timetve, base, time, tve)
  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAVE ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dta0 %>% nrow() ==  dta_out_train %>% nrow()
rm(dta_out_train)

step3 = c(step3, list(errors = errors3))
pipe  = tibble::lst(step1a, step1b, step1c, step2, step3, step4)
model = eval(parse(text = paste(bestform_mod, "dta0)", sep = ", ")))

pipeline <- list(pipe = pipe, model = model)
saveRDS(pipeline, here::here(results1))

rm(step1, step1a, step1b, step1c, step2, step3, step4, pred2, order_h, order_o)
rm(bestform, errors3, pred)
rm(pipe, model)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# END ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
