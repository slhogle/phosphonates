Modelling prochlorococcus/SAR11 pepm abundance using environmental
covariates
================

  - [INTRO](#intro)
  - [READ DATA](#read-data)
  - [GENERAL FORMATTING](#general-formatting)
  - [DATA SPLITTING & RESAMPLING](#data-splitting-resampling)
  - [CROSS FOLD SET](#cross-fold-set)
  - [RANDOM FOREST MODEL (RANGER)](#random-forest-model-ranger)
  - [RECIPES](#recipes)
      - [CAR SCORE](#car-score)
      - [BORUTA](#boruta)
  - [WORKFLOW](#workflow)
      - [CARSCORE](#carscore)
      - [BORUTA](#boruta-1)
  - [TUNE HYPERPARAMETERS](#tune-hyperparameters)
      - [CARSCORE](#carscore-1)
      - [BORUTA](#boruta-2)
  - [FINAL FIT](#final-fit)
      - [MODEL](#model)
      - [RECIPE](#recipe)
      - [WORKFLOW](#workflow-1)
      - [EXECUTE](#execute)

# INTRO

useful tutorials:
<http://www.rebeccabarter.com/blog/2020-03-25_machine_learning>
<https://juliasilge.com/blog/sf-trees-random-tuning/>
<https://www.tidymodels.org/start/case-study/>
<https://bradleyboehmke.github.io/HOML/random-forest.html>

feature selection:
<http://blog.datadive.net/selecting-good-features-part-iii-random-forests/>
<https://academic.oup.com/bib/article/20/2/492/4554516>

SO questions on colinearity:
<https://stats.stackexchange.com/questions/168622/why-is-multicollinearity-not-checked-in-modern-statistics-machine-learning>
<https://stats.stackexchange.com/questions/141619/wont-highly-correlated-variables-in-random-forest-distort-accuracy-and-feature>

Evaluating performance
<https://bookdown.org/max/FES/measuring-performance.html>

``` r
library(here)
library(tidyverse)

library(tidymodels)
library(recipeselectors)
library(vip)

#library(parallel)
#library(foreach)
#library(doParallel)

library(ranger)
library(Boruta)
library(care)

library(ggforce)
```

# READ DATA

Preformatted dataset

``` r
pepm <- readRDS(here::here("input", "pepm.global.final")) %>%
  group_by(sampleID, group) %>% 
  slice(1) %>%
  ungroup() %>%
  filter(sampleID != "S0468")
```

# GENERAL FORMATTING

``` r
a <- left_join(filter(pepm, group=="bactarc") %>% select(sampleID, bactarc_pepm=RA),
               filter(pepm, group=="sar11") %>% select(sampleID, sar11_pepm=RA)) %>%
  mutate(sar11_pepm=ifelse(is.na(sar11_pepm), 1e-4, sar11_pepm)) %>%
  mutate(sar11_pepm=ifelse(sar11_pepm > 1, 0.99, sar11_pepm)) %>%
  mutate(bactarc_pepm=ifelse(is.na(bactarc_pepm), 1e-4, bactarc_pepm)) %>%
  mutate(bactarc_pepm=ifelse(bactarc_pepm > 1, 0.99, bactarc_pepm))

PEPMreads2forest <- pepm %>%
  select_if(!str_detect(names(.), "imputed_")) %>%
  select_if(!str_detect(names(.), "pel_")) %>%
  select_if(!str_detect(names(.), "pro_LLVII")) %>% 
  select_if(!str_detect(names(.), "Long")) %>% 
  filter(group=="prochlorococcus") %>%
  select(-section, -lat, -lon, -group, -W, -M,
         -synDarwin_umolC.kg, -proDarwin_umolC.kg,
         -DOFe_darwin_clim, -POFe_darwin_clim) %>% 
  drop_na() %>%
  mutate(RA=ifelse(RA > 1, 0.99, RA)) %>%
  mutate(RA=ifelse(is.na(RA), 1e-4, RA)) %>%
  left_join(., a) %>%
  select(-sampleID,
         -bacteria,                       # redundant with archaea
         -pro_LLII_LLIII,                 # redundant with pro_LLIV
         -DOironDarwin_dissolved_nmol.kg, # redundant with DOPDarwin_dissolved_umol.kg
         -copper_dissolved_nmol.kg,       # nitrateDarwin_dissolved_umol.kg
         -ironDarwin_dissolved_nmol.kg)   # nitrateDarwin_dissolved_umol.kg
```

# DATA SPLITTING & RESAMPLING

``` r
splits      <- initial_split(PEPMreads2forest, strata = RA, breaks=10)

pepm_other <- training(splits)
pepm_test  <- testing(splits)
```

# CROSS FOLD SET

``` r
folds <- vfold_cv(pepm_other,
                  v = 10,
                  strata = RA, 
                  breaks = 10)
```

# RANDOM FOREST MODEL (RANGER)

Perform with parallel processing

``` r
cores <- parallel::detectCores()
cores
```

``` r
rf_mod <- 
  rand_forest() %>%
  set_args(   mtry = tune(), 
              min_n = tune(), 
              trees = 1000) %>% 
  set_engine("ranger", 
             num.threads = cores,
             importance = "permutation", 
             splitrule = "beta",
             replace = FALSE) %>% 
  set_mode("regression")

rf_mod %>%    
  parameters()
```

# RECIPES

## CAR SCORE

This version contains the variable filter step using the CAR score
algorithm from the `care` package. See here:
<http://www.strimmerlab.org/software/care/>

The “care” package implements the regression approach of Zuber and
Strimmer (2011). CAR scores measure the correlation between the response
and the Mahalanobis-decorrelated predictors. The squared CAR score is a
natural measure of variable importance and provides a canonical ordering
of variables. This package provides functions for estimating CAR scores,
for variable selection using CAR scores, and for estimating
corresponding regression coefficients. Both shrinkage as well as
empirical estimators are available.

``` r
rf_recipe_carscore <- recipe(RA ~ ., data = pepm_other) %>%
  step_normalize(all_predictors(), -all_nominal()) %>%
  step_corr(all_predictors(), -all_nominal(), threshold = 0.9) %>%
  step_zv(all_predictors()) %>%
  step_select_carscore(all_predictors(), 
                       outcome = "RA", 
                       threshold = 0.5, # Retain features only within the top 50th percentile
                       top_p = tune())  # can be tuned
```

## BORUTA

``` r
rf_recipe_boruta <- recipe(RA ~ ., data = pepm_other) %>%
  step_normalize(all_predictors(), -all_nominal()) %>%
  step_corr(all_predictors(), -all_nominal(), threshold = 0.9) %>%
  step_zv(all_predictors()) %>%
  step_select_boruta(all_predictors(), outcome = "RA", res="boruta_tune")
```

# WORKFLOW

## CARSCORE

``` r
rf_workflow_carscore <- 
  workflow() %>% 
  add_recipe(rf_recipe_carscore) %>% 
  add_model(rf_mod)
```

### FEATURE SELECTION PARAMETERS

In case we are tuning top\_p or threshold we need to update the tuned
hyperparameters using upper bounds

``` r
p_info <- 
  rf_workflow_carscore %>% 
  parameters() %>% 
  update(threshold = threshold(c(0.5, 0.9))) %>%  # Pick an upper bound for the threshold
  update(top_p = top_p(c(1, 20))) %>%             # Pick an upper bound for number of variables 
  update(mtry = mtry(c(1, 10)))                   # Pick an upper bound for mtry 

ctrl <- control_grid(extract = identity,
                     verbose = TRUE)
```

## BORUTA

``` r
rf_workflow_boruta <- 
  workflow() %>% 
  add_recipe(rf_recipe_boruta) %>% 
  add_model(rf_mod)
```

# TUNE HYPERPARAMETERS

## CARSCORE

``` r
rf_res_carscore <-
  rf_mod %>% 
  tune_grid(
    rf_recipe,
    resamples = folds,
    grid = 20,
    param_info = p_info, # use only with the CARscore recipe
    control = ctrl,
    metrics = metric_set(rmse, rsq)
  )

saveRDS(rf_res_carscore, here::here("data", "rf_tune.carscore.pro"))
```

``` r
rf_res_carscore %>% 
  show_best(metric = "rsq")
```

``` r
rf_res_carscore %>% 
  show_best(metric = "rmse")
```

``` r
autoplot(rf_res_carscore)
```

``` r
rf_best_carscore <- 
  rf_res_carscore %>% 
  select_best(metric = "rsq")
rf_best_carscore
```

## BORUTA

``` r
rf_res_boruta <-
  rf_mod %>% 
  tune_grid(
    rf_recipe,
    resamples = folds,
    grid = 20,
    control = ctrl,
    metrics = metric_set(rmse, rsq)
  )

saveRDS(rf_res_boruta, here::here("data", "rf_tune.boruta.pro"))
```

``` r
rf_res_boruta %>% 
  show_best(metric = "rsq")
```

``` r
rf_res_boruta %>% 
  show_best(metric = "rmse")
```

``` r
autoplot(rf_res_boruta)
```

``` r
rf_best_boruta <- 
  rf_res_boruta %>% 
  select_best(metric = "rsq")
rf_best_boruta
```

# FINAL FIT

After selecting our best model and hyperparameter values, our last step
is to fit the final model on all the rows of data not originally held
out for testing (both the training and the validation sets combined),
and then evaluate the model performance one last time with the held-out
test set. We’ll start by building our parsnip model object again from
scratch. We take our best hyperparameter values from our random forest
model.

## MODEL

Will run with Boruta

``` r
final_rf_mod <- 
  rand_forest() %>%
  set_args(mtry = 9, 
              min_n = 4, 
              trees = 1000) %>% 
  set_engine("ranger", 
             num.threads = cores,
             importance = "permutation", 
             splitrule = "beta",
             replace = FALSE) %>% 
  set_mode("regression")
```

## RECIPE

``` r
final_rf_recipe <- recipe(RA ~ ., data = pepm_other) %>%
  step_normalize(all_predictors(), -all_nominal()) %>%
  step_corr(all_predictors(), -all_nominal(), threshold = 0.9) %>%
  step_zv(all_predictors()) %>%
  step_select_boruta(all_predictors(), outcome = "RA", res="boruta_tune")
  #step_select_carscore(all_predictors(), -all_nominal(), outcome = "RA", top_p = 16)
```

## WORKFLOW

``` r
final_rf_workflow <- 
  workflow() %>% 
  add_recipe(final_rf_recipe) %>% 
  add_model(final_rf_mod)
```

## EXECUTE

``` r
final_rf_res <-
  final_rf_workflow %>% 
  last_fit(splits)

saveRDS(final_rf_res, here::here("data", "final_rf_res.pro"))
```

This fitted workflow contains everything, including our final metrics
based on the test set. So, how did this model do on the test set? Was
the validation set a good estimate of future performance?

``` r
final_rf_res %>% 
  collect_metrics()
```

``` r
final_rf_res %>% 
  pluck(".workflow", 1) %>%   
  pull_workflow_fit() %>% 
  vip(num_features = 59)
```

``` r
final_rf_prep <- prep(final_rf_recipe)
boruta_obj <- final_rf_prep$steps[[4]]$res
```

``` r
plot(boruta_obj)
```

``` r
boruta_obj$ImpHistory
```
