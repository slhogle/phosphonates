.libPaths(c("/projappl/project_2001175/project_rpackages", .libPaths()))
library(tidyverse)
library(tidymodels)
library(ranger)
library(parallel)
library(doParallel)

# cores
ncores <- 20 # number of cores requested 
doParallel::registerDoParallel(cores=20)

groups <- c("pro", "sar11", "bactarc")

for (group in groups) {
  myinfile <- paste("PEPMreads2forest.", group, sep="")
  myoutfile <- paste("rand.forest.mantune.par.", group, sep="")
  
  # set up train/test partition
  PEPMreads2forest <- readRDS(myinfile)
  pepm.split <- initial_split(PEPMreads2forest, prop = 0.8)
  pepm.split.train <- training(pepm.split)
  pepm.split.test <- testing(pepm.split)
  pepm.split.train.crossval <- vfold_cv(pepm.split.train)
  
  # define recipe for data preprocessing
  pepm.recipe <- recipe(RA ~ ., data = PEPMreads2forest) %>% 
    update_role(sampleID, new_role = "ID") %>% 
    step_normalize(all_predictors(), -dcm.layer, -dcm.max, -hdb_cl, -ocean)
  
  # define model
  pepm.rand.forest.model <- rand_forest() %>%
    set_args(mtry = tune(), 
             min_n = tune(), 
             trees = 2000) %>%
    set_engine("ranger",
               importance = "permutation",
               splitrule = "beta",
               replace = FALSE) %>%
    set_mode("regression")
  
  # define model workflow
  pepm.rand.forest.workflow <- workflow() %>%
    add_recipe(pepm.recipe) %>%
    add_model(pepm.rand.forest.model)
  
  # create tuning grid for parameters mtry and min_n
  pepm.rand.forest.grid <- grid_regular(
    mtry(range = c(10, 25)),
    min_n(range = c(1, 3)),
    levels = 7)
  
  # run manual tuning of model
  pepm.rand.forest.mantune <- pepm.rand.forest.workflow %>%
    tune_grid(resamples = pepm.split.train.crossval,
              grid = pepm.rand.forest.grid,
              metrics = metric_set(rmse, mae, ccc, rsq),
              control = control_grid(verbose = TRUE,
                                     allow_par = TRUE,
                                     extract = NULL,
                                     save_pred = FALSE,
                                     pkgs = NULL))
  
  # save results
  saveRDS(pepm.rand.forest.mantune, myoutfile)
  
}