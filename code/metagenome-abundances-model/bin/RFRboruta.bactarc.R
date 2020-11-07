.libPaths(c("/projappl/project_2001175/project_rpackages", .libPaths()))
library(tidyverse)
library(tidymodels)
library(ranger)
library(Boruta)

input <- "PEPMreads2forest.bactarc"
trees <- 750
m.try <- 10
min.n <- 15
output <- "PEPMboruta.bactarc"

PEPMreads2forest <- readRDS(input)

# define recipe for data preprocessing
pepm.recipe <- recipe(RA ~ ., data = PEPMreads2forest) %>% 
  step_normalize(all_predictors(), -all_nominal())

# extract the normalized values
pepm.preprocessed <- pepm.recipe %>%
  prep() %>%
  juice()

# custom function for getting variable importance
getImpBetaZ <- function(x, y, ntree=trees, num.trees=ntree, ...){
  x$shadow.Boruta.decision <- y
  ranger::ranger(data=x, 
                 dependent.variable.name = "shadow.Boruta.decision",
                 num.trees = num.trees, 
                 mtry = m.try,             # values determined from tune step
                 min.node.size = min.n,     # values determined from tune step
                 importance = "permutation",
                 splitrule = "beta",
                 scale.permutation.importance = TRUE,
                 write.forest = FALSE,...)$variable.importance
}

borutaPEPM <- Boruta(RA ~ .,
                     getImp=getImpBetaZ,
                     data=pepm.preprocessed,
                     doTrace=2)

saveRDS(borutaPEPM, output)