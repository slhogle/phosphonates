sim_pa <- function(x, B){
  # init output
  mod <- x
  M <- mod$M
  W <- mod$W
  ymin <- ymax <- rep(NA, length(M))
  sims <- matrix(NA, nrow = B, ncol = length(W))
  newdat <- mod$dat
  
  simiter <- function(mod, W, newdat){
    library(corncob)
    # simulate init W by drawing from model
    sim <- simulate(mod, nsim = length(W))
    # replace old W with simulated W
    newdat$W <- sim
    # refit model using simulate W
    refit <-
      suppressWarnings(try(corncob::bbdml(mod$formula, phi.formula = mod$phi.formula,
          link = mod$link, phi.link = mod$phi.link, inits = mod$inits, data = newdat), silent = TRUE))
    
    if ("try-error" %in% c(class(refit))) {
      return(list(NA, NA, NA))
    }
    
    # pull simulated observations from model
    sim_obs <- simulate(refit, nsim = length(W))
    
    # pull new coefficients
    sim_mu <- refit$b.mu
    sim_phi <- refit$b.phi
    
    # output as list
    list(sim_obs, sim_mu, sim_phi)
  }
  
  # init parallelization framework
  cores=detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  
  # run parallel simulations
  sims <- foreach(i=1:B, .combine=rbind) %dopar% {
    tmpmat <- simiter(mod, W, newdat)
    # this is the same to doing an rbind with sims and the tmpmat
    tmpmat
  }
  
  # close parallelization framework
  stopCluster(cl)
  
  # output the list
  return(sims)
}

########################################

sim_summary <- function(x, sims){
  
  # original model para
  mod <- x
  M   <- mod$M
  W   <- mod$W
  mu  <- mod$b.mu
  phi <- mod$b.phi
  
  # formatting obs
  obs.predint <- apply(do.call(rbind, sims[[1]]), 2, stats::quantile, c(.025, .975), TRUE)
  
  # formatting mu
  mu.predint <- apply(do.call(rbind, sims[[2]]), 2, stats::quantile, c(.025, .975), TRUE)
  
  # formatting phi
  phi.predint <- apply(do.call(rbind, sims[[3]]), 2, stats::quantile, c(.025, .975), TRUE)
  
  df.obs <- data.frame(
    sample = mod$dat[,1],
    tot = M,
    W = W,
    yminW = obs.predint[1, ],
    ymaxW = obs.predint[2, ], 
    RA = W / M,
    yminRA = obs.predint[1, ] / M,
    ymaxRA = obs.predint[2, ] / M
  )
  
  df.mu <- data.frame(
    Mu = mu, 
    yminMu = mu.predint[1, ],
    ymaxMu = mu.predint[2, ]
  )
  
  df.phi <- data.frame(
    Phi = phi,
    yminPhi = phi.predint[1, ],
    ymaxPhi = phi.predint[2, ]
  )
  
  list(df.obs, df.mu, df.phi)
}

########################################

contrastmod <- function(invcont, mod_input, levelvector){
  cmat <- solve(invcont)[,-1]
  
  X <- model.matrix(~region,
                    data=mod_input,
                    contrasts=list(region=cmat))
  
  mod_input_1 <- with(mod_input, data.frame(W, M, region, X[,-1]))
  
  bbdml(formula     = as.formula(paste("cbind(W, M)", paste(paste0("region", levelvector), collapse=" + "), sep=" ~ ")),
        phi.formula = as.formula(paste("cbind(W, M)", paste(paste0("region", levelvector), collapse=" + "), sep=" ~ ")),
        link        = "logit",
        phi.link    = "logit",
        data        = mod_input_1)
}

########################################

bootLRT <- function(mod, mod_null) {
  # Simulate n samples from the model fit under the null
  library(corncob)
  newW <- simulate(object = mod_null, nsim = nrow(mod_null$dat))
  
  newdf <- mod$dat
  newdf$W <- newW
  
  # Refit models using simulated data
  newout_null <- suppressWarnings(try(bbdml(formula = mod_null$formula,
                                            phi.formula = mod_null$phi.formula,
                                            link = mod_null$link, phi.link = mod_null$phi.link,
                                            data = newdf, inits = mod_null$inits), silent = TRUE))
  newout_alt <- suppressWarnings(try(bbdml(formula = mod$formula,
                                           phi.formula = mod$phi.formula,
                                           link = mod$link, phi.link = mod$phi.link,
                                           data = newdf, inits = mod$inits), silent = TRUE))
  
  if ("try-error" %in% c(class(newout_null), class(newout_alt))) {
    return(NA)
  }
  
  test.stat <- 2 * abs(newout_alt$logL - newout_null$logL)
  return(test.stat)
}

########################################

PBLRTfun <- function(mod, mod_null, B, bootfun){

  sims <- matrix(NA, nrow = B, ncol = 1)
  
  # init parallelization framework
  cores=detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  
  # run parallel simulations
  sims <- foreach(i=1:B, .combine=rbind) %dopar% {
    tmpmat <- bootfun(mod, mod_null)
    # this is the same to doing an rbind with sims and the tmpmat
    tmpmat
  }
  
  # close parallelization framework
  stopCluster(cl)
  
  # output the list
  sims
}

########################################

makenull <- function(null_form, full_mod, input_data){
  
  null.m <- bbdml(formula = null_form, phi.formula = null_form, link="logit",
                  phi.link="logit", data = input_data)
  
  null.m
}

########################################

mypbLRT <- function(mod, mod_null, sims){ #, input_data
  
  t.observed <- 2 * abs(mod$logL - mod_null$logL)

  pboot.m <- mean(sims >= t.observed, na.rm=T)
  pboot.o <- (1 + sum(stats::na.omit(sims) >= t.observed)) / (length(stats::na.omit(sims)) + 1)
  
  lrt <- lrtest(mod_null = mod_null, mod = mod)
  
  out <- list(pboot.m=pboot.m, pboot.o=pboot.o, lrt=lrt)
  out
}

########################################

remove_terms <- function(form, term) {
  fterms <- terms(form)
  fac <- attr(fterms, "factors")
  idx <- which(as.logical(fac[term, ]))
  new_fterms <- drop.terms(fterms, dropx = idx, keep.response = TRUE)
  return(formula(new_fterms))
}


