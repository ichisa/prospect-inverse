
##FUNCTIONS-needed

#' A)Check convergence of BayesianTools output using Gelman convergence
#' #########################################
converge_test <- function(samples, threshold = 1.1, use_CI = TRUE) {
  
  i <- ifelse(use_CI, 2, 1)
  gelman <- try(BayesianTools::gelmanDiagnostics(samples))
  if (class(gelman) == 'try-error') {
    message('Error trying to calculate gelman diagnostic. Assuming no convergence')
    return(FALSE)
  }
  gelman_vec <- gelman$psrf[,i]
  exceeds <- gelman_vec > threshold
  
  if (any(exceeds)) {
    exceeds_vec <- gelman_vec[exceeds]
    exceeds_char <- sprintf('%s: %.2f', names(exceeds_vec), exceeds_vec)
    exceeds_str <- paste(exceeds_char, collapse = '; ')
    message('The following parameters exceed threshold: ', exceeds_str)
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#'
#'B) create.bprior####
############################

##' Create priors for BayesianTools
##' 
##' prior.sel must contain the following columns:
##' `distn` -- String describing a distribution; e.g. `norm` for `dnorm`, `rnorm`, etc.
##' `parama`, `paramb` -- First and second parameters, respectively, of the corresponding distribution
##'
##' Optionally, `prior.sel` may also contain the following columns:
##'   * `param_name` -- Parameter name, which will be carried through to the prior object and sampler
##'   * `lower`, `upper` -- Lower and upper bounds, respectively. These can be leveraged by the BayesianTools samplers.
##'   * `best` -- Best guess for a parameter estimate. BayesianTools can also use this, though 
##'    
create.bprior <- function(prior.sel) {
  
  # Returns a function that calculates the density of the specified 
  # distribution given the parameters
  ddist_generator <- function(distn, a, b) {
    fun_string <- paste0('d', distn)
    f <- match.fun(fun_string)
    out <- function(x) f(x, a, b, log = TRUE)
    return(out)
  }
  
  # Returns a function that draws from the specified distribution with the 
  # specified parameters
  rdist_generator <- function(distn, a, b) {
    fun_string <- paste0('r', distn)
    f <- match.fun(fun_string)
    out <- function(n = 1) f(n, a, b)
    return(out)
  }
  
  # Create a list of density and random draw functions
  ddist_funs <- with(prior.sel, mapply(ddist_generator, distn, parama, paramb))
  rdist_funs <- with(prior.sel, mapply(rdist_generator, distn, parama, paramb))
  if ('param_name' %in% names(prior.sel)) {
    names(ddist_funs) <- names(rdist_funs) <- prior.sel[['param_name']]
  }
  
  # `mapply` statement returns
  density <- function(params) {
    dens_vec <- mapply(function(f, x) f(x), ddist_funs, params)   # Returns vector of log densities
    out <- sum(dens_vec)
    return(out)
  }
  
  # Returns vector of random draws
  sampler <- function() {
    out <- vapply(rdist_funs, function(f) f(), numeric(1))
    return(out)
  }
  
  # BayesianTools lower and upper bounds and best guess, if specified in data.frame
  lower <- NULL
  if ('lower' %in% names(prior.sel)) {
    lower <- prior.sel[['lower']]
  }
  upper <- NULL
  if ('upper' %in% names(prior.sel)) {
    upper <- prior.sel[['upper']]
  }
  best <- NULL
  if ('best' %in% names(prior.sel)) {
    best <- prior.sel[['best']]
  }
  
  out <- BayesianTools::createPrior(density = density, sampler = sampler, 
                                    lower = lower, upper = upper, best = best)
  return(out)
} #
