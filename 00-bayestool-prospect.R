
library(Rprospect)
library(BayesianTools)
library(readr)


#####################################

#Loads needed functions
source('bayes_functions.R')
rm(lopex1993metadata, lopex1993reflectance)

nparams<-5
model<-prospect5 # defines model from library Rprospect

#Define spectrum to invert
observed<-as.numeric(angers2003reflectance[2,])
#observed[2052:2101]<-observed[2051]
default_parameters <- c(N=1.5, Cab=40, Car=8, Cw=0.01, Cm=0.009) #Valores que vienen de feret 2008biblio

#'#######################
#'1) Likelihood Function#
#'#######################


fail_ll <- -1e10
n_obs <- length(observed)
llfunction <- function(x) {
  rtm_params <- x[seq_len(nparams)]
  rsd <- x[nparams + 1]
  mod <- model(rtm_params[1],rtm_params[2],rtm_params[3],rtm_params[4],rtm_params[5])[1:2051,2] ##Run prospect- selects indices up to 2051 because measured data goes to 2051
  
  err <- mod - observed 
  ss <- sum(err * err) 
  sigma2 <- rsd * rsd #residual
  
  ll <- -0.5 * (n_obs * log(sigma2) + ss / sigma2)
  
  return(ll)
}

#'####################################
#'2)default parameters to create prior
#'#####################################  


defpram <- default_parameters 
col_names <- c('param_name', 'distn', 'parama', 'paramb', 'lower', 'best', 'upper')
prior_list <- list(#values from Feret 2008
  N = list('N', 'norm', 1.4, 0.8, 1, 1.4, 4),
  Cab = list('Cab', 'lnorm', log(40), 0.9, 0, 40, 110),
  Car = list('Car', 'lnorm', log(10), 1.1, 0, 10, 30),
  Cw = list('Cw', 'lnorm', log(0.01), 1, 0, 0.01, 0.5),
  Cm = list('Cm', 'lnorm', log(0.009), 1, 0, 0.009, 0.16),
  residual = list('residual', 'lnorm', log(0.001), 2.5, 0, 0.000001, 10000)
)

prior_df_all <- do.call(rbind.data.frame, prior_list)
colnames(prior_df_all) <- col_names
prior <- create.bprior(prior_df_all) # funtions parsed at the beginning
rm(prior_df_all, prior_list)
##########################################################################
#'3) Bayesian inversion using BayesianTools package
#'

#Define the settings
loglike<-llfunction
settings <- list(
  common = list(),
  init = list(iterations = 10000),
  loop = list(iterations = 2000),
  other = list(
    sampler = 'DEzs',
    min_samp = 5000,
    max_iter = 1e6,
    threshold = 1.1
  )
)

###Function

invert_prospect<-function(observed, loglike, settings){
  
  min_samp <- settings[['other']][['min_samp']]
  max_iter <- settings[['other']][['max_iter']]
  threshold <- settings[['other']][['threshold']]
  
  test_samp <- prior$sampler()
  param_names <- names(test_samp)
  
  
  
  
  setup <- BayesianTools::createBayesianSetup(
    likelihood = loglike,
    prior = prior,
    names = param_names
  )
  
  
  init_settings <- modifyList(settings[['common']], settings[['init']])
  stop_iter <- init_settings[["iterations"]]
  
  message('Running initial ', stop_iter, ' iterations.')
  samples <- BayesianTools::runMCMC(
    bayesianSetup = setup,
    sampler = settings[['other']][['sampler']],
    settings = init_settings
  )
  
  #
  #'4)Check convergence and go on if not converged
  
  converged <- converge_test(samples = samples, threshold = threshold)
  
  loop_settings <- modifyList(settings[['common']], settings[['loop']])
  
  next_iter <- loop_settings[['iterations']]
  
  
  
  while (!(converged)) {
    start_iter <- stop_iter + 1
    stop_iter <- stop_iter + next_iter
    if (start_iter > max_iter) {
      warning('Next start iteration (', start_iter, ') greater than maximum iteration count (', max_iter, ') ',
              'but convergence has not been achieved. ',
              'Terminating sampling and returning results as is.')
      break
    }
    message('Running ', next_iter, ' more iterations (', start_iter, ' to ', stop_iter, ').')
    samples <- BayesianTools::runMCMC(samples, sampler = sampler, settings = loop_settings)
    converged <- converge_test(samples = samples, threshold = threshold)
    
    if (converged) {
      coda_samples <- BayesianTools::getSample(samples, coda = TRUE)
    }
  }
  #return(list(samples, coda_samples))
  return(samples)
}#end functions invert prospect


#Test the funcions 

inversion <- invert_prospect( observed = observed, loglike = llfunction, settings = settings) 


pdf("check-plots.pdf")
tracePlot(samples, parametersOnly = TRUE, start = 1, whichParameters = 1:5)
marginalPlot(samples, scale = T, best = T, start = 5000)
correlationPlot(samples, parametersOnly = TRUE, start = 2000)
dev.off()

save(inversion,samples, observed, model, file="C:/Users/ichas/Documents/master thesis/current topic/code-bayesian-inversion-prospect/inversion_example.RData")



