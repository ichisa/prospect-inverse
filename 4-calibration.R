
require(Rprospect)
library(lhs)
library(coda)
library(BayesianTools)
library(ExtDist)


data_dir <- ("input_data/")
load(paste0(data_dir, "angers_data_param.RData"))
load(file=paste0(data_dir, "prospect.RData"))#



par_ind<-ref_ind <- seq(from=1, to=length(angers2003reflectance)) 
ndat_ref <- length(ref_ind)

#Calculate sd of reflectance to include in the ll calculation
ref_ind <- 400 # index use to optimice
samp <- 1 #sample which parameters are used
sd_ref<- sd(angers2003reflectance[ref_ind,])
parind<-c(1:5)

#1- Likelihood function
likelihood <- function(pValues) {
  p <- param$def
  p[parind] <- pValues # new parameter values
  predicted<- prospect5(p[1], p[2],p[3],p[4],p[5])[,2]
  diff <- predicted[ref_ind]-angers2003reflectance[samp,ref_ind]  #llvalues <- sum(dnorm(predicted$GPP, mean = s1$GPPobs, sd = p[31], log=T)) ###   llvalues <- sum(dnorm(diff_GPP, sd = p[31], log=T))
  llvalues <- sum(dexp(abs(diff),rate = 1/(sd_ref*predicted[ref_ind]),log=T))
  return(llvalues)
}



#2- Prior
prior <- createUniformPrior(lower = param$min, upper = param$max)

#=Bayesian set up=#
param$names<-row.names(param)
BSpreles <- createBayesianSetup(likelihood, prior, best = param$def, names = param$names, parallel = F)

#=Run the MCMC with three chains=#

settings <- list(iterations = 3e3, optimize=F, nrChains = 2)
chainDE <- runMCMC(BSpreles, sampler="DEzs", settings = settings)
par.opt<-MAP(chainDE) #gets the optimized maximum value for the parameters

#=Check convergence -- Trace plots=#

pdf("check-plots.pdf")
tracePlot(chainDE, parametersOnly = TRUE, start = 1, whichParameters = 1:5)

marginalPlot(chainDE, scale = T, best = T, start = 5000)
correlationPlot(chainDE, parametersOnly = TRUE, start = 2000)
dev.off()

