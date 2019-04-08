##CALCULATION of RMSE for parameters

set.seed(123)

n.samples<-20
optimised.parameters<-matrix(data=NA, nrow=n.samples, ncol=5)
sample.id<-round(runif(n.samples, 1 , nrow(lopex1993metadata)))
e<-(as.matrix(lopex1993reflectance[sample.id, ]))
observed<-t(e)
counter<-0

extract.optim.param<-function(observed){
  counter<-counter+1
  obs<-as.numeric(observed)
  # print("couunt: ", counter)
  out.inversion<-invert_prospect(observed=obs, loglike = llfunction, settings = settings)
  par.opt<-MAP(out.inversion)
  p<-par.opt$parametersMAP[1:5]
  return(p)
}

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

ppangers<-pp
pplopex<-apply(observed, MARGIN=2, FUN=extract.optim.param)

#save(pplopex, sample.id, file="calibration of lopex.RData")
#save(pp, sample.id, file="calibration of angers.RData")

###
#Calculate RMSE for parameters
##
load ("C:/Users/ichas/Documents/master thesis/current topic/code-PROSPECT-sensitivity analysis and bayesian attempt/input_data/prospect.RData")
load(file="calibration of angers.RData")
ppangers<-pp
xi <- t(pp[,1:20])
x0 <- angers2003metadata[sample.id[1:20],c(1:5)]
RMSEper.a<-sqrt((colSums(((xi-x0)/x0)^2))/length(sample.id))*100
RMSE.a<-sqrt((colSums(((xi-x0))^2))/length(sample.id))
BIAS.a<-(colSums((xi-x0)))/length(sample.id)
bias.a<-matrix(BIAS, ncol=5, nrow = 20, byrow = TRUE)
SPEC.a<-sqrt((colSums((xi-x0-bias)^2))/20)



load(file="calibration of lopex.RData")
pp<-pplopex
xi <- t(pp[,1:20])
x0 <- lopex1993metadata[sample.id[1:20],c(1:5)]
RMSEper.l<-sqrt((colSums(((xi-x0)/x0)^2))/length(sample.id))*100
RMSE.l<-sqrt((colSums(((xi-x0))^2))/length(sample.id))
BIAS.l<-(colSums((xi-x0)))/length(sample.id)
bias.l<-matrix(BIAS, ncol=5, nrow = 20, byrow = TRUE)
SPEC.l<-sqrt((colSums((xi-x0-bias)^2))/20)

#
# Wavelengths RMSE---Plots of wavelength
#

lopex.index <- as.numeric(colnames(pplopex))
lopex.observed<-t(lopex1993reflectance[lopex.index,])
parameters<-t(pplopex)

inverse <- matrix(data=NA, nrow = 2051, ncol = length(lopex.index))
for (i in 1:length(lopex.index)){
  inverse[,i]<-prospect5(parameters[i, 1],parameters[i, 2],parameters[i, 3],parameters[i, 4],parameters[i,5])[1:2051,2]
}

diff.lopex<-inverse-lopex.observed

RMSE.wave.l <- sqrt((sum(diff.lopex^2))/(dim(diff.lopex)[1]*dim(diff.lopex)[2]))
BIAS.wave.l <- (sum(diff.lopex))/(dim(diff.lopex)[1]*dim(diff.lopex)[2])
sqrt((sum((diff.lopex-BIAS.wave.l)^2))/(dim(diff.lopex)[1]*dim(diff.lopex)[2]))


q<-apply(diff.lopex, 1, quantile)


plot(c(400:2450),q[3,], type="l", ylim = c(-0.05, 0.07), lwd=2, ylab="Observed-Modeled", xlab ="Wavelength", main="LOPEX")
lines(c(400:2450),q[2,], col="darkgray")
lines(c(400:2450),q[4,], col="darkgray")
abline(h=0, lty=2)
legend("topright", legend=c("median", "50% confidence interval"), col = c("black", "gray"), lty=1)

#########



angers.index <- as.numeric(colnames(ppangers))
angers.observed<-t(angers2003reflectance[angers.index,])
parameters<-t(ppangers)

inverse.a <- matrix(data=NA, nrow = 2051, ncol = length(angers.index))
for (i in 1:length(angers.index)){
  inverse.a[,i]<-prospect5(parameters[i, 1],parameters[i, 2],parameters[i, 3],parameters[i, 4],parameters[i,5])[1:2051,2]
}

diff.angers<-inverse.a-angers.observed
q.a<-apply(diff.angers, 1, quantile)

RMSE.wave.a <- sqrt((sum(diff.angers^2))/(dim(diff.angers)[1]*dim(diff.angers)[2]))
BIAS.wave.a <- (sum(diff.angers))/(dim(diff.angers)[1]*dim(diff.angers)[2])
sqrt((sum((diff.angers-BIAS.wave.a)^2))/(dim(diff.angers)[1]*dim(diff.angers)[2]))


plot(c(400:2450),q.a[3,], type="l", ylim = c(-0.04, 0.1), lwd=2, ylab="Observed-Modeled", xlab ="Wavelength", main="ANGERS")
lines(c(400:2450),q.a[2,], col="darkgray")
lines(c(400:2450),q.a[4,], col="darkgray")
abline(h=0, lty=2)
legend("topright", legend=c("median", "50% confidence interval"), col = c("black", "gray"), lty=1)
#########

