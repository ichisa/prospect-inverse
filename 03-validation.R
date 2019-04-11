

#loads prospect data
load ("prospect.RData")
observed<-as.numeric(angers2003reflectance[2,])

#Test the funcions 
inversion <- invert_prospect( observed = observed, loglike = llfunction, settings = settings) 


pdf("check-plots.pdf")
tracePlot(samples, parametersOnly = TRUE, start = 1, whichParameters = 1:5)
marginalPlot(samples, scale = T, best = T, start = 5000)
correlationPlot(samples, parametersOnly = TRUE, start = 2000)
dev.off()

save(inversion,samples, observed, model, file="C:/Users/ichas/Documents/master thesis/current topic/code-bayesian-inversion-prospect/inversion_example.RData")



#Run 20 different inversions to validate model
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


##Save the inversions
#ppangers<-pp
#pplopex<-apply(observed, MARGIN=2, FUN=extract.optim.param)
###################################################################################################################
#save(pplopex, sample.id, file="calibration of lopex.RData")
#save(pp, sample.id, file="calibration of angers.RData")

###
#Calculate RMSE for parameters
##
#load observations and resuts of inversion

load ("prospect.RData")
load(file="calibration of angers.RData")
ppangers<-pp
sample.id.angers<-sample.id
rm(pp, sample.id)
load(file="calibration of lopex.RData")

xi <- t(ppangers[,1:20])
x0 <- angers2003metadata[sample.id[1:20],c(1:5)]

RMSEper.a<-sqrt((colSums(((xi-x0)/x0)^2))/length(sample.id))*100
RMSE.a<-sqrt((colSums(((xi-x0))^2))/length(sample.id))
BIAS.a<-(colSums((xi-x0)))/length(sample.id)
bias.a<-matrix(BIAS.a, ncol=5, nrow = 20, byrow = TRUE)
SPEC.a<-sqrt((colSums((xi-x0-bias.a)^2))/20)



load(file="calibration of lopex.RData")
lopex.index<-sample.id
pp<-pplopex
xi <- t(pp[,1:20])
x0 <- lopex1993metadata[sample.id[1:20],c(1:5)]
RMSE.l<-sqrt((colSums(((xi-x0))^2))/length(sample.id))
RMSEper.l<-sqrt((colSums(((xi-x0)/x0)^2))/length(sample.id))*100
RMSE.l<-sqrt((colSums(((xi-x0))^2))/length(sample.id))
BIAS.l<-(colSums((xi-x0)))/length(sample.id)
bias.l<-matrix(BIAS.l, ncol=5, nrow = 20, byrow = TRUE)
SPEC.l<-sqrt((colSums((xi-x0-bias.l)^2))/20)
lopex.validation<-data.frame(parameter=colnames(xi),RMSE=RMSE.l, RMSEper=RMSEper.l, BIAS=BIAS.l, SPEC=SPEC.l)
lopex.validation


# Calculate errors for the wavelengths
# Wavelengths RMSE---Plots of wavelength
#

lopex.index <- as.numeric(colnames(pplopex))

lopex.observed<-t(lopex1993reflectance[lopex.index,])
parameters<-t(pplopex)

#run the model for the result of the inversion
inverse <- matrix(data=NA, nrow = 2051, ncol = length(lopex.index))
for (i in 1:length(lopex.index)){
  inverse[,i]<-prospect5(parameters[i, 1],parameters[i, 2],parameters[i, 3],parameters[i, 4],parameters[i,5])[1:2051,2]
}

diff.lopex<-inverse-lopex.observed


RMSE.wave.l.vis <- sqrt((sum(diff.lopex[1:401, ]^2))/(dim(diff.lopex[1:401, ])[1]*dim(diff.lopex[1:401, ])[2]))
BIAS.wave.l.vis <- (sum(diff.lopex[1:401, ]))/(dim(diff.lopex[1:401, ])[1]*dim(diff.lopex[1:401, ])[2])
SPEC<-sqrt((sum((diff.lopex-BIAS.wave.l)^2))/(dim(diff.lopex)[1]*dim(diff.lopex)[2]))

RMSE.wave.l.ir <- sqrt((sum(diff.lopex[401:2051, ]^2))/(dim(diff.lopex[401:2051, ])[1]*dim(diff.lopex[401:2051, ])[2]))
BIAS.wave.l.ir <- (sum(diff.lopex[401:2051, ]))/(dim(diff.lopex[401:2051, ])[1]*dim(diff.lopex[401:2051, ])[2])
SPEC<-sqrt((sum((diff.lopex-BIAS.wave.l.ir)^2))/(dim(diff.lopex)[1]*dim(diff.lopex)[2]))


#calculate conf interval
q<-apply(diff.lopex, 1, quantile, probs=c(0.1, 0.25, 0.5, 0.75, 0.9))

#plot
plot(c(400:2450),q[3,], type="l", ylim = c(-0.1, 0.1), lwd=2, ylab="Observed-Inversion", xlab ="Wavelength", main="LOPEX")
lines(c(400:2450),q[2,], col="darkgray", lwd=1.5)
lines(c(400:2450),q[4,], col="darkgray", lwd=1.5)
lines(c(400:2450),q[1,], col="lightgray", lwd=1.5, lty=3)
lines(c(400:2450),q[5,], col="lightgray",lwd=1.5, lty=3)

abline(h=0, lty=2)
legend("topright", legend=c("median", "50% confidence interval", "90% confidence interval"), col = c("black", "gray", "lightgray"), lty=1)







######### Same analysis for angers

angers.index <- as.numeric(colnames(ppangers))
angers.observed<-t(angers2003reflectance[angers.index,])
parameters<-t(ppangers)

inverse.a <- matrix(data=NA, nrow = 2051, ncol = length(angers.index))
for (i in 1:length(angers.index)){
  inverse.a[,i]<-prospect5(parameters[i, 1],parameters[i, 2],parameters[i, 3],parameters[i, 4],parameters[i,5])[1:2051,2]
}

diff.angers<-inverse.a-angers.observed
q.a<-apply(diff.angers, 1, quantile, probs=c(0.1, 0.25, 0.5, 0.75, 0.9))

RMSE.wave.a <- sqrt((sum(diff.angers^2))/(dim(diff.angers)[1]*dim(diff.angers)[2]))
BIAS.wave.a <- (sum(diff.angers))/(dim(diff.angers)[1]*dim(diff.angers)[2])
sqrt((sum((diff.angers-BIAS.wave.a)^2))/(dim(diff.angers)[1]*dim(diff.angers)[2]))


plot(c(400:2450),q.a[3,], type="l", ylim = c(-0.1, 0.1), lwd=2, ylab="Observed-Invevrse", xlab ="Wavelength", main="ANGERS")
lines(c(400:2450),q.a[2,], col="darkgray")
lines(c(400:2450),q.a[4,], col="darkgray")
lines(c(400:2450),q.a[1,], col="lightgray", lwd=1.5, lty=3)
lines(c(400:2450),q.a[5,], col="lightgray",lwd=1.5, lty=3)

abline(h=0, lty=2)
legend("topright", legend=c("median", "50% confidence interval"), col = c("black", "gray"), lty=1)
#########


################plot observed - model but for prospect 5
angers.index<-sample.id.angers
predicted.a <- matrix(data=NA, nrow = 2051, ncol = length(angers.index))

dat <- as.matrix(angers2003metadata[angers.index,1:5])

for (i in 1:length(angers.index)){
 
  predicted.a[,i]<-prospect5(dat[i, 1],dat[i, 2],dat[i, 3],dat[i, 4],dat[i,5])[1:2051,2]
}

diff.angers.2<-predicted.a-angers.observed


#calculate conf interval
a2<-apply(diff.angers.2, 1, quantile, probs=c(0.1, 0.25, 0.5, 0.75, 0.9))

#plot
plot(c(400:2450),a2[3,], type="l", ylim = c(-0.04, 0.04), lwd=2, ylab="Observed-Modeled", xlab ="Wavelength", main="ANGERS")
lines(c(400:2450),a2[2,], col="darkgray", lwd=1.5)
lines(c(400:2450),a2[4,], col="darkgray", lwd=1.5)
lines(c(400:2450),a2[1,], col="lightgray", lwd=1.5, lty=3)
lines(c(400:2450),a2[5,], col="lightgray",lwd=1.5, lty=3)

abline(h=0, lty=2)
legend("topright", legend=c("median", "50% confidence interval", "90% confidence interval"), col = c("black", "gray", "lightgray"), lty=1)

######################same for LOPEX

lopex.index<-sample.id
predicted.l <- matrix(data=NA, nrow = 2051, ncol = length(lopex.index))

dat <- as.matrix(lopex1993metadata[lopex.index,1:5])

for (i in 1:length(lopex.index)){
  
  predicted.l[,i]<-prospect5(dat[i, 1],dat[i, 2],dat[i, 3],dat[i, 4],dat[i,5])[1:2051,2]
}

diff.lopex.2<-predicted.l-lopex.observed


#calculate conf interval
l2<-apply(diff.lopex.2, 1, quantile, probs=c(0.1, 0.25, 0.5, 0.75, 0.9))

#plot
plot(c(400:2450),l2[3,], type="l", ylim = c(-0.06, 0.06), lwd=2, ylab="Observed-Modeled", xlab ="Wavelength", main="LOPEX")
lines(c(400:2450),l2[2,], col="darkgray", lwd=1.5)
lines(c(400:2450),l2[4,], col="darkgray", lwd=1.5)
lines(c(400:2450),l2[1,], col="lightgray", lwd=1.5, lty=3)
lines(c(400:2450),l2[5,], col="lightgray",lwd=1.5, lty=3)

abline(h=0, lty=2)
legend("topright", legend=c("median", "50% confidence interval", "90% confidence interval"), col = c("black", "gray", "lightgray"), lty=1)


#Plots of parameters measured vs predicted###########################################################################


observed.lopex<-lopex1993metadata[lopex.index, 1:5]
x<- as.vector(observed.lopex[,1])
y<- as.numeric(pplopex[1,])
plot(x~y)

