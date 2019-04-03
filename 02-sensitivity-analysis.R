###
#This code loads the input data for the model
###

require(readr)
require (dplyr)
require(Rprospect)
library(ggplot2)
require(sensitivity)
#load prospectd data
data_dir <- ("input_data/")
load(file=paste0(data_dir, "prospect.RData"))#

#Find wavelength that hast highest variability to work with it
#calculates variability of all wavelengths
ref.var <- sapply(angers2003reflectance[,],FUN=var, simplify = "array")
#select wavelength in the spectrum with higher variability
reflectance.out.index <- which(ref.var==max(ref.var))

#Values from literature Feret 2008
N<-c(1.6, 1, 4)
Cab<-c(45, 6, 100)
Car<-c(10, 2, 60)
Cw<-c(0.01,  0.001, 0.05)
Cm<-c(0.005, 0.001, 0.015)
param<-data.frame(t(cbind(N, Cab, Car, Cw, Cm)))
names(param)<- c("def", "min", "max")

#Return one value for each parameter
#Sensitivity tfunction
sensitivityTarget <- function(par.new, parSel){
  refPars$def[parSel] <- par.new
  predicted <- prospect5(refPars$def[1], refPars$def[2],refPars$def[3],refPars$def[4],refPars$def[5])[,2]
  return(predicted)
}

#Calculates output changing one parameter at the time for maximum and minimum parameters
out <-array(data=NA, c(5,3,2101))
refPars <- param
for (i in 1:nrow(out)) {
  for (j in 1:3) {
    par <- param$def
    par.new <- refPars[i, j] 
    parSel <- i
    e <- sensitivityTarget(par.new , i)
    out[i,j,] <- e#saves
  }
}

#Loop for wavelength values 
wvl<-seq(1, 2101, 10)
abs.change.low<-abs.change.upper<-matrix(NA,length(wvl),5)
per.change.low<-per.change.upper<-matrix(NA,length(wvl),5)
WL.out<-param

for (w in c(1:length(wvl))){

WL.out[,c(1:3)] <- out[,,wvl[w]]
#Absolute change
abs.change.low[w,] <- abs((WL.out$def-WL.out$min)/(param$def-param$min))
abs.change.upper[w,] <- abs((WL.out$def-WL.out$max)/(param$def-param$max))
#Elasticity
per.change.low[w,]<- abs(((WL.out$def-WL.out$min)/WL.out$def)/((param$def-param$min)/param$def))*100
per.change.upper[w,] <- abs(((WL.out$def-WL.out$max)/WL.out$def)/((param$def-param$max)/param$def))*100

}#end loop local sensitivity


#Plots
jpeg("rplot.png", width = 700, height = 550)

par(mfrow=c(3,2), mar=c(3, 3, 1.5, 0.6), cex=0.9, cex.axis=0.85, mgp=c(2,1,0))
for (p in 1:5){
plot(out[p,3,], type = "l", col="darkgreen",xlab="Reflectance", ylab="Wavelength", lwd=3, ylim = c(0, 0.8))
lines(out[p,2,], col="lightgreen", lwd=3)
lines(out[p,1,], col="limegreen", lwd=3)
text(50, 0.65, rownames(param)[p], cex=1.3)
}
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("bottomleft",legend =c("max"," def"," mim"), col=c("darkgreen", "limegreen","lightgreen" ), lwd=3.5, cex=1.2, bty="n")
dev.off()
#evaluate parameters FROM min to max values


##
#Local sensitivity
##
reflectance.out.index<-300
r<-seq(51, 2051, 500)
for (reflectance.index in r ){
local.df <- data.frame(result=numeric(20*5), par.value=numeric(20*5), parameter=rep(NA, 20*5) )
count <- 1
for (k in c(1:5)){ #Look through pars
  results<-numeric(20)
  seq1<-seq(param$min[k], param$max[k], length.out=20)
  
  for (i in 1:20){#loop through the sequence
    parSel<-seq1[i]
    local.df$parameter[count] <- param[k,1]
    out.e <- sensitivityTarget(parSel, k)
    local.df$result[count] <- out.e[reflectance.index]
    local.df$par.value[count] <- seq1[i]
    count=count+1
  }
}

print(
local.df %>% ggplot(aes(x=par.value, y=result, group=parameter))+
  geom_line( size=1.2)+
  facet_wrap(~parameter, scales="free")+
  xlab("Parameter value") + ylab("Predicted")+ ggtitle(reflectance.out.index)+
  theme_bw() 
)
}

#############################################################################
##
#Global sensitivity-morris screning
##

par <- param

#par_sel <- 3
sensitivityTarget2 <- function(mat){
  result = numeric(nrow(mat))
  for (i in 1:nrow(mat)){
    par_temp <-param$def
    par_temp = as.vector(mat[i,])
    predicted <- prospect5(par_temp[1], par_temp[2],par_temp[3],par_temp[4],par_temp[5])[,2]
    result[i] = predicted[1000]
  }
  return(result)
}

param$names <- c("N", "Cab", "Car", "Cw", "Cm")
res_morris <-morris(model=sensitivityTarget2, factors = param$name,
                    r=500, design=list(type="oat", levels=10, grid.jump=1), binf=par$min, bsup=par$max, scale=T)

par(mfrow=c(1,1))
plot(res_morris)


####Morris sensitivity analysisi for each wavelength###########################
 
param$names <- c("N", "Cab", "Car", "Cw", "Cm")

ws <- seq(from=1, to=2101, by=1000)

morris_mu<-matrix(NA, nrow = length(ws), ncol = 5)
morris_mustar<-matrix(NA, nrow = length(ws), ncol = 5)
morris_sigma <- matrix(NA, nrow = length(ws), ncol = 5)


#par_sel <- 3
sensitivityTarget2 <- function(mat){
  result = numeric(nrow(mat))
  for (i in 1:nrow(mat)){
    par_temp <-param$def
    par_temp = as.vector(mat[i,])
    predicted <- prospect5(par_temp[1], par_temp[2],par_temp[3],par_temp[4],par_temp[5])[wvl,2]
    result[i] = predicted
  }
  return(result)
}


run_morris<-function(wvl){
  res_morris <-morris(model=sensitivityTarget2, 
                    factors = param$name,
                    r=500, design=list(type="oat", levels=10, grid.jump=1), 
                    binf=par$min, bsup=par$max, scale=T)

  mu<- apply(res_morris$ee, 2, mean)
  mustar <- apply(res_morris$ee, 2, function(x) mean(abs(x)))
  sigma <- apply(res_morris$ee, 2, sd)
return(list(mu, musigma, mustar))
}

ptm <- proc.time()
mustar<-sapply(1:200, run_morris)
t2<-proc.time()
t<-ptm$user-t2$user

par(mfrow=c(1,1))
plot(res_morris)

save(mu, file="101-200 mu sensitivity morris")

