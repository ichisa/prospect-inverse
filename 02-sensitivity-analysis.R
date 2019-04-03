#Sensitivity analysis for PROSPECT5

require (dplyr)
require(Rprospect)
library(ggplot2)
require(sensitivity)

#############################
###'()Local sensitivity()'###
#############################

#Parameter Values from literature Feret 2008
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

}#  end loop local sensitivity


#Plots#################

jpeg("rplot.png", width = 700, height = 550)

par(mfrow=c(3,2), mar=c(3, 3, 1.5, 0.6), cex=0.9, cex.axis=0.85, mgp=c(2,1,0))
for (p in 1:5){
  plot(out[p,3,], type = "l", col="darkgreen",xlab="Wavelength", ylab="Reflectance", 
       lwd=3, ylim = c(0, 0.8))
  lines(out[p,2,], col="lightgreen", lwd=3)
  lines(out[p,1,], col="limegreen", lwd=3)
  text(50, 0.65, rownames(param)[p], cex=1.3)
}
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("bottomleft",legend =c("max"," def"," mim"), col=c("darkgreen", "limegreen","lightgreen" ), lwd=3.5, cex=1.2, bty="n")
dev.off()
dev.off()
###################
par(mfrow=c(3,2), mar=c(3, 3, 1.5, 0.6), cex=0.9, cex.axis=0.85, mgp=c(2,1,0))
for (p in 1:5){
  plot(wvl, per.change.low[,p], type="l", col="darkred", xlab="Wavelength", ylab="Elasticity(%)", 
       lwd=3)
  lines(wvl, per.change.upper[,p], type="l", lwd=3, col="orange")
}
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("bottomleft",legend =c("Change from defined to lower","Change from defined to upper" ), col=c("orange", "darkred"), lwd=3.5, cex=1.2, bty="n")
dev.off

# PArameter importance for different wavelength

r<-seq(51, 2051, 200)
n<-20

local.df <- data.frame(result=numeric(20*5*length(r)), 
                       par.value=numeric(20*5*length(r)), 
                       parameter=rep(NA, 20*5*length(r)),
                       Wavelength=numeric(20*5*length(r)))

count<-1
for (ri in 1:length(r) ){
  wvl<-r[ri]
  for (k in c(1:5)){ #Look through pars
    results<-numeric(20)
    seq1<-seq(param$min[k], param$max[k], length.out=20)
    for (i in 1:20){#loop through the sequence
      parSel<-seq1[i]
      local.df$parameter[count] <- rownames(param)[k]
      out.e <- sensitivityTarget(parSel, k)
      local.df$result[count] <- out.e[wvl]
      local.df$par.value[count] <- seq1[i]
      local.df$Wavelength[count] <- wvl
      count<-count+1
      
    }
  }
}

local.df$Wavelength<-local.df$Wavelength+399
local.df$parameter<-ordered(local.df$parameter, levels=c("N", "Cab", "Car", "Cw", "Cm"))#get right oreder in plot

 local.df %>% ggplot(aes(x=par.value, y=result, group=Wavelength, col=Wavelength))+
    geom_line( size=1)+
    facet_wrap(~parameter, scales="free")+
    xlab("Parameter value") + ylab("Predicted Reflectance")+
    scale_color_gradient(low="green", high="blue")+
    theme_bw()


#########################################
#'()Global sensitivity-morris screning()#
#########################################

par <- param

####Morris sensitivity analysisi for each wavelength###########################
 
param$names <- c("N", "Cab", "Car", "Cw", "Cm")
ws <- seq(from=1, to=2101, by=200)

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
return(sigma)
}

ptm <- proc.time()
sigma<-sapply(ws, run_morris)
sigma
t2<-proc.time()
t<-ptm$user-t2$user

par(mfrow=c(1,1))
plot(res_morris)

save(musar,sigma ,ws,file="morris-sensitivity-analysis.RData")

