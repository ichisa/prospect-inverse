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


spectrum2 <- prospect5(N=angers2003metadata$N[1], 
                      Cab=angers2003metadata$Cab[1], 
                      Car = angers2003metadata$Car[1], 
                      Cw=angers2003metadata$Cw[1], 
                      Cm=angers2003metadata$Cm[1])
plot(spectrum2$Wavelength, spectrum2$Reflectance)

#Sensitivity analysis

param <- data.frame(name=colnames(lopex1993metadata[1:5]), 
                  def=colMeans(lopex1993metadata[,1:5]), 
                  min=sapply(lopex1993metadata[,1:5], FUN=min, simplify = "array"), 
                  max=sapply(lopex1993metadata[,1:5],FUN=max, simplify = "array"))

param <- param[,2:4]

save(spectrum2, param, file=paste0(data_dir, "angers_data_param.RData"))

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

out[]

WL.out<-param
WL.out[,c(1:3)] <- out[,,reflectance.out.index]
WL.out
#Absolute change
abs.change.low <- abs((WL.out$def-WL.out$min)/(param$def-param$min))
abs.change.upper <- abs((WL.out$def-WL.out$max)/(param$def-param$max))

abs.change.low
abs.change.upper

#Elasticity
per.change.low <- abs(((WL.out$def-WL.out$min)/WL.out$def)/((param$def-param$min)/param$def))*100
per.change.upper <- abs(((WL.out$def-WL.out$max)/WL.out$def)/((param$def-param$max)/param$def))*100

per.change.low
per.change.upper
#evaluate parameters FROM min to max values
seq1<-seq(refPars$min[2], refPars$max[2], length.out=20)

##
#Local sensitivity
##
r<-seq(51, 2051, 500)
r<-reflectance.out.index
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
