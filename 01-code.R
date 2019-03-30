#########################
#This code only loads the data, converts it into the specific hsdar library format and plot the 
#wavelengths using hsdar library
#########################
# Not anymore used because not the HSDAR libray will be used at the beginning


library(readr)
library (dplyr)
library(hsdar)

#Loads data Angers2003
data_dir <- ("input_data/")
angers2003reflectance <- read_csv(paste0(data_dir, "angers2003reflectance.csv"))
head(angers2003reflectance)
#Transform data into the necessary format to put it into a spec library
matrix <- t(as.matrix(angers2003reflectance[,2:ncol(angers2003reflectance)]))
wavel <- as.numeric(as.matrix(angers2003reflectance %>%  dplyr::select(`Sample_#`)))
wavel
#crecreate Speclibdata and plots data
nSpeclib <- speclib(matrix, wavel)
plot(nSpeclib, main="angers2003")

##
#The same for the data Lopex 1993
#Lopex1993
lopex1993 <- read_csv(paste0(data_dir,"lopex1993reflectance.csv"))
head(lopex1993)
matrix2 <- t(as.matrix(lopex1993[,2:ncol(lopex1993)]))
wavel2 <- as.numeric(as.matrix(lopex1993 %>%  dplyr::select(`Sample_#`)))
wavel2
#recreate S?eclibdata
nSpeclib2 <- speclib(matrix2, wavel2)
plot(nSpeclib2, main="lopex93")


#### Wavelength with maximum standard deviation to choose for calibration
which(apply(matrix,2,FUN=sd) == max(apply(matrix,2,FUN=sd)))
which(apply(matrix2,2,FUN=sd) == max(apply(matrix2,2,FUN=sd)))

