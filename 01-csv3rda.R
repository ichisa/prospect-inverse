
library(readr)
library(dplyr)

data_dir <- ("input_data/")

###################################################
#Loads metadata - leaf characteristics - paramters#
###################################################


angers2003metadata <- read_csv(paste0(data_dir,"angers2003metadata.csv"), 
                               col_types = cols(`English Name` = col_skip(), 
                                                X19 = col_skip(), X20 = col_skip(), 
                                                X21 = col_skip(), X22 = col_skip()))
head(angers2003metadata, 20)
colnames(angers2003metadata)
angers2003metadata <- angers2003metadata %>% rename( N=`Structure parameter`, 
                                                     Cab=`Chlorophyll_a+b (µg/cm²)`, 
                                                     Car=`Carotenoid (µg/cm²)`, 
                                                     Cw=`Equivalent Water Thickness (g/cm²)`,
                                                     Cm=`Leaf mass per area (g/cm²)`,
                                                     sample=`Sample_#`, 
                                                     species=`Latin Name`) %>% 
  dplyr::select(N, Cab, Car, Cw, Cm, sample, species) %>% 
  na.omit()

#Lopex
lopex1993metadata <- read_csv("input_data/lopex1993metadata.csv")

head(lopex1993metadata, 20)
tail(lopex1993metadata, 20)
colnames(lopex1993metadata)
lopex1993metadata <- lopex1993metadata %>% rename( N=`Structure parameter`, 
                                                   Cab=`Chlorophyll_a+b (µg/cm²)`, 
                                                   Car=`Carotenoid (µg/cm²)`, 
                                                   Cw=`Equivalent Water Thickness (g/cm²)`,
                                                   Cm=`Leaf mass per area (g/cm²)`,
                                                   sample=`Sample_#`, 
                                                   species=`Latin Name`) %>% 
  dplyr::select(N, Cab, Car, Cw, Cm, sample, species) %>% 
  na.omit()

lopex1993NA<- unique(which(lopex1993metadata == -999, arr.ind = TRUE)[,1])
lopex1993metadata <- lopex1993metadata[-lopex1993NA,]# remove columns with no measurement



#####################
#Load spectrum data
#####################

#Loads data Angers2003
data_dir <- ("input_data/")
angers2003reflectance <- read_csv(paste0(data_dir, "angers2003reflectance.csv"))

#Wavelength
matrix <- t(as.matrix(angers2003reflectance[,2:ncol(angers2003reflectance)]))
angers2003wavel <- as.numeric(as.matrix(angers2003reflectance %>%  dplyr::select(`Sample_#`)))

#Spectra
angers2003reflectance <- t(angers2003reflectance)
angers2003reflectance[1:10, 1:10]
cn <- as.character(angers2003reflectance[1,])
cn <- paste0("W", cn)#change name of variables 

angers2003reflectance <- angers2003reflectance [2:nrow(angers2003reflectance),]
colnames(angers2003reflectance) <- cn
angers2003reflectance <- as.data.frame(angers2003reflectance)
angers2003reflectance[1:10, 1:10]

rm(matrix, cn)

#Lopex1993

lopex1993 <- read_csv(paste0(data_dir,"lopex1993reflectance.csv"))
#wavel
matrix <- t(as.matrix(lopex1993[,2:ncol(lopex1993)]))
lopex1993wavel <- as.numeric(as.matrix(lopex1993 %>%  dplyr::select(`Sample_#`)))

#Spectra
lopex1993 <- t(lopex1993)
cn <- as.character(lopex1993[1,])
lopex1993 <- lopex1993 [2:nrow(lopex1993),]
colnames(lopex1993) <- paste0("W", cn)#change name of variables 
lopex1993[1:10, 1:10]
plot(lopex1993[1,])

#remove samples that have no info about parameters

rm(matrix, cn)

#discard wavelength hat are not in both datasets
#checks data goes one by one
unique(diff(lopex1993wavel))
unique(diff(angers2003wavel))
lopex1993reflectance <- lopex1993[,1:2051]
lopex1993reflectance <- as.data.frame(lopex1993reflectance)
lopex1993wavel <- lopex1993wavel[1:2051]
lopex1993reflectance <- lopex1993reflectance[-lopex1993NA,]
rm(lopex1993)


save(angers2003metadata, angers2003reflectance, lopex1993metadata, lopex1993reflectance, file=paste0(data_dir, "prospect.RData"))
