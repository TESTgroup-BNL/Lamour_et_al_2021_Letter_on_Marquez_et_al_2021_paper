###############################################################
### The aim of this code is to import and check the quality ###
### of the ACi curves                                       ###
###############################################################

library(LeafGasExchange) #https://github.com/TESTgroup-BNL/LeafGasExchange
library(here)
setwd(here())
source('0_LICOR_Recalculations_functions_M2021.R')



Aci_data_Panama=f.import_licor6800(file = '~/GitHub/Marquez_et_al_2021_New_Gasex_theory/Datasets/20200122_ACi_Lindsey_BNL17167.xlsx')
colnames(Aci_data_Panama)[colnames(Aci_data_Panama)=='X.U.0394.Pcham']='X.Pcham'

plot(x=Aci_data_Panama$Ci,Aci_data_Panama$A,cex=2)
text(x=Aci_data_Panama$Ci,Aci_data_Panama$A, labels=Aci_data_Panama$obs,cex=0.7)
remove=c(19,9:11)
Aci_data_Panama=Aci_data_Panama[!Aci_data_Panama$obs%in%remove,]
points(x=Aci_data_Panama$Ci,y=Aci_data_Panama$A,cex=2,col='red')
'Guatteria dumetorum'

#####################
#### Arctic data  ###
#####################
Aci_data_arctic=read.csv('Datasets/NGEE-Arctic_A-Ci_curves_2012-2015.csv',skip = 8,na.strings = '-9999',header = TRUE)
Aci_data_arctic=Aci_data_arctic[Aci_data_arctic$Sample_Barcode=='1557',]
Aci_data_arctic$obs=1:18
plot(x=Aci_data_arctic$Ci,Aci_data_arctic$Photo,cex=2)
text(x=Aci_data_arctic$Ci,Aci_data_arctic$Photo, labels=Aci_data_arctic$obs,cex=0.7)
remove=c(7,8:10,18)
Aci_data_arctic=Aci_data_arctic[!Aci_data_arctic$obs%in%remove,]
points(x=Aci_data_arctic$Ci,y=Aci_data_arctic$Photo,cex=2,col='red')
'Petasites frigidus'

################
### Oak data ###
################
Aci_data_oak=read.csv('Datasets/raw_gas_exchange_li6400_oak.csv',header = TRUE)
for (i in 1:length(unique(Aci_data_oak$Sample_ID))){
plot(x=Aci_data_oak[Aci_data_oak$Sample_ID==unique(Aci_data_oak$Sample_ID)[i],'Ci'],
     y=Aci_data_oak[Aci_data_oak$Sample_ID==unique(Aci_data_oak$Sample_ID)[i],'Photo'],cex=2, main= i)
text(x=Aci_data_oak[Aci_data_oak$Sample_ID==unique(Aci_data_oak$Sample_ID)[i],'Ci'],
     y=Aci_data_oak[Aci_data_oak$Sample_ID==unique(Aci_data_oak$Sample_ID)[i],'Photo'],
     labels = Aci_data_oak[Aci_data_oak$Sample_ID==unique(Aci_data_oak$Sample_ID)[i],'Obs'],cex=0.7)  
}
#unique(Aci_data_oak$Sample_ID)[14]
Aci_data_oak=Aci_data_oak[Aci_data_oak$Sample_ID=='15072',]



save(Aci_data_oak,Aci_data_arctic,Aci_data_Panama,file='1_Aci_data_QaQc.Rdata')
