####################################
### Function to compute the LICOR6800 gasex variables according to the new theory by Marquez et al. 2021
#' Title
#'
#' @param LICOR6800_data Data.frame object corresponding to the excel file raw output from a LICOR6800. Note that the column 'Delta'.Pcham in the excel file is important but is often not correctly imported in R due to the Delta grec symbol. In this data.frame, the column name has to names X.Pcham for the function to work.
#' @param gcw Cuticular conductance in mol m-2 s-1
#' @param Beta Ratio of cuticular conductance to CO2 and cuticular conductance to H20, see 
#' @references Márquez DA, Stuart-Williams H, Farquhar GD. 2021. An improved theory for calculating leaf gas exchange more precisely accounting for small fluxes. Nature Plants.
#' @return Return a list with: 
#' Cs the leaf surface CO2 concentration (ppm)
#' Ci the leaf internal CO2 concentration (ppm)
#' gsw the leaf stomatal conductance to water (mol m-2 s-1)
#' glw the leaf conductance to water (mol m-2 s-1)
#' ws the leaf surface water concentration
#' wi the leaf internal water concentration
#' Es the transpiration through the stomata (mol m-2 s-1)
#' Ec the transpiration through the cuticule (mol m-2 s-1)
#' @export
#'
#' @examples
f.Comput_GasEx_6800<-function(LICOR6800_data,gcw,Beta=0.05){
  if(is.null(LICOR6800_data$X.Pcham)){print('Please, be sure to import the data from the excel output of the LICOR6800 and that the column name Delta.Pcham which uses the grec letter Delta is called X.Pcham in your R data.frame')}
  ## Correspondance between lICOR6800 column names and unit with Marquez et al. 2021 variables
  ET=LICOR6800_data$E
  AT=LICOR6800_data$A
  wa=LICOR6800_data$H2O_s/1000
  gbw=2*LICOR6800_data$gbw # Licor gives the one sided boundary layer, not the total boundary layer, thus the multiplication by 2.
  Ca = LICOR6800_data$CO2_s
  Pa=LICOR6800_data$Pa
  DeltaPa=LICOR6800_data$X.Pcham
  Tleaf=LICOR6800_data$Tleaf

  ## Eq. 9
  ws=(wa+(1-wa/2)*ET/gbw)/(1+ET/(2*gbw))
  
  ## Eq 30 and 40 in LICOR 6800 manual
  ## https://www.licor.com/env/support/LI-6800/topics/additional-calculations.html?Highlight=stomatal%20ratio
  wi=0.61365*exp(17.502*Tleaf/(240.97+Tleaf))/(Pa+DeltaPa)
  
  ## Eq.15
  gtw=1/((wi-ws)/(ET-(ET-gcw*(wi-ws))*((wi+ws)/2))+1/gbw)
  
  ## Eq 66 in LICOR 6800 manual
  ## https://www.licor.com/env/support/LI-6800/topics/additional-calculations.html?Highlight=stomatal%20ratio
  gbc=gbw/1.37
  
  ## Eq. 12
  Cs=(gbc*Ca-AT-ET/2*Ca)/(gbc+ET/2)
  
  ## Eq. 14
  D=wi-ws
  wbars=(wi+ws)/2
  Gi=ET/(1.6*(wi-ws))*(1-wbars)+ET/2
  alpha=(1-wbars)/1.6+D/2-Beta
  
  Ci=Cs-(AT+Cs*ET-Cs*D*gcw)/(Gi-alpha*gcw)
  
  Ec = (gcw)*(wi-ws) ## Eq 10
  Es = ET - Ec
  
  gsw = (Es-Es*((wi+ws)/2))/(wi-ws)
  
  return(list(Cs=Cs,Ci=Ci,gsw=gsw,glw=gsw+gcw,ws=ws,wi=wi,Es=Es,Ec=Ec)) 
}

## The function is compared to the result Marquez et al., (2021) obtained on their Figure 4 
data_LICOR6800=read.csv('Figure4.csv') ## data from figure 4, with only the first experiment 

test=f.Comput_GasEx_6800(LICOR6800_data = data_LICOR6800,gcw = 21*10^-3,Beta = 0.05) ## The gcw value was communicated by Farquhar et Marquez through email
test0=f.Comput_GasEx_6800(LICOR6800_data = data_LICOR6800,gcw = 0,Beta = 0.05) ## For comparison with the LICOR6800 data. This simulation should give a very similar result as with the output of the LICOR6800. The result should be identical for the parameter K in the LICOR6800 = 0.5 

par(mfrow=c(1,2),mar=c(4,4.5,2,1))
plot(y=data_LICOR6800$A,x=data_LICOR6800$Ci,xlim=c(0,300),ylim=c(0,22),col='red',ylab=expression(A~(micro*mol~m^-2~s^-1)),xlab=expression(C[i]~(ppm)))
points(y=data_LICOR6800$A,x=test$Ci,col='blue')
points(y=data_LICOR6800$A,x=test0$Ci,col='green')


plot(y=data_LICOR6800$Ci,x=data_LICOR6800$Air.VP.def,xlim=c(1,3),ylim=c(0,300),col='red',ylab=expression(C[i]~(ppm)),xlab = expression(ASD~(kPa)))
points(y=test$Ci,x=data_LICOR6800$Air.VP.def,col='blue')
points(y=test0$Ci,x=data_LICOR6800$Air.VP.def,col='green')

####################################
### Function to compute the LICOR6400 gasex variables according to the new theory by Marquez et al. 2021
#' Title
#'
#' @param LICOR6400_data Data.frame object corresponding to the csv file
#' @param gcw Cuticular conductance in mol m-2 s-1
#' @param Beta Ratio of cuticular conductance to CO2 and cuticular conductance to H20, see 
#' @references Márquez DA, Stuart-Williams H, Farquhar GD. 2021. An improved theory for calculating leaf gas exchange more precisely accounting for small fluxes. Nature Plants.
#' @return Return a list with: 
#' Cs the leaf surface CO2 concentration (ppm)
#' Ci the leaf internal CO2 concentration (ppm)
#' gsw the leaf stomatal conductance to water (mol m-2 s-1)
#' glw the leaf conductance to water (mol m-2 s-1)
#' ws the leaf surface water concentration
#' wi the leaf internal water concentration
#' @export
#'
#' @examples
f.Comput_GasEx_6400<-function(LICOR6400_data,gcw,Beta=0.05){
  ## Correspondance between lICOR6400 column names and unit with Marquez et al. 2021 variables
  ET=LICOR6400_data$Trmmol/1000
  AT=LICOR6400_data$Photo
  wa=LICOR6400_data$H2OS/1000
  gbw=2*LICOR6400_data$BLCond # Licor gives the one sided boundary layer, not the total boundary layer, thus the multiplication by 2.
  Ca = LICOR6400_data$CO2S
  Pa=LICOR6400_data$Press
  DeltaPa=0
  Tleaf=LICOR6400_data$Tleaf
  
  ## Eq. 9
  ws=(wa+(1-wa/2)*ET/gbw)/(1+ET/(2*gbw))
  
  ## Eq 30 and 40 in LICOR 6400 manual
  ## https://www.licor.com/env/support/LI-6400/topics/additional-calculations.html?Highlight=stomatal%20ratio
  wi=0.61365*exp(17.502*Tleaf/(240.97+Tleaf))/(Pa+DeltaPa)
  
  ## Eq.15
  gtw=1/((wi-ws)/(ET-(ET-gcw*(wi-ws))*((wi+ws)/2))+1/gbw)
  
  ## Eq 66 in LICOR 6400 manual
  ## https://www.licor.com/env/support/LI-6400/topics/additional-calculations.html?Highlight=stomatal%20ratio
  gbc=gbw/1.37
  
  ## Eq. 12
  Cs=(gbc*Ca-AT-ET/2*Ca)/(gbc+ET/2)
  
  ## Eq. 14
  D=wi-ws
  wbars=(wi+ws)/2
  Gi=ET/(1.6*(wi-ws))*(1-wbars)+ET/2
  alpha=(1-wbars)/1.6+D/2-Beta
  
  Ci=Cs-(AT+Cs*ET-Cs*D*gcw)/(Gi-alpha*gcw)
  
  Ec = (gcw)*(wi-ws)
  Es = ET - Ec
  
  gsw = (Es-Es*((wi+ws)/2))/(wi-ws)
  
  return(list(Cs=Cs,Ci=Ci,gsw=gsw,glw=gsw+gcw,ws=ws,wi=wi,Ec=Ec,Es=Es)) 
}

