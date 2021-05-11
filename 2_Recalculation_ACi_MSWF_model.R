## The aim of this code is to test the effect of the new theory on the ACi measurement
## and fitting

setwd("~/GitHub/Marquez_et_al_2021_New_Gasex_theory")
library(LeafGasExchange)
library(cowplot)
source('0_NewTheoryGasex.R')
load('1_Aci_data_QaQc.Rdata')


#############################################################
#### Recalculation of the ACi values from the Panama data ###
#############################################################

colnames(Aci_data_Panama)[colnames(Aci_data_Panama)=='X.U.0394.Pcham']='X.Pcham'

Aci=Aci_data_Panama[order(Aci_data_Panama$Ci),]
Recomp_Aci_Panama=Aci[,c('A','Tleaf','Ci','Qin','gsw')]
Recomp_Aci_Panama$Recomp='Reference'

## I recompute this file using Marquez et al. 2021 theory using different values 
## of cuticular conductance
for(gcw in seq(0,25*10^-3,1*10^-3)){
  Recomp=f.Comput_GasEx_6800(LICOR6800_data = Aci,gcw = gcw)
  Recomp_Aci_Panama=rbind.data.frame(Recomp_Aci_Panama,data.frame(A=Aci$A,
                                                    Tleaf=Aci$Tleaf,
                                                    Ci=Recomp$Ci,
                                                    Qin=Aci$Qin,
                                                    gsw=Recomp$gsw,
                                                    Recomp=gcw))

}

Recomp_Aci_Panama$Species='Guatteria dumetorum'
Recomp_Aci_Panama$mean_gbw=mean(2*Aci$gbw)

#############################################################
#### Recalculation of the ACi values from the Arctic data ###
#############################################################
Aci=Aci_data_arctic[order(Aci_data_arctic$Ci),]
Recomp_Aci_arctic=Aci[,c('Photo','Tleaf','Ci','PARi','Cond')]
Recomp_Aci_arctic$Recomp='Reference'

## I recompute this file using Marquez et al. 2021 theory using different values 
## of cuticular conductance
for(gcw in seq(0,25*10^-3,1*10^-3)){
  Recomp=f.Comput_GasEx_6400(LICOR6400_data = Aci,gcw = gcw)
  Recomp_Aci_arctic=rbind.data.frame(Recomp_Aci_arctic,data.frame(Photo=Aci$Photo,
                                                    Tleaf=Aci$Tleaf,
                                                    Ci=Recomp$Ci,
                                                    PARi=Aci$PARi,
                                                    Cond=Recomp$gsw,
                                                    Recomp=gcw))
  
}
Recomp_Aci_arctic$Species='Petasites frigidus'
Recomp_Aci_arctic$mean_gbw=mean(2*Aci_data_arctic$BLCond)
colnames(Recomp_Aci_arctic)=colnames(Recomp_Aci_Panama)

## Effect of different gcw on an Aci curve
ggplot(data=Recomp_Aci_arctic[Recomp_Aci_arctic$Recomp!='Reference',],aes(x=Ci,y=A,color=Recomp))+geom_point()+xlab(expression(italic(C)[i]~calculated~using~Marquez~et~al.~2021))+labs(color = expression(italic(g)[cw]~mol~m^-2~s^-1))

##########################################################
#### Recalculation of the ACi values from the Oak data ###
##########################################################

Aci=Aci_data_oak[order(Aci_data_oak$Ci),]
Recomp_Aci_oak=Aci[,c('Photo','Tleaf','Ci','PARi','Cond')]
Recomp_Aci_oak$Recomp='Reference'

## I recompute this file using Marquez et al. 2021 theory using different values 
## of cuticular conductance
for(gcw in seq(0,25*10^-3,1*10^-3)){
  Recomp=f.Comput_GasEx_6400(LICOR6400_data = Aci,gcw = gcw)
  Recomp_Aci_oak=rbind.data.frame(Recomp_Aci_oak,data.frame(Photo=Aci$Photo,
                                                                  Tleaf=Aci$Tleaf,
                                                                  Ci=Recomp$Ci,
                                                                  PARi=Aci$PARi,
                                                                  Cond=Recomp$gsw,
                                                                  Recomp=gcw))
  
}
Recomp_Aci_oak$Species='Quercus coccinea Münchh'
Recomp_Aci_oak$mean_gbw=mean(2*Aci_data_oak$BLCond)
colnames(Recomp_Aci_oak)=colnames(Recomp_Aci_Panama)

Recomp_Aci=rbind.data.frame(Recomp_Aci_Panama,Recomp_Aci_arctic,Recomp_Aci_oak)

#######################################
## Figure effect of cc on Aci curves ##
#######################################

# First, I check if the A-Ci curves measured by the LICOR correspond to the recomputed values for gcw = 0

plot(x=Recomp_Aci[Recomp_Aci$Recomp=='Reference','Ci'],y=Recomp_Aci[Recomp_Aci$Recomp=='0','Ci'])
abline(c(0,1))
       
# Figure of the effect of cc
data_fig=Recomp_Aci[Recomp_Aci$Recomp!='Reference',]
data_fig$Recomp=as.numeric(data_fig$Recomp)  
data_fig$Species=factor(as.character(data_fig$Species),levels = c("Quercus coccinea Münchh","Petasites frigidus","Guatteria dumetorum"),ordered = TRUE)
library(viridis)
jpeg(filename = 'Figure1.jpeg',width = 150,height = 92,units = 'mm',res=300)
(ggplot(data=data_fig,
       aes(x=Ci,y=A,color=Recomp,shape=Species))+geom_point()
      +ylab(expression(italic(A)[n]~mu*mol~m^-2~s^-1))
      +xlab(expression(italic(C)[i]~ppm))
      +labs(color = expression(italic(g)[cw]~mol~m^-2~s^-1))
  + scale_color_gradientn(colours = c('grey92','grey60','grey20'))
      +scale_shape_discrete(labels=
          c(expression(italic(Quercus~cocinear)~Münchh),
            expression(italic(Petasides~frigidus)),
            expression(italic(Guatteria~dumetorum))
            ))
      +theme_bw()
  +theme(legend.text.align = 0,panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
dev.off()

###########################################################################
## Effect of the differnt Aci curves on the estimation of the parameters ##
## of the A-Ci curves                                                    ##
###########################################################################

## Correcting the units of the temperature
Recomp_Aci_Panama$Tleaf=Recomp_Aci_Panama$Tleaf+273.16
Recomp_Aci_arctic$Tleaf=Recomp_Aci_arctic$Tleaf+273.16
Recomp_Aci_oak$Tleaf=Recomp_Aci_oak$Tleaf+273.16


##### Panama

Without_TPU=f.fitting(measures = Recomp_Aci_Panama[Recomp_Aci_Panama$Recomp=='Reference',],id.name = 'Species',
                      Start =list(JmaxRef=65,RdRef=1,VcmaxRef=45),
                      param = f.make.param(TpRef=100),type = 'Aci')

With_TPU=f.fitting(measures = Recomp_Aci_Panama[Recomp_Aci_Panama$Recomp=='Reference',],id.name = 'Species',
                   Start =list(JmaxRef=Without_TPU[[2]]@coef['JmaxRef'],
                               RdRef=Without_TPU[[2]]@coef['RdRef'],
                               VcmaxRef=Without_TPU[[2]]@coef['VcmaxRef'],
                               TpRef=Without_TPU[[2]]@coef['VcmaxRef']/10),
                   param = f.make.param(),type = 'Aci')

AIC(Without_TPU[[2]])
AIC(With_TPU[[2]])

pdf(file = '1_Aci_fitting_Panama.pdf')
result_various_gcw_Panama=by(data = Recomp_Aci_Panama,INDICES = list(Recomp_Aci_Panama$Recomp),
                            FUN = function(x){
                              f.fitting(measures = x,
                                        Start = list(JmaxRef=With_TPU[[2]]@coef['JmaxRef'],
                                                     RdRef=With_TPU[[2]]@coef['RdRef'],
                                                     VcmaxRef=With_TPU[[2]]@coef['VcmaxRef'],
                                                     TpRef=With_TPU[[2]]@coef['TpRef']),
                                        param=f.make.param(),
                                        id.name = 'Recomp',modify.init = FALSE)})
dev.off()


res_Panama=as.data.frame(t(sapply(result_various_gcw_Panama,
                                             FUN = function(x) c(x[[2]]@coef,
                                                                 AIC=AIC(x[[2]])))))
res_Panama$Recomp=row.names(res_Panama)


##### arctic

Without_TPU=f.fitting(measures = Recomp_Aci_arctic[Recomp_Aci_arctic$Recomp=='Reference',],id.name = 'Species',
                      Start =list(JmaxRef=65,RdRef=1,VcmaxRef=45),
                      param = f.make.param(TpRef=100),type = 'Aci')

With_TPU=f.fitting(measures = Recomp_Aci_arctic[Recomp_Aci_arctic$Recomp=='Reference',],id.name = 'Species',
                   Start =list(JmaxRef=Without_TPU[[2]]@coef['JmaxRef'],
                               RdRef=Without_TPU[[2]]@coef['RdRef'],
                               VcmaxRef=Without_TPU[[2]]@coef['VcmaxRef'],
                               TpRef=Without_TPU[[2]]@coef['VcmaxRef']/10),
                   param = f.make.param(),type = 'Aci')

AIC(Without_TPU[[2]])
AIC(With_TPU[[2]])

pdf(file = '1_Aci_fitting_arctic.pdf')
result_various_gcw_arctic=by(data = Recomp_Aci_arctic,INDICES = list(Recomp_Aci_arctic$Recomp),
                             FUN = function(x){
                               f.fitting(measures = x,
                                         Start = list(JmaxRef=Without_TPU[[2]]@coef['JmaxRef'],
                                                      RdRef=Without_TPU[[2]]@coef['RdRef'],
                                                      VcmaxRef=Without_TPU[[2]]@coef['VcmaxRef']),
                                         param=f.make.param(TpRef = 100),
                                         id.name = 'Recomp',modify.init = FALSE)})
dev.off()


res_arctic=as.data.frame(t(sapply(result_various_gcw_arctic,
                                  FUN = function(x) c(x[[2]]@coef,
                                                      AIC=AIC(x[[2]])))))
res_arctic$Recomp=row.names(res_arctic)


######## oak 


##### oak

Without_TPU=f.fitting(measures = Recomp_Aci_oak[Recomp_Aci_oak$Recomp=='Reference',],id.name = 'Species',
                      Start =list(JmaxRef=65,RdRef=1,VcmaxRef=45),
                      param = f.make.param(TpRef=100),type = 'Aci')

With_TPU=f.fitting(measures = Recomp_Aci_oak[Recomp_Aci_oak$Recomp=='Reference',],id.name = 'Species',
                   Start =list(JmaxRef=Without_TPU[[2]]@coef['JmaxRef'],
                               RdRef=Without_TPU[[2]]@coef['RdRef'],
                               VcmaxRef=Without_TPU[[2]]@coef['VcmaxRef'],
                               TpRef=Without_TPU[[2]]@coef['VcmaxRef']/10),
                   param = f.make.param(),type = 'Aci')

AIC(Without_TPU[[2]])
AIC(With_TPU[[2]])

pdf(file = '1_Aci_fitting_oak.pdf')
result_various_gcw_oak=by(data = Recomp_Aci_oak,INDICES = list(Recomp_Aci_oak$Recomp),
                             FUN = function(x){
                               f.fitting(measures = x,
                                         Start = list(JmaxRef=Without_TPU[[2]]@coef['JmaxRef'],
                                                      RdRef=Without_TPU[[2]]@coef['RdRef'],
                                                      VcmaxRef=Without_TPU[[2]]@coef['VcmaxRef']),
                                         param=f.make.param(TpRef = 100),
                                         id.name = 'Recomp',modify.init = FALSE)})
dev.off()


res_oak=as.data.frame(t(sapply(result_various_gcw_oak,
                                  FUN = function(x) c(x[[2]]@coef,
                                                      AIC=AIC(x[[2]])))))
res_oak$Recomp=row.names(res_oak)



#########################
### Figure parameters ###
#########################
res_oak$TpRef=NA
res_arctic$TpRef=NA
res_oak$Species="Quercus coccinea Münchh"
res_Panama$Species="Guatteria dumetorum"
res_arctic$Species="Petasites frigidus"
res_all=rbind.data.frame(res_oak[,colnames(res_Panama)],res_arctic[,colnames(res_Panama)],res_Panama)
res_all$Species=factor(as.character(res_all$Species),levels = c("Quercus coccinea Münchh","Petasites frigidus","Guatteria dumetorum"),ordered = TRUE)

data_fig=res_all[res_all$Recomp!='Reference',]
for(species in unique(data_fig$Species)){
  data_fig[data_fig$Species==species,'Vcmax25']=data_fig[data_fig$Species==species,'VcmaxRef']/data_fig[data_fig$Species==species&data_fig$Recomp==0,'VcmaxRef']
  data_fig[data_fig$Species==species,'Jmax25']=data_fig[data_fig$Species==species,'JmaxRef']/data_fig[data_fig$Species==species&data_fig$Recomp==0,'JmaxRef']
  data_fig[data_fig$Species==species,'Tp25']=data_fig[data_fig$Species==species,'TpRef']/data_fig[data_fig$Species==species&data_fig$Recomp==0,'TpRef']
  data_fig[data_fig$Species==species,'Rd25']=data_fig[data_fig$Species==species,'RdRef']/data_fig[data_fig$Species==species&data_fig$Recomp==0,'RdRef']
 }
data_fig$Recomp=as.numeric(data_fig$Recomp) 
size_p=1
a=(ggplot(data=data_fig,
          aes(x=Recomp,y=Vcmax25,shape=Species))+geom_point(size=size_p)
   +ylab(expression(Standardized~italic(V)[cmax]))
   +xlab(expression(italic(g)[cw]~mol~m^-2~s^-1))
   +labs(color = expression(italic(g)[cw]~mol~m^-2~s^-1))
   +ylim(c(0.5,2))
   #+ scale_color_gradientn(colours = c('grey92','grey60','grey20'))
   +scale_shape_discrete(labels=
                           c(expression(italic(Quercus~cocinear)~Münchh),
                             expression(italic(Petasides~frigidus)),
                             expression(italic(Guatteria~dumetorum))
                           ))
   +theme_bw()
   +theme(legend.text.align = 0,panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


b=(ggplot(data=data_fig,
          aes(x=Recomp,y=Jmax25,shape=Species))+geom_point(size=size_p)
   +ylab(expression(Standardized~italic(J)[max]))
   +xlab(expression(italic(g)[cw]~mol~m^-2~s^-1))
   +labs(color = expression(italic(g)[cw]~mol~m^-2~s^-1))
   +ylim(c(0.5,2))
   #+ scale_color_gradientn(colours = c('grey92','grey60','grey20'))
   +scale_shape_discrete(labels=
                           c(expression(italic(Quercus~cocinear)~Münchh),
                             expression(italic(Petasides~frigidus)),
                             expression(italic(Guatteria~dumetorum))
                           ))
   +theme_bw()
   +theme(legend.text.align = 0,panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position='none'))



c=(ggplot(data=data_fig,
          aes(x=Recomp,y=Tp25,shape=Species))+geom_point(size=size_p)
   +ylab(expression(Standardized~italic(TPU)))
   +xlab(expression(italic(g)[cw]~mol~m^-2~s^-1))
   +labs(color = expression(italic(g)[cw]~mol~m^-2~s^-1))
   +ylim(c(0.5,2))
   #+ scale_color_gradientn(colours = c('grey92','grey60','grey20'))
   +scale_shape_discrete(labels=
                           c(expression(italic(Quercus~cocinear)~Münchh),
                             expression(italic(Petasides~frigidus)),
                             expression(italic(Guatteria~dumetorum))
                           ))
   +theme_bw()
   +theme(legend.text.align = 0,panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position='none'))




d=(ggplot(data=data_fig,
          aes(x=Recomp,y=Rd25,shape=Species))+geom_point(size=size_p)
   +ylab(expression(Standardized~italic(R)[day]))
   +xlab(expression(italic(g)[cw]~mol~m^-2~s^-1))
   +labs(color = expression(italic(g)[cw]~mol~m^-2~s^-1))
   +ylim(c(0.5,2))
   #+ scale_color_gradientn(colours = c('grey92','grey60','grey20'))
   +scale_shape_discrete(labels=
                           c(expression(italic(Quercus~cocinear)~Münchh),
                             expression(italic(Petasides~frigidus)),
                             expression(italic(Guatteria~dumetorum))
                           ))
   +theme_bw()
   +theme(legend.text.align = 0,panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position='none'))

jpeg(filename = 'Figure2.jpeg',width = 210,height = 145,units = 'mm',res=300)
half=plot_grid(a+theme(legend.position = 'none'),b,c,d,align='hv',ncol=2,labels = 'auto')
leg=get_legend(a)
plot_grid(half,leg,ncol = 2,rel_widths = c(0.75,0.25))
dev.off()

save(res_all,file='2_Aci_parameters.Rdata')