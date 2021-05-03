####
## The aim of this code is to test the effect of the new theory on the fitting 
## of ACi curves

setwd("~/GitHub/NewTheoryGasEx")
library(LeafGasExchange)

load('1_Aci_data_QaQc.Rdata')
colnames(Aci_data)[colnames(Aci_data)=='X.U.0394.Pcham']='X.Pcham'
## I first work on the curve BNL17167 which is quite nice

Aci=Aci_data[Aci_data$Barcode=='BNL17167',]
Aci$Tleaf=Aci$Tleaf
Aci=Aci[order(Aci$Ci),]
Recomp_Aci=Aci[,c('A','Tleaf','Ci','Qin','gsw')]
Recomp_Aci$Recomp='Reference'


## I recompute this file using Marquez et al. 2021 theory using different values 
## of cuticular conductance
for(gcw in seq(0,20*10^-3,1*10^-3)){
  Recomp=f.Comput_GasEx(LICOR6800_data = Aci,gcw = gcw)
  Recomp_Aci=rbind.data.frame(Recomp_Aci,data.frame(A=Aci$A,
                                                    Tleaf=Aci$Tleaf,
                                                    Ci=Recomp$Ci,
                                                    Qin=Aci$Qin,
                                                    gsw=Recomp$gsw,
                                                    Recomp=gcw))

}

## Here, I check that the model from Marquez et al. 2021 for a gcw of 0 correponds
## to the outputs given by the LICOR

plot(x=Recomp_Aci[Recomp_Aci$Recomp=='Reference','Ci'],
     y=Recomp_Aci[Recomp_Aci$Recomp=='0','Ci'],
     xlab='Ci LICOR6800',ylab='Ci Marquez et al. 2021, gcw=0')
abline(c(0,1))

## Effect of different gcw on an Aci curve
ggplot(data=Recomp_Aci[Recomp_Aci$Recomp!='Reference',],aes(x=Ci,y=A,color=Recomp))+geom_point()+xlab(expression(italic(C)[i]~calculated~using~Marquez~et~al.~2021))+labs(color = expression(italic(g)[cw]~mol~m^-2~s^-1))


## Here I show the effects on the fitting of Vcmax25
Recomp_Aci$Tleaf=Recomp_Aci$Tleaf+273.16
Aci$Tleaf=Aci$Tleaf+273.16
## Computing Vcmax25, Jmax25 and Tp25 for the different curves
Without_TPU=f.fitting(measures = Aci,id.name = 'Barcode',
                      Start =list(JmaxRef=65,RdRef=1,VcmaxRef=45),
                      param = f.make.param(TpRef=100),type = 'Aci')

With_TPU=f.fitting(measures = Aci,id.name = 'Barcode',
                   Start =list(JmaxRef=Without_TPU[[2]]@coef['JmaxRef'],
                               RdRef=Without_TPU[[2]]@coef['RdRef'],
                               VcmaxRef=Without_TPU[[2]]@coef['VcmaxRef'],
                               TpRef=Without_TPU[[2]]@coef['VcmaxRef']/10),
                   param = f.make.param(),type = 'Aci')

AIC(Without_TPU[[2]])
AIC(With_TPU[[2]])

pdf(file = paste('1_Aci_fitting','.pdf'))
result_various_gcw=by(data = Recomp_Aci,INDICES = list(Recomp_Aci$Recomp),
                            FUN = function(x){
                              f.fitting(measures = x,
                                        Start = list(JmaxRef=With_TPU[[2]]@coef['JmaxRef'],
                                                     RdRef=With_TPU[[2]]@coef['RdRef'],
                                                     VcmaxRef=With_TPU[[2]]@coef['VcmaxRef'],
                                                     TpRef=With_TPU[[2]]@coef['TpRef']),
                                        param=f.make.param(),
                                        id.name = 'Recomp',modify.init = FALSE)})
dev.off()


res=as.data.frame(t(sapply(result_various_gcw,
                                             FUN = function(x) c(x[[2]]@coef,
                                                                 sqrt(diag(vcov(x[[2]]))),
                                                                 AIC=AIC(x[[2]])))))
res$Recomp=row.names(res)
plot(y=res$VcmaxRef,x=as.numeric(res$Recomp)*1000,ylab=expression(italic(V)[cmax25]),xlab=expression(italic(g)[cw]))
plot(y=res$JmaxRef,x=as.numeric(res$Recomp)*1000,ylab=expression(italic(J)[max25]),xlab=expression(italic(g)[cw]))
plot(y=res$TpRef,x=as.numeric(res$Recomp)*1000,ylab=expression(italic(TPU)[25]),xlab=expression(italic(g)[cw]))
plot(y=res$RdRef,x=as.numeric(res$Recomp)*1000,ylab=expression(italic(R)[day25]),xlab=expression(italic(g)[cw]))
