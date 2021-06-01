setwd("~/GitHub/Marquez_et_al_2021_New_Gasex_theory")
source('Coupled_FvCF_USO_MSWF_models.R')
load('2_Aci_parameters.Rdata',verbose=TRUE)
library(cowplot)

### Variables and parameters for the simulations
##Trash part colours = c('#133831','white','#3CA4A7') c('grey92','grey60','grey20')
PFD=1000;cs=400;Tleaf=Tair=25+273.15;RH=60

size_p=1
#########################
###   Qin variation   ###
#########################
param=f.make.param(g0=0.03,model.gs = 'USO',TBM = 'FATES',g1=1.89,VcmaxRef = 36.7,JmaxRef=66.9,RdRef=0.33,TpRef=4.09)
Qin_Fick=f.A(PFD=0.01:2000,cs=cs,Tleaf=Tleaf,Tair=Tair,RH=RH,model_diff='Fick',param=param)
simu_Qin=data.frame()
for(gcw in seq(0,25*10^-3,0.1*10^-3)){
  res=f.A(PFD=0.01:2000,cs=cs,Tleaf=Tleaf,Tair=Tair,RH=RH,model_diff='MSWF',
          gcw=gcw,param=param)
  simu_Qin=rbind.data.frame(simu_Qin,cbind.data.frame(as.data.frame(res),gcw=gcw,Qin=0.01:2000))
}
## An figure
a=(ggplot(data=simu_Qin,
        aes(x=Qin,y=A,color=gcw))+geom_point(size=size_p)+ylim(c(-0.5,12))
  +ylab(expression(italic(A)[n]~mu*mol~m^-2~s^-1))
  +xlab(expression(Irradiance~mu*mol~m^-2~s^-1))
  +labs(color = expression(italic(g)[cw]~mol~m^-2~s^-1))
  + scale_color_gradientn(colours = c('#133831','white','#3CA4A7'))
  +theme_bw()+geom_line(data=cbind.data.frame(as.data.frame(Qin_Fick),Qin=0.01:2000),aes(linetype='Fick'),color='red')
  +theme(legend.text.align = 0,panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank()))

## ET figure
e=(ggplot(data=simu_Qin,
          aes(x=Qin,y=ET*1000,color=gcw))+geom_point(size=size_p)+ylim(c(0,3.2))
   +ylab(expression(italic(E)[T]~mmol~m^-2~s^-1))
   +xlab(expression(Irradiance~mu*mol~m^-2~s^-1))
   +labs(color = expression(italic(g)[cw]~mol~m^-2~s^-1))
   + scale_color_gradientn(colours = c('#133831','white','#3CA4A7'))
   +theme_bw()+geom_line(data=cbind.data.frame(as.data.frame(Qin_Fick),Qin=0.01:2000),aes(linetype='Fick'),color='red')
   +theme(legend.text.align = 0,panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))

#####################
### cs variation  ###
#####################


cs_Fick=f.A(PFD=PFD,cs=seq(200,800,1),Tleaf=Tleaf,Tair=Tair,RH=RH,model_diff='Fick',param=param)
simu_cs=data.frame()
for(gcw in seq(0,25*10^-3,0.1*10^-3)){
  res=f.A(PFD=PFD,cs=seq(200,800,1),Tleaf=Tleaf,Tair=Tair,RH=RH,model_diff='MSWF',
          gcw=gcw,param=param)
  simu_cs=rbind.data.frame(simu_cs,cbind.data.frame(as.data.frame(res),gcw=gcw,cs=seq(200,800,1)))
}
b=(ggplot(data=simu_cs,
          aes(x=cs,y=A,color=gcw))+geom_point(size=size_p)+ylim(c(-0.5,12))+xlim(c(200,800))
   +ylab(expression(italic(A)[n]~mu*mol~m^-2~s^-1))
   +xlab(expression(italic(C)[s]~ppm))
   +labs(color = expression(italic(g)[cw]~mol~m^-2~s^-1))
   + scale_color_gradientn(colours = c('#133831','white','#3CA4A7'))
   +theme_bw()
   +geom_line(data=cbind.data.frame(as.data.frame(cs_Fick),cs=seq(200,800,1)),aes(linetype='Fick'),color='red')
   +theme(legend.text.align = 0,panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))

f=(ggplot(data=simu_cs,
          aes(x=cs,y=ET*1000,color=gcw))+geom_point(size=size_p)+ylim(c(0,3.2))+xlim(c(200,800))
   +ylab(expression(italic(E)[T]~mmol~m^-2~s^-1))
   +xlab(expression(italic(C)[s]~ppm))
   +labs(color = expression(italic(g)[cw]~mol~m^-2~s^-1))
   + scale_color_gradientn(colours = c('#133831','white','#3CA4A7'))
   +theme_bw()
   +geom_line(data=cbind.data.frame(as.data.frame(cs_Fick),cs=seq(200,800,1)),aes(linetype='Fick'),color='red')
   +theme(legend.text.align = 0,panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))

########################
### Tleaf variation  ###
########################

Tleaf_Fick=f.A(PFD=PFD,cs=cs,Tleaf=seq(5,45,0.1)+273.15,Tair=seq(5,45,0.1)+273.15,RH=RH,model_diff='Fick',param=param)
simu_Tleaf=data.frame()
for(gcw in seq(0,25*10^-3,0.1*10^-3)){
  res=f.A(PFD=PFD,cs=cs,Tleaf=seq(5,45,0.1)+273.15,Tair=seq(5,45,0.1)+273.15,RH=RH,model_diff='MSWF',
          gcw=gcw,param=param)
  simu_Tleaf=rbind.data.frame(simu_Tleaf,cbind.data.frame(as.data.frame(res),gcw=gcw,Tleaf=seq(5,45,0.1)+273.15))
}
c=(ggplot(data=simu_Tleaf,
          aes(x=Tleaf-273.16,y=A,color=gcw))+geom_point(size=size_p)+ylim(c(-0.5,12))
   +ylab(expression(italic(A)[n]~mu*mol~m^-2~s^-1))
   +xlab(expression(Leaf~Temperature~degree*C))
   +labs(color = expression(italic(g)[cw]~mol~m^-2~s^-1))
   + scale_color_gradientn(colours = c('#133831','white','#3CA4A7'))
   +theme_bw()+geom_line(data=cbind.data.frame(as.data.frame(Tleaf_Fick),Tleaf=seq(5,45,0.1)+273.15),aes(linetype='Fick'),color='red')
   +theme(legend.text.align = 0,panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))

g=(ggplot(data=simu_Tleaf,
          aes(x=Tleaf-273.16,y=ET*1000,color=gcw))+geom_point(size=size_p)+ylim(c(0,3.2))
   +ylab(expression(italic(E)[T]~mmol~m^-2~s^-1))
   +xlab(expression(Leaf~Temperature~degree*C))
   +labs(color = expression(italic(g)[cw]~mol~m^-2~s^-1))
   + scale_color_gradientn(colours = c('#133831','white','#3CA4A7'))
   +theme_bw()+geom_line(data=cbind.data.frame(as.data.frame(Tleaf_Fick),Tleaf=seq(5,45,0.1)+273.15),aes(linetype='Fick'),color='red')
   +theme(legend.text.align = 0,panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))

#####################
### RH variation  ###
#####################

RH_Fick=f.A(PFD=PFD,cs=cs,Tleaf=Tleaf,Tair=Tair,RH=seq(10,95,0.1),model_diff='Fick',param=param)
simu_RH=data.frame()
for(gcw in seq(0,25*10^-3,0.1*10^-3)){
  res=f.A(PFD=PFD,cs=cs,Tleaf=Tleaf,Tair=Tair,RH=seq(10,95,0.1),model_diff='MSWF',
          gcw=gcw,param=param)
  simu_RH=rbind.data.frame(simu_RH,cbind.data.frame(as.data.frame(res),gcw=gcw,RH=seq(10,95,0.1)))
}
d=(ggplot(data=simu_RH,
          aes(x=RH,y=A,color=gcw*1000))+geom_point(size=size_p)+ylim(c(-0.5,12))
   +ylab(expression(italic(A)[n]~mol~m^-2~s^-1))
   +xlab(expression(Relative~humidity~'%'))+xlim(c(0,100))
   +labs(color = expression(italic(g)[cw]~mmol~m^-2~s^-1),linetype='')
   + scale_color_gradientn(colours = c('#133831','white','#3CA4A7'))
   +theme_bw()+geom_line(data=cbind.data.frame(as.data.frame(RH_Fick),RH=seq(10,95,0.1)),color='red')
   +theme(legend.text.align = 0,panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))







h=(ggplot(data=simu_RH,
          aes(x=RH,y=ET*1000,color=gcw))+geom_point(size=size_p)+ylim(c(0,3.2))
   +ylab(expression(italic(E)[T]~mmol~m^-2~s^-1))
   +xlab(expression(Relative~Humidity~'%'))+xlim(c(0,100))
   +labs(color = expression(MSWF:~italic(g)[cw]~mol~m^-2~s^-1),linetype='')
   + scale_color_gradientn(colours = c('#133831','white','#3CA4A7'))
   +theme_bw()+geom_line(data=cbind.data.frame(as.data.frame(RH_Fick),RH=seq(10,95,0.1)),aes(linetype='Fick'),color='red')
   +theme(legend.text.align = 0,panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))


jpeg(filename = '3_diagnostic_plots.jpeg',width = 200,height = 170,units = 'mm',res=300)

fig=plot_grid(a+theme(legend.position = 'none'),
              b+theme(legend.position = 'none'),
              c+theme(legend.position = 'none'),
              d+theme(legend.position = 'none'),
              align = "hv",labels = "auto")
leg=get_legend(d)
plot_grid(fig,leg,ncol=2,rel_widths = c(0.75,0.25))
dev.off()

jpeg(filename = '3_diagnostic_plots_ET.jpeg',width = 200,height = 170,units = 'mm',res=300)

fig=plot_grid(e+theme(legend.position = 'none'),
              f+theme(legend.position = 'none'),
              g+theme(legend.position = 'none'),
              h+theme(legend.position = 'none'),
              align = "hv",labels = "auto")
leg=get_legend(d)
plot_grid(fig,leg,ncol=2,rel_widths = c(0.75,0.25))
dev.off()

jpeg(filename = '3_diagnostic_plots_combined.jpeg',width = 200,height = 140,units = 'mm',res=300)

fig=plot_grid(a+theme(legend.position = 'none'),
              b+theme(legend.position = 'none'),
              c+theme(legend.position = 'none'),
              d+theme(legend.position = 'none'),
              e+theme(legend.position = 'none'),
              f+theme(legend.position = 'none'),
              g+theme(legend.position = 'none'),
              h+theme(legend.position = 'none'),
              align = "hv",labels = "auto",ncol=4)
leg=get_legend(d+theme(legend.position = 'bottom'))
plot_grid(fig,leg,ncol=1,rel_heights = c(0.8,0.2))
dev.off()