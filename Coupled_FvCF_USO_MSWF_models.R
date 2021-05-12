library(LeafGasExchange)
f.A<-function(PFD,cs,Tleaf,Tair,RH,gcw=10*10^-3,param=f.make.param(),model_diff='Fick'){
  #Calculation of temperature dependence of the parameters
  Kc=f.arrhenius(param[['KcRef']],param[['KcHa']],Tleaf)
  Ko=f.arrhenius(param[['KoRef']],param[['KoHa']],Tleaf)
  Gstar=f.arrhenius(param[['GstarRef']],param[['GstarHa']],Tleaf)
  Rd=f.modified.arrhenius(PRef=param[['RdRef']],param[['RdHa']],param[['RdHd']],param[['RdS']],Tleaf)
  Vcmax=f.modified.arrhenius(PRef=param[['VcmaxRef']],param[['VcmaxHa']],param[['VcmaxHd']],param[['VcmaxS']],Tleaf)
  Jmax=f.modified.arrhenius(PRef=param[['JmaxRef']],param[['JmaxHa']],param[['JmaxHd']],param[['JmaxS']],Tleaf)
  
  
  I2=PFD*param[['abso']]*(param[['aQY']])
  J=(I2+Jmax-((I2+Jmax)^2-4*(param[['Theta']])*I2*Jmax)^0.5)/(2*(param[['Theta']]))
  
  ds=f.ds(Tleaf,Tair,RH)
  cc=NA
  
  wi=0.61365*exp(17.502*(Tleaf-273.16)/(240.97+(Tleaf-273.16)))/(param[['Patm']])
  ws=0.61365*exp(17.502*(Tair-273.16)/(240.97+(Tair-273.16)))/(param[['Patm']])*RH/100
  
  
  #if(model_diff=='Fick'){k=0;l=0}
 # if(model_diff=='vCF'){k=k=ds/(param[['Patm']]*1000)/2*1.6;l=0} ## *1.6 due to a mistake in f.solv
  #if(model_diff=='MSWF'){k=ds/(param[['Patm']]*1000)/2*1.6;l=gcw*0.05;param[['g0']]=param[['g0']]-gcw}
  if(model_diff=='Fick'){k=0;l=0}
  if(model_diff=='vCF'){k=(wi-ws)/(2-(wi+ws))*1.6;l=0} ## *1.6 due to a mistake in f.solv
  if(model_diff=='MSWF'){k=(wi-ws)/(2-(wi+ws))*1.6;l=gcw*0.05;param[['g0']]=param[['g0']]-gcw}
  
   (wi-ws)/(2-(wi+ws))
  #Resolution for CLM4.5 and FATES
  if(param[['TBM']]%in%c(0,2)){
    
    # Analytical solution of the system of equations {E1 : A=f(ci), E2 : gs=f(A,cs) and ci=f(cs)}
    cic=f.solv(x=Vcmax,y=Kc*(1+param[['O2']]/Ko),cs=cs,Rd=Rd,Gstar=Gstar,g0=param[['g0']],g1=param[['g1']],power=param[['power']],ds=ds,RH=RH,k=k,l=l)
    Wc=(cic-Gstar)*Vcmax/(cic+Kc*(1+param[['O2']]/Ko))
    
    cij=f.solv(x=J/4,y=2*Gstar,cs=cs,Rd=Rd,Gstar=Gstar,g0=param[['g0']],g1=param[['g1']],param[['power']],ds=ds,RH=RH,k=k,l=l)
    Wj=(cij-Gstar)*J/(4*cij+8*Gstar)
    
    Tp=f.modified.arrhenius(PRef=param[['TpRef']],param[['TpHa']],param[['TpHd']],param[['TpS']],Tleaf)
    Wp=3*Tp
    cip=f.solv(x=3*Tp,y=-Gstar,cs=cs,Rd=Rd,Gstar=Gstar,g0=param[['g0']],g1=param[['g1']],power=param[['power']],ds=ds,RH=RH,k=k,l=l)
    
    ci=cij
    if(!is.null(which(Wc<Wj))&length(cic)==length(cij)){ci[which(Wc<Wj)]=cic[which(Wc<Wj)]}
    if(!is.null(which(Wc<Wj))&length(cic)!=length(cij)){ci[which(Wc<Wj)]=cic}
    W=pmin(Wc,Wj)
    if(!is.null(which(Wp<W))){ci[which(Wp<W)]=cip[which(Wp<W)]}
    Ai=f.smooth(A1 = Wc,A2 = Wj,theta=param[['thetacj']])
    A=f.smooth(A1=Ai,A2=Wp,theta=param[['thetaip']])-Rd
    gs=f.gs(A=A,cs=cs,ds=ds,Rd=Rd,Gstar=Gstar,RH=RH,g0=param[['g0']],g1=param[['g1']],power=param[['power']],model =param[['model.gs']])
    output=list(A=A,Ac=Wc-Rd,Aj=Wj-Rd,Ap=Wp-Rd,Ag=A+Rd,gs=gs,ci=ci,ds=ds,Transp=gs*ds/(param[['Patm']]*1000)*18)
    return(output)
  }
}


## Apply only for USO model
f.solv<-function(x,y,cs,Rd,Gstar,g0,g1,power,ds,RH,k,l){
  a=(8*Rd*g1*k +8*sqrt(ds)*Rd*k -8*g1*x -8*sqrt(ds)*x -8*g1*k*x
     -5*g0*sqrt(ds)*cs -8*sqrt(ds)*cs*l-8*sqrt(ds)*k*x -5*sqrt(ds)*cs*g0*k +8*sqrt(ds)*Rd +8*Rd*g1)
  
  b=(-5*sqrt(ds)*cs^2*g0*k -5*sqrt(ds)*cs*g0*k*y +8*Gstar*sqrt(ds)*k*x +8*Gstar*g1*k*x
     +8*sqrt(ds)*Rd*cs*k +8*sqrt(ds)*Rd*k*y +5*sqrt(ds)*cs^2*g0 +8*sqrt(ds)*cs^2*l -5*sqrt(ds)*cs*g0*y
     -8*sqrt(ds)*cs*k*x -8*sqrt(ds)*cs*l*y + 8*Rd*cs*g1*k +8*Rd*g1*k*y -8*cs*g1*k*x
     +8*Gstar*sqrt(ds)*x +8*Gstar*g1*x +8*sqrt(ds)*Rd*y -8*Rd*cs*g1 +8*Rd*g1*y +8*cs*g1*x)
  
  c=(-5*sqrt(ds)*cs^2*g0*k*y +8*Gstar*sqrt(ds)*cs*k*x +8*Gstar*cs*g1*k*x +8*sqrt(ds)*Rd*cs*k*y
     +5*sqrt(ds)*cs^2*g0*y +8*sqrt(ds)*cs^2*l*y +8*Rd*cs*g1*k*y -8*Gstar*cs*g1*x -8*Rd*cs*g1*y)
  
  ci2=(-b-(b^2-4*a*c)^0.5)/(2*a)
  ci1=(-b+(b^2-4*a*c)^0.5)/(2*a)
  return(pmax(ci1,ci2)) ## Be careful that the max does not always correspond to the more sounding solution at low light
}

