#######################################################
### The aim of this code is to modify the functions ###
### of the package LeafGasExchange to consider the  ###
### Marquez et al. 2021 theoryLeagGasExchange       ###
#######################################################

library(LeafGasExchange) #https://github.com/TESTgroup-BNL/LeafGasExchange
## See the funciton f.A of the package LeafGasExchange for more info
f.A<-function(PFD,cs,Tleaf,Tair,RH,gcw=10*10^-3,param=f.make.param(),model_diff='Fick'){
  if(!param[['TBM']]%in%c(0,2)){print('Please, use FATES or CLM4.5 TBMs')}
  #Calculation of temperature dependence of the parameters
  Kc=f.arrhenius(param[['KcRef']],param[['KcHa']],Tleaf)
  Ko=f.arrhenius(param[['KoRef']],param[['KoHa']],Tleaf)
  Gstar=f.arrhenius(param[['GstarRef']],param[['GstarHa']],Tleaf)
  Rd=f.modified.arrhenius(PRef=param[['RdRef']],param[['RdHa']],param[['RdHd']],param[['RdS']],Tleaf)
  Vcmax=f.modified.arrhenius(PRef=param[['VcmaxRef']],param[['VcmaxHa']],param[['VcmaxHd']],param[['VcmaxS']],Tleaf)
  Jmax=f.modified.arrhenius(PRef=param[['JmaxRef']],param[['JmaxHa']],param[['JmaxHd']],param[['JmaxS']],Tleaf)
  
  
  I2=PFD*param[['abso']]*(param[['aQY']])
  J=(I2+Jmax-((I2+Jmax)^2-4*(param[['Theta']])*I2*Jmax)^0.5)/(2*(param[['Theta']]))
  
  wi=0.61365*exp(17.502*(Tleaf-273.16)/(240.97+(Tleaf-273.16)))/(param[['Patm']])
  ws=0.61365*exp(17.502*(Tair-273.16)/(240.97+(Tair-273.16)))/(param[['Patm']])*RH/100
  
  ### VPDleaf for USO model in kPa
  ds=(wi-ws)*param[['Patm']]
  
  if(model_diff=='Fick'){k=0;l=0;gcw=0;q=param[['g0']]}
  if(model_diff=='vCF'){k=(wi-ws)/(2-(wi+ws));l=0;gcw=0;q=param[['g0']]} 
  if(model_diff=='MSWF'){k=(wi-ws)/(2-(wi+ws));l=gcw/20;q=param[['g0']]-gcw}
  

  #Resolution for CLM4.5 and FATES
  if(param[['TBM']]%in%c(0,2)){
    
    # Analytical solution of the system of equations {E1 : A=f(ci), E2 : gs=f(A,cs) and ci=f(cs)}
    cic=f.solv(x=Vcmax,y=Kc*(1+param[['O2']]/Ko),cs=cs,Rd=Rd,Gstar=Gstar,q=q,g1=param[['g1']],power=param[['power']],ds=ds,RH=RH,k=k,l=l,model=param[['model.gs']])
    Wc=(cic-Gstar)*Vcmax/(cic+Kc*(1+param[['O2']]/Ko))
    
    cij=f.solv(x=J/4,y=2*Gstar,cs=cs,Rd=Rd,Gstar=Gstar,q=q,g1=param[['g1']],param[['power']],ds=ds,RH=RH,k=k,l=l,model=param[['model.gs']])
    Wj=(cij-Gstar)*J/(4*cij+8*Gstar)
    
    Tp=f.modified.arrhenius(PRef=param[['TpRef']],param[['TpHa']],param[['TpHd']],param[['TpS']],Tleaf)
    Wp=3*Tp
    cip=f.solv(x=3*Tp,y=-Gstar,cs=cs,Rd=Rd,Gstar=Gstar,q=q-gcw,g1=param[['g1']],power=param[['power']],ds=ds,RH=RH,k=k,l=l,model=param[['model.gs']])
    
    ci=cij
    if(!is.null(which(Wc<Wj))&length(cic)==length(cij)){ci[which(Wc<Wj)]=cic[which(Wc<Wj)]}
    if(!is.null(which(Wc<Wj))&length(cic)!=length(cij)){ci[which(Wc<Wj)]=cic}
    W=pmin(Wc,Wj)
    if(!is.null(which(Wp<W))){ci[which(Wp<W)]=cip[which(Wp<W)]}
    Ai=f.smooth(A1 = Wc,A2 = Wj,theta=param[['thetacj']])
    A=f.smooth(A1=Ai,A2=Wp,theta=param[['thetaip']])-Rd
    gs=f.gs(A=A,cs=cs,ds=ds*1000,Rd=Rd,Gstar=Gstar,RH=RH,g0=param[['g0']],g1=param[['g1']],power=param[['power']],model =param[['model.gs']])
    
    Ec=gcw*(wi-ws)
    Es=gs*(wi-ws)/(1-(wi+ws)/2)
    ET=Ec+Es
    output=list(A=A,Ac=Wc-Rd,Aj=Wj-Rd,Ap=Wp-Rd,Ag=A+Rd,gs=gs,ci=ci,ds=ds,ET=ET,Ec=Ec,Es=Es)
    return(output)
  }
}


## Apply only for USO models and BWB model of conductance
f.solv<-function(x,y,cs,Rd,Gstar,g1,power,ds,RH,k,l,q,model){
  if(model=="USO"|model==0){
    m=1.6*(1+g1/(ds)^power)*1/cs
  }else if(model=="USO_simpl"|model==1){
    m=1.6*(g1/(ds)^power)*1/cs
  }else if(model=="BWB"|model==2){
    m=(g1*RH/100)*1/cs
  }else{print('Please, use BWB or USO or USO_simpl conductance model')}
  a=1/1.6*Rd*m+Rd*k*m-1/1.6*m*x-k*m*x-1/1.6*q-k*q-l 
  b=1/1.6*Gstar*m*x+Gstar*k*m*x-1/1.6*Rd*cs*m+Rd*cs*k*m+1/1.6*Rd*m*y+ Rd*k*m*y+1/1.6*cs*m*x-cs*k*m*x+1/1.6*cs*q-cs*k*q-1/1.6*q*y-k*q*y+cs*l- l*y+ Rd-x
  c=-1/1.6*Gstar*cs*m*x+Gstar*cs*k*m*x-1/1.6*Rd*cs*m*y+Rd*cs*k*m*y+1/1.6*cs*q*y-cs*k*q*y+cs*l*y+x*Gstar+Rd*y  
  
  ci2=(-b-(b^2-4*a*c)^0.5)/(2*a)
  ci1=(-b+(b^2-4*a*c)^0.5)/(2*a)
  return(pmax(ci1,ci2)) 
}

