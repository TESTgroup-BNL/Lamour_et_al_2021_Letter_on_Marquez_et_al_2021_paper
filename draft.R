### Test:
## simulation of an Aci curve using the MSWF model (gcw = 10) and the vCF models
## using the same set of parameters for Vcmax Jmax etc

cs=10:2000
Aci_MSWF=f.A(PFD = 1800,cs = cs,Tleaf = 25+273.15,Tair = 25+273.15,RH = 70,gcw = 0.01,param = f.make.param(g1=3,thetacj = 1,thetaip = 1,TpRef=100),model_diff = 'MSWF')
Aci_vCF=f.A(PFD = 1800,cs = cs,Tleaf = 25+273.15,Tair = 25+273.15,RH = 70,param = f.make.param(g1=3,thetacj = 1,thetaip = 1,TpRef=100),model_diff = 'vCF')

plot(x=Aci_MSWF$ci,y=Aci_MSWF$A,type='l')
lines(x=Aci_vCF$ci,y=Aci_vCF$A,col='red') ## They overlap for ci

plot(x=cs,y=Aci_MSWF$A,type='l')
lines(x=cs,y=Aci_vCF$A,col='red') ## They overlap for ci

## fitting of the simulated Aci curves

reg_MSWF=f.fitting(measures = data.frame(Qin=1800,cs=cs,Tleaf = 25+273.15,A=Aci_MSWF$A,Ci=Aci_MSWF$ci,Species='Test'),id.name = 'Species',
          Start =list(JmaxRef=65,RdRef=1,VcmaxRef=45),
          param = f.make.param(thetacj = 1,thetaip = 1),type = 'Aci')

reg_vCF=f.fitting(measures = data.frame(Qin=1800,cs=cs,Tleaf = 25+273.15,A=Aci_vCF$A,Ci=Aci_vCF$ci,Species='Test'),id.name = 'Species',
          Start =list(JmaxRef=65,RdRef=1,VcmaxRef=45),
          param = f.make.param(thetacj = 1,thetaip = 1),type = 'Aci')


## simulation of an Acs curve using the MSWF and vCF model using the two set of parameters
