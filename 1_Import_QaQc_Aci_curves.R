

library(LeafGasExchange)
source('0_NewTheoryGasex.R')
setwd("~/GitHub/NewTheoryGasEx")

ls.files=dir('~/GitHub/NewTheoryGasEx/Datasets/ACi_Panama_2020/',recursive = TRUE,full.names = TRUE)

remove_Aci_dupplicate<-function(Aci_curve){
  Aci_curve=round(Aci_curve,-1)
  remove=rep('NO',length(Aci_curve))
  remove[length(Aci_curve)]='YES'
  remove[which(duplicated(Aci_curve[1:(length(Aci_curve)-1)],fromLast = TRUE))]='YES'
  remove[1:5]='NO'
  remove[which(duplicated(Aci_curve[1:5],fromLast = TRUE))]='YES'
  return(remove)
}

Aci_data=do.call("rbind", apply(X = as.matrix(ls.files),FUN = f.import_licor6800,MARGIN = 1))
colnames(Aci_data)['X.U.0394.Pcham']='X.Pcham'

## A binary column is added to the file to detect the points that should be removed for the fitting procedures
Aci_data$remove='NO'
## All the Aci_data with Ci<0 or Ci>2000 are considered wrong
Aci_data[Aci_data$Ci<0|Aci_data$Ci>2000,'remove']='YES'

## The duplicated points are removed using this function
remove=ave(Aci_data$CO2_r,Aci_data$file,FUN=remove_Aci_dupplicate)
Aci_data[which(remove=='YES'),'remove']='YES'

data$remove=factor(data$remove)
pdf(file='QA_QC_Aci_Panama_2020.pdf',)
par(mar = c(5,5,2,5))
by(data = Aci_data,INDICES = Aci_data$Barcode,FUN = function(x){
  
  mean_temp=mean(x$Tleaf)
  max_A=max(x$A)  
  b=28
  a=(32-28)/max_A
  
  #(ggplot(data=x,aes(x=Ci,y=A,color=remove,label=Obs))+ggtitle(unique(x$Leaf_Barcode))
  # +geom_point(size=4,alpha=0.5)+geom_text(aes(label=Obs),col='black')
  # +geom_line(data=x,aes(x=Ci,y=(Tleaf-b)/a),type='dotted',color='grey')+theme_bw() 
  # +scale_y_continuous(sec.axis = sec_axis(~ . *a+b, name = "Leaf Temp")
  # ))
  max_ci=max(x$Ci)
  if(max(x$Ci)>2000){max_ci=1800}
  plot(x=x$Ci,y=x$A,main=unique(x$Barcode),xlab="Ci",ylab='A',cex=2,xlim=c(0,max_ci))
  points(x=x[x$remove=='NO','Ci'],y=x[x$remove=='NO','A'],col='red',cex=2)
  text(x=x$Ci,y=x$A, labels=x$obs,cex=0.7)
  try(text(x=x[x$remove=='NO','Ci'],y=x[x$remove=='NO','A'],labels=x[x$remove=='NO','Obs'],col='red',cex=0.7))
  par(new=TRUE)
  plot(x=x$Ci,y=x$Tleaf,col='grey',axes=F, xlab=NA, ylab=NA,ylim=c(27,32),xlim=c(0,max_ci))
  axis(side = 4)
  mtext(side = 4, line = 3, 'Tleaf')
})
dev.off()

ls_remove=c('BNL17516','BNL17565','BNL16587','BNL17659','BNL17662','BNL18443')
Aci_data=Aci_data[!Aci_data$Barcode%in%ls_remove,]
Aci_data[Aci_data$Barcode=='BNL17195'&Aci_data$obs%in%c(13),'Remove']='YES'
Aci_data[Aci_data$Barcode=='BNL17197'&Aci_data$obs%in%c(2,5),'Remove']='YES'
Aci_data[Aci_data$Barcode=='BNL17521'&Aci_data$obs%in%c(18),'Remove']='YES'
Aci_data[Aci_data$Barcode=='BNL17534'&Aci_data$obs%in%c(6),'Remove']='YES'
Aci_data[Aci_data$Barcode=='BNL17540'&Aci_data$obs%in%c(6),'Remove']='YES'
Aci_data[Aci_data$Barcode=='BNL17639'&Aci_data$obs%in%c(6),'Remove']='YES'
Aci_data[Aci_data$Barcode=='BNL17542'&Aci_data$obs%in%c(5),'Remove']='YES'
Aci_data[Aci_data$Barcode=='BNL17741'&Aci_data$obs%in%c(6),'Remove']='YES'

Aci_data=Aci_data[Aci_data$remove=='NO',]

save(Aci_data,file='1_Aci_data_QaQc.Rdata')