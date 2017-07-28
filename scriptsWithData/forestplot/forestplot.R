data=read.csv("forest.csv")

library(meta)


data$Pathology_xingyang=as.numeric(Pathology_xingyang)
data$Surgery_xingyang=as.numeric(data$Surgery_xingyang)
levels(data$T_xingyang)[5]
data$T_xingyang=as.numeric(data$T_xingyang)
data$TNM=as.numeric(data$TNM)
data$ERCC1=as.numeric(data$ERCC1)
data$MSH2=as.numeric(data$MSH2)
data$KI67=as.numeric(data$KI67)
data$Site=as.numeric(data$Site)

attach(data)

temp=t(data[,name])
temp=cbind(name,temp)
colnames(temp)[2:ncol(temp)]=c(2:ncol(temp))

detach(data)
temp=data.frame(temp)
attach(temp)

mn<-metabin(ee,en,ce,cn,studlab=paste(data[,1]),method="MH",sm="RR",data=data)
mn

forest(mn)
