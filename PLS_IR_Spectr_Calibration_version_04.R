#Load all the packages required
library(soil.spec)
library(pls)
library(chemometrics)
library(e1071)
#Start by specifying the folder to keep results. For me I chose "/Regional training/Nigeria/" and declare it as:

#Note for windows OS ensure the complete address path is defined inclusing the computer driver letter

results.folder<-"~/Ir-training/calibrations"

#Read into R both IR and chem data

raw<-read.table(file.choose(),sep=",",header=T)

#Remove sample with ssn "icr075956"; this is AfSIS sample with high CN 

nraw<-names(raw)
w2<-substr(colnames(raw[,200:250]),1,1)[20]

g<-which(nraw==paste(w2,2379.8,sep=""))
h<-which(nraw==paste(w2,2350.8,sep=""))

ifelse(length(g)>0,"The CO2 region 2379.8 to 2350.8 cm-1 will be removed","No CO2 bands found in this IR data")

ifelse(length(g)>0,raw<-raw[,-c(g:h)],raw<-raw)


chem<-read.table(file.choose(),sep=",",header=T)
summary(chem)
##################################################################
#Create  folders for storing PLS plots and predicted results     #
	
dir.create(paste(results.folder,"PLS plots",sep="/"),showWarnings=FALSE)

dir.create(paste(results.folder,"PLS residual plots",sep="/"),showWarnings=FALSE)
	
dir.create(paste(results.folder,"Predictions",sep="/"),showWarnings=FALSE)

dir.create(paste(results.folder,"Calibration models",sep="/"),showWarnings=FALSE)


pls.d<-paste(results.folder,"PLS plots",sep="/")
resid.d<-paste(results.folder,"PLS residual plots",sep="/")
mvr.tn.d<-paste(results.folder,"Predictions",sep="/")
calib.d<-paste(results.folder,"Calibration models",sep="/")

#Remove the first column from raw spectra table which contains SSN/SSN and also the first and last columns of the spectra

rh<-colnames(raw)
rhw<-substr(rh,1,1)
l<-which(rhw==w2)[2]
z<-ncol(raw)-1

raw.cutoff<-raw[,l:z]

#Set the spectra table as a matrix
mim<-as.matrix(raw.cutoff)
dim(mim)
#Set the colnames of mim matrix as numeric
colnames(mim)<-as.numeric(substr(colnames(mim),2,15))

#msc preprocessing
raw.msc<-msc(mim)
dim(raw.msc)

#Combine msc with the ssn
mscssn<-cbind(as.vector(raw[,1]),raw.msc)
#View the first four columns and rows
mscssn[1:4,1:4]

#give label to the ssn column
colnames(mscssn)<-c("SSN",colnames(raw.msc))

mscssn<-as.data.frame(mscssn)
setwd(results.folder)
write.table(mscssn,"Mulplicative Scatter corrected spectra.csv",sep=",",row.names=F)


#Preprocess the raw spectra by first derivative; use library soil.spec
de<-trans(mim,tr="derivative",order=1,gap=21)
der<-as.matrix(de$trans)

dc<-ncol(der)

#Combine derivative with the ssn
derssn<-cbind(as.vector(raw[,1]),der)
#View the first four columns and rows
derssn[1:4,1:4]

#give label to the ssn column
colnames(derssn)<-c("SSN",colnames(der))

der<-as.data.frame(derssn)
setwd(results.folder)
write.table(derssn,"First derivative.csv",sep=",",row.names=F)

#Select either the first derivative or MSC file to use for calibration models
files<-c("First derivative.csv","Mulplicative Scatter corrected spectra.csv")
spe<-menu(files,graphics=TRUE,title="Select IR data")

#Read the file selected
der<-read.csv(paste(results.folder,files[spe],sep="/"),sep=",",header=T)
###----Merge first derivative with the chemical data

#First explore chem data and decide on the best transformation method

wtc<-colnames(chem)

#################################################
#Split the plot region into 4
par(mfrow=c(2,2))
plot.var<-menu(wtc,graphics=TRUE,title="Select variable to check distribution")
plotv<-as.vector(subset(chem,select=wtc[plot.var]))
hist(na.omit(plotv[-1,1]),main=wtc[plot.var])
plot(density(sqrt(na.omit(plotv[-1,1]))),main=wtc[plot.var])
################################################
chem[1:4,1:4]
#der$SSN<-paste("icr0",substr(der$SSN,4,13),sep="")
calibset<-merge(chem,der,by.y="SSN",by.x="SSN")
#calibset<-merge(chem,der,by.y="SSN",by.x="oSSN")

#To know the cutoff regions get the columns of the calibset
colnames(calibset)

c<-ncol(chem)
wavenumbers<-as.numeric(substr(colnames(calibset[,-c(1:c)]),2,14))
y<-select.list(as.character(wavenumbers),graphics=TRUE,title="Select region where MIR begins and ends",multiple=TRUE)

ynum<-as.numeric(y)
ynumg<-max(ynum)
ynumn<-min(ynum)


a<-which(wavenumbers==ynumg)+c
b<-which(wavenumbers==ynumn)+c

#cutspec should lie in the regions 4000 to 600
cutspec<-as.matrix(calibset[,a:b])

y<-ncol(cutspec)

#Quick view of merged data
cutspec[1:4,1:4];cutspec[1:4,(y-4):y]

#Obtain from the calibset, the part with chem data
calibset[1:5,1:c]

#Where does the reference chem start?
refh<-colnames(calibset[,1:c])

z<-select.list(as.character(refh),graphics=TRUE,title="Select the first column where reference data starts",multiple=FALSE)
refh.z<-which(refh==z)

cutchem<-as.data.frame(calibset[,refh.z:c])

colnames(cutspec)<-as.numeric(substr(colnames(cutspec),2,12))

#Create the prediction set from the first derivative table; ensuring the same regions in the cutspec are in the prediction set
der[1:4,1:4]

ps<-as.matrix(der[,(a-(c-1)):(b-(c-1))])

par(mfrow=c(1,1))
colnames(ps)<-colnames(cutspec)
dim(cutspec)
#############
#Natural log
############
summ.all2<-c()
for(q in 1:ncol(cutchem)){
setwd(pls.d)

hd<-colnames(cutchem)

png(paste(hd[q],"calibration.png"),height=480,width=480,units="px")


cutspe<-as.data.frame(cutspec)
#Exclude constants
ifelse(min(na.omit(cutchem[,q]))==max(na.omit(cutchem[,q])),q<-q+1,q<-q)
#Determine the max pcs to be tested btn 1:20
cutchem[,q]<-ifelse(cutchem[,q]>0,cutchem[,q],NA)
dvcl<-nrow(cutchem)+1
dvcl<-length(na.omit(cutchem[,q]))
pc<-c(1:20)
jj<-length(na.omit(cutchem[,q]))
	te<-c()
 for (u in 1:length(pc)) {
 	te[u]<-jj/pc[u]}
 
te.sel<-which(te>2)
p<-te.sel[length(te.sel)]

ifelse(jj<32,p<-10,p<-13)

cpe<-na.omit(cbind(as.vector(cutchem[,q]),cutspec))
colnames(cpe)<-c(hd[q],colnames(cutspec))
cpe[1:4,1:4]

ref<-as.matrix(cpe[,1])#response
cutspe<-as.data.frame(cpe[,-1])#predictors
#Perform double cross-validation if calibration dataset<501 spectra
ifelse(jj>0,res.pls<-mvr_dcv(log(ref)~.,data=cutspe,plot.opt=FALSE,ncomp=p,method="widekernelpls",repl=3,selstrat="hastie",segments0=4,segment0.type="random",na.action=na.omit,segments=4,segment.type="random"),"Do nothing")

af.min<-3
#Ensure we do not fit a model with less than 3 PCs.
nc<-ifelse(res.pls$afinal<af.min,res.pls$afinal<-af.min,res.pls$afinal)

#ifelse(jj<dvcl,af<-res.pls$afinal,af<-fx)

#Run the best model
mvr.tn<-plsr(na.omit(log(as.numeric(cutchem[,q]))~cutspec,ncomp=nc,validation ="CV"))

rmse<-RMSEP(mvr.tn)$val[nc]
se<-round(rmse,2)

#mvr.tn<-mvr.tn

#a<-min(exp(mvr.tn$fitted),na.omit(as.numeric(cutchem[,q])))

a<-round(min(na.omit(as.numeric(cutchem[,q])),na.omit(exp(mvr.tn$fitted.values[,,nc]))),1)

am<-round(max(na.omit(as.numeric(cutchem[,q])),na.omit(exp(mvr.tn$fitted.values[,,nc]))),1)

#am<-round(max(exp(mvr.tn$fitted),na.omit(as.numeric(cutchem[,q])),1))

b<-am+(0.05*am)
#b<-am
    #y<-na.omit(as.numeric(cutchem[,q]))
    #x<-na.omit(exp(mvr.tn$fitted.values[,,nc]))
    #xx<-na.omit(exp(mvr.tn$fitted.values[,,nc]))
    #xxx<-exp(mvr.tn$fitted.values[,,nc])
plot(na.omit(as.numeric(cutchem[,q]))~na.omit(exp(mvr.tn$fitted.values[,,nc])),pch=19,col="blue",ylab="",xlab="",main="",ylim=c(a,b),xlim=c(a,b))
bl<-lm(na.omit(as.numeric(cutchem[,q]))~na.omit(exp(mvr.tn$fitted.values[,,nc])))
fl.a<-coef(bl)[1]; fl.b<-coef(bl)[2]
abline(a=0,b=1,col="red",lwd=1,lty=1)
#abline(a=fl.a,b=fl.b,col="red",lwd=1)
r<-round(summary(bl)$r.squared,2)

#Another method for getting r squared-values
#require(stats) 
#r<-round(cor(exp(predicted$fitted.values[,,nc]),na.omit(cutchem[,q]),method="pearson")^2,2)
mtext(paste(hd[q]," (using ",nc," PCs); ","n=",jj,sep=""),side=3,line=0.7,cex=1)
mtext(paste("predicted",sep=" "),side=1,line=2.4,cex=0.8)
mtext(paste("Measured ",sep=" "),side=2,line=2.4,cex=0.8)
#legend("topleft",c(paste("r.squared",r,sep="="),paste("rmse=",se,sep=""),"Best line of fit","1 to 1 line"),lty=c(1,1,1,2),col=c("white","white",2,2),bty="n",cex=0.7)
legend("topleft",c(paste0("rsq = ", r),paste0("rmse=",se),"1 to 1 line"),lty=c(1,1,1),col=c("white","white",2),bty="n",cex=0.7)

dev.off()

#Calculate residuals
resid<-na.omit(as.numeric(cutchem[,q]))-na.omit(exp(mvr.tn$fitted.values[,,nc]))

#Look at the residuals
#a<-min(exp(mvr.tn$fitted),na.omit(as.numeric(cutchem[,q])))

aa<-round(min(na.omit(as.numeric(cutchem[,q])),na.omit(resid)),1)

amm<-round(max(na.omit(as.numeric(cutchem[,q])),na.omit(resid)),1)
m.r<-mean(na.omit(resid))
m.d<-mean(na.omit(cutchem[,q]))

bb<-amm+(0.05*amm)

#Set working directory for the residuals
setwd(resid.d)
png(paste(hd[q],"residuals.png"),height=480,width=480,units="px")

plot(na.omit(as.numeric(cutchem[,q]))~resid,pch=19,col="grey",ylab="",xlab="",main="",ylim=c(aa,bb),xlim=c(aa,bb))
abline(v=m.r,lty=2,col="grey")
abline(h=m.d,lty=2,col="grey")

require(stats)
mtext(paste("Residuals",sep=" "),side=1,line=2.4,cex=0.8)
mtext(paste("Measured ",sep=" "),side=2,line=2.4,cex=0.8)

dev.off()
#######################################################

pds<-c()
summ2<-as.vector(c(paste(hd[q]," (",jj,")",sep=""),r,nc,se))
summ.all2<-rbind(summ.all2,summ2)
#To predict write the general syntax as: predict("model","new spectra"," optimal no. of PCs  selected in the model")
colnames(ps)<-colnames(cutspec)
pds<-predict(mvr.tn,ps,nc)
ssn<-as.matrix(der[,1])
pdss<-cbind(ssn,pds)
colnames(pdss)<-c("SSN",hd[q])

#Get SSN for residuals from the calibset
sn<-as.matrix(as.matrix(calibset[,1]))
resid.ssn<-cbind(sn,resid)
colnames(resid.ssn)<-c("SSN",paste(hd[q],"residuals",sep="."))

#Save loadings plot
#Create wd

setwd(pls.d)
lod<-as.matrix(mvr.tn$loadings)

wavl<-as.numeric(colnames(cutspec))

png(filename=paste(hd[q],"Loadings for PC ",nc,".png",sep=""),height=480,width=480,units="px")

ro<-round((nc/2+0.4),0)
ror<-round(sqrt(ro),0)+1
par(mfrow=c(ror,ror))

for (z in 1:nc){
plot(wavl,lod[,z],type="l",col=z,xlim=c(max(wavl),min(wavl)),ylab=paste("PC ",z," loading",sep=""),xlab=expression("Wavenumbers cm"^-1),main=hd[q])}
dev.off()

setwd(mvr.tn.d)
write.table(pdss,file=paste(hd[q],".csv"),sep=",",row.names=F)

setwd(resid.d)
write.table(resid.ssn,file=paste(hd[q],".csv"),sep=",",row.names=F)

#Save model into the calibration model folder for future predictions
setwd(calib.d)
save(mvr.tn,file=paste(hd[q],".rda",sep=""))
}

#Write model summary into a file
setwd(pls.d)
colnames(summ.all2)<-c("Property (n)","r-squared","PCs","rmse")
write.table(summ.all2,file="Model summary.csv",sep=",",row.names=F)


######################################################################
#Get the prediction and output in one file
#setwd(mvr.tn.d)
setwd(mvr.tn.d)
pd<-list.files(pattern=".csv");
al.p<-which(pd=="All predictions.csv")
ifelse(al.p>0,pd<-pd[-al.p],pd<-pd)

tt<-as.vector(pd)
ttt<-as.vector(strsplit(tt, " .csv"))
#tt<-substr(tt,11,40)
pdm<-read.table(pd[1],sep=",",header=T)
		pdm<-as.matrix(pdm[,1])
for(k in 1:length(pd)){
	pda<-read.table(pd[k],sep=",",header=T)
		pdm<-cbind(pdm,exp(pda[,2]))}
		colnames(pdm)<-c("SSN",ttt)
write.table(pdm,file="All predictions.csv",sep=",",row.names=F)

#Get the residuals and output in one file
setwd(resid.d)
pd<-list.files(pattern=".csv");
al.r<-which(pd=="All residuals.csv")
ifelse(al.r>0,pd<-pd[-al.r],pd<-pd)
tt<-as.vector(pd)
ttt<-as.vector(strsplit(tt, " .csv"))
#tt<-substr(tt,11,40)
pdm<-read.table(pd[1],sep=",",header=T)
		pdm<-as.matrix(pdm[,1])
for(k in 1:length(pd)){
	pda<-read.table(pd[k],sep=",",header=T)
		pdm<-cbind(pdm,pda[,2])}
		colnames(pdm)<-c("SSN",ttt)
write.table(pdm,file="All residuals.csv",sep=",",row.names=F)
