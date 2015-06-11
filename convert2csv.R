#######################
##Converting opus files
#######################

## Read me
##This scripts works with spectral files named using a prefix 3-letter prefix followed
##by a six digit number

## 1. To use the script:
## 2.source the read.opus function
source.path<-source('~/Training/India/Spectral_processing/scripts', chdir = TRUE)#On the same computer remains unchanged
opus.file.path<-'~/HTS-xt_spectra/Kiberashi'
output.path<-'~/HTS-xt_spectra'
file<-"Raw spectra.csv"

source('~/Training/India/Spectral_processing/scripts/read.opus.R', chdir = TRUE,echo=TRUE)
## 3.Change path to point to the folder with OPUS files

lst <- as.list(list.files(path=opus.file.path, pattern=".[0-9]$", full.names=TRUE))
spectra<-c()
for ( i in 1:length(lst)){
spec <- opus(lst[[i]], speclib="ICRAF",plot.spectra=TRUE,print.progress=TRUE)
spectra<-rbind(spectra,spec)
}
#View part of spectra
spectra[,1:8]
setwd(output.path) #set working directory where to store converted folder
write.table(spectra,file=file,sep=",",row.names=FALSE)