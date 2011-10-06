#----this r script get the best fit plane.

library(scatterplot3d) 


#-----upperleft ------------
b=read.table('/home/jghao/research/ccd/imager/flatness_7_29_11/upperleft/flatness_upperleft.txt')
x=b$V1
y=b$V2
z=b$V3

fit=lm(z~x+y)
s3d<-scatterplot3d(x,y,z,pch=16,highlight.3d=T,type="h",angle=-60,box=F,main='left side best fit plane')
s3d$plane3d(fit)




#-----upperright------------
b=read.table('/home/jghao/research/ccd/imager/flatness_7_29_11/upperright/flatness_upperright.txt')

x=b$V1
y=b$V2
z=b$V3
fit=lm(z~x+y)
s3d<-scatterplot3d(x,y,z,pch=16,highlight.3d=T,type="h",angle=-60,box=F,main='left side best fit plane')
s3d$plane3d(fit)

#-----lowerleft------------
b=read.table('/home/jghao/research/ccd/imager/flatness_7_29_11/lowerleft/flatness_lowerleft.txt')

x=b$V1
y=b$V2
z=b$V3
fit=lm(z~x+y)
s3d<-scatterplot3d(x,y,z,pch=16,highlight.3d=T,type="h",angle=-60,box=F,main='left side best fit plane')
s3d$plane3d(fit)


#-----lowerright------------
b=read.table('/home/jghao/research/ccd/imager/flatness_7_29_11/lowerright/flatness_lowerright.txt')

x=b$V1
y=b$V2
z=b$V3
fit=lm(z~x+y)
s3d<-scatterplot3d(x,y,z,pch=16,highlight.3d=T,type="h",angle=-60,box=F,main='left side best fit plane')
s3d$plane3d(fit)


