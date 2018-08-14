## Realdata

data1<-read.table("C:\\Users\\auz5836\\Desktop\\papers\\realdata\\MLN0002-C13006\\C13006.0_001.txt",head=T)
data2<-read.table("C:\\Users\\auz5836\\Desktop\\papers\\realdata\\MLN0002-C13006\\effpat6.pca.txt",head=T)
data2<-data2[sapply(data1$FID,function(x) grep(x,data2$FID)),]


loc<-intersect(which(is.na(data2$CMYCHG6)==F), which(is.na(data2$CALPRO)==F))
data1<-data1[loc,]
data2<-data2[loc,]
data1<-droplevels(data1)
data2<-droplevels(data2)

x0<-cbind(data2$SEXCD,data2$AGE,data2$CALPRO,data2$BASECMY)
w<-as.integer(data2$IARM)
y<-data2$CMYCHG6
x<-data1[,-c(1:6)]


## Input Missing Data
x<-apply(x,2, function(x) {x[which(is.na(x)==T)]<-mean(x, na.rm=T); return(x)})

## Save original data before scaling
xx0<-x0
ww<-w
yy<-y
xx<-x

## Scale
x<-apply(x,2,scale)
x0<-apply(x0,2,scale)
y<-scale(y)
w<-2*(w-1.5)

## Generate final data
m_X<-dim(x0)[2]
m_W<-1
m_G<-dim(x)[2]
m_I<-m_G
x<-cbind(x0,w,x,w*x)

## Run
x0<-rep(0,dim(x)[2])
lamb_opt_glasso<-3
lamb_opt2_glasso<-2
solg<-FASTA(x,y,f, gradf, g, proxg, x0, tau1, max_iters = 300, w = 10, 
            backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
            eps_n = 1e-15,m_X,m_W,m_G,m_G,lamb_opt_glasso,lamb_opt2_glasso,restart=TRUE)
beta.glasso<-split_beta(solg$x,m_X,m_W,m_G,m_G )
beta.glasso$G_gene<-colnames(xx)[which(beta.glasso$G!=0)]
beta.glasso$I_gene<-colnames(xx)[which(beta.glasso$I!=0)]
