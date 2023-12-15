#Aim 2 binary outcome
#generate data
gen_data<-function(nF,n,beta_0,beta_1,beta_2){
  for(i in 1:nF){
    df=subset(kins_all,FID==UFID[i],) #process one family at one time
    df2=subset(shs.var,FID==UFID[i],)
    df2<-df2[order(df2$IDNO),]#this makes sure the bi was added accordingly
    Fsize=sqrt(length(df$KINCOEF))
    kinmat=matrix(df$KINCOEF,nrow=Fsize)
    G1=2*sig2*kinmat #calculate G, the cov matrix for each family
    A=as.matrix(G1)#this will change G1 as a matrix A
    r=eigen(A) #spectral decomposition, special case of svd
    vectors=r$vectors #eigen vectors
    values=r$values #eigen values  
    #G is the square root of A, i.e, G%*%t(G)=A
    G=vectors%*%diag(sqrt(values))
    ui=rnorm(Fsize)
    bi=G%*%ui
    family.data=data.frame(df2$IDNO,df2$FID,df2$sex,df2$age,bi)
    whole.data=rbind(whole.data,family.data)
  }
  names(whole.data)[names(whole.data) == 'df2.IDNO'] <- 'IDNO'
  names(whole.data)[names(whole.data) == 'df2.age'] <- 'age'
  names(whole.data)[names(whole.data) == 'df2.sex'] <- 'sex'
  names(whole.data)[names(whole.data) == 'df2.FID'] <- 'fid'
  logit.p<-beta_0 + beta_1*whole.data$age + beta_2*whole.data$sex +whole.data$bi
  odds.p <- exp(logit.p)
  p=odds.p/(1+odds.p)
  whole.data$y<-rbinom(n,1,p)
  whole.data
}
library(Matrix) #kinship2 loading requires package Matrix and quadprog
library(quadprog)
library(MASS)
library(geeM)
#library(pROC) #package to calculate ROC
kins_all <- read.table("S:/Data/Data created/kins_noid.csv", header=TRUE,sep=",")
shs.var<- read.table("S:/Data/Data created/shs_noid.csv", header=TRUE,sep=",")
#define other constants
UFID=unique(kins_all$FID)
nF=length(UFID)
n=length(shs.var$IDNO)
sig2=1
beta_0=1 
beta_1=-0.1
beta_2=3
beta=c(beta_0,beta_1,beta_2)

###################################################
#calculate the marginal probabilities
#Inverse logit function to calculate the probability from log odds)
expit<-function(rs) {1/(1+exp(-rs))}
pm=NULL
pj_all=NULL
whole.data=NULL
set.seed(1001)
whole.data=gen_data(nF,n,beta_0,beta_1,beta_2)
summary(whole.data)
for(i in 1:nF){
  df=subset(whole.data,fid==UFID[i],) #process one family at one time
  x=data.matrix(cbind(1,df[,c('age','sex')])) #x is an Nxp design matrix 
  # MC Integrate over the distribution of the random effect for a given patient by family
  # and output population average predictions
  #Sample from the distribution of the random effects by family
  #get the cov matrix A for each family
  df=subset(kins_all,FID==UFID[i],) #process one family at one time
  df2=subset(shs.var,FID==UFID[i],)
  df2<-df2[order(df2$IDNO),]#this makes sure the bi was added accordingly
  Fsize=sqrt(length(df$KINCOEF))
  kinmat=matrix(df$KINCOEF,nrow=Fsize)
  G1=2*sig2*kinmat #calculate G, the cov matrix for each family
  A=as.matrix(G1)#this will change G1 as a matrix A
  #get the marginal p for one family 
  x_beta=c(x%*%(as.matrix(beta))) 
  pj_f=NULL
  for (j in 1: 1000){
    set.seed(j)
    bj=mvrnorm(n=1,mu=rep(0, nrow(df2)),Sigma=A) #draw a random sample of b (vector of random effects)
    rsj=x_beta+bj
    pj=expit(rsj)
    pj_f=cbind(pj_f,pj)
  }
  pj_all=rbind(pj_all,pj_f)
}
pm=rowMeans(pj_all)
#not much difference between individual and marginal mean
whole.data$pm=pm
diff=whole.data$pm-whole.data$y
summary(diff) 
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.990303 -0.191716  0.044846  0.003403  0.237737  0.893802 
#calculate the coefficient of the marginal effect 
#prepare the design matrix X and the logit of pm
xd=data.matrix(cbind(1,whole.data[,c('age','sex')])) #x is an Nxp design matrix 
odds.pm=whole.data$pm/(1-whole.data$pm)
whole.data$logitpm=log(odds.pm)
lpm=data.matrix(whole.data$logitpm)
#coefficient of the marginal effect 
betam=ginv(xd)%*%lpm
#betam
# [,1]
# [1,]  0.89832757
# [2,] -0.08858366
# [3,]  2.64219429
beta_t=c(betam)
##################################################
#GEE performance
sum.coef=NULL
sum.cover=NULL
whole.data=NULL
# Start the clock!
ptm <- proc.time()
for (i in 1:1000)
{set.seed(i)
  whole.data=NULL
  whole.data=gen_data(nF,n,beta_0,beta_1,beta_2)
  #gee.es=geem(y ~ age+ sex,id=fid,corstr="exchangeable", data=whole.data,family = binomial)
  #gee.es=geem(y ~ age+ sex,id=fid,corstr="independence", data=whole.data,family = binomial)
  #gee.es=geem(y ~ age+ sex,id=fid,corstr="unstructured", data=whole.data,family = binomial)
  gee.es=geem(y ~ age+ sex,id=fid,corstr="ar1", data=whole.data,family = binomial)
  model.coef=coef(gee.es)
  sum.coef=rbind(sum.coef,model.coef)
  model.sde=summary(gee.es)$se.robust
  LCI95=model.coef-qnorm(.975)*model.sde
  UCI95=model.coef+qnorm(.975)*model.sde
  #judge whether mean is within the 95%CI
  model.cover= beta_t< UCI95 & beta_t>LCI95
  sum.cover=rbind(sum.cover,model.cover)
}
summary(whole.data)
# Stop the clock
proc.time() - ptm

bias=beta_t-colMeans(sum.coef,na.rm=T)
cover.p=colMeans(sum.cover,na.rm=T)
relative.bias= (beta_t-colMeans(sum.coef,na.rm=T)) / abs(beta_t)

bias
cover.p
relative.bias

#exchangeable
# > proc.time() - ptm
# user  system elapsed 
# 931.48    9.69  942.60 
# > bias
# (Intercept)          age          sex 
# 0.030282647 -0.002580399  0.072580766 
# > cover.p
# (Intercept)       age         sex 
# 0.931       0.888       0.898 
# > relative.bias
# (Intercept)         age         sex 
# 0.03371003 -0.02912951  0.02746988 
#independent
# user  system elapsed 
# 461.90    2.95  465.84 
# bias
# (Intercept)          age          sex 
# 0.031530250 -0.002535435  0.071791981 
# > cover.p
# (Intercept)         age         sex 
# 0.933       0.884       0.900 
# > relative.bias
# (Intercept)         age         sex 
# 0.03509883 -0.02862193  0.02717135 

#AR1   user  system elapsed 
# 594.19   32.85  629.86 
# > 
#   > bias=beta_t-colMeans(sum.coef,na.rm=T)
# > cover.p=colMeans(sum.cover,na.rm=T)
# > relative.bias= (beta_t-colMeans(sum.coef,na.rm=T)) / abs(beta_t)
# > 
#   > bias
# (Intercept)         age         sex 
# 0.03219806 -0.00256051  0.07168580 
# > cover.p
# (Intercept)         age         sex 
# 0.932       0.884       0.905 
# > relative.bias
# (Intercept)         age         sex 
# 0.03584223 -0.02890499  0.02713116 
# > 

