#binary outcome 
source("S:/AIM 2/Programs/combination of family simulate binary outcome.R")
library(kinship2)
library(Matrix) #kinship2 loading requires package Matrix and quadprog
library(quadprog)
library(MASS)
library(geeM)
#define other constants
beta_0=1 
beta_1=-0.1
beta_2=3
beta=c(beta_0,beta_1,beta_2)
sig2=3
#calculate beta for the marginal model
expit<-function(rs) {1/(1+exp(-rs))}
pheno.file=gen_data(n1,n2,k2,n3,k3,n4,k4,n5,k5,nf1,nf2,nf3,nf4,nf5)
whole.data=pheno.file
family.list=unique(pheno.file$fid)
names(pheno.file)[names(pheno.file) == 'id'] <- 'subject'
names(pheno.file)[names(pheno.file) == 'fid'] <- 'gpedid'
names(pheno.file)[names(pheno.file) == 'dadid'] <- 'dadsubj'
names(pheno.file)[names(pheno.file) == 'momid'] <- 'momsubj'
pm=NULL
pj_all=NULL
set.seed(1001)
for(i in 1:nF){
  p <-pheno.file[which(pheno.file$gpedid
                       == family.list[i]),]
  mped <-pedigree(id=p$subject,dadid=p$dadsubj,
                  momid=p$momsubj,sex=p$sex,famid
                  =p$gpedid,missid=0)
  kmat <-kinship(mped)
  G1=2*sig2*kmat #calculate G, the cov matrix for each family
  A=as.matrix(G1)#this will change G1 as a matrix A
  df2=p
  x=data.matrix(cbind(1,df2[,c('age','sex1')])) #x is an Nxp design matrix 
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
#calculate the coefficient of the marginal effect 
#prepare the design matrix X and the logit of pm
xd=data.matrix(cbind(1,whole.data[,c('age','sex1')])) #x is an Nxp design matrix 
odds.pm=whole.data$pm/(1-whole.data$pm)
whole.data$logitpm=log(odds.pm)
lpm=data.matrix(whole.data$logitpm)
#coefficient of the marginal effect 
beta_t=as.numeric(ginv(xd)%*%lpm)
# 
# sig2=1 
#> beta_t
#[1]  0.87746321 -0.08720622  2.59138052
#sig2=5
# beta_t
# [1]  0.59836399 -0.05945772  1.76432325
# sig2=10
# > beta_t
# [1]  0.45549178 -0.04543398  1.35462405
sum.coef=NULL
sum.cover=NULL
# Start the clock!
ptm <- proc.time()
for (i in 1:1000)
{set.seed(i)
  whole.data=gen_data(n1,n2,k2,n3,k3,n4,k4,n5,k5,nf1,nf2,nf3,nf4,nf5)
  #gee.es=geem(y ~ age+ sex1,id=fid,corstr="exchangeable", data=whole.data,family = binomial)
  #gee.es=geem(y ~ age+ sex1,id=fid,corstr="independence", data=whole.data,family = binomial)
  gee.es=geem(y ~ age+ sex1,id=fid,corstr="ar1", data=whole.data,family = binomial)
  model.coef=coef(gee.es)
  sum.coef=rbind(sum.coef,model.coef)
  model.sde=summary(gee.es)$se.robust
  LCI95=model.coef-qnorm(.975)*model.sde
  UCI95=model.coef+qnorm(.975)*model.sde
  #judge whether mean is within the 95%CI
  model.cover= beta_t< UCI95 & beta_t>LCI95
  sum.cover=rbind(sum.cover,model.cover)
}

# Stop the clock
proc.time() - ptm

bias=beta_t-colMeans(sum.coef,na.rm=T)
cover.p=colMeans(sum.cover,na.rm=T)
relative.bias= (beta_t-colMeans(sum.coef,na.rm=T)) / abs(beta_t)

bias
cover.p
relative.bias

#02062022 update
# with geeM corr=exchangable
# user  system elapsed 
# 299.53    0.14  300.66 
# > 
#   > bias=beta-colMeans(sum.coef)
# > cover.p=colMeans(sum.cover)
# > relative.bias= (beta-colMeans(sum.coef)) / abs(beta)
# > 
# > bias
# (Intercept)          age         sex1 
# -0.007317379  0.001349522 -0.080095129 
# > cover.p
# (Intercept)         age        sex1 
# 0.937       0.949       0.937 
# > relative.bias
# (Intercept)          age         sex1 
# -0.008339243  0.015475066 -0.030908285 

# 
#with geeM corr=independence
# > bias
# (Intercept)          age         sex1 
# -0.007435813  0.001404280 -0.082539079 
# > cover.p
# (Intercept)         age        sex1 
# 0.940       0.949       0.939 
# > relative.bias
# (Intercept)          age         sex1 
# -0.008474216  0.016102981 -0.031851393 

#AR1 sig=1
#  
# bias
# (Intercept)          age         sex1 
# -0.007430836  0.001365230 -0.079660920 
# > cover.p
# (Intercept)         age        sex1 
# 0.938       0.950       0.939 
# > relative.bias
# (Intercept)         age        sex1 
# -0.00847270  0.01565936 -0.03073604 

# AR1 sig=2
# bias
# (Intercept)          age         sex1 
# -0.015781848  0.002662589 -0.140268985 
# > cover.p
# (Intercept)         age        sex1 
# 0.955       0.950       0.902 
# > relative.bias
# (Intercept)         age        sex1 
# -0.02025461  0.03441864 -0.06117918 
# sig2=3
# > bias
# (Intercept)          age         sex1 
# -0.022183501  0.003841484 -0.187230020 
# > cover.p
# (Intercept)         age        sex1 
# 0.939       0.904       0.840 
# > relative.bias
# (Intercept)         age        sex1 
# -0.03181163  0.05518652 -0.09052599 