library(Matrix) 
library(quadprog)
library(MASS)
library(geeM)
library(pROC) #package to calculate ROC
library(coda)
library(rjags)
library(Epi)
library(lme4)
datacv <- read.table("S:/Data/Data created/aim2.bin.csv", header=TRUE,sep=",")
#str(aim2.bin)
#5CV 
############################GEE MODEL##################################
#loop to calculate the average of c statistics for GEE model
expit<-function(rs) {1/(1+exp(-rs))}
c_all=NULL
for(i in 1:5){
  datatest <- datacv[which(datacv$group==i), ]
  datamodel<- datacv[which(datacv$group!=i), ]
  gee.es=geem(ANYCHDS4 ~ age+sbp+ldl+hdl+sex+diab+smoke+hptr+micro+macro,
              id=FID,corstr="ar1", data=datamodel,family = binomial)
  #gee.es=geem(ANYCHDS4 ~ age+sbp+ldl+hdl+sex+diab+smoke+hptr+micro+macro,
             # id=FID,corstr="independence", data=datamodel,family = binomial)
  #using the estimated coefficient to calculate the pm_gee on test data
  model.coef=as.matrix(coef(gee.es))
  x=data.matrix(cbind(1,datatest[,c('age','sbp','ldl','hdl','sex','diab','smoke',
                                     'hptr','micro','macro')])) #x is an Nxp design matrix 
  rsj=c(x%*%(model.coef)) 
  datatest$pp=expit(rsj)
  #calculate the ROC 
  model.roc=roc(response = datatest$ANYCHDS4, predictor = datatest$pp)
  model.c=model.roc$auc
  c_all=rbind(c_all,model.c)
}
gee.c=colMeans(c_all)
# > gee.c #str=independence
# [1] 0.7939483
# > gee.c #str=exchangeable
# [1] 0.7944331
# > gee.c #str=AR1
# [1] 0.79389
#create a auc curver
# rocobj <- plot.roc(datatest$ANYCHDS4, predictor = datatest$pm_gee,
#                    main = "Confidence intervals", 
#                    percent=TRUE,
#                    ci = TRUE,                  # compute AUC (of AUC by default)
#                    print.auc = TRUE)           # print the AUC (will contain the CI)

############################BAYESIAN MODEL###############
kins_all<-read.table("S:\\Data\\Data created\\aim2.kin.csv", header=TRUE,sep=",")
setwd("S:/AIM 2/Programs/JAGS")
#create a vector offset using current order of data
#The JAGS code
###################jags file#####################################
## this cat function will save these code in jags file
cat( "model
 {
 for( i in 1:n.fam)
 {                    ## loop over individuals within each family
    for(j in offset[i]:(offset[i+1]- 1))
    {
      v[j]<-inprod(G[j,offset[i]:(offset[i+1]-1)], u[offset[i]:(offset[i+1]-1)])  ## G *%* u
      y[j] ~dbern(mu[j])
      logit(mu[j]) <-b0+b.age*age[j] +b.sex*sex[j] +b.sbp*sbp[j] +b.smoke*smoke[j]+ b.ldl*ldl[j]+ b.hdl*hdl[j]
                    +b.diab*diab[j] + b.hptr*hptr[j]+ b.macro*macro[j] +b.micro*micro[j] +v[j]

    }
 }
  ## Model random effects as univariate normal
  for( t in 1:N){u[t]~dnorm(0,tau.g)}
  ## priors for fixed effects
  b0 ~dnorm(0,0.001)
  b.age ~dnorm(0,0.001)
  b.sex ~dnorm(0,0.001)
  b.sbp ~dnorm(0,0.001)
  b.smoke ~dnorm(0,0.001)
  b.ldl ~dnorm(0,0.001)
  b.hdl ~dnorm(0,0.001)
  b.diab ~dnorm(0,0.001)
  b.hptr ~dnorm(0,0.001)
  b.macro ~dnorm(0,0.001)
  b.micro ~dnorm(0,0.001)
  ## varance components
  tau.g ~dgamma(1,1)
  sigma.g2 <-1/tau.g
 }",
     file="m5.jag"
)

###########################c statistic using marginal probability
ptm <- proc.time()
c_bay_all=NULL
for(i in 1:5){
  datatest <- datacv[which(datacv$group==i), ]
  datamodel<- datacv[which(datacv$group!=i), ]
  UFID=unique(datamodel$FID)
  nF=length(UFID)
  mydata=datamodel
  N=nrow(mydata)
  mydata$IDoffset=1:N
  Ranks <- with(mydata, ave(IDNO, FID, FUN = function(x) 
    rank(x, ties.method="first")))
  mydata1=mydata[Ranks == 1, ]
  offset=mydata1$IDoffset
  offset[nF+1]=nrow(mydata)+1
  family.list=UFID
  pheno.file=datamodel
  names(pheno.file)[names(pheno.file) == 'IDNO'] <- 'subject'
  names(pheno.file)[names(pheno.file) == 'FID'] <- 'gpedid'
  my.U <-vector("list")
  my.S <-vector("list")
  my.kmat <-vector("list")
  my.G <-vector("list")
  for (i in 1:length(family.list)){
    df=subset(kins_all,FID==family.list[i],) #process one family at one time
    df2=subset(pheno.file,gpedid==family.list[i],)
    FamIDs=df2$subject
    Fsize=sqrt(length(df$KINCOEF))
    kinmat=matrix(df$KINCOEF,nrow=Fsize)
    kmat=2*kinmat #calculate G, the cov matrix for each family
    row.names(kmat)<-FamIDs
    ## Now get kinship submatrix for subjects
    ind.subj <-match(pheno.file$subject
                     [which(pheno.file$gpedid ==
                              family.list[i])], row.names(kmat))
    test <-svd(kmat[ind.subj[which(is.na (ind.subj)==F)],
                    ind.subj[ which(is.na(ind.subj)==F)]])
    my.U[[i]] <-test$u
    my.S[[i]] <-test$d
    my.G[[i]] <-test$u%*% diag(sqrt(test$d))
    my.kmat[[i]] <-kmat[ind.subj
                        [ which(is.na(ind.subj)==F)],
                        ind.subj[ which(is.na(ind.subj)==F)]]
  }
  ## This is the matrix to be loaded to bugs
  G.mat <-my.G[[1]]
  for(i in 2:length(my.G)){
    ith.block.1 <-rep(0,nrow(G.mat)*ncol(my.G[[i]]))
    dim(ith.block.1) <-c(nrow(G.mat),ncol(my.G[[i]]))
    ith.block.2 <-rep(0,ncol(G.mat)*nrow(my.G[[i]]))
    dim(ith.block.2) <-c(nrow(my.G[[i]]),ncol(G.mat))
    G.mat <-rbind(cbind(G.mat,ith.block.1),
                  cbind(ith.block.2, my.G[[i]]))
  }
  glmm.es = glmer(ANYCHDS4 ~ age+sbp+ldl+hdl+sex+diab+smoke+hptr+micro+macro+(1|FID), 
                  data = datamodel,family=binomial)
  b=summary(glmm.es)
  beta_hat=coef(b)
  c=as.data.frame(b$varcor)#get the variance of random effects (c[1,4]) 
  reg.dat <- list( y=datamodel$ANYCHDS4,
                   age=datamodel$age,
                   sbp=datamodel$sbp,
                   ldl=datamodel$ldl,
                   hdl=datamodel$hdl,
                   sex=datamodel$sex,
                   diab=datamodel$diab,
                   smoke=datamodel$smoke,
                   hptr=datamodel$hptr,
                   micro=datamodel$micro,
                   macro=datamodel$macro,
                   n.fam=nF,
                   N=N,
                   offset=offset,
                   G=G.mat) 
  reg.par <- c("b0","b.age","b.sbp","b.ldl","b.hdl","b.sex","b.diab","b.smoke",
               "b.hptr","b.micro","b.macro")
  if (c[1,4]==0){
    reg.ini <- list( b0=beta_hat[1], 
                     b.age=beta_hat[2], 
                     b.sbp=beta_hat[3],
                     b.ldl=beta_hat[4],
                     b.hdl=beta_hat[5],
                     b.sex=beta_hat[6],
                     b.diab=beta_hat[7],
                     b.smoke=beta_hat[8],
                     b.hptr=beta_hat[9],
                     b.micro=beta_hat[10],
                     b.macro=beta_hat[11])
  }else{reg.ini <- list(b0=beta_hat[1], 
                        b.age=beta_hat[2], 
                        b.sbp=beta_hat[3],
                        b.ldl=beta_hat[4],
                        b.hdl=beta_hat[5],
                        b.sex=beta_hat[6],
                        b.diab=beta_hat[7],
                        b.smoke=beta_hat[8],
                        b.hptr=beta_hat[9],
                        b.micro=beta_hat[10],
                        b.macro=beta_hat[11]
                        ,tau.g=1/c[1,4])
  }

  reg.mod <- jags.model( file = "m5.jag",
                         data = reg.dat,
                         n.chains = 1,
                         inits = reg.ini,
                         n.adapt = 1000)
  reg.res <- coda.samples( reg.mod,
                           var = reg.par,
                           n.iter = 2000,
                           thin = 2)
  a=summary( reg.res )

  #using the estimated coefficient to calculate the predicted probability on test data
  model.coef=as.matrix(a$statistics[,1])
  x=data.matrix(cbind(datatest[,c('age','diab','hdl','hptr','ldl','macro','micro','sbp','sex','smoke')],1)) #x is an Nxp design matrix 
  rsj=c(x%*%(model.coef)) 
  datatest$pp=expit(rsj)
  #calculate the ROC 
  model.roc=roc(response = datatest$ANYCHDS4, predictor = datatest$pp)
  model.c=model.roc$auc
  c_bay_all=rbind(c_bay_all,model.c)
}
Bayesian.c=colMeans(c_bay_all)
# Stop the clock
proc.time() - ptm
Bayesian.c
# > Bayesian.c
# [1] 0.7940024
# 





