library(Matrix) 
library(quadprog)
library(MASS)
library(geeM)
library(pROC) #package to calculate ROC
library(coda)
library(rjags)
library(Epi)
library(lme4)
kins_all<-read.table("S:\\Data\\Data created\\aim2.kin.csv", header=TRUE,sep=",")
setwd("S:/AIM 2/Programs/JAGS")
mydata <- read.table("S:/Data/Data created/aim2.bin.csv", header=TRUE,sep=",")
#GEE 
gee.exes=geem(ANYCHDS4 ~ age+sbp+ldl+hdl+sex+diab+smoke+hptr+micro+macro,
            id=FID,corstr="ar1", data=mydata,family = binomial)
summary(gee.exes)
gee.indes=geem(ANYCHDS4 ~ age+sbp+ldl+hdl+sex+diab+smoke+hptr+micro+macro,
            id=FID,corstr="independence", data=mydata,family = binomial)
summary(gee.indes)
#Bayesian model
UFID=unique(mydata$FID)
nF=length(UFID)
N=nrow(mydata)
mydata$IDoffset=1:N
Ranks <- with(mydata, ave(IDNO, FID, FUN = function(x) 
  rank(x, ties.method="first")))
mydata1=mydata[Ranks == 1, ]
offset=mydata1$IDoffset
offset[nF+1]=nrow(mydata)+1
family.list=UFID
pheno.file=mydata
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
datamodel=mydata
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
             "b.hptr","b.micro","b.macro","sigma.g2")
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
a
#GEE exchangeable
# Estimates Model SE Robust SE   wald         p
# (Intercept) -5.820000 0.646400  0.644100 -9.036 0.0000000
# age          0.045910 0.005322  0.005326  8.620 0.0000000
# sbp          0.004995 0.004300  0.004859  1.028 0.3040000
# ldl          0.007912 0.002311  0.002128  3.717 0.0002012
# hdl         -0.012090 0.005558  0.005833 -2.072 0.0382600
# sex         -0.384000 0.154700  0.146200 -2.626 0.0086520
# diab         0.545500 0.166000  0.163800  3.330 0.0008690
# smoke        0.311400 0.159400  0.154200  2.020 0.0433900
# hptr         0.399200 0.172600  0.170600  2.340 0.0192600
# micro        0.445200 0.188700  0.201800  2.206 0.0273600
# macro        0.781000 0.280700  0.297600  2.624 0.0086830

#GEE Independent
#              Estimates Model SE Robust SE    wald         p
# (Intercept) -5.694000 0.620200  0.638900 -8.9120 0.0000000
# age          0.044970 0.005108  0.005192  8.6610 0.0000000
# sbp          0.004745 0.004199  0.004801  0.9883 0.3230000
# ldl          0.007773 0.002248  0.002029  3.8310 0.0001275
# hdl         -0.012670 0.005429  0.005791 -2.1870 0.0287200
# sex         -0.363400 0.151100  0.143600 -2.5310 0.0113700
# diab         0.547200 0.161400  0.159200  3.4370 0.0005873
# smoke        0.329000 0.155200  0.152500  2.1570 0.0309800
# hptr         0.400600 0.168900  0.169000  2.3700 0.0177800
# micro        0.456000 0.184500  0.198900  2.2930 0.0218500
# macro        0.829900 0.273500  0.289700  2.8650 0.0041680

#Bayesian
          #Mean       SD  Naive SE Time-series SE
# b.age    0.047656 0.004779 1.511e-04      0.0006629
# b.diab   0.557945 0.175520 5.550e-03      0.0114009
# b.hdl   -0.013500 0.006308 1.995e-04      0.0009961
# b.hptr   0.413009 0.175690 5.556e-03      0.0087623
# b.ldl    0.007975 0.002436 7.704e-05      0.0003155
# b.macro  0.851276 0.294079 9.300e-03      0.0119162
# b.micro  0.466903 0.192151 6.076e-03      0.0075502
# b.sbp    0.004140 0.003171 1.003e-04      0.0006460
# b.sex   -0.394096 0.158882 5.024e-03      0.0086578
# b.smoke  0.331696 0.164566 5.204e-03      0.0073402
# b0      -5.865389 0.566831 1.792e-02      0.1477336
# sigma.g2  0.283158 0.081778 2.586e-03      0.0323311
#           2.5%       25%       50%       75%     97.5%
# b.age    0.039598  0.044403  0.047285  0.050509  0.058809
# b.diab   0.238921  0.438833  0.544822  0.669540  0.936757
# b.hdl   -0.026603 -0.018039 -0.012893 -0.009078 -0.002492
# b.hptr   0.067523  0.288626  0.409931  0.533384  0.738089
# b.ldl    0.002778  0.006458  0.008171  0.009681  0.012453
# b.macro  0.273160  0.662342  0.852016  1.046601  1.386308
# b.micro  0.088977  0.332357  0.465161  0.598079  0.831770
# b.sbp   -0.002296  0.002028  0.004519  0.006074  0.010653
# b.sex   -0.720245 -0.497181 -0.395875 -0.295060 -0.065812
# b.smoke  0.020084  0.218020  0.326945  0.442528  0.661610
# b0      -6.856890 -6.247966 -5.899856 -5.517183 -4.710120
#sigma.g2  0.165586  0.223829  0.272026  0.325771  0.487883

# str=AR1
# > summary(gee.exes)
# Estimates Model SE Robust SE   wald         p
# (Intercept) -5.699000 0.620200  0.637300 -8.943 0.0000000
# age          0.045090 0.005094  0.005193  8.683 0.0000000
# sbp          0.004700 0.004202  0.004791  0.981 0.3266000
# ldl          0.007797 0.002248  0.002015  3.870 0.0001089
# hdl         -0.012740 0.005433  0.005817 -2.190 0.0285500
# sex         -0.359200 0.151400  0.143100 -2.509 0.0120900
# diab         0.548700 0.161500  0.159400  3.443 0.0005743
# smoke        0.331800 0.155300  0.152200  2.180 0.0292700
# hptr         0.402700 0.169000  0.169000  2.383 0.0171900
# micro        0.455300 0.184600  0.198300  2.296 0.0216500
# macro        0.836300 0.273500  0.288800  2.896 0.0037820
# 
# Estimated Correlation Parameter:  -0.01161 
# Correlation Structure:  ar1 
# Est. Scale Parameter:  0.9822 
# 
# Number of GEE iterations: 5 
# Number of Clusters:  91    Maximum Cluster Size:  102 
# Number of observations with nonzero weight:  2674 