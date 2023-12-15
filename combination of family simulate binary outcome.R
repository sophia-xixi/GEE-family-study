
gen_sing<-function(n1){
  id=(1:n1)
  fid=(1:n1)
  pid=rep(1,n1)
  momid=rep(0,n1)
  dadid=rep(0,n1)
  sex1=rbinom(n1,1,0.5)
  age=runif(n1,7,40)
  sig = matrix(0.5, n1, n1)+0.5*diag(n1)
  #bi=mvrnorm(n=1, mu = rep(0, n1), Sigma = sig2*diag(n1)) #no random effect in singlton family
  bi=0
  epsi=rnorm(n1)
  whole.data1=data.frame(id,fid,pid,momid,dadid,sex1,age,bi,epsi)
  whole.data1$sex[whole.data1$sex1==1]="M"
  whole.data1$sex[whole.data1$sex1==0]="F"
  whole.data1
}
#generate nuclear data
gen_nuclear<-function(n2,k2,nf2){
  #the first few lines are used to build the covariance matrix
  df=data.frame(id=c(1,2,3,4),dadid=c(0,0,1,1),momid=c(0,0,2,2),sex=c("M","F","M","F"),famid=c(1,1,1,1))
  pedlist = with(df,pedigree(id=id,dadid=dadid,momid=momid,sex=sex,famid=famid))
  kinmat=kinship(pedlist)
  G1=2*sig2*kinmat #calculate G, the cov matrix for each family
  A=as.matrix(G1)#this will change G1 as a matrix A
  r=eigen(A) #spectral decomposition, special case of svd
  vectors=r$vectors #eigen vectors
  values=r$values #eigen values
  G=vectors%*%diag(sqrt(values))
  for (i in (nf1+1):(nf1+nf2)) #a total of n/k families(loops),also considering the n1 families in above
  {#create a series of id
    id=c( (n1+(i-nf1-1)*k2+1):(n1+(i-nf1)*k2))
    fid=rep(i,k2)
    pid=c(1:k2)
    momid=c(0,0,n1+(i-nf1-1)*k2+2,n1+(i-nf1-1)*k2+2)
    dadid=c(0,0,n1+1+(i-(1+nf1))*k2,n1+1+(i-(1+nf1))*k2)
    #sex and age, sex is fixed, age is randomly selected from unif [30,40] and unif[7,11] for parents and kid respectly
    sex=c("M","F","M","F") 
    sex1=c(1,0,1,0)
    age=c(runif(1,30,40), runif(1,30,40), runif(1,7,11),runif(1,7,11))
    #random effect
    ui=rnorm(k2)
    bi=G%*%ui
    epsi=rnorm(k2)
    family.data2=data.frame(id,fid,pid,momid,dadid,sex,sex1,age,bi,epsi)
    whole.data2=rbind(whole.data2,family.data2)
  }
  whole.data2
}
#generate one-trio data
gen_data_onetrio<-function(n3,k3,nf3){
  df=data.frame(id=c(1:5),dadid=c(0,0,1,0,3),momid=c(0,0,2,0,4),sex=c("M","F","M","F","F"),famid=rep(1,5))
  pedlist = with(df,pedigree(id=id,dadid=dadid,momid=momid,sex=sex,famid=famid))
  kinmat=kinship(pedlist)
  G1=2*sig2*kinmat #calculate G, the cov matrix for each family
  A=as.matrix(G1)#this will change G1 as a matrix A
  r=eigen(A) #spectral decomposition, special case of svd
  vectors=r$vectors #eigen vectors
  values=r$values #eigen values  #G is the square root of A, i.e, G%*%t(G)=A
  G=vectors%*%diag(sqrt(values))
  for (i in (1+nf1+nf2):(nf1+nf2+nf3)) #a total of n/k families(loops)
  {#create a series of id
    id=c( (n1+n2+1+(i-(1+nf1+nf2))*k3):(n1+n2+k3+(i-(1+nf1+nf2))*k3))
    fid=rep(i,k3)
    pid=c(1:k3)
    dadid=c(0,0,n1+n2+1+(i-(1+nf1+nf2))*k3,0,n1+n2+3+(i-(1+nf1+nf2))*k3)
    momid=c(0,0,n1+n2+2+(i-(1+nf1+nf2))*k3,0,n1+n2+4+(i-(1+nf1+nf2))*k3)
    #sex and age, sex is fixed, age is randomly selected from unif [30,40] and unif[7,11] for parents and kid respectly
    sex=c("M","F","M","F","F")
    sex1=c(1,0,1,0,0)
    age=c(runif(1,50,60), runif(1,50,60), runif(1,25,35),runif(1,25,35),runif(1,7,11))
    #random effect 
    ui=rnorm(k3)
    bi=G%*%ui
    epsi=rnorm(k3)
    family.data3=data.frame(id,fid,pid,momid,dadid,sex,sex1,age,bi,epsi)
    whole.data3=rbind(whole.data3,family.data3)
  }
  whole.data3
}
#two-trio
gen_data_twotrio<-function(n4,k4,nf4){
  df=data.frame(id=c(1:8),dadid=c(0,0,1,0,3,1,0,6),momid=c(0,0,2,0,4,2,0,7),sex=c("M","F","M","F","F","M","F","F"),famid=rep(1,8))
  pedlist = with(df,pedigree(id=id,dadid=dadid,momid=momid,sex=sex,famid=famid))
  kinmat=kinship(pedlist)
  G1=2*sig2*kinmat #calculate G, the cov matrix for each family
  A=as.matrix(G1)#this will change G1 as a matrix A
  r=eigen(A) #spectral decomposition, special case of svd
  vectors=r$vectors #eigen vectors
  values=r$values #eigen values  #G is the square root of A, i.e, G%*%t(G)=A
  G=vectors%*%diag(sqrt(values))
  for (i in (1+nf1+nf2+nf3):(nf1+nf2+nf3+nf4)) #a total of n/k families(loops)
  {#create a series of id
    id=c((1+n1+n2+n3+(i-(1+nf1+nf2+nf3))*k4):(k4+n1+n2+n3+(i-(1+nf1+nf2+nf3))*k4)) #create id by k 
    fid=rep(i,k4)
    pid=c(1:k4)
    dadid=c(0,0,1+n1+n2+n3+(i-(1+nf1+nf2+nf3))*k4,0,3+n1+n2+n3+(i-(1+nf1+nf2+nf3))*k4,1+n1+n2+n3+(i-(1+nf1+nf2+nf3))*k4,0,6+n1+n2+n3+(i-(1+nf1+nf2+nf3))*k4)
    momid=c(0,0,2+n1+n2+n3+(i-(1+nf1+nf2+nf3))*k4,0,4+n1+n2+n3+(i-(1+nf1+nf2+nf3))*k4,2+n1+n2+n3+(i-(1+nf1+nf2+nf3))*k4,0,7+n1+n2+n3+(i-(1+nf1+nf2+nf3))*k4)
    sex=c("M","F","M","F","F","M","F","F")
    sex1=c(1,0,1,0,0,1,0,0)
    age=c(runif(1,50,60), runif(1,50,60), runif(1,25,35),runif(1,25,35),
          runif(1,7,11), runif(1,25,35),runif(1,25,35),runif(1,7,11))
    ui=rnorm(k4)
    bi=G%*%ui
    epsi=rnorm(k4)
    family.data4=data.frame(id,fid,pid,momid,dadid,sex,sex1,age,bi,epsi)
    whole.data4=rbind(whole.data4,family.data4)
  }
  whole.data4
}
gen_data_threetrio<-function(n5,k5,nf5){
  df=data.frame(id=c(1:11),dadid=c(0,0,1,0,3,1,0,6,1,0,9),
                momid=c(0,0,2,0,4,2,0,7,2,0,10),
                sex=c("M","F","M","F","F","M","F","F","M","F","F"),famid=rep(1,11))
  pedlist = with(df,pedigree(id=id,dadid=dadid,momid=momid,sex=sex,famid=famid))
  kinmat=kinship(pedlist)
  G1=2*sig2*kinmat #calculate G, the cov matrix for each family
  A=as.matrix(G1)#this will change G1 as a matrix A
  r=eigen(A) #spectral decomposition, special case of svd
  vectors=r$vectors #eigen vectors
  values=r$values #eigen values  #G is the square root of A, i.e, G%*%t(G)=A
  G=vectors%*%diag(sqrt(values))
  for (i in (nf1+nf2+nf3+nf4+1):(nf1+nf2+nf3+nf4+nf5)) #a total of n/k families(loops)
  {#create a series of id
    id=c((1+n1+n2+n3+n4+(i-(1+nf1+nf2+nf3+nf4))*k5):(k5+n1+n2+n3+n4+(i-(1+nf1+nf2+nf3+nf4))*k5)) #create id by k 
    fid=rep(i,k5)
    pid=c(1:k5)
    dadid=c(0,0,1+n1+n2+n3+n4+(i-(1+nf1+nf2+nf3+nf4))*k5,0,3+n1+n2+n3+n4+(i-(1+nf1+nf2+nf3+nf4))*k5,
            1+n1+n2+n3+n4+(i-(1+nf1+nf2+nf3+nf4))*k5,0,6+n1+n2+n3+n4+(i-(1+nf1+nf2+nf3+nf4))*k5,
            1+n1+n2+n3+n4+(i-(1+nf1+nf2+nf3+nf4))*k5,0,9+n1+n2+n3+n4+(i-(1+nf1+nf2+nf3+nf4))*k5)
    momid=c(0,0,2+n1+n2+n3+n4+(i-(1+nf1+nf2+nf3+nf4))*k5,0,4+n1+n2+n3+n4+(i-(1+nf1+nf2+nf3+nf4))*k5,
            2+n1+n2+n3+n4+(i-(1+nf1+nf2+nf3+nf4))*k5,0,7+n1+n2+n3+n4+(i-(1+nf1+nf2+nf3+nf4))*k5,
            2+n1+n2+n3+n4+(i-(1+nf1+nf2+nf3+nf4))*k5,0,10+n1+n2+n3+n4+(i-(1+nf1+nf2+nf3+nf4))*k5)
    sex=c("M","F","M","F","F","M","F","F","M","F","F")
    sex1=c(1,0,1,0,0,1,0,0,1,0,0)
    age=c(runif(1,50,60),runif(1,50,60),runif(1,25,35),runif(1,25,35),runif(1,7,11),
          runif(1,25,35),runif(1,25,35),runif(1,7,11),
          runif(1,25,35),runif(1,25,35),runif(1,7,11))
    ui=rnorm(k5)
    bi=G%*%ui
    epsi=rnorm(k5)
    family.data5=data.frame(id,fid,pid,momid,dadid,sex,sex1,age,bi,epsi)
    whole.data5=rbind(whole.data5,family.data5)
  }
  whole.data5
}
gen_data<-function(n1,n2,k2,n3,k3,n4,k4,n5,k5,nf1,nf2,nf3,nf4,nf5){
  data1=gen_sing(n1)
  data2=gen_nuclear(n2,k2,nf2)
  data3=gen_data_onetrio(n3,k3,nf3)
  data4=gen_data_twotrio(n4,k4,nf4)
  data5=gen_data_threetrio(n5,k5,nf5)
  whole.data=rbind(data1,data2,data3,data4,data5)
  logit.p<-beta_0 + beta_1*whole.data$age + beta_2*whole.data$sex1 +whole.data$bi
  odds.p <- exp(logit.p)
  p=odds.p/(1+odds.p)
  whole.data$y<-rbinom(n,1,p)
  whole.data
}
whole.data=NULL
whole.data1=NULL
whole.data2=NULL
whole.data3=NULL
whole.data4=NULL
whole.data5=NULL
#define other constants
n1=200 #sample size for data1
n2=200 #sample size for data2 nuclear family
k2=4 #family size for data2 nuclear family
n3=200 #sample size for data3 one-trio family
k3=5 #family size for data3 one-trio family
n4=200 #sample size for data2 two-trio family
k4=8 #family size for data2 two-trio family
n5=220 #sample size for data2 two-trio family
k5=11 #family size for data2 three-trio family
n=n1+n2+n3+n4+n5
nf1=n1/1#number of family in data1
nf2=n2/k2
nf3=n3/k3
nf4=n4/k4
nf5=n5/k5
nF=nf1+nf2+nf3+nf4+nf5