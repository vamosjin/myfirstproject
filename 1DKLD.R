#KLD but still have some problems in KL_d function(Update needed)
############################################################
############################################################
rm(list=ls())
dev.off()
############################################################
############################################################
#(STEP1) Making 1D KL divergence function
KL_d=function(q,p){
  f=function(x, q , p) {q(x)*(q(x, T)-p(x, T))}
  t<-integrate(f,lower=-Inf,upper=Inf, q=q, p=p)
  t$value  #Forgot this $value and cause problems in optimize function
}

q=function(x,log=FALSE) dnorm(x,mean=5,sd=1,log=log)
p=function(x,log=FALSE) dnorm(x,mean=2,sd=1,log=log)
#Sanity Check 
KL_d(q,p)
KL_d(p,q)
############################################################
############################################################
#(STEP2) Changing a q function 
q_optim=function(x,qmean,log=FALSE) dnorm(x,mean=qmean,sd=1,log=log)
q_optim(1,qmean = 2)
KL_d(q_optim(5,qmean = 2),p) # Whatever the value x is, it has the same value 4.5 

############################################################
############################################################
#(STEP3) ((qmean))make a objective function which is going to be optimized 
# p~N(2,1) : The target Distribution
# q~N(qmean,1)

#For qmean                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
func=function(qmean){
  q_optim=function(x,log=FALSE){dnorm(x,mean=qmean,sd=1,log=log)
}
     KL_d(q_optim,p)
}

#Sanity Check
func(5)

#Grid Searching from Advisor (Highly Recommended)
cand<-seq(-10,10, length.out=10000)
res<-vector("numeric")
for(k in 1:length(cand)){
  res[k]<-func(qmean=cand[k])}
plot(x=cand, y=res, typ="l")
abline(v=cand[which.min(res)], col="Red",main="Grid searching")
text(x=cand[which.min(res)], y=60, labels = cand[which.min(res)]) # The answer is 1.999

func(qmean=5) #checking #4.5 with absolute error < 0.00019
optimize(func, interval = c(0,10) ) #I got the same value with the grid searching 


############################################################
############################################################
#(STEP 4) ((qmean, qsd)) make a objective function which is going to be optimized 
# p~N(2,1)
# q~N(qmean,qsd)

##For qmean, qsd
func2=function(par){
  qmean<-par[1]
  qsd<-par[2]
  q_optim=function(x,log=FALSE){dnorm(x,mean=qmean,sd=qsd,log=log)
  }
  print(par)
  return(KL_d(q_optim,p))
}
func2(par=c(2,0.02233233)) #Checking #4.5 with absolute error < 0.00019
fit=optim(par=c(50,10), fn=func2, method = "L-BFGS-B",
          lower=c(-Inf,0.000001), upper=c(Inf,Inf) ) #should check the error at some points 
fit

################################################




