#################################################################
#################################################################
#################################################################
#################################################################
#From youtube CODE FOR MFVB, CAVI (2D KL DIVERGENCE)
#################################################################
#################################################################
#################################################################
#################################################################
#https://www.youtube.com/watch?v=4MBNfT45Lps&list=PLJ71tqAZr196GJ5G36s64xifr1tURUCSJ&index=8


# #Cholesky Decomposition
# A = as.matrix(data.frame(c(3,4,3),c(4,8,6),c(3,6,9)))
# colnames(A)<-NULL
# A
# A.chol<-chol(A)
# A.chol
# t(A.chol)%*%A.chol

#################################################################
#################################################################
rm(list=ls())
dev.off()
#Minimizing KL-divergence between a correlated 2D Gaussian
#and a q1~N(0,s1), q2~N(0,s2)

library(mvtnorm) #for multivariate normal density
#install.packages("pracma")
library(pracma) #for 2d-integral, sqrtm
library(ggplot2) #for plotting 

#2D KL divergence
KL=function(q,p){
  f=function(x,y) q(x,y)*(q(x,y,T)-p(x,y,T))
  integral2(f,-20,20,-20,20)$Q
}


# KL_a=function(q,p){
#   f=function(x,y){q(x,y)*(log(q(x,y))-log(p(x,y)))}
#   integral2(f,-20,20,-20,20)$Q
# }

# ?integral2
# log(10)-log(2)
# KL(qq,post)
# KL_a(qq,post)




#covariance matrix of true posterior
Sigma=matrix(c(1,0.9,0.9,1),ncol=2)


#(True) Posterior PDF
post=function(x,y,log=F){
  dmvnorm(matrix(c(x,y),ncol=2),mean=rep(0,2),sigma=Sigma, log=F)
}





#Sanity checks
post(0,0)
integral2(post,-10,10,-10,10)$Q






#Mean Field variational family distribution (product of 1d's Gaussians)
s1=2 # variance of 1st Gaussians
s2=2 # variance of 2nd Gaussians

qq=function(x,y,lg=F){
  if(lg){
    dnorm(x,0,s1,log=lg)+dnorm(y,0,s2,log=lg)
  } else {
    dnorm(x,0,s1)*dnorm(y,0,s2)
  }
}

# qq_j=function(x,y,log=F){
#   if(log){
#     dnorm(x,0,s1,log=F)+dnorm(y,0,s2,log=F)
#   } else{
#     dnorm(x,0,s1)*dnorm(y,0,s2)
#   }
#   }
# 
# 
# 

#Sanity Check
qq(0,0)
integral2(qq,-10,10,-10,10)$Q

#KL divergence between q and posterior
KL(qq,post)
KL.divergence(qq,post)














#Wrap the q function so we can give it different values of sigma's and optimize it
func=function(t){
  qt=function(x,y,lg=F){ 
    if (lg) {
      dnorm(x,0,t[1],log=T)+dnorm(y,0,t[2],log=T)
    } else {
      dnorm(x,0,t[1])*dnorm(y,0,t[2])
    }
  }
  KL(qt,post)
}

func(c(0.4,0.4))







#plot the two distributions 
circleFun=function(center=c(0,0), diameter=1, npoints=100){
  r=diameter/2
  tt=seq(0,2*pi,length.out=npoints)
  xx=center[1]+r*cos(tt)
  yy=center[2]+r*sin(tt)
  return(data.frame(x=xx, y=yy))
}


plotCircle=function(s){
  dat=circleFun(c(0,0),4,npoints=100) #Represents the circle of 2*sigma of MVN(c(0,0),diag(1,1))
  dat1=as.matrix(dat) %*% sqrtm(Sigma)$B #Represents the posterior MVN(c(0,0), Sigma)
  dat2=as.matrix(dat) %*% diag(s) #Represents the Variational Family MVN(c(0,0),diag(s1,s2) )
  ggplot(as.data.frame(dat1),aes(x=V1,y=V2))+theme_light()+
    ylim(-5,5)+xlim(-5,5)+xlab("x")+ylab("y")+
    geom_polygon(color="blue", fill="blue", alpha=0.3)+
    geom_polygon(data=as.data.frame(dat2), aes(x=V1, y=V2, color="#ff7d9d"), 
                 fill="#ff7d9d", alpha=0.3, show.legend = F)
}


#starting value - s1=1, s2=1
plotCircle(c(1,1))

#Minimize KL(q,post)
fit=optim(par=c(1,1), fn=func, method = "L-BFGS-B",lower=0 )
fit$par










#Wait! we are cheating, we know the posterior exactly!
#so let's instead, suppose we dont know the posterior exactly, but up to some N.C

NC=0.5
post.unk=function(x,y,lg=F){
  if (lg){
    log(NC)+post(x,y,T)
  } else {
    NC*post(x,y)
  }
}

integral2(post.unk,-10,10,-10,10)$Q
post.unk(0,0)
KL(qq,post.unk)


#Wrap the q function so we can give it different values of sigma's and optimize it. 

func2=function(t){
  qt=function(x,y,lg=F){ 
    if (lg) {
      dnorm(x,0,t[1],log=T)+dnorm(y,0,t[2],log=T)
    } else {
      dnorm(x,0,t[1])*dnorm(y,0,t[2])
    }
  }
  KL(qt,post.unk) #note that since we are giving it only the joint(unnormalized posterior)
  #this is equivalent to -ELBO
}

#optim finds the minimum, so using this finds the maximum of the ELBO
fit=optim(par=c(1,1),fn=func2 ,method="L-BFGS-B", lower=0)
fit$par
plotCircle(fit$par)










#################################################################
#################################################################
#################################################################
#################################################################


#Using the CAVI algorithm on a toy example

rm(list=ls())
dev.off()

library(mvtnorm) #for multivariate normal density
library(pracma) #for 2-D integral, sqrtm
#install.packages("matlib")
library(matlib) #for solving linear equations
library(ggplot2) #for plotting 

#Prior Mu~N((0,0), I)
mu0=c(0,0)
I=diag(c(1,1))

#Likelihood x|mu~N(mu,Sigma)
Sigma=matrix(c(1,0.9,0.9,1),ncol=2)

set.seed(42)
#Assuming for simplicity that 
mu=c(0.5,-0.5)
x=rmvnorm(1,mean=mu, sigma=I) #draw 1 observation of x 

#(True) Posterior mu|X~N(mu.p, Sigma.p)
Sig.p=solve(I+solve(Sigma))
mu.p=Sig.p%*%solve(Sigma)%*%t(x)


circleFun=function(center=c(0,0),diameter=1,npoints=100){
  r=diameter/2
  tt=seq(0,2*pi, length.out=npoints )
  xx=center[1]+r*cos(tt)
  yy=center[2]+r*sin(tt)
  return(data.frame(x=xx,y=yy))
}




#solve the CAVI Linear equations
A=matrix(c(1,-4/73/6.26,-4/73/6.26,1),ncol=2)
x_=c(x[2],x[1])
b=t(x_*(5.26/6.26)-x_*(4.73/6.26))
Solve(A, b) #######This is not working porperly!!!!!!!!!!! 
mu.cavi=c(1.47, -1.19)


#The derived CAVI variance 
s=sqrt(1/6.26)


#Plotting the Posterior and the CAVI solution 
dat=circleFun(c(0,0),4,npoints=100) #Represents the cicle of 2*sigma of MVN(c(0,0),diag(1,1))
dat1=as.matrix(dat)%*%sqrtm(Sig.p)$B #Represents the posterior MVN(c(0,0),Sig.p)
dat1=t(t(dat1)+as.vector(mu.p))    #Represents the posterior MVN(mu.p,Sig.p)
dat2=as.matrix(dat)%*%diag(c(s,s)) #Represents the CAVI MVN(c(0,0),diag(s,s))
dat2=t(t(dat2)+as.vector(mu.cavi))  #Represents the CAVI MVN(mu.cavi,daig(s,s))



ggplot(as.data.frame(dat1),aes(x=V1, y=V2))+theme_light()+
  ylim(-3,3)+xlim(-3,3)+xlab("mu1")+ylab("mu2")+
  geom_polygon(color="#2E9FDF",fill="#2E9FDF",alpha=0.3)+
  geom_point(data=as.data.frame(t(mu.p)),aes=(x=V1, y=V2),color="#2E9FDF")+
  geom_polygon(data=as.data.frame(dat2), aes(x=V1,y=V2,color="#ff7d9d"),
               fill="#ff7d9d", alpha=0.3, show.legend=F)+
  geom_point(data=as.data.frame(t(mu.cavi)), aes(x=V1, y=V2), color="#ff7d9d")




ggplot(as.data.frame(dat1))+theme_light()+
  ylim(-3,3)+xlim(-3,3)+xlab("mu1")+ylab("mu2")+
  geom_polygon(color="#2E9FDF",fill="#2E9FDF",alpha=0.3)+
  geom_point(data=as.data.frame(t(mu.p)),color="#2E9FDF")+
  geom_polygon(data=as.data.frame(dat2), fill="#ff7d9d", alpha=0.3, show.legend=F)+
  geom_point(data=as.data.frame(t(mu.cavi)), color="#ff7d9d")



#n data points
n=10
x=rmvnorm(10,mean=mu, sigma=I)

#(True) Posterior mu|x~ N(mu.p, Sigma.p)
Sig.p=solve(I+n*solve(Sigma))
mu.p=n*Sig.p%*%solve(Sigma)%*%colMeans(x)

#Solve the CAVI linear equations
A=matrix(c(1,-4.73*n/(1+n*5.26),-4.73*n/(1+n*5.26),1),ncol=2)
x1=colSums(x)
x2=c(x1[2], x1[1])
b=t(t(x1*(5.26/(1+n*5.26)) - x2*(4.73/(1+n*5.26))))
Solve(A,b)
mu.cavi=c(0.591, -0.512)

#the derived CAVI variance
s=sqrt(1/(1+n*5.26))

#Plotting the Posterior and the CAVI solution 
dat=circleFun(c(0,0),4,npoints=100) #Represents the cicle of 2*sigma of MVN(c(0,0),diag(1,1))
dat1=as.matrix(dat)%*%sqrtm(Sig.p)$B #Represents the posterior MVN(c(0,0),Sig.p)
dat1=t(t(dat1)+as.vector(mu.p))    #Represents the posterior MVN(mu.p,Sig.p)
dat2=as.matrix(dat)%*%diag(c(s,s)) #Represents the CAVI MVN(c(0,0),diag(s,s))
dat2=t(t(dat2)+as.vector(mu.cavi))  #Represents the CAVI MVN(mu.cavi,daig(s,s))


ggplot(as.data.frame(dat1),aes(x=V1, y=V2))+theme_light()+
  ylim(-3,3)+xlim(-3,3)+xlab("mu1")+ylab("mu2")+
  geom_polygon(color="#2E9FDF",fill="#2E9FDF",alpha=0.3)+
  geom_point(data=as.data.frame(t(mu.p)),aes=(x=V1, y=V2),color="#2E9FDF")+
  geom_polygon(data=as.data.frame(dat2), aes(x=V1,y=V2,color="#ff7d9d"),
               fill="#ff7d9d", alpha=0.3, show.legend=F)+
  geom_point(data=as.data.frame(t(mu.cavi)), aes(x=V1, y=V2), color="#ff7d9d")


#D.Raaeli
