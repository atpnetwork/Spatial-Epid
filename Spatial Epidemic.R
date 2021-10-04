# author: Akhil_Thomas_Panicker
library(spatstat)
library(spatgraphs)
library(igraph)
#3x3 layout
l=c(
  1, 1, 2, 2, 3, 3, 
  1, 1, 2, 2, 3, 3, 
  4, 4, 5, 5, 6, 6, 
  4, 4, 5, 5, 6, 6,
  7, 7, 8, 8, 9, 9,
  7, 7, 8, 8, 9, 9
)
m1 <- matrix(l, nrow = 6, ncol = 6, byrow = TRUE)
layout(m1)

#################################################
#prevalence over time
#################################################
#N = number of individuals
#r1 = the characteristic connection distance
#tmax = Maximum time of the simulation
#l = length of the square patch
#f factor by which Social Distancing is done
#th is the threshold "c" that decides the behaviour conditions 
SIRspatial1=function(N,r1,tmax,l,f,th){
  beta=0.1
  gamma=0.05
  xa=runif(N,0,l) #initial x location of N agents between 0 and l
  ya=runif(N,0,l) #initial y location of N agents between 0 and l
  dfa=data.frame(xa,ya)
  g1=spatgraph(dfa,"geometric",par=r1) #converting to random geometric graph with radius r1
  p1=g1[]
  gr1=graph_from_adj_list(p1$edges)
  I=matrix(rep(0,N),nrow=N,ncol=1) 
  S=matrix(rep(1,N),nrow=N,ncol=1) 
  R=matrix(rep(0,N),nrow=N,ncol=1) 
  I1=sample(1:N, size=1)
  I[I1,1]=1
  S[I1,1]=0
  t=1
  sus=c()#list to add susceptibles
  inf=c()#list to add infected
  rec=c()#list to add recovered
  gr=gr1
  while (t<=tmax) {
    t=t+1
    infneigh=gr[]%*%I[,t-1]
    #-------------------------------------------------------------------------
    # this part of code adapted from Epidemics: Models and Data using R by Ottar N Bjornstad
    pinf=1-(1-beta)**infneigh 
    newI=c() 
    for (i in 1:N) {
      new=rbinom(1,S[,t-1][i],pinf[i])
      newI=c(newI,new)
    }
    newR=rbinom(N,I[,t-1],gamma)
    #-------------------------------------------------------------------------
    nextS=S[,t-1]-newI 
    nextI=I[,t-1]+newI-newR 
    nextR=R[,t-1]+newR 
    inf=c(inf,sum(nextI)) 
    sus=c(sus,sum(nextS)) 
    rec=c(rec,sum(nextR)) 
    I=cbind(I, nextI)
    S=cbind(S, nextS)
    R=cbind(R, nextR)
    if(sum((nextI)/N)<th){ #if prevalence fraction less than threshold
      xa=runif(N,0,l) #no social distancing
      ya=runif(N,0,l) 
      dfa=data.frame(xa,ya)
      g1=spatgraph(dfa,"geometric",par=r1)
      p1=g1[]
      gr1=graph_from_adj_list(p1$edges)
      gr=gr1
    }else{ #if prevalence fraction more than threshold
      xb=runif(N,0,l*f) #follow social distancing, 
      yb=runif(N,0,l*f) 
      dfb=data.frame(xb,yb)
      g2=spatgraph(dfb,"geometric",par=r1)
      p2=g2[]
      gr2=graph_from_adj_list(p2$edges)
      gr=gr2
    }
  }
  infe=(c(1,inf)[1:tmax])#list of infected people over tmax time
  return(infe)
}
#####################################################
#numerical soultion of the difference equation
SIR1=function(N,tmax,r,l,f,th){
  S=c()
  I=c()
  R=c()
  dS=c()
  dI=c()
  t <- seq(0,tmax,1)
  S[1] <- N-1
  I[1] <- 1
  R[1] <- 0
  rho1=N/(l**2) #density without social distancing
  rho2=N/((l*f)**2)#density with social distancing
  inf=c()
  b=0.1
  g=0.05
  rho=rho1
  for (h in seq_len(length(t)-1)) {
    dS[h+1] <- S[h]*(1-(1-b)**(pi*r*r*rho*I[h]/N))
    S[h+1]  <- S[h]-dS[h+1]
    I[h+1]  <- I[h] + dS[h+1] - g*I[h]
    R[h+1]  <- R[h]+g*I[h]
    if((I[h+1]/N)>th){ #if prevalence fraction greater than thresold
      rho=rho2# density during social distancing 
    }else{    #if prevalence fraction less than threshold
      rho=rho1#density when not following social distancing 
    }
    if(I[h+1]<0){
      break
    }
    inf=c(inf,I[h+1])
  }
  return(c(1,inf)[1:tmax])
}
layout(m1)
################################################################
#prevalence over time for threshold 0.08, r =0.05
plot(c(0,350),c(0,0.15),type="n",xlab = "Time",ylab = "Prevalence",
     main = "c = 0.08, r = 0.05",cex.main = 1.8,cex.axis=1.1,cex.lab=1.3)
#10 simulations
for (i in 1:10) {
  s=SIRspatial1(100,0.05,350,0.5,4,0.08)/100
  lines(1:350,s,col="green",lwd=0.35)
}
lines(1:350,SIR1(100,350,0.05,0.5,4,0.08)/100,col="blue",lwd=1.1)
abline(h=0.08,col="brown",lty=2)
text(300, 0.1, expression(c == 0.08),cex=1.4)
#prevalence over time for threshold 0.2, r =0.05
plot(c(0,250),c(0,0.35),type="n",xlab = "Time",ylab = "Prevalence",
     main = "c = 0.2, r = 0.05",cex.main = 1.8,cex.axis=1.1,cex.lab=1.3)
#10 simulations
for (i in 1:10) {
  s=SIRspatial1(100,0.05,250,0.5,4,0.2)/100
  lines(1:250,s,col="green",lwd=0.35)
}
lines(1:250,SIR1(100,250,0.05,0.5,4,0.2)/100,col="blue",lwd=1)
abline(h=0.2,col="brown",lty=2)
text(200, 0.25, expression(c == 0.2),cex=1.4)
plot(c(0,200),c(0,1),axes=FALSE,xlab="",ylab = "",type="n",
     frame.plot=TRUE)
legend(10, 0.75, legend=c("Mean field Solution", "Individual Runs"),
       col=c("blue","green"), lty=1, cex=1.4)
#prevalence over time for threshold 0.08, r =0.01
plot(c(0,100),c(0,0.04),type="n",xlab = "Time",ylab = "Prevalence",
     main = "c = 0.08, r = 0.01",cex.main = 1.8,cex.axis=1.1,cex.lab=1.3)
#10 simulations
for (i in 1:10) {
  s=SIRspatial1(100,0.01,100,0.5,4,0.08)/100
  lines(1:100,s,col="green",lwd=0.35)
}
lines(1:100,SIR1(100,100,0.01,0.5,4,0.08)/100,col="blue",lwd=1)
#prevalence over time for threshold 0.02, r =0.01
plot(c(0,100),c(0,0.05),type="n",xlab = "Time",ylab = "Prevalence",
     main = "c = 0.02, r = 0.01",cex.main = 1.8,cex.axis=1.1,cex.lab=1.3)
for (i in 1:10) {
  s=SIRspatial1(100,0.01,100,0.5,4,0.02)/100
  lines(1:100,s,col="green",lwd=0.35)
}
lines(1:100,SIR1(100,100,0.01,0.5,4,0.02)/100,col="blue",lwd=1)
####################################################################
#code for hospitalisation of infected people
SIRspatialHosp=function(N,r1,tmax,l,f,th){
  beta=0.1
  gamma=0.05
  xa=runif(N,0,l) 
  ya=runif(N,0,l) 
  dfa=data.frame(xa,ya)
  g1=spatgraph(dfa,"geometric",par=r1) 
  p1=g1[]
  gr1=graph_from_adj_list(p1$edges)
  I=matrix(rep(0,N),nrow=N,ncol=1) 
  S=matrix(rep(1,N),nrow=N,ncol=1) 
  R=matrix(rep(0,N),nrow=N,ncol=1) 
  I1=sample(1:N, size=1)
  I[I1,1]=1
  S[I1,1]=0
  t=1
  sus=c()
  inf=c()
  rec=c()
  gr=gr1
  perI=c()
  while (t<=tmax) {
    t=t+1
    infneigh=gr[]%*%I[,t-1]
    pinf=1-(1-beta)**infneigh #effective probability for each agent 
    newI=c() #new infected people
    for (i in 1:N) {
      new=rbinom(1,S[,t-1][i],pinf[i])
      newI=c(newI,new)
    }
    newR=rbinom(N,I[,t-1],gamma)
    ne=newI-I[,t-1]
    ne=replace(ne,ne<0,0)
    perI=c(perI,sum(ne))#per day infections
    nextS=S[,t-1]-newI #updated S
    nextI=I[,t-1]+newI-newR #updated I
    nextR=R[,t-1]+newR #updated R
    inf=c(inf,sum(nextI)) #appending S
    sus=c(sus,sum(nextS)) #appending I
    rec=c(rec,sum(nextR)) #appending R
    I=cbind(I, nextI)
    S=cbind(S, nextS)
    R=cbind(R, nextR)
    if(sum((nextI)/N)<th){ #if prevalence fraction less than threshold
      xa=runif(N,0,l) #no social distancing
      ya=runif(N,0,l) 
      dfa=data.frame(xa,ya)
      g1=spatgraph(dfa,"geometric",par=r1)
      p1=g1[]
      gr1=graph_from_adj_list(p1$edges)
      gr=gr1
    }else{ #if prevalence fraction more than threshold
      xb=runif(N,0,l*f) #follow social distancing, 
      yb=runif(N,0,l*f) #x and y location of N agents between 0 and l*f
      dfb=data.frame(xb,yb)
      g2=spatgraph(dfb,"geometric",par=r1)
      p2=g2[]
      gr2=graph_from_adj_list(p2$edges)
      gr=gr2
    }
  }
  return(perI)
}
########################################################
#Hospitalised case 
#we cumulatively add the perday infection for hmax days so that a fraction
#of which can be assumed to ne hospitalised
HC=function(x,hmax){
  n=length(x)
  s = rep(0, n)
  # hmax = 3
  h = 1
  s[1] = x[1]
  for (i in 2:n){
    j = i - h
    #print(c(j,i, h))
    s[i] = sum(x[j:i])
    h = min(h+1, hmax)
  }
  return(s)
}
#numerical  solution for calculating hoslpitalised fraction
SIRHos=function(N,t,r,l,f,th){
  tmax=t
  t <- seq(0,t,1)
  S[1] <- N-1
  I[1] <- 1
  R[1] <- 0
  rho1=N/(l**2)
  rho2=N/((l*f)**2)
  su=c()
  inf=c()
  re=c()
  b=0.1
  g=0.05
  rho=rho1
  for (h in seq_len(length(t)-1)) {
    dS[h+1] <- S[h]*(1-(1-b)**(pi*r*r*rho*I[h]/N))
    S[h+1]  <- S[h]-dS[h+1]
    I[h+1]  <- I[h] + dS[h+1] - I[h]*g
    R[h+1]  <- R[h]+I[h]*g
    if((I[h+1]/N)>th){
      rho=rho2
    }else{
      rho=rho1
    }
    if(I[h+1]<0){
      break
    }
    su=c(su,dS[h+1])
    inf=c(inf,I[h+1])
    re=c(re,R[h+1])
  }
  return(su)
}
#10 percentage of the per day cumulative upto hmax =15 hospitalisation is assumed to 
#be hospitalised
#r=0.05, c = 0.08
plot(c(0,300),c(0,0.025),col="red",type="n",
     xlab = "Time",ylab = "Hospitalised Fraction",
     main="c = 0.08, r = 0.05",cex.main=2,cex.lab=1.2)
for (i in 1:10) {
  s=SIRspatialHosp(100,0.05,300,0.5,4,0.08)/100
  lines(1:300,HC(s,15)*0.1,col="green",lwd=0.35)
}
s=SIRHos(100,300,0.05,0.5,4,0.08)/100
lines(1:300,HC(s,15)*0.1,col="blue",lwd=1.3)
#r=0.05, c = 0.2
plot(c(0,200),c(0,0.04),col="red",type="n",
     xlab = "Time",ylab = "Hospitalised Fraction",
     main="c = 0.2, r = 0.05",cex.main=2,cex.lab=1.2)
for (i in 1:10) {
  s=SIRspatialHosp(100,0.05,200,0.5,4,0.2)/100
  lines(1:200,HC(s,15)*0.1,col="green",lwd=0.35)
}
s=SIRHos(100,200,0.05,0.5,4,0.2)/100
lines(1:200,HC(s,15)*0.1,col="blue",lwd=1.3)
#r=0.01, c = 0.08
plot(c(0,100),c(0,0.004),col="red",type="n",
     xlab = "Time",ylab = "Hospitalised Fraction",
     main="c = 0.08, r = 0.01",cex.main=2,cex.lab=1.2)
for (i in 1:10) {
  s=SIRspatialHosp(100,0.01,100,0.5,4,0.08)/100
  lines(1:100,HC(s,15)*0.1,col="green",lwd=1.2)
}
s=SIRHos(100,100,0.01,0.5,4,0.08)/100
lines(1:100,HC(s,15)*0.1,col="blue",lwd=1.3)
#r=0.01, c = 0.2
plot(c(0,100),c(0,0.004),col="red",type="n",
     xlab = "Time",ylab = "Hospitalised Fraction",
     main="c = 0.2, r = 0.01",cex.main=2,cex.lab=1.2)
for (i in 1:10) {
  s=SIRspatialHosp(100,0.01,100,0.5,4,0.2)/100
  lines(1:100,HC(s,15)*0.1,col="green",lwd=1.35)
}
s=SIRHos(100,100,0.01,0.5,4,0.2)/100
lines(1:100,HC(s,15)*0.1,col="blue",lwd=1.1)
################################################################
###########Delay################################################
################################################################
m1 <- matrix(c(
  1, 1, 2, 2, 
  1, 1, 2, 2, 
  3, 3, 4, 4, 
  3, 3, 4, 4), nrow = 4, ncol = 4, byrow = TRUE)
layout(m1)
#peak infection vs delay numerical
#whenever the prevalence exceeds the threshold, 
#social distancing is implemented after d days
#whenever it comes below threshols, SD is removed after d days
SIRDelay=function(N,t,r,l,f,th,d){
  tmax=t
  t <- seq(0,t,1)
  S[1] <- N-1
  I[1] <- 1
  R[1] <- 0
  rho1=N/(l**2)
  rho2=N/((l*f)**2)
  su=c()
  inf=c()
  re=c() 
  b=0.1
  g=0.05
  rho=rho1
  for (h in seq_len(length(t)-1)) {
    dS[h+1] <- S[h]*(1-(1-b)**(pi*r*r*rho*I[h]/N))
    S[h+1]  <- S[h]-dS[h+1]
    I[h+1]  <- I[h] + dS[h+1] - I[h]*g
    R[h+1]  <- R[h]+I[h]*g
    su=c(su,S[h+1])
    inf=c(inf,I[h+1])
    re=c(re,R[h+1])
    if(d==0){
      if((inf[length(inf)]/N)>th){
        rho=rho2
      }else{
        rho=rho1
      }
    }else{
      if(h+1<(d+1)){
        rho=rho1
      }else{
        v=(inf[h+1-d]/100)
        if(v>th){
          rho=rho2
        }else{
          rho=rho1
        }
      }
    }
  }
  return(max(inf))
}
######################################################################
#radius0.02
plottwoa=function(N){
  l=c()
  k=c()
  for (i in 0:8) {
    l=c(l,SIRDelay(N,200,0.02,0.5,4,0.010,i))
    k=c(k,SIRDelay(N,200,0.02,0.5,4,0.011,i))
  }
  plot(c(0,8),c(0.009,0.015),type = "n",xlab = "Delay (Days)",ylab = "Peak Infection",
       cex.lab=1.5, cex.axis=1.5,
       cex.main=2, main = "r = 0.02")
  lines(0:8,l/100,col="blue",lwd=1.5)
  points(0:8,l/100,col="blue",lwd=1.5,pch=16)
  lines(0:8,k/100,col="green",lwd=1.5)
  points(0:8,k/100,col="green",lwd=1.5,pch=16)
  legend(5, 0.0145, legend=c("0.010","0.011"),
         col=c("blue","green"),
         lty=1,cex=1.3,title = "c")
}
plottwoa(100)
##############################################################
#radius0.05
plottwob=function(N){
  l=c()
  k=c()
  x=c()
  for (i in 0:20) {
    l=c(l,SIRDelay(N,200,0.05,0.5,4,0.1,i))
    k=c(k,SIRDelay(N,200,0.05,0.5,4,0.2,i))
    x=c(x,SIRDelay(N,200,0.05,0.5,4,0.05,i))
  }
  plot(c(0,20),c(0,0.65),type = "n",xlab = "Delay (Days)",ylab = "Peak Infection",
       cex.lab=1.5, cex.axis=1.5,
       cex.main=2, main = "r = 0.05")
  lines(0:20,l/100,col="blue",lwd=1.5)
  points(0:20,l/100,col="blue",lwd=1.5,pch=16)
  lines(0:20,k/100,col="orange",lwd=1.5)
  points(0:20,k/100,col="orange",lwd=1.5,pch=16)
  lines(0:20,x/100,col="green",lwd=1.5)
  points(0:20,x/100,col="green",lwd=1.5,pch=16)
  legend(15, 0.2, legend=c("0.2","0.1","0.05"),
         col=c("orange","blue","green"),
         lty=1,cex=1.2,title = "c")
}
plottwob(100)

######################################################################
#radius0.07
plottwoc=function(N){
  l=c()
  k=c()
  x=c()
  for (i in 0:12) {
    l=c(l,SIRDelay(N,200,0.07,0.5,4,0.1,i))
    k=c(k,SIRDelay(N,200,0.07,0.5,4,0.2,i))
    x=c(x,SIRDelay(N,200,0.07,0.5,4,0.05,i))
  }
  plot(c(0,12),c(0,0.8),type = "n",xlab = "Delay (Days)",ylab = "Peak Infection",
       cex.lab=1.5, cex.axis=1.5,
       cex.main=2, main = "r = 0.07")
  lines(0:12,l/100,col="blue",lwd=1.5)
  points(0:12,l/100,col="blue",lwd=1.5,pch=16)
  lines(0:12,k/100,col="orange",lwd=1.5)
  points(0:12,k/100,col="orange",lwd=1.5,pch=16)
  lines(0:12,x/100,col="green",lwd=1.5)
  points(0:12,x/100,col="green",lwd=1.5,pch=16)
  legend(8, 0.4, legend=c("0.2","0.1","0.05"),
         col=c("orange","blue","green"),
         lty=1,cex=1.2,title = "c")
}
plottwoc(100)
######################################################################
#radius 0.1
plottwod=function(N){
  l=c()
  k=c()
  x=c()
  for (i in 0:8) {
    l=c(l,SIRDelay(N,200,0.1,0.5,4,0.1,i))
    k=c(k,SIRDelay(N,200,0.1,0.5,4,0.2,i))
    x=c(x,SIRDelay(N,200,0.1,0.5,4,0.05,i))
  }
  plot(c(0,8),c(0,0.9),type = "n",xlab = "Delay (Days)",ylab = "Peak Infection",
       cex.lab=1.5, cex.axis=1.5,
       cex.main=2, main = "r = 0.1")
  lines(0:8,l/100,col="blue",lwd=1.5)
  points(0:8,l/100,col="blue",lwd=1.5,pch=16)
  lines(0:8,k/100,col="orange",lwd=1.5)
  points(0:8,k/100,col="orange",lwd=1.5,pch=16)
  lines(0:8,x/100,col="green",lwd=1.5)
  points(0:8,x/100,col="green",lwd=1.5,pch=16)
  legend(5, 0.4, legend=c("0.2","0.1","0.05"),
         col=c("orange","blue","green"),
         lty=1,cex=1.2,title = "c")
  
}
#######################################################################
plottwod(100)

