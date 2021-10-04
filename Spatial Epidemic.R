author: Akhil Thomas Panicker
library(spatstat)
library(spatgraphs)
library(igraph)
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
l=c()
#10 simulations
for (i in 1:10) {
  s=SIRspatial1(100,0.05,350,0.5,4,0.08)/100
  l=cbind(l,s)
}
file=write.csv(l,"plotprevalence1simulation")
#numerical
Num1=SIR1(100,250,0.05,0.5,4,0.08)
file=write.csv(Num1,"plotprevalence1numerical")
#prevalence over time for threshold 0.2, r =0.05
l=c()
#10 simulations
for (i in 1:10) {
  s=SIRspatial1(100,0.05,350,0.5,4,0.2)/100
  l=cbind(l,s)
}
file=write.csv(l,"plotprevalence2simulation")
#numerical
Num2=SIR1(100,250,0.05,0.5,4,0.08)
file=write.csv(Num2,"plotprevalence2numerical")
#prevalence over time for threshold 0.08, r =0.01
l=c()
#10 simulations
for (i in 1:10) {
  s=SIRspatial1(100,0.01,350,0.5,4,0.08)/100
  l=cbind(l,s)
}
file=write.csv(l,"plotprevalence3simulation")
#numerical
Num3=SIR1(100,250,0.01,0.5,4,0.08)
file=write.csv(Num3,"plotprevalence3numerical")
#prevalence over time for threshold 0.2, r =0.01
l=c()
#10 simulations
for (i in 1:10) {
  s=SIRspatial1(100,0.01,350,0.5,4,0.2)/100
  l=cbind(l,s)
}
file=write.csv(l,"plotprevalence4simulation")
#numerical
Num4=SIR1(100,250,0.01,0.5,4,0.2)
file=write.csv(Num4,"plotprevalence4numerical")
####################################################################
#code for hospitalisation of infected people
SIRspatialHosp=function(N,r1,tmax,l,f,th){
  beta=0.1
  gamma=0.05
  xa=runif(N,0,l) #initial x location of N agents between 0 and l
  ya=runif(N,0,l) #initial y location of N agents between 0 and l
  dfa=data.frame(xa,ya)
  g1=spatgraph(dfa,"geometric",par=r1) #converting to random geometric graph with radius r1
  p1=g1[]
  gr1=graph_from_adj_list(p1$edges)
  I=matrix(rep(0,N),nrow=N,ncol=1) #First infecteds
  S=matrix(rep(1,N),nrow=N,ncol=1) #First susceptibles
  R=matrix(rep(0,N),nrow=N,ncol=1) #First removed
  I1=sample(1:N, size=1)#Pick first random infected
  I[I1,1]=1
  S[I1,1]=0
  t=1
  sus=c()#list to add susceptibles
  inf=c()#list to add infected
  rec=c()#list to add recovered
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
    perI=c(perI,sum(ne))
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
HC=function(x,hmax){
  n=length(x)
  s = rep(0, n)
  # hmax = 3
  h = 1
  s[1] = x[1]
  for (i in 2:n){
    j = i - h
    #print(c(j,i, h))
    s[i] = sum(x[j:i]) # cumulative sum of at most `hmax` previuous days (based on data availability) counts
    h = min(h+1, hmax) # period : at least 1 and at most `hmax`
  }
  return(s)
}
#10 percent of this fraction is assumed to be hospitalised at each day 
#numerical calculation of perday infection to estimate hospitalised fraction
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
#assuming 15 day hospitalisation

#r = 0.05, c = 0.08
l=c()
for (i in 1:10) {
  s=SIRspatialHosp(100,0.05,300,0.5,4,0.08)/100
  l=cbind(l,(HC(s,15)*0.1))
}
file=write.csv(l,"plothosp1simulation")
#numerical
s=SIRHos(100,300,0.05,0.5,4,0.08)
Hp1=HC(s,15)*0.1
file=write.csv(Hp1,"plothosp1numerical")
#r = 0.05, c = 0.2
l=c()
for (i in 1:10) {
  s=SIRspatialHosp(100,0.05,300,0.5,4,0.2)/100
  l=cbind(l,(HC(s,15)*0.1))
}
file=write.csv(l,"plothosp2simulation")
#numerical
s=SIRHos(100,300,0.05,0.5,4,0.2)
Hp1=HC(s,15)*0.1
file=write.csv(Hp1,"plothosp2numerical")
#r = 0.05, c = 0.08
l=c()
for (i in 1:10) {
  s=SIRspatialHosp(100,0.01,300,0.5,4,0.08)/100
  l=cbind(l,(HC(s,15)*0.1))
}
file=write.csv(l,"plothosp3simulation")
#numerical
s=SIRHos(100,300,0.01,0.5,4,0.08)
Hp1=HC(s,15)*0.1
file=write.csv(Hp1,"plothosp3numerical")

#hospitalised fraction, r = 0.05, c = 0.08
l=c()
for (i in 1:10) {
  s=SIRspatialHosp(100,0.01,300,0.5,4,0.2)/100
  l=cbind(l,(HC(s,15)*0.1))
}
file=write.csv(l,"plothosp4simulation")
#numerical
s=SIRHos(100,300,0.01,0.5,4,0.2)
Hp1=HC(s,15)*0.1
file=write.csv(Hp1,"plothosp4numerical")
################################################################
###########Delay################################################
################################################################
m1 <- matrix(c(
  1, 1, 2, 2,
  1, 1, 2, 2,
  3, 3, 4, 4,
  3, 3, 4, 4), nrow = 4, ncol = 4, byrow = TRUE)
layout(m1)
#peak infection vs delay numerical solution
#implement delay "d" days for introducing and removing social distancind
SIRDelay=function(N,t,r,l,f,th,d){
  tmax=t
  S=c()
  I=c()
  R=c()
  dI=c()
  dS=c()
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
      if((inf[length(inf)]/N)>th){ #no delay
        rho=rho2
      }else{
        rho=rho1
      }
    }else{#delay
      if(h+1<(d+2)){
        rho=rho1
      }else{
        v=(inf[h-d]/100)
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
  #Peak infection for various delays, cutoff 0.10 and  0.11, r - 0.02
  l=c()
  k=c()
  for (i in 0:8) {
    l=c(l,SIRDelay(N,200,0.02,0.5,4,0.010,i))
    k=c(k,SIRDelay(N,200,0.02,0.5,4,0.011,i))
  }
  d=data.frame(l,k)
  write.csv(d,"Delayplt1")
#######################################################################
 #Peak infection for various delays, cutoff 0.05, 0.1 and 0.2, r = 0.05
  l=c()
  k=c()
  x=c()
  for (i in 0:20) {
    l=c(l,SIRDelay(N,200,0.05,0.5,4,0.1,i))
    k=c(k,SIRDelay(N,200,0.05,0.5,4,0.2,i))
    x=c(x,SIRDelay(N,200,0.05,0.5,4,0.05,i))
  }
  d=data.frame(l,k,x)
  write.csv(d,"Delayplt2")

########################################################################
  #Peak infection for various delays, cutoff 0.05, 0.1 and 0.2, r = 0.07
  l=c()
  k=c()
  x=c()
  for (i in 0:12) {
    l=c(l,SIRDelay(N,200,0.07,0.5,4,0.1,i))
    k=c(k,SIRDelay(N,200,0.07,0.5,4,0.2,i))
    x=c(x,SIRDelay(N,200,0.07,0.5,4,0.05,i))
  }
  d=data.frame(l,k,x)
  write.csv(d,"Delayplt3")
  
######################################################################
#Peak infection for various delays, cutoff 0.05, 0.1 and 0.2, r = 0.1
  l=c()
  k=c()
  x=c()
  for (i in 0:8) {
    l=c(l,SIRDelay(N,200,0.1,0.5,4,0.1,i))
    k=c(k,SIRDelay(N,200,0.1,0.5,4,0.2,i))
    x=c(x,SIRDelay(N,200,0.1,0.5,4,0.05,i))
  }
  d=data.frame(l,k,x)
  write.csv(d,"Delayplt4")
