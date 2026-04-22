OpenBUGS code
#ns: number of studies
#nc: number of treatments
#nt: number of doses
#t: dose code per study arm
#na: number of doses per study
#nat: number of treatments per study arm
#r: number of events per study arm
#n: sample size per study arm
#ref: reference treatment (codes assume ref=1 is placebo [with a single dose] and require
modification if this is not the case)
#tau.sq=between-study variance within dose-level
#sigma.sq=between-dose variance within treatment-level
##### The following part of the code is common in all models #####
model{
  for(i in 1:ns) {
    w[i,1] <-0
    delta[i,t[i,1]]<-0
    ##prior distribution for log-odds in baseline arm of study i
    u[i] ~ dnorm(0,.001)
    ##binomial likelihood of number of events for each arm k of study i
    for (k in 1:na[i]) {
      r[i,t[i,k]] ~ dbin(p[i,t[i,k]],n[i,t[i,k]])
    }
    ##parameterization of the 'true' effect of each comparison
    ##of arm k vs. baseline arm (1) of study i
    logit(p[i,t[i,1]])<-u[i]
    for (k in 2:na[i]) {
      logit(p[i,t[i,k]])<-u[i] + delta[i,t[i,k]]
      ##distribution of random effects
      delta[i,t[i,k]] ~ dnorm(md[i,t[i,k]],precd[i,t[i,k]])
      precd[i,t[i,k]] <- prec *2*(k-1)/k
      ##assumption of consistency on the dose level (mu in File 1)
      md[i,t[i,k]]<- d[t[i,k]] - d[t[i,1]] + sw[i,k]
      w[i,k]<- delta[i,t[i,k]] - d[t[i,k]] + d[t[i,1]]
      sw[i,k]<- sum(w[i,1:k-1])/(k-1)
    }
  }
  ##prior for between-study variance
  #tau~dnorm(0,1)I(0,) #minimally informative prior
  #tau.sq<-pow(tau,2)
  tau.sq~dlnorm(-3.02,0.29) #stroke, nausea, headache (informative prior)
  prec<-1/tau.sq
  var <-1/prec
  ##### For equal dose effects use the following part of the code #####
  ##dose effect and treatment effect are zero for reference treatment (here we assume ref=1
  [placebo])
d[ref]<-0 # estimated LOR for each dose-comparison vs. reference treatment-dose (mu
in File 1)
D[ref]<-0 # estimated LOR for each treatment-comparison vs. reference treatment
(lambda in File 1)
for(k in 1:(ref-1)){
  for(j in class[k]:(class[k+1]-1)){
    d[j]<- D[k]
  }
}
for(k in (ref+1):nc){
  for(j in class[k]:(class[k+1]-1)){
    d[j]<- D[k]
  }
}
## Vague priors for treatment effects (basic parameters)
for(k in 1:(ref-1)){
  D[k] ~ dnorm(0,.0001)
}
for(k in (ref+1):nc){
  D[k] ~ dnorm(0,.0001)
}
##estimated & predictive OR for each treatment-comparison
for(i in 1:(nc-1)) {
  for (j in (i+1):nc) {
    OR[j,i]<- exp(D[j] - D[i])
    LOR[j,i]<- D[j] - D[i]
    predLOR[j,i] ~ dnorm(LOR[j,i],prec)
    predOR[j,i]<- exp(predLOR[j,i])
  }
}
##treatment ranking
for(k in 1:nc) {
  order[k]<- rank(D[],k) #harmful outcome <- t+1-rank (d[],k) when beneficial
  outcome
  most.effective[k]<-equals(order[k],1)
  for(j in 1:nc) {effectiveness[k,j]<- equals(order[k],j)
  cumeffectiveness[k,j]<- sum(effectiveness[k,1:j])}}
##SUCRA
for(k in 1:nc) {
  SUCRA[k]<- sum(cumeffectiveness[k,1:(nc-1)]) /(nc-1)
}
###########################################################
##### For separate dose effects use the following part of the code #####
##Independent dose-effects (not related with their parent treatments) (mu in File 1)
d[ref]<-0
for (k in 1:(ref-1)){
  d[k] ~ dnorm(0,.0001)
}
for (k in (ref+1):nt){
  d[k] ~ dnorm(0,.0001)
}
##estimated & predicted OR for each treatment-comparison
for(i in 1:(nc-1)) {
  for (j in (i+1):nc) {
    OR[j,i]<- exp(mean(d[class[j]:(class[j+1]-1)]) -
                    mean(d[class[i]:(class[i+1]-1)]))
    LOR[j,i]<- mean(d[class[j]:(class[j+1]-1)]) - mean(d[class[i]:(class[i+1]-
                                                                     1)])
    predLOR[j,i] ~ dnorm(LOR[j,i],prec)
    predOR[j,i]<- exp(predLOR[j,i])
  }
}
##estimated & predicted OR for each dose-comparison
for(i in 1:(nt-1)) {
  for (j in (i+1):nt) {
    OR.dose[j,i]<- exp(d[j] - d[i])
    LOR.dose[j,i]<- d[j] - d[i]
    predLOR.dose[j,i] ~ dnorm(LOR.dose[j,i],prec)
    predOR.dose[j,i]<- exp(predLOR.dose[j,i])
  }
}
#dose ranking
for(k in 1:nt) {
  order[k]<- rank(d[],k) #harmful outcome <- t+1-rank (d[],k) when beneficial
  outcome
  most.effective[k]<-equals(order[k],1)
  for(j in 1:nt) {effectiveness[k,j]<- equals(order[k],j)
  cumeffectiveness[k,j]<- sum(effectiveness[k,1:j])}}
##SUCRA
for(k in 1:nt) {
  SUCRA[k]<- sum(cumeffectiveness[k,1:(nt-1)]) /(nt-1)
}
###########################################################
##### For exchangeable dose effects use the following part of the code #####
##dose effect and treatment effect are zero for reference treatment (here we assume ref=1
[placebo])
d[ref]<-0 # estimated LOR for each dose-comparison vs. reference treatment-dose (mu in
File 1)
D[ref]<-0 # estimated LOR for each treatment-comparison vs. reference treatment
(lambda in File 1)
for(k in 1:(ref-1)){
  for(j in class[k]:(class[k+1]-1)){
    d[j]<- D[k]+var[j]
  }
}
for(k in (ref+1):nc){
  for(j in class[k]:(class[k+1]-1)){
    d[j]<- D[k]+var[j]
  }
}
##no between-dose variance for the reference treatment
var[ref]<-0
for(k in 2:nt){
  var[k] ~ dnorm(0,sigma.prec)
}
##between-dose variance
sigma~dnorm(0,1)I(0,)
sigma.sq<-pow(sigma,2)
sigma.prec<-1/sigma.sq
sigma.var<-1/sigma.prec
##vague priors for treatment effects (basic parameters)
for(k in 1:(ref-1)){
  D[k] ~ dnorm(0,.0001)
}
for(k in (ref+1):nc){
  D[k] ~ dnorm(0,.0001)
}
##estimated & predictive OR for each treatment-comparison
for(i in 1:(nc-1)) {
  for (j in (i+1):nc) {
    OR[j,i]<- exp(D[j] - D[i])
    LOR[j,i]<- D[j] - D[i]
    predLOR[j,i] ~ dnorm(LOR[j,i],prec)
    predOR[j,i]<- exp(predLOR[j,i])
  }
}
##estimated & predicted OR for each dose-comparison
for(i in 1:(nt-1)) {
  for (j in (i+1):nt) {
    OR.dose[j,i]<- exp(d[j] - d[i])
    LOR.dose[j,i]<- d[j] - d[i]
    predLOR.dose[j,i] ~ dnorm(LOR.dose[j,i],prec)
    predOR.dose[j,i]<- exp(predLOR.dose[j,i])
  }
}
#Ranking of treatment-doses#
for(k in 1:nt) {
  order[k]<- rank(d[],k) #harmful outcome <- t+1-rank (d[],k) when beneficial
  outcome
  most.effective[k]<-equals(order[k],1)
  for(j in 1:nt) {effectiveness[k,j]<- equals(order[k],j)
  cumeffectiveness[k,j]<- sum(effectiveness[k,1:j])}}
##SUCRA
for(k in 1:nt) {
  SUCRA[k]<- sum(cumeffectiveness[k,1:(nt-1)]) /(nt-1)
  ###########################################################
  ##### The following part of the code is common in all models #####
  ##model fit
  for(i in 1:ns) {
    for (k in 1:na[i]) {
      Darm[i,k]<- -2*( r[i,t[i,k]] *log(n[i,t[i,k]]*p[i,t[i,k]]/ r[i,t[i,k]])+(n[i,t[i,k]]
                                                                               - r[i,t[i,k]])*log((n[i,t[i,k]]-n[i,t[i,k]]* p[i,t[i,k]])/(n[i,t[i,k]]- r[i,t[i,k]])))
    }
    Dsumarm[i]<- sum(Darm[i,1:na[i]])
  }
  D.bar<- sum(Dsumarm[])
}