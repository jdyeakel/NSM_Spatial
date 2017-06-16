using StatsBase
using Distributions
using RCall

include("$(homedir())/Dropbox/PostDoc/2017_NSM_Spatial/NSM_Spatial/simulation/src/meanfield_event.jl")


L=1000;
dim=2;
t_term=1000000000;
#For M=100g
lambda= 1.22149*10^-7.0;
sigma=0.0000427472;
rho=7.91067*10^-7.0;
mu=4.12274*10^-6.0;
delta=7.76415*10^-9.0;
beta=2.3468*10^-7.0;
alpha=9.45*10^-9.0;

FStar = -((alpha*lambda*mu^2*(mu + 2*rho)*(lambda - sigma))/((2*lambda*rho + mu*sigma)*(lambda*(2*delta*lambda + 2*beta*mu - lambda*mu)*rho + mu*(beta*mu + lambda*(delta + rho))* sigma)))
HStar = -((alpha*lambda^2*mu*(mu + 2*rho)*(lambda - sigma))/((2*lambda*rho + mu*sigma)*(lambda*(2*delta*lambda + 2*beta*mu - lambda*mu)*rho + mu*(beta*mu + lambda*(delta + rho))*sigma)))
RStar = (mu*(-lambda + sigma))/(2*lambda*rho + mu*sigma)

mult = 1.2;
initialdensity = [FStar*mult,HStar*mult,RStar];


time_out,prop_out,N_out = meanfield_event(L,dim,initialdensity,t_term,lambda,sigma,rho,mu,beta,delta,alpha);


R"""
plot($time_out,$(prop_out[1,:]),type='l',ylim=c(min($(prop_out[[1,2],:])),max($(prop_out[[1,2],:]))),col='green')
lines($time_out,$(prop_out[2,:]),col='orange')
lines($time_out,$(prop_out[3,:]),col='blue')
"""
