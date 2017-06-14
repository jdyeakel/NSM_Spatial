
using StatsBase
using Gadfly
using Cairo


#include("/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/src/starvingforager_event.jl")
include("$(homedir())/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/src/starvingforager_event_rate.jl")
# include("$(homedir())/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/src/starvingforager_event_spatial.jl")




sigmavec = collect(0.4:0.1:1.0);
Fstar = zeros(length(sigmavec));
Hstar = zeros(length(sigmavec));
Rstar = zeros(length(sigmavec));
Fsimstar = zeros(length(sigmavec));
Hsimstar = zeros(length(sigmavec));
Rsimstar = zeros(length(sigmavec));
for i = 1:length(sigmavec)
  println("sigma = ",sigmavec[i])
  L = 50;
  dim = 2;
  prop_fill = 0.5
  initsize = convert(Int64,round(((L-2)^dim)*prop_fill));
  t_term = 100;
  alpha = 0.5;
  K = 1;
  sigma = sigmavec[i];
  rho = 0.2;
  m = 0.8;
  lambda = 0.2;
  mu = 0.1;

  Fstar[i] = (alpha*lambda*mu*(mu + rho))/((lambda*rho + m*mu)*(lambda*rho + sigma*mu));
  Hstar[i] = (alpha*lambda^2*(mu + rho))/((lambda*rho + m*mu)*(lambda*rho + sigma*mu));
  Rstar[i] = (mu*(-lambda+sigma))/(lambda*rho + mu*sigma);

  #The simulation
  time_out, prop_out, N_out = starvingforager_event_rate(L,dim,initsize,t_term,alpha,K,sigma,rho,m,lambda,mu);
  F = prop_out[1,:];
  H = prop_out[2,:];
  R = prop_out[3,:];
  last = length(F);
  burnin = Int(round(last*0.75));
  Fsimstar[i] = mean(F[burnin:last]);
  Hsimstar[i] = mean(H[burnin:last]);
  Rsimstar[i] = mean(R[burnin:last]);

end

comparison = plot(
layer(x=sigmavec,y=Fstar,Geom.line,Theme(default_color=colorant"green")),
layer(x=sigmavec,y=Fsimstar,Geom.point,Theme(default_color=colorant"green",default_point_size=3pt)),
layer(x=sigmavec,y=Hstar,Geom.line,Theme(default_color=colorant"orange")),
layer(x=sigmavec,y=Hsimstar,Geom.point,Theme(default_color=colorant"orange",default_point_size=3pt)),
layer(x=sigmavec,y=Rstar,Geom.line,Theme(default_color=colorant"blue")),
layer(x=sigmavec,y=Rsimstar,Geom.point,Theme(default_color=colorant"blue",default_point_size=3pt)),
Guide.XLabel("sigma"),Guide.YLabel("Steady state; F=green; H=orange; R=blue"));

draw(PDF("$(homedir())/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/figs/fig_comparison.pdf", 8inch, 5inch), comparison)

F = prop_out[1,:];
H = prop_out[2,:];
R = prop_out[3,:];
timeseries = plot(layer(x=time_out,y=F,Geom.line,Theme(default_color=colorant"green")),
layer(x=time_out,y=H,Geom.line,Theme(default_color=colorant"orange")),
layer(x=time_out,y=R,Geom.line,Theme(default_color=colorant"blue")));

draw(PNG("/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/figs/fig_timeseries.png", 8inch, 5inch), timeseries)

ts_l = length(time_out);
burnin = convert(Int64,round(ts_l/2,0));

HvsF = plot(x=H[burnin:ts_l],y=F[burnin:ts_l],Geom.point,Theme(default_point_size=0.8pt,default_color=colorant"black",highlight_width = 0pt));

RvsHF = plot(x=R[burnin:ts_l],y=H[burnin:ts_l]+F[burnin:ts_l],Geom.point,Theme(default_point_size=0.8pt,default_color=colorant"black",highlight_width = 0pt));

draw(PNG("/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/figs/fig_HvsF.png", 8inch, 5inch), HvsF)
draw(PNG("/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/DiffusingForager/figs/fig_RvsHF.png", 8inch, 5inch), RvsHF)


time_out, prop_out, N_out = starvingforager_event_spatial(L,dim,initsize,t_term,alpha,K,sigma,rho,lambda,mu,DF,DH);
