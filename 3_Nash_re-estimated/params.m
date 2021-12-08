function [ydiffdataPercent,idiffdataPercent,alfa,rho,u,Q,sigma,pibar,betta,lambda,deltak,tau,nuf,...
    etag,recSHAREpercent,b,xi,DSHARE,deltapercent,Spp,sigmaa,rhoR,rpi,ry,sig_epsR,kappaf,...
    iota,doNash,doEHL,sigmaL,xiw,lambdaw,kappaw,rhomupsi,rhomuz,sig_mupsi,sig_muz,...
    profy,dolagZBAR,varkappaw,pibreve,thetaw,varkappaf,M,thetaG, thetaPHI,... 
    thetaKAP, thetaD, thetaGAM, do_given_bets,bet1,bet2,bet3,searchSHAREpercent]=params
 
%params.m

%model switches

doNash=1; %if=0, AOB; 
          %if=1, Nash bargaining
doEHL=0; %if=0, AOB or Nash bargaining; 
         %if=1, EHL labor market

do_given_bets=0; %if =1, given values for beta's in sharing rule
bet1 = 0.0907; %AOB sharing rule coefficients
bet2 = 28.9219;
bet3 = 0.4562;

%Labor market parameters

u=0.055; %unemployment rate
rho=0.9; %job survival rate
sigma=0.5570; %matching function share of unemployment.
recSHAREpercent=0.5; %hiring cost parameter
searchSHAREpercent=0.05; %search cost parameter

DSHARE=0.368861073612840; %replacement ratio
Q=0.7; %vacancy filling rate

deltapercent=0.3022; %probability of bargaining session determination
M=60; %Maximum bargaining rounds per quarter 
 
if M<2, error('M must be at least equal to 2');end
if mod(M,2)==1, error('M must be an even number');end 

%prices  

xi=0.5841; %Calvo prices
pibar=1.00625; %inflation rate, gross, quarterly
kappaf=0; %price indexation to past inflation
varkappaf=1; %price indexation parameter 
pibreve=1; %price indexation to value pibreve
%if kappaf=0 and varkappaf=1 and pibreve=1, no indexation at all

lambda=1.4256; %steady state price markup
nuf=1; %working capital fraction
tau=1; %steady state markup shifter

%technology and adjustment cost

alfa=0.2366; %capital share
deltak=0.025; %depreciation rate of capital
Spp=13.6684; %second derivative of investment adjustment costs
sigmaa=0.0760; %capacity utilization costs

%growth rates 

idiffdataPercent=2.9; %growth rates of investment
ydiffdataPercent=1.7; %growth rates of real GDP

%Preferences

betta=0.996783170280770; %discount factor households
b=0.8320; %habit formation in consumption

%Monetary Policy by Taylor rule

rhoR=0.8555; %Interest rate smoothing
rpi=1.3552; %Coefficient on inflation
ry=0.0359; %Coefficient on GDP

%Government

etag=0.2; %Government spending-to-output in steady state

%technology diffusion into gov. spending and other parameters of the model

thetaG =0.0136; %government spending
thetaPHI=thetaG; %fixed cost of production
thetaKAP=thetaG; %hiring cost
thetaGAM=thetaG; %cost of counteroffer
thetaD=0.7416; %replacement ratio
dolagZBAR=1; %if =1, zbar_t=(zplus_t-1)^thetaj*(zbar_t-1)^(1-thetaj); 
             %if =0, zbar_t=(zplus_t)^thetaj*(zbar_t-1)^(1-thetaj) for j=G,GAM,KAP,BU,PHI

% steady state profits as share of output

profy=0;

%EHL parameter

lambdaw=1.2; %wage markup 
xiw=0.75; %Calvo wage stickiness 
sigmaL=1; %inverse labor supply elasticcity
kappaw=0; %wage indexation to past inflation  
varkappaw=1; %wage indexation parameter
thetaw=0; %wage indexation to 
%if kappaw=0 and varkappaw=1 and pibreve=1 and thetaw=0, no indexation at all

% standard deviation of shocks

sig_muz=0.1362; %unit root neutral technology
sig_mupsi=0.1100; %unit root investment technology
sig_epsR=0.6028; %monetary policy

sig_epsilon=0; %stationary neutral technology
sig_upsilon=0; %stationary investment technology
sig_g=0; %government spending
sig_taud=0; %price markup
sig_zetac=0; %consumption preference

%AR(1)

rhomupsi=0.7279; %unit root investment technology
rhomuz=0; %unit root neutral technology

rhoepsilon=0; %stationary neutral technology 
rhoupsilon=0; %stationary investment technology
rhozetac=0; %consumsption preference
rhog=0; %government spending
rhotaud=0; %price markup
rhosigmam=0; %matching function shock

iota=1;      
 
end


