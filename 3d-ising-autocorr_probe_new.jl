#=

THis uses the normal analysis
=#

# this version loops in phi rather than beta.
using Dates
using Plots
using DelimitedFiles

include("3d-ising-engine.jl");

println("Starting again:");



#Nsteps = 10000000;
#Nsteps = 5000000;
##Nsteps = 20;

#Nsteps = 10000;
# this is almost certainly not enough steps, but oh well.
# and this one:
#sizeLat = 8;
sizeLat = 8;
#sample_rate = 5000;
sample_rate = sizeLat*sizeLat*sizeLat;

#sample_rate = 10;
#Nsteps = 1000;
#Nsteps = 5000*(sizeLat*sizeLat*sizeLat);
Nsteps = 3000*sample_rate;

#Nsteps = 5*sample_rate;
#Nsteps = 20000;
#Nsteps = 10000
s = config(sizeLat,sizeLat,sizeLat);
randomizeSpins(s);

# okay, I think this is now working!
phi = 0;


#beta = 0.5;
#beta = 0.22;

#betas = [0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5];
#betas = [0.5];
#beta = 0.5;
beta = 0.4;

#taus = zeros(length(betas));
#beta = +0.5;

numSamples = convert(Int32,floor(Nsteps/sample_rate))+1

i = 1;
# first, equilibrate for a few steps
mcInit(s,beta,phi,1000*sample_rate);

# then do the data taking load
@time (Marray,Earray,sMarray,o3array,o4array,rejectionRate) = mcRunEfficient(s,beta,phi,Nsteps,sample_rate);
#tau = sum(autocorrM[1:convert(Int32,floor(length(autocorrM)/2))])/autocorrM[1];
#println("tau ", tau);
#display(plot(autocorrM[1:convert(Int32,floor(length(autocorrM)/2))]/autocorrM[1],xlabel="time",ylabel="autocorr", title =string("N = ", sizeLat)));
#display(plot(autocorrM/autocorrM[1],xlabel="time",ylabel="autocorr", title =string("N = ", sizeLat)));
#Plots.savefig(string("autocorr-L-",sizeLat));
#plot(1:length(Marray), Marray);

# okay; now let's see how to do this. for each of them, construct the
filename = string(pwd(), "/data/symmetric/" , today(), "-L-", sizeLat, "-beta-", beta, "-phi-",phi);
writedlm(string(filename,"-Marray.txt"),Marray);
binning(Marray);

#binning(Earray.^2);
#@time (Marray,Earray,sMarray,o3array,o4array,chiArray,rejectionRate,autocorrM) = mcRunErrors(s,beta,phi,Nsteps,sample_rate);
#plot(autocorrM);
