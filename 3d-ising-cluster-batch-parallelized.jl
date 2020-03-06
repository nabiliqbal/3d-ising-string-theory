

# this is a version of the code that is parallelized


using Dates,Plots,Distributed

println("Starting again:");

# figure out how many cores we have, add a few more:
n_PROCS=Sys.CPU_THREADS;
addprocs(n_PROCS-nprocs());

println("Number of workers: ", n_PROCS);

# doing it now with the branch point engine
@everywhere include(string(pwd(),"/nabil code/better topological parallelized/3d-ising-engine.jl"));


# the usual range of beta scan.
#=beta_i = -0.4;
beta_f = 0.4;
nbetas = 9;
betas = range(beta_i, stop=beta_f, length = nbetas);

# recall phiScan syntax:
# phiScan(beta,phi_i,phi_f,nphis,sizeLat,Nsweeps)

# note that I separate it into two bits, one starting at zero and moving up
# the other starting slightly less than zero, moving down
# this is because of issues with equilibration if I start in an extremely
# ordered region and move down.

# doing a run on a 10 x 10 x 10 lattice.
# now try running it with more steps.
@time pmap(x->phiScan(x,0,2,6,12,10000),betas);
@time pmap(x->phiScan(x,-0.4,-2,5,12,10000),betas);=#

# this is the scan that I do for probing the precise binder coeff.
sizes = [6,8,10];
@time pmap(x->phiScan(0,-0.5,-0.3,21,x,20000),sizes)

#sizes = [6,8];
#@time pmap(x->phiScan(0,-0.5,-0.3,11,x,20),sizes)
