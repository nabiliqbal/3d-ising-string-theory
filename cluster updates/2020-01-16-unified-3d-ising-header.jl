
using ProgressMeter
using Distributed

n_PROCS=Sys.CPU_THREADS - 1;
#n_PROCS=1;
addprocs(n_PROCS)
#wp = CachingPool(workers())

#@everywhere include("2018-10-30-binder-parallel02.jl");
#@everywhere include("2018-11-09-3d-ising02.jl");
#@everywhere include("/Users/mcgreevy/Dropbox/READ-directed/3d-ising/2018-11-02-3d-ising01.jl");
@everywhere include("2020-01-08-unified-3d-ising-engine.jl");
@everywhere include("2019-06-21-wall-clusters-no-size-distribution.jl");


SAMPLENUM=n_PROCS;
println("hello.  number of workers = ", n_PROCS)

println("inputted value of gs = ", ARGS[1])
gs = parse(Float32,ARGS[1])
#arglist=map(x->(v = tryparse(Float32,x); isnull(v) ? 1.0 : v), ARGS)
#gs=arglist[1]
#gs=1.0
#println(gs*2)

# Initialize parameters on all threads.
@everywhere begin
    using Nullables
    using Dates
    using DelimitedFiles
    using Statistics



#     NMC0 = 100;
#    NMC0=10000;

#    CUTOFF = 0;
#    gs = 0.9;
## i am being lazy and calling the number of gs values NBETA.

#    BETA = -0.5;
#    phis  = range(-1.5, stop=1.5, length = NBETA)
#    gss= exp.(phis)

#    fil=open("/Users/mcgreevy/Dropbox/READ-directed/3d-ising/lookuptable.dat", "r")
fil=open("lookuptable_sym.dat", "r")
Finput = readlines(fil)
close(fil)
#    lookup=map(x->(v = tryparse(Int8,x); isnull(v) ? 0 : get(v)),input)
#lookup1=map(x->(v = tryparse(Int64,x); isnull(v) ? 0 : get(v)),input)
const lookup=map(x->(v = tryparse(Int64,x); isnull(v) ? 0 : v),Finput)
#    lookup[1]=0;

# this we don't need
fil=open("lookup_connections.dat", "r")
Finput = readlines(fil)
close(fil)
lookup_c=map(x->(v = tryparse(Int64,x); isnull(v) ? 0 : v),Finput)
lookup_c = reshape(lookup_c, (6, 4096))
const lookup_connections = lookup_c'

#    fil=open("/Users/mcgreevy/Dropbox/READ-directed/3d-ising/lookuptable.dat", "r")
#    lookup = lookup1;



end


scheme = "branch-point";

println("regularization scheme = ", scheme);
#scheme = "no-touching";
#NMC0=2^20+50;
#NMC0=2^14+50;
#NMC0=2^11+50;
#NMC0=2^19+50;
NMC0=2^16+50;
NBETA=5;
#   NBETA = 3;
BETAmin = 0.05;
BETAmax = 0.4;
BETAs = range(BETAmin, stop=BETAmax, length=NBETA);

#gs = ARGS[1]

#gs = 1.0;

# these are the linear system sizes
LLIST = [6];
#LLIST = [8];

function jackknife(dataSet)
    n = length(dataSet);
    #assume for now that the estimator is linear
    avg=sum(dataSet)/n;
    subsetVictim = zeros(n-1);
    sesq=0;
    for i = 1:n
        subsetVictim = deepcopy(dataSet);
        deleteat!(subsetVictim, i);
        #println(subsetVictim)
        sesq += (sum(subsetVictim)/(n-1))^2;
#        sesq += sum([(subsetVictim[i] - avg)^2 for i in 1:n-1]);
    end
    buffer = 0.00000000001;
    return (avg, sqrt( (n-1)/n * (sesq - n*avg^2 + buffer) ))
end


# here we actually do the calculating

a=1:SAMPLENUM
binder = zeros(length(LLIST), NBETA);
MAV = zeros(length(LLIST), NBETA);
M2AV = zeros(length(LLIST), NBETA);
EAV = zeros(length(LLIST), NBETA);
E2AV = zeros(length(LLIST), NBETA);
AFMAV = zeros(length(LLIST),NBETA);
chiAV = zeros(length(LLIST),NBETA)
chi2AV = zeros(length(LLIST),NBETA)

#O2AV =  zeros(length(LLIST),NBETA)
O3AV =  zeros(length(LLIST),NBETA)
O4AV =  zeros(length(LLIST),NBETA)
acceptanceRates = zeros(length(LLIST), NBETA);
totalClustersAV = zeros(length(LLIST), NBETA);
#Mhistory = zeros(NMC0-50, length(LLIST), NBETA);


ERRORbinder = zeros(length(LLIST), NBETA);
ERRORMAV = zeros(length(LLIST), NBETA);
ERRORM2AV = zeros(length(LLIST), NBETA);
ERROREAV = zeros(length(LLIST), NBETA);
ERRORE2AV = zeros(length(LLIST), NBETA);
ERRORAFMAV = zeros(length(LLIST),NBETA);
ERRORchiAV = zeros(length(LLIST),NBETA)
ERRORchi2AV = zeros(length(LLIST),NBETA)

#O2AV =  zeros(length(LLIST),NBETA)
ERRORO3AV =  zeros(length(LLIST),NBETA)
ERRORO4AV =  zeros(length(LLIST),NBETA)
ERRORacceptanceRates = zeros(length(LLIST), NBETA);
ERRORtotalClustersAV = zeros(length(LLIST), NBETA);

eqTime=50;
historyLength =NMC0 - eqTime;
#convert(Int64,ceil(NMC0/VOL)-eqTime)
number_of_decimations=(convert(Int64,  floor(log(2,historyLength)))-5)
binnedErrors = zeros(number_of_decimations,length(LLIST), NBETA);


@time  @showprogress 1 "Computing..." for (Liter, L) in enumerate(LLIST)
    NX=convert(Int, L);
    NY=NX;
    NZ=NX;
#    NMCLIST = fill(NMC0*L^3, length(a))
    NMCLIST = fill(NMC0, length(a))
    NXLIST=fill(NX, length(a))
    NYLIST=fill(NY, length(a))
    NZLIST=fill(NZ, length(a))
    gsLIST = fill(gs, length(a))

#    println("  ", NMC0*L^3, "  ");
    println("  ", NMC0, "  ");

    for (Titer, BETA) in enumerate(BETAs)
#    binder[LLIST, :]=sum(pmap( x-> runTheNumbers(NX,NY, Ts, NMC, x),1:SAMPLENUM))/SAMPLENUM;

        BETALIST = fill(BETA, length(a))

#    binder[Liter, Titer]=sum(pmap(runTheNumbers, wp, a, NXLIST,NYLIST,NZLIST, TLIST, NMCLIST))/SAMPLENUM;
        if scheme=="branch-point"
            ploop = pmap(runTheNumbersBP, a, NXLIST,NYLIST,NZLIST, BETALIST, gsLIST, NMCLIST);
        else
            ploop = pmap(runTheNumbersNT, a, NXLIST,NYLIST,NZLIST, BETALIST, gsLIST, NMCLIST);
        end

        (binder[Liter, Titer], ERRORbinder[Liter,Titer]) =  jackknife([ploop[i][1] for i = 1:SAMPLENUM]);
        (MAV[Liter, Titer], ERRORMAV[Liter, Titer]) =  jackknife([ploop[i][2] for i = 1:SAMPLENUM]);
        (EAV[Liter, Titer], ERROREAV[Liter,Titer]) =  jackknife([ploop[i][3] for i = 1:SAMPLENUM]);
        (M2AV[Liter, Titer], ERRORM2AV[Liter, Titer]) =  jackknife([ploop[i][4] for i = 1:SAMPLENUM]);
        (AFMAV[Liter, Titer], ERRORAFMAV[Liter, Titer]) =  jackknife([ploop[i][6] for i = 1:SAMPLENUM]);
        (chiAV[Liter, Titer], ERRORchiAV[Liter,Titer]) =  jackknife([ploop[i][7] for i = 1:SAMPLENUM]);

        (totalClustersAV[Liter, Titer], ERRORtotalClustersAV[Liter,Titer]) =  jackknife([ploop[i][8] for i = 1:SAMPLENUM]);
        (O3AV[Liter, Titer], ERRORO3AV[Liter,Titer]) =  jackknife([ploop[i][9] for i = 1:SAMPLENUM]);
        (O4AV[Liter, Titer], ERRORO4AV[Liter,Titer]) =  jackknife([ploop[i][10] for i = 1:SAMPLENUM]);
        (acceptanceRates[Liter, Titer], ERRORacceptanceRates[Liter,Titer]) =  jackknife([ploop[i][11] for i = 1:SAMPLENUM]);
        (E2AV[Liter, Titer], ERRORE2AV[Liter, Titer]) =  jackknife([ploop[i][12] for i = 1:SAMPLENUM]);
        (chi2AV[Liter, Titer], ERRORchi2AV[Liter,Titer]) =  jackknife([ploop[i][13] for i = 1:SAMPLENUM]);

#        Mhistory[:, Liter, Titer] = sum(ploop[i][13] for i = 1:SAMPLENUM)/SAMPLENUM
        binnedErrors[:,Liter,Titer] = [sqrt(sum((ploop[i][14][j])^2 for i=1:SAMPLENUM)/SAMPLENUM) for j =1:number_of_decimations]


#(binder, MAV, EAV, M2AV, M4AV, AFMAV, chiAV, totalClustersAV, acceptance/NMC)

#    (binder[Liter, Titer], MAV[Liter,Titer])=(sum(ploop[1])/SAMPLENUM,sum(ploop[2])/SAMPLENUM);
# binder, MAV, EAV, M2AV, M4AV, ETAV, AFMAV,acceptance/NMC
    end
end

# save to file
using Dates
using DelimitedFiles

#fil= open(join(["output/", join([Dates.format(Dates.today(), "yyyy-mm-dd"), "binders", reverse(chop(reverse(chop(string([gs, CUTOFF, NMC0, NBETA, LLIST]))))), ".dat"], "-", "")],""), "w")
fil= open(join(["output/", join([Dates.now(), "sign-BETA-$scheme-point-method-vary-BETA", reverse(chop(reverse(chop(string([gs, NMC0, NBETA, LLIST]))))), ".dat"], "-", "")],""), "w")

DelimitedFiles.writedlm(fil, [BETAs, LLIST, binder, MAV, EAV, acceptanceRates, AFMAV, M2AV, chiAV, totalClustersAV, O3AV, O4AV,E2AV, chi2AV, binnedErrors, ERRORbinder, ERRORMAV, ERROREAV, ERRORacceptanceRates, ERRORAFMAV, ERRORM2AV, ERRORchiAV, ERRORtotalClustersAV, ERRORO3AV, ERRORO4AV,E2AV, ERRORchi2AV], "\t")
#	writedlm(fil, [BETA, EAV, MAV, hammings], "\t")
close(fil)


#=
##   make pictures

using PyPlot
#ioff() #turns off interactive plotting
#clf()
pygui(false)  # or else it will crash when running from terminal.

#fig = figure("Binders-crossing",figsize=(15,5)) # Create a new blank figure
#fig[:set_figheight](7) # Doesn't work
#fig[:set_figwidth](3) # Doesn't work

for (Liter, L) in enumerate(LLIST)
    plot(gss, chiAV[Liter, :],"o")
end
legend(LLIST)
xlabel(L"g_s"); ylabel(L"\langle \chi \rangle")
text(2, 0.2, join([L"T = ", string(T)]))

savefig(join(["output/", join([Dates.now(), "average-euler",reverse(chop(reverse(chop(string([T,NMC0, NBETA, LLIST]))))), ".pdf"], "-", "")],""))
#close(fig)
















=#
