




using PyPlot
using DelimitedFiles
using Statistics
using Dates
# here is the script to rename files.
#for file in *; do mv "$file" `echo $file | tr -cd 'A-Za-z0-9_.-'` ; done

## here is the command to change the fonts in the plots.
#matplotlib.rc("font", size=20.0)
#matplotlib.rc("mathtext", fontset= "cm")
#matplotlib.rc("mathtext", rm= "serif")

cd("output-pictures")
#fil=open("/Users/mcgreevy/GoogleDrive/READ-directed/fractons/output/2019-11-13T215523.691-belt-model-2118752-30-6.dat",  "r")
#fil=open("/Users/mcgreevy/GoogleDrive/READ-directed/fractons/output/2019-11-13T221258.108-belt-model-8389408-30-2.dat",  "r")
#fil=open("/Users/mcgreevy/GoogleDrive/READ-directed/fractons/output/2019-11-14T101531.675-belt-model-536970912-30-10.dat",  "r")
#fil=open("/Users/mcgreevy/GoogleDrive/READ-directed/fractons/output/2019-11-14T174756.453-belt-model-67130464-30-6.dat",  "r")
#fil=open("/Users/mcgreevy/GoogleDrive/READ-directed/fractons/output/2019-11-16T090742.788-belt-model-268535456-30-6810.dat",  "r")
#fil=open("/Users/mcgreevy/GoogleDrive/READ-directed/3d-ising/output/2019-11-27T180040.405-sign-BETA-no-touching-method-vary-BETA-ny1.0824256.dat",  "r")
#fil=open("/Users/mcgreevy/GoogleDrive/READ-directed/3d-ising/output/2019-11-28T133839.848-sign-BETA-no-touching-method-vary-BETA-ny1.0104862656.dat",  "r")
#fil=open("/Users/mcgreevy/GoogleDrive/READ-directed/3d-ising/output/2019-12-03T105103.932-sign-BETA-no-touching-method-vary-BETA-ny0.9104862656.dat",  "r")
#fil=open("/Users/mcgreevy/GoogleDrive/READ-directed/3d-ising/output/2019-12-05T141800.366-sign-BETA-no-touching-method-vary-BETA-ny0.8104862656.dat",  "r")
#fil=open("/Users/mcgreevy/GoogleDrive/READ-directed/3d-ising/output/2020-01-14T125324.423-sign-BETA-branch-point-point-method-vary-BETA-ny1.01643456810.dat",  "r")
#fil=open("/Users/mcgreevy/GoogleDrive/READ-directed/3d-ising/output/2020-01-16T190708.096-sign-BETA-branch-point-point-method-vary-BETA-ny1.0209836.dat",  "r")
fil=open("/Users/mcgreevy/GoogleDrive/READ-directed/3d-ising/output/2020-01-24T182733.596-sign-BETA-branch-point-point-method-vary-BETA-ny0.778801524338306810.dat",  "r")
#fil=open("/Users/mcgreevy/GoogleDrive/READ-directed/3d-ising/output/2020-02-02T234109.508-sign-BETA-branch-point-point-method-vary-BETA-ny1.2841524338306810.dat",  "r")
#fil=open("/Users/mcgreevy/GoogleDrive/READ-directed/3d-ising/output/2020-02-06T202117.814-sign-BETA-branch-point-point-method-vary-BETA-ny1.0524338306810.dat",  "r")
#fil=open("/Users/mcgreevy/GoogleDrive/READ-directed/3d-ising/output/2020-02-09T090328.529-sign-BETA-no-touching-method-vary-BETA-ny0.778801524338306810.dat",  "r")
#fil=open("/Users/mcgreevy/GoogleDrive/READ-directed/3d-ising/output/2020-02-12T082247.067-sign-BETA-no-touching-method-vary-BETA-ny1.28403524338306810.dat",  "r")
#fil=open("/Users/mcgreevy/GoogleDrive/READ-directed/3d-ising/output/2020-02-15T075856.767-sign-BETA-no-touching-method-vary-BETA-ny1.0524338306810.dat",  "r")
stuff = readdlm(fil)
close(fil)
NBETA=30;
LLIST = [6,8,10]

#BETAs, LLIST, binder, MAV, EAV, acceptanceRates, AFMAV, M2AV, chiAV, totalClustersAV, O3AV, O4AV,E2AV, chi2AV, binnedErrors, ERRORbinder, ERRORMAV, ERROREAV, ERRORacceptanceRates, ERRORAFMAV, ERRORM2AV, ERRORchiAV, ERRORtotalClustersAV, ERRORO3AV, ERRORO4AV,E2AV, ERRORchi2AV


nDATS = 15+12



(size1,size2)  = size(stuff)
raw_dats = zeros(nDATS, NBETA*length(LLIST))
#raw_dats = zeros(size1, size2)
dats = zeros(nDATS, length(LLIST), NBETA)
for j = 3:nDATS-1
    for k = 1:length(LLIST)*NBETA
        raw_dats[j,k] = stuff[j, k]
        dats[j,:,:] = reshape(raw_dats[j,:], (length(LLIST), NBETA))
    end
end

#NMC=convert(Int64, raw_dats[nDATS])
#NMC = 2^19+50
NMC = 2^19+50



#VOL = LLIST[length(LLIST)]^3
eqTime=50;
#NMC=2^26+eqTime*VOL;

#M = length(mhistory[:, 1])
historyLength = NMC - eqTime #convert(Int64,ceil(NMC)-eqTime)

number_of_decimations=(convert(Int64,  floor(log(2,historyLength)))-5)

#for k = 1:length(LLIST)*NBETA
    binnedErrors = stuff[15, :]

#    raw_dats[j,k] = stuff[j, k]
    binnedErrors = reshape(binnedErrors, (number_of_decimations, length(LLIST), NBETA))
#end

#mhistory = reshape(mhistory, convert(Int64,ceil(NMC/VOL)-eqTime), length(LLIST), NBETA)
##binnedErrors = reshape(binnedErrors, number_of_decimations, length(LLIST), NBETA)


binder = dats[3,:,:]
MAV = dats[4,:,:]
EAV = dats[5,:,:]
#ETAV = dats[6,:]
acceptanceRates = dats[6,:,:]
AFMAV = dats[7,:,:]
M2AV = dats[8,:,:]
chiAV = dats[9,:,:]
clusterAV = dats[10,:,:]
#O2AV = dats[11,:]
O3AV= dats[11,:,:]
O4AV=dats[12,:,:]
E2AV=dats[13,:,:]
chi2AV = dats[14,:,:]


ERRORbinder = dats[16,:,:]
ERRORMAV = dats[17,:,:]
ERRORM2AV= dats[21,:,:]


#chi2AV = dats[14,:,:]

BETAmin = 0.05;
BETAmax = 0.4;
BETAs = range(BETAmin, stop=BETAmax, length=NBETA);

#CUTOFF= 0;
#NMC = 5000*55;
#TC = 4.51;
BETAc = .22;
nu=0.6299
phi = -0.25
#phi=0;
#BETAc
#gs=exp(phi)
gs=.7788
#gs = 1.284
#gs=1.0;

scheme="branch-point"
#scheme="no-touching"







#
#BETAc= .225 # phi = 0 branch point
BETAc = .164; # phi= -.25 branch point
#BETAc = .253;# phi = +.25 branch point
#.22 + 1/2*.25
#.22 - 1/2*.25
#BETAc=.207  # phi = +.25 no touching
#BETAc=.178# phi = 0 no touching
#BETAc=.139   # phi = -.25 no touching
TC = 1/BETAc;
clf()
for (Liter, L) in enumerate(LLIST)
    abscissa = [];
    ordinate =[];
    for (Titer, BETA) in enumerate(BETAs)
        T = 1/BETAs[Titer];
        append!(abscissa, (T-TC)*L^(1/nu));
#        append!(abscissa, (BETA-BETAc)*L^(-1/nu));
        append!(ordinate, binder[Liter, Titer]);
    end
    plot(abscissa, ordinate, ".")
end
xlim([-100,55])
ylim([-.1,1.1])
#text(3, 0.2, join([L"T_c = ", string(TC)]))
text(-90, -.05, join([L"\beta_c = ", string(BETAc), L", \nu = ", string(nu)]), fontsize="20")
#title(join(["Data collapse, ", L"g_s=", string(gs)]), fontsize="20")
title(join(["Data collapse, ", L"\phi = ", phi, ", ", scheme]), fontsize="20")
legend(LLIST)
xlabel(L"(T-T_c)L^{1/\nu}", fontsize="20"); ylabel("binder cumulant", fontsize="20")
#xlabel(L"(\beta-\beta_c)L^{-1/\nu}"); ylabel("binder cumulant")
gcf()

clf()
title(join(["Euler number, ", L"\phi = ", phi, ", ", scheme]), fontsize="20")

#ax2[:set_ylim]([-2,2])  # set plotrange
for (Liter, L) in enumerate(LLIST)
    plot(BETAs, chiAV[Liter,:]/L^3, "o");
end
legend(LLIST)
xlabel(L"\beta", fontsize="20"); ylabel(L"\langle \chi \rangle/L^3", fontsize="20")

blip=5;
#x = fill(BETAc, blip)
x = range(BETAmin,stop=BETAmax,length=blip)
y = fill(0, blip)
plot(x,y)
x = fill(BETAc, blip)
y = range(-.02,stop=.01,length=blip)
plot(x,y)
gcf()
savefig("$(Dates.today())-$scheme-euler-number-crossing-gs=$gs-L=$LLIST-$BETAmin-$BETAmax.pdf", bbox_inches = "tight")
