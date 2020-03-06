#=

THis uses the normal analysis
=#

# this version loops in phi rather than beta.
using Dates
using Plots
using DelimitedFiles

include("3d-ising-engine.jl");

println("Starting again:");

println("Starting again:");

# create all possible 2 x 2 x 2 unit cells
s = config(2,2,2);

chis = zeros(Int64,128);
walls = zeros(Int64,128);
safeChis = zeros(Int64,128);
safeWalls = zeros(Int64,128);

spinConfs = zeros(2,2,2,128);

o3s = zeros(Float64, 128);
o4s = zeros(Float64,128);
ms = zeros(Float64,128);
sms = zeros(Float64,128);

# runs from 1 to 128 in most general case
for conf = 1:128
	# now set the spins equal to the right things here
	ds = digits(conf-1,base =2, pad = 8);
	for i = 0:1
		for j = 0:1
			for k =0:1
				spinSet(s,i+1,j+1,k+1,(-1)^(ds[fromDigits([i j k])+1]));
			end
		end
	end
	chis[conf] = chi(s);

	walls[conf] = isingEnergy(s)[1];

	if(abs(chis[conf])<100);
		safeChis[conf]=chis[conf];
		safeWalls[conf]=walls[conf];

		ms[conf]=abs(sum(s.spins))/8;
		sms[conf] = abs(staggeredM(s))/8;
		o3s[conf]= orderParameters(s)[1];
		o4s[conf] = orderParameters(s)[2];

		spinConfs[:,:,:,conf]=s.spins;
	end
	#println("Conf: ",ds, "  Chi: ", chi(s), " Walls: ", isingEnergy(s)[1]);
end

#println("Staggered M ", staggeredM(s), "Droplet M ", dropletM(s));

# now find all of them:

maxChilist = findall(x->x==maximum(safeChis),safeChis);
minChilist = findall(x->x==minimum(safeChis),safeChis);
maxWalllist = findall(x->x==maximum(safeWalls),safeWalls);
minWallllist = findall(x->x==minimum(safeWalls),safeWalls);

println("Max chi no-touching: " , maximum(safeChis), " config #: ",findall(x->x==maximum(safeChis),safeChis));
println("Min chi no-touching: config #: ", minimum(safeChis), " config #: ",findall(x->x==minimum(safeChis),safeChis));
println("Max walls no-touching: " , maximum(safeWalls), " config #: ", findall(x->x==maximum(safeWalls),safeWalls));
println("Min walls no-touching: config #: ", minimum(safeWalls), " config #: ",findall(x->x==minimum(safeWalls),safeWalls));

# now print the order parameters for all of them.

for i in maxChilist
	println("conf: ", i, " FM ", ms[i], " AFM ", sms[i], " o3 ", o3s[i]/8, " o4 ", o4s[i]/8);
end
