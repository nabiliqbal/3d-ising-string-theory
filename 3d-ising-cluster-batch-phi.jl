#=

# this is the 3d topological Ising model, v2.0: it now includes the branch points
# when the resolution of the face reconnections are not done correctly.
# it is however extremely inefficient. now try to get this to work faster -- update energy
# by only updating each individual spin rather than all of them.

# this is faster, but still very very slow compared to the usual 3d ising.

# this version has some extra stuff to count clusters.
=#

# this version loops in phi rather than beta.
using Dates
using Plots

include("3d-ising-engine.jl");

println("Starting again:");



sizeLat = 8;
sample_rate = sizeLat*sizeLat*sizeLat;

#Nsteps = 2000*sizeLat^3;
# the below one is the number I usually use.
#Nsteps = 1000000;

#try again with twice as many and see if the answer is different!

#Nsteps = 2000000

# run for 5000 sweeps. this is pretty damn long for a 10x10x10 lattice.
Nsteps = 5000*(sizeLat*sizeLat*sizeLat);
#Nsteps = 10000



s = config(sizeLat,sizeLat,sizeLat);
randomizeSpins(s);


#println("Initial value of droplet order: ",dropletM(s)[2]/(s.NX*s.NY*s.NZ));
# set up a table of phis:
phi_i=2;
phi_f=-2;
nphis =11;

# making the number 11 makes each step rather nice.

#phis = range(phi_i,stop=phi_f,length=nphis);
phis = range(phi_i, stop=phi_f, length = nphis);
#beta = 0.5;
#beta = 0.22;

beta = 0;

# set up a table of things to store the final values
Ms = zeros(nphis);
Es = zeros(nphis);
sMs = zeros(nphis);
#chis = zeros(nphis);
#clusters = zeros(nphis);
o3s = zeros(nphis);
o4s = zeros(nphis);
rejectionRates = zeros(nphis);

to_store = zeros(nphis, 8);

# just equilibrate by running the whole thing at this value of phi
mcInit(s,beta,phi_i,1000*(sizeLat*sizeLat*sizeLat));

# now if it doesnt exist, just construct a directory for this batch, labeled by today's date

if (!isdir(string(pwd(),"/data/",today())))
	mkdir(string(pwd(),"/data/",today()))
end

@time for i=1:length(phis);
	println("On phi ", phis[i])
	@time (Marray,Earray,sMarray,o3array,o4array,rejectionRate) = mcRunEfficient(s,beta,phis[i],Nsteps,sample_rate);


	Ms[i] = mean(Marray);
	Es[i] = mean(Earray);
	sMs[i] = mean(sMarray);
	#chis[i] = mean;
	o3s[i] = mean(o3array);
	o4s[i] = mean(o4array);
	#clusters[i] = sum(clustersArray)/length(clustersArray);
	rejectionRates[i] = rejectionRate;

	# okay; i want to store all the stuff. make the prefix
	filename = string(pwd(), "/data/" , today(), "/",today(), "-L-", sizeLat, "-beta-", beta, "-phi-", phis[i]);

	# now write the Ms files and so all.
	writedlm(string(filename,"-Marray.txt"),Marray);
	writedlm(string(filename,"-sMarray.txt"),sMarray);
	writedlm(string(filename,"-o3array.txt"),o3array);
	writedlm(string(filename,"-o4array.txt"),o4array);
	writedlm(string(filename,"-Earray.txt"),Earray);

	# this is the format in which we output everything
	to_store[i,:]=[beta phis[i] Ms[i] sMs[i] o3s[i] o4s[i] Es[i] rejectionRates[i]];



	println("Rejection rate:", rejectionRate)
	# output something at the very end (towards the interesting large-cluster part)
	if i == length(phis)
		outputData(s);
	end
end

# okay, now I think its time to learn how to write everything to a file!

filename_summary = string(pwd(), "/data/" , today(), "/",today(), "-aaa-L-", sizeLat, "-beta-", beta, "-phi-", phi);
#output the summary to a file
writedlm(string(filename_summary,"-summary.txt"),to_store);


plot_sMs = plot(phis, sMs,xlabel="phi",ylabel="sM");
plot_Ms= plot(phis, Ms,xlabel="phi",ylabel="M");
plot_reject = plot(phis, rejectionRates,xlabel="phi",ylabel="rejection rate");
plot_o3 = plot(phis, o3s,xlabel="phi",ylabel="o3");
plot_o4 = plot(phis, o4s,xlabel="phi",ylabel="o4");
plot_Es = plot(phis, Es, xlabel="phi",ylabel="E");

plot(plot_sMs, plot_Ms, plot_o3, plot_o4, plot_reject, plot_Es, layout = (3,2))
#title(string("L = ", sizeLat, " beta = ", beta));
savefig(string(filename,"-plots.pdf"));
