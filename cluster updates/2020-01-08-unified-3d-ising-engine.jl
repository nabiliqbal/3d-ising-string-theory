



function initSpinsHot(NX, NY,NZ)
	s = zeros(Int,NX,NY,NZ);
# choose initial spin config.
  	for i=1:NX, j = 1:NY, k = 1:NZ
    	s[i,j,k] = rand([-1,1]);
	end
  	return s
end


function initSpinsCold(NX, NY,NZ)
	s = zeros(Int, NX,NY,NZ);
# choose initial spin configuration iid randomly.
  	for i=1:NX, j = 1:NY, k = 1:NZ
	     s[i,j,k] = 1;
    end
    return s
end

#list = Array{Int64}(2)
#print(list[1,1])

function nbrs_cubic(i,j,k,Nx,Ny,Nz)
  #list = [zeros(2,2)]
  nbrs = Array{Int}[];
  if i > 1 push!(nbrs, [i-1, j,k])
  else push!(nbrs, [Nx, j,k])
  end
  if j > 1 push!(nbrs, [i, j-1,k])
  else push!(nbrs, [i, Ny,k])
  end
  if k > 1 push!(nbrs, [i,j,k-1])
  else push!(nbrs, [i,j,Nz])
  end
  if i < Nx push!(nbrs, [i+1,j,k])
  else push!(nbrs, [1, j,k])
  end
  if j < Ny push!(nbrs, [i, j+1,k])
  else push!(nbrs, [i, 1,k])
  end
  if k < Nz push!(nbrs, [i, j,k+1])
  else push!(nbrs, [i, j,1])
  end
    return nbrs
end


function shift(i,j,k,a, NX,NY,NZ)
	return a==0 ? [periodicAdd(i,NX),j,k] : ( a==1 ? [i,periodicAdd(j,NY),k] : [i,j,periodicAdd(k,NZ)] )
end

function periodicAdd(i, N)
	if i<N
		i = i+1
    else
		i=1
	end
	return i
end

function periodicSubtract(i, N)
	if i>1
        i = i-1
	else
		i=N
	end
	return i
end



function shiftReverse(i,j,k,a, NX,NY,NZ)
	return a==0 ? [periodicSubtract(i,NX),j,k] : ( a==1 ? [i,periodicSubtract(j,NY),k] : [i,j,periodicSubtract(k,NZ)] )
end

## if BETA< 0, try to make AFM clusters
function makeClusterSigned(sig, BETA, NX::Int,NY::Int,NZ::Int)
    prob = 1.0 - exp(- 2*abs(BETA));
    signBETA= BETA/abs(BETA);
	nCLUSTERS = 1;   # number of clusters tried.
	# pick a random site to start the cluster
	flag = zeros(Int64, NX,NY,NZ);
	kx=convert(Int, rand(1:NX));
	ky=convert(Int, rand(1:NY));
	kz=convert(Int,rand(1:NZ));
	pocket = [[kx,ky,kz]]
	cluster =[[kx,ky,kz]];
	flag[kx,ky,kz] = nCLUSTERS;
	while length(pocket) > 0
		pock=rand(1:length(pocket))  # pick a random element of the pocket.
		i=pocket[pock][1]; j = pocket[pock][2]; k = pocket[pock][3];
		bers = nbrs_cubic(i,j,k, NX,NY,NZ);
		for sites in bers
			# print(sites)
			if sig[sites[1],sites[2],sites[3]]== signBETA*sig[i,j,k]
				if flag[sites[1],sites[2],sites[3]] != nCLUSTERS
					if rand() < prob
						push!(pocket, sites);
						push!(cluster, sites);
						flag[sites[1],sites[2],sites[3]] = nCLUSTERS;
					end
				end
			end
		end
		deleteat!(pocket, pock)  # remove from the pocket
	end
    return (cluster, flag)
end


## given the spin config on the primal lattice, s, output the domain wall configuration
## p[i,j,k,a+1] = 0,1.  a=0,1,2 =x,y,z is the normal to the plaquette in the unit cell
## of the dual lattice at i,j,k.
## we should only need to run this at the beginning of the simulation.
function domainWallConfig(s::Array{Int,3},NX::Int,NY::Int,NZ::Int)
	p = zeros(Int8, NX,NY,NZ,3);
	for i=1:NX, j=1:NY, k = 1:NZ, a=0:2
	    (i1,j1,k1) = shiftReverse(i,j,k,a,NX,NY,NZ);
	    p[i,j,k,a+1] = convert(Int8, (1 - s[i,j,k]*s[i1,j1,k1])/2);
#		if p[i,j,k,a+1] >0 println(i,j,k,a+1) end
	end
	return p
end

## this function outputs the ps configuration associated with the dual-lattice site ijk
function getSiteConfig(i,j,k, p,NX,NY,NZ)
# (where a=0,1,2 labels the direction in which the link points)
	ps = zeros(Int8, 12);
	for a = 0:2
		b = mod(a + 1,3); c = mod(a+2,3);
		ps[a+1]= p[i,j,k,a+1];
		(i1,j1,k1) = shiftReverse(i,j,k,a,NX,NY,NZ);
		ps[a+4] = p[i1,j1,k1,b+1];
		ps[a+7] = p[i1,j1,k1,c+1];
		(i1,j1,k1) = shiftReverse(i1,j1,k1,b,NX,NY,NZ);
		ps[a+10] = p[i1,j1,k1,c+1];
	end
    return ps
#	return ESite(ps);
#	siteenergy = lookup[fromDigits(ps)+1]
#	if siteenergy % 2 == 1 println("site ", i,j,k, " is in configuration ", fromDigits(ps)+1, " which gives ", siteenergy) end
#    return siteenergy
end





function getSiteEnergy(p,NX,NY,NZ)
    collision=0;
    #nubbinData = zeros(Int8, NX,NY,NZ,6);
    siteEnergy = zeros(Int8,NX,NY,NZ);
    for i=1:NX,j=1:NY,k=1:NZ
        ps = getSiteConfig(i,j,k, p,NX,NY,NZ);
#        for b=1:6
#            nubbinData[i,j,k,b] =lookup_connections[fromDigits(ps)+1,b];
#        end
#        nubbinData[i,j,k,b] =lookup_connections[fromDigits(ps)+1,b];
        siteEnergy[i,j,k] = lookupCUTOFF[fromDigits(ps)+1];
        if siteEnergy[i,j,k] < 0
            collision = 1;
        end
    end
    return (siteEnergy, collision)
end




function ELink(ps)
    return ps[1]*ps[2]*(1-ps[3])*(1-ps[4])+ps[3]*ps[2]*(1-ps[1])*(1-ps[4])+
        ps[4]*ps[2]*(1-ps[3])*(1-ps[1]) + ps[1]*ps[3]*(1-ps[2])*(1-ps[4])+
        ps[1]*ps[4]*(1-ps[3])*(1-ps[2]) +ps[3]*ps[4]*(1-ps[1])*(1-ps[2]) + 2*ps[1]*ps[2]*ps[3]*ps[4];
end
## extracts from p the four faces in which the link (ijka) sits
## looks up the energy of this configuration.
function linkCounting(i,j,k,a, p,NX,NY,NZ)
# (where a=0,1,2 labels the direction in which the link points)
	ps = zeros(Int8, 4);
	b = mod(a + 1,3); c = mod(a+2,3);
	ps[1]= p[i,j,k,b+1];
	ps[2]= p[i,j,k,c+1];
	(i1,j1,k1) = shiftReverse(i,j,k,b,NX,NY,NZ); # THIS IS CORRECT AND IT MAKES A BIG DIFFERENCE!
	ps[3] = p[i1,j1,k1, c+1];
	(i1,j1,k1) = shiftReverse(i,j,k,c,NX,NY,NZ);
	ps[4] = p[i1,j1,k1, b+1];
	#return ELink[p1,p2,p3,p4] ;
	return ELink(ps)
end



function fromDigits(p)
    len=length(p);
    sum([p[k]*2^(len-k) for k=1:len]);
end



function latticeToCell(x,y,z)  ## give the label of the site in the 2x2x2 unit cell
    X = mod.([z,y,x],2)  # beware flipped convention!!
    return fromDigits(X)
end

function orderParameters(s, NX,NY,NZ)
# fixed 2019-12-16


	sqrt3 = sqrt(3)
	sqrt6 = sqrt(6)
	sqrt8 = sqrt(8)


  vplus = [ .5 0 0 -.5 -.5 0 0 .5;
  	  -.5/sqrt3 1/sqrt3 0 -.5/sqrt3  -.5/sqrt3  0 1/sqrt3 -.5/sqrt3;
	  -.5/sqrt6 -.5/sqrt6 sqrt3/sqrt8  -.5/sqrt6 -.5/sqrt6 sqrt3/sqrt8 -.5/sqrt6  .5/sqrt6]

 vminus = [-.5 0 0 .5 -.5 0 0 .5;
 	   -.5/sqrt3 -1/sqrt3 0 -.5/sqrt3 .5/sqrt3 0 1/sqrt3 .5/sqrt3;
	  -.5/sqrt6 .5/sqrt6 -sqrt3/sqrt8  -.5/sqrt6 .5/sqrt6 sqrt3/sqrt8 -.5/sqrt6  .5/sqrt6]

# better
#    v3=[ 1 0 0 -1 -1 0 0 1;
#        0 1 0 -1 -1 0 1 0 ;
#        0 0 1 -1 -1 1 0 0]

#    v4= [-1 0 0 1 -1 0 0 1;
#        0 -1 0 -1 1 0 1 0;
#        0 0 -1 -1 1 1 0 0]

    #vAFM = [1 -1 -1 1 -1 1 1 -1]

    return (
	# this is AFM, as a check.
# 	sum(vAFM[latticeToCell(x,y,z)+1]*s[x,y,z] for x=1:NX, y=1:NY, z=1:NZ)^2,
    #This is Oplus:
	sum(sum(vplus[a,latticeToCell(x,y,z)+1]*s[x,y,z] for x=1:NX, y=1:NY, z=1:NZ)^2 for a=1:3),
    #This is Ominus:
    sum(sum(vminus[a,latticeToCell(x,y,z)+1]*s[x,y,z] for x=1:NX, y=1:NY, z=1:NZ)^2 for a=1:3)
	)
end




function decimate(timeSeries)
    M = length(timeSeries)
    halfM=convert(Int64, floor(M/2))
    binned = zeros(halfM)
    for i =1:halfM
        binned[i] = (timeSeries[2*i-1] + timeSeries[2*i])/2;
    end
    return (binned)
end

function binning(timeSeries, number_of_decimations)
    M = length(timeSeries)
#    number_of_decimations=(convert(Int64,  floor(log(2,M)))-3)
    Deltas = zeros(number_of_decimations)
    for decimation_step = 1:number_of_decimations
        Ml = length(timeSeries)
        Abar = mean(timeSeries)
        Deltas[decimation_step] = sqrt(sum( (timeSeries[i] - Abar )^2  for i = 1:Ml)/(Ml*(Ml-1)))
        (timeSeries) = decimate(timeSeries);
    end

    return Deltas
end




### above here is shared between the two methods.


function getNubbinData(p,NX,NY,NZ)
    nubbinData = zeros(Int8, NX,NY,NZ,6);
    siteEnergy = zeros(Int8,NX,NY,NZ);
    for i=1:NX,j=1:NY,k=1:NZ
        ps = getSiteConfig(i,j,k, p,NX,NY,NZ);
        for b=1:6
            nubbinData[i,j,k,b] =lookup_connections[fromDigits(ps)+1,b];
        end
#        nubbinData[i,j,k,b] =lookup_connections[fromDigits(ps)+1,b];
        siteEnergy[i,j,k] = lookup[fromDigits(ps)+1];
    end
    return (nubbinData, siteEnergy)
end



function linkEnergy(i,j,k, a, p, nubbinData,NX,NY,NZ)
## the conflicted edge is (5,4,2, 2).
#    (i1,j1,k1) = shiftReverse(i,j,k, a,NX,NY,NZ)
    (i1,j1,k1) = shift(i,j,k, a,NX,NY,NZ)
    nubbin1=nubbinData[i,j,k,2*(a+1)]  # a+
    nubbin2=nubbinData[i1,j1,k1, 2*(a+1)-1] # a-
#    nubbin1=nubbinData[i,j,k,2*(a+1)-1]  # a-
#    nubbin2=nubbinData[i1,j1,k1, 2*(a+1)] # a+

#    (nubbin1 == 0 ? 0 : (nubbin1 == 1 ? 1 : (nubbin1 != nubbin2 ? 3 : 2)))
    return (nubbin1 == 0 ? 0 : (nubbin1 == 1 ? 1 : (nubbin1 != nubbin2 ? 3 : 2)))
end

function energyIsingBP(s,p, nubbinData, siteEnergy,NX, NY, NZ)
    #EE = 0;	MM =0; #	hh=0;
	AFM=0;
	edges =0;	faces =0; vertices=0
    for i = 1:NX, j = 1:NY, k = 1:NZ
#	    hh=0;
#		MM+= ss[i,j,k];
		AFM += (-1)^(i+j+k)*s[i,j,k];
#		hh += ss[periodicAdd(i,NX),j,k];
#		hh += ss[i,periodicAdd(j,NY),k];
#        hh += ss[i,j,periodicAdd(k,NZ)];
#    	EE -= ss[i,j,k]*hh;
#        EE += siteEnergy(p, i,j,NX,NY);
        for a =0:2
#			chicken = linkCounting(i,j,k,a, p,NX,NY,NZ)
            chicken = linkEnergy(i,j,k,a,p, nubbinData, NX,NY,NZ)
#			if chicken >0 println("edge ", i,j,k, " ", a, " counted for ", chicken);end
#            ET += chicken;
			edges += chicken
        end
		chicken = siteEnergy[i,j,k]
        #chicken = siteCounting(i,j,k, p,NX,NY,NZ)
#		if chicken >0 println("site ", i,j,k, " counted for ", chicken);end
#        ET += gV*chicken;
		vertices += chicken
	 end
#     ET += gF*sum(p);
	 faces = sum(p)
#	 println("edges = ", edges, " faces = ", faces, " vertices = ", vertices)
	 #(2*sum(p)- 3*NX*NY*NZ);  # sum(p)
	 EE= 2*sum(p)- 3*NX*NY*NZ
	 MM = sum(s)
 ## i know this is redundant, but for some reason i want to keep it.
	 return (EE,faces-edges+vertices,MM,AFM)
end

function siteEnergyAndNubbinData(i,j,k,p,NX,NY,NZ)
    ps = getSiteConfig(i,j,k, p,NX,NY,NZ);
    nubbinData = zeros(Int, 6)
    for b=1:6
        nubbinData[b] =lookup_connections[fromDigits(ps)+1,b];
    end
#        nubbinData[i,j,k,b] =lookup_connections[fromDigits(ps)+1,b];
    siteEnergy = lookup[fromDigits(ps)+1];
    return (siteEnergy, nubbinData)
end


function mcStepClusterShortBP(s,p,nubbinData, siteEnergy, E, M, AFM, chi,gs,BETA, NX::Int, NY::Int,NZ::Int)
	(cluster, flag)=makeClusterSigned(s,BETA,NX,NY,NZ);
	return updateClusterBruteForceBP(s,p,nubbinData, siteEnergy, cluster,flag, E,M,AFM,chi,gs,BETA,NX,NY,NZ)
end



function updateClusterBruteForceBP(s,p,nubbinData, siteEnergy, cluster,flag, E,M,AFM,chi,gs, BETA, NX,NY,NZ)
	accept =0;
   stemp = deepcopy(s)
   for flips in cluster
	   (i,j,k) = flips
	   stemp[i,j,k] = - stemp[i,j,k]
#	   println("flipping ", i,j,k)
#	   M -= 2*s[flips[1], flips[2],flips[3]];
   end  # update domain wall configuration
# TOTAL BRUTE FORCE METHOD
   ptemp = domainWallConfig(stemp, NX,NY,NZ)
#   println("the updated magnetization will be ", sum(stemp))
#println("p is      ", listOfps(p))
#println("temp p is ", listOfps(ptemp))

## brute force method:
    (nubbinDataTemp,siteEnergyTemp) = getNubbinData(ptemp,NX,NY,NZ)
	(E1, chi1, M1,AFM1)=energyIsingBP(stemp,ptemp,nubbinDataTemp, siteEnergyTemp, NX,NY,NZ)
	DeltaChi = chi1-chi;

#	println(DeltaET);
    ### here is where we decide to flip the cluster or not.
#    Ups = exp(DeltaET) # = gs^{\Delta \chi};  CHECK THE SIGN.  ok the sign was wrong.
	Ups  = gs^(-DeltaChi)
#	println("probability to flip was ", Ups)
    if rand() < Ups
		accept = 1;
#        println("flipping!")
		s = stemp; p = ptemp; nubbinData = nubbinDataTemp; siteEnergy = siteEnergyTemp;
		E= E1; chi=chi1; M = M1;AFM=AFM1;
    end


	return (s, p, nubbinData, siteEnergy, E, M,AFM, chi, accept)
end




# s, p, nubbinData, siteEnergy, E, ET, M,AFM,ETOTAL,
#                     s, p, nubbinData, siteEnergy, E, ET, M,AFM,ETOTAL, ETTOTAL,MTOTAL, MTOTAL2, MTOTAL4, AFMTOTAL,chiTOTAL,BETA, gF,gE,gV,NX, NY, NZ, NMC
function mcRUNCLUSTERBP(s, p, nubbinData, siteEnergy, E, chi,M,AFM, ETOTAL,E2TOTAL,MTOTAL, M2TOTAL, MTOTAL4,AFMTOTAL, chiTOTAL, chi2TOTAL, totalClustersTOTAL,O3TOTAL,O4TOTAL, BETA, gs,NX::Int, NY::Int, NZ::Int, NMC) # s = initial config.
	 	      	  # could be the last one of the previous temp. E, M = its energy..
	measurementCounter=0;
	chiMeasurementCounter=0;
	OMeasurementCounter=0;
	VOL = (NX*NY*NZ);
	acceptance=0;
	equilibrateWait = 50;
	Mhistory = zeros(convert(Int64, floor(NMC-equilibrateWait)));
	for itMC = 1:NMC
#		(s, p, nubbinData, siteEnergy, E, M,AFM, chi, accept)
		(s, p, nubbinData, siteEnergy, E, M,AFM, chi, accept)=
		mcStepClusterShortBP(s, p,nubbinData, siteEnergy, E, M, AFM, chi,gs,BETA,NX, NY,NZ)
		if accept == 1 acceptance += 1; end
		# measure  (here the measurement is easy).
		if itMC > 50  ## equilibrate
	  		ETOTAL += E/VOL;
			E2TOTAL += (E/VOL)^2
	  		MTOTAL  += abs(M)/VOL;
			M2TOTAL += (M/VOL)^2;
			MTOTAL4 += (M/VOL)^4;
#			ETTOTAL += ET/VOL;
			AFMTOTAL+= (AFM/VOL)^2;  # AFM is extensive.  we square it so it doesn't average to zero.
#			chiTOTAL += avgChi
			measurementCounter +=1;
			Mhistory[measurementCounter] = abs(M)/VOL

			if itMC % 5 == 0  # how often do we want to do this?
				(o3,o4) = orderParameters(s,NX,NY,NZ)
#				O2TOTAL += o2/VOL^2;
				O3TOTAL += o3/VOL^2;  # O3,4 are squared
				O4TOTAL += o4/VOL^2;
				OMeasurementCounter +=1;
			end
		end
		if itMC % 50 == 49
			(label, totalClusters)= makeClusterLabels(p,nubbinData,NX,NY,NZ)
			avgChi = (totalClusters!=0) ?  chi/totalClusters : 0;
			if avgChi > 2
				println("something bad happened, average chi > 2")
				fil= open(join(["output/", join([Dates.now(), "problem-config-for-branch-point", reverse(chop(reverse(chop(string([NX,NY,NZ]))))), ".dat"], "-", "")],""), "w")

				DelimitedFiles.writedlm(fil, [s], "\t")
				#	writedlm(fil, [BETA, EAV, MAV, hammings], "\t")
				close(fil)

			end
			chiTOTAL += avgChi
			chi2TOTAL += avgChi^2
			totalClustersTOTAL += totalClusters
			chiMeasurementCounter+=1;
		end

		#println(s)
#		if mod(itMC, 100)== 0 push!(history, copy(s)) end
	end # itmc
	 return (s, p, nubbinData, siteEnergy, E,chi,M,AFM,ETOTAL, E2TOTAL, MTOTAL, M2TOTAL, MTOTAL4, AFMTOTAL,chiTOTAL, chi2TOTAL, totalClustersTOTAL, O3TOTAL,O4TOTAL,OMeasurementCounter, chiMeasurementCounter, measurementCounter, acceptance, Mhistory)
end



function energyIsingNT(s,p, siteEnergy,NX, NY, NZ)
    #EE = 0;	MM =0; #	hh=0;
	AFM=0;
	edges =0;	faces =0; vertices=0
    for i = 1:NX, j = 1:NY, k = 1:NZ
#	    hh=0;
#		MM+= ss[i,j,k];
		AFM += (-1)^(i+j+k)*s[i,j,k];
#		hh += ss[periodicAdd(i,NX),j,k];
#		hh += ss[i,periodicAdd(j,NY),k];
#        hh += ss[i,j,periodicAdd(k,NZ)];
#    	EE -= ss[i,j,k]*hh;
#        EE += siteEnergy(p, i,j,NX,NY);
        for a =0:2
			chicken = linkCounting(i,j,k,a, p,NX,NY,NZ)
#            chicken = linkEnergy(i,j,k,a,p, nubbinData, NX,NY,NZ)
#			if chicken >0 println("edge ", i,j,k, " ", a, " counted for ", chicken);end
#            ET += chicken;
			edges += chicken
        end
		chicken = siteEnergy[i,j,k]
        #chicken = siteCounting(i,j,k, p,NX,NY,NZ)
#		if chicken >0 println("site ", i,j,k, " counted for ", chicken);end
#        ET += gV*chicken;
		vertices += chicken
	 end
#     ET += gF*sum(p);
	 faces = sum(p)
#	 println("edges = ", edges, " faces = ", faces, " vertices = ", vertices)
	 #(2*sum(p)- 3*NX*NY*NZ);  # sum(p)
	 EE= 2*sum(p)- 3*NX*NY*NZ
	 MM = sum(s)
 ## i know this is redundant, but for some reason i want to keep it.
	 return (EE,faces-edges+vertices,MM,AFM)
end

function mcStepClusterShortNT(s,p,siteEnergy, E, M, AFM, chi,gs,BETA, NX::Int, NY::Int,NZ::Int)
	(cluster, flag)=makeClusterSigned(s,BETA,NX,NY,NZ);
	return updateClusterBruteForce(s,p,siteEnergy, cluster,flag, E,M,AFM,chi,gs,BETA,NX,NY,NZ)
end

function updateClusterBruteForceNT(s,p,siteEnergy, cluster,flag, E,M,AFM,chi,gs, BETA, NX,NY,NZ)

	# this actually does not depend on BETA at all.
	accept =0;
   stemp = deepcopy(s)
   for flips in cluster
	   (i,j,k) = flips
	   stemp[i,j,k] = - stemp[i,j,k]
#	   println("flipping ", i,j,k)
#	   M -= 2*s[flips[1], flips[2],flips[3]];
   end  # update domain wall configuration
# TOTAL BRUTE FORCE METHOD
   ptemp = domainWallConfig(stemp, NX,NY,NZ)
#   println("the updated magnetization will be ", sum(stemp))
#println("p is      ", listOfps(p))
#println("temp p is ", listOfps(ptemp))

## brute force method:
    (siteEnergyTemp, collision) = getSiteEnergy(ptemp,NX,NY,NZ)
	(E1, chi1, M1,AFM1)=energyIsingNT(stemp,ptemp,siteEnergyTemp, NX,NY,NZ)
	DeltaChi = chi1-chi;

#	println(DeltaET);
    ### here is where we decide to flip the cluster or not.
	Ups  = gs^(DeltaChi)

    ## modify this to include a penalty if siteEnergy contains collisions.
#	println("probability to flip was ", Ups)
    if rand() < Ups
        if collision == 0
		      accept = 1;
#        println("flipping!")
		      s = stemp; p = ptemp; siteEnergy = siteEnergyTemp;
		      E= E1; chi=chi1; M = M1;AFM=AFM1;
        else
#              println("Collision! No touching!")
        end
    end


	return (s, p, siteEnergy, E, M,AFM, chi, accept)
end




function mcRUNCLUSTERNT(s, p, siteEnergy, E, chi,M,AFM, ETOTAL,E2TOTAL,MTOTAL, MTOTAL2, MTOTAL4,AFMTOTAL, chiTOTAL, totalClustersTOTAL, O3TOTAL,O4TOTAL, BETA, gs,NX::Int, NY::Int, NZ::Int, NMC) # s = initial config.
	 	      	  # could be the last one of the previous temp. E, M = its energy..
	measurementCounter=0;
	chiMeasurementCounter=0;
	OMeasurementCounter=0;
	VOL = (NX*NY*NZ);
	acceptance=0;
	equilibrateWait = 50;
	Mhistory = zeros(convert(Int64, floor(NMC-equilibrateWait)));
	for itMC = 1:NMC
		(s, p, siteEnergy, E, M,AFM, chi, accept)=
		mcStepClusterShortNT(s, p,siteEnergy, E, M, AFM, chi,gs,BETA,NX, NY,NZ)
		if accept == 1 acceptance += 1; end
		# measure  (here the measurement is easy).
		if itMC > equilibrateWait  ## equilibrate
	  		ETOTAL += E/VOL;
			E2TOTAL += (E/VOL)^2;
	  		MTOTAL  += abs(M)/VOL;
			MTOTAL2 += (M/VOL)^2;
			MTOTAL4 += (M/VOL)^4;
#			ETTOTAL += ET/VOL;
			AFMTOTAL+= (AFM/VOL)^2;  #
#			chiTOTAL += avgChi
			measurementCounter +=1;
			Mhistory[measurementCounter] = abs(M)/VOL
			if itMC % 5 == 0  # how often do we want to do this?
				(o3,o4) = orderParameters(s,NX,NY,NZ)
#				O2TOTAL += o2/VOL^2;
				O3TOTAL += o3/VOL^2;  # O3,4 are squared
				O4TOTAL += o4/VOL^2;
				OMeasurementCounter +=1;
			end


		end
		if itMC % 50 == 49
			(label, correctTotalClusters) = makeClusterLabelsNoTouching(p,NX,NY,NZ)
#			correctTotalClusters = max(totalClusters, totalUpClusters)
#			avgChi = (totalClusters!=0) ?  chi/totalClusters : 0;
			avgChi = (correctTotalClusters != 0 ? chi/correctTotalClusters : 0)
			if avgChi > 2
				println("something bad happened, average chi > 2: ", "chitotal = ", chi,
				", correct total clusters=", correctTotalClusters,
				", correct chi =", avgChi)
				fil= open(join(["output/", join([Dates.now(), "problem-config-for-no-touching", reverse(chop(reverse(chop(string([NX,NY,NZ]))))), ".dat"], "-", "")],""), "w")

				DelimitedFiles.writedlm(fil, [s], "\t")
				#	writedlm(fil, [BETA, EAV, MAV, hammings], "\t")
				close(fil)

			end
			chiTOTAL += avgChi
			totalClustersTOTAL += correctTotalClusters
			chiMeasurementCounter+=1;
		end

		#println(s)
#		if mod(itMC, 100)== 0 push!(history, copy(s)) end
	end # itmc
	 return (s, p, siteEnergy, E,chi,M,AFM,ETOTAL, E2TOTAL, MTOTAL, MTOTAL2, MTOTAL4, AFMTOTAL,chiTOTAL,totalClustersTOTAL, O3TOTAL,O4TOTAL,OMeasurementCounter, chiMeasurementCounter, measurementCounter, acceptance, Mhistory)
end


function runTheNumbersNT(j::Int64, NX::Int,NY::Int,NZ::Int, BETA, gs, NMC)
	#function runTheNumbers(NX,NY, Ts, NMC, j::Int64)
#	  BETA = 1/T;
#	  gF = log(gs)
#	  gE = -log(gs)
#	  gV = log(gs)
	  s = initSpinsCold(NX,NY,NZ);  # initial spin config.
	  p = domainWallConfig(s, NX,NY,NZ)
      (siteEnergy, collision) = getSiteEnergy(p,NX,NY,NZ)
      if collision == 1
          print("There is a collision!  No touching!")
      end

## here we must get nubbinData and siteEnergies.

	  (E,chi,M,AFM) = energyIsingNT(s, p,siteEnergy,NX,NY,NZ) ;
#	  history = [copy(s)]
	  ETOTAL=0; E2TOTAL=0;MTOTAL=0; MTOTAL2=0; MTOTAL4=0;AFMTOTAL=0; chiTOTAL=0;
	  totalClustersTOTAL=0;O3TOTAL=0; O4TOTAL=0;
	  #																													(s, p, E, ET,M,AFM, ETOTAL, ETTOTAL,MTOTAL, MTOTAL2, MTOTAL4,AFMTOTAL, BETA, NX::Int, NY::Int, NZ::Int, NMC, history) # s = initial con
	  (s, p, siteEnergy,
      E,ET,M,AFM,ETOTAL, E2TOTAL,MTOTAL, MTOTAL2, MTOTAL4, AFMTOTAL,chiTOTAL,totalClustersTOTAL, O3TOTAL,O4TOTAL,OMeasurementCounter,chiMeasurementCounter, measurementCounter,
      acceptance, Mhistory)=mcRUNCLUSTERNT(s, p, siteEnergy,
       E, chi, M,AFM,ETOTAL, E2TOTAL, MTOTAL, MTOTAL2, MTOTAL4, AFMTOTAL, chiTOTAL,totalClustersTOTAL, O3TOTAL,O4TOTAL, BETA, gs,NX, NY, NZ, NMC);
	  EAV = ETOTAL/measurementCounter;	MAV = MTOTAL/measurementCounter;	M2AV = MTOTAL2/measurementCounter; M4AV = MTOTAL4/measurementCounter;
	  E2AV = E2TOTAL/measurementCounter;
	  # this is the TOTAL euler character of all components.
#	  ETAV= ETTOTAL/measurementCounter;
	  AFMAV = AFMTOTAL/measurementCounter;
	  chiAV = chiTOTAL/chiMeasurementCounter;
	  totalClustersAV = totalClustersTOTAL/chiMeasurementCounter;
	  O3 = O3TOTAL/OMeasurementCounter
	  O4 = O4TOTAL/OMeasurementCounter
	  binder = 0.5*(3 - M4AV/M2AV^2);

	  ## BINNING
	  	M = length(Mhistory)
	  	number_of_decimations=(convert(Int64,  floor(log(2,M)))-5)
	  	#DeltaCollection = zeros(number_of_decimations)
	  	binnedError = binning(Mhistory, number_of_decimations)



	  return (binder, MAV, EAV, M2AV, M4AV, AFMAV, chiAV, totalClustersAV, O3, O4, acceptance/NMC,E2AV, binnedError)
end




function runTheNumbersBP(j::Int64, NX::Int,NY::Int,NZ::Int, BETA, gs, NMC)
	#function runTheNumbers(NX,NY, Ts, NMC, j::Int64)
#	  BETA = 1/T;
#	  gF = log(gs)
#	  gE = -log(gs)
#	  gV = log(gs)
s = initSpinsHot(NX,NY,NZ);  # initial spin config.
p = domainWallConfig(s, NX,NY,NZ)
(nubbinData,siteEnergy) = getNubbinData(p,NX,NY,NZ)

## here we must get nubbinData and siteEnergies.

	(E,chi,M,AFM) = energyIsingBP(s, p,nubbinData, siteEnergy,NX,NY,NZ) ;
	#	  history = [copy(s)]
	ETOTAL=0; E2TOTAL=0;MTOTAL=0; M2TOTAL=0; MTOTAL4=0;AFMTOTAL=0; chiTOTAL=0;chi2TOTAL=0;
	totalClustersTOTAL=0; O3TOTAL=0; O4TOTAL=0;
	#																													(s, p, E, ET,M,AFM, ETOTAL, ETTOTAL,MTOTAL, MTOTAL2, MTOTAL4,AFMTOTAL, BETA, NX::Int, NY::Int, NZ::Int, NMC, history) # s = initial con
	(s, p, nubbinData, siteEnergy,
	E,ET,M,AFM,ETOTAL, E2TOTAL,MTOTAL, M2TOTAL, MTOTAL4, AFMTOTAL,chiTOTAL,chi2TOTAL, totalClustersTOTAL, O3TOTAL,O4TOTAL,OMeasurementCounter, chiMeasurementCounter, measurementCounter,
	acceptance, Mhistory)=mcRUNCLUSTERBP(s, p, nubbinData, siteEnergy,
	 E, chi, M,AFM,ETOTAL, E2TOTAL,MTOTAL, M2TOTAL, MTOTAL4, AFMTOTAL, chiTOTAL,chi2TOTAL, totalClustersTOTAL,  O3TOTAL,O4TOTAL, BETA, gs,NX, NY, NZ, NMC);
	  EAV = ETOTAL/measurementCounter;	MAV = MTOTAL/measurementCounter;	M2AV = M2TOTAL/measurementCounter; M4AV = MTOTAL4/measurementCounter;
	  E2AV = E2TOTAL/measurementCounter;
	  # this is the TOTAL euler character of all components.
#	  ETAV= ETTOTAL/measurementCounter;
	  AFMAV = AFMTOTAL/measurementCounter;
	  chiAV = chiTOTAL/chiMeasurementCounter;
	  chi2AV = chi2TOTAL/chiMeasurementCounter;

	  totalClustersAV = totalClustersTOTAL/chiMeasurementCounter;
	  O3 = O3TOTAL/OMeasurementCounter
	  O4 = O4TOTAL/OMeasurementCounter
	  binder = 0.5*(3 - M4AV/M2AV^2);

	  ## BINNING
	  	M = length(Mhistory)
	  	number_of_decimations=(convert(Int64,  floor(log(2,M)))-5)
	  	#DeltaCollection = zeros(number_of_decimations)
	  	binnedError = binning(Mhistory, number_of_decimations)



	  return (binder, MAV, EAV, M2AV, M4AV, AFMAV, chiAV, totalClustersAV, O3, O4, acceptance/NMC,E2AV, chi2AV, binnedError)
end
















#
