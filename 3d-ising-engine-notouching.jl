# this is the engine with all the functions in it -- from now on, I think it is
# good to just edit this one.

# this is the 3d topological Ising model, v2.0: it now includes the branch points
# when the resolution of the face reconnections are not done correctly.
# it is however extremely inefficient. now try to get this to work faster -- update energy
# by only updating each individual spin rather than all of them.

# this is faster, but still very very slow compared to the usual 3d ising.

# this version has some extra stuff to count clusters.

# note the "otherDir" function here is defined differently (preserves cyclic permutations
#); I think it makes no difference to anything.

# finally this version has been edited to use the new more symmetric lookup table!

# okay and this version will also use the no-touching algorithm, i.e. it will exclude
# any configuration that has an edge number equal to 2.


using DelimitedFiles
using Plots
using Statistics
using LinearAlgebra
using Dates


# this is a configuration of walls and spins, store the sizes of the config here
# as well.
struct config
    spins::AbstractArray{Int8,3}
    walls::AbstractArray{Int8,4}
    NX::Int16
    NY::Int16
    NZ::Int16
	# constructor just makes a new thing of all spins up, no domain walls.
    config(NX,NY,NZ) = new(ones(Int8, NX, NY,NZ),zeros(Int8, NX, NY,NZ,3),NX,NY,NZ);
	# note that the wall object is labeled by a point in the dual lattice plus a normal direction
end

# this thing flips a spin on a config (and updates all the important walls!)
function spinSet(c::config,i,j,k,new_val)
	# flip the spin!

    c.spins[i,j,k] = new_val;

    # next, write the code to update the wall data. there are 6 walls total to update.
	# the wall array is parametrized by a point in the dual lattice as well as a direction
	for dir=1:3
		# now first, there are three wall sites that have the same index.
		shifted = shiftAdd(i,j,k,dir,c.NX,c.NY,c.NZ);
		c.walls[i,j,k,dir] = (1-c.spins[i,j,k]*c.spins[shifted[1],shifted[2],shifted[3]])/2;

		# now there are three wall sites that differ by 1.
		shifted = shiftSubtract(i,j,k,dir,c.NX,c.NY,c.NZ);
		c.walls[shifted[1],shifted[2],shifted[3],dir]=(1-c.spins[i,j,k]*c.spins[shifted[1],shifted[2],shifted[3]])/2;
	end


	return c
end


# this adds 1 mod N
function periodicAdd(i,N)
    if i < N
        i = i +1;
    else
        i = 1;
    end
    return i;
end

# subtracts 1 mod N
function periodicSubtract(i, N)
    if i > 1
        i = i - 1;
    else
        i = N;
    end
    return i;
end

# this is a copy of john's shift code -- a difference is that my a runs over 1, 2, 3
# as I wanted this to agree more with the usual array index conventions in julia.
function shiftAdd(i,j,k,a, NX,NY,NZ)
	return a==1 ? [periodicAdd(i,NX),j,k] : ( a==2 ? [i,periodicAdd(j,NY),k] : [i,j,periodicAdd(k,NZ)] )
end

function shiftSubtract(i,j,k,a, NX,NY,NZ)
	return a==1 ? [periodicSubtract(i,NX),j,k] : ( a==2 ? [i,periodicSubtract(j,NY),k] : [i,j,periodicSubtract(k,NZ)] )
end

function isingEnergy(c::config)
    # This actually simply returns the number of: walls, edges, vertices, the total Euler character and the total magnetization.
    walls = 0;
	vertices = 0;
	edges = 0;
    MM = 0;
    for i = 1:c.NX
        for j = 1:c.NY
            for k = 1:c.NZ

				# simply count the number of walls and edges to get the energy
                for dir = 1:3
					walls += c.walls[i,j,k,dir];
					# let me modify it to use the new edge counting
					edges += new_edge_notouch(c,i,j,k,dir);
					#edges += edge(c,i,j,k,dir);
				end

				# alco count the number of vertices
				vertices += vertex(c,i,j,k);

				# add the spins themselves
                MM += c.spins[i,j,k];
            end
        end
    end

	# return everything -- the last line is the total Euler character!
    return (walls,edges,vertices,MM);
end

function otherDirs(dir)
	# given a direction, returns the other two directions
	# see if its okay to flip this
	return dir==1 ? [2,3] : ( dir==2 ? [3,1] : [1,2] )
end

function edge(c,i,j,k,dir)
	# returns the effective number of edges at that point and direction, from the wall configs
	# note that i,j,k parametrize a point in the dual lattice, and dir is a direction off that point.
	# the conventions are chosen so that the three edges at the point (i,j,k) are each in between the
	# walls that are stored at (i,j,k) (i.e. edge(i,j,k,1) is in between wall(i,j,k,2) and wall(i,j,k,3))

	od = otherDirs(dir);

	# the two shifted vertex points
	shift1 = shiftAdd(i,j,k,od[1],c.NX,c.NY,c.NZ);
	shift2 = shiftAdd(i,j,k,od[2],c.NY,c.NY,c.NZ);

	w = [c.walls[i,j,k,od[1]], c.walls[i,j,k,od[2]],c.walls[shift1[1],shift1[2],shift1[3],od[2]],c.walls[shift2[1],shift2[2],shift2[3],od[1]]];


	# the number of possible active walls is zero, 2, or 4. (anything else is a wall with end)
	# thus a simple mapping to number of edges is:
	return convert(Int16,sum(w)/2);
end

function new_edge(c,i,j,k,dir)
	# this is a new edge function; it uses the vertex reconnection idea and includes branch points.
	# this is the one that should be used!
	# note, this direction is (x =1, y = 2, z = 3)

	# note: the edge sits between the vertex(i,j,k) and (i-1, j, k) (or so on, depending on dir)

	odir1 = dir*2; # this is for the vertex at (i,j,k)
	odir2 = dir*2 - 1; # this is for the vertex at (i-1, j, k)
	shift2 = shiftSubtract(i,j,k,dir,c.NX,c.NY,c.NZ); #these are the coordinates (i-1, j, k)

	nubbin1 = vertex_connection(c,i,j,k,odir1);
	nubbin2 = vertex_connection(c,shift2[1],shift2[2],shift2[3],odir2);

	# now, if either of them is 0 or 1 then its clear how this is working.
	# if they disagree then we are introducing a branch point, which is like an edge
	# of weight 3
	# if they agree and they're not 0 or 1 then we should return 2
	return (nubbin1 == 0 ? 0 : (nubbin1 == 1 ? 1 : (nubbin1 != nubbin2 ? 3 : 2)));
end

function new_edge_notouch(c,i,j,k,dir)
	# this is a new edge function; it returns a 1000 for the edge number if there are two
	# domain walls touching at an edge.

	# note, this direction is (x =1, y = 2, z = 3)

	# note: the edge sits between the vertex(i,j,k) and (i-1, j, k) (or so on, depending on dir)

	odir1 = dir*2; # this is for the vertex at (i,j,k)
	odir2 = dir*2 - 1; # this is for the vertex at (i-1, j, k)
	shift2 = shiftSubtract(i,j,k,dir,c.NX,c.NY,c.NZ); #these are the coordinates (i-1, j, k)

	nubbin1 = vertex_connection(c,i,j,k,odir1);
	nubbin2 = vertex_connection(c,shift2[1],shift2[2],shift2[3],odir2);

	# now, if either of them is 0 or 1 then its clear how this is working.
	# if they disagree then we are introducing a branch point, which is like an edge
	# of weight 3
	# if they agree and they're not 0 or 1 then we should return 2
	return (nubbin1 == 0 ? 0 : (nubbin1 == 1 ? 1 : 1000));
end



# shamelessly copied from john:
function fromDigits(p)
    len=length(p);
    sum([p[k]*2^(len-k) for k=1:len]);
end

function basisChange(q)
	return convert(Int8,(1-q)/2);
end

function vertex(c,i,j,k)
	# this is a bit complicated. okay, here we go. Sites is a bit string.
	sites = zeros(Int8,8);

	# ugly but I find this the easiest way to remember which site is associated
	# with which point. basis change just maps (-1,1) -> (0,1).
	# (worry about overheads in function calling? eh, not for now.)

	sites[1]=basisChange(c.spins[i,j,k]);
	sites[2]=basisChange(c.spins[i,j,periodicAdd(k,c.NZ)]);
	sites[3]=basisChange(c.spins[i,periodicAdd(j,c.NY),k]);
	sites[4]=basisChange(c.spins[i,periodicAdd(j,c.NX),periodicAdd(k,c.NY)]);
	sites[5]=basisChange(c.spins[periodicAdd(i,c.NX),j,k]);
	sites[6]=basisChange(c.spins[periodicAdd(i,c.NX),j,periodicAdd(k,c.NZ)]);
	sites[7]=basisChange(c.spins[periodicAdd(i,c.NX),periodicAdd(j,c.NY),k]);
	sites[8]=basisChange(c.spins[periodicAdd(i,c.NX),periodicAdd(j,c.NY),periodicAdd(k,c.NZ)]);

	# now just return the thing from the previously computed lookup table!
	return lookup[fromDigits(sites)+1];

end

function vertex_debug(c,i,j,k)
	# this is the same as the vertex function, but it also returns the config #
	sites = zeros(Int8,8);

	# ugly but I find this the easiest way to remember which site is associated
	# with which point. basis change just maps (-1,1) -> (0,1).
	# (worry about overheads in function calling? eh, not for now.)

	sites[1]=basisChange(c.spins[i,j,k]);
	sites[2]=basisChange(c.spins[i,j,periodicAdd(k,c.NZ)]);
	sites[3]=basisChange(c.spins[i,periodicAdd(j,c.NY),k]);
	sites[4]=basisChange(c.spins[i,periodicAdd(j,c.NX),periodicAdd(k,c.NY)]);
	sites[5]=basisChange(c.spins[periodicAdd(i,c.NX),j,k]);
	sites[6]=basisChange(c.spins[periodicAdd(i,c.NX),j,periodicAdd(k,c.NZ)]);
	sites[7]=basisChange(c.spins[periodicAdd(i,c.NX),periodicAdd(j,c.NY),k]);
	sites[8]=basisChange(c.spins[periodicAdd(i,c.NX),periodicAdd(j,c.NY),periodicAdd(k,c.NZ)]);

	# now just return the thing from the previously computed lookup table!
	return (lookup[fromDigits(sites)+1],fromDigits(sites)+1);

end

function vertex_connection(c,i,j,k,ordir)
	# this is a function which returns the way in which the faces are connected
	# off of a particular "oriented direction" at vertex (i,j,k).
	# recall: ordir: (1,2,3,4,5,6) = (x+, x-, y+, y-, z+, z-)

	# this is just copied from above: should make this more efficient; store it, I think!
	# this is a bit complicated. okay, here we go. Sites is a bit string.
	sites = zeros(Int8,8);


	sites[1]=basisChange(c.spins[i,j,k]);
	sites[2]=basisChange(c.spins[i,j,periodicAdd(k,c.NZ)]);
	sites[3]=basisChange(c.spins[i,periodicAdd(j,c.NY),k]);
	sites[4]=basisChange(c.spins[i,periodicAdd(j,c.NX),periodicAdd(k,c.NY)]);
	sites[5]=basisChange(c.spins[periodicAdd(i,c.NX),j,k]);
	sites[6]=basisChange(c.spins[periodicAdd(i,c.NX),j,periodicAdd(k,c.NZ)]);
	sites[7]=basisChange(c.spins[periodicAdd(i,c.NX),periodicAdd(j,c.NY),k]);
	sites[8]=basisChange(c.spins[periodicAdd(i,c.NX),periodicAdd(j,c.NY),periodicAdd(k,c.NZ)]);

	# now just return it.
	return lookup_connections[fromDigits(sites)+1,ordir];
end


# randomizes the spins
function randomizeSpins(c::config)

    for i = 1:c.NX
        for j = 1:c.NY
            for k = 1:c.NZ
                spinSet(c,i,j,k,rand([-1 1]));
            end
        end
    end
end

function mcStepSlow(c::config,EE,MM, beta, phi) # phi is log of string coupling (i.e. dilaton).
	# old fashioned slow code (sums over all sites at each step)
    # this is code to perform the Monte Carlo step! Wat leuk

    # pick a site at random
    x = rand(1:c.NX);
    y = rand(1:c.NY);
    z = rand(1:c.NZ);

    # now construct the energy difference from flipping the spin
    # for now, do this in a terrible manner -- compute the full thing.


	# flip the spin!
    spinSet(c,x,y,z,-c.spins[x,y,z]);
	(EEnew, MMnew) = boltzWeight(c,beta,phi);

	DeltaE = EEnew - EE;
	#println("Delta E is: ", DeltaE);


    Ups = exp(-DeltaE);


    if (rand() > Ups)
		# I should not accept it -- so, unflip the spin and undo the assignments
		spinSet(c,x,y,z,-c.spins[x,y,z]);
		EEnew = EE;
		MMnew = MM;
		#println("Rejected");
    end

	return(EEnew,MMnew)
	# note that the spin configuration is updated automatically through the miracle
	# of programming, so I don't pass it back. I like this better as EE contains the ugly
	# physicsy things beta and phi and not just beautiful 1s and 0s.
end

function mcStep(c::config,EE,MM, beta, phi) # more efficient mc step

    # pick a site at random
    x = rand(1:c.NX);
    y = rand(1:c.NY);
    z = rand(1:c.NZ);

    # now find the local geometry about the point, first before:

	(walls_i, edges_i, vertices_i) = new_localGeometry(c,x,y,z);



	# flip spin:
    spinSet(c,x,y,z,-c.spins[x,y,z]);

	# still a diagnostic:

	#(EEnew2, MMnew2) = boltzWeight(c,beta,phi);
	# then after:
	(walls_f, edges_f, vertices_f) = new_localGeometry(c,x,y,z);

	DeltaE = 2*beta*(walls_f - walls_i) + phi*(walls_f-walls_i-edges_f+edges_i+vertices_f-vertices_i);
	MMnew = MM + 2*c.spins[x,y,z];

	EEnew = EE + DeltaE;


#=	if (abs(EEnew-EEnew2)>0.001)
		println("DISASTER DISASTER!");
		println(EEnew, " ",EEnew2, " ");
		println(abs(EEnew-EEnew2));
		println("Hello");
		quit();
	end=#

	# now compare:
    Ups = exp(-DeltaE);

	rejected = false;

	# note this new thing:
	# reject it if the final config is forbidden
	# reject it energetically only if the *initial* config is allowed
	# thus it should really relax away from the initial state.
    if (edges_f > 500 || (rand() > Ups  && edges_i < 500))
		# I should not accept it -- so, unflip the spin and undo the assignments
		spinSet(c,x,y,z,-c.spins[x,y,z]);
		EEnew = EE;
		MMnew = MM;
		rejected = true;
		#println("Rejected");
    end

	return(EEnew,MMnew,rejected)
	# note that the spin configuration is updated automatically through the miracle
	# of programming, so I don't pass it back. I like this better as EE contains the ugly
	# physicsy things beta and phi and not just beautiful 1s and 0s.
end

function localGeometry(c::config,x,y,z)
	# returns the number of walls, edges, and vertices surrounding a particular point
	walls = 0;
	for dir=1:3
		# now first, there are three wall sites that have the same index.
		shifted = shiftAdd(x,y,z,dir,c.NX,c.NY,c.NZ);
		walls += c.walls[x,y,z,dir];

		# now there are three wall sites that differ by 1.
		shifted = shiftSubtract(x,y,z,dir,c.NX,c.NY,c.NZ);
		walls += c.walls[shifted[1],shifted[2],shifted[3],dir];
	end

	edges = 0;

	# need to finish writing the edge and vertex code.
	for dir = 1:3
		# first the easy 3 that naturally bound the thing.
		edges += new_edge_notouch(c,x,y,z,dir);

		# next, the other 3 -- move in the opposite direction to find them!

		od = otherDirs(dir);
		shifted_1 = shiftSubtract(x,y,z,od[1],c.NX,c.NY,c.NZ);
		shifted_2 = shiftSubtract(x,y,z,od[2],c.NX,c.NY,c.NZ);
		shifted_both = shiftSubtract(shifted_1[1],shifted_1[2],shifted_1[3],od[2],c.NX,c.NY,c.NZ);
		edges+=new_edge_notouch(c,shifted_1[1],shifted_1[2],shifted_1[3],dir);
		edges+=new_edge_notouch(c,shifted_2[1],shifted_2[2],shifted_2[3],dir);
		edges+=new_edge_notouch(c,shifted_both[1],shifted_both[2],shifted_both[3],dir);
	end

	# finally, need to find the vertices.
	vertices = 0;

	for k=0:7
		shifted = [x,y,z];
		# make a binary number telling me whether or not to shift each one.
		to_shift = digits(k,base=2,pad=3);
		for dir = 1:3
			if (to_shift[dir]==1)
				shifted = shiftSubtract(shifted[1],shifted[2],shifted[3],dir,c.NX,c.NY,c.NZ);
			end
		end
		# okay, now I have the point. now just add the vertex contribution:
		vertices += vertex(c,shifted[1],shifted[2],shifted[3]);
	end

	# send it all back.
	return (walls, edges, vertices);
end

function new_localGeometry(c::config,x,y,z)
	# returns the number of walls, edges, and vertices surrounding a particular point
	# this one actually does it for all edges that can be affected by this -- with the new
	# branch point code this is quite a few (36 altogether). walls and vertices are the same.
	walls = 0;
	for dir=1:3
		# now first, there are three wall sites that have the same index.
		shifted = shiftAdd(x,y,z,dir,c.NX,c.NY,c.NZ);
		walls += c.walls[x,y,z,dir];

		# now there are three wall sites that differ by 1.
		shifted = shiftSubtract(x,y,z,dir,c.NX,c.NY,c.NZ);
		walls += c.walls[shifted[1],shifted[2],shifted[3],dir];
	end

	edges = 0;

	# new edge code: checks all 36 edges that can possibly be affected by the flipping
	# of a spin.
	for dir = 1:3



		# the direction specifies which way the edge points. there are 4 different axes,
		# generate these by moving in the 4 directions transverse to the edge:

		od = otherDirs(dir);
		shifted_1 = shiftSubtract(x,y,z,od[1],c.NX,c.NY,c.NZ);
		shifted_2 = shiftSubtract(x,y,z,od[2],c.NX,c.NY,c.NZ);
		shifted_both = shiftSubtract(shifted_1[1],shifted_1[2],shifted_1[3],od[2],c.NX,c.NY,c.NZ);

		# form all the axes into an array .
		allaxes = [[x,y,z],shifted_1,shifted_2,shifted_both];

		for n=1:4
			# now each axis contributes either + or -
			# the middle one
			edges += new_edge_notouch(c,allaxes[n][1],allaxes[n][2],allaxes[n][3],dir);

			# the ones off the edge
			shift_plus = shiftAdd(allaxes[n][1],allaxes[n][2],allaxes[n][3],dir,c.NX,c.NY,c.NZ);
			edges += new_edge_notouch(c,shift_plus[1],shift_plus[2],shift_plus[3],dir)

			shift_minus = shiftSubtract(allaxes[n][1],allaxes[n][2],allaxes[n][3],dir,c.NX,c.NY,c.NZ);
			edges += new_edge_notouch(c,shift_minus[1],shift_minus[2],shift_minus[3],dir)
		end



	end

	# finally, need to find the vertices.
	vertices = 0;

	for k=0:7
		shifted = [x,y,z];
		# make a binary number telling me whether or not to shift each one.
		to_shift = digits(k,base=2,pad=3);
		for dir = 1:3
			if (to_shift[dir]==1)
				shifted = shiftSubtract(shifted[1],shifted[2],shifted[3],dir,c.NX,c.NY,c.NZ);
			end
		end
		# okay, now I have the point. now just add the vertex contribution:
		vertices += vertex(c,shifted[1],shifted[2],shifted[3]);
	end

	# send it all back.
	return (walls, edges, vertices);
end

# computes the Boltzman weight and the  magnetization.
function boltzWeight(c::config,beta,phi)
	(walls, edges, vertices, MM) = isingEnergy(c);
	EE = 2*beta*walls + phi*(walls-edges+vertices);
	return (convert(Float32,EE),convert(Float32,MM));
end

function staggeredM(c::config)
	sMM = 0;
	for i = 1:c.NX
		for j = 1:c.NY
			for k =1:c.NZ
				sMM += (-1)^(i+j+k)*c.spins[i,j,k];
			end
		end
	end
	return convert(Float32,sMM);
end

#=
function dropletM(c::config)
	dMM = zeros(Float32,4);
	# the order parameter was constructed in the mathematica code.
	for i = 1:c.NX
		for j = 1:c.NY
			for k =1:c.NZ
				#println(i, " ", j, " ", k, " ");
				dMM[min(fromDigits([i%2 j%2 k%2])+1, 8-fromDigits([i%2 j%2 k%2]))]+=(c.spins[i,j,k])/2;
			end
		end
	end
	# send back both the order parameter vector and the "magnitude" that is invariant.
	return (dMM, dMM[1]-dMM[2]-dMM[3]+dMM[4]);
end=#

# debugging function that immediately returns chi:
function chi(c::config)
	(walls,edges,vertices,MM) = isingEnergy(c);
	return walls - edges + vertices;
end


# this is an older version of mcRun
function mcRunOld(c::config,beta,phi,Nsteps,sample_rate)


    # initialize the energy etc. from the thing.
	(EE,MM) = boltzWeight(c,beta,phi);

	# initialize the Marray and Earray
    sites = convert(Float32,c.NX*c.NY*c.NZ);
    Marray = [MM/sites];
    Earray = [EE/sites];
	sMarray = [staggeredM(c)/sites];

	(o3, o4) = orderParameters(c);
	o3array = [o3/(sites*sites)];
	o4array = [o4/(sites*sites)];



	chiArray =[chi(c)];
	clustersArray = [length(findClusters(c)[2])];

	rejected_steps = 0;

    for i=1:Nsteps
		#println("Hello", i);
        (EE,MM,rejected) = mcStep(c,EE,MM,beta,phi);

		if rejected
			rejected_steps +=1
		end
        #Etotal += EE/(NX*NY);
        #Mtotal += MM/(NX*NY);

        # store the magnetization and energy every few steps.
        if(mod(i,sample_rate)==0)
            Marray = [Marray;MM^2/(sites*sites)];
            Earray = [Earray;EE/sites];
			sMarray = [sMarray;staggeredM(c)^2/(sites*sites)];

			# compute the funny order parameters
			(o3, o4) = orderParameters(c);
			o3array = [o3array; o3/(sites*sites)];
			o4array = [o4array; o4/(sites*sites)];

			chiArray = [chiArray; chi(c)/sites];
			clustersArray = [clustersArray; length(findClusters(c)[2])];
			println("Sampling ", i);
        end
    end
    return (Marray,Earray,sMarray,o3array,o4array,chiArray,clustersArray,convert(Float32,rejected_steps/Nsteps));
end


function mcInit(c::config,beta,phi,Nsteps)
	# this just runs a few steps on c to initialize; do this at the beginning
	(EE,MM) = boltzWeight(c,beta,phi);
	for i=1:Nsteps
		#println("Hello", i);
        (EE,MM,rejected) = mcStep(c,EE,MM,beta,phi);
		if (mod(i,1000)==0)
			println("At beta ", beta, " phi ", phi, " Equilibrating ", (Float32)(i/Nsteps));
		end
	end
end
# this is a new version of mcRun that is (in principle) more efficient, and returns the autocorrelation
# it also in principle lets you compute the errors.


function mcRunErrors(c::config,beta,phi,Nsteps,sample_rate)


    # initialize the energy etc. from the thing.
	(EE,MM) = boltzWeight(c,beta,phi);

	# initialize the Marray and Earray
	numSamples = convert(Int32,floor(Nsteps/sample_rate))+1;
	whichSample = 1;
	sites = convert(Float32,c.NX*c.NY*c.NZ);
	Marray = zeros(numSamples);
	Earray = zeros(numSamples);
	sMarray = zeros(numSamples);
	o3array = zeros(numSamples);
	o4array = zeros(numSamples);
	chiArray = zeros(numSamples);

	# in this version let me not keep track of clusters. delete that.
    Marray[1] = MM/sites;
    Earray[1] = EE/sites;
	sMarray[1] = staggeredM(c)/sites;

	(o3, o4) = orderParameters(c);
	o3array[1] = o3/(sites*sites);
	o4array[1] = o4/(sites*sites);



	chiArray[1] =chi(c);
	#clustersArray = [length(findClusters(c)[2])];

	rejected_steps = 0;

    for i=1:Nsteps
		#println("Hello", i);
        (EE,MM,rejected) = mcStep(c,EE,MM,beta,phi);

		if rejected
			rejected_steps +=1;
		end
        #Etotal += EE/(NX*NY);
        #Mtotal += MM/(NX*NY);

        # store the magnetization and energy every few steps.
        if(mod(i,sample_rate)==0)
			whichSample += 1;
            Marray[whichSample] = MM/sites;
            Earray[whichSample] = EE/sites;
			sMarray[whichSample] = staggeredM(c)/sites;

			# compute the funny order parameters
			(o3, o4) = orderParameters(c);
			o3array[whichSample] = o3/(sites*sites);
			o4array[whichSample] = o4/(sites*sites);

			chiArray[whichSample] = chi(c)/sites;
			#clustersArray = [clustersArray; length(findClusters(c)[2])];
			println("Sampling ", i/Nsteps, " MM ", Marray[whichSample]);
        end
    end

	# now, let's compute the autocorrelation. The first thing to do is take only the data
	# from a certain ratio of things at the *end* (so after its correlated)

	# okay, I think now the point is -- we drop a certain fraction. let's do it;
	# start with half.

	dropFrac = 0.5;
	start = convert(Int32,floor(numSamples*dropFrac));

	# lowest value it can have is 0.
	usedSamples = numSamples - start;
	autocorrM = zeros(usedSamples);
	for t=1:usedSamples
		# this is using the formula in Newman+Barkema, p63
		Mbar1 = sum(Marray[1+start:(numSamples-t+1)])/(convert(Float64,usedSamples-t+1));
		Mbar2 = sum(Marray[t+start:numSamples])/(convert(Float64,usedSamples-t+1));
		for tprime = 1:(usedSamples-t+1) # so that I never go outside the domain
			autocorrM[t] += ((Marray[start+tprime]*Marray[start+tprime+t-1])/convert(Float64,usedSamples-t+1));
		end
		autocorrM[t] -= Mbar1*Mbar2;
		#println(Mbar1);
	end

	# now determine integrated correlation time
	#tau = sum(autocorrM)/autocorrM[1];
	#println(tau);

    return (Marray,Earray,sMarray,o3array,o4array,chiArray,convert(Float32,rejected_steps/Nsteps),autocorrM);
end


# this is going to be the workhorse version.
# try to make it as useful and efficient as possible
function mcRunEfficient(c::config,beta,phi,Nsteps,sample_rate)


    # initialize the energy etc. from the thing.
	(EE,MM) = boltzWeight(c,beta,phi);

	# initialize the Marray and Earray
	numSamples = convert(Int32,floor(Nsteps/sample_rate))+1;
	whichSample = 1;
	sites = convert(Float32,c.NX*c.NY*c.NZ);
	Marray = zeros(numSamples);
	Earray = zeros(numSamples);
	sMarray = zeros(numSamples);
	o3array = zeros(numSamples);
	o4array = zeros(numSamples);
	#chiArray = zeros(numSamples);

	# in this version let me not keep track of clusters. delete that. Note
	# that each array contains a *square-rooted* version.
    Marray[1] = abs(MM)/sites;
    Earray[1] = EE/sites;
	sMarray[1] = abs(staggeredM(c))/sites;

	(o3, o4) = orderParameters(c);
	o3array[1] = o3/sites;
	o4array[1] = o4/sites;



	#chiArray[1] =chi(c);
	#clustersArray = [length(findClusters(c)[2])];

	rejected_steps = 0;

    for i=1:Nsteps
		#println("Hello", i);
        (EE,MM,rejected) = mcStep(c,EE,MM,beta,phi);

		if rejected
			rejected_steps +=1;
		end
        #Etotal += EE/(NX*NY);
        #Mtotal += MM/(NX*NY);

        # store the magnetization and energy every few steps.
        if(mod(i,sample_rate)==0)
			whichSample += 1;
            Marray[whichSample] = abs(MM)/sites;
            Earray[whichSample] = EE/sites;
			sMarray[whichSample] = abs(staggeredM(c))/sites;

			# compute the funny order parameters
			(o3, o4) = orderParameters(c);
			o3array[whichSample] = o3/sites;
			o4array[whichSample] = o4/sites;

			# okay, let me just not compute the chis for now (since I don't use that!)
			#chiArray[whichSample] = chi(c);
			#clustersArray = [clustersArray; length(findClusters(c)[2])];
			#println("Sampling ", i/Nsteps, " MM ", Marray[whichSample]);
        end
		# output a progress bar every couple of steps
		if (mod(i,4000)==0)
			println("Sampling ", i/Nsteps);
		end
    end


	# return the full arrays. don't do any computations here.
    return (Marray,Earray,sMarray,o3array,o4array,convert(Float32,rejected_steps/Nsteps));
end

function outputData(c::config)
	# this outputs all the data. also cluster decompositions.

	(wallCluster, uniqueClusters) = findClusters(s);

	println("Number of clusters: ", length(uniqueClusters));

	# keep track of walls
	to_store_walls = zeros(Int16,0,5);

	# similarly, edges.
	to_store_edges = zeros(Int16,0,5);

	# similarly, vertices.
	to_store_vertices = zeros(Int16,0,5);
	# find the nonzero ones
	for i=1:c.NX
		for j=1:c.NY
			for k =1:c.NZ

				# for each point, check the vertices:
				(v1, val1) = vertex_debug(s,i,j,k)

				if (v1!=0)
					to_store_vertices = [to_store_vertices; i j k v1 val1];
				end

				for dir=1:3

					# first the walls
					if (c.walls[i,j,k,dir] != 0)
						to_store_walls = [to_store_walls; i j k dir wallCluster[i,j,k,dir]];
					end

					# then the edges
					e1 = new_edge_notouch(c,i,j,k,dir)

					if (e1 != 0 )
						to_store_edges = [to_store_edges; i j k dir e1];
					end

				end

			end
		end
	end


	open(string(pwd(),"/mathematica/new code (summer 2019)/walls.txt"), "w") do io
		writedlm(io, to_store_walls);
	end

	writedlm(string(pwd(),"/mathematica/new code (summer 2019)/edges.txt"), to_store_edges);
	writedlm(string(pwd(),"/mathematica/new code (summer 2019)/vertices.txt"), to_store_vertices);
end

# I will need union find etc. to run the cluster finding algorithm. let's get this to
# work first

# this is a function that returns the four walls that are adjacent to an edge.

function adjoinedWalls(c::config, i,j,k,dir)
	od = otherDirs(dir);

	# so the way to do this, is to first use this bit copied from the new_edge function
	odir1 = dir*2; # this is for the vertex at (i,j,k)
	odir2 = dir*2 - 1; # this is for the vertex at (i-1, j, k)
	shift_a = shiftSubtract(i,j,k,dir,c.NX,c.NY,c.NZ); #these are the coordinates (i-1, j, k)

	# these two are the nubbins on each vertex on the two sides of the edge.
	nubbins = zeros(Int8,2);
	nubbins[1] = vertex_connection(c,i,j,k,odir1);
	nubbins[2] = vertex_connection(c,shift_a[1],shift_a[2],shift_a[3],odir2);

	# now build the table of connections: note this thing copied from the emails
	#= tells me what the nubbin does
	0: no faces
	1: one DW passing through, no ambiguity
	2: connections are: 1-2, 3-4
	3: connections are 1-3, 2-4
	4: connections are 1-4, 2-3.  =#

	shift1 = shiftAdd(i,j,k,od[1],c.NX,c.NY,c.NZ);
	shift2 = shiftAdd(i,j,k,od[2],c.NX,c.NY,c.NZ);

	walls_out = [shift1' od[2] c.walls[shift1[1],shift1[2],shift1[3],od[2]];
		shift2' od[1] c.walls[shift2[1],shift2[2],shift2[3],od[1]];
		i j k od[2] c.walls[i,j,k,od[2]];
		i j k od[1] c.walls[i,j,k,od[1]]];
	# okay, now then just return the right wall info.
	# there is a mathematica notebook that visualizes this and it is indeed correct.

	connections = zeros(Int8,4,4);

	# now we build the connections table. this takes into account the info from above.
	# I tried to work hard to get the walls conventions etc. to line up, so this m
	# makes it not so difficult to populate the connections table. note that if
	# there is a branch point it really populates the convention table with both
	# sets of connections.
	for i = 1:2
		if nubbins[i]==1
			# then this is a single DW. ironically we then have to work harder.
			for m = 1:4
				for n = 1:4
					if (walls_out[m,5]==1 && walls_out[n,5]==1)
						connections[m,n]=1;
					end
				end
			end
		elseif nubbins[i]==2
			connections[1,2]=1;
			connections[2,1]=1;
			connections[3,4]=1;
			connections[4,3]=1;
		elseif nubbins[i]==3
			connections[1,3]=1;
			connections[3,1]=1;
			connections[2,4]=1;
			connections[4,2]=1;
		elseif nubbins[i]==4
			connections[1,4]=1;
			connections[4,1]=1;
			connections[2,3]=1;
			connections[3,2]=1;
		end
		#println(nubbins[1]," ",nubbins[2]);

	end
	# returns the following:
	# a 4-element array containing the coords of the wall in question, its direction, and also whether
	# or not it is occupied.

	# also returns a 4x4 matrix telling whether each wall is connected to each other wall.

	return (walls_out, connections);
end

# this implementation is copied from: https://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html
function ufInitLabels(N)
	labels = zeros(Int32,N);
	for x = 1:N
		labels[x] = x;
	end
	return labels;
end

function ufFind(labels,x)
	while (labels[x]!=x)
		x = labels[x];
	end
	return x;
end

function ufUnion(labels,x,y)
	labels[ufFind(labels,x)]=ufFind(labels,y);
end

#function
# okay, now to decompose into clusters -- here we go!

function findClusters(c::config)
	# init the uf stuff.
	labels = ufInitLabels(c.NX*c.NY*c.NZ*3);
	wallCluster = zeros(Int64,c.NX,c.NY,c.NZ,3);
	maxCluster = 1;

	# loop through all edges
	for i = 1:c.NX
		for j = 1:c.NY
			for k = 1:c.NZ
				for dir = 1:3
					# get the info
					(tryWalls,connected) = adjoinedWalls(c,i,j,k,dir);
					for wall_i = 1:4
						# okay now here we go:
						if tryWalls[wall_i,5]==1
							# wall is occupied
							if wallCluster[tryWalls[wall_i,1],tryWalls[wall_i,2],tryWalls[wall_i,3],tryWalls[wall_i,4]]==0
								# previously unclustered
								wallCluster[tryWalls[wall_i,1],tryWalls[wall_i,2],tryWalls[wall_i,3],tryWalls[wall_i,4]]=maxCluster;
								maxCluster+=1;
							end

							# loop through the *other* walls, adjoining this, up to this one.
							for wall_j = 1:(wall_i-1)
								# also check that wall is occupied:
								if tryWalls[wall_j,5]==1

									if connected[wall_i,wall_j]==1
										# union the clusters
										ufUnion(labels,wallCluster[tryWalls[wall_i,1],tryWalls[wall_i,2],tryWalls[wall_i,3],tryWalls[wall_i,4]],
										wallCluster[tryWalls[wall_j,1],tryWalls[wall_j,2],tryWalls[wall_j,3],tryWalls[wall_j,4]]);

									end
								end
							end
						end
					end

				end
			end
		end
	end

	# okay, now that its done, make all of the clusters simpler
	# i.e. just run through and run find on everything to collapse the equivalence classes

	uniqueClusters = zeros(Int64,0);
	for i =1:c.NX
		for j = 1:c.NY
			for k =1:c.NZ
				for dir =1:3
					if wallCluster[i,j,k,dir]!=0
						wallCluster[i,j,k,dir]=ufFind(labels,wallCluster[i,j,k,dir]);

						# also keep track of the list of unique clusters so you can return it
						if (!(wallCluster[i,j,k,dir] in uniqueClusters))
							uniqueClusters = [uniqueClusters; wallCluster[i,j,k,dir]];
						end
					end

				end
			end
		end
	end
	return (wallCluster, uniqueClusters);
end


function latticeToCell(x,y,z)  ## give the label of the site in the 2x2x2 unit cell
    X = mod.([z,y,x],2)  # beware flipped convention!!
    return fromDigits(X)
end


function orderParameters(c::config)

	# this is taken from John's code; in principle, I should probably rewrite this myself to make sure its okay.
    #=v3=[ 1 0 0 -1 -1 0 0 1;
        0 1 0 -1 -1 0 1 0 ;
        0 0 1 -1 -1 1 0 0]

    v4= [-1 0 0 1 -1 0 0 1;
        0 -1 0 -1 1 0 1 0;
        0 0 -1 -1 1 1 0 0]
=#


    return (
	# this is AFM, as a check.
# 	sum(vAFM[latticeToCell(x,y,z)+1]*s[x,y,z] for x=1:NX, y=1:NY, z=1:NZ)^2,
    #This is O3:
	sqrt(sum(sum(v3[a,latticeToCell(x,y,z)+1]*c.spins[x,y,z] for x=1:c.NX, y=1:c.NY, z=1:c.NZ)^2 for a=1:3)),
    #This is O4:
    sqrt(sum(sum(v4[a,latticeToCell(x,y,z)+1]*c.spins[x,y,z] for x=1:c.NX, y=1:c.NY, z=1:c.NZ)^2 for a=1:3))
	)
end

function bootstrap(data,reps)
	# this figures out the errors in X via a bootstrap procedure. Following Newman+Barkema
	# Section 3.4.3.

	# reps is the number of resamples to do.
	estimates = zeros(reps);

	N = length(data);
	println(N);
	for k = 1:reps


		# randomly sample N times from this set of N things
		sampled = zeros(N);
		for i = 1:N
			sampled[i] = data[rand(1:N)];
		end

		# now calculate the mean of the sample
		estimates[k] = sum(sampled)/N;

	end
	# now I have a lot of things in this distribution. what do I do with them?
	mean = sum(estimates)/reps;
	meansq = sum(estimates.^2)/reps;

	# calculate the obvious mean and deviation as well
	obvmean = sum(data)/N;
	obvdev = sqrt((sum(data.^2)/N - obvmean^2)/convert(Float32,N-1));
	display(histogram(estimates));

	# return all of this stuff.
	return (sqrt(meansq - mean^2),obvdev,mean,obvmean);
end

function decimate(data)
	# this returns an array of half the size
	# with each entry being an average of two consecutive ones
	N = convert(Int32, floor(length(data)/2));
	decimated = zeros(N);
	for i = 1:N
		decimated[i] = (data[2*i-1]+data[2*i])/2;
	end
	return(decimated);
end


function meandev(data)
	# this just returns the mean and the deviation of the sample mean of the data
	return (mean(data), std(data)/sqrt(length(data)));
end

function binning(data)
	# this does a "binning" analysis. check how many times I can do it while keeping
	# 30 bins
	num_steps = convert(Int32,floor(log2(length(data)/30)));
	println(num_steps);
	devs = zeros(num_steps);
	for i = 1:num_steps
		data = decimate(data);
		devs[i]=meandev(data)[2]
		println(meandev(data)," ",devs[i]);
	end
	display(plot(devs));
	return devs;
end

# this function does a scan through phi values for a fixed beta

function phiScan(beta,phi_i,phi_f,nphis,sizeLat,Nsweeps)

	#sizeLat = 8;
	sample_rate = sizeLat*sizeLat*sizeLat;

	#Nsteps = 2000*sizeLat^3;
	# the below one is the number I usually use.
	#Nsteps = 1000000;

	#try again with twice as many and see if the answer is different!

	#Nsteps = 2000000

	# run for 5000 sweeps. this is pretty damn long for a 10x10x10 lattice.
	#Nsteps = 5000*(sizeLat*sizeLat*sizeLat);
	Nsteps = Nsweeps*(sizeLat*sizeLat*sizeLat);
	#Nsteps = 10000



	global s = config(sizeLat,sizeLat,sizeLat);
	#randomizeSpins(s);
	# let me try to initialize with spins up; not clear this is good, but let's do it.

	#println("Initial value of droplet order: ",dropletM(s)[2]/(s.NX*s.NY*s.NZ));
	# set up a table of phis:
	#phi_i=2;
	#phi_f=-2;
	#nphis =11;

	# making the number 11 makes each step rather nice.

	#phis = range(phi_i,stop=phi_f,length=nphis);
	phis = range(phi_i, stop=phi_f, length = nphis);
	#beta = 0.5;
	#beta = 0.22;

	#beta = 0;

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

	#mcInit(s,beta,phi_i,10*(sizeLat*sizeLat*sizeLat));
	# now if it doesnt exist, just construct a directory for this batch, labeled by today's date

	# store this in a variable, because it is possible for the date to change in the course
	# of the run!
	global datestart = today();

	if (!isdir(string(pwd(),"/data/",datestart)))
		mkdir(string(pwd(),"/data/",datestart))
	end

	#global filename = "";

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
		filename = string(pwd(), "/data/" , datestart, "/",datestart, "-L-", sizeLat, "-beta-", beta, "-phi-", phis[i]);

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

	filename_summary = string(pwd(), "/data/" , datestart, "/",datestart, "-aaa-L-", sizeLat, "-beta-", beta,"-phif-",phi_f);
	#output the summary to a file
	writedlm(string(filename_summary,"-summary.txt"),to_store);


	plot_sMs = plot(phis, sMs,xlabel="phi",ylabel="sM",legend=false);
	plot_Ms= plot(phis, Ms,xlabel="phi",ylabel="M",legend=false);
	plot_reject = plot(phis, rejectionRates,xlabel="phi",ylabel="rejection rate",legend=false);
	plot_o3 = plot(phis, o3s,xlabel="phi",ylabel="o3",legend=false);
	plot_o4 = plot(phis, o4s,xlabel="phi",ylabel="o4",legend=false);
	plot_Es = plot(phis, Es, xlabel="phi",ylabel="E",legend=false);

	plot(plot_sMs, plot_Ms, plot_o3, plot_o4, plot_reject, plot_Es, layout = (3,2))
	#title(string("L = ", sizeLat, " beta = ", beta));
	savefig(string(filename_summary,"-plots.pdf"));

end


string(pwd(),"/nabil code/better topological symmetric/lookuptable_spins_sym.txt")
const lookup = readdlm(string(pwd(),"/nabil code/better topological symmetric/lookuptable_spins_sym.txt"),Int8);
const lookup_connections = readdlm(string(pwd(),"/nabil code/better topological symmetric/lookuptable_connections_sym.txt"),Int8);

# new vectors taken from paper. v3 is called O+ in the paper
v3 = [1/2 0 0 -1/2 -1/2 0 0 1/2;
	-1/(2*sqrt(3)) 1/sqrt(3) 0 -1/(2*sqrt(3)) -1/(2*sqrt(3)) 0 1/sqrt(3) -1/(2*sqrt(3));
	-1/(2*sqrt(6)) -1/(2*sqrt(6)) sqrt(3/8) -1/(2*sqrt(6)) -1/(2*sqrt(6)) sqrt(3/8) -1/(2*sqrt(6)) -1/(2*sqrt(6))];
#vAFM = [1 -1 -1 1 -1 1 1 -1]

v4 = [-1/2 0 0 1/2 -1/2 0 0 1/2;
	-1/(2*sqrt(3)) -1/sqrt(3) 0 -1/(2*sqrt(3)) 1/(2*sqrt(3)) 0 1/sqrt(3) 1/(2*sqrt(3));
	-1/(2*sqrt(6)) 1/(2*sqrt(6)) -sqrt(3/8) -1/(2*sqrt(6)) 1/(2*sqrt(6)) sqrt(3/8) -1/(2*sqrt(6)) 1/(2*sqrt(6))];
#=
# lowest value it can have is 1.
usedSamples = numSamples - start + 1;
autocorrM = zeros(usedSamples);
for t=1:numSamples
	# this is using the formula in Newman+Barkema, p63
	Mbar1 = sum(Marray[1:(numSamples-t+1)])/(convert(Float64,numSamples-t+1));
	Mbar2 = sum(Marray[t:numSamples])/(convert(Float64,numSamples-t+1));
	for tprime = 1:(numSamples-t+1) # so that I never go outside the domain
		autocorrM[t] += ((Marray[tprime]*Marray[tprime+t-1])/convert(Float64,numSamples-t+1));
	end
	autocorrM[t] -= Mbar1*Mbar2;
	#println(Mbar1);
end=#
