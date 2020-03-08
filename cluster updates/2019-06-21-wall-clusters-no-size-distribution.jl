

## from NI
function otherDirs(dir)
	# given a direction, returns the other two directions
	# see if its okay to flip this
	return dir==1 ? [2,3] : ( dir==2 ? [3,1] : [1,2] )
end


function adjoinedWalls(i,j,k,a,p,nubbinData,NX,NY,NZ)

## this is the site of the dual lattice at the other end of the dual-lattice link ijka:
    (i1,j1,k1) = shift(i,j,k, a-1,NX,NY,NZ)
    nubbins=zeros(Int8,2)
    nubbins[1]=nubbinData[i,j,k,2*(a)]  # a+  ## here i've followed nabil's 1:3 convention
    nubbins[2]=nubbinData[i1,j1,k1, 2*(a)-1] # a-
    # now build the table of connections: note this thing copied from the emails
    #= tells me what the nubbin does
    0: no faces
    1: one DW passing through, no ambiguity
    2: connections are: 1-2, 3-4
    3: connections are 1-3, 2-4
    4: connections are 1-4, 2-3.  =#

    #tryWalls = zeros(Int8,5,4);
    #tryWalls[1]

# we want to enumerate the faces which touch the link ijkla
# and say whether they are occupied in the 5th entry
    od = otherDirs(a);
    shift1 = shiftReverse(i,j,k,od[1]-1,NX,NY,NZ);
	shift2 = shiftReverse(i,j,k,od[2]-1,NX,NY,NZ);
    tryWalls = [shift1' od[2] p[shift1[1],shift1[2],shift1[3],od[2]];
		shift2' od[1] p[shift2[1],shift2[2],shift2[3],od[1]];
		i j k od[2] p[i,j,k,od[2]];
		i j k od[1] p[i,j,k,od[1]]];



    connections = zeros(Int8,4,4);

    for i = 1:2
		if nubbins[i]==1
			# then this is a single DW. ironically we then have to work harder.
			for m = 1:4
				for n = 1:4
					if (tryWalls[m,5]==1 && tryWalls[n,5]==1)
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



    return (tryWalls, connections)
end



function makeClusterLabels(p, nubbinData, NX,NY,NZ)

    labels = zeros(Int64, NX*NY*NZ*3)
    wallCluster = zeros(Int64, NX,NY,NZ,3)  # this was called label before

    currentLabel = 0;
    for i=1:NX*NY*NZ*3
        labels[i]=i
    end

    for x=1:NX, y=1:NY, z=1:NZ, dir = 1:3
            (tryWalls, connections) = adjoinedWalls(x,y,z,dir,p,nubbinData,NX,NY,NZ)
            for wall_i =1:4
                if tryWalls[wall_i,5]==1 # if occupied
                    (i,j,k,a) = tryWalls[wall_i, 1:4]
                    if wallCluster[i,j,k,a]==0 # if not yet visitied
                        currentLabel+=1;
                        wallCluster[i,j,k,a] = currentLabel
                    end
                    for wall_j = 1:(wall_i-1)  # we want to hit each *connnection* once.
                        if tryWalls[wall_j, 5]==1  # if that one is occupied too
                            if connections[wall_i, wall_j]==1
                                (i1,j1,k1,a1) = tryWalls[wall_j,1:4];
                                unionify(wallCluster[i,j,k,a], wallCluster[i1,j1,k1,a1],labels)
                            end
                        end
                    end
                end

            end
    end

    (label, totalClusters)=makeSequential(wallCluster, labels, NX,NY,NZ)
    return (label, totalClusters)
end





function adjoinedWallsNoTouching(i,j,k,a,p,NX,NY,NZ)

## this is the site of the dual lattice at the other end of the dual-lattice link ijka:
#    (i1,j1,k1) = shift(i,j,k, a-1,NX,NY,NZ)
##    nubbins=zeros(Int8,2)
#    nubbins[1]=nubbinData[i,j,k,2*(a)]  # a+  ## here i've followed nabil's 1:3 convention
#    nubbins[2]=nubbinData[i1,j1,k1, 2*(a)-1] # a-


# we want to enumerate the faces which touch the link ijkla,
# and say whether they are occupied in the 5th entry
    od = otherDirs(a);
    shift1 = shiftReverse(i,j,k,od[1]-1,NX,NY,NZ);
	shift2 = shiftReverse(i,j,k,od[2]-1,NX,NY,NZ);
    tryWalls = [shift1' od[2] p[shift1[1],shift1[2],shift1[3],od[2]];
		shift2' od[1] p[shift2[1],shift2[2],shift2[3],od[1]];
		i j k od[2] p[i,j,k,od[2]];
		i j k od[1] p[i,j,k,od[1]]];

    return tryWalls
end

function makeClusterLabelsNoTouching(p, NX,NY,NZ)

    labels = zeros(Int64, NX*NY*NZ*3)
    wallCluster = zeros(Int64, NX,NY,NZ,3)  # this was called label before

    currentLabel = 0;
    for i=1:NX*NY*NZ*3
        labels[i]=i
    end

    for x=1:NX, y=1:NY, z=1:NZ, dir = 1:3
            tryWalls = adjoinedWallsNoTouching(x,y,z,dir,p,NX,NY,NZ)
            for wall_i =1:4
                if tryWalls[wall_i,5]==1 # if occupied
                    (i,j,k,a) = tryWalls[wall_i, 1:4]
                    if wallCluster[i,j,k,a]==0 # if not yet visitied
                        currentLabel+=1;
                        wallCluster[i,j,k,a] = currentLabel
                    end
                    for wall_j = 1:(wall_i-1)  # we want to hit each *connnection* once.
                        if tryWalls[wall_j, 5]==1
                            # if that one is occupied too, then they must be in the same cluster, given the no touching rules.
#                            if connections[wall_i, wall_j]==1
                                (i1,j1,k1,a1) = tryWalls[wall_j,1:4];
                                unionify(wallCluster[i,j,k,a], wallCluster[i1,j1,k1,a1],labels)
#                            end
                        end
                    end
                end

            end
    end

    (label, totalClusters) =makeSequential(wallCluster, labels, NX,NY,NZ)
    return (label, totalClusters)
	#, sizeOfCluster)
end




function findCanonical(x::Int64, labels)
	y = x;
	while (labels[y]!=y)
		y = labels[y];
	end

#	while (labels[x]!=x)
#		z = labels[x]
#		labels[x] = y;
#		x=z;
#	end
	return y
end


# this function makes a new equivalence relation
function unionify(x::Int64, y::Int64,labels)
	if x>y
		labels[findCanonical(x, labels)] = labels[findCanonical(y, labels)];
	else
		labels[findCanonical(y, labels)] = labels[findCanonical(x, labels)];
	end
    return labels
end



function makeSequential(label, labels, NX,NY,NZ)
    newLabels = zeros(Int64, NX*NY*NZ*3)
    labelIterator= 0;
#    sizeOfCluster = zeros(Int64, NX*NY*NZ*3)
    for x=1:NX, y =1:NY, z = 1:NZ, dir=1:3
        i = label[x,y,z,dir]
        if i > 0
            j = findCanonical(i, labels);
            if newLabels[j]==0
                labelIterator +=1;
                newLabels[j] = labelIterator;
#                sizeOfCluster[labelIterator]=1;
#            else
#                sizeOfCluster[newLabels[j]] += 1;
            end
            label[x,y,z,dir] = newLabels[j]
        end
    end
    totalClusters = labelIterator
    return (label, totalClusters)
	#, sizeOfCluster)
end







##
