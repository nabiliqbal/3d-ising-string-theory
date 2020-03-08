# this is all the data analysis subroutines


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

# this performs a jacknife operation on whatever you pass to it.

function jackknife(data,operation::Function)
	# this performs a jacknife analysis to get the error (following "Statistical
	# Analysis of Simulations: Data Correlations and Error Estimation") by
	# Wolfhard Janke, section 3.6

	# it is assumed that the operation acts on an array and gives the observable


	N = length(data);

	# create a table of jacknifey values
	estimate = zeros(N);

	# remove the n-th guy.
	for n = 1:N
		#removed = vcat(data[1:n-1], data[n+1:N]);
		#println(removed);
		estimate[n] = operation(vcat(data[1:n-1], data[n+1:N]));
	end

	# this will return the standard deviation and the estimate
	return (mean(estimate),sqrt((N-1)/N * sum((estimate .- mean(estimate)).^2)));
end

function test(data)
	for i =1:4
       data = decimate(data);
       println(jackknife(data, calcBinder), " ", calcBinder(data));
   end
end

function contBinning(data, binSize, operation)
	# this is a bit different, it lets us do a continuous variation of the bin size
	# rather than the logarithmic thing that Troyer advocated.
	# lets give it a try

	nBins = convert(Int32,floor(length(data)/binSize));
	binned = zeros(nBins);
	for i = 1:nBins
		binned[i] = mean(data[(i-1)*binSize+1:i*binSize]);
	end
	return jackknife(binned,operation);
	#return meandev(binned);
end

# this generates a binning summary with a certain maxBinSize
function binningSummary(data,maxBinSize,operation)
	vals = zeros(maxBinSize);
	devs = zeros(maxBinSize);
	for k = 1:maxBinSize
		(vals[k],devs[k])=contBinning(data,k,operation);
	end
	display(plot(devs));
	return (vals, devs);
end

function binning(data)
	# this does a "binning" analysis. check how many times I can do it while keeping
	# 30 bins
	num_steps = convert(Int32,floor(log2(length(data)/30)));
	println(num_steps);
	devs_norm = zeros(num_steps);
	#devs_jack = zeros(num_steps);
	for i = 1:num_steps
		data = decimate(data);
		devs_norm[i]=meandev(data)[2];
		#devs_jack[i] = jacknife(data);
		println(meandev(data)," ",devs_norm[i]);
	end
	display(plot(devs_norm));
	return devs_norm;
end

function binningJack(data, operation::Function)
	# this is a different binning operation! it does the same thing as the
	# other one, but for an arbitrary operation (called operation). it also does
	# it with the jacknife method.

	# this does a "binning" analysis. check how many times I can do it while keeping
	# 30 bins
	num_steps = convert(Int32,floor(log2(length(data)/30)));
	println(num_steps);
	means_jack = zeros(num_steps+1);
	devs_jack = zeros(num_steps+1);
	(means_jack[1],devs_jack[1]) = jackknife(data,operation);

	for i = 2:num_steps+1;
		data = decimate(data);
		(means_jack[i],devs_jack[i]) = jackknife(data,operation);
		println(means_jack[i]," ",devs_jack[i]);
	end
	display(plot(devs_jack));
	return (means_jack, devs_jack);
end

function calcBinder(data)
	# just calculates the binder coefficient:
	return (3-mean(data.^4)/(mean(data.^2)^2))/2;
end
