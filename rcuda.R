library(RCUDA)
library(coda)
cuGetContext(TRUE)

truncnormCPU<-function(N,mu,sigma,a,b,maxRejection)
{
	# A counter. Stops rejection-sampling when sampleCount==maxRejection
	sampleCount=0
	x=rnorm(N, mu, sigma) # randomly sample from Normal
	# resample has indexes for which x is not in [a,b].
	# Then, x[resample] needs to be resampled.
	resample=which(x<a | x>b)

	while(length(resample)>0 && sampleCount<maxRejection)
	{
		x[resample]=rnorm(length(resample), mu[resample], sigma[resample])
		# Below: Find the index for which x is not in [a,b] again,
		# but we only need to check the ones that were resampled.
		resample.tmp=which(x[resample]<a[resample] | x[resample]>b[resample])
		# resample.tmp lists indexes of x.subset that needs to be
		# changed. Now, we get the actual index needed to be changed for x.
		resample=resample[resample.tmp]
		sampleCount=sampleCount+1
	}

	# If the max number of tries failed, then just sample from a uniform.
	if(length(resample)>0)
	{
		# runif does not work on Inf, so replace -Inf/Inf with min/max number for R
		a[resample][which(is.infinite(a[resample])==TRUE)]=.Machine$double.xmin
		b[resample][which(is.infinite(b[resample])==TRUE)]=.Machine$double.xmax
		x[resample]=runif(length(resample), a[resample], b[resample])
	}

	invisible(x) # Returns x
}

truncnormCPU_Robert<-function(N,mu,sigma,a,b,maxRejection)
{
	sampleCount=0 # A counter. Stops rejection-sampling when sampleCount==maxRejection
	sign.a=ifelse(a>=0, 1, -1)
	a=abs(a)
	sign.b=ifelse(b>=0, 1, -1)
	b=abs(b)
	alpha=a+sqrt(a^2+4)/2

	z=rexp(N, alpha)
	g=exp(-(z-alpha)^2/2)
	u=runif(N)
	resample=which(g>u) # Has indexes for which x is not in [a,b]. Then, x[resample] needs to be resampled.

	while(length(resample)>0 && sampleCount<maxRejection)
	{
		x[resample]=rexp(length(resample), alpha[resample])
# Below: Find the index for which x is not in [a,b] again, but we only need to check the ones that were resampled.
		g=exp(-(x[resample]-alpha[resample])^2/2)
		u=runif(length(resample))
		resample.tmp=which(g>u)
		resample=resample[resample.tmp] # resample.tmp lists indexes of x.subset that needs to be changed. Now, we get the actual index needed to be changed for x.
		sampleCount=sampleCount+1
	}

# If the max number of tries failed, then just sample from a uniform.
	if(length(resample)>0)
	{
		x[resample]=runif(length(resample), a[resample], b[resample])
	}

	invisible(x) # Returns x
}

truncnormGPU<-function(file, N, mu, sigma, a, b, maxRejection, grid_dims, block_dims, returnTimes=FALSE, randNumGenId)
{
	# Get the kernel function from a ptx file.
	m = loadModule(file)
	k = m$rtruncnorm_kernel

	t1<-system.time({
			x_device<-copyToDevice(numeric(N));
			mu_device<-copyToDevice(mu);
			sigma_device<-copyToDevice(sigma);
			a_device<-copyToDevice(a);
			b_device<-copyToDevice(b);
			})['user.self']

	t2<-system.time(.cuda(k, x_device, N, mu_device, sigma_device,
			 a_device, b_device, length(mu), length(sigma),
			 length(a), length(b), maxRejection,
			 gridDim=grid_dims, blockDim=block_dims, randNumGenId))['user.self']

	t3<-system.time(x<-copyFromDevice(obj=x_device, nels=x_device@nels, type="float"))['user.self']
	rm(x_device, mu_device, sigma_device, a_device, b_device)

	if(returnTimes==TRUE)
		invisible(c(t1,t2,t3))
	else
		invisible(x)
}

