source("utility.R")
source("rcuda.R")
library(MASS)

probit_mcmc<-function(
		y,           # vector of length n 
		X,           # (n x p) design matrix
		beta_0,      # (p x 1) prior mean
		Sigma_0_inv, # (p x p) prior precision 
		niter,       # number of post burnin iterations
		burnin,      # number of burnin iterations
		maxRejection,
		useCPU=TRUE,
		dims,
		ptxfile
		)
{
	N=length(y)
	samp=rep(0,length(beta_0)) # initiate start state will be 0's
	postVar=solve(Sigma_0_inv + t(X)%*%X) # Posterior variance.
	constant=Sigma_0_inv%*%beta_0 # Used in posterior mean calculation.
	samples=matrix(0, ncol=ncol(X), nrow=(niter+burnin))

	upper=which(y>0) # these rows sample from TruncNormal in [0, Inf)
	lower=which(y<=0) # these rows sample from TruncNormal in (-Inf, 0]
	for(i in 1:(niter+burnin))
	{
		# Sample from upper TruncNorm
		if(useCPU==TRUE)
			z.upper=truncnormCPU(length(upper), X[upper,]%*%beta_0, rep(1,length(upper)), rep(0,length(upper)), rep(Inf,length(upper)), maxRejections)
		else
			z.upper=truncnormGPU(ptxfile, length(upper), X[upper,]%*%beta_0, rep(1,length(upper)), rep(0,length(upper)), rep(Inf,length(upper)), maxRejections, block_dims=dims$block_dims, grid_dims=dims$grid_dims, randNumGenId=i)

		# Sample from lower TruncNorm
		if(useCPU==TRUE)
			z.lower=truncnormCPU(length(lower), X[lower,]%*%beta_0, rep(1,length(lower)), rep(-Inf,length(lower)), rep(0,length(lower)), maxRejections)
		else
			z.lower=truncnormGPU(ptxfile, length(lower), X[lower,]%*%beta_0, rep(1,length(lower)), rep(-Inf,length(lower)), rep(0,length(lower)), maxRejections, block_dims=dims$block_dims, grid_dims=dims$grid_dims, randNumGenId=i)

		# Put the sampled TruncNorm values into a vector
		samp[upper]=z.upper
		samp[lower]=z.lower

		# Recalculate posterior mean every 100 iterations
		if(i%%100==1 && i<burnin) betaMean=postVar %*% (constant + t(X)%*%samp)
			beta_0=mvrnorm(1, betaMean, postVar)

		# Put the beta estimates in a matrix
		samples[i,]=beta_0
		}

	invisible(mcmc(samples[(burnin+1):nrow(samples),]))
}

files=paste0("data_0", 1:5, ".txt") # list of files: data_01, ..., data_05
niter=2000
burnin=500
maxRejections=2000
# CPU times for data_01, data_02, ..., data_05
# Note: CPU and GPU loops are separate in case something happens.
for(i in 1:length(files))
{
	# read data and set up values for the probit mcmc
	dat=read.table(files[i], colClasses="numeric", header=TRUE)
	N=nrow(dat)
	X=as.matrix(dat[,2:ncol(dat)])
	y=as.vector(dat[,1])
	beta_0=rep(0,ncol(X))
	Sigma_0_inv=matrix(0,ncol=ncol(X),nrow=ncol(X))
	dims=compute_grid(N)

	# write means and compute times to files
	fileconnection=paste0("results_timesCPU_", i, ".txt")
	fileconnection_means=paste0("results_meansCPU_", i, ".txt")
	fileconnection2=paste0("results_timesGPU_", i, ".txt")
	fileconnection_means2=paste0("results_meansGPU_", i, ".txt")
	CPUtime=system.time(sample<-probit_mcmc(y,X,rep(0,ncol(X)),Sigma_0_inv,niter,burnin,maxRejection=maxRejections, useCPU=TRUE))['user.self']
	write(CPUtime, file=fileconnection, append=TRUE)
	write(colMeans(sample), file=fileconnection_means, append=TRUE)
}

# GPU times for data_01, data_02, ..., data_05
for(i in 1:length(files))
{
	# read data and set up values for the probit mcmc
	dat=read.table(files[i], colClasses="numeric", header=TRUE)
	N=nrow(dat)
	X=as.matrix(dat[,2:ncol(dat)])
	y=as.vector(dat[,1])
	beta_0=rep(0,ncol(X))
	Sigma_0_inv=matrix(0,ncol=ncol(X),nrow=ncol(X))
	dims=compute_grid(N)

	# write means and compute times to files
	fileconnection2=paste0("results_timesGPU_", i, ".txt")
	fileconnection_means2=paste0("results_meansGPU_", i, ".txt")
	GPUtime=system.time(sample<-probit_mcmc(y,X,rep(0,ncol(X)),Sigma_0_inv,niter,burnin,maxRejection=maxRejections, useCPU=FALSE, dims, ptxfile="tn.ptx"))['user.self']
	write(GPUtime, file=fileconnection2, append=TRUE)
	write(colMeans(sample), file=fileconnection_means2, append=TRUE)
}

