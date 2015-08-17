source("utility.R")
source("rcuda.R")

timeTruncNormSample<-function(N, mu, sigma, a, b, maxRejection, pureR=TRUE, dims)
{
if(pureR==TRUE){
	return(system.time(truncnormCPU(N,mu,sigma,a,b,maxRejection))['user.self'])
}else{ # else, use GPU
	return((truncnormGPU(file="tn.ptx", N, mu, sigma, a, b, maxRejection,
			     dims$grid_dims, dims$block_dims,returnTimes=TRUE)))
}
}

#######################################
# ~~==| Make plots for question 1 |==~~
#######################################

k=1:8 # Will use samples of size 10^k (i.e. 10^1 10^2 ... 10^8)
# The mu, sigma, a, and b vectors have the same values, so don't need
# to make the entire vector. Instead, set them as a value and use recycling.
mu = 2
sigma = 1
a=0 # lower bound
b=1.5 # upper bound
maxRejection=2000

# Time it takes CPU to sample for N=10^k, where k=1,2,...8.
fileconnection="results_timesCPU.txt" # File to save the CPU times to
fileconnection2="results_timesGPU.txt" # "" "" "" "" GPU times to
for(i in k)
{
	# check if number of threads is valid
	size=10^i
	dims=compute_grid(size)
	nthreads <- prod(dims$grid_dims)*prod(dims$block_dims)
	if (nthreads < size){
	    stop("Grid is not large enough...!")
	}

	# CPU time
	timed=timeTruncNormSample(N=size, rep(mu, length.out=size), 
				  rep(sigma, length.out=size),
				  rep(a,length.out=size),
				  rep(b,length.out=size),
				  maxRejection=maxRejection, pureR=TRUE)
	write(as.numeric(timed), file=fileconnection, append=TRUE)

	# GPU time
	timed=timeTruncNormSample(N=size, mu, sigma, a, b,
				  maxRejection=maxRejection, pureR=FALSE,
				  dims=dims)
	write(as.numeric(timed), file=fileconnection2, append=TRUE)
}

 CPUtime=read.table("results_timesCPU.txt", colClasses="numeric")
 CPUtime=as.numeric(CPUtime[,1])
 GPUtime.mat=read.table("results_timesGPU.txt", colClasses="numeric")
 GPUtime=colSums(GPUtime.mat)

xrange=range(log(10^k))
yrange=range(c(CPUtime, GPUtime))
pdf("plot_runtime.pdf")
plot(NA, ylim=yrange, xlim=xrange, xlab="log(n)", ylab="time (s)",
     main="Runtimes of CPU and GPU")
lines(log(10^k), CPUtime, lty=2)
lines(log(10^(1:length(GPUtime))), GPUtime)
legend("topleft", legend=c("CPU", "GPU"), lty=c(2,1))
dev.off()

# CPUtime=read.table("results_timesCPU.txt", colClasses="numeric")
# CPUtime=as.numeric(CPUtime[,1])
# GPUtime.mat=read.table("results_timesGPU.txt", colClasses="numeric")
# GPUtime=colSums(GPUtime.mat)

