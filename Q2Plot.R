# Done after obtaining results for data_01, ..., data_05

#######################################
# ~~==| Make plots for question 2 |==~~
#######################################

CPUtime=read.table("results_timesCPUall.txt", colClasses="numeric")
CPUtime=as.numeric(CPUtime[,1])
GPUtime=read.table("results_timesGPUall.txt", colClasses="numeric")
GPUtime=as.numeric(GPUtime[,1])

k=3:7
xrange=range(log(10^k))
yrange=range(c(CPUtime, GPUtime))
pdf("plot_runtimeDatas.pdf")
plot(NA, ylim=yrange, xlim=xrange, xlab="log(n)", ylab="time (s)",
     main="Runtimes of CPU and GPU - Probit MCMC")
lines(log(10^(k[1:length(CPUtime)])), CPUtime, lty=2)
lines(log(10^(k[1:length(GPUtime)])), GPUtime)
legend("topleft", legend=c("CPU", "GPU"), lty=c(2,1))
dev.off()

