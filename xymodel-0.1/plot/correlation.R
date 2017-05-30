## change
pdf("corr.2.pdf")
tmp <- read.table("corr.2.dat", col.names=c("T", "dist", "cor", "err"))
#########

corr <- split(tmp,tmp$T)


poly_fit <- numeric()
poly_exp <- numeric()
poly_exp_err <- numeric()
exp_fit <- numeric()
exp_c <- numeric()
exp_c_err <- numeric()
T <- as.numeric(names(corr))
for (i in seq(1,length(corr)))
{
	data <- corr[[i]]
	#data <- data[2:dim(data)[1],]
	data <- data[which(data$dist != 0),]
	if (0==var(data$cor))
	{
		# fit (nls) does not work for non noisy data
		poly_exp[i] <- 0
		poly_exp_err[i] <- 0
		poly_fit[i] <- 0

		exp_fit[i] <- 0
		exp_c[i] <- Inf
		exp_c_err[i] <- 0
	}
	else
	{
		p_fit <- nls(cor ~ a*(dist)^-eta, data=data, start=list(eta=2,a=1))
		poly_exp[i] <- coefficients(p_fit)[1]
		poly_exp_err[i] <- summary(p_fit)$coefficients['eta','Std. Error']
		poly_fit[i] <- summary(p_fit)$sigma

		e_fit <- nls(cor ~ b*exp(a*dist), data=data, start=list(a=-1, b=1))
		exp_fit[i] <- summary(e_fit)$sigma
		exp_c[i] <- -1/coefficients(e_fit)[1]
		exp_c_err[i] <- exp_c[i]^2*summary(e_fit)$coefficients['a', 'Std. Error']
	}



	plot(cor ~ dist, data=data, main=paste("T =", T[i]), ylim=c(0:1), xlab="r", ylab="correlation")
	if (0==var(data$cor))
	{
		abline(h=1, col="red")
		abline(h=1, col="blue")
	}
	else
	{
		lines(spline(data$dist, predict(p_fit)), col="blue")
		lines(spline(data$dist, predict(e_fit)), col="red")
	}
}

T_C <- 0.892
T <- as.numeric(names(corr))

plot(T, poly_fit, type='b', ylab="error G(r) ~ |r|^-eta(T)")
abline(v=T_C, lty=2)

plot(T, exp_fit, type='b', ylab="error G(r) = exp(c*r)")
abline(v=T_C, lty=2)

plot(T, poly_fit-exp_fit, type='b', ylab="delta error [poly - exp]")
abline(v=T_C, lty=2)
abline(0,0)


plot(T, poly_exp, type="b", ylab="eta(T)")
abline(v=T_C, lty=2)
abline(h=1/4, lty=2)
arrows(T, poly_exp + poly_exp_err, T, poly_exp -poly_exp_err, angle=90, code=3, length=0.1)

exp_c_right = exp_c[which(T >= 1)]
T_right = T[which(T>=1)]
plot(T, exp_c, type="b", ylab = "correlation length", ylim=c(0,10))
arrows(T, exp_c + exp_c_err, T, exp_c -exp_c_err, angle=90, code=3, length=0.1)
#exp_c_fit <- nls(exp_c_right ~ a*1/sqrt(Tc/(T_right-Tc))+b, start=list(a=-1.4, Tc = T_C, b=-0.1),
#algorithm="port", upper=c(Inf, 1))
#exp_c_fit <- nls(exp_c_right ~ a*log(T_right/2)+b, start=list(a=-1.4, b=-1))
#lines(spline(T_right, predict(exp_c_fit)))
abline(v=T_C, lty=2)
