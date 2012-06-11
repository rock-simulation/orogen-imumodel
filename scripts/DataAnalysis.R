#Rgui.exe --sdi --max-mem-size=2047Mb

#libraries
library(stats)

setwd("C:/Documents and Settings/Javier Hidalgo/My Documents/TEC-MMA/ExoMars/ExoMaDeR/GuidanceSystem/AllanVariable")

values <- read.table ("Data/22122009-150533.data", sep="\t")
values <- read.table ("Data/23122009-103116.data", sep="\t")
values <- read.table ("Data/07012010-201649.data", sep="\t")
values <- read.table ("Data/15012010-976herz.data", sep="\t")
values <- read.table ("Data/15012010-9herz.data", sep="\t")
values <- read.table ("Data/15012010-195700.data", sep="\t", nrows = 1000000)
values <- read.table ("Data/19012010-185825-6hrs-588Hz.data", sep="\t",  nrows = 2000000)
values <- read.table ("Data/TestCenter/28012010-153920.data", sep="\t")

/* Test center data acquisition */
values <- read.table ("Data/TestCenter/28012010-153920.data", sep="\t")
values <- read.table ("Data/TestCenter/Inclinometers/28012010-163846.data", sep="\t")
values <- read.table ("Data/TestCenter/Accelerometers/01022010-110510.data", sep="\t")








names(values) = c("time","gyrox", "gyroy", "gyroz", "accx", "accy", "accz", "clock")

names(values) = c("time","inclx", "incly", "inclz", "clock")

names(values) = c("time","accx", "accy", "accz", "clock")

names(values) = c("time","gyrox", "gyroy", "gyroz", "clock")





delta=NULL
for (i in 2:length(values$time))
{
	delta[i-1] = values$time[i]-values$time[i-1]
}

imt30Gyros <- ts (data.frame(values$gyrox, values$gyroy, values$gyroz), start=c(values$time[1]), delta=mean(delta))

imt30Acc <- ts (data.frame(values$accx, values$accy, values$accz), start=c(values$time[1]), delta=mean(delta))

imt30Incl <- ts (data.frame(values$inclx, values$incly, values$inclz), start=c(values$time[1]), delta=mean(delta))




windows()
plot(imt30Gyros , lag(imt30Gyros , 1), cex = .8, col="blue", main = substitute(Gyros~values~(rad/s)~-~Sampling~rate~(F~Hertz),list(F=frequency (imt30Gyros))))

windows()
plot(imt30Acc, lag(imt30Acc , 1), cex = .8, col="blue", main = "Acc values (m/s)")

windows()
plot(imt30Incl, lag(imt30Incl, 1), cex = .8, col="blue", main = "Incl values (m/s)")



#structure of time serie (ts)
# imt30Gyros@.Data[1,1]
# imt30Gyros@.Data[1,2]
# imt30Gyros@.Data[1,3]
# first component is the value and second component is the gyro (x, y or z)
# imt30Gyros@.Data[1:10,1]

#javi <- ts (c(1,1,1,1,1,1,1,1,1,1), delta=1)
#javi <- ts (c(1,2,3,4,5,6,7,8,9,10,11,12), delta=1)
#javi <- ts (rep(1, 10000), delta=1)
#N = length(javi)
#tau = 1/frequency(javi)
#n = ceiling((N-1)/2)
#2^floor (log10(n)/log10(2))

#####################################################################
######## Verification of the AV using the IEEE std 1554-2005 ########
#####################################################################

## White noise simulation ##
sequ<- rnorm (2^16, mean=0, sd=1)#unit Gaussian distribution (sigma = 1)
sigma <- (6.04E-05/sqrt(0.01))
sequ <- sequ * sigma
wnoise <- ts (sequ, delta=0.001024)

avwnoise <- avar (wnoise@.Data, frequency (wnoise))
avwnoise <- avari (wnoise@.Data, frequency (wnoise))

#integrate the value to process the avar2
sequi = rep(0, length(sequ))
sequi[1] = sequ[1]
for (k in 2: (length(sequ)))
{
	sequi[k] = sequi[k-1]+sequ[k]
}

windows()
plot (wnoise)
wnoisei <- ts (sequi, delta=1)
avwnoise2 <- avari (wnoisei@.Data, frequency (wnoisei))


windows()
plot(avwnoise$time,sqrt(avwnoise$av),log= "xy", xaxt="n" , yaxt="n", col="blue", xlab="", ylab="")
lines(avwnoise$time, sqrt(avwnoise$av),type="l", col="green4", lwd=4)
lines(avwnoise1$time, sqrt(avwnoise1$av),type="l", col="green")
lines(avwnoise2$time, sqrt(avwnoise2$av),type="l", col="orange")
axis(1, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
axis(2, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
grid(equilogs=TRUE, lwd=1, col="orange")
title(main = "Allan variance Analysis", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation")
legend(100, 3e-01, c("avar (logN)", "avar (N)", "avari"),  fill = c("blue", "green", "orange"))

#Check value at sigma*sqrt(delta)

#Find the T=1 sec
approx (x=c(avwnoise$time[7], avwnoise$time[8]), y= c(sqrt(avwnoise$av[7]), sqrt(avwnoise$av[8])), n=100)



## Random walk simulation ##
sequ <- rep (0, 2^16)
At <- 0.01024
m <- 10
b <- 3.31e-06

sequ[1] = 0
for (i in 2: (2^16))
{
	y <- rnorm (m, mean=0, sd=1)
	sequ[i] = sequ[i-1]+ b*sum(y[1:m])*sqrt (At/m)
	#print(sequ[i])
}

rwnoise <- ts (sequ, delta=At)

avrwnoise <- avar (rwnoise@.Data, frequency (rwnoise))
avrwnoise1 <- avar1 (rwnoise@.Data, frequency (rwnoise))

#integrate the value to process the avar2
sequi = rep(0, length(sequ))
sequi[1] = sequ[1]
for (k in 2: (length(sequ)))
{
	sequi[k] = sequi[k-1]+sequ[k]
}

rwnoisei <- ts (sequi, delta=At)

windows()
plot (rwnoise)

avrwnoise2 <- avari (rwnoisei@.Data, frequency (rwnoisei))


windows()
plot(avrwnoise$time,sqrt(avrwnoise$av),log= "xy", xaxt="n" , yaxt="n", type="l", col="blue", xlab="", ylab="")
lines(avrwnoise$time, sqrt(avrwnoise$av), type="l", col="green")
lines(avrwnoise1$time, sqrt(avrwnoise1$av), type="l", col="green")
lines(avrwnoise2$time, sqrt(avrwnoise2$av), type="l", col="orange")
axis(1, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
axis(2, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
grid(equilogs=TRUE, lwd=1, col="black")
title(main = "Allan variance Analysis", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation")
legend(10, 100, c("avar (logN)", "avar (N)", "avari"),  fill = c("blue", "green", "orange"))

#check value at b/sqrt(delta)

############################################
# Estimation Quality of the Allan Variance #
############################################
#In practice, estimation of the Allan variance is based on a finite number of independent
#clusters that can be formed from any finite length of data. The Allan variance of any
#noise terms is estimated using the total number of clusters of a given length that can be
#created. The confidence of the estimation improves as the number of independent clusters
#is increased.

#delta_avar as the percentage error in estimating the Allan standard
#deviation of the cluster due to the finiteness of the number of clusters gives

#Based on Papoulis 1991

sigma_delta_avar = 1/sqrt(2((N/n)-1))

#Where N is the total number of data points in the entire data set
#and n is the number of the points contained in the cluster

#Example
N <- 100000
n <- seq(1,ceiling((N-1)/2))

sigma_delta_avar = 1/sqrt(2*((N/n)-1))

windows()
plot (n,1-sigma_delta_avar, type="l", col = "blue")
grid(10,10, col="black")
legend (N/4, (1-sigma_delta_avar[2]), c("Quality of Allan Variance"), fill = c("blue"))


##############################################
# Calculate the Allan Variance of imt30Gyros #
##############################################
# N = is the number of total samples
# tau = is the rate= 1/frequency(x)
# n = max. size of clusters.

N = length(imt30GyroX)
tau = 1/frequency(imt30GyroX)
n = ceiling((N-1)/2)
p = floor (log10(n)/log10(2))

#Allan variance array
av <- rep (0,p+1)
#Time array
time <- rep(0,p+1)

# Minimal cluster size is 1 and max is n
# in time would be 2*tau and max time would be n*tau
for (i in 0:(p))
{
	print("i")
	print(i)
	print (2^i)
	print (floor(N/(2^i)))
	omega = rep(0,floor(N/(2^i)))
	T = (2^i)*tau

	l <- 1
	k <- 1
	while (k <= floor(N/(2^i)))
	{
		print ("l")
		print(l)
		print (l+((2^i)-1))
		omega[k] = sum (javi@.Data[l:(l+((2^i)-1))])/T
		l <- l + (2^i)
		k <- k + 1
		print (omega)
	}
	sumvalue <- 0

	for (k in 1: (length(omega)-1))
	{
		sumvalue = sumvalue + (omega[k+1]-omega[k])^2
		print(sumvalue)
	}

	av[i+1] = sumvalue/(2*(length(omega)-1))
	time[i+1] = T
	print("Allan variance")
	print (av[i+1])
}

?plot.default

windows()
plot (time,sqrt(av), log= "xy", yaxt="n")
lines (time,sqrt(av))
axis(2, c(0.00001, 0.0001, 0.001, 0.01, 0.1))
grid(equilogs=TRUE)

windows()
plot (time[1:25000],sqrt(av[1:25000]), log= "xy", yaxt="n")
lines (time[1:25000],sqrt(av[1:25000]), col="red")
axis(2, c(0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 10))
grid(equilogs=TRUE)

lines (lm(time ~ sqrt(av)), col="green")


allan2 <- function (values, freq)
{

	N = length(values)
	tau = 1/freq
	n = ceiling((N-1)/2)

	print (N)
	print(tau)
	print(n)

	#Allan variance array
	av <- rep (0,n)
	#Time array
	time <- rep(0,n)

	# Minimal cluster size is 1 and max is n
	# in time would be 1*tau and max time would be n*tau
	for (i in 1:(n))
	{
		print(i)
		print (floor(N/i))
		omega = rep(0,floor(N/i))
		T = i*tau

		l <- 1
		k <- 1
		while (k <= floor(N/i))
		{
			omega[k] = sum (values[l:(l+(i-1))])/T
			l <- l + i
			k <- k + 1
		}
		sumvalue <- 0
	
		for (k in 1: (length(omega)-1))
		{
			sumvalue = sumvalue + (omega[k+1]-omega[k])^2
		}
	
		av[i] = sumvalue/(2*(length(omega)-1))
		time[i] = T
		print("Allan variance")
		print (av[i])
	}

	return (data.frame(av, time))
}
###################################################
############# Allan with group of 2^p #############
###################################################
allan <- function (values, freq)
{

	N = length(values)
	tau = 1/freq
	n = ceiling((N-1)/2)
	p = floor (log10(n)/log10(2))

	print (N)
	print(tau)
	print(n)

	#Allan variance array
	av <- rep (0,p+1)
	#Time array
	time <- rep(0,p+1)

	# Minimal cluster size is 1 and max is 2^p
	# in time would be 1*tau and max time would be (2^p)*tau
	for (i in 0:(p))
	{
		print(i)
		print (floor(N/(2^i)))
		omega = rep(0,floor(N/(2^i)))
		T = (2^i)*tau

		l <- 1
		k <- 1
		while (k <= floor(N/(2^i)))
		{
			omega[k] = sum (values[l:(l+((2^i)-1))])/T
			l <- l + (2^i)
			k <- k + 1
		}
		sumvalue <- 0
	
		for (k in 1: (length(omega)-1))
		{
			sumvalue = sumvalue + (omega[k+1]-omega[k])^2
		}
	
		av[i+1] = sumvalue/(2*(length(omega)-1))
		time[i+1] = T
		print("Allan variance")
		print (av[i+1])
	}

	return (data.frame(av, time))
}



#### Calculating the Allan Variance for the gyros ####
length (imt30Gyros@.Data[,1])
frequency (imt30Gyros)

#Gyro X
avgyrox <- allan (imt30Gyros@.Data[,1], frequency (imt30Gyros))

#Gyro Y
avgyroy <- allan (imt30Gyros@.Data[,2], frequency (imt30Gyros))

#Gyro Z
avgyroz <- allan (imt30Gyros@.Data[,3], frequency (imt30Gyros))

#### Calculating the Allan Variance for the accelerometers ####
imt30Acc@.Data[1,1]
frequency (imt30Acc)

#Acc X
avaccx <- allan (imt30Acc@.Data[,1], frequency (imt30Acc))

#Acc Y
avaccy <- allan (imt30Acc@.Data[,2], frequency (imt30Acc))

#Acc Z
avaccz <- allan (imt30Acc@.Data[,3], frequency (imt30Acc))

#### Calculating the Allan Variance for the inclinometers ####
imt30Incl@.Data[1,1]
frequency (imt30Incl)

#Incl X
avinclx <- allan (imt30Incl@.Data[,1], frequency (imt30Incl))

#Incl Y
avincly <- allan (imt30Incl@.Data[,2], frequency (imt30Incl))

#Incl Z
avinclz <- allan (imt30Incl@.Data[,3], frequency (imt30Incl))



#### Plotting the results ####


windows()
plot (avgyrox$time,sqrt(avgyrox$av),log= "xy", xaxt="n" , yaxt="n", type="l", col="blue", xlab="", ylab="")
lines (avgyroy$time,sqrt(avgyroy$av), type="l", col="green")
lines (avgyroz$time,sqrt(avgyroz$av), type="l", col="red")
axis(1, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
axis(2, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
grid(equilogs=TRUE, lwd=1, col="orange")
title(main = "Allan variance Analysis", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation (rad/s)")



windows()
plot (avaccz$time,sqrt(avaccz$av),log= "xy", xaxt="n" , yaxt="n", type="l", col="red", xlab="", ylab="")
lines (avaccy$time,sqrt(avaccy$av), col="green")
lines (avaccx$time,sqrt(avaccx$av), col="blue")
axis(1, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
axis(2, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
grid(equilogs=TRUE, lwd=1, col="orange")
title(main = "Allan variance Analysis", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation (m/s)")

legend(10, 50, c("AccX", "AccY", "AccZ"),  fill = c("blue", "green", "red"))


windows()
plot (avinclx$time,sqrt(avinclx$av),log= "xy", xaxt="n" , yaxt="n", type="l", col="blue", xlab="", ylab="")
lines (avincly$time,sqrt(avincly$av), col="green")
lines (avinclz$time,sqrt(avinclz$av), col="red")
axis(1, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
axis(2, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
grid(equilogs=TRUE, lwd=1, col="orange")
title(main = "Allan variance Analysis", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation (m/s)")


######## creating cluster whitout superpousing 

#Minimal cluster size is 2 and max is n
for (i in 2:(n))
{
	print(i)
	print (floor(N/i))
	omega = rep(0,floor(N/i))
	T = (i-1)*tau

	k <- 1
	l <- 1
	while (k <= floor(N/i))
	{
		print ("l")
		print(l)
		print (l+(i-1))
		print (javi@.Data[l:(l+(i-1))])
		omega[k] = sum (javi@.Data[l:(l+(i-1))])/T
		l <- l + i
		k <- k + 1
	}
}

######## creating cluster superpousing values
N = length(imt30GyroX)
tau = 1/frequency(imt30GyroX)
n = ceiling((N-1)/2)

#Allan variance array
av <- rep (0,n)
#Time array
time <- rep(0,n)

# Minimal cluster size is 2 and max is n
# in time would be 2*tau and max time would be n*tau
for (i in 1:(n))
{
	print(i)
	print (N-(i-1))
	omega = rep(0,(N-(i-1)))
	T = i*tau

	l <- 1
	for (k in 1: (N-(i-1)))
	{
		print ("l")
		print(l)
		print (l+(i-1))
		omega[k] = sum (imt30GyroX@.Data[l:(l+(i-1))])/T
		l = l+1
		print (omega)
	}
	sumvalue <- 0
	for (k in 1: (length(omega)-1))
	{
		sumvalue = sumvalue + (omega[k+1]-omega[k])^2
		print(sumvalue)
	}

	av[i] = sumvalue/(2*(N - (i-1)))
	time[i] = T
	print("Allan variance")
	print (av[i])
}


for (i in 1:(n-1))
{
	print(i)
}

