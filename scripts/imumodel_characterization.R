###############################################
# Post Processing R Script to compute         
# Allan Variance of a IMU simulated data
# Javier Hidalgo Carrió
# DFKI-RIC March 2012
###############################################

#libraries
library(allanvar)

################
# ASQUARG Test 1
################
# You may change this path
setwd ("~/iMoby/iMoby-dev/simulation/orogen/imumodel/scripts")
load ("imumodel.0.log-imumodel.imuout.RData")


ls()
time <- data.frame(time)
names(time) <- seq(1, length(time))

acc <- as.data.frame(acc)
names(acc) <- seq(1, length(acc))

gyro <- as.data.frame(gyro)
names(gyro) <- seq(1, length(gyro))


#Data time adaptation
sec <- as.integer(sapply(time,function(x) x[1]))
usec <- as.integer(sapply(time,function(x) x[2]))

t1 <- as.numeric(paste(sec, usec, sep=''))

delta1 = NULL
for (i in 1:(length(t1)-1))
{
  delta1[i] <- (t1[i+1] - t1[i])/1000000
}

delta <- mean(delta1)


values <- NULL
values$accx <- unlist(acc [1, ])
values$accy <- unlist(acc [2, ])
values$accz <- unlist(acc [3, ])
values$gyrox <- unlist(gyro [1, ])
values$gyroy <- unlist(gyro [2, ])
values$gyroz <- unlist(gyro [3, ])

########################## X Axis ##########################
xsensAcc <- ts (data.frame(values$accx), start=c(t1[1]), delta=mean(delta))

#Frequency
frequency (xsensAcc)
mean (values$accx)
sd (values$accx)

#### Calculating the Allan Variance for the accelerometers ####
avacc1x <- avar (xsensAcc@.Data, frequency (xsensAcc))

########################## Y Axis ##########################
xsensAcc <- ts (data.frame(values$accy), start=c(t1[1]), delta=mean(delta))

#Frequency
frequency (xsensAcc)
mean (values$accy)
sd (values$accy)

#### Calculating the Allan Variance for the accelerometers ####
avacc1y <- avar (xsensAcc@.Data, frequency (xsensAcc))


########################## Z Axis ##########################
xsensAcc <- ts (data.frame(values$accz), start=c(t1[1]), delta=mean(delta))

#Frequency
frequency (xsensAcc)
mean (values$accz)
sd (values$accz)

#### Calculating the Allan Variance for the accelerometers ####
avacc1z <- avar (xsensAcc@.Data, frequency (xsensAcc))

#### Ploting Acc X
plot (xsensAcc)

#### Ploting the results ####
plot (avacc1x$time,sqrt(avacc1x$av),log= "xy", xaxt="n" , yaxt="n", type="l", col="blue", xlab="", ylab="")
lines (avacc1y$time,sqrt(avacc1y$av), col="green")
lines (avacc1z$time,sqrt(avacc1z$av), col="red")
axis(1, c(0.001, 0.01, 0.1, 0, 1, 10, 100, 1000, 10000))
axis(2, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
grid(equilogs=TRUE, lwd=1, col="orange")
title(main = "Allan variance Analysis", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation (m/s^2)")

legend(10, 1e-02, c("AccelerometerX", "AccelerometerY", "AccelerometerZ"),  fill = c("blue", "green", "red"))


########################## X Axis ##########################
xsensGyro <- ts (data.frame(values$gyrox), start=c(t1[1]), delta=mean(delta))

#Frequency
frequency (xsensGyro)
mean (values$gyrox)
sd (values$gyrox)

#### Calculating the Allan Variance for the gyroscope ####
avgyro1x <- avar (xsensGyro@.Data, frequency (xsensGyro))

########################## Y Axis ##########################
xsensGyro <- ts (data.frame(values$gyroy), start=c(t1[1]), delta=mean(delta))

#Frequency
frequency (xsensGyro)
mean (values$gyroy)
sd (values$gyroy)

#### Calculating the Allan Variance for the gyroscope ####
avgyro1y <- avar (xsensGyro@.Data, frequency (xsensGyro))


########################## Z Axis ##########################
xsensGyro <- ts (data.frame(values$gyroz), start=c(t1[1]), delta=mean(delta))

#Frequency
frequency (xsensGyro)
mean (values$gyroz)
sd (values$gyroz)

#### Calculating the Allan Variance for the gyroscope ####
avgyro1z <- avar (xsensGyro@.Data, frequency (xsensGyro))

#### Plotting the results ####
x11()
plot (avgyro1x$time,sqrt(avgyro1x$av),log= "xy", xaxt="n" , yaxt="n", type="l", col="blue", xlab="", ylab="")
lines (avgyro1y$time,sqrt(avgyro1y$av), col="green", lwd=1)
lines (avgyro1z$time,sqrt(avgyro1z$av), col="red")
axis(1, c(0.001, 0.01, 0.1, 0, 1, 10, 100, 1000, 10000, 100000))
axis(2, c(0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000, 10000))
grid(equilogs=TRUE, lwd=1, col="orange")
title(main = "Allan variance Analysis", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation (rad/s)")

legend(10, 5e-03, c("GyroscopeX", "GyroscopeY", "GyroscopeZ"),  fill = c("blue", "green", "red"))


plotav (avgyro1x)
##
#Random Walk, can be directly obtained by reading the slope line at T = 1
##
#Test 1 #  100Hz look at 57
approx (x=c(avgyro1x$time[7], avgyro1x$time[8]), y= c(sqrt(avgyro1x$av[7]), sqrt(avgyro1x$av[8])), n=100)
6.549113e-05

plotav (avgyro1y)
approx (x=c(avgyro1y$time[7], avgyro1y$time[8]), y= c(sqrt(avgyro1y$av[7]), sqrt(avgyro1y$av[8])), n=100)
8.019877e-05

plotav (avgyro1z)
approx (x=c(avgyro1z$time[7], avgyro1z$time[8]), y= c(sqrt(avgyro1z$av[7]), sqrt(avgyro1z$av[8])), n=100)
7.288081e-05


##
#Bias Instability, can be directly obtained by reading the 0 slope line 
##
#Bias instability Coefficient
sqrt(avgyro1x$av[9])/(sqrt(2*log(2)/pi))


#Rate Random Walk can be obtained by reading the allan variance value
#by a slope at +1/2. K=sqrt(3)*allandeviation(t)/sqrt(t)
# More or less at 20

sqrt(avgyro1x$av[12])/(sqrt(avgyro1x$time[12]/3))
sqrt(avgyro1y$av[12])/(sqrt(avgyro1y$time[12]/3))
sqrt(avgyro1z$av[12])/(sqrt(avgyro1z$time[12]/3))

## White noise simulation ##
sequ<- rnorm (2^16, mean=0, sd=1)#unit Gaussian distribution (sigma = 1)
sigma <- (6.04E-05/sqrt(0.01))
sequ <- sequ * sigma
wnoise <- ts (sequ, delta=0.01)

avwnoise <- avar (wnoise@.Data, frequency (wnoise))

plot (wnoise)

x11()
plot (avgyro1x$time,sqrt(avgyro1x$av),log= "xy", xaxt="n" , yaxt="n", type="l", col="blue", xlab="", ylab="")
lines(avwnoise$time, sqrt(avwnoise$av),type="l", col="green4", lwd=2)
axis(1, c(0.001, 0.01, 0.1, 0, 1, 10, 100, 1000, 10000, 100000))
axis(2, c(0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000, 10000))
grid(equilogs=TRUE, lwd=1, col="orange")

## Random walk simulation ##
sequ <- rep (0, 2^32)
At <- 0.01
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

x11()
plot (rwnoise)

plot(avrwnoise$time,sqrt(avrwnoise$av),log= "xy", xaxt="n" , yaxt="n", type="l", col="blue", xlab="", ylab="")
lines(avrwnoise$time, sqrt(avrwnoise$av), type="l", col="green")

plotav (avrwnoise)

sqrt(avrwnoise$av[7])/(sqrt(avrwnoise$time[3]/3))
sqrt(avrwnoise$av[8])/(sqrt(avrwnoise$time[4]/3))
sqrt(avrwnoise$av[7])/(sqrt(avrwnoise$time[5]/3))
sqrt(avrwnoise$av[8])/(sqrt(avrwnoise$time[6]/3))
sqrt(avrwnoise$av[7])/(sqrt(avrwnoise$time[7]/3))
sqrt(avrwnoise$av[8])/(sqrt(avrwnoise$time[8]/3))

## Addind both noises ##
noise <- wnoise + rwnoise
noise <- ts (noise, delta = 0.01)

avnoise <- avar (noise@.Data, frequency (noise))
x11()
plot(avnoise$time,sqrt(avnoise$av),log= "xy", xaxt="n" , yaxt="n", type="l", col="blue", xlab="", ylab="")

x11()
plot (avgyro1x$time,sqrt(avgyro1x$av),log= "xy", xaxt="n" , yaxt="n", type="l", col="blue", xlab="", ylab="")
lines(avnoise$time, sqrt(avnoise$av),type="l", col="green4", lwd=2)
axis(1, c(0.001, 0.01, 0.1, 0, 1, 10, 100, 1000, 10000, 100000))
axis(2, c(0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000, 10000))
grid(equilogs=TRUE, lwd=1, col="orange")
