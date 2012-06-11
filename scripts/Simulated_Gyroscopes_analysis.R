#libraries
library(stats)


#########
# TEST 1#
#########
rm(values)
values <- read.csv ("imt30_model.csv")



names(values) = c("time", "gyroz")
names(values) = c("time", "gyrox", "gyroy", "gyroz")
#names(values) = c("time","inclz", "inclx", "incly")

delta = values$time[2]- values$time[1]

#### Calculating the Allan Variance for the inclinometers ####

#Incl X
imt30Gyros <- ts (data.frame(values$gyrox), start=c(values$time[1]), delta=mean(delta))
avsimugyro1x <- avar (imt30Gyros@.Data, frequency (imt30Gyros))

#frequency (imt30Gyros)

#Incl Y
imt30Gyros <- ts (data.frame(values$gyroy), start=c(values$time[1]), delta=mean(delta))
avsimugyro1y <- avar (imt30Gyros@.Data, frequency (imt30Gyros))

#Incl Z
imt30Gyros <- ts (data.frame(values$gyroz), start=c(values$time[1]), delta=mean(delta))
avsimugyro1z <- avar (imt30Gyros@.Data, frequency (imt30Gyros))


windows()
plot (imt30Gyros, xlab="", ylab="")
title(main = "Gyroscopes X Values Position 1", ylab = "Angular velocity (rad/s)", xlab = "Time(s)")
lines (rep(mean(imt30Gyros), length(imt30Gyros)), xlab="", ylab="", col="red")
lines (rep(mean(imt30Gyros)-sd(imt30Gyros), length(imt30Gyros)), xlab="", ylab="", col="blue")
lines (rep(mean(imt30Gyros)+sd(imt30Gyros), length(imt30Gyros)), xlab="", ylab="", col="blue")


#### Plotting the results ####
windows()
lines(avsimugyro1x$time,sqrt(avsimugyro1x$av),lwd=2,log= "xy", xaxt="n" , yaxt="n", type="l", col="black", xlab="", ylab="")
lines (avsimugyro1y$time,sqrt(avsimugyro1y$av), col="black",lwd=2)
lines (avsimugyro1z$time,sqrt(avsimugyro1z$av), col="black", lwd=2)
axis(1, c(0.00000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000, 10000, 100000))
axis(2, c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
grid(equilogs=TRUE, lwd=1, col="orange")
title(main = "Allan variance Analysis", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation (m/s)")

legend(100, 3e-01, c("GyroscopeX", "GyroscopeY", "GyroscopeZ"),  fill = c("blue", "green", "red"))

#Comparation of the results
windows()
plot (avgyro3x$time[1:length(avgyro3x$time)-1],sqrt(avgyro3x$av[1:length(avgyro3x$av)-1]),log= "xy", xaxt="n" , yaxt="n", col="blue", xlab="", ylab="")
lines (avgyro3x$time[1:length(avgyro3x$time)-1],sqrt(avgyro3x$av[1:length(avgyro3x$av)-1]), col="blue")
lines(avsimugyro1x$time,sqrt(avsimugyro1x$av),lwd=3,log= "xy", xaxt="n" , yaxt="n", type="l", col="black", xlab="", ylab="")
axis(1, c(0.00000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000, 10000, 100000))
axis(2, c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
grid(equilogs=TRUE, lwd=1, col="orange")

title(main = "Allan variance Analysis", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation (rad/s)")

legend(0.01, 1e-03, c("Gyroscope(real values)", "Simulated Gyroscope (whitenoise+randomwalk+biasinst)"),  fill = c("blue", "black"))


#min value
min(sqrt(avsimugyro1x$av[1:length(avsimugyro1x$av)-1]))

min(sqrt(avgyro5x$av[1:length(avgyro5x$av)-1]))-min(sqrt(avsimugyro1x$av[1:length(avsimugyro1x$av)-1]))


#Delete the variables for a new test
rm (imt30Incl, values)





