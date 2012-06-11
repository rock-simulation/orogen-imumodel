###############################################
# Post Processing R Script to compute         #
# the Allan variance over the Gyroscopes      #
# data recorded in the Test Center            #
# Javier Hidalgo Carrió                       #
# TEC-MMA  February 2010                      #
###############################################

#libraries
#Need packages(tseries, stats, ast and lmtest)
library(stats)

# You have to load the Allan Variance function (avar)

# You may change this path
#setwd("Z:/jhidalgo/AllanVariable/Data/TestCenter/Gyroscopes")

setwd("/Users/javi/Misdatos/esa/Traineeship/TEC-MMA/jhidalgo/Characterization_Inertial_Sensors/Allan_Variance/Characterization&Modeling_IMT30/Data/TestCenter/Gyroscopes")

#############
# TEST 1#
#########
values <- read.table ("24062010-194634-x.data", sep="\t")

names(values) = c("time", "gyrox")

delta = 1/976

imt30Gyros <- ts (data.frame(values$gyrox), start=c(values$time[1]), delta=mean(delta))

#### Calculating the Allan Variance for the gyroscopes ####
#imt30Gyros@.Data[1,1]
frequency (imt30Gyros)

#Incl X
avgyro1x <- avar (imt30Gyros@.Data, frequency (imt30Gyros))
#Delete previous values
rm (imt30Gyros, values)

#Incl Y
values <- read.table ("24062010-194634-y.data", sep="\t")

names(values) = c("time", "gyroy")

imt30Gyros <- ts (data.frame(values$gyroy), start=c(values$time[1]), delta=mean(delta))
avgyro1y <- avar (imt30Gyros@.Data, frequency (imt30Gyros))
#Delete previous values
rm (imt30Gyros, values)

#Incl Z
values <- read.table ("24062010-194634-z.data", sep="\t")

names(values) = c("time", "gyroz")

imt30Gyros<- ts (data.frame(values$gyroz), start=c(values$time[1]), delta=mean(delta))
avgyro1z <- avar (imt30Gyros@.Data, frequency (imt30Gyros))
#Delete previous values
rm (imt30Gyros, values)


#### Plotting the results ####
windows()
plot (avgyro1x$time,sqrt(avgyro1x$av),log= "xy", xaxt="n" , yaxt="n", type="l", col="blue", xlab="", ylab="")
lines (avgyro1y$time,sqrt(avgyro1y$av), col="green", lwd=1)
lines (avgyro1z$time,sqrt(avgyro1z$av), col="red")
axis(1, c(0.001, 0.01, 0.1, 0, 1, 10, 100, 1000, 10000, 100000))
axis(2, c(0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000, 10000))
grid(equilogs=TRUE, lwd=1, col="orange")
title(main = "Allan variance Analysis", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation (rad/s)")

legend(1, 1e-03, c("GyroscopeX", "GyroscopeY", "GyroscopeZ"),  fill = c("blue", "green", "red"))


#########
# TEST 2#
#########
values <- read.table ("26062010-191649-x.data", sep="\t")

names(values) = c("time", "gyrox")

delta = 1/976

imt30Gyros <- ts (data.frame(values$gyrox), start=c(values$time[1]), delta=mean(delta))

#### Calculating the Allan Variance for the gyroscopes ####
#imt30Gyros@.Data[1,1]
frequency (imt30Gyros)

#Incl X
avgyro2x <- avar (imt30Gyros@.Data, frequency (imt30Gyros))
#Delete previous values
rm (imt30Gyros, values)

#Incl Y
values <- read.table ("26062010-191649-y.data", sep="\t")

names(values) = c("time", "gyroy")

imt30Gyros <- ts (data.frame(values$gyroy), start=c(values$time[1]), delta=mean(delta))
avgyro2y <- avar (imt30Gyros@.Data, frequency (imt30Gyros))
#Delete previous values
rm (imt30Gyros, values)

#Incl Z
values <- read.table ("26062010-191649-z.data", sep="\t")

names(values) = c("time", "gyroz")

imt30Gyros<- ts (data.frame(values$gyroz), start=c(values$time[1]), delta=mean(delta))
avgyro2z <- avar (imt30Gyros@.Data, frequency (imt30Gyros))
#Delete previous values
rm (imt30Gyros, values)


#### Plotting the results ####
windows()
plot (avgyro2x$time,sqrt(avgyro2x$av),log= "xy", xaxt="n" , yaxt="n", type="l", col="blue", xlab="", ylab="")
lines (avgyro2y$time,sqrt(avgyro2y$av), col="green", lwd=2)
lines (avgyro2z$time,sqrt(avgyro2z$av), col="red")
axis(1, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
axis(2, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
grid(equilogs=TRUE, lwd=1, col="orange")
title(main = "Allan variance Analysis", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation (rad/s)")

legend(1, 1e-03, c("GyroscopeX", "GyroscopeY", "GyroscopeZ"),  fill = c("blue", "green", "red"))


#########
# TEST 3#
#########
values <- read.table ("29062010-202049-x.data", sep="\t")

names(values) = c("time", "gyrox")

delta = 1/976

imt30Gyros <- ts (data.frame(values$gyrox), start=c(values$time[1]), delta=mean(delta))

#### Calculating the Allan Variance for the gyroscopes ####
#imt30Gyros@.Data[1,1]
frequency (imt30Gyros)

#Incl X
avgyro3x <- avar (imt30Gyros@.Data, frequency (imt30Gyros))
#Delete previous values
rm (imt30Gyros, values)

#Incl Y
values <- read.table ("29062010-202049-y.data", sep="\t")

names(values) = c("time", "gyroy")

imt30Gyros <- ts (data.frame(values$gyroy), start=c(values$time[1]), delta=mean(delta))
avgyro3y <- avar (imt30Gyros@.Data, frequency (imt30Gyros))
#Delete previous values
rm (imt30Gyros, values)

#Incl Z
values <- read.table ("29062010-202049-z.data", sep="\t")

names(values) = c("time", "gyroz")

imt30Gyros<- ts (data.frame(values$gyroz), start=c(values$time[1]), delta=mean(delta))
avgyro3z <- avar (imt30Gyros@.Data, frequency (imt30Gyros))
#Delete previous values
rm (imt30Gyros, values)


#### Plotting the results ####
windows()
plot (avgyro3x$time,sqrt(avgyro3x$av),log= "xy", xaxt="n" , yaxt="n", type="l", col="blue", xlab="", ylab="")
lines (avgyro3y$time,sqrt(avgyro3y$av), col="green", lwd=2)
lines (avgyro3z$time,sqrt(avgyro3z$av), col="red")
axis(1, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
axis(2, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
grid(equilogs=TRUE, lwd=1, col="orange")
title(main = "Allan variance Analysis", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation (rad/s)")

legend(1, 1e-03, c("GyroscopeX", "GyroscopeY", "GyroscopeZ"),  fill = c("blue", "green", "red"))


#########
# TEST 4#
#########
values <- read.table ("30062010-211701-x.data", sep="\t")

names(values) = c("time", "gyrox")

delta = 1/976

imt30Gyros <- ts (data.frame(values$gyrox), start=c(values$time[1]), delta=mean(delta))

#### Calculating the Allan Variance for the gyroscopes ####
#imt30Gyros@.Data[1,1]
frequency (imt30Gyros)

#Incl X
avgyro4x <- avar (imt30Gyros@.Data, frequency (imt30Gyros))
#Delete previous values
rm (imt30Gyros, values)

#Incl Y
values <- read.table ("30062010-211701-y.data", sep="\t")

names(values) = c("time", "gyroy")

imt30Gyros <- ts (data.frame(values$gyroy), start=c(values$time[1]), delta=mean(delta))
avgyro4y <- avar (imt30Gyros@.Data, frequency (imt30Gyros))
#Delete previous values
rm (imt30Gyros, values)

#Incl Z
values <- read.table ("30062010-211701-z.data", sep="\t")

names(values) = c("time", "gyroz")

imt30Gyros<- ts (data.frame(values$gyroz), start=c(values$time[1]), delta=mean(delta))
avgyro4z <- avar (imt30Gyros@.Data, frequency (imt30Gyros))
#Delete previous values
rm (imt30Gyros, values)


#### Plotting the results ####
windows()
plot (avgyro4x$time,sqrt(avgyro4x$av),log= "xy", xaxt="n" , yaxt="n", type="l", col="blue", xlab="", ylab="")
lines (avgyro4y$time,sqrt(avgyro4y$av), col="green", lwd=2)
lines (avgyro4z$time,sqrt(avgyro4z$av), col="red")
axis(1, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
axis(2, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
grid(equilogs=TRUE, lwd=1, col="orange")
title(main = "Allan variance Analysis", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation (rad/s)")

legend(1, 1e-03, c("GyroscopeX", "GyroscopeY", "GyroscopeZ"),  fill = c("blue", "green", "red"))


#########
# TEST 5#
#########
values <- read.table ("01072010-193849-x.data", sep="\t")

names(values) = c("time", "gyrox")

delta = 1/976

imt30Gyros <- ts (data.frame(values$gyrox), start=c(values$time[1]), delta=mean(delta))

#### Calculating the Allan Variance for the gyroscopes ####
#imt30Gyros@.Data[1,1]
frequency (imt30Gyros)

#Incl X
avgyro5x <- avar (imt30Gyros@.Data, frequency (imt30Gyros))
#Delete previous values
rm (imt30Gyros, values)

#Incl Y
values <- read.table ("01072010-193849-y.data", sep="\t")

names(values) = c("time", "gyroy")

imt30Gyros <- ts (data.frame(values$gyroy), start=c(values$time[1]), delta=mean(delta))
avgyro5y <- avar (imt30Gyros@.Data, frequency (imt30Gyros))
#Delete previous values
rm (imt30Gyros, values)

#Incl Z
values <- read.table ("01072010-193849-z.data", sep="\t")

names(values) = c("time", "gyroz")

imt30Gyros<- ts (data.frame(values$gyroz), start=c(values$time[1]), delta=mean(delta))
avgyro5z <- avar (imt30Gyros@.Data, frequency (imt30Gyros))
#Delete previous values
rm (imt30Gyros, values)


#### Plotting the results ####
windows()
plot (avgyro5x$time,sqrt(avgyro5x$av),log= "xy", xaxt="n" , yaxt="n", type="l", col="blue", xlab="", ylab="")
lines (avgyro5y$time,sqrt(avgyro5y$av), col="green", lwd=2)
lines (avgyro5z$time,sqrt(avgyro5z$av), col="red")
axis(1, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
axis(2, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
grid(equilogs=TRUE, lwd=1, col="orange")
title(main = "Allan variance Analysis", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation (rad/s)")

legend(1, 1e-03, c("GyroscopeX", "GyroscopeY", "GyroscopeZ"),  fill = c("blue", "green", "red"))




###############################
#Calculate the values (Axis X)#
###############################

windows()
plot (avgyro1x$time[1:length(avgyro1x$time)-1],sqrt(avgyro1x$av[1:length(avgyro1x$av)-1]),log= "xy", xaxt="n" , yaxt="n", col="blue", xlab="", ylab="")
lines (avgyro1x$time[1:length(avgyro1x$time)-1],sqrt(avgyro1x$av[1:length(avgyro1x$av)-1]), col="blue")
lines (avgyro2x$time[1:length(avgyro2x$time)-1],sqrt(avgyro2x$av[1:length(avgyro2x$av)-1]), col="green")
lines (avgyro3x$time[1:length(avgyro3x$time)-1],sqrt(avgyro3x$av[1:length(avgyro3x$av)-1]), col="red")
lines (avgyro4x$time[1:length(avgyro4x$time)-1],sqrt(avgyro4x$av[1:length(avgyro4x$av)-1]), col="violet")
lines (avgyro5x$time[1:length(avgyro5x$time)-1],sqrt(avgyro5x$av[1:length(avgyro5x$av)-1]), col="orange")

axis(1, c(0.01, 0.1, 0, 1, 10, 100, 1000, 10000))
axis(2, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
grid(equilogs=TRUE, lwd=1, col="orange")
title(main = "Allan variance Analysis for Gyroscopes", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation (rad/s^2)")


#legend(10, 5e-04, c("Axis X (Test1)", "Axis X (Test2)", "Axis X (Test3)", "Axis X (Test4)", "Axis X (Test5)", "Simulated X axis"),  fill = c("blue", "green", "red", "violet", "orange", "black"))
legend(10, 5e-04, c("Axis X (Test1)", "Axis X (Test2)", "Axis X (Test3)", "Axis X (Test4)", "Axis X (Test5)"),  fill = c("blue", "green", "red", "violet", "orange"))



##
#Random Walk, can be directly obtained by reading the slope line at T = 1
##
#Test 1 # at 977Hz look at 91
approx (x=c(avgyro1x$time[10], avgyro1x$time[11]), y= c(sqrt(avgyro1x$av[10]), sqrt(avgyro1x$av[11])), n=100)

#Test 2 # at 977Hz look at 91
approx (x=c(avgyro2x$time[10], avgyro2x$time[11]), y= c(sqrt(avgyro2x$av[10]), sqrt(avgyro2x$av[11])), n=100)

#Test 3 # at 977Hz look at 91
approx (x=c(avgyro3x$time[10], avgyro3x$time[11]), y= c(sqrt(avgyro3x$av[10]), sqrt(avgyro3x$av[11])), n=100)

#Test 4 # at 977Hz look at 91
approx (x=c(avgyro4x$time[10], avgyro4x$time[11]), y= c(sqrt(avgyro4x$av[10]), sqrt(avgyro4x$av[11])), n=100)

#Test 5 # at 977Hz look at 91
approx (x=c(avgyro5x$time[10], avgyro5x$time[11]), y= c(sqrt(avgyro5x$av[10]), sqrt(avgyro5x$av[11])), n=100)

##
#Bias Instability, can be directly obtained by reading the 0 slope line 
##

#Bias instability Coefficient
sqrt(avgyro1x$av[16])/(sqrt(2*log(2)/pi))
sqrt(avgyro2x$av[16])/(sqrt(2*log(2)/pi))
sqrt(avgyro3x$av[16])/(sqrt(2*log(2)/pi))
sqrt(avgyro4x$av[16])/(sqrt(2*log(2)/pi))
sqrt(avgyro5x$av[16])/(sqrt(2*log(2)/pi))

#Bias instability correlation time
1/avgyro1x$time[18]
avgyro2x$time[16]/1.89
avgyro3x$time[16]/1.89
avgyro4x$time[16]/1.89
avgyro5x$time[16]/1.89


#Rate Random Walk can be obtained by reading the allan variance value
#by a slope at +1/2. K=sqrt(3)*allandeviation(t)/sqrt(t)
# More or less at 20

sqrt(avgyro1x$av[18])/(sqrt(avgyro1x$time[18]/3))
sqrt(avgyro2x$av[18])/(sqrt(avgyro2x$time[18]/3))
sqrt(avgyro3x$av[18])/(sqrt(avgyro3x$time[18]/3))
sqrt(avgyro4x$av[18])/(sqrt(avgyro4x$time[18]/3))
sqrt(avgyro5x$av[18])/(sqrt(avgyro5x$time[18]/3))

#Mean value of the allanvariance at the coefficiente K
(sqrt(avgyro1x$av[21])+sqrt(avgyro2x$av[21])+
sqrt(avgyro3x$av[21])+sqrt(avgyro4x$av[21])+sqrt(avgyro5x$av[22]))/5


###############################
#Calculate the values (Axis Y)#
###############################

windows()
plot (avgyro1y$time[1:length(avgyro1y$time)-1],sqrt(avgyro1y$av[1:length(avgyro1y$av)-1]),log= "xy", xaxt="n" , yaxt="n", col="blue", xlab="", ylab="")
lines (avgyro1y$time[1:length(avgyro1y$time)-1],sqrt(avgyro1y$av[1:length(avgyro1y$av)-1]), col="blue")
lines (avgyro2y$time[1:length(avgyro2y$time)-1],sqrt(avgyro2y$av[1:length(avgyro2y$av)-1]), col="green")
lines (avgyro3y$time[1:length(avgyro3y$time)-1],sqrt(avgyro3y$av[1:length(avgyro3y$av)-1]), col="red")
lines (avgyro4y$time[1:length(avgyro4y$time)-1],sqrt(avgyro4y$av[1:length(avgyro4y$av)-1]), col="violet")
lines (avgyro5y$time[1:length(avgyro5y$time)-1],sqrt(avgyro5y$av[1:length(avgyro5y$av)-1]), col="orange")
axis(1, c(0.01, 0.1, 0, 1, 10, 100, 1000, 10000))
axis(2, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
grid(equilogs=TRUE, lwd=1, col="orange")
title(main = "Allan variance Analysis for Gyroscopes", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation (rad/s^2)")


#legend(10, 5e-04, c("Axis Y (Test1)", "Axis Y (Test2)", "Axis Y (Test3)", "Axis Y (Test4)", "Axis Y (Test5)", "Simulated Y axis"),  fill = c("blue", "green", "red", "violet", "orange", "black"))

legend(10, 5e-04, c("Axis Y (Test1)", "Axis Y (Test2)", "Axis Y (Test3)", "Axis Y (Test4)", "Axis Y (Test5)"),  fill = c("blue", "green", "red", "violet", "orange"))



##
#Random Walk, can be directly obtained by reading the slope line at T = 1
##

#Test 1 # at 500Hz look at 91
approx (x=c(avgyro1y$time[10], avgyro1y$time[11]), y= c(sqrt(avgyro1y$av[10]), sqrt(avgyro1y$av[11])), n=100)

#Test 2 # at 500Hz look at 91
approx (x=c(avgyro2y$time[10], avgyro2y$time[11]), y= c(sqrt(avgyro2y$av[10]), sqrt(avgyro2y$av[11])), n=100)

#Test 3 # at 500Hz look at 91
approx (x=c(avgyro3y$time[10], avgyro3y$time[11]), y= c(sqrt(avgyro3y$av[10]), sqrt(avgyro3y$av[11])), n=100)

#Test 4 # at 977Hz look at 91
approx (x=c(avgyro4y$time[10], avgyro4y$time[11]), y= c(sqrt(avgyro4y$av[10]), sqrt(avgyro4y$av[11])), n=100)

#Test 5 # at 977Hz look at 91
approx (x=c(avgyro5y$time[10], avgyro5y$time[11]), y= c(sqrt(avgyro5y$av[10]), sqrt(avgyro5y$av[11])), n=100)


##
#Bias Instability, can be directly obtained by reading the 0 slope line 
##
#Bias instability

#Bias instability Coefficient
sqrt(avgyro1y$av[16])/(sqrt(2*log(2)/pi))
sqrt(avgyro2y$av[16])/(sqrt(2*log(2)/pi))
sqrt(avgyro3y$av[16])/(sqrt(2*log(2)/pi))
sqrt(avgyro4y$av[16])/(sqrt(2*log(2)/pi))
sqrt(avgyro5y$av[16])/(sqrt(2*log(2)/pi))


#Rate Random Walk can be obtained by reading the allan variance value
#by a slope at +1/2. K=sqrt(3)*allandeviation(t)/sqrt(t)
# More or less at 21

(sqrt(3)*sqrt(avgyro1y$av[18]))/sqrt(avgyro1y$time[18])
(sqrt(3)*sqrt(avgyro2y$av[18]))/sqrt(avgyro2y$time[18])
(sqrt(3)*sqrt(avgyro3y$av[18]))/sqrt(avgyro3y$time[18])
(sqrt(3)*sqrt(avgyro4y$av[18]))/sqrt(avgyro4y$time[18])
(sqrt(3)*sqrt(avgyro5y$av[18]))/sqrt(avgyro5y$time[18])

#Mean value of the allanvariance at the coefficiente K
(sqrt(avgyro1y$av[21])+sqrt(avgyro2y$av[21])+
sqrt(avgyro3y$av[21])+sqrt(avgyro4y$av[21])+sqrt(avgyro5y$av[22]))/5



###############################
#Calculate the values (Axis Z)#
###############################

windows()
plot (avgyro1z$time[1:length(avgyro1z$time)-1],sqrt(avgyro1z$av[1:length(avgyro1z$av)-1]),log= "xy", xaxt="n" , yaxt="n", col="blue", xlab="", ylab="")
lines (avgyro1z$time[1:length(avgyro1z$time)-1],sqrt(avgyro1z$av[1:length(avgyro1z$av)-1]), col="blue")
lines (avgyro2z$time[1:length(avgyro2z$time)-1],sqrt(avgyro2z$av[1:length(avgyro2z$av)-1]), col="green")
lines (avgyro3z$time[1:length(avgyro3z$time)-1],sqrt(avgyro3z$av[1:length(avgyro3z$av)-1]), col="red")
lines (avgyro4z$time[1:length(avgyro4z$time)-1],sqrt(avgyro4z$av[1:length(avgyro4z$av)-1]), col="violet")
lines (avgyro5z$time[1:length(avgyro5z$time)-1],sqrt(avgyro5z$av[1:length(avgyro5z$av)-1]), col="orange")
axis(1, c(0.01, 0.1, 0, 1, 10, 100, 1000, 10000))
axis(2, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
grid(equilogs=TRUE, lwd=1, col="orange")
title(main = "Allan variance Analysis for Gyroscopes", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation (rad/s^2)")

#legend(10, 5e-04, c("Axis Z (Test1)", "Axis Z (Test2)", "Axis Z (Test3)", "Axis Z (Test4)", "Axis Z (Test5)", "Simulated Z axis"),  fill = c("blue", "green", "red", "violet", "orange", "black"))

legend(10, 5e-04, c("Axis Z (Test1)", "Axis Z (Test2)", "Axis Z (Test3)", "Axis Z (Test4)", "Axis Z (Test5)"),  fill = c("blue", "green", "red", "violet", "orange"))


##
#Random Walk, can be directly obtained by reading the slope line at T = 1
##

#Test 1 # at 500Hz look at 91
approx (x=c(avgyro1z$time[10], avgyro1z$time[11]), y= c(sqrt(avgyro1z$av[10]), sqrt(avgyro1z$av[11])), n=100)

#Test 2 # at 500Hz look at 91
approx (x=c(avgyro2z$time[10], avgyro2z$time[11]), y= c(sqrt(avgyro2z$av[10]), sqrt(avgyro2z$av[11])), n=100)

#Test 3 # at 500Hz look at 91
approx (x=c(avgyro3z$time[10], avgyro3z$time[11]), y= c(sqrt(avgyro3z$av[10]), sqrt(avgyro3z$av[11])), n=100)

#Test 4 # at 977Hz look at 91
approx (x=c(avgyro4z$time[10], avgyro4z$time[11]), y= c(sqrt(avgyro4z$av[10]), sqrt(avgyro4z$av[11])), n=100)

#Test 5 # at 977Hz look at 91
approx (x=c(avgyro5z$time[10], avgyro5z$time[11]), y= c(sqrt(avgyro5z$av[10]), sqrt(avgyro5z$av[11])), n=100)

#Bias instability Coefficient
sqrt(avgyro1z$av[16])/(sqrt(2*log(2)/pi))
sqrt(avgyro2z$av[16])/(sqrt(2*log(2)/pi))
sqrt(avgyro3z$av[16])/(sqrt(2*log(2)/pi))
sqrt(avgyro4z$av[16])/(sqrt(2*log(2)/pi))
sqrt(avgyro5z$av[16])/(sqrt(2*log(2)/pi))


#Rate Random Walk can be obtained by reading the allan variance value
#by a slope at +1/2. K=sqrt(3)*allandeviation(t)/sqrt(t)
# More or less at 21

(sqrt(3)*sqrt(avgyro1z$av[19]))/sqrt(avgyro1z$time[19])
(sqrt(3)*sqrt(avgyro2z$av[19]))/sqrt(avgyro2z$time[19])
(sqrt(3)*sqrt(avgyro3z$av[19]))/sqrt(avgyro3z$time[19])
(sqrt(3)*sqrt(avgyro4z$av[19]))/sqrt(avgyro4z$time[19])
(sqrt(3)*sqrt(avgyro5z$av[19]))/sqrt(avgyro5z$time[19])

#Mean value of the allanvariance at the coefficiente K
(sqrt(avgyro1z$av[21])+sqrt(avgyro2z$av[21])+
sqrt(avgyro3z$av[21])+sqrt(avgyro4z$av[21])+sqrt(avgyro5z$av[22]))/5

##################
# Filtering
# Low-pass at 10Hz
##################
values <- read.table ("977Hz/moreinstables/30062010-211701-z.data", sep="\t", nrows=1576800)
names(values) = c("time","gyrosz")

delta = 1/976

#First allan variance of the original one(bandwidth 100Hz, sampling rate 976Hz)
imt30Gyros <- ts (data.frame(values$gyrosz), start=c(values$time[1]), delta=mean(delta))
avgyro4z <- avar (imt30Gyros@.Data, frequency (imt30Gyros))


#Low-pass filter 10Hz
#gyrosz = filtervalues$gyrosz
gyrosz = values$gyrosz
#time = filtervalues$time
timez = values$time

Tc=1/10#10Hz

filtervalues=NULL
filtervalues$gyrosz[1] = gyrosz[1];
for (i in 2: length(gyrosz))
{
	filtervalues$gyrosz[i]=filtervalues$gyrosz[i-1]+(delta/Tc)*(gyrosz[i]-filtervalues$gyrosz[i-1])
}

filtervalues$time = timez

#Allan Variance of filtered values
imt30Gyros <- ts (data.frame(filtervalues$gyrosz), start=c(filtervalues$time[1]), delta=mean(delta))
avgyro6z <- avar (imt30Gyros@.Data, frequency (imt30Gyros))

#Taking 30Hz
gyrosz = filtervalues$gyrosz
timez = filtervalues$time

filtervalues=NULL
time = 0
i = 1
k = 1
while (i<=(length(gyrosz)-(33)))
{
	filtervalues$gyrosz[k] = gyrosz[(i)]
	filtervalues$time[k] = timez[i]

	i <- i + 33
	k <- k + 1
}
delta = 1/30

imt30Gyros <- ts (data.frame(filtervalues$gyrosz), start=c(filtervalues$time[1]), delta=mean(delta))
avgyro7z <- avar (imt30Gyros@.Data, frequency (imt30Gyros))



#Plot the results
windows()
plot (avgyro4z$time[1:length(avgyro4z$time)-1],sqrt(avgyro4z$av[1:length(avgyro4z$av)-1]),log= "xy", xaxt="n" , yaxt="n", col="blue", xlab="", ylab="")
lines (avgyro4z$time[1:length(avgyro4z$time)-1],sqrt(avgyro4z$av[1:length(avgyro4z$av)-1]), col="blue")
points (avgyro6z$time[1:length(avgyro6z$time)-1],sqrt(avgyro6z$av[1:length(avgyro6z$av)-1]), col="red")
lines (avgyro6z$time[1:length(avgyro6z$time)-1],sqrt(avgyro6z$av[1:length(avgyro6z$av)-1]), col="red")
points (avgyro7z$time[1:length(avgyro7z$time)-1],sqrt(avgyro7z$av[1:length(avgyro7z$av)-1]), col="green")
lines (avgyro7z$time[1:length(avgyro7z$time)-1],sqrt(avgyro7z$av[1:length(avgyro7z$av)-1]), col="green")

axis(1, c(0.01, 0.1, 0, 1, 10, 100, 1000, 10000))
axis(2, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
grid(equilogs=TRUE, lwd=1, col="orange")

legend(0.03, 1e-03, c("Gyros(100HzBandwidth)(977Hz sample rate)", "Filtered Gyros(10HzBandwicth)(977Hz sample rate)", "Filtered Gyros(10HzBandwicth)(30Hz sample rate)"),  fill = c("blue", "red", "green"))
title(main = "Allan variance Analysis for Gyroscopes", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation (rad/s^2)")


#Averaging 30Hz(Optional)
filtervalues=NULL
time = 0
i = 1
k = 1
while (i<=(length(values$gyrosz)-(33)))
{
	filtervalues$gyrosz[k] = sum(values$gyrosz[i:(i+(33))])/33
	filtervalues$time[k] = values$time[i]

	i <- i + 33
	k <- k + 1
}

#Allan Variance of averaging values
delta = 1/30
imt30Gyros <- ts (data.frame(filtervalues$gyrosz), start=c(filtervalues$time[1]), delta=mean(delta))
avgyro4zbis <- avar (imt30Gyros@.Data, frequency (imt30Gyros))

points (avgyro4zbis$time[1:length(avgyro4zbis$time)-1],sqrt(avgyro4zbis$av[1:length(avgyro4zbis$av)-1]),log= "xy", xaxt="n" , yaxt="n", col="green", xlab="", ylab="")
lines (avgyro4zbis$time[1:length(avgyro4zbis$time)-1],sqrt(avgyro4zbis$av[1:length(avgyro4zbis$av)-1]), col="green")


##################
# Filtering
# Averaging at 10Hz
##################

values <- read.table ("977Hz/moreinstables/30062010-211701-z.data", sep="\t")
#values <- read.table ("976Hz/20100612-142840-imegogyros.data", sep="\t")

names(values) = c("time","gyrosz")
#names(values) = c("gyrosx","gyrosy", "gyrosz")


delta = 1/977

filtervalues=NULL
time = 0
i = 1
k = 1
while (i<=(length(values$gyrosz)-(98)))
{
	filtervalues$gyrosz[k] = sum(values$gyrosz[i:(i+(98))])/98
	filtervalues$time[k] = values$time[i]

	i <- i + 98
	k <- k + 1
}

windows()
plot (values$time[1:(10*3400)], values$gyrosz[1:(10*3400)], type="l", xlab="", ylab="")
lines (filtervalues$time[1:((10*3400))], filtervalues$gyrosz[1:((10*3400))], col="green")
title(main = "Gyroscopes values - average filtering", xlab = "Times (Sec)", ylab = "Angular velocity(rad/s^2)")
lines (rep(mean(values$gyrosz[1:(10*3400)]), length(values$gyrosz[1:(10*3400)])), xlab="", ylab="", col="red")
lines (rep(mean(values$gyrosz[1:(10*3400)])-sd(values$gyrosz[1:(10*3400)]), length(values$gyrosz[1:(10*3400)])), xlab="", ylab="", col="blue")
lines (rep(mean(values$gyrosz[1:(10*3400)])+sd(values$gyrosz[1:(10*3400)]), length(values$gyrosz[1:(10*3400)])), xlab="", ylab="", col="blue")
lines (rep(mean(filtervalues$gyrosz[1:((10*3400)/98)]), length(filtervalues$gyrosz[1:(10*3400)])), xlab="", ylab="", col="red")
lines (rep(mean(filtervalues$gyrosz[1:((10*3400)/98)])-sd(filtervalues$gyrosz[1:((10*3400)/98)]), length(filtervalues$gyrosz[1:((10*3400)/98)])), xlab="", ylab="", col="blue")
lines (rep(mean(filtervalues$gyrosz[1:((10*3400)/98)])+sd(filtervalues$gyrosz[1:((10*3400)/98)]), length(filtervalues$gyrosz[1:((10*3400)/98)])), xlab="", ylab="", col="blue")

windows()
plot (filtervalues$time, filtervalues$gyrosz, type="l", xlab="", ylab="")




windows()
plot (values$time, values$gyrosx, type="l", )
lines (rep(mean(values$gyrosx), length(values$gyrosx)), xlab="", ylab="", col="red")
lines (rep(mean(values$gyrosx)-sd(values$gyrosx), length(values$gyrosx)), xlab="", ylab="", col="blue")
lines (rep(mean(values$gyrosx)+sd(values$gyrosx), length(values$gyrosx)), xlab="", ylab="", col="blue")
lines (filtervalues$time, filtervalues$gyrosx, col="blue")


#compute Allan Variance
delta = 1/10 #1/10Hz

imt30Gyros <- ts (data.frame(filtervalues$gyrosz), start=c(filtervalues$time[1]), delta=mean(delta))

#### Calculating the Allan Variance for the gyroscopes ####
#imt30Gyros@.Data[1,1]
frequency (imt30Gyros)

#Gyros filter Z Axis
avgyro6z <- avar (imt30Gyros@.Data, frequency (imt30Gyros))
#Delete previous values
rm (imt30Gyros, values)

windows()
plot (avgyro4z$time[1:length(avgyro4z$time)-1],sqrt(avgyro4z$av[1:length(avgyro4z$av)-1]),log= "xy", xaxt="n" , yaxt="n", col="blue", xlab="", ylab="")
lines (avgyro4z$time[1:length(avgyro4z$time)-1],sqrt(avgyro4z$av[1:length(avgyro4z$av)-1]), col="blue")
points (avgyro6z$time[1:length(avgyro6z$time)-1],sqrt(avgyro6z$av[1:length(avgyro6z$av)-1]), col="red")
lines (avgyro6z$time[1:length(avgyro6z$time)-1],sqrt(avgyro6z$av[1:length(avgyro6z$av)-1]), col="red")
axis(1, c(0.01, 0.1, 0, 1, 10, 100, 1000, 10000))
axis(2, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
grid(equilogs=TRUE, lwd=1, col="orange")

legend(1, 5e-04, c("GyrosValues 977Hz", "GyrosValues Aver at 10Hz"),  fill = c("blue", "red"))
title(main = "Allan variance Analysis for Gyroscopes", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation (rad/s^2)")

#Random Walk, can be directly obtained by reading the slope line at T = 1
##
#Test 4 # at 977Hz look at 91
approx (x=c(avgyro4z$time[10], avgyro4z$time[11]), y= c(sqrt(avgyro4z$av[10]), sqrt(avgyro4z$av[11])), n=100)
#Test 6 # at 10Hz look at 26
approx (x=c(avgyro6z$time[4], avgyro6z$time[5]), y= c(sqrt(avgyro6z$av[4]), sqrt(avgyro6z$av[5])), n=100)

whitenoise = 5.950260e-05

#Rate Random Walk can be obtained by reading the allan variance value
#by a slope at +1/2. K=sqrt(3)*allandeviation(t)/sqrt(t)
# More or less at 11 because is at 10Hz
sqrt(avgyro6z$av[11])/(sqrt(avgyro6z$time[11]/3))
sqrt(avgyro4z$av[18])/(sqrt(avgyro4z$time[18]/3))

rrwalk = 4.704003e-06 

#Simulation of the measurement model of the Kalman filtering
#and simple attitude computation.
drift[1] = 0.00
for (k in 2: ((10*3400)/98))
{
	noiseterm <- rnorm (1, mean=0, sd=(whitenoise))
	driftterm <- rnorm (1, mean=0, sd=rrwalk)
	drift[k] = drift[k-1] + delta * driftterm
	modelgyrosz[k] = (filtervalues$gyrosz[k]) -drift[k] - (noiseterm)
}

windows()
plot (filtervalues$time[1:((10*3400)/98)], filtervalues$gyrosz[1:((10*3400)/98)], type="l", xlab="", ylab="")
lines (filtervalues$time[1:((10*3400)/98)], modelgyrosz, col="blue")
lines (rep(mean(filtervalues$gyrosz[1:((10*3400)/98)]), length(filtervalues$gyrosz[1:(10*3400)])), xlab="", ylab="", col="red")
lines (rep(mean(filtervalues$gyrosz[1:((10*3400)/98)])-sd(filtervalues$gyrosz[1:((10*3400)/98)]), length(filtervalues$gyrosz[1:((10*3400)/98)])), xlab="", ylab="", col="blue")
lines (rep(mean(filtervalues$gyrosz[1:((10*3400)/98)])+sd(filtervalues$gyrosz[1:((10*3400)/98)]), length(filtervalues$gyrosz[1:((10*3400)/98)])), xlab="", ylab="", col="blue")

sd(filtervalues$gyrosz[1:((10*3400)/98)])
sd (modelgyrosz)
mean(filtervalues$gyrosz[1:((10*3400)/98)])
mean (modelgyrosz)
plot (drift, type="l")


#Gyros of Xaxis
windows()
plot (filtervalues$time, filtervalues$gyrosx, col="black", type="l", xlab="", ylab="")
lines (rep(mean(filtervalues$gyrosx), length(filtervalues$gyrosx)), xlab="", ylab="", col="red")
lines (rep(mean(filtervalues$gyrosx)-sd(filtervalues$gyrosx), length(filtervalues$gyrosx)), xlab="", ylab="", col="blue")
lines (rep(mean(filtervalues$gyrosx)+sd(filtervalues$gyrosx), length(filtervalues$gyrosx)), xlab="", ylab="", col="blue")

title(main = "Gyroscopes values", xlab = "Times (Sec)", ylab = "Angular velocity(rad/s^2)")

#Gyros of Z axis
windows()
plot (filtervalues$time, filtervalues$gyrosz, col="black", type="l", xlab="", ylab="")
lines (rep(mean(filtervalues$gyrosz), length(filtervalues$gyrosz)), xlab="", ylab="", col="red")
lines (rep(mean(filtervalues$gyrosz)-sd(filtervalues$gyrosz), length(filtervalues$gyrosz)), xlab="", ylab="", col="blue")
lines (rep(mean(filtervalues$gyrosz)+sd(filtervalues$gyrosz), length(filtervalues$gyrosz)), xlab="", ylab="", col="blue")

title(main = "Gyroscopes values", xlab = "Times (Sec)", ylab = "Angular velocity(rad/s^2)")



delta = 1

#Incl X
imt30Gyros <- ts (data.frame(filtervalues$gyrosx), start=c(filtervalues$time[1]), delta=mean(delta))
avgyrofiltx <- avar (imt30Gyros@.Data, frequency (imt30Gyros ))

windows()
plot (avgyro1x$time[1:length(avgyro1x$time)-1],sqrt(avgyro1x$av[1:length(avgyro1x$av)-1]),log= "xy", xaxt="n" , yaxt="n", col="blue", xlab="", ylab="")
lines (avgyro1x$time[1:length(avgyro1x$time)-1],sqrt(avgyro1x$av[1:length(avgyro1x$av)-1]), col="blue")
lines (avgyrofiltx$time[1:length(avgyrofiltx$time)-1],sqrt(avgyrofiltx$av[1:length(avgyrofiltx$av)-1]), col="black")

plotAV(avgyrofiltx)




###########################
# Autocorelation Function #
###########################

library(tseries)
library(bspec)


values <- read.table ("977Hz/moreinstables/24062010-194634-x.data", sep="\t")

names(values) = c("time", "gyrox")

delta = 1/976

imt30Gyros <- ts (data.frame(values$gyrox), start=c(values$time[1]), delta=mean(delta))
windows()
acf (imt30Gyros,two.sided=TRUE, main = "IMT30 Gyro X Autocorrelation function")

acf(imt30Gyros, type = "covariance")
pacf (imt30Gyros)


require(graphics)

## Examples from Venables & Ripley
acf(lh)
acf(lh, type = "covariance")
pacf(lh)




