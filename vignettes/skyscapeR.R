## ---- echo = F-----------------------------------------------------------
knitr::opts_chunk$set(fig.width = 3.5, fig.height = 3.5, fig.align='center')

## ------------------------------------------------------------------------
library(skyscapeR)

## ------------------------------------------------------------------------
data(RugglesRSC)
hist <- histogram(RugglesRSC$Dec, 2)

## ---- fig.show='hold'----------------------------------------------------
plot(hist)

## ---- fig.show='hold'----------------------------------------------------
data(RugglesRSC)
unc <- runif(length(RugglesRSC$Dec),1,10)  # random values between 1 and 10
curv <- curvigram(RugglesRSC$Dec, unc)
plot(curv)

## ------------------------------------------------------------------------
lunar <- sky.objects('moon', epoch=-2000, col='red', lty=2)

## ---- fig.show='hold'----------------------------------------------------
plot(curv, obj=lunar)

## ------------------------------------------------------------------------
georef <- rbind( c(35.1, -7.1),     # GPS data
                 c(35.1, -7),
                 c(35.2, -7.1),
                 c(35.1, -7.3) )
azimuths <- c(93, 108, 105, 98)    # Compass measurements
altitudes <- c(2, 1.5, 0.5, 1)    # Clinometer measurements

data <- reduct.compass(loc=georef, mag.az=azimuths, date="2017/06/13", alt=altitudes)
data

## ------------------------------------------------------------------------
data <- reduct.compass(loc=georef, mag.az=azimuths, date="2017/06/13")
data

## ------------------------------------------------------------------------
georef <- c( ten(35,50,37.8),     # GPS data
             ten(14,34,6.4) )
az <- c( ten(298,24,10),     # Theodolite H measurements
         ten(302,20,40) )
alt <- c( ten(1,32,10),     # Theodolite V measurements
          ten(0,2,27) )
az.sun <- ten(327,29,50)    # The azimuth of the sun as measured at time
limb <- "right"   # Which limb of the sun was targeted
date <- "2016/02/20"
time <- "11:07:17"    # Time the sun measurement was taken
timezone <- "Europe/Malta"    # Timezone corresponding to time above

data <- reduct.theodolite(loc=georef, az, date, time, timezone, az.sun, limb, alt)
data

## ---- fig.show='hold'----------------------------------------------------
az <- rnorm(30, 85, 20)    # This creates 30 random azimuths
plotAz(az)

## ---- fig.show='hold'----------------------------------------------------
sunandmoon <- sky.objects(c('sun','moon'), epoch=-4000, col=c('blue','red'), lty=c(2,3))
plotAz(az, obj=sunandmoon, loc=c(52,0))

## ------------------------------------------------------------------------
dec <- az2dec(az, loc= c(35,-7), alt=0)
dec

## ------------------------------------------------------------------------
az <- c(0,90,180,270,360)
alt <- c(0,5,5,0,0)
georef = c(40.1, -8)
hor <- createHor(az, alt, loc=georef, name= 'Horizon Profile 1')

## ---- fig.width = 7, fig.height = 3--------------------------------------
plot(hor)

## ---- fig.width = 7, fig.height = 3--------------------------------------
hor <- downloadHWT('NML6GMSX')
plot(hor)

## ------------------------------------------------------------------------
hor2alt(hor, az=90)
hor2alt(hor, az=110)

## ------------------------------------------------------------------------
data <- reduct.compass(loc=hor, mag.az=az, date="2017/06/13")
data

## ---- fig.width = 7, fig.height = 3--------------------------------------
plot(hor, obj=lunar)

## ---- fig.width = 7, fig.height = 3--------------------------------------
aux <- sky.objects(names=c('sun', 'Aldebaran'), epoch=c(-4300,-3700), col=c('blue', 'red'))
plot(hor, obj=aux)

## ------------------------------------------------------------------------
ss <- star('Sirius')
ss

## ------------------------------------------------------------------------
ss <- star('Sirius', year=-4000)
ss$dec

## ------------------------------------------------------------------------
# ncores forced to 2 for production of this document
sp <- star.phases('Sirius', -3000, loc=c(30.0,31.2), alt.hor=2, ncores=2)  

## ------------------------------------------------------------------------
sp$phase
sp$events
sp$seasons

## ---- fig.show='hold', fig.width = 5, fig.height = 1---------------------
plot(sp)

## ------------------------------------------------------------------------
data <- c()
data$az <- RugglesRSC$CL_Rec_C
data$az.unc <- rep(5, length(data$az))
data$lat <- RugglesRSC$Latitude
data$alt <- RugglesRSC$CL_Mean_Alt
# ncores forced to 2 and nsims to 100 for production of this document
sg <- sigTest(data, ncores=2, nsims=100)   

## ---- fig.show='hold'----------------------------------------------------
plot(sg)

## ---- fig.show='hold'----------------------------------------------------
plot(sg, show.local=T)

## ------------------------------------------------------------------------
print(sg)

