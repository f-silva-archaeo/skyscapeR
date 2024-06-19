## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(skyscapeR)

## -----------------------------------------------------------------------------
BC.AD(1)
BC.AD(0)
BC.AD(-1)
BC.AD(-99)
BC.AD(-100)

## -----------------------------------------------------------------------------
skyscapeR.vars() # shows current values
skyscapeR.vars(timezone='CET') # to set timezone to Central European Time
skyscapeR.vars(calendar='J') # to set calendar to Julian
skyscapeR.vars(refraction=F) # to switch refraction calculations off

## -----------------------------------------------------------------------------
skyscapeR.vars(timezone='GMT') # to set timezone to Central European Time
skyscapeR.vars(calendar='G') # to set calendar to Julian
skyscapeR.vars(refraction=T) # to switch refraction calculations off


## -----------------------------------------------------------------------------
az2dec(az=92, loc=35.125, alt=2)

## -----------------------------------------------------------------------------
# geocentric declination for june solstice at 2501 BC
jS(-2500)

# topocentric declination for december solstice at 2501 BC from Stonehenge
dS(-2500, loc=c(51.17889,1.826111,10))

# geocentric declination for northern minor lunar extreme at 3001 BC
nmnLX(-3000)

# topocentric declination for northern minor lunar extreme at 3001 BC from Stonehenge
sMjLX(-3000, loc=c(51.17889,1.826111,10))

# zenith sun for Teotihuacan
zenith(loc=c(19.6925,98.84389,200))

## -----------------------------------------------------------------------------
mag.dec(loc=c(35,-8), date="2020/06/07")

## -----------------------------------------------------------------------------
loc <- c(35,-8,100)
mag.az <- c(89.5, 105, 109.5)
data <- reduct.compass(loc, mag.az, "2016/04/02", alt=c(1,2,0))
data

## -----------------------------------------------------------------------------
loc <- c(52,-3,100)
sunAz(loc, '2017-10-04 12:32:14', 'Europe/London')

## -----------------------------------------------------------------------------
lat <- ten(35,50,37.8) # latitude
lon <- ten(14,34,6.4) # longitude
elev <- 100
az <- c(ten(298,24,10), ten(302,20,40))
alt <- c(2, 5)
az.sun <- ten(327,29,50)
date <- "2016/02/20"
time <- "11:07:17"
data <- reduct.theodolite(c(lat,lon,elev), az, date , time, tz= "Europe/Malta", az.sun, alt=alt)
data

## ---- fig.width = 6, fig.height=4---------------------------------------------
az <- c(0,90,180,270,360)
alt <- c(0,5,5,0,0)
loc <- c(51.17889,1.826111,10) 
hor <- createHor(az, alt, alt.unc=0.5, loc=loc, name='Stonehenge Test')
plot(hor)

## ---- fig.width = 6, fig.height=4---------------------------------------------
hor <- createHWT(lat=ten(51,30,45), lon=ten(0,5,26.1), name='London Mithraeum')
plot(hor)

## ---- fig.width = 5, fig.height=5---------------------------------------------
az <- c(120, 100, 93, 97, 88, 115, 112, 67)
plotAzimuth(az)

## ---- fig.width = 5, fig.height=5---------------------------------------------
tt <- sky.objects('solar extremes', epoch=-2000, loc=c(35,-8), col='red')
plotAzimuth(az, obj=tt)

## ---- fig.width = 5, fig.height=5---------------------------------------------
tt <- sky.objects(c('solar extremes', 'lunar extremes'), epoch=-2000, loc=c(35,-8), col=c('red','blue'))
plotAzimuth(obj=tt)

## ---- fig.width = 5, fig.height=5---------------------------------------------
plotAzimuth(az, obj=tt, lwd=2, col='black')

## ---- fig.width = 6, fig.height=4---------------------------------------------
az <- c(120, 100, 93, 97, 88, 115, 112, 67)
hor <- createHWT(51.17889, 1.826111, name='Stonehenge')
dec <- az2dec(az, loc=hor)
kde <- density(dec)
plot(kde)

## ---- fig.width = 6, fig.height=4---------------------------------------------
kde <- density(dec, bw=2) # forces uncertainty to be a given value, in this case 2 degrees
plot(kde)

kde <- density(dec, bw=2, kernel='epanechnikov') # changes kernel to epanechnikov
plot(kde)

## -----------------------------------------------------------------------------
bernoulli.trial(n=8, p=12/180, r=4)

## ---- fig.width = 6, fig.height=4---------------------------------------------
az <- c(120, 100, 93, 97, 88, 115, 112, 67) # same azimuths as before
unc <- c(2, 3, 2, 2.5, 1.5, 3, 4, 3.5)
az.prob <- az.pdf(az=az, unc=unc)
plot(az.prob)

## ---- fig.width = 6-----------------------------------------------------------
dec.prob <- coordtrans(az.prob, hor) # same horizon as before
plot(dec.prob)

## ---- fig.width = 6, fig.height=4---------------------------------------------
ss <- spd(dec.prob)
plot(ss)

## ---- fig.width = 6, fig.height=4---------------------------------------------
# on azimuth data
st.az <- randomTest(az.prob, nsims=10, ncores=1) 
plot(st.az)

# on declination data
st.dec <- randomTest(dec.prob, nsims=10, ncores=1)
plot(st.dec)

## ---- fig.width = 6, fig.height=4---------------------------------------------
plot(st.dec, show.local=T)
st.dec$metadata$local.pval

## -----------------------------------------------------------------------------
findTargets(c(-25,-18), c(-2499,-1749))

## -----------------------------------------------------------------------------
findTargets(c(-25,-18), c(-2499,-1749), max.mag=3)

## ---- results='hide'----------------------------------------------------------
sp <- star.phases('Sirius', -2999, loc=c(30.0, 31.2, 25), alt.hor=2)

## -----------------------------------------------------------------------------
sp$metadata$type
sp$metadata$events
sp$metadata$seasons

## ---- fig.width = 6-----------------------------------------------------------
plot(sp)

## ---- fig.width = 6-----------------------------------------------------------
sky.sketch(time='2019/01/07 09:30', loc=c(35,-8,100))

