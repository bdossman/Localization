
##################################################################################################
#
# 				Localization Project Version 2 Update: 9/18/13
#					Bryant Dossman, Phil Taylor, John B.
#
##################################################################################################
##Libraries

library(rgl); library(circular);library(geometry); library(lattice); library(latticeExtra)


##################################################################################################
#
# 				Antenna layout file : Antenna Positions
#
##################################################################################################

##############         Creating Circular  Antenna Layout				###################

r = 200		# in meters
angular.steps =as.numeric(c( pi/4, pi/2, 3/4*pi, pi, 5/4*pi, 3/2*pi, 7/4*pi,0))

Easting <- r*sin(angular.steps)
Northing <- r*cos(angular.steps)
Elevation <- rep(10, length(angular.steps))	# Height of Towers

ID <- letters[1:length(angular.steps)]
Type <- rep("5", length(angular.steps))
Location <- paste("Tower",1:length(angular.steps), sep=" ")
Orientation <- rep(0, length(angular.steps))

circle_ant_layout <- as.data.frame(cbind(ID,Type,Location,Easting,Northing,Elevation,Orientation))

write.table(circle_ant_layout,"circle_ant_layout.csv",sep=",")

##############         Creating Football Field Antenna Layout			###################

n.towers = 9
Easting <- c(rep(-50,3),rep(0,3),rep(50,3))
Northing <- rep(c(-25,0,25),3)
Elevation <- rep(10, length(n.towers))	# Height of Towers
ID <- letters[1:n.towers]
Type <- rep("5", n.towers)
Location <- paste("Tower",1:n.towers, sep=" ")
Orientation <- rep(0, n.towers)

football_ant_layout <- cbind(ID,Type,Location,Easting,Northing,Elevation,Orientation)

write.table(football_ant_layout,"football_ant_layout.csv",sep=",")

# Loading antenna layout files

ant.layout.filename <- "circle_ant_layout.csv"


ant.pattern.filenames <- list (
                               "5" = "yagi_5_pattern.txt",
                               "9" = "yagi_9_pattern.txt"
                               )

skip.ant <- c("OC-h", "Manual")

#################################################################################
#
#  Functions
#
################################################################################
#	Using this simulated data create a dataset of radio hits on a collection of
# 	towers. Implementing portions of John's code (## svn: $Id: mapper.R 11
#	2010-01-14 04:46:02Z john $)
################################################################################

## Read the antenna layout file

read.antenna <- function(filename) {
  f <- file(filename, "r")
  apd <- scan(f, n=2, quiet=TRUE)
  ## The gain pattern matrix has theta.Orfan changing by row, phi.Orfan changing by column.
  pat <- matrix(scan(f, quiet=TRUE), nrow=apd[1], ncol=apd[2], byrow=TRUE)
  close(f)

  ## angular resolution:
  th.step <- 360 / (apd[1] - 1)
  phi.step <- 90 / (apd[2] - 1)

  ## remap to our coordinates; create a matrix with phi changing by column, theta by row
  g <- NA * pat

  ## Our equivalents of Orfanidis' coordinates
  th <- c(rad(phi.step * (col(pat)-1)))
  phi<- c(rad(90 - th.step * (row(pat)-1)))

  ## This laborious conversion works but there is a cleaner way that I wasn't
  ## able to figure out.

  x <- cos(phi)*cos(th)
  y <- cos(phi)*sin(th)
  z <- sin(phi)

  phi2 <- deg(abs(asin(z)))
  th2 <- deg(atan2(y, x)) %% 360

  g[cbind(1 + round(th2 / th.step), 1 + round(phi2 / phi.step))] <- c(pat)

  ## fill in the other half by symmetry
  th2 <- 360 - th2
  g[cbind(1 + round(th2 / th.step), 1 + round(phi2 / phi.step))] <- c(pat)

  return(structure(g, phi.step=phi.step, th.step=th.step, phi.step.rad=rad(phi.step), th.step.rad=rad(th.step)))
}

read.layout <- function(filename) {
  x <- read.csv(filename, header=TRUE, as.is=TRUE)
  rownames(x) <- x$ID
  x$Type <- as.character(x$Type)
  x
}

## lookup antenna gain by spherical angles
## with phi and theta in radians

gain <- function(pat, phi, theta) {
  pat[round((c(theta) %% (2*pi)) / attr(pat, "th.step.rad")) + 1, round(pmin(abs(c(phi)), pi/2) / attr(pat, "phi.step.rad")) + 1]
}

## predicted relative power from target given range and spherical
## angles from antenna axis.  This is returned in dbm.

## pat is an antenna pattern matrix
## r is in metres; phi, theta are in radians
## P is transmitter power

## We square both 1/(r+1) and voltage pattern to get power.  The
## "+1" is to avoid the singularity at 0.

pow.sph <- function(pat, r, phi, theta, P=1) {
  ## Power in dbm
  10 * log10(P * (gain(pat, phi, theta)/(r+1))^2)

  ## Power in absolute units (use this instead if the Lotek 
  ## power measurement is not in log units)
  
  ## P * (gain(pat, phi, theta)/(r+1))^2
}

## predict power from target given its cartesian coordinates
## and a row from the antenna matrix

pow <- function(ant, pats, x, y, z, P=1) {
  ax <- ant$Easting
  ay <- ant$Northing
  az <- ant$Elevation
  r <- sqrt((x-ax)^2 + (y-ay)^2 + (z-az)^2)
  phi <- asin((z-az) / max(1, r))  ## use at least r = 1
  theta <- atan2(y - ay, x - ax) - rad(90 - ant$Orientation)
  
  pow.sph(pats[[ant$Type]], r, phi, theta, P)
}

## Read the antenna layout file

read.layout <- function(filename) {
  x <- read.csv(filename, header=TRUE, as.is=TRUE)
  rownames(x) <- x$ID
  x$Type <- as.character(x$Type)
  x
}

## transform points in 3d:

## 1: rotate around the origin in the xz-plane by theta degrees
## 2: scale by "scale" w.r.t. the origin
## 3: shift by offs along x, y, z axes

## l: list with components x, y, and z
## returns new list with components x, y, and z

tx3d <- function(l, offs=c(0, 0, 0), scale = 1, theta = 0) {
  theta <- rad(90 - theta)
  ct <- cos(theta)
  st <- sin(theta)
  list (x = offs[1] + scale * (l$x * ct - l$z * st),
        y = offs[3] + scale * l$y,
        z = offs[2] + scale * (l$z * ct + l$x * st))
}

map <- function(ants, pats, dat) {
  ## antenna coordinates
  
  ax <- ants[dat$antenna, "Easting"]
  ay <- ants[dat$antenna, "Northing"]
  az <- ants[dat$antenna, "Elevation"]
  aazi <- rad(90 - ants[dat$antenna, "Orientation"])  ## convert from compass to math angle
  atype <- ants[dat$antenna, "Type"]

  dat$dbm <- dat$power
  
  ## given estimated target position and transmitter power, get
  ## deviance between observed and predicted power
  
  dev <- function(x) {
    ## x[1:3]: target x, y, and z position

    ## compute offset from each antenna to target in spherical coordinates
    ## KLUDGE: force r to have a minimum value so that phi is defined

    ## do assignments in the parent frame so that if debugging is enabled
    ## by uncommenting the browser() call below, the variables set by
    ## dev can be examined.
    
    r <<- pmax(1, sqrt((x[1] - ax)^2 + (x[2] - ay)^2 + (x[3] - az)^2))

    theta <<- atan2(x[2] - ay, x[1] - ax)

    phi <<- asin((x[3] - az) / r)

    ## For each antenna, predict power and get deviance
    pp <<- numeric(length(ax))
    for (i in seq(along=ax))
      ## predict power
      pp[i] <<- pow.sph(pats[[atype[i]]], r[i], phi[i], theta[i] - aazi[i])
    
    rv <- sum((dat$dbm - pp)^2)

    if (is.finite(rv))
      return(rv)

    return(1e6)  ## bogus large value
  }

  ## initial target position estimate is centroid of antenna positions
  cent <<- c(mean(ax), mean(ay), mean(az))
  
  mod <- optim(cent, dev)
##  browser()
  return(mod)
}


map.set <- function(ants, pats, dat, win.time=60) {

  ## break up time into slots
  t <- as.numeric(dat$ts)
  tlen <- diff(range(t))
  nslot <- ceiling(tlen / win.time)
  t.cuts <- min(t) + tlen * (0:nslot) / nslot
  ## weird: why is this necessary for getting the top couple of values included?
  t.cuts[nslot+1] <- t.cuts[nslot+1] * 1.1
  
  win <- unclass(cut(t, t.cuts, include.lowest=TRUE))

  x <- y <- z <- code <- fit <- iter <- rep(NA, nslot)
  
  for (i in 1:nslot) {
    ## catch errors so that we fit the model whenever possible
    tryCatch({
      mod <- map(ants, pats, dat[win==i,])
      x[i] <- mod$par[1]
      y[i] <- mod$par[2]
      z[i] <- mod$par[3]
      code[i] <- mod$convergence
      fit[i] <- mod$value
      iter[i] <- mod$counts[1]
    }, error=function(e) print(paste("map.set error: ", as.character(e))))
  }
  return(data.frame(Easting=x, Northing=y, Elevation=z, t=min(t) + tlen * (1:nslot - 0.5) / nslot, code=code, fit=fit, iter=iter))
}

plot.pat <- function(ant, loc=c(0, 0, 0), scale=1, orient = 90, legend=TRUE) {
 
  ## convert gain to log scale and set floor at -20 db
  g <- 10*log10(.Machine$double.eps+ant)
  g <- g - max(g)
  g <- 20 + pmax(g, -20)

  ## get spherical coordinates
  th <- rad(attr(ant, "th.step") * (col(g)-1))
  phi<- rad(attr(ant, "phi.step") * (row(g)-1))
  
  ## colour palette
  hc <- heat.colors(201)
  
  ## plot points
  rgl.points(tx3d(list(x=g*cos(phi)*cos(th), y=g*sin(phi), z=g*cos(phi)*sin(th)), loc, scale, orient), col=hc[1+c(round(10*g))])

  ## reflect in the xz plane by negating y (table only holds values for phi >= 0)
  rgl.points(tx3d(list(x=g*cos(phi)*cos(th), y=g*sin(phi), z=-g*cos(phi)*sin(th)), loc, scale, orient), col=hc[1+c(round(10*g))])

  if (legend) {
    
    ## plot a gain legend
    
    rgl.texts(tx3d(list(x=20, y=2.5*-1:10, z=20), loc, scale, orient), text=c(paste(-20 + 2*(0:10), "db Max"), "Gain"), col=c(hc[1+20*(0:10)], "#FFFFFF"))
    
    ## show positive axes
    
    rgl.lines(tx3d(list(x=c(0,30), y=c(0,0), z=c(0,0)), loc, scale, orient),col="green")
    rgl.texts(tx3d(list(x=30, y=0, z=0), loc, scale, orient),text="X")
    rgl.lines(tx3d(list(x=c(0,0), y=c(0,30), z=c(0,0)), loc, scale, orient),col="green")
    rgl.texts(tx3d(list(x=0, y=30, z=0), loc, scale, orient),text="Z")
    rgl.lines(tx3d(list(x=c(0,0), y=c(0,0), z=c(0,30)), loc, scale, orient),col="green")
    rgl.texts(tx3d(list(x=0, y=0, z=30), loc, scale, orient),text="Y")
  }
}

do.test <- function(antfile, pats, np = 500, win=10, jitter=FALSE, stationary=TRUE, init=c(0,0,0),model=FALSE, move.model=NULL, surpress=FALSE) {

  ants <- read.layout(antfile)
  
  xdim <- diff(range(ants$Easting))
  ydim <- diff(range(ants$Northing))
  zdim <- diff(range(ants$Elevation))

  scale <- max(xdim, ydim, zdim)
  
  xc <- init[1]
  yc <- init[2]
  zc <- init[3]
  
  if (stationary==TRUE) {		## simulates a stationary tag point amongst towers centered around init=c(x,y,z)
  
  	x <- rep(init[1], np)
  	y <- rep(init[2], np)
  	z <- rep(init[3], np)
  
  }
  if(stationary==FALSE){

  	v <- 2    ## mean speed in metres/second
  	vth <- 0   ## mean angular speed in degrees/second
  	hmax <- 100 ## maximum height, in m
  	hstep <- 0 ## mean change in height, in metres / second
  
  	## generate motion steps, using an exponential distribution of
  	## step lengths and turning angles

  	step.r <- rexp(np, 1 / v)

  	## generate turning angles, with random signs
  	turn.th <- c(pi/4,rep(0,np-1))

  dh <- rexp(np, 1/hstep)
  
  ## positions:

  theta <- rad(cumsum(turn.th))
  dz <- pmin(hmax, pmax(0, cumsum(dh)))  ## constrain within [0, hmax]

  x <- xc + cumsum(step.r * cos(theta))
  y <- yc + cumsum(step.r * sin(theta))
  z <- zc + dz

  }
  
  ## if you give the function constant.move...it will replace the location
  ## simulation with what you give it
  
if (model==TRUE){
  x <- move.model[,1]
  y <- move.model[,2]  
  z <- move.model[,3]
} 
  ## plot true path and antennas
	
  if(surpress==FALSE){
  rgl.open()
  points3d(x, z, y, col="green")
  lines3d(x, z, y, col="green")
}
  nant <- dim(ants)[1]
  
 if(surpress==FALSE){
  for(i in 1:dim(ants)[1]){
    plot.pat(pats[[ants$Type[i]]], as.numeric(ants[i,4:6]), 2, as.numeric(ants[i, 7]), i == 1)
}
}
  ## for each time step, pick an antenna at random and determine
  ## the power it would receive from the target at its current location

  ant <- rep(1:nant, length=np)
  pow <- numeric(np)

  for (i in 1:np){
    pow[i] <- pow(ants[ant[i],], pats, x[i], y[i], z[i])
  }

  if (jitter)
    pow <- jitter(pow)
  
  ## create a fake data frame from the path and antenna power:

  fake <- data.frame(ts = 1:np,
                     id = "fake",
                     antenna=as.character(ants$ID[ant]),
                     power = pow,                          ## FIXME: add some noise with jitter()
                     Easting = x,
                     Northing = y,
                     Elevation = z,
                     stringsAsFactors = FALSE
                     )


  ## try estimate locations from the fake data, in 10 second groups
  m <- map.set(ants, pats, fake, win)


  ## modifying code to surpress the plotting function
  
  if(surpress==FALSE){
 
  ## map the reconstructed path
  points3d(m$Easting, m$Elevation, m$Northing, color="red")

  lines3d(m$Easting, m$Elevation, m$Northing, color="red")

  bbox3d()
  }
  ## calculating a position accuracy estimate XY direction
  
  groups 	  = sort(rep(1:nrow(m),win))
  east        = tapply(fake$Easting, groups, mean)
  north       = tapply(fake$Northing, groups, mean)
  elevation	  = tapply(fake$Elevation, groups, mean)
  acc.xy      = sqrt((east-m$Easting)^2+(north-m$Northing)^2)
  acc.xyz	  = sqrt((east-m$Easting)^2+(north-m$Northing)^2+(elevation-m$Elevation)^2)
  
  acc 		  = data.frame (acc.xy=acc.xy, acc.xyz=acc.xyz)
  
  return(list(fake=fake, fit=m, acc=acc))
}

#################################################################################################
#
#			Simulating movements from a given movement model and then translating those xyz 
#			coordinates into a dataset of signal strengths from multiple antennas usings
#           john's modified code (above)
#
################################################################################################

## Reading in Pattern Layout Files
setwd("/Users/Dossman/Desktop/Current Projects/Localization Project/August Work/lotek_mapper")
pats <- lapply(ant.pattern.filenames, read.antenna)

## Reading Circle Antenna Layout Files
setwd("~/Desktop/Current Projects/Localization Project")

ants <- read.layout("circle_ant_layout.csv")

## Making omnidirectional antennas/ isotropic

iso.pats <- pats
iso.pats[["5"]][] <- 1  ## uniform gain of 1 in all angular directions

## Simulating dataset... set stationary = FALSE for straight line movement

do.test("circle_ant_layout.csv", iso.pats,np=500, win=20, stationary=FALSE,init=c(-200,0,0))

#############################################################################################
#
#			Visualizing Parameter space to asses eficacy of the optimize function
#
#############################################################################################
#############################################################################################
#
#			Testing # putting together code to plot parameter space.
#
##############################################################################################
## TODO: Wrap below code into a function. It probably can be made faster by simply vecotorize
## alot of these loops. Also a spatial resolution of a meter is probably not neccsary for this
## figure.


n=8 # number of antennas
x=seq(-500,500)
y=seq(-500,500)
ant <- rep(1:n, length=n)


## Picking a starting location. You need to plot deviance from a known location. So must start
## with a given location first. e.g 0,0,0

loc <- NULL
for (i in 1:n){
	loc[i] <- jitter(pow(ants[ant[i],], iso.pats, 0,0,0), amount=15)
}

## levelplot works well with lists.

dev.raster <- list("x"=x,"y"=y, z=rep(NA, length(x)*length(y))) # we will use this list to plot deviance

area<-expand.grid(x=x,y=y)

## building a n (number of antennas) by X by Y dimensional arrays. Predicited powers from each
## antenna for any combination of x and y coordinates.

k <- array(data=NA,dim=c(length(x),length(y),n), dimnames=c("x", "y", "ant"))

for (i in 1:n){
	k[,,i] <- pow(ants[ant[i],], iso.pats, area$x,area$y,0)
}

## plotting received power from a given tower

# test <- list(x=x,y=y, z=k[,,1])
# image(test)
# for (i in 1:n){
# points(ants[ant[i],4],ants[ant[i],5], pch=19, col="green")
# }

## Calculating deviance of all the predicted powers from a given known location. 

tmp=NULL
dev.raster$z <- matrix(rep(0, length(x)*length(y)), nrow=length(x), ncol=length(y))

for (i in 1:length(y)){
	for (j in 1:length(x)){
		dev.raster$z[j,i] <- sum((loc-k[j,i,])^2)
	}
}

require(lattice)

#image(dev.raster, col = topo.colors(20), xlab="X", ylab="Y")

levelplot(dev.raster$z,row.values=dev.raster$x, column.values=dev.raster$y, xlim=range(x), ylim=range(y), xlab="X", ylab="Y", main="Parameter Space")

#####################################################################################
#
#			Testing Localization VS Triangulation
#
#####################################################################################
## Simple movement model we will use for each experiment

steps = 500
start = -100
end   =  100

x <- seq(start, end, length.out=steps)
y <- seq(start, end, length.out=steps)
z <- rep(0, steps)

move.model <- cbind(x,y,z)

## Modified do.test code to predict locations with this movement model and provide the accuracy in
## location estimates.

do.test("circle_ant_layout.csv", iso.pats,np=500, win=20, stationary=FALSE, model=TRUE, move.model=move.model)

## Creating function to create a random antenna layout

ant.layout <- function(nants,elev.fixed=TRUE, elev=10, max.elev=10, D=100, circular=TRUE, r=0){
	nants=nants
	
	Easting  = rnorm(nants,mean=0, sd=D)
	Northing = rnorm(nants,mean=0, sd=D)
		
	## creating the option to randomly assign an elevation to each ant as well
	if (elev.fixed==FALSE) Elevation <- runif(nants,0,max.elev)
						   Elevation <- rep(elev, nants)
	
	## if circular==TRUE
	if(circular==TRUE){
			sp <- 360/nants
			Easting <- r*sin(rad(seq(1:nants) * sp))
			Northing <- r*cos(rad(seq(1:nants) * sp))
	}
	
	# assuming we use omni-directional ants and orientation doesn't matter. Same goes for IDs and type of 
	# ants used.

	
	ID <- 1:nants
	Type <- rep("5", nants)
	Location <- paste("Tower",1:nants, sep=" ")
	Orientation <- rep(0, nants) 
		
	## Putting it all together in a list
	
	ant_layout <- data.frame(ID,
							Type,
							Location,
							Easting,
							Northing,
							Elevation,
							Orientation
							)

   return(ant_layout)
}

### Testing the random array generator. I am interested in looking at how position accuracy varies
### according to the number of ants present and the dispersion of the antennas i.e extent of coverage


test<-ant.layout(5,D=150, circular=FALSE)

write.table(test,"test.csv", sep=",")

plot(test$Northing~test$Easting, ylim=c(-500,500), xlim=c(-500,500), pch=10, ylab="Northing", xlab="Easting")


do.test("test.csv", iso.pats,np=500, win=20, stationary=FALSE, model=TRUE, move.model=move.model, surpress=TRUE)


## Testing the effect of nants and D on localization accuracy

# # n      <- sort(rep(seq(5,30,by=5), 20))
# D      <- rep(c(100,400, 800),40)

# test   <- data.frame(n = n,
					 # D = D,
					 # stringsAsFactors = TRUE
					 # )
# start <- Sys.time()					 

# for (i in 1:length(n)){
	# layout <- ant.layout(n[i],D=D[i], circular=FALSE)
	# write.table(layout,"test.csv", sep=",")
	
	# trial = do.test("test.csv", 
					# iso.pats,np=500, 
					# win=20, 
					# stationary=FALSE, 
					# model=TRUE, 
					# move.model=move.model,
					# surpress=TRUE
					# )
	
	# test$mean.xy[i] <- mean(trial$acc$acc.xy)
	# test$var.xy[i]  <- var(trial$acc$acc.xy)
	
	# test$mean.xyz[i] <- mean(trial$acc$acc.xyz)
	# test$var.xyz[i]  <- var(trial$acc$acc.xyz)
	
# }

# end <- Sys.time()

# difftime(end,start)

#####################################################################
#		Preforming a simulation with two different layouts and Analyzing the effect that number of antennas 
#   	and the distances between antennas has on localization accuracy.
#    
#####################################################################

test<-NULL

n      <- sort(rep(seq(1,50,by=2), 12))
D      <- rep(c(100,400,1200),100)

test   <- data.frame(n = n,
					 D = D,
					 stringsAsFactors = TRUE
					 )
start <- Sys.time()						 
for (i in 1:length(n)){

	layout <- ant.layout(n[i],r=D[i], circular=TRUE) # modified for circular layout
	write.table(layout,"test.csv", sep=",")
	
	trial = do.test("test.csv", 
					iso.pats,np=500, 
					win=25, 
					stationary=FALSE, 
					model=TRUE, 
					move.model=move.model,
					surpress=TRUE
					)
	
	test$mean.xy[i] <- mean(trial$acc$acc.xy)
	test$median.xy[i] <- median(trial$acc$acc.xy)
	test$var.xy[i]  <- var(trial$acc$acc.xy)
	
	test$mean.xyz[i] <- mean(trial$acc$acc.xyz)
	test$median.xyz[i] <- median(trial$acc$acc.xyz)
	test$var.xyz[i]  <- var(trial$acc$acc.xyz)
	
}

end <- Sys.time()
print(difftime(end,start))

circle.layout.data <- test

write.table(circle.layout.data, "circle_layout_sim5.csv", sep=",")

#### For random layout now


start <- Sys.time()						 
for (i in 1:length(n)){

	layout <- ant.layout(n[i],D=D[i], circular=FALSE) # modified for circular layout
	write.table(layout,"test.csv", sep=",")
	
	trial = do.test("test.csv", 
					iso.pats,np=500, 
					win=25, 
					stationary=FALSE, 
					model=TRUE, 
					move.model=move.model,
					surpress=TRUE
					)
	
	test$mean.xy[i] <- mean(trial$acc$acc.xy)
	test$median.xy[i] <- median(trial$acc$acc.xy)
	test$var.xy[i]  <- var(trial$acc$acc.xy)
	
	test$mean.xyz[i] <- mean(trial$acc$acc.xyz)
	test$median.xyz[i] <- median(trial$acc$acc.xyz)
	test$var.xyz[i]  <- var(trial$acc$acc.xyz)	
}

end <- Sys.time()
print(difftime(end,start))

random.layout.data <- test

write.table(random.layout.data, "random_layout_sim5.csv", sep=",")

#####################################################################
#  Analysis of Data Simulation and mean XY localization accuracy
#####################################################################
# Restart Analysis here by reading in datasets created above
random.layout.data <- read.csv("random_layout_sim3.csv", header=TRUE)
circle.layout.data <- read.csv("circle_layout_sim3.csv", header=TRUE)

## merging both datasets

random.layout.data$method<-"random"
circle.layout.data$method<-"circle"

sim.data <- as.data.frame(rbind(random.layout.data,circle.layout.data))

circular.model <- lm(mean.xy ~ log(D) * n, data = subset(sim.data,method=="circle"))
summary(circular.model)
capture.output(summary(circular.model), file="test.doc")

random.model <- lm(mean.xy ~ D * n, data = subset(sim.data,method=="random"))
summary(random.model)
capture.output(summary(random.model), file="test.doc")


## plotting localization accuracy as a function of antenna dispersion

xyplot(mean.xy ~ log(D), data = sim.data,
	   type   	= c("p","r"),
	   scales 	= list(log=F),
	   groups= method,
	   pch    	= 1:2,
	   xlab   	= c("Log (Random: D or Circular: r)"),	
	   ylab   	= c("Accuracy (m)"),
	   auto.key =  list( 
           			 lines = TRUE,
				 	 columns = 2,
				 	 corner=c(0.01,0.98)),
	 main = "Localization Accuracy as a Function of Antenna Dispersion"
	 )
	 
	 

## plotting localization accuracy as a function of the log(number of antennas)

xyplot(median.xy ~ log(n), data = sim.data,
	   type   	= c("p","smooth"),
	   scales 	= list(log=F),
	   groups= method,
	   pch    	= 1:2,
	   xlab   	= c("Log Number of Antennas log(n)"),	
	   ylab   	= c("Accuracy (m)"),
	   auto.key =  list( 
           			 lines = TRUE,
				 columns = 2,
				 corner=c(0.01,0.98)),
	 main = "Localization Accuracy as a Function of the Number of Antennas"
	 )

## plotting localization accuracy as a function of the number of antennas
	  
xyplot(median.xy ~ n, data = sim.data,
	   type   	= c("p","smooth"),
	   scales 	= list(log=T),
	   groups= method,
	   pch    	= 1:2,
	   xlab   	= c("Number of Antennas (n)"),	
	   ylab   	= c("Accuracy (m)"),
	   auto.key =  list( 
           			 lines = TRUE,
				 columns = 2,
				 corner=c(0.01,0.98)),
	 main = "Localization Accuracy as a Function of the Number of Antennas"
	 )






