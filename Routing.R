# Optimal Routing Code (shortest distance to each mark)
  
# Load packages
suppressPackageStartupMessages(require(data.table)) # data.table and fread
suppressPackageStartupMessages(require(rgdal)) # GRIB file parser
suppressPackageStartupMessages(require(raster)) # GRIB file parser
suppressPackageStartupMessages(require(geosphere)) # Geo distance calculations
suppressPackageStartupMessages(require(aspace)) # Trig angle calculations 
suppressPackageStartupMessages(require(ggmap)) # Maps
suppressPackageStartupMessages(require(OpenStreetMap)) # Maps
suppressPackageStartupMessages(require(R.utils)) # Arduino serial monitor read
suppressPackageStartupMessages(require(plotrix)) # Polar plots

# Load Polars 
# (using J111 polars at wind speed of 8 - will incorporate full polars in next build!)
# polars.raw <- fread("J111_Polars_WS8.csv", header = TRUE,data.table = TRUE)
polars.raw <- data.table("Angle" = c(0,47,55,60,70,80,90,100,110,120,130,137,180),
                         "BS" = c(0,6,6.95,7.2,7.35,7.60,7.80,7.92,7.94,7.75,7.43,6.75,4.47))
polars <- spline(polars.raw$Angle, polars.raw$BS, n = 181) # spline interpolarion of polars
polars <- data.table("Angle" = polars$x, "BS" = polars$y) 
polars <- append(polars$BS[1:180],polars$BS[181:2])
polars <- data.table("Angle" = 0:359, "BS" = polars/8) # dividing polars by 8 then multiplying by wind speed (will correct later)

# # Load GRIB data
url <- "http://www.glerl.noaa.gov/ftp/EMF/glcfs/grib/NCAST/glcfs-win.grb" # NCAST
# url <- "http://www.glerl.noaa.gov/ftp/EMF/glcfs/grib/FCAST/glcfs-win.grb" # FCAST
GRIBfile <- "glcfs-win.grb"
download.file(url, GRIBfile, mode="wb")

# Create Map for plot
myLocation <- c(-88, 41, -86, 43)
myMap <- get_map(location=myLocation, source="google", maptype="satellite", crop=FALSE)
map <- ggmap(myMap)

# Location of starting point and waypoints (Longitude, Latitude)
CHI <- c(-87.562880, 41.881406) # Race start
MAN <- c(-86.222431, 44.936899) # Manituos
# MAN <- c(-84.616287, 45.842111) # Mackinac
# MAN <- c(-86.508403, 42.097567) # St Joe

wind.grib <- readGDAL(GRIBfile)
wind <- spTransform(wind.grib,CRS("+init=epsg:4326"))
wind.data <- data.table(wind@coords,wind@data)[,.(x,y,band2,band8)]

## Start Routing

# location must be null each time you run a route loop
LOC <- NULL 

# Set current location to Chicago 
if ( is.null(LOC) ) { LOC = CHI }

# Set observation interval (in hours)
obs.interval = .1

# Set counter
counter = 1

destination = NULL

# Begin calculating optimal route
while ( distGeo(LOC,MAN)>10000 ) {

# Set prior location for plot
LOC.prior = LOC

# Find the location of the nearest weather observation point
wind.nearest <- wind.data[which.min(distGeo(wind.data[,.(x,y)],LOC))]

# Define wind speed and direction for polars
wind.nearest.dir = 155 #as.numeric(round(wind.nearest$band8,0)) 
wind.nearest.spd = 9 # as.numeric(wind.nearest$band2)

# Target waypoint direction
bearing = round(bearing(LOC,MAN),0)
heading = bearing + -180:179
heading = ifelse(heading<0,heading+360,heading)

# True wind angles (dir of wind - dir of boat)
TWA = (wind.nearest.dir - heading)
TWA = ifelse(TWA<0,TWA+360,TWA)

# Calculate VMC
VMC <- polars[Angle %in% TWA, BS]*cos_d(heading)

polars$Angle <- as.numeric(polars$Angle)
d = polars[Angle %in% TWA, BS]*wind.nearest.spd*obs.interval

cos.a = cos_d(heading) # = cos(pi/6 radians) #= Sqrt(3)/2 = 0.866025
sin.a = sin_d(heading) # = sin(pi/6 radians) #= 1/2 = 0.5
cos.lat = cos_d(LOC[2]) # = cos(-0.00548016 radian) = 0.999984984
r = d/0.000539957 # Nautical mile to meter conversion (knots/(meters/knots))

east.disp = r * sin.a / cos.lat / 111111 # east displacement
north.disp = r * cos.a / 111111 # north displacement

if ( counter == 1 ) {
LOC.range = data.table("Lon" = (LOC[1] + east.disp),"Lat" = (LOC[2] + north.disp), "VMC" = VMC)#[which(VMC>0)] 
} else {
LOC.range = data.table("Lon" = (LOC[1] + east.disp),"Lat" = (LOC[2] + north.disp), "VMC" = VMC)#[which(distGeo(LOC.prior,]
}
#distance Sailed to new coordinate (lat and lon)
LOC = as.double(LOC.range[which.min(distGeo(LOC.range,MAN))])[1:2] # distGeo(LOC.prior,MAN)<
LOC.range$VMC <- NULL

# Generate plot and points
if ( counter == 1 ) {
plot(LOC.range, col = 'red', ylim = c(41.8,45), xlim = c(-87.8,-86.2), type = 'l', main = "Optimal VMC with Final Bearing", ylab = "Lat", xlab = "Lon") }
points(LOC.range, type ='l')
points(LOC.prior[1],LOC.prior[2], col = 'blue', pch= 17)
points(LOC, col = 'blue', pch= 4)

destination = rbind(destination, LOC.range[which.max(distGeo(LOC.prior,MAN)>distGeo(LOC.range,MAN))])

counter = counter + 1

}

map + geom_point(aes(x = Lon, y= Lat, size = 10), data = destination, alpha = .5)

counter*obs.interval # hours
counter*obs.interval/24 # days


