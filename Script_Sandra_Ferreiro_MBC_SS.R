#Package used
install.packages("spatstats")
library(spatstats)

#Data inspection
data("nztrees") 
plot(nztrees, main = NULL) # to plot the ppp (planar point pattern)
str(nztrees) # to see data (coordinates points, spatial region description (window))
summary(nztrees) # points´ number, average intensity per unit area, window dimension, area and type and units
nz<-rescale(nztrees,2.87356,"meter") #rescaling the data. The original data is divided by the conversion factor (2.87356)
summary(nz) # points´ number, average intensity per unit area, window dimension, area and type and units of the rescaled sample


#Is the data homogenous?
# intensity (homogeneous process?)
intensity(nz)
# standard error of the previous measurement (points per sq meter)
lam <- intensity(nz)
sqrt(lam/area(Window(nz)))

# quadrat counting to check homogeinity 
Q<- quadratcount(nz,nx=2, ny=2) # x and y are the coordinates where the points are being evaluated
Q #similar number of points per region
Q2<- quadratcount(nz,nx=4,ny=4)
Q2 #in this case the variance between regions is higher

intensity(Q) # average intensity in each quadrat (4 quadrants)
plot(intensity(Q, image = TRUE)) # average intensity per quadrat plot
intensity(Q2)# average intensity in each quadrat (16 quadrants)
plot(intensity(Q2,image = TRUE))

quadrat.test(nz)#to check if our point pattern is an homogeneous Poisson process (CRS) or inhomogeneus (alternative hypothesis)

#To visualize the observed vs expected counts for 25 quadrats
plot(quadrat.test(nz), add = TRUE) #per quadrat: the upper left number is the observed counts  and the upper right number is the expected counts. The number of the bottom is the Pearson residual

#To avoid warning messages, use 4 quadrants (2x2 grid)
Q4 <-quadrat.test(nz, nx = 2, ny = 2) 
Q4

#To visualize the intensity per quadrat + observed vs expected counts for 4 quadrats
plot(intensity(Q, image = TRUE))
plot(Q4, add = TRUE)
# p-value = 0.2093, Ho still being rejected. No homogenous

# Gaussian density
d <- density(nz, edge = TRUE, kernel="gaussian") 
plot(d, main = "Density") # to visualize density


#IS THE PATTERN REGULAR?
fryplot(nz, main = NULL)
fryplot(nz, main = NULL, width=10) #to zoom the plot. A rectangular width of 10 is being observed
#When zoom it, the point pattern seems regular. Subjetive interpretation anyway

# The K FUNCTION is the cumulative average number of data points lying within a distance r of a typical data point
K <- Kest(nz)
plot(K, main=NULL, legend = FALSE) #To visualize the K function
K #To obtain more info abount the different estimators of the emperical K function

#Pointwise envelopes for K
a <- capture.output(plot (envelope(nz, Kest), main = NULL, legend = FALSE))
a <- capture.output(e <- envelope(nz, Kest, nsim = 39)) # to capture the pointwise critical envelopes for K(r)
e #to see details of the test

# Global envelopes test for K
set.seed(0)
a <- capture.output(plot(envelope(nz, Kest, nsim = 19, global = TRUE), main = NULL, legend = FALSE))
# the black line is inside the envelope

# L FUNCTION 
plot(Lest(nz), main = NULL)
# global envelopes for L
set.seed(0)
a <- capture.output(plot(envelope(nz, Lest, nsim = 19, global = TRUE), main = NULL, legend = FALSE)) 

# F FUNCTION (empty-space function)
Fc <- Fest(nz)
plot(Fc, main = NULL) # the experimental curves differ from the theoretical CSR curve (blue dolled-line), and suggest a regular pattern
#Pointwise envelope for F
a <- capture.output(plot(envelope(nz, Fest, nsim = 39), main = NULL, legend = FALSE)) # simulation with envelopes
# the black line only falls within the envelope of F for values of r < 4

# G FUNCTION 
Gc <- Gest(nz) 
plot(Gc, main = NULL)

# J FUNCTION 
Jc <- Jest(nz)
plot(Jc, main = NULL, legend = FALSE) 

