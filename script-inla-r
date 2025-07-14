################################################################################
################# Spatial point pattern analysis using inlabru #################
################################################################################

# In this script the Chorley-Ribble dataset is analysed by performing
# a multivariate point pattern analysis. The lung cancer cases are 
# assumed as the controls and the exposure to an old inicinerator is
# included in the model as a LINEAR effect of the distance from every
# point to the location of this old incinerator.

#install.packages('INLA')
#install.packages("INLA", repos = INLA:::getOption("inla.binary.repo"))
install.packages("INLA", repos = "https://inla.r-inla-download.org/R/stable")
install.packages("sn")

#Load necessary packages
library(sp)
library(INLA)
library(inlabru)
library(spatstat)
library(RColorBrewer)
library(ggmap)
library(splancs)
library(fmesher)
library(sn)
library(sf)

# load the Rdata

load("Mult_Exp_Linear_code_artigo.RData")

#Step 0: Load, visualize and prepare the data. 

#Load data
data(southlancs)

#Plot the data
#Display a raw plot of the data
plot(southlancs.bdy, col="darkgreen")
points(southlancs.pts, col="black")
points(old.incinerator, col="red3", pch=17)


#Change the scale from m to km dividing by 1000 the coordinates
sl.pts <- as.data.frame(southlancs.pts/1000)
sl.bdy <- as.data.frame(southlancs.bdy/1000)
old.inc <- as.data.frame(old.incinerator/1000)

#Create a polygon with the boundary points
bdy.orig <- cbind(rev(sl.bdy[, 1]), rev(sl.bdy[, 2])) 
bdy.SP <- as(sf::st_polygon(list(bdy.orig[c(seq_len(NROW(bdy.orig)), 1), ])),
             "Spatial")

#Separate the two different cancer type points
sl.lun <- sl.pts[1:917, ]
sl.lar <- sl.pts[918:974, ]

#Build spatial object
coordinates(sl.pts) <- ~ arx + ary
coordinates(sl.lun) <- ~ arx + ary
coordinates(sl.lar) <- ~ arx + ary
coordinates(sl.bdy) <- ~ V1 + V2
coordinates(old.inc) <- ~ V1 + V2

#Step 1: Build the mesh
#Do not what to do? A useful tool is meshbuilder(), check it.

#Define the mesh

mesh <- fm_mesh_2d_inla(
  boundary = list(bdy.SP, NULL),
  cutoff = 0.1,
  max.edge = c(0.4, 0.8),
  min.angle = 27,
  offset = c(0.5, 2),
  n=c(16,16))   ### O Q SIGNIFICA ESSE NUM INICIAL DE NÓS????

#Number of mesh nodes
nv <- mesh$n #5556

#Plot together boundary, mesh and all the points



# ggplot() + gg(mesh) + gg(sl.lun, col="dodgerblue") +
#   gg(sl.lar, col="orangered")  +
#   gg(bdy.SP, col="darkgreen") +
#   gg(old.inc, col="red3", size=3) + coord_fixed(ratio = 1) ##### DEU ERRO!!!!


ggplot() +
  gg(mesh) +
  gg(sl.lun, col="dodgerblue") +
  gg(sl.lar, col="orangered") +
  gg(bdy.SP, col="green") +
  gg(old.inc, col="red3", size=3) +
  coord_sf()



#Step 2: Compute the distance to the old incinerator

#Calculate distance between the old inicinerator and the mesh points
grid.pts <- mesh$loc[, 1:2]
grid.dist <- spDistsN1(grid.pts, old.inc); head(grid.dist)

#Change the mesh to spatialpixels object
pix.cov.grid <- fm_pixels(mesh, format = "sp")

#Interpolate linearly on each triangle of the mesh
pix.cov <- SpatialPixelsDataFrame(
  pix.cov.grid@coords,
  data = data.frame(cov1 = fmesher::fm_evaluate(
    mesh = mesh,
    loc = pix.cov.grid,
    field = grid.dist)))

#Define a functions which inlabru uses to evaluate the covariate
f.cov1 <- function(where) {
  # Extract the values
  v <- eval_spatial(pix.cov, where, layer = "cov1")
  # Fill in missing values
  if (any(is.na(v))) {
    v <- bru_fill_missing(pix.cov, where, v)
  }
  return(v)
}



#Step 3: Create the set of basis for the spde approximation

#Create spde basis functions using PC-priors
pcmatern <- inla.spde2.pcmatern(
  # Mesh and smoothness parameter
  mesh = mesh, alpha = 2,
  # P(practic.range < 5) = 0.95
  prior.range = c(5, 0.95),
  # P(sigma > ln(10)) = 0.01
  prior.sigma = c(log(10), 0.01))

#Step 4: Specify the components of the model and the likelihoods.

#Define the components: A common (conSpde) and the effect of the 
# old inicinerator as a linear effect are added here.

cmp <- coordinates ~ -1 + Inter.con(1) + Inter.lar(1) +
  sharedspde(coordinates, model = pcmatern) +
  larspde(coordinates, model = pcmatern) +
  cov.lar(f.cov1(.data.), model="linear")


#Specify the stracutre of the likelihoods
#One likelihood for each point pattern 

# Likelihood for the controls (lung cancer)
con.lik <- like(
  family = "cp",
  formula = coordinates ~ Inter.con + sharedspde,
  samplers = bdy.SP,
  data = sl.lun,
  domain = list(coordinates = mesh)
)

# Likelihood for the cases (larynx cancer)
lar.lik <- like(
  family = "cp",
  formula = coordinates ~ Inter.lar + sharedspde + larspde + cov.lar,
  samplers = bdy.SP,
  data = sl.lar,
  domain = list(coordinates = mesh)
)

#Step 5: Fit the model

#Fit the model
t1 <- Sys.time();t1
fit <- bru(cmp,
  con.lik, lar.lik,
  options = list(
    control.inla = list(int.strategy = "eb"),
    verbose = FALSE
  ))
t2 <- Sys.time(); t2
fit.time<- t2 - t1; fit.time


summary(fit)

# Fixed effects:
#   mean    sd 0.025quant 0.5quant 0.975quant   mode kld
# Inter.con -1.502 0.819     -3.107   -1.502      0.104 -1.502   0
# Inter.lar -4.261 0.883     -5.991   -4.261     -2.531 -4.261   0
# cov.lar   -0.005 0.035     -0.074   -0.005      0.064 -0.005   0
# 
# Random effects:
#   Name	  Model
# sharedspde SPDE2 model
# larspde SPDE2 model
# 
# Model hyperparameters:
#   mean    sd 0.025quant 0.5quant 0.975quant  mode
# Range for sharedspde 6.075 1.092      4.282    5.956      8.566 5.684
# Stdev for sharedspde 2.421 0.387      1.775    2.383      3.293 2.292
# Range for larspde    0.358 0.472      0.038    0.218      1.555 0.094
# Stdev for larspde    0.229 0.225      0.023    0.162      0.827 0.065
# 
# Deviance Information Criterion (DIC) ...............: -4200.10
# Deviance Information Criterion (DIC, saturated) ....: -4206.60
# Effective number of parameters .....................: -6370.91
# 
# Watanabe-Akaike information criterion (WAIC) ...: 2524.09
# Effective number of parameters .................: 177.00
# 
# Marginal log-Likelihood:  -4611.53

#Define where to estimate
where <- fm_pixels(mesh, mask = bdy.SP, format = "sp")

#Estimate the intensity of the lung cancer cases
t3 <- Sys.time();t3
lambda.con <- predict(fit, where, ~exp(Inter.con + sharedspde))
t4 <- Sys.time(); p.time<- t4 - t3; p.time

#Estimate the intensity of the larynx cancer cases
t3 <- Sys.time();t3
lambda.lar <- predict(fit, where, ~exp(Inter.lar + sharedspde + larspde + cov.lar ))
t4 <- Sys.time(); p.time.lun <- t4 - t3; p.time.lun

#Estimate the common spatial effect
t3 <- Sys.time();t3
sp.eff.con <- predict(fit, where, ~sharedspde)
t4 <- Sys.time(); p.time<- t4 - t3; p.time

#Estimate the linear exposure effect of dist 
t3 <- Sys.time();t3
cov.eff.lar <- predict(fit, where, ~ cov.lar)
t4 <- Sys.time(); p.time<- t4 - t3; p.time

#Quick multiplot of the results
multiplot(
  ggplot() + gg(lambda.con, mask = bdy.SP) + coord_equal() + gg(bdy.SP),
  ggplot() + gg(lambda.lar, mask = bdy.SP) + coord_equal() + gg(bdy.SP),
  ggplot() + gg(sp.eff.con, mask = bdy.SP) + coord_equal() + gg(bdy.SP),
  ggplot() + gg(cov.eff.lar, mask = bdy.SP) + coord_equal() + gg(bdy.SP),
  cols=2)

## fazendo os graficos em separado

 ggplot() + gg(lambda.con, mask = bdy.SP) + coord_equal() + gg(bdy.SP)
# 
 class(lambda.con) # arq sp
# 
 ggplot() + gg(lambda.lar, mask = bdy.SP) + coord_equal() + gg(bdy.SP)
# 
 class(lambda.lar) # "bru_prediction" "data.frame"  
# 
 ggplot() + gg(sp.eff.con, mask = bdy.SP) + coord_equal() + gg(bdy.SP)
# 
 class(sp.eff.con) # arq sp
# 
 ggplot() + gg(cov.eff.lar, mask = bdy.SP) + coord_equal() + gg(bdy.SP) 
# 
 #class(cov.eff.lar) "bru_prediction" "data.frame" 

#Plot the estimated trend of the covariate effect in 1D
ggplot(
  data = cbind(data.frame(dist = over(where, pix.cov)[, 1, drop = TRUE]),
    cov.eff.lar@data)) +
  geom_line(aes(dist, mean)) +
  geom_ribbon(aes(dist, ymin = q0.025, ymax = q0.975,
    alpha = 0.5))

## não achou o 'pix.cov'

#Save the results
save.image("Mult_Exp_Linear_code_artigo.RData")

rm(list=ls())
