mu_init = 30; sigma_squared_init = 21

x = sample_data; alpha = 1; epsilon = 1e-5

mu_new <- mu_old <- mu_init  # Initialize variables
sigma_squared_new <- sigma_squared_old <- sigma_squared_init
convergence_flag <- TRUE; count <- 0; backtracks <- 0

mu_old <- mu_new; sigma_squared_old <- sigma_squared_new

theta_test
mu_new = theta_test[1]
sigma_squared_new = theta_test[2]
# (1) 30.04762 21.00227

step_old <- g0(x, mu_old, sigma_squared_old) # Calculate steps
step_new <- g0(x, mu_new, sigma_squared_new)

step_new - step_old
# -0.002542405
# -0.001203784
# -0.0005850358
# -0.0002883037
# -0.0001430983
# -7.128576e-05
# -3.557703e-05
# -1.777205e-05

# -2.219706e-06 (10x)
# -2.117417e-12 (30x)


z <- 2 * volcano        # Exaggerate the relief
x <- 10 * (1:nrow(z))   # 10 meter spacing (S to N)
y <- 10 * (1:ncol(z))   # 10 meter spacing (E to W)

# draw volcano and store transformation matrix
pmat <- persp(x, y, z, theta = 35, phi = 40, col = 'green4', scale = FALSE,
              ltheta = -120, shade = 0.75, border = NA, box = TRUE)

# take some xyz values from the matrix
s = sample(1:prod(dim(z)), size=500)
xx = x[row(z)[s] ]
yy = y[col(z)[s]]
zz = z[s] + 10

# depth calculation function (adapted from Duncan Murdoch at https://stat.ethz.ch/pipermail/r-help/2005-September/079241.html)
depth3d <- function(x,y,z, pmat, minsize=0.2, maxsize=2) {
  
  # determine depth of each point from xyz and transformation matrix pmat
  tr <- as.matrix(cbind(x, y, z, 1)) %*% pmat
  tr <- tr[,3]/tr[,4]
  
  # scale depth to point sizes between minsize and maxsize
  psize <- ((tr-min(tr) ) * (maxsize-minsize)) / (max(tr)-min(tr)) + minsize
  return(psize)
}

# determine distance to eye
psize = depth3d(xx,yy,zz,pmat,minsize=0.1, maxsize = 1)

# from 3D to 2D coordinates
mypoints <- trans3d(xx, yy, zz, pmat=pmat)

# plot in 2D space with pointsize related to distance
points(mypoints, pch=8, cex=psize, col=4)


points()