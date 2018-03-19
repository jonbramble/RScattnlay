library(plotrix)
S <- Scatterer()  # create a scatterer object

np <- Layer()     # and a layer representing the np
lipid <- Layer()  # and other layers

na(S) <- 1.33     # set the ambient refractive index eg water 
r(np) <- 20       # and the size (radius) of the np in nm

r(lipid) <- 21        # and the size (radius) of the other layers
m(lipid) <- 1.45+0i   # fixed value as a complex number

# load up silver data from package or own source
palik_ag_vis <- read.table("data/palik_ag_vis_hb.csv",sep=",", header=TRUE)
colnames(palik_ag_vis) <- c("lambda","n","k","eps_real","eps_imag")

# extract the n and k values
lambda_palik <- palik_ag_vis$lambda #convert to nm
n_palik <- palik_ag_vis$n 
k_palik <- palik_ag_vis$k

# interpolation as palik data is sparse using approx
n = 1024 #number of points
spl_n <- approx(lambda_palik,n_palik,n=n)
spl_k <- approx(lambda_palik,k_palik,n=n)
lambda = spl_n$x

k <- 200
m(np) <- spl_n$y[k]+spl_k$y[k]*(0+1i)
lambda(S) <- lambda[k]

nt(S) <- 360
tf(S) <- 360

St <- S+np+lipid
Q <- amplitudes(St)

I <- abs(Q$S1)^2 + abs(Q$S2)^2
angles<-Q$Theta
radial.plot(I,angles,main="Mie Angular Scattering",show.grid.labels=4,show.grid=TRUE,labels=NA,rp.type = "p")

