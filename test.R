#scattnlay()

S <- Scatterer()  # create a scatterer object

np <- Layer()     # and a layer representing the np
lipid <- Layer()  # and other layers

na(S) <- 1.33     # set the ambient index

d(np) <- 40       # and the size of the np in nm

d(lipid) <- 42    # and the size of the other layers
m(lipid) <- 1.45+0i   #fixed value

# load up silver data
palik_ag_vis <- read.csv("palik_ag_vis.csv", header=FALSE)
colnames(palik_ag_vis) <- c("lambda","n","k")

lambda_palik <- palik_ag_vis$lambda*1000  #convert to nm
n_palik <- palik_ag_vis$n 
k_palik <- palik_ag_vis$k

#interpolation as palik data is sparse
n = 1024 #number of points
spl_n <- approx(lambda_palik,n_palik,n=n)
spl_k <- approx(lambda_palik,k_palik,n=n)

lambda = spl_n$x

Qf = matrix(nrow = n, ncol = 9, byrow=TRUE)
colnames(Qf) <- c("Qext", "Qsca", "Qabs", "Qbk", "Qpr", "g", "Albedo", "nmax","Lambda")

# array based for now
for(k in 1:n) {
 m(np) <- spl_n$y[k]+spl_k$y[k]*(0+1i);
 lambda(S) <- lambda[k];
 St <- S+np;
 Q <- scattnlay(St);
 Q[9] <- lambda[k];
 Qf[k,] <- Q;
}


Qf_l = matrix(nrow = n, ncol = 9, byrow=TRUE)
colnames(Qf_l) <- c("Qext", "Qsca", "Qabs", "Qbk", "Qpr", "g", "Albedo", "nmax","Lambda")

# array based for now
for(k in 1:n) {
  m(np) <- spl_n$y[k]+spl_k$y[k]*(0+1i);
  lambda(S) <- lambda[k];
  St <- S+np+lipid;
  Q <- scattnlay(St);
  Q[9] <- lambda[k];
  Qf_l[k,] <- Q;
}

plot(Qf[,"Lambda"],Qf[,"Qsca"],type="l",col="black")
lines(Qf_l[,"Lambda"],Qf_l[,"Qsca"],type="l",col="red")


