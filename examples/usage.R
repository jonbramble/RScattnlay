library(Rscattnlay)

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

# write a function to loop over lambda
scatnlay_ar_lipid <-function(k){
  m(np) <- spl_n$y[k]+spl_k$y[k]*(0+1i);
  lambda(S) <- lambda[k];
  St <- S+np+lipid;              # PUT THE STACK HERE
  Q <- scattnlay(St);
}

output_lipid <- sapply(1:n,scatnlay_ar_lipid)

data_mat_lipid <- data.frame(t(output_lipid),lambda,'lipid')
colnames(data_mat_lipid) <- c("Qext", "Qsca", "Qabs", "Qbk", "Qpr", "g", "Albedo", "nmax","Lambda","Layers")