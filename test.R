#scattnlay()
library(Rscattnlay)
library(ggplot2)

S <- Scatterer()  # create a scatterer object

np <- Layer()     # and a layer representing the np
lipid <- Layer()  # and other layers
na(S) <- 1.33     # set the ambient index
d(np) <- 40       # and the size of the np in nm
d(lipid) <- 45    # and the size of the other layers
m(lipid) <- 1.45+0i   #fixed value

# load up silver data
palik_ag_vis <- read.csv("palik_ag_vis.csv", header=FALSE)
colnames(palik_ag_vis) <- c("lambda","n","k")

lambda_palik <- palik_ag_vis$lambda*1000  #convert to nm
n_palik <- palik_ag_vis$n 
k_palik <- palik_ag_vis$k

#interpolation as palik data is sparse
n = 256 #number of points
spl_n <- approx(lambda_palik,n_palik,n=n)
spl_k <- approx(lambda_palik,k_palik,n=n)
lambda = spl_n$x

scatnlay_ar <-function(k){
  m(np) <- spl_n$y[k]+spl_k$y[k]*(0+1i);
  lambda(S) <- lambda[k];
  St <- S+np;              # PUT THE STACK HERE
  Q <- scattnlay(St);
}

scatnlay_ar_lipid <-function(k){
  m(np) <- spl_n$y[k]+spl_k$y[k]*(0+1i);
  lambda(S) <- lambda[k];
  St <- S+np+lipid;              # PUT THE STACK HERE
  Q <- scattnlay(St);
}

output <- sapply(1:n,scatnlay_ar)
output_lipid <- sapply(1:n,scatnlay_ar_lipid)
data_mat_np <- data.frame(t(output),lambda,1)
data_mat_lipid <- data.frame(t(output_lipid),lambda,2)

colnames(data_mat_np) <- c("Qext", "Qsca", "Qabs", "Qbk", "Qpr", "g", "Albedo", "nmax","Lambda","Layers")
colnames(data_mat_lipid) <- c("Qext", "Qsca", "Qabs", "Qbk", "Qpr", "g", "Albedo", "nmax","Lambda","Layers")

data_mat <- rbind(data_mat_np,data_mat_lipid)

scat_plot <- ggplot(data=data_mat, aes(x=Lambda,y=Qsca, group=as.factor(Layers), color=as.factor(Layers)))
scat_plot + geom_line() + theme_minimal()
