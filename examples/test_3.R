#scattnlay()
library(Rscattnlay)
library(ggplot2)
library(reshape2)

S <- Scatterer()  # create a scatterer object

np <- Layer()     # and a layer representing the np
na(S) <- 1.52   # set the ambient index
r(np) <- 20     # and the size (radius) of the np in nm

xsec <- pi*20^2

# load up silver data
palik_ag_vis <- data.frame(read.table("data/johnson_ag_vis.csv",sep=",", header=TRUE))
colnames(palik_ag_vis) <- c("lambda","n","k","eps_real","eps_imag")

lambda_palik <- palik_ag_vis$lambda #convert to nm
n_palik <- palik_ag_vis$n 
k_palik <- palik_ag_vis$k

#interpolation as palik data is sparse
n = 500#number of points
spl_n <- approx(lambda_palik,n_palik,n=n)
spl_k <- approx(lambda_palik,k_palik,n=n)
lambda = spl_n$x

scatnlay_ar <-function(k){
  m(np) <- spl_n$y[k]+spl_k$y[k]*(0+1i);
  lambda(S) <- lambda[k];
  St <- S+np;              # PUT THE STACK HERE
  Q <- scattnlay(St);
}

output <- sapply(1:n,scatnlay_ar)
data_mat_np <- data.frame(t(output),lambda,'water')

#Data1 <- data.frame(read.table("/media/mbajb/data/Linux/Comsol Projects/Nottingham Simulations/Data8.txt",skip=5))
#colnames(Data1) <- c("Lambda","Frequency","Qscat")
#Data1$scal <- Data1$Qscat/max(Data1$Qscat)
#Data1$Lambda <- Data1$Lambda*1e9

colnames(data_mat_np) <- c("Qext", "Qsca", "Qabs", "Qbk", "Qpr", "g", "Albedo", "nmax","Lambda","Layers")

data_mat <- subset(data_mat_np,Lambda<700 & Lambda >300)
#data_mat$Qsca <- data_mat$Qsca/max(data_mat$Qsca)


palik_ag_melt <- melt(subset(palik_ag_vis,select = c("lambda","n","k")), id=c("lambda"),variable.name = "type", 
                      value.name = "Index")

scat_plot <- ggplot(data=data_mat, aes(x=Lambda,y=Qsca*xsec/1000)) 
scat_plot + geom_line() + theme_minimal(base_size=22) + xlim(c(300,550)) + xlab("Wavelength (nm)") + ylab("Scattering Cross section")

scat_plot2 <- ggplot(data=data_mat, aes(x=Lambda,y=Qabs*xsec/1000)) 
scat_plot2 + geom_line() + theme_minimal(base_size=22) + xlim(c(300,550)) + xlab("Wavelength (nm)") + ylab("Absorption Cross section")


#scat_plot2 <- ggplot(data=Data1, aes(x=Lambda,y=scal)) 
#scat_plot2 + geom_line() + theme_minimal(base_size=22) + xlim(c(350,700))+ ylim(c(0,1.2)) + xlab("Wavelength (nm)")

