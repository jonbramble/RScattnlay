library(Rscattnlay)
library(ggplot2)
library(reshape2)

S <- Scatterer()  # create a scatterer object

np <- Layer()     # and a layer representing the np
lipid <- Layer()  # and other layers
inner_layer <- Layer()
water <- Layer()

na(S) <- 1.33     # set the ambient index
r(np) <- 40       # and the size (radius) of the np in nm

# load up silver data
palik_ag_vis <- read.table("data//palik_ag_vis_hb.csv",sep=",", header=TRUE, nrows=49)
colnames(palik_ag_vis) <- c("lambda","n","k","eps_real","eps_imag")

lambda_palik <- palik_ag_vis$lambda #convert to nm
n_palik <- palik_ag_vis$n 
k_palik <- palik_ag_vis$k

#interpolation as palik data is sparse
n = 1024 #number of points
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
colnames(data_mat_np) <- c("Qext", "Qsca", "Qabs", "Qbk", "Qpr", "g", "Albedo", "nmax","Lambda","Layers")

data_mat_np$Lambda[which.max(data_mat_np$Qext)]

ext_plot <- ggplot(data=data_mat_np, aes(x=Lambda,y=Qext, group=as.factor(Layers), color=as.factor(Layers))) + 
  geom_line() + theme_minimal(base_size=22) + xlim(c(350,730))+ ylim(c(0,15)) + labs(colour = "Experiment")  + xlab("Wavelength (nm)")

ext_plot
#ggsave(file="test.svg", plot=ext_plot, width=10, height=8)
