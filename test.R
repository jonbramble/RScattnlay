#scattnlay()
library(Rscattnlay)
library(ggplot2)

S <- Scatterer()  # create a scatterer object

np <- Layer()     # and a layer representing the np
lipid <- Layer()  # and other layers
water <- Layer()
na(S) <- 1.33     # set the ambient index
d(np) <- 40       # and the size (radius) of the np in nm

d(water) <- 41
m(water) <- 1.33+0i

d(lipid) <- 41    # and the size (radius) of the other layers
m(lipid) <- 1.40+0i   #fixed value

# load up silver data
palik_ag_vis <- read.table("palik_ag_vis_hb.csv", header=TRUE)
colnames(palik_ag_vis) <- c("lambda","n","k","eps_real","eps_imag")

lambda_palik <- palik_ag_vis$lambda #convert to nm
n_palik <- palik_ag_vis$n 
k_palik <- palik_ag_vis$k

#interpolation as palik data is sparse
n = 512 #number of points
spl_n <- approx(lambda_palik,n_palik,n=n)
spl_k <- approx(lambda_palik,k_palik,n=n)
lambda = spl_n$x

scatnlay_ar <-function(k){
  m(np) <- spl_n$y[k]+spl_k$y[k]*(0+1i);
  lambda(S) <- lambda[k];
  St <- S+np+water;              # PUT THE STACK HERE
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
data_mat_np <- data.frame(t(output),lambda,'water')
data_mat_lipid <- data.frame(t(output_lipid),lambda,'lipid')

colnames(data_mat_np) <- c("Qext", "Qsca", "Qabs", "Qbk", "Qpr", "g", "Albedo", "nmax","Lambda","Layers")
colnames(data_mat_lipid) <- c("Qext", "Qsca", "Qabs", "Qbk", "Qpr", "g", "Albedo", "nmax","Lambda","Layers")

data_mat <- rbind(data_mat_np,data_mat_lipid)

ext <- subset(data_mat_lipid, Lambda < 700, select= c("Qext","Lambda"))
ext$Lambda[which.max(ext$Qext)]


palik_ag_melt <- melt(subset(palik_ag_vis,select = c("lambda","n","k")), id=c("lambda"),variable.name = "type", 
                      value.name = "Index")

ag_plot <- ggplot(data=palik_ag_melt, aes(x=lambda,y=Index, group=type, color=type))
ag_plot + geom_line() + theme_minimal(base_size=22) + xlim(c(300,700)) + ylim(c(0,6)) + xlab("Wavelength (nm)")

scat_plot <- ggplot(data=data_mat, aes(x=Lambda,y=Qext, group=as.factor(Layers), color=as.factor(Layers))) 
scat_plot + geom_line() + theme_minimal(base_size=22) + xlim(c(300,700)) + labs(colour = "Layers") 
