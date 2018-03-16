library(Rscattnlay)
library(ggplot2)
library(reshape2)

S <- Scatterer()  # create a scatterer object

dr <- Layer()     # and a layer representing the droplet
na(S) <- 1.33    # set the ambient index
lambda(S) <- 600
nt(S) <- 1


m(dr) <- 1.431+0i;

scatnlay_ar <- function(r){
  r(dr) <- r
  St <- S+dr;              # PUT THE STACK HERE
  Q <- scattnlay(St);
}

dr_size <- c(1:100)

output <- sapply(dr_size,scatnlay_ar)

data_mat <- data.frame(t(output),dr_size,'oil')
colnames(data_mat) <- c("Qext", "Qsca", "Qabs", "Qbk", "Qpr", "g", "Albedo", "nmax","Size","Layers")

scat_plot <- ggplot(data=data_mat, aes(x=Size,y=Qsca)) 
scat_plot + geom_line() + theme_minimal(base_size=22) + xlim(c(0,100)) + xlab("Radius (nm)")
