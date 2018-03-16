library(Rscattnlay)

scatdat <- function(rp,na,title) {
  S <- Scatterer()  # create a scatterer object
  
  np <- Layer()     # and a layer representing the np
  na(S) <- na   # set the ambient index
  r(np) <- rp     # and the size (radius) of the np in nm
  
  xsec <- pi*rp^2
  
  # load up silver data
  data_ag_vis <- data.frame(read.table("data/johnson_ag_vis.csv",sep=",", header=TRUE))
  colnames(data_ag_vis) <- c("lambda","n","k","eps_real","eps_imag")
  
  lambda_data <- data_ag_vis$lambda #convert to nm
  n_data <- data_ag_vis$n 
  k_data <- data_ag_vis$k
  
  scatnlay_ar <-function(k){
    m(np) <- spl_n$y[k]+spl_k$y[k]*(0+1i);
    lambda(S) <- lambda[k];
    St <- S+np;              # PUT THE STACK HERE
    Q <- scattnlay(St);
  }
  
  #interpolation as palik data is sparse
  n = 500 #number of points
  spl_n <- approx(lambda_data,n_data,n=n)
  spl_k <- approx(lambda_data,k_data,n=n)
  lambda = spl_n$x
  
  output <- sapply(1:n,scatnlay_ar)
  data_mat_np <- data.frame(t(output),lambda,title)
  
  colnames(data_mat_np) <- c("Qext", "Qsca", "Qabs", "Qbk", "Qpr", "g", "Albedo", "nmax","Lambda","Layers")
}


