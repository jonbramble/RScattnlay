#scattnlay()

S <- Scatterer()

l <- Layer()

lambda(S) <- 400

d(l) <- 40
m(l) <- 1.33+0i;

St <- S+l

scattnlay(St)
