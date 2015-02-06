#scattnlay()

S <- Scatterer()

np <- Layer()
lipid <- Layer()

lambda(S) <- 400
na(S) <- 1.00

d(np) <- 40
m(np) <- 1.33+0i

d(lipid) <- 42
m(lipid) <- 1.45+0i

St <- S+np+lipid

scattnlay(St)
