colnames(dat) <- c("Qext", "Qsca", "Qabs", "Qbk", "Qpr", "g", "Albedo", "nmax","Lambda","Layers")
data_mat <- subset(dat,Lambda<700 & Lambda >300)

xsec <- (pi*20*20)

ylabsca <- expression('Scattering Cross Section nm'^2 )
ylababs <- expression('Absorption Cross Section nm'^2 )

scaplot <- ggplot(data=data_mat, aes(x=Lambda,y=Qsca*xsec,color=Layers))
scaplot + geom_line() + theme_minimal() + ylab(ylabsca) + xlim(c(300,550)) + xlab("Wavelength (nm)")

absplot <- ggplot(data=data_mat, aes(x=Lambda,y=Qabs*xsec,color=Layers))
absplot + geom_line() + theme_minimal() + ylab(ylababs) + xlim(c(300,550)) + xlab("Wavelength (nm)")
