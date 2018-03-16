qtable <- read.table("/home/mbzjpb/Programming/R/ddscat/qtable", header=TRUE, skip=15, quote="\"")
#qtable2 <- read.table("/home/mbzjpb/Programming/R/ddscat/qtable_ext", header=TRUE, skip=15, quote="\"")
cols <- c("aeff","wavelength","Qext","Q_abs","Q_sca","g(1)=<cos>","<cos^2>","Q_bk","Nsca")
df <- data.frame(qtable)
colnames(df) <- cols

#df2 <- data.frame(qtable2)
#colnames(df2) <- cols

#df<-rbind(df,df2)
df$lambda = df$wavelength*1000

p<-ggplot(data=df,aes(x=lambda,y=Qext))
p + geom_point()+ geom_line(data=data_mat_np) + theme_minimal(base_size=22) + xlim(350,650)+ xlab("Wavelength (nm)")
