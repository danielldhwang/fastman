# Create a test dataset
data <- data.frame(CHR=c(rep(c(1:22),1000)),
                   BP=c(round(runif(22000,1,4*10^7))),
                   P=c(runif(10000,0.1,1),runif(8000,0.01,0.1),runif(2000,0.001,0.01),runif(1000,0.0001,0.001),runif(900,1e-4,1e-3),runif(100,1e-8,1e-4)))
fastman(data)
