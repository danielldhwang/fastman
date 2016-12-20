#' fastqq
#'
#' Creates Q-Q plots from GWAS summaries.
#' @param data A GWAS summary with at least a column of P-value. P-values have to be numeric.
#' @param p A string for the header for the column of P-value in data.
#' @param lambda Set FALSE to exclude lambda from the plot.
#' @param main A string for the title of the plot.
#' @importFrom stats median qchisq ppoints
#' @export
fastqq <- function(data, p="P", lambda = T, main = "Q-Q plot"){
    # Check for sensible dataset ** these steps are adapted from qqman **
    ## Make sure you have p column.
    if (!(p %in% names(data))) stop(paste("Column", p, "not found!"))
    ## make sure p column is numeric.
    if (!is.numeric(data[[p]])) stop(paste(p, "column should be numeric."))

    # Calculate lambda based on the median of P-values.
    if(lambda){
        lambda_to_print<-qchisq((1-median(data[,p])),1)/qchisq(0.5,1)
    }
    # Generate expected P-values.
    data$exp <- NA
    data <- data[order(data[,p]),]
    data$exp <- ppoints(nrow(data) )
    # Restrict the number of SNPs with P > 0.1 to ~20000.
    nrow_0.1 <- nrow(data[which(data[,p]>0.1),])
    i1 <- 1
    n_0.1 <- 1
    while(i1<10000){
        if(nrow_0.1/i1 <= 200000){
            n_0.1 <- i1
            break()
        }
        i1 <- i1+1
    }
    temp1 <- rep(c(1:n_0.1), nrow(data)/n_0.1)
    temp1 <- temp1[1:nrow(data)]
    data$gp_0.1 <- temp1
    data <- data[which(data[,p]<=0.1 | data[,"gp_0.1"]==1),]

    # Restrict the number of SNPs with P <=0.1 and P > 0.01 to ~20000.
    nrow_0.01 <- nrow(data[which(data[,p]<=0.1 & data[,p]>0.01),])
    if(nrow_0.01 > 200000){
        i2 <- 1
        n_0.01 <- 1
        while(i2<10000){
            if(nrow_0.01/i2 <= 200000){
                n_0.01 <- i2
                break()
            }
            i2 <- i2+1
        }
        temp2 <- rep(c(1:n_0.01), nrow(data)/n_0.01)
        temp2 <- temp2[1:nrow(data)]
        data$gp_0.01 <- temp2
        data <- data[which(data[,p]>0.1 | data[,p]<=0.01 | data[,"gp_0.01"]==1),]
    }
    # Restrict the number of SNPs with P <=0.01 and P > 0.001 to ~20000.
    nrow_0.001 <- nrow(data[which(data[,p]<=0.01 & data[,p]>0.001),])
    if(nrow_0.001 > 200000){
        i3 <- 1
        n_0.001 <- 1
        while(i2<10000){
            if(nrow_0.001/i3 <= 200000){
                n_0.001 <- i3
                break()
            }
            i3 <- i3+1
        }
        temp3 <- rep(c(1:n_0.001), nrow(data)/n_0.001)
        temp3 <- temp3[1:nrow(data)]
        data$gp_0.001 <- temp3
        data <- data[which(data[,p]>0.01 | data[,p]<=0.001 | data[,"gp_0.001"]==1),]
    }
    # Calculate minus log P and plot.
    minus_log_P <- -log10(data[,p])
    exp_minus_log_P <- -log10(data$exp)
    plot(exp_minus_log_P, minus_log_P, pch = 16, cex = 0.7, xaxs="i", yaxs="i",
         xlim=c(0,ceiling(max(exp_minus_log_P))), ylim=c(0,ceiling(max(minus_log_P))),
         xlab = expression(Expected -log[10](P-value)), ylab = expression(Observed -log[10](P-value)))
    abline(0,1,col="red")
    mtext(main, side = 3, line = 1)
    # Add lamdba information.
    if(lambda){
        legend("topleft", legend = bquote(lambda == .(lambda_to_print)), bty = "n")
    }
}
