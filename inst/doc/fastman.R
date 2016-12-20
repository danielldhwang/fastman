#' fastman
#'
#' Creates Manhattan plots from GWAS summaries.
#' @param data A GWAS summary with columns of chromosome number, position (base pair) information, and P-value.
#' @param chr A string for the header of the column of chromosome number in data.
#' @param ps A string for the header of the column of position (base pair) information in data.
#' @param p A string for the header for the column of P-value in data.
#' @param main A string for the Title of the plot.
#' @param suggest_line Set FALSE to remove the suggestive threshold line.
#' @param gws_line Set FALSE to remove the genome-wide significant line.
#' @param color A string vector containg the colors to be used.
#' @param chr_build A string describing what Genetic Reference Consortium build is used. It has to be either "GRCh37" or "CRCh38".
#' @param yscale A numeric value for y scale.
#' @param xlab_all Set TRUE to show labels for all chromosomes.
#' @param turbo Set TRUE to remove all rows with P <= 0.1.
#' @importFrom graphics abline axis mtext plot legend
#' @export
fastman <- function(data, chr="CHR", ps="BP", p="P", main="Manhattan plot", suggest_line=-log10(1e-5), gws_line=-log10(5e-8), color=c("grey","black"), chr_build="GRCh37", yscale=NA, xlab_all=F, turbo=F){
    # Check for sensible dataset ** these steps are adapted from qqman **
    ## Make sure you have chr, bp and p columns.
    if (!(chr %in% names(data))) stop(paste("Column", chr, "not found!"))
    if (!(ps %in% names(data))) stop(paste("Column", ps, "not found!"))
    if (!(p %in% names(data))) stop(paste("Column", p, "not found!"))
    ## make sure chr, bp, and p columns are numeric.
    if (!is.numeric(data[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(data[[ps]])) stop(paste(ps, "column should be numeric."))
    if (!is.numeric(data[[p]])) stop(paste(p, "column should be numeric."))

        # When turbo is TRUE, remove rows with P > 0.1, otherwise restrict the number of SNPs with P > 0.1 to ~20000.
    if(turbo){
        data <- data[which(data[,p]<=0.1),]
    }else{
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
    }
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
    # Select whether the SNPs position info is based on GRCh38 or GRCh38
    if(chr_build=="GRCh37"){
        chr_base <- c(0, 249250621, 492449994, 690472424, 881626700, 1062541960, 1233657027, 1392795690, 1539159712, 1680373143, 1815907890, 1950914406, 2084766301, 2199936179, 2307285719, 2409817111, 2500171864, 2581367074, 2659444322, 2718573305, 2781598825, 2829728720, 2881033286, 3036303846, 3095677412)
        chr_label <- c(124625310.5, 370850307.5, 591461209, 786049562, 972084330, 1148099494, 1313226359, 1465977701, 1609766428, 1748140517, 1883411148, 2017840354, 2142351240, 2253610949, 2358551415, 2454994488, 2540769469, 2620405698, 2689008814, 2750086065, 2815663773, 2865381003, 2958668566, 3065990629)
    }
    if(chr_build=="GRCh38"){
        chr_base <- c(0, 248956422, 491149951, 689445510, 879660065, 1061198324, 1232004303, 1391350276, 1536488912, 1674883629, 1808681051, 1943767673, 2077042982, 2191407310, 2298451028, 2400442217, 2490780562, 2574038003, 2654411288, 2713028904, 2777473071, 2824183054, 2875001522, 3031042417, 3088269832)
        chr_label <- c(124478211, 370053186.5, 590297730.5, 784552787.5, 970429194.5, 1146601314, 1311677290, 1463919594, 1605686271, 1741782340, 1876224362, 2010405328, 2134225146, 2244929169, 2349446623, 2445611390, 2532409283, 2614224646, 2683720096, 2745250988, 2800828063, 2849592288, 2953021970, 3059656125)
    }

    # Assign chromosome labels and their positions
    chr_n <- c(1:22)
    xlabel <- c(1:22)
    x_at <- chr_label[1:22]

    # Check if any SNP is on sex chromosomes.
    if(24 %in% data[,chr]){
        chr_n <- c(chr_n,23,24)
        xlabel <- c(xlabel,"X","Y")
        x_at <- chr_label[1:24]
    } else if(23 %in% data[,chr]){
        chr_n <- c(chr_n,23)
        xlabel <- c(xlabel,"X")
        x_at <- chr_label[1:23]
    }

    # Calculate the positions to plot for SNPs on chromosome 2 and after.
    data$ps2 <- NA
    for (i in chr_n){
        data[which(data[,chr]==i),"ps2"] <- data[which(data[,chr]==i),ps] + chr_base[i]
    }
    # Set colors
    data$colors <- NA
    for (i in chr_n){
        temp_n <- i%%length(color)
        if(temp_n==0){
            col<-color[length(color)]
        }else col<-color[temp_n]
        data[which(data[,chr]==i),"colors"] <- col
    }
    # If yscale is not provided, set the yscale to default of 10 when the strongest association < 1e-10, otherwise set it to match the strongest p.
    if(!is.na(yscale)){
        y_max <- yscale
        if(y_max<=10){
            y_at <- c(0:10)
        }else if(y_max<=20){
            y_at <- c(seq(0,y_max,by=2))
        }else if(y_max<=50){
            y_at <- c(seq(0,y_max,by=5))
        }else if(y_max<=100){
            y_at <- c(seq(0,y_max,by=10))
        }else {y_at <- c(seq(0,y_max,by=20))}
    }else if(is.na(yscale) & -log10(min(data[,p])) > 10){
        y_max <- ceiling(-log10(min(data[,p])))
        if(y_max<=10){
            y_at <- c(0:10)
        }else if(y_max<=20){
            y_at <- c(seq(0,y_max,by=2))
        }else if(y_max<=50){
            y_at <- c(seq(0,y_max,by=5))
        }else if(y_max<=100){
            y_at <- c(seq(0,y_max,by=10))
        }else {y_at <- c(seq(0,y_max,by=50))}
    }else {
        y_max <- 10
        y_at <- c(0:10)
    }
    if(y_at[length(y_at)]!=y_max){
        y_at <- c(y_at, y_max)
    }
    # Assign y scale labels based on turbo option
    y_min <- 0
    if(turbo){
        y_min <- 1
        y_at <- y_at[-1]
    }
    # Calculate the -logP and plot.
    data$minus_log10_p <- -log10(as.numeric(as.character(data[,p])))
    plot(minus_log10_p~ps2, data, axes = F, ann = F, col=data$colors, pch = 16, cex = 0.7,
         xlim=c(0,chr_base[length(chr_n+1)]), ylim=c(y_min,y_max))
    axis(side = 2, at=y_at, labels=y_at)
    mtext(expression(-log[10](P-value)), side = 2, line = 2)
    if(xlab_all){
        xlas=2
    }else{xlas=1}
    axis(side=1, line = -1, at=x_at, labels=xlabel, tick=F, las=xlas)
    mtext("Chromosome", side = 1, line = 1)
    mtext(main, side = 3, line = 1)
    # Plot genome-wide significant line
    if(gws_line){
        abline(h=gws_line, col="red")
    }
    # Plot suggestive significant line
    if(suggest_line){
        abline(h=suggest_line, col="blue", lty=2)
    }
}
