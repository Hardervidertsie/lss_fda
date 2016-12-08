

celloscillate <- function(x, Frame = NULL, value = NULL, locationID = NULL, timeInterval = NULL, aboveRegr = 5,
                           plotFit = FALSE, plotDomains = FALSE, 
                           nbasis = 19, lambda = NULL, THpeaks = 0.1, f = 0.2) {
  
  if(!is.numeric(Frame)) {
    stop("input Frame must be numeric")
  }
  if(!is.numeric(value)){
    stop("input value must be numeric")
  }
  if(length(Frame) != length(value)){
    stop("unequal input vectors")
  }
  if(length(unique(Frame)) != length(value)){
    stop("Duplicate frame entries")
  }
  if(length(unique(locationID)) != 1 && !is.null(locationID)){
    stop("Duplicate locationID found for current track")
  }
  if(nbasis %% 2 == 0 ){
    stop("Provide uneven number of nbasis")
  }
  if(!THpeaks > 0) {
    stop('THpeaks not > 0')
  }
  if(!is.character(timeInterval)){
    stop("provide timeInterval as \"HH:MM:SS\"")
  }
  if(any(
    class(
      strptime(timeInterval, format = "%H:%M:%S")
      ) %in% "POSIXlt"
  ) && nchar(timeInterval) != 8 ){
    stop("provide timeInterval as \"HH:MM:SS\"")
  }
    
  if(any(as.integer(Frame) != Frame)){
  stop('Frame numbers must be integers or conversion to integer should be possible')
    }

  
      Frame <- as.integer(Frame)
    
  timeBetweenFrames <- round(as.integer(strftime(strptime(timeInterval, format = "%H:%M:%S"), "%H")) + 
                               1/60 * as.integer(strftime(strptime(timeInterval, 
                                                                   format = "%H:%M:%S"), "%M")) +
                               1/3600 * as.integer(strftime(strptime(timeInterval, 
                                                                     format = "%H:%M:%S"), "%S"))
                             , digit = 4 )
  
  
  Frame <- Frame * timeBetweenFrames
  
  
  
  
  # estimation of not-provided parameters
  if(is.null(lambda)){
    lambda <- max(Frame)/50000
  }
  output = alist()
  
  # fit tracks
  
  rangval <- range(Frame)
  period = max(rangval)
  basisfd = create.fourier.basis(rangval,nbasis, period = 2 * max(Frame) )
  datafdPar <- fdPar(basisfd, lambda = lambda) # smoothing
  datalist <- smooth.basis(Frame, value, datafdPar) ## data
  predictdata <- predict(newdata = Frame, datalist)
  
  inputArgs <- seq(from = range(Frame)[1], to = range(Frame)[2], length.out = 500)
  
  fitData <- eval.fd( evalarg = inputArgs, fdobj = datalist$fd, Lfdobj = 0)
  firstDeriv <- eval.fd( evalarg = inputArgs, fdobj = datalist$fd, Lfdobj = 1)
  secondDeriv <- eval.fd( evalarg = inputArgs, fdobj = datalist$fd, Lfdobj = 2)
  
  residSD <- round( sum((predictdata - value)^2) / datalist$df, digits = 3)
  
  # plot fits
  
  if(plotFit){

    if(!file.exists('celloscillate')){
      dir.create('celloscillate')
    }
    
    if(!file.exists('celloscillate/plots')){
      dir.create('celloscillate/plots')
    }
    
    locationID <- unique(locationID)
    pdf( file = paste0('celloscillate/plots/fit_', locationID, '.pdf') , height = 6, width = 6 )
    
      plot(Frame, value, main = "raw data + fourier fit")
      lines(datalist$fd)
      text( x = range(Frame)[2] - round(range(Frame)[2]/10, digits = 1)  , 
            y = range(value)[2] - round(range(value)[2]/10, digits = 1), labels = 
              paste('resiSD:', residSD) )
    
    dev.off()
      
  }
  
   # determine domains  
  
  inputArgs <- seq(from = range(Frame)[1], to = range(Frame)[2], length.out = 500)
  
  fitData <- eval.fd( evalarg = inputArgs, fdobj = datalist$fd, Lfdobj = 0)
  firstDeriv <- eval.fd( evalarg = inputArgs, fdobj = datalist$fd, Lfdobj = 1)
  secondDeriv <- eval.fd( evalarg = inputArgs, fdobj = datalist$fd, Lfdobj = 2)
  
  myData <- data.frame( inputArgs = inputArgs, fitData = fitData, firstDeriv = firstDeriv, secondDeriv = secondDeriv)
  myData <- as.data.table(myData)
  myData[, lowessfit := lowess(inputArgs, fitData, f = f)$y]
  
  myData[ firstDeriv < 0, sign:="neg"]
  myData[ firstDeriv >= 0, sign:="pos"]
  
  # check if switch sign (local max/min)
  myData [ , signlagged := c(NA, sign[-length(sign)]) ]
  myData[ sign != signlagged , localmaxmin := TRUE ]
  myData[ sign == signlagged , localmaxmin := FALSE ]
  
  myData[ , localmax := FALSE]
  myData[ , localmin := FALSE]
  
  myData[ localmaxmin & secondDeriv < -0.1, localmax := TRUE] # peaks should be clear, for minima this is tricky
  myData[ localmaxmin & secondDeriv > 0.005, localmin := TRUE] 
  
  # filter based on peak-height criterion: more than THpeaks-fraction difference between the lowess regression and peak of median raw values
  
  myData[ ((fitData - lowessfit)/ median(value)) < THpeaks  , localmax:= FALSE]
  myData[ ((lowessfit - fitData)/ median(value)) < THpeaks  , localmin:= FALSE]
  
  # calculate peak width: from peak to next minima, and then determine when slope is low enough to define end of peak
  slopeTH <- 0.02
  myData[, lowSlope := FALSE]
  myData[ abs(firstDeriv) > slopeTH, lowSlope := TRUE]
  
  
  # define the peak regions around local maxima; ( not sure this is the way to go: too sensitive to setting this par )
  
  myData[ fitData > lowessfit, peakArea := TRUE]

  
  
  # gather measurements based on all peaks
  
  output$npeaks = sum(myData$localmax)
  output$peakheights = myData$fitData[myData$localmax == TRUE] 
  output$peaktimes = myData$inputArgs[myData$localmax == TRUE] 
  if(length(output$peakheights) > 1){
  output$peak12damp = output$peakheights[2] / output$peakheights[1]
  } else {
    output$peak12damp <- NA
  }
  
  if(length(output$peakheights) > 2){
    output$peak23damp = output$peakheights[2] / output$peakheights[3] 
  } else {
    output$peak34damp <- NA
  }
  output$maxslope <- max(myData$firstDeriv)
  output$minslope <- min(myData$firstDeriv)
  output$maxvalue <- max(value)
  output$minvalue <- min(value)
  output$maxfitvalue <- max(myData$fitData)
  output$minfitvalue <- min(myData$fitData)
  
  output$residSD <- residSD
  output$df <- datalist$df
  
  # split data in seperate peaks to determine peak specific data
  # slopes of both sides of each peak & peak width
  
  # strategy: fill the gaps of the lowSlope around the local maxima
  # then use the remaining gaps together with lowess line (peakArea) to identify the peaks
  # note that peakArea ensures regions being a peak
  
  indmax <- which(myData[ ,localmax] )
  indmin <- which(myData[ ,localmin] )
  fp <- 10
  while(any((indmax - fp) < 0) || any((indmax - fp) >= 500) ){ # no negative indexes and no out of boundary
    fp <- fp - 1
    if(fp == 1) break
  }
  indfill = vector()
  for( i in seq_along(indmax)){
    indfill <- c(indfill, (indmax[i] - fp) : (indmax[i] + fp) )
  }
  fp_old <- fp
  while(any(indmin %in% indfill)) { # and no filling minima
    fp <- fp - 1
    if(fp == 1) break
    indfill = vector()
    for( i in seq_along(indm)){
      indfill <- c(indfill, (indmax[i] - fp) : (indmax[i] + fp) )
    }
  }
  
  
    
  myData[ indfill  , lowSlope := TRUE ] # actual filling of peak lowSlope holes
  
  myData[ lowSlope == TRUE & shift(lowSlope, type = 'lead', n = 1, fill = TRUE) == FALSE   , truefalse := 1] 
  myData[is.na(truefalse), truefalse := 0]
  myData[ , peakids := cumsum(truefalse)]
  
  myDataPeaks <- split(myData, myData$peakids)
  # do the peakid regions contain peakArea > length 5? do they contain a localmax?
  indKeep1 <- sapply(myDataPeaks, function(x) sum(x$peakArea, na.rm = TRUE)) >= aboveRegr # # of  points above lowess regression line)
  indKeep2 <- sapply(myDataPeaks, function(x) sum(x$localmax, na.rm = TRUE)) > 0
  indKeep <- indKeep1 & indKeep2
  myDataPeaks <- myDataPeaks[indKeep]
  
  output$fp <- fp # handy to check for out of ordinary values
  output$nsplitpeaks <- length(myDataPeaks)
  
  output$perpeak_max <- sapply(myDataPeaks, function(x) max(x$fitData))
  output$perpeak_maxslope <- sapply(myDataPeaks, function(x) max(x$firstDeriv))
  output$perpeak_minslope <- sapply(myDataPeaks, function(x) min(x$firstDeriv))
  output$perpeak_peakwidth <- sapply(myDataPeaks, function(x) sum(x$lowSlope))
  output$perpeak_skew <- sapply(myDataPeaks, function(x) {
    left_dT <- x$inputArgs[x$localmax] - min(x$inputArgs)
    right_dT <- max(x$inputArgs) -  x$inputArgs[x$localmax]
    right_dT / left_dT # equal = 1, right tailed( left slope steeper than right slope ) then > 1, left tailed then < 1
  }
  )
  output$locationID <- unique(locationID)
  
  output$decreasing <- ifelse( all( unlist(output$perpeak_max) == rev(sort(unlist(output$perpeak_max)) )), 'yes', 'no')

  
  if(plotDomains){
    pdf( file = paste0('celloscillate/plots/domains_', locationID, '.pdf') , height = 6, width = 6 )
    
      plot(x = Frame, value, xlab = 'time', ylim = (c(0, max(value))))
      lines(x = inputArgs, fitData)
      points( myData[localmax == TRUE, inputArgs],   myData[localmax == TRUE, fitData], cex = 5, col = 'blue')
      points( myData[localmin == TRUE, inputArgs],   myData[localmin == TRUE, fitData], cex = 5, col = 'red')
      points( myData[lowSlope== TRUE ,inputArgs], rep(0.06, sum(myData$lowSlope,na.rm = TRUE)) , col ="orange", cex = 0.5)
      points( myData[peakArea== TRUE ,inputArgs], myData[peakArea == TRUE, rep(0.05, sum(myData$peakArea, na.rm=TRUE)) ], cex= 0.5)
      lines(myData$inputArgs, myData$lowessfit, lty = 2)
      abline( v = myData[truefalse != 0, inputArgs], lty = 4)
    
    dev.off()
  }
  
  
  
  return(output)
  
}



