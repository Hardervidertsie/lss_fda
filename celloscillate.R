

require(reshape2)
require(ggplot2)
require(fda)
require(data.table)

celloscillate <- function(x, Frame = NULL, keepFrameTime = FALSE, value = NULL, locationID = NULL, timeInterval = NULL, aboveRegr = 5, suppresMaxFrames = 5,
                          slopeTH = NULL, basisType = NULL, slopeDomain = NULL,slopeDomain2 = NULL, timeToTH = NULL, maxValueTime = TRUE, plotFit = FALSE, plotDomains = FALSE, 
                           nbasis = 19, lambda = NULL, THpeaks = 0.1, f = 0.2, signal1 = NULL, signal1TH = NULL, signal2 = NULL,
                          xCoord = NULL, yCoord = NULL, pix.x = NULL, pix.y = NULL, ...) {
  
 
  
  if(any(is.na(c(value, Frame)))){
    stop("NA in value or Frame not implemented")
  }
  if(!is.logical(keepFrameTime) ){
    stop('keepFrameTime is not logical')
  }
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
  if(!is.null(suppresMaxFrames)){
    if(!is.numeric(suppresMaxFrames)){
      stop('suppresMaxFrames not numeric')
    }
  }
  if( !keepFrameTime & is.null(timeInterval) ){
    stop("provide timeInterval when KeepFrameTime is FALSE")
  }
  
  if(!is.character(timeInterval) & !is.null(timeInterval)) {
    stop("provide timeInterval as \"HH:MM:SS\"")
  }
  if(!is.null(timeInterval)){
    if(any(
      class(
        strptime(timeInterval, format = "%H:%M:%S")
        ) %in% "POSIXlt"
    ) && nchar(timeInterval) != 8 ){
      stop("provide timeInterval as \"HH:MM:SS\"")
    }
  }
  if(!is.logical(maxValueTime)){
    stop("maxValueTime not logical")
  }
  
  if( is.null(basisType)){
    stop('basisType must be fourier or spline')
  }
  if( !basisType %in% c('fourier', 'spline')) {
    stop('basisType must be fourier or spline')
  } 
  
  if(!keepFrameTime){
  if(any(as.integer(Frame) != Frame)){
  stop('Frame numbers must be integers or conversion to integer should be possible')
  }
  }
  
  if(!is.null(slopeDomain) ) {
   if(length(slopeDomain) != 2) {
     stop('slopeDomain must be vector of length 2 with start and end time')
   } 
  }
  if(!is.null(slopeDomain) ) { 
     if(!any(
       class(
         strptime(slopeDomain, format = "%H:%M:%S")
       ) %in% "POSIXlt"
     ) || any(nchar(slopeDomain) != 8) ){
       stop("provide slopeDomain as \"HH:MM:SS\"")
     } 
  }
  
  if(!is.null(timeToTH)) {
    if(length(timeToTH) != 2) {
      stop('timeToTH must be vector of length 2')
    }
    if(!is.numeric(as.numeric(timeToTH[1])) || as.numeric(timeToTH[1]) != timeToTH[1] ){
      stop('First entry in timeToTH must be numeric')
    }
    if(!timeToTH[2] %in% c('relative', 'absoluut')){
      stop('Second entry in timeToTH must be one of \"relative\", \"absoluut\"')
    }
    if(timeToTH[2] == 'relative' && (timeToTH[1] >= 1 | timeToTH[1] <= 0 )){
      stop('relative timeToTH requires input between 0 and 1 ')
    }
    }     
  
  if(!is.null(signal1)){
    if(length(signal1) != length(value)){
      stop('signal1 should be equal length to value')
    }
  if(!is.numeric(signal1)){
    stop('signal1 should be numeric')
  }    
  }
  if(!is.null(signal2)){
    if(length(signal2) != length(value)){
      stop('signal2 should be equal length to value')
    }
    if(!is.numeric(signal2)){
      stop('signal2 should be numeric')
    }
  }
  if(!is.null( xCoord )& !is.null( yCoord) ) {
   if(!is.numeric(xCoord) | !is.numeric(yCoord)){
     stop('x and yCoord must be numeric')
   }
  
    if(length(xCoord) != length(yCoord) | length(xCoord) != length(Frame)){
      stop('Please provide equal x, yCoord and Frame lengths')
    }  
  }
  
  if(!is.null(c(pix.x, pix.y))){
    if(!is.numeric(c(pix.x, pix.y))){
      stop("pix.x and pix.y not numeric")
    }
  }
   
  if(!is.null(signal1TH)){
    if(!is.numeric(signal1TH)){
      stop("signal1TH not numeric")
    }
  } 
  if(is.null(signal1) & !is.null(signal1TH)){
    stop("signal1TH provided but not signal1")
  }
  
  
  
  
  convertToHours <- function(x) {
    round(as.integer(strftime(strptime(x, format = "%H:%M:%S"), "%H")) + 
            1/60 * as.integer(strftime(strptime(x, 
                                                format = "%H:%M:%S"), "%M")) +
            1/3600 * as.integer(strftime(strptime(x, 
                                                  format = "%H:%M:%S"), "%S"))
          , digit = 4 )
    
  }
  
  
  
  if(!keepFrameTime){
      Frame <- as.integer(Frame)
      timeBetweenFrames <- convertToHours(timeInterval)
      Frame <- Frame * timeBetweenFrames    
    
  }
  
  
  if(!is.null(slopeDomain)){
  slopeDomain <- sapply(slopeDomain, convertToHours)
  }
  
  
  # estimation of not-provided parameters
  if(is.null(lambda)){
    lambda <- max(Frame)/50000
  }
  output = alist()
  
  # fit tracks
  
  rangval <- range(Frame)
  
  if(basisType == 'fourier'){
  period = max(rangval)
  basisfd = create.fourier.basis(rangval,nbasis, period = 2 * max(Frame) )
  } else {
    basisfd = create.bspline.basis(rangval, nbasis )
  }
  datafdPar <- fdPar(basisfd, lambda = lambda) # smoothing
  datalist <- smooth.basis(Frame, value, datafdPar) ## data
  predictdata <- predict(newdata = Frame, datalist)
  
  inputArgs <- seq(from = range(Frame)[1], to = range(Frame)[2], length.out = 500)
  
  fitData <- eval.fd( evalarg = inputArgs, fdobj = datalist$fd, Lfdobj = 0)
  firstDeriv <- eval.fd( evalarg = inputArgs, fdobj = datalist$fd, Lfdobj = 1)
  secondDeriv <- eval.fd( evalarg = inputArgs, fdobj = datalist$fd, Lfdobj = 2)
  
  residSD <-  sum((predictdata - value)^2) / (nbasis+1)
  
  # plot fits
  
  if(plotFit){

    
    locationID <- unique(locationID)
   # pdf( file = paste0('celloscillate/plots/fit_', locationID, '.pdf') , height = 6, width = 6 )
    
      plot(Frame, value, main = locationID, cex.main = 0.5,  ...)
      lines(datalist$fd)
      text( x = range(Frame)[2] - round(range(Frame)[2]/7, digits = 1)  , 
            y = range(value)[2] - round(range(value)[2]/10, digits = 1), labels = 
              paste('resiSD:', residSD), cex = 0.5, xlab = "time[h]" )
    
   # dev.off()
      
  }
  
   # determine domains  
  
  inputArgs <- seq(from = range(Frame)[1], to = range(Frame)[2], length.out = 1000)
  
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
  
  myData[, lowSlope := FALSE]
  myData[ abs(firstDeriv) > slopeTH, lowSlope := TRUE]
  
  
  # define the peak regions around local maxima; ( not sure this is the way to go: too sensitive to setting this par )
  
  myData[ fitData > lowessfit, peakArea := TRUE]

  # gather measurements based on all peaks
  output$timedomain <- range(Frame)
  output$starttime <- min(Frame)
  output$timelength <- max(Frame) - min(Frame)
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
  # max slope, return slope of max(abs(slope)
  if( max(myData$firstDeriv) >= max(abs(myData$firstDeriv)) ){
  output$maxslope <- max(myData$firstDeriv)
  } else {
    output$maxslope <- min(myData$firstDeriv)
  }
  
  
  output$minslope <- min(myData$firstDeriv)
  output$maxvalue <- max(value)
  indmax <- min(which(value == max(value)))
  output$maxValueTime <- Frame[indmax]
  output$minvalue <- min(value)
  output$maxfitvalue <- max(myData$fitData)
  output$minfitvalue <- min(myData$fitData)
  output$valueStart <- value[1]
  output$valueEnd <- value[length(value)]
  
  output$residSD <- residSD
  output$df <- datalist$df
  
  if(!is.null(slopeDomain)){
    output$slopeDomain <- myData[ inputArgs > slopeDomain[1]  & inputArgs < slopeDomain[2] , mean(firstDeriv)]
  }
  if(!is.null(slopeDomain2)){
    output$slopeDomain2 <- myData[ inputArgs > slopeDomain2[1]  & inputArgs < slopeDomain2[2] , mean(firstDeriv)]
  }
  
  if(!is.null(timeToTH)){
    if(timeToTH[2] == 'relative' ){
      relTH <- as.numeric(timeToTH[1]) * max(myData[,fitData])
      ind <- min(which(myData$fitData >= relTH ))
      output$timeToTH <- myData$inputArgs[ind]
    } else {
      if(any(myData$fitData >= as.numeric(timeToTH[1]))){
        ind <- min(which(abs(myData$fitData) >= as.numeric(timeToTH[1]) ))
        output$timeToTH <- myData$inputArgs[ind]
      } else{
        output$timeToTH <- NA
      }
      
    }
  }
  
  
  if(!is.null(signal1) & !is.null(signal1TH) ){
      if(any(signal1 > as.numeric(signal1TH))){
        ind <- min(which(signal1 > as.numeric(signal1TH)))
        output$timesignal1TH <- Frame[ind] # when is threshold reached.
      }
    } else{
      output$timesignal1TH <- NA
    }
  
  
  
  if(!is.null(signal1)){
  output$meansignal1 <- mean(signal1, na.rm = TRUE)
  } 
  if(!is.null(signal2)){
  output$meansignal2 <- mean(signal2, na.rm = TRUE)
  }

  # split data in seperate peaks to determine peak specific data
  # slopes of both sides of each peak & peak width
  
  # strategy: fill the gaps of the lowSlope around the local maxima (because of higher certainty of peaks as compared to minima )
  # then use the remaining gaps in lowSlope together with lowess line (peakArea) to identify the peaks (gaps are minima)
  # note that peakArea ensures regions being a peak
  # regions to fill gaps are determined in enfollowing while blocks
 
  indmax <- which(myData[ ,localmax] )
  indmin <- which(myData[ ,localmin] )
  if(length(indmax) !=0 & length(indmin) != 0){
  fp <- suppresMaxFrames 
  while(all((indmax - fp) > 0) && all((indmax + fp) < 1000) ){ # grow fill-vector untill out of boundary
    fp <- fp + 1
    }
  
  indfill = vector()
  
  for( i in seq_along(indmax)){
    indfill <- c(indfill, (indmax[i] - fp) : (indmax[i] + fp) )
  }
  
  while(any(indmin %in% indfill)) { # no filling of minima allowed
    fp <- fp - 1
    if(fp == suppresMaxFrames) break
    indfill = vector()
    for( i in seq_along(indmax)){
      
      indfill <- c(indfill, (indmax[i] - fp) : (indmax[i] + fp) )
    
      }
  }
  
  while(sum(indmax %in% indfill) > length(indmax)) { # only one indmax filled per indmax entry
    fp <- fp - 1
    if(fp == suppresMaxFrames) break
    indfill = vector()
    for( i in seq_along(indmin)){
      
      indfill <- c(indfill, (indmax[i] - fp) : (indmax[i] + fp) )
      
    }
  }
  myData[ indfill  , lowSlope := TRUE ] # actual filling of peak lowSlope holes
  output$fp <- fp # handy to check for out of ordinary values
  } # end if statement  
    
  
  #myData[ indmin, ]
  myData[ lowSlope == TRUE & shift(lowSlope, type = 'lead', n = 1, fill = TRUE) == FALSE   , truefalse := 1] 
  myData[is.na(truefalse), truefalse := 0]
  myData[ , peakids := cumsum(truefalse)]
  
  myDataPeaks <- split(myData, myData$peakids)
  # do the peakid regions contain peakArea > length 5? do they contain a localmax?
  indKeep1 <- sapply(myDataPeaks, function(x) sum(x$peakArea, na.rm = TRUE)) >= aboveRegr # # of  points above lowess regression line)
  indKeep2 <- sapply(myDataPeaks, function(x) sum(x$localmax, na.rm = TRUE)) > 0
  indKeep <- indKeep1 & indKeep2
  myDataPeaks <- myDataPeaks[indKeep]
  
  output$locationID <- unique(locationID)
  
  output$nsplitpeaks <- length(myDataPeaks)
  
  if(output$nsplitpeaks != 0 ){
  
  output$perpeak_max <- sapply(myDataPeaks, function(x) max(x$fitData))
  output$perpeak_maxslope <- sapply(myDataPeaks, function(x) max(x$firstDeriv))
  output$perpeak_minslope <- sapply(myDataPeaks, function(x) min(x$firstDeriv))
  output$perpeak_peakwidth <- sapply(myDataPeaks, function(x) sum(x$lowSlope))
  output$perpeak_peakwidth <- output$perpeak_peakwidth * output$timedomain[2]/ 1000
  output$perpeak_skew <- sapply(myDataPeaks, function(x) {
    left_dT <- x$inputArgs[x$localmax] - min(x$inputArgs)
    right_dT <- max(x$inputArgs) -  x$inputArgs[x$localmax]
    right_dT / left_dT # equal = 1, right tailed( left slope steeper than right slope ) then > 1, left tailed then < 1
  }
  )
  
  
  
  
  if( length(output$perpeak_max) != 0 ) {
  output$decreasing <- ifelse( all( unlist(output$perpeak_max) == rev(sort(unlist(output$perpeak_max)) )), 'yes', 'no')
  } else{
    output$decreasing <- ifelse( mean(myData$firstDeriv) < 0 , 'yes', 'no')
}
  
  } else{ # if no peaks..
  
    output$perpeak_max <- NA
    output$perpeak_maxslope <- NA
    output$perpeak_minslope <- NA
    
    output$perpeak_peakwidth <- NA
    output$perpeak_skew <- NA
    output$decreasing <- NA
    
  }

  
  
  
  if(plotDomains){
    
    if(!is.null(xCoord)){
   # pdf( file = paste0('celloscillate/plots/domains_', locationID, '.pdf') , height = 12, width = 6 )
    } else{
    #  pdf( file = paste0('celloscillate/plots/domains_', locationID, '.pdf') , height = 6, width = 6 )
    }
    if(!is.null(xCoord)){
    par(mfrow=c(2,1))
    }
    if(!is.null(xCoord) & !is.null(signal1)){
      par(mfrow = c(3,1))
    }
    plot(x = Frame, value, xlab = 'time', ylim = c(0, max(value)), main = locationID)
      lines(x = inputArgs, fitData)
      points( myData[localmax == TRUE, inputArgs],   myData[localmax == TRUE, fitData], cex = 5, col = 'blue')
      points( myData[localmin == TRUE, inputArgs],   myData[localmin == TRUE, fitData], cex = 5, col = 'red')
      points( myData[lowSlope== TRUE ,inputArgs], rep(0.06, sum(myData$lowSlope,na.rm = TRUE)) , col ="orange", cex = 0.5)
      points( myData[peakArea== TRUE ,inputArgs], myData[peakArea == TRUE, rep(0.05, sum(myData$peakArea, na.rm=TRUE)) ], cex= 0.5)
      
      
      lines(myData$inputArgs, myData$lowessfit, lty = 2)
      
      if(exists('truefalse', where = myData)){
      abline( v = myData[truefalse != 0, inputArgs], lty = 4)
      }
      if(!is.null(signal1)){
      text( x = range(Frame)[2] - round(range(Frame)[2]/5, digits = 1)  , 
            y = range(value)[2] - round(range(value)[2]/8, digits = 2), labels = 
              paste('mean_signal1:', round(output$meansignal1, digits = 3) ) )
    }
      if(!is.null(signal2)){
      text( x = range(Frame)[2] - round(range(Frame)[2]/5, digits = 1)  , 
            y = range(value)[2] - round(range(value)[2]/4, digits = 2), labels = 
              paste('mean_signal2:', round(output$meansignal2, digits = 3) ) )
      }
      
      if(!is.null(xCoord)){
      plot( x = xCoord, y = pix.y - yCoord, xlim = c(0, pix.x), ylim = c(0, pix.y), cex = 0.2, pch = '.')
      text(x = xCoord + 5, y = pix.y - yCoord + 5, round(Frame, digits = 0 ), cex = 0.2)  
      
      }
      
      if(!is.null(signal1)){
      plot(Frame, signal1, main = "signal1", col = "red", ylim = c(0, c(max(signal1, signal2))))
        legend( legend = "signal1", x = Frame[1], y = max(signal1) - round(max(signal2)/10, digits= 0), pch= 1 , col = "red")
      }
      
      if(!is.null(signal2)){
      points(Frame, signal2, main = "signal2", col = "blue")
        legend( legend = "signal2", x = Frame[1], y = max(signal2) - round(max(signal2)/12, digits= 0), pch= 1 , col = "blue")
      }
    
  #  dev.off()
  }

  return(output)
  
}




