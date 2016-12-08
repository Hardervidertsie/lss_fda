
ratioDampener <- function(nom = NULL, denom = NULL, dampingFactor = c(0, 0.1), n = 20, ...) { 
  # dampingFactor should be start and end values of length 2, between 0 and 1
  # n number of factors between dampingFactor range to evaluate function at
  require(fda)
  require(data.table)
  
  if(any(is.null(c(nom, denom, dampingFactor)))){
    stop('Provide nom, denom and dampingFactor[0,1]')
  }
  
  if(!length(nom) == length(denom)){
    stop('unequal vector length')
  }
  
  if(any(!is.numeric(c(nom, denom, dampingFactor)))){
    stop('input must be numeric')
  }
  
  if(length(dampingFactor) != 2) {
    stop('dampingFactor must be bector of length 2' )
  }
  
  if( max(dampingFactor) >= 1 || min(dampingFactor) < 0) {
    stop('dampingFactor should live in domain [0,1]')
  }
  
  OutputList = alist()
  Damplist = alist()
  Evalrange <- seq( from = dampingFactor[1], to = dampingFactor[2], length.out = n)
  
  for( i in seq_along(Evalrange)){
    
    Damp <- Evalrange[i] * mean(denom, na.rm = TRUE)
    OutputList[[i]] <- quantile( nom / (denom + Damp), na.rm = TRUE)
    
  }
  Output <- do.call('cbind',OutputList)
  colnames(Output) <- round(Evalrange, digits = 4)
  boxplot( x = Output, main = 'distribution of ratio', ...)
  
}