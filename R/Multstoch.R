
#' Create a matrix with output from multiple simulations for a single parameter combination
#'
#' This function creates a matrix for output from multiple simulations for a single parameter combination across nseasons
#'
#' Updated 2018-08-21

#' @param Kx total number of plants
#' @param betax maximum seasonal transmission rate
#' @param wx environmental effect on transmission rate
#' @param hx host effect on transmission rate
#' @param mx vector management effect on transmission rate
#' @param ax roguing effect in terms of decreased DP
#' @param maY maximum attainable yield, end of season, in the absence of disease
#' @param miY minimum yield when all plants are diseased (useable yield despite disease)
#' @param thetax rate of decline of Yld with increasing disease incidence
#' @param Ex amount of external inoculum around field
#' @param rx reversion rate
#' @param zx proportional selection against diseased plants
#' @param gx seed production rate in healthy plants
#' @param cx proportional seed production rate in diseased plants
#' @param phix proportion clean seed purchased
#' @param nseasons number of seasons
#' @keywords seed health
#' @export 
#' @examples
#' Multstoch()
#' Multstoch()

# to do - GENERAL TESTING 
# to do - check whether parameter list is correct

# Columns of output matrix
# col 1 - timestep (initial time step is season 0)
# col 2 - HP healthy plant number
# col 3 - DP diseased plant number (after roguing)
# col 4 - HS healthy seed number
# col 5 - DS diseased seed number
# col 6 - pHS proportion healthy seed
# col 7 - pDS proportion diseased seed
# col 8 - Yld end of season yield
# col 9 - YL end of season yield loss
# col 10 - DPbr (diseased plants before roguing) 

# Step 1B. Create matrix for output from one simulation (stochastic model)

# Weather (wx), vector management (mx), positive selection (zx) and roguing (zx) are stochastic
# Each have a mean and associated standard deviation 

# set.seed(1234)  


Multstoch <- function(pHSinit2=0.8, Kx2 = 100, betax2=0.02, wxtnormm2=0.8, wxtnormsd2=0.3, hx2=1, mxtnormm2=1, mxtnormsd2=0.1, axtnormm2=1, axtnormsd2=0.1, rx2=0.1, zxtnormm2=1, zxtnormsd2= 0.1, gx2=4, cx2=0.9, phix2=0, nseasons2=10, nsim2=10, HPcut2=0.5,pHScut2=0.5,maY2=100, miY2=0, thetax2=0.2, Ex2=0) {

# nsim - number of simulations

outmf <- as.data.frame(matrix(data=-999, nrow=nsim2, ncol=12, dimnames = list(1:nsim2,c('fHP', 'fDP', 'fHS', 'fDS', 'fpHS', 'fpDS', 'HPtrans', 'pHStrans', 'HPpseas', 'pHSpseas', 'fYld', 'fYL'))))

for(si in 1:nsim2) {
temp <- SeasModstoch(pHSinit=pHSinit2, Kx = Kx2, betax=betax2, wxtnormm=wxtnormm2, wxtnormsd=wxtnormsd2, hx=hx2, mxtnormm=mxtnormm2, mxtnormsd=mxtnormsd2, axtnormm=axtnormm2, axtnormsd=axtnormsd2, rx=rx2,zxtnormm=zxtnormm2, zxtnormsd= zxtnormsd2, gx=gx2, cx=cx2, phix=phix2, nseasons=nseasons2, HPcut=HPcut2, pHScut=pHScut2, maY=maY2, miY=miY2, thetax=thetax2, Ex=Ex2)$outfin   
   outmf$fHP[si] <- temp$HP
   outmf$fDP[si] <- temp$DP
   outmf$fHS[si] <- temp$HS
   outmf$fDS[si] <- temp$DS
   outmf$fpHS[si] <- temp$pHS
   outmf$fpDS[si] <- temp$pDS
   outmf$HPtrans[si] <- temp$HPtrans
   outmf$pHStrans[si] <- temp$pHStrans
   outmf$HPpseas[si] <- temp$HPpseas
   outmf$pHSpseas[si] <- temp$pHSpseas
   outmf$fYld[si]<-temp$Yld
   outmf$fYL[si]<-temp$YL
 }

quantile0.05 <- function(x){quantile(x,probs=0.05,na.rm=T)}
quantile0.95 <- function(x){quantile(x,probs=0.95,na.rm=T)}

outfsum <- as.data.frame(matrix(data=-999, nrow=5, ncol=12, dimnames = list(c('mean','median','var','q0.05','q0.95'), c('fHP', 'fDP', 'fHS', 'fDS', 'fpHS', 'fpDS',  'HPtrans', 'pHStrans', 'HPpseas', 'pHSpseas', 'fYld', 'fYL'))))

# the first row gives the mean for each of the responses, the second the median, etc.
outfsum[1,] <- apply(outmf,MARGIN=2,FUN=mean)
outfsum[2,] <- apply(outmf,MARGIN=2,FUN=median)
outfsum[3,] <- apply(outmf,MARGIN=2,FUN=var)
outfsum[4,] <- apply(outmf,MARGIN=2,FUN=quantile0.05)
outfsum[5,] <- apply(outmf,MARGIN=2,FUN=quantile0.95)

list(outmf=outmf,outfsum=outfsum)
}
