
#' Create a matrix with output for multiple parameter combinations
#'
#' This function creates a matrix with output for multiple parameter combinations across nseasons
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

# Weather (wx), vector management (mx), positive selection (zx) and roguing (zx) are stochastic
# Each have a mean and associated standard deviation 

# set.seed(1234)  


mpc <- function(vpHSinit=c(0.2,0.5,0.8), vKx = 100, vbetax=c(0.02,0.04), vwxtnormm=c(0.3,0.7), vwxtnormsd=c(0.3, 0.1), vhx=1, vmxtnormm=1, vmxtnormsd=0.1, vaxtnormm=c(1,0.5), vaxtnormsd= c(0.3, 0.1), vrx=c(0.1,0.3), vzxtnormm=c(1,0.5), vzxtnormsd= c(0.3, 0.1), vgx=4, vcx=0.9, vphix=c(0,0.5), vnseasons=10, vnsim=1, vHPcut=0.5, vpHScut=0.5,vmaY=100, vmiY=0,  vthetax=c(-0.5,0,0.5), vEx=0.02){

ncomb <- length(vpHSinit) * length(vKx) * length(vbetax) * length(vwxtnormm) * length(vwxtnormsd)*length(vhx) * length(vmxtnormm)*length(vmxtnormsd) * length(vaxtnormm) * length(vaxtnormsd) * length(vrx) * length(vzxtnormm) * length(vzxtnormsd) * length(vgx) * length(vcx) * length(vphix) * length(vnseasons) * length(vnsim)*length(vHPcut)*length(vpHScut)*length(vmaY)*length(vmiY)*length(vthetax)*length(vEx)

basn <- c('fHP', 'fDP', 'fHS', 'fDS', 'fpHS', 'fpDS', 'HPtrans', 'pHStrans', 'HPpseas', 'pHSpseas', 'fYld', 'fYL')

outmpc <- as.data.frame(matrix(data=-999, nrow=ncomb, ncol=84, dimnames = list(1:ncomb,c('pHSinit','Kx','betax','wxtnormm','wxtnormsd','hx','mxtnormm','mxtnormsd','axtnormm','axtnormsd','rx','zxtnormm','zxtnormsd','gx','cx','phix','nseasons','nsim', 'HPcut', 'pHScut', 'maY','miY','thetax','Ex', paste(basn,'mean',sep=''),paste(basn,'median',sep=''),paste(basn,'var',sep=''),paste(basn,'q0.05',sep=''),paste(basn,'q0.95',sep='')))))

icomb <- 1 # indicates number of current parameter combination, and corresponding row of outmpc

for(i1 in 1:length(vpHSinit)){ tpHSinit <- vpHSinit[i1]
for(i2 in 1:length(vKx)){tKx <- vKx[i2]
for(i3 in 1:length(vbetax)){tbetax <- vbetax[i3]
for(i4 in 1:length(vwxtnormm)){twxtnormm <- vwxtnormm[i4]
for(i5 in 1:length(vwxtnormsd)){twxtnormsd<- vwxtnormsd[i5]
for(i6 in 1:length(vhx)){thx <- vhx[i6]
for(i7 in 1:length(vmxtnormm)){tmxtnormm <- vmxtnormm[i7]
for(i8 in 1:length(vmxtnormsd)){tmxtnormsd <- vmxtnormsd[i8]
for(i9 in 1:length(vaxtnormm)){taxtnormm <-vaxtnormm[i9]
for(i10 in 1:length(vaxtnormsd)){taxtnormsd <- vaxtnormsd[i10]
for(i11 in 1:length(vrx)){trx <- vrx[i11]
for(i12 in 1:length(vzxtnormm)){tzxtnormm <-vzxtnormm[i12]
for(i13 in 1:length(vzxtnormsd)){tzxtnormsd <- vzxtnormsd[i13]
for(i14 in 1:length(vgx)){tgx <- vgx[i14]
for(i15 in 1:length(vcx)){tcx <- vcx[i15]
for(i16 in 1:length(vphix)){tphix <- vphix[i16]
for(i17 in 1:length(vnseasons)){tnseasons <- vnseasons[i17]
for(i18 in 1:length(vnsim)){tnsim <- vnsim[i18]
for(i19 in 1:length(vHPcut)){tHPcut<- vHPcut[i19]
for(i20 in 1:length(vpHScut)){tpHScut<- vpHScut[i20]
for(i21 in 1:length(vmaY)){tmaY<- vmaY[i21]
for(i22 in 1:length(vmiY)){tmiY<- vmiY[i22]
for(i23 in 1:length(vthetax)){tthetax<- vthetax[i23]
for(i24 in 1:length(vEx)){tEx<-vEx[i24]
   

temp <- Multstoch(pHSinit2=tpHSinit, Kx2 = tKx, betax2=tbetax, wxtnormm2=twxtnormm, wxtnormsd2=twxtnormsd, hx2=thx, mxtnormm2=tmxtnormm, mxtnormsd2=tmxtnormsd, axtnormm2=taxtnormm, axtnormsd2=taxtnormsd, rx2=trx, zxtnormm2=tzxtnormm, zxtnormsd2=tzxtnormsd, gx2=tgx, cx2=tcx, phix2=tphix, nseasons2=tnseasons, nsim2=tnsim, HPcut2=tHPcut, pHScut2=tpHScut,maY2=tmaY,miY2=tmiY, thetax2=tthetax, Ex2=tEx)$outfsum 

outmpc[icomb,1:24] <- c(tpHSinit,tKx,tbetax,twxtnormm,twxtnormsd,thx,tmxtnormm, tmxtnormsd,taxtnormm,taxtnormsd,trx,tzxtnormm,tzxtnormsd,tgx,tcx, tphix,tnseasons,tnsim,tHPcut,tpHScut,tmaY,tmiY,tthetax, tEx)
outmpc[icomb,25:36] <- temp[1,]
outmpc[icomb,37:48] <- temp[2,]
outmpc[icomb,49:60] <- temp[3,]
outmpc[icomb,61:72] <- temp[4,]
outmpc[icomb,73:84] <- temp[5,]

icomb <- icomb + 1
}}}}}}}}}}}}}}}}}}}}}}}}
outmpc
}
