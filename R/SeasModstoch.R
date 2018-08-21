
#' Generate an output matrix for one simulation
#'
#' This function simulates one parameter combination across nseasons once 
#'
#' Updated 2018-08-14

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
#' sim1()
#' sim1()

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

# This function simulates nseasons for one parameter combination once 

SeasModstoch <- function(pHSinit=0.8, Kx = 100, betax=0.02, wxtnormm=0.8, wxtnormsd=0.3, hx=1, mxtnormm=1, mxtnormsd=0.1, axtnormm=1, axtnormsd=0.1, rx=0.1, zxtnormm=1, zxtnormsd= 0.1, gx=4, cx=0.9, phix=0, nseasons=10, HPcut=0.5, pHScut=0.5, maY=100, miY=0, thetax=0.2,  Ex=0) {

outm <- as.data.frame(matrix(data=-999, nrow=(nseasons+1), ncol=14, dimnames = list(1:(nseasons+1),c('season','HP', 'DP','HS','DS','pHS','pDS','mx','zx','ax','wx', 'Yld','YL','DPbr'))))

outm[1,] <- NA # row one gives initial conditions
outm$season <- 0:nseasons

# Initial value for state variables
outm$pHS[1] <- pHSinit # initial proportion healthy seed (for nseasons=0)
outm$pDS[1] <- 1 - pHSinit

# seasons 2 and higher

for(si in 2:(nseasons+1)) {
#Generating stochastic  varaibles 
   outm$wx[si] <-ifelse(wxtnormsd==0, wxtnormm, (altrtruncnorm(1,a=0,b=1,mean=wxtnormm,sd=wxtnormsd))) 
  outm$mx[si] <- ifelse(mxtnormsd==0, mxtnormm, (altrtruncnorm(1,a=0,b=1,mean=mxtnormm,sd=mxtnormsd))) 
   outm$ax[si] <- ifelse(axtnormsd==0, axtnormm, (altrtruncnorm(1,a=0,b=1,mean=axtnormm,sd=axtnormsd))) 
   outm$zx[si] <-ifelse(zxtnormsd==0, zxtnormm, (altrtruncnorm(1,a=0,b=1,mean=zxtnormm,sd=zxtnormsd))) 
#Calculating rate of disease transmission
  tempnewinf<-  betax * outm$wx[si] * hx * outm$mx[si] * ((Kx * outm$pDS[si-1]) * (Kx * outm$pHS[si-1])+ Ex*(Kx * outm$pHS[si-1]))
#Equation A1
  outm$HP[si] <- min(max(0,Kx * outm$pHS[si-1] - tempnewinf),100)
# Equation A2
   outm$DPbr[si] <- min(max(0, Kx * outm$pDS[si-1] + tempnewinf),100)
# Equation A3
   # Equation A3
   outm$DP[si] <- outm$ax[si]*outm$DPbr[si]
# Equation B1
   outm$HS[si] <- max(0,gx * (outm$HP[si] + rx * outm$DP[si]))
# Equation B2
    outm$DS[si] <- max(0,outm$zx[si] * cx * gx * (1 - rx) * outm$DP[si])
#Equation C1
   outm$Yld[si]<- ((outm$HP[si]+outm$DP[si])/(outm$HP[si]+outm$DPbr[si]))*(miY+(maY-miY)*((1-(outm$DP[si]/(outm$DP[si]+outm$HP[si])))/((1-thetax)+ thetax*(1-(outm$DP[si]/(outm$DP[si]+outm$HP[si]))))^2))
# Equation C2
   outm$YL[si]<- maY-outm$Yld[si]
# Equation D1
   outm$pHS[si] <- phix + (1-phix) * outm$HS[si]/(outm$HS[si] + outm$DS[si])
# Equation D2
   outm$pDS[si] <- 1 - outm$pHS[si]
    }

outfin <- outm[(nseasons+1),]

# Season in which HP first transitions below HPcut*Kx
HPtrans <- outm$season[which(outm$HP < HPcut*Kx)][1]
# if HPtrans is NA i.e., HP never less than HPcut, set to max seasons tested
HPtrans[is.na(HPtrans)]<-nseasons
# Season in which pHS first transitions below pHScut
pHStrans <- outm$season[which(outm$pHS < pHScut)][1]
# if pHStrans is NA i.e., pHS never less than pHScut, set to max seasons tested
pHStrans[is.na(pHStrans)]<-nseasons
# Proportion seasons with HP below HPcut*Kx
# Healthy plants are calculated from season 1 onwards so +1 not needed in denominator
 HPpseas <- length(outm$season[which(outm$HP < HPcut*Kx)])/nseasons
# Proportion seasons with pHS below pHScut
# if +1 is not included in denominator, if in season 0 pHS<pHScut, the proportion is >1
 pHSpseas <- length(outm$season[2:nseasons][which(outm$pHS < pHScut)])/(nseasons+1)

# ptempa <- matrix(c(HPtrans,pHStrans,HPpseas,pHSpseas), dimnames=list(c('HPtrans','pHStrans','HPpseas','pHSpseas'),'1'))
# ptempb <- as.data.frame(t(ptempa))

ptemp <- as.data.frame(t(matrix(c(HPtrans,pHStrans,HPpseas,pHSpseas), dimnames=list(c('HPtrans','pHStrans','HPpseas','pHSpseas'),'1'))))

outfin <- cbind(outfin,ptemp)

list(outm=outm,outfin=outfin)
}
