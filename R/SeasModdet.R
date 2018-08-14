
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


SeasModdet <- function(pHSinit=0.8, Kx = 100, betax=0.02, wx=0.8, hx=1, mx=1, ax=1, rx=0.1,zx=0.9,gx=4, cx=0.9, phix=0, maY=100, miY=0, thetax=0.2, nseasons=10, Ex=0) {

outm <- as.data.frame(matrix(data=-999, nrow=(nseasons+1), ncol=10, dimnames = list(1:(nseasons+1),c('season','HP', 'DP','HS','DS','pHS','pDS','Yld','YL', 'DPbr'))))

outm[1,] <- NA# row one gives initial conditions
outm$season <- 0:nseasons

# Initial value for state variables
outm$pHS[1] <- pHSinit # initial proportion healthy seed (for nseasons=0)
outm$pDS[1] <- 1 - pHSinit

# seasons 2 and higher

for(si in 2:(nseasons+1)) {
# Equation A1
   outm$HP[si] <- min(max(0,Kx * outm$pHS[si-1] - betax * wx * hx * mx * ((Kx * outm$pDS[si-1]) * (Kx * outm$pHS[si-1])+ Ex*(Kx * outm$pHS[si-1]))), 100) 
# max prevents <0 plants & min prevents >K plants
# Equation A2
   outm$DPbr[si] <- min(max(0, Kx * outm$pDS[si-1] + betax * wx * hx * mx * ((Kx * outm$pDS[si-1]) * (Kx * outm$pHS[si-1]) + Ex*(Kx * outm$pHS[si-1]))), 100)
# Equation A3
   outm$DP[si] <- ax*outm$DPbr[si]
# Equation B1
   outm$HS[si] <- max(0,gx * (outm$HP[si] + rx * outm$DP[si]))
# Equation B2
   outm$DS[si] <- max(0,zx * cx * gx * (1 - rx) * outm$DP[si])
# Equation C1 
   outm$Yld[si]<- ((outm$HP[si]+outm$DP[si])/(outm$HP[si]+outm$DPbr[si]))*(miY+(maY-miY)*((1-(outm$DP[si]/(outm$DP[si]+outm$HP[si])))/((1-thetax)+ thetax*(1-(outm$DP[si]/(outm$DP[si]+outm$HP[si]))))^2))
#Equation C2
   outm$YL[si]<- maY-outm$Yld[si]
# Equation D1
   outm$pHS[si] <- phix + (1-phix) * outm$HS[si]/(outm$HS[si] + outm$DS[si])
# Equation D2
   outm$pDS[si] <- 1 - outm$pHS[si]
 }

outm[(nseasons+1),]
}
