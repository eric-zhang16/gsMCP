## Design 1 ##
rm(list=ls())
library(gsMCP)
G <- matrix(0, nr=4, nc = 4)                                  # G is transition matrix
G[1,2] <- G[1,3] <- G[2,1] <- G[2,4] <- 1/2
G[3,2] <- G[4,1] <- 1
w.start <- c(0.5,0.5,0,0)                                        # w.start is the starting weights
t<- seq(1,2)                                                             # t are the IA/FA indicators
h <- seq(1,4)                                                           # h are endpoint indicators
p <- matrix(NA, nr=length(t), nc = length(h))   # p are the p-value matrix
p[1,]<-c(0.0339,0.1411,0.1533,0.0313)
p[2,]<-c(0.0108,0.0183,0.0007,0.00002)
timing=0.7                                                            # timing is the information fraction for IA




# sf is the spending function names in ('OF','Pocock','WT','HSD')
# sfpar is the spending function parameter, no need to specify for sf=OF/Pocock
# debug=1 to print some processing details 
gMCPgSD_IW(G,w.start,t,h,p,alpha=0.025,timing,sf='OF',debug=0)      
gMCPgSD_BF(G,w.start,t,h,p,alpha=0.025,timing,sf='HSD',sfpar=-2,debug=0)



           
