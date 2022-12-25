#library(gsMCP)
# H1: OS PDL1 positve
# H2: PFS PDL1 positve
# H3: OS Overall
# H4: PFS Overall

G <- matrix(0, nr=4, nc = 4)                                  # G is transition matrix
G[1,2] <- G[1,3] <- G[2,1] <- G[2,4] <- 1/2
G[3,2] <- G[4,1] <- 1
w.start <- c(0.5,0.5,0,0)                                        # w.start is the starting weights

t<- seq(1,4)                                    # t are the total IA/FA indicators, e.g. if OS will be tested at IA1-3 and PFS tested at IA2-4, then t=c(1:4)
h <- seq(1,4)                                   # h are endpoint indicators
p <- matrix(NA, nr=length(t), nc = length(h))   # p are the p-value matrix
#p[1,]<-c(0.0339,0.1411,0.1533,0.0313)
#p[2,]<-c(0.0108,0.0183,0.0007,0.00002)
p[1,]<-c(0.13,NA,0.23,NA)
p[2,]<-c(0.001,0.12,0.01,0.2)
p[3,]<-c(0.001,0.002,0.002,0.04)
p[4,]<-c(NA,0.03,NA,0.02)
timing<-matrix(0, nr=length(t), nc = length(h))                                                             # timing is the information fraction for IA
timing[1,]<-c(0.6,NA,0.6,NA)
timing[2,]<-c(0.8,0.5,0.8,0.5)
timing[3,]<-c(1,0.71,1,0.71)
timing[4,]<-c(NA,1,NA,1)
sfpar<-c(-4,-4,-15,-15)



# spending function is set to be  Hwang-Shih-DeCani spending function, sfpar is the spending function parameter vector for each H
# debug=1 to print some processing details
gMCPgSD_MS(g=G,w.start=w.start,t=t,h=h,p=p,alpha=0.025,timing=timing,sfpar=sfpar,debug=1)
gMCPgSD_BF(g=G,w.start=w.start,t=t,h=h,p=p,alpha=0.025,timing=timing,sfpar=sfpar,debug=0)

