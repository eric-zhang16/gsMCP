# gsMCP
This is an example for how to uapply a modified Simes test to a group sequential trial. 

## Installation

Install the released version of gsMCP from GitHub with:

``` r
devtools::install_github("eric-zhang16/gsMCP")
library(gsMCP)
```
## Set up initial alpha weights and transition weights per the below design

![Testing Strategy](https://github.com/eric-zhang16/gsMCP/blob/main/design.PNG?raw=true)

G is transition matrix and w.start is the starting weights
``` r
G <- matrix(0, nr=4, nc = 4)                                  # G is transition matrix
G[1,2] <- G[1,3] <- G[2,1] <- G[2,4] <- 1/2
G[3,2] <- G[4,1] <- 1
w.start <- c(0.5,0.5,0,0)       
```
## Set up hypotheses and IA timeline
h is a vector of hypothesis indicators, t is a vector of IA/FA indicators. Assume 4 hypotheses with IA and one FA
``` r
h <- seq(1,4)
t - seq(1,2)
```
Assume information fraction is 0.76 for H1 and H2, and 0.7 for H3 and H4
``` r
timing<-matrix(0, nr=length(t), nc = length(h))                                                            
timing[,1]<-c(0.76,1)
timing[,2]<-c(0.76,1)
timing[,3]<-c(0.7,1)
timing[,4]<-c(0.7,1)
```
## specify spending functions
gsMCP uses Hwang-Shih-DeCani spending function. Users need to specify the gamma parameter
sfpar<- c(-2,-2,-4,-4)


