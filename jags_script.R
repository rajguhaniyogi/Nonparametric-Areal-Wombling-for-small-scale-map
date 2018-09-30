library('rjags')
 
jags <- jags.model('sbn1.txt',
                   data = list('k' = k,
                               't' = time,
                               'Y' = Y,
                               'x' = x,
                               'm0' = m0),
                   n.chains = 1,
                   n.adapt = 100)                                                  
 
update(jags, 10000)                                                                 ## execute the sampler for 10000 iterations
 
dd <- jags.samples(jags,
       c('M1', 'beta0', 'beta1', 'Z', 'theta', 'tau.c', 'sigma.c'),10000)                                                                   

ff <- coda.samples(jags,c('M1', 'beta0', 'beta1', 'Z', 'theta', 'tau.c', 'sigma.c'),
             10000, thin=2)                                                         ## collect post burn-in samples

dd.new <- ff[[1]]
theta.temp <- dd.new[,(5+k+1):(5+2*k)]                                              ## collect post burn-in theta samples
Z.temp     <- dd.new[,(1+1):(1+k)]                                                  ## collect post burn-in phi indices
phi.rec    <- matrix(NA, 10000,k)                                                   

for(ll in 1:10000){
      phi.rec[ll,] <- theta.temp[ll,Z.temp[ll,]]                                    ## collect post burn-in phi samples
}

## functions 

cvv<-function(l,v){
      return(sum(ifelse(arrange_x[,l]==v,1,0)))
}

long.size <- c(1:length(arrange_x[1,]))

## Post processing

for(i in 1:10000){
      a11           <- rep(0,k)
      distinct.thet <- as.numeric(names(table(phi.rec[i,])))
      flag1         <- 0
      lengt         <- as.list(as.numeric(names(table(phi.rec[i,]))))
      r.thet        <- lapply(round(distinct.thet,6),cv,v=round(phi.rec[i,],6))

      for( j in 1:k){
           l <- which(round(distinct.thet,6)==round(phi.rec[i,j],6))
           if(a11[r.thet[[l]][1]]==0){
                 flag1            <- flag1+1
                 a11[r.thet[[l]]] <- rep(flag1,length(l))
           }else{
                 a11[r.thet[[l]]] <- a11[r.thet[[l]]]
                 flag1            <- flag1
            }
      }
      entry  <- as.numeric(lapply(long.size,cvv,v=a11))
      hyp[i] <- which(as.numeric(entry)==length(arrange_x[,1]))
}


connecticut <- read.table("connecticut_hypo_accept.txt")                              ## hypotheses permissible by Connecticut map
conn        <- connecticut$x
conn        <- as.numeric(conn)
prob_hy     <- vector() 

for(hhp in 1:length(conn)){
      prob_hy[hhp] <- length(which(hyp[1001:it]==conn[hhp]))
}
prob_hy     <- prob_hy/sum(prob_hy)
table       <- cbind(conn[(sort(prob_hy,index.return=TRUE)$ix)[511:610]],
                     prob_hy[(sort(prob_hy,index.return=TRUE)$ix)[511:610]])          ## top 100 DP hypothesis corresponding to the Conn. map
arr_x_int   <- arrange_x[,conn]
boundary    <- matrix(c(1,2,2,3,3,4,4,5,5,6,6,7,7,8,1,8,1,7,2,5,2,6,2,7,3,5),
                        13,2,byrow=T)                                                 ## All posiible geographical boundaries in Conn. map
prob_bound  <- vector()

for( boun in 1:nrow(boundary)){
      diff_bound <- which(arr_x_int[boundary[boun,1],]==arr_x_int[boundary[boun,2],])
      prob_bound[boun] <- sum(prob_hy[diff_bound])                                    ## probability of boundary hypotheses
}



### Conditional Autoregressive Normal Competitor

Y       <- t(z)
x       <- t(covariat)
adj     <- c(2, 7, 8,
            1, 3, 5, 6, 7,
            2, 4, 5,
            3, 5,
            2, 3, 4, 6,
            2, 5, 7,
            1, 2, 6, 8,
            1, 7)
weights <- c(1,1,1,
            1, 1, 1, 1, 1,
            1, 1, 1,
            1, 1,
            1, 1, 1, 1,
            1, 1, 1,
            1, 1, 1, 1,
            1, 1)
num     <- c(3,5,3,2,4,3,4,2)
data    <- list("k", "time", "Y", "x", "adj", "weights", "num")
inits   <- list(tau.c=1, tauC=1, beta0=0, beta1=0, 
                phi=c(0,0,0,0,0,0,0,0))
inits   <- list(inits)
parameters <- c("tau.c", "tauC", "beta0", "beta1", "phi")

## Not run: 
## both write access in the working directory and package BRugs required:
schools.sim <- bugs(data, inits, parameters, model.file="winbug_lu_carlin.txt",
                    n.chains = 1, n.iter = 40000, n.burnin=10000, n.thin=1, 
                    bugs.directory="C:/Program Files (x86)/WinBUGS14/", program="OpenBUGS",debug=T)

dd             <- schools.sim$sims.matrix
sp.eff         <- dd[,5:12]
bound_lucarlin <- numeric()

for(boun in 1:nrow(boundary)){
     bound_lucarlin[boun] <- mean(abs(sp.eff[,boundary[boun,1]]-sp.eff[,boundary[boun,2]]))
}

cutoff        <- quantile(bound_lucarlin,c(0.05,0.20,0.30))
prob_2_bound  <- matrix(NA,nrow(boundary),length(cutoff))
uncertain_2   <- matrix(NA,nrow(boundary),length(cutoff))

for(i in 1:length(cutoff)){
    for(boun in 1:nrow(boundary)){
       prob_2_bound[boun,i] <- length(which(abs(phi.rec[,boundary[boun,1]]-phi.rec[,boundary[boun,2]])>cutoff[i]))/
                                      nrow(phi.rec)
       uncertain_2[boun,i]  <- sqrt(prob_2_bound[boun,i]*(1-prob_2_bound[boun,i])/nrow(phi.rec))
    }
}
cbind(boundary,prob_bound,bound_lucarlin,prob_2_bound)


########### boundary detection

cutoff        <- quantile(bound_lucarlin,c(0.05,0.20,0.30))  ## boundary detection cut-offs for Lu and Carlin procedure
prob_lc_bound <- matrix(NA,nrow(boundary),length(cutoff))
uncertain_lc  <- matrix(NA,nrow(boundary),length(cutoff))

for(i in 1:length(cutoff)){
    for(boun in 1:nrow(boundary)){
         prob_lc_bound[boun,i] <- length(which(abs(sp.eff[,boundary[boun,1]]-sp.eff[,boundary[boun,2]])>cutoff[i]))/
                                         nrow(sp.eff)        ## Boundary detection prob. for various cut-offs in Lu and Carlin
         uncertain_lc[boun,i]  <- sqrt(prob_lc_bound[boun,i]*(1-prob_lc_bound[boun,i])/nrow(sp.eff)) 
                                                             ## uncertainties for Lu and Carlin boundary detection procedure
    }
}


