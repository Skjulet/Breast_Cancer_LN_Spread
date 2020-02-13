#Reads the simulated data set
sizeNodesHRT <- read.csv("Z:/'''insert path'''/simDataSet.csv",header=TRUE,sep=",")
tumsize <- sizeNodesHRT[,1]
nodes <- sizeNodesHRT[,2]
hrt <- sizeNodesHRT[,3]

#Negative log-likelihood function for number of lymph nodes conditional on tumour size and HRT use.
c0 <- (0.5^3)*pi/6
fNatural <- function(par){
	-sum(dnbinom(nodes, exp(par[1]), 
		exp(par[2])/(exp(par[2])+exp(par[3]*(hrt==1))*(log(((tumsize^3)*pi/6)/c0))^5),
		log = TRUE))
}

#Optimizes the negative log likelihood function to obtain parameter estimates.
opt <- optim(c(0.1,0.1,1),fNatural,method="L-BFGS-B")


#Bootstraps 95% confidence intervals
nodes1<-nodes
hrt1<-hrt
tumsize1<-tumsize
numberBoots <- 1000
bootNatural <- matrix(0,numberBoots,3)
for(n in 1:numberBoots){
	s <- sample(length(tumsize1), length(tumsize1), replace = TRUE)
	tumsize<-tumsize1[s]
	nodes<-nodes1[s]
	hrt<-hrt1[s]
	fNatural <- function(par){
		-sum(dnbinom(nodes, exp(par[1]), 
				exp(par[2])/(exp(par[2])+exp(par[3]*(hrt==1))*(log(((tumsize^3)*pi/6)/c0))^5), 
				log = TRUE))
	}
	bootNatural[n,] <- optim(opt$par,fNatural,method="L-BFGS-B")$par
}

#print parameter estimates and 95% confidence intervals of the parameters.
opt$par
sort(bootNatural[,1])[c(51,950)]
sort(bootNatural[,2])[c(51,950)]
sort(bootNatural[,3])[c(51,950)]





























