# simulation 
# all tissues sharing snps contribute to trait, but the effect size is random, 



#  type 1 err under NULL (hz=0)
# parallel --dryrun Rscript evalType1Err_args_diffeff_alltraits_011421.R  --m 200 --p 300 --n 50000 --Ti 10 --Nsim 200  --hc 0.05 --hz 0 --s 0.05 --shares 0.3 --shareTi {1} --rhoX 0.1 --rhoE 0.5 --batch {2} ::: 1 2 4 6 8 10 ::: $(seq 1 1 10) > tmpjob
#   power (hz > 0 )
# parallel --dryrun Rscript evalType1Err_args_diffeff_alltraits_011421.R  --m 200 --p 300 --n 50000 --Ti 10 --Nsim 200  --hc 0.05 --hz 0.05 --s 0.05 --shares 0.3 --shareTi {1} --rhoX 0.1 --rhoE 0.5 --batch {2} ::: 1 2 4 6 8 10 ::: $(seq 1 1 10) > tmpjob

# use R 3.6
.libPaths("R/rlib-3.6.0")
library(mvtnorm)
library(GBJ)
library(Rcpp)
library(TisCoMM)
library(optparse)
source("functions/functions.r")
source("functions/multiEstB.R")
sourceCpp("functions/PXem_ss.cpp")
sourceCpp("functions/PXem.cpp")

# a file with cca script
source("functions/cca_092930.r")



option_list <- list(
    
    make_option(c("--m"), action="store", # training omic sample number
                default=200,
                type='integer', help="number of omic samples [default: %default]"),
    make_option(c("--batch"), action="store", # training omic sample number
                default=1,
                type='integer', help="batch [default: %default]"),
    make_option(c("--n"), action="store", # trait sample number
                default=5000,
                type='integer', help="number of GWAS trait samples [default: %default]"),
    #make_option(c("--nref"), action="store", # trait sample number
    #            default=5000,
    #            type='integer', help="number of ref panel samples [default: %default]"),
                            
    make_option(c("--p"), action="store", #total snp number
                default=300,
                type='integer', help="number of dimension [default: %default]"),
    make_option(c("--Ti"), action="store", # nonzero snp number
                default=10,
                type='integer', help="number of tissues [default: %default]"),
   
	make_option(c("--Nsim"), action="store", # number of simulations
                default=100,
                type='integer', help="number of simulations expression[default: %default]"),		
    make_option(c("--hc"), action="store", # signal-to-noise ratio/heritibility of splicing
                default=0.05, 
                type='double', help="cell level h2: signal (fixed part) to noise ratio [default: %default]"),
    make_option(c("--hz"), action="store", # signal-to-noise ratio/heritibility of trait
                default=0.01, 
                type='double', help="trait level h2: signal (fixed part) to noise ratio [default: %default]"),
    make_option(c("--s"), action="store", # fraction of nonzero
                default=0.1, 
                type='double', help="sparsity (fraction of snps nonzero effects on splicing/expression)[default: %default]"),
     make_option(c("--shares"), action="store", # fraction of nonzero
                default=0.5, 
                type='double', help="homogeneity (fraction of nonzero snps to work in all tissues)[default: %default]"),         
                
     make_option(c("--shareTi"), action="store", # nonzero snp number
                default=2,
                type='integer', help="number of tissues with shared effects [default: %default]"),
                
    make_option(c("--rhoX"), action="store", # correlation between x (ld) 
                default=0, 
                type='double', help=" correlation between snps genotype [default: %default]"),
    make_option(c("--rhoE"), action="store", # correlation between errors 
                default=0.5, 
                type='double', help=" correlation between tissues in error [default: %default]")

)


######################################################
#######        parse arguments ########################
#######################################################

opt <- parse_args(OptionParser(option_list=option_list))

# training sample
m = opt$m
# snp count
p = opt$p
# gwas sample
n = opt$n
# reference panel sample
#nref = opt$nref
# tissue 
Ti = opt$Ti
nc = opt$nc
# number of sim
Rep = opt$Nsim
# sparsity
s = opt$s
# shared fraction of s
shares = opt$shares
# tissues with shared effects
shareTi = opt$shareTi

# expression h2
h_c = opt$hc
# trait h2
h_z = opt$hz

# correlation between tissues
rhoX = opt$rhoX
# correlation in error
rhoE = opt$rhoE

# batch 
batch = opt$batch

print(c(h_z,h_c,rhoX,rhoE,s,shares,Ti,shareTi,Rep,batch))

print(paste0("pvalue_hz",h_z*10, "_hc", h_c*10, "_rhoX", rhoX*10,"_rhoE",rhoE*10, "_s", s*10, "_shares",shares*10, "_Ti",Ti,"_shareTi",shareTi, "_rep",Rep, "_batch-", batch,".txt"))


evalType1Err <- function(m,n,p, h_z, h_c, s,shares,shareTi, rhoX, rhoE, Ti, batch, Rep = 1000)
{
	PX <- T
	tic <- Sys.time()
	#pvalue <- matrix(0, ncol=5, nrow=Rep)
	pvalue <- matrix(0,ncol=8,nrow=Rep)
	
	# number of true snps found
	# strue <- mtrue <- mccatrue <- list()
# 	sest <- mest <- mccaest <- list()
# 	
 	# mammot, mammot_ss, multixcan, multixcan_ss, utmost.
	for(i in 1:Rep)	{
#		m = 400
     #   m = 200
      #  n = 5000
		#p = 300

# h_z <- c(0, 0.01)[jobs[1, number]]
# cellular h2: snp to expression
# 	h_c <- c(0, 0.05, 0.1)[jobs[2, number]]
# sparsity: fraction of snp with nonzero effects on expression
# 	s = c(0.1, 0.3, 0.5)[jobs[3, number]]
# # correlation between tissues
# 	rhoX = c(0.2, 0.5, 0.8)[jobs[4, number]]
#  	batch <- (1:10)[jobs[5, number]]
# correlation between errors in expression
# representing the shared environmental or other non-genetic factors among over- lap samples in 10 tissues
# rhoE = 0.5
# 	Ti = 10
# 	Rep = 10
# evalType1Err(h_z, h_c, s, rhoX, rhoE, Ti, batch, Rep)

		# genotypes in eQTL
		x1 <- x2snp(rmvnorm(m, mean=rep(0, p), sigma = AR(rhoX, p)))
		x1mean <- matrix(rep(colMeans(x1), m), ncol=p, byrow = T)
		x1 <- (x1 - x1mean)

		# genotypes in GWAS
		x2 <- x2snp(rmvnorm(n, mean=rep(0, p), sigma = AR(rhoX, p)))
		x2mean <- matrix(rep(colMeans(x2), n), ncol=p, byrow = T)
		x2 <- (x2 - x2mean)
 
 		# eQTL coefficient matrix, across Ti tissues.
		b <- rnorm(p) # just 300 numbers, shared effects
		W <- matrix(0, nrow = p, ncol=Ti) # tissue specific, so have effect =1 
		
		# shared count
		sharen = round(s*shares*p,0)
		print(paste("share n:",sharen))
		
		nshn = round( s*(1-shares)*p,0)
		print(paste("share Ti:",shareTi))
		
  		
  		W[sample(p, sharen, replace=FALSE),1:shareTi] <- runif(shareTi,  min=-1, max=1)
		# random fraction of nonzero effects in tissues
		for (j in 1:Ti)	{ # sample some entires to be 1
			W[sample(p, nshn, replace=FALSE),j] <- runif(1,min=-1,max=1)
		}
		# not shared tissues:
		if (shareTi+1 < Ti){
			for (j in (shareTi+1):Ti){	 # make up for not shared tissues: make sure each tissue get same fraction of 1s
				W[sample(p, sharen, replace=FALSE),j] <- runif(1,min=-1,max=1)
			}
		}
  		
#   		
		B <-  diag(b) %*% W
		lp <- x1 %*% B # linear part

		# error matrix in eQTL.
		if (Ti != 1){
		# h_c : expression h2
			sd.lp = diag(sqrt(diag(var(lp)) * (1/h_c - 1)))
			Ve = sd.lp %*% AR(rhoE, Ti) %*% sd.lp
			E <- rmvnorm(m, mean = rep(0, Ti), sigma = Ve)
		} else {
			sd.lp = sqrt(diag(var(lp)* (1/h_c - 1)))
			E <- rnorm(m, sd = sd.lp)
		}
		y <- lp + E

# null: h_z = 0, z is just random
		if (h_z == 0) z <- rnorm(n, 0, sqrt(3)) else {
			a <- runif(Ti, -1, 1) # alpha: exp on trait, 10 uniform numbers -1 to 1
			#a <- rep(0,Ti)
			# nc is the number of causal splicing isoforms (number of splicing have effects on trait)
			#a[1:shareTi] = 1
			
			lp_z <- x2 %*% B %*% a #linear part of z
			#h_z <- .01
			sig2 <- c(var(lp_z)) * ((1-h_z)/h_z) 
			z <- lp_z + rnorm(n, 0, sqrt(sig2))
  		}

 		maxIter = 3000
 		
  		n_p <- 400
		x3 <- rmvnorm(n_p, mean=rep(0, p), sigma = AR(rhoX, p))
		x3mean <- matrix(rep(colMeans(x3), n_p), ncol=p, byrow = T)
		#x3sd <- matrix(rep(apply(x3, 2, sd), n_p), ncol=p, byrow = T)
		x3 <- (x3 - x3mean)#/x3sd/sqrt(n_p)

		 

		lam = 0.95
 		sumx3 = apply(x3*x3, 2, sum)
		RR = matrix(0, p, p);
		for (i1 in 1:p){
    		for (j1 in 1:p){
        		RR[i1, j1] = t(x3[,i1])%*%x3[,j1]/sqrt(sumx3[i1]*sumx3[j1])                        
    		}
		}
 
		R = RR*lam + (1 - lam)*diag(p)

	
		# summary statisitcs for GWAS	
 		hatmu = matrix(0, p, 1)
		hats  = matrix(0, p, 1)
		hatp = matrix(0,p,1)
		for (j in 1:p){
    		fm <- lm(z ~ 1 + x2[, j]);
    		hatmu[j] <- summary(fm)$coefficients[2,1]
    		hats[j]  <- summary(fm)$coefficients[2,2]
    		hatp[j] <- summary(fm)$coefficients[2,4]
		}


		
		
 
		

		# multiXcan (individual data)
		cond <- 30
		print(MulTiXcan(x1, y, x2, z, cond)$pval)
		print(dim(pvalue))
		print(pvalue[i,])
 		pvalue[i,1] <- MulTiXcan(x1, y, x2, z, cond)$pval	


	
		
		pvalue[i,4] <- cca(x1,y,x2,z)$pval

		# cca: summary level
		ccao = ccapa(x1,y)

    	B_hat1 <- ccao[[1]] # cca weight


		if (sum(abs(B_hat1)) > 1e-6) {
			B_hat1 <- as.matrix(B_hat1)
			Ge_impute1 <- x3%*%B_hat1
			eta <- apply(Ge_impute1, 2, sd)
			#sig_z <- sd(z)
			#sig_x <- sqrt(diag(var(x3)))
			sig_x <- sqrt(diag(var(x3)))

			Ztilde <- hatmu/hats
			#Ztilde <- hatmu/hats
			Lambda1 <-  diag(sig_x) %*% sweep(B_hat1, 2, eta, "/") # 300 x 1
			Lambda_sub1 <- Lambda1[, which(eta != 0)]
			# num 1:300
			test_stats1 <- t(Lambda_sub1) %*% Ztilde
			
			
			rhoGE1 <- cor(Ge_impute1)
			rhoGE_svd1 <- svd(rhoGE1)		
			ind_top <- 1

			
			u <- rhoGE_svd1$u
			v <- rhoGE_svd1$v
			us <- as.matrix(u[, ind_top])
			vs <- as.matrix(v[, ind_top])
			d <- rhoGE_svd1$d[ind_top]
			if(length(d) > 1) ds <- diag(1/d) else ds = matrix(1/d)
			rhoGE_ginv <- vs %*% ds %*% t(us)
			chisq <- c(t(test_stats1) %*% rhoGE_ginv %*% test_stats1)
			dof = length(ind_top)
			ccamp <- pchisq(chisq,length(ind_top),lower.tail=F) 

			pvalue[i,5] <- ccamp
			
		} else{
		
			pvalue[i,5] <- NA
			
		}


		

		# mcca individual data
		
		mccaind <- function(x1, y, x2, z, cond = 30){
			n <- nrow(x2)
			Ti <- ncol(y)
			#yhat2 <- matrix(0, nrow = n, ncol = Ti)
		
			mod = CCA(x1,y,K=ncol(y))
			yhat2 = x2%*%as.matrix(mod$u)
			
			yhat2_var <- apply(yhat2, 2, var)
			if (all(yhat2_var == 0))	{
				alpha <- 0
				pval <- 1
			}	else  if (sum(yhat2_var != 0) == 1) { # only 1 tissue nonzero
				twas <- summary(lm(z ~ yhat2))
				alpha <- twas$coefficients[2,1]
				pval <- twas$coefficients[2,4]
			} 	else  { #more than 1 tissue nonzero
				yhat2 <- yhat2[, yhat2_var != 0]

				yhat2_pca <- summary(prcomp(yhat2, center = TRUE, scale. = TRUE))
				ind_top <- which(yhat2_pca$importance[2, 1]/yhat2_pca$importance[2, ] < cond)
				pc_top <- yhat2_pca$x[,ind_top]
				twas <- summary(lm(z ~ pc_top))
				alpha <- twas$coefficients[-1]
			
		 
				pval <- pf(twas$fstatistic[1],twas$fstatistic[2],twas$fstatistic[3],lower.tail=FALSE)
			} # f-test, 
			list(alpha = alpha, 
				 pval  =  unname(pval))
		}
      
      	pvalue[i, 6] <- mccaind(x1, y, x2, z, cond)$pval
		
		

		# mcca summary data
		
		mod = CCA(x1,y,K=ncol(y))

		B_hat1 = mod$u
			
		mccasum <- function(B_hat1,nref,hatmu,hats){
			
		
			n_p <- nref
			x3 <- rmvnorm(n_p, mean=rep(0, p), sigma = AR(rhoX, p))
			x3mean <- matrix(rep(colMeans(x3), n_p), ncol=p, byrow = T)
			#x3sd <- matrix(rep(apply(x3, 2, sd), n_p), ncol=p, byrow = T)
			x3 <- (x3 - x3mean)#/x3sd/sqrt(n_p)

		
			if (sum(abs(B_hat1)) > 1e-6) {
				B_hat1 <- as.matrix(B_hat1[, apply(B_hat1, 2, sd)!= 0])
			
			
				Ge_impute1 <- x3%*%B_hat1
				eta <- apply(Ge_impute1, 2, sd)
				#sig_z <- sd(z)
				#sig_x <- sqrt(diag(var(x3)))
				sig_x <- sqrt(diag(var(x3)))

				Ztilde <- hatmu/hats
				#Ztilde <- hatmu/hats
				Lambda1 <-  diag(sig_x) %*% sweep(B_hat1, 2, eta, "/")
				Lambda_sub1 <- Lambda1[, which(eta != 0)]
				# 333 6
				test_stats1 <- t(Lambda_sub1) %*% Ztilde
				rhoGE1 <- cor(Ge_impute1)
				rhoGE_svd1 <- svd(rhoGE1)		
				ind_top <- which(rhoGE_svd1$d[1]/rhoGE_svd1$d < cond)

				if(length(ind_top) == 0) ind_top <- 1
				u <- rhoGE_svd1$u
				v <- rhoGE_svd1$v
				us <- as.matrix(u[, ind_top])
				vs <- as.matrix(v[, ind_top])
				d <- rhoGE_svd1$d[ind_top]
				if(length(d) > 1) ds <- diag(1/d) else ds = matrix(1/d)
				rhoGE_ginv <- vs %*% ds %*% t(us)
				mcca.chisq <- c(t(test_stats1) %*% rhoGE_ginv %*% test_stats1)
				mcca.dof = length(ind_top)
				mcca.p <- pchisq(mcca.chisq,length(ind_top),lower.tail=F) 

				return(mcca.p)
			
			} else{
		
				return(NA)
			
			}
	}

		pvalue[i,7] = mccasum(B_hat1=B_hat1,nref=400,hatmu=hatmu,hats=hats)
		
		pvalue[i,8] = mccasum(B_hat1 =B_hat1,nref=5000,hatmu=hatmu,hats=hats)
		
	
		# utmost (summary data)
		B_hat = multiEstB(x1, y)
		
		
		
		if (sum(abs(B_hat)) > 1e-6) {
			B_hat <- as.matrix(B_hat[, apply(B_hat, 2, sd)!= 0])
			Ge_impute <- x3 %*% B_hat
			eta <- apply(Ge_impute, 2, sd)
			#sig_z <- sd(z)
			sig_x <- sqrt(diag(var(x3)))
			Ztilde <- hatmu/hats
			Lambda <-  diag(sig_x) %*% sweep(B_hat, 2, eta, "/")
			Lambda_sub <- Lambda[, which(eta != 0)]
			test_stats <- t(Lambda_sub) %*% Ztilde
			cov_Z <- t(Lambda_sub) %*% R %*% Lambda_sub
			pvalue[i, 3] <- GBJ(test_stats=test_stats, cor_mat=cov_Z)$GBJ_pvalue

			# multiXcan (summary data)
 			rhoGE <- cor(Ge_impute)
			rhoGE_svd <- svd(rhoGE)
			ind_top <- which(rhoGE_svd$d[1]/rhoGE_svd$d < cond)
			if(length(ind_top) == 0) ind_top <- 1
			u <- rhoGE_svd$u
			v <- rhoGE_svd$v
			us <- as.matrix(u[, ind_top])
			vs <- as.matrix(v[, ind_top])
			d <- rhoGE_svd$d[ind_top]
			if(length(d) > 1) ds <- diag(1/d) else ds = matrix(1/d)
			rhoGE_ginv <- vs %*% ds %*% t(us)
			chisq <- c(t(test_stats) %*% rhoGE_ginv %*% test_stats)
			pvalue[i, 2] <- pchisq(chisq,  length(ind_top), lower.tail=F)
			
		} else { 
			pvalue[i, 3] <- NA
			pvalue[i, 2] <- NA
		}	
		

		
		
 		
 	}
 	toc <- Sys.time()
 	print(toc - tic)

 	colnames(pvalue) <- c( "multixan", "multixcan_ss", "utmost","cca","cca.sum","mcca.ind","mcca.s.400","mcca.s.5000")
 	print(pvalue)
 	
 
 	
 	
 	
	write.table(pvalue, paste0("pvalue_hz", h_z*10,"_Nsim",Rep, "_hc", h_c*10, "_rhoX", rhoX*10,"_rhoE",rhoE*10, "_s", s*10, "_shares",shares*10, "_Ti",Ti,"_shareTi",shareTi, "_rep",Rep, "_batch-", batch,".txt"), row.names = F, col.names = T)



}

evalType1Err(m,n,p,h_z, h_c, s,shares,shareTi, rhoX, rhoE, Ti, batch, Rep)

