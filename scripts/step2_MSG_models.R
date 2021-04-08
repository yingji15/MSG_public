
# goal: take processed x_all, y_all from gtex, compute multi-omics results
# models: S-MultiXcan, UTMOST, MSG (sCCA) 

# splicing only 
# permute y rows (so person id mixed up)



# IF: generate LD matrix from reference: 
# awk '{print "Rscript step2_MSG_models.R --x" , "data/"$1"."$8".x_all",  "--y","data/"$1"."$8".y_all --ref_ld_chr  r5000/chr  --sumstats MDD2018_2.sumstats.gz --generate_db_and_cov"}'  geneid.nomiss.output  > myjob.120120


# IF: model training used computed LD reference:
# awk '{print "Rscript step2_MSG_models.R --x" , "data/"$1"."$8".x_all",  "--y","data/"$1"."$8".y_all   --sumstats MDD2018_2.sumstats.gz --model_training"}'  geneid.nomiss.output  > myjob.train.120120


start_time <- Sys.time()

# use R 3.6
.libPaths("R/rlib-3.6.0")
library(mvtnorm)
library(GBJ)
library(glmnet)
library(Rcpp)
library(TisCoMM)
library(optparse)
library(plink2R)
# these scripts are adapted from TisCoMM (https://github.com/XingjieShi/TisCoMM)
source("functions/functions.r")
source("functions/multiEstB.R")
sourceCpp("functions/PXem_ss.cpp")
sourceCpp("functions/PXem.cpp")
# this script is adapted from https://github.com/gusevlab/fusion_twas
source("functions/allele_qc.r")
# a script of applying sCCA
source("functions/cca_092930.r")



option_list <- list(
    ## Type : logical, integer, souble complex, character
	make_option(c("--x"), action="store", #dimension
                type='character', help="x matrix genotype"),
     make_option(c("--y"),action="store", #z share with b2
                type='character', help="y phenotype"),
     make_option(c("--sumstats"), action="store", default="clozukscz.sumstats.gz", type='character',help="Path to summary statistics (must have SNP and Z column headers) [required]"),
    make_option(c("--ref_ld_chr"), action="store", 
    #default="1000G.EUR.",
    default=NA,
    type='character',help="Prefix to reference LD files in binary PLINK format by chromosome [required]"),
    make_option(c("--main_dir"), action="store", 
    default="BioVU_5000REF",
    type='character',help="Prefix to reference LD files in RData format by gene id [required]"),
    make_option("--model_training", action="store_true", default=FALSE),
    make_option("--generate_db_and_cov", action="store_true", default=FALSE),
    make_option("--save_model",action="store_true",default=FALSE)
)



######################################################
#######        parse arguments ########################
#######################################################


opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

if(opt$model_training & opt$generate_db_and_cov){
  stop('You can only do model training or generating db and cov at a time')
}
if(!(opt$model_training | opt$generate_db_and_cov)){
  stop('Please use "--model_training" for model training or use "--generate_db_and_cov" for db and cov generating')
}



# x.file="ENSG00000100197.CYP2D6.x_all"
# y.file="ENSG00000100197.CYP2D6.y_all"


x.file = opt$x
print(x.file)
y.file = opt$y
print(y.file)
#output = opt$output




######################
# generate output loc
######################

## Provide the dir name(i.e sub dir) that you want to create under main dir:
output_dir <- "results/"

if (!dir.exists(output_dir)){
    dir.create(output_dir)
} else {
    print("results Dir already exists!")
}

setwd(output_dir)

getwd()

odcmd <- "mkdir all"

system(odcmd,wait=T)

if (opt$save_model){
	odcmd <- "mkdir models"
	system(odcmd,wait=T)
}


######################################################################################
##loading input and filter
######################################################################################
####x_all
if (is.null(x.file) || x.file == "None" || x.file == "NULL") {
	stop("x file is required");
} else {
	x <- read.table(x.file,as.is=T, head = F);
	
	x1<-x[c(1,2,3,4,5,6)] # basic info, ref , alt
	x2<-as.matrix(x[c(-1,-2,-3,-4,-5,-6,-7)])
	x2 <- t(x2)	# 205 429 person x snps
	rm(x)
	
	p <- ncol(x2)
	
	x.na <- sapply(1:p, function(j) {any(is.na(x2[, j]))})
	
    if (any(x.na)) {
      x2 <- x2[, !x.na]
	  p <- ncol(x2)	  
	  x1 <- x1[, !x.na]    
    }
}
############
#y_all
if (is.null(y.file) || y.file == "None" || y.file == "NULL") {
	stop("y file is required");
} else {
	#y <- read.table(y.file)
	y <- read.table(y.file,as.is=T,head=F)
    y3 <- y[!(y[,1] == "Expression"),] # delete expression 

	y1<- y3[1]	# name of this
    #print(y1)
    # splicing only
	y2<- y3[-1]	
    # permute rows of this 
    set.seed(42)
    rows <- sample(nrow(y2))
    y2 <- y2[rows,]
	rm(y)
	
#	y2 <- as.matrix(y2);	
    y2 <- as.matrix(t(y2))
   # print(colnames(y2))
    # row scale: so each omic become mean 0 sd 1
#	y_all <- lapply(1:nrow(y2), function(r) { y2 <- scale(as.matrix(y2[r,])) })
	y_all <- scale(y2)
    #y_all <- lapply(1:col(y2), function(r) { y2 <- scale(as.matrix(y2[,r])) })
	rm(y2)
    
    
}

#n= length(y_all)

x_all <- x2

#x_all <- lapply(1:n,function(k){x2})
rm(x2)



##############
# extract meta info
#################



snps = x1[,c(1:6)]

gidpre = sapply( strsplit(x.file,"/"),'[', length(unlist(strsplit(x.file,"/")) ))

gid = sapply( strsplit(gidpre,".x_all"),'[', 1)

N.tot = nrow(x_all)

chr = unique(snps$V1)





########################
# make cov files
########################

if(opt$generate_db_and_cov){
  print(paste0('INFO generating db and cov files: ',gid))


  # Load in reference data: 1kg
  #e.g.  genos = read_plink("1000G.EUR.22",impute="avg")
	genos = read_plink(paste(opt$ref_ld_chr,chr,sep=''),impute="avg")
# 
  # need to subset to certain snps
  	m = match( snps[,2] , genos$bim[,2] )
	m.keep = !is.na(m)
	snps = snps[m.keep,]
#	wgt.matrix = wgt.matrix[m.keep,,drop=F]
	cur.genos = scale(genos$bed[,m[m.keep]])
	print(dim(cur.genos))
	#cur.genos = scale(genos$bed)
	LD = var(cur.genos)
	if(!dir.exists(paste0(opt$main_dir,'/cov/'))){
    	dir.create(paste0(opt$main_dir,'/cov/'))
  	}
  	outname1 <- paste0(opt$main_dir,'/cov/',paste(gid,"cov.RData",sep="."))
	print(outname1)
	print(str(genos))
	bim = genos$bim
  	save(bim,LD,file=outname1)
  	
}

######################
# sum stats
#######################


if(opt$model_training){
  print(paste0('INFO model training: ',gid))
	s.out2<-data.frame(GID=gid,CHR=chr,stringsAsFactors=FALSE)
  
 # load in LD file
 	load(paste0(opt$main_dir,'/cov/',paste(gid,"cov.RData",sep=".")))
 	
 
# Load in summary stats
	sumstat = read.table(opt$sumstats,head=T,as.is=T)

	# Match summary data to input, record NA where summary data is missing
	# this is just limit to the snps in this chr in ldref
	m = match( bim[,2] , sumstat$SNP )
	sum.missing = is.na(m)
	sumstat = sumstat[m,]
	sumstat$SNP = bim[,2]
	sumstat$A1[ sum.missing ] = bim[sum.missing,5]
	sumstat$A2[ sum.missing ] = bim[sum.missing,6]

	# QC / allele-flip the input and output
	qc = allele.qc( sumstat$A1 , sumstat$A2 , bim[,5] , bim[,6] )

	# Flip Z-scores for mismatching alleles in sumstats
	sumstat$Z[ qc$flip ] = -1 * sumstat$Z[ qc$flip ]
	sumstat$A1[ qc$flip ] = bim[qc$flip,5]
	sumstat$A2[ qc$flip ] = bim[qc$flip,6]

	# Remove strand ambiguous SNPs (if any)
	if ( sum(!qc$keep) > 0 ) {

		bim = bim[qc$keep,]
		#genos$bed = genos$bed[,qc$keep]
		sumstat = sumstat[qc$keep,]
		
	}

# Match up the LDREF SNPs and weights (betwee weight matrix and 1kg genotype)
	# Match up the LDREF SNPs and weights (betwee weight matrix and 1kg genotype)
	m = match( snps[,2] , bim[,2] )
	m.keep = !is.na(m)
	snps = snps[m.keep,]
	x_all=x_all[,m.keep]

	m = match( snps[,2] , bim[,2] )
	m.keep = !is.na(m)

	B_hat = multiEstB(x_all, y_all)

	wgt.matrix = B_hat

	# Remove NAs 
	
	wgt.matrix[is.na(wgt.matrix)] = 0
	print(dim(wgt.matrix))
	print(length(m.keep))
	wgt.matrix = wgt.matrix[m.keep,,drop=F]
	cur.bim = bim[m[m.keep],]
	# Flip WEIGHTS for mismatching alleles
	qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
	wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]
	#rm(snps)

	
	# Match up the LDREF SNPs and the summary stats
	m = match(cur.bim[,2] , sumstat$SNP)

	cur.Z = sumstat$Z[m]

	cur.miss = is.na(cur.Z)


	cur.LD = LD[!cur.miss,!cur.miss]


	##########################
	# MSG: use sCCA 
	############################

	mod = CCA(x_all,y_all,K=ncol(y_all))


	B_hat1 = mod$u
	print(dim(B_hat1))
	print(dim(snps))
	#rownames(B_hat1) = snps[,2]
	print(head(B_hat1))
	sel = which(apply(B_hat1, 2, sd)!= 0)
	#B_hat1 <- as.matrix(B_hat1[!cur.miss, apply(B_hat1, 2, sd)!= 0])
	B_hat1 <- as.matrix(B_hat1[!cur.miss, sel])

	

	# 333 333					
	B_hats = B_hat1
	eta2 = diag(sqrt( t(B_hats) %*% cur.LD %*% B_hats ))

	sig_x2 <- sqrt(diag(cur.LD))
	
	Ztilde <- cur.Z[!cur.miss] 
	Lambda2 <-  diag(sig_x2) %*% sweep(B_hats, 2, eta2, "/")
	Lambda_sub2 <- Lambda2[, which(eta2 != 0)]
	test_stats2 <- t(Lambda_sub2) %*% Ztilde


	cur.Zs = cur.Z[!cur.miss]	

	# now get rhoGE2
	
	rhoGE2 = matrix(0,nrow=ncol(B_hats),ncol=ncol(B_hats))
	for (i in 1:ncol(B_hats)){
		for (j in 1:ncol(B_hats)){
		
			nom = t(B_hats[,i])%*%cur.LD%*%B_hats[,j]
			#print(nom)
			denom = sqrt((t(B_hats[,i])%*%cur.LD%*%B_hats[,i])*(t(B_hats[,j])%*%cur.LD%*%B_hats[,j]))
			rhoGE2[i,j] = nom/denom
	
		}
	
	}


	rhoGE_svd2 = svd(rhoGE2)

	cond = 30

	ind_top <- which(rhoGE_svd2$d[1]/rhoGE_svd2$d < cond)

	if(length(ind_top) == 0) ind_top <- 1
	u <- rhoGE_svd2$u
	v <- rhoGE_svd2$v
	us <- as.matrix(u[, ind_top])
	vs <- as.matrix(v[, ind_top])
	d <- rhoGE_svd2$d[ind_top]
	if(length(d) > 1) ds <- diag(1/d) else ds = matrix(1/d)
	rhoGE_ginv <- vs %*% ds %*% t(us)
	mcca.chisq2 <- c(t(test_stats2) %*% rhoGE_ginv %*% test_stats2)

	mcca.dof = length(ind_top)
	mcca.p <- pchisq(mcca.chisq2,length(ind_top),lower.tail=F) 




	s.out2[,"mcca.p"] = rep(mcca.p,nrow(s.out2))
	s.out2[,"mcca.chisq"] = rep(mcca.chisq2,nrow(s.out2))
	s.out2[,"mcca.dof"] = rep(mcca.dof,nrow(s.out2))

	print(s.out2)
	
	outname1 <- paste0("all/",paste(gid,"results","ccaonly",sep="."))
    print(outname1)
    write.table( s.out2 , quote=F , row.names=F , sep='\t' , file=outname1 )

	


	#######################
	# S-MultiXcan and UTMOST (adapted from TisCoMM implementation)
	#######################

	B_hat.org = B_hat
	B_hat = wgt.matrix # use the flipped wgt

	# summarize results 
	pvalue <- matrix(0,ncol=3,nrow=1)

		if (sum(abs(B_hat)) > 1e-6) {
			B_hat <- as.matrix(B_hat[!cur.miss, apply(B_hat, 2, sd)!= 0])
			
		
			B_hats = B_hat
			eta2 = diag(sqrt( t(B_hats) %*% cur.LD %*% B_hats ))

			sig_x2 <- sqrt(diag(cur.LD))


			Ztilde <- cur.Z[!cur.miss] 
			Lambda2 <-  diag(sig_x2) %*% sweep(B_hats, 2, eta2, "/")
			Lambda_sub <- Lambda2[, which(eta2 != 0)]
			test_stats <- t(Lambda_sub) %*% Ztilde

			# R: LD matrix
			lam = 0.95
			RR = cur.LD
			p = sum(!cur.miss)
	
			R = RR*lam + (1 - lam)*diag(p)

	
			cov_Z <- t(Lambda_sub) %*% R %*% Lambda_sub
			tc.gbj <- GBJ(test_stats=test_stats, cor_mat=cov_Z)

				
		
		
			s.out2[,"tc.gbj.p"] = rep(tc.gbj$GBJ_pvalue,nrow(s.out2))
			s.out2[,"tc.gbj"] = rep(tc.gbj$GBJ,nrow(s.out2))
			s.out2[,"tc.gbj.dof"] = rep(length(test_stats),nrow(s.out2))
		

			# multiXcan (summary data)
			cond <- 30
		
			rhoGE2 = matrix(0,nrow=ncol(B_hats),ncol=ncol(B_hats))
			for (i in 1:ncol(B_hats)){
				for (j in 1:ncol(B_hats)){
		
					nom = t(B_hats[,i])%*%cur.LD%*%B_hats[,j]
					print(nom)
					denom = sqrt((t(B_hats[,i])%*%cur.LD%*%B_hats[,i])*(t(B_hats[,j])%*%cur.LD%*%B_hats[,j]))
					rhoGE2[i,j] = nom/denom
	
				}
	
			}


			rhoGE_svd = svd(rhoGE2)
		
	
			ind_top <- which(rhoGE_svd$d[1]/rhoGE_svd$d < cond)
			if(length(ind_top) == 0) ind_top <- 1
			u <- rhoGE_svd$u
			v <- rhoGE_svd$v
			us <- as.matrix(u[, ind_top])
			vs <- as.matrix(v[, ind_top])
			d <- rhoGE_svd$d[ind_top]
			if(length(d) > 1) ds <- diag(1/d) else ds = matrix(1/d)
			rhoGE_ginv <- vs %*% ds %*% t(us)
			tc.multixcan.chisq <- c(t(test_stats) %*% rhoGE_ginv %*% test_stats)
			tc.multixcan.dof = length(ind_top)
			tc.multixcan.p <- pchisq(tc.multixcan.chisq,length(ind_top),lower.tail=F) 
		
		
			s.out2[,"tc.multixcan.p"] = rep(tc.multixcan.p,nrow(s.out2))
			s.out2[,"tc.multixcan.chisq"] = rep(tc.multixcan.chisq,nrow(s.out2))
			s.out2[,"tc.multixcan.dof"] = rep(tc.multixcan.dof,nrow(s.out2))
		
		
		} else { 
		
		
			s.out2[,"tc.gbj.p"] = rep(NA,nrow(s.out2))
			s.out2[,"tc.gbj"] = rep(NA,nrow(s.out2))
			s.out2[,"tc.gbj.dof"] = rep(NA,nrow(s.out2))
		
		
		
			s.out2[,"tc.multixcan.p"] = rep(NA,nrow(s.out2))
			s.out2[,"tc.multixcan.chisq"] = rep(NA,nrow(s.out2))
			s.out2[,"tc.multixcan.dof"] = rep(NA,nrow(s.out2))
		
		
		}	

outname1 <- paste0("all/",paste(gid,"results","txt",sep="."))
print(outname1)
write.table( s.out2 , quote=F , row.names=F , sep='\t' , file=outname1 )

if (opt$save_model) {
	outname2 <- paste0("models/",paste(gid,"models","RData",sep="."))
	print(outname2)
	save(wgt.matrix, B_hats, rhoGE_svd2, test_stats2, mcca.dof, rhoGE_ginv, ind_top, mcca.chisq2, snps, file=outname2)
	}

}



end_time <- Sys.time()
end_time - start_time


# test one: scz summstats
# Rscript step2_MSG_models.R --x ENSG00000137601.NEK1.x_all --y ENSG00000137601.NEK1.y_all    --model_training --save_model