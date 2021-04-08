# /fs0/jiy1/CCA/TisCoMM/simulation/cca_092930.R
# usage: pvalue[i,4] <- cca(x1,y,x2,z)$pval
# sample data loc: /fs0/jiy1/ASD_paper/sim/simgwas_0929/sdata_x1_y_x2_z.092930.RData

# dim(x1)
# [1] 400 300
# > dim(y)
# [1] 400  10
# > dim(x2)
# [1] 5000  300
# > dim(z)
# [1] 5000    1
.libPaths("/home/jiy1/R/rlib-3.6.0")
library(PMA)

ccapa <- function( x1,y ){

tryCatch( {
    perm.out <- CCA.permute(x = x1, z = as.matrix(y),typex="standard",typez="standard", nperms=20,standardize=FALSE)
            # use The z penalty that resulted in the highest z-statistic.
    cout <- CCA(x = x1, z = as.matrix(y) ,typex="standard",typez="standard",K=1,penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,v=perm.out$v.init,standardize=FALSE)

            print(cout)
        #}
 #   str(cout)
        mod <- cout$u
        vec <- cout$v
        #selcor <- sout[sout$V1 == selpx,"V2"]
        selcor <-  cout$cors
        wgt.matrix = matrix(0,nrow=ncol(x1),ncol=1)
        selpx = perm.out$bestpenaltyx

        outpy <- perm.out$bestpenaltyz
        print(outpy)

        wgt.matrix[,1] <- mod
        
        return(list(wgt.matrix,vec, selpx,outpy, selcor))
        },
        error = function(error_message){
        message("error")
        message(error_message)
        return(list(matrix(0,nrow=ncol(x1),ncol=1),c(0,0))) #if not work,just save 0
        }
       )
    }
    
  
cca <- function(x1,y,x2,z){  
	#z <- genz(fixphe, h  = h , design = design)

    ccao = ccapa(x1,y)

    cca.wgt <- ccao[[1]]

    pred = x2 %*% cca.wgt

    model <- lm(z ~ pred)

# beta est, t, p, se
    results <- coef(summary(model))[c(2,6,8,4)]

	pval = results[3]
	beta = results[1]
	t = results[2]
	se = results[4]
	list = list("pval"=pval,"beta"=beta,"t"=t,"se"=se,"cca.wgt"=cca.wgt,"cca.v"=ccao[[2]],"cca.px"=ccao[[3]],"cca.py"=ccao[[4]],"cca.cor"=ccao[[5]])
    return(list)
}



