# Rscript step3_summarize_and_plots.R --list scz.0221.all.txt --silver scz.silver.txt --manplot

# summary results
# first get a file with all results, one gene a row
#awk 'NR==FNR||FNR!=1' all/*txt  > scz.0221.all.txt

.libPaths("R/rlib-3.6.0")
library(optparse)
library(ggvenn)
library(dplyr)
library(ggrepel)
library(ggplot2)


option_list <- list(
    ## Type : logical, integer, souble complex, character
    make_option(c("--list"), action="store", #dimension
                type='character', help="gene scoring list: 1st column gene name, 2nd column score"),
    make_option(c("--silver"), action="store", #dimension
    			default="none",
                type='character', help="silver gene list: 1st column gene name"),
    make_option("--manplot",action="store_true",default=FALSE)
)




opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

tab.file = opt$list


# tab.file = "scz.0221.all.txt"




tab = read.table(tab.file,as.is=T,head=T)


ts = tab
ts[is.na(ts)] <- 1



checkn <- function(tab,col){

        e= nrow(tab[tab[{{col}}] < 0.05/nrow(tab) & tab[{{col}}] >0 & !is.na(tab[{{col}}]),] )

        return(e)

}

tab[,c("HGNC")] = sapply( strsplit(tab$GID,"[.]"),'[',2)

ccag = tab[tab$mcca.p < 0.05/nrow(tab) & tab$mcca.p >0,"HGNC"]
gbjg = tab[tab$tc.gbj.p < 0.05/nrow(tab) & tab$tc.gbj.p >0 & !is.na(tab$tc.gbj.p),"HGNC"]
multg = tab[tab$tc.multixcan.p < 0.05/nrow(tab) & tab$tc.multixcan.p >0 & !is.na(tab$tc.multixcan.p),"HGNC"]
x= list("MSG"=ccag,"UTMOST"=gbjg,"S-MultiXcan"=multg)

# number
ngbj = checkn(tab=ts,col="tc.gbj.p")

print(ngbj)


print("tc.multixcan")
nmulti = checkn(tab=ts,col="tc.multixcan.p")
print(nmulti)


print("mcca")
nmcca = checkn(tab=ts,col="mcca.p")
print(nmcca)

fname = sapply( strsplit(tab.file,".txt"),'[',1)
print( fname)



cnt <- data.frame(method=c("MSG", "UTMOST", "S-MultiXcan"), number=c(nmcca,ngbj,nmulti))
print(str(cnt)) 

pdf(paste0(fname,".barplot.pdf"))
ggplot(cnt, aes(x=method, y=number, fill=method)) + geom_bar(stat="identity")+theme_minimal()+scale_fill_manual(values=c("#0073C2FF", "#EFC000FF", "#868686FF"))
dev.off()

pdf(paste0(fname,".ggvenn.pdf"))
ggvenn( x,  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"), stroke_size = 0.5, set_name_size = 3)
dev.off()





# 26




source("functions/qqunif_092930.r")

pvalue = ts[ ,c("tc.gbj.p","tc.multixcan.p","mcca.p")]


pvalue[is.na(pvalue)] <- 1
print(str(pvalue))

pv <- list( MultiXcan_S = pvalue$tc.multixcan.p,
			  		  UTMOST = pvalue$tc.gbj.p,
			  		  mCCA = pvalue$mcca.p )
	


pdf(paste0(fname,Sys.Date(),".qq.pdf"))
qqunif.plot(pv, xlab="Expected: -log10(p-value)", ylab = "Observed: -log10(p-value)", auto.key=list(corner=c(.95,.05)), main = " ")
dev.off()


# basic man plot

library(qqman)
loc = read.table("geneid.nomiss.output",as.is=T,head=F)

locs <- loc[,c("V5","V8")]

tab2 <- merge(tab, locs, by.x="HGNC",by.y="V8")
tab3=tab2[,c("CHR","V5","mcca.p","HGNC")]
colnames(tab3) = c("CHR","BP","P","HGNC")
pdf(paste0(fname,"basicman.pdf"))
manhattan(tab3, p = "P",snp="HGNC" ,main = "",  col = c("blue4", "orange3"), annotatePval=0.05/nrow(tab),annotateTop = FALSE, suggestiveline = F, genomewideline = T)
dev.off()


save(tab,ccag,gbjg,multg,x,cnt,tab3,file=paste0(fname,Sys.Date(),"siglist_cnt_tab.RData"))



#######################
# man plot
# part of the code is adapted from https://github.com/gamazonlab/MR-JTI
########################

if(opt$manplot){
  	print(paste0('generating manhatten plot: ',fname))
	silver.file = opt$silver
	silver = read.table(silver.file,as.is=T,head=F)
	
	
	tabcomb = tab2
	
	tab2 = tab3
	
	tab2$BP<-as.double(tab2$BP)
	print(head(tab2))
	# 1) compute chromosome length
	chr_len <- tab2 %>% group_by(CHR) %>% summarise(chr_len=max(BP))
	head(chr_len)
	# 2） start location of each chromosome
	chr_pos <- chr_len  %>% mutate(total = cumsum(chr_len) - chr_len) %>% select(-chr_len)
	head(chr_pos) 
	#3) cumulative loc of each chrom
	Snp_pos <- chr_pos %>% left_join(tab2, ., by="CHR") %>% arrange(CHR, BP) %>% mutate( BPcum = BP + total)

	#take a look 
	head(Snp_pos,2)   

	gstar <- intersect(ccag,silver$V1)
	print(gstar)
	
	print(intersect(ccag,gstar))
	tab3 <- Snp_pos %>% mutate( is_highlight=ifelse(HGNC %in% gstar, "yes", "no"))

	
	# add genes want to highlight in plot 

	print(sum(tab3[,"is_highlight"]=="yes"))
	
	X_axis <-  tab3 %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

	pvalue.extreme2 <- function(p) {
	  
	  log10.pvalue <- log10(p)
	
	  mantissa <- 10^(log10.pvalue %% 1)
	  exponent <- log10.pvalue %/% 1
  
	  logp<-log10.pvalue*-1
	  return(logp)
	  ##return(sprintf("p value is %1.2f times 10^(%d)",mantissa,exponent))
	}


	tab3$plog = pvalue.extreme2(tab3$P)

	tab4<-tab3[tab3$plog>log(0.05,10)*-1,]

	print(head(tab4))
	
	p1 <- ggplot(tab4, aes(x=BPcum, y=plog)) +
		   geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
		  # scale_color_manual(values = rep(c("blue4", "orange3"), 22 )) +
		   scale_color_manual(values = rep("blue4", 22 )) +
		   scale_x_continuous( labels = c(as.character(seq(1,15)),' ','17','','19','','21',''), breaks= X_axis$center) + 
		   #scale_x_continuous( label = X_axis$CHR, breaks= X_axis$center ) +
		scale_y_continuous(breaks=c(1,2,3,5,10,20,30,50,100,150,200,250),trans='log10')+
		#   scale_y_continuous(expand = c(0, 0) ) +  
		   # highlight
		   geom_point(data=subset(tab4, is_highlight=="yes"), color="red", size=2) +
		   # don't forget the * between expressions
		   labs(x="chromosome",y= expression("-log"[10]*"(P)"),title="")+
		   # add labels to highlights
	   geom_label_repel( data=subset(tab4, is_highlight=="yes"), aes(label=HGNC), size=2) +
		   theme_bw() + 
		theme(
			 legend.position="none",
			 panel.border = element_blank(),
			 panel.grid.major.x = element_blank(),
			 panel.grid.minor.x = element_blank()
		  )   

	pdf( paste0(fname,"man.pdf"), width = 7,height = 7*0.5)
	print(p1)
	dev.off()

	save(tab4, X_axis,tabcomb,tab2, file=paste0(fname,"man.RData"))

#ts[,"symbol"] = sapply( strsplit(ts$GID,'[.]'),'[',2)

} 
#ts[order(ts$mcca.p),][1:10,c("symbol","tc.gbj.p","tc.multixcan.p","mcca.p")]




