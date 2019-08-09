# analyze boruta permutations bootstraps
library(UpSetR)
library(grid)

info.bbs <- function(savefn, destpath, summplotname ){
  # get boruta bootstraps info 
  # Arguments
  # * savefn : filename or path to save robject
  # * destpath : path to read in boruta r-results iterables from
  # * summplotname : name or path of summary plot to save
  
  lfli <- list.files(destpath)
  dlbs <- list()
  nns <- seq(1,length(lfli),1)
  for(i in 1:length(nns)){
    load(paste0(destpath,"/",
                lfli[grepl(paste0("_",i,".rda"),lfli)]))
    dlbs[[paste0("iter",i)]] <- bdat
    message(i)
  }
  
  # bsm bootstraps matrix
  genesi <- names(dlbs$iter1$finalDecision)
  bsm <- matrix(nrow=0, ncol=length(genesi))
  for(g in 1:length(dlbs)){
    bsm <- rbind(bsm, as.character(dlbs[[g]]$finalDecision))
  }
  colnames(bsm) <- genesi
  
  # get most frequently confirmed genes
  dff <- data.frame(gid=genesi, stringsAsFactors = F)
  dff$perc.confirmed <- dff$num.confirmed <- dff$num.tent <- dff$num.rej <- 0
  for(i in 1:ncol(bsm)){
    gii <- bsm[,i]
    # counts
    dff$num.conf[i] <- length(gii[gii=="Confirmed"])
    dff$num.tent[i] <- length(gii[gii=="Tentative"])
    dff$num.rej[i] <- length(gii[gii=="Rejected"])
    # percents
    dff$perc.conf[i] <- 100*(length(gii[gii=="Confirmed"])/length(gii))
    dff$perc.tent[i] <- 100*(length(gii[gii=="Tentative"])/length(gii))
    dff$perc.rej[i] <- 100*(length(gii[gii=="Rejected"])/length(gii))
    
    message(i)
  }
  dff <- dff[rev(order(dff$perc.conf)),]
  
  # save r object
  rl <- list("bdat.iter"=dlbs,
             "feature.statfreq"=dff,
             "results.matrix"=bsm)
  save(rl, file=savefn)
  
  # save composite summary plot
  jpeg(summplotname,8,8,units="in",res=400)
  par(mfrow=c(2,2), oma=c(4,2,2,2))
  plot(dff$perc.conf, dff$num.conf,
       xlab=paste0("% Confirmed\n(Num. Bootstrap Iter. = ",nrow(bsm),")"), 
       ylab="Num. Confirmed (N Iter.)")
  plot(dff$perc.rej, dff$perc.conf,
       xlab="% Rejected (N Iter.)",
       ylab="% Confirmed (N Iter.)")
  plot(dff$perc.tent, dff$perc.conf,
       xlab="% Tentative (N Iter.)",
       ylab="% Confirmed (N Iter.)")
  boxplot(dff[,c("perc.conf","perc.tent","perc.rej")],
          ylab="% Iterations", main="Status Freq. by Feature",
          xlab="", cex.axis=0.8,
          las=2)
  dev.off()
}

# nrank
savefn <- "bbs-nrank_rlist.rda"
destpath <- "borutaiter/nrank"
summplotname <- "bbs-nrank_compplot.jpg"
info.bbs(savefn, destpath, summplotname)

# random forest
savefn <- "bbs-rf_rlist.rda"
destpath <- "borutaiter/rf"
summplotname <- "bbs-rf_compplot.jpg"
info.bbs(savefn, destpath, summplotname)

# svm
savefn <- "bbs-svm_rlist.rda"
destpath <- "borutaiter/svm"
summplotname <- "bbs-svm_compplot.jpg"
info.bbs(savefn, destpath, summplotname)

# xgb
savefn <- "bbs-xgb_rlist.rda"
destpath <- "borutaiter/xgb"
summplotname <- "bbs-xgb_compplot.jpg"
info.bbs(savefn, destpath, summplotname)

# lasso
savefn <- "bbs-lasso_rlist.rda"
destpath <- "borutaiter/lasso"
summplotname <- "bbs-lasso_compplot.jpg"
info.bbs(savefn, destpath, summplotname)

#-------------
# upset plots
#-------------
fl <- c("bbs-nrank_rlist.rda",
        "bbs-svm_rlist.rda",
        "bbs-rf_rlist.rda",
        "bbs-xgb_rlist.rda",
        "bbs-lasso_rlist.rda")

# no cutoff, >1x confirmed
# genesi <- names(rl$bdat.iter$iter1$finalDecision)
load("~/scratch/consensusML/boruta_supplement/bbs-nrank_rlist.rda")
genesi <- names(rl$bdat.iter$iter1$finalDecision)
um <- data.frame(geneid=genesi, stringsAsFactors = F)
um$nrank <- um$svm <- um$rf  <- um$xgb <- um$lasso <- 0
for(f in 1:length(fl)){
  load(fl[f])
  gff <- rl$feature.statfreq
  gffid <- gff[gff$num.conf>0,]$gid
  um[um$geneid %in% gffid, f+1] <- 1
  message(f)
}

# get total consensus genes by algo
table(um$lasso) # lasso = 501
table(um$xgb) # xgb = 47
table(um$rf) # rf = 376
table(um$svm) # svm = 485
table(um$nrank) # nrank = 106

jpeg("nrank-4algo_gconfall_upsetr.jpg", 7,5, units="in", res=400)

upset(um, order.by = "freq")
grid.text("Feature Overlap (All, >1 Confirmed)", x=0.6, 
          y=0.95, gp=gpar(fontsize=10))

dev.off()

# cutoff at 20% of iterations confirmed
um20 <- data.frame(geneid=genesi, stringsAsFactors = F)
um20$nrank <- um20$svm <- um20$rf  <- um20$xgb <- um20$lasso <- 0
for(f in 1:length(fl)){
  load(fl[f])
  gff <- rl$feature.statfreq
  gffid <- gff[gff$perc.conf>20,]$gid
  um20[um20$geneid %in% gffid, f+1] <- 1
  message(f)
}

jpeg("nrank-4algo_gconf20perc_upsetr.jpg", 7,5, units="in", res=400)
upset(um20, order.by = "freq")
grid.text("Feature Overlap (>20% Iter. Confirmed)", x=0.6, 
          y=0.95, gp=gpar(fontsize=10))
dev.off()

# cutoff at 5% of iterations confirmed
um5 <- data.frame(geneid=genesi, stringsAsFactors = F)
um5$nrank <- um5$svm <- um5$rf  <- um5$xgb <- um5$lasso <- 0
for(f in 1:length(fl)){
  load(fl[f])
  gff <- rl$feature.statfreq
  gffid <- gff[gff$perc.conf>5,]$gid
  um5[um5$geneid %in% gffid, f+1] <- 1
  message(f)
}

jpeg("nrank-4algo_gconf5perc_upsetr.jpg", 7,5, units="in", res=400)
upset(um5, order.by = "freq")
grid.text("Feature Overlap (>5% Iter. Confirmed)", x=0.6, 
          y=0.95, gp=gpar(fontsize=10))
dev.off()

#--------------------------
# consensus genes summary
#--------------------------
library(SummarizedExperiment)
load("sesetfilt_degseahack_targetaml.rda")

# annotate times confirmed
um.all.nrank <- um[um$nrank==1,]
genes.cany <- um.all.nrank[,1]
gdatall <- as.data.frame(rowData(degfilt.se))
cgdat <- as.data.frame(gdatall[genes.cany,])
cgdat$nrank.conf.n1 <- 1
cgdat$nrank.conf.5perc <- 0
cgdat$nrank.conf.20perc <- 0
cgdat[rownames(cgdat) %in% um5[um5$nrank==1,]$geneid,]$nrank.conf.5perc <- 1
cgdat[rownames(cgdat) %in% um20[um20$nrank==1,]$geneid,]$nrank.conf.20perc <- 1
write.csv(cgdat,file="consensus-genes_borutaboot-nrank-results.csv")

#--------------------
# imp hist composite
#--------------------

lfl = list.files()
lfl.filt <- lfl[grepl("rlist",lfl)]

# select random iters to plot
fhl <- list()
for(i in 1:length(lfl.filt)){
  fni = lfl.filt[i]
  load(fni)
  iterl <- rl$bdat.iter
  isel <- sample(length(iterl), 4)
  iname <- gsub("bbs-|_.*","",fni)
  fhl[[iname]] <- iterl[c(isel)]
  message(i)
}

jpeg("feathist_composite_5algo.jpg", 12,12,units="in",res=400)
par(mfcol=c(4,5))
for(l in 1:length(fhl)){
  fii <- fhl[[l]]
  titlei <- names(fhl)[l]
  for(i in 1:length(fii)){
    if(i==1){
      plotImpHistory(fii[[i]], main=titlei)
    } else{
      plotImpHistory(fii[[i]])
    }
  }
  message(l)
}
dev.off()


