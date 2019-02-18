# Make summarized experiment objects for TARGET AML data
# author: Sean Maden
# contact: maden@ohsu.edu

# Description and Details
# * Read in data tables, annotate gene assay information, 
#   and create summarized experiment (SE) and multi assay experiment (MAE) objects
# * Use '?SummarizedExperiment' for info on the SE object class and methods
# * Objects assume hg19 build alignment for expr data, and so Ens75 ranges were used
# * For details about MAE objects see reference: <http://bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/>

# dependencies
require(SummarizedExperiment)
require(GenomicRanges)
require(edgeR)
require(limma)
require(EnsDb.Hsapiens.v75)
require(MultiAssayExperiment)

# define globals
sys.sep = "/"
data.dir = "data"
seobj.dir = "seobjects"

gene.counts.tablename = "TARGET_NBL_AML_RT_WT_HTSeq_Counts.csv" # to load
clinical.tablename = "AML_assay_clinical.csv" # to load
tmm.counts.tablename = "TARGET_NBL_AML_WT_RT_TMMCPM_log2_Norm_Counts_17k.csv" # to be written
testset.prior.name = "TARGET_AML_Testing_Samples.csv"
geneanno.name = "edb_filt_anno.rda"
deg.tablename = "TARGET_AML_High.Std.Risk_vs_LowRisk_DEGs.csv"

countsseset.name <- "seset_genecounts_targetaml.rda"
tmmseset.name <- "seset_genetmmfilt_targetaml.rda"
degseset.name <- "seset_degseahack_targetaml.rda"
maeobj.name <- "mae_targetaml.rda"

# load data 
clinical <- read.csv(paste0(data.dir,sys.sep,clinical.tablename), 
                     row.names = 1, stringsAsFactors = F)

counts <- read.csv(paste0(data.dir, sys.sep, gene.counts.tablename), 
                   row.names = 1, stringsAsFactors = F)
tmm <- read.csv(paste0(data.dir, sys.sep, tmm.counts.tablename), 
                row.names = 1, stringsAsFactors = F)
degtable <- read.csv(paste0(data.dir,sys.sep,deg.tablename), 
                     row.names = 1, stringsAsFactors = F)

testsamp.prior = read.csv(paste0(data.dir,sys.sep,testset.prior.name), row.names=1,
                          stringsAsFactors = F)

#============
# Preprocess
#============
# filter clinical variables
clinical.filt = clinical[,c(1:84)]

# filter AML clinical samples on tissue type (retain primary only)
clinical.filt = clinical.filt[clinical.filt$Diagnostic.ID %in% c("03A","09A"),]
dim(clinical.filt)
# check for repeated patient ids
summary(as.data.frame(table(clinical.filt$TARGET.USI))[,2]) # max = 1, no repeated ids

# apply sample filters to expr data
cnames.counts = colnames(counts)
counts.filt = counts[,grepl(".*\\.20\\..*",cnames.counts)] # project id filt
patidfilt = substr(colnames(counts.filt),11,16) %in% substr(clinical.filt$TARGET.USI,11,16) # patient id filt
counts.filt = counts.filt[,patidfilt]
tisstypefilt = substr(colnames(counts.filt),18,20) %in% clinical.filt$Diagnostic.ID # tissue sample type filt
counts.filt = counts.filt[,tisstypefilt]
# match clinical and counts data
counts.filt = counts.filt[,order(match(substr(colnames(counts.filt),11,16),
                                       substr(clinical.filt$TARGET.USI,11,16)
)
)
]
identical(substr(colnames(counts.filt),11,16),substr(clinical.filt$TARGET.USI,11,16)) # true
identical(substr(colnames(counts.filt),18,20),clinical.filt$Diagnostic.ID) # true

# add train/test data info
clinical.filt$exptset.seahack <- ifelse(clinical.filt$TARGET.USI %in% testsamp.prior[,1],
                                        "test", "train")
table(clinical.filt$exptset.seahack)

#======================
# Normalize Expression
#======================
counts.sub <- counts.filt
dge <- DGEList(counts = counts.sub)
samp <- ncol(counts.sub)
#Note: used a minium # of samples as 5 to ensure that normalized values will include all DEGs identified with the training set counts. Higher thresholds lead to genes included in DEGs but excluded in the "master" TMM normalized counts. 
keep.dge <- rowSums(cpm(dge) >= 1) >= 5
dge <- dge[keep.dge,] #subset for those genes with cmp >= 1 per gene in AML samples
dge <- calcNormFactors(dge) #Do TMM normalization
dim(dge) # 18243 genes meet these criteria in AML only
cpm <- cpm(dge,log = TRUE, prior.count = 1) # all expression as counts per million at filtered genes

#=================================
# Get Gene Annotations and Ranges
#=================================
# gene annotations and granges objects
edb <- EnsDb.Hsapiens.v75 # columns(edb) to check available annotation info
# example:
# select(edb, keys="TP53", columns=colnames(edb),keytype="SYMBOL")
# head(rownames(counts.filt))
# [1] "ENSG00000000003.13" "ENSG00000000005.5"  "ENSG00000000419.11" "ENSG00000000457.12"
# [5] "ENSG00000000460.15" "ENSG00000000938.11"
# NOTE: ids are transcript ids

# simply use gene-level information for now
genes.edb <- genes(edb)
counts.genes.grdf <- as.data.frame(matrix(ncol=6,nrow=0))
colnames(counts.genes.grdf) <- c("gene.id","gene.symbol","countsdf.id","chr.seqname","start","end")
for(i in 1:nrow(counts.filt)){
  gene.info.i = as.data.frame(genes.edb[gsub("\\..*","",rownames(counts.filt)[i])])
  counts.genes.grdf <- rbind(counts.genes.grdf,data.frame(gene.id=rownames(gene.info.i)[1],
                                                          gene.symbol=gene.info.i$gene_name[1],
                                                          countsdf.id=rownames(counts.filt)[i],
                                                          chr.seqname=gene.info.i$seqnames,
                                                          start=gene.info.i$start,
                                                          end=gene.info.i$end,
                                                          stringsAsFactors = F))
  message(i," perc. complete = ",round(100*(i/nrow(counts.filt)),4),"%")
}
save(counts.genes.grdf, file=paste0(data.dir, sys.sep, geneanno.name))

# se experiments using filtered genes
length(intersect(rownames(counts.filt), counts.genes.grdf$countsdf.id)) # 54713
length(intersect(counts.genes.grdf$countsdf.id, rownames(dge))) # 17637

counts.se <- counts.filt[rownames(counts.filt) %in% counts.genes.grdf$countsdf.id,]
cpm.se <- cpm[rownames(cpm) %in% counts.genes.grdf$countsdf.id,]
dim(cpm.se)
# [1] 17637   145

# order genes for counts se
ganno.counts <- counts.genes.grdf[order(match(counts.genes.grdf$countsdf.id,
                                              rownames(counts.se))),]
identical(counts.genes.grdf$countsdf.id, rownames(counts.se))
ganno.tmm <- counts.genes.grdf[counts.genes.grdf$countsdf.id %in% rownames(cpm.se),]
ganno.tmm <- ganno.tmm[order(match(ganno.tmm$countsdf.id,rownames(cpm.se))),]
identical(ganno.tmm$countsdf.id,rownames(cpm.se))
colnames(ganno.counts) <- colnames(ganno.tmm) <- c("gene.id","gene.symbol","countsdf.id","seqnames","start","end")

ggr.counts <- makeGRangesFromDataFrame(ganno.counts, 
                                       keep.extra.columns = T, 
                                       ignore.strand = T)
names(ggr.counts) <- ggr.counts$countsdf.id
ggr.tmm <- makeGRangesFromDataFrame(ganno.tmm,
                                    keep.extra.columns = T,
                                    ignore.strand = T)
names(ggr.tmm) <- ggr.tmm$countsdf.id

#====================================
# Make Summarized Experiment Objects
#====================================
# Gene Expr Counts SE object
identical(ggr.counts$countsdf.id, rownames(counts.se)) # TRUE
identical(names(ggr.counts), rownames(counts.se)) # TRUE
identical(substr(colnames(counts.se),11,16), 
          substr(clinical.filt$TARGET.USI,11,16)) # TRUE
counts.seset <- SummarizedExperiment(assays = as.matrix(counts.se),
                                     rowRanges = ggr.counts, 
                                     colData = DataFrame(clinical.filt, 
                                                         row.names = colnames(counts.se)
                                     ),
                                     metadata = list(dataset = "TARGET_AML", 
                                                     assay_source = "GDC",
                                                     genome_build = "hg19")
)

# Gene TMM SE object
identical(ggr.tmm$countsdf.id, rownames(cpm.se)) # TRUE
identical(names(ggr.tmm), rownames(cpm.se))
identical(substr(colnames(cpm.se),11,16), 
          substr(clinical.filt$TARGET.USI,11,16)) # TRUE
tmm.seset <- SummarizedExperiment(assays = as.matrix(cpm.se),
                                  rowRanges = ggr.tmm,
                                  colData = DataFrame(clinical.filt,
                                                      row.names = colnames(cpm.se)
                                  ),
                                  metadata = list(dataset = "TARGET_AML",
                                                  assay_source = "GDC",
                                                  genome_build = "hg19",
                                                  normalization_strategy = "TMM, log_cpm, limma, edgeR"))

# DEG TMM SE object
deglist = rownames(degtable)
length(intersect(deglist, counts.genes.grdf$countsdf.id)) # 1937 of 1998
degfilt = deglist[deglist %in% counts.genes.grdf$countsdf.id]
ggr.deg = ggr.counts[names(ggr.counts) %in% degfilt]
# deg.assay <- counts.se[rownames(counts.se) %in% degfilt,]
deg.assay <- cpm.se[rownames(cpm.se) %in% degfilt,]
ggr.deg <- ggr.deg[order(match(names(ggr.deg), rownames(deg.assay)))]
identical(names(ggr.deg), rownames(deg.assay)) # TRUE
identical(substr(colnames(cpm.se),11,16), 
          substr(clinical.filt$TARGET.USI,11,16)) # TRUE
# add the deg statistics to gene annotation
degstats = degtable[rownames(degtable) %in% degfilt,]
degstats = degstats[order(match(rownames(degstats), names(ggr.deg))),]
identical(rownames(degstats), names(ggr.deg)) # TRUE
ggr.deg$logFC <- degstats$logFC
ggr.deg$AveExpr <- degstats$AveExpr
ggr.deg$t <- degstats$t
ggr.deg$p.unadj <- degstats$P.Value
ggr.deg$p.adj.bh <- degstats$adj.P.Val
ggr.deg$b <- degstats$B
# make the se object
deg.seset <- SummarizedExperiment(assays = as.matrix(deg.assay),
                                  rowRanges = ggr.deg,
                                  colData = DataFrame(clinical.filt,
                                                      row.names = colnames(deg.assay)
                                  ),
                                  metadata = list(dataset = "TARGET_AML",
                                                  assay_source = "GDC",
                                                  genome_build = "hg19",
                                                  normalization_strategy = "DEGs_trainset_binaryrisk, Low=0 notLow=1, reference: Low, tmm_log_cpm, voom_DE function"))

# Save new SE objects
save(counts.seset, file=paste0(seobj.dir, sys.sep, countsseset.name))
save(tmm.seset, file=paste0(seobj.dir, sys.sep, tmmseset.name))
save(deg.seset, file=paste0(seobj.dir, sys.sep, degseset.name))

#===============================
# Multi Assay Experiment class
#===============================
counts.map <- data.frame(primary = colnames(counts.seset),
                         colname = colnames(counts.seset),
                         stringsAsFactors = F)
tmm.map <- data.frame(primary = colnames(tmm.seset),
                      colname = colnames(tmm.seset),
                      stringsAsFactors = F)
deg.map <- data.frame(primary = colnames(deg.seset),
                      colname = colnames(deg.seset),
                      stringsAsFactors = F)
listmap <- list(counts.map, tmm.map, deg.map)
names(listmap) <- c("gene_counts", "tmm_log2norm_filtcounts", "deg_tmm_counts")
dfmap <- listToMap(listmap)
objlist = list("gene_counts" = counts.seset,
               "tmm_log2norm_filtcounts" = tmm.seset,
               "deg_tmm_counts" = deg.seset)
patient.data <- clinical.filt
rownames(patient.data) <- colnames(counts.se)
targetaml.mae <- MultiAssayExperiment(objlist, patient.data, dfmap)
save(targetaml.mae, file=paste0(seobj.dir, sys.sep, maeobj.name))