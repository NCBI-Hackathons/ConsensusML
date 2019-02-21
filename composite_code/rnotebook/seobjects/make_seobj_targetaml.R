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
require(biomaRt)
require(MultiAssayExperiment)
# require(EnsDb.Hsapiens.v86)

# define globals
sys.sep = "/"
data.dir = "data"
seobj.dir = "seobjects"

gene.counts.tablename = "TARGET_NBL_AML_RT_WT_HTSeq_Counts.csv" # to load
clinical.tablename = "AML_assay_clinical.csv" # to load
tmm.counts.tablename = "TARGET_NBL_AML_WT_RT_TMMCPM_log2_Norm_Counts_17k.csv" # to be written
testset.prior.name = "TARGET_AML_Testing_Samples.csv"
geneanno.name = "edb86_filt_anno.rda"
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
require(biomaRt)

eids <- gsub("\\..*","",rownames(cpm))
dfens <- getBM(attributes = c("chromosome_name","start_position","end_position","strand",
                              "hgnc_id", "hgnc_symbol", 
                              "ensembl_gene_id"),
               filters = 'ensembl_gene_id', 
               values = eids, 
               mart = ensembl)
dfens <- dfens[dfens$ensembl_gene_id %in% eids,]
dfens <- dfens[!is.na(dfens$ensembl_gene_id),]
dim(dfens)
# [1] 18101     7
dfens <- dfens[dfens$ensembl_gene_id %in% gsub("\\..*","",rownames(counts.se)),]
dfens <- dfens[!duplicated(dfens$ensembl_gene_id),]
colnames(dfens) <- c("seqnames","start","end","strand","hgnc_id","hgnc_symbol","ensembl_gene_id")
save(dfens, file="./data/dfens_v95hg38_bmart.rda")

# se experiments using filtered genes
counts.se <- counts.filt[gsub("\\..*","",rownames(counts.filt)) %in% dfens$ensembl_gene_id,]
counts.se <- counts.se[order(match(gsub("\\..*","",rownames(counts.se)), dfens$ensembl_gene_id)),]
identical(as.character(gsub("\\..*","",rownames(counts.se))), 
          as.character(dfens$ensembl_gene_id))
cpm.se <- cpm[gsub("\\..*","",rownames(cpm)) %in% dfens$ensembl_gene_id,]
cpm.se <- cpm.se[order(match(gsub("\\..*","",rownames(cpm.se)), dfens$ensembl_gene_id)),]
identical(gsub("\\..*","",rownames(cpm.se)), dfens$ensembl_gene_id)

# counts.se <- counts.filt[rownames(counts.filt) %in% counts.genes.grdf$countsdf.id,]
# cpm.se <- cpm[rownames(cpm) %in% counts.genes.grdf$countsdf.id,]
# dim(cpm.se)
# [1] 18166   145

# order genes for counts se
#ganno.counts <- counts.genes.grdf[order(match(counts.genes.grdf$countsdf.id,
#                                              rownames(counts.se))),]
#identical(counts.genes.grdf$countsdf.id, rownames(counts.se))
#ganno.tmm <- counts.genes.grdf[counts.genes.grdf$countsdf.id %in% rownames(cpm.se),]
#ganno.tmm <- ganno.tmm[order(match(ganno.tmm$countsdf.id,rownames(cpm.se))),]
#identical(ganno.tmm$countsdf.id,rownames(cpm.se))
#colnames(ganno.counts) <- colnames(ganno.tmm) <- c("gene.id","gene.symbol","countsdf.id","seqnames","start","end")
#ggr.counts <- makeGRangesFromDataFrame(ganno.counts, 
#                                       keep.extra.columns = T, 
#                                       ignore.strand = T)
#ggr.tmm <- makeGRangesFromDataFrame(ganno.tmm,
#                                    keep.extra.columns = T,
#                                    ignore.strand = T)

dfens$strand <- ifelse(dfens$strand==-1,"-","+")
ggr.counts <- makeGRangesFromDataFrame(dfens, 
                                       keep.extra.columns = T, 
                                       ignore.strand = F)
names(ggr.counts) <- rownames(counts.se)

ggr.tmm <- makeGRangesFromDataFrame(dfens,
                                    keep.extra.columns = T,
                                    ignore.strand = F)
names(ggr.tmm) <- rownames(cpm.se)

identical(gsub("\\..*","",names(ggr.counts)),ggr.counts$ensembl_gene_id)
identical(gsub("\\..*","",names(ggr.tmm)),ggr.tmm$ensembl_gene_id)

#====================================
# Make Summarized Experiment Objects
#====================================
# Gene Expr Counts SE object
identical(substr(colnames(counts.se),11,16), 
          substr(clinical.filt$TARGET.USI,11,16)) # TRUE
counts.seset <- SummarizedExperiment(assays = as.matrix(counts.se),
                                     rowRanges = ggr.counts, 
                                     colData = DataFrame(clinical.filt, 
                                                         row.names = colnames(counts.se)
                                     ),
                                     metadata = list(dataset = "TARGET_AML", 
                                                     assay_source = "GDC",
                                                     genome_build = "hg38")
)

# Gene TMM SE object
identical(ggr.tmm$ensembl_gene_id, gsub("\\..*","",rownames(cpm.se))) # TRUE
identical(substr(colnames(cpm.se),11,16), 
          substr(clinical.filt$TARGET.USI,11,16)) # TRUE
tmm.seset <- SummarizedExperiment(assays = as.matrix(cpm.se),
                                  rowRanges = ggr.tmm,
                                  colData = DataFrame(clinical.filt,
                                                      row.names = colnames(cpm.se)
                                  ),
                                  metadata = list(dataset = "TARGET_AML",
                                                  assay_source = "GDC",
                                                  genome_build = "hg38",
                                                  normalization_strategy = "TMM, log_cpm, limma, edgeR"))

# DEG TMM SE object
deglist = gsub("\\..*","",rownames(degtable))
length(intersect(deglist, dfens$ensembl_gene_id)) # 1984 of 1998
degfilt = deglist[deglist %in% dfens$ensembl_gene_id]
ggr.deg = ggr.counts[ggr.counts$ensembl_gene_id %in% degfilt]
# deg.assay <- counts.se[rownames(counts.se) %in% degfilt,]
deg.assay <- cpm.se[gsub("\\..*","",rownames(cpm.se)) %in% degfilt,]
deg.assay <- deg.assay[order(match(gsub("\\..*","",rownames(deg.assay)), degfilt)),]
identical(gsub("\\..*","",rownames(deg.assay)), degfilt)
ggr.deg <- ggr.deg[order(match(names(ggr.deg), rownames(deg.assay)))]
identical(names(ggr.deg), rownames(deg.assay)) # TRUE
identical(substr(colnames(cpm.se),11,16), 
          substr(clinical.filt$TARGET.USI,11,16)) # TRUE
# add the deg statistics to gene annotation
degstats = degtable[gsub("\\..*","",rownames(degtable)) %in% degfilt,]
degstats = degstats[order(match(gsub("\\..*","",rownames(degstats)), names(ggr.deg))),]
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
                                                  genome_build = "hg38",
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