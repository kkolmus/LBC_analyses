suppressPackageStartupMessages({
  library(futile.logger)
  library(GEOquery)
  library(R.utils)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(rlist)
})

flog.threshold(DEBUG)

flog.debug("Set directories and analysis parameters")

GEOid = "GSE120268"
input.dir = "~/Desktop/HT projects/LCB_analyses/DZB/GSE120268/data/input"
output.dir = "~/Desktop/HT projects/LCB_analyses/DZB/GSE120268/data/output"
UP = 0.6
DOWN = -0.6
pval = 0.05

flog.debug("Get metadata")
gse <- GEOquery::getGEO(GEO = GEOid, destdir = input.dir, GSEMatrix = TRUE)
metadata = list.files(path = input.dir, pattern = ".gz", full.names = TRUE)
R.utils::gunzip(metadata, remove = TRUE)
metadata = list.files(path = input.dir, pattern = ".txt", full.names = TRUE)

flog.debug("Get count data")
GEOquery::getGEOSuppFiles(GEO = GEOid, baseDir = input.dir)
data = list.files(path = file.path(input.dir, GEOid), pattern = ".gz", full.names = TRUE)
R.utils::gunzip(data, remove = TRUE)
data = list.files(path = file.path(input.dir, GEOid), pattern = ".txt", full.names = TRUE)

data <- read.table(file = data, header = TRUE, fill = TRUE)
metadata <- pData(phenoData(gse[[1]]))


flog.debug("Select upregulated genes")
genes_UP <- filter(data, data$log2FoldChange >= UP, data$padj < pval)
genes_DOWN <- filter(data, data$log2FoldChange <= DOWN, data$padj < pval)
  
id_genes_UP <- as.character(genes_UP[[1]])
id_genes_DOWN <- as.character(genes_DOWN[[1]])
id_genes_background <- as.character(data$ensembl_gene_id)



### Function for converting gene names between
convert.IDs <- function(dataset, OrgDb = "org.Hs.eg.db",
                        fromType = "ENSEMBL", toType = "ENTREZID") {
  bitr(geneID = dataset, fromType = fromType, toType = toType,  OrgDb = OrgDb)
}

id_genes_UP_converted <- convert.IDs(id_genes_UP)[[2]]
id_genes_DOWN_converted <- convert.IDs(id_genes_DOWN)[[2]]
id_genes_background_converted <- convert.IDs(id_genes_background)[[2]]


flog.debug("GO analysis of BP")
### Functions for GO analysis of Biological processes
GOanalyzer <- function(dataset, background,
                       ont = "BP", pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05, qvalueCutoff  = 0.05,
                       minGSSize = 10, maxGSSize = 500) {
  enrichGO(gene = dataset, universe = background, OrgDb = org.Hs.eg.db,
           ont = ont, pAdjustMethod = pAdjustMethod,
           pvalueCutoff  = pvalueCutoff, qvalueCutoff  = qvalueCutoff,
           minGSSize = minGSSize, maxGSSize = maxGSSize, readable = TRUE)
}

GOsimplifier <- function(dataset, cutoff, by = "p.adjust") {
  clusterProfiler::simplify(x = dataset, cutoff = cutoff, by = by)
}

# first analysis
BP_UP <- GOanalyzer(dataset = id_genes_UP_converted, background = id_genes_background_converted)
# BP_UP <- as.data.frame(BP_UP@result)
# BP_UP <- filter(BP_UP, BP_UP$Count > 10)

BP_DOWN <- GOanalyzer(dataset = id_genes_DOWN_converted, background = id_genes_background_converted)
# BP_DOWN <- as.data.frame(BP_DOWN@result)
# BP_DOWN <- filter(BP_DOWN, BP_DOWN$Count > 10)

BP <- list(BP_UP, BP_DOWN)
names(BP) <- c("UP", "DOWN")


# tuning analysis parameters
BP.up.sim.cutoff <- list()
BP.sim.condition <- NULL

# for loop
for (i in seq(from = 0.4, to = 0.70, by = 0.05)) {
  print(i)
  for (ds in names(BP)) {
    print(ds)
    if (nrow(as.data.frame(BP[[paste0(ds)]])) == 0) {
      next}
    else {
      temp = as.data.frame(GOsimplifier(dataset = BP[[paste0(ds)]], cutoff = i))
      temp = mutate(temp, Type = factor(paste0(ds)))
      a = as.numeric(gsub("[0-9]+[:/:]", "", temp$GeneRatio))
      temp = mutate(temp, GeneRatio = Count/a)
      temp = filter(temp, Count >= 10) }
    BP.sim.condition = rbind(BP.sim.condition, temp)
    BP.sim.condition = BP.sim.condition[!duplicated(BP.sim.condition), ]
  }
  BP.sim.cutoff = list.append(BP.sim.cutoff, BP.sim.condition)
}

names(BP.sim.cutoff) <- paste0("cutoff_", seq(from = 0.40, to = 0.70, by = 0.05))


flog.debug("Analysis of signaling pathways")
### Functions for analysis of signaling pathway

PathwayAnalyzer <- function(dataset, background,
                            pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05,
                            minGSSize = 10, maxGSSize = 500, organism = "human") {
  enrichPathway(gene = dataset, 
                universe = background, 
                organism = organism,
                pAdjustMethod = pAdjustMethod,
                pvalueCutoff  = pvalueCutoff, 
                qvalueCutoff  = qvalueCutoff,
                minGSSize = minGSSize, 
                maxGSSize = maxGSSize, 
                readable = TRUE)
}

pathways_UP <- PathwayAnalyzer(dataset = id_genes_UP_converted, background = id_genes_background_converted)
pathways_UP <- as.data.frame(pathways_UP@result)

pathways_DOWN <- PathwayAnalyzer(dataset = id_genes_DOWN_converted, background = id_genes_background_converted)
pathways_DOWN <- as.data.frame(pathways_DOWN@result)

sessionInfo()