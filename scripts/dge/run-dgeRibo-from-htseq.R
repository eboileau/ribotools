#!/biosw/R/3.5.1/bin/Rscript
# #!/usr/bin/Rscript

## Usage: ./run-degRibo-from-htseq.R [1/2/3] [config] [num] [denom] [dirloc]

## 1/2/3: 1 mouse, 2 human, 3 rat
## config: same configuration file used when calling run-htseq-workflow
## num:    condition level to compare, e.g. treatment
## denom:  reference condition level, e.g. control
## dirloc: output directory (full path), also used as tmp/

## NOTE: num, denom must match the names used by run-htseq-workflow
##       in constructing the file names, given by ribo_ and rnaseq_sample_name_map

##       We assume that gene ids are ensembl ids ** TODO generalise, pass keytype as args

# ---------------------------------------------------------

## Ribo- and RNA-seq are both filtered independently and ONLY genes where
## both have non-zero values for ALL replicates are kept.  
## Joint normalisation is performed within DESeq, see below for separate normalisation.
## ** TODO test separate normalisation, relax filter (non-zero ALL), to allow 
##         Ribo- and/or RNA-"specific" genes

## Independent filtering using row median (mean is default)

## Ribo- and RNA-seq LFC are determined using pairwise comparisons (Wald test).
## We set a LFC threshold for the results, however we cannot do this for the shrunken FC
## values. For the latter, we filter the results using the "un-shrunken p-values" and LFC.

## Ratio of ratios (interaction term) is used to determined if genes are regulated differently
## in ribo vs. rna: (ribo/rna)_Trt / (ribo/rna)_Ctrl = (ribo_Trt/ribo_Ctrl) / (rna_Trt/rna_Ctrl).
## In the latter, we use LRT to compare the full model to the reduced model to identify significant genes.
## The p-values are determined by the difference in deviance between the full and reduced model and 
## not the LFC. Thus for this term, we only filter based on the "un-shrunken p-values".
## There are no LFC threshold set. 

## Shrinking done using "ashr".

# ---------------------------------------------------------

## Default LFC and FDR threshold

lfcThreshold.set <- log2(1.2)
altHypothesis.set <- "greaterAbs"
alpha.set <- 0.05

# ---------------------------------------------------------

library(DESeq2)
library(IHW)
library(ashr)

library("Glimma")
library("genefilter")

library("org.Mm.eg.db")
library("org.Hs.eg.db")
library("org.Rn.eg.db")

library(yaml)

library(dplyr)
library(tibble)
library(purrr)

library(openxlsx)

# ---------------------------------------------------------
## Functions

write_results <- function(dds, inter, ribo, rna, num, denom, shrunken, genome) {

    ## filter based on the "un-shrunken p-values"
    
    if (shrunken == "shrunken") {
    
        ## interaction
        res.inter.tib <- inter %>%
            data.frame() %>%
            select(log2FoldChange, padj, svalue) %>%
            rename(log2FC.inter = log2FoldChange, 
                padj.inter = padj,
                svalue.inter = svalue) %>%
            rownames_to_column(var="gene") %>% 
            as_tibble()
   
        ## ribo
        res.ribo.tib <- ribo %>%
            data.frame() %>%
            select(log2FoldChange, padj, svalue) %>%
            rename(log2FC.ribo = log2FoldChange, 
                padj.ribo = padj,
                svalue.ribo = svalue) %>%
            rownames_to_column(var="gene") %>% 
            as_tibble()

        ## rna
        res.rna.tib <- rna %>%
            data.frame() %>%
            select(log2FoldChange, padj, svalue) %>%
            rename(log2FC.rna = log2FoldChange, 
                padj.rna = padj,
                svalue.rna = svalue) %>%
            rownames_to_column(var="gene") %>% 
            as_tibble()
  
    } else{
    
        ## interaction
        res.inter.tib <- inter %>%
            data.frame() %>%
            select(log2FoldChange, padj) %>%
            rename(log2FC.inter = log2FoldChange, 
                padj.inter = padj) %>%
            rownames_to_column(var="gene") %>% 
            as_tibble()

        ## ribo
        res.ribo.tib <- ribo %>%
            data.frame() %>%
            select(log2FoldChange, padj) %>%
            rename(log2FC.ribo = log2FoldChange, 
                padj.ribo = padj) %>%
            rownames_to_column(var="gene") %>% 
            as_tibble()

        ## rna
        res.rna.tib <- rna %>%
            data.frame() %>%
            select(log2FoldChange, padj) %>%
            rename(log2FC.rna = log2FoldChange, 
                padj.rna = padj) %>%
            rownames_to_column(var="gene") %>% 
            as_tibble()

    }
    
    ## add raw counts, annotations and re-order columns
    cts <- counts(dds, normalized=FALSE) %>%
        data.frame() %>%    
        rownames_to_column(var="gene") %>% 
        as_tibble()
    
    # # nromalised and/or cpm?
    # c <- counts(dds.simple, normalized=F)
    # c <- cpm(c)

    # any baseMean, they are all the same
    means <- ribo %>%
        data.frame() %>%
        select(baseMean) %>%
        rownames_to_column(var="gene") %>% 
        as_tibble()
        
    res <- res.inter.tib %>% 
        full_join(res.ribo.tib, by="gene") %>% 
        full_join(res.rna.tib, by="gene") %>% 
        left_join(means, by="gene")
        
    # <- Reduce(function(x,y) merge(x, y, by = "gene", all.x = TRUE, all.y = TRUE),
    # list(res.inter.tib, res.ribo.tib, res.rna.tib))
    # <- merge(tata, by = "gene", all.x = TRUE, all.y = FALSE)
    
    res <- res %>%
      filter((padj.inter < alpha.set) | (padj.ribo < alpha.set & abs(log2FC.ribo) > lfcThreshold.set) | (padj.rna < alpha.set & abs(log2FC.rna) > lfcThreshold.set))
    
    res$symbol <- map_ids(genome, res$gene)
    
    if (shrunken == "shrunken") {
    
        res <- res %>% 
            select(gene, symbol, baseMean, log2FC.inter, padj.inter, svalue.inter, log2FC.ribo, padj.ribo, svalue.ribo, log2FC.rna, padj.rna, svalue.rna)
            
    } else {
        res <- res %>% 
            select(gene, symbol, baseMean, log2FC.inter, padj.inter, log2FC.ribo, padj.ribo, log2FC.rna, padj.rna)
            
    }
    
    # add counts
    res <- res %>% left_join(cts, by="gene")
    
    ## write to disk, add size factors for reference
    wb <- createWorkbook()

    addWorksheet(wb, sheetName=paste0(num, "_vs_", denom, sep=""))
    writeDataTable(wb, sheet=1, x=res)

    addWorksheet(wb, sheetName="sizeFactors")
    
    sf <- dds$sizeFactor %>%
        data.frame() %>%    
        rownames_to_column(var="sample") %>% 
        rename(sizeFactor = ".") %>%
        as_tibble()
            
    writeDataTable(wb, sheet=2, x= sf)
    
    if (shrunken == "shrunken") {
        filen <- paste0("condition_", num, "_vs_", denom, "_shrunken.xlsx", sep="")
        filen <- file.path(dirloc.out, filen, fsep=.Platform$file.sep)
    } else {
        filen <- paste0("condition_", num, "_vs_", denom, ".xlsx", sep="")
        filen <- file.path(dirloc.out, filen, fsep=.Platform$file.sep)
    }
    saveWorkbook(wb, filen, overwrite = TRUE)

}


map_ids <- function(genome, keys, keytype="ENSEMBL", column="SYMBOL") {

    genome <- switch(genome, org.Mm.eg.db, org.Hs.eg.db, org.Rn.eg.db)
    mapIds(genome, keys=keys, column=column, keytype=keytype, multiVals="first")
}


# ---------------------------------------------------------
## Call

# config file
args <- commandArgs(trailingOnly=TRUE)

genome <- args[1]
if (!genome %in% c(1,2,3)) { stop("Genome 1=mouse, 2=human, 3=rat") }

params.file <- args[2]
params <- yaml::read_yaml(params.file)

num <- args[3]
denom <- args[4]
dirloc.out <- args[5]

# default from pipeline
base.loc <- "count-tables"

# write tmp files to
tmp.loc <- "tmp"
dirloc.tmp <- file.path(dirloc.out, tmp.loc, fsep=.Platform$file.sep)
dir.create(dirloc.tmp)

# first use separate locations, this is easier, then 
# link all files to base.loc, remove afterwards.

dirloc.ribo <- params$riboseq_data
dirloc.ribo <- file.path(dirloc.ribo, base.loc, fsep=.Platform$file.sep)
ribo.files <- list.files(dirloc.ribo)

dirloc.rna <- params$rnaseq_data
dirloc.rna <- file.path(dirloc.rna, base.loc, fsep=.Platform$file.sep)
rna.files <- list.files(dirloc.rna)

## construct sample table (coldata)

# Currently one contrast at a time, but yet construct full table.
# We could do all in one run, and subset the results using contrast...
## TODO do not call each contrast separately...

rna.table <- params$rnaseq_sample_name_map %>% 
    data.frame() %>% t %>% data.frame(stringsAsFactors=FALSE) %>% 
    rename(sampleName = ".")

rna.table <- rna.table %>% 
    rowwise() %>% mutate(fileName=rna.files[grep(sampleName, rna.files)])

rna.table$assay <- 'rna' 
rna.table$condition <- NA
rna.table$condition[grep(num, rna.table$sampleName, fixed=TRUE)] <- num
rna.table$condition[grep(denom, rna.table$sampleName, fixed=TRUE)] <- denom


ribo.table <- params$riboseq_sample_name_map %>% 
    data.frame() %>% t %>% data.frame(stringsAsFactors=FALSE) %>% 
    rename(sampleName = ".")

ribo.table <- ribo.table %>% 
    rowwise() %>% mutate(fileName=ribo.files[grep(sampleName, ribo.files)])

ribo.table$assay <- 'ribo' 
ribo.table$condition <- NA
ribo.table$condition[grep(num, ribo.table$sampleName, fixed=TRUE)] <- num
ribo.table$condition[grep(denom, ribo.table$sampleName, fixed=TRUE)] <- denom


sampleTable <- bind_rows(ribo.table, rna.table)
sampleTable <- as.data.frame(sampleTable) # need dataframe for DESEq
sampleTable$condition <- factor(sampleTable$condition)
sampleTable$assay <- factor(sampleTable$assay)

## reduce table to one contrast
sampleTable <- sampleTable %>% filter(!is.na(condition))

# write for reference
filen <- paste("condition_", num, "_vs_", denom, "_sampleTable.csv", sep="")
write.csv(sampleTable, file=file.path(dirloc.out, filen, fsep=.Platform$file.sep))

ribo.files <- sampleTable$fileName[sampleTable$assay=="ribo"]
rna.files <- sampleTable$fileName[sampleTable$assay=="rna"]

# link all files for DESeq
for (f in ribo.files) {
    from <- file.path(dirloc.ribo, f, fsep=.Platform$file.sep)
    to <- file.path(dirloc.tmp, f, fsep=.Platform$file.sep)
    file.symlink(from, to)
}

for (f in rna.files) {
    from <- file.path(dirloc.rna, f, fsep=.Platform$file.sep)
    to <- file.path(dirloc.tmp, f, fsep=.Platform$file.sep)
    file.symlink(from, to)
}

## create DESEq data object
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                       directory=dirloc.tmp,
                                       design=~assay+condition+assay:condition)

stopifnot(all(sampleTable$sampleName == colnames(ddsHTSeq)))

# always rna as reference level for assay
ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref=denom)
ddsHTSeq$assay <- relevel(ddsHTSeq$assay, ref="rna")

# remove files
for (f in list.files(dirloc.tmp)) {
    to <- file.path(dirloc.tmp, f, fsep=.Platform$file.sep)
    file.remove(to)
}

# separate normalisation ???
# sf <- numeric(ncols(ddsHTSeq))
# idx.ribo <- ddsHTSeq$assay == 'ribo'
# sf[idx.ribo] <- estimateSizeFactorsFromMatrix(counts(ddsHTSeq)[,idx.ribo])
# idx.rna <- ddsHTSeq$assay == 'rna'
# sf[idx.rna] <- estimateSizeFactorsFromMatrix(counts(ddsHTSeq)[,idx.rna])
# sizeFactors(ddsHTSeq) <- sf

## filter

ddsHTSeq.cts <- counts(ddsHTSeq)

m_ribo <- ddsHTSeq.cts[,which(colData((ddsHTSeq))$assay == "ribo")]
m_rna <- ddsHTSeq.cts[,which(colData((ddsHTSeq))$assay == "rna")]

idx.ribo <- which(apply(m_ribo, 1, function(x){length(which(x>0))})==ncol(m_ribo))
idx.rna <- which(apply(m_rna, 1, function(x){length(which(x>0))})==ncol(m_rna))

idx <- intersect(idx.ribo, idx.rna)
ddsHTSeq <- ddsHTSeq[idx,]


## fit the model
ddsHTSeq.reduced <- DESeq(ddsHTSeq, test="LRT", reduced=~assay+condition)

## main effect of condition on rna counts (interaction terms for base levels are absorbed by the intercept)
    
res.rna <- results(ddsHTSeq.reduced,
                   contrast=c("condition", num, denom),
                   lfcThreshold=lfcThreshold.set, 
                   altHypothesis=altHypothesis.set,
                   alpha=alpha.set,
                   filterFun=ihw,
                   filter=rowMedians(counts(ddsHTSeq.reduced)), # comment for default independent filtering (mean)
                   test="Wald")
                        
res.rna$padj[is.na(res.rna$padj)] <- 1

## effect of condition on ribo counts: condition_Trt_vs_Ctrl + assayribo.conditionTrt
                
res.ribo <- results(ddsHTSeq.reduced,
                    contrast=list(c(resultsNames(ddsHTSeq.reduced)[3],resultsNames(ddsHTSeq.reduced)[4])),
                    lfcThreshold=lfcThreshold.set, 
                    altHypothesis=altHypothesis.set,
                    alpha=alpha.set,
                    filterFun=ihw,
                    filter=rowMedians(counts(ddsHTSeq.reduced)),
                    test="Wald")

res.ribo$padj[is.na(res.ribo$padj)] <- 1

## interaction: condition effect across assays, i.e. how different is the response 
## from ribo vs. rna counts with treatment: (ribo/rna)_Trt / (ribo/rna)_Ctrl = (ribo_Trt/ribo_Ctrl) / (rna_Trt/rna_Ctrl)

## we cannot set lfc threshold 
res.inter <- results(ddsHTSeq.reduced,
                     name=resultsNames(ddsHTSeq.reduced)[4],
                     alpha=alpha.set,
                     filterFun=ihw,
                     filter=rowMedians(counts(ddsHTSeq.reduced)))

res.inter$padj[is.na(res.inter$padj)] <- 1

## assay ribo vs rna effect under control
# use: contrast=c("assay", "ribo", "rna"), i.e resultsNames(ddsHTSeq.reduced)[2] or assay_ribo_vs_rna, test="Wald"

## assay ribo vs rna effect for treatment (condition)
# use: contrast=list(c(resultsNames(ddsHTSeq.reduced)[2],resultsNames(ddsHTSeq.reduced)[4])), test="Wald", i.e.
# assay_ribo_vs_rna + assayribo.conditionTrt

## use lfcShrink for downstream analysis
## use "ashr" instead of "apeglm", which only takes coef (unless we rewrite the design...)
## with "ashr", we cannot use lfcThreshold
## write s-values, but only use p-values

res.rna.shrunken <- lfcShrink(ddsHTSeq.reduced, 
                              contrast=c("condition", num, denom), # coef=resultsNames(ddsHTSeq.reduced)[2]
                              res=res.rna,
                              type="ashr",
                              svalue=TRUE)

res.rna.shrunken$padj <- res.rna$padj
res.rna.shrunken$pvalue <- res.rna$pvalue

res.ribo.shrunken <- lfcShrink(ddsHTSeq.reduced, 
                               contrast=list(c(resultsNames(ddsHTSeq.reduced)[3],resultsNames(ddsHTSeq.reduced)[4])), 
                               res=res.ribo,
                               type="ashr",
                               svalue=TRUE)
                                
res.ribo.shrunken$padj <- res.ribo$padj
res.ribo.shrunken$pvalue <- res.ribo$pvalue

res.inter.shrunken <- lfcShrink(ddsHTSeq.reduced, 
                                coef=resultsNames(ddsHTSeq.reduced)[4], 
                                res=res.inter,
                                type="ashr",
                                svalue=TRUE)

res.inter.shrunken$padj <- res.inter$padj
res.inter.shrunken$pvalue <- res.inter$pvalue

## Glimma, only for the interaction, but write all results to disk     

res.inter.shrunken$symbol <- map_ids(genome, rownames(res.inter.shrunken))
is.de <- as.numeric(res.inter.shrunken$padj < alpha.set)
anno <- data.frame(GeneID=rownames(res.inter.shrunken), symbol=res.inter.shrunken$symbol)

## log2FC on y-axis, log mean expression on x-axis [log(res.df$baseMean + 0.5)]
## side plot: average of the counts normalized by size factor on y-axis [counts(dds ,normalized=TRUE)]
## if transform = TRUE, as.matrix(edgeR::cpm(counts, log=TRUE))

glMDPlot(res.inter.shrunken, 
         counts=counts(ddsHTSeq.reduced ,normalized=TRUE),
         anno, # anno
         ddsHTSeq.reduced$condition, # groups
         samples=colnames(ddsHTSeq.reduced), 
         status=is.de, 
         transform = FALSE,
         xlab = "logMeanExpr",
         ylab = "log2FoldChange",
         side.ylab = "NormalizedCount",
         path=dirloc.out, 
         folder=paste("glimma-plots_", num, "_vs_", denom, "_", resultsNames(ddsHTSeq.reduced)[4], sep=""), 
         launch=FALSE)


## xlsx file: both shrunken and not
write_results(ddsHTSeq.reduced,
              res.inter, 
              res.ribo, 
              res.rna,
              num, 
              denom, 
              "",
              genome)
                
write_results(ddsHTSeq.reduced, 
              res.inter.shrunken, 
              res.ribo.shrunken, 
              res.rna.shrunken,
              num, 
              denom, 
              "shrunken",
              genome)
              
