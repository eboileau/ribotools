#!/biosw/R/3.5.1/bin/Rscript
# #!/usr/bin/Rscript

## Usage: ./run-degRibo-from-htseq.R [mm/rn/hs] [config] [num] [denom] [dirloc]

## mm = mouse, rn = rat, hs = human
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


run_analysis_fit_sub_lrt <- function(contrast, coldata, dds, model) {

    # fit one given contrast
    # NOTE: LFC shrinking is commented out

    # subset the DESeqDataSet for the selected contrast
    num.condition <- as.character(contrast[1])
    denom.condition <- as.character(contrast[2]) # reference

    coldata.sub <- coldata %>% dplyr::filter(condition %in% c(num.condition, denom.condition))
    dds.sub <- dds[,colnames(dds) %in% coldata.sub$sampleName]
    dds.sub$condition <- relevel(dds.sub$condition, ref=denom.condition)
    dds.sub$condition <- droplevels(dds.sub$condition)

    stopifnot(all(rownames(colData(dds.sub)) == colnames(dds.sub)))

    # fit the model
    dds.sub <- DESeq(dds.sub, test="LRT", reduced=model)

    coefs <- tail(resultsNames(dds.sub), 2)
    print(resultsNames(dds.sub))
    print(paste("Testing for: ", coefs, sep=""))
    print(paste("@ lfcThreshold: ", lfcThreshold.set, sep=""))
    print(paste("@ alpha: ", alpha.set, sep=""))


    ## main effect of condition on rna counts (interaction terms for base levels are absorbed by the intercept)
    res.rna <- results(dds.sub,
                       name=resultsNames(dds.sub)[3], #contrast=c("condition", num, denom),
                       lfcThreshold=lfcThreshold.set,
                       altHypothesis=altHypothesis.set,
                       alpha=alpha.set,
                       filterFun=ihw,
                       filter=rowMedians(counts(dds.sub)), # comment for default independent filtering (mean)
                       test="Wald")
    res.rna$padj[is.na(res.rna$padj)] <- 1

    ## effect of condition on ribo counts: condition_Trt_vs_Ctrl + assayribo.conditionTrt
    res.ribo <- results(dds.sub,
                        contrast=list(c(resultsNames(dds.sub)[3],resultsNames(dds.sub)[4])),
                        lfcThreshold=lfcThreshold.set,
                        altHypothesis=altHypothesis.set,
                        alpha=alpha.set,
                        filterFun=ihw,
                        filter=rowMedians(counts(dds.sub)),
                        test="Wald")
    res.ribo$padj[is.na(res.ribo$padj)] <- 1

    ## interaction: condition effect across assays, i.e. how different is the response
    ## from ribo vs. rna counts with treatment: (ribo/rna)_Trt / (ribo/rna)_Ctrl = (ribo_Trt/ribo_Ctrl) / (rna_Trt/rna_Ctrl)

    ## we cannot set lfc threshold
    res.inter <- results(dds.sub,
                        name=resultsNames(dds.sub)[4],
                        alpha=alpha.set,
                        filterFun=ihw,
                        filter=rowMedians(counts(dds.sub)))
    res.inter$padj[is.na(res.inter$padj)] <- 1

    ## assay ribo vs rna effect under control
    # use: contrast=c("assay", "ribo", "rna"), i.e resultsNames(dds.sub)[2] or assay_ribo_vs_rna, test="Wald"

    ## assay ribo vs rna effect for treatment (condition)
    # use: contrast=list(c(resultsNames(dds.sub)[2],resultsNames(dds.sub)[4])), test="Wald", i.e.
    # assay_ribo_vs_rna + assayribo.conditionTrt

    ## use lfcShrink for downstream analysis
    ## use "ashr" instead of "apeglm", which only takes coef (unless we rewrite the design...)
    ## with "ashr", we cannot use lfcThreshold
    ## write s-values, but only use p-values

    # res.rna.shrunken <- lfcShrink(dds.sub,
    #                             coef=resultsNames(dds.sub)[3], #contrast=c("condition", num, denom)
    #                             res=res.rna,
    #                             type="ashr",
    #                             svalue=TRUE)
    # res.rna.shrunken$padj <- res.rna$padj
    # res.rna.shrunken$pvalue <- res.rna$pvalue
    #
    #  res.ribo.shrunken <- lfcShrink(dds.sub,
    #                             contrast=list(c(resultsNames(dds.sub)[3],resultsNames(dds.sub)[4])),
    #                             res=res.ribo,
    #                             type="ashr",
    #                             svalue=TRUE)
    # res.ribo.shrunken$padj <- res.ribo$padj
    # res.ribo.shrunken$pvalue <- res.ribo$pvalue
    #
    #  res.inter.shrunken <- lfcShrink(dds.sub,
    #                                 coef=resultsNames(dds.sub)[4],
    #                                 res=res.inter,
    #                                 type="ashr",
    #                                 svalue=TRUE)
    # res.inter.shrunken$padj <- res.inter$padj
    # res.inter.shrunken$pvalue <- res.inter$pvalue

    ## write to xlsx
    # all.res <- list("dte"=res.inter,
    #                 "sdte"=res.inter.shrunken,
    #                 "rna"=res.rna,
    #                 "srna"=res.rna.shrunken,
    #                 "ribo"=res.ribo,
    #                 "sribo"=res.ribo.shrunken)
    all.res <- list("dte"=res.inter,
                    "rna"=res.rna,
                    "ribo"=res.ribo)
    write_results(dds.sub, all.res, num.condition, denom.condition)
}


run_analysis_fit_sub_deltaTE <- function(contrast, coldata, dds, model) {

    # use the deltaTE method
    # NOTE: in the orginal method, they do not apply alpha, filterFun
    # and use default filter when calling results
    # if there are more than 2 conditions, they fit the whole data, and
    # extract the relevant contrast in results, here we subset the object
    # and fit one contrast at a time

    # subset the DESeqDataSet for the selected contrast
    num.condition <- as.character(contrast[1])
    denom.condition <- as.character(contrast[2]) # reference

    coldata.sub <- coldata %>% dplyr::filter(condition %in% c(num.condition, denom.condition))
    dds.sub <- dds[,colnames(dds) %in% coldata.sub$sampleName]
    dds.sub$condition <- relevel(dds.sub$condition, ref=denom.condition)
    dds.sub$condition <- droplevels(dds.sub$condition)

    stopifnot(all(rownames(colData(dds.sub)) == colnames(dds.sub)))

    # fit the model
    dds.sub <- DESeq(dds.sub)

    coefs <- tail(resultsNames(dds.sub), 2)
    print(resultsNames(dds.sub))
    print(paste("Testing for: ", coefs, sep=""))
    print(paste("@ lfcThreshold: ", lfcThreshold.set, sep=""))
    print(paste("@ alpha: ", alpha.set, sep=""))

    ## assayribo.conditionTrt = changes in Ribo in conditon Trt vs Ctrl accounting for changes in RNA in Trt vs Ctrl
    res.inter <- results(dds.sub,
                        name=resultsNames(dds.sub)[4],
                        alpha=alpha.set,
                        filterFun=ihw,
                        filter=rowMedians(counts(dds.sub)))
    res.inter$padj[is.na(res.inter$padj)] <- 1

    # recreate the object to fit ribo and rna separately w/o the assay variable
    idx.ribo <- which(coldata.sub$assay == "ribo")
    coldata.ribo <- coldata.sub[idx.ribo,]
    coldata.ribo$assay <- NULL
    idx.rna <- which(coldata.sub$assay == "rna")
    coldata.rna <- coldata.sub[idx.rna,]
    coldata.rna$assay <- NULL

    cts <- counts(dds.sub)
    cts.ribo <- cts[,idx.ribo]
    cts.rna <- cts[,idx.rna]

    dds.ribo <- DESeqDataSetFromMatrix(countData=cts.ribo,
                                       colData=coldata.ribo,
                                       design=model)
    dds.ribo <- DESeq(dds.ribo)
    res.ribo <- results(dds.ribo,
                        contrast=c("condition", num.condition, denom.condition),
                        lfcThreshold=lfcThreshold.set,
                        altHypothesis=altHypothesis.set,
                        alpha=alpha.set,
                        filterFun=ihw,
                        filter=rowMedians(counts(dds.ribo)))
    res.ribo$padj[is.na(res.ribo$padj)] <- 1
    res.ribo <- lfcShrink(dds.ribo, coef=2, type="apeglm", res=res.ribo)

    dds.rna <- DESeqDataSetFromMatrix(countData=cts.rna,
                                       colData=coldata.rna,
                                       design=model)
    dds.rna <- DESeq(dds.rna)
    res.rna <- results(dds.rna,
                       contrast=c("condition", num.condition, denom.condition),
                       lfcThreshold=lfcThreshold.set,
                       altHypothesis=altHypothesis.set,
                       alpha=alpha.set,
                       filterFun=ihw,
                       filter=rowMedians(counts(dds.rna)), # comment for default independent filtering (mean)
                       test="Wald")
    res.rna$padj[is.na(res.rna$padj)] <- 1
    res.rna <- lfcShrink(dds.rna, coef=2, type="apeglm", res=res.rna)

    all.res <- list("dte"=res.inter,
                    "rna"=res.rna,
                    "ribo"=res.ribo)
    write_results(dds.sub, all.res, num.condition, denom.condition)
}


write_results <- function(dds, all.res, num, denom) {

    # NOTE: LFC shrinking is commented out

    # first get baseMean, any
    means <- all.res[["ribo"]] %>%
        data.frame() %>%
        dplyr::select(baseMean) %>%
        rownames_to_column(var="gene") %>%
        as_tibble()

    all.res[["dte"]] <- all.res[["dte"]] %>%
        data.frame() %>%
        dplyr::select(log2FoldChange, pvalue, padj) %>%
        dplyr::rename(log2FC.dte = log2FoldChange,
            pvalue.dte = pvalue,
            padj.dte = padj) %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
    # all.res[["sinter"]] <- all.res[["sinter"]] %>%
    #     data.frame() %>%
    #     dplyr::select(log2FoldChange, svalue) %>%
    #     dplyr::rename(shrunken.log2FC.dte = log2FoldChange,
    #         svalue.dte = svalue) %>%
    #     rownames_to_column(var="gene") %>%
    #     as_tibble()

    all.res[["rna"]] <- all.res[["rna"]] %>%
        data.frame() %>%
        dplyr::select(log2FoldChange, pvalue, padj) %>%
        dplyr::rename(log2FC.rna = log2FoldChange,
            pvalue.rna = pvalue,
            padj.rna = padj) %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
    # all.res[["srna"]] <- all.res[["srna"]] %>%
    #     data.frame() %>%
    #     dplyr::select(log2FoldChange, svalue) %>%
    #     dplyr::rename(shrunken.log2FC.rna = log2FoldChange,
    #         svalue.rna = svalue) %>%
    #     rownames_to_column(var="gene") %>%
    #     as_tibble()

    all.res[["ribo"]] <- all.res[["ribo"]] %>%
        data.frame() %>%
        dplyr::select(log2FoldChange, pvalue, padj) %>%
        dplyr::rename(log2FC.ribo = log2FoldChange,
            pvalue.ribo = pvalue,
            padj.ribo = padj) %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
    # all.res[["sribo"]] <- all.res[["sribo"]] %>%
    #     data.frame() %>%
    #     dplyr::select(log2FoldChange, svalue) %>%
    #     dplyr::rename(shrunken.log2FC.ribo = log2FoldChange,
    #         svalue.ribo = svalue) %>%
    #     rownames_to_column(var="gene") %>%
    #     as_tibble()

    # add raw counts, annotations and re-order columns
    cts <- counts(dds, normalized=FALSE) %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()

    res <- all.res[["dte"]] %>%
    #     full_join(all.res[["sinter"]], by="gene") %>%
        full_join(all.res[["rna"]], by="gene") %>%
    #     full_join(all.res[["srna"]], by="gene") %>%
        full_join(all.res[["ribo"]], by="gene") %>%
    #     full_join(all.res[["sribo"]], by="gene") %>%
        left_join(means, by="gene")

    res$symbol <- id.mapping

    res <- res %>%
            dplyr::select(gene, symbol, baseMean,
                    log2FC.dte, pvalue.dte, padj.dte,
    #                shrunken.log2FC.dte, svalue.dte,
                    log2FC.rna, pvalue.rna, padj.rna,
    #                 shrunken.log2FC.rna, svalue.rna,
                    log2FC.ribo, pvalue.ribo, padj.ribo)
    #                  shrunken.log2FC.ribo, svalue.ribo)


    res <- res %>% left_join(cts, by="gene")

    ## write to disk
    wb <- createWorkbook()

    # if sheetName is too long this won't work...
    sheetName <- paste0(num, "_vs_", denom, sep="")
    addWorksheet(wb, sheetName=sheetName)
    writeDataTable(wb, sheet=1, x=res)

    filen <- paste0("contrast_", num, "_vs_", denom, sep="")
    dirloc.sub <- file.path(dirloc.out, filen, fsep=.Platform$file.sep)
    created <- ifelse(!dir.exists(dirloc.sub), dir.create(dirloc.sub, recursive=TRUE), FALSE)
    if (!created) { print(paste("Using existing directory ", dirloc.sub, " ...", sep="")) }
    filen <- paste0(filen, ".xlsx", sep="")
    filen <- file.path(dirloc.sub, filen, fsep=.Platform$file.sep)
    print(paste("Writing results to: ", filen, sep=""))
    saveWorkbook(wb, filen, overwrite = TRUE)
    get_reg_layers(res, dirloc.sub)
}


read_symbols <- function(sampleTable, dirloc) {
    # read gene symbols from HTSeq output tables
    l <- lapply( as.character( sampleTable[,2] ), function(fn) read.table( file.path( dirloc, fn ), fill=TRUE ))
    if( ! all( sapply( l, function(a) all( head(a$V2, -5) == head(l[[1]]$V2, -5) ) ) ) )
        stop( "Gene names (second column) differ between files." )
    tbl <- sapply( l, function(a) a[,2] )
    colnames(tbl) <- sampleTable[,1]
    rownames(tbl) <- l[[1]]$V1
    oldSpecialNames <- c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")
    specialRows <- (substr(rownames(tbl),1,1) == "_") | rownames(tbl) %in% oldSpecialNames
    tbl <- tbl[ !specialRows, , drop=FALSE ]
    # in fact we just need the first columnn - format as mapIds
    tbl[,1]
}

get_reg_layers <- function(res, dirloc) {
    # use same terminology as Chothani et al. (deltaTE)
    # use same global alpha.set

    # Translationally forwarded genes have a significant change in mRNA and RPF at the same rate, with no significant change in TE
    forwarded <- res[which(res$padj.dte > alpha.set & res$padj.ribo < alpha.set & res$padj.rna < alpha.set),]
    write.table(forwarded[,1:2], file.path(dirloc, "forwarded.txt", fsep=.Platform$file.sep), quote=F, sep="\t", col.names = F, row.names = F)

    # Exclusive genes are DTEGs that have a significant change in RPF, with no change in mRNA leading to a significant change in TE
    exclusive <- res[which(res$padj.dte < alpha.set & res$padj.ribo < alpha.set & res$padj.rna > alpha.set),]
    write.table(exclusive[,1:2], file.path(dirloc, "exclusive.txt", fsep=.Platform$file.sep), quote=F, sep="\t", col.names = F, row.names = F)

    # direction of change is taken into account
    both < which(res$padj.dte < alpha.set & res$padj.ribo < alpha.set & res$padj.rna < alpha.set)
    # Translationally intensified genes have a significant change in TE that acts with the effect of transcription
    intensified <- res[both[which(res[both,]$log2FC.dte*res[both,]$log2FC.rna > 0)],]
    write.table(intensified[,1:2], file.path(dirloc, "intensified.txt", fsep=.Platform$file.sep), quote=F, sep="\t", col.names = F, row.names = F)

    # Translationally buffered genes have a significant change in TE that counteracts the change in RNA (buffering transcription)
    buffered <- res[both[which(res[both,]$log2FC.dte*res[both,]$log2FC.rna < 0)],]
    buffered <- rbind(buffered, res[which(res$padj.dte < alpha.set & res$padj.ribo > alpha.set & res$padj.rna < alpha.set),])
    write.table(buffered[,1:2], file.path(dirloc, "buffered.txt", fsep=.Platform$file.sep), quote=F, sep="\t", col.names = F, row.names = F)

    filen <- file.path(dirloc, "scatter.pdf", fsep=.Platform$file.sep)
    pdf(filen, useDingbats = F)
    max_val = max(res$log2FC.ribo,res$log2FC.rna, na.rm = T)
    plot(y=res$log2FC.ribo, x=res$log2FC.rna, xlab="RNA-seq log2FC",ylab = "Ribo-seq log2FC", asp=1, pch=16, col=rgb(128/255,128/255,128/255,0.1), ylim=c(-max_val,max_val), xlim=c(-max_val,max_val), cex=0.4)
    abline(a=0, b=1, col="gray")
    abline(h=0, v=0, col="gray")
    points(y=res[res$gene %in% forwarded$gene,]$log2FC.ribo, x=res[res$gene %in% forwarded$gene,]$log2FC.rna, pch=16, col=rgb(0,0,1,1))
    points(y=res[res$gene %in% exclusive$gene,]$log2FC.ribo, x=res[res$gene %in% exclusive$gene,]$log2FC.rna, pch=16, col=rgb(0,1,0,1))
    points(y=res[res$gene %in% intensified$gene,]$log2FC.ribo, x=res[res$gene %in% intensified$gene,]$log2FC.rna, pch=16, col=rgb(1,0,0,1))
    points(y=res[res$gene %in% buffered$gene,]$log2FC.ribo, x=res[res$gene %in% buffered$gene,]$log2FC.rna, pch=16, col=rgb(1,0,0,1))
    dev.off()
}

# ---------------------------------------------------------
## Call

# defaults
base.loc <- "count-tables"

# arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<4) { stop("./run-tea.R [YAML] [CONTRASTS] [OUTPUT] [METHOD] <BATCH> <FILTER>\n", call.=FALSE) }
# config file
params <- yaml::read_yaml(args[1])
# contrasts file
contrasts <- read.table(args[2], header=T)
# output directory
dirloc.out <- args[3]
created <- ifelse(!dir.exists(dirloc.out), dir.create(dirloc.out, recursive=TRUE), FALSE)
if (!created) { print(paste("Using existing directory ", dirloc.out, " ...", sep="")) }
# method 0=LRT, 1=deltaTE
method <- "LRT"
if (as.integer(args[4])==1) {
    method <- "deltaTE"
}
print(paste("Using the following method ", method, " ...", sep=""))
# batch effect
batch <- NULL
if (length(args)>4) {
    batch <- read.table(args[5], header=T)
    colnames(batch) <- c("sampleName", "batch")
}
# filter 0 counts 0=No 1=Yes
filter <- FALSE
if (length(args)>5 & as.integer(args[6])==1) {
    filter <- TRUE
}


# write tmp files to
tmp.loc <- uuid::UUIDgenerate()
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

## construct sample table
# rna
rna.table <- params$rnaseq_sample_name_map %>%
    data.frame() %>% t %>% data.frame(stringsAsFactors=FALSE) %>%
    dplyr::rename(sampleName = ".")

rna.table <- rna.table %>%
    rowwise() %>% mutate(fileName=rna.files[grep(sampleName, rna.files)])

rna.table$assay <- 'rna'
rna.table$condition <- NA
conditions <- unique(as.vector(as.matrix(contrasts)))
used <- lapply(conditions, function(c) {rna.table$condition[grep(c, rna.table$sampleName, fixed=TRUE)] <<- c})

# ribo
ribo.table <- params$riboseq_sample_name_map %>%
    data.frame() %>% t %>% data.frame(stringsAsFactors=FALSE) %>%
    dplyr::rename(sampleName = ".")

ribo.table <- ribo.table %>%
    rowwise() %>% mutate(fileName=ribo.files[grep(sampleName, ribo.files)])

ribo.table$assay <- 'ribo'
ribo.table$condition <- NA
used <- lapply(conditions, function(c) {ribo.table$condition[grep(c, ribo.table$sampleName, fixed=TRUE)] <<- c})

sampleTable <- bind_rows(ribo.table, rna.table)
sampleTable <- as.data.frame(sampleTable) # need dataframe for DESEq
sampleTable$condition <- factor(sampleTable$condition)
sampleTable$assay <- factor(sampleTable$assay)

if (!is.null(batch)) {
    sampleTable <- cbind(sampleTable, batch[match(sampleTable$sampleName, batch$sampleName),]$batch)
    colnames(sampleTable)[ncol(sampleTable)] <- 'batch'
}

# write for reference
write.csv(sampleTable, file=file.path(dirloc.out, "sampleTable.csv", fsep=.Platform$file.sep))

# no filtering - should be the same as above
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
model <- ~assay+condition+assay:condition
if (method == "LRT") {
  model.reduced <- ~assay+condition
} else if (method == "deltaTE") {
    model.reduced <- ~condition
    if (!is.null(batch)) { model.reduced <- ~batch+condition }
}
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                       directory=dirloc.tmp,
                                       design=model)

stopifnot(all(sampleTable$sampleName == colnames(ddsHTSeq)))

# always rna as reference level for assay
ddsHTSeq$assay <- relevel(ddsHTSeq$assay, ref="rna")

## annotation
# the order should be the same...
id.mapping <- read_symbols(sampleTable, dirloc.tmp)

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

if (filter) {
    print("Removing genes with 0 counts in each assay separately, using the intersection...")

    ddsHTSeq.cts <- counts(ddsHTSeq)

    m_ribo <- ddsHTSeq.cts[,which(colData((ddsHTSeq))$assay == "ribo")]
    m_rna <- ddsHTSeq.cts[,which(colData((ddsHTSeq))$assay == "rna")]

    idx.ribo <- which(apply(m_ribo, 1, function(x){length(which(x>0))})==ncol(m_ribo))
    idx.rna <- which(apply(m_rna, 1, function(x){length(which(x>0))})==ncol(m_rna))

    idx <- intersect(idx.ribo, idx.rna)
    ddsHTSeq <- ddsHTSeq[idx,]
}

## fit contrasts
if (method == "LRT") {
    apply(contrasts, 1, run_analysis_fit_sub_lrt, coldata=sampleTable, dds=ddsHTSeq, model=model.reduced)
} else if (method == "deltaTE") {
    apply(contrasts, 1, run_analysis_fit_sub_deltaTE, coldata=sampleTable, dds=ddsHTSeq, model=model.reduced)
}
