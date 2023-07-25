#! /usr/bin/env Rscript

# Usage: ./pre-sample-table.R [-sample SAMPLE] [-contrast CONTRAST] [-output OUTPUT] <-method LRT/deltaTE> <--batch> <--filter>
# 1: [-sample SAMPLE]      Sample table (CSV) - or output from prep-sample-table.R
# 2. [-contrast CONTRAST]  CSV file with contrasts to test on each line in 2 columns w/o header - condition,reference
# 3. [-output OUTPUT]      Output directory
# 4. <-method LRT/deltaTE> Method to use (default LRT)
# 5. <--batch>             If present, uses "batch" from sample table
# 6. <--filter>            If present, filter genes


library(DESeq2)
library(IHW)
library(ashr)

library(dplyr)
library(tibble)
library(purrr)

library(openxlsx)

# ---------------------------------------------------------

## Default LFC and FDR threshold

lfcThreshold.set <- log2(1.2)
altHypothesis.set <- "greaterAbs"
alpha.set <- 0.05

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
                        contrast=list(c(resultsNames(dds.sub)[3], resultsNames(dds.sub)[4])),
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
    #res.ribo <- lfcShrink(dds.ribo, coef=2, type="apeglm", res=res.ribo)

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
    #res.rna <- lfcShrink(dds.rna, coef=2, type="apeglm", res=res.rna)

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
    dirloc.sub <- file.path(args$output, filen, fsep=.Platform$file.sep)
    created <- ifelse(!dir.exists(dirloc.sub), dir.create(dirloc.sub, recursive=TRUE), FALSE)
    if (!created) { print(paste("Using existing directory ", dirloc.sub, " ...", sep="")) }
    filen <- paste0(filen, ".xlsx", sep="")
    filen <- file.path(dirloc.sub, filen, fsep=.Platform$file.sep)
    print(paste("Writing results to: ", filen, sep=""))
    saveWorkbook(wb, filen, overwrite = TRUE)
    get_reg_layers(res, dirloc.sub)
}


read_symbols <- function(samples, dirloc) {
    # read gene symbols from HTSeq output tables
    l <- lapply( as.character( samples[,2] ), function(fn) read.table( file.path( dirloc, fn ), fill=TRUE ))
    if( ! all( sapply( l, function(a) all( head(a$V2, -5) == head(l[[1]]$V2, -5) ) ) ) )
        stop( "Gene names (second column) differ between files." )
    tbl <- sapply( l, function(a) a[,2] )
    colnames(tbl) <- samples[,1]
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
    both <- which(res$padj.dte < alpha.set & res$padj.ribo < alpha.set & res$padj.rna < alpha.set)
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
    points(y=res[res$gene %in% forwarded$gene,]$log2FC.ribo, x=res[res$gene %in% forwarded$gene,]$log2FC.rna, pch=21, col=rgb(0,0,1,1), cex=0.6)
    points(y=res[res$gene %in% exclusive$gene,]$log2FC.ribo, x=res[res$gene %in% exclusive$gene,]$log2FC.rna, pch=21, col=rgb(0,1,0,1), cex=0.6)
    points(y=res[res$gene %in% intensified$gene,]$log2FC.ribo, x=res[res$gene %in% intensified$gene,]$log2FC.rna, pch=21, col=rgb(1,0,0,1), cex=0.6)
    points(y=res[res$gene %in% buffered$gene,]$log2FC.ribo, x=res[res$gene %in% buffered$gene,]$log2FC.rna, pch=22, col=rgb(1,0,0,1), cex=0.6)
    legend(x="topleft",
           legend=c("Forwarded", "Exclusive (DTEG)", "Intensified", "Buffered"),
           col=c(rgb(0,0,1,1), rgb(0,1,0,1), rgb(1,0,0,1), rgb(1,0,0,1)), lwd=1, lty=c(0,0),
           pch=c(21,21,21,22), bty='n')
    dev.off()
}

# ---------------------------------------------------------
## Call

# arguments
defaults <- list(method="LRT", batch=FALSE, filter=FALSE)
args <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE, defaults=defaults)
if (length(args)<3) {
    stop("./run-tea.R [-sample SAMPLE] [-contrast CONTRAST] [-output OUTPUT] <-method LRT/deltaTE> <--batch> <--filter>\n",
         call.=FALSE)
}
# sample table file
samples <- read.delim(args$sample, sep=",")
# contrast file
contrasts <- read.delim(args$contrast, sep=",", header=FALSE)
# output directory
created <- ifelse(!dir.exists(args$output), dir.create(args$output, recursive=TRUE), FALSE)
if (!created) { print(paste0("Using existing directory ", args$output, "...")) }
# method
method <- "LRT"
if (args$method=="deltaTE") { method <- "deltaTE" }
print(paste0("Using the following method ", method, "..."))

# write tmp files to
tmp.loc <- uuid::UUIDgenerate()
dirloc.tmp <- file.path(args$output, tmp.loc, fsep=.Platform$file.sep)
dir.create(dirloc.tmp)

# wrangle
samples <- as.data.frame(samples) # need dataframe for DESEq
samples$condition <- factor(samples$condition)
samples$assay <- factor(samples$assay)

# check if we use batch
if (args$batch) {
    if (nlevels(sampleTable$batch)<2) {
        args$batch <- FALSE
        print("Ignoring --batch! Use 'batch' as column header, and at least 2 levels!")
    } else {
        samples$batch <- factor(samples$batch)
    }
}

# link all files for DESeq
for (from in samples$fileName) {
    to <- file.path(dirloc.tmp, basename(from), fsep=.Platform$file.sep)
    file.symlink(from, to)
}

# create DESEq data object
samples$fileName <- basename(samples$fileName)
model <- ~assay+condition+assay:condition
if (method == "LRT") {
  model.reduced <- ~assay+condition
} else if (method == "deltaTE") {
    model.reduced <- ~condition
    if (args$batch) { model.reduced <- ~batch+condition }
}
ddsHTSeq <- DESeqDataSetFromHTSeqCount(samples,
                                       directory=dirloc.tmp,
                                       design=model)

stopifnot(all(samples$sampleName == colnames(ddsHTSeq)))

# always rna as reference level for assay
ddsHTSeq$assay <- relevel(ddsHTSeq$assay, ref="rna")

# annotation
# the order should be the same...
id.mapping <- read_symbols(samples, dirloc.tmp)

# remove files
unlink(dirloc.tmp, recursive=TRUE)

# separate normalisation ???
# sf <- numeric(ncols(ddsHTSeq))
# idx.ribo <- ddsHTSeq$assay == 'ribo'
# sf[idx.ribo] <- estimateSizeFactorsFromMatrix(counts(ddsHTSeq)[,idx.ribo])
# idx.rna <- ddsHTSeq$assay == 'rna'
# sf[idx.rna] <- estimateSizeFactorsFromMatrix(counts(ddsHTSeq)[,idx.rna])
# sizeFactors(ddsHTSeq) <- sf

# filter

if (args$filter) {
    print("Removing genes with 0 counts in each assay separately, using the intersection...")

    ddsHTSeq.cts <- counts(ddsHTSeq)

    m_ribo <- ddsHTSeq.cts[,which(colData((ddsHTSeq))$assay == "ribo")]
    m_rna <- ddsHTSeq.cts[,which(colData((ddsHTSeq))$assay == "rna")]

    idx.ribo <- which(apply(m_ribo, 1, function(x){length(which(x>0))})==ncol(m_ribo))
    idx.rna <- which(apply(m_rna, 1, function(x){length(which(x>0))})==ncol(m_rna))

    idx <- intersect(idx.ribo, idx.rna)
    ddsHTSeq <- ddsHTSeq[idx,]
    id.mapping <- id.mapping[idx]
}

## fit contrasts
if (method == "LRT") {
    apply(contrasts, 1, run_analysis_fit_sub_lrt, coldata=samples, dds=ddsHTSeq, model=model.reduced)
} else if (method == "deltaTE") {
    apply(contrasts, 1, run_analysis_fit_sub_deltaTE, coldata=samples, dds=ddsHTSeq, model=model.reduced)
}
