######################################################################
### Find pausing sites in genes using a method similar to the one
### published in Stein et al. Nature 2022

library(tidyverse)
library(DESeq2)
library(org.Mm.eg.db)
library(cowplot)
library(Logolas)

ggColorHue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


readCorr <- function(file){
    c.1 <- read.table(file)
    c.cod <- tapply(c.1[1:(dim(c.1)[1]-2),2],
                    rep(seq(-99, 100-3, 3), each = 3), sum)
    cData <- data.frame(position = as.numeric(names(c.cod)) / 3,
                        value = c.cod / median(c.cod),
                        raw = c.cod)
    return(cData)
}

readChipExtract <- function(file){
    m1 <- read.table(file)
    colnames(m1) <- -99:100
    pos <- seq(-99, 100-3, 3)
    r1 <- matrix(ncol = length(pos), nrow = dim(m1)[1])
    colnames(r1) <- pos
    for(I in pos){
        r1[,as.character(I)] <- rowSums(m1[,as.character(c(I,I+1, I+2))])
    }
    return(r1)
}

condense.matrix <- function(mat, rows=20){
    groups <- sort(rep(1:ceiling(dim(mat)[1]/rows), rows))[1:dim(mat)[1]]
    c.mat <- matrix(ncol=dim(mat)[2], nrow=max(groups))
    for (I in 1:dim(mat)[2]){
        c.mat[,I] <- tapply(mat[,I], groups, mean, na.rm=T)
    }
    colnames(c.mat) = colnames(mat)
    return(c.mat)
}

compareCodons <- function(codonData, cdsData, sampleAnno, contrasts,
                          cmin = 0){
    id1 <- which(sampleAnno[,contrasts[1]] == contrasts[2])
    id2 <- which(sampleAnno[,contrasts[1]] == contrasts[3])
    cdsCounts1 <- round(rowSums(cdsData[,rownames(sampleAnno)[id1]]))
    cdsCounts2 <- round(rowSums(cdsData[,rownames(sampleAnno)[id2]]))
    codonCounts1 <- round(rowSums(codonData[,rownames(sampleAnno)[id1]]))
    codonCounts2 <- round(rowSums(codonData[,rownames(sampleAnno)[id2]]))
    codCounts.mean <- data.frame(codonData[,1:3],
                                 chow = codonCounts1,
                                 hfd = codonCounts2,
                                 chow.tot = cdsCounts1[codonData[,1]],
                                 hfd.tot = cdsCounts2[codonData[,1]]
                                 )
    codCounts.kept <- codCounts.mean[
        which(codCounts.mean$chow > cmin &
              codCounts.mean$hfd > cmin), ]
    oddRatio <- vector()
    pValues <- vector()
    for(I in 1:dim(codCounts.kept)[1]){
        x1 <- matrix(ncol = 2, nrow = 2)
        x1[1,1] <- codCounts.kept[I,4]
        x1[2,1] <- codCounts.kept[I,5]
        x1[1,2] <- codCounts.kept[I,6] - codCounts.kept[I,4]
        x1[2,2] <- codCounts.kept[I,7] - codCounts.kept[I,5]
        ft <- fisher.test(x=x1)
        pValues <- c(pValues, ft$p.value)
        oddRatio <- c(oddRatio, ft$estimate)
    }
    padj <- p.adjust(pValues)
    codCounts.kept <- cbind(codCounts.kept,
                            oddRatios = log2(oddRatio),
                            pValues = pValues, padj = padj)
    colnames(codCounts.kept) <- c(colnames(codonData)[1:3],
                                  contrasts[2:3],
                                  paste(contrasts[2:3], "cds", sep="."),
                                  "oddRatios", "pValue", "padj")
    rownames(codCounts.kept) <- paste(codCounts.kept$GID,
                                      codCounts.kept$Position,
                                      sep = ".")
    codCounts.kept$oddRatios[which(abs(codCounts.kept$oddRatios) == Inf)] <- NA
    return(codCounts.kept)
}


pauseScore <- function(codonData, cdsLengths){
    codonScore <- codonData
    for(I in 4:dim(codonData)[2]){
        codonSum <- tapply(codonData[,I], codonData[,1], sum)
        codonDens <- codonSum / (cdsLengths[names(codonSum)] / 3)
        codonScore[,I] <- codonData[,I] / codonDens[codonData[,1]]
    }
    return(codonScore)
}


plotPauseScores <- function(id, scores, cds, sampleAnno, contrasts){
    id1 <- which(sampleAnno[,contrasts[1]] == contrasts[2])
    id2 <- which(sampleAnno[,contrasts[1]] == contrasts[3])
    cdsLen <- cds[,4] - cds[,3] + 1
    cdsStart <- cds[,3]
    cdsStop <- cds[,4]
    names(cdsLen) <- names(cdsStart) <- names(cdsStop) <- cds[,1]
    d1 <- subset(codonsScore, GID == id)
    allPos <- seq(cdsStart[id], cdsStop[id]+1, 3)
    codPos <- (allPos - cdsStart[id]) / 3
    d2 <- cbind(codPos = rep(0, dim(d1)[1]), d1)
    for(I in 1:length(allPos)){
        p1 <- allPos[I]
        if (p1 %in% d1$Position){
            d2[which(d2$Position == p1), "codPos"] <- codPos[I]
        }else{
            newrow <- cbind(codPos = codPos[I], d2[1,2:dim(d2)[2]])
            newrow[,5:dim(d2)[2]] <- 0
            d2 <- rbind(d2, newrow)
        }
    }
    d1lt <- pivot_longer(d2, cols=5:dim(d2)[2],
                         names_to = "Sample", values_to = "Score")
    d1lt$contrast <- sampleAnno[d1lt$Sample, contrasts[1]]
    df1 <- subset(d1lt, contrast == contrasts[2] | contrast == contrasts[3])
    mdata <- aggregate(Score ~ codPos:contrast, df1, FUN=mean)
    p <- ggplot(df1, aes(x = codPos, y = Score, col = contrast)) +
        geom_point() +
        geom_line(data = mdata, aes())
    return(p)
}


######################################################################
### Read in codon data (calculated in 3 codons window)

codons <- read.table("../codonCounts/2022-05-09.mouse_cDNA.codonCountsMod.dat",
                     header = FALSE)

sampleOrder <- read.table("../2022-05-06_samples.order",
                          header = FALSE)

colnames(codons) <- c("GID", "TID", "Position",
                      paste(sampleOrder[,1],
                            sampleOrder[,2], sep = "."))

head(codons)

samplesAnnoRaw <- read.table("../2022-05-06_samples.anno",
                             header = TRUE)
rownames(samplesAnnoRaw) <- paste(samplesAnnoRaw$Library,
                                  samplesAnnoRaw$Adapter , sep = ".")

######################################################################
### Read in genes RPF numbers

sampleAnno <- data.frame(Sample   = paste(samplesAnnoRaw$Library,
                                          samplesAnnoRaw$Adapter,
                                          sep='.'),
                         FileName = paste(samplesAnnoRaw$Library,
                                          samplesAnnoRaw$Adapter,
                                          'readCounts', 'dat',
                                          sep = '.'),
                         Library = samplesAnnoRaw$Library,
                         Barcode = samplesAnnoRaw$Adapter,
                         Genotype = factor(samplesAnnoRaw$Genotype),
                         Treatment = factor(samplesAnnoRaw$Treatment,
                                            levels = c("SD", "Leu")),
                         Type = factor(paste(samplesAnnoRaw$Genotype,
                                             samplesAnnoRaw$Treatment,
                                             sep = ".")),
                         Replica = factor(samplesAnnoRaw$Replica)
                         )

head(sampleAnno)
str(sampleAnno)
rownames(sampleAnno) <- sampleAnno$Sample

directory <- '../readCounts'


ddsData <- DESeqDataSetFromHTSeqCount(sampleTable = sampleAnno,
                                      directory = directory,
                                      design = ~ Treatment + Genotype)

## keep <- rowMeans(counts(ddsData, norm=F)) > 10
## ddsData <- ddsData[keep,]
ddsData <- DESeq(ddsData)
resultsNames(ddsData)

rawCounts <- counts(ddsData, norm = FALSE)
head(rawCounts)
colSums(rawCounts)
cbind(sampleAnno, colSums(rawCounts))

######################################################################
### Start with analyzing treatment effect in the WT

codonCounts <- compareCodons(codons, rawCounts, sampleAnno,
                             c("Type", "WT.Leu", "WT.SD"),
                             cmin = 20)
head(codonCounts)
codonCounts.sign <- codonCounts[which(codonCounts$padj < 0.05),]


enrichedCod <- which(codonCounts.sign$oddRatios > 0)
depletedCod <- which(codonCounts.sign$oddRatios < 0)

sum(codonCounts.sign$WT.Leu[enrichedCod])
sum(codonCounts.sign$WT.SD[enrichedCod])

sum(codonCounts.sign$WT.Leu[depletedCod])
sum(codonCounts.sign$WT.SD[depletedCod])

length(enrichedCod)
length(depletedCod)

### Volcano plot:
vp1 <- ggplot(codonCounts, aes(x = oddRatios, y = -log10(padj))) +
    geom_point(col = "lightgray", alpha = 0.5, size = 0.8) +
    geom_point(data = codonCounts.sign[enrichedCod,],
               aes(x = oddRatios, y = -log10(padj)),
               col = "orange", size = 0.8) +
    geom_point(data = codonCounts.sign[depletedCod,],
               aes(x = oddRatios, y = -log10(padj)),
               col = "skyblue", size = 0.8) +
    labs(x = "odd Ratios") +
    theme_cowplot(14)
ggsave("figures/2022-05-09/volcanoPlotWt.pdf", vp1)

######################################################################
### Annotate genes:
affGenes <- unique(codonCounts.sign$GID)

gene.symbols<-mapIds(org.Mm.eg.db,
                     keys=affGenes,
                     column=c("SYMBOL"),
                     keytype="ENSEMBL",
                     multiVals='first')
gene.descr<-mapIds(org.Mm.eg.db,
                   keys=affGenes,
                   column=c("GENENAME"),
                   keytype="ENSEMBL",
                   multiVals='first')

head(gene.symbols)

dataOut<-data.frame(codonCounts.sign,
                    GBid = paste(codonCounts.sign$GID,
                                 codonCounts.sign$TID, sep = "|"),
                    Name = gene.symbols[codonCounts.sign$GID],
                    Description = gene.descr[codonCounts.sign$GID])

head(dataOut)

write.csv(dataOut,
          file = "results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons.csv")


######################################################################
######################################################################
### Do the treatment in the mutant
######################################################################
######################################################################

codCountsZak <- compareCodons(codons, rawCounts, sampleAnno,
                              c("Type", "ZAK.Leu", "ZAK.SD"),
                              cmin = 20)
codCountsZak$oddRatios[which(codCountsZak$oddRatios == -Inf)] <- NA
head(codCountsZak)

codCountsZak.sign <- codCountsZak[which(codCountsZak$padj < 0.05),]

as.numeric(summary(codCountsZak.sign$oddRatios > 0))

### Volcano plot:
ggplot(codCountsZak, aes(x = oddRatios, y = -log10(padj))) +
    geom_point()

######################################################################
### Annotate genes:
affGenes <- unique(codCountsZak.sign$GID)

gene.symbols<-mapIds(org.Mm.eg.db,
                     keys=affGenes,
                     column=c("SYMBOL"),
                     keytype="ENSEMBL",
                     multiVals='first')
gene.descr<-mapIds(org.Mm.eg.db,
                   keys=affGenes,
                   column=c("GENENAME"),
                   keytype="ENSEMBL",
                   multiVals='first')

dataOut<-data.frame(codCountsZak.sign,
                    GBid = paste(codCountsZak.sign$GID,
                                 codCountsZak.sign$TID, sep = "|"),
                    Name = gene.symbols[codCountsZak.sign$GID],
                    Description = gene.descr[codCountsZak.sign$GID])

head(dataOut)

write.csv(dataOut,
          file = "results/2022-05-09/ZAKdeficientVsStandardDiet_differentialCodon.csv")


######################################################################
######################################################################
### Compare the two treatments:
######################################################################
######################################################################

rownames(codCountsZak.sign) <- paste(codCountsZak.sign$GID,
                                     codCountsZak.sign$Position,
                                     sep = ".")
rownames(codonCounts.sign) <- paste(codonCounts.sign$GID,
                                    codonCounts.sign$Position,
                                    sep = ".")

commonSites <- rownames(codCountsZak.sign)[which(rownames(codCountsZak.sign) %in%
                                                 rownames(codonCounts.sign))]
allSite <- unique(c(rownames(codCountsZak.sign),
                    rownames(codonCounts.sign)))


length(commonSites)
length(allSite)

summary(codCountsZak[allSite, "oddRatios"])

plot(codCountsZak[allSite, "oddRatios"],
     codonCounts[allSite, "oddRatios"], pch = ".")
points(codCountsZak[commonSites, "oddRatios"],
       codonCounts[commonSites, "oddRatios"], pch = 20, col = "red")

### Compare volcano plots:

noTreatSites <- rownames(codonCounts)[which(codonCounts$padj < 0.05)]
zakSites <- rownames(codCountsZak)[which(codCountsZak$padj < 0.05)]

length(noTreatSites)
length(zakSites)

commonSites <- noTreatSites[which(noTreatSites %in% zakSites)]

### Plot Non treated and ZAK positive sites
set2col <- RColorBrewer::brewer.pal(8, "Set2")
set1col <- RColorBrewer::brewer.pal(9, "Set1")

defCol <- ggColorHue(2)

p1 <- ggplot(codonCounts, aes(x = oddRatios, y = -log10(padj))) +
    geom_point(alpha = 1, col = "lightgray", size = 0.5) +
    geom_point(data = subset(codonCounts[noTreatSites,]),
               col = defCol[2], size = 0.5) +
    ylim(0, 150) +
    xlim(-3, 3) +
    labs(title = "WT",
         subtitle = "Deficient vs. standard diet",
         x = "log2[oddRatio]",
         y = "log10[adjusted P]") +
    theme_cowplot(12)
p2 <- ggplot(codCountsZak, aes(x = oddRatios, y = -log10(padj))) +
    geom_point(alpha = 1, col = "lightgray", size = 0.5) +
    geom_point(data = subset(codCountsZak[zakSites,]),
               col = defCol[1], size = 0.5) +
    ylim(0, 150) +
    xlim(-3, 3) +
    labs(title = "ZAK",
         subtitle = "Deficient vs. standard diet",
         x = "log2[oddRatio]",
         y = "log10[adjusted P]") +
    theme_cowplot(12)
p3 <- plot_grid(p1, p2)
ggsave("figures/2022-05-09/volcanoPlots.png", p3, width = 8, height = 4)


### Plot Non treated site on ZAK:
p1 <- ggplot(codonCounts, aes(x = oddRatios, y = -log10(padj))) +
    geom_point(alpha = 1, col = "lightgray", size = 0.5) +
    geom_point(data = subset(codonCounts[noTreatSites,]),
               col = defCol[2], size = 0.5) +
    ylim(0, 150) +
    xlim(-3, 3) +
    labs(title = "WT",
         subtitle = "Deficient vs. standard diet",
         x = "log2[oddRatio]",
         y = "log10[adjusted P]") +
    theme_cowplot(12)
p2 <- ggplot(codCountsZak, aes(x = oddRatios, y = -log10(padj))) +
    geom_point(alpha = 1, col = "lightgray", size = 0.5) +
    geom_point(data = subset(codCountsZak[noTreatSites,]),
               col = defCol[2], size = 0.5) +
    ylim(0, 150) +
    xlim(-3, 3) +
    labs(title = "ZAK",
         subtitle = "Deficient vs. standard diet",
         x = "log2[oddRatio]",
         y = "log10[adjusted P]") +
    theme_cowplot(12)
p3 <- plot_grid(p1, p2)
ggsave("figures/2022-05-09/volcanoPlots_WTsites.png", p3,
       width = 8, height = 4)


### Plot ZAK sites on non treated sites
p1 <- ggplot(codonCounts, aes(x = oddRatios, y = -log10(padj))) +
    geom_point(alpha = 1, col = "lightgray", size = 0.5) +
    geom_point(data = subset(codonCounts[zakSites,]),
               col = defCol[1], size = 0.5) +
    ylim(0, 150) +
    xlim(-3, 3) +
    labs(title = "WT",
         subtitle = "Deficient vs. standard diet",
         x = "log2[oddRatio]",
         y = "log10[adjusted P]") +
    theme_cowplot(12)
p2 <- ggplot(codCountsZak, aes(x = oddRatios, y = -log10(padj))) +
    geom_point(alpha = 1, col = "lightgray", size = 0.5) +
    geom_point(data = subset(codCountsZak[zakSites,]),
               col = defCol[1], size = 0.5) +
    ylim(0, 150) +
    xlim(-3, 3) +
    labs(title = "ZAK",
         subtitle = "Deficient vs. standard diet",
         x = "log2[oddRatio]",
         y = "log10[adjusted P]") +
    theme_cowplot(12)
p3 <- plot_grid(p1, p2)
ggsave("figures/2022-05-09/volcanoPlots_zakSites.png", p3,
       width = 8, height = 4)

### Plot common
p1 <- ggplot(codonCounts, aes(x = oddRatios, y = -log10(padj))) +
    geom_point(alpha = 1, col = "lightgray", size = 0.5) +
    geom_point(data = subset(codonCounts[commonSites,]),
               col = set1col[4], size = 0.5) +
    ylim(0, 150) +
    xlim(-3, 3) +
    labs(title = "WT",
         subtitle = "Deficient vs. standard diet",
         x = "log2[oddRatio]",
         y = "log10[adjusted P]") +
    theme_cowplot(12)
p2 <- ggplot(codCountsZak, aes(x = oddRatios, y = -log10(padj))) +
    geom_point(alpha = 1, col = "lightgray", size = 0.5) +
    geom_point(data = subset(codCountsZak[commonSites,]),
               col = set1col[4], size = 0.5) +
    ylim(0, 150) +
    xlim(-3, 3) +
    labs(title = "ZAK",
         subtitle = "Deficient vs. standard diet",
         x = "log2[oddRatio]",
         y = "log10[adjusted P]") +
    theme_cowplot(12)
p3 <- plot_grid(p1, p2)
ggsave("figures/2022-05-09/volcanoPlots_commonSites.png", p3,
       width = 8, height = 4)


eData <- data.frame(counts = c(as.numeric(summary(codCountsZak.sign$oddRatios > 0))[2:3],
                               as.numeric(summary(codonCounts.sign$oddRatios > 0))[2:3]),
                    genotype = rep(c("ZAK", "WT"), each = 2),
                    type = rep(c("Depleted", "Enriched"), 2)
                    )

p4 <- ggplot(eData, aes(x = counts, y = genotype, fill = type)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_brewer(palette="Set1") +
    labs(x = "",
         y = "",
         fill = "") +
    theme_minimal(16)
ggsave("figures/2022-05-09/differentialSites.pdf", p4,
       width = 7, height = 2)



######################################################################
### Compare WT.Leu vs ZAK.Leu

codCountsLeu <- compareCodons(codons, rawCounts, sampleAnno,
                              c("Type", "WT.Leu", "ZAK.Leu"),
                              cmin = 20)
codCountsLeu.sign <- codCountsLeu[which(codCountsLeu$padj < 0.05),]

summary(codCountsLeu.sign$oddRatios > 0)

### Volcano plot:
p5 <- ggplot(codCountsLeu, aes(x = oddRatios, y = -log10(padj))) +
    geom_point(alpha = 1, col = "lightgray", size = 0.5) +
    geom_point(data = codCountsLeu.sign,
               col = set1col[3], size = 0.8) +
    ylim(0, 150) +
    xlim(-3, 3) +
    labs(title = "WT vs ZAK",
         subtitle = "Deficient diet",
         x = "log2[oddRatio]",
         y = "log10[adjusted P]") +
    theme_cowplot(12)
ggsave("figures/2022-05-09/volcano_WTvsZAKdeficientDiet.png", p5,
       width = 7, height = 7)



######################################################################
### Annotate genes:
affGenes <- unique(codCountsLeu.sign$GID)
length(affGenes)

gene.symbols<-mapIds(org.Mm.eg.db,
                     keys=affGenes,
                     column=c("SYMBOL"),
                     keytype="ENSEMBL",
                     multiVals='first')
gene.descr<-mapIds(org.Mm.eg.db,
                   keys=affGenes,
                   column=c("GENENAME"),
                   keytype="ENSEMBL",
                   multiVals='first')

dataOut<-data.frame(codCountsLeu.sign,
                    GBid = paste(codCountsLeu.sign$GID,
                                 codCountsLeu.sign$TID, sep = "|"),
                    Name = gene.symbols[codCountsLeu.sign$GID],
                    Description = gene.descr[codCountsLeu.sign$GID])

head(dataOut)

write.csv(dataOut,
          file = "results/2022-05-09/WTvsZAKdeficientDiet_differentialCodon.csv")

######################################################################
### Compare WT vs ZAK standard diet

codCountsSD <- compareCodons(codons, rawCounts, sampleAnno,
                               c("Type", "WT.SD", "ZAK.SD"),
                               cmin = 20)
codCountsSD.sign <- codCountsSD[which(codCountsSD$padj < 0.05),]

summary(codCountsSD.sign$oddRatios > 0)

### Volcano plot:
p6 <- ggplot(codCountsSD, aes(x = oddRatios, y = -log10(padj))) +
    geom_point(alpha = 1, col = "lightgray", size = 0.5) +
    geom_point(data = codCountsSD.sign,
               col = set1col[5], size = 0.8) +
    ylim(0, 150) +
    xlim(-3, 3) +
    labs(title = "WT vs ZAK",
         subtitle = "Standard diet",
         x = "log2[oddRatio]",
         y = "log10[adjusted P]") +
    theme_cowplot(12)
ggsave("figures/2022-05-09/volcano_WTvsZAKstandardDiet.png", p6,
       width = 7, height = 7)

### Volcano plot:
p7 <- plot_grid(p5, p6)
ggsave("figures/2022-05-09/volcanoPlots_otherComparisons.png", p7,
       width = 8, height = 4)


######################################################################
### Annotate genes:
affGenes <- unique(codCountsSD.sign$GID)
length(affGenes)

gene.symbols<-mapIds(org.Mm.eg.db,
                     keys=affGenes,
                     column=c("SYMBOL"),
                     keytype="ENSEMBL",
                     multiVals='first')
gene.descr<-mapIds(org.Mm.eg.db,
                   keys=affGenes,
                   column=c("GENENAME"),
                   keytype="ENSEMBL",
                   multiVals='first')

dataOut<-data.frame(codCountsSD.sign,
                    GBid = paste(codCountsSD.sign$GID,
                                 codCountsSD.sign$TID, sep = "|"),
                    Name = gene.symbols[codCountsSD.sign$GID],
                    Description = gene.descr[codCountsSD.sign$GID])

dim(dataOut)

write.csv(dataOut,
          file = "results/2022-05-09/WTvsZAKstandardDiet_differentialCodon.csv")

######################################################################
######################################################################
### Plot coverage around enriched sites:
######################################################################
######################################################################

normFact <- colSums(rawCounts) / mean(colSums(rawCounts))

wtSd1.d <- readCorr('../chipCor/2022-05-06/PF35.ATCGT_vs_WT_DD.Enriched.dat')
wtSd2.d <- readCorr('../chipCor/2022-05-06/PF35.AGCTA_vs_WT_DD.Enriched.dat')

wtDd1.d <- readCorr('../chipCor/2022-05-06/PF35.CGTAA_vs_WT_DD.Enriched.dat')
wtDd2.d <- readCorr('../chipCor/2022-05-06/PF35.CTAGA_vs_WT_DD.Enriched.dat')

normData <- rbind(wtSd1.d, wtDd1.d)
normData$type <- rep(c("SD", "DD"), each = dim(wtSd1.d)[1])
normData$value <- c((wtSd1.d$value + wtSd2.d$value) / 2, (wtDd1.d$value + wtDd2.d$value) / 2)
normData$raw <- c((wtSd1.d$raw/normFact['PF35.ATCGT'] + wtSd2.d$raw/normFact['PF35.AGCTA']) / 2, (wtDd1.d$raw/normFact['PF35.CGTAA'] + wtDd2.d$raw/normFact['PF35.CTAGA']) / 2)

head(normData)

ggplot(normData, aes(x = position, y = value, col = type)) +
    geom_vline(xintercept = 0) +
    geom_line() +
    labs(x = "Distance from pause site (codon)",
         y = "Normalised ribosome occupancy",
         col = "") +
    xlim(-20, 20) +
    ## ylim(0.8, 1.6) +
    theme_cowplot(16) +
    theme(legend.position = c(0.1, 0.9))



#### Now the NAK in the same regions:
zakSd1.d <- readCorr('../chipCor/2022-05-06/PF35.GATCA_vs_WT_DD.Enriched.dat')
zakSd2.d <- readCorr('../chipCor/2022-05-06/PF35.GCATA_vs_WT_DD.Enriched.dat')

zakDd1.d <- readCorr('../chipCor/2022-05-06/PF35.TAGAC_vs_WT_DD.Enriched.dat')
zakDd2.d <- readCorr('../chipCor/2022-05-06/PF35.TCTAG_vs_WT_DD.Enriched.dat')

normZak <- rbind(zakSd1.d, zakDd1.d)
normZak$type <- rep(c("SD", "DD"), each = dim(zakSd1.d)[1])
normZak$value <- c((zakSd1.d$value + zakSd2.d$value) / 2, (zakDd1.d$value + zakDd2.d$value) / 2)
normZak$raw <- c((zakSd1.d$raw/normFact['PF35.GATCA'] + zakSd2.d$raw/normFact['PF35.GCATA']) / 2, (zakDd1.d$raw/normFact['PF35.TAGAC'] + zakDd2.d$raw/normFact['PF35.TCTAG']) / 2)

head(normZak)

ggplot(normZak, aes(x = position, y = value, col = type)) +
    geom_vline(xintercept = 0) +
    geom_line() +
    labs(x = "Distance from pause site (codon)",
         y = "Normalised ribosome occupancy",
         col = "") +
    xlim(-20, 20) +
    ## ylim(0.8, 1.6) +
    theme_cowplot(16) +
    theme(legend.position = c(0.1, 0.9))





#### Depleted:
wtSd1.dd <- readCorr('../chipCor/2022-05-06/PF35.ATCGT_vs_WT_DD.Depleted.dat')
wtSd2.dd <- readCorr('../chipCor/2022-05-06/PF35.AGCTA_vs_WT_DD.Depleted.dat')

wtDd1.dd <- readCorr('../chipCor/2022-05-06/PF35.CGTAA_vs_WT_DD.Depleted.dat')
wtDd2.dd <- readCorr('../chipCor/2022-05-06/PF35.CTAGA_vs_WT_DD.Depleted.dat')

normDd <- rbind(wtSd1.dd, wtDd1.dd)
normDd$type <- rep(c("SD", "DD"), each = dim(wtSd1.dd)[1])
normDd$value <- c((wtSd1.dd$value + wtSd2.dd$value) / 2, (wtDd1.dd$value + wtDd2.dd$value) / 2)
normDd$raw <- c((wtSd1.dd$raw/normFact['PF35.ATCGT'] + wtSd2.dd$raw/normFact['PF35.AGCTA']) / 2, (wtDd1.dd$raw/normFact['PF35.CGTAA'] + wtDd2.dd$raw/normFact['PF35.CTAGA']) / 2)

head(normDd)

ggplot(normDd, aes(x = position, y = log(raw), col = type)) +
    geom_vline(xintercept = 0) +
    geom_line() +
    labs(x = "Distance from pause site (codon)",
         y = "Normalised ribosome occupancy",
         col = "") +
    xlim(-20, 20) +
    ## ylim(0.8, 1.6) +
    theme_cowplot(16) +
    theme(legend.position = c(0.1, 0.9))


#### Depleted:
zakSd1.dd <- readCorr('../chipCor/2022-05-06/PF35.GATCA_vs_WT_DD.Depleted.dat')
zakSd2.dd <- readCorr('../chipCor/2022-05-06/PF35.GCATA_vs_WT_DD.Depleted.dat')

zakDd1.dd <- readCorr('../chipCor/2022-05-06/PF35.TAGAC_vs_WT_DD.Depleted.dat')
zakDd2.dd <- readCorr('../chipCor/2022-05-06/PF35.TCTAG_vs_WT_DD.Depleted.dat')

normZakDd <- rbind(zakSd1.dd, zakDd1.dd)
normZakDd$type <- rep(c("SD", "DD"), each = dim(zakSd1.dd)[1])
normZakDd$value <- c((zakSd1.dd$value + zakSd2.dd$value) / 2, (zakDd1.dd$value + zakDd2.dd$value) / 2)
normZakDd$raw <- c((zakSd1.dd$raw/normFact['PF35.GATCA'] + zakSd2.dd$raw/normFact['PF35.GCATA']) / 2, (zakDd1.dd$raw/normFact['PF35.TAGAC'] + zakDd2.dd$raw/normFact['PF35.TCTAG']) / 2)

head(normZakDd)

ggplot(normZakDd, aes(x = position, y = log(raw), col = type)) +
    geom_vline(xintercept = 0) +
    geom_line() +
    labs(x = "Distance from pause site (codon)",
         y = "Normalised ribosome occupancy",
         col = "") +
    xlim(-20, 20) +
    ## ylim(0.8, 1.6) +
    theme_cowplot(16) +
    theme(legend.position = c(0.1, 0.9))


### Plot all together
allCorr <- rbind(normData, normZak, normDd, normZakDd)
allCorr$diet <- factor(rep(c("Enriched Sites", "Depleted sites"),
                           each = dim(normZakDd)[1] * 2))
allCorr$geno <- factor(rep(rep(c("WT", "ZAK"),
                               each = dim(normZakDd)[1]), 2))

pg1 <- ggplot(allCorr, aes(x = position, y = log(raw), col = type)) +
    geom_vline(xintercept = 0) +
    geom_line(size = 1.1) +
    labs(title = "Footprints distribution around WT deficient diet sites",
         x = "Distance from pause site (codon)",
         y = "Normalised ribosome occupancy",
         col = "") +
    xlim(-20, 20) +
    ylim(9.5, 11.1) +
    facet_grid(geno ~ diet) +
    theme_cowplot(16) +
    theme(legend.position = c(0.05, 0.9))
ggsave("figures/2022-05-09/footprintDistributionAroundWTsites.pdf",
       pg1, width = 10, height = 10)


######################################################################
### Make a matrix plot for them:

## Plot enriched sites in depleted diet
wtSd1.m <- readChipExtract('../chipCor/2022-05-06/PF35.ATCGT_vs_WT_DD.Enriched.mat')
wtSd2.m <- readChipExtract('../chipCor/2022-05-06/PF35.AGCTA_vs_WT_DD.Enriched.mat')
wtSd.m <- (wtSd1.m + wtSd2.m)

wtDd1.m <- readChipExtract('../chipCor/2022-05-06/PF35.CGTAA_vs_WT_DD.Enriched.mat')
wtDd2.m <- readChipExtract('../chipCor/2022-05-06/PF35.CTAGA_vs_WT_DD.Enriched.mat')
wtDd.m <- (wtDd1.m + wtDd2.m)


m = colMeans(wtDd.m)
## order = order(cov(t(wtDd.m), m))
order <- order(rowSums(wtDd.m), decreasing = T)
w <- ceiling(dim(wtDd.m)[1]/1000)
x1 <- condense.matrix(wtDd.m[order,], rows=5)
max <- sort(x1, decreasing=T)[0.003*dim(x1)[1]*dim(x1)[2]]
for(i in 1:dim(x1)[1]) {for (j in 1:dim(x1)[2]) {x1[i,j]=min(x1[i,j], max)}}
### 
x2 <- condense.matrix(wtSd.m[order,], rows=5)
max <- sort(x2, decreasing=T)[0.003*dim(x2)[1]*dim(x2)[2]]
for(i in 1:dim(x2)[1]) {for (j in 1:dim(x2)[2]) {x2[i,j]=min(x2[i,j], max)}}


pdf("figures/2022-05-09/depletedDietEnrichedSites.pdf", width = 4,
    height = 3)
color <- colorRampPalette(c("white", "black"), space = "rgb")(100)
layout(matrix(c(1,2), nrow=1, ncol=2), widths = c(1, 1))
par(mar=c(1,1,1,1))
image(t(x1), col=color, xaxt="n", yaxt="n", bty="n")
par(mar=c(1,1,1,1))
image(t(x2), col=color, xaxt="n", yaxt="n", bty="n")
dev.off()

### Plot the depleted sites in depleted diet
wtSdD1.m <- readChipExtract('../chipCor/2022-05-06/PF35.ATCGT_vs_WT_DD.Depleted.mat')
wtSdD2.m <- readChipExtract('../chipCor/2022-05-06/PF35.AGCTA_vs_WT_DD.Depleted.mat')
wtSdD.m <- (wtSdD1.m + wtSdD2.m)

wtDdD1.m <- readChipExtract('../chipCor/2022-05-06/PF35.CGTAA_vs_WT_DD.Depleted.mat')
wtDdD2.m <- readChipExtract('../chipCor/2022-05-06/PF35.CTAGA_vs_WT_DD.Depleted.mat')
wtDdD.m <- (wtDdD1.m + wtDdD2.m)


m = colMeans(wtDdD.m)
## order = order(cov(t(wtDdD.m), m))
order <- order(rowSums(wtDdD.m), decreasing = T)
w <- ceiling(dim(wtDdD.m)[1]/1000)
x1 <- condense.matrix(wtDdD.m[order,], rows=5)
max <- sort(x1, decreasing=T)[0.003*dim(x1)[1]*dim(x1)[2]]
for(i in 1:dim(x1)[1]) {for (j in 1:dim(x1)[2]) {x1[i,j]=min(x1[i,j], max)}}
###
x2 <- condense.matrix(wtSdD.m[order,], rows=5)
max <- sort(x2, decreasing=T)[0.003*dim(x2)[1]*dim(x2)[2]]
for(i in 1:dim(x2)[1]) {for (j in 1:dim(x2)[2]) {x2[i,j]=min(x2[i,j], max)}}



pdf("figures/2022-05-09/depletedDietDepletedSites.pdf", width = 4,
    height = 3)
color <- colorRampPalette(c("white", "black"), space = "rgb")(100)
layout(matrix(c(1,2), nrow=1, ncol=2), widths = c(1, 1))
par(mar=c(1,1,1,1))
image(t(x1), col=color, xaxt="n", yaxt="n", bty="n")
par(mar=c(1,1,1,1))
image(t(x2), col=color, xaxt="n", yaxt="n", bty="n")
dev.off()




######################################################################
### Check distribution of up and down sites:
ratio <- sum(codonCounts.sign$oddRatios < 0)[1] / sum(codCountsZak.sign$oddRatios < 0)[1]


wtDd <- read.table("../chipCor/2022-05-06/WT_DD_UP_vs_DOWN.dat")
colnames(wtDd) <- c("Position", "Counts")
wtDd$Counts <- wtDd$Counts / sum(codonCounts.sign$oddRatios > 0)

zakDd <- read.table("../chipCor/2022-05-06/ZAK_DD_UP_vs_DOWN.dat")
colnames(zakDd) <- c("Position", "Counts")
zakDd$Counts <- zakDd$Counts / sum(codCountsZak.sign$oddRatios > 0)

dataDd <- rbind(wtDd, zakDd)
dataDd$type <- rep(c("WT", "ZAK"), each = dim(wtDd)[1])

ggplot(dataDd, aes(x = Position, y = Counts, col = type)) +
    geom_vline(alpha = 0.2, xintercept = 0, linetype = "dashed") +
    geom_line(size = 1.3, alpha = 0.5) +
    labs(title = "",
         x = "Distance from enriched sites (nt)",
         y = "Depleted sites counts",
         col = "Genotype") +
    xlim(-2000, 2000) +
    theme_cowplot(16) +
    theme(legend.position = c(0.1, 0.9))
ggsave("figures/2022-05-09/upVsDownDistribution.pdf",
       width = 7, height = 7)


ggplot(dataDd, aes(x = Position, y = Counts, col = type, fill = type)) +
    geom_vline(alpha = 0.2, xintercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.6, size = 2) +
    geom_line(alpha = 1, linetype = "dotted") +
    geom_smooth(method = "loess", se = TRUE, span = 0.25, alpha = 0.15) +
    labs(title = "WT, deficient diet",
         x = "Distance from enriched sites (nt)",
         y = "Depleted sites enrichment",
         col = "Treatment",
         fill = "Treatment") +
    xlim(-2000, 2000) +
    theme_cowplot(16) +
    theme(legend.position = c(0.1, 0.9))
ggsave("figures/2022-05-09/upVsDownDistributionSmooth.pdf",
       width = 7, height = 7)


######################################################################
### Make an autocorrelation plot for enriched sites to see if the re
### is any evidence of collisions:

wtDdUpAuto <- read.table('../chipCor/2022-05-06/WT_DD_UP_autocorrelation.dat')
colnames(wtDdUpAuto) <- c("position", "counts")
wtDdUpAuto$position <- wtDdUpAuto$position / 3
## wtDdUpAuto$counts <- wtDdUpAuto$counts / ratio

wtDdDownAuto <- read.table('../chipCor/2022-05-06/WT_DD_DOWN_autocorrelation.dat')
colnames(wtDdDownAuto) <- c("position", "counts")
wtDdDownAuto$position <- wtDdDownAuto$position / 3

wtDdAuto <- rbind(wtDdUpAuto, wtDdDownAuto)
wtDdAuto$type <- rep(c("Enriched", "Depleted"), each = dim(wtDdUpAuto)[1])


ap1 <- ggplot(wtDdAuto, aes(x = position, y = counts, col = type)) +
    geom_line(size = 1.3) +
    ## geom_vline(xintercept = seq(-40, 40, 10), linetype = "dashed",
    ##            alpha = 0.5) +
    labs(title = "WT, deficient diet",
         x = "Dinstance from Sites (codon)",
         y = "Sites counts",
         col = "Sites") +
    xlim(-40 , 40) +
    theme_half_open(16) +
    background_grid()
ggsave("figures/2022-05-09/wtDdEnrichedSitesAutocorrelation.pdf", ap1,
       width = 6, height = 7)

######################################################################
### Draw logo for pausing sites:
bgRaw <- read.table('results/2022-03-03/HFD_Chow_differentialCodons_UP.aaFreq',
                 header = FALSE, sep = "\t")

bg <- bgRaw[,2] / sum(bgRaw[,2])
names(bg) <- bgRaw[,1]

aaSeqNone <- read.table('results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_UP_AminoA.pwm',
                        header = FALSE, row.names = 1)
colnames(aaSeqNone) <- -10:10

aaSeqZak <- read.table('results/2022-05-09/HFD.ZAK_Chow.ZAK_differentialCodons_UP_AminoA.pwm',
                        header = FALSE, row.names = 1)
colnames(aaSeqZak) <- -10:10



noCounts <- which(rowSums(aaSeqNone) == 0)
aaSeqNone[-noCounts,]
aaSeqZak[-noCounts,]


pdf('figures/2022-05-09/aaSeqNoneLogo.pdf', width = 10)
logomaker(aaSeqNone[-noCounts,], type = "EDLogo", bg = bg[-noCounts])
dev.off()


pdf('figures/2022-05-09/aaSeqZakLogo.pdf', width = 10)
logomaker(aaSeqZak[-noCounts,], type = "EDLogo", bg = bg[-noCounts])
dev.off()


######################################################################
### Calculate amino acid frequency near the pausing site
bgRaw <- read.table('../config/Mus_musculus.GRCm38.100.startStop.APPRIS_expressed.aaFreq',
                 header = FALSE, sep = " ")

bg <- bgRaw[,2] / sum(bgRaw[,2])
names(bg) <- bgRaw[,1]


aaUpSite <- read.table('results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_UP.aaFreq',
                       header = FALSE, sep = " ")
aaUp <- aaUpSite[,2] / sum(aaUpSite[,2])
names(aaUp) <- aaUpSite[,1]
totAa <- sum(aaUpSite[,2])

aaDwSite <- read.table('results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_DOWN.aaFreq',
                       header = FALSE, sep = " ")
aaDw <- aaDwSite[,2] / sum(aaDwSite[,2])
names(aaDw) <- aaDwSite[,1]

aaData <- data.frame(counts = round(c(bg * totAa, aaUp * totAa,
                                      aaDw * totAa)),
                     amino  = c(names(bg), names(aaUp), names(aaDw)),
                     type   = factor(c(rep("Expected", length(bg)),
                                       rep("Enriched", length(aaUp)),
                                       rep("Depleted", length(aaDw))),
                                     levels = c("Enriched",
                                                "Depleted",
                                                "Expected")
                                     )
                     )

subset(aaData, amino == "Y")

aa1 <- ggplot(aaData, aes(x = amino, y = counts, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_y_continuous(
        expand = expansion(mult = c(0, 0.05))
    ) +
    labs(x = "Amino acid",
         y = "Normalised counts",
         fill = "") +
    ## theme_half_open(14) +
    theme_minimal_hgrid(14) + 
    theme(legend.position = c(0.85, 0.95))
ggsave('figures/2022-05-09/aaEnrichment.pdf', aa1, width = 10,
       height = 6)


aaData2 <- data.frame(enrichment = log2(c(aaUp / bg[names(aaUp)],
                                          aaDw / bg[names(aaDw)])),
                      amino  = c(names(aaUp), names(aaDw)),
                      type   = factor(c(rep("Enriched", length(aaUp)),
                                        rep("Depleted", length(aaDw))),
                                      levels = c("Depleted",
                                                 "Enriched")
                                      )
                      )

sd(aaData2$enrichment)
mean(aaData2$enrichment)

aa1 <- ggplot(aaData2, aes(x = amino, y = enrichment, fill = type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.6) +
    geom_hline(yintercept = c(sd(aaData2$enrichment),
                              -1*sd(aaData2$enrichment)),
               linetype = "dashed", alpha = 0.6) +
    ## scale_y_continuous(
    ##     expand = expansion(mult = c(0, 0.05))
    ## ) +
    scale_fill_manual(values=c("skyblue", "orange")) +
    labs(x = "Amino acid",
         y = "Enrichment",
         fill = "") +
    ## theme_half_open(14) +
    theme_minimal_hgrid(14) + 
    theme(legend.position = "none")
ggsave('figures/2022-05-09/aaEnrichment_2.pdf', aa1, width = 10,
       height = 5)


######################################################################
### Make an autocorrelation plot for enriched sites that have a Y in
### the center
wtDdYUpAuto <- read.table('../chipCor/2022-05-06/WT_DD_UP_Y_autocorrelation.dat')
colnames(wtDdYUpAuto) <- c("position", "counts")
wtDdYUpAuto$position <- wtDdYUpAuto$position / 3
wtDdYUpAuto$counts[wtDdYUpAuto$position == 0] <- 0
wtDdYUpAuto$enrichment <- wtDdYUpAuto$counts / max(wtDdYUpAuto$counts)

wtDdYUpCorr <- read.table('../chipCor/2022-05-06/WT_DD_UP_Y_corr.dat')
colnames(wtDdYUpCorr) <- c("position", "counts")
wtDdYUpCorr$position <- wtDdYUpCorr$position / 3
wtDdYUpCorr$counts[wtDdYUpCorr$position == 0] <- 0
wtDdYUpCorr$enrichment <- wtDdYUpCorr$counts / max(wtDdYUpCorr$counts)

wtDdYDownAuto <- read.table('../chipCor/2022-05-06/WT_DD_DOWN_Y_autocorrelation.dat')
colnames(wtDdYDownAuto) <- c("position", "counts")
wtDdYDownAuto$position <- wtDdYDownAuto$position / 3
wtDdYDownAuto$enrichment <- wtDdYDownAuto$counts / max(wtDdYDownAuto$counts)

wtDdUpAuto$enrichment <- wtDdUpAuto$counts / max(wtDdUpAuto$counts)
wtDdDownAuto$enrichment <- wtDdDownAuto$counts / max(wtDdDownAuto$counts)

wtDdYAuto <- rbind(wtDdUpAuto, wtDdYUpAuto, wtDdYUpCorr)
wtDdYAuto$type <- rep(c("Enriched"), each = dim(wtDdYUpAuto)[1]*3)
wtDdYAuto$site <- rep(c("All vs All", "Y vs Y", "Y vs All"),
                      each = dim(wtDdYUpAuto)[1])
wtDdYAuto$class <- paste(wtDdYAuto$type, wtDdYAuto$site, sep = ", ")

ap1 <- ggplot(subset(wtDdYAuto, type == "Enriched"),
              aes(x = position, y = enrichment, col = site)) +
    geom_line(size = 1.2) +
    ## geom_vline(xintercept = seq(-40, 40, 10), linetype = "dashed",
    ##            alpha = 0.5) +
    labs(title = "WT, deficient diet",
         x = "Dinstance from Sites (codon)",
         y = "Normlised enrichment",
         col = "Sites") +
    scale_x_continuous(breaks = seq(-40, 40, 10), limits = c(-40, 40)) +
    theme_half_open(16) +
    background_grid()
ggsave("figures/2022-05-09/wtDdEnrichedYSitesAutocorrelation.pdf", ap1,
       width = 6, height = 7)

ap1 <- ggplot(subset(wtDdYAuto, type == "Enriched"),
              aes(x = position, y = counts, col = site)) +
    geom_line(size = 1.2) +
    ## geom_vline(xintercept = seq(-40, 40, 10), linetype = "dashed",
    ##            alpha = 0.5) +
    labs(title = "WT, deficient diet",
         x = "Dinstance from Sites (codon)",
         y = "Counts",
         col = "Sites") +
    scale_x_continuous(breaks = seq(-40, 40, 10), limits = c(-40, 40)) +
    theme_half_open(16) +
    background_grid()
ggsave("figures/2022-05-09/wtDdEnrichedYSitesAutocorrelationCounts.pdf", ap1,
       width = 6, height = 7)


######################################################################
### Plot Pause score of some interesting genes
cdsLenData <- read.table("../config/Mus_musculus.GRCm38.100.APPRISuniq.cdsLength.txt",
                         header = FALSE)
cdsLen <- cdsLenData[,4] - cdsLenData[,3] + 1
names(cdsLen) <- cdsLenData[,1]

codonsScore <- pauseScore(codons, cdsLen)

id <- "ENSMUSG00000054422"
id <- "ENSMUSG00000064356"


## Hpx
sp1 <- plotPauseScores("ENSMUSG00000030895", codonsScore, cdsLenData,
                sampleAnno, contrasts = c("Type", "WT.Leu", "WT.SD")) +
    scale_color_brewer(palette="Set1") +
    labs(title = "Hpx",
         x = "Codon position",
         y = "Pause score",
         col = "") +
    theme_cowplot(16)
ggsave("figures/2022-05-09/ENSMUSG00000030895_Hpx_pauseScores.pdf",
       sp1, width = 15, height = 6)

sp1 <- plotPauseScores("ENSMUSG00000030895", codonsScore, cdsLenData,
                sampleAnno, contrasts = c("Type", "WT.Leu", "WT.SD")) +
    scale_color_brewer(palette="Set1") +
    labs(title = "Hpx",
         x = "Codon position",
         y = "Pause score",
         col = "") +
    xlim(300, 470) +
    theme_cowplot(16)
ggsave("figures/2022-05-09/ENSMUSG00000030895_Hpx_pauseScores_zoom.pdf",
       sp1, width = 15, height = 6)


## 
sp1 <- plotPauseScores("ENSMUSG00000074207", codonsScore, cdsLenData,
                sampleAnno, contrasts = c("Type", "WT.Leu", "WT.SD")) +
    scale_color_brewer(palette="Set1") +
    ## ylim(0, 7) +
    labs(title = "ENSMUSG00000074207",
         x = "Codon position",
         y = "Pause score",
         col = "") +
    theme_cowplot(16)
ggsave("figures/2022-05-09/ENSMUSG00000074207_NA_pauseScores.pdf",
       sp1, width = 15, height = 6)


## Alb
sp1 <- plotPauseScores("ENSMUSG00000029368", codonsScore, cdsLenData,
                sampleAnno, contrasts = c("Type", "WT.Leu", "WT.SD")) +
    scale_color_brewer(palette="Set1") +
    labs(title = "Alb",
         x = "Codon position",
         y = "Pause score",
         col = "") +
    theme_cowplot(16)
ggsave("figures/2022-05-09/ENSMUSG00000029368_Alb_pauseScores.pdf",
       sp1, width = 15, height = 6)

## Cat
sp1 <- plotPauseScores("ENSMUSG00000027187", codonsScore, cdsLenData,
                sampleAnno, contrasts = c("Type", "WT.Leu", "WT.SD")) +
    scale_color_brewer(palette="Set1") +
    labs(title = "Cat",
         x = "Codon position",
         y = "Pause score",
         col = "") +
    theme_cowplot(16)
ggsave("figures/2022-05-09/ENSMUSG00000027187_Cat_pauseScores.pdf",
       sp1, width = 15, height = 6)

## COX1
sp1 <- plotPauseScores("ENSMUSG00000064351", codonsScore, cdsLenData,
                sampleAnno, contrasts = c("Type", "WT.Leu", "WT.SD")) +
    scale_color_brewer(palette="Set1") +
    labs(title = "COX1",
         x = "Codon position",
         y = "Pause score",
         col = "") +
    theme_cowplot(16)
ggsave("figures/2022-05-09/ENSMUSG00000064351_COX1_pauseScores.pdf",
       sp1, width = 15, height = 6)

## Calr
sp1 <- plotPauseScores("ENSMUSG00000003814", codonsScore, cdsLenData,
                sampleAnno, contrasts = c("Type", "WT.Leu", "WT.SD")) +
    scale_color_brewer(palette="Set1") +
    labs(title = "Calr",
         x = "Codon position",
         y = "Pause score",
         col = "") +
    theme_cowplot(16)
ggsave("figures/2022-05-09/ENSMUSG00000003814_Calr_pauseScores.pdf",
       sp1, width = 15, height = 6)

## Gnmt
sp1 <- plotPauseScores("ENSMUSG00000002769", codonsScore, cdsLenData,
                sampleAnno, contrasts = c("Type", "WT.Leu", "WT.SD")) +
    scale_color_brewer(palette="Set1") +
    labs(title = "Gnmt",
         x = "Codon position",
         y = "Pause score",
         col = "") +
    theme_cowplot(16)
ggsave("figures/2022-05-09/ENSMUSG00000002769_Gnmt_pauseScores.pdf",
       sp1, width = 15, height = 6)

sp1 <- plotPauseScores("ENSMUSG00000002769", codonsScore, cdsLenData,
                sampleAnno, contrasts = c("Type", "WT.Leu", "WT.SD")) +
    scale_color_brewer(palette="Set1") +
    labs(title = "Gnmt",
         x = "Codon position",
         y = "Pause score",
         col = "") +
    xlim(195, 240) +
    theme_cowplot(16) +
    theme(legend.position = "none")
ggsave("figures/2022-05-09/ENSMUSG00000002769_Gnmt_pauseScores_zoom.pdf",
       sp1, width = 6, height = 6)


## Aldob
sp1 <- plotPauseScores("ENSMUSG00000028307", codonsScore, cdsLenData,
                sampleAnno, contrasts = c("Type", "WT.Leu", "WT.SD")) +
    scale_color_brewer(palette="Set1") +
    labs(title = "Aldob",
         x = "Codon position",
         y = "Pause score",
         col = "") +
    theme_cowplot(16)
ggsave("figures/2022-05-09/ENSMUSG00000028307_Aldob_pauseScores.pdf",
       sp1, width = 15, height = 6)


## Ambp
sp1 <- plotPauseScores("ENSMUSG00000028356", codonsScore, cdsLenData,
                sampleAnno, contrasts = c("Type", "WT.Leu", "WT.SD")) +
    scale_color_brewer(palette="Set1") +
    labs(title = "Ambp",
         x = "Codon position",
         y = "Pause score",
         col = "") +
    theme_cowplot(16)
ggsave("figures/2022-05-09/ENSMUSG00000028356_Ambp_pauseScores.pdf",
       sp1, width = 15, height = 6)

######################################################################
### Save image
## save.image("2022-05-09_findPausingSites.RData")
load("2022-05-09_findPausingSites.RData")
