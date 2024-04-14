
## packages are required for differential methylation analysis. 
## Additionally, I attached tutorial "A cross-package Bioconductorworkfl ow for analysingmethylation array data"
## Authors: Jovana Maksimovic,Belinda Phipson and AliciaOshlack
## which I follow for methylation data analysis:

library(tidyverse)
library(lumi)  
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(compEpiTools)
library(org.Hs.eg.db)

load("van_Methylation.Rda")

methylation_sample_annotation  <- read.csv("methylation_sample_annotation.csv")

vanSampleInfo_PDX2   <- tail(methylation_sample_annotation,8)
vanSampleInfo_PDX2_f2   <- vanSampleInfo_PDX2[-c(4),] ## Only contain control and f2_aza treated


vanMethData_PDX2  <- vanMethData %>% select(4,8,9,11:15)
vanMethData_PDX2_f2  <- vanMethData %>% select(4,8,9,12:15)   ## Only contain control and f2_aza treated



beta_value_Matrix_vanMethData_PDX2 <- as.matrix(vanMethData_PDX2)
typeof(vanMethData_PDX2)
head(vanMethData_PDX2)
typeof(beta_value_Matrix_vanMethData_PDX2)
head(beta_value_Matrix_vanMethData_PDX2)
M_value_Matrix_vanMethData_PDX2 <- beta2m(beta_value_Matrix_vanMethData_PDX2)


##### only contain f2_aza and control
beta_value_Matrix_vanMethData_PDX2_f2 <- as.matrix(vanMethData_PDX2_f2)
typeof(vanMethData_PDX2_f2)
head(vanMethData_PDX2_f2)
typeof(beta_value_Matrix_vanMethData_PDX2_f2)
head(beta_value_Matrix_vanMethData_PDX2_f2)
M_value_Matrix_vanMethData_PDX2_f2 <- beta2m(beta_value_Matrix_vanMethData_PDX2_f2)

write.table(vanSampleInfo_PDX2_f2, file="vanSampleInfo_PDX2_f2.csv", sep=",", row.names=FALSE)
write.table(beta_value_Matrix_vanMethData_PDX2_f2, file="beta_value_Matrix_vanMethData_PDX2_f2.csv", sep=",", row.names=TRUE)
write.table(M_value_Matrix_vanMethData_PDX2_f2, file="M_value_Matrix_vanMethData_PDX2_f2.csv", sep=",", row.names=TRUE)

#####


# picture of beta_value and M_value

pdf("PDX2_study.pdf") 
par(mfrow=c(1,2))
densityPlot(beta_value_Matrix_vanMethData_PDX2, sampGroups=vanSampleInfo_PDX2$Treatment, main="Beta values", legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(vanSampleInfo_PDX2$Treatment)),text.col=brewer.pal(8,"Dark2"))
densityPlot(M_value_Matrix_vanMethData_PDX2, sampGroups=vanSampleInfo_PDX2$Treatment, main="M-values",legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(vanSampleInfo_PDX2$Treatment)), text.col=brewer.pal(8,"Dark2"))
dev.off() 


##### only contain f2_aza and control
pdf("PDX2_f2_study.pdf") 
par(mfrow=c(1,2))
densityPlot(beta_value_Matrix_vanMethData_PDX2_f2, sampGroups=vanSampleInfo_PDX2_f2$Treatment, main="Beta values", legend=FALSE, xlab="Beta values")
legend("topright", legend = levels(factor(vanSampleInfo_PDX2_f2$Treatment)),text.col=brewer.pal(8,"Dark2"))
densityPlot(M_value_Matrix_vanMethData_PDX2_f2, sampGroups=vanSampleInfo_PDX2_f2$Treatment, main="M-values",legend=FALSE, xlab="M values")
legend("topright", legend = levels(factor(vanSampleInfo_PDX2_f2$Treatment)), text.col=brewer.pal(8,"Dark2"))
dev.off() 
####



cellType <- factor(vanSampleInfo_PDX2$Treatment)
individual <- factor(vanSampleInfo_PDX2$Model)

design <- model.matrix(~0+cellType, data=vanSampleInfo_PDX2) 
colnames(design) <- c(levels(cellType),levels(individual)[-1])

fit <- lmFit(M_value_Matrix_vanMethData_PDX2, design)
contMatrix <- makeContrasts(control-f1_aza,control-f2_aza,levels=design)
contMatrix

fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)


# 1 for upregulated (as a control), -1 for downregulated (act as case) and 0 for not significant corresponding to probes.
write.table(decideTests(fit2), file="PDX2_result_Significant_differential_expression.csv", sep="\t", row.names=TRUE)

summary(decideTests(fit2))

fullannotSub <- fullannot[match(rownames(M_value_Matrix_vanMethData_PDX2),fullannot[, "NAME"]),c(1:ncol(fullannot))]

# "control - f2_aza" can be saved as a data.frame by setting coef=2. Here, I took forward only "control versus f2_aza" for analysis.


DMPs <- topTable(fit2, num=Inf, coef=2, genelist=fullannotSub)

head(DMPs)

write.table(DMPs, file="result_control-f2_aza_DMPs_PDX2.csv", sep=",", row.names=FALSE)


pdf("result_control-f2_aza_PDX2.pdf") 
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){plotCpg(vanMethData_PDX2, cpg=cpg, pheno=vanSampleInfo_PDX2$Treatment, ylab = "Beta values")})
dev.off() 


## we need to remove those probes which contains "NA" values.
modified_M_value_Matrix_vanMethData_PDX2 <- na.omit(M_value_Matrix_vanMethData_PDX2)
modified_beta_value_Matrix_vanMethData_PDX2 <- na.omit(beta_value_Matrix_vanMethData_PDX2)


myAnnotation <- cpg.annotate(object = modified_M_value_Matrix_vanMethData_PDX2, datatype = "array", what = "M", analysis.type = "differential", design = design, contrasts = TRUE, cont.matrix = contMatrix, coef ="control - f2_aza", arraytype = "450K")

str(myAnnotation)

DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
results.ranges


write.table(results.ranges, file="result_control-f2_aza_DM_region_PDX2.csv", sep=",", row.names=FALSE)


########## Here, below I deal with only those probes which has lost methylation after treatment.


system("head PDX2_result_Significant_differential_expression.csv")


system("awk 'NR>1 && $3==1{print $1}' PDX2_result_Significant_differential_expression.csv > keep_PDX2_f2_aza")


keep_PDX2_f2_aza  <- read.csv("keep_PDX2_f2_aza",header=FALSE)


# filter probes having methylation loss after treatment:

result_control_f2_aza_DMPs_PDX2_methylation_loss  <- subset(DMPs, rownames(DMPs) %in% keep_PDX2_f2_aza$V1)

write.table(result_control_f2_aza_DMPs_PDX2_methylation_loss, file="result_control_f2_aza_DMPs_PDX2_methylation_loss.csv", sep=",", row.names=FALSE)


modified_M_value_Matrix_vanMethData_PDX2_methylation_loss  <- subset(modified_M_value_Matrix_vanMethData_PDX2, rownames(modified_M_value_Matrix_vanMethData_PDX2) %in% keep_PDX2_f2_aza$V1)


modified_beta_value_Matrix_vanMethData_PDX2_methylation_loss  <- subset(modified_beta_value_Matrix_vanMethData_PDX2, rownames(modified_beta_value_Matrix_vanMethData_PDX2) %in% keep_PDX2_f2_aza$V1)


fullannotSub_methylation_loss <- fullannot[match(rownames(modified_M_value_Matrix_vanMethData_PDX2_methylation_loss),fullannot[, "NAME"]),c(1:ncol(fullannot))]


myAnnotation_methylation_loss <- cpg.annotate(object = modified_M_value_Matrix_vanMethData_PDX2_methylation_loss, datatype = "array", what = "M", analysis.type = "differential", design = design, contrasts = TRUE, cont.matrix = contMatrix, coef ="control - f2_aza", arraytype = "450K")

str(myAnnotation_methylation_loss)

DMRs_methylation_loss <- dmrcate(myAnnotation_methylation_loss, lambda=1000, C=2)
results.ranges_methylation_loss <- extractRanges(DMRs_methylation_loss)
results.ranges_methylation_loss

write.table(results.ranges_methylation_loss, file="result_control-f2_aza_DM_region_PDX2_methylation_loss.csv", sep=",", row.names=FALSE)


##-- Genomic location of DMRs

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
pdf("ratio_of_promoter_intergeneic_interagenic_region.pdf")
par(mfrow=c(1,1))
DMRs_Simpleanno <- GRannotateSimple(results.ranges,txdb=txdb, upstream=2000, downstream=1000)
dev.off()


pdf("ratio_of_promoter_intergeneic_interagenic_region_after_methylation_loss.pdf")
par(mfrow=c(1,1))
DMRs_Simpleanno <- GRannotateSimple(results.ranges_methylation_loss,txdb=txdb, upstream=2000, downstream=1000)
dev.off()



