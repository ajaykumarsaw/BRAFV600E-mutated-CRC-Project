
library("gplots")
load("van_Methylation.Rda")

annotation = read.csv(file = "methylation_sample_annotation.csv",header = TRUE)
annotation1 <- annotation[-c(1:6,10),]
current_PDX2_f2 <- vanMethData[,c(4,8,9,12:15)]
colnames(current_PDX2_f2) <- c("a_Control_1", "a_Control_2", "a_Control_3", "b_f2_Aza_1", "b_f2_Aza_2", "b_f2_Aza_3", "b_f2_Aza_4" )
PDX_data <- current_PDX2_f2
## sort with respect to "b_f2_Aza_1" sample
sort_PDX_data <- PDX_data[order(PDX_data$b_f2_Aza_1, decreasing = TRUE),]
dataNum <- matrix(data = NA, nrow = dim(sort_PDX_data)[1], ncol = dim(sort_PDX_data)[2])
colnames(dataNum) <- colnames(dataNum, do.NULL = F, prefix = "col_")
rownames(dataNum) <- rownames(dataNum, do.NULL = F, prefix = "row_")
for (i in 1:dim(sort_PDX_data)[2]){dataNum[,i] <- c(as.numeric(sort_PDX_data[[i]]))}
colnames(dataNum) <- colnames(sort_PDX_data)
rownames(dataNum) <- rownames(sort_PDX_data)
top_50000 <- head(dataNum, 50000)
pdf("PDX2_control_aza_Heatmap_50K_decreasing_order.pdf")
heatmap.2(x=top_50000,
    dendrogram="col",
    scale="col",
    col="bluered",
    trace="none",
    labRow=FALSE,
    main="Heatmap")
    dev.off()
