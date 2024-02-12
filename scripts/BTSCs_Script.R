# Loading packages.
library(fusca)
library(dplyr)
library(Matrix)

#  Loading GSC raw data.
#  >69k unfiltered GSCs.
GSC_RawData <- read.csv(
    "data/Richards_NatureCancer_GSC_scRNAseq_counts.csv",
    header = TRUE,
    row.names = 1,
    check.names = FALSE
)

# Here we have GSCs classified as "Adherent" or "Sphere".
GSC_Metadata <- read.csv(
    "data/Richards_NatureCancer_GSC_scRNAseq_meta.csv",
    sep = ";",
    row.names = 1,
    header = TRUE
)
# Subset adherent cells.
GSC_Adherent <- GSC_Metadata[GSC_Metadata$CultureMethod == "Adherent",]
# Subset cells cultured as spheres.
GSC_Sphere <- GSC_Metadata[GSC_Metadata$CultureMethod == "Sphere",]
 
# Here we have >65k GSCs post-QC and filtering (done by the authors of the original publication).
GSC_Tumor_Metadata <- read.csv(
    "data/GSC_Tumor_Meta.csv",
    header = TRUE,
    row.names = 1,
    sep = ";",
    check.names = FALSE
)
# Keep only GSCs, excluding whole-tumour cells.
GSC_Metadata_Filtered <- GSC_Metadata_Filtered[GSC_Metadata_Filtered$Sample.Type == "GSC",] 
# Remove "BTSC_" from cell names so that they can be matched with the names in the other table ("GSC_RawData").
rownames(GSC_Tumor_Metadata) <- gsub('BTSC_', '', rownames(GSC_Tumor_Metadata)) 
 
# Filter raw count matrix, keeping only cells with good quality (as identified by the authors of the original publication).
GSC_RawData_Filtered <- GSC_RawData %>% select(rownames(GSC_Metadata_Filtered))

# Get count matrices for adherent cells and spheres separately.
GSC_RawData_Adherent <- GSC_RawData_Filtered %>% select(intersect(rownames(GSC_Adherent), colnames(GSC_RawData_Filtered)))
GSC_RawData_Sphere <- GSC_RawData_Filtered %>% select(intersect(rownames(GSC_Sphere), colnames(GSC_RawData_Filtered)))

# Creating CellRouter object for adherent cells.
adherent <- CreateCellRouter(
    as.matrix(GSC_RawData_Adherent),
    min.genes = 0,
    min.cells = 0,
    is.expr = 0
)
adherent@assays$RNA@rawdata <- adherent@assays$RNA@rawdata[rownames(adherent@assays$RNA@ndata), colnames(adherent@assays$RNA@ndata)]

# Normalize & scale data according to all genes.
adherent <- Normalize(adherent)
adherent <- scaleData(
    adherent, 
    blocksize = nrow(adherent@assays$RNA@ndata)
)
# Compute PCA, UMAP & tSNE.
adherent <- computePCA(adherent, num.pcs = 50, seed = 42)
adherent <- computeUMAP(adherent, num.pcs = 15, seed = 42)
adherent <- computeTSNE(adherent, num.pcs = 15, seed = 42)

# Find clusters.
adherent <- findClusters(
    object = adherent,
    assay.type = "RNA",
    method = "graph.clustering",
    k = 20,
    num.pcs = 15,
    sim.type = "jaccard",
    nn.type = "knn"
)

dir.create("results/adherent/")
save(adherent, file = "results/adherent/Adherent_Normalized_Scaled_PCA_UMAP_tSNE_CellRouter_object_with_Clusters.R")


plotReducedDimension(
    adherent,
    assay.type = "RNA",
    reduction.type = "umap",
    dims.use = c(1,2),
    annotation = "population",
    annotation.color = "colors",
    showlabels = TRUE,
    width = 4.5,
    height = 3.5,
    convex = FALSE,
    dotsize = 0.01,
    filename = "results/adherent/UMAP_CellRouter_Clusters.pdf"
)
dev.off()
plotReducedDimension(
    adherent,
    assay.type = "RNA",
    reduction.type = "tsne",
    dims.use = c(1,2),
    annotation = "population",
    annotation.color = "colors",
    showlabels = TRUE,
    width = 4.5,
    height = 3.5,
    convex = FALSE,
    dotsize = 0.01,
    filename = "results/adherent/tSNE_CellRouter_Clusters.pdf"
)
dev.off()
 
adherent_markers <- findSignatures(
    adherent,
    assay.type = "RNA",
    column = "population",
    test.use = "wilcox",
    pos.only = "TRUE"
)
write.csv(adherent_markers, file = "results/adherent/Adherent_markers_CellRouter_clusters.csv")
  
prnp_exp <- adherent@assays$RNA@ndata["PRNP",]
prnp_exp <- as.data.frame(prnp_exp)
colnames(prnp_exp) <- "PRNP"
prnp_exp <- prnp_exp %>% mutate(classification = case_when(PRNP == 0 ~ "PRNP_negative_cells", PRNP > 0 ~ "PRNP_positive_cells"))
 
adherent <- addInfo(
    adherent, 
    metadata = prnp_exp, 
    colname = "PRNP_status", 
    metadata.column = "classification"
)
print("Added PRNP classification to CellRouter object.")
 
plotDRExpression(adherent,
                  genelist = c("PRNP"),
                  title = "",
                 width = 4.5,
                 height = 3.5,
                  threshold = 3,
                  reduction.type = "umap",
                  dotsize = 0.01,
                  filename = "BTSCs/adherent/Adherent_UMAP_DREexpression_PRNP.pdf")
dev.off()
print("Plotted UMAP displaying PRNP expression across cells.")
 
plotDRExpression(adherent,
                  genelist = c("PRNP"),
                  title = "",
                  width = 2,
                  height = 2.2,
                  threshold = 3,
                  reduction.type = "tsne",
                  dotsize = 0.01,
                  filename = "BTSCs/adherent/Adherent_tSNE_DREexpression_PRNP.pdf")

plotViolin(adherent,
           assay.type = "RNA",
           geneList = "PRNP",
           column = "population",
           column.color = "colors",
           cols = 1,
           width = 8,
           height = 3.5,
           filename = "BTSCs/adherent/Adherent_tSNE_VlnPlot_PRNP.pdf")

dev.off()
print("Plotted tSNE displaying PRNP expression across cells.")
 
adherent.markers.PRNP.levels <- findSignatures(object = adherent, 
                                                assay.type = "RNA",
                                                column = "PRNP_status", 
                                                test.use = "wilcox")
print("Found PRNP+ and PRNP- marker gene signatures.")
 
write.csv(adherent.markers.PRNP.levels, 
         "BTSCs/adherent/Adherent_PRNP+_vs_PRNP-_signature_markers.csv")

adherent <- scoreGeneSets(adherent, assay.type = "RNA", genes.list = LIST)

plotSignatureScore(adherent,
                   assay.type = "RNA",
                   names(LIST), 
                   reduction.type='umap', 
                   columns = 2,
                   width=7, 
                   height=7, 
                   dotsize=0.5, 
                   filename='BTSCs/adherent/signatures_umap.pdf')
 
#print("Starting correlation analysis.")
# df_adherent <- data.frame(matrix(0, ncol=4))
#  colnames(df_adherent) <- c('gene1', 'gene2', 'correlation', 'pvalue')
#  i=1
#  for(g in rownames(adherent@assays$RNA@ndata)){
#    c <- cor.test(as.numeric(adherent@assays$RNA@ndata[g,]), as.numeric(adherent@assays$RNA@ndata['PRNP',]))
#    df_adherent[i,'gene1'] <- 'PRNP'
#    df_adherent[i,'gene2'] <- g
#    df_adherent[i,'correlation'] <- c$estimate stores the correlation
#    df_adherent[i,'pvalue'] <- c$p.value stores the correlation
#    i <- i + 1
#  }
#  print("Computed correlation between PRNP and all genes.")
#  
#  write.csv(df_adherent, "BTSCs/adherent/Adherent_Correlation_PRNPvsAllGenes.csv")
#  
print("Finished the analysis of Adherent cultures.")

################################################################################

print("Starting the analysis of Sphere cultures.")
 
# Creating CellRouter object.
Sphere <- CreateCellRouter(GSC_RawData_Sphere_Sparse,
                            min.genes = 0,
                            min.cells = 0,
                            is.expr = 0)
 
print("Created CellRouter object.")
 
Sphere@assays$RNA@rawdata <- Sphere@assays$RNA@rawdata[rownames(Sphere@assays$RNA@ndata), colnames(Sphere@assays$RNA@ndata)]
 
# NOTE: no QC filtering will be performed.
 
Sphere <- Normalize(Sphere)
print("Normalized data.")
 
Sphere <- scaleData(Sphere, 
                   blocksize = nrow(Sphere@assays$RNA@ndata))
print("Scaled data according to all genes.")
 
Sphere <- computePCA(Sphere, num.pcs = 50, seed = 42)
pdf(file = 'BTSCs/Sphere/Sphere_elbow_plot.pdf', width = 5, height = 5)
plot(Sphere@pca$sdev)
dev.off()
print("Computed PCA.")
 
Sphere <- computeUMAP(Sphere, num.pcs = 15, seed = 42)
print("Computed UMAP.")
 
Sphere <- computeTSNE(Sphere, num.pcs = 15, seed = 42)
print("Computed tSNE.")
 
save(Sphere, 
    file = "BTSCs/Sphere/Sphere_Normalized_Scaled_PCA_UMAP_tSNE_CellRouter_object.R")

Sphere <- findClusters(Sphere,
                        assay.type = "RNA",
                        method = "graph.clustering",
                        k = 20,
                        num.pcs = 15,
                        sim.type = "jaccard",
                        nn.type = "knn")

save(Sphere, 
    file = "BTSCs/Sphere/Sphere_Normalized_Scaled_PCA_UMAP_tSNE_CellRouter_object_with_Clusters.R")

plotReducedDimension(Sphere,
                      assay.type = "RNA",
                      reduction.type = "umap",
                      dims.use = c(1,2),
                      annotation = "population",
                      annotation.color = "colors",
                      showlabels = TRUE,
                      width = 4.5,
                      height = 3.5,
                      convex = FALSE,
                      dotsize = 0.01,
                      filename = "BTSCs/Sphere/UMAP_CellRouter_Clusters.pdf")
dev.off()
plotReducedDimension(Sphere,
                      assay.type = "RNA",
                      reduction.type = "tsne",
                      dims.use = c(1,2),
                      annotation = "population",
                      annotation.color = "colors",
                      showlabels = TRUE,
                      width = 4.5,
                      height = 3.5,
                      convex = FALSE,
                      dotsize = 0.01,
                      filename = "BTSCs/Sphere/tSNE_CellRouter_Clusters.pdf")
dev.off()

plotViolin(sphere,
           assay.type = "RNA",
           geneList = "PRNP",
           column = "population",
           column.color = "colors",
           cols = 1,
           width = 8,
           height = 3.5,
           filename = "BTSCs/Sphere/Sphere_VlnPlot_PRNP.pdf")

Sphere_markers <- findSignatures(Sphere,
                                  assay.type = "RNA",
                                  column = "population",
                                  test.use = "wilcox",
                                  pos.only = "TRUE")
write.csv(Sphere_markers,
         file = "BTSCs/Sphere/Sphere_markers_CellRouter_clusters.csv")
 
prnp_exp <- Sphere@assays$RNA@ndata["PRNP",]
prnp_exp <- as.data.frame(prnp_exp)
colnames(prnp_exp) <- "PRNP"
prnp_exp <- prnp_exp %>% mutate(classification = case_when(PRNP == 0 ~ "PRNP_negative_cells",
                                                           PRNP > 0 ~ "PRNP_positive_cells"))

write.csv(prnp_exp,
          file = "BTSCs/Sphere/Sphere_PRNP_counts.csv")

Sphere <- addInfo(Sphere, 
                    metadata = prnp_exp, 
                    colname = "PRNP_status", 
                    metadata.column = "classification")
print("Added PRNP classification to CellRouter object.")

plotDRExpression(Sphere,
                 genelist = c("PRNP"),
                 title = "",
                 width = 2,
                 height = 2.2,
                 threshold = 2,
                 reduction.type = "umap",
                 dotsize = 0.01,
                 filename = "BTSCs/Sphere/Sphere_UMAP_DREexpression_PRNP.pdf")
dev.off()
print("Plotted UMAP displaying PRNP expression across cells.")

plotDRExpression(Sphere,
                 genelist = c("PRNP"),
                 title = "",
                 width = 2,
                 height = 2.2,
                 threshold = 2,
                 reduction.type = "tsne",
                 dotsize = 0.01,
                 filename = "BTSCs/Sphere/Sphere_tSNE_DREexpression_PRNP.pdf")
dev.off()
print("Plotted tSNE displaying PRNP expression across cells.")

Sphere.markers.PRNP.levels <- findSignatures(object = Sphere, 
                                               assay.type = "RNA",
                                               column = "PRNP_status", 
                                               test.use = "wilcox")
print("Found PRNP+ and PRNP- marker gene signatures.")

write.csv(Sphere.markers.PRNP.levels, 
          "BTSCs/Sphere/Sphere_PRNP+_vs_PRNP-_signature_markers.csv")

# print("Starting correlation analysis.")
# df_Sphere <- data.frame(matrix(0, ncol=4))
# colnames(df_Sphere) <- c('gene1', 'gene2', 'correlation', 'pvalue')
# i=1
# for(g in rownames(Sphere@assays$RNA@ndata)){
#   c <- cor.test(as.numeric(Sphere@assays$RNA@ndata[g,]), as.numeric(Sphere@assays$RNA@ndata['PRNP',]))
#   df_Sphere[i,'gene1'] <- 'PRNP'
#   df_Sphere[i,'gene2'] <- g
#   df_Sphere[i,'correlation'] <- c$estimate stores the correlation
#   df_Sphere[i,'pvalue'] <- c$p.value stores the correlation
#   i <- i + 1
# }
# print("Computed correlation between PRNP and all genes.")

#write.csv(df_Sphere, "BTSCs/Sphere/Sphere_Correlation_PRNPvsAllGenes.csv")

print("Finished the analysis of Sphere cultures.")

sphere <- scoreGeneSets(sphere, assay.type = "RNA", genes.list = LIST)

plotSignatureScore(sphere,
                   assay.type = "RNA",
                   names(LIST), 
                   reduction.type='umap', 
                   columns = 2,
                   width=7, 
                   height=7, 
                   dotsize=0.5, 
                   filename='BTSCs/Sphere/signatures_umap.pdf')

plotSignaturesBoxplot(sphere,
                      assay.type = "RNA",
                      names(LIST), 
                      width=7, 
                      height=7, 
                      filename='BTSCs/Sphere/signatures_boxplot.pdf')
write.csv(sphere@assays$RNA@sampTab, "BTSCs/Sphere/Sphere_metadata_signatures.csv")

