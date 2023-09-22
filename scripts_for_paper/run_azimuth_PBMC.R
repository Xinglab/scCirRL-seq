library(ggplot2)
library(Seurat)
library(Azimuth)
library(Seurat)


argv <- commandArgs(trailingOnly = TRUE)

data_dir <- argv[1]
name <- argv[2]
out_dir <- data_dir

make_azimuth_seurat_obj <- function(in_mtx_dir, name) {
    if (file.exists(paste(in_mtx_dir, "/pbmc_gene.RData", sep=""))) {
        load(paste(in_mtx_dir, "/pbmc_gene.RData", sep=""))
        print("pbmc_gene.RData loaded.")
    } else {
        print("pbmc_gene.RData not found, creating one from count matrix now.")

        # load data
        dt <- Read10X(in_mtx_dir)
        # create seurat object
        pbmc_gene <- CreateSeuratObject(counts = dt, project = "PBMC_Gene", min.cells = 3, min.features = 200)
        pbmc_gene[["percent.mt"]] <- PercentageFeatureSet(pbmc_gene, pattern = "^MT-")
        pbmc_gene <- subset(pbmc_gene, subset = percent.mt < 25)

        # run azimuth
        pbmc_gene <- RunAzimuth(pbmc_gene, reference = "pbmcref")

        # write cell cluster assignments to file
        # bc_to_clu_file <- paste(in_mtx_dir, "/", name, "_azimuth_bc_to_cluster.tsv", sep = "")
        # pbmc_gene$predicted.celltype.l1<-gsub(" ", "", pbmc_gene$predicted.celltype.l1)
        # write.table(pbmc_gene$predicted.celltype.l1, bc_to_clu_file, sep = "\t", quote=F, col.names=F, row.names=T)
        # save(pbmc_gene, file = paste(in_mtx_dir, "/pbmc_gene.RData", sep=""))
    }
    return(pbmc_gene)
}

# load data
pbmc <- make_azimuth_seurat_obj(data_dir, name)

# umap plot
l1_umap_plot_file <- paste(out_dir, "/", name, "_azimuth_umap.l1.pdf", sep = "")
pdf(l1_umap_plot_file, width = 7, height = 7)
DimPlot(pbmc, group.by = "predicted.celltype.l1", label = TRUE, label.size = 6, pt.size = 0.5) + 
        NoLegend() +
        theme(text = element_text(size = 18, family = "sans"),
            axis.text=element_text(size = 18, family = "sans"),
            plot.title = element_blank())
dev.off()
l2_umap_plot_file <- paste(out_dir, "/", name, "_azimuth_umap.l2.pdf", sep = "")
pdf(l2_umap_plot_file, width = 7, height = 7)
DimPlot(pbmc, group.by = "predicted.celltype.l2", label = TRUE, label.size = 6, pt.size = 0.5) + 
        NoLegend() +
        theme(text = element_text(size = 18, family = "sans"),
            axis.text=element_text(size = 18, family = "sans"),
            plot.title = element_blank())
dev.off()

# write cell cluster assignments to file
l1_bc_to_clu_file <- paste(out_dir, "/", name, "_azimuth_bc_to_cluster_l1.tsv", sep = "")
# pbmc$predicted.celltype.l1<-gsub(" ", "", pbmc$predicted.celltype.l1)
write.table(pbmc$predicted.celltype.l1, l1_bc_to_clu_file, sep = "\t", quote=F, col.names=F, row.names=T)

l2_bc_to_clu_file <- paste(out_dir, "/", name, "_azimuth_bc_to_cluster_l2.tsv", sep = "")
# pbmc$predicted.celltype.l2<-gsub(" ", "", pbmc$predicted.celltype.l2)
write.table(pbmc$predicted.celltype.l2, l2_bc_to_clu_file, sep = "\t", quote=F, col.names=F, row.names=T)
