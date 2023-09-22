library(Seurat)
library(Azimuth)
library(ggplot2)



argv <- commandArgs(trailingOnly = TRUE)

data_dir <- argv[1]
name <- argv[2]
out_dir <- argv[3]

make_azimuth_seurat_obj <- function(in_mtx_dir, name) {
    rdata_fn=paste(in_mtx_dir, "/", name, "_pbmc_gene.RData", sep="")
    if (file.exists(rdata_fn)) {
        load(rdata_fn)
        print(paste0(rdata_fn, " loaded."))
    } else {
        print(paste0(rdata_fn, " not found, creating one from count matrix now."))

        # load data
        dt <- Read10X(in_mtx_dir)
        # create seurat object
        pbmc_gene <- CreateSeuratObject(counts = dt, project = paste0(name, "_PBMC_Gene", sep=""), min.cells = 3, min.features = 200)
        pbmc_gene[["percent.mt"]] <- PercentageFeatureSet(pbmc_gene, pattern = "^MT-")
        pbmc_gene <- subset(pbmc_gene, subset = percent.mt < 25)

        # run azimuth
        pbmc_gene <- RunAzimuth(pbmc_gene, reference = "pbmcref")

        # write cell cluster assignments to file
        # bc_to_clu_file <- paste(in_mtx_dir, "/", name, "_azimuth_bc_to_cluster.tsv", sep = "")
        # pbmc_gene$predicted.celltype.l1<-gsub(" ", "", pbmc_gene$predicted.celltype.l1)
        # write.table(pbmc_gene$predicted.celltype.l1, bc_to_clu_file, sep = "\t", quote=F, col.names=F, row.names=T)
        save(pbmc_gene, file = rdata_fn)
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
l1_bc_to_clu_file <- paste(out_dir, "/", name, "_azimuth_bc_to_cell_type_l1.tsv", sep = "")
# pbmc$predicted.celltype.l1<-gsub(" ", "", pbmc$predicted.celltype.l1)
write.table(pbmc$predicted.celltype.l1, l1_bc_to_clu_file, sep = "\t", quote=F, col.names=F, row.names=T)

l2_bc_to_clu_file <- paste(out_dir, "/", name, "_azimuth_bc_to_cell_type_l2.tsv", sep = "")
# pbmc$predicted.celltype.l2<-gsub(" ", "", pbmc$predicted.celltype.l2)
write.table(pbmc$predicted.celltype.l2, l2_bc_to_clu_file, sep = "\t", quote=F, col.names=F, row.names=T)
