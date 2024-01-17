library(ggplot2)
library(Seurat)
# library(Azimuth)
library(ggpubr)
library(gridExtra)
library(grid)
library(ggtext)
library(scales)

umap_dot_red_cols <- c("#eeecec", "#FB6A4A", "#67000D")
umap_dot_blue_cols <- c("#eeecec", "#4292C6", "#084594")

append_list <- function(x, y) {
    if (length(x) == 0) {
        x <- c(y)
    } else {
        x <- append(x, y)
    }
    return(x)
}

make_seurat_obj <- function(in_mtx_dir) {
    if (is.object(in_mtx_dir)) {
        return(in_mtx_dir)
    }
    if (!is.character(in_mtx_dir) || !dir.exists(in_mtx_dir)) {
        stop("gene_mtx_dir must be a Seurat object or a folder.")
    }
    rdata = paste(in_mtx_dir, "/seurat.RData", sep="")
    if (file.exists(rdata)) {
        load(rdata)
        print(paste(rdata, "loaded.", sep = " "))
    } else {
        print(paste(rdata, "not found, creating one from count matrix now.", sep = " "))
        # load data
        dt <- Read10X(in_mtx_dir)
        # create seurat object
        set.seed(123)
        gene_obj <- CreateSeuratObject(counts = dt, project = name, min.cells = 3, min.features = 200)
        gene_obj[["percent.mt"]] <- PercentageFeatureSet(gene_obj, pattern = "^MT-")
        gene_obj <- subset(gene_obj, subset = percent.mt < 25)

        gene_obj <- NormalizeData(gene_obj)
        gene_obj <- FindVariableFeatures(gene_obj, selection.method = "vst", nfeatures = 2000)
        gene_obj <- ScaleData(gene_obj)
        gene_obj <- RunPCA(gene_obj)
        gene_obj <- FindNeighbors(gene_obj, dims = 1:10)
        gene_obj <- FindClusters(gene_obj, resolution = 0.03, graph.name = "RNA_snn")  # change resolution for different data
        gene_obj <- RunUMAP(gene_obj, dims = 1:10, graph.name = "RNA_snn")

        save(gene_obj, file = rdata)
    }
    return(gene_obj)
}

feature_plot <- function(seurat_vec,
                         pt_size, cells, 
                         feature_list, feature_label_list,
                         cols, legend_pos, min1, max1, min2, max2,
                         full_idx) {
    pls <- list()
    vjust=0
    for (i in seq_along(feature_list)) {
        if (i == length(feature_list)) {
            legened_pos = legend_pos
            if (length(feature_list) > 3) {
                margin = unit(c(0.1, 0.8, 0, 0.2), "cm") # margin: t, r, b, l 
                
            } else {
                margin = unit(c(0.1, 1.3, 0, 0.2), "cm") # margin: t, r, b, l
            }
        } else if (i %in% full_idx) {
            legened_pos = "none"
            margin = unit(c(0.1, 0.3, 0, 0.2), "cm")
        } else {
            legened_pos = "none"
            margin = unit(c(0.1, 0.3, 0, 0.2), "cm")
        }
        seurat_obj = seurat_vec[[i]]
        feature = feature_list[[i]]
        label = feature_label_list[[i]]

        if (i %in% full_idx) {
             pl <- FeaturePlot(seurat_obj, pt.size=pt_size, cells=cells, 
                          features = feature,
                          keep.scale = NULL, #"all", #feature"
                          order = TRUE,
                          cols = cols) + 
                          labs(title=label) +
                          scale_colour_gradientn(colours = cols,
                                                breaks = c(1, 3),
                                                labels = c("Low", "High")) +
                          guides(colour = guide_colorbar(title.position = "left", 
                                                title.hjust = 0.5,
                                                direction = "vertical",
                                                barwidth = 0.5,
                                                barheight = 6)) +
                          theme(text=element_text(size=font_size, family=font_family), 
                          legend.position = legened_pos,
                          # legend.title = element_text(size=font_size, family=font_family, vjust=0.5, angle = 90),
                          plot.title = element_markdown(hjust = 0.5, vjust = vjust, size=font_size, face="plain"),
                          # plot.title = element_text(size=font_size, face="plain", vjust = vjust),
                          plot.margin = margin,
                          axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                          axis.text.y=element_blank(), axis.ticks.y=element_blank())
        } else {
            pl <- FeaturePlot(seurat_obj, pt.size=pt_size, cells=cells, 
                              features = feature,
                              keep.scale = NULL, #"all", #feature"
                              order = TRUE,
                              cols = cols) + 
                labs(title=label, y=NULL) +
                scale_colour_gradientn(colours = cols,
                                       breaks = c(1, 3),
                                       labels = c("Low", "High")) +
                guides(colour = guide_colorbar(title.position = "left", 
                                               title.hjust = 0.5,
                                               direction = "vertical",
                                               barwidth = 0.5,
                                               barheight = 6)) +
                # scale_x_continuous(limits=c(min1, max1)) + scale_y_continuous(limits = c(min2, max2)) +
                theme(text=element_text(size=font_size, family=font_family), 
                      legend.position = legened_pos,
                      plot.title = element_markdown(hjust = 0.5, vjust = vjust, size=font_size, face="plain"),
                      # plot.title = element_text(size=font_size, face="plain", vjust = vjust),
                      plot.margin = margin,
                      axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                      axis.text.y=element_blank(), axis.title.y=element_blank(),
                      axis.ticks.y=element_blank(), axis.line.y=element_blank())
        }
        pls[[length(pls)+1]] <- pl
    }
    return(pls)
}

# gene_mtx_dir: gene matrix folder generated by NanoHunter
# feature_mtx_dir: gene/transcript matrix folder generated by NanoHunter, either a folder or vector of folders with same length as feature_list/feature_label_list
# feature_list: list of gene name/transcript id, e.g. c("VIM", "ENST00000433892", "ENST00000263398")
# feature_label_list: list of gene/transcript label, e.g. c("Gene", "Trans1", "Trans2")
# nrow: number of rows in the combined UMAP plot (default: 1)
nanohunter_umap_plot <- function(gene_mtx_dir, 
                                 feature_mtx_dir, 
                                 feature_list, 
                                 feature_label_list=list(), 
                                 umap_col=umap_dot_red_cols, 
                                 text_pos="left", 
                                 text_rot=0, 
                                 nrow=1) {
    if (length(feature_label_list) == 0) {
        feature_label_list = feature_list
    } else if (length(feature_list) != length(feature_label_list)) {
        stop("feature_list must have same length as feature_label_list.")
    }
    feature_mtx_dirs <- list()
    if (length(feature_list) > 1 && length(feature_mtx_dir) == 1) {
        for (i in seq_along(feature_list)) {
            feature_mtx_dirs <- c(feature_mtx_dirs, feature_mtx_dir)
        }
    } else {
        feature_mtx_dirs <- c(feature_mtx_dir)
    }
    gene_obj <- make_seurat_obj(gene_mtx_dir)
    feature_obj_vec <- list()

    for (i in seq_along(feature_mtx_dirs)) {
        feature_obj <- make_seurat_obj(feature_mtx_dirs[[i]])
        common_cells = intersect(colnames(gene_obj), colnames(feature_obj))
        feature_obj@reductions$umap@cell.embeddings[common_cells, ] <- gene_obj@reductions$umap@cell.embeddings[common_cells, ]
        feature_obj_vec <- append_list(feature_obj_vec, feature_obj)
    }
    # min_exp=0
    # g <- FetchData(object = gene_obj, vars = gene)
    # common_cells <- intersect(common_cells, rownames(g)[g > min_exp])

    min1 <- min(gene_obj@reductions$umap@cell.embeddings[, 1]) - 0.1
    min2 <- min(gene_obj@reductions$umap@cell.embeddings[, 2]) - 0.1
    max1 <- max(gene_obj@reductions$umap@cell.embeddings[, 1]) + 0.1
    max2 <- max(gene_obj@reductions$umap@cell.embeddings[, 2]) + 0.1

    umap_legend_pos <- c(0.9, 0.4)
    full_idx <- c(1)
    if (nrow > 1) {
        s = 1
        unit_len = ceiling(length(feature_list)/nrow)
        for (i in 2:nrow) {
            full_idx <- c(full_idx, s + unit_len*(i-1))
        }
    }
    features_plts <- feature_plot(feature_obj_vec,
                                  umap_pt_size, common_cells, 
                                  feature_list, feature_label_list,
                                  umap_col,
                                  #c("lightgrey", "red"), #,  #eeecec
                                  umap_legend_pos, min1, max1, min2, max2,
                                  full_idx)
    # isoform plot
    pls = list()# (gene_plt)
    for (i in seq_along(features_plts)) {
        pls[[1+length(pls)]] <- features_plts[[i]]
    }
    if (length(pls) > 3) {
        if (nrow == 1) {
            gt <- arrangeGrob(grobs = pls, nrow = 1, 
                              widths = c(4.5, rep(4, length(feature_list)-2), 4.5))
        } else {
            gt <- arrangeGrob(grobs = pls, nrow=nrow)
        }
    } else if (length(pls) == 1) {
        gt <- arrangeGrob(grobs = pls, nrow=nrow)
    } else {
        gt <- arrangeGrob(grobs = pls, nrow = 1, widths = c(4.5, rep(4, length(feature_list)-2), 5))
    }
    
    p <- as_ggplot(gt)
    return(p)
}