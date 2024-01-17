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
        x <- y
    } else {
        x <- append(x, y)
    }
    return(x)
}

make_gene_seurat_obj <- function(in_mtx_dir) {
    if (is.object(in_mtx_dir)) {
        return(in_mtx_dir)
    }
    if (!is.character(in_mtx_dir) || !dir.exists(in_mtx_dir)) {
        stop("gene_mtx_dir must be a Seurat object or a folder.")
    }
    rdata = paste(in_mtx_dir, "/gene.RData", sep="")
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

make_trans_seurat_obj <- function(in_mtx_dir) {
    print(in_mtx_dir)
    if (is.object(in_mtx_dir)) {
        return(in_mtx_dir)
    }
    if (!is.character(in_mtx_dir) || !dir.exists(in_mtx_dir)) {
        stop("trans_mtx_dir must be a Seurat object or a folder.")
    }
    rdata = paste(in_mtx_dir, "/transcript.RData", sep="")
    if (file.exists(rdata)) {
        load(rdata)
        print(paste(rdata, "loaded.", sep=" "))
    } else {
        print(paste(rdata, "not found, creating one from count matrix now.", sep=" "))
        min_features <- 0 #500 for old plots
        # load the matrix
        counts_dt <- Read10X(data.dir = in_mtx_dir)
        # create seurat object
        set.seed(123)
        trans_obj <- CreateSeuratObject(counts = counts_dt, project = "transripts")
        trans_obj[["percent.mt"]] <- PercentageFeatureSet(trans_obj, pattern = "^MT-")

        trans_obj <- subset(trans_obj, subset = nFeature_RNA > min_features & percent.mt < 25)
        trans_obj <- NormalizeData(trans_obj)
        trans_obj <- FindVariableFeatures(trans_obj, selection.method = "vst", nfeatures = 2000)
        all.trans <- rownames(trans_obj)
        trans_obj <- ScaleData(trans_obj, features = all.trans)
        trans_obj <- RunPCA(trans_obj)
        trans_obj <- FindNeighbors(trans_obj, dims = 1:10)
        trans_obj <- FindClusters(trans_obj, resolution = 0.03, graph.name = "RNA_snn")  # change resolution for different data
        trans_obj <- RunUMAP(trans_obj, dims = 1:10, graph.name = "RNA_snn")

        save(trans_obj, file = rdata)
    }
    return(trans_obj)
}

gene_feature_plot <- function(Gene_seurat, pt_size, gene, cols, umap_legend_pos, min1, max1, min2, max2) {
    margin = unit(c(0,0.3,0,0), "cm")
    FeaturePlot(Gene_seurat, pt.size=pt_size, 
                #cells=cells, 
                features = gene, 
                # keep.scale = "all", #NULL,
                order = TRUE,
                label = FALSE, label.size = umap_label_size-2, label.color = "black", #label.color = labels_cols,
                repel=TRUE,
                cols = cols) + 
                # scale_x_continuous(limits=c(min1, max1)) + scale_y_continuous(limits = c(min2, max2)) +
                # labs(y=NULL) +
                theme(text=element_text(size=font_size, family=font_family), 
                      # legend.position = umap_legend_pos,
                      legend.position = "none",
                      plot.margin = margin,
                      # plot.title = element_text(size=font_size, face="plain"),
                      plot.title = element_markdown(hjust = 0.5, vjust = 0, size=font_size, face="italic"),
                      axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                      axis.text.y=element_blank(), axis.ticks.y=element_blank())
}

trans_feature_plot <- function(Trans_seurat_vec,
                               pt_size, cells, 
                               trans_list, trans_label_list,
                               cols, legend_pos, min1, max1, min2, max2,
                               full_idx) {
    pls <- list()
    vjust=0
    for (i in seq_along(trans_list)) {
        if (i == length(trans_list)) {
            legened_pos = legend_pos
            if (length(trans_list) > 3) {
                margin = unit(c(0, 0.8, 0, 0), "cm") # margin: t, r, b, l 
                
            } else {
                margin = unit(c(0, 1.3, 0, 0), "cm") # margin: t, r, b, l
            }
        } else if (i %in% full_idx) {
            legened_pos = "none"
            margin = unit(c(0, 0.3, 0, -0.5), "cm")
        } else {
            legened_pos = "none"
            margin = unit(c(0, 0.3, 0, 0), "cm")
        }
        Trans_seurat = Trans_seurat_vec[[i]]
        trans = trans_list[i]
        label = trans_label_list[i]
        
        if (i %in% full_idx) {
             pl <- FeaturePlot(Trans_seurat, pt.size=pt_size, cells=cells, 
                          features = trans,
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
            pl <- FeaturePlot(Trans_seurat, pt.size=pt_size, cells=cells, 
                              features = trans,
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
# trans_mtx_dir: transcript matrix folder generated by NanoHunter, either a folder or vector of folders with same length as trans_list/trans_label_list
# gene: gene name
# sample_name: sample name
# trans_list: list of transcript id, e.g. c("ENST00000373020.8", "ENST00000373031.8", "ENST00000373034.8")
# trans_label_list: list of transcript label, e.g. c("Trans1", "Trans2", "Trans3")
# nrow: number of rows in the combined UMAP plot (default: 1)
trans_umap_plot <- function(gene_mtx_dir, 
                            trans_mtx_dir, 
                            sample_name, 
                            trans_list, 
                            trans_label_list=c(), 
                            umap_col=umap_dot_red_cols, 
                            text_pos="left", 
                            text_rot=0, 
                            nrow=1) {
    if (length(trans_label_list) == 0) {
        trans_label_list = trans_list
    } else if (length(trans_list) != length(trans_label_list)) {
        stop("trans_list must have same length as trans_label_list.")
    }
    trans_mtx_dirs <- c()
    if (length(trans_list) > 1 && length(trans_mtx_dir) == 1) {
        # trans_mtx_dir <- rep(trans_mtx_dir, length(trans_list))
        for (i in seq_along(trans_list)) {
            trans_mtx_dirs <- c(trans_mtx_dirs, trans_mtx_dir)
        }
    } else {
        trans_mtx_dirs <- c(trans_mtx_dir)
    }
    gene_obj <- make_gene_seurat_obj(gene_mtx_dir)
    trans_obj_vec <- list()

    for (i in seq_along(trans_mtx_dirs)) {
        trans_obj <- make_trans_seurat_obj(trans_mtx_dirs[[i]])
        common_cells = intersect(colnames(gene_obj), colnames(trans_obj))
        trans_obj@reductions$umap@cell.embeddings[common_cells, ] <- gene_obj@reductions$umap@cell.embeddings[common_cells, ]
        trans_obj_vec <- append_list(trans_obj_vec, trans_obj)
    }
    # min_exp=0
    # g <- FetchData(object = gene_obj, vars = gene)
    # common_cells <- intersect(common_cells, rownames(g)[g > min_exp])

    min1 <- min(gene_obj@reductions$umap@cell.embeddings[, 1]) - 0.1
    min2 <- min(gene_obj@reductions$umap@cell.embeddings[, 2]) - 0.1
    max1 <- max(gene_obj@reductions$umap@cell.embeddings[, 1]) + 0.1
    max2 <- max(gene_obj@reductions$umap@cell.embeddings[, 2]) + 0.1

    umap_legend_pos <- c(1, 0.4)
    full_idx <- c(1)
    if (nrow > 1) {
        s = 1
        unit_len = ceiling(length(trans_list)/nrow)
        for (i in 2:nrow) {
            full_idx <- c(full_idx, s + unit_len*(i-1))
        }
    }
    trans_plts <- trans_feature_plot(trans_obj_vec,
                                     umap_pt_size, common_cells, 
                                     trans_list, trans_label_list,
                                     umap_col,
                                     #c("lightgrey", "red"), #,  #eeecec
                                     umap_legend_pos, min1, max1, min2, max2,
                                     full_idx)
    # isoform plot
    pls = list()# (gene_plt)
    for (i in seq_along(trans_plts)) {
        pls[[1+length(pls)]] <- trans_plts[[i]]
    }
    if (length(pls) > 3) {
        if (nrow == 1) {
            gt <- arrangeGrob(grobs = pls, nrow = 1, 
                              widths = c(4.5, rep(4, length(trans_list)-2), 4.5))
        } else {
            gt <- arrangeGrob(grobs = pls, nrow=nrow)
        }
    } else if (length(pls) == 1) {
        gt <- arrangeGrob(grobs = pls, nrow=nrow)
    } else {
        gt <- arrangeGrob(grobs = pls, nrow = 1, widths = c(4.5, rep(4, length(trans_list)-2), 5))
    }
    
    p <- as_ggplot(gt)
    
    if (text_pos == "left") {
        p<- as_ggplot(
            grid.arrange(
                p, 
                left = grid.text(sample_name, rot=text_rot, 
                                 gp=gpar(fontsize=title_font_size, 
                                         lineheight=text_line_height), 
                                 y=0.5, x=-0.5),
                nrow = 1)
        )
    } else if (text_pos == "top") {
        p<- as_ggplot(
            grid.arrange(
                p, 
                top = grid.text(sample_name, rot=text_rot, 
                                 gp=gpar(fontsize=title_font_size, 
                                         lineheight=text_line_height), 
                                 y=1, x=0.5),
                nrow = 1)
        )
    } else if (text_pos == "bottom") {
        p<- as_ggplot(
            grid.arrange(
                p, 
                bottom = grid.text(sample_name, rot=text_rot, 
                                gp=gpar(fontsize=title_font_size, 
                                        lineheight=text_line_height), 
                                y=0, x=0.5),
                nrow = 1)
        )
    } else if (text_pos == "right") {
        p<- as_ggplot(
            grid.arrange(
                p, 
                right = grid.text(sample_name, rot=text_rot, 
                                   gp=gpar(fontsize=title_font_size, 
                                           lineheight=text_line_height), 
                                   y=0.5, x=1),
                nrow = 1)
        )
    }
    return(p)
}

gene_umap_plot <- function(gene_mtx_dir, gene, sample_name, umap_col=umap_dot_red_cols) {
    if (is.object(gene_mtx_dir)) {
        gene_obj <- gene_mtx_dir
    } else {
        if (is.character(gene_mtx_dir) && dir.exists(gene_mtx_dir)) {
            gene_obj <- make_gene_seurat_obj(gene_mtx_dir)
        } else {
            stop("gene_mtx_dir must be a Seurat object or a folder.")
        }
    }
    
    g <- FetchData(object = gene_obj, vars = gene)

    min1 <- min(gene_obj@reductions$umap@cell.embeddings[, 1]) - 0.1
    min2 <- min(gene_obj@reductions$umap@cell.embeddings[, 2]) - 0.1
    max1 <- max(gene_obj@reductions$umap@cell.embeddings[, 1]) + 0.1
    max2 <- max(gene_obj@reductions$umap@cell.embeddings[, 2]) + 0.1
    
    umap_legend_pos <- c(1, 0.4)
    # gene plot
    gene_plt <- gene_feature_plot(gene_obj, umap_pt_size, gene,
                                  umap_col,
                                  umap_legend_pos, min1, max1, min2, max2)
    
    p<- as_ggplot(
            grid.arrange(
                p, 
                left = grid.text(sample_name, rot=0, 
                                 gp=gpar(fontsize=title_font_size, 
                                         lineheight=text_line_height), 
                                 y=0.5, x=-0.5),
                nrow = 1)
    )
    return(p)
}