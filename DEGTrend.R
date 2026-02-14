# DEGTrend - Differentially Expressed Genes Trending Analysis Program 
# Author: Michelle Kojekine 
# Language: R
# version of R: must be 4.4 or higher

# Set CRAN mirror so install.packages() works (e.g. on locked-down or custom setups)
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", dependencies = TRUE)
}
cran_packages <- c(
  "shiny", "ggVennDiagram", "plotly", "DT", "miniUI", "sf", "ggupset",
  "readxl", "ggplot2", "RColorBrewer", "circlize", "viridis", "reshape2",
  "openxlsx", "VennDiagram", "dplyr", "tidyr", "stringr", "matrixStats"
)
bioc_packages <- c(
  "DESeq2", "RUVSeq", "edgeR", "rtracklayer", "ComplexHeatmap",
  "Biobase", "SummarizedExperiment"
)
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, dependencies = TRUE, update = TRUE, ask = FALSE)
  }
}
# Uncomment if install fails for RCurl or GenomeInfoDb:
# install.packages("RCurl", type = "source")
# BiocManager::install("GenomeInfoDb", force = TRUE)

# Load necessary libraries
suppressPackageStartupMessages({
    library(shiny)
    library(ggVennDiagram)
    library(plotly)
    library(DT)
    library(miniUI)
    library(sf)
    library(ggupset)
    library(DESeq2)
    library(RUVSeq)
    library(edgeR)
    library(rtracklayer)
    library(readxl)
    library(ggplot2)
    library(RColorBrewer)
    library(ComplexHeatmap)
    library(circlize)
    library(viridis)
    library(reshape2)
    library(openxlsx)
    library(VennDiagram)
    library(grid)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(matrixStats)   # for rowVars() for pca plot function
    library(Biobase)
    library(SummarizedExperiment)

})
# Suppress repeated ComplexHeatmap messages (use_raster, magick)
tryCatch(
  { ComplexHeatmap::ht_opt$message <- FALSE },
  error = function(e) { }
)
# Only place in the script that may call stop() - used for user-initiated quit
do_quit <- function() {
  cat("User exited the script.\n")
  stop("Script terminated by user.")
}
# Function to handle quit command
handle_quit <- function(input) {
  input <- trimws(as.character(input))
  if (tolower(input) %in% c("quit", "exit")) {
    do_quit()
  }
  if (tolower(input) %in% c("undo", "back")) {
    return("back")
  }
  out <- gsub("^['\"]|['\"]$", "", input)
  return(trimws(out))
}

# Function to construct full paths for output files
get_output_path <- function(filename) {
    file.path(output_dir, filename)
}
# PCA/RLE plots: saved under output_dir/PCA_RLE/<raw|normalised|ruv|final>/
get_pca_rle_path <- function(level, filename) {
    dir_path <- file.path(output_dir, "PCA_RLE", level)
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    file.path(dir_path, filename)
}
# Results (DEG, normalized counts, coldata, mega_results): saved under output_dir/results/
get_results_path <- function(filename) {
    results_dir <- file.path(output_dir, "results")
    dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    file.path(results_dir, filename)
}

# Central color palette for customized colors (gene expression, Venn/UpSet, trending)
DEGTREND_PALETTE <- unique(c(
    "#80b1d3", "#fb8072", "#fdb462", "#b3de69", "#bc80bd", "#c2b975",
    "#ffb3b3", "#bebada", "#8dd3c7", "#27775c", "#31678c", "#574e92",
    "#a80a5b", "#624d0e", "#ffb580", "#99CC66", "#FFCC99", "#FF9999", "#99CCFF",
    "#87CEEB", "#F08080", "#FFCC99", "#98FB98", "#DDA0DD", "#F0E68C"
))
get_custom_colors <- function() {
    DEGTREND_PALETTE
}

# Per-feature output paths (folders created only when the feature is used)
get_genes_of_interest_path <- function(filename) {
    dir_path <- file.path(output_dir, "genes_of_interest")
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    file.path(dir_path, filename)
}
get_heatmaps_path <- function(scaling, type, filename) {
    dir_path <- file.path(output_dir, "heatmaps", scaling, type)
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    file.path(dir_path, filename)
}
get_venn_upset_path <- function(filename, subdir = "") {
    dir_path <- file.path(output_dir, "venn_upset", subdir)
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    file.path(dir_path, filename)
}
get_trending_path <- function(filename) {
    dir_path <- file.path(output_dir, "trending")
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    file.path(dir_path, filename)
}

# Genomic pipeline: outputs under output_dir/genomic/
get_genomic_genes_of_interest_path <- function(filename) {
    dir_path <- file.path(output_dir, "genomic", "genes_of_interest")
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    file.path(dir_path, filename)
}
get_genomic_heatmap_path <- function(filename, subdir = "") {
    dir_path <- if (nzchar(subdir)) file.path(output_dir, "genomic", "heatmaps", subdir) else file.path(output_dir, "genomic", "heatmaps")
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    file.path(dir_path, filename)
}

# For a comparison string "X vs Y", return log2((counts+1)/(ref_avg+1)) with ref = Y
get_logged_for_comparison <- function(comp_str, counts_data, coldata, group_col) {
  parts <- trimws(strsplit(comp_str, " vs ")[[1]])
  if (length(parts) < 2) return(log2(as.matrix(counts_data) + 1))
  ref_condition <- parts[2]
  ref_samps <- rownames(coldata)[as.character(coldata[[group_col]]) == as.character(ref_condition)]
  ref_samples_in_counts <- intersect(ref_samps, colnames(counts_data))
  if (length(ref_samples_in_counts) == 0) return(log2(as.matrix(counts_data) + 1))
  ref_avg <- rowMeans(as.matrix(counts_data)[, ref_samples_in_counts, drop = FALSE], na.rm = TRUE)
  log2((as.matrix(counts_data) + 1) / (ref_avg + 1))
}

# Function to create a heatmap with z-scores
heatmap_vd_z <- function(mat, plot_incl, cols = c("red","black","blue"), col_bias = 1, kmeans_split = NULL, legend_range = 0.1, lim1 = NULL, lim2 = NULL) {
    m <- mat[, plot_incl, drop=FALSE]
    split_vec <- kmeans_split
    rvar <- matrixStats::rowSds(as.matrix(m), na.rm=TRUE)
    keep <- rvar > 0
    m   <- m[keep, , drop=FALSE]
    m_z <- t(scale(t(m)))
    split_vec <- split_vec[keep]
  if (is.null(lim1) || is.null(lim2)) {
    lim1 <- quantile(m_z, probs = legend_range/20, na.rm = TRUE)
    lim2 <- quantile(m_z, probs = 1 - legend_range/20, na.rm = TRUE)
  }
  breaks <- seq(lim1, lim2, length.out = 101)
  mycols <- colorRampPalette(cols, bias = col_bias)(length(breaks))
  plot_obj <- ComplexHeatmap::Heatmap(
    m_z,
    split            = split_vec,
    col              = circlize::colorRamp2(breaks, mycols),
    show_row_dend    = TRUE,
    show_column_dend = FALSE,
    show_column_names= TRUE,
    show_row_names   = FALSE,
    cluster_columns  = FALSE,
    cluster_rows     = TRUE,
    row_title_rot    = 0,
    row_title_side   = "left",
    heatmap_legend_param = list(
      title      = "z-score",
      color_bar  = "continuous"
    ),
    gap = grid::unit(3, "mm")
  )
  
  return(plot_obj)
}


heatmap_vd= function (mat,plot_incl,cols=c("red","black","blue"),col_bias=1,kmeans_split=NULL,legend_range=0.1,lim1=NULL,lim2=NULL){
  if(is.null(lim1)){
    m <- mat[, plot_incl, drop = FALSE]
    n <- length(m)
    if (n > 100000) {
      set.seed(1)
      sidx <- sample(n, min(50000L, n))
      lim1 <- quantile(m[sidx], legend_range/20, na.rm = TRUE)
      lim2 <- quantile(m[sidx], 1 - (legend_range/20), na.rm = TRUE)
    } else {
      lim1 <- quantile(m, legend_range/20, na.rm = TRUE)
      lim2 <- quantile(m, 1 - (legend_range/20), na.rm = TRUE)
    }
  }
    use_raster_val <- nrow(mat) > 2000
    if(any(is.infinite(c(lim1,lim2)))){
        mycols <- colorRampPalette(cols, col_bias)( length(seq(lim1, lim2, length.out = 101)) )
        plot_obj = Heatmap(mat[,plot_incl], split = kmeans_split, col = colorRamp2(seq(lim1, lim2, length.out = 101), mycols),
                    show_row_dend=T, show_column_dend = F, show_column_names=T,show_row_names=F,
                    cluster_columns = F,cluster_rows = T,row_title_rot = 0,row_title_side = "left",
                    heatmap_legend_param=list(title = "logFC", color_bar = "continuous"),gap = unit(3, "mm"),
                    use_raster = use_raster_val)
        cat("Lims: ",lim1," ",lim2)
    } else{
        mycols=colorRampPalette(cols,col_bias)( length(seq(lim1,lim2,(lim2-lim1)/100) ))
        plot_obj= Heatmap(mat[,plot_incl], split = kmeans_split,col = colorRamp2(seq(lim1,lim2,(lim2-lim1)/100),mycols),
                    show_row_dend=T, show_column_dend = F, show_column_names=T,show_row_names=F,
                    cluster_columns = F,cluster_rows = T,row_title_rot = 0,row_title_side = "left",
                    heatmap_legend_param=list(title = "logFC", color_bar = "continuous"),gap = unit(3, "mm"),
                    use_raster = use_raster_val)
        cat("Lims: ",lim1," ",lim2)
    }	
  
  return(plot_obj)

}
heatmap_z <- function(mat, plot_incl, cols = c("red","black","blue"), col_bias = 1,
                        kmeans_split = NULL, legend_range = 0.1, lim1 = NULL, lim2 = NULL,
                        row_split = NULL, palette = NULL, sum_style = FALSE) {
  m <- mat[, plot_incl, drop=FALSE]; m_z <- t(scale(t(m))); m_z[!is.finite(m_z)] <- 0
  if (is.null(lim1)||is.null(lim2)) {
    lim1 <- quantile(m_z, legend_range/20, na.rm=TRUE); lim2 <- quantile(m_z, 1-legend_range/20, na.rm=TRUE)
  }
  breaks <- seq(lim1, lim2, length.out=101)
  mycols <- if (!is.null(palette) && length(palette) >= 101) palette[seq_len(101)] else viridis::inferno(101)
  split_factor <- if (!is.null(row_split)) {
    if (length(row_split)!=nrow(m_z)) return(invisible(NULL))
    factor(row_split)
  } else kmeans_split
  ComplexHeatmap::Heatmap(
    m_z, split=split_factor,
    col=circlize::colorRamp2(breaks, mycols),
    show_row_dend=sum_style, show_column_dend=FALSE,
    show_column_names=TRUE, show_row_names=FALSE,
    cluster_columns=FALSE, cluster_rows=TRUE,
    row_title_rot=0, row_title_side="left",
    heatmap_legend_param=list(
      title="z-score",
      title_gp=grid::gpar(fontsize=14),
      labels_gp=grid::gpar(fontsize=12)
    ),
    gap = if (sum_style) grid::unit(3, "mm") else grid::unit(0, "mm")
  )
}

# Heatmap Function (log2 / logFC)
heatmap_v <- function(mat, plot_incl, cols = c("red", "black", "blue"), col_bias = 1, kmeans_split = NULL, legend_range = 0.1, lim1 = NULL, lim2 = NULL, row_split = NULL, palette = NULL, sum_style = FALSE) {
  if (is.null(lim1) || is.null(lim2)) {
    m <- mat[, plot_incl, drop = FALSE]
    n <- length(m)
    if (n > 100000) {
      set.seed(1)
      sidx <- sample(n, min(50000L, n))
      lim1 <- quantile(m[sidx], probs = legend_range / 20, na.rm = TRUE)
      lim2 <- quantile(m[sidx], probs = 1 - (legend_range / 20), na.rm = TRUE)
    } else {
      lim1 <- quantile(m, probs = legend_range / 20, na.rm = TRUE)
      lim2 <- quantile(m, probs = 1 - (legend_range / 20), na.rm = TRUE)
    }
  }
  if (!is.finite(lim1) || !is.finite(lim2)) return(invisible(NULL))
  breaks <- seq(lim1, lim2, length.out = 101)
  mycols <- if (!is.null(palette) && length(palette) >= 101) palette[seq_len(101)] else viridis::inferno(101)
  if (!is.null(row_split)) {
    if (length(row_split) != nrow(mat)) return(invisible(NULL))
    split <- factor(row_split)
  } else {
    split <- kmeans_split
  }
  use_raster_val <- nrow(mat) > 2000
  Heatmap(mat[, plot_incl],
          split = split,
          col = colorRamp2(breaks, mycols),
          show_row_dend = sum_style,
          show_column_dend = FALSE,
          show_column_names = TRUE,
          show_row_names = FALSE,
          cluster_columns = FALSE,
          cluster_rows = TRUE,
          heatmap_legend_param = list(
            title = "logFC",
            title_gp = gpar(fontsize = 14),
            labels_gp = gpar(fontsize = 12)
          ),
          gap = if (sum_style) grid::unit(3, "mm") else grid::unit(0, "mm"),
          use_raster = use_raster_val)
}

# helper function for venn diagrams
generate_gene_table <- function(selected_genes, category_names, output_file_name = "gene_table.csv") {
    # Ensure `output_file_name` is just the file name, not a full path
    output_file_name <- basename(output_file_name)

    # Initialize an empty data frame to store results
    gene_table <- data.frame(Gene = character(), Segment = character(), stringsAsFactors = FALSE)

    # Iterate over each category to assign segment labels
    for (i in seq_along(selected_genes)) {
        if (length(selected_genes[[i]]) > 0) {
            current_genes <- selected_genes[[i]]
            segment_label <- category_names[i]
            gene_table <- rbind(gene_table, data.frame(Gene = current_genes, Segment = segment_label, stringsAsFactors = FALSE))
        }
    }

    # Adjust segment labels for overlaps
    all_genes <- unique(gene_table$Gene)
    gene_segments <- sapply(all_genes, function(gene) {
        belonging_sets <- names(selected_genes)[sapply(selected_genes, function(set) gene %in% set)]
        if (length(belonging_sets) > 1) {
            paste(belonging_sets, collapse = " intersect ")
        } else {
            belonging_sets
        }
    })
    gene_table <- data.frame(Gene = all_genes, Segment = gene_segments, stringsAsFactors = FALSE)

    # Define the output file path using `get_output_path`
    output_file <- get_output_path(output_file_name)

    # Debug: Print output file path
    cat("Writing to file:", output_file, "\n")

    # Try writing to CSV and handle errors
    tryCatch({
        write.csv(gene_table, file = output_file, row.names = FALSE)
        cat("Gene table successfully exported to:", output_file, "\n")
    }, error = function(e) {
        cat("Error writing file:", e$message, "\n")
    })
}
# Function to determine significance stars
get_significance_stars <- function(padj_values) {
    sapply(padj_values, function(p) {
        if (is.na(p)) return("")  # No significance if NA
        if (p < 0.001) return("***")  # Highly significant
        if (p < 0.01) return("**")  # Very significant
        if (p < 0.05) return("*")  # Significant
        return("")  # Not significant
    })
}

# Bar width for gene expression / trending plots (by number of levels)
gene_expr_bar_width <- function(n_levels) {
    if (n_levels == 3L) return(0.65)
    if (n_levels %in% c(4L, 5L)) return(0.55)
    if (n_levels == 6L) return(0.45)
    max(0.25, min(0.45, 0.9 / max(1L, n_levels)))
}
# Stroke width for bar outline
gene_expr_stroke <- function(n_levels) {
    if (n_levels > 6L) 0.45 else 0.9
}
# Export dimensions for gene expression plots
gene_expr_width_in  <- 2.853
gene_expr_height_in <- 3.963
gene_expr_dpi       <- 300
# Bar/circle colors from central palette
get_gene_expression_bar_colors <- function(k) {
    pal <- get_custom_colors()
    pal[((seq_len(k) - 1L) %% length(pal)) + 1L]
}

# Function to plot gene expression with significance brackets (ggplot2, SVG + PNG)
plot_gene_expression <- function(gene, expr_data, coldata, group_col, matches) {
    if (missing(gene) || is.null(gene) || length(gene) == 0) {
        cat("Error: Gene name must be provided\n")
        return(invisible(NULL))
    }
    if (missing(expr_data) || is.null(expr_data) || length(expr_data) == 0) {
        cat("Error: Expression data must be provided\n")
        return(invisible(NULL))
    }
    if (missing(coldata) || is.null(coldata)) {
        cat("Error: Coldata must be provided\n")
        return(invisible(NULL))
    }
    if (!group_col %in% colnames(coldata)) {
        cat(sprintf("Error: Column '%s' not found in coldata\n", group_col))
        return(invisible(NULL))
    }
    if (length(expr_data) != nrow(coldata)) {
        cat(sprintf("Error: expr_data length (%d) must equal coldata rows (%d).\n", length(expr_data), nrow(coldata)))
        return(invisible(NULL))
    }

    group_vec <- factor(trimws(as.character(coldata[[group_col]])))
    avg_expr_per_condition <- tapply(expr_data, group_vec, mean, na.rm = TRUE)
    if (is.null(avg_expr_per_condition) || length(avg_expr_per_condition) == 0L) {
        cat(sprintf("⚠️ No conditions in group column '%s'. Skipping plot...\n", group_col))
        return(invisible(NULL))
    }

    valid_values <- avg_expr_per_condition[is.finite(avg_expr_per_condition)]
    if (length(valid_values) == 0) {
        cat(sprintf("⚠️ No finite values available for gene %s. Skipping plot...\n", gene))
        return(invisible(NULL))
    }

    k <- length(avg_expr_per_condition)
    bar_width <- gene_expr_bar_width(k)
    max_y <- max(valid_values, na.rm = TRUE)
    y_limit <- max_y * 1.15
    tip_length <- 0.02 * y_limit
    step_between <- 0.04 * y_limit

    bar_values <- as.numeric(avg_expr_per_condition)
    bar_values[!is.finite(bar_values) | bar_values < 0] <- 0
    conditions <- names(avg_expr_per_condition)
    df_bars <- data.frame(
        condition = factor(conditions, levels = conditions),
        value     = bar_values,
        stringsAsFactors = FALSE
    )

    significance_stars <- get_significance_stars(matches$padj)
    significant_matches <- matches[significance_stars != "", , drop = FALSE]
    significant_stars <- significance_stars[significance_stars != ""]

    segs_h <- list()
    segs_v <- list()
    star_df <- list()
    y_positions <- numeric(0)

    for (i in seq_len(nrow(significant_matches))) {
        comp <- significant_matches$comparison[i]
        comp_conditions <- unlist(strsplit(comp, " vs "))
        if (length(comp_conditions) != 2L) next
        idx1 <- match(comp_conditions[1L], conditions)
        idx2 <- match(comp_conditions[2L], conditions)
        if (is.na(idx1) || is.na(idx2) || idx1 == idx2) next
        val1 <- bar_values[idx1]
        val2 <- bar_values[idx2]
        y_bar_max <- max(val1, val2)
        y_line <- if (length(y_positions) == 0) {
            y_bar_max + tip_length
        } else {
            max(max(y_positions) + step_between, y_bar_max + tip_length)
        }
        y_positions <- c(y_positions, y_line)
        segs_h[[i]] <- data.frame(x = idx1, xend = idx2, y = y_line, yend = y_line)
        segs_v[[i]] <- data.frame(
            x = c(idx1, idx2), xend = c(idx1, idx2),
            y = c(val1, val2), yend = c(y_line, y_line)
        )
        star_df[[i]] <- data.frame(
            x = (idx1 + idx2) / 2,
            y = y_line + 0.02 * y_limit,
            label = significant_stars[i]
        )
    }

    y_max_plot <- if (length(y_positions) > 0) max(y_positions) + 0.08 * y_limit else y_limit

    condition_levels <- levels(df_bars$condition)
    fill_values <- setNames(get_gene_expression_bar_colors(k), condition_levels)
    bar_stroke <- gene_expr_stroke(k)

    p <- ggplot2::ggplot(df_bars, ggplot2::aes(x = condition, y = value, fill = condition)) +
        ggplot2::geom_col(width = bar_width, colour = "black", linewidth = bar_stroke) +
        ggplot2::scale_fill_manual(values = fill_values, limits = condition_levels, guide = "none") +
        ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.12, 0.08))) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::coord_cartesian(ylim = c(0, y_max_plot), clip = "off") +
        ggplot2::labs(title = gene, y = "Normalized counts", x = NULL) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            plot.title       = ggplot2::element_text(size = 24, face = "bold", hjust = 0.5),
            axis.title.y     = ggplot2::element_text(size = 18, face = "bold"),
            axis.text.x      = ggplot2::element_text(size = 14, face = "bold", angle = 45, hjust = 1, vjust = 1),
            axis.text.y      = ggplot2::element_text(size = 14, face = "bold"),
            axis.line        = ggplot2::element_line(linewidth = 1.1),
            axis.ticks       = ggplot2::element_line(linewidth = 1.1),
            axis.ticks.length = grid::unit(8, "pt")
        )

    for (i in seq_along(segs_h)) {
        p <- p + ggplot2::geom_segment(
            data = segs_h[[i]],
            ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
            inherit.aes = FALSE, linewidth = 1.3
        )
    }
    for (i in seq_along(segs_v)) {
        p <- p + ggplot2::geom_segment(
            data = segs_v[[i]],
            ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
            inherit.aes = FALSE, linewidth = 1.3
        )
    }
    for (i in seq_along(star_df)) {
        p <- p + ggplot2::geom_text(
            data = star_df[[i]],
            ggplot2::aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 5, fontface = "bold"
        )
    }

    plot_file_svg <- get_genes_of_interest_path(paste0(gene, "_expression_plot.svg"))
    plot_file_png <- get_genes_of_interest_path(paste0(gene, "_expression_plot.png"))
    ggplot2::ggsave(
        plot_file_svg, plot = p,
        width = gene_expr_width_in, height = gene_expr_height_in,
        dpi = gene_expr_dpi, device = "svg", units = "in"
    )
    ggplot2::ggsave(
        plot_file_png, plot = p,
        width = gene_expr_width_in, height = gene_expr_height_in,
        dpi = gene_expr_dpi, device = "png", units = "in"
    )
    cat(sprintf("✅ Plot for %s saved to: %s, %s\n", gene, plot_file_svg, plot_file_png))
}

# Genomic pipeline: plot log2FC by comparison (line + points, significance stars)
plot_gene_expression_from_genomic <- function(gene, matches) {
    if (missing(gene) || is.null(gene) || length(gene) == 0) {
        cat("Error: Gene name must be provided\n")
        return(invisible(NULL))
    }
    if (missing(matches) || is.null(matches) || nrow(matches) == 0) {
        cat("Error: matches (genomic rows) must be provided\n")
        return(invisible(NULL))
    }
    if (!"log2FoldChange" %in% colnames(matches)) {
        cat("Error: Column 'log2FoldChange' not found in genomic data. Cannot plot.\n")
        return(invisible(NULL))
    }
    if (!"comparison" %in% colnames(matches)) {
        cat("Error: Column 'comparison' not found in genomic data. Cannot plot.\n")
        return(invisible(NULL))
    }

    lfc_vals <- as.numeric(matches$log2FoldChange)
    comp_names <- as.character(matches$comparison)
    valid <- is.finite(lfc_vals)
    if (!any(valid)) {
        cat(sprintf("⚠️ No finite log2FoldChange for gene %s. Skipping genomic plot...\n", gene))
        return(invisible(NULL))
    }
    lfc_vals <- lfc_vals[valid]
    comp_names <- comp_names[valid]
    per_bar_stars <- if ("padj" %in% colnames(matches)) {
        get_significance_stars(as.numeric(matches$padj)[valid])
    } else {
        rep("", length(lfc_vals))
    }

    k <- length(lfc_vals)
    max_y <- max(abs(lfc_vals), na.rm = TRUE)
    df_bars <- data.frame(
        condition = factor(comp_names, levels = comp_names),
        value     = lfc_vals,
        star      = per_bar_stars
    )
    y_max_plot <- max_y * 1.3

    p <- ggplot2::ggplot(df_bars, ggplot2::aes(x = condition, y = value)) +
        ggplot2::geom_line(ggplot2::aes(group = 1), linewidth = 2, colour = "black") +
        ggplot2::geom_point(size = 4, colour = "black") +
        ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.22, 0.08))) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::coord_cartesian(ylim = c(-y_max_plot, y_max_plot), clip = "off") +
        ggplot2::labs(title = gene, y = "log2 Fold Change", x = NULL) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            plot.title       = ggplot2::element_text(size = 24, face = "bold", hjust = 0.5),
            axis.title.y     = ggplot2::element_text(size = 18, face = "bold"),
            axis.text.x      = ggplot2::element_text(size = 14, face = "bold", angle = 45, hjust = 1, vjust = 1),
            axis.text.y      = ggplot2::element_text(size = 14, face = "bold"),
            axis.line        = ggplot2::element_line(linewidth = 1.1),
            axis.ticks       = ggplot2::element_line(linewidth = 1.1),
            axis.ticks.length = grid::unit(8, "pt"),
            plot.margin      = ggplot2::margin(l = 12, unit = "pt")
        )

    for (i in seq_len(k)) {
        if (nzchar(per_bar_stars[i])) {
            p <- p + ggplot2::geom_text(
                data = data.frame(x = i, y = lfc_vals[i] + 0.12 * max_y, label = per_bar_stars[i]),
                ggplot2::aes(x = x, y = y, label = label),
                inherit.aes = FALSE, size = 5, fontface = "bold", vjust = -0.2
            )
        }
    }

    plot_file_svg <- get_genomic_genes_of_interest_path(paste0(gene, "_expression_plot.svg"))
    plot_file_png <- get_genomic_genes_of_interest_path(paste0(gene, "_expression_plot.png"))
    ggplot2::ggsave(
        plot_file_svg, plot = p,
        width = gene_expr_width_in, height = gene_expr_height_in,
        dpi = gene_expr_dpi, device = "svg", units = "in"
    )
    ggplot2::ggsave(
        plot_file_png, plot = p,
        width = gene_expr_width_in, height = gene_expr_height_in,
        dpi = gene_expr_dpi, device = "png", units = "in"
    )
    cat(sprintf("✅ Genomic expression plot for %s saved to: %s, %s\n", gene, plot_file_svg, plot_file_png))
}

# Function to plot with significance brackets
plot_deseq2_patterns2 <- function(de_results_list, logged_df, sample_conditions, main_titles, degs_data) {
    op <- par(no.readonly = TRUE)  # Save current graphical parameters
    par(mfrow = c(1, length(de_results_list)), mar = c(6, 6, 8, 4), mgp = c(3, 1, 0), lend  = 1)  
    
    excel_grid  <- function(ylim) {
    y_grid <- pretty(c(0, ylim), n = 6)
    abline(h = y_grid, col = "grey80", lwd = .7)
    }

    axis_excel  <- function(ylim) {
        axis(2, las = 1, col = "grey40", col.axis = "grey20", tck = -0.02)
        axis(1, labels = FALSE, tick = FALSE)       # we’ll draw labels via barplot()
        box(col = "grey40")
    }

    custom_colors <- c("#80b1d3", "#fb8072", "#fdb462", "#b3de69", "#bc80bd", "#c2b975", 
                       "#ffb3b3", "#bebada", "#8dd3c7", "#27775c", "#31678c", "#574e92", 
                       "#a80a5b", "#624d0e", "#ffb580")

    for (i in seq_along(de_results_list)) {
        # expr_data <- logged_df[de_results_list[[i]], , drop = FALSE]

        ## --- keep only genes that exist in the matrix ----
        genes_present <- intersect(de_results_list[[i]], rownames(logged_df))
        if (length(genes_present) == 0) {
            message(sprintf("⚠️  No matching genes for '%s' – skipping panel.", main_titles[i]))
            next
        }

        expr_data <- logged_df[genes_present, , drop = FALSE]
            # Compute average expression per condition
            avg_expr_per_condition <- sapply(sample_conditions, function(samples) {
                rowMeans(expr_data[, samples, drop = FALSE], na.rm = TRUE)
            })
        
        overall_mean_expr <- colMeans(avg_expr_per_condition, na.rm = TRUE)
        # plot_x_positions <- 1:length(overall_mean_expr)

        # # Assign colors based on the number of bars
        # bar_colors <- custom_colors[seq_len(length(overall_mean_expr))]
        # Scale the y-axis properly
        max_y <- max(overall_mean_expr, na.rm = TRUE)
        y_limit <- max_y * 1.3  # Extra space for significance brackets

        ## pick one width / spacing and reuse it everywhere
        nbars   <- length(overall_mean_expr)      # how many conditions?
        bar_w   <- min(0.00000000001, 0.00000000001 / nbars)            # width shrinks as nbars grows
        bar_sp  <- 5

        ## ----– 1. work out bar centres BEFORE opening the window –------------
        bar_centres <- barplot(overall_mean_expr,
                            width = bar_w,
                            space = bar_sp,
                            plot  = FALSE)                 # <─ calculates x coords

        xlim <- c(bar_centres[1] - bar_w/2, tail(bar_centres,1) + bar_w/2)

        ## -------- blank canvas, add grid & axes ----------
        plot.new()
        par(xaxs = "r")
        plot.window(xlim = xlim, ylim = c(0, y_limit))
        excel_grid(y_limit)
        axis_excel(y_limit)

        ## -------- bars  -----------------------------------
        bar_x <- barplot(overall_mean_expr,
                        width      = bar_w,
                        space      = bar_sp,             # <─ keeps clear gaps between bars
                        col        = custom_colors[seq_along(overall_mean_expr)],
                        names.arg  = names(sample_conditions),
                        las        = 1,
                        ylim       = c(0, y_limit),
                        border     = "black",
                        add        = TRUE)          # draw inside existing canvas

        title(main = main_titles[i], line = 2.5)

        ## ADDED: Parse the user’s chosen comparisons from `main_titles[i]`
        selected_comps <- unlist(strsplit(main_titles[i], " INT "))
        # For each selected comparison, see if it matches `degs_data$comparison`
        # or if removing _up/_down helps find a match.
        actual_comps <- character(0)

        for (sc in selected_comps) {
            sc <- sub(" \\((up|down)\\)$", "_\\1", sc)
            if (sc %in% degs_data$comparison) {
                # Perfect match
                actual_comps <- c(actual_comps, sc)
            } else {
                # Try removing trailing _up or _down
                sc_no_dir <- sub("_down$", "", sub("_up$", "", sc))
                if (sc_no_dir %in% degs_data$comparison) {
                    actual_comps <- c(actual_comps, sc_no_dir)
                }
            }
        }
        # Now `actual_comps` has only the comparisons that truly exist in degs_data$comparison

        # We'll draw brackets only for these actual_comps
        # y_offset <- max_y * 1.1
        y_offset <- max_y * 1.05
        y_positions <- c()

        # for (comp in actual_comps) {
        #     comp_conditions <- unlist(strsplit(comp, " vs "))
        #     if (length(comp_conditions) == 2) {
        #             c1 <- which(names(sample_conditions) == comp_conditions[1])
        #             c2 <- which(names(sample_conditions) == comp_conditions[2])
        #             if (length(c1) == 1 && length(c2) == 1 && c1 != c2) {
        #             bracket_y  <- ifelse(length(y_positions) == 0,
        #                                 y_offset,
        #                                 max(y_positions) + 0.08 * max_y)
        #             y_positions <- c(y_positions, bracket_y)
        #             segments(bar_x[c1], bracket_y, bar_x[c2], bracket_y, lwd = 2)
        #             segments(bar_x[c1], bracket_y, bar_x[c1], bracket_y - 0.03 * max_y, lwd = 2)
        #             segments(bar_x[c2], bracket_y, bar_x[c2], bracket_y - 0.03 * max_y, lwd = 2)

        #             star <- get_significance_stars(degs_data$padj[degs_data$comparison == comp])[1]
        #             text(mean(c(bar_x[c1], bar_x[c2])), bracket_y + 0.03 * max_y,
        #                 labels = star, cex = 1.5, font = 2)
        #         }
        #     }
        # }

        norm_name <- function(x) {
            x <- trimws(tolower(x))          # lower-case + trim
            x <- sub("[ _]+$", "", x)        # drop trailing spaces/underscores
            x <- sub("[ _]+",  "_", x)       # collapse runs to single “_”
            x
        }

        names_norm <- norm_name(names(sample_conditions))

        for (comp in actual_comps) {
            parts <- norm_name(unlist(strsplit(comp, " vs ")))

            c1 <- match(parts[1], names_norm)
            c2 <- match(parts[2], names_norm)

            if (!is.na(c1) && !is.na(c2) && c1 != c2) {
                bracket_y  <- ifelse(length(y_positions) == 0,
                                    y_offset,
                                    max(y_positions) + 0.08 * max_y)
                y_positions <- c(y_positions, bracket_y)

                segments(bar_x[c1], bracket_y, bar_x[c2], bracket_y, lwd = 2)
                segments(bar_x[c1], bracket_y, bar_x[c1], bracket_y - 0.03 * max_y, lwd = 2)
                segments(bar_x[c2], bracket_y, bar_x[c2], bracket_y - 0.03 * max_y, lwd = 2)

                star <- get_significance_stars(
                            degs_data$padj[degs_data$comparison == comp])[1]
                text(mean(c(bar_x[c1], bar_x[c2])), bracket_y + 0.03 * max_y,
                    labels = star, cex = 1.5, font = 2)
            }
        }

    }
    par(op)
}

# Helpers for trending pattern logic
approx_eq <- function(a, b, x) {
  r <- diff(range(x, na.rm = TRUE))
  tol <- if (r > 0) 0.05 * r else 0
  abs(a - b) <= tol
}
strong_inc <- function(a, b, x) {
  r <- diff(range(x, na.rm = TRUE))
  if (r <= 0) return(b > a)
  (b - a) >= 0.1 * r
}
strong_dec <- function(a, b, x) {
  r <- diff(range(x, na.rm = TRUE))
  if (r <= 0) return(b < a)
  (a - b) >= 0.1 * r
}

get_trending_patterns <- function(n_conditions) {
  if (n_conditions == 3L) {
    list(
      "Gradual Uptrend" = function(x) x[1] < x[2] && x[2] < x[3],
      "Gradual Downtrend" = function(x) x[1] > x[2] && x[2] > x[3],
      "Linear Uptrend with Plateau" = function(x) x[1] < x[2] && approx_eq(x[2], x[3], x),
      "Linear Downtrend with Plateau" = function(x) x[1] > x[2] && approx_eq(x[2], x[3], x),
      "Gentle Uptrend with Acceleration" = function(x) (x[2] - x[1]) > 0 && (x[3] - x[2]) > (x[2] - x[1]),
      "Gentle Downtrend with Acceleration" = function(x) (x[1] - x[2]) > 0 && (x[2] - x[3]) > (x[1] - x[2]),
      "Immediate Rise" = function(x) approx_eq(x[1], x[2], x) && x[2] < x[3],
      "Immediate Drop" = function(x) approx_eq(x[1], x[2], x) && x[2] > x[3],
      "Spike - Peak" = function(x) x[1] < x[2] && x[2] > x[3],
      "Spike - Trough" = function(x) x[1] > x[2] && x[2] < x[3],
      "Sharp Linear Increase" = function(x) strong_inc(x[1], x[2], x) && strong_inc(x[2], x[3], x),
      "Sharp Linear Decrease" = function(x) strong_dec(x[1], x[2], x) && strong_dec(x[2], x[3], x)
    )
  } else if (n_conditions == 4L) {
    list(
      "Gradual Uptrend" = function(x) x[1] < x[2] && x[2] < x[3] && x[3] < x[4],
      "Gradual Downtrend" = function(x) x[1] > x[2] && x[2] > x[3] && x[3] > x[4],
      "Plateau Uptrend" = function(x) x[1] < x[2] && approx_eq(x[2], x[3], x) && approx_eq(x[3], x[4], x),
      "Plateau Downtrend" = function(x) x[1] > x[2] && approx_eq(x[2], x[3], x) && approx_eq(x[3], x[4], x),
      "Gentle Uptrend" = function(x) x[1] < x[2] && x[2] < x[3] && x[3] <= x[4] && (x[3] - x[2]) >= (x[2] - x[1]),
      "Gentle Downtrend" = function(x) x[1] > x[2] && x[2] > x[3] && x[3] >= x[4] && (x[2] - x[3]) >= (x[1] - x[2]),
      "Immediate Rise" = function(x) approx_eq(x[1], x[2], x) && approx_eq(x[2], x[3], x) && x[3] < x[4],
      "Immediate Drop" = function(x) approx_eq(x[1], x[2], x) && approx_eq(x[2], x[3], x) && x[3] > x[4],
      "Spike - Peak" = function(x) {
        imax <- which.max(x)
        imax %in% c(2L, 3L) && (imax == 2L && x[2] >= x[4] && x[1] < x[2] && x[2] > x[3] ||
                                imax == 3L && x[3] >= x[1] && x[2] < x[3] && x[3] > x[4])
      },
      "Spike - Trough" = function(x) {
        imin <- which.min(x)
        imin %in% c(2L, 3L) && (imin == 2L && x[2] <= x[4] && x[1] > x[2] && x[2] < x[3] ||
                                imin == 3L && x[3] <= x[1] && x[2] > x[3] && x[3] < x[4])
      },
      "Sharp Linear Increase" = function(x) strong_inc(x[1], x[2], x) && strong_inc(x[2], x[3], x) && strong_inc(x[3], x[4], x),
      "Sharp Linear Decrease" = function(x) strong_dec(x[1], x[2], x) && strong_dec(x[2], x[3], x) && strong_dec(x[3], x[4], x),
      "Sudden Jump at Start" = function(x) x[1] < x[2] && approx_eq(x[2], x[3], x) && approx_eq(x[3], x[4], x),
      "Sudden Drop at Start" = function(x) x[1] > x[2] && approx_eq(x[2], x[3], x) && approx_eq(x[3], x[4], x),
      "Oscillatory (Rise-Fall-Rise)" = function(x) x[1] < x[2] && x[2] > x[3] && x[3] < x[4],
      "Oscillatory (Fall-Rise-Fall)" = function(x) x[1] > x[2] && x[2] < x[3] && x[3] > x[4]
    )
  } else {
    list()
  }
}

plot_trending_ggplot <- function(condition_means, condition_names, title,
                                 degs_data = NULL, comparisons = NULL, show_title = TRUE) {
  if (length(condition_means) != length(condition_names) || length(condition_means) == 0L) {
    cat("Error: condition_means and condition_names must have same positive length\n")
    return(invisible(NULL))
  }
  k <- length(condition_means)
  bar_width <- gene_expr_bar_width(k)
  max_y <- max(condition_means, na.rm = TRUE)
  y_limit <- max_y * 1.15
  tip_length <- 0.02 * y_limit
  step_between <- 0.04 * y_limit

  bar_values <- as.numeric(condition_means)
  bar_values[!is.finite(bar_values) | bar_values < 0] <- 0
  conditions <- condition_names
  df_bars <- data.frame(
    condition = factor(conditions, levels = conditions),
    value     = bar_values,
    stringsAsFactors = FALSE
  )

  y_max_plot <- y_limit
  segs_h <- list()
  segs_v <- list()
  star_df <- list()
  y_positions <- numeric(0)

  if (!is.null(degs_data) && !is.null(comparisons) && length(comparisons) > 0L &&
      "comparison" %in% colnames(degs_data) && "padj" %in% colnames(degs_data)) {
    norm_name <- function(x) {
      x <- trimws(tolower(x))
      x <- sub("[ _]+$", "", x)
      x <- sub("[ _]+", "_", x)
      x
    }
    names_norm <- norm_name(conditions)
    for (comp in comparisons) {
      comp_clean <- sub("_down$", "", sub("_up$", "", comp))
      if (!comp_clean %in% degs_data$comparison) next
      parts <- norm_name(unlist(strsplit(comp_clean, " vs ")))
      if (length(parts) != 2L) next
      idx1 <- match(parts[1L], names_norm)
      idx2 <- match(parts[2L], names_norm)
      if (is.na(idx1) || is.na(idx2) || idx1 == idx2) next
      padj_vals <- degs_data$padj[degs_data$comparison == comp_clean]
      if (length(padj_vals) == 0L) next
      star <- get_significance_stars(as.numeric(padj_vals[1L]))[1L]
      if (is.na(star) || star == "") next
      val1 <- bar_values[idx1]
      val2 <- bar_values[idx2]
      y_line <- if (length(y_positions) == 0) {
        max(val1, val2) + tip_length
      } else {
        max(max(y_positions) + step_between, max(val1, val2) + tip_length)
      }
      y_positions <- c(y_positions, y_line)
      segs_h <- c(segs_h, list(data.frame(x = idx1, xend = idx2, y = y_line, yend = y_line)))
      segs_v <- c(segs_v, list(data.frame(
        x = c(idx1, idx2), xend = c(idx1, idx2),
        y = c(val1, val2), yend = c(y_line, y_line)
      )))
      star_df <- c(star_df, list(data.frame(
        x = (idx1 + idx2) / 2,
        y = y_line + 0.02 * y_limit,
        label = star
      )))
    }
    if (length(y_positions) > 0) y_max_plot <- max(y_positions) + 0.08 * y_limit
  }

  condition_levels <- levels(df_bars$condition)
  fill_values <- setNames(get_gene_expression_bar_colors(k), condition_levels)
  bar_stroke <- gene_expr_stroke(k)
  p <- ggplot2::ggplot(df_bars, ggplot2::aes(x = condition, y = value, fill = condition)) +
    ggplot2::geom_col(width = bar_width, colour = "black", linewidth = bar_stroke) +
    ggplot2::scale_fill_manual(values = fill_values, limits = condition_levels, guide = "none") +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.12, 0.08))) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = c(0, y_max_plot), clip = "off") +
    ggplot2::labs(title = if (show_title) title else NULL, y = "Normalized counts", x = NULL) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = 24, face = "bold", hjust = 0.5),
      axis.title.y     = ggplot2::element_text(size = 18, face = "bold"),
      axis.text.x      = ggplot2::element_text(size = 14, face = "bold", angle = 45, hjust = 1, vjust = 1),
      axis.text.y      = ggplot2::element_text(size = 14, face = "bold"),
      axis.line        = ggplot2::element_line(linewidth = 1.1),
      axis.ticks       = ggplot2::element_line(linewidth = 1.1),
      axis.ticks.length = grid::unit(8, "pt")
    )

  for (i in seq_along(segs_h)) {
    p <- p + ggplot2::geom_segment(
      data = segs_h[[i]],
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE, linewidth = 1.3
    )
  }
  for (i in seq_along(segs_v)) {
    p <- p + ggplot2::geom_segment(
      data = segs_v[[i]],
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE, linewidth = 1.3
    )
  }
  for (i in seq_along(star_df)) {
    p <- p + ggplot2::geom_text(
      data = star_df[[i]],
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE, size = 5, fontface = "bold"
    )
  }
  p
}

# Function to calculate and plot default trending patterns.
# Uses counts_data only (not logged_df). X-axis shows condition names (not comparisons).
generate_default_trending_plots <- function(counts_data, coldata, significant_genes, group_col, padj_threshold) {
    significant_genes <- significant_genes[significant_genes$padj < padj_threshold, , drop = FALSE]
    significant_genes <- significant_genes[complete.cases(significant_genes), ]
    if (nrow(significant_genes) == 0) {
        cat(sprintf("No significant genes found for padj < %.2f. Skipping default trending plots.\n", padj_threshold))
        return()
    }

    comparisons_all <- unique(as.character(significant_genes$comparison))
    genes_sig_all <- if (length(comparisons_all) > 0L) {
        sig_by_comp <- split(significant_genes$Gene, significant_genes$comparison)
        Reduce(intersect, sig_by_comp)
    } else character(0)
    if (length(genes_sig_all) == 0) {
        cat("No genes significant in all comparisons. Skipping default trending plots.\n")
        return()
    }

    counts_data_filtered <- counts_data[rownames(counts_data) %in% genes_sig_all, , drop = FALSE]
    sample_conditions <- split(rownames(coldata), coldata[[group_col]])
    trending_conditions <- unique(coldata[[group_col]])
    n_cond <- length(trending_conditions)

    tryCatch({
        avg_expr_per_condition <- t(apply(counts_data_filtered, 1, function(row) {
            tapply(row, coldata[[group_col]], mean, na.rm = TRUE)
        }))
    }, error = function(e) {
        cat("Error calculating average expressions:", e$message, "\n")
        return()
    })

    patterns <- get_trending_patterns(n_cond)
    if (length(patterns) == 0) {
        cat("Default trending plots are only available for 3 or 4 conditions.\n")
        return()
    }

    condition_names <- colnames(avg_expr_per_condition)
    for (pattern_name in names(patterns)) {
        tryCatch({
            matching_genes <- rownames(avg_expr_per_condition)[apply(avg_expr_per_condition, 1, patterns[[pattern_name]])]
            if (length(matching_genes) == 0) next

            condition_means <- colMeans(avg_expr_per_condition[matching_genes, , drop = FALSE], na.rm = TRUE)
            p <- plot_trending_ggplot(
                as.numeric(condition_means),
                condition_names,
                pattern_name,
                degs_data = significant_genes,
                comparisons = comparisons_all
            )
            if (is.null(p)) next

            plot_name_png <- get_trending_path(paste0("trending_", gsub(" ", "_", pattern_name), ".png"))
            plot_name_svg <- get_trending_path(paste0("trending_", gsub(" ", "_", pattern_name), ".svg"))
            output_file <- get_trending_path(paste0("trending_data_", gsub(" ", "_", pattern_name), ".csv"))

            ggplot2::ggsave(plot_name_png, plot = p, width = gene_expr_width_in, height = gene_expr_height_in, dpi = gene_expr_dpi, device = "png", units = "in")
            ggplot2::ggsave(plot_name_svg, plot = p, width = gene_expr_width_in, height = gene_expr_height_in, dpi = gene_expr_dpi, device = "svg", units = "in")

            export_data <- data.frame(
                Gene = rep(matching_genes, each = ncol(avg_expr_per_condition)),
                Condition = rep(condition_names, length(matching_genes)),
                Expression = as.vector(avg_expr_per_condition[matching_genes, , drop = FALSE])
            )
            write.csv(export_data, output_file, row.names = FALSE)
            cat(sprintf("Trending pattern '%s': %d genes. Plot: %s, %s. Data: %s\n", pattern_name, length(matching_genes), plot_name_png, plot_name_svg, output_file))
        }, error = function(e) {
            cat(sprintf("Error for pattern %s: %s\n", pattern_name, e$message))
        })
    }
}
# --- Venn/UpSet helpers (shared by export_venn_upset and export_venn_upset_split) ---
build_upset_data <- function(sets_named_list) {
  Nsets <- length(sets_named_list)
  set_names <- names(sets_named_list)
  all_genes <- sort(unique(unlist(sets_named_list)))
  incidence_list <- lapply(all_genes, function(gene) {
    member_sets <- set_names[vapply(sets_named_list, function(s) gene %in% s, logical(1))]
    member_sorted <- sort(member_sets)
    key <- paste(member_sorted, collapse = "&")
    data.frame(Gene = gene, IntersectKey = key, stringsAsFactors = FALSE)
  })
  incidence_df <- do.call(rbind, incidence_list)
  region_ids_present <- sort(unique(incidence_df$IntersectKey))
  region_to_genes <- setNames(
    lapply(region_ids_present, function(rid) sort(unique(incidence_df$Gene[incidence_df$IntersectKey == rid]))),
    region_ids_present
  )
  counts_table <- as.integer(table(incidence_df$IntersectKey)[region_ids_present])
  df_counts <- data.frame(region = region_ids_present, count = counts_table, stringsAsFactors = FALSE)
  human_labels <- setNames(
    sapply(region_ids_present, function(rid) {
      parts <- strsplit(rid, "&")[[1]]
      if (length(parts) == 1) paste0(parts, " only") else paste(parts, collapse = "&")
    }),
    region_ids_present
  )
  df_counts$label <- human_labels[df_counts$region]
  list(counts_df = df_counts, region_to_genes = region_to_genes, region_ids_present = region_ids_present, human_labels = human_labels)
}

get_circle_fills <- function(n_sets) get_gene_expression_bar_colors(n_sets)

blend_hex_colors <- function(hex_vec) {
  if (length(hex_vec) == 0) return("#808080")
  if (length(hex_vec) == 1) return(hex_vec[1])
  rgb_mat <- grDevices::col2rgb(hex_vec)
  rgb_avg <- rowMeans(rgb_mat)
  grDevices::rgb(rgb_avg[1], rgb_avg[2], rgb_avg[3], maxColorValue = 255)
}

get_region_fill_color <- function(rid, set_names, circle_fills) {
  parts <- strsplit(rid, "&", fixed = TRUE)[[1]]
  hex_vec <- circle_fills[parts]
  hex_vec <- hex_vec[!is.na(hex_vec)]
  blend_hex_colors(as.character(hex_vec))
}

get_clicked_region <- function(click, df_counts, region_ids_present) {
  if (is.null(click)) return(NULL)
  rid <- if (is.data.frame(click) && "customdata" %in% names(click)) click$customdata[[1]] else click[["customdata"]]
  rid <- as.character(rid)
  if (length(rid) > 1) rid <- rid[1]
  if (!is.na(rid) && rid %in% region_ids_present) return(rid)
  pn <- click[["pointNumber"]]
  if (!is.null(pn) && !is.na(pn)) {
    idx <- as.integer(pn) + 1
    if (idx >= 1 && idx <= nrow(df_counts)) return(df_counts$region[idx])
  }
  NULL
}

create_venn_diagram <- function(sets_named_list, set_names, circle_fills, highlight_region = NULL) {
  n_sets <- length(sets_named_list)
  if (n_sets < 2L || n_sets > 5L) return(NULL)
  fill_vec <- if (!is.null(circle_fills) && length(circle_fills) == n_sets) unname(circle_fills[set_names]) else get_circle_fills(n_sets)
  base_venn <- tryCatch({
    VennDiagram::venn.diagram(
      x = sets_named_list,
      category.names = set_names,
      filename = NULL,
      lwd = 2, lty = "blank", col = "black",
      fill = fill_vec, alpha = 0.5,
      cex = 1.2, fontfamily = "sans",
      cat.cex = 1.2, cat.fontface = "bold", cat.fontfamily = "sans",
      margin = 0.3,
      disable.logging = TRUE
    )
  }, error = function(e) NULL)
  if (!is.null(base_venn)) {
    grid::grid.draw(base_venn)
    if (!is.null(highlight_region) && nzchar(highlight_region)) {
      parts <- strsplit(highlight_region, "&")[[1]]
      bitmask <- paste0(sapply(set_names, function(snm) if (snm %in% parts) "1" else "0"), collapse = "")
      for (gname in names(base_venn)) {
        if (grepl(paste0("region\\.", bitmask, "$"), gname)) {
          grob_region <- base_venn[[gname]]
          grid::grid.draw(grid::editGrob(grob_region, gp = grid::gpar(fill = "#FF6666", alpha = 0.5, col = "red", lwd = 2)))
        }
      }
    }
  }
  base_venn
}

build_upset_ggplot <- function(counts_df, region_colors, x_lab, y_lab) {
  region_levels <- counts_df[["region"]]
  ggplot2::ggplot(counts_df, ggplot2::aes(x = factor(region, levels = region_levels), y = count, fill = region)) +
    ggplot2::geom_col(colour = "black", linewidth = 0.9) +
    ggplot2::scale_fill_manual(values = region_colors, guide = "none") +
    ggplot2::geom_text(ggplot2::aes(label = count), vjust = -0.5, size = 3.5) +
    ggplot2::scale_x_discrete(labels = counts_df[["label"]]) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      axis.title.y = ggplot2::element_text(size = 11),
      axis.title.x = ggplot2::element_text(size = 11),
      plot.margin = ggplot2::margin(55, 45, 45, 45)
    ) +
    ggplot2::labs(x = x_lab, y = y_lab)
}

# Export Venn and UpSet to files (no Shiny). Saves Venn.png, Venn.svg, UpSet.png, UpSet.svg, UpSet_counts.csv, UpSet_genes_by_region.csv.
export_venn_upset <- function(sets_named_list, subdir = "") {
  set_names <- names(sets_named_list)
  N <- length(set_names)
  if (N < 2) {
    cat("You must supply at least 2 sets.\n")
    return(invisible(NULL))
  }
  if (N > 5) warning("VennDiagram supports up to 5 sets. Diagrams for N>5 may not render correctly.")
  upset_data <- build_upset_data(sets_named_list)
  df_counts <- upset_data$counts_df
  region_to_genes <- upset_data$region_to_genes
  region_ids_present <- upset_data$region_ids_present
  circle_fills <- get_circle_fills(N)
  names(circle_fills) <- set_names
  region_colors <- setNames(
    vapply(df_counts$region, function(rid) get_region_fill_color(rid, set_names, circle_fills), character(1)),
    df_counts$region
  )
  grDevices::png(get_venn_upset_path("Venn.png", subdir), width = 800, height = 800, bg = "white")
  tryCatch({ create_venn_diagram(sets_named_list, set_names, circle_fills, NULL) }, finally = grDevices::dev.off())
  grDevices::svg(get_venn_upset_path("Venn.svg", subdir), width = 10, height = 10, bg = "white")
  tryCatch({ create_venn_diagram(sets_named_list, set_names, circle_fills, NULL) }, finally = grDevices::dev.off())
  g <- build_upset_ggplot(df_counts, region_colors, "Intersection of Sets", "Number of Genes")
  ggplot2::ggsave(get_venn_upset_path("UpSet.png", subdir), plot = g, width = 10, height = 6.5, units = "in", dpi = 150, limitsize = FALSE)
  ggplot2::ggsave(get_venn_upset_path("UpSet.svg", subdir), plot = g, device = "svg", width = 10, height = 6.5, units = "in", limitsize = FALSE)
  write.csv(df_counts, get_venn_upset_path("UpSet_counts.csv", subdir), row.names = FALSE)
  long <- do.call(rbind, lapply(region_ids_present, function(rid) {
    genes <- region_to_genes[[rid]]
    if (length(genes) == 0) return(NULL)
    data.frame(Region = rid, Gene = sort(genes), stringsAsFactors = FALSE)
  }))
  if (is.null(long) || nrow(long) == 0) long <- data.frame(Region = character(0), Gene = character(0), stringsAsFactors = FALSE)
  write.csv(long, get_venn_upset_path("UpSet_genes_by_region.csv", subdir), row.names = FALSE)
  cat("Venn and UpSet exported to", file.path(output_dir, "venn_upset", subdir), "\n")
  invisible(NULL)
}

export_venn_upset_split <- function(up_list, down_list, subdir_prefix = "") {
  subdir_up   <- if (nzchar(subdir_prefix)) file.path(subdir_prefix, "up")   else "up"
  subdir_down <- if (nzchar(subdir_prefix)) file.path(subdir_prefix, "down") else "down"
  if (length(up_list) >= 2) export_venn_upset(up_list, subdir = subdir_up)
  if (length(down_list) >= 2) export_venn_upset(down_list, subdir = subdir_down)
  if (length(up_list) < 2 && length(down_list) < 2) cat("Need at least two UP or two DOWN sets to export.\n")
  invisible(NULL)
}

# Function to plot PCA with custom colors and save as PNG and SVG
plotPCA_custom <- function(set,
                           group,
                           colors,
                           prefix = "PCA_plot",
                           width  = 6,
                           height = 6,
                           dpi    = 300){

  # ——————————————————————
  # 1) Extract & transform exactly like plotPCA(SeqExpressionSet)
  # ——————————————————————
  if (inherits(set, "SeqExpressionSet") ||
      inherits(set, "ExpressionSet")) {
    # EDASeq/limma‐style SeqExpressionSet
    mat0 <- Biobase::exprs(set)
    mat  <- log2(mat0 + 1)
  } else if (inherits(set, "SummarizedExperiment")) {
    # DESeqTransform or SE
    mat <- SummarizedExperiment::assay(set)
  } else {
    cat("Unsupported object class for PCA: ", class(set)[1], ". Skipping PCA plot.\n")
    return(invisible(NULL))
  }

  # ——————————————————————
  # 2) PCA on *all* rows, same centering/scaling as plotPCA()
  # ——————————————————————
  pca <- prcomp(t(mat), center=TRUE, scale.=FALSE)
  pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)[1:2]

  # ——————————————————————
  # 3) Build DF & ggplot styling (unchanged)
  # ——————————————————————
  df <- data.frame(
    PC1   = pca$x[,1],
    PC2   = pca$x[,2],
    name  = colnames(mat),
    group = factor(group, levels = names(colors))
  )

  p_base <- ggplot(df, aes(PC1, PC2, fill=group)) +
    geom_point(shape=21, colour="black", stroke=0.5, size=3) +
    scale_fill_manual(values=colors) +
    labs(
      x = sprintf("PC1 (%.1f%%)", pct[1]),
      y = sprintf("PC2 (%.1f%%)", pct[2])
    ) +
  theme_classic() +
  theme(
    panel.border      = element_rect(colour="black", fill=NA, size=1),
    plot.margin       = unit(c(1, 4, 1, 1), "lines"),      # top, right, bottom, left
    legend.position     = c(1.02, 0.9),
    legend.justification= c(0, 1),
    legend.margin     = margin(5, 5, 5, 5),
    legend.background    = element_rect(fill="white", colour="black", size=0.5)
    # legend.box.background = element_rect(fill="white", colour="black", size=0.5),
    # legend.box.margin    = margin(5,5,5,5)
    ) +
    guides(fill = guide_legend(
        keywidth  = unit(0.4, "lines"),
        keyheight = unit(0.4, "lines"),
        override.aes = list(shape=21, size=3, stroke=0.5),
        label.theme  = element_text(size = 6)
    )) +
    theme(
        legend.key.size = unit(0.4, "lines"),
        legend.text     = element_text(size = 6),
        legend.title    = element_text(size = 8)
    )


  # ——————————————————————
  # 4) Save dots‐only & dots+labels
  # ——————————————————————
  ggsave(paste0(prefix, "_dots.png"), p_base,
         width=width, height=height,
         dpi=dpi, units="in", limitsize=FALSE)
  ggsave(paste0(prefix, "_dots.svg"), p_base,
         device="svg", width=width, height=height,
         units="in", limitsize=FALSE)

  p_lab <- p_base +
    geom_text(aes(label=name),
              fontface="bold", colour="black",
              nudge_y=0.02)

  ggsave(paste0(prefix, "_labels.png"), p_lab,
         width=width, height=height,
         dpi=dpi, units="in", limitsize=FALSE)
  ggsave(paste0(prefix, "_labels.svg"), p_lab,
         device="svg", width=width, height=height,
         units="in", limitsize=FALSE)
}


###################################################################################################
# **Main Program Start**
# Start the main loop
repeat {
    repeat {
        # Prompt user for output directory
        cat("Please provide a directory to save all exported files (or type 'quit' to exit):\n")
        output_dir <- handle_quit(readline())
        
        # Handle "back" command
        if (output_dir == "back") {
            cat("Already at the first step. Cannot go back further.\n")
            next
        }

        # Ensure the directory is valid
        dir_error <- NULL
        tryCatch({
            # Convert to absolute path
            full_path <- normalizePath(output_dir, mustWork = FALSE)
            parent_dir <- dirname(full_path)
            
            if (!dir.exists(parent_dir)) {
                dir_error <<- sprintf("The parent directory '%s' does not exist. Please provide a valid path.", parent_dir)
            } else if (!dir.exists(full_path)) {
                if (file.access(parent_dir, 2) != 0) {
                    dir_error <<- "Parent directory is not writable. Please choose another path."
                } else {
                    dir.create(full_path, recursive = TRUE)
                    cat(sprintf("Output directory created at: %s\n", full_path))
                    break
                }
            } else {
                if (file.access(full_path, 2) != 0) {
                    dir_error <<- "Directory exists but is not writable. Please choose another directory."
                } else {
                    cat(sprintf("Using existing directory: %s\n", full_path))
                    break
                }
            }
        }, error = function(e) {
            dir_error <<- e$message
        })
        if (!is.null(dir_error)) {
            cat(sprintf("Error: %s\n", dir_error))
            cat("Please try again with a valid directory path.\n")
            next
        }
    } # End of directory input loop

    # **Menu-Based Analysis Start**
    repeat {
        cat("\nMain Menu\n")
        cat("[1] Load DEGs tables and analyze\n")
        cat("[2] Process count files and perform differential expression analysis\n")
        cat("[3] Load genomic data tables and analyze\n")
        cat("[4] Quit\n")
        cat("Please select an option by entering the corresponding number (or type 'back' to return to the previous step):\n")
        
        menu_choice <- handle_quit(readline())
        
        # Handle "back" command
        if (menu_choice == "back") {
            cat("Returning to the previous step.\n")
            break # Exit to go back to Step 1
        }

        # Validate menu choice
        menu_choice <- as.numeric(menu_choice)
        
        if (is.na(menu_choice) || menu_choice < 1 || menu_choice > 4) {
            cat("Invalid choice. Please enter a valid number between 1 and 3.\n")
            next
        }
        
        if (menu_choice == 4) {
            cat("Exiting the program. Goodbye!\n")
            do_quit()
        }

        if (menu_choice == 1) {
            repeat {
                # **Option 1: Load DEGs Tables**
                cat("Please provide the path to the DEGs file (or type 'quit' to exit or type 'back' to return to the main menu):\n")
                degs_file_path <- handle_quit(readline())

                if (degs_file_path == "back") {
                    cat("Returning to the main menu.\n")
                    break # Return to the main menu
                }    

                if (!file.exists(degs_file_path)) {
                    cat("Error: The DEGs file does not exist. Returning to the menu.\n")
                    next
                }
                tryCatch({
                    # Determine the file type based on extension
                    if (grepl("\\.csv$", degs_file_path, ignore.case = TRUE)) {
                        degs_data <- read.csv(degs_file_path, header = TRUE, stringsAsFactors = FALSE)
                    } else if (grepl("\\.txt$|\\.tsv$", degs_file_path, ignore.case = TRUE)) {
                        degs_data <- read.delim(degs_file_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t") 
                    } else if (grepl("\\.xlsx$", degs_file_path, ignore.case = TRUE)) {
                        degs_data <- read_excel(degs_file_path, col_names = TRUE)
                        # Ensure consistency with CSV reading
                        degs_data <- as.data.frame(degs_data, stringsAsFactors = FALSE)  # Convert to data.frame
                        colnames(degs_data) <- trimws(colnames(degs_data))               # Trim whitespace from headers
                        degs_data[] <- lapply(degs_data, function(x) {                   # Trim whitespace from data
                            if (is.character(x)) trimws(x) else x
                        })
                        degs_data[degs_data == ""] <- NA                                # Handle empty cells as NA
                    } else {
                        cat("Unsupported file format. Please provide a file in CSV, TXT, TSV, or XLSX format.\n")
                        next # Retry the input
                    }
                    
                    # Ensure the structure of the loaded DEGs data
                    required_columns <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "comparison", "Gene")
                    missing_columns <- setdiff(required_columns, colnames(degs_data))
                    
                    if (length(missing_columns) > 0) {
                        cat(sprintf("Error: The DEGs file is missing the following required columns: %s\n", paste(missing_columns, collapse = ", ")))
                        next # Retry the input
                    }
                    
                    # Check for duplicate rows (optional informational message)
                    if (anyDuplicated(degs_data)) {
                        cat("Warning: Duplicate rows detected in the DEGs file. Proceeding with all entries.\n")
                    }
                    
                    cat("DEGs table loaded successfully:\n")
                    print(head(degs_data))
                    break # Exit loop after successful load
                }, error = function(e) {
                    cat("Error reading the DEGs file: ", e$message, "\nReturning to the menu.\n")
                    next
                })
            }

            repeat {
                # Step 2: Ask for counts data and coldata files
                cat("Please provide the path to the counts data file (or type 'quit' to exit):\n")
                counts_data_path <- handle_quit(readline())

                # Handle "back" command
                if (counts_data_path == "back") {
                    cat("Returning to the DEGs file step.\n")
                    break # Return to the previous step
                }

                # Load counts data 
                tryCatch({
                    if (grepl("\\.csv$", counts_data_path, ignore.case = TRUE)) {
                        counts_data <- read.csv(counts_data_path, header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
                    } else if (grepl("\\.txt$|\\.tsv$", counts_data_path, ignore.case = TRUE)) {
                        counts_data <- read.delim(counts_data_path, header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = "\t")
                    } else if (grepl("\\.xlsx$", counts_data_path, ignore.case = TRUE)) {
                        counts_data <- read_excel(counts_data_path, col_names = TRUE)
                        # Ensure consistency with CSV reading
                        counts_data <- as.data.frame(counts_data, stringsAsFactors = FALSE)  #
                        colnames(counts_data) <- trimws(colnames(counts_data))               # Trim whitespace from headers
                        counts_data[] <- lapply(counts_data, function(x) {                   # Trim whitespace
                            if (is.character(x)) trimws(x) else x
                        })
                        counts_data[counts_data == ""] <- NA                                # Handle empty cells as NA
                    } else {
                        cat("Unsupported counts data file format. Please provide a file in CSV, TXT, TSV, or XLSX format.\n")
                        next
                    }
                    cat("Counts data loaded successfully:\n")
                    print(head(counts_data))
                    break # Exit the loop on success

                }, error = function(e) {
                    cat("Error reading the counts data file: ", e$message, "\n")
                    next
                })
            }

            repeat{
                cat("Please provide the path to the coldata file (or type 'quit' to exit):\n")
                coldata_path <- handle_quit(readline())

                # Handle "back" command
                if (coldata_path == "back") {
                    cat("Returning to the counts data step.\n")
                    break # Return to the previous step
                }
                tryCatch({
                    if (grepl("\\.csv$", coldata_path, ignore.case = TRUE)) {
                        coldata <- read.csv(coldata_path, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
                    } else if (grepl("\\.txt$|\\.tsv$", coldata_path, ignore.case = TRUE)) {
                        coldata <- read.delim(coldata_path, header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = "\t")
                    } else if (grepl("\\.xlsx$", coldata_path, ignore.case = TRUE)) {
                        coldata <- read_excel(coldata_path, col_names = TRUE)
                    } else {
                        cat("Unsupported coldata file format. Please provide a file in CSV, TXT, TSV, or XLSX format.")
                        next
                    }
                    cat("Coldata loaded successfully:\n")
                    print(head(coldata))
                    break # Exit the loop on success
                }, error = function(e) {
                    cat("Error reading coldata: ", e$message, "\n")
                    next
                })
            }

            # Proceed to downstream analysis menu
            # **Downstream Analysis Menu**
            # User chooses automatic or manual analysis
            cat("\nSelect analysis mode:\n")
            cat("[1] Automatic Streamlined Analysis (Heatmap, Venn, and Trending Plots)\n")
            cat("[2] Manual Analysis\n")
            analysis_mode_raw <- handle_quit(readline("Enter your choice (1 or 2), or type 'back' to return to main menu: "))
            if (analysis_mode_raw == "back") break
            analysis_mode <- as.numeric(analysis_mode_raw)

            if (analysis_mode == 1) {
                # **Automatic Analysis Pipeline**
                cat("\nRunning Automatic Streamlined Analysis...\n")
                # Ensure DEGs data is available
                if (nrow(degs_data) == 0) {
                    cat("No DEGs data available for generating heatmaps.\n")
                    next
                }
                padj_threshold <- 0.05
                cat(sprintf("Using padj threshold: %.2f.\n", padj_threshold))
                # Filter significant genes once and reuse
                significant_genes <- degs_data[degs_data$padj < padj_threshold, , drop = FALSE]
                significant_genes <- significant_genes[complete.cases(significant_genes), ]

                if (nrow(significant_genes) == 0) {
                    cat("No significant genes found for the selected padj threshold.\n")
                    next
                }

                # Check alignment between coldata row names and counts_data column names
                missing_samples <- setdiff(rownames(coldata), colnames(counts_data))
                if (length(missing_samples) > 0) {
                    cat("The following samples are present in coldata but missing in counts_data:\n")
                    print(missing_samples)
                    cat("Please ensure the coldata row names match the counts_data column names.\n")
                    next
                }

                aligned_samples <- intersect(rownames(coldata), colnames(counts_data))
                cat("Aligned samples between coldata and counts_data:\n")
                print(aligned_samples)

                if (length(aligned_samples) == 0) {
                    cat("No aligned samples found between coldata and counts_data. Please verify the input files.\n")
                    break
                }
                # Ask user to select a grouping column
                cat("Available columns in coldata:\n")
                for (i in seq_along(colnames(coldata))) {
                    cat(sprintf("[%d] %s\n", i, colnames(coldata)[i]))
                }

                repeat {
                    cat("Please select the column to group samples by (enter column name or index, or 'back' to return to main menu):\n")
                    group_input <- handle_quit(readline())
                    if (group_input == "back") break
                    if (suppressWarnings(!is.na(as.numeric(group_input)))) {
                        col_index <- as.numeric(group_input)
                    if (col_index >= 1 && col_index <= ncol(coldata)) {
                        group_col <- colnames(coldata)[col_index]
                        break
                    }
                    } else if (group_input %in% colnames(coldata)) {
                        group_col <- group_input
                        break
                    }
                    cat("Invalid input. Please enter a valid column name or index.\n")
                }
                if (group_input == "back") next

                cat(sprintf("Grouping samples by column: '%s'\n", group_col))
                # Extract unique conditions from coldata dynamically
                unique_conditions <- unique(as.vector(unlist(coldata)))
                cat("Unique conditions extracted from coldata:\n")
                print(unique_conditions)

                if (length(unique_conditions) == 0) {
                    cat("No unique conditions found in coldata. Please verify your coldata file.\n")
                    next
                }

                # Display unique conditions with indices
                cat("Available unique conditions in your experiment:\n")
                for (i in seq_along(unique_conditions)) {
                    cat(sprintf("[%d] %s\n", i, unique_conditions[i]))
                }

                # Prompt user to select the control condition by index
                control_index_input <- NULL
                repeat {
                    cat("Select the index of the control condition (or type 'back' to return to group selection):\n")
                    control_index_input <- handle_quit(readline())
                    if (control_index_input == "back") break
                    control_index <- as.numeric(control_index_input)
                    if (!is.na(control_index) && control_index >= 1 && control_index <= length(unique_conditions)) {
                        control_condition <- unique_conditions[control_index]
                        cat(sprintf("Selected control condition: %s\n", control_condition))
                        break
                    } else {
                        cat("Invalid index. Please select a valid index from the list.\n")
                    }
                }
                if (control_index_input == "back") next

                # Create logged data matrix
                tryCatch({
                    control_samples <- rownames(coldata)[apply(coldata, 1, function(row) control_condition %in% row)]
                    cat("Control samples identified for the selected control condition:\n")
                    print(control_samples)

                    if (length(control_samples) == 0) {
                        cat(sprintf("No samples found for the selected control condition: %s", control_condition))
                        next
                    }

                    control_samples_in_counts <- intersect(control_samples, colnames(counts_data))
                    if (length(control_samples_in_counts) == 0) {
                        cat(sprintf("None of the control samples are present in counts_data. Please verify the input files."))
                        next
                    }
                    if (!is.matrix(counts_data)) counts_data <- data.matrix(counts_data)
                    
                    cont_avg <- rowMeans(counts_data[, control_samples_in_counts, drop = FALSE], na.rm = TRUE)
                    cat("Control condition averages computed:\n")
                    print(head(cont_avg))

                    # div_df <- counts_data / cont_avg
                    logged_df <- log2((as.matrix(counts_data) + 1) / (cont_avg + 1)) # Add 1 to avoid division by zero or log of zero

                    # Remove rows with NA in logged_df
                    logged_df <- logged_df[complete.cases(logged_df), ]
                    cat("Final logged data dimensions after removing NA rows:\n")
                    print(dim(logged_df))
                    # # Construct the full path for the logged data file
                    # logged_df_file <- get_output_path("logged_df.csv")
                    # # Save the logged data matrix
                    # write.csv(logged_df, file = logged_df_file, row.names = TRUE)
                    # # Inform the user about the file's location
                    # cat("Logged data matrix saved to:", logged_df_file, "\n")

                    ################## for haetmap_vd ##################
                    degs_data <- degs_data %>%
                    mutate(
                        Gene       = as.character(Gene),
                        comparison = as.character(comparison)
                    )

                    best_per_gene <- degs_data %>%
                        filter(!is.na(padj), padj < padj_threshold) %>%
                        group_by(Gene) %>%
                        slice_max(abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
                        ungroup() %>%
                        mutate(
                            comparison = paste0(
                            comparison,
                            ifelse(log2FoldChange > 0, "_up", "_down")
                            )
                        ) %>%
                        select(Gene, comparison)

                    # 3. Bring your logged_df rownames into a “Gene” column
                    logged_df_comp <- data.frame(
                        Gene = rownames(logged_df),
                        logged_df,
                        stringsAsFactors = FALSE,
                        check.names = FALSE
                    )

                    # 4. Join on the best comparison per gene
                    logged_df_comp <- left_join(logged_df_comp, best_per_gene, by = "Gene")

                    # 5. (Optional) restore Gene to rownames
                    rownames(logged_df_comp) <- logged_df_comp$Gene
                    logged_df_comp$Gene <- NULL

                    # remove any rows with NA in the comparison column
                    logged_df_comp <- logged_df_comp[!is.na(logged_df_comp$comparison), ]

                    # Inspect
                    head(logged_df_comp)
                    write.csv(logged_df_comp, file = get_output_path("logged_df_comp.csv"), row.names = TRUE)
                    # print how many genes in each group in comp
                    cat("Comparison groups in logged_df_comp:\n")
                    print(table(logged_df_comp$comparison))


                    scored_df <- as.matrix(counts_data)
                    scored_df_comp <- data.frame(
                        Gene = rownames(scored_df),
                        scored_df,
                        stringsAsFactors = FALSE,
                        check.names = FALSE
                    )
                    scored_df_comp <- left_join(scored_df_comp, best_per_gene, by = "Gene")

                    rownames(scored_df_comp) <- scored_df_comp$Gene
                    scored_df_comp$Gene <- NULL

                    scored_df_comp <- scored_df_comp[!is.na(scored_df_comp$comparison), ]

                }, error = function(e) {
                    cat("Error creating logged data matrix. Please ensure the control condition is valid.\n")
                    print(e)
                    next
                })

                # Generate heatmap
                cat("Generating Heatmap...\n")
                png(get_heatmaps_path("log2", "combined", "sum_heatmap_auto_vd.png"), width = 400, height = 900)

                ht <- heatmap_vd(
                        logged_df_comp,
                        plot_incl=c(1:(ncol(logged_df_comp)-1)),
                        kmeans_split=logged_df_comp$comparison, 
                        col=mako(100),legend_range=1,col_bias=1,lim1=NULL,lim2=NULL)
                ComplexHeatmap::draw(ht)
                dev.off()

                png(get_heatmaps_path("z_score", "combined", "sum_heatmap_auto_z.png"), width = 400, height = 900)

                ht <- heatmap_vd_z(
                        scored_df_comp,
                        plot_incl=c(1:(ncol(logged_df_comp)-1)),
                        kmeans_split=logged_df_comp$comparison, 
                        col=mako(100),legend_range=1,col_bias=1,lim1=NULL,lim2=NULL)
                ComplexHeatmap::draw(ht)
                dev.off()

                tryCatch({
                    # Per-comparison and sum up/down heatmaps (manual-aligned: get_logged_for_comparison, both scales, get_heatmaps_path, sum_style, CSV)
                    cat("Filtering upregulated and downregulated genes...\n")
                    heatmap_palette <- viridis::mako(101)
                    comparison_groups <- unique(significant_genes$comparison)
                    upregulated_heatmaps <- list()
                    downregulated_heatmaps <- list()
                    for (gr in comparison_groups) {
                        logged_gr <- get_logged_for_comparison(gr, counts_data, coldata, group_col)
                        logged_gr <- logged_gr[complete.cases(logged_gr), , drop = FALSE]
                        group_genes_gr <- significant_genes[significant_genes$comparison == gr, "Gene"]
                        group_heatmap_data <- logged_gr[rownames(logged_gr) %in% group_genes_gr, , drop = FALSE]
                        group_heatmap_data <- as.matrix(group_heatmap_data)
                        if (nrow(group_heatmap_data) == 0) next
                        if (any(!is.finite(group_heatmap_data))) {
                            group_heatmap_data <- apply(group_heatmap_data, 2, function(col) {
                                col[!is.finite(col)] <- median(col, na.rm = TRUE); return(col)
                            })
                        }
                        upregulated_genes <- significant_genes[significant_genes$comparison == gr & significant_genes$log2FoldChange > 0, "Gene"]
                        downregulated_genes <- significant_genes[significant_genes$comparison == gr & significant_genes$log2FoldChange < 0, "Gene"]
                        upregulated_overlap <- intersect(upregulated_genes, rownames(group_heatmap_data))
                        downregulated_overlap <- intersect(downregulated_genes, rownames(group_heatmap_data))
                        if (length(upregulated_overlap) > 0) {
                            upregulated_data <- group_heatmap_data[upregulated_overlap, , drop = FALSE]
                            if (nrow(upregulated_data) > 0) {
                                upregulated_heatmaps[[length(upregulated_heatmaps) + 1]] <- list(heatmap_data = upregulated_data, title = gr)
                            }
                        }
                        if (length(downregulated_overlap) > 0) {
                            downregulated_data <- group_heatmap_data[downregulated_overlap, , drop = FALSE]
                            if (nrow(downregulated_data) > 0) {
                                downregulated_heatmaps[[length(downregulated_heatmaps) + 1]] <- list(heatmap_data = downregulated_data, title = gr)
                            }
                        }
                    }
                    for (heatmap_scale in c("log2", "z_score")) {
                        for (gr in comparison_groups) {
                            logged_gr <- get_logged_for_comparison(gr, counts_data, coldata, group_col)
                            logged_gr <- logged_gr[complete.cases(logged_gr), , drop = FALSE]
                            group_genes_gr <- significant_genes[significant_genes$comparison == gr, "Gene"]
                            group_heatmap_data <- logged_gr[rownames(logged_gr) %in% group_genes_gr, , drop = FALSE]
                            if (nrow(group_heatmap_data) == 0) next
                            group_heatmap_data <- as.matrix(group_heatmap_data)
                            if (any(!is.finite(group_heatmap_data))) {
                                group_heatmap_data <- apply(group_heatmap_data, 2, function(col) {
                                    col[!is.finite(col)] <- median(col, na.rm = TRUE); return(col)
                                })
                            }
                            upregulated_genes <- significant_genes[significant_genes$comparison == gr & significant_genes$log2FoldChange > 0, "Gene"]
                            downregulated_genes <- significant_genes[significant_genes$comparison == gr & significant_genes$log2FoldChange < 0, "Gene"]
                            upregulated_overlap <- intersect(upregulated_genes, rownames(group_heatmap_data))
                            downregulated_overlap <- intersect(downregulated_genes, rownames(group_heatmap_data))
                            if (length(upregulated_overlap) > 0) {
                                upregulated_data <- group_heatmap_data[upregulated_overlap, , drop = FALSE]
                                if (nrow(upregulated_data) > 0) {
                                    upregulated_data_with_names <- data.frame(Gene = rownames(upregulated_data), upregulated_data, check.names = FALSE)
                                    heatmap_obj_up <- if (heatmap_scale == "z_score") heatmap_z(upregulated_data, plot_incl = seq_len(ncol(upregulated_data)), kmeans_split = NULL, palette = heatmap_palette) else heatmap_v(upregulated_data, plot_incl = seq_len(ncol(upregulated_data)), kmeans_split = NULL, palette = heatmap_palette)
                                    write.csv(upregulated_data_with_names, file = get_heatmaps_path(heatmap_scale, "up", sprintf("heatmap_upregulated_%s.csv", gsub(" ", "_", gr))), row.names = FALSE)
                                    upregulated_png_file <- get_heatmaps_path(heatmap_scale, "up", sprintf("heatmap_upregulated_%s.png", gsub(" ", "_", gr)))
                                    if (!is.null(heatmap_obj_up)) { png(upregulated_png_file, width = 900, height = 500); ComplexHeatmap::draw(heatmap_obj_up, heatmap_legend_side = "right"); dev.off() }
                                }
                            }
                            if (length(downregulated_overlap) > 0) {
                                downregulated_data <- group_heatmap_data[downregulated_overlap, , drop = FALSE]
                                if (nrow(downregulated_data) > 0) {
                                    downregulated_data_with_names <- data.frame(Gene = rownames(downregulated_data), downregulated_data, check.names = FALSE)
                                    heatmap_obj_down <- if (heatmap_scale == "z_score") heatmap_z(downregulated_data, plot_incl = seq_len(ncol(downregulated_data)), kmeans_split = NULL, palette = heatmap_palette) else heatmap_v(downregulated_data, plot_incl = seq_len(ncol(downregulated_data)), kmeans_split = NULL, palette = heatmap_palette)
                                    write.csv(downregulated_data_with_names, file = get_heatmaps_path(heatmap_scale, "down", sprintf("heatmap_downregulated_%s.csv", gsub(" ", "_", gr))), row.names = FALSE)
                                    downregulated_png_file <- get_heatmaps_path(heatmap_scale, "down", sprintf("heatmap_downregulated_%s.png", gsub(" ", "_", gr)))
                                    if (!is.null(heatmap_obj_down)) { png(downregulated_png_file, width = 900, height = 500); ComplexHeatmap::draw(heatmap_obj_down, heatmap_legend_side = "right"); dev.off() }
                                }
                            }
                        }
                        if (length(upregulated_heatmaps) > 0) {
                            combined_upregulated_data <- do.call(rbind, lapply(upregulated_heatmaps, function(hm) hm$heatmap_data))
                            row_ann_up <- unlist(lapply(upregulated_heatmaps, function(hm) rep(hm$title, nrow(hm$heatmap_data))))
                            row_split_upregulated <- factor(row_ann_up, levels = unique(row_ann_up))
                            heatmap_obj_upregulated <- if (heatmap_scale == "z_score") heatmap_z(combined_upregulated_data, plot_incl = seq_len(ncol(combined_upregulated_data)), row_split = row_split_upregulated, palette = heatmap_palette, sum_style = TRUE) else heatmap_v(combined_upregulated_data, plot_incl = seq_len(ncol(combined_upregulated_data)), row_split = row_split_upregulated, palette = heatmap_palette, sum_style = TRUE)
                            png(get_heatmaps_path(heatmap_scale, "up", "sum_upregulated_heatmap.png"), width = 1600, height = 2400)
                            if (!is.null(heatmap_obj_upregulated)) ComplexHeatmap::draw(heatmap_obj_upregulated, heatmap_legend_side = "right")
                            dev.off()
                            write.csv(data.frame(Gene = rownames(combined_upregulated_data), combined_upregulated_data, check.names = FALSE), get_heatmaps_path(heatmap_scale, "up", "sum_upregulated_heatmap.csv"), row.names = FALSE)
                        }
                        if (length(downregulated_heatmaps) > 0) {
                            combined_downregulated_data <- do.call(rbind, lapply(downregulated_heatmaps, function(hm) hm$heatmap_data))
                            row_ann_down <- unlist(lapply(downregulated_heatmaps, function(hm) rep(hm$title, nrow(hm$heatmap_data))))
                            row_split_downregulated <- factor(row_ann_down, levels = unique(row_ann_down))
                            heatmap_obj_downregulated <- if (heatmap_scale == "z_score") heatmap_z(combined_downregulated_data, plot_incl = seq_len(ncol(combined_downregulated_data)), row_split = row_split_downregulated, palette = heatmap_palette, sum_style = TRUE) else heatmap_v(combined_downregulated_data, plot_incl = seq_len(ncol(combined_downregulated_data)), row_split = row_split_downregulated, palette = heatmap_palette, sum_style = TRUE)
                            png(get_heatmaps_path(heatmap_scale, "down", "sum_downregulated_heatmap.png"), width = 1600, height = 2400)
                            ComplexHeatmap::draw(heatmap_obj_downregulated, heatmap_legend_side = "right")
                            dev.off()
                            write.csv(data.frame(Gene = rownames(combined_downregulated_data), combined_downregulated_data, check.names = FALSE), get_heatmaps_path(heatmap_scale, "down", "sum_downregulated_heatmap.csv"), row.names = FALSE)
                        }
                    }
                    cat("Heatmap generation completed.\n")
                }, error = function(e) {
                    cat(sprintf("Error generating heatmaps: %s\n", e$message))
                })

                # Generate Venn diagrams
                cat("Generating Venn Diagrams...\n")
                # Filter DEGs based on padj
                degs_data <- degs_data[!is.na(degs_data$padj) & degs_data$padj < padj_threshold, ]

                # Check if any significant genes exist
                if (nrow(degs_data) == 0) {
                    cat("No significant genes found for the given padj threshold. Returning to menu.\n")
                    next
                }
                # Split into upregulated and downregulated genes
                degs_data_up <- subset(degs_data, log2FoldChange > 0)
                degs_data_down <- subset(degs_data, log2FoldChange < 0)
                # Generate indexed list of comparisons
                comparisons <- unique(degs_data$comparison)
                options_list <- list()
                for (i in seq_along(comparisons)) {
                    options_list[[paste0(i)]] <- comparisons[i]
                    options_list[[paste0(i, "_up")]] <- paste0(comparisons[i], "_up")
                    options_list[[paste0(i, "_down")]] <- paste0(comparisons[i], "_down")
                }
                venn_data <- list()
                venn_data[["Upregulated"]] <- degs_data_up
                venn_data[["Downregulated"]] <- degs_data_down
                for (category in names(venn_data)) {
                    category_data <- venn_data[[category]]
                    if (nrow(category_data) == 0) {
                        cat(sprintf("No genes in %s category. Skipping Venn diagrams for this category.\n", category))
                        next
                    }
                    # Subset comparisons dynamically
                    subset_genes <- split(category_data$Gene, category_data$comparison)
                    # Add total gene counts to the category names
                    category_names_with_counts <- sapply(names(subset_genes), function(name) {
                        paste0(name, " (", length(subset_genes[[name]]), ")")
                    })
                    # Save under venn_upset/up/ or venn_upset/down/
                    venn_subdir <- if (category == "Upregulated") "up" else "down"
                    venn_filename <- get_venn_upset_path(paste0("default_venn_", category, ".png"), venn_subdir)
                    table_filename <- get_venn_upset_path(paste0("default_venn_table_", category, ".csv"), venn_subdir)
                    
                    tryCatch({
                        # Save the Venn diagram
                        png(venn_filename, width = 600, height = 600)
                        venn.plot <- venn.diagram(
                            x = subset_genes,
                            category.names = category_names_with_counts,
                            # category.names = names(subset_genes),
                            filename = NULL,
                            lwd = 2,
                            lty = 'blank',
                            col = "transparent",
                            fill = colorRampPalette(brewer.pal(8, "Pastel2"))(length(subset_genes)),
                            # alpha = 0.5,
                            cex = 1.5,
                            fontface = "bold",
                            fontfamily = "sans",
                            cat.cex = 1,
                            cat.default.pos = "outer",
                            cat.fontface = "bold",
                            cat.fontfamily = "sans",
                            # rotation = 1,
                            main = paste("Default Venn Diagram:", category),
                            disable.logging = TRUE
                        )
                        grid.draw(venn.plot)
                        dev.off()
                        cat("Default Venn diagram saved as:", venn_filename, "\n")

                        # Generate and export detailed gene table
                        generate_gene_table(subset_genes, names(subset_genes), table_filename)
                    }, error = function(e) {
                        cat("Error generating default Venn diagram:", e$message, "\n")
                    })
                }

                # Generate trending plots
                cat("Generating Trending Plots...\n")
                generate_default_trending_plots(
                    counts_data = counts_data,
                    coldata = coldata,
                    significant_genes = degs_data,
                    group_col = group_col,
                    padj_threshold = padj_threshold
                )
                cat("\nDefault trending plots completed.\n")
                cat("\nAutomatic Analysis Completed!\n")
            } # End of automatic analysis mode
            else if (analysis_mode == 2) {
                repeat {
                    cat("\nInteractive Analysis Menu\n")
                    cat("[1] Check if there are some genes of interest\n")
                    cat("[2] Generate a heatmap\n")
                    cat("[3] Create Venn and UpSet diagrams\n")
                    cat("[4] Create trending plots\n")
                    cat("[5] Quit\n")
                    cat("Please select an option by entering the corresponding number (or type 'back' to return to main menu):\n")
                    
                    user_choice_raw <- handle_quit(readline())
                    if (user_choice_raw == "back") break
                    user_choice <- as.numeric(user_choice_raw)
                    
                    if (is.na(user_choice) || user_choice < 1 || user_choice > 5) {
                    cat("Invalid choice. Please enter a valid number between 1 and 5.\n")
                    next
                    }
                    
                    if (user_choice == 5) {
                    cat("Exiting the program. Goodbye!\n")
                    break
                    }
                    
                    # **Analysis Options**
                    if (user_choice == 1) {
                        # **Option 1: Check Genes of Interest**
                        repeat {
                            cat("Insert the gene names (comma-separated, or provide a file path, or type 'quit' to exit):\n")
                            genes_input <- handle_quit(readline())
                            # genes_in_interest <- unlist(strsplit(gsub("\\s+", "", genes_input), ",")) 
                            # Check if the input is a valid file path
                            if (file.exists(genes_input)) {
                                cat("Detected file input. Extracting genes from the first column...\n")
                                if (grepl("\\.xlsx$", genes_input, ignore.case = TRUE)) {
                                    genes_df <- read_excel(genes_input, col_names = TRUE)
                                } else if (grepl("\\.(txt|tsv)$", matrix_path, ignore.case = TRUE)) {
                                    genes_df <- read.table(genes_input, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
                                } else if (grepl("\\.csv$", matrix_path, ignore.case = TRUE)) {
                                    genes_df <- read.csv(genes_input, header = TRUE, stringsAsFactors = FALSE)
                                } else {
                                    cat("Error: Unsupported file format. Please provide a .csv, .tsv, .txt, or .xlsx file.\n")
                                    next # Return to the start of the loop to retry
                                }

                                if (ncol(genes_df) < 1) {
                                    cat("Error: The file does not contain any columns. Please try again.\n")
                                    next
                                }
                                
                                genes_in_interest <- unique(na.omit(genes_df[[1]]))  # Extract first column and remove NAs
                            } else {
                                # Process manual input (comma-separated values)
                                genes_in_interest <- unlist(strsplit(gsub("\\s+", "", genes_input), ","))
                            }                   
                            if (length(genes_in_interest) == 0) {
                                cat("No valid genes entered. Please try again.\n")
                            } else {
                                break
                            }
                        } # End of repeat loop for gene input

                        cat("Searching for genes in the DEGs data and generating plots...\n")

                        # **Step 1: Display Available Columns in coldata**
                        cat("Available columns in coldata:\n")
                        for (i in seq_along(colnames(coldata))) {
                            cat(sprintf("[%d] %s\n", i, colnames(coldata)[i]))
                        }

                        # **Step 2: Ask user to select a grouping column**
                        repeat {
                            cat("Please select the column to group samples by (enter column name or index, or 'back' to cancel):\n")
                            group_input <- handle_quit(readline())
                            if (group_input == "back") break
                            if (suppressWarnings(!is.na(as.numeric(group_input)))) {
                                col_index <- as.numeric(group_input)
                                if (col_index >= 1 && col_index <= ncol(coldata)) {
                                    group_col <- colnames(coldata)[col_index]
                                    break
                                }
                            } else if (group_input %in% colnames(coldata)) {
                                group_col <- group_input
                                break
                            }
                            cat("Invalid input. Please enter a valid column name or index.\n")
                        } # End of repeat loop for grouping column selection
                        if (group_input == "back") next

                        cat(sprintf("Grouping samples by column: '%s'\n", group_col))

                        # **Step 3: Group samples by the selected condition column**
                        sample_conditions <- split(rownames(coldata), coldata[[group_col]])

                        # **Step 4: Ensure user input is valid**
                        if (length(sample_conditions) < 2) {
                            cat("Error: The selected grouping column does not contain enough unique conditions for plotting. Please select another column or fix the data.\n")
                            next
                        }

                        cat("Sample grouping completed.\n")
                        # Initialize a results list to store summary information
                        results_list <- list()

                        # Iterate over each gene of interest
                        for (gene in genes_in_interest) {
                            # Match all rows where the `Gene` column contains the gene of interest correctly
                            matched_rows <- degs_data[grepl(paste0("^", gene, "$|^", gene, "\\."), degs_data$Gene, ignore.case = TRUE), , drop = FALSE]

                            if (nrow(matched_rows) > 0) {
                                cat(sprintf("Gene %s found in %d comparisons:\n", gene, nrow(matched_rows)))

                                for (i in seq_len(nrow(matched_rows))) {
                                    row <- matched_rows[i, ]
                                    cat(sprintf("  Comparison: %s, padj = %.4f, log2FoldChange = %.4f\n",
                                                row["comparison"], as.numeric(row["padj"]), as.numeric(row["log2FoldChange"])))
                                }

                                # Store all matches
                                results_list[[gene]] <- matched_rows

                                # ---- Extract & Plot Expression Data for Each Unique Matched Gene ----
                                matched_gene_names <- unique(matched_rows$Gene)  # Get unique gene names found

                                for (matched_gene in matched_gene_names) {
                                    gene_lower <- tolower(matched_gene)
                                    rownames_lower <- tolower(rownames(counts_data))

                                    if (gene_lower %in% rownames_lower) {
                                        expr_data <- as.numeric(unlist(counts_data[rownames_lower == gene_lower, , drop = FALSE]))
                                        
                                        # Prevent duplicate plots
                                        if (!file.exists(get_genes_of_interest_path(paste0(matched_gene, "_expression_plot.png")))) {
                                            plot_gene_expression(matched_gene, expr_data, coldata, group_col, matched_rows)
                                            cat(sprintf("✅ Plot for %s saved.\n", matched_gene))
                                        } else {
                                            cat(sprintf("🚀 Skipping duplicate plot for %s\n", matched_gene))
                                        }
                                    } else {
                                        cat(sprintf("⚠️ Gene %s not found in counts data. Skipping plot...\n", matched_gene))
                                    }
                                }

                            } else {
                                cat(sprintf("Gene %s not found in the DEGs data.\n", gene))
                                # Create a placeholder row for unmatched genes
                                placeholder_row <- data.frame(
                                    Gene = gene,
                                    baseMean = NA,
                                    log2FoldChange = NA,
                                    lfcSE = NA,
                                    stat = NA,
                                    pvalue = NA,
                                    padj = NA,
                                    comparison = NA,
                                    stringsAsFactors = FALSE
                                )
                                results_list[[gene]] <- placeholder_row
                            }
                        }

                        # Combine all results into a single summary table
                        summary_table <- do.call(rbind, results_list)

                        # Ensure numeric consistency for columns
                        numeric_columns <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
                        summary_table[numeric_columns] <- lapply(summary_table[numeric_columns], as.numeric)

                        # Save the summary table
                        output_file <- get_genes_of_interest_path("gene_search_summary.csv")
                        write.csv(summary_table, file = output_file, row.names = FALSE)
                        cat("Summary table saved to:", output_file, "\n")
                    } # End of Option 1: Check Genes of Interest
                    if (user_choice == 2) {
                        # **Option 2: Generate a Heatmap**
                        # Step 1: Ask user whether to generate default or customized heatmap
                        cat("Do you want to generate default or customized heatmaps? ([1] Default / [2] Customized):\n")
                        heatmap_mode_raw <- handle_quit(readline())
                        if (heatmap_mode_raw == "back") next
                        heatmap_mode <- as.numeric(heatmap_mode_raw)
                        if (!(heatmap_mode %in% c(1, 2))) {
                            cat("Invalid input. Please enter 1 or 2. Returning to menu.\n")
                            next
                        }

                        # Step 2: Ask user to select a grouping column
                        cat("Available columns in coldata:\n")
                        for (i in seq_along(colnames(coldata))) {
                            cat(sprintf("[%d] %s\n", i, colnames(coldata)[i]))
                        }
                        repeat {
                            cat("Please select the column to group samples by (enter column name or index, or 'back' to cancel):\n")
                            group_input <- handle_quit(readline())
                            if (group_input == "back") break
                            if (suppressWarnings(!is.na(as.numeric(group_input)))) {
                                col_index <- as.numeric(group_input)
                                if (col_index >= 1 && col_index <= ncol(coldata)) {
                                    group_col <- colnames(coldata)[col_index]
                                    break
                                }
                            } else if (group_input %in% colnames(coldata)) {
                                group_col <- group_input
                                break
                            }
                            cat("Invalid input. Please enter a valid column name or index.\n")
                        }
                        if (group_input == "back") next
                        cat(sprintf("Grouping samples by column: '%s'\n", group_col))

                        # Step 3: Scaling for heatmap
                        cat("Heatmap scaling: [1] Log2 / [2] z-score:\n")
                        scale_raw <- handle_quit(readline())
                        if (scale_raw == "back") next
                        scale_choice <- as.numeric(scale_raw)
                        if (!scale_choice %in% c(1, 2)) {
                            cat("Invalid input. Please enter 1 or 2. Returning to menu.\n")
                            next
                        }
                        heatmap_scale <- if (scale_choice == 1L) "log2" else "z_score"

                        # Step 4: Color palette for heatmap
                        cat("Color palette: [1] viridis / [2] mako / [3] inferno:\n")
                        pal_raw <- handle_quit(readline())
                        if (pal_raw == "back") next
                        pal_choice <- as.numeric(pal_raw)
                        if (!pal_choice %in% c(1, 2, 3)) {
                            cat("Invalid input. Please enter 1, 2, or 3. Returning to menu.\n")
                            next
                        }
                        heatmap_palette <- switch(pal_choice, viridis::viridis(101), viridis::mako(101), viridis::inferno(101))

                        # Step 5: Ask whether to split into upregulated and downregulated genes
                        cat("Do you want to split the results into upregulated and downregulated genes? (y/n):\n")
                        split_choice <- tolower(handle_quit(readline()))
                        split_genes <- split_choice == "y"

                        # Step 6: Set padj threshold for significance filtering
                        repeat {
                            cat("Enter the maximum padj value for filtering significant genes (or type 'back' to cancel):\n")
                            padj_threshold_raw <- handle_quit(readline())
                            if (padj_threshold_raw == "back") break
                            padj_threshold <- as.numeric(padj_threshold_raw)
                            if (!is.na(padj_threshold) && padj_threshold > 0 && padj_threshold <= 1) break
                            cat("Invalid input. Please enter a value between 0 and 1.\n")
                        }
                        if (exists("padj_threshold_raw") && padj_threshold_raw == "back") next

                        # Step 5: Filter DEGs based on padj threshold
                        significant_genes <- degs_data[degs_data$padj < padj_threshold, , drop = FALSE]
                        significant_genes <- significant_genes[complete.cases(significant_genes), ]

                        if (nrow(significant_genes) == 0) {
                            cat("No significant genes found for the selected padj threshold.\n")
                            next
                        }

                        # Step 6: Select comparisons for the heatmap
                        comparisons <- unique(as.character(na.omit(degs_data$comparison)))
                        repeat {
                            cat("Available comparisons:\n")
                            comparison_list <- list()
                            for (i in seq_along(comparisons)) {
                                cat(sprintf("[%d] %s\n", i, comparisons[i]))
                                comparison_list[[as.character(i)]] <- comparisons[i]
                                if (split_genes) {
                                    cat(sprintf("[%d_up] %s_up\n", i, comparisons[i]))
                                    cat(sprintf("[%d_down] %s_down\n", i, comparisons[i]))
                                }
                            }
                            if (heatmap_mode == 1) {
                                # Default heatmap with all comparisons
                                comparison_groups <- unique(significant_genes$comparison)

                                for (group in comparison_groups) {
                                    group_genes <- significant_genes[significant_genes$comparison == group, "Gene"]
                                    if (length(group_genes) == 0) next

                                    if (split_genes){
                                        upregulated_heatmaps <- list()
                                        downregulated_heatmaps <- list()
                                        for (gr in comparison_groups) {
                                            logged_gr <- get_logged_for_comparison(gr, counts_data, coldata, group_col)
                                            logged_gr <- logged_gr[complete.cases(logged_gr), , drop = FALSE]
                                            group_genes_gr <- significant_genes[significant_genes$comparison == gr, "Gene"]
                                            group_heatmap_data <- logged_gr[rownames(logged_gr) %in% group_genes_gr, , drop = FALSE]
                                            group_heatmap_data <- as.matrix(group_heatmap_data)
                                            if (nrow(group_heatmap_data) == 0) next
                                            if (any(!is.finite(group_heatmap_data))) {
                                                group_heatmap_data <- apply(group_heatmap_data, 2, function(col) {
                                                    col[!is.finite(col)] <- median(col, na.rm = TRUE); return(col)
                                                })
                                            }
                                            upregulated_genes <- significant_genes[significant_genes$comparison == gr & significant_genes$log2FoldChange > 0, "Gene"]
                                            downregulated_genes <- significant_genes[significant_genes$comparison == gr & significant_genes$log2FoldChange < 0, "Gene"]
                                            upregulated_overlap <- intersect(upregulated_genes, rownames(group_heatmap_data))
                                            downregulated_overlap <- intersect(downregulated_genes, rownames(group_heatmap_data))
                                            if (length(upregulated_overlap) > 0) {
                                                upregulated_data <- group_heatmap_data[upregulated_overlap, , drop = FALSE]
                                                if (nrow(upregulated_data) > 0) {
                                                    upregulated_heatmaps[[length(upregulated_heatmaps) + 1]] <- list(heatmap_data = upregulated_data, title = gr)
                                                    upregulated_data_with_names <- data.frame(Gene = rownames(upregulated_data), upregulated_data, check.names = FALSE)
                                                    heatmap_obj_up <- if (heatmap_scale == "z_score") heatmap_z(upregulated_data, plot_incl = seq_len(ncol(upregulated_data)), kmeans_split = NULL, palette = heatmap_palette) else heatmap_v(upregulated_data, plot_incl = seq_len(ncol(upregulated_data)), kmeans_split = NULL, palette = heatmap_palette)
                                                    upregulated_csv_file <- get_heatmaps_path(heatmap_scale, "up", sprintf("heatmap_upregulated_%s.csv", gsub(" ", "_", gr)))
                                                    write.csv(upregulated_data_with_names, file = upregulated_csv_file, row.names = FALSE)
                                                    upregulated_png_file <- get_heatmaps_path(heatmap_scale, "up", sprintf("heatmap_upregulated_%s.png", gsub(" ", "_", gr)))
                                                    if (!is.null(heatmap_obj_up)) { png(upregulated_png_file, width = 900, height = 500); ComplexHeatmap::draw(heatmap_obj_up, heatmap_legend_side = "right"); dev.off() }
                                                }
                                            }
                                            if (length(downregulated_overlap) > 0) {
                                                downregulated_data <- group_heatmap_data[downregulated_overlap, , drop = FALSE]
                                                if (nrow(downregulated_data) > 0) {
                                                    downregulated_heatmaps[[length(downregulated_heatmaps) + 1]] <- list(heatmap_data = downregulated_data, title = gr)
                                                    downregulated_data_with_names <- data.frame(Gene = rownames(downregulated_data), downregulated_data, check.names = FALSE)
                                                    heatmap_obj_down <- if (heatmap_scale == "z_score") heatmap_z(downregulated_data, plot_incl = seq_len(ncol(downregulated_data)), kmeans_split = NULL, palette = heatmap_palette) else heatmap_v(downregulated_data, plot_incl = seq_len(ncol(downregulated_data)), kmeans_split = NULL, palette = heatmap_palette)
                                                    downregulated_csv_file <- get_heatmaps_path(heatmap_scale, "down", sprintf("heatmap_downregulated_%s.csv", gsub(" ", "_", gr)))
                                                    write.csv(downregulated_data_with_names, file = downregulated_csv_file, row.names = FALSE)
                                                    downregulated_png_file <- get_heatmaps_path(heatmap_scale, "down", sprintf("heatmap_downregulated_%s.png", gsub(" ", "_", gr)))
                                                    if (!is.null(heatmap_obj_down)) { png(downregulated_png_file, width = 900, height = 500); ComplexHeatmap::draw(heatmap_obj_down, heatmap_legend_side = "right"); dev.off() }
                                                }
                                            }
                                        }
                                        if (length(upregulated_heatmaps) > 0) {
                                            combined_upregulated_data <- do.call(rbind, lapply(upregulated_heatmaps, function(hm) hm$heatmap_data))
                                            row_ann_up <- unlist(lapply(upregulated_heatmaps, function(hm) rep(hm$title, nrow(hm$heatmap_data))))
                                            row_split_upregulated <- factor(row_ann_up, levels = unique(row_ann_up))
                                            heatmap_obj_upregulated <- if (heatmap_scale == "z_score") heatmap_z(combined_upregulated_data, plot_incl = seq_len(ncol(combined_upregulated_data)), row_split = row_split_upregulated, palette = heatmap_palette, sum_style = TRUE) else heatmap_v(combined_upregulated_data, plot_incl = seq_len(ncol(combined_upregulated_data)), row_split = row_split_upregulated, palette = heatmap_palette, sum_style = TRUE)
                                            png(get_heatmaps_path(heatmap_scale, "up", "sum_upregulated_heatmap.png"), width = 1600, height = 2400)
                                            if (!is.null(heatmap_obj_upregulated)) ComplexHeatmap::draw(heatmap_obj_upregulated, heatmap_legend_side = "right")
                                            dev.off()
                                            write.csv(data.frame(Gene = rownames(combined_upregulated_data), combined_upregulated_data, check.names = FALSE), get_heatmaps_path(heatmap_scale, "up", "sum_upregulated_heatmap.csv"), row.names = FALSE)
                                        }
                                        if (length(downregulated_heatmaps) > 0) {
                                            combined_downregulated_data <- do.call(rbind, lapply(downregulated_heatmaps, function(hm) hm$heatmap_data))
                                            row_ann_down <- unlist(lapply(downregulated_heatmaps, function(hm) rep(hm$title, nrow(hm$heatmap_data))))
                                            row_split_downregulated <- factor(row_ann_down, levels = unique(row_ann_down))
                                            heatmap_obj_downregulated <- if (heatmap_scale == "z_score") heatmap_z(combined_downregulated_data, plot_incl = seq_len(ncol(combined_downregulated_data)), row_split = row_split_downregulated, palette = heatmap_palette, sum_style = TRUE) else heatmap_v(combined_downregulated_data, plot_incl = seq_len(ncol(combined_downregulated_data)), row_split = row_split_downregulated, palette = heatmap_palette, sum_style = TRUE)
                                            png(get_heatmaps_path(heatmap_scale, "down", "sum_downregulated_heatmap.png"), width = 1600, height = 2400)
                                            ComplexHeatmap::draw(heatmap_obj_downregulated, heatmap_legend_side = "right")
                                            dev.off()
                                            write.csv(data.frame(Gene = rownames(combined_downregulated_data), combined_downregulated_data, check.names = FALSE), get_heatmaps_path(heatmap_scale, "down", "sum_downregulated_heatmap.csv"), row.names = FALSE)
                                        }
                                    } else if (!split_genes) {
                                        individual_heatmaps <- list()
                                        for (gr in comparison_groups) {
                                            logged_gr <- get_logged_for_comparison(gr, counts_data, coldata, group_col)
                                            logged_gr <- logged_gr[complete.cases(logged_gr), , drop = FALSE]
                                            group_genes <- significant_genes[significant_genes$comparison == gr, "Gene"]
                                            group_heatmap_data <- logged_gr[rownames(logged_gr) %in% group_genes, , drop = FALSE]
                                            if (nrow(group_heatmap_data) == 0) next
                                            group_heatmap_data <- as.matrix(group_heatmap_data)
                                            if (any(!is.finite(group_heatmap_data))) {
                                                group_heatmap_data <- apply(group_heatmap_data, 2, function(col) {
                                                    col[!is.finite(col)] <- median(col, na.rm = TRUE); return(col)
                                                })
                                            }
                                            group_heatmap_data_with_names <- data.frame(Gene = rownames(group_heatmap_data), group_heatmap_data, check.names = FALSE)
                                            heatmap_obj_combined <- if (heatmap_scale == "z_score") heatmap_z(group_heatmap_data, plot_incl = seq_len(ncol(group_heatmap_data)), kmeans_split = NULL, palette = heatmap_palette) else heatmap_v(group_heatmap_data, plot_incl = seq_len(ncol(group_heatmap_data)), kmeans_split = NULL, palette = heatmap_palette)
                                            individual_heatmaps[[length(individual_heatmaps) + 1]] <- list(heatmap_data = group_heatmap_data, title = gr)
                                            combined_csv_file <- get_heatmaps_path(heatmap_scale, "combined", sprintf("heatmap_combined_%s.csv", gsub(" ", "_", gr)))
                                            write.csv(group_heatmap_data_with_names, file = combined_csv_file, row.names = FALSE)
                                            combined_png_file <- get_heatmaps_path(heatmap_scale, "combined", sprintf("heatmap_combined_%s.png", gsub(" ", "_", gr)))
                                            if (!is.null(heatmap_obj_combined)) { png(combined_png_file, width = 900, height = 500); ComplexHeatmap::draw(heatmap_obj_combined, heatmap_legend_side = "right"); dev.off() }
                                        }
                                        if (length(individual_heatmaps) > 0) {
                                            combined_data <- do.call(rbind, lapply(individual_heatmaps, function(hm) hm$heatmap_data))
                                            row_annotations <- unlist(lapply(individual_heatmaps, function(hm) rep(hm$title, nrow(hm$heatmap_data))))
                                            row_split <- factor(row_annotations, levels = unique(row_annotations))
                                            heatmap_obj_combined <- if (heatmap_scale == "z_score") heatmap_z(combined_data, plot_incl = seq_len(ncol(combined_data)), row_split = row_split, palette = heatmap_palette, sum_style = TRUE) else heatmap_v(combined_data, plot_incl = seq_len(ncol(combined_data)), row_split = row_split, palette = heatmap_palette, sum_style = TRUE)
                                            png(get_heatmaps_path(heatmap_scale, "combined", "sum_heatmap.png"), width = 1600, height = 2400)
                                            if (!is.null(heatmap_obj_combined)) ComplexHeatmap::draw(heatmap_obj_combined, heatmap_legend_side = "right")
                                            dev.off()
                                            write.csv(data.frame(Gene = rownames(combined_data), combined_data, check.names = FALSE), get_heatmaps_path(heatmap_scale, "combined", "sum_heatmap.csv"), row.names = FALSE)
                                        }
                                    }
                                }
                                break
                            } else if (heatmap_mode == 2) {
                                # User selects comparisons
                                cat("Select indices (e.g., 1,2_up,3_down) for comparisons to intersect (or type 'quit' to exit):\n")
                                selected_indices <- unlist(strsplit(handle_quit(readline()), ","))

                                # Step 6: Filter selected results
                                selected_results <- list()
                                for (index in selected_indices) {
                                    if (grepl("_up", index)) {
                                        comp <- gsub("_up", "", index)
                                        selected_results[[index]] <- degs_data$Gene[degs_data$comparison == comparison_list[[comp]] & 
                                                                                    degs_data$log2FoldChange > 0 & 
                                                                                    degs_data$padj < padj_threshold]
                                    } else if (grepl("_down", index)) {
                                        comp <- gsub("_down", "", index)
                                        selected_results[[index]] <- degs_data$Gene[degs_data$comparison == comparison_list[[comp]] & 
                                                                                    degs_data$log2FoldChange < 0 & 
                                                                                    degs_data$padj < padj_threshold]
                                    } else {
                                        selected_results[[index]] <- degs_data$Gene[degs_data$comparison == comparison_list[[index]] & 
                                                                                    degs_data$padj < padj_threshold]
                                    }
                                    selected_results[[index]] <- unique(na.omit(selected_results[[index]]))
                                    cat(sprintf("Number of genes for %s: %d\n", index, length(selected_results[[index]])))
                                }

                                # Step 7: Intersect genes
                                intersected_genes <- Reduce(intersect, selected_results)
                                cat(sprintf("Number of genes in intersection: %d\n", length(intersected_genes)))
                                if (length(intersected_genes) == 0) {
                                    cat("No genes found in the intersection. Please try other comparisons.\n")
                                    next
                                }
                                # Map indices to actual comparison names
                                comparison_names <- sapply(selected_indices, function(index) {
                                    if (grepl("_up", index)) {
                                        comp <- gsub("_up", "", index)
                                        paste0(comparison_list[[comp]], "_up")
                                    } else if (grepl("_down", index)) {
                                        comp <- gsub("_down", "", index)
                                        paste0(comparison_list[[comp]], "_down")
                                    } else {
                                        comparison_list[[index]]
                                    }
                                })

                                group_genes <- significant_genes[significant_genes$Gene %in% intersected_genes, "Gene"]
                                comparison_groups <- comparison_names
                                individual_heatmaps <- list()
                                file_safe_name <- gsub("[^A-Za-z0-9_]", "_", paste(comparison_names, collapse = "_intersect_"))

                                for (group in comparison_groups) {
                                    comp_base <- gsub("_down$", "", gsub("_up$", "", group))
                                    logged_gr <- get_logged_for_comparison(comp_base, counts_data, coldata, group_col)
                                    logged_gr <- logged_gr[complete.cases(logged_gr), , drop = FALSE]
                                    group_heatmap_data <- logged_gr[rownames(logged_gr) %in% group_genes, , drop = FALSE]
                                    if (nrow(group_heatmap_data) == 0) next
                                    group_heatmap_data <- as.matrix(group_heatmap_data)
                                    if (any(!is.finite(group_heatmap_data))) {
                                        group_heatmap_data <- apply(group_heatmap_data, 2, function(col) {
                                            col[!is.finite(col)] <- median(col, na.rm = TRUE); return(col)
                                        })
                                    }
                                    group_heatmap_data_with_names <- data.frame(Gene = rownames(group_heatmap_data), group_heatmap_data, check.names = FALSE)
                                    heatmap_obj_combined <- if (heatmap_scale == "z_score") heatmap_z(group_heatmap_data, plot_incl = seq_len(ncol(group_heatmap_data)), kmeans_split = NULL, palette = heatmap_palette) else heatmap_v(group_heatmap_data, plot_incl = seq_len(ncol(group_heatmap_data)), kmeans_split = NULL, palette = heatmap_palette)
                                    individual_heatmaps[[length(individual_heatmaps) + 1]] <- list(heatmap_data = group_heatmap_data, title = group)
                                    combined_csv_file <- get_heatmaps_path(heatmap_scale, "combined", sprintf("heatmap_%s.csv", file_safe_name))
                                    write.csv(group_heatmap_data_with_names, file = combined_csv_file, row.names = FALSE)
                                    combined_png_file <- get_heatmaps_path(heatmap_scale, "combined", sprintf("heatmap_%s.png", file_safe_name))
                                    if (!is.null(heatmap_obj_combined)) { png(combined_png_file, width = 900, height = 500); ComplexHeatmap::draw(heatmap_obj_combined, heatmap_legend_side = "right"); dev.off() }
                                }

                                cat("Do you want to create another heatmap? (y/n):\n")
                                repeat_choice <- tolower(handle_quit(readline()))
                                if (repeat_choice == "n"){
                                    if (length(individual_heatmaps) > 1) {
                                        combined_data <- do.call(rbind, lapply(individual_heatmaps, function(hm) hm$heatmap_data))
                                        row_annotations <- unlist(lapply(individual_heatmaps, function(hm) rep(hm$title, nrow(hm$heatmap_data))))
                                        row_split <- factor(row_annotations, levels = unique(row_annotations))
                                        heatmap_obj_combined <- if (heatmap_scale == "z_score") heatmap_z(combined_data, plot_incl = seq_len(ncol(combined_data)), row_split = row_split, palette = heatmap_palette, sum_style = TRUE) else heatmap_v(combined_data, plot_incl = seq_len(ncol(combined_data)), row_split = row_split, palette = heatmap_palette, sum_style = TRUE)
                                        png(get_heatmaps_path(heatmap_scale, "combined", "sum_intersections_heatmap.png"), width = 1600, height = 2400)
                                        if (!is.null(heatmap_obj_combined)) ComplexHeatmap::draw(heatmap_obj_combined, heatmap_legend_side = "right")
                                        dev.off()
                                        write.csv(data.frame(Gene = rownames(combined_data), combined_data, check.names = FALSE), get_heatmaps_path(heatmap_scale, "combined", "sum_intersections_heatmap.csv"), row.names = FALSE)
                                    }
                                    break
                                }
                            } else {
                                cat("Invalid choice. Please enter 1 or 2.\n")
                            }
                        } # repeat for heatmap generation
                    } # End of heatmap generation
                    if (user_choice == 3) {
                        # **Option 3: Create Venn and UpSet diagrams**
                        # (A) Ask user about splitting into up/down
                        cat("Do you want to split results into upregulated and downregulated genes? (y/n):\n")
                        split_choice <- tolower(handle_quit(readline()))
                        if (! split_choice %in% c("y","n")) {
                            cat("Please enter 'y' or 'n'.\n")
                            next
                        }

                        # (B) Ask for padj threshold
                        cat("Enter the padj threshold for significant genes (or type 'back' to cancel):\n")
                        padj_threshold_raw <- handle_quit(readline())
                        if (padj_threshold_raw == "back") next
                        padj_threshold <- as.numeric(padj_threshold_raw)
                        if (is.na(padj_threshold) || padj_threshold <= 0 || padj_threshold > 1) {
                            cat("Invalid padj threshold. Please enter a value between 0 and 1.\n")
                            next
                        }

                        # (C) Filter DEGs
                        degs_data <- degs_data[!is.na(degs_data$padj) & degs_data$padj < padj_threshold, ]
                        if (nrow(degs_data) == 0) {
                            cat("No significant genes found for the given padj threshold. Returning to menu.\n")
                            next
                        }

                        # (D) Split into up/down if requested
                        if (split_choice == "y") {
                            degs_data_up   <- subset(degs_data, log2FoldChange > 0)
                            degs_data_down <- subset(degs_data, log2FoldChange < 0)
                        } else {
                            degs_data_up   <- NULL
                            degs_data_down <- NULL
                        }

                        # (E) Build indexed list of comparisons (+ “_up” and “_down” if split)
                        comparisons   <- unique(degs_data$comparison)
                        options_list  <- list()
                        for (i in seq_along(comparisons)) {
                            cmp <- comparisons[i]
                            options_list[[paste0(i)]] <- cmp
                            if (split_choice == "y") {
                                options_list[[paste0(i, "_(up)")]]   <- paste0(cmp, "_(up)")
                                options_list[[paste0(i, "_(down)")]] <- paste0(cmp, "_(down)")
                            }
                        }
                        # Ask the user if they want default or custom Venn diagrams
                        cat("Do you want to create default Venn diagrams or custom? ([1] default/[2] custom), or type 'back' to cancel:\n")
                        venn_choice_raw <- handle_quit(readline())
                        if (venn_choice_raw == "back") next
                        venn_choice <- as.numeric(venn_choice_raw)
                        if (venn_choice == 1) {
                            repeat {
                                selected_up   <- list()
                                selected_down <- list()
                                for (cmp in comparisons) {
                                    up_vec   <- subset(degs_data,   comparison == cmp & log2FoldChange > 0)$Gene
                                    down_vec <- subset(degs_data,   comparison == cmp & log2FoldChange < 0)$Gene
                                    if (length(up_vec)   > 0) selected_up[[paste0(cmp, "_up")]]   <- unique(up_vec)
                                    if (length(down_vec) > 0) selected_down[[paste0(cmp, "_down")]] <- unique(down_vec)
                                }
                                selected_up   <- selected_up  [sapply(selected_up,   length) > 0]
                                selected_down <- selected_down[sapply(selected_down, length) > 0]

                                if (split_choice == "y") {
                                    # Launch the parallel gadget if at least one side has ≥2 sets
                                    if (length(selected_up) >= 2 || length(selected_down) >= 2) {
                                        export_venn_upset_split(selected_up, selected_down)
                                        cat("Run another interactive UpSet? (y/n):\n")
                                        again <- tolower(handle_quit(readline()))
                                        if (again != "y") break

                                    } else {
                                        cat("Need at least two UP or two DOWN sets to launch. Returning to menu.\n")
                                    }
                                } else {
                                    # No splitting: just run one UpSet on the combined list
                                    selected_all <- list()
                                    for (cmp in comparisons) {
                                        all_vec <- subset(degs_data, comparison == cmp)$Gene
                                        if (length(all_vec) > 0) selected_all[[cmp]] <- unique(all_vec)
                                    }
                                    selected_all <- selected_all[sapply(selected_all, length) > 0]

                                    if (length(selected_all) >= 2) {
                                        cat("Launching interactive UpSet (no split)…\n")
                                        export_venn_upset(selected_all)
                                        # After user closes the Shiny gadget, ask if they want to run again:
                                        cat("Run another default interactive UpSet? (y/n):\n")
                                        again <- tolower(handle_quit(readline()))
                                        if (again != "y") break
                                    } else {
                                        cat("Need at least two sets to launch interactive UpSet. Returning to menu.\n")
                                    }
                                }
                            } # repeat for default Venn diagrams
                        } # Default Venn diagrams
                        if (venn_choice == 2) {
                            repeat{
                                cat("Available comparisons:\n")
                                for (opt in names(options_list)) {
                                    if (split_choice == "n" && grepl("_up|_down", opt)) next
                                    cat(sprintf("[%s] %s\n", opt, options_list[[opt]]))
                                }

                                cat("Select the indices (e.g., 1,2,3) for comparisons to include in the Venn diagram (or type 'back' to cancel):\n")
                                selection <- handle_quit(readline())
                                if (selection == "back") break

                                selected_indices <- unlist(strsplit(selection, ","))
                                selected_genes   <- list()

                                for (idx in selected_indices) {
                                    comparison_name <- options_list[[idx]]
                                    if (grepl("_up$", idx)) {
                                        comp <- gsub("_(up)$", "", comparison_name)
                                        selected_genes[[comparison_name]] <- unique(degs_data_up$Gene[degs_data_up$comparison == comp])
                                    } else if (grepl("_down$", idx)) {
                                        comp <- gsub("_(down)$", "", comparison_name)
                                        selected_genes[[comparison_name]] <- unique(degs_data_down$Gene[degs_data_down$comparison == comp])
                                    } else {
                                        selected_genes[[comparison_name]] <- unique(degs_data$Gene[degs_data$comparison == comparison_name])
                                    }
                                }
                                # Drop any empty sets
                                selected_genes <- selected_genes[sapply(selected_genes, length) > 0]
                                if (length(selected_genes) < 2) {
                                    cat("At least two non‐empty sets are required. Try again.\n")
                                    next
                                }

                                # Instead of building a static Venn, launch the interactive UpSet+Venn gadget for these custom sets:
                                cat("Launching interactive UpSet+Venn for your custom selection…\n")
                                export_venn_upset(selected_genes)

                                # After the user closes the Shiny gadget, ask if they want to generate another custom Venn
                                cat("Do you want to generate another custom Venn? (y/n):\n")
                                again <- tolower(handle_quit(readline()))
                                if (again != "y") break
                            } # repeat for custom Venn diagrams 
                        } # Custom Venn diagrams
                    } # Venn Diagrams

                    if (user_choice == 4) {
                        # **Option 4: Create Trending Plots**
                        # Step 1: Ask user to select a grouping column
                        cat("Available columns in coldata:\n")
                        for (i in seq_along(colnames(coldata))) {
                            cat(sprintf("[%d] %s\n", i, colnames(coldata)[i]))
                        }

                        repeat {
                            cat("Please select the column to group samples by (enter column name or index, or 'back' to cancel):\n")
                            group_input <- handle_quit(readline())
                            if (group_input == "back") break
                            if (suppressWarnings(!is.na(as.numeric(group_input)))) {
                                col_index <- as.numeric(group_input)
                                if (col_index >= 1 && col_index <= ncol(coldata)) {
                                    group_col <- colnames(coldata)[col_index]
                                    break
                                }
                            } else if (group_input %in% colnames(coldata)) {
                                group_col <- group_input
                                break
                            }
                            cat("Invalid input. Please enter a valid column name or index.\n")
                        }
                        if (group_input == "back") next

                        cat(sprintf("Grouping samples by column: '%s'\n", group_col))

                        # Step 3: Ask whether to split into upregulated and downregulated genes
                        cat("Do you want to split the results into upregulated and downregulated genes? (y/n):\n")
                        split_choice <- tolower(handle_quit(readline()))
                        split_genes <- split_choice == "y"

                        # Step 4: Set padj threshold
                        cat("Enter the padj threshold for significant genes (or type 'back' to cancel):\n")
                        padj_threshold_raw <- handle_quit(readline())
                        if (padj_threshold_raw == "back") next
                        padj_threshold <- as.numeric(padj_threshold_raw)
                        if (is.na(padj_threshold) || padj_threshold <= 0 || padj_threshold > 1) {
                            cat("Invalid padj threshold. Returning to menu.\n")
                            next
                        }
                        # Step 5: Check for default or custom trending plots
                        cat("Do you want to create default trending plots or custom? ([1] default/[2] custom), or type 'back' to cancel:\n")
                        plot_choice_raw <- handle_quit(readline())
                        if (plot_choice_raw == "back") next
                        plot_choice <- as.numeric(plot_choice_raw)

                        # Step 5: Prepare comparisons and splits
                        comparisons <- unique(na.omit(degs_data$comparison))
                        repeat {
                            cat("Available comparisons:\n")
                            comparison_list <- list()
                            for (i in seq_along(comparisons)) {
                                cat(sprintf("[%d] %s\n", i, comparisons[i]))
                                comparison_list[[as.character(i)]] <- comparisons[i]
                                if (split_genes) {
                                    cat(sprintf("[%d_up] %s_up\n", i, comparisons[i]))
                                    cat(sprintf("[%d_down] %s_down\n", i, comparisons[i]))
                                }
                            }
                            if (plot_choice == 1) {
                                generate_default_trending_plots(
                                    counts_data = counts_data,
                                    coldata = coldata,
                                    significant_genes = degs_data,
                                    group_col = group_col,
                                    padj_threshold = padj_threshold
                                )
                                cat("\nDefault trending plots completed. Returning to the main menu...\n")
                                break
                            } else if (plot_choice == 2) {
                                # User selects comparisons
                                cat("Select indices (e.g., 1,2_up,3_down) for comparisons to intersect (or type 'quit' to exit):\n")
                                selected_indices <- unlist(strsplit(handle_quit(readline()), ","))

                                # Step 6: Filter selected results
                                selected_results <- list()
                                for (index in selected_indices) {

                                    if (grepl("_up$", index)) {
                                        comp_index <- gsub("_up$", "", index)
                                        comp_name <- comparison_list[[comp_index]]
                                        genes <- degs_data$Gene[
                                            degs_data$comparison == comp_name &
                                            degs_data$log2FoldChange > 0 &
                                            degs_data$padj < padj_threshold
                                        ]
                                    } else if (grepl("_down$", index)) {
                                        comp_index <- gsub("_down$", "", index)
                                        comp_name <- comparison_list[[comp_index]]
                                        genes <- degs_data$Gene[
                                            degs_data$comparison == comp_name &
                                            degs_data$log2FoldChange < 0 &
                                            degs_data$padj < padj_threshold
                                        ]
                                    } else {
                                        comp_name <- comparison_list[[index]]
                                        genes <- degs_data$Gene[
                                            degs_data$comparison == comp_name &
                                            degs_data$padj < padj_threshold
                                        ]
                                    }

                                selected_results[[index]] <- unique(na.omit(genes))
                                cat(sprintf("Genes for %s: %d\n", index, length(selected_results[[index]])))
                                }

                                # Step 7: Intersect genes
                                intersected_genes <- Reduce(intersect, selected_results)
                                cat(sprintf("Number of genes in intersection: %d\n", length(intersected_genes)))
                                if (length(intersected_genes) == 0) {
                                    cat("No genes found in the intersection. Please try other comparisons.\n")
                                    next
                                }

                                base_selected <- gsub("_up$|_down$", "", selected_indices)        # numeric indices w/o direction
                                other_idx     <- setdiff(seq_along(comparison_list), as.numeric(base_selected))
                                other_genes   <- unique(na.omit(
                                degs_data$Gene[
                                    degs_data$comparison %in% comparison_list[other_idx] &
                                    degs_data$padj < padj_threshold
                                ]))
                                # other_genes <- unique(na.omit(unlist(
                                #     lapply(selected_indices, function(idx) {
                                #         if (grepl("_up$", idx)) {                     # keep direction
                                #             degs_data$Gene[
                                #                 degs_data$comparison %in% comparison_list[other_idx] &
                                #                 degs_data$log2FoldChange > 0 &            # same sign only
                                #                 degs_data$padj < padj_threshold
                                #             ]
                                #         } else if (grepl("_down$", idx)) {
                                #             degs_data$Gene[
                                #                 degs_data$comparison %in% comparison_list[other_idx] &
                                #                 degs_data$log2FoldChange < 0 &            # same sign only
                                #                 degs_data$padj < padj_threshold
                                #             ]
                                #         } else {                                      # no-direction set → keep any sign
                                #             degs_data$Gene[
                                #                 degs_data$comparison %in% comparison_list[other_idx] &
                                #                 degs_data$padj < padj_threshold
                                #             ]
                                #         }
                                #     })
                                # )))

                                intersected_genes <- setdiff(intersected_genes, other_genes)
                                cat(sprintf("Unique genes after exclusion: %d\n", length(intersected_genes)))
                                if (!length(intersected_genes)) {
                                    cat("No genes left after uniqueness filter. Try different comparisons.\n")
                                    next
                                }
                                # Step 8: Prepare the plot (trending uses counts_data only; x-axis = conditions)
                                trending_data <- counts_data[rownames(counts_data) %in% intersected_genes, , drop = FALSE]
                                group_vec      <- coldata[[group_col]]
                                group_levels   <- unique(group_vec)                 # first appearance order
                                coldata[[group_col]] <- factor(group_vec, levels = group_levels, ordered = TRUE)

                                sample_conditions <- split(rownames(coldata), coldata[[group_col]])

                                comparison_names <- sapply(selected_indices, function(idx) {
                                    if (grepl("_up$", idx))   paste0(comparison_list[[gsub("_up$",   "", idx)]], "_up")
                                    else if (grepl("_down$", idx)) paste0(comparison_list[[gsub("_down$", "", idx)]], "_down")
                                    else comparison_list[[idx]]
                                })

                                file_safe_name <- gsub("[^A-Za-z0-9_]", "_", paste(comparison_names, collapse = "_intersect_"))
                                plot_name_png <- get_trending_path(paste0("trending_plot_", file_safe_name, ".png"))
                                plot_name_svg <- get_trending_path(paste0("trending_plot_", file_safe_name, ".svg"))
                                output_file <- get_trending_path(paste0("trending_data_", file_safe_name, ".csv"))

                                tryCatch({
                                    temp_degs_data <- degs_data
                                    temp_degs_data$comparison <- gsub("_down$", "", gsub("_up$", "", temp_degs_data$comparison))
                                    display_names <- gsub("_(up|down)$", " (\\1)", comparison_names)
                                    title <- paste(display_names, collapse = " INT ")

                                    avg_expr <- t(apply(trending_data, 1, function(x)
                                        tapply(x, coldata[[group_col]], mean, na.rm = TRUE)))
                                    condition_means <- colMeans(avg_expr, na.rm = TRUE)
                                    condition_names_plot <- colnames(avg_expr)
                                    comps_for_brackets <- unique(unlist(comparison_list[base_selected]))

                                    p <- plot_trending_ggplot(
                                        as.numeric(condition_means),
                                        condition_names_plot,
                                        title,
                                        degs_data = temp_degs_data,
                                        comparisons = comps_for_brackets,
                                        show_title = FALSE
                                    )
                                    if (!is.null(p)) {
                                        ggplot2::ggsave(plot_name_png, plot = p, width = gene_expr_width_in, height = gene_expr_height_in, dpi = gene_expr_dpi, device = "png", units = "in")
                                        ggplot2::ggsave(plot_name_svg, plot = p, width = gene_expr_width_in, height = gene_expr_height_in, dpi = gene_expr_dpi, device = "svg", units = "in")
                                    }

                                    export_data <- data.frame(
                                        Gene      = rep(rownames(avg_expr), times = ncol(avg_expr)),
                                        Condition = rep(colnames(avg_expr), each = nrow(avg_expr)),
                                        Expression = as.vector(avg_expr)
                                    )
                                    write.csv(export_data, output_file, row.names = FALSE)
                                    cat("Trending plot saved as:", plot_name_png, "and", plot_name_svg, "\n")
                                    cat("Trending data exported as:", output_file, "\n")
                                }, error = function(e) {
                                    cat(sprintf("Error generating trending plot for comparisons: %s. Error: %s\n",
                                                paste(selected_indices, collapse = ", "), e$message))
                                })

                                # Ask if the user wants another plot
                                cat("Do you want to create another trending plot? (y/n):\n")
                                repeat_choice <- tolower(handle_quit(readline()))
                                if (repeat_choice == "n") break
                            }
                            else {
                                cat("Invalid choice. Please enter 1 or 2.\n")
                            }
                        } # repeat for trending plots
                    } # Trending Plots
                } # repeat for manual analysis menu
            } # End of manual analysis menu
        } # Load DEGS tables and analyze

        if (menu_choice == 2) {
            # **Option 2: Process count files and perform differential expression analysis**
            repeat {
                # ask the user to choose if to include [1] to process the count files or directly[2] to perform the analysis menu
                cat("Please select an option by entering the corresponding number:\n")
                cat("[1] Process count files and perform differential expression analysis\n")
                cat("[2] Perform differential expression analysis without processing count files\n")
                cat("[3] Quit\n")
                cat("Please select an option by entering the corresponding number(or type 'back' to return to the main menu):\n")

                sub_menu_choice <- handle_quit(readline())

                if (sub_menu_choice == "back"){
                    cat("Returning to the main menu.\n")
                    break # Return to the main menu
                }

                sub_menu_choice <- as.numeric(sub_menu_choice)

                if (is.na(sub_menu_choice) || sub_menu_choice < 1 || sub_menu_choice > 3) {
                    cat("Invalid choice. Please enter a valid number between 1 and 3.\n")
                    next
                }

                if (sub_menu_choice == 3) {
                    cat("Exiting the program. Goodbye!\n")
                    do_quit()
                }

                if (sub_menu_choice == 1) {
                    # Ask the user if they have a count matrix or a folder of count files
                    cat("Do you have a:\n")
                    cat("[1] Count matrix\n")
                    cat("[2] Folder of count files\n")
                    repeat {
                        data_input_choice <- handle_quit(readline("Enter your choice (1 or 2): "))
                        if (data_input_choice == "back") break
                        data_input_choice <- as.numeric(data_input_choice)
                        if (data_input_choice %in% c(1, 2)) break
                        cat("Invalid choice. Please enter 1 or 2.\n")
                    }

                    if (data_input_choice == "back") {
                        cat("Returning to the previous step.\n")
                        next # Return to the previous step
                    }

                    if (data_input_choice == 1) {
                        # Initialize loop control variables
                        # navigate_back <- FALSE
                        # User has a count matrix
                        repeat {
                            cat("Please provide the path to the count matrix file (or type 'back' to return to the previous step): \n")
                            matrix_path <- handle_quit(readline())
                            if (matrix_path == "back") {
                                cat("Returning to the data input choice.\n")
                                # navigate_back <- TRUE
                                break # Return to the data input choice
                            }

                            if (file.exists(matrix_path)) {
                                tryCatch({
                                    cat("Reading the count matrix...\n")
                                    if (grepl("\\.xlsx$", matrix_path, ignore.case = TRUE)) {
                                        combined_data <- read_excel(matrix_path, col_names = TRUE)
                                        # Ensure consistency with CSV reading
                                        combined_data <- as.data.frame(combined_data, stringsAsFactors = FALSE)  # Convert to data.frame
                                        colnames(combined_data) <- trimws(colnames(combined_data))               # Trim whitespace from headers
                                        combined_data[] <- lapply(combined_data, function(x) {                   # Trim whitespace from data
                                            if (is.character(x)) trimws(x) else x
                                        })
                                        combined_data[combined_data == ""] <- NA                                # Handle empty cells as NA
                                    } else if (grepl("\\.(txt|tsv)$", matrix_path, ignore.case = TRUE)) {
                                        combined_data <- read.delim(matrix_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
                                    } else if (grepl("\\.csv$", matrix_path, ignore.case = TRUE)) {
                                        combined_data <- read.csv(matrix_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
                                    } else {
                                        cat("Error: Unsupported file format. Please provide a .csv, .tsv, .txt, or .xlsx file.\n")
                                        next # Return to the start of the loop to retry
                                    }
                                    cat("Count matrix loaded successfully.\n")
                                    # navigate_back <- FALSE
                                    # break # Exit the loop if the matrix is loaded successfully
                                }, error = function(e) {
                                    cat("Error reading the count matrix: ", e$message, "\n")
                                    next # Retry the input if there’s an error
                                })
                                break
                            } else {
                                cat("Invalid file path. Please try again.\n")
                            }
                        }
                        
                        # if (navigate_back) break # # Return to main menu if "back" was selected
                        # Ask if there is a sample info table
                        repeat {
                            cat("Do you have a sample info table? (y/n, or type 'back' to return to the matrix input, 'quit' to exit): \n")
                            info_table_choice <- handle_quit(readline())
                            # Handle "back" command to return to the matrix input step
                            if (info_table_choice == "back") {
                                cat("Returning to the matrix input step.\n")
                                # Reset variables to ensure clean navigation
                                combined_data <- NULL
                                final_conditions <- NULL
                                # navigate_back <- TRUE
                                break # Exit to return to the matrix input loop. TOFIX: breaks to the main menu
                            }
                            if (tolower(info_table_choice) == "y") {
                                repeat {
                                    cat("Please provide the path to the sample info table (or type 'quit' to exit): \n")
                                    info_table_path <- handle_quit(readline())
                                    # Handle "back" command to return to the previous step
                                    if (info_table_path == "back") {
                                        cat("Returning to the sample info table choice.\n")
                                        final_conditions <- NULL
                                        break # Exit to return to the sample info table choice loop
                                    }

                                    if (file.exists(info_table_path)) {
                                        tryCatch({
                                            # Function to read sample info table
                                            read_sample_info <- function(file_path) {
                                                if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
                                                    read.csv(file_path, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
                                                } else if (grepl("\\.(txt|tsv)$", file_path, ignore.case = TRUE)) {
                                                    read.delim(file_path, header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = "\t")
                                                } else if (grepl("\\.xlsx$", file_path, ignore.case = TRUE)) {
                                                    read_excel(file_path, col_names = TRUE)
                                                } else {
                                                    cat("Error: Unsupported file format. Please provide a .csv, .txt, .tsv, or .xlsx file.\n")
                                                    return(NULL)
                                                }
                                            }

                                            # Read the sample info table
                                            sample_info <- read_sample_info(info_table_path)
                                            if (is.null(sample_info)) next # Retry if unsupported format

                                            # Debugging: Print parsed sample info table
                                            cat("Sample info table loaded successfully:\n")
                                            print(sample_info)

                                            # Exit the loop on successful read
                                            break
                                        }, error = function(e) {
                                            cat("Error reading the sample info table: ", e$message, "\n")
                                            next # Retry on error
                                        })
                                    } else {
                                        cat("Invalid file path. Please try again.\n")
                                    }
                                }
                                # # Align sample info with count matrix
                                tryCatch({
                                    # Align sample info with count matrix
                                    coldata <- sample_info
                                    # Align with the count matrix columns except the first which is assumed to be geneID
                                    coldata <- coldata[match(colnames(combined_data)[-1], rownames(coldata)), , drop = FALSE]

                                    if (any(is.na(coldata))) {
                                        cat("Error: Mismatch between sample info and count matrix. Please ensure sample names align.\n")
                                        next # Retry alignment
                                    }

                                    # Assign conditions from the sample info table
                                    conditions_list <- list()
                                    for (i in seq_len(ncol(coldata))) {
                                        conditions_list[[paste0("Condition", i)]] <- coldata[, i]
                                    }

                                    # Combine sample names and all conditions into a final dataframe
                                    final_conditions <- data.frame(
                                        Sample = colnames(combined_data)[-1],
                                        conditions_list,
                                        stringsAsFactors = FALSE
                                    )

                                    colnames(combined_data)[-1] <- final_conditions$Sample  # Update column names
                                    rownames(combined_data) <- combined_data$geneID
                                    combined_data$geneID <- NULL

                                    # Debugging: Print the final data structure
                                    cat("\nFinal preview of assigned conditions:\n")
                                    print(final_conditions)

                                    cat("\nFinal head of data structure with the assigned conditions:\n")
                                    print(head(combined_data))
                                    # navigate_back <- FALSE
                                    break # Success, exit loop        
                                }, error = function(e) {
                                    cat("Error aligning sample info with count matrix: ", e$message, "\n")
                                    next # Retry alignment
                                })
                            } else if (tolower(info_table_choice) == "n") {
                                cat("\nNo sample info table provided. Extracting sample names from the matrix.\n")
                                tryCatch({
                                    # Extract sample names from the count matrix
                                    repeat {
                                        extracted_sample_names <- colnames(combined_data)[-1]  # Assuming the first column is geneID
                                        cat("Extracted sample names:\n")
                                        print(extracted_sample_names)
                                        
                                        # Validate extracted names with the user
                                        repeat {
                                            cat("\nAre these sample names correct? (y/n): ")
                                            validation_response <- handle_quit(readline())

                                            if (validation_response == "back") {
                                                cat("Returning to the matrix input step.\n")
                                                break # Exit to return to the previous step
                                            }
                                            
                                            if (tolower(validation_response) == "y") {
                                                break
                                            } else if (tolower(validation_response) == "n") {
                                                cat("\nPlease provide the correct sample names. Enter the names separated by commas:\n")
                                                new_sample_names <- handle_quit(readline())

                                                if (new_sample_names == "back") {
                                                cat("Returning to extracted sample names validation.\n")
                                                next # Retry sample name validation
                                                }                                    
                                                
                                                extracted_sample_names <- unlist(strsplit(new_sample_names, ","))
                                                
                                                # Validate length of new sample names
                                                if (length(extracted_sample_names) == ncol(combined_data) - 1) {
                                                    colnames(combined_data)[-1] <- extracted_sample_names
                                                    break
                                                } else {
                                                    cat("The number of sample names does not match the number of columns in the count matrix. Please try again.\n")
                                                }
                                            } else {
                                                cat("Invalid input. Please enter 'y' or 'n'.\n")
                                            }
                                        }
                                        break # Exit sample name extraction loop if valid names are confirmed
                                    }

                                    # Assign conditions
                                    assign_conditions <- function(sample_names) {
                                        conditions <- list()  # List to store multiple conditions
                                        condition_index <- 1  # Counter for condition number
                                        
                                        repeat {
                                        cat(paste("\nAssigning Condition", condition_index, "...\n"))
                                        condition_values <- vector("character", length(sample_names))  # Initialize condition values
                                        
                                        cat("Please assign a condition value for each sample:\n")
                                        for (i in seq_along(sample_names)) {
                                            cat(sprintf("[%d] %s: ", i, sample_names[i]))
                                            condition_values[i] <- handle_quit(readline())
                                        }
                                        
                                        # Validation Step
                                        repeat {
                                            cat("\nCondition assignment preview:\n")
                                            preview <- data.frame(Sample = sample_names, Condition = condition_values, stringsAsFactors = FALSE)
                                            print(preview)
                                            
                                            cat("\nAre these conditions correct? (y/n): ")
                                            validation_response <- handle_quit(readline())
                                            if (tolower(validation_response) == "y") {
                                                break
                                            } else if (tolower(validation_response) == "n") {
                                                cat("Please reassign the conditions.\n")
                                                for (i in seq_along(sample_names)) {
                                                    cat(sprintf("[%d] %s: ", i, sample_names[i]))
                                                    condition_values[i] <- handle_quit(readline())
                                                }
                                            } else {
                                                cat("Invalid input. Please enter 'y' or 'n'.\n")
                                            }
                                        }
                                        
                                        # Add the new condition to the conditions list
                                        conditions[[paste0("Condition", condition_index)]] <- condition_values
                                        
                                        cat("\nPreview of assigned conditions:\n")
                                        preview <- data.frame(Sample = sample_names, conditions, stringsAsFactors = FALSE)
                                        print(preview)
                                        
                                        # Ask if more conditions need to be defined
                                        cat("\nDo you want to define another condition? (y/n): ")
                                        user_response <- handle_quit(readline())
                                        if (tolower(user_response) == "n") {
                                            break
                                        }
                                        
                                        # Increment condition index
                                        condition_index <- condition_index + 1
                                        }
                                        
                                        return(conditions)
                                    }

                                    # Assign conditions
                                    conditions_list <- assign_conditions(extracted_sample_names)
                                    
                                    # Combine sample names and all conditions into a final dataframe
                                    final_conditions <- data.frame(
                                        Sample = extracted_sample_names,
                                        conditions_list,
                                        stringsAsFactors = FALSE
                                    )
                                    
                                    cat("\nFinal preview of assigned conditions:\n")
                                    print(final_conditions)

                                    # Align with the combined_data structure
                                    colnames(combined_data)[-1] <- extracted_sample_names  # Update column names
                                    rownames(combined_data) <- combined_data$geneID
                                    combined_data$geneID <- NULL

                                    # Debugging: Preview final structure
                                    cat("\nFinal head of data structure with the assigned conditions:\n")
                                    print(head(combined_data))
                                    cat(sprintf("\nNumber of genes before filtering: %d\n", nrow(combined_data)))
                                    # break # Exit the loop after successful extraction
                                }, error = function(e) {
                                    cat("An error occurred while processing the sample names or conditions: ", e$message, "\n")
                                    cat("Retrying the step...\n")
                                    next # Retry the whole process
                                })
                                # navigate_back <- FALSE
                                # break
                            } else {
                                cat("Invalid input. Please enter 'y' or 'n'.\n")
                                next
                            }
                        } 
                        # if (navigate_back) next # Return to the count matrix input loop
                        # break
                    } else if (data_input_choice == 2) {
                        # User has a folder of count files
                        # 1. Ask for the folder path
                        repeat {
                            cat("Please insert the folder path (or type 'back' to return to the previous step): \n")
                            folder_path <- handle_quit(readline())
                            
                            if (folder_path == "back") {
                                cat("Returning to the data input choice.\n")
                                break # Return to the data input choice
                            }

                            if (dir.exists(folder_path)) {
                                setwd(folder_path)
                                break
                            } else {
                                cat("Invalid folder path. Please try again.\n")
                            }
                        }

                        # 2. Check for optional sample info table
                        cat("Do you have a sample info table? (y/n, or type 'quit' to exit): \n")
                        info_table_choice <- handle_quit(readline())
                        if (tolower(info_table_choice) == "y") {
                            repeat {
                                cat("Please provide the path to the sample info table (or type 'quit' to exit): \n")
                                info_table_path <- handle_quit(readline())
                                # Handle "back" command to return to the previous step
                                if (info_table_path == "back") {
                                    cat("Returning to the sample info table choice.\n")
                                    final_conditions <- NULL
                                    break # Exit to return to the sample info table choice loop
                                }

                                if (file.exists(info_table_path)) {
                                    tryCatch({
                                        # Function to read sample info table
                                        read_sample_info <- function(file_path) {
                                            if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
                                                read.csv(file_path, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
                                            } else if (grepl("\\.(txt|tsv)$", file_path, ignore.case = TRUE)) {
                                                read.delim(file_path, header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = "\t", fileEncoding = "UTF-8")
                                            } else if (grepl("\\.xlsx$", file_path, ignore.case = TRUE)) {
                                                read_excel(file_path, col_names = TRUE)
                                            } else {
                                                cat("Error: Unsupported file format. Please provide a .csv, .txt, .tsv, or .xlsx file.\n")
                                                return(NULL)
                                            }
                                        }

                                        # Read the sample info table
                                        sample_info <- read_sample_info(info_table_path)
                                        if (is.null(sample_info)) next # Retry if unsupported format

                                        # Debugging: Print parsed sample info table
                                        cat("Sample info table loaded successfully:\n")
                                        print(sample_info)

                                        # Exit the loop on successful read
                                        break
                                    }, error = function(e) {
                                        cat("Error reading the sample info table: ", e$message, "\n")
                                        next # Retry on error
                                    })
                                } else {
                                    cat("Invalid file path. Please try again.\n")
                                }
                            }
                            # Extract count file names from the last column
                            count_files <- sample_info[, ncol(sample_info)]

                            # Debugging: Print extracted file paths
                            cat("Debugging Info: Extracted file paths from the sample info table:\n")
                            print(count_files)

                            # Check for missing or invalid file paths
                            if (any(is.na(count_files) | count_files == "")) {
                                cat("Error: Some file paths in the sample info table are missing or invalid.")
                                next
                            }

                            # Concatenate folder path with filenames
                            # count_files <- file.path(folder_path, count_files)
                            count_files <- normalizePath(sample_info$path, winslash = "/", mustWork = FALSE)


                            # Validate that the count files exist
                            if (!all(file.exists(count_files))) {
                                missing_files <- count_files[!file.exists(count_files)]
                                cat("The following count files from the sample info table are missing: ", paste(missing_files, collapse = ", "))
                                next
                            }

                            # Use all other columns (excluding the last column) as metadata for coldata
                            coldata <- sample_info[, -ncol(sample_info), drop = FALSE]
                            colnames(coldata) <- make.names(colnames(coldata))  # Ensure column names are valid

                            cat("Sample info table loaded successfully. Using the following coldata:\n")
                            print(coldata)

                            # Read and process count files (supports .txt, .tsv, .csv, .xlsx)
                            cat("Reading gene count data...\n")
                            tryCatch({
                                read_one_count_file <- function(file) {
                                    if (grepl("\\.xlsx$", file, ignore.case = TRUE)) {
                                        data <- read_excel(file, col_names = FALSE)
                                    } else if (grepl("\\.csv$", file, ignore.case = TRUE)) {
                                        data <- read.csv(file, header = FALSE, stringsAsFactors = FALSE)
                                    } else {
                                        data <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)
                                    }
                                    colnames(data) <- c("geneID", paste0("Sample_", seq_len(ncol(data) - 1)))
                                    return(data)
                                }
                                Genes <- lapply(count_files, read_one_count_file)
                            }, error = function(e) {
                                cat("Error reading files. Please ensure they are properly formatted:", e$message, "\n")
                                next
                            })

                            # Combine data into a single matrix
                            combined_data <- Reduce(function(x, y) merge(x, y, by = "geneID", all = TRUE), Genes)
                            # Rename columns based on sample info table
                            colnames(combined_data)[-1] <- rownames(sample_info)  # Use sample names from the sample info table
                            cat("Gene count matrix created successfully.\n")

                            # # save the combined data to a file
                            # combined_data_filename <- get_output_path("combined_data.csv")
                            # write.csv(combined_data, combined_data_filename, row.names = FALSE)

                            cat("The col names of the combined data are: ") #DEBUG
                            print(colnames(combined_data)[-1])

                            # Validate that the sample names match the coldata
                            if (!all(colnames(combined_data)[-1] %in% rownames(coldata))) {
                                cat("Mismatch between sample names in the count files and the provided sample info table. Please correct the sample info table or paths.\n")
                                next
                            }

                            # Assign conditions from the sample info table
                            conditions_list <- list()
                            for (i in seq_len(ncol(coldata))) {
                                conditions_list[[paste0("Condition", i)]] <- coldata[, i]
                            }

                            # Combine sample names and all conditions into a final dataframe
                            final_conditions <- data.frame(
                                Sample = colnames(combined_data)[-1],
                                conditions_list,
                                stringsAsFactors = FALSE
                            )
                            sample_names <- colnames(combined_data)[-1]
                            cat("\nThe sample names are: ") #DEBUG
                            print(sample_names)

                            cat("\nFinal preview of assigned conditions:\n")
                            print(final_conditions)

                            # Preview the final data structure
                            cat("\nThe head of data structure with the assigned conditions:\n")
                            print(head(combined_data))
                            write.csv(combined_data, get_output_path("combined_data.csv"), row.names = TRUE)
                            cat("Combined data saved as 'combined_data.csv'\n")
                        } else {
                            # 2. Read the files into a list (supports txt, csv, tsv, xlsx)
                            cfiles <- list.files(folder_path, pattern = "\\.(txt|csv|tsv|xlsx)$", full.names = TRUE)
                            if (length(cfiles) == 0) {
                                cat("No valid files found in the specified folder. Please check the folder and try again.\n")
                                next
                            }

                            # Normalize file paths to use `/` consistently
                            cfiles <- gsub("\\\\", "/", cfiles)
                            sample_names <- sub("\\.(txt|csv|tsv|xlsx)$", "", basename(cfiles), ignore.case = TRUE)

                            # Read a single count file by extension (supports .txt, .tsv, .csv, .xlsx)
                            read_count_file <- function(file) {
                                if (grepl("\\.xlsx$", file, ignore.case = TRUE)) {
                                    read_excel(file, col_names = FALSE)
                                } else if (grepl("\\.csv$", file, ignore.case = TRUE)) {
                                    read.csv(file, header = FALSE, stringsAsFactors = FALSE)
                                } else {
                                    read.delim(file, header = FALSE, stringsAsFactors = FALSE)
                                }
                            }

                            # Validate the structure of input files
                            validate_file_format <- function(file) {
                                tryCatch({
                                    data <- read_count_file(file)
                                    if (ncol(data) < 2) return(FALSE)
                                    return(TRUE)
                                }, error = function(e) {
                                    cat("Error in file:", file, "-", e$message, "\n")
                                    return(FALSE)
                                })
                            }

                            # Filter valid files
                            valid_files <- sapply(cfiles, validate_file_format)
                            cfiles <- cfiles[valid_files]
                            sample_names <- sample_names[valid_files]

                            if (length(cfiles) == 0) {
                                cat("No valid input files found. Please ensure files are correctly formatted and try again.\n")
                                next
                            }

                            cat("Loaded valid files:\n")
                            print(cfiles)

                            # 3. Read and process files
                            cat("Reading gene count data...\n")
                            tryCatch({
                                # Genes <- lapply(cfiles, function(file) {
                                #     data <- if (grepl("\\.xlsx$", file)) {
                                #         read_excel(file, col_names = FALSE)
                                #     } else {
                                #         read.delim(file, header = FALSE)
                                #     }
                                #     # colnames(data) <- c("geneID", paste0("Sample_", seq_len(ncol(data) - 1)))
                                #     return(data)
                                ## TODO: Check if this works:
                                gene_list <- lapply(seq_along(cfiles), function(i) {
                                file <- cfiles[i]
                                sample_name <- sample_names[i]
                                data <- read_count_file(file)
                                colnames(data) <- c("geneID", sample_name)
                                return(data)
                                })
                            }, error = function(e) {
                                cat("Error reading files. Please ensure they are properly formatted:", e$message, "\n")
                                next
                            })

                            # Combine data into a single matrix
                            # combined_data <- Reduce(function(x, y) merge(x, y, by = "geneID", all = TRUE), Genes)
                            combined_data <- Reduce(function(x, y) merge(x, y, by = "geneID", all = TRUE), gene_list)## TODO: Check if this works:
                            cat("Gene count matrix created successfully.\n")

                            # 4. Extract sample information
                            cat("Please insert from where to where to read the file name for relevant information (or type 'quit' to exit):\n")

                            repeat {
                                cat("From (leave empty for no delimiter): ")
                                from <- handle_quit(readline())
                                cat("To (leave empty for no delimiter): ")
                                to <- handle_quit(readline())

                                # Escape special characters in delimiters
                                from <- ifelse(from != "", gsub("([][\\^$.|?*+(){}])", "\\\\\\1", from), "")
                                to <- ifelse(to != "", gsub("([][\\^$.|?*+(){}])", "\\\\\\1", to), "")

                                # Extract file parts based on delimiters
                                tryCatch({
                                    dat.names <- basename(cfiles)  # Default to file names
                                    if (from != "") {
                                        dat.names <- sub(from, "", dat.names)
                                    }
                                    if (to != "") {
                                        dat.names <- sub(to, "", dat.names)
                                    }
                                    if (length(dat.names) > 0) {
                                        break
                                    } else {
                                        cat("No matches found for the specified delimiters. Please try again.\n")
                                    }
                                }, error = function(e) {
                                    cat("Error during processing: ", e$message, "\nPlease try again.\n")
                                })
                            }

                            cat("Extracted file parts:\n")
                            print(dat.names)

                            # Function to generate combinations from file names
                            generate_combinations <- function(name) {
                                delimiters <- "[/_\\.-]+"  # Matches '/', '_', '\', '.' or '-'
                                parts <- unlist(strsplit(name, delimiters))
                                unlist(lapply(seq_along(parts), function(i) {
                                    sapply(seq(i, length(parts)), function(j) {
                                        paste(parts[i:j], collapse = "_")
                                    })
                                }))
                            }

                            # 5. Extract sample names
                            extract_sample_names <- function(dat.names) {
                                remaining_files <- dat.names
                                selected_samples <- vector("character", length(dat.names))  # Preallocate with default empty
                                names(selected_samples) <- dat.names  # Map files to names

                                repeat {
                                    if (length(remaining_files) == 0) {
                                        cat("\nAll sample names have been defined.\n")
                                        break
                                    }
                                    
                                    file <- remaining_files[1]
                                    cat("\nProcessing file:", file, "\n")
                                    
                                    # Generate all possible combinations for this file
                                    sample_combinations <- generate_combinations(file)
                                    names(sample_combinations) <- seq_along(sample_combinations)  # Add indices for clarity

                                    # Display all combinations
                                    cat("Select the part of the file that represents the sample name:\n")
                                    for (i in seq_along(sample_combinations)) {
                                        cat(sprintf("[%d] %s\n", i, sample_combinations[i]))
                                    }
                                    cat("[manual] Enter manually\n[done] Finish selection\n[quit] Exit program\n")

                                    # User input for selection
                                    user_input <- handle_quit(readline("Enter the index for sample names (or 'manual', 'done', 'quit'): "))
                                    
                                    if (tolower(user_input) == "done") {
                                        break
                                    }
                                    if (tolower(user_input) == "manual") {
                                        # Allow manual entry
                                        manual_name <- handle_quit(readline("Enter the sample name manually: "))
                                        selected_samples[file] <- manual_name
                                        remaining_files <- remaining_files[-1]  # Remove the current file
                                        next
                                    }

                                    # Handle index-based selection
                                    selected_index <- as.integer(user_input)
                                    if (!is.na(selected_index) && selected_index > 0 && selected_index <= length(sample_combinations)) {
                                        selected_sample <- sample_combinations[selected_index]
                                        selected_samples[file] <- selected_sample

                                        # Remove files that match this selection
                                        remaining_files <- remaining_files[!sapply(remaining_files, function(name) {
                                            selected_parts <- unlist(strsplit(selected_sample, "_"))
                                            parts <- generate_combinations(name)
                                            any(selected_parts %in% parts)
                                        })]

                                        cat("\nPreview of selected sample names:\n")
                                        print(selected_samples[!selected_samples == ""])
                                    } else {
                                        cat("Invalid input. Please enter a valid index, 'manual', 'done', or 'quit'.\n")
                                    }
                                }

                                # Handle unassigned files
                                unassigned_files <- names(selected_samples)[selected_samples == ""]
                                if (length(unassigned_files) > 0) {
                                    cat("\nSome files were not assigned sample names. Define names for these files:\n")
                                    for (file in unassigned_files) {
                                        cat("\nProcessing unassigned file:", file, "\n")
                                        sample_combinations <- generate_combinations(file)
                                        names(sample_combinations) <- seq_along(sample_combinations)

                                        # Display all combinations for the unassigned file
                                        cat("Select the part of the file that represents the sample name:\n")
                                        for (i in seq_along(sample_combinations)) {
                                            cat(sprintf("[%d] %s\n", i, sample_combinations[i]))
                                        }
                                        cat("[manual] Enter manually\n[skip] Skip this file\n")

                                        # Get user input for the unassigned file
                                        user_input <- handle_quit(readline("Enter the index for sample names (or 'manual', 'skip'): "))
                                        if (tolower(user_input) == "manual") {
                                            manual_name <- handle_quit(readline("Enter the sample name manually: "))
                                            selected_samples[file] <- manual_name
                                        } else if (tolower(user_input) == "skip") {
                                            selected_samples[file] <- "unassigned"
                                        } else {
                                            selected_index <- as.integer(user_input)
                                            if (!is.na(selected_index) && selected_index > 0 && selected_index <= length(sample_combinations)) {
                                                selected_samples[file] <- sample_combinations[selected_index]
                                            } else {
                                                cat("Invalid input. Skipping this file.\n")
                                                selected_samples[file] <- "unassigned"
                                            }
                                        }
                                    }
                                }
                                
                                return(selected_samples)
                            }

                            # Apply the function to extract sample names
                            sample_names <- extract_sample_names(dat.names)

                            cat("\nFinal list of extracted sample names:\n")
                            print(sample_names)

                            # 6. Condition Assignment
                            assign_conditions <- function(dat.names, samples) {
                                conditions <- list()  # List to store multiple conditions
                                condition_index <- 1  # Counter for condition number

                                repeat {
                                    cat(paste("\nAssigning Condition", condition_index, "...\n"))
                                    condition_values <- rep(NA, length(dat.names))  # Initialize condition values
                                    remaining_files <- dat.names[is.na(condition_values)]  # Only process files without current condition

                                    while (length(remaining_files) > 0) {
                                        file <- remaining_files[1]
                                        cat("\nProcessing file:", file, "\n")

                                        # Generate possible combinations
                                        sample_combinations <- generate_combinations(file)
                                        names(sample_combinations) <- seq_along(sample_combinations)

                                        # Display options
                                        cat("Select the part of the file that represents Condition", condition_index, ":\n")
                                        for (i in seq_along(sample_combinations)) {
                                            cat(sprintf("[%d] %s\n", i, sample_combinations[i]))
                                        }
                                        cat("[manual] Enter manually\n[skip] Skip this file\n[done] Finish assignment\n[quit] Exit program\n")

                                        user_input <- handle_quit(readline("Enter the index for Condition (or 'manual', 'skip', 'done', 'quit'): "))

                                        if (tolower(user_input) == "done") {
                                            break
                                        }
                                        if (tolower(user_input) == "manual") {
                                            manual_condition <- handle_quit(readline("Enter the condition manually: "))
                                            condition_values[dat.names == file] <- manual_condition
                                            remaining_files <- remaining_files[-1]
                                            next
                                        }
                                        if (tolower(user_input) == "skip") {
                                            remaining_files <- remaining_files[-1]
                                            next
                                        }

                                        selected_index <- as.integer(user_input)
                                        if (!is.na(selected_index) && selected_index > 0 && selected_index <= length(sample_combinations)) {
                                            selected_condition <- sample_combinations[selected_index]
                                            matched_files <- sapply(dat.names, function(name) {
                                                parts <- generate_combinations(name)
                                                any(selected_condition %in% parts)
                                            })

                                            # Assign selected condition to matching files
                                            condition_values[matched_files & is.na(condition_values)] <- selected_condition

                                            # Preview currently assigned conditions
                                            cat("\nPreview of assigned conditions:\n")
                                            preview <- data.frame(Sample = dat.names, Condition = condition_values, stringsAsFactors = FALSE)
                                            print(preview[!is.na(preview$Condition), ])

                                            # Remove matched files from the remaining list
                                            remaining_files <- remaining_files[is.na(condition_values)]
                                        } else {
                                            cat("Invalid input. Please select a valid index, 'manual', 'skip', or 'done'.\n")
                                        }
                                    }

                                    # Validation Step: Review and make corrections
                                    cat("\nValidation: Review assigned conditions for Condition", condition_index, ":\n")
                                    preview <- data.frame(Sample = dat.names, Condition = condition_values, stringsAsFactors = FALSE)
                                    print(preview)

                                    repeat {
                                        cat("\nDo you want to make any changes to the assigned conditions? (y/n): ")
                                        validation_response <- handle_quit(readline())
                                        if (tolower(validation_response) == "y") {
                                            cat("\nSelect a sample to modify (enter the row number or 'done' to finish):\n")
                                            for (i in seq_along(dat.names)) {
                                                cat(sprintf("[%d] %s -> %s\n", i, dat.names[i], condition_values[i]))
                                            }

                                            modify_input <- handle_quit(readline("Enter the row number to modify (or 'done'): "))
                                            if (tolower(modify_input) == "done") {
                                                break
                                            }
                                            row_index <- as.integer(modify_input)
                                            if (!is.na(row_index) && row_index > 0 && row_index <= length(dat.names)) {
                                                cat("\nCurrent condition for sample:", dat.names[row_index], "is", condition_values[row_index])
                                                new_condition <- handle_quit(readline("\nEnter the new condition value: "))
                                                condition_values[row_index] <- new_condition
                                                cat("\nCondition updated. Preview of assigned conditions:\n")
                                                preview <- data.frame(Sample = dat.names, Condition = condition_values, stringsAsFactors = FALSE)
                                                print(preview)
                                            } else {
                                                cat("Invalid row number. Please select a valid sample.\n")
                                            }
                                        } else if (tolower(validation_response) == "n") {
                                            break
                                        } else {
                                            cat("Invalid input. Please enter 'y' or 'n'.\n")
                                        }
                                    }

                                    # Add the new condition to the conditions list
                                    conditions[[paste0("Condition", condition_index)]] <- condition_values

                                    cat("\nPreview of assigned conditions:\n")
                                    preview <- data.frame(Sample = dat.names, conditions, stringsAsFactors = FALSE)
                                    print(preview)

                                    # Ask if more conditions need to be defined
                                    cat("\nDo you want to define another condition? (y/n): ")
                                    user_response <- handle_quit(readline())
                                    if (tolower(user_response) == "n") {
                                        break
                                    }

                                    # Increment condition index
                                    condition_index <- condition_index + 1
                                }

                                return(conditions)
                            }

                            # Call the function to assign multiple conditions
                            conditions_list <- assign_conditions(dat.names, sample_names)

                            # Combine sample names and all conditions into a final dataframe
                            final_conditions <- data.frame(
                                Sample = dat.names,
                                conditions_list,
                                stringsAsFactors = FALSE
                            )

                            cat("\nFinal preview of assigned conditions:\n")
                            print(final_conditions)

                            # Preview the final data structure
                            cat("\nThe head of data structure with the assigned conditions:\n")
                            print(head(combined_data))
                            write.csv(combined_data, get_output_path("combined_data.csv"), row.names = TRUE)
                            cat("Combined data saved as 'combined_data.csv'\n")

                        }
                    }

                    names(combined_data) <- c("geneID", final_conditions$Sample)
                    row.names(combined_data) <- combined_data$geneID
                    combined_data$geneID <- NULL

                    # Preview the final data structure
                    cat("\nFinal head of data structure with the assigned conditions:\n")
                    print(head(combined_data))

                    cat(sprintf("\nNumber of genes before filtering: %d\n", nrow(combined_data)))


                    # 6. Filter non-expressed genes
                    cat("🔍 Filtering Explanation:\n")
                    cat("- min_reads: Minimum count a gene must have to be considered expressed.\n")
                    cat("- min_samples: Minimum number of samples where the gene must be expressed.\n")
                    cat("- Example: Enter '0,1' to remove all genes with only zero counts in all samples.\n\n")

                    cat("Enter the minimum reads and samples as 'min_reads,min_samples' (or type 'quit' to exit): \n")

                    repeat {
                        filter_input <- handle_quit(readline())
                        filter_values <- unlist(strsplit(filter_input, ","))
                        
                        if (length(filter_values) == 2 && all(!is.na(as.numeric(filter_values)))) {
                            min_reads <- as.numeric(filter_values[1])
                            min_samples <- as.numeric(filter_values[2])

                            # Ensure data is numeric
                            combined_data <- as.matrix(combined_data)
                            mode(combined_data) <- "numeric"

                            cat("\n🔎 Applying filtering criteria:\n")
                            cat(sprintf("- Genes must have at least %d read(s) in at least %d sample(s) to be kept.\n", min_reads, min_samples))

                            # Apply filtering logic
                            filter <- apply(combined_data, 1, function(x) {
                                if (min_reads == 0 && min_samples == 0) {
                                    TRUE  # Include all rows
                                } else {
                                    length(x[x > min_reads]) >= min_samples
                                }
                            })

                            filtered <- combined_data[filter, ]

                            if (nrow(filtered) == 0) {
                                cat("No genes meet the filtering criteria. Please adjust your input.\n")
                            } else {
                                cat("Filtered genes successfully.\n")
                                break
                            }
                        } else {
                            cat("Invalid input. Please enter two numeric values separated by a comma.\n")
                        }
                    }

                    # Print summary of the filtered data
                    cat("\nSummary of filtered data:\n")
                    print(head(filtered))
                    cat(sprintf("\nNumber of genes after filtering: %d\n", nrow(filtered)))


                    # 7. Create expression set and assign conditions
                    if (nrow(filtered) == 0) {
                        cat("No data available after filtering. Please adjust your filtering criteria and try again.\n")
                        next
                    }

                    filtered <- filtered[!grepl("__", rownames(filtered)),] #remove the rows with __ in the name
                    genes <- rownames(filtered)  ##names of genes##


                    # 8. Optional: Convert gene IDs to gene names
                    cat("To convert gene IDs to their names? (y/n): \n")
                    convert_genes_choice <- handle_quit(readline())

                    if (tolower(convert_genes_choice) == "y") {
                        repeat {
                                
                            cat("Please insert the annotation file path (or type 'quit' to exit): \n")
                            annotation_file_path <- handle_quit(readline())
                            if (!file.exists(annotation_file_path)) {
                                cat("The annotation file does not exist. Please check the path.")
                                next
                            }

                            # Determine file type and process accordingly
                            if (grepl("\\.gtf$", annotation_file_path, ignore.case = TRUE)) {
                                cat("Loading annotation from GTF file...\n")
                                gtf <- tryCatch({
                                    import(annotation_file_path)
                                }, error = function(e) {
                                    cat("Error reading the GTF file: ", e$message)
                                    next
                                })
                                # Performance: subset to gene-level rows if type column exists
                                mc <- mcols(gtf)
                                if ("type" %in% names(mc)) {
                                    gtf <- gtf[as.character(mc$type) == "gene"]
                                    mc <- mcols(gtf)
                                }
                                # Support different GTF annotation formats: detect gene ID and name columns
                                id_candidates <- c("gene_id", "GeneID", "geneID", "gene")
                                name_candidates <- c("gene_name", "gene_symbol", "symbol", "Name", "gene_name")
                                id_col <- NULL
                                for (cand in id_candidates) {
                                    if (cand %in% names(mc)) { id_col <- cand; break }
                                }
                                name_col <- NULL
                                for (cand in name_candidates) {
                                    if (cand %in% names(mc)) { name_col <- cand; break }
                                }
                                if (is.null(id_col)) {
                                    cat("GTF does not contain a recognized gene ID column (tried: ", paste(id_candidates, collapse = ", "), "). Available: ", paste(names(mc), collapse = ", "), "\n")
                                    next
                                }
                                if (is.null(name_col)) name_col <- id_col
                                # Build annotation table with canonical names for merge
                                annot <- as.data.frame(mc[, c(id_col, name_col)], stringsAsFactors = FALSE)
                                colnames(annot) <- c("gene_id", "gene_name")
                                annot$gene_id <- as.character(annot$gene_id)
                                annot$gene_name <- as.character(annot$gene_name)
                                annot$gene_name[is.na(annot$gene_name) | annot$gene_name == ""] <- annot$gene_id
                                annot <- annot[!duplicated(annot$gene_id), ]

                                # Merge filtered genes with annotation (single merge)
                                filtered <- merge(filtered, annot, by.x = "row.names", by.y = "gene_id", all.x = TRUE)

                                filtered$final_gene_name <- ifelse(is.na(filtered$gene_name) | filtered$gene_name == "",
                                                                filtered$Row.names,
                                                                filtered$gene_name)
                                rownames(filtered) <- make.unique(filtered$final_gene_name)
                                filtered$gene_name <- NULL
                                filtered$final_gene_name <- NULL
                                filtered$Row.names <- NULL

                                cat("Gene IDs successfully converted with proper alignment.\n")
                                write.csv(filtered, get_output_path("filtered_data.csv"), row.names = TRUE)
                                cat("Filtered data saved as 'filtered_data.csv'\n")

                                # -------------------------------
                                # -------------------------------
                                # Memory-Efficient Approach
                                # -------------------------------

                                # # 0) Ensure 'filtered' is a data frame
                                # filtered <- as.data.frame(filtered)

                                # # 1) Extract only the necessary annotation columns and remove duplicates
                                # annot <- as.data.frame(mcols(gtf)[, c("exon_id", "gene_name")])
                                # annot <- annot[!duplicated(annot$exon_id), ]
                                # annot$exon_id   <- as.character(annot$exon_id)
                                # annot$gene_name <- as.character(annot$gene_name)

                                # # Replace missing gene names with exon_id
                                # missing_idx <- is.na(annot$gene_name) | annot$gene_name == ""
                                # annot$gene_name[missing_idx] <- annot$exon_id[missing_idx]

                                # # 2) Use match() to add gene names to filtered data without a heavy merge.
                                # # Assume the row names of filtered are the exon IDs.
                                # exon_ids <- rownames(filtered)
                                # idx <- match(exon_ids, annot$exon_id)
                                # gene_names <- annot$gene_name[idx]

                                # # Remove trailing ".number" from gene names (e.g., "Gad1.6" becomes "Gad1")
                                # gene_names <- sub("\\.[0-9]+$", "", gene_names)

                                # # Add gene_names as a new column
                                # filtered$gene_name <- gene_names

                                # # 3) Compute row sums over expression columns.
                                # # Identify expression columns: all except the new "gene_name" column.
                                # expr_cols <- setdiff(colnames(filtered), "gene_name")

                                # # Ensure the expression columns are numeric
                                # filtered[, expr_cols] <- lapply(filtered[, expr_cols, drop = FALSE], as.numeric)

                                # # Compute the sum of each row across these columns
                                # filtered$RowSum <- rowSums(filtered[, expr_cols, drop = FALSE])

                                # # 4) For genes with multiple exons, keep only the row with the highest total expression.
                                # # Split the data frame by gene name.
                                # split_by_gene <- split(filtered, filtered$gene_name)

                                # # For each gene, select the row with the maximum RowSum (if ties, the first is kept)
                                # max_rows <- lapply(split_by_gene, function(df) df[which.max(df$RowSum), ])

                                # # Reassemble the filtered data
                                # filtered <- do.call(rbind, max_rows)

                                # # 5) Clean up: set row names and remove extra columns.
                                # rownames(filtered) <- filtered$gene_name
                                # filtered$gene_name <- NULL
                                # filtered$RowSum <- NULL

                                # cat("Gene IDs successfully converted, trailing exon numbers removed, and most-expressed rows retained.\n")


                                # write.csv(filtered, get_output_path("filtered_data.csv"), row.names = TRUE)
                                # cat("Filtered data saved as 'filtered_data.csv'\n")
                                # -------------------------------
                                # -------------------------------
                                # Exons Approach
                                # -------------------------------

                                # # Extract relevant annotation columns
                                # annot <- unique(as.data.frame(mcols(gtf)[, c("exon_id", "gene_name")]))

                                # # Ensure no factor conversion issues
                                # annot$exon_id <- as.character(annot$exon_id)
                                # annot$gene_name <- as.character(annot$gene_name)

                                # # Handle missing gene names: Replace NA values with gene_id
                                # annot$gene_name[is.na(annot$gene_name) | annot$gene_name == ""] <- annot$exon_id[is.na(annot$gene_name) | annot$gene_name == ""]
                                
                                # # Ensure gene_name is unique by appending gene_id where needed
                                # annot$gene_name <- make.unique(annot$gene_name)

                                # # Merge with filtered genes
                                # filtered <- merge(filtered, annot, by.x = "row.names", by.y = "exon_id", all.x = TRUE)

                                # # Ensure row names are unique and maintain alignment
                                # rownames(filtered) <- filtered$gene_name
                                # filtered$gene_name <- NULL
                                # filtered$Row.names <- NULL
                                # cat("Gene IDs successfully converted with proper alignment.\n")
                                # write.csv(filtered, get_output_path("filtered_data.csv"), row.names = TRUE)
                                # cat("Filtered data saved as 'filtered_data.csv'\n")
                            } else if (grepl("\\.xlsx$", annotation_file_path, ignore.case = TRUE)) {
                                cat("Loading annotation from XLSX file...\n")
                                annotation <- tryCatch({
                                    read_excel(annotation_file_path)
                                }, error = function(e) {
                                    cat("Error reading the XLSX file:", e$message, "\n")
                                    return(NULL)
                                })
                                if (is.null(annotation)) next

                                # Extract gene ID and gene name columns
                                annot <- annotation[, c("Gene_stable_ID", "Gene_name")]
                                rownames(annot) <- annot$Gene_stable_ID

                                # Merge with filtered genes and remove duplicates
                                filtered <- merge(filtered, annot, by.x = "row.names", by.y = "Gene_stable_ID", all.x = TRUE)
                                filtered <- filtered[!duplicated(filtered$Gene_name), ]
                                rownames(filtered) <- filtered$Gene_name
                                filtered$Gene_name <- NULL
                                filtered$Row.names <- NULL

                                # Extract gene names
                                genes <- rownames(filtered)

                                cat("Gene IDs successfully converted to gene names using XLSX annotation.\n")


                            } else {
                                cat("Unsupported file format. Please provide either a GTF or XLSX file.\n")
                            }
                            break
                        }
                    } else {
                        cat("Skipping gene ID conversion step.\n")
                    }

                    # Create a sample info table for DESeq2
                    coldata <- final_conditions[final_conditions$Sample %in% colnames(filtered), ]
                    rownames(coldata) <- coldata$Sample
                    coldata$Sample <- NULL

                    # make combined conditions vector from the conditions list if there is more than 1 condition
                    if (ncol(coldata) >= 1) {
                        if (ncol(coldata) > 1) {
                            # Dynamically identify columns to exclude
                            exclude_cols <- sapply(seq_along(coldata), function(idx) {
                                col <- coldata[[idx]]
                                # Extract all other columns
                                other_columns <- coldata[, -idx, drop = FALSE]
                                # Generate all possible combinations of unique words from other columns
                                other_combinations <- apply(other_columns, 1, paste, collapse = "_")
                                # Check if more than half of the column values match the combinations
                                mean(col %in% other_combinations) > 0.5
                            })
                            # Select columns to combine (exclude the identified ones)
                            columns_to_combine <- coldata[, !exclude_cols, drop = FALSE]
                            # Avoid duplicating "conditions_combined"
                            if ("conditions_combined" %in% colnames(coldata) || ncol(columns_to_combine) < ncol(coldata)) {
                                conditions_combined <- coldata[, ncol(coldata)]
                                # conditions_combined <- apply(coldata[, 1:(ncol(coldata) - 1), drop = FALSE], 1, paste, collapse = "_")          
                                # conditions_combined <- apply(coldata[, -ncol(coldata), drop = FALSE], 1, paste, collapse = "_")
                            } else {
                                conditions_combined <- apply(columns_to_combine, 1, paste, collapse = "_")
                            }
                        } else {
                            # Handle single-column case
                            conditions_combined <- coldata[, 1]
                        }
                        conditions_combined <- factor(conditions_combined)
                        coldata$conditions_combined <- conditions_combined
                    } else {
                        cat("Error: coldata has no valid columns remaining. Please check the sample info table and try again.\n")
                        next
                    }
                    # User option: choose which column (or combined) to use for PCA, RLE, and DESeq design
                    if (ncol(coldata) > 1) {
                        design_cols <- setdiff(colnames(coldata), "conditions_combined")
                        if (length(design_cols) == 0) design_cols <- colnames(coldata)[1]
                        cat("\nAvailable columns in sample info (for grouping in PCA, RLE, and comparisons):\n")
                        for (i in seq_along(design_cols)) cat(sprintf("  [%d] %s\n", i, design_cols[i]))
                        cat(sprintf("  [%d] Combined (paste all columns)\n", length(design_cols) + 1))
                        repeat {
                            choice_input <- handle_quit(readline("Enter the number of your choice (or 'back' to keep current): "))
                            if (choice_input == "back") break
                            idx <- as.numeric(choice_input)
                            if (!is.na(idx) && idx >= 1 && idx <= length(design_cols) + 1) {
                                if (idx <= length(design_cols)) {
                                    conditions_combined <- coldata[[design_cols[idx]]]
                                } else {
                                    conditions_combined <- apply(coldata[, design_cols, drop = FALSE], 1, paste, collapse = "_")
                                }
                                conditions_combined <- factor(conditions_combined)
                                coldata$conditions_combined <- conditions_combined
                                cat("Using", if (idx <= length(design_cols)) design_cols[idx] else "Combined", "for grouping.\n")
                                break
                            }
                            cat("Invalid choice. Please enter a number between 1 and", length(design_cols) + 1, ".\n")
                        }
                    }
                    conditions_combined <- factor(conditions_combined)
                    # Preview of chosen condition column
                    cat("\nCondition column preview (sample -> grouping):\n")
                    print(data.frame(Sample = rownames(coldata), Condition = coldata$conditions_combined, row.names = NULL))
                    cat("\nCounts per condition:\n")
                    print(table(coldata$conditions_combined))
                
                    # Order coldata by the combined conditions and reorder the filtered data accordingly
                    coldata <- coldata[order(coldata$conditions_combined, decreasing = TRUE), ]
                    filtered <- filtered[, rownames(coldata)]

                    # Clean up the filtered data
                    colnames(filtered) <- make.unique(colnames(filtered))
                    filtered <- as.matrix(filtered)
                    filtered[is.na(filtered)] <- 0
                    filtered[is.nan(filtered)] <- 0
                    filtered[is.infinite(filtered)] <- 0

                    conditions_combined <- coldata$conditions_combined

                    i <- 1
                    while (i < ncol(coldata)) {  # Use while to dynamically adjust column count
                        for (j in (i + 1):ncol(coldata)) {
                            if (all(coldata[, i] == coldata[, j])) {  # Check if column values are identical
                                cat(sprintf("The column %s is the same as the column %s. Removing the first duplicated column.\n",
                                            colnames(coldata)[i], colnames(coldata)[j]))
                                coldata <- coldata[, -i, drop = FALSE]  # Remove first duplicated column
                                i <- i - 1  # Adjust index to stay at the correct position
                                break  # Restart comparison for the new structure
                            }
                        }
                        i <- i + 1  # Move to the next column
                    }

                    # Print the coldata for debugging
                    cat("\nSample info table for DESeq2:\n")
                    print(coldata)

                    # Create the expression set for further analysis
                    set <- newSeqExpressionSet(as.matrix(filtered),
                                                phenoData = data.frame(conditions_combined, row.names = colnames(filtered)))

                    cat("\nExpression set created successfully.\n")

                    # Verify counts in set
                    cat("Expression set counts preview:\n")
                    print(head(counts(set)))

                    # Create a color palette for conditions
                    # if there is less than 3 conditions, make error and exit
                    if (length(unique(conditions_combined)) < 3) {
                        cat("Less than 3 unique conditions detected. Please ensure there are at least 3 unique conditions for analysis and try again.\n")
                        next
                    }

                    # print the unique conditions
                    cat("\nUnique conditions detected:\n")
                    print(unique(conditions_combined))
                    par(mfrow = c(1, 2))  # Create a 1x2 layout for plots

                    custom_colors <- c("#80b1d3", "#fb8072", "#fdb462", "#b3de69", "#bc80bd", "#c2b975", "#ffb3b3", "#bebada", "#8dd3c7", "#27775c", "#31678c", "#574e92", "#a80a5b", "#624d0e", "#ffb580")
                    
                    # Assign colors based on the number of levels in conditions_combined
                    if (length(custom_colors) < nlevels(conditions_combined)) {
                        cat("Not enough custom colors defined for all conditions!")
                        # restart the loop
                        next
                    }

                    # Map the custom colors to condition levels
                    colors <- setNames(custom_colors[1:nlevels(conditions_combined)], levels(conditions_combined))

                    # # Generate color palette
                    # colors <- brewer.pal(max(3, nlevels(conditions_combined)), "Set3")[1:nlevels(conditions_combined)]
                    # names(colors) <- levels(conditions_combined)

                    # Plot RLE and PCA with error handling
                    cat("Creating plots...\n")
                    tryCatch({
                        plotRLE(set, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])              
                        png(get_pca_rle_path("raw", "RLE_plot.png"), width = 900, height = 500)
                        plotRLE(set, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                        dev.off()
                        svg(get_pca_rle_path("raw", "RLE_plot.svg"), width = 900/72, height = 500/72)
                        plotRLE(set, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                        dev.off()
                    }, error = function(e) {
                        cat("Error in plotRLE:", e$message, "\n")
                    })
                    tryCatch({
                        plotPCA(set, col = colors[conditions_combined], cex = 1.2)
                        # Adding group labels (levels of conditions_combined)
                        legend("topright", legend = levels(conditions_combined), col = colors, pch = 16, cex = 0.8)

                        png(get_pca_rle_path("raw", "PCA_plot.png"), width = 900, height = 900)
                        plotPCA(set, col = colors[conditions_combined], cex = 1.2)
                        # Adding group labels (levels of conditions_combined)
                        legend("topright", legend = levels(conditions_combined), col = colors, pch = 16, cex = 0.8)
                        dev.off()
                        svg(get_pca_rle_path("raw", "PCA_plot.svg"), width = 900/72, height = 900/72)
                        plotPCA(set, col = colors[conditions_combined], cex = 1.2)
                        legend("topright", legend = levels(conditions_combined), col = colors, pch = 16, cex = 0.8)
                        dev.off()

                        plotPCA_custom(set,
                                    group      = conditions_combined,
                                    colors        = colors,
                                    prefix        = get_pca_rle_path("raw", "newPCA_plot"),
                                    width         = 6,      # in inches
                                    height        = 6,
                                    dpi           = 300)
                                    
                    }, error = function(e) {
                        cat("Error in plotPCA:", e$message, "\n")
                    })

                    # Option to change condition column after raw PCA/RLE
                    repeat {
                        cat("Would you like to use a different condition column for grouping? (y/n):\n")
                        change_cond_ans <- handle_quit(readline())
                        if (tolower(change_cond_ans) != "y") break
                        if (ncol(coldata) <= 1) { cat("Only one column in sample info. No other choice available.\n"); break }
                        design_cols <- setdiff(colnames(coldata), "conditions_combined")
                        if (length(design_cols) == 0) design_cols <- colnames(coldata)[1]
                        cat("\nAvailable columns in sample info (for grouping in PCA, RLE, and comparisons):\n")
                        for (i in seq_along(design_cols)) cat(sprintf("  [%d] %s\n", i, design_cols[i]))
                        cat(sprintf("  [%d] Combined (paste all columns)\n", length(design_cols) + 1))
                        updated_cond <- FALSE
                        repeat {
                            choice_input <- handle_quit(readline("Enter the number of your choice (or 'back' to keep current): "))
                            if (choice_input == "back") break
                            idx <- as.numeric(choice_input)
                            if (!is.na(idx) && idx >= 1 && idx <= length(design_cols) + 1) {
                                if (idx <= length(design_cols)) {
                                    conditions_combined <- coldata[[design_cols[idx]]]
                                } else {
                                    conditions_combined <- apply(coldata[, design_cols, drop = FALSE], 1, paste, collapse = "_")
                                }
                                conditions_combined <- factor(conditions_combined)
                                coldata$conditions_combined <- conditions_combined
                                if (exists("set")) pData(set)$conditions_combined <- conditions_combined
                                if (exists("set4")) pData(set4)$conditions_combined <- conditions_combined
                                if (length(custom_colors) < nlevels(conditions_combined)) {
                                    cat("Not enough custom colors defined for all conditions. Keeping current colors.\n")
                                } else {
                                    colors <- setNames(custom_colors[1:nlevels(conditions_combined)], levels(conditions_combined))
                                }
                                cat("Using", if (idx <= length(design_cols)) design_cols[idx] else "Combined", "for grouping.\n")
                                updated_cond <- TRUE
                                break
                            }
                            cat("Invalid choice. Please enter a number between 1 and", length(design_cols) + 1, ".\n")
                        }
                        if (!updated_cond) break
                        tryCatch({
                            plotRLE(set, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                            png(get_pca_rle_path("raw", "RLE_plot.png"), width = 900, height = 500)
                            plotRLE(set, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                            dev.off()
                            svg(get_pca_rle_path("raw", "RLE_plot.svg"), width = 900/72, height = 500/72)
                            plotRLE(set, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                            dev.off()
                            plotPCA(set, col = colors[conditions_combined], cex = 1.2)
                            legend("topright", legend = levels(conditions_combined), col = colors, pch = 16, cex = 0.8)
                            png(get_pca_rle_path("raw", "PCA_plot.png"), width = 900, height = 900)
                            plotPCA(set, col = colors[conditions_combined], cex = 1.2)
                            legend("topright", legend = levels(conditions_combined), col = colors, pch = 16, cex = 0.8)
                            dev.off()
                            svg(get_pca_rle_path("raw", "PCA_plot.svg"), width = 900/72, height = 900/72)
                            plotPCA(set, col = colors[conditions_combined], cex = 1.2)
                            legend("topright", legend = levels(conditions_combined), col = colors, pch = 16, cex = 0.8)
                            dev.off()
                            plotPCA_custom(set, group = conditions_combined, colors = colors,
                                prefix = get_pca_rle_path("raw", "newPCA_plot"), width = 6, height = 6, dpi = 300)
                        }, error = function(e) cat("Error re-plotting:", e$message, "\n"))
                    }

                    # 8. Normalization and RUV adjustment
                    cat("To normalize? y/n (or type 'quit' to exit):\n")
                    norm_choice <- handle_quit(readline())
                    par(mfrow = c(1, 2))  # Create a 1x2 layout for plots
                    if (tolower(norm_choice) == "y") {
                        set1 <- betweenLaneNormalization(set, which = "upper")
                        plotRLE(set1, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                        plotPCA(set1, col = colors[conditions_combined], cex = 1.2)
                        # Adding group labels (levels of conditions_combined)
                        legend("topright", legend = levels(conditions_combined), col = colors, pch = 16, cex = 0.8)

                        png(get_pca_rle_path("normalised", "RLE_plot.png"), width = 900, height = 500)
                        plotRLE(set1, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                        dev.off()
                        svg(get_pca_rle_path("normalised", "RLE_plot.svg"), width = 900/72, height = 500/72)
                        plotRLE(set1, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                        dev.off()

                        png(get_pca_rle_path("normalised", "PCA_plot.png"), width = 900, height = 900)
                        plotPCA(set1, col = colors[conditions_combined], cex = 1.2)
                        # Adding group labels (levels of conditions_combined)
                        legend("topright", legend = levels(conditions_combined), col = colors, pch = 16, cex = 0.8)

                        dev.off()
                        svg(get_pca_rle_path("normalised", "PCA_plot.svg"), width = 900/72, height = 900/72)
                        plotPCA(set1, col = colors[conditions_combined], cex = 1.2)
                        legend("topright", legend = levels(conditions_combined), col = colors, pch = 16, cex = 0.8)
                        dev.off()
                        plotPCA_custom(set1,
                                    group      = conditions_combined,
                                    colors        = colors,
                                    prefix        = get_pca_rle_path("normalised", "newPCA_plot"),
                                    width         = 6,      # in inches
                                    height        = 6,
                                    dpi           = 300)

                        # Option to change condition column after normalised PCA/RLE
                        repeat {
                            cat("Would you like to use a different condition column for grouping? (y/n):\n")
                            change_cond_ans <- handle_quit(readline())
                            if (tolower(change_cond_ans) != "y") break
                            if (ncol(coldata) <= 1) { cat("Only one column in sample info. No other choice available.\n"); break }
                            design_cols <- setdiff(colnames(coldata), "conditions_combined")
                            if (length(design_cols) == 0) design_cols <- colnames(coldata)[1]
                            cat("\nAvailable columns in sample info (for grouping in PCA, RLE, and comparisons):\n")
                            for (i in seq_along(design_cols)) cat(sprintf("  [%d] %s\n", i, design_cols[i]))
                            cat(sprintf("  [%d] Combined (paste all columns)\n", length(design_cols) + 1))
                            updated_cond <- FALSE
                            repeat {
                                choice_input <- handle_quit(readline("Enter the number of your choice (or 'back' to keep current): "))
                                if (choice_input == "back") break
                                idx <- as.numeric(choice_input)
                                if (!is.na(idx) && idx >= 1 && idx <= length(design_cols) + 1) {
                                    if (idx <= length(design_cols)) {
                                        conditions_combined <- coldata[[design_cols[idx]]]
                                    } else {
                                        conditions_combined <- apply(coldata[, design_cols, drop = FALSE], 1, paste, collapse = "_")
                                    }
                                    conditions_combined <- factor(conditions_combined)
                                    coldata$conditions_combined <- conditions_combined
                                    if (exists("set")) pData(set)$conditions_combined <- conditions_combined
                                    if (exists("set4")) pData(set4)$conditions_combined <- conditions_combined
                                    if (length(custom_colors) < nlevels(conditions_combined)) {
                                        cat("Not enough custom colors defined for all conditions. Keeping current colors.\n")
                                    } else {
                                        colors <- setNames(custom_colors[1:nlevels(conditions_combined)], levels(conditions_combined))
                                    }
                                    cat("Using", if (idx <= length(design_cols)) design_cols[idx] else "Combined", "for grouping.\n")
                                    updated_cond <- TRUE
                                    break
                                }
                                cat("Invalid choice. Please enter a number between 1 and", length(design_cols) + 1, ".\n")
                            }
                            if (!updated_cond) break
                            tryCatch({
                                plotRLE(set1, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                                png(get_pca_rle_path("normalised", "RLE_plot.png"), width = 900, height = 500)
                                plotRLE(set1, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                                dev.off()
                                svg(get_pca_rle_path("normalised", "RLE_plot.svg"), width = 900/72, height = 500/72)
                                plotRLE(set1, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                                dev.off()
                                png(get_pca_rle_path("normalised", "PCA_plot.png"), width = 900, height = 900)
                                plotPCA(set1, col = colors[conditions_combined], cex = 1.2)
                                legend("topright", legend = levels(conditions_combined), col = colors, pch = 16, cex = 0.8)
                                dev.off()
                                svg(get_pca_rle_path("normalised", "PCA_plot.svg"), width = 900/72, height = 900/72)
                                plotPCA(set1, col = colors[conditions_combined], cex = 1.2)
                                legend("topright", legend = levels(conditions_combined), col = colors, pch = 16, cex = 0.8)
                                dev.off()
                                plotPCA_custom(set1, group = conditions_combined, colors = colors,
                                    prefix = get_pca_rle_path("normalised", "newPCA_plot"), width = 6, height = 6, dpi = 300)
                            }, error = function(e) cat("Error re-plotting:", e$message, "\n"))
                        }

                        cat("To remove the unwanted variation? y/n (or type 'quit' to exit):\n")
                        ruv_choice <- handle_quit(readline())
                        if (tolower(ruv_choice) == "y") {
                            repeat {
                                cat("Please insert the number of unwanted variation vectors (k) (or type 'back' to skip):\n")
                                k_raw <- handle_quit(readline())
                                if (k_raw == "back") break
                                k <- as.numeric(k_raw)
                                if (is.na(k) || k <= 0) {
                                    cat("Invalid input. Please enter a positive numeric value for k.\n")
                                } else {
                                    # Warn if k is large relative to smallest condition group (risk of overfitting)
                                    min_group_n <- min(table(conditions_combined))
                                    while (k >= min_group_n / 2) {
                                        cat("Warning: k =", k, "is >= half the sample size (n/2) in the smallest condition group (n =", min_group_n, ").\n")
                                        cat("This may lead to overfitting or reduced power. Suggest using k < n/2 in the smallest group.\n")
                                        cat("Change k? (y/n):\n")
                                        change_k <- handle_quit(readline())
                                        if (tolower(change_k) == "n") break
                                        if (tolower(change_k) == "y") {
                                            cat("Please enter a new value for k (or type 'back' to skip RUV):\n")
                                            k_raw <- handle_quit(readline())
                                            if (k_raw == "back") { k <- 0; break }
                                            k <- as.numeric(k_raw)
                                            if (is.na(k) || k <= 0) {
                                                cat("Invalid input. Please enter a positive numeric value for k.\n")
                                            } else {
                                                min_group_n <- min(table(conditions_combined))
                                            }
                                        }
                                    }
                                    if (k > 0) tryCatch({
                                        design <- model.matrix(~conditions_combined, data = pData(set1))
                                        y <- DGEList(counts = counts(set1), group = conditions_combined)
                                        y <- calcNormFactors(y, method = "upperquartile")
                                        y <- estimateGLMCommonDisp(y, design)
                                        y <- estimateGLMTagwiseDisp(y, design)
                                        fit <- glmFit(y, design)
                                        # res <- residuals(fit, type = "deviance")
                                        # Try using deviance residuals first
                                        tryCatch({
                                            res <- residuals(fit, type = "deviance")  # First attempt with deviance
                                            set4 <- RUVr(set1, rownames(filtered), k = k, res)
                                            cat("✅ Using deviance residuals for RUV normalization.\n")
                                        }, error = function(e) {
                                            cat("⚠️ Deviance residuals failed, switching to Pearson residuals.\n")
                                            res <<- residuals(fit, type = "pearson")  # Second attempt with Pearson
                                            set4 <- RUVr(set1, rownames(filtered), k = k, res)
                                        })

                                        # set4 <- RUVr(set1, rownames(filtered), k = k, res)
                                        plotRLE(set4, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                                        plotPCA(set4, col = colors[conditions_combined], cex = 1.2)
                                        # Adding group labels (levels of conditions_combined)
                                        legend("topright", legend = levels(conditions_combined), col = colors, pch = 16, cex = 0.8)

                                        png(get_pca_rle_path("ruv", paste0("RLE_plot_RUV_", k, ".png")), width = 900, height = 500)
                                        plotRLE(set4, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                                        dev.off()
                                        svg(get_pca_rle_path("ruv", paste0("RLE_plot_RUV_", k, ".svg")), width = 900/72, height = 500/72)
                                        plotRLE(set4, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                                        dev.off()

                                        png(get_pca_rle_path("ruv", paste0("PCA_plot_RUV_", k, ".png")), width = 900, height = 900)
                                        plotPCA(set4, col = colors[conditions_combined], cex = 1.2)
                                        # Adding group labels (levels of conditions_combined)
                                        legend("topright", legend = levels(conditions_combined), col = colors, pch = 16, cex = 0.8)

                                        dev.off()
                                        svg(get_pca_rle_path("ruv", paste0("PCA_plot_RUV_", k, ".svg")), width = 900/72, height = 900/72)
                                        plotPCA(set4, col = colors[conditions_combined], cex = 1.2)
                                        legend("topright", legend = levels(conditions_combined), col = colors, pch = 16, cex = 0.8)
                                        dev.off()
                                        plotPCA_custom(set4,
                                                    group      = conditions_combined,
                                                    colors        = colors,
                                                    prefix        = get_pca_rle_path("ruv", paste0("newPCA_plot_RUV_", k)),
                                                    width         = 6,      # in inches
                                                    height        = 6,
                                                    dpi           = 300)

                                        # Option to change condition column after RUV PCA/RLE
                                        repeat {
                                            cat("Would you like to use a different condition column for grouping? (y/n):\n")
                                            change_cond_ans <- handle_quit(readline())
                                            if (tolower(change_cond_ans) != "y") break
                                            if (ncol(coldata) <= 1) { cat("Only one column in sample info. No other choice available.\n"); break }
                                            design_cols <- setdiff(colnames(coldata), "conditions_combined")
                                            if (length(design_cols) == 0) design_cols <- colnames(coldata)[1]
                                            cat("\nAvailable columns in sample info (for grouping in PCA, RLE, and comparisons):\n")
                                            for (i in seq_along(design_cols)) cat(sprintf("  [%d] %s\n", i, design_cols[i]))
                                            cat(sprintf("  [%d] Combined (paste all columns)\n", length(design_cols) + 1))
                                            updated_cond <- FALSE
                                            repeat {
                                                choice_input <- handle_quit(readline("Enter the number of your choice (or 'back' to keep current): "))
                                                if (choice_input == "back") break
                                                idx <- as.numeric(choice_input)
                                                if (!is.na(idx) && idx >= 1 && idx <= length(design_cols) + 1) {
                                                    if (idx <= length(design_cols)) {
                                                        conditions_combined <- coldata[[design_cols[idx]]]
                                                    } else {
                                                        conditions_combined <- apply(coldata[, design_cols, drop = FALSE], 1, paste, collapse = "_")
                                                    }
                                                    conditions_combined <- factor(conditions_combined)
                                                    coldata$conditions_combined <- conditions_combined
                                                    if (exists("set")) pData(set)$conditions_combined <- conditions_combined
                                                    if (exists("set4")) pData(set4)$conditions_combined <- conditions_combined
                                                    if (length(custom_colors) < nlevels(conditions_combined)) {
                                                        cat("Not enough custom colors defined for all conditions. Keeping current colors.\n")
                                                    } else {
                                                        colors <- setNames(custom_colors[1:nlevels(conditions_combined)], levels(conditions_combined))
                                                    }
                                                    cat("Using", if (idx <= length(design_cols)) design_cols[idx] else "Combined", "for grouping.\n")
                                                    updated_cond <- TRUE
                                                    break
                                                }
                                                cat("Invalid choice. Please enter a number between 1 and", length(design_cols) + 1, ".\n")
                                            }
                                            if (!updated_cond) break
                                            tryCatch({
                                                plotRLE(set4, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                                                png(get_pca_rle_path("ruv", paste0("RLE_plot_RUV_", k, ".png")), width = 900, height = 500)
                                                plotRLE(set4, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                                                dev.off()
                                                svg(get_pca_rle_path("ruv", paste0("RLE_plot_RUV_", k, ".svg")), width = 900/72, height = 500/72)
                                                plotRLE(set4, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                                                dev.off()
                                                png(get_pca_rle_path("ruv", paste0("PCA_plot_RUV_", k, ".png")), width = 900, height = 900)
                                                plotPCA(set4, col = colors[conditions_combined], cex = 1.2)
                                                legend("topright", legend = levels(conditions_combined), col = colors, pch = 16, cex = 0.8)
                                                dev.off()
                                                svg(get_pca_rle_path("ruv", paste0("PCA_plot_RUV_", k, ".svg")), width = 900/72, height = 900/72)
                                                plotPCA(set4, col = colors[conditions_combined], cex = 1.2)
                                                legend("topright", legend = levels(conditions_combined), col = colors, pch = 16, cex = 0.8)
                                                dev.off()
                                                plotPCA_custom(set4, group = conditions_combined, colors = colors,
                                                    prefix = get_pca_rle_path("ruv", paste0("newPCA_plot_RUV_", k)), width = 6, height = 6, dpi = 300)
                                            }, error = function(e) cat("Error re-plotting:", e$message, "\n"))
                                        }

                                        ruv_done <- TRUE  # Mark RUV as completed
                                    }, error = function(e) {
                                        # cat("Error during RUV adjustment with k =", k, ". Please try again.\n")
                                        # next
                                        cat("⚠️ RUV Normalization Failed (Both Residual Methods). Skipping RUV adjustment.\n")
                                        set4 <<- set1  # Use normalized set instead of stopping
                                        k <<- 0  # No RUV adjustment
                                        ruv_done <<- FALSE  # Mark RUV as not performed
                                    })
                                    if (k == 0) {
                                        set4 <<- set1
                                        ruv_done <<- FALSE
                                        break
                                    }
                                    # Allow iterative adjustment of k
                                    cat("To adjust the unwanted variation vector number? y/n (or type 'quit' to exit):\n")
                                    adjust_k_choice <- handle_quit(readline())
                                    if (tolower(adjust_k_choice) == "n") {
                                        break
                                    }
                                }
                            }
                        } else {
                            set4 <- set1
                            k <- 0  # No RUV adjustment
                            ruv_done <- FALSE # Mark RUV as not performed
                        }
                    } else {
                        set4 <- set
                        k <- 0  # No normalization or RUV adjustment
                        ruv_done <- FALSE # Mark RUV as not performed
                    }

                    # 9. Clustering and outlier detection
                    repeat {
                        # Perform hierarchical clustering
                        d4 <- tryCatch({
                        if (ruv_done == TRUE) {
                            dist(t(normCounts(set4)))  # Use normalized counts if RUV was performed
                        } else {
                            dist(t(counts(set4)))  # Use raw counts otherwise
                        }
                        }, error = function(e) {
                            cat("Error in distance calculation. Please check the data:", e$message, "\n")
                            next
                        })

                        hc4 <- tryCatch({
                        hclust(d4)
                        }, error = function(e) {
                            cat("Error in hierarchical clustering. Please check the data:", e$message, "\n")
                            next
                        })

                        # Function to plot all clustering and PCA
                        plot_clustering <- function() {
                        par(mfrow = c(1, 3))  # Create a 1x3 layout for plots

                        # Plot RLE
                        tryCatch({
                            plotRLE(set4, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                            title("RLE Plot", cex.main = 1.5)  # Adjust title size

                            # Save RLE plot
                            png(get_pca_rle_path("final", "RLE_plot_final.png"), width = 900, height = 500)
                            plotRLE(set4, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                            dev.off()
                            svg(get_pca_rle_path("final", "RLE_plot_final.svg"), width = 900/72, height = 500/72)
                            plotRLE(set4, outline = FALSE, ylim = c(-1.5, 1.5), col = colors[conditions_combined])
                            dev.off()
                        }, error = function(e) {
                            cat("Error in RLE plot:", e$message, "\n")
                        })

                        # Plot PCA
                        tryCatch({
                            plotPCA(set4, col = colors[conditions_combined], cex = 1.2)
                            # Adding group labels (levels of conditions_combined)
                            legend("topright", legend = levels(conditions_combined), col = colors, pch = 16, cex = 0.8)
                            title("PCA Plot", cex.main = 1.5)  # Adjust title size

                            # Save PCA plot
                            png(get_pca_rle_path("final", "PCA_plot_final.png"), width = 900, height = 900)
                            plotPCA(set4, col = colors[conditions_combined], cex = 1.2)
                            # Adding group labels (levels of conditions_combined)
                            legend("topright", legend = levels(conditions_combined), col = colors, pch = 16, cex = 0.8)
                            dev.off()
                            svg(get_pca_rle_path("final", "PCA_plot_final.svg"), width = 900/72, height = 900/72)
                            plotPCA(set4, col = colors[conditions_combined], cex = 1.2)
                            legend("topright", legend = levels(conditions_combined), col = colors, pch = 16, cex = 0.8)
                            dev.off()
                            plotPCA_custom(set4,
                                        group      = conditions_combined,
                                        colors        = colors,
                                        prefix        = get_pca_rle_path("final", "PCA_plot_final"),
                                        width         = 6,      # in inches
                                        height        = 6,
                                        dpi           = 300)

                            }, error = function(e) {
                                cat("Error in PCA plot:", e$message, "\n")
                            })

                        # Plot hierarchical clustering
                        plot(hc4, main = "Hierarchical Clustering", cex.main = 1.5)  # Adjust title size

                        # save the plot
                        png(get_pca_rle_path("final", "Hierarchical_Clustering.png"), width = 900, height = 900)
                        plot(hc4, main = "Hierarchical Clustering", cex.main = 1.5)
                        dev.off()
                        svg(get_pca_rle_path("final", "Hierarchical_Clustering.svg"), width = 900/72, height = 900/72)
                        plot(hc4, main = "Hierarchical Clustering", cex.main = 1.5)
                        dev.off()

                        # Reset plotting cache
                        par(mfrow = c(1, 1))  # Clear plotting layout
                        }

                        # Initial plot
                        plot_clustering()

                        # Prompt for outlier detection
                        cat("To detect outliers? y/n (or type 'quit' to exit):\n")
                        detect_choice <- handle_quit(readline())
                        if (tolower(detect_choice) == "n" || tolower(detect_choice) == "quit") {
                            break
                        }

                        # Ask for detection method
                        repeat {
                            cat("\nChoose the outlier detection method:\n")
                            cat("[1] Automatic\n")
                            cat("[2] Manual\n")
                            cat("[3] Quit outlier detection\n")
                            detect_method_raw <- handle_quit(readline("Enter your choice (1/2/3), or type 'back' to skip: "))
                            if (detect_method_raw == "back") break
                            detect_method <- as.numeric(detect_method_raw)

                            if (detect_method == 3) {
                                break
                            } else if (detect_method == 1) {
                                # Automatic Outlier Detection
                                cat("\nPerforming automatic outlier detection based on hierarchical clustering...\n")

                                # Define automatic outlier detection logic
                                cluster_cutoff <- 2  # Define a threshold for outlier cluster size
                                cluster_assignments <- cutree(hc4, h = max(d4) * 0.7)  # Adjust the height threshold as needed
                                cluster_sizes <- table(cluster_assignments)

                                # Identify clusters smaller than the cutoff size
                                outlier_clusters <- names(cluster_sizes[cluster_sizes <= cluster_cutoff])
                                outlier_samples <- colnames(set4)[cluster_assignments %in% as.numeric(outlier_clusters)]

                                if (length(outlier_samples) == 0) {
                                    cat("No outliers detected using automatic method.\n")
                                    break
                                }

                                # Remove outliers automatically
                                cat("The following outlier samples were detected and removed:\n")
                                print(outlier_samples)
                                set4 <- set4[, !(colnames(set4) %in% outlier_samples)]
                                # DEBUG: check if set4 updated correctly print TRUE if updated
                                print(paste0("set4 still contains outliers: ", any(colnames(set4) %in% outlier_samples)))
                                filtered <- filtered[, !(colnames(filtered) %in% outlier_samples)]
                                print(paste0("filtered still contains outliers: ", any(colnames(filtered) %in% outlier_samples)))
                                coldata <- coldata[!(rownames(coldata) %in% outlier_samples), , drop = FALSE]  
                                print(paste0("coldata updated: ", any(rownames(coldata) %in% outlier_samples)))
                                # update conditions_combined
                                conditions_combined <- factor(coldata$conditions_combined)

                                cat("Outliers automatically removed. Updating clustering plots...\n")
                                # Recompute clustering and re-plot
                                d4 <- dist(t(if (ruv_done == TRUE) normCounts(set4) else counts(set4)))
                                hc4 <- hclust(d4)
                                plot_clustering()
                                # break

                            } else if (detect_method == 2) {
                                # Manual Outlier Detection
                                repeat {
                                    cat("\nSwitching to manual outlier detection...\n")

                                    # List samples with indices
                                    cat("\nAvailable samples:\n")
                                    sample_names <- colnames(set4)
                                    for (i in seq_along(sample_names)) {
                                        cat(sprintf("[%d] %s\n", i, sample_names[i]))
                                    }

                                    # Validate sample index for removal
                                    repeat {
                                        cat("\nPlease select the index of the sample to remove (or type 'quit' to exit manual mode):\n")
                                        sample_index_input <- handle_quit(readline())
                                        sample_index <- as.numeric(sample_index_input)

                                        if (tolower(sample_index_input) == "quit") {
                                            break
                                        } else if (!is.na(sample_index) && sample_index >= 1 && sample_index <= length(sample_names)) {
                                            sample_to_remove <- sample_names[sample_index]

                                            # Remove the selected sample
                                            set4 <- set4[, colnames(set4) != sample_to_remove]
                                            filtered <- filtered[, colnames(filtered) != sample_to_remove]
                                            # coldata <- coldata[rownames(coldata) != sample_to_remove, ]
                                            coldata <- coldata[rownames(coldata) != sample_to_remove, , drop = FALSE]  # ✅ Fix applied
                                            # update conditions_combined
                                            conditions_combined <- factor(coldata$conditions_combined)

                                            cat(sprintf("Sample '%s' removed. Updating clustering plots...\n", sample_to_remove))

                                            # Recompute clustering and re-plot
                                            d4 <- dist(t(if (ruv_done == TRUE) normCounts(set4) else counts(set4)))
                                            hc4 <- hclust(d4)
                                            plot_clustering()
                                            # break
                                        } else {
                                            cat("Invalid index. Please try again.\n")
                                        }
                                    }

                                    # Prompt if user wants to continue removing samples manually
                                    cat("Do you want to remove more samples manually? y/n (or type 'quit' to exit):\n")
                                    remove_more <- handle_quit(readline())
                                    if (tolower(remove_more) == "n" || tolower(remove_more) == "quit") {
                                        break
                                    }
                                }
                                break
                            } else {
                                cat("Invalid choice. Please enter 1, 2, or 3.\n")
                            }
                        }
                        # Update `filtered` and `coldata` once after all samples are removed
                        filtered <- filtered[, rownames(coldata)]
                        if (ncol(coldata) >= 1) {
                            if (ncol(coldata) > 1) {
                                # Avoid duplicating "conditions_combined"
                                if ("conditions_combined" %in% colnames(coldata)) {
                                    conditions_combined <- apply(coldata[, 1:(ncol(coldata) - 1), drop = FALSE], 1, paste, collapse = "_")
                                } else {
                                    conditions_combined <- apply(coldata, 1, paste, collapse = "_")
                                }
                            } else {
                                # Handle single-column case
                                conditions_combined <- coldata[, 1]
                            }
                            conditions_combined <- factor(conditions_combined)
                            coldata$conditions_combined <- conditions_combined
                        } else {
                            cat("Error: coldata has no valid columns remaining. Please check the sample info table and try again.\n")
                            next
                        }
                        conditions_combined <- factor(conditions_combined)
                        coldata$conditions_combined <- conditions_combined
                    }

                    # Option to change condition column after outlier removal, before DESeq
                    repeat {
                        cat("Would you like to use a different condition column for grouping? (y/n):\n")
                        change_cond_ans <- handle_quit(readline())
                        if (tolower(change_cond_ans) != "y") break
                        if (ncol(coldata) <= 1) { cat("Only one column in sample info. No other choice available.\n"); break }
                        design_cols <- setdiff(colnames(coldata), "conditions_combined")
                        if (length(design_cols) == 0) design_cols <- colnames(coldata)[1]
                        cat("\nAvailable columns in sample info (for grouping in PCA, RLE, and comparisons):\n")
                        for (i in seq_along(design_cols)) cat(sprintf("  [%d] %s\n", i, design_cols[i]))
                        cat(sprintf("  [%d] Combined (paste all columns)\n", length(design_cols) + 1))
                        updated_cond <- FALSE
                        repeat {
                            choice_input <- handle_quit(readline("Enter the number of your choice (or 'back' to keep current): "))
                            if (choice_input == "back") break
                            idx <- as.numeric(choice_input)
                            if (!is.na(idx) && idx >= 1 && idx <= length(design_cols) + 1) {
                                if (idx <= length(design_cols)) {
                                    conditions_combined <- coldata[[design_cols[idx]]]
                                } else {
                                    conditions_combined <- apply(coldata[, design_cols, drop = FALSE], 1, paste, collapse = "_")
                                }
                                conditions_combined <- factor(conditions_combined)
                                coldata$conditions_combined <- conditions_combined
                                if (exists("set")) pData(set)$conditions_combined <- conditions_combined
                                if (exists("set4")) pData(set4)$conditions_combined <- conditions_combined
                                if (length(custom_colors) < nlevels(conditions_combined)) {
                                    cat("Not enough custom colors defined for all conditions. Keeping current colors.\n")
                                } else {
                                    colors <- setNames(custom_colors[1:nlevels(conditions_combined)], levels(conditions_combined))
                                }
                                cat("Using", if (idx <= length(design_cols)) design_cols[idx] else "Combined", "for grouping.\n")
                                updated_cond <- TRUE
                                break
                            }
                            cat("Invalid choice. Please enter a number between 1 and", length(design_cols) + 1, ".\n")
                        }
                        if (!updated_cond) break
                        tryCatch({
                            d4 <- dist(t(if (ruv_done == TRUE) normCounts(set4) else counts(set4)))
                            hc4 <- hclust(d4)
                            plot_clustering()
                        }, error = function(e) cat("Error re-plotting clustering:", e$message, "\n"))
                    }

                    # Save the filtered data, normalized data, and coldata
                    cat("\nData saved to the following files:\n")

                    # Save the filtered data
                    filtered_file <- get_output_path("filtered_data.csv")
                    write.csv(filtered, file = filtered_file)
                    print(filtered_file)

                    # Save the normalized data if ruv_done is TRUE
                    if (ruv_done == TRUE) {
                        normalized_file <- get_results_path("normalized_data.csv")
                        write.csv(normCounts(set4), file = normalized_file)
                        print(normalized_file)
                    }

                    # Save the coldata
                    if ((tolower(info_table_choice) == "n")||(tolower(detect_choice) == "y")) {
                        coldata_file <- get_results_path("coldata.csv")
                        write.csv(coldata, file = coldata_file)
                        print(coldata_file)
                    }

                    # Design choice: single combined condition vs multiple condition columns
                    design_single <- TRUE
                    design_rhs_multi <- NULL
                    design_cols_selected <- NULL
                    cat("\nUse (1) single combined condition for design, or (2) select multiple condition columns for design (e.g. ~cond_i + cond_j)?\n")
                    repeat {
                        design_mode_input <- handle_quit(readline("Enter 1 or 2 (or 'back' to use single): "))
                        if (design_mode_input == "back") break
                        design_mode <- as.numeric(design_mode_input)
                        if (!is.na(design_mode) && design_mode %in% c(1, 2)) {
                            if (design_mode == 1) break
                            design_cols_avail <- setdiff(colnames(coldata), "conditions_combined")
                            if (length(design_cols_avail) == 0) design_cols_avail <- colnames(coldata)[1]
                            if (length(design_cols_avail) < 2) {
                                cat("Only one design column available. Using single combined condition.\n")
                                break
                            }
                            cat("\nAvailable columns for design (enter indices separated by commas, e.g. 1,2):\n")
                            for (i in seq_along(design_cols_avail)) cat(sprintf("  [%d] %s\n", i, design_cols_avail[i]))
                            design_indices_input <- handle_quit(readline("Enter column indices to include in design: "))
                            if (design_indices_input == "back") break
                            design_indices <- as.numeric(trimws(unlist(strsplit(design_indices_input, ","))))
                            design_indices <- design_indices[!is.na(design_indices) & design_indices >= 1 & design_indices <= length(design_cols_avail)]
                            if (length(design_indices) == 0) {
                                cat("No valid indices. Using single combined condition.\n")
                                break
                            }
                            design_cols_selected <- design_cols_avail[unique(sort(design_indices))]
                            design_rhs_multi <- paste(make.names(design_cols_selected), collapse = " + ")
                            design_single <- FALSE
                            cat("Design will use:", design_rhs_multi, "\n")
                            break
                        }
                        cat("Invalid choice. Please enter 1 or 2.\n")
                    }

                    # 10. DESeq2 Dataset Creation
                    if (design_single) {
                        if (k == 0) {
                            # No RUV adjustment; use raw counts
                            cat("\nCreating DESeq2 dataset with raw counts and ~conditions_combined...\n")
                            tryCatch({
                            dds <- DESeqDataSetFromMatrix(countData = counts(set4), 
                                                            colData = pData(set4), 
                                                            design = ~ conditions_combined)
                            }, error = function(e) {
                                cat("Error creating DESeq2 dataset. Please ensure the input data is correctly formatted:", e$message, "\n")
                                next
                            })
                        } else {
                            # RUV adjustment applied
                            dds_choice <- "n"
                            tryCatch({
                                if (tolower(dds_choice) == "y") {
                                    cat("\nCreating DESeq2 dataset with normalized data and ~conditions_combined...\n")
                                    dds <- DESeqDataSetFromMatrix(countData = normCounts(set4), 
                                                                colData = pData(set4), 
                                                                design = ~ conditions_combined)
                                } else {
                                    cat("\nCreating DESeq2 dataset with raw counts and design formula including W factors...\n")
                                    w_factors <- paste0("W_", seq_len(k), collapse = " + ")
                                    design_formula <- as.formula(paste("~", w_factors, "+ conditions_combined"))
                                    dds <- DESeqDataSetFromMatrix(countData = counts(set4), 
                                                                colData = pData(set4), 
                                                                design = design_formula)
                                }
                            }, error = function(e) {
                                cat("Error creating DESeq2 dataset. Please ensure the input data is correctly formatted:", e$message, "\n")
                                next
                            })
                        }
                    } else {
                        # Multi-factor design: ~cond_i + ... + cond_j (optionally + W_1+...+W_k)
                        colData_dds <- as.data.frame(pData(set4))
                        for (cn in design_cols_selected) {
                            cname <- make.names(cn)
                            colData_dds[[cname]] <- coldata[rownames(colData_dds), cn]
                        }
                        if (k == 0) {
                            design_formula <- as.formula(paste("~", design_rhs_multi))
                            cat("\nCreating DESeq2 dataset with raw counts and design", paste("~", design_rhs_multi), "\n")
                            tryCatch({
                                dds <- DESeqDataSetFromMatrix(countData = counts(set4), 
                                                            colData = colData_dds, 
                                                            design = design_formula)
                            }, error = function(e) {
                                cat("Error creating DESeq2 dataset. Please ensure the input data is correctly formatted:", e$message, "\n")
                                next
                            })
                        } else {
                            w_factors <- paste0("W_", seq_len(k), collapse = " + ")
                            design_formula <- as.formula(paste("~", w_factors, "+", design_rhs_multi))
                            cat("\nCreating DESeq2 dataset with raw counts and design", paste("~", w_factors, "+", design_rhs_multi), "\n")
                            tryCatch({
                                dds <- DESeqDataSetFromMatrix(countData = counts(set4), 
                                                            colData = colData_dds, 
                                                            design = design_formula)
                            }, error = function(e) {
                                cat("Error creating DESeq2 dataset. Please ensure the input data is correctly formatted:", e$message, "\n")
                                next
                            })
                        }
                    }
                    # 11. Differential Expression Analysis
                    cat("\nRunning DESeq analysis...\n")
                    tryCatch({
                        # Debugging: Ensure counts are integer mode
                        cat("\nDebugging Info: Count data preview:\n")
                        print(head(counts(dds)))

                        # Debugging: Check colData structure
                        cat("\nDebugging Info: colData preview:\n")
                        print(pData(set4))

                        # Debugging: Check samples by condition
                        cat("\nDebugging Info: Samples by condition:\n")
                        print(table(pData(set4)$conditions_combined))

                        # Run DESeq
                        dds <- DESeq(dds)

                        # check sizeFactors(dds)
                        cat("\nDebugging Info: Size factors of DESeq2 dataset:\n")
                        print(sizeFactors(dds))

                        # write the normalized counts
                        normalized_counts <- counts(dds, normalized = TRUE)
                        normalized_counts_file <- get_results_path("normalized_dds_counts.csv")
                        write.csv(normalized_counts, file = normalized_counts_file)

                    }, error = function(e) {
                        cat("\nDebugging Info: Error during DESeq analysis.\n")
                        cat("Possible reasons could include:\n")
                        cat("- Mismatch in input data format\n")
                        cat("- Issues with the design formula\n")
                        cat("- Factor levels not appropriately set\n")
                        cat("\nDetailed Error Message:\n")
                        print(e)
                        cat("Error during DESeq analysis. Please check the input data and design formula and try again.\n")
                        next
                    })

                    # 12. User-defined comparisons
                    # Initialize storage for formatted results
                    mega_results <- list()
                    gene_names <- c()
                    counts_data_final <- NULL  # To store counts data for the final section

                    degs_data <- data.frame()  # Initialize an empty data frame to store consolidated results

                    repeat {
                        cat("\n📊 Input Format:\n")
                        cat("Enter two condition indices as: X,Y\n")
                        cat("- X = Treatment/Condition of Interest\n")
                        cat("- Y = Reference/Control Condition\n")
                        cat("log2FoldChange shows changes in X relative to Y:\n")
                        cat("  - Positive log2FC → Upregulated in X vs Y\n")
                        cat("  - Negative log2FC → Downregulated in X vs Y\n\n")
                        cat("--------------------------------------------------\n")

                        cat("\nAvailable conditions for comparisons:\n")
                        
                        # Display unique conditions with indices
                        conditions <- levels(pData(set4)$conditions_combined)
                        for (i in seq_along(conditions)) {
                            cat(sprintf("[%d] %s\n", i, conditions[i]))
                        }
                        
                        # Ask user to select conditions for comparisons
                        cat("\nPlease specify the indices of the conditions you want to compare (comma-separated, or type 'quit' to exit):\n")
                        user_input <- handle_quit(readline())
                        
                        if (tolower(user_input) == "quit") {
                            break
                        }
                        
                        # Parse and validate user input
                        selected_indices <- as.numeric(trimws(unlist(strsplit(user_input, ","))))
                        if (length(selected_indices) != 2 || any(is.na(selected_indices)) || any(selected_indices < 1) || any(selected_indices > length(conditions))) {
                            cat("Invalid input. Please enter exactly two valid indices.\n")
                            next
                        }
                        
                        # Map indices to condition levels
                        selected_conditions <- conditions[selected_indices]
                        
                        # Debugging: Display selected conditions
                        cat("\nDebugging Info: Selected conditions for comparison:\n")
                        print(selected_conditions)
                        
                        # Construct manual contrast
                        contrast_vector <- c("conditions_combined", selected_conditions[1], selected_conditions[2])
                        
                        # Perform comparison
                        cat(sprintf("\nPerforming comparison: %s vs %s\n", selected_conditions[1], selected_conditions[2]))

                        tryCatch({
                            # Use the manual contrast
                            res <- results(dds, contrast = contrast_vector)
                            
                            # Check if `res` contains valid data
                            if (all(is.na(res$padj)) || nrow(res) == 0) {
                                cat(sprintf("No valid results for comparison: %s vs %s. Skipping this comparison.\n",
                                            selected_conditions[1], selected_conditions[2]))
                                next
                            }
                            
                            # Print summary of differentially expressed genes
                            cat("\nSummary of differentially expressed genes:\n")
                            valid_padj <- !is.na(res$padj)
                            cat(sprintf("Number of genes with adjusted p-value < 0.1: %d\n", sum(res$padj < 0.1, na.rm = TRUE)))
                            cat(sprintf("Number of genes with adjusted p-value < 0.05: %d\n", sum(res$padj < 0.05, na.rm = TRUE)))
                            cat(sprintf("Number of genes with adjusted p-value < 0.01: %d\n", sum(res$padj < 0.01, na.rm = TRUE)))
                            cat(sprintf("Number of genes with NA adjusted p-value: %d\n", sum(is.na(res$padj))))
                            
                            # Remove rows where all values are NA
                            res <- res[complete.cases(res), ]
                            if (nrow(res) == 0) {
                                cat(sprintf("No significant genes detected for %s vs %s. Skipping this comparison.\n", 
                                            selected_conditions[2], selected_conditions[1]))
                                next  # Skip adding empty results
                            }
                            
                            # Order the results by adjusted p-value
                            res <- res[order(res$padj, na.last = TRUE), ]

                            # Save DEG results to the specified output directory
                            output_file <- get_results_path(paste0("DEG_results_", selected_conditions[1], "_vs_", selected_conditions[2], ".csv"))
                            write.csv(as.data.frame(res), file = output_file, row.names = TRUE)
                            cat("Results saved to", output_file, "\n")

                            # Add comparison column and append to `degs_data`
                            res$comparison <- paste(selected_conditions[1], "vs", selected_conditions[2])
                            
                            # Include gene names explicitly as a column instead of relying on row names
                            res$Gene <- rownames(res)
                            degs_data <- rbind(degs_data, as.data.frame(res))
                                    
                            # Store gene names
                            gene_names <- unique(c(gene_names, rownames(res)))

                            # Store results in formatted structure
                            col_name <- paste(selected_conditions[1], "vs", selected_conditions[2])
                            mega_results[[col_name]] <- data.frame(
                                Gene = rownames(res),
                                baseMean = res$baseMean,
                                log2FoldChange = res$log2FoldChange,
                                lfcSE = res$lfcSE,
                                stat = res$stat,
                                pvalue = res$pvalue,
                                padj = res$padj
                            )
                        }, error = function(e) {
                                cat("Error generating results for comparison:", selected_conditions[1], "vs", selected_conditions[2], "\n")
                                cat("Detailed Error Message:\n")
                                print(e)
                        })

                        # Prompt to perform more comparisons
                        cat("\nDo you want to perform more comparisons? (y/n, or type 'quit' to exit):\n")
                        more_comparisons <- handle_quit(readline())
                        if (tolower(more_comparisons) == "n" || tolower(more_comparisons) == "quit") {
                            break
                        }
                    }

                    # Save consolidated DEGs data if any comparisons were performed
                    if (nrow(degs_data) > 0) {
                        # Ensure rows are not overwritten by using explicit unique identifiers
                        output_file <- get_results_path("all_DEG_results.csv")
                        write.csv(degs_data, file = output_file, row.names = FALSE)
                        cat(sprintf("Consolidated DEGs data saved as '%s'.\n", output_file))
                    } else {
                        cat("No differential expression results generated.\n")
                    }

                    # Create and save the final data frame
                    if (length(mega_results) > 0) {
                        # Convert list to data frame
                        all_genes <- sort(unique(gene_names))
                        
                        # Create a unified data frame
                        mega_df <- data.frame()

                        # List to track the starting column of each comparison block
                        comparison_col_indices <- c()

                        # Merge all comparisons into one table with missing values filled as NA
                        for (comp_name in names(mega_results)) {
                            df_comp <- mega_results[[comp_name]]
                            
                            # Skip empty results (shouldn't happen, but extra check)
                            if (nrow(df_comp) == 0) {
                                cat(sprintf("Skipping empty comparison: %s\n", comp_name))
                                next
                            }

                            # Ensure df_comp contains all genes in all_genes, filling missing ones with NA
                            aligned_df <- data.frame(Gene = all_genes)  # Create a full gene list
                            df_comp <- merge(aligned_df, df_comp, by = "Gene", all.x = TRUE)  # Fill missing values with NA

                            # Rename the Gene column to include the comparison name
                            colnames(df_comp)[colnames(df_comp) == "Gene"] <- comp_name

                            # colnames(df_comp)[1] <- comp_name  # Rename Gene column to match comparison name
                            df_comp <- df_comp[order(df_comp[[comp_name]]), ]  # Ensure sorting by gene name

                            # Track the start of this comparison block
                            if (nrow(mega_df) == 0) {
                                comparison_col_indices <- c(comparison_col_indices, 1)  # First block starts at column 1
                                mega_df <- df_comp
                            } else {
                                comparison_col_indices <- c(comparison_col_indices, ncol(mega_df) + 2)  # Add 2 to account for the empty column
                                mega_df <- cbind(mega_df, NA, df_comp)  # Insert empty column after each block
                            }
                        }

                        # Insert exactly 7 empty columns before adding counts data
                        for (i in 1:7) {
                            mega_df[[" "]] <- NA
                        }
                        counts_data <- counts(dds, normalized = TRUE)

                        # Extract and align the raw counts data
                        counts_data_final <- counts_data[rownames(counts_data) %in% all_genes, , drop = FALSE]

                        counts_data_final <- counts_data_final[order(rownames(counts_data_final)), ]  # Sort alphabetically
                        counts_data_final <- as.data.frame(counts_data_final)
                        counts_data_final$Gene <- rownames(counts_data_final)
                        # Reorder columns to move Gene to the first position
                        counts_data_final <- counts_data_final[, c("Gene", setdiff(colnames(counts_data_final), "Gene"))]

                        # Append counts data at the right end
                        mega_df <- cbind(mega_df, counts_data_final)


                        # Save as Excel file with pastel colors
                        wb <- createWorkbook()
                        addWorksheet(wb, "DEG Results")

                        # Write data with styling
                        writeData(wb, "DEG Results", mega_df, rowNames = FALSE)

                        light_colors <- c("#FFDDC1", "#C1E1FF", "#C1FFC1", "#FFD1E1", "#E1C1FF", "#FFFFC1")
                        
                        for (i in seq_along(comparison_col_indices)) {
                            col_start <- comparison_col_indices[i]
                            col_end <- col_start + 6  # Assuming 7 columns per comparison block
                            style <- createStyle(fgFill = light_colors[i %% length(light_colors) + 1])
                            addStyle(wb, "DEG Results", style, cols = col_start:col_end, rows = 1:(nrow(mega_df) + 1), gridExpand = TRUE)
                        }

                        # Save the Excel file1
                        output_xlsx <- get_results_path("mega_results.xlsx")
                        # saveWorkbook(wb, output_xlsx, overwrite = TRUE)
                        # cat("Mega results saved as", output_xlsx, "\n")
                        # Step 1: Detect or create a writable temporary directory
                        safe_tmp_dir <- tempdir()

                        if (!dir.exists(safe_tmp_dir)) {
                            cat("⚠️ Temporary directory does not exist. Creating a new one...\n")
                            safe_tmp_dir <- file.path(Sys.getenv("USERPROFILE", unset = Sys.getenv("HOME")), "R_openxlsx_temp")
                            dir.create(safe_tmp_dir, recursive = TRUE, showWarnings = FALSE)
                        }

                        # Step 2: Ensure openxlsx uses the writable temp directory
                        options(openxlsx.tempdir = safe_tmp_dir)

                        # Step 3: Save the workbook safely
                        tryCatch({
                            saveWorkbook(wb, output_xlsx, overwrite = TRUE)
                            cat("✅ Mega results saved as", output_xlsx, "\n")
                        }, error = function(e) {
                            cat("❌ Error saving workbook:", e$message, "\n")
                            
                            # Retry with a completely fresh temp directory
                            retry_tmp_dir <- file.path(tempdir(), paste0("R_openxlsx_retry_", Sys.getpid()))
                            if (!dir.exists(retry_tmp_dir)) {
                                dir.create(retry_tmp_dir, recursive = TRUE, showWarnings = FALSE)
                                options(openxlsx.tempdir = retry_tmp_dir)
                            }

                            tryCatch({
                                saveWorkbook(wb, output_xlsx, overwrite = TRUE)
                                cat("✅ Mega results saved successfully after retrying!\n")
                            }, error = function(e) {
                                cat("❌ Second attempt failed. Please check write permissions.\n")
                            })
                        })
                    } else {
                        cat("No comparisons were performed.\n")
                    }
                    # return to the main menu
                    break
                } else if (sub_menu_choice == 2) {
                    # **Option 2: Perform Differential Expression Analysis Without Processing Count Files**
                
                    repeat {
                        cat("Please provide the path to the count matrix file (or type 'back' to return to the previous step): \n")
                        matrix_path <- handle_quit(readline())
                        if (matrix_path == "back") {
                            cat("Returning to the data input choice.\n")
                            # navigate_back <- TRUE
                            break # Return to the data input choice
                        }

                        if (file.exists(matrix_path)) {
                            tryCatch({
                                cat("Reading the count matrix...\n")
                                if (grepl("\\.xlsx$", matrix_path, ignore.case = TRUE)) {
                                    combined_data <- read_excel(matrix_path, col_names = TRUE)
                                    # Ensure consistency with CSV reading
                                    combined_data <- as.data.frame(combined_data, stringsAsFactors = FALSE)  # Convert to data.frame
                                    colnames(combined_data) <- trimws(colnames(combined_data))               # Trim whitespace from headers
                                    combined_data[] <- lapply(combined_data, function(x) {                   # Trim whitespace from data
                                        if (is.character(x)) trimws(x) else x
                                    })
                                    combined_data[combined_data == ""] <- NA                                # Handle empty cells as NA
                                } else if (grepl("\\.(txt|tsv)$", matrix_path, ignore.case = TRUE)) {
                                    combined_data <- read.delim(matrix_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
                                } else if (grepl("\\.csv$", matrix_path, ignore.case = TRUE)) {
                                    combined_data <- read.csv(matrix_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
                                } else {
                                    cat("Error: Unsupported file format. Please provide a .csv, .tsv, .txt, or .xlsx file.\n")
                                    next # Return to the start of the loop to retry
                                }
                                cat("Count matrix loaded successfully.\n")
                                # navigate_back <- FALSE
                                # break # Exit the loop if the matrix is loaded successfully
                            }, error = function(e) {
                                cat("Error reading the count matrix: ", e$message, "\n")
                                next # Retry the input if there’s an error
                            })
                            break
                        } else {
                            cat("Invalid file path. Please try again.\n")
                        }
                    }                    

                    # Prompt for coldata file path
                    repeat {
                        cat("Please provide the path to the sample info table (or type 'quit' to exit): \n")
                        info_table_path <- handle_quit(readline())
                        # Handle "back" command to return to the previous step
                        if (info_table_path == "back") {
                            cat("Returning to the sample info table choice.\n")
                            final_conditions <- NULL
                            break # Exit to return to the sample info table choice loop
                        }

                        if (file.exists(info_table_path)) {
                            tryCatch({
                                # Function to read sample info table
                                read_sample_info <- function(file_path) {
                                    if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
                                        read.csv(file_path, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
                                    } else if (grepl("\\.(txt|tsv)$", file_path, ignore.case = TRUE)) {
                                        read.delim(file_path, header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = "\t")
                                    } else if (grepl("\\.xlsx$", file_path, ignore.case = TRUE)) {
                                        read_excel(file_path, col_names = TRUE)
                                    } else {
                                        cat("Error: Unsupported file format. Please provide a .csv, .txt, .tsv, or .xlsx file.\n")
                                        return(NULL)
                                    }
                                }

                                # Read the sample info table
                                sample_info <- read_sample_info(info_table_path)
                                # if (is.null(sample_info)) next # Retry if unsupported format

                                # Debugging: Print parsed sample info table
                                cat("Sample info table loaded successfully:\n")
                                print(sample_info)

                                # # Exit the loop on successful read
                                # break
                            }, error = function(e) {
                                cat("Error reading the sample info table: ", e$message, "\n")
                                # next # Retry on error
                            })
                            # Exit the loop on successful read
                            break
                        } else {
                            cat("Invalid file path. Please try again.\n")
                        }
                    }
                    # # Align sample info with count matrix
                    tryCatch({
                        # Align sample info with count matrix
                        coldata <- sample_info
                        # Align with the count matrix columns except the first which is assumed to be geneID
                        coldata <- coldata[match(colnames(combined_data)[-1], rownames(coldata)), , drop = FALSE]

                        if (any(is.na(coldata))) {
                            cat("Error: Mismatch between sample info and count matrix. Please ensure sample names align.\n")
                            next # Retry alignment
                        }

                        # Assign conditions from the sample info table
                        conditions_list <- list()
                        for (i in seq_len(ncol(coldata))) {
                            conditions_list[[paste0("Condition", i)]] <- coldata[, i]
                        }

                        # Combine sample names and all conditions into a final dataframe
                        final_conditions <- data.frame(
                            Sample = colnames(combined_data)[-1],
                            conditions_list,
                            stringsAsFactors = FALSE
                        )

                        colnames(combined_data)[-1] <- final_conditions$Sample  # Update column names
                        rownames(combined_data) <- combined_data$geneID
                        combined_data$geneID <- NULL

                        # Debugging: Print the final data structure
                        cat("\nFinal preview of assigned conditions:\n")
                        print(final_conditions)

                        cat("\nFinal head of data structure with the assigned conditions:\n")
                        print(head(combined_data))
                        # navigate_back <- FALSE
                        # break # Success, exit loop        
                    }, error = function(e) {
                        cat("Error aligning sample info with count matrix: ", e$message, "\n")
                        next # Retry alignment
                    })

                    names(combined_data) <- c("geneID", final_conditions$Sample)
                    row.names(combined_data) <- combined_data$geneID
                    combined_data$geneID <- NULL
                    filtered <- combined_data[!grepl("__", rownames(combined_data)),] #remove the rows with __ in the name

                    # Create a sample info table for DESeq2
                    coldata <- final_conditions[final_conditions$Sample %in% colnames(filtered), ]
                    rownames(coldata) <- coldata$Sample
                    coldata$Sample <- NULL

                    # make combined conditions vector from the conditions list if there is more than 1 condition
                    if (ncol(coldata) >= 1) {
                        if (ncol(coldata) > 1) {
                            # Dynamically identify columns to exclude
                            exclude_cols <- sapply(seq_along(coldata), function(idx) {
                                col <- coldata[[idx]]
                                # Extract all other columns
                                other_columns <- coldata[, -idx, drop = FALSE]
                                # Generate all possible combinations of unique words from other columns
                                other_combinations <- apply(other_columns, 1, paste, collapse = "_")
                                # Check if more than half of the column values match the combinations
                                mean(col %in% other_combinations) > 0.5
                            })
                            # Select columns to combine (exclude the identified ones)
                            columns_to_combine <- coldata[, !exclude_cols, drop = FALSE]
                            # Avoid duplicating "conditions_combined"
                            if ("conditions_combined" %in% colnames(coldata) || ncol(columns_to_combine) < ncol(coldata)) {
                                conditions_combined <- coldata[, ncol(coldata)]
                                # conditions_combined <- apply(coldata[, 1:(ncol(coldata) - 1), drop = FALSE], 1, paste, collapse = "_")          
                                # conditions_combined <- apply(coldata[, -ncol(coldata), drop = FALSE], 1, paste, collapse = "_")
                            } else {
                                conditions_combined <- apply(columns_to_combine, 1, paste, collapse = "_")
                            }
                        } else {
                            # Handle single-column case
                            conditions_combined <- coldata[, 1]
                        }
                        conditions_combined <- factor(conditions_combined)
                        coldata$conditions_combined <- conditions_combined
                    } else {
                        cat("Error: coldata has no valid columns remaining. Please check the sample info table and try again.\n")
                        next
                    }
                    # User option: choose which column (or combined) to use for PCA, RLE, and DESeq design
                    if (ncol(coldata) > 1) {
                        design_cols <- setdiff(colnames(coldata), "conditions_combined")
                        if (length(design_cols) == 0) design_cols <- colnames(coldata)[1]
                        cat("\nAvailable columns in sample info (for grouping in PCA, RLE, and comparisons):\n")
                        for (i in seq_along(design_cols)) cat(sprintf("  [%d] %s\n", i, design_cols[i]))
                        cat(sprintf("  [%d] Combined (paste all columns)\n", length(design_cols) + 1))
                        repeat {
                            choice_input <- handle_quit(readline("Enter the number of your choice (or 'back' to keep current): "))
                            if (choice_input == "back") break
                            idx <- as.numeric(choice_input)
                            if (!is.na(idx) && idx >= 1 && idx <= length(design_cols) + 1) {
                                if (idx <= length(design_cols)) {
                                    conditions_combined <- coldata[[design_cols[idx]]]
                                } else {
                                    conditions_combined <- apply(coldata[, design_cols, drop = FALSE], 1, paste, collapse = "_")
                                }
                                conditions_combined <- factor(conditions_combined)
                                coldata$conditions_combined <- conditions_combined
                                cat("Using", if (idx <= length(design_cols)) design_cols[idx] else "Combined", "for grouping.\n")
                                break
                            }
                            cat("Invalid choice. Please enter a number between 1 and", length(design_cols) + 1, ".\n")
                        }
                    }
                    conditions_combined <- factor(conditions_combined)
                
                    # Order coldata by the combined conditions and reorder the filtered data accordingly
                    coldata <- coldata[order(coldata$conditions_combined, decreasing = TRUE), ]
                    filtered <- filtered[, rownames(coldata)]

                    # Clean up the filtered data
                    colnames(filtered) <- make.unique(colnames(filtered))
                    filtered <- as.matrix(filtered)
                    filtered[is.na(filtered)] <- 0
                    filtered[is.nan(filtered)] <- 0
                    filtered[is.infinite(filtered)] <- 0

                    conditions_combined <- coldata$conditions_combined

                    i <- 1
                    while (i < ncol(coldata)) {  # Use while to dynamically adjust column count
                        for (j in (i + 1):ncol(coldata)) {
                            if (all(coldata[, i] == coldata[, j])) {  # Check if column values are identical
                                cat(sprintf("The column %s is the same as the column %s. Removing the first duplicated column.\n",
                                            colnames(coldata)[i], colnames(coldata)[j]))
                                coldata <- coldata[, -i, drop = FALSE]  # Remove first duplicated column
                                i <- i - 1  # Adjust index to stay at the correct position
                                break  # Restart comparison for the new structure
                            }
                        }
                        i <- i + 1  # Move to the next column
                    }

                    # Print the coldata for debugging
                    cat("\nSample info table for DESeq2:\n")
                    print(coldata)

                    # Create DESeq2 dataset
                    cat("Creating DESeq2 dataset...\n")
                    tryCatch({
                        # dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = coldata, design = ~ conditions_combined)
                        dds <- DESeqDataSetFromMatrix(countData = as.matrix(filtered), colData = coldata, design = ~ conditions_combined)
                        cat("DESeq2 dataset created successfully.\n")
                    }, error = function(e) {
                        cat("Error creating DESeq2 dataset: ", e$message, "\n")
                        next
                    })

                    # Perform DESeq analysis
                    cat("Running DESeq analysis...\n")
                    tryCatch({
                        dds <- DESeq(dds)
                        cat("DESeq analysis completed successfully.\n")
                    }, error = function(e) {
                        cat("Error during DESeq analysis: ", e$message, "\n")
                        next
                    })

                    # write the normalized counts
                    normalized_counts <- counts(dds, normalized = TRUE)
                    normalized_counts_file <- get_results_path("normalized_dds_counts.csv")
                    write.csv(normalized_counts, file = normalized_counts_file)

                    # Initialize storage for formatted results
                    mega_results <- list()
                    gene_names <- c()
                    counts_data_final <- NULL  # To store counts data for the final section
                    degs_data <- data.frame()  # Initialize an empty data frame for results

                    repeat {
                        cat("\n📊 Input Format:\n")
                        cat("Enter two condition indices as: X,Y\n")
                        cat("- X = Treatment/Condition of Interest\n")
                        cat("- Y = Reference/Control Condition\n")
                        cat("log2FoldChange shows changes in X relative to Y:\n")
                        cat("  - Positive log2FC → Upregulated in X vs Y\n")
                        cat("  - Negative log2FC → Downregulated in X vs Y\n\n")
                        cat("--------------------------------------------------\n")

                        cat("\nAvailable conditions for comparisons:\n")
                        
                        # Display unique conditions with indices
                        conditions <- levels(conditions_combined)
                        for (i in seq_along(conditions)) {
                            cat(sprintf("[%d] %s\n", i, conditions[i]))
                        }
                        
                        # Ask user to select conditions for comparisons
                        cat("\nPlease specify the indices of the conditions you want to compare (comma-separated, or type 'quit' to exit):\n")
                        user_input <- handle_quit(readline())
                        
                        if (tolower(user_input) == "quit") {
                            break
                        }
                        
                        # Parse and validate user input
                        selected_indices <- as.numeric(unlist(strsplit(user_input, ",")))
                        if (length(selected_indices) != 2 || any(is.na(selected_indices)) || any(selected_indices < 1) || any(selected_indices > length(conditions))) {
                            cat("Invalid input. Please enter exactly two valid indices.\n")
                            next
                        }
                        
                        # Map indices to condition levels
                        selected_conditions <- conditions[selected_indices]
                        
                        # Perform comparison
                        cat(sprintf("\nPerforming comparison: %s vs %s\n", selected_conditions[1], selected_conditions[2]))
                        contrast_vector <- c("conditions_combined", selected_conditions[1], selected_conditions[2])

                         tryCatch({
                            # Use the manual contrast
                            res <- results(dds, contrast = contrast_vector)
                            
                            # Check if `res` contains valid data
                            if (all(is.na(res$padj)) || nrow(res) == 0) {
                                cat(sprintf("No valid results for comparison: %s vs %s. Skipping this comparison.\n",
                                            selected_conditions[1], selected_conditions[2]))
                                next
                            }
                            
                            # Print summary of differentially expressed genes
                            cat("\nSummary of differentially expressed genes:\n")
                            valid_padj <- !is.na(res$padj)
                            cat(sprintf("Number of genes with adjusted p-value < 0.1: %d\n", sum(res$padj < 0.1, na.rm = TRUE)))
                            cat(sprintf("Number of genes with adjusted p-value < 0.05: %d\n", sum(res$padj < 0.05, na.rm = TRUE)))
                            cat(sprintf("Number of genes with adjusted p-value < 0.01: %d\n", sum(res$padj < 0.01, na.rm = TRUE)))
                            cat(sprintf("Number of genes with NA adjusted p-value: %d\n", sum(is.na(res$padj))))
                            
                            # Remove rows where all values are NA
                            res <- res[complete.cases(res), ]
                            if (nrow(res) == 0) {
                                cat(sprintf("No significant genes detected for %s vs %s. Skipping this comparison.\n", 
                                            selected_conditions[2], selected_conditions[1]))
                                next  # Skip adding empty results
                            }
                            
                            # Order the results by adjusted p-value
                            res <- res[order(res$padj, na.last = TRUE), ]

                            # Save DEG results to the specified output directory
                            output_file <- get_results_path(paste0("DEG_results_", selected_conditions[1], "_vs_", selected_conditions[2], ".csv"))
                            write.csv(as.data.frame(res), file = output_file, row.names = TRUE)
                            cat("Results saved to", output_file, "\n")

                            # Add comparison column and append to `degs_data`
                            res$comparison <- paste(selected_conditions[1], "vs", selected_conditions[2])
                            
                            # Include gene names explicitly as a column instead of relying on row names
                            res$Gene <- rownames(res)
                            degs_data <- rbind(degs_data, as.data.frame(res))
                                    
                            # Store gene names
                            gene_names <- unique(c(gene_names, rownames(res)))

                            # Store results in formatted structure
                            col_name <- paste(selected_conditions[1], "vs", selected_conditions[2])
                            mega_results[[col_name]] <- data.frame(
                                Gene = rownames(res),
                                baseMean = res$baseMean,
                                log2FoldChange = res$log2FoldChange,
                                lfcSE = res$lfcSE,
                                stat = res$stat,
                                pvalue = res$pvalue,
                                padj = res$padj
                            )
                        }, error = function(e) {
                                cat("Error generating results for comparison:", selected_conditions[1], "vs", selected_conditions[2], "\n")
                                cat("Detailed Error Message:\n")
                                print(e)
                        })

                        # Prompt for more comparisons
                        cat("\nDo you want to perform more comparisons? (y/n, or type 'quit' to exit):\n")
                        more_comparisons <- handle_quit(readline())
                        if (tolower(more_comparisons) == "n" || tolower(more_comparisons) == "quit") {
                            break
                        }
                    }

                    # Save consolidated results
                    if (nrow(degs_data) > 0) {
                        # Construct the full path using the user-defined output directory
                        output_file <- get_results_path("all_DEG_results.csv")
                        # Write the consolidated DEGs data to the file
                        write.csv(degs_data, file = output_file, row.names = FALSE)
                        # Inform the user where the file was saved
                        cat("Consolidated DEGs data saved to:", output_file, "\n")
                    } else {
                        cat("No differential expression results generated.\n")
                    }

                    # Create and save the final data frame
                    if (length(mega_results) > 0) {
                        # Convert list to data frame
                        all_genes <- sort(unique(gene_names))
                        
                        # Create a unified data frame
                        mega_df <- data.frame()

                        # List to track the starting column of each comparison block
                        comparison_col_indices <- c()

                        # Merge all comparisons into one table with missing values filled as NA
                        for (comp_name in names(mega_results)) {
                            df_comp <- mega_results[[comp_name]]
                            
                            # Skip empty results (shouldn't happen, but extra check)
                            if (nrow(df_comp) == 0) {
                                cat(sprintf("Skipping empty comparison: %s\n", comp_name))
                                next
                            }

                            # Ensure df_comp contains all genes in all_genes, filling missing ones with NA
                            aligned_df <- data.frame(Gene = all_genes)  # Create a full gene list
                            df_comp <- merge(aligned_df, df_comp, by = "Gene", all.x = TRUE)  # Fill missing values with NA

                            # Rename the Gene column to include the comparison name
                            colnames(df_comp)[colnames(df_comp) == "Gene"] <- comp_name

                            # colnames(df_comp)[1] <- comp_name  # Rename Gene column to match comparison name
                            df_comp <- df_comp[order(df_comp[[comp_name]]), ]  # Ensure sorting by gene name

                            # Track the start of this comparison block
                            if (nrow(mega_df) == 0) {
                                comparison_col_indices <- c(comparison_col_indices, 1)  # First block starts at column 1
                                mega_df <- df_comp
                            } else {
                                comparison_col_indices <- c(comparison_col_indices, ncol(mega_df) + 2)  # Add 2 to account for the empty column
                                mega_df <- cbind(mega_df, NA, df_comp)  # Insert empty column after each block
                            }
                        }

                        # Insert exactly 7 empty columns before adding counts data
                        for (i in 1:7) {
                            mega_df[[" "]] <- NA
                        }
                        counts_data <- counts(dds, normalized = TRUE)

                        # Extract and align the raw counts data
                        counts_data_final <- counts_data[rownames(counts_data) %in% all_genes, , drop = FALSE]

                        counts_data_final <- counts_data_final[order(rownames(counts_data_final)), ]  # Sort alphabetically
                        counts_data_final <- as.data.frame(counts_data_final)
                        counts_data_final$Gene <- rownames(counts_data_final)
                        # Reorder columns to move Gene to the first position
                        counts_data_final <- counts_data_final[, c("Gene", setdiff(colnames(counts_data_final), "Gene"))]

                        # Append counts data at the right end
                        mega_df <- cbind(mega_df, counts_data_final)


                        # Save as Excel file with pastel colors
                        wb <- createWorkbook()
                        addWorksheet(wb, "DEG Results")

                        # Write data with styling
                        writeData(wb, "DEG Results", mega_df, rowNames = FALSE)

                        light_colors <- c("#FFDDC1", "#C1E1FF", "#C1FFC1", "#FFD1E1", "#E1C1FF", "#FFFFC1")
                        
                        for (i in seq_along(comparison_col_indices)) {
                            col_start <- comparison_col_indices[i]
                            col_end <- col_start + 6  # Assuming 7 columns per comparison block
                            style <- createStyle(fgFill = light_colors[i %% length(light_colors) + 1])
                            addStyle(wb, "DEG Results", style, cols = col_start:col_end, rows = 1:(nrow(mega_df) + 1), gridExpand = TRUE)
                        }

                        # Save the Excel file1
                        output_xlsx <- get_results_path("mega_results.xlsx")
                        # saveWorkbook(wb, output_xlsx, overwrite = TRUE)
                        # cat("Mega results saved as", output_xlsx, "\n")
                        # Step 1: Detect or create a writable temporary directory
                        safe_tmp_dir <- tempdir()

                        if (!dir.exists(safe_tmp_dir)) {
                            cat("⚠️ Temporary directory does not exist. Creating a new one...\n")
                            safe_tmp_dir <- file.path(Sys.getenv("USERPROFILE", unset = Sys.getenv("HOME")), "R_openxlsx_temp")
                            dir.create(safe_tmp_dir, recursive = TRUE, showWarnings = FALSE)
                        }

                        # Step 2: Ensure openxlsx uses the writable temp directory
                        options(openxlsx.tempdir = safe_tmp_dir)

                        # Step 3: Save the workbook safely
                        tryCatch({
                            saveWorkbook(wb, output_xlsx, overwrite = TRUE)
                            cat("✅ Mega results saved as", output_xlsx, "\n")
                        }, error = function(e) {
                            cat("❌ Error saving workbook:", e$message, "\n")
                            
                            # Retry with a completely fresh temp directory
                            retry_tmp_dir <- file.path(tempdir(), paste0("R_openxlsx_retry_", Sys.getpid()))
                            if (!dir.exists(retry_tmp_dir)) {
                                dir.create(retry_tmp_dir, recursive = TRUE, showWarnings = FALSE)
                                options(openxlsx.tempdir = retry_tmp_dir)
                            }

                            tryCatch({
                                saveWorkbook(wb, output_xlsx, overwrite = TRUE)
                                cat("✅ Mega results saved successfully after retrying!\n")
                            }, error = function(e) {
                                cat("❌ Second attempt failed. Please check write permissions.\n")
                            })
                        })
                    } else {
                        cat("No comparisons were performed.\n")
                    }
                } else {
                    cat("Invalid choice. Please enter a valid number between 1 and 2.\n")
                }
            }
        }
        if (menu_choice == 3) {
            # **Option 3: Load DEGs & genomic results from Files**
            repeat {
                # **Option 1: Load DEGs Tables**
                cat("Please provide the path to the genomic data file (or type 'quit' to exit or type 'back' to return to the main menu):\n")
                degs_file_path <- handle_quit(readline())

                if (degs_file_path == "back") {
                    cat("Returning to the main menu.\n")
                    break # Return to the main menu
                }    

                if (!file.exists(degs_file_path)) {
                    cat("Error: The DEGs file does not exist. Returning to the menu.\n")
                    next
                }
                tryCatch({
                    # Determine the file type based on extension
                    if (grepl("\\.csv$", degs_file_path, ignore.case = TRUE)) {
                        degs_data <- read.csv(degs_file_path, header = TRUE, stringsAsFactors = FALSE)
                    } else if (grepl("\\.txt$|\\.tsv$", degs_file_path, ignore.case = TRUE)) {
                        degs_data <- read.delim(degs_file_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t") 
                    } else if (grepl("\\.xlsx$", degs_file_path, ignore.case = TRUE)) {
                        degs_data <- read_excel(degs_file_path, col_names = TRUE)
                        # Ensure consistency with CSV reading
                        degs_data <- as.data.frame(degs_data, stringsAsFactors = FALSE)  # Convert to data.frame
                        colnames(degs_data) <- trimws(colnames(degs_data))               # Trim whitespace from headers
                        degs_data[] <- lapply(degs_data, function(x) {                   # Trim whitespace from data
                            if (is.character(x)) trimws(x) else x
                        })
                        degs_data[degs_data == ""] <- NA                                # Handle empty cells as NA
                    } else {
                        cat("Unsupported file format. Please provide a file in CSV, TXT, TSV, or XLSX format.\n")
                        next # Retry the input
                    }
                    
                    # Ensure the structure of the loaded DEGs data
                    required_columns <- c("log2FoldChange", "padj", "comparison", "Gene")
                    missing_columns <- setdiff(required_columns, colnames(degs_data))
                    
                    if (length(missing_columns) > 0) {
                        cat(sprintf("Error: The DEGs file is missing the following required columns: %s\n", paste(missing_columns, collapse = ", ")))
                        next # Retry the input
                    }
                    
                    # Check for duplicate rows (optional informational message)
                    if (anyDuplicated(degs_data)) {
                        cat("Warning: Duplicate rows detected in the DEGs file. Proceeding with all entries.\n")
                    }
                    
                    cat("DEGs table loaded successfully:\n")
                    print(head(degs_data))
                    break # Exit loop after successful load
                }, error = function(e) {
                    cat("Error reading the DEGs file: ", e$message, "\nReturning to the menu.\n")
                    next
                })
            } # repeat for loading DEGs and genomic file
            repeat {
                cat("\nInteractive Analysis Menu\n")
                cat("[1] Check if there are some genes of interest\n")
                cat("[2] Generate a heatmap\n")
                cat("[3] Create Venn diagrams\n")
                cat("[4] Quit\n")
                cat("Please select an option by entering the corresponding number (or type 'back' to return to main menu):\n")
                
                user_choice_raw <- handle_quit(readline())
                if (user_choice_raw == "back") break
                user_choice <- as.numeric(user_choice_raw)
                
                if (is.na(user_choice) || user_choice < 1 || user_choice > 4) {
                    cat("Invalid choice. Please enter a valid number between 1 and 4.\n")
                    next
                }
                
                if (user_choice == 4) {
                    cat("Exiting the program. Goodbye!\n")
                    break
                }
                
                # **Analysis Options**
                if (user_choice == 1) {
                    # **Option 1: Check Genes of Interest**
                    repeat {
                        cat("Insert the gene names (comma-separated, or provide a file path, or type 'quit' to exit):\n")
                        genes_input <- handle_quit(readline())
                        # Check if the input is a valid file path
                        if (file.exists(genes_input)) {
                            cat("Detected file input. Extracting genes from the first column...\n")
                            if (grepl("\\.xlsx$", genes_input, ignore.case = TRUE)) {
                                genes_df <- read_excel(genes_input, col_names = TRUE)
                                # Ensure consistency with CSV reading
                            } else if (grepl("\\.(txt|tsv)$", matrix_path, ignore.case = TRUE)) {
                                genes_df <- read.table(genes_input, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
                            } else if (grepl("\\.csv$", matrix_path, ignore.case = TRUE)) {
                                genes_df <- read.csv(genes_input, header = TRUE, stringsAsFactors = FALSE)
                            } else {
                                cat("Error: Unsupported file format. Please provide a .csv, .tsv, .txt, or .xlsx file.\n")
                                next # Return to the start of the loop to retry
                            }

                            if (ncol(genes_df) < 1) {
                                cat("Error: The file does not contain any columns. Please try again.\n")
                                next
                            }
                            
                            genes_in_interest <- unique(na.omit(genes_df[[1]]))  # Extract first column and remove NAs
                        } else {
                            # Process manual input (comma-separated values)
                            genes_in_interest <- unlist(strsplit(gsub("\\s+", "", genes_input), ","))
                        }                   
                        if (length(genes_in_interest) == 0) {
                            cat("No valid genes entered. Please try again.\n")
                        } else {
                            break
                        }
                    }

                    cat("Searching for genes in the DEGs data ...\n")


                    # Initialize a results list to store summary information
                    results_list <- list()

                    # Iterate over each gene of interest
                    for (gene in genes_in_interest) {
                        # Match all rows where the `Gene` column contains the gene of interest correctly
                        matched_rows <- degs_data[grepl(paste0("^", gene, "$|^", gene, "\\."), degs_data$Gene, ignore.case = TRUE), , drop = FALSE]

                        if (nrow(matched_rows) > 0) {
                            cat(sprintf("Gene %s found in %d comparisons:\n", gene, nrow(matched_rows)))

                            for (i in seq_len(nrow(matched_rows))) {
                                row <- matched_rows[i, ]
                                cat(sprintf("  Comparison: %s, padj = %.4f, log2FoldChange = %.4f\n",
                                            row["comparison"], as.numeric(row["padj"]), as.numeric(row["log2FoldChange"])))
                            }

                            # Store all matches
                            results_list[[gene]] <- matched_rows
                            # Plot log2FC by comparison (genomic line plot)
                            plot_gene_expression_from_genomic(gene, matched_rows)

                        } else {
                            cat(sprintf("Gene %s not found in the DEGs data.\n", gene))
                            # Create a placeholder row for unmatched genes
                            placeholder_row <- data.frame(
                                Gene = gene,
                                log2FoldChange = NA,
                                padj = NA,
                                comparison = NA,
                                stringsAsFactors = FALSE
                            )
                            results_list[[gene]] <- placeholder_row
                        }
                    }

                    # Combine all results into a single summary table
                    summary_table <- do.call(rbind, results_list)

                    # Ensure numeric consistency for columns
                    numeric_columns <- c("log2FoldChange", "padj")
                    summary_table[numeric_columns] <- lapply(summary_table[numeric_columns], as.numeric)

                    # Save the summary table under genomic folder
                    output_file <- get_genomic_genes_of_interest_path("gene_search_summary.csv")
                    write.csv(summary_table, file = output_file, row.names = FALSE)
                    cat("Summary table saved to:", output_file, "\n")
                }
                if (user_choice == 2) {
                    # **Option 2: Generate Heatmap**
                    # 1) Prompt for padj threshold
                    repeat {
                        padj_input <- handle_quit(
                            readline(prompt = "Enter the padj threshold (e.g., 0.05), or type 'back' to return to menu:\n")
                        )
                        if (padj_input == "back") {
                            cat("Returning to main menu.\n")
                            break  # exit the heatmap block
                        }
                        padj_cutoff <- suppressWarnings(as.numeric(padj_input))
                        if (!is.na(padj_cutoff) && padj_cutoff >= 0 && padj_cutoff <= 1) {
                            break
                        }
                        cat("Invalid padj. Please enter a number between 0 and 1.\n")
                    }
                    if (!exists("padj_cutoff")) {
                        # User quit
                        break
                    }

                    # 2) Prompt for absolute log2FoldChange threshold
                    repeat {
                        log2fc_input <- handle_quit(
                            readline(prompt = "Enter the absolute log2FoldChange threshold (e.g., 1 or 2), or type 'back' to return to menu:\n")
                        )
                        if (log2fc_input == "back") {
                            cat("Returning to main menu.\n")
                            break
                        }
                        log2fc_cutoff <- suppressWarnings(as.numeric(log2fc_input))
                        if (!is.na(log2fc_cutoff) && log2fc_cutoff >= 0) {
                            break
                        }
                        cat("Invalid log2FoldChange. Please enter a non-negative number.\n")
                    }
                    if (!exists("log2fc_cutoff")) {
                        break
                    }
                    # if str_detect(comparison, "rna")
                    if (degs_data %>% 
                        filter(str_detect(comparison, "RNA")) %>% nrow() > 0) {
                        # 3) Split DEGs into RNA vs genomic
                        rna_df  <- degs_data %>% filter(str_detect(comparison, "RNA"))
                        genomic_df <- degs_data %>% filter(!str_detect(comparison, "RNA"))
                        cat("Detected comparisons:\n")
                        rna_comparisons <- unique(rna_df$comparison)
                        genomic_comparisons <- unique(genomic_df$comparison)
                        cat("RNA comparisons:\n")
                        for (cmp in rna_comparisons) {
                            cat("  -", cmp, "\n")
                        }
                        cat("Genomic comparisons:\n")
                        for (cmp in genomic_comparisons) {
                            cat("  -", cmp, "\n")
                        }
                    } else {
                        # 3) Split DEGs into RNA vs genomic
                        rna_df  <- degs_data %>% filter(str_detect(comparison, "vs"))
                        genomic_df <- degs_data %>% filter(!str_detect(comparison, "vs"))
                    }
                    # 4) Filter by padj < cutoff and |log2FoldChange| ≥ cutoff
                    # rna_sig  <- rna_df  %>% filter(padj < padj_cutoff & abs(log2FoldChange) >= log2fc_cutoff)
                    rna_sig <- rna_df %>%
                        filter(padj < padj_cutoff & abs(log2FoldChange) >= log2fc_cutoff) %>%
                        arrange(desc(abs(log2FoldChange))) %>%
                        distinct(Gene, comparison, .keep_all = TRUE)


                    # genomic_sig <- genomic_df %>% filter(padj < padj_cutoff & abs(log2FoldChange) >= log2fc_cutoff)

                    # 5) In genomic, if a (Gene, comparison) appears multiple times, keep the row with largest |log2FoldChange|
                    # genomic_sig <- genomic_sig %>%
                    # group_by(Gene, comparison) %>%
                    # slice_max(order_by = abs(log2FoldChange), n = 1) %>%
                    # ungroup()
                    genomic_sig <- genomic_df %>%
                        filter(padj < padj_cutoff & abs(log2FoldChange) >= log2fc_cutoff) %>%
                        arrange(desc(abs(log2FoldChange))) %>%
                        distinct(Gene, comparison, .keep_all = TRUE)


                    # # 6) Remove exact duplicates for RNA (just keep the first if any)
                    # rna_sig  <- rna_sig  %>% distinct(Gene, comparison, .keep_all = TRUE)

                    # 7) Pivot RNA to wide
                    rna_wide_tib <- rna_sig %>%
                        select(Gene, comparison, log2FoldChange) %>%
                        pivot_wider(names_from = comparison, values_from = log2FoldChange, values_fill = NA)
                    if (nrow(rna_wide_tib) == 0) {
                        cat("No RNA-seq genes passed thresholds. Returning to main menu.\n")
                        break
                    }
                    rna_wide_df <- as.data.frame(rna_wide_tib)
                    rownames(rna_wide_df) <- rna_wide_df$Gene
                    rna_wide_df$Gene <- NULL
                    rna_mat_all <- as.matrix(rna_wide_df)

                    # 8) Pivot genomic to wide
                    genomic_mat_all <- matrix(nrow = 0, ncol = 0)
                    if (nrow(genomic_sig) > 0) {
                        genomic_wide_tib <- genomic_sig %>%
                            select(Gene, comparison, log2FoldChange) %>%
                            pivot_wider(names_from = comparison, values_from = log2FoldChange, values_fill = NA)
                        if (nrow(genomic_wide_tib) > 0) {
                            genomic_wide_df <- as.data.frame(genomic_wide_tib)
                            rownames(genomic_wide_df) <- genomic_wide_df$Gene
                            genomic_wide_df$Gene <- NULL
                            genomic_mat_all <- as.matrix(genomic_wide_df)
                        }
                    }
                    if (ncol(genomic_mat_all) == 0) {
                        cat("No Genomic comparisons passed thresholds. Returning to main menu.\n")
                        break
                    }
                    # Only one genomic comparison assumed; pipeline supports single RNA and single genomic (e.g. example RNAseq + RRBS).
                    genomic_cmp <- colnames(genomic_mat_all)[1]
                    if (ncol(rna_mat_all) == 1L && ncol(genomic_mat_all) == 1L) {
                        cat("Using 1 RNA comparison and 1 genomic comparison.\n")
                    }

                    # 9) Prompt for Default vs Customized mode (no UNION/INTERSECTION choice)
                    repeat {
                        cat("\nChoose heatmap mode:\n")
                        cat("[1] Default  (all 4 direction combinations, per-comparison & combined)\n")
                        cat("[2] Customized  (select genomic data direction(s) & RNA comparisons/directions)\n")
                        hm_input <- handle_quit(
                            readline(prompt = "Enter 1 or 2, or 'back' to return to menu:\n")
                        )
                        if (hm_input == "back") {
                            cat("Returning to main menu.\n")
                            break
                        }
                        hm_type <- suppressWarnings(as.integer(hm_input))
                        if (!is.na(hm_type) && hm_type %in% c(1, 2)) {
                            break
                        }
                        cat("Invalid choice. Please enter 1 or 2.\n")
                    }
                    if (!exists("hm_type")) {
                        break
                    }

                    # Helper: Filter gene names of a numeric vector by direction
                    filter_genes_by_direction <- function(vec, direction) {
                        if (direction == "up") {
                            return(names(vec)[which(vec >= 0 & !is.na(vec))])
                        }
                        if (direction == "down") {
                            return(names(vec)[which(vec <= 0 & !is.na(vec))])
                        }
                        return(character(0))
                    }

                    draw_and_save <- function(mat, title_text, base_name) {
                    if (is.null(mat) || nrow(mat) == 0) {
                        cat(paste0("Skipping \"", title_text, "\" (no genes remain).\n"))
                        return()
                    }
                    mat[is.na(mat)] <- 0
                    my_palette <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
                    ht <- Heatmap(
                        mat,
                        name                = "log2FC",
                        heatmap_legend_param= list(
                            title    = "log2FC",
                            title_gp = gpar(fontfamily = "Arial", fontsize = 14),
                            labels_gp= gpar(fontfamily = "Arial", fontsize = 12)
                        ),
                        col                 = my_palette,
                        cluster_rows        = TRUE,
                        cluster_columns     = TRUE,
                        show_row_dend       = FALSE,
                        show_column_dend    = FALSE,
                        show_row_names      = FALSE,
                        row_names_gp        = gpar(fontfamily = "Arial", fontsize = 10),
                        show_column_names   = TRUE,
                        column_names_rot    = 45,
                        column_names_gp     = gpar(fontfamily = "Arial", fontsize = 10),
                        column_title        = title_text,
                        column_title_gp     = gpar(fontfamily = "Arial", fontsize = 16, fontface = "bold"),
                        rect_gp             = gpar(col = "black", lwd = 1),
                        width               = unit(ncol(mat) * 6, "mm")
                    )
                    png_path <- get_genomic_heatmap_path(paste0(base_name, ".png"))
                    png(filename = png_path, width = 1400, height = 1000, res = 150, family = "Arial")
                    draw(ht)
                    dev.off()
                    cat("Saved heatmap PNG: ", png_path, "\n")
                    svg_path <- get_genomic_heatmap_path(paste0(base_name, ".svg"))
                    grDevices::svg(svg_path, width = 1400/150, height = 1000/150, family = "Arial")
                    draw(ht)
                    dev.off()
                    cat("Saved heatmap SVG: ", svg_path, "\n")
                    csv_path <- get_genomic_heatmap_path(paste0(base_name, ".csv"))
                    out_df <- as.data.frame(mat)
                    out_df <- cbind(Gene = rownames(out_df), out_df)
                    write.csv(out_df, file = csv_path, row.names = FALSE)
                    cat("Saved data CSV:  ", csv_path, "\n")
                    }

                    #--------------------------------------
                    # 10) DEFAULT MODE: always use INTERSECTION logic
                    #--------------------------------------
                    if (hm_type == 1) {
                        cat("\n--- DEFAULT MODE: Generating 4 direction combinations ---\n")
                        
                        rna_cols <- colnames(rna_mat_all)
                        
                        # Per-comparison heatmaps
                        for (cmp in rna_cols) {
                            cat(sprintf("\n--- RNA comparison: %s ---\n", cmp))
                            rna_vec  <- rna_mat_all[, cmp]
                            genomic_vec <- genomic_mat_all[, genomic_cmp]
                            
                            rna_up    <- filter_genes_by_direction(rna_vec,   "up")
                            rna_down  <- filter_genes_by_direction(rna_vec, "down")
                            genomic_up   <- filter_genes_by_direction(genomic_vec,  "up")
                            genomic_down <- filter_genes_by_direction(genomic_vec,"down")
                            
                            # 1) genomic_up vs RNA_up (intersection)
                            genes_uu <- intersect(genomic_up, rna_up)
                            if (length(genes_uu) > 0) {
                                # mat_uu <- cbind(rna_val = rna_vec[genes_uu], genomic_val = genomic_vec[genes_uu])
                                valid_genes <- genes_uu[!is.na(rna_vec[genes_uu]) & !is.na(genomic_vec[genes_uu])]
                                mat_uu <- cbind(rna_val = rna_vec[valid_genes], genomic_val = genomic_vec[valid_genes])

                                colnames(mat_uu) <- c(cmp, genomic_cmp)
                                rownames(mat_uu) <- genes_uu
                                draw_and_save(
                                    mat = mat_uu,
                                    title_text = paste0("Hyper: ", cmp),
                                    base_name  = paste0(gsub("[^A-Za-z0-9]", "_", cmp), "_hyper")
                                )
                            } else {
                                cat(sprintf("Skipping Hyper for %s (no genes).\n", cmp))
                            }
                            
                            # 2) genomic_up vs RNA_down (intersection)
                            genes_ud <- intersect(genomic_up, rna_down)
                            if (length(genes_ud) > 0) {
                                mat_ud <- cbind(rna_val = rna_vec[genes_ud], genomic_val = genomic_vec[genes_ud])
                                colnames(mat_ud) <- c(cmp, genomic_cmp)
                                rownames(mat_ud) <- genes_ud
                                draw_and_save(
                                    mat = mat_ud,
                                    title_text = paste0("genomic_up vs RNA_down: ", cmp),
                                    base_name  = paste0(gsub("[^A-Za-z0-9]", "_", cmp), "_genomicup_vs_RNAdown")
                                )
                            } else {
                                cat(sprintf("Skipping genomic_up vs RNA_down for %s (no genes).\n", cmp))
                            }
                            
                            # 3) genomic_down vs RNA_up (intersection)
                            genes_du <- intersect(genomic_down, rna_up)
                            if (length(genes_du) > 0) {
                                mat_du <- cbind(rna_val = rna_vec[genes_du], genomic_val = genomic_vec[genes_du])
                                colnames(mat_du) <- c(cmp, genomic_cmp)
                                rownames(mat_du) <- genes_du
                                draw_and_save(
                                    mat = mat_du,
                                    title_text = paste0("genomic_down vs RNA_up: ", cmp),
                                    base_name  = paste0(gsub("[^A-Za-z0-9]", "_", cmp), "_genomicdown_vs_RNAup")
                                )
                            } else {
                                cat(sprintf("Skipping genomic_down vs RNA_up for %s (no genes).\n", cmp))
                            }
                            
                            # 4) genomic_down vs RNA_down (intersection)
                            genes_dd <- intersect(genomic_down, rna_down)
                            if (length(genes_dd) > 0) {
                                mat_dd <- cbind(rna_val = rna_vec[genes_dd], genomic_val = genomic_vec[genes_dd])
                                colnames(mat_dd) <- c(cmp, genomic_cmp)
                                rownames(mat_dd) <- genes_dd
                                draw_and_save(
                                    mat = mat_dd,
                                    title_text = paste0("Hypo: ", cmp),
                                    base_name  = paste0(gsub("[^A-Za-z0-9]", "_", cmp), "_hypo")
                                )
                            } else {
                                cat(sprintf("Skipping genomic_down vs RNA_down for %s (no genes).\n", cmp))
                            }
                        }
                        
                        # Combined across *all* RNA comparisons (intersection)
                        cat("\n--- Combined across all RNA comparisons ---\n")
                        for (dir_pair in c("up_up", "up_down", "down_up", "down_down")) {
                            parts    <- str_split(dir_pair, "_")[[1]]
genomic_dir <- parts[1]
                            rna_dir  <- parts[2]

                            # For RNA intersection across all columns:
                            if (rna_dir == "up") {
                                rna_genes <- rownames(rna_mat_all)[rowSums(rna_mat_all > 0, na.rm = TRUE) == ncol(rna_mat_all)]
                            } else {
                                rna_genes <- rownames(rna_mat_all)[rowSums(rna_mat_all < 0, na.rm = TRUE) == ncol(rna_mat_all)]
                            }
                            # For genomic direction:
                            if (genomic_dir == "up") {
                                genomic_genes <- rownames(genomic_mat_all)[which(genomic_mat_all[, genomic_cmp] > 0)]
                            } else {
                                genomic_genes <- rownames(genomic_mat_all)[which(genomic_mat_all[, genomic_cmp] < 0)]
                            }

                            genes_combined <- intersect(genomic_genes, rna_genes)
                            if (length(genes_combined) == 0) {
                                cat(sprintf("Skipping combined %s (no genes).\n", dir_pair))
                                next
                            }
                            mat_all_rna  <- rna_mat_all[genes_combined, , drop = FALSE]
                            genomic_col_sub <- genomic_mat_all[genes_combined, genomic_cmp, drop = FALSE]
                            combined_mat <- cbind(mat_all_rna, genomic_col_sub)
                            draw_and_save(
                                mat = combined_mat,
                                title_text = paste0("Genomic", genomic_dir, " vs RNA_", rna_dir),
                                base_name  = paste0("Genomic", genomic_dir, "_RNA", rna_dir)
                            )
                        }
                    }

                    #--------------------------------------
                    # 11) CUSTOMIZED MODE: loop until user says no
                    #--------------------------------------
                    if (hm_type == 2) {
                        repeat {
                            cat("\n--- CUSTOMIZED MODE: Pick genomic direction(s) and RNA comparisons/directions ---\n")
                            cat("genomic directions available:\n")
                            cat("[1] genomic_up   (hyper-methylated genes)\n")
                            cat("[2] genomic_down (hypo-methylated genes)\n")
                            cat("[3] genomic (all significant genes)\n")
                            repeat {
                                genomic_input <- handle_quit(
                                readline(prompt = "Enter 1, 2, 3, or combinations (e.g. 1,2); or 'back' to return to menu:\n")
                                )
                                if (genomic_input == "back") {
                                    cat("Returning to main menu.\n")
                                    break
                                }
                                tokens <- strsplit(genomic_input, ",")[[1]] %>% str_trim()
                                valid <- all(tokens %in% c("1","2","3","1","2","3"))
                                if (!valid) {
                                    cat("Invalid entry. Choose 1, 2, 3, or combination.\n")
                                    next
                                }
                                break
                            }
                            if (!exists("genomic_input")) break
                            genomic_idxs <- unique(as.integer(tokens))
                            chosen_genomic_dirs <- c()
                            if (3 %in% genomic_idxs) {
                                chosen_genomic_dirs <- c("combined")
                            } else {
                                if (1 %in% genomic_idxs) chosen_genomic_dirs <- c(chosen_genomic_dirs, "up")
                                if (2 %in% genomic_idxs) chosen_genomic_dirs <- c(chosen_genomic_dirs, "down")
                            }

                            all_rna_cols <- colnames(rna_mat_all)
                            cat("\nRNA comparisons available (base index => comparison name):\n")
                            for (i in seq_along(all_rna_cols)) {
                                base_name <- all_rna_cols[i]
                                cat(sprintf("[%d]        %s\n", i, base_name))
                                cat(sprintf("[%d_up]   %s_up   [%d_down]   %s_down\n", i, base_name, i, base_name))
                            }
                            repeat {
                                rna_input <- handle_quit(
                                readline(prompt = "Enter indices (e.g. \"1_up,2_down,3\"), or 'back' to return to menu:\n")
                                )
                                if (rna_input == "back") {
                                    cat("Returning to main menu.\n")
                                    break
                                }
                                tokens2 <- strsplit(rna_input, ",")[[1]] %>% str_trim()
                                chosen_rna_list <- list()
                                valid2 <- TRUE
                                for (tok in tokens2) {
                                    if (grepl("^[0-9]+$", tok)) {
                                        idx <- as.integer(tok)
                                        chosen_rna_list <- c(chosen_rna_list, list(list(idx = idx, dir = "both")))
                                    } else if (grepl("^[0-9]+_up$", tok)) {
                                        idx <- as.integer(sub("_up$", "", tok))
                                        chosen_rna_list <- c(chosen_rna_list, list(list(idx = idx, dir = "up")))
                                    } else if (grepl("^[0-9]+_down$", tok)) {
                                        idx <- as.integer(sub("_down$", "", tok))
                                        chosen_rna_list <- c(chosen_rna_list, list(list(idx = idx, dir = "down")))
                                    } else {
                                        valid2 <- FALSE
                                        break
                                    }
                                }
                                if (!valid2) {
                                    cat("Invalid token detected. Use format: 1, 2_up, 3_down, etc.\n")
                                    next
                                }
                                ok_idx <- sapply(chosen_rna_list, function(x) x$idx >= 1 && x$idx <= length(all_rna_cols))
                                if (!all(ok_idx)) {
                                    cat("One or more indices out of range. Try again.\n")
                                    next
                                }
                                break
                            }
                            if (!exists("chosen_rna_list")) break

                            if ("combined" %in% chosen_genomic_dirs) {
                                genomic_dir_genes <- rownames(genomic_mat_all)
                            } else {
                                genomic_dir_genes <- c()
                                if ("up" %in% chosen_genomic_dirs) genomic_dir_genes <- union(genomic_dir_genes, filter_genes_by_direction(genomic_mat_all[, genomic_cmp], "up"))
                                if ("down" %in% chosen_genomic_dirs) genomic_dir_genes <- union(genomic_dir_genes, filter_genes_by_direction(genomic_mat_all[, genomic_cmp], "down"))
                            }

                            rna_selected_matrices <- list()
                            col_names <- c()
                            for (entry in chosen_rna_list) {
                                idx <- entry$idx
                                dir <- entry$dir
                                cmp <- all_rna_cols[idx]
                                vec <- rna_mat_all[, cmp]
                                if (dir == "both") {
                                    col_vec <- vec
                                    col_names <- c(col_names, cmp)
                                    rna_selected_matrices[[cmp]] <- col_vec
                                } else {
                                    genes_dir <- filter_genes_by_direction(vec, dir)
                                    col_vec <- rep(NA, length(vec))
                                    names(col_vec) <- names(vec)
                                    col_vec[genes_dir] <- vec[genes_dir]
                                    col_name <- paste0(cmp, "_", dir)
                                    col_names <- c(col_names, col_name)
                                    rna_selected_matrices[[col_name]] <- col_vec
                                }
                            }
                            if (length(rna_selected_matrices) > 0) {
                                mat_rna_custom <- do.call(cbind, rna_selected_matrices)
                                colnames(mat_rna_custom) <- col_names
                            } else {
                                mat_rna_custom <- matrix(nrow = 0, ncol = 0)
                            }

                            if (ncol(mat_rna_custom) == 0) {
                                genes_custom <- genomic_dir_genes
                            } else {
                                rna_nonNA <- rowSums(!is.na(mat_rna_custom)) >= 1
                                rna_genes_all <- rownames(mat_rna_custom)[rna_nonNA]
                                genes_custom <- intersect(genomic_dir_genes, rna_genes_all)
                            }

                            if (length(genes_custom) == 0) {
                                cat("No genes remain after filtering for this custom selection.\n")
                            } else {
                                final_rna_mat <- mat_rna_custom[genes_custom, , drop = FALSE]
                                genomic_col_sub  <- genomic_mat_all[genes_custom, genomic_cmp, drop = FALSE]
                                combined_custom <- cbind(final_rna_mat, genomic_col_sub)
                                colnames(combined_custom)[ncol(combined_custom)] <- genomic_cmp
                                safe_genomic_part <- paste0(
                                    ifelse("combined" %in% chosen_genomic_dirs, "genomic_combined",
                                        paste0(
                                        ifelse("up" %in% chosen_genomic_dirs, "genomic_up", ""),
                                        ifelse(all(c("up","down") %in% chosen_genomic_dirs), "_", ""),
                                        ifelse("down" %in% chosen_genomic_dirs, "genomic_down", "")
                                        )
                                    )
                                )
                                safe_rna_part <- paste0(sapply(col_names, function(x) gsub("[^A-Za-z0-9]", "_", x)), collapse = "_")
                                base_name_custom <- paste0("heatmap_custom_", safe_genomic_part, "_vs_RNA_", safe_rna_part, "_combined")
                                title_custom <- paste0("Custom: genomic (", paste(chosen_genomic_dirs, collapse = ", "), 
                                                    ") vs RNA(", paste(col_names, collapse = ", "), ")")
                                draw_and_save(
                                    mat = combined_custom,
                                    title_text = title_custom,
                                    base_name  = base_name_custom
                                )
                            }

                            repeat {
                                more_input <- handle_quit(
                                readline(prompt = "Generate another customized heatmap? (yes/no), or 'back' to exit:\n")
                                )
                                if (more_input == "back") {
                                    cat("Exiting customized heatmap mode.\n")
                                    break
                                }
                                mi <- tolower(more_input)
                                if (mi %in% c("yes", "y")) {
                                    break
                                }
                                if (mi %in% c("no", "n")) {
                                    cat("Exiting customized heatmap mode.\n")
                                    break
                                }
                                cat("Please answer 'yes' or 'no'.\n")
                            }
                            if (exists("more_input") && (tolower(more_input) %in% c("no", "n") || more_input == "back")) {
                                break
                            }
                        }
                    }
                } # Heatmap generation
                if (user_choice == 3) {
                    # **Option 3: Create Venn Diagrams**
                    # (A) Ask user about splitting into up/down
                    cat("Do you want to split results into upregulated and downregulated genes? (y/n):\n")
                    split_choice <- tolower(handle_quit(readline()))
                    if (! split_choice %in% c("y","n")) {
                        cat("Please enter 'y' or 'n'.\n")
                        next
                    }

                    # (B) Ask for padj threshold
                    cat("Enter the padj threshold for significant genes (or type 'back' to cancel):\n")
                    padj_threshold_raw <- handle_quit(readline())
                    if (padj_threshold_raw == "back") next
                    padj_threshold <- as.numeric(padj_threshold_raw)
                    if (is.na(padj_threshold) || padj_threshold <= 0 || padj_threshold > 1) {
                        cat("Invalid padj threshold. Please enter a value between 0 and 1.\n")
                        next
                    }

                    # (C) Filter DEGs
                    degs_data <- degs_data[!is.na(degs_data$padj) & degs_data$padj < padj_threshold, ]
                    if (nrow(degs_data) == 0) {
                        cat("No significant genes found for the given padj threshold. Returning to menu.\n")
                        next
                    }

                    # (D) Split into up/down if requested
                    if (split_choice == "y") {
                        degs_data_up   <- subset(degs_data, log2FoldChange > 0)
                        degs_data_down <- subset(degs_data, log2FoldChange < 0)
                    } else {
                        degs_data_up   <- NULL
                        degs_data_down <- NULL
                    }

                    # (E) Build indexed list of comparisons (+ “_up” and “_down” if split)
                    comparisons   <- unique(degs_data$comparison)
                    options_list  <- list()
                    for (i in seq_along(comparisons)) {
                        cmp <- comparisons[i]
                        options_list[[paste0(i)]] <- cmp
                        if (split_choice == "y") {
                            options_list[[paste0(i, "_up")]]   <- paste0(cmp, "_up")
                            options_list[[paste0(i, "_down")]] <- paste0(cmp, "_down")
                        }
                    }
                    # Ask the user if they want default or custom Venn diagrams
                    cat("Do you want to create default Venn diagrams or custom? ([1] default/[2] custom), or type 'back' to cancel:\n")
                    venn_choice_raw <- handle_quit(readline())
                    if (venn_choice_raw == "back") next
                    venn_choice <- as.numeric(venn_choice_raw)
                    if (venn_choice == 1) {
                        repeat {
                            selected_up   <- list()
                            selected_down <- list()
                            for (cmp in comparisons) {
                                up_vec   <- subset(degs_data,   comparison == cmp & log2FoldChange > 0)$Gene
                                down_vec <- subset(degs_data,   comparison == cmp & log2FoldChange < 0)$Gene
                                if (length(up_vec)   > 0) selected_up[[paste0(cmp, "_up")]]   <- unique(up_vec)
                                if (length(down_vec) > 0) selected_down[[paste0(cmp, "_down")]] <- unique(down_vec)
                            }
                            selected_up   <- selected_up  [sapply(selected_up,   length) > 0]
                            selected_down <- selected_down[sapply(selected_down, length) > 0]

                            if (split_choice == "y") {
                                # Launch the parallel gadget if at least one side has ≥2 sets
                                if (length(selected_up) >= 2 || length(selected_down) >= 2) {
                                    export_venn_upset_split(selected_up, selected_down, subdir_prefix = "genomic")
                                    cat("Run another interactive UpSet? (y/n):\n")
                                    again <- tolower(handle_quit(readline()))
                                    if (again != "y") break

                                } else {
                                    cat("Need at least two UP or two DOWN sets to launch. Returning to menu.\n")
                                }
                            } else {
                                # No splitting: just run one UpSet on the combined list
                                selected_all <- list()
                                for (cmp in comparisons) {
                                    all_vec <- subset(degs_data, comparison == cmp)$Gene
                                    if (length(all_vec) > 0) selected_all[[cmp]] <- unique(all_vec)
                                }
                                selected_all <- selected_all[sapply(selected_all, length) > 0]

                                if (length(selected_all) >= 2) {
                                    cat("Launching interactive UpSet (no split)…\n")
                                    export_venn_upset(selected_all, subdir = "genomic")
                                    # After user closes the Shiny gadget, ask if they want to run again:
                                    cat("Run another default interactive UpSet? (y/n):\n")
                                    again <- tolower(handle_quit(readline()))
                                    if (again != "y") break
                                } else {
                                    cat("Need at least two sets to launch interactive UpSet. Returning to menu.\n")
                                }
                            }
                        } # repeat for default Venn diagrams
                    } # Default Venn diagrams
                    if (venn_choice == 2) {
                        repeat{
                            cat("Available comparisons:\n")
                            for (opt in names(options_list)) {
                                if (split_choice == "n" && grepl("_up|_down", opt)) next
                                cat(sprintf("[%s] %s\n", opt, options_list[[opt]]))
                            }

                            cat("Select the indices (e.g., 1,2,3) for comparisons to include in the Venn diagram (or type 'back' to cancel):\n")
                            selection <- handle_quit(readline())
                            if (selection == "back") break

                            selected_indices <- unlist(strsplit(selection, ","))
                            selected_genes   <- list()

                            for (idx in selected_indices) {
                                comparison_name <- options_list[[idx]]
                                if (grepl("_up$", idx)) {
                                    comp <- gsub("_up$", "", comparison_name)
                                    selected_genes[[comparison_name]] <- unique(degs_data_up$Gene[degs_data_up$comparison == comp])
                                } else if (grepl("_down$", idx)) {
                                    comp <- gsub("_down$", "", comparison_name)
                                    selected_genes[[comparison_name]] <- unique(degs_data_down$Gene[degs_data_down$comparison == comp])
                                } else {
                                    selected_genes[[comparison_name]] <- unique(degs_data$Gene[degs_data$comparison == comparison_name])
                                }
                            }
                            # Drop any empty sets
                            selected_genes <- selected_genes[sapply(selected_genes, length) > 0]
                            if (length(selected_genes) < 2) {
                                cat("At least two non‐empty sets are required. Try again.\n")
                                next
                            }

                            # Instead of building a static Venn, launch the interactive UpSet+Venn gadget for these custom sets:
                            cat("Launching interactive UpSet+Venn for your custom selection…\n")
                            export_venn_upset(selected_genes, subdir = "genomic")

                            # After the user closes the Shiny gadget, ask if they want to generate another custom Venn
                            cat("Do you want to generate another custom Venn? (y/n):\n")
                            again <- tolower(handle_quit(readline()))
                            if (again != "y") break
                        } # repeat for custom Venn diagrams 
                    } # Custom Venn diagrams
                } # venn diagrams
            } # End of interactive analysis menu
        } # End of loading DEGs and genomic results 
    } # End of main menu loop
} # End of main function
