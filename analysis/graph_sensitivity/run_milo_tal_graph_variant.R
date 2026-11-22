suppressPackageStartupMessages({
    library(Matrix)
    library(SingleCellExperiment)
    library(miloR)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(script_arg) > 0) {
    sub("^--file=", "", script_arg[[1]])
} else {
    "work/run_milo_tal_graph_variant.R"
}
script_dir <- dirname(normalizePath(script_path))

metadata_path <- Sys.getenv(
    "MILO_METADATA_PATH",
    unset = file.path(
        "D:/CodexData/RAMP3/milo_full_tal_input",
        "KPMP_milo_full_tal_cells.csv"
    )
)
output_dir <- Sys.getenv("MILO_OUTPUT_DIR")
if (!nzchar(output_dir)) {
    stop("MILO_OUTPUT_DIR must be set")
}
output_prefix <- Sys.getenv(
    "MILO_OUTPUT_PREFIX",
    unset = "KPMP_milo_graph"
)
k <- as.integer(Sys.getenv("MILO_K", unset = "30"))
prop <- as.numeric(Sys.getenv("MILO_PROP", unset = "0.1"))
seed_text <- Sys.getenv("MILO_SEED", unset = "20261122")
norm_method <- Sys.getenv("MILO_NORM_METHOD", unset = "logMS")
reduced_dim_name <- Sys.getenv(
    "MILO_REDUCED_DIM",
    unset = "UMAP"
)
dim_columns_text <- Sys.getenv(
    "MILO_DIM_COLUMNS",
    unset = ""
)
if (!nzchar(dim_columns_text)) {
    if (reduced_dim_name == "UMAP") {
        dim_columns_text <- "UMAP1,UMAP2"
    } else {
        stop("MILO_DIM_COLUMNS is required for a non-UMAP embedding")
    }
}
dim_columns <- trimws(
    strsplit(dim_columns_text, ",", fixed = TRUE)[[1]]
)
if (length(dim_columns) < 2) {
    stop("At least two embedding columns are required")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cell_meta <- read.csv(
    metadata_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM"
)
rownames(cell_meta) <- cell_meta$cell_id
missing_dim_columns <- setdiff(dim_columns, colnames(cell_meta))
if (length(missing_dim_columns) > 0) {
    stop(
        "Missing embedding columns: ",
        paste(missing_dim_columns, collapse = ", ")
    )
}

build_milo <- function(meta, k_value) {
    counts <- Matrix(
        0,
        nrow = 1,
        ncol = nrow(meta),
        sparse = TRUE
    )
    rownames(counts) <- "dummy"
    colnames(counts) <- meta$cell_id
    reduced_dim <- as.matrix(meta[, dim_columns, drop = FALSE])
    rownames(reduced_dim) <- meta$cell_id
    reduced_dims <- SimpleList(reduced_dim)
    names(reduced_dims) <- reduced_dim_name

    sce <- SingleCellExperiment(
        assays = list(logcounts = counts),
        colData = DataFrame(meta),
        reducedDims = reduced_dims
    )
    milo <- Milo(sce)
    milo <- buildGraph(
        milo,
        k = k_value,
        d = length(dim_columns),
        reduced.dim = reduced_dim_name
    )
    milo <- makeNhoods(
        milo,
        prop = prop,
        k = k_value,
        d = length(dim_columns),
        refined = TRUE,
        reduced_dims = reduced_dim_name
    )
    milo <- calcNhoodDistance(
        milo,
        d = length(dim_columns),
        reduced.dim = reduced_dim_name
    )
    milo <- countCells(
        milo,
        samples = "donor_id",
        meta.data = meta
    )
    milo
}

write_membership <- function(milo, meta, output_path) {
    membership <- summary(nhoods(milo))
    stopifnot(nrow(nhoods(milo)) == nrow(meta))
    cell_index <- membership$i
    membership_table <- data.frame(
        cell_id = meta$cell_id[cell_index],
        subclass = meta$subclass[cell_index],
        donor_id = meta$donor_id[cell_index],
        condition = meta$condition[cell_index],
        Nhood = membership$j - 1L,
        weight = membership$x,
        stringsAsFactors = FALSE
    )
    connection <- gzfile(output_path, open = "wt")
    on.exit(close(connection), add = TRUE)
    write.csv(
        membership_table,
        connection,
        row.names = FALSE
    )
}

run_comparison <- function(label, conditions, seed_offset) {
    message("Running Milo graph variant: ", label, " k=", k)
    set.seed(as.integer(seed_text) + seed_offset)
    comparison <- cell_meta[
        cell_meta$condition %in% conditions,
    ]
    comparison$group <- factor(
        ifelse(
            comparison$condition == "Reference",
            "Reference",
            "Disease"
        ),
        levels = c("Reference", "Disease")
    )

    milo <- build_milo(comparison, k)
    nhood_count_matrix <- nhoodCounts(milo)
    sample_totals <- Matrix::colSums(nhood_count_matrix)
    keep_samples <- names(sample_totals)[sample_totals > 0]
    dropped_samples <- setdiff(
        colnames(nhood_count_matrix),
        keep_samples
    )
    if (length(dropped_samples) > 0) {
        message(
            "Dropping ",
            length(dropped_samples),
            " zero-neighborhood donor(s): ",
            paste(dropped_samples, collapse = ", ")
        )
    }
    sample_info <- unique(
        comparison[, c("donor_id", "group")]
    )
    rownames(sample_info) <- sample_info$donor_id
    sample_info <- sample_info[
        keep_samples,
        ,
        drop = FALSE
    ]
    nhood_totals <- Matrix::rowSums(
        nhood_count_matrix[
            ,
            keep_samples,
            drop = FALSE
        ]
    )
    keep_nhoods <- which(nhood_totals > 0)
    if (length(keep_nhoods) < nrow(nhood_count_matrix)) {
        message(
            "Dropping ",
            nrow(nhood_count_matrix) - length(keep_nhoods),
            " zero-total neighborhood(s)"
        )
    }

    da_results <- testNhoods(
        milo,
        design = ~ group,
        design.df = sample_info,
        reduced.dim = reduced_dim_name,
        norm.method = norm_method,
        subset.nhoods = keep_nhoods
    )
    da_results <- annotateNhoods(
        milo,
        da_results,
        coldata_col = "subclass"
    )
    da_results$nhood_purity <- da_results$subclass_fraction
    da_results$comparison <- label
    da_results <- da_results[
        order(
            da_results$SpatialFDR,
            -abs(da_results$logFC)
        ),
    ]

    result_path <- file.path(
        output_dir,
        paste0(
            output_prefix,
            "_",
            label,
            "_nhood_results.csv"
        )
    )
    membership_path <- file.path(
        output_dir,
        paste0(
            output_prefix,
            "_",
            label,
            "_nhood_membership.csv.gz"
        )
    )
    write.csv(da_results, result_path, row.names = FALSE)
    write_membership(milo, comparison, membership_path)
    message(
        "Saved ",
        result_path,
        "; significant nhoods: ",
        sum(da_results$SpatialFDR < 0.05, na.rm = TRUE)
    )
}

run_comparison(
    "all_disease",
    c("Reference", "AKI", "CKD"),
    1L
)
run_comparison(
    "AKI_vs_Reference",
    c("Reference", "AKI"),
    2L
)
run_comparison(
    "CKD_vs_Reference",
    c("Reference", "CKD"),
    3L
)
