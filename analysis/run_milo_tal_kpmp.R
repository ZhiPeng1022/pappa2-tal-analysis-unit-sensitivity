suppressPackageStartupMessages({
    library(Matrix)
    library(SingleCellExperiment)
    library(miloR)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
    sub("^--file=", "", file_arg[[1]])
} else {
    "scripts/run_milo_tal_kpmp.R"
}
script_dir <- dirname(normalizePath(script_path))
project_root <- normalizePath(file.path(script_dir, "..", ".."))
data_root <- Sys.getenv(
    "PAPPA2_DATA_ROOT",
    unset = file.path(project_root, "data", "RAMP3")
)
metadata_path <- Sys.getenv(
    "MILO_METADATA_PATH",
    unset = file.path(
        data_root,
        "milo_full_tal_input",
        "KPMP_milo_full_tal_cells.csv"
    )
)
output_dir <- Sys.getenv(
    "MILO_OUTPUT_DIR",
    unset = file.path(project_root, "outputs")
)
output_prefix <- Sys.getenv(
    "MILO_OUTPUT_PREFIX",
    unset = "KPMP_milo_full_TAL_final"
)
seed_text <- Sys.getenv("MILO_SEED", unset = "")
export_membership <- identical(
    tolower(Sys.getenv("MILO_EXPORT_MEMBERSHIP", unset = "0")),
    "1"
)
norm_method <- Sys.getenv(
    "MILO_NORM_METHOD",
    unset = "logMS"
)

cell_meta <- read.csv(
    metadata_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM"
)
rownames(cell_meta) <- cell_meta$cell_id


build_milo <- function(meta) {
    counts <- Matrix(
        0,
        nrow = 1,
        ncol = nrow(meta),
        sparse = TRUE
    )
    rownames(counts) <- "dummy"
    colnames(counts) <- meta$cell_id
    reduced_dim <- as.matrix(meta[, c("UMAP1", "UMAP2")])
    rownames(reduced_dim) <- meta$cell_id

    sce <- SingleCellExperiment(
        assays = list(logcounts = counts),
        colData = DataFrame(meta),
        reducedDims = SimpleList(UMAP = reduced_dim)
    )
    milo <- Milo(sce)
    milo <- buildGraph(
        milo,
        k = 30,
        d = 2,
        reduced.dim = "UMAP"
    )
    milo <- makeNhoods(
        milo,
        prop = 0.1,
        k = 30,
        d = 2,
        refined = TRUE,
        reduced_dims = "UMAP"
    )
    milo <- calcNhoodDistance(
        milo,
        d = 2,
        reduced.dim = "UMAP"
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
    message("Running TAL Milo comparison: ", label)
    if (nzchar(seed_text)) {
        set.seed(as.integer(seed_text) + seed_offset)
    }
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

    milo <- build_milo(comparison)
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
        reduced.dim = "UMAP",
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
    output_path <- file.path(
        output_dir,
        paste0(
            output_prefix,
            "_",
            label,
            "_nhood_results.csv"
        )
    )
    write.csv(
        da_results,
        output_path,
        row.names = FALSE
    )
    if (export_membership) {
        membership_path <- file.path(
            output_dir,
            paste0(
                output_prefix,
                "_",
                label,
                "_nhood_membership.csv.gz"
            )
        )
        write_membership(milo, comparison, membership_path)
        message("Saved membership: ", membership_path)
    }
    message(
        "Saved ", output_path,
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
