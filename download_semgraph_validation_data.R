args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript download_semgraph_validation_data.R --genes-file path/to/genes.csv [--column gene_label] [--out-dir data]",
      "",
      "Accepted gene files:",
      "  - .txt with one gene symbol per line",
      "  - .csv / .tsv containing a column like gene, gene_label, hgnc_symbol, or symbol",
      "",
      "Outputs:",
      "  - STRING_edges.tsv",
      "  - CellRegulonDB_query_genes.txt",
      "  - CellRegulonDB_download_instructions.txt",
      sep = "\n"
    )
  )
}

parse_args <- function(x) {
  out <- list(
    genes_file = NULL,
    column = NULL,
    out_dir = file.path(getwd(), "data")
  )
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (identical(key, "--genes-file") && i < length(x)) {
      out$genes_file <- x[[i + 1L]]
      i <- i + 2L
    } else if (identical(key, "--column") && i < length(x)) {
      out$column <- x[[i + 1L]]
      i <- i + 2L
    } else if (identical(key, "--out-dir") && i < length(x)) {
      out$out_dir <- x[[i + 1L]]
      i <- i + 2L
    } else if (key %in% c("-h", "--help")) {
      usage()
      quit(save = "no", status = 0)
    } else {
      stop(sprintf("Unknown or incomplete argument: %s", key), call. = FALSE)
    }
  }
  if (is.null(out$genes_file) || !nzchar(out$genes_file)) {
    usage()
    stop("Missing required argument --genes-file", call. = FALSE)
  }
  out
}

read_gene_file <- function(path, column = NULL) {
  if (!file.exists(path)) stop(sprintf("Gene file not found: %s", path), call. = FALSE)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("txt", "list")) {
    genes <- readLines(path, warn = FALSE)
    return(unique(trimws(genes[nzchar(trimws(genes))])))
  }
  sep <- if (ext == "tsv") "\t" else ","
  tbl <- utils::read.table(path, header = TRUE, sep = sep, quote = "\"", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(tbl)) stop("Gene file is empty.", call. = FALSE)
  if (!is.null(column) && nzchar(column) && column %in% names(tbl)) {
    genes <- tbl[[column]]
  } else {
    lower_names <- tolower(names(tbl))
    idx <- match(TRUE, lower_names %in% c("gene", "gene_label", "hgnc_symbol", "symbol", "gene_symbol"), nomatch = 0L)
    if (idx == 0L) {
      stop("Could not detect a gene-symbol column. Use --column with the correct column name.", call. = FALSE)
    }
    genes <- tbl[[idx]]
  }
  genes <- unique(trimws(as.character(genes)))
  genes[nzchar(genes)]
}

download_string_edges <- function(genes) {
  if (length(genes) < 2L) {
    stop("Need at least two genes to query STRING.", call. = FALSE)
  }
  endpoint <- paste0(
    "https://string-db.org/api/tsv/network?identifiers=",
    utils::URLencode(paste(genes, collapse = "\r"), reserved = TRUE),
    "&species=9606&required_score=0&network_type=functional"
  )
  tbl <- utils::read.delim(endpoint, stringsAsFactors = FALSE, check.names = FALSE)
  if (is.null(tbl) || !nrow(tbl)) {
    warning("STRING returned no interactions for the supplied genes.")
    return(data.frame())
  }
  nms <- tolower(names(tbl))
  a_col <- names(tbl)[match(TRUE, nms %in% c("preferredname_a", "preferrednamea"), nomatch = 0L)]
  b_col <- names(tbl)[match(TRUE, nms %in% c("preferredname_b", "preferrednameb"), nomatch = 0L)]
  score_col <- names(tbl)[match(TRUE, nms %in% c("score", "combined_score"), nomatch = 0L)]
  if (!length(a_col) || !length(b_col) || !length(score_col)) {
    stop("STRING response did not contain expected columns.", call. = FALSE)
  }
  out <- data.frame(
    preferredName_A = as.character(tbl[[a_col]]),
    preferredName_B = as.character(tbl[[b_col]]),
    score = suppressWarnings(as.numeric(tbl[[score_col]])),
    stringsAsFactors = FALSE
  )
  out <- out[out$preferredName_A %in% genes & out$preferredName_B %in% genes, , drop = FALSE]
  out$score[is.na(out$score)] <- 0
  out
}

write_cellregulondb_instructions <- function(genes, out_dir) {
  gene_file <- file.path(out_dir, "CellRegulonDB_query_genes.txt")
  writeLines(genes, gene_file, useBytes = TRUE)
  instr_file <- file.path(out_dir, "CellRegulonDB_download_instructions.txt")
  instructions <- c(
    "CellRegulonDB export guide",
    "",
    "1. Open: https://www.cellregulondb.org/network.html",
    "2. Paste or enter the genes from CellRegulonDB_query_genes.txt into the Genes filter.",
    "3. If you know the relevant tissue, lineage, or cell type, add those filters to make the export more specific.",
    "4. Run the search and download the resulting interaction table as CSV.",
    "5. Save the downloaded file into this data directory as one of:",
    "   - CellRegulonDB_edges.csv",
    "   - CellRegulonDB_regulons.csv",
    "",
    "Expected columns for the app:",
    "   - a source TF/regulator column such as tf, regulator, source, or tf_symbol",
    "   - a target column such as target, target_gene, gene, or target_symbol",
    "Optional columns that improve the annotations:",
    "   - tissue / organ",
    "   - cell type",
    "   - regulon / network / module",
    "",
    "The app will then use that file for directed TF-target validation in the SEMgraph section."
  )
  writeLines(instructions, instr_file, useBytes = TRUE)
  list(gene_file = gene_file, instr_file = instr_file)
}

opt <- parse_args(args)
if (!dir.exists(opt$out_dir)) dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

genes <- read_gene_file(opt$genes_file, opt$column)
cat(sprintf("Loaded %d unique genes.\n", length(genes)))

string_out <- file.path(opt$out_dir, "STRING_edges.tsv")
string_tbl <- download_string_edges(genes)
utils::write.table(string_tbl, file = string_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Saved STRING interactions to %s\n", string_out))

cell_files <- write_cellregulondb_instructions(genes, opt$out_dir)
cat(sprintf("Saved CellRegulonDB gene list to %s\n", cell_files$gene_file))
cat(sprintf("Saved CellRegulonDB instructions to %s\n", cell_files$instr_file))
