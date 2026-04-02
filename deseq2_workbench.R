#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x

detect_worker_threads <- function() {
  detected <- NA_integer_
  detected <- tryCatch(parallel::detectCores(logical = TRUE), error = function(...) NA_integer_)
  if (!is.finite(detected) || is.na(detected) || detected < 1L) {
    env_threads <- suppressWarnings(as.integer(Sys.getenv("NUMBER_OF_PROCESSORS", unset = "1")))
    if (is.finite(env_threads) && !is.na(env_threads) && env_threads >= 1L) {
      detected <- env_threads
    } else {
      detected <- 1L
    }
  }
  max(1L, floor(detected * 0.8))
}

WORKER_THREADS <- detect_worker_threads()

configure_parallelism <- function(threads = WORKER_THREADS) {
  options(mc.cores = max(1L, threads))
  if (threads >= 2L && requireNamespace("WGCNA", quietly = TRUE)) {
    tryCatch(
      WGCNA::allowWGCNAThreads(nThreads = threads),
      error = function(...) NULL
    )
  }
  invisible(threads)
}

pkgs <- c("shiny", "plotly", "DT")
all_required_packages <- c(
  "shiny", "plotly", "DT", "DESeq2", "clusterProfiler", "WGCNA", "SEMgraph",
  "AnnotationDbi", "org.Hs.eg.db", "GO.db", "jsonlite", "igraph"
)
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  message("Missing CRAN packages: ", paste(missing, collapse = ", "))
  message("Install with: R -e \"install.packages(c(", paste(sprintf("'%s'", missing), collapse = ", "), "))\"")
  quit(status = 1L)
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  message("Missing Bioconductor package: DESeq2")
  message("Install with: R -e \"if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('DESeq2')\"")
  quit(status = 1L)
}
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  message("Missing Bioconductor package: clusterProfiler")
  message("Install with: R -e \"if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('clusterProfiler')\"")
  quit(status = 1L)
}
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  message("Missing CRAN package: WGCNA")
  message("Install with: R -e \"install.packages('WGCNA')\"")
  quit(status = 1L)
}
configure_parallelism()
if (!requireNamespace("SEMgraph", quietly = TRUE)) {
  message("Missing package: SEMgraph")
  message("Install with: R -e \"install.packages('SEMgraph')\"")
  quit(status = 1L)
}
if (!requireNamespace("igraph", quietly = TRUE)) {
  message("Missing package: igraph")
  message("Install with: R -e \"install.packages('igraph')\"")
  quit(status = 1L)
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  message("Missing CRAN package: jsonlite")
  message("Install with: R -e \"install.packages('jsonlite')\"")
  quit(status = 1L)
}
if (!requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  message("Missing Bioconductor annotation packages: AnnotationDbi and/or org.Hs.eg.db")
  message("Install with: R -e \"if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c('AnnotationDbi','org.Hs.eg.db'))\"")
  quit(status = 1L)
}
if (!requireNamespace("GO.db", quietly = TRUE)) {
  message("Missing Bioconductor annotation package: GO.db")
  message("Install with: R -e \"if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('GO.db')\"")
  quit(status = 1L)
}

data_dir <- file.path(getwd(), "data")
cache_dir <- file.path(getwd(), "cache")
if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
reactome_cache_dir <- file.path(cache_dir, "reactome")
if (!dir.exists(reactome_cache_dir)) dir.create(reactome_cache_dir, recursive = TRUE, showWarnings = FALSE)
default_counts <- file.path(data_dir, "GSE267213_counts.csv")
default_meta <- file.path(data_dir, "SraRunTable.csv")

read_inputs <- function(counts_path, meta_path, min_count, min_samples) {
  counts_df <- read.csv(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
  meta <- read.csv(meta_path, check.names = FALSE, stringsAsFactors = FALSE)
  sample_key <- if ("run_id" %in% names(meta)) {
    "run_id"
  } else if ("sample_id" %in% names(meta)) {
    "sample_id"
  } else {
    stop("Metadata file must contain either a 'run_id' column or a 'sample_id' column.")
  }
  if (sample_key != "run_id") {
    names(meta)[names(meta) == sample_key] <- "run_id"
  }
  stopifnot(ncol(counts_df) >= 2L)
  counts <- as.matrix(counts_df[, -1L, drop = FALSE])
  storage.mode(counts) <- "numeric"
  rownames(counts) <- counts_df[[1L]]
  samples <- colnames(counts)
  meta <- merge(data.frame(run_id = samples, ord = seq_along(samples), stringsAsFactors = FALSE), meta, by = "run_id", all.x = TRUE, sort = FALSE)
  meta <- meta[order(meta$ord), , drop = FALSE]
  meta$ord <- NULL
  for (nm in setdiff(names(meta), "run_id")) {
    meta[[nm]][is.na(meta[[nm]]) | meta[[nm]] == ""] <- "Unannotated"
    if (is.character(meta[[nm]])) meta[[nm]] <- factor(meta[[nm]])
  }
  keep <- rowSums(counts >= min_count) >= min_samples
  filtered <- counts[keep, , drop = FALSE]
  if (nrow(filtered) == 0L) stop("Filtering removed all genes.")
  lib <- colSums(filtered)
  log_cpm <- log2(sweep(filtered, 2L, lib / 1e6, "/") + 1)
  pca <- prcomp(t(log_cpm), center = TRUE, scale. = TRUE)
  gene_map <- map_gene_symbols(rownames(counts))
  list(counts = counts, filtered = filtered, meta = meta, log_cpm = log_cpm, pca = pca, dist = as.matrix(dist(t(log_cpm))), gene_map = gene_map)
}

resolve_input_path <- function(upload_info, typed_path) {
  if (!is.null(upload_info) && nrow(upload_info) > 0L && "datapath" %in% names(upload_info)) {
    return(upload_info$datapath[[1L]])
  }
  typed_path
}

map_gene_symbols <- function(gene_ids) {
  clean_ids <- sub("\\..*$", "", gene_ids)
  symbol_map <- suppressMessages(
    AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = unique(clean_ids),
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  )
  entrez_map <- suppressMessages(
    AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = unique(clean_ids),
      column = "ENTREZID",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  )
  fullname_map <- suppressMessages(
    AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = unique(clean_ids),
      column = "GENENAME",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  )
  symbols <- unname(symbol_map[clean_ids])
  entrez <- unname(entrez_map[clean_ids])
  full_names <- unname(fullname_map[clean_ids])
  symbols[is.na(symbols) | !nzchar(symbols)] <- gene_ids[is.na(symbols) | !nzchar(symbols)]
  full_names[is.na(full_names) | !nzchar(full_names)] <- "Full gene name unavailable"
  display <- ifelse(duplicated(symbols) | duplicated(symbols, fromLast = TRUE), paste0(symbols, " (", gene_ids, ")"), symbols)
  data.frame(
    gene_id = gene_ids,
    ensembl_id = clean_ids,
    entrez_id = entrez,
    hgnc_symbol = symbols,
    full_gene_name = full_names,
    gene_label = display,
    stringsAsFactors = FALSE
  )
}

attach_gene_labels <- function(df, gene_map) {
  merge(df, gene_map, by = "gene_id", all.x = TRUE, sort = FALSE)
}

annotate_tf_status <- function(gene_ids, gene_map) {
  tf_root_terms <- c("GO:0003700", "GO:0140110")
  available_bp <- intersect(tf_root_terms, names(GO.db::GOBPOFFSPRING))
  available_mf <- intersect(tf_root_terms, names(GO.db::GOMFOFFSPRING))
  tf_terms <- unique(unlist(c(
    tf_root_terms,
    if (length(available_bp) > 0L) unname(GO.db::GOBPOFFSPRING[available_bp]) else NULL,
    if (length(available_mf) > 0L) unname(GO.db::GOMFOFFSPRING[available_mf]) else NULL
  )))
  tf_terms <- tf_terms[!is.na(tf_terms)]

  ann <- gene_map[match(gene_ids, gene_map$gene_id), c("gene_id", "ensembl_id", "hgnc_symbol", "full_gene_name", "gene_label"), drop = FALSE]
  valid_keys <- unique(ann$ensembl_id[!is.na(ann$ensembl_id) & nzchar(ann$ensembl_id)])
  valid_keys <- valid_keys[grepl("^ENSG", valid_keys)]
  if (length(valid_keys) == 0L) {
    ann$is_tf <- FALSE
    ann$tf_source <- "No valid Ensembl IDs available for TF annotation"
    return(ann)
  }
  go_tbl <- tryCatch(
    suppressMessages(
      AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = valid_keys,
        keytype = "ENSEMBL",
        columns = c("GO", "ONTOLOGY")
      )
    ),
    error = function(...) data.frame()
  )
  if (is.null(go_tbl) || nrow(go_tbl) == 0L) {
    ann$is_tf <- FALSE
    ann$tf_source <- "No TF annotation found in GO.db/org.Hs.eg.db for the available Ensembl IDs"
    return(ann)
  }
  go_tbl <- go_tbl[go_tbl$GO %in% tf_terms, , drop = FALSE]
  tf_ens <- unique(go_tbl$ENSEMBL)
  ann$is_tf <- ann$ensembl_id %in% tf_ens
  ann$tf_source <- ifelse(
    ann$is_tf,
    "GO.db + org.Hs.eg.db transcription-factor GO terms",
    "No TF annotation found in GO.db/org.Hs.eg.db"
  )
  ann
}

extract_gene_ratio <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) == 2L) as.numeric(p[1L]) / as.numeric(p[2L]) else NA_real_
  }, numeric(1))
}

format_enrichment_df <- function(result_df, method_label) {
  if (is.null(result_df) || nrow(result_df) == 0L) {
    return(data.frame())
  }
  out <- as.data.frame(result_df, stringsAsFactors = FALSE)
  out$method <- method_label
  if ("Description" %in% names(out)) {
    out$pathway <- out$Description
  } else if ("ID" %in% names(out)) {
    out$pathway <- out$ID
  } else {
    out$pathway <- paste(method_label, seq_len(nrow(out)))
  }
  if ("GeneRatio" %in% names(out)) {
    out$gene_ratio_value <- extract_gene_ratio(out$GeneRatio)
  } else if ("setSize" %in% names(out) && "Count" %in% names(out)) {
    out$gene_ratio_value <- suppressWarnings(as.numeric(out$Count) / as.numeric(out$setSize))
  } else if ("NES" %in% names(out)) {
    out$gene_ratio_value <- abs(out$NES)
  } else {
    out$gene_ratio_value <- 1
  }
  out
}

bind_rows_fill <- function(dfs) {
  if (!length(dfs)) {
    return(data.frame())
  }
  all_cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  aligned <- lapply(dfs, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    if (length(missing_cols) > 0L) {
      for (col in missing_cols) {
        df[[col]] <- NA
      }
    }
    df[, all_cols, drop = FALSE]
  })
  do.call(rbind, aligned)
}

color_text_input <- function(id, label, value) {
  shiny::tags$div(
    style = "display:flex; align-items:flex-end; gap:8px;",
    shiny::tags$div(
      style = "flex:1 1 auto;",
      shiny::textInput(id, label, value = value)
    ),
    shiny::tags$div(
      id = paste0(id, "_preview"),
      class = "live-color-preview",
      `data-color-input` = id,
      title = "Current selected color"
    )
  )
}

.drug_cache <- new.env(parent = emptyenv())

make_cache_key <- function(symbols) {
  txt <- paste(sort(unique(symbols[!is.na(symbols) & nzchar(symbols)])), collapse = "|")
  if (!nzchar(txt)) return("empty")
  ints <- utf8ToInt(txt)
  hash <- 0
  mod <- 2147483647
  for (val in ints) {
    hash <- (hash * 131 + val) %% mod
  }
  paste0("k_", length(ints), "_", hash)
}

cache_file_path <- function(prefix, key_parts) {
  key <- make_cache_key(as.character(key_parts))
  file.path(cache_dir, paste0(prefix, "_", key, ".rds"))
}

read_cached_rds <- function(prefix, key_parts, refresh = FALSE) {
  path <- cache_file_path(prefix, key_parts)
  if (!refresh && file.exists(path)) {
    return(tryCatch(readRDS(path), error = function(...) NULL))
  }
  NULL
}

write_cached_rds <- function(prefix, key_parts, value) {
  path <- cache_file_path(prefix, key_parts)
  tryCatch(saveRDS(value, path), error = function(...) NULL)
  invisible(value)
}

reactome_mapping_file <- file.path(reactome_cache_dir, "NCBI2Reactome.txt")
reactome_mapping_urls <- c(
  "https://download.reactome.org/current/NCBI2Reactome.txt",
  "https://reactome.org/download/current/NCBI2Reactome.txt"
)
reactome_url_stats_file <- file.path(reactome_cache_dir, "reactome_url_stats.rds")

load_reactome_url_stats <- function() {
  if (!file.exists(reactome_url_stats_file)) {
    return(data.frame(url = reactome_mapping_urls, last_elapsed = Inf, success = FALSE, stringsAsFactors = FALSE))
  }
  stats <- tryCatch(readRDS(reactome_url_stats_file), error = function(...) NULL)
  if (is.null(stats) || !all(c("url", "last_elapsed", "success") %in% names(stats))) {
    return(data.frame(url = reactome_mapping_urls, last_elapsed = Inf, success = FALSE, stringsAsFactors = FALSE))
  }
  missing <- setdiff(reactome_mapping_urls, stats$url)
  if (length(missing)) {
    stats <- rbind(stats, data.frame(url = missing, last_elapsed = Inf, success = FALSE, stringsAsFactors = FALSE))
  }
  stats
}

save_reactome_url_stats <- function(stats) {
  tryCatch(saveRDS(stats, reactome_url_stats_file), error = function(...) NULL)
  invisible(stats)
}

rank_reactome_urls <- function() {
  stats <- load_reactome_url_stats()
  stats$pref_success <- ifelse(isTRUE(stats$success), 0, 1)
  stats$pref_elapsed <- ifelse(is.finite(stats$last_elapsed), stats$last_elapsed, Inf)
  stats <- stats[order(stats$pref_success, stats$pref_elapsed, match(stats$url, reactome_mapping_urls)), , drop = FALSE]
  stats$url
}

download_reactome_mapping <- function(force = FALSE, progress = NULL) {
  if (!force && file.exists(reactome_mapping_file)) {
    return(reactome_mapping_file)
  }
  old_timeout <- getOption("timeout")
  options(timeout = max(300, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  urls <- rank_reactome_urls()
  stats <- load_reactome_url_stats()
  for (i in seq_along(urls)) {
    url <- urls[i]
    if (!is.null(progress)) progress(0.02 + 0.01 * i, sprintf("Downloading Reactome mapping file (source %s/%s)", i, length(reactome_mapping_urls)))
    started <- proc.time()[["elapsed"]]
    ok <- tryCatch({
      utils::download.file(url, destfile = reactome_mapping_file, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(...) FALSE, warning = function(w) FALSE)
    elapsed <- proc.time()[["elapsed"]] - started
    stats$last_elapsed[stats$url == url] <- elapsed
    stats$success[stats$url == url] <- isTRUE(ok) && file.exists(reactome_mapping_file) && file.info(reactome_mapping_file)$size > 0
    save_reactome_url_stats(stats)
    if (isTRUE(ok) && file.exists(reactome_mapping_file) && file.info(reactome_mapping_file)$size > 0) {
      return(reactome_mapping_file)
    }
  }
  stop("Could not download Reactome mapping file from the available official sources within the timeout window.")
}

load_reactome_mapping <- function(force_refresh = FALSE, progress = NULL) {
  cached <- read_cached_rds("reactome_mapping_parsed", "ncbi2reactome", refresh = force_refresh)
  if (!is.null(cached)) return(cached)
  path <- download_reactome_mapping(force = force_refresh, progress = progress)
  if (!file.exists(path)) return(NULL)
  if (!is.null(progress)) progress(0.04, "Parsing local Reactome mapping cache")
  tbl <- tryCatch(
    utils::read.delim(
      path,
      header = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      fill = TRUE
    ),
    error = function(...) NULL
  )
  if (is.null(tbl) || ncol(tbl) < 6L) return(NULL)
  names(tbl)[1:6] <- c("entrez_id", "reactome_id", "url", "event_name", "evidence_code", "species")
  tbl <- tbl[tbl$species == "Homo sapiens", c("entrez_id", "reactome_id", "event_name", "species"), drop = FALSE]
  tbl$entrez_id <- as.character(tbl$entrez_id)
  write_cached_rds("reactome_mapping_parsed", "ncbi2reactome", tbl)
  tbl
}

classify_drug_action <- function(interaction_types) {
  txt <- paste(tolower(interaction_types %||% character()), collapse = " | ")
  if (!nzchar(txt)) return("unspecified")
  inhib <- grepl("inhib|antagon|block|suppress|decreas|negative|inverse agon|downreg", txt)
  activ <- grepl("activ|agon|stimulat|induc|potentiat|positive|upreg", txt)
  if (inhib && activ) return("mixed")
  if (inhib) return("inhibiting")
  if (activ) return("activating")
  "other or unclear"
}

empty_drug_annotation <- function(symbols) {
  data.frame(
    hgnc_symbol = symbols,
    is_druggable = FALSE,
    drug_mode = "no interaction found",
    drug_summary = "No external drug-target interaction found",
    drug_count = 0L,
    drug_source = "None",
    drug_source_details = "No DGIdb, HCDT, or Reactome drug context found",
    reactome_context = "No Reactome drug-mediated pathway context found",
    stringsAsFactors = FALSE
  )
}

empty_edge_validation <- function(edges = NULL) {
  if (is.null(edges) || !nrow(edges)) {
    return(data.frame(
      source = character(),
      target = character(),
      string_supported = logical(),
      string_score = numeric(),
      string_source = character(),
      string_action_supported = logical(),
      string_action_mode = character(),
      string_action_effect = character(),
      string_action_source = character(),
      string_direction_available = logical(),
      string_direction_matches = logical(),
      validation_tier = character(),
      validation_score = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    source = edges$source,
    target = edges$target,
    string_supported = FALSE,
    string_score = 0,
    string_source = "STRING lookup unavailable",
    string_action_supported = FALSE,
    string_action_mode = "No STRING action annotation",
    string_action_effect = "No STRING effect annotation",
    string_action_source = "STRING action lookup unavailable",
    string_direction_available = FALSE,
    string_direction_matches = FALSE,
    validation_tier = "No external support",
    validation_score = 0,
    stringsAsFactors = FALSE
  )
}

normalize_edge_key <- function(a, b) {
  ifelse(a <= b, paste(a, b, sep = "||"), paste(b, a, sep = "||"))
}

query_dgidb_targets <- function(symbols, refresh = FALSE) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (!length(symbols)) {
    return(empty_drug_annotation(character()))
  }

  cache_key <- make_cache_key(symbols)
  if (!refresh && exists(cache_key, envir = .drug_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .drug_cache, inherits = FALSE))
  }
  disk_cached <- read_cached_rds("dgidb", symbols, refresh = refresh)
  if (!is.null(disk_cached)) {
    assign(cache_key, disk_cached, envir = .drug_cache)
    return(disk_cached)
  }

  endpoint <- paste0(
    "https://dgidb.org/api/v2/interactions.json?genes=",
    utils::URLencode(paste(symbols, collapse = ","))
  )

  out <- tryCatch({
    raw <- jsonlite::fromJSON(endpoint, simplifyVector = FALSE)
    matched <- raw$matchedTerms %||% list()
    rows <- lapply(matched, function(term) {
      gene <- term$searchTerm %||% term$geneName %||% NA_character_
      interactions <- term$interactions %||% list()
      if (!length(interactions)) {
        return(data.frame(
          hgnc_symbol = gene,
          is_druggable = FALSE,
          drug_mode = "no interaction found",
          drug_summary = "No DGIdb interaction found",
          drug_count = 0L,
          drug_source = "DGIdb",
          drug_source_details = "DGIdb: no interaction found",
          reactome_context = NA_character_,
          stringsAsFactors = FALSE
        ))
      }
      drugs <- vapply(interactions, function(x) {
        x$drugName %||% x$drug_claim_name %||% x$drugConceptId %||% "Unnamed drug"
      }, character(1))
      interaction_types <- lapply(interactions, function(x) x$interactionTypes %||% character())
      modes <- vapply(interaction_types, classify_drug_action, character(1))
      unique_drugs <- unique(drugs[!is.na(drugs) & nzchar(drugs)])
      mode_counts <- sort(table(modes), decreasing = TRUE)
      top_mode <- if (length(mode_counts)) names(mode_counts)[1L] else "unspecified"
      top_pairs <- unique(paste0(drugs, " (", modes, ")"))
      top_pairs <- top_pairs[!is.na(top_pairs) & nzchar(top_pairs)]
      summary_txt <- if (length(top_pairs)) {
        paste(head(top_pairs, 6L), collapse = "; ")
      } else {
        "DGIdb interaction present, but no readable drug names were returned"
      }
      data.frame(
        hgnc_symbol = gene,
        is_druggable = TRUE,
        drug_mode = top_mode,
        drug_summary = summary_txt,
        drug_count = length(unique_drugs),
        drug_source = "DGIdb drug-gene interaction knowledgebase",
        drug_source_details = paste0("DGIdb: ", summary_txt),
        reactome_context = NA_character_,
        stringsAsFactors = FALSE
      )
    })
    if (!length(rows)) {
      data.frame(
        hgnc_symbol = symbols,
        is_druggable = FALSE,
        drug_mode = "no interaction found",
        drug_summary = "No DGIdb interaction found",
        drug_count = 0L,
        drug_source = "DGIdb",
        drug_source_details = "DGIdb: no interaction found",
        reactome_context = NA_character_,
        stringsAsFactors = FALSE
      )
    } else {
      res <- do.call(rbind, rows)
      missing_syms <- setdiff(symbols, unique(res$hgnc_symbol))
      if (length(missing_syms)) {
        res <- rbind(
          res,
          data.frame(
            hgnc_symbol = missing_syms,
            is_druggable = FALSE,
            drug_mode = "no interaction found",
            drug_summary = "No DGIdb interaction found",
            drug_count = 0L,
            drug_source = "DGIdb",
            drug_source_details = "DGIdb: no interaction found",
            reactome_context = NA_character_,
            stringsAsFactors = FALSE
          )
        )
      }
      res
    }
  }, error = function(...) {
    data.frame(
      hgnc_symbol = symbols,
      is_druggable = FALSE,
      drug_mode = "lookup unavailable",
      drug_summary = "DGIdb lookup could not be completed on this run",
      drug_count = 0L,
      drug_source = "DGIdb",
      drug_source_details = "DGIdb lookup unavailable on this run",
      reactome_context = NA_character_,
      stringsAsFactors = FALSE
    )
  })

  out <- out[!duplicated(out$hgnc_symbol), , drop = FALSE]
  assign(cache_key, out, envir = .drug_cache)
  write_cached_rds("dgidb", symbols, out)
  out
}

locate_hcdt_file <- function() {
  explicit <- c(
    file.path(data_dir, "HCDT_drug_gene.tsv"),
    file.path(data_dir, "HCDT_drug_gene.csv"),
    file.path(data_dir, "HCDT.tsv"),
    file.path(data_dir, "HCDT.csv")
  )
  hits <- explicit[file.exists(explicit)]
  if (length(hits)) return(hits[1L])
  fuzzy <- list.files(data_dir, pattern = "HCDT.*\\.(tsv|csv|txt)$", full.names = TRUE, ignore.case = TRUE)
  if (length(fuzzy)) fuzzy[1L] else ""
}

read_hcdt_table <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  reader <- if (ext %in% c("tsv", "txt")) utils::read.delim else utils::read.csv
  tryCatch(reader(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(...) NULL)
}

query_hcdt_targets <- function(symbols, refresh = FALSE) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (!length(symbols)) return(empty_drug_annotation(character()))
  disk_cached <- read_cached_rds("hcdt", symbols, refresh = refresh)
  if (!is.null(disk_cached)) return(disk_cached)
  hcdt_path <- locate_hcdt_file()
  if (!nzchar(hcdt_path)) {
    out <- empty_drug_annotation(symbols)
    out$drug_source <- "HCDT"
    out$drug_source_details <- "HCDT local export not found in data/. Add an HCDT drug-gene export to enable this source."
    return(out)
  }
  tbl <- read_hcdt_table(hcdt_path)
  if (is.null(tbl) || nrow(tbl) == 0L) {
    out <- empty_drug_annotation(symbols)
    out$drug_source <- "HCDT"
    out$drug_source_details <- "HCDT file was found but could not be parsed"
    return(out)
  }
  nms <- tolower(names(tbl))
  gene_col <- names(tbl)[match(TRUE, nms %in% c("gene_symbol", "hgnc_symbol", "symbol", "gene"), nomatch = 0L)]
  drug_col <- names(tbl)[match(TRUE, nms %in% c("drug_name", "drug", "drugname", "name"), nomatch = 0L)]
  action_col <- names(tbl)[match(TRUE, nms %in% c("action", "interaction_type", "interaction", "mechanism", "effect"), nomatch = 0L)]
  if (length(gene_col) == 0L || length(drug_col) == 0L) {
    out <- empty_drug_annotation(symbols)
    out$drug_source <- "HCDT"
    out$drug_source_details <- "HCDT file is missing recognizable gene/drug columns"
    return(out)
  }
  keep <- tbl[[gene_col]] %in% symbols
  sub <- tbl[keep, , drop = FALSE]
  if (!nrow(sub)) {
    out <- empty_drug_annotation(symbols)
    out$drug_source <- "HCDT"
    out$drug_source_details <- "HCDT: no interaction found"
    return(out)
  }
  rows <- lapply(symbols, function(sym) {
    cur <- sub[sub[[gene_col]] == sym, , drop = FALSE]
    if (!nrow(cur)) {
      return(data.frame(
        hgnc_symbol = sym,
        is_druggable = FALSE,
        drug_mode = "no interaction found",
        drug_summary = "No HCDT interaction found",
        drug_count = 0L,
        drug_source = "HCDT",
        drug_source_details = "HCDT: no interaction found",
        reactome_context = NA_character_,
        stringsAsFactors = FALSE
      ))
    }
    drugs <- unique(as.character(cur[[drug_col]]))
    drugs <- drugs[!is.na(drugs) & nzchar(drugs)]
    action_vals <- if (length(action_col)) as.character(cur[[action_col]]) else character()
    mode <- classify_drug_action(action_vals)
    summary_txt <- if (length(drugs)) {
      paste(head(unique(paste0(drugs, if (length(action_vals) && any(nzchar(action_vals))) paste0(" (", head(action_vals, length(drugs)), ")") else "")), 6L), collapse = "; ")
    } else {
      "HCDT interaction found, but no readable drug name column value was available"
    }
    data.frame(
      hgnc_symbol = sym,
      is_druggable = TRUE,
      drug_mode = mode,
      drug_summary = summary_txt,
      drug_count = length(unique(drugs)),
      drug_source = "HCDT",
      drug_source_details = paste0("HCDT: ", summary_txt),
      reactome_context = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  write_cached_rds("hcdt", symbols, out)
  out
}

locate_cellregulondb_file <- function() {
  explicit <- c(
    file.path(data_dir, "CellRegulonDB_edges.tsv"),
    file.path(data_dir, "CellRegulonDB_edges.csv"),
    file.path(data_dir, "CellRegulonDB_regulons.tsv"),
    file.path(data_dir, "CellRegulonDB_regulons.csv"),
    file.path(data_dir, "cellregulondb.tsv"),
    file.path(data_dir, "cellregulondb.csv")
  )
  hits <- explicit[file.exists(explicit)]
  if (length(hits)) return(hits[1L])
  fuzzy <- list.files(data_dir, pattern = "CellRegulonDB.*\\.(tsv|csv|txt)$", full.names = TRUE, ignore.case = TRUE)
  if (length(fuzzy)) fuzzy[1L] else ""
}

locate_string_file <- function() {
  explicit <- c(
    file.path(data_dir, "STRING_edges.tsv"),
    file.path(data_dir, "STRING_edges.csv"),
    file.path(data_dir, "string_edges.tsv"),
    file.path(data_dir, "string_edges.csv")
  )
  hits <- explicit[file.exists(explicit)]
  if (length(hits)) return(hits[1L])
  fuzzy <- list.files(data_dir, pattern = "STRING.*\\.(tsv|csv|txt)$", full.names = TRUE, ignore.case = TRUE)
  if (length(fuzzy)) fuzzy[1L] else ""
}

locate_string_action_file <- function() {
  explicit <- c(
    file.path(data_dir, "STRING_actions.tsv"),
    file.path(data_dir, "STRING_actions.csv"),
    file.path(data_dir, "string_actions.tsv"),
    file.path(data_dir, "string_actions.csv"),
    file.path(data_dir, "STRING_action.tsv"),
    file.path(data_dir, "STRING_action.csv")
  )
  hits <- explicit[file.exists(explicit)]
  if (length(hits)) return(hits[1L])
  fuzzy <- list.files(data_dir, pattern = "STRING.*action.*\\.(tsv|csv|txt)$", full.names = TRUE, ignore.case = TRUE)
  if (length(fuzzy)) fuzzy[1L] else ""
}

validation_source_status <- function(symbols = character(), allow_live_string = FALSE) {
  string_file <- locate_string_file()
  string_cache <- read_cached_rds("string_network", c(symbols, if (nzchar(string_file) && file.exists(string_file)) paste(basename(string_file), as.numeric(file.info(string_file)$mtime), file.info(string_file)$size) else "no_local_string"), refresh = FALSE)
  list(
    string_local = nzchar(string_file),
    string_local_path = string_file,
    string_cache = !is.null(string_cache) && nrow(string_cache) > 0L,
    string_live_enabled = isTRUE(allow_live_string)
  )
}

string_cache_token <- function(path = "") {
  if (!nzchar(path) || !file.exists(path)) return("no_local_string")
  info <- file.info(path)
  paste(basename(path), as.numeric(info$mtime), info$size)
}

format_bytes <- function(x) {
  if (!is.finite(x) || is.na(x)) return("unknown size")
  units <- c("B", "KB", "MB", "GB", "TB")
  idx <- 1L
  while (x >= 1024 && idx < length(units)) {
    x <- x / 1024
    idx <- idx + 1L
  }
  sprintf("%.1f %s", x, units[idx])
}

format_duration_short <- function(seconds) {
  if (!is.finite(seconds) || is.na(seconds) || seconds < 0) return("ETA unknown")
  if (seconds < 60) return(sprintf("ETA %.0fs", seconds))
  mins <- floor(seconds / 60)
  secs <- round(seconds %% 60)
  sprintf("ETA %dm %02ds", mins, secs)
}

download_url_with_progress <- function(url, destfile, progress = NULL, chunk_bytes = 1024L * 256L) {
  headers <- tryCatch(utils::curlGetHeaders(url, redirect = TRUE), error = function(...) "")
  content_length <- NA_real_
  if (nzchar(headers)) {
    header_lines <- unlist(strsplit(headers, "\r?\n"))
    cl_line <- header_lines[grepl("^content-length\\s*:", header_lines, ignore.case = TRUE)]
    if (length(cl_line)) {
      content_length <- suppressWarnings(as.numeric(trimws(sub("^[^:]+:", "", cl_line[[length(cl_line)]]))))
    }
  }
  con <- url(url, open = "rb")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  out <- file(destfile, open = "wb")
  on.exit(try(close(out), silent = TRUE), add = TRUE)
  start_time <- Sys.time()
  downloaded <- 0
  last_report <- start_time
  repeat {
    buf <- readBin(con, what = "raw", n = chunk_bytes)
    if (!length(buf)) break
    writeBin(buf, out)
    downloaded <- downloaded + length(buf)
    now <- Sys.time()
    elapsed <- as.numeric(difftime(now, start_time, units = "secs"))
    speed <- if (elapsed > 0) downloaded / elapsed else NA_real_
    eta <- if (is.finite(content_length) && is.finite(speed) && speed > 0) (content_length - downloaded) / speed else NA_real_
    if (!is.null(progress) && as.numeric(difftime(now, last_report, units = "secs")) >= 0.2) {
      frac <- if (is.finite(content_length) && content_length > 0) min(downloaded / content_length, 1) else 0.5
      detail <- paste(
        "Downloading STRING network",
        sprintf("Downloaded %s", format_bytes(downloaded)),
        if (is.finite(content_length)) sprintf("of %s", format_bytes(content_length)) else "of unknown total size",
        if (is.finite(speed)) sprintf("at %s/s", format_bytes(speed)) else "at unknown speed",
        format_duration_short(eta),
        sep = " | "
      )
      progress(frac, detail)
      last_report <- now
    }
  }
  invisible(destfile)
}

write_cellregulondb_bundle <- function(symbols, out_dir = data_dir) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  symbols <- sort(unique(symbols[!is.na(symbols) & nzchar(symbols)]))
  gene_file <- file.path(out_dir, "CellRegulonDB_query_genes.txt")
  instr_file <- file.path(out_dir, "CellRegulonDB_download_instructions.txt")
  writeLines(symbols, gene_file, useBytes = TRUE)
  writeLines(
    c(
      "CellRegulonDB export guide",
      "",
      "1. Open https://www.cellregulondb.org/network.html",
      "2. Paste the genes from CellRegulonDB_query_genes.txt into the gene filter.",
      "3. Add tissue or cell-type filters if you know the most relevant biological context.",
      "4. Run the query and download the resulting CSV.",
      "5. Save that file into the data/ folder as CellRegulonDB_edges.csv or CellRegulonDB_regulons.csv.",
      "",
      "The app expects a source TF/regulator column and a target-gene column.",
      "Optional tissue, cell-type, and regulon columns will be shown in the SEMgraph validation context."
    ),
    instr_file,
    useBytes = TRUE
  )
  list(gene_file = gene_file, instruction_file = instr_file)
}

query_cellregulondb_edges <- function(symbols, refresh = FALSE) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (!length(symbols)) {
    return(data.frame(
      source = character(),
      target = character(),
      cellregulondb_context = character(),
      stringsAsFactors = FALSE
    ))
  }
  disk_cached <- read_cached_rds("cellregulondb", symbols, refresh = refresh)
  if (!is.null(disk_cached)) return(disk_cached)
  path <- locate_cellregulondb_file()
  if (!nzchar(path)) {
    return(data.frame(
      source = character(),
      target = character(),
      cellregulondb_context = character(),
      stringsAsFactors = FALSE
    ))
  }
  tbl <- read_hcdt_table(path)
  if (is.null(tbl) || nrow(tbl) == 0L) {
    return(data.frame(
      source = character(),
      target = character(),
      cellregulondb_context = character(),
      stringsAsFactors = FALSE
    ))
  }
  nms <- tolower(names(tbl))
  src_col <- names(tbl)[match(TRUE, nms %in% c("tf", "regulator", "source", "tf_symbol", "regulator_symbol"), nomatch = 0L)]
  tgt_col <- names(tbl)[match(TRUE, nms %in% c("target", "target_gene", "gene", "target_symbol"), nomatch = 0L)]
  tissue_col <- names(tbl)[match(TRUE, grepl("tissue|organ", nms), nomatch = 0L)]
  celltype_col <- names(tbl)[match(TRUE, grepl("cell.?type", nms), nomatch = 0L)]
  regulon_col <- names(tbl)[match(TRUE, grepl("regulon|network|module", nms), nomatch = 0L)]
  if (length(src_col) == 0L || length(tgt_col) == 0L) {
    return(data.frame(
      source = character(),
      target = character(),
      cellregulondb_context = character(),
      stringsAsFactors = FALSE
    ))
  }
  keep <- tbl[[src_col]] %in% symbols & tbl[[tgt_col]] %in% symbols
  sub <- tbl[keep, , drop = FALSE]
  if (!nrow(sub)) {
    return(data.frame(
      source = character(),
      target = character(),
      cellregulondb_context = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- data.frame(
    source = as.character(sub[[src_col]]),
    target = as.character(sub[[tgt_col]]),
    cellregulondb_context = trimws(paste(
      if (length(tissue_col)) paste0("tissue=", as.character(sub[[tissue_col]])) else NULL,
      if (length(celltype_col)) paste0("cell_type=", as.character(sub[[celltype_col]])) else NULL,
      if (length(regulon_col)) paste0("regulon=", as.character(sub[[regulon_col]])) else NULL,
      sep = "; "
    )),
    stringsAsFactors = FALSE
  )
  out$cellregulondb_context[!nzchar(out$cellregulondb_context)] <- "CellRegulonDB local TF-target support"
  out <- out[!duplicated(out[, c("source", "target")]), , drop = FALSE]
  write_cached_rds("cellregulondb", symbols, out)
  out
}

query_string_network <- function(symbols, refresh = FALSE, allow_live = FALSE) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (length(symbols) < 2L) {
    return(data.frame(
      gene_a = character(),
      gene_b = character(),
      string_score = numeric(),
      string_source = character(),
      stringsAsFactors = FALSE
    ))
  }
  local_files <- c(
    file.path(data_dir, "STRING_edges.tsv"),
    file.path(data_dir, "STRING_edges.csv"),
    file.path(data_dir, "string_edges.tsv"),
    file.path(data_dir, "string_edges.csv")
  )
  local_hits <- local_files[file.exists(local_files)]
  local_path <- if (length(local_hits)) local_hits[1L] else locate_string_file()
  cache_key <- c(symbols, string_cache_token(local_path), if (allow_live) "live_allowed" else "live_disabled")
  disk_cached <- read_cached_rds("string_network", cache_key, refresh = refresh)
  if (!is.null(disk_cached)) return(disk_cached)
  if (nzchar(local_path)) {
    tbl <- read_hcdt_table(local_path)
    if (!is.null(tbl) && nrow(tbl) > 0L) {
      nms <- tolower(names(tbl))
      a_col <- names(tbl)[match(TRUE, nms %in% c("preferredname_a", "protein1", "gene_a", "source"), nomatch = 0L)]
      b_col <- names(tbl)[match(TRUE, nms %in% c("preferredname_b", "protein2", "gene_b", "target"), nomatch = 0L)]
      score_col <- names(tbl)[match(TRUE, nms %in% c("score", "combined_score", "string_score"), nomatch = 0L)]
      if (length(a_col) && length(b_col) && length(score_col)) {
        sub <- tbl[tbl[[a_col]] %in% symbols & tbl[[b_col]] %in% symbols, , drop = FALSE]
        out <- data.frame(
          gene_a = as.character(sub[[a_col]]),
          gene_b = as.character(sub[[b_col]]),
          string_score = suppressWarnings(as.numeric(sub[[score_col]])),
          string_source = rep("STRING local network export", nrow(sub)),
          stringsAsFactors = FALSE
        )
        out$string_score[is.na(out$string_score)] <- 0
        write_cached_rds("string_network", cache_key, out)
        return(out)
      }
    }
  }
  if (!allow_live) {
    return(data.frame(
      gene_a = character(),
      gene_b = character(),
      string_score = numeric(),
      string_source = character(),
      stringsAsFactors = FALSE
    ))
  }
  endpoint <- paste0(
    "https://string-db.org/api/tsv/network?identifiers=",
    utils::URLencode(paste(symbols, collapse = "\r"), reserved = TRUE),
    "&species=9606&required_score=0&network_type=functional"
  )
  out <- tryCatch({
    tbl <- utils::read.delim(endpoint, stringsAsFactors = FALSE, check.names = FALSE)
    if (is.null(tbl) || !nrow(tbl)) {
      data.frame(gene_a = character(), gene_b = character(), string_score = numeric(), string_source = character(), stringsAsFactors = FALSE)
    } else {
      a_col <- names(tbl)[match(TRUE, tolower(names(tbl)) %in% c("preferredname_a", "preferrednamea"), nomatch = 0L)]
      b_col <- names(tbl)[match(TRUE, tolower(names(tbl)) %in% c("preferredname_b", "preferrednameb"), nomatch = 0L)]
      score_col <- names(tbl)[match(TRUE, tolower(names(tbl)) %in% c("score", "combined_score"), nomatch = 0L)]
      if (length(a_col) && length(b_col) && length(score_col)) {
        sub <- tbl[tbl[[a_col]] %in% symbols & tbl[[b_col]] %in% symbols, , drop = FALSE]
        data.frame(
          gene_a = as.character(sub[[a_col]]),
          gene_b = as.character(sub[[b_col]]),
          string_score = suppressWarnings(as.numeric(sub[[score_col]])),
          string_source = rep("STRING live network API", nrow(sub)),
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(gene_a = character(), gene_b = character(), string_score = numeric(), string_source = character(), stringsAsFactors = FALSE)
      }
    }
  }, error = function(...) {
    data.frame(gene_a = character(), gene_b = character(), string_score = numeric(), string_source = character(), stringsAsFactors = FALSE)
  })
  if (nrow(out)) {
    out$string_score[is.na(out$string_score)] <- 0
    write_cached_rds("string_network", cache_key, out)
  }
  out
}

query_string_actions <- function(symbols, refresh = FALSE) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (length(symbols) < 2L) {
    return(data.frame(
      gene_a = character(),
      gene_b = character(),
      string_action_mode = character(),
      string_action_effect = character(),
      a_is_acting = logical(),
      string_action_source = character(),
      stringsAsFactors = FALSE
    ))
  }
  candidate_paths <- unique(c(locate_string_action_file(), locate_string_file()))
  candidate_paths <- candidate_paths[nzchar(candidate_paths) & file.exists(candidate_paths)]
  cache_key <- c(symbols, if (length(candidate_paths)) vapply(candidate_paths, string_cache_token, character(1)) else "no_action_file")
  disk_cached <- read_cached_rds("string_actions", cache_key, refresh = refresh)
  if (!is.null(disk_cached)) return(disk_cached)
  for (local_path in candidate_paths) {
    tbl <- read_hcdt_table(local_path)
    if (is.null(tbl) || !nrow(tbl)) next
    nms <- tolower(names(tbl))
    a_col <- names(tbl)[match(TRUE, nms %in% c("preferredname_a", "preferrednamea", "item_id_a", "protein1", "gene_a", "source"), nomatch = 0L)]
    b_col <- names(tbl)[match(TRUE, nms %in% c("preferredname_b", "preferrednameb", "item_id_b", "protein2", "gene_b", "target"), nomatch = 0L)]
    mode_col <- names(tbl)[match(TRUE, nms %in% c("mode", "interaction_mode", "string_mode"), nomatch = 0L)]
    effect_col <- names(tbl)[match(TRUE, nms %in% c("action", "effect", "string_action"), nomatch = 0L)]
    acting_col <- names(tbl)[match(TRUE, nms %in% c("a_is_acting", "is_directed", "direction", "acting"), nomatch = 0L)]
    if (!length(a_col) || !length(b_col) || (!length(mode_col) && !length(effect_col))) next
    sub <- tbl[tbl[[a_col]] %in% symbols & tbl[[b_col]] %in% symbols, , drop = FALSE]
    if (!nrow(sub)) next
    acting_vals <- if (length(acting_col)) {
      raw <- tolower(trimws(as.character(sub[[acting_col]])))
      raw %in% c("t", "true", "1", "yes", "a", "a->b", "forward")
    } else {
      rep(NA, nrow(sub))
    }
    out <- data.frame(
      gene_a = as.character(sub[[a_col]]),
      gene_b = as.character(sub[[b_col]]),
      string_action_mode = if (length(mode_col)) as.character(sub[[mode_col]]) else "unspecified",
      string_action_effect = if (length(effect_col)) as.character(sub[[effect_col]]) else "unspecified",
      a_is_acting = acting_vals,
      string_action_source = rep(sprintf("STRING action export: %s", basename(local_path)), nrow(sub)),
      stringsAsFactors = FALSE
    )
    out$string_action_mode[is.na(out$string_action_mode) | !nzchar(out$string_action_mode)] <- "unspecified"
    out$string_action_effect[is.na(out$string_action_effect) | !nzchar(out$string_action_effect)] <- "unspecified"
    out <- out[!duplicated(out), , drop = FALSE]
    write_cached_rds("string_actions", cache_key, out)
    return(out)
  }

  data.frame(
    gene_a = character(),
    gene_b = character(),
    string_action_mode = character(),
    string_action_effect = character(),
    a_is_acting = logical(),
    string_action_source = character(),
    stringsAsFactors = FALSE
  )
}

build_edge_validation_table <- function(graph, refresh = FALSE, use_live_string = FALSE) {
  edge_df <- igraph::as_data_frame(graph, what = "edges")
  if (is.null(edge_df) || !nrow(edge_df)) {
    return(empty_edge_validation())
  }
  edges <- data.frame(source = edge_df$from, target = edge_df$to, stringsAsFactors = FALSE)
  symbols <- unique(c(edges$source, edges$target))
  string_df <- query_string_network(symbols, refresh = refresh, allow_live = use_live_string)
  action_df <- query_string_actions(symbols, refresh = refresh)

  out <- empty_edge_validation(edges)
  if (nrow(string_df)) {
    string_df$edge_key <- normalize_edge_key(string_df$gene_a, string_df$gene_b)
    best_string <- stats::aggregate(string_score ~ edge_key, data = string_df, FUN = max)
    src_string <- stats::aggregate(string_source ~ edge_key, data = string_df, FUN = function(x) paste(unique(x), collapse = "; "))
    best_string <- merge(best_string, src_string, by = "edge_key", all.x = TRUE, sort = FALSE)
    out$edge_key <- normalize_edge_key(out$source, out$target)
    out$string_score <- best_string$string_score[match(out$edge_key, best_string$edge_key)]
    out$string_supported <- !is.na(out$string_score) & out$string_score > 0
    out$string_score[is.na(out$string_score)] <- 0
    out$string_source <- best_string$string_source[match(out$edge_key, best_string$edge_key)]
  } else {
    out$string_source <- "STRING unavailable"
  }

  if (nrow(action_df)) {
    action_df$edge_key <- normalize_edge_key(action_df$gene_a, action_df$gene_b)
    out$edge_key <- normalize_edge_key(out$source, out$target)
    action_split <- split(action_df, action_df$edge_key)
    for (i in seq_len(nrow(out))) {
      cur <- action_split[[out$edge_key[i]]]
      if (is.null(cur) || !nrow(cur)) next
      out$string_action_supported[i] <- TRUE
      out$string_action_mode[i] <- paste(unique(cur$string_action_mode[nzchar(cur$string_action_mode)]), collapse = "; ")
      if (!nzchar(out$string_action_mode[i])) out$string_action_mode[i] <- "unspecified"
      out$string_action_effect[i] <- paste(unique(cur$string_action_effect[nzchar(cur$string_action_effect)]), collapse = "; ")
      if (!nzchar(out$string_action_effect[i])) out$string_action_effect[i] <- "unspecified"
      out$string_action_source[i] <- paste(unique(cur$string_action_source[nzchar(cur$string_action_source)]), collapse = "; ")
      directed_idx <- which(!is.na(cur$a_is_acting))
      if (length(directed_idx)) {
        out$string_direction_available[i] <- TRUE
        dir_match <- vapply(directed_idx, function(j) {
          if (isTRUE(cur$a_is_acting[j])) {
            identical(cur$gene_a[j], out$source[i]) && identical(cur$gene_b[j], out$target[i])
          } else {
            identical(cur$gene_b[j], out$source[i]) && identical(cur$gene_a[j], out$target[i])
          }
        }, logical(1))
        out$string_direction_matches[i] <- any(dir_match, na.rm = TRUE)
      }
    }
  } else {
    out$string_action_source <- "STRING action data unavailable"
  }

  out$validation_tier <- ifelse(
    out$string_direction_matches,
    "STRING action direction agrees",
    ifelse(
      out$string_action_supported & out$string_direction_available,
      "STRING action available but direction differs/unclear",
      ifelse(out$string_supported, "STRING pair supported", "No external support")
    )
  )
  out$validation_score <- out$string_score
  out$edge_key <- NULL
  out
}

filter_validated_edges <- function(edge_df, mode = "all", string_cutoff = 0.7) {
  if (is.null(edge_df) || !nrow(edge_df)) return(edge_df)
  keep <- switch(
    mode,
    any_support = edge_df$string_supported,
    string_only = edge_df$string_supported & edge_df$string_score >= string_cutoff,
    rep(TRUE, nrow(edge_df))
  )
  edge_df[keep, , drop = FALSE]
}

annotate_semgraph_edges <- function(graph, res) {
  edge_df <- igraph::as_data_frame(graph, what = "edges")
  if (is.null(edge_df) || !nrow(edge_df)) return(data.frame())
  names(edge_df)[1:2] <- c("source", "target")
  valid <- res$edge_validation %||% empty_edge_validation()
  if (nrow(valid)) {
    edge_df <- merge(edge_df, valid, by = c("source", "target"), all.x = TRUE, sort = FALSE)
  }
  edge_df$string_supported <- ifelse(is.na(edge_df$string_supported), FALSE, edge_df$string_supported)
  edge_df$string_score <- ifelse(is.na(edge_df$string_score), 0, edge_df$string_score)
  edge_df$string_source <- ifelse(is.na(edge_df$string_source), "STRING unavailable", edge_df$string_source)
  edge_df$string_action_supported <- ifelse(is.na(edge_df$string_action_supported), FALSE, edge_df$string_action_supported)
  edge_df$string_action_mode <- ifelse(is.na(edge_df$string_action_mode), "No STRING action annotation", edge_df$string_action_mode)
  edge_df$string_action_effect <- ifelse(is.na(edge_df$string_action_effect), "No STRING effect annotation", edge_df$string_action_effect)
  edge_df$string_action_source <- ifelse(is.na(edge_df$string_action_source), "STRING action data unavailable", edge_df$string_action_source)
  edge_df$string_direction_available <- ifelse(is.na(edge_df$string_direction_available), FALSE, edge_df$string_direction_available)
  edge_df$string_direction_matches <- ifelse(is.na(edge_df$string_direction_matches), FALSE, edge_df$string_direction_matches)
  edge_df$validation_tier <- ifelse(is.na(edge_df$validation_tier), "No external support", edge_df$validation_tier)
  edge_df$validation_score <- ifelse(is.na(edge_df$validation_score), 0, edge_df$validation_score)
  edge_df
}

summarize_node_validation <- function(edge_df) {
  if (is.null(edge_df) || !nrow(edge_df)) {
    return(data.frame(
      gene = character(),
      externally_supported_edges = numeric(),
      string_supported_edges = numeric(),
      string_direction_supported_edges = numeric(),
      best_string_score = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  support_rows <- lapply(seq_len(nrow(edge_df)), function(i) {
    data.frame(
      gene = c(edge_df$source[i], edge_df$target[i]),
      externally_supported_edges = c(as.numeric(edge_df$string_supported[i]), as.numeric(edge_df$string_supported[i])),
      string_supported_edges = c(as.numeric(edge_df$string_supported[i]), as.numeric(edge_df$string_supported[i])),
      string_direction_supported_edges = c(as.numeric(edge_df$string_direction_matches[i]), as.numeric(edge_df$string_direction_matches[i])),
      best_string_score = c(edge_df$string_score[i], edge_df$string_score[i]),
      stringsAsFactors = FALSE
    )
  })
  support_df <- do.call(rbind, support_rows)
  genes <- unique(support_df$gene)
  out <- data.frame(
    gene = genes,
    externally_supported_edges = numeric(length(genes)),
    string_supported_edges = numeric(length(genes)),
    string_direction_supported_edges = numeric(length(genes)),
    best_string_score = numeric(length(genes)),
    stringsAsFactors = FALSE
  )
  for (i in seq_along(genes)) {
    sub <- support_df[support_df$gene == genes[i], , drop = FALSE]
    out$externally_supported_edges[i] <- sum(sub$externally_supported_edges, na.rm = TRUE)
    out$string_supported_edges[i] <- sum(sub$string_supported_edges, na.rm = TRUE)
    out$string_direction_supported_edges[i] <- sum(sub$string_direction_supported_edges, na.rm = TRUE)
    out$best_string_score[i] <- if (all(is.na(sub$best_string_score))) 0 else max(sub$best_string_score, na.rm = TRUE)
  }
  out
}

query_reactome_drug_context <- function(symbols, gene_map = NULL, refresh = FALSE, progress = NULL, allow_live = FALSE) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (!length(symbols)) return(empty_drug_annotation(character()))
  disk_cached <- read_cached_rds("reactome", symbols, refresh = refresh)
  if (!is.null(disk_cached)) return(disk_cached)
  if (!is.null(progress)) progress(0.83, "Loading local Reactome mapping cache")
  reactome_tbl <- load_reactome_mapping(force_refresh = refresh, progress = progress)
  if (!is.null(reactome_tbl) && !is.null(gene_map)) {
    gm <- gene_map[gene_map$hgnc_symbol %in% symbols, c("hgnc_symbol", "entrez_id"), drop = FALSE]
    gm <- gm[!is.na(gm$entrez_id) & nzchar(gm$entrez_id), , drop = FALSE]
    gm <- gm[!duplicated(gm$hgnc_symbol), , drop = FALSE]
    merged <- merge(gm, reactome_tbl, by = "entrez_id", all.x = TRUE, sort = FALSE)
    rows <- lapply(symbols, function(sym) {
      cur <- merged[merged$hgnc_symbol == sym, , drop = FALSE]
      pathway_names <- unique(cur$event_name[!is.na(cur$event_name) & nzchar(cur$event_name)])
      pathway_names <- pathway_names[grepl("drug|inhibit|agon|antagon|compound", pathway_names, ignore.case = TRUE)]
      if (!length(pathway_names)) {
        return(data.frame(
          hgnc_symbol = sym,
          is_druggable = FALSE,
          drug_mode = "context only",
          drug_summary = "No Reactome drug-mediated pathway context found",
          drug_count = 0L,
          drug_source = "Reactome",
          drug_source_details = "Reactome local mapping: no drug-mediated pathway context found",
          reactome_context = "No Reactome drug-mediated pathway context found",
          stringsAsFactors = FALSE
        ))
      }
      context_txt <- paste(head(pathway_names, 4L), collapse = "; ")
      data.frame(
        hgnc_symbol = sym,
        is_druggable = TRUE,
        drug_mode = "context only",
        drug_summary = context_txt,
        drug_count = length(pathway_names),
        drug_source = "Reactome",
        drug_source_details = paste0("Reactome local mapping context: ", context_txt),
        reactome_context = context_txt,
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, rows)
    write_cached_rds("reactome", symbols, out)
    return(out)
  }
  if (!allow_live) {
    out <- empty_drug_annotation(symbols)
    out$drug_source <- "Reactome"
    out$drug_mode <- "cache only"
    out$drug_source_details <- "Reactome local mapping unavailable and live lookup skipped"
    out$reactome_context <- "Reactome local mapping unavailable and live lookup skipped"
    return(out)
  }
  total <- length(symbols)
  rows <- lapply(seq_along(symbols), function(i) {
    sym <- symbols[i]
    if (!is.null(progress) && (i == 1L || i %% 10L == 0L || i == total)) {
      progress(0.85, sprintf("Loading Reactome pathway context (%s/%s genes)", i, total))
    }
    endpoint <- paste0(
      "https://reactome.org/ContentService/search/query?query=",
      utils::URLencode(sym, reserved = TRUE),
      "&species=Homo%20sapiens&types=Pathway&cluster=true"
    )
    raw <- suppressWarnings(tryCatch(jsonlite::fromJSON(endpoint, simplifyVector = FALSE), error = function(...) NULL))
    pathway_names <- character()
    recurse_names <- function(x) {
      if (is.list(x)) {
        if (!is.null(x$name) && is.character(x$name) && length(x$name) == 1L) {
          pathway_names <<- c(pathway_names, x$name)
        }
        lapply(x, recurse_names)
      }
      NULL
    }
    recurse_names(raw)
    pathway_names <- unique(pathway_names[grepl("drug|inhibit|agon|antagon|compound", pathway_names, ignore.case = TRUE)])
    if (!length(pathway_names)) {
      return(data.frame(
        hgnc_symbol = sym,
        is_druggable = FALSE,
        drug_mode = "context only",
        drug_summary = "No Reactome drug-mediated pathway context found",
        drug_count = 0L,
        drug_source = "Reactome",
        drug_source_details = "Reactome: no drug-mediated pathway context found",
        reactome_context = "No Reactome drug-mediated pathway context found",
        stringsAsFactors = FALSE
      ))
    }
    context_txt <- paste(head(pathway_names, 4L), collapse = "; ")
    data.frame(
      hgnc_symbol = sym,
      is_druggable = TRUE,
      drug_mode = "context only",
      drug_summary = context_txt,
      drug_count = length(pathway_names),
      drug_source = "Reactome",
      drug_source_details = paste0("Reactome drug-mediated pathway context: ", context_txt),
      reactome_context = context_txt,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  write_cached_rds("reactome", symbols, out)
  out
}

aggregate_drug_annotations <- function(symbols, gene_map = NULL, refresh = FALSE, progress = NULL, reactome_live = FALSE) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (!length(symbols)) return(empty_drug_annotation(character()))
  disk_cached <- read_cached_rds("drug_agg", symbols, refresh = refresh)
  if (!is.null(disk_cached)) return(disk_cached)
  if (!is.null(progress)) progress(0.76, sprintf("Loading DGIdb drug-target annotations for %s genes", length(symbols)))
  dgidb <- query_dgidb_targets(symbols, refresh = refresh)
  if (!is.null(progress)) progress(0.80, sprintf("Loading HCDT drug-target annotations for %s genes", length(symbols)))
  hcdt <- query_hcdt_targets(symbols, refresh = refresh)
  if (!is.null(progress)) progress(0.82, sprintf("Loading Reactome pathway context for %s genes", length(symbols)))
  reactome <- query_reactome_drug_context(symbols, gene_map = gene_map, refresh = refresh, progress = progress, allow_live = reactome_live)
  rows <- lapply(symbols, function(sym) {
    srcs <- list(
      DGIdb = dgidb[dgidb$hgnc_symbol == sym, , drop = FALSE],
      HCDT = hcdt[hcdt$hgnc_symbol == sym, , drop = FALSE],
      Reactome = reactome[reactome$hgnc_symbol == sym, , drop = FALSE]
    )
    is_hit <- vapply(srcs, function(df) nrow(df) > 0L && isTRUE(df$is_druggable[1L]), logical(1))
    active_sources <- names(srcs)[is_hit]
    mode_candidates <- unlist(lapply(srcs[is_hit], function(df) df$drug_mode[1L]), use.names = FALSE)
    summary_candidates <- unlist(lapply(names(srcs)[is_hit], function(nm) {
      df <- srcs[[nm]]
      if (!nrow(df)) return(character())
      paste0(nm, ": ", df$drug_summary[1L])
    }), use.names = FALSE)
    detail_candidates <- unlist(lapply(names(srcs), function(nm) {
      df <- srcs[[nm]]
      if (!nrow(df)) return(character())
      df$drug_source_details[1L]
    }), use.names = FALSE)
    reactome_context <- if (nrow(srcs$Reactome)) srcs$Reactome$reactome_context[1L] else "No Reactome drug-mediated pathway context found"
    chosen_mode <- {
      preferred <- mode_candidates[mode_candidates %in% c("inhibiting", "activating", "mixed", "other or unclear")]
      if (length(preferred)) preferred[1L] else if (length(mode_candidates)) mode_candidates[1L] else "no interaction found"
    }
    drug_counts <- sum(vapply(srcs, function(df) if (nrow(df)) df$drug_count[1L] else 0L, integer(1)), na.rm = TRUE)
    data.frame(
      hgnc_symbol = sym,
      is_druggable = any(is_hit),
      drug_mode = chosen_mode,
      drug_summary = if (length(summary_candidates)) paste(summary_candidates, collapse = " | ") else "No external drug-target interaction found",
      drug_count = drug_counts,
      drug_source = if (length(active_sources)) paste(active_sources, collapse = ", ") else "None",
      drug_source_details = if (length(detail_candidates)) paste(detail_candidates, collapse = " | ") else "No DGIdb, HCDT, or Reactome drug context found",
      reactome_context = reactome_context,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  write_cached_rds("drug_agg", symbols, out)
  out
}

run_enrichment <- function(de_res, gene_map, methods, padj_cutoff, lfc_cutoff, min_gs, max_gs, p_adjust_method, top_n, refresh = FALSE) {
  cache_key <- c(
    de_res$label %||% "de",
    sort(methods),
    padj_cutoff, lfc_cutoff, min_gs, max_gs, p_adjust_method, top_n,
    paste(head(de_res$df$gene_id, 200L), collapse = "|")
  )
  disk_cached <- read_cached_rds("pathway_enrichment", cache_key, refresh = refresh)
  if (!is.null(disk_cached)) return(disk_cached)
  joined <- de_res$df
  if (!"entrez_id" %in% names(joined)) {
    joined <- attach_gene_labels(joined, gene_map)
  }
  sig_genes <- joined$entrez_id[!is.na(joined$entrez_id) & !is.na(joined$padj) & joined$padj <= padj_cutoff & abs(joined$log2FoldChange) >= lfc_cutoff]
  sig_genes <- unique(as.character(sig_genes))

  ranked <- joined[!is.na(joined$entrez_id) & !is.na(joined$stat), c("entrez_id", "stat")]
  ranked <- ranked[!duplicated(ranked$entrez_id), , drop = FALSE]
  ranked_vec <- ranked$stat
  names(ranked_vec) <- ranked$entrez_id
  ranked_vec <- sort(ranked_vec, decreasing = TRUE)

  results <- list()
  tables <- list()

  if ("GO ORA" %in% methods && length(sig_genes) > 0L) {
    go_ora <- clusterProfiler::enrichGO(gene = sig_genes, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = p_adjust_method, pvalueCutoff = padj_cutoff, qvalueCutoff = padj_cutoff, minGSSize = min_gs, maxGSSize = max_gs, readable = TRUE)
    tables[["GO ORA"]] <- format_enrichment_df(as.data.frame(go_ora), "GO ORA")
    results[["GO ORA"]] <- go_ora
  }
  if ("KEGG ORA" %in% methods && length(sig_genes) > 0L) {
    kegg_ora <- clusterProfiler::enrichKEGG(gene = sig_genes, organism = "hsa", keyType = "kegg", pAdjustMethod = p_adjust_method, pvalueCutoff = padj_cutoff, qvalueCutoff = padj_cutoff, minGSSize = min_gs, maxGSSize = max_gs)
    tables[["KEGG ORA"]] <- format_enrichment_df(as.data.frame(kegg_ora), "KEGG ORA")
    results[["KEGG ORA"]] <- kegg_ora
  }
  if ("GO GSEA" %in% methods && length(ranked_vec) > 0L) {
    go_gsea <- clusterProfiler::gseGO(geneList = ranked_vec, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = p_adjust_method, minGSSize = min_gs, maxGSSize = max_gs, verbose = FALSE)
    tables[["GO GSEA"]] <- format_enrichment_df(as.data.frame(go_gsea), "GO GSEA")
    results[["GO GSEA"]] <- go_gsea
  }
  if ("KEGG GSEA" %in% methods && length(ranked_vec) > 0L) {
    kegg_gsea <- clusterProfiler::gseKEGG(geneList = ranked_vec, organism = "hsa", keyType = "ncbi-geneid", pAdjustMethod = p_adjust_method, minGSSize = min_gs, maxGSSize = max_gs, verbose = FALSE)
    tables[["KEGG GSEA"]] <- format_enrichment_df(as.data.frame(kegg_gsea), "KEGG GSEA")
    results[["KEGG GSEA"]] <- kegg_gsea
  }

  tables <- tables[vapply(tables, nrow, integer(1)) > 0L]
  combined <- if (length(tables)) bind_rows_fill(tables) else data.frame()
  if (nrow(combined) > 0L) {
    if (!"p.adjust" %in% names(combined)) {
      combined$p.adjust <- NA_real_
    }
    if (!"pathway" %in% names(combined)) {
      combined$pathway <- if ("Description" %in% names(combined)) combined$Description else paste("Pathway", seq_len(nrow(combined)))
    }
    combined$consensus_hits <- ave(combined$pathway, combined$pathway, FUN = length)
    combined$p.adjust <- suppressWarnings(as.numeric(combined$p.adjust))
    combined$consensus_hits <- suppressWarnings(as.numeric(combined$consensus_hits))
    combined <- combined[order(combined$p.adjust, na.last = TRUE), , drop = FALSE]
  }
  consensus <- if (nrow(combined) > 0L) {
    consensus_order <- order(-combined$consensus_hits, combined$p.adjust, na.last = TRUE)
    ordered <- combined[consensus_order, , drop = FALSE]
    consensus_df <- data.frame(
      pathway = ordered$pathway,
      consensus_hits = ordered$consensus_hits,
      p.adjust = ordered$p.adjust,
      stringsAsFactors = FALSE
    )
    consensus_df[!duplicated(consensus_df$pathway), , drop = FALSE]
  } else {
    data.frame()
  }
  method_tables <- lapply(split(combined, combined$method), function(df) {
    df[order(df$p.adjust, na.last = TRUE), , drop = FALSE]
  })

  out <- list(
    sig_gene_count = length(sig_genes),
    ranked_gene_count = length(ranked_vec),
    combined = combined,
    method_tables = method_tables,
    consensus = head(consensus, top_n),
    ranked = ranked_vec
  )
  write_cached_rds("pathway_enrichment", cache_key, out)
  out
}

prepare_wgcna_traits <- function(metadata, selected_cols) {
  if (length(selected_cols) == 0L) {
    return(NULL)
  }
  trait_parts <- list()
  for (col in selected_cols) {
    values <- metadata[[col]]
    if (is.factor(values) || is.character(values)) {
      mm <- stats::model.matrix(~ factor(values) - 1)
      colnames(mm) <- paste0(col, ": ", sub("^factor\\(values\\)", "", colnames(mm)))
      trait_parts[[col]] <- mm
    } else {
      trait_parts[[col]] <- matrix(as.numeric(values), ncol = 1L, dimnames = list(metadata$run_id, col))
    }
  }
  out <- do.call(cbind, trait_parts)
  rownames(out) <- metadata$run_id
  out
}

run_wgcna_analysis <- function(expr_matrix, metadata, trait_cols, power_mode, manual_power, min_module_size, deep_split, merge_cut_height, top_var_genes) {
  dat_expr <- t(expr_matrix)
  gene_var <- apply(dat_expr, 2L, stats::var, na.rm = TRUE)
  keep_order <- order(gene_var, decreasing = TRUE, na.last = NA)
  keep_n <- min(length(keep_order), top_var_genes)
  dat_expr <- dat_expr[, keep_order[seq_len(keep_n)], drop = FALSE]

  gsg <- WGCNA::goodSamplesGenes(dat_expr, verbose = 0)
  if (!gsg$allOK) {
    dat_expr <- dat_expr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
    metadata <- metadata[match(rownames(dat_expr), metadata$run_id), , drop = FALSE]
  } else {
    metadata <- metadata[match(rownames(dat_expr), metadata$run_id), , drop = FALSE]
  }

  candidate_powers <- c(1:10, seq(12, 20, by = 2))
  sft <- WGCNA::pickSoftThreshold(dat_expr, powerVector = candidate_powers, verbose = 0)
  fit_tbl <- sft$fitIndices
  auto_power <- fit_tbl$Power[which(fit_tbl$SFT.R.sq >= 0.8)[1L]]
  if (is.na(auto_power) || length(auto_power) == 0L) {
    auto_power <- fit_tbl$Power[which.max(fit_tbl$SFT.R.sq)]
  }
  selected_power <- if (identical(power_mode, "auto")) auto_power else manual_power

  gene_cor <- stats::cor(dat_expr, use = "pairwise.complete.obs", method = "pearson")
  adjacency <- ((1 + gene_cor) / 2) ^ selected_power
  diag(adjacency) <- 0
  tom <- WGCNA::TOMsimilarity(adjacency, TOMType = "signed")
  diss_tom <- 1 - tom
  gene_tree <- stats::hclust(stats::as.dist(diss_tom), method = "average")

  dynamic_labels <- dynamicTreeCut::cutreeDynamic(
    dendro = gene_tree,
    distM = diss_tom,
    deepSplit = deep_split,
    pamRespectsDendro = FALSE,
    minClusterSize = min_module_size
  )
  dynamic_colors <- WGCNA::labels2colors(dynamic_labels)
  names(dynamic_colors) <- colnames(dat_expr)
  merged <- WGCNA::mergeCloseModules(
    exprData = dat_expr,
    colors = dynamic_colors,
    cutHeight = merge_cut_height,
    verbose = 0
  )

  module_colors <- merged$colors
  names(module_colors) <- colnames(dat_expr)
  eigengenes <- WGCNA::orderMEs(merged$newMEs)
  trait_mat <- prepare_wgcna_traits(metadata, trait_cols)
  module_trait_cor <- NULL
  module_trait_p <- NULL
  if (!is.null(trait_mat) && ncol(trait_mat) > 0L) {
    module_trait_cor <- WGCNA::cor(eigengenes, trait_mat, use = "p")
    module_trait_p <- WGCNA::corPvalueStudent(module_trait_cor, nrow(dat_expr))
  }

  module_sizes <- sort(table(module_colors), decreasing = TRUE)
  gene_modules <- data.frame(
    gene_id = colnames(dat_expr),
    module = module_colors,
    stringsAsFactors = FALSE
  )
  intramod_connectivity <- numeric(length(module_colors))
  names(intramod_connectivity) <- colnames(dat_expr)
  for (mod in unique(module_colors)) {
    mod_genes <- names(module_colors)[module_colors == mod]
    if (length(mod_genes) == 0L) next
    mod_adj <- adjacency[mod_genes, mod_genes, drop = FALSE]
    intramod_connectivity[mod_genes] <- rowSums(mod_adj, na.rm = TRUE)
  }
  gene_modules$kWithin <- intramod_connectivity[gene_modules$gene_id]

  list(
    dat_expr = dat_expr,
    metadata = metadata,
    fit_tbl = fit_tbl,
    selected_power = selected_power,
    auto_power = auto_power,
    module_colors = module_colors,
    module_eigengenes = eigengenes,
    module_trait_cor = module_trait_cor,
    module_trait_p = module_trait_p,
    module_sizes = module_sizes,
    gene_modules = gene_modules,
    adjacency = adjacency,
    dendrogram = gene_tree,
    unmerged_colors = dynamic_colors,
    params = list(min_module_size = min_module_size, deep_split = deep_split, merge_cut_height = merge_cut_height, top_var_genes = top_var_genes)
  )
}

run_module_enrichment <- function(module_genes, gene_map, padj_cutoff = 0.05, min_gs = 10, max_gs = 500, p_adjust_method = "BH", top_n = 15, refresh = FALSE, cache_label = "module") {
  cache_key <- c(cache_label, sort(module_genes), padj_cutoff, min_gs, max_gs, p_adjust_method, top_n)
  disk_cached <- read_cached_rds("module_enrichment", cache_key, refresh = refresh)
  if (!is.null(disk_cached)) return(disk_cached)
  mod_map <- gene_map[match(module_genes, gene_map$gene_id), , drop = FALSE]
  entrez <- unique(as.character(mod_map$entrez_id[!is.na(mod_map$entrez_id) & nzchar(mod_map$entrez_id)]))
  if (length(entrez) == 0L) {
    return(list(go = data.frame(), kegg = data.frame(), combined = data.frame(), consensus = data.frame(), gene_count = 0L))
  }

  go <- tryCatch(
    clusterProfiler::enrichGO(
      gene = entrez,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = p_adjust_method,
      pvalueCutoff = padj_cutoff,
      qvalueCutoff = padj_cutoff,
      minGSSize = min_gs,
      maxGSSize = max_gs,
      readable = TRUE
    ),
    error = function(...) NULL
  )
  kegg <- tryCatch(
    clusterProfiler::enrichKEGG(
      gene = entrez,
      organism = "hsa",
      keyType = "kegg",
      pAdjustMethod = p_adjust_method,
      pvalueCutoff = padj_cutoff,
      qvalueCutoff = padj_cutoff,
      minGSSize = min_gs,
      maxGSSize = max_gs
    ),
    error = function(...) NULL
  )

  go_df <- format_enrichment_df(as.data.frame(go), "Module GO ORA")
  kegg_df <- format_enrichment_df(as.data.frame(kegg), "Module KEGG ORA")
  dfs <- Filter(function(x) nrow(x) > 0L, list(go_df, kegg_df))
  combined <- if (length(dfs)) bind_rows_fill(dfs) else data.frame()
  if (nrow(combined) > 0L) {
    if (!"p.adjust" %in% names(combined)) combined$p.adjust <- NA_real_
    if (!"pathway" %in% names(combined)) {
      combined$pathway <- if ("Description" %in% names(combined)) combined$Description else paste("Pathway", seq_len(nrow(combined)))
    }
    combined$p.adjust <- suppressWarnings(as.numeric(combined$p.adjust))
    combined$consensus_hits <- ave(combined$pathway, combined$pathway, FUN = length)
    combined$consensus_hits <- suppressWarnings(as.numeric(combined$consensus_hits))
    combined <- combined[order(combined$p.adjust, na.last = TRUE), , drop = FALSE]
  }
  consensus <- if (nrow(combined) > 0L) {
    ord <- order(-combined$consensus_hits, combined$p.adjust, na.last = TRUE)
    tmp <- combined[ord, c("pathway", "consensus_hits", "p.adjust"), drop = FALSE]
    tmp[!duplicated(tmp$pathway), , drop = FALSE]
  } else {
    data.frame()
  }

  out <- list(
    go = go_df,
    kegg = kegg_df,
    combined = combined,
    consensus = head(consensus, top_n),
    gene_count = length(entrez),
    has_signal = nrow(combined) > 0L
  )
  write_cached_rds("module_enrichment", cache_key, out)
  out
}

infer_module_de_direction <- function(wgcna_res, de_df = NULL, alpha = 0.05, lfc_cutoff = 0) {
  mods <- sort(unique(wgcna_res$gene_modules$module))
  mods <- mods[mods != "grey"]
  if (is.null(de_df) || !nrow(de_df)) {
    return(data.frame(module = mods, mean_log2fc = NA_real_, up_fraction = NA_real_, down_fraction = NA_real_, direction = "unknown", stringsAsFactors = FALSE))
  }
  do.call(rbind, lapply(mods, function(mod) {
    genes <- wgcna_res$gene_modules$gene_id[wgcna_res$gene_modules$module == mod]
    dd <- de_df[match(genes, de_df$gene_id), , drop = FALSE]
    mean_fc <- mean(dd$log2FoldChange, na.rm = TRUE)
    up_fraction <- mean(deg_mask(dd, alpha = alpha, lfc_cutoff = lfc_cutoff) & dd$log2FoldChange > 0, na.rm = TRUE)
    down_fraction <- mean(deg_mask(dd, alpha = alpha, lfc_cutoff = lfc_cutoff) & dd$log2FoldChange < 0, na.rm = TRUE)
    direction <- if (!is.finite(mean_fc)) "unknown" else if (mean_fc > 0) "up" else if (mean_fc < 0) "down" else "flat"
    data.frame(module = mod, mean_log2fc = mean_fc, up_fraction = up_fraction, down_fraction = down_fraction, direction = direction, stringsAsFactors = FALSE)
  }))
}

suggest_module_programs <- function(wgcna_res, de_df = NULL, alpha = 0.05, lfc_cutoff = 0) {
  me <- wgcna_res$module_eigengenes
  if (is.null(me) || ncol(me) < 2L) return(data.frame())
  me_cols <- colnames(me)
  modules <- sub("^ME", "", me_cols)
  keep <- modules != "grey"
  me <- me[, keep, drop = FALSE]
  modules <- modules[keep]
  if (length(modules) < 2L) return(data.frame())
  cor_mat <- stats::cor(me, use = "pairwise.complete.obs")
  de_dir <- infer_module_de_direction(wgcna_res, de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
  top_trait <- NULL
  top_trait_sign <- NULL
  if (!is.null(wgcna_res$module_trait_cor) && ncol(wgcna_res$module_trait_cor) > 0L) {
    top_idx <- apply(abs(wgcna_res$module_trait_cor), 1, which.max)
    top_trait <- colnames(wgcna_res$module_trait_cor)[top_idx]
    top_trait_sign <- sign(wgcna_res$module_trait_cor[cbind(seq_len(nrow(wgcna_res$module_trait_cor)), top_idx)])
    names(top_trait) <- sub("^ME", "", rownames(wgcna_res$module_trait_cor))
    names(top_trait_sign) <- sub("^ME", "", rownames(wgcna_res$module_trait_cor))
  }
  pairs <- utils::combn(modules, 2L, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    m1 <- pair[1L]; m2 <- pair[2L]
    anti_corr <- -cor_mat[paste0("ME", m1), paste0("ME", m2)]
    if (!is.finite(anti_corr) || anti_corr < 0.45) return(NULL)
    same_trait <- !is.null(top_trait) && identical(top_trait[[m1]], top_trait[[m2]])
    opposite_trait_sign <- !is.null(top_trait_sign) && isTRUE(top_trait_sign[[m1]] * top_trait_sign[[m2]] < 0)
    d1 <- de_dir[match(m1, de_dir$module), , drop = FALSE]
    d2 <- de_dir[match(m2, de_dir$module), , drop = FALSE]
    opposite_de <- nrow(d1) && nrow(d2) && is.finite(d1$mean_log2fc) && is.finite(d2$mean_log2fc) && sign(d1$mean_log2fc) * sign(d2$mean_log2fc) < 0
    score <- anti_corr + if (same_trait) 0.2 else 0 + if (opposite_trait_sign) 0.25 else 0 + if (opposite_de) 0.35 else 0
    rationale <- paste(
      sprintf("Eigengenes are strongly inverse (r=%.2f).", cor_mat[paste0('ME', m1), paste0('ME', m2)]),
      if (same_trait) sprintf("They point to the same dominant trait (%s).", top_trait[[m1]]) else "They do not share the same dominant trait.",
      if (opposite_trait_sign) "That trait is associated in opposite directions across the two modules." else NULL,
      if (opposite_de) sprintf("Their mean DE directions are opposite (%s vs %s).", d1$direction, d2$direction) else "Current DE results do not clearly support opposite module directions."
    )
    data.frame(
      modules = paste(pair, collapse = " + "),
      module_a = m1,
      module_b = m2,
      eigengene_correlation = cor_mat[paste0("ME", m1), paste0("ME", m2)],
      dominant_trait = if (same_trait) top_trait[[m1]] else "different or unavailable",
      opposite_trait_sign = opposite_trait_sign,
      opposite_de_direction = opposite_de,
      suggestion_score = score,
      rationale = rationale,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out) || !nrow(out)) return(data.frame())
  out[order(out$suggestion_score, decreasing = TRUE), , drop = FALSE]
}

build_module_program <- function(wgcna_res, gene_map, modules, name = NULL, de_df = NULL, trait_cols = NULL, alpha = 0.05, lfc_cutoff = 0, include_enrichment = TRUE, include_traits = TRUE) {
  modules <- unique(modules[modules %in% unique(wgcna_res$gene_modules$module)])
  modules <- modules[modules != "grey"]
  if (length(modules) < 2L) return(NULL)
  me_cols <- paste0("ME", modules)
  eigs <- wgcna_res$module_eigengenes[, me_cols, drop = FALSE]
  anchor <- me_cols[1L]
  signs <- rep(1, length(me_cols)); names(signs) <- me_cols
  for (nm in me_cols[-1L]) {
    cc <- suppressWarnings(stats::cor(wgcna_res$module_eigengenes[, anchor], wgcna_res$module_eigengenes[, nm], use = "pairwise.complete.obs"))
    if (is.finite(cc) && cc < 0) signs[nm] <- -1
  }
  aligned_me <- sweep(eigs, 2L, signs[colnames(eigs)], "*")
  program_score <- rowMeans(aligned_me, na.rm = TRUE)
  genes <- unique(wgcna_res$gene_modules$gene_id[wgcna_res$gene_modules$module %in% modules])
  gene_tbl <- wgcna_res$gene_modules[wgcna_res$gene_modules$gene_id %in% genes, , drop = FALSE]
  gene_tbl <- attach_gene_labels(gene_tbl, gene_map)
  gene_tbl$program_name <- name %||% paste(modules, collapse = " + ")
  gene_tbl$module_orientation <- ifelse(signs[paste0("ME", gene_tbl$module)] < 0, "inverse arm", "aligned arm")
  if (!is.null(de_df) && nrow(de_df) > 0L) {
    de_match <- de_df[match(gene_tbl$gene_id, de_df$gene_id), , drop = FALSE]
    gene_tbl$log2FoldChange <- de_match$log2FoldChange
    gene_tbl$padj <- de_match$padj
    gene_tbl$de_supported <- deg_mask(de_match, alpha = alpha, lfc_cutoff = lfc_cutoff)
  } else {
    gene_tbl$log2FoldChange <- NA_real_
    gene_tbl$padj <- NA_real_
    gene_tbl$de_supported <- FALSE
  }
  trait_cor <- NULL
  trait_p <- NULL
  if (isTRUE(include_traits) && !is.null(wgcna_res$module_trait_cor) && !is.null(wgcna_res$metadata)) {
    traits <- prepare_wgcna_traits(wgcna_res$metadata, trait_cols %||% setdiff(names(wgcna_res$metadata), "run_id"))
    if (!is.null(traits) && ncol(traits) > 0L) {
      trait_cor <- stats::cor(program_score, traits, use = "pairwise.complete.obs")
      trait_p <- WGCNA::corPvalueStudent(matrix(trait_cor, nrow = 1L), nrow(wgcna_res$dat_expr))
      colnames(trait_p) <- colnames(traits)
    }
  }
  list(
    name = name %||% paste(modules, collapse = " + "),
    modules = modules,
    signs = signs,
    score = program_score,
    aligned_eigengenes = aligned_me,
    gene_table = gene_tbl,
    enrichment = if (isTRUE(include_enrichment)) run_module_enrichment(genes, gene_map, top_n = 15) else list(combined = data.frame(), summary = data.frame(), has_signal = FALSE),
    trait_cor = trait_cor,
    trait_p = trait_p
  )
}

prioritize_semgraph_rows <- function(unit_rows, de_df = NULL, max_genes = 500L, alpha = 0.05, lfc_cutoff = 0) {
  rows <- unit_rows
  if (nrow(rows) <= max_genes) return(rows)
  rows$priority_score <- rows$kWithin %||% 0
  rows$abs_log2fc <- 0
  rows$de_supported <- FALSE
  if (!is.null(de_df) && nrow(de_df) > 0L && "gene_id" %in% names(rows)) {
    de_match <- de_df[match(rows$gene_id, de_df$gene_id), , drop = FALSE]
    rows$abs_log2fc <- ifelse(is.na(de_match$log2FoldChange), 0, abs(de_match$log2FoldChange))
    rows$de_supported <- deg_mask(de_match, alpha = alpha, lfc_cutoff = lfc_cutoff)
  }
  kw <- rows$kWithin
  kw[!is.finite(kw)] <- 0
  if (diff(range(kw, na.rm = TRUE)) > 0) {
    kw_scaled <- (kw - min(kw, na.rm = TRUE)) / diff(range(kw, na.rm = TRUE))
  } else {
    kw_scaled <- rep(0, length(kw))
  }
  rows$priority_score <- 5 * rows$de_supported + rows$abs_log2fc + 2 * kw_scaled
  rows <- rows[order(rows$priority_score, decreasing = TRUE), , drop = FALSE]
  rows <- rows[!duplicated(rows$hgnc_symbol) & nzchar(rows$hgnc_symbol), , drop = FALSE]
  head(rows, max_genes)
}

pick_binary_group <- function(metadata, col_name) {
  if (is.null(col_name) || !nzchar(col_name) || !col_name %in% names(metadata)) {
    return(NULL)
  }
  vals <- metadata[[col_name]]
  if (is.factor(vals) || is.character(vals)) {
    vals <- factor(vals)
    if (nlevels(vals) != 2L) {
      return(NULL)
    }
    return(as.integer(vals == levels(vals)[2L]))
  }
  uniq <- sort(unique(stats::na.omit(vals)))
  if (length(uniq) == 2L) {
    return(as.integer(vals == uniq[2L]))
  }
  NULL
}

resolve_semgraph_unit <- function(unit_id, wgcna_res, gene_map, program_tbl = NULL, de_df = NULL, trait_cols = NULL, alpha = 0.05, lfc_cutoff = 0) {
  if (!nzchar(unit_id)) return(NULL)
  if (startsWith(unit_id, "module:")) {
    module_color <- sub("^module:", "", unit_id)
    gm <- wgcna_res$gene_modules
    rows <- gm[gm$module == module_color, , drop = FALSE]
    rows <- attach_gene_labels(rows, gene_map)
    rows <- rows[!duplicated(rows$hgnc_symbol) & nzchar(rows$hgnc_symbol), , drop = FALSE]
    return(list(
      rows = rows,
      label = module_color,
      type = "module",
      source_modules = module_color
    ))
  }
  if (startsWith(unit_id, "program:")) {
    program_name <- sub("^program:", "", unit_id)
    row <- program_tbl[program_tbl$name == program_name, , drop = FALSE]
    if (!nrow(row)) return(NULL)
    modules <- trimws(unlist(strsplit(row$modules[1L], "\\+", perl = TRUE)))
    modules <- modules[nzchar(modules)]
    prog <- build_module_program(
      wgcna_res,
      gene_map,
      modules = modules,
      name = row$name[1L],
      de_df = de_df,
      trait_cols = trait_cols,
      alpha = alpha,
      lfc_cutoff = lfc_cutoff,
      include_enrichment = FALSE,
      include_traits = FALSE
    )
    if (is.null(prog)) return(NULL)
    rows <- prog$gene_table
    rows <- rows[!duplicated(rows$hgnc_symbol) & nzchar(rows$hgnc_symbol), , drop = FALSE]
    return(list(
      rows = rows,
      label = row$name[1L],
      type = "program",
      source_modules = paste(modules, collapse = " + ")
    ))
  }
  NULL
}

build_semgraph_input <- function(unit_id, wgcna_res, gene_map, program_tbl = NULL, de_df = NULL, trait_cols = NULL, max_genes = 500L, alpha = 0.05, lfc_cutoff = 0) {
  unit <- resolve_semgraph_unit(unit_id, wgcna_res, gene_map, program_tbl = program_tbl, de_df = de_df, trait_cols = trait_cols, alpha = alpha, lfc_cutoff = lfc_cutoff)
  if (is.null(unit)) {
    return(NULL)
  }
  unit_rows <- prioritize_semgraph_rows(unit$rows, de_df = de_df, max_genes = max_genes, alpha = alpha, lfc_cutoff = lfc_cutoff)
  if (nrow(unit_rows) < 4L) {
    return(NULL)
  }
  expr <- wgcna_res$dat_expr[, unit_rows$gene_id, drop = FALSE]
  colnames(expr) <- unit_rows$hgnc_symbol
  corr <- stats::cor(expr, use = "pairwise.complete.obs")
  diag(corr) <- 0
  thr <- stats::quantile(abs(corr[upper.tri(corr)]), probs = 0.75, na.rm = TRUE)
  adj <- abs(corr) >= thr
  diag(adj) <- 0
  graph <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  graph <- igraph::delete_vertices(graph, igraph::V(graph)[igraph::degree(graph) == 0])
  if (igraph::vcount(graph) < 4L) {
    return(NULL)
  }

  expr <- expr[, igraph::V(graph)$name, drop = FALSE]
  transformed <- SEMgraph::transformData(expr)$data
  list(graph = graph, data = transformed, unit = unit, input_gene_count = nrow(unit_rows), original_gene_count = nrow(unit$rows))
}

run_semgraph_analysis <- function(unit_id, wgcna_res, gene_map, metadata, program_tbl = NULL, de_df = NULL, trait_cols = NULL, refresh_external = FALSE, refresh_fit = FALSE, max_model_genes = 500L, annotation_max_genes = 200L, reactome_live = FALSE, use_live_string = FALSE, group_col = NULL, alpha = 0.05, lfc_cutoff = 0, progress = NULL) {
  string_file <- locate_string_file()
  string_stamp <- if (nzchar(string_file) && file.exists(string_file)) paste(basename(string_file), as.numeric(file.info(string_file)$mtime)) else "no_string_file"
  fit_cache_key <- c(
    "semgraph_fit_v3",
    unit_id,
    max_model_genes,
    group_col %||% "",
    paste(sort(trait_cols %||% character()), collapse = "|"),
    ifelse(use_live_string, "live_string", "cached_string"),
    string_stamp
  )
  if (!refresh_fit) {
    cached <- read_cached_rds("semgraph_fit", fit_cache_key, refresh = FALSE)
    if (!is.null(cached)) {
      if (!is.null(progress)) progress(0.85, "Refreshing external edge validation on cached SEMgraph fit")
      cached$edge_validation <- build_edge_validation_table(cached$graph, refresh = refresh_external, use_live_string = use_live_string)
      if (!is.null(progress)) progress(0.95, "Loaded cached SEMgraph fit")
      return(cached)
    }
  }
  if (!is.null(progress)) progress(0.05, "Selecting genes for SEMgraph fit")
  base_input <- build_semgraph_input(unit_id, wgcna_res, gene_map, program_tbl = program_tbl, de_df = de_df, trait_cols = trait_cols, max_genes = max_model_genes, alpha = alpha, lfc_cutoff = lfc_cutoff)
  if (is.null(base_input)) {
    return(NULL)
  }

  group_vec <- pick_binary_group(metadata[match(rownames(base_input$data), metadata$run_id), , drop = FALSE], group_col)
  if (!is.null(progress)) progress(0.1, "Converting prior graph to a directed graph")
  directed_graph <- igraph::as_directed(base_input$graph, mode = "arbitrary")
  if (!is.null(progress)) progress(0.12, "Converting directed graph to a DAG")
  prior_dag <- SEMgraph::graph2dag(graph = directed_graph, data = base_input$data)
  if (!is.null(progress)) progress(0.18, sprintf("Fitting SEMgraph DAG on %s prioritized genes", ncol(base_input$data)))
  dag_fit <- withCallingHandlers(
    SEMgraph::SEMdag(graph = prior_dag, data = base_input$data, LO = "TO", beta = 0.1),
    message = function(m) {
      msg <- conditionMessage(m)
      if (!is.null(progress)) {
        if (grepl("Node Linear Ordering", msg, fixed = TRUE)) {
          progress(0.32, sprintf("SEMgraph node ordering in progress on %s genes. This is the slowest step and may take a few minutes.", ncol(base_input$data)))
        } else if (grepl("Conducting the nonparanormal transformation", msg, fixed = TRUE)) {
          progress(0.22, "Applying SEMgraph data transformation")
        }
      }
    }
  )
  if (!is.null(progress)) progress(0.58, "Extracting directed regulator structure")
  dag_graph <- dag_fit$dag

  node_names <- igraph::V(dag_graph)$name
  expr_cols <- match(node_names, colnames(base_input$data))
  cor_to_all <- stats::cor(base_input$data[, expr_cols, drop = FALSE], use = "pairwise.complete.obs")
  diag(cor_to_all) <- 0
  out_deg <- igraph::degree(dag_graph, mode = "out")
  in_deg <- igraph::degree(dag_graph, mode = "in")
  total_deg <- out_deg + in_deg
  focal_gene <- names(sort(total_deg, decreasing = TRUE))[1L]
  upstream_nodes <- SEMgraph::ancestors(dag_graph, focal_gene)
  downstream_nodes <- SEMgraph::descendants(dag_graph, focal_gene)
  parent_nodes <- SEMgraph::parents(dag_graph, focal_gene)
  neighbor_nodes <- igraph::neighbors(dag_graph, focal_gene, mode = "all")
  neighbor_nodes <- igraph::V(dag_graph)$name[as.integer(neighbor_nodes)]
  local_targets <- unique(c(
    focal_gene,
    head(parent_nodes, 8L),
    head(upstream_nodes, 12L),
    head(downstream_nodes, 12L),
    head(neighbor_nodes, 12L)
  ))
  local_targets <- local_targets[local_targets %in% node_names]
  if (length(local_targets) < 3L) {
    ranked_nodes <- names(sort(total_deg, decreasing = TRUE))
    local_targets <- unique(c(focal_gene, head(ranked_nodes, 20L)))
    local_targets <- local_targets[local_targets %in% node_names]
  }
  if (!is.null(progress)) progress(0.68, "Preparing ranked genes for annotation")
  gene_lookup <- gene_map[, c("gene_id", "hgnc_symbol", "full_gene_name", "gene_label"), drop = FALSE]
  gene_lookup <- gene_lookup[!is.na(gene_lookup$hgnc_symbol) & nzchar(gene_lookup$hgnc_symbol), , drop = FALSE]
  gene_lookup <- gene_lookup[!duplicated(gene_lookup$hgnc_symbol), , drop = FALSE]
  all_gene_ids <- gene_lookup$gene_id[match(node_names, gene_lookup$hgnc_symbol)]
  local_gene_modules <- wgcna_res$gene_modules[match(all_gene_ids, wgcna_res$gene_modules$gene_id), , drop = FALSE]
  local_gene_modules$graph_gene <- node_names
  local_gene_modules <- local_gene_modules[!is.na(local_gene_modules$gene_id) & nzchar(local_gene_modules$gene_id), , drop = FALSE]
  hub_cut <- stats::quantile(local_gene_modules$kWithin, probs = 0.9, na.rm = TRUE)
  local_gene_modules$is_hub <- local_gene_modules$kWithin >= hub_cut
  tf_ann <- annotate_tf_status(local_gene_modules$gene_id, gene_map)
  local_annotations <- merge(local_gene_modules, tf_ann[, c("gene_id", "is_tf", "tf_source", "gene_label", "full_gene_name", "hgnc_symbol")], by = "gene_id", all.x = TRUE, sort = FALSE)
  local_annotations <- local_annotations[!is.na(local_annotations$gene_id) & nzchar(local_annotations$gene_id), , drop = FALSE]
  local_annotations <- local_annotations[!duplicated(local_annotations$gene_id), , drop = FALSE]
  annotate_rows <- prioritize_annotation_rows(local_annotations, de_df = de_df, max_genes = annotation_max_genes, alpha = alpha, lfc_cutoff = lfc_cutoff)
  if (!is.null(progress)) progress(0.74, sprintf("Loading pharmacology and pathway context for %s ranked genes", nrow(annotate_rows)))
  drug_ann <- aggregate_drug_annotations(annotate_rows$hgnc_symbol, gene_map = gene_map, refresh = refresh_external, progress = progress, reactome_live = reactome_live)
  local_annotations <- merge(local_annotations, drug_ann, by = "hgnc_symbol", all.x = TRUE, sort = FALSE)
  local_annotations$is_druggable <- ifelse(is.na(local_annotations$is_druggable), FALSE, local_annotations$is_druggable)
  local_annotations$drug_mode <- ifelse(is.na(local_annotations$drug_mode), "no interaction found", local_annotations$drug_mode)
  local_annotations$drug_summary <- ifelse(is.na(local_annotations$drug_summary), "No external drug-target interaction found", local_annotations$drug_summary)
  local_annotations$drug_count <- ifelse(is.na(local_annotations$drug_count), 0L, local_annotations$drug_count)
  local_annotations$drug_source <- ifelse(is.na(local_annotations$drug_source), "None", local_annotations$drug_source)
  local_annotations$drug_source_details <- ifelse(is.na(local_annotations$drug_source_details), "No DGIdb, HCDT, or Reactome drug context found", local_annotations$drug_source_details)
  local_annotations$reactome_context <- ifelse(is.na(local_annotations$reactome_context), "No Reactome drug-mediated pathway context found", local_annotations$reactome_context)
  if (nrow(local_annotations) > 0L) {
    rownames(local_annotations) <- local_annotations$hgnc_symbol
  }
  if (!is.null(progress)) progress(0.80, "Validating SEMgraph edges with STRING")
  edge_validation <- build_edge_validation_table(dag_graph, refresh = refresh_external, use_live_string = use_live_string)
  ace_df <- data.frame()
  ace_note <- "ACE summary not computed."
  if (!is.null(group_vec) && length(local_targets) >= 2L && length(local_targets) <= 25L) {
    if (!is.null(progress)) progress(0.86, "Estimating local causal effects")
    ace_df <- tryCatch({
      subg <- igraph::induced_subgraph(dag_graph, vids = local_targets)
      subdata <- base_input$data[, local_targets, drop = FALSE]
      out <- utils::capture.output(
        ace <- SEMgraph::SEMace(graph = subg, data = subdata, group = group_vec),
        type = "message"
      )
      ace_tab <- as.data.frame(ace, stringsAsFactors = FALSE)
      ace_tab
    }, error = function(...) data.frame())
    ace_note <- if (nrow(ace_df) > 0L) {
      "ACE summary computed on the focal regulator neighborhood only, to avoid unstable all-by-all effect estimation."
    } else {
      "ACE summary could not be estimated robustly for this module and group on the current data."
    }
  } else if (is.null(group_vec)) {
    ace_note <- "ACE summary requires a binary group selection."
  } else {
    ace_note <- "ACE summary was skipped because the local causal neighborhood was too small or too large for a stable quick estimate."
  }

  out <- list(
    graph = dag_graph,
    local_nodes = local_targets,
    data = base_input$data,
    focal_gene = focal_gene,
    degree_df = data.frame(
      gene = node_names,
      in_degree = in_deg[node_names],
      out_degree = out_deg[node_names],
      total_degree = total_deg[node_names],
      stringsAsFactors = FALSE
    ),
    local_annotations = local_annotations,
    edge_validation = edge_validation,
    upstream = upstream_nodes,
    downstream = downstream_nodes,
    parents = parent_nodes,
    ace = ace_df,
    ace_note = ace_note,
    group_col = group_col %||% "",
    module_color = base_input$unit$label,
    analysis_unit = base_input$unit$label,
    analysis_unit_type = base_input$unit$type,
    source_modules = base_input$unit$source_modules,
    semgraph_input_gene_count = base_input$input_gene_count,
    semgraph_original_gene_count = base_input$original_gene_count
  )
  if (!is.null(progress)) progress(0.97, "Saving SEMgraph fit to cache")
  write_cached_rds("semgraph_fit", fit_cache_key, out)
  out
}

prioritize_annotation_rows <- function(local_gene_modules, de_df = NULL, max_genes = 200L, alpha = 0.05, lfc_cutoff = 0) {
  rows <- local_gene_modules
  if (nrow(rows) <= max_genes) return(rows)
  rows$priority_score <- rows$kWithin
  if (!is.null(de_df) && nrow(de_df) > 0L) {
    de_match <- de_df[match(rows$gene_id, de_df$gene_id), , drop = FALSE]
    rows$de_supported <- deg_mask(de_match, alpha = alpha, lfc_cutoff = lfc_cutoff)
    rows$abs_log2fc <- ifelse(is.na(de_match$log2FoldChange), 0, abs(de_match$log2FoldChange))
    rows$priority_score <- 5 * rows$de_supported + rows$abs_log2fc + rows$kWithin / max(1, max(rows$kWithin, na.rm = TRUE))
  }
  rows <- rows[order(rows$priority_score, decreasing = TRUE), , drop = FALSE]
  head(rows, max_genes)
}

rank_semgraph_nodes <- function(res, de_df = NULL, alpha = 0.05, lfc_cutoff = 0) {
  graph_nodes <- igraph::V(res$graph)$name
  ann <- res$local_annotations[match(graph_nodes, rownames(res$local_annotations)), , drop = FALSE]
  deg <- res$degree_df[match(graph_nodes, res$degree_df$gene), , drop = FALSE]
  de_supported <- rep(FALSE, length(graph_nodes))
  if (!is.null(de_df) && nrow(de_df) > 0L) {
    de_match <- de_df[match(ann$gene_id, de_df$gene_id), , drop = FALSE]
    de_supported <- deg_mask(de_match, alpha = alpha, lfc_cutoff = lfc_cutoff)
  }
  graph_ud <- igraph::as_undirected(res$graph, mode = "collapse")
  dist_vec <- tryCatch(as.numeric(igraph::distances(graph_ud, v = res$focal_gene, to = graph_nodes)), error = function(...) rep(Inf, length(graph_nodes)))
  data.frame(
    gene = graph_nodes,
    total_degree = deg$total_degree,
    in_degree = deg$in_degree,
    out_degree = deg$out_degree,
    is_hub = ifelse(is.na(ann$is_hub), FALSE, ann$is_hub),
    is_tf = ifelse(is.na(ann$is_tf), FALSE, ann$is_tf),
    is_druggable = ifelse(is.na(ann$is_druggable), FALSE, ann$is_druggable),
    de_supported = de_supported,
    distance_to_focal = dist_vec,
    stringsAsFactors = FALSE
  )
}

select_semgraph_display_nodes <- function(res, de_df, mode = "de_connectors", min_nodes = 6L, max_context = 10L, max_nodes = 60L, alpha = 0.05, lfc_cutoff = 0) {
  graph_nodes <- igraph::V(res$graph)$name
  ann <- res$local_annotations
  if (is.null(ann) || nrow(ann) == 0L) {
    return(unique(res$local_nodes[res$local_nodes %in% graph_nodes]))
  }
  if (length(max_nodes) == 0L || is.na(max_nodes) || !is.finite(max_nodes)) {
    max_nodes <- 60L
  }
  max_nodes <- max(1L, as.integer(max_nodes[[1L]]))
  if (length(min_nodes) == 0L || is.na(min_nodes) || !is.finite(min_nodes)) {
    min_nodes <- 6L
  }
  min_nodes <- max(1L, as.integer(min_nodes[[1L]]))
  if (length(max_context) == 0L || is.na(max_context) || !is.finite(max_context)) {
    max_context <- 10L
  }
  max_context <- max(0L, as.integer(max_context[[1L]]))
  target_nodes <- max(min_nodes, min(max_nodes, max(20L, floor(max_nodes * 0.6))))
  max_context <- max(max_context, target_nodes)

  supported_symbols <- character()
  if (!is.null(de_df) && nrow(de_df) > 0L) {
    de_match <- de_df[match(ann$gene_id, de_df$gene_id), , drop = FALSE]
    supported_idx <- deg_mask(de_match, alpha = alpha, lfc_cutoff = lfc_cutoff)
    supported_symbols <- unique(ann$hgnc_symbol[supported_idx & !is.na(ann$hgnc_symbol) & nzchar(ann$hgnc_symbol)])
  }

  if (identical(mode, "full_local")) {
    ranked <- rank_semgraph_nodes(res, de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
    ranked$score <- 5 * ranked$de_supported + 2 * ranked$is_hub + 2 * ranked$is_tf + 2 * ranked$is_druggable + ranked$total_degree / max(1, max(ranked$total_degree, na.rm = TRUE)) - pmin(ranked$distance_to_focal, 6)
    ranked <- ranked[order(ranked$score, decreasing = TRUE), , drop = FALSE]
    keep <- unique(c(res$focal_gene, head(ranked$gene, max_nodes)))
    return(keep[keep %in% graph_nodes])
  }

  display_nodes <- unique(c(res$focal_gene, supported_symbols))
  display_nodes <- display_nodes[display_nodes %in% graph_nodes]

  if (identical(mode, "de_connectors") && length(display_nodes) >= 2L) {
    connector_nodes <- character()
    graph_ud <- igraph::as_undirected(res$graph, mode = "collapse")
    for (target in setdiff(display_nodes, res$focal_gene)) {
      path <- tryCatch(
        igraph::shortest_paths(graph_ud, from = res$focal_gene, to = target, output = "vpath")$vpath[[1]],
        error = function(...) NULL
      )
      if (!is.null(path) && length(path) > 0L) {
        connector_nodes <- c(connector_nodes, igraph::V(graph_ud)$name[as.integer(path)])
      }
    }
    display_nodes <- unique(c(display_nodes, connector_nodes))
  }

  if (!identical(mode, "de_only") && length(display_nodes) < target_nodes) {
    ranked_tbl <- rank_semgraph_nodes(res, de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
    ranked_tbl$score <- 3 * ranked_tbl$is_hub + 2 * ranked_tbl$is_tf + 2 * ranked_tbl$is_druggable + ranked_tbl$total_degree / max(1, max(ranked_tbl$total_degree, na.rm = TRUE)) - pmin(ranked_tbl$distance_to_focal, 6)
    ranked_context <- ranked_tbl$gene[order(ranked_tbl$score, decreasing = TRUE)]
    ranked_context <- setdiff(ranked_context, display_nodes)
    display_nodes <- unique(c(display_nodes, head(ranked_context, max_context)))
  }
  ranked_tbl <- rank_semgraph_nodes(res, de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
  ranked_tbl$score <- 6 * ranked_tbl$de_supported + 3 * ranked_tbl$is_hub + 2 * ranked_tbl$is_tf + 2 * ranked_tbl$is_druggable + ranked_tbl$total_degree / max(1, max(ranked_tbl$total_degree, na.rm = TRUE)) - pmin(ranked_tbl$distance_to_focal, 8) / 2
  ranked_genes <- ranked_tbl$gene[order(ranked_tbl$score, decreasing = TRUE)]
  prioritized <- unique(c(display_nodes, ranked_genes))
  prioritized <- prioritized[prioritized %in% graph_nodes]
  head(prioritized, max_nodes)
}

select_semgraph_table_nodes <- function(res, de_df, scope = "evidence_context", max_rows = 200L, alpha = 0.05, lfc_cutoff = 0) {
  ranked <- rank_semgraph_nodes(res, de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
  ranked$score <- 6 * ranked$de_supported + 3 * ranked$is_hub + 2 * ranked$is_tf + 2 * ranked$is_druggable + ranked$total_degree / max(1, max(ranked$total_degree, na.rm = TRUE)) - pmin(ranked$distance_to_focal, 10) / 3
  if (identical(scope, "de_only")) {
    keep <- ranked$gene[ranked$de_supported | ranked$gene == res$focal_gene]
  } else if (identical(scope, "de_plus_regulators")) {
    keep <- ranked$gene[ranked$de_supported | ranked$is_hub | ranked$is_tf | ranked$is_druggable | ranked$gene == res$focal_gene]
  } else {
    keep <- ranked$gene[ranked$de_supported | ranked$is_hub | ranked$is_tf | ranked$is_druggable | ranked$distance_to_focal <= 5 | ranked$gene == res$focal_gene]
  }
  keep_tbl <- ranked[ranked$gene %in% keep, , drop = FALSE]
  if (nrow(keep_tbl) < min(max_rows, 40L)) {
    ranked_extra <- ranked[!ranked$gene %in% keep_tbl$gene, , drop = FALSE]
    keep_tbl <- rbind(keep_tbl, head(ranked_extra, min(max_rows, 40L) - nrow(keep_tbl)))
  }
  keep_tbl <- keep_tbl[order(keep_tbl$score, decreasing = TRUE), , drop = FALSE]
  head(keep_tbl$gene, max_rows)
}

compute_semgraph_layout <- function(graph, method = "fr") {
  if (igraph::vcount(graph) <= 1L) {
    return(matrix(c(0, 0), ncol = 2))
  }
  layout_key <- paste(method, paste(sort(igraph::V(graph)$name), collapse = "||"), sep = "::")
  seed_val <- sum(utf8ToInt(layout_key)) %% .Machine$integer.max
  if (!is.finite(seed_val) || seed_val <= 0) seed_val <- 1L
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    had_seed <- TRUE
  } else {
    old_seed <- NULL
    had_seed <- FALSE
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed_val)
  switch(
    method,
    kk = igraph::layout_with_kk(graph),
    graphopt = igraph::layout_with_graphopt(graph),
    lgl = igraph::layout_with_lgl(graph),
    circle = igraph::layout_in_circle(graph),
    igraph::layout_with_fr(graph)
  )
}

find_semgraph_paths <- function(graph, from, to, top_n = 5L, max_steps = 6L) {
  if (!from %in% igraph::V(graph)$name || !to %in% igraph::V(graph)$name) {
    return(list(paths = list(), nodes = character()))
  }
  raw_paths <- tryCatch(
    igraph::all_simple_paths(graph, from = from, to = to, mode = "out", cutoff = max_steps),
    error = function(...) list()
  )
  if (!length(raw_paths)) {
    return(list(paths = list(), nodes = unique(c(from, to))))
  }
  path_nodes <- lapply(raw_paths, function(vs) igraph::V(graph)$name[as.integer(vs)])
  ord <- order(vapply(path_nodes, length, integer(1)))
  path_nodes <- path_nodes[head(ord, top_n)]
  list(paths = path_nodes, nodes = unique(unlist(path_nodes, use.names = FALSE)))
}

select_connector_nodes <- function(graph, selected_nodes, max_hops = 3L) {
  selected_nodes <- unique(selected_nodes[selected_nodes %in% igraph::V(graph)$name])
  if (!length(selected_nodes)) return(character())
  if (max_hops <= 0L) {
    return(selected_nodes)
  }
  graph_ud <- igraph::as_undirected(graph, mode = "collapse")
  if (length(selected_nodes) == 1L) {
    ego_nodes <- tryCatch(
      igraph::ego(graph_ud, order = max_hops, nodes = selected_nodes[1L], mode = "all")[[1L]],
      error = function(...) integer()
    )
    return(unique(igraph::V(graph_ud)$name[as.integer(ego_nodes)]))
  }
  keep <- selected_nodes
  pairs <- utils::combn(selected_nodes, 2L, simplify = FALSE)
  for (pair in pairs) {
    sp <- tryCatch(
      igraph::shortest_paths(graph_ud, from = pair[1L], to = pair[2L], output = "vpath")$vpath[[1L]],
      error = function(...) NULL
    )
    if (!is.null(sp) && length(sp) > 0L && (length(sp) - 1L) <= max_hops) {
      keep <- c(keep, igraph::V(graph_ud)$name[as.integer(sp)])
    }
  }
  unique(keep)
}

collect_semgraph_focus_nodes <- function(res, ranked, ann, selected_nodes = character(),
                                         include_tf = FALSE, include_druggable = FALSE,
                                         hub_mode = "off", focal_top_n = 0L) {
  focus <- unique(selected_nodes[selected_nodes %in% ranked$gene])
  if (isTRUE(include_tf)) {
    focus <- unique(c(focus, ann$hgnc_symbol[ifelse(is.na(ann$is_tf), FALSE, ann$is_tf)]))
  }
  if (isTRUE(include_druggable)) {
    focus <- unique(c(focus, ann$hgnc_symbol[ifelse(is.na(ann$is_druggable), FALSE, ann$is_druggable)]))
  }
  if (identical(hub_mode, "wgcna_hubs")) {
    focus <- unique(c(focus, ann$hgnc_symbol[ifelse(is.na(ann$is_hub), FALSE, ann$is_hub)]))
  }
  if (isTRUE(focal_top_n > 0L)) {
    regulator_score <- ranked$out_degree - ranked$in_degree + ranked$total_degree / max(1, max(ranked$total_degree, na.rm = TRUE))
    ord <- order(regulator_score, decreasing = TRUE, na.last = NA)
    focus <- unique(c(focus, head(ranked$gene[ord], focal_top_n)))
  }
  unique(focus[nzchar(focus)])
}

select_semgraph_graph_nodes <- function(res, de_df, display_mode = "de_connectors", graph_scope = "overview",
                                        max_nodes = 60L, start_gene = NULL, end_gene = NULL,
                                        hops = 2L, path_n = 5L, path_max_steps = 6L,
                                        selected_nodes = character(), selected_hops = 3L,
                                        selected_only = FALSE, filter_tf = FALSE,
                                        filter_druggable = FALSE, filter_hub_mode = "off",
                                        focal_top_n = 0L, restrict_to_focus = FALSE,
                                        alpha = 0.05, lfc_cutoff = 0) {
  if (length(max_nodes) == 0L || is.na(max_nodes) || !is.finite(max_nodes)) {
    max_nodes <- 60L
  }
  max_nodes <- max(1L, as.integer(max_nodes[[1L]]))
  if (length(hops) == 0L || is.na(hops) || !is.finite(hops)) {
    hops <- 2L
  }
  hops <- max(0L, as.integer(hops[[1L]]))
  if (length(path_n) == 0L || is.na(path_n) || !is.finite(path_n)) {
    path_n <- 5L
  }
  path_n <- max(1L, as.integer(path_n[[1L]]))
  if (length(path_max_steps) == 0L || is.na(path_max_steps) || !is.finite(path_max_steps)) {
    path_max_steps <- 6L
  }
  path_max_steps <- max(1L, as.integer(path_max_steps[[1L]]))
  if (length(selected_hops) == 0L || is.na(selected_hops) || !is.finite(selected_hops)) {
    selected_hops <- 3L
  }
  selected_hops <- max(0L, as.integer(selected_hops[[1L]]))
  if (length(focal_top_n) == 0L || is.na(focal_top_n) || !is.finite(focal_top_n)) {
    focal_top_n <- 0L
  }
  focal_top_n <- max(0L, as.integer(focal_top_n[[1L]]))
  graph_nodes <- igraph::V(res$graph)$name
  base_nodes <- select_semgraph_display_nodes(res, de_df, mode = display_mode, max_nodes = max(max_nodes, 20L), alpha = alpha, lfc_cutoff = lfc_cutoff)
  base_nodes <- unique(base_nodes[base_nodes %in% graph_nodes])
  if (!length(base_nodes)) {
    return(list(
      nodes = unique(res$local_nodes[res$local_nodes %in% graph_nodes]),
      paths = list(),
      mode = graph_scope,
      start_gene = res$focal_gene,
      end_gene = end_gene %||% res$focal_gene
    ))
  }

  ranked <- rank_semgraph_nodes(res, de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
  ann <- res$local_annotations[match(ranked$gene, rownames(res$local_annotations)), , drop = FALSE]
  ranked$score <- 6 * ranked$de_supported + 3 * ranked$is_hub + 2 * ranked$is_tf + 2 * ranked$is_druggable +
    ranked$total_degree / max(1, max(ranked$total_degree, na.rm = TRUE)) - pmin(ranked$distance_to_focal, 8) / 2
  focus_anchor_nodes <- collect_semgraph_focus_nodes(
    res,
    ranked,
    ann,
    selected_nodes = selected_nodes,
    include_tf = filter_tf,
    include_druggable = filter_druggable,
    hub_mode = filter_hub_mode,
    focal_top_n = focal_top_n
  )
  candidate_nodes <- unique(c(focus_anchor_nodes, selected_nodes, start_gene, end_gene, base_nodes, head(ranked$gene[order(ranked$score, decreasing = TRUE)], max(max_nodes * 4L, 120L))))
  candidate_nodes <- candidate_nodes[candidate_nodes %in% graph_nodes]
  candidate_graph <- igraph::induced_subgraph(res$graph, vids = candidate_nodes)

  start_gene <- start_gene %||% res$focal_gene
  if (!start_gene %in% igraph::V(candidate_graph)$name) {
    start_gene <- res$focal_gene
  }

  if (identical(graph_scope, "neighborhood")) {
    ego_nodes <- tryCatch(
      igraph::ego(candidate_graph, order = hops, nodes = start_gene, mode = "all")[[1L]],
      error = function(...) integer()
    )
    hop_nodes <- igraph::V(candidate_graph)$name[as.integer(ego_nodes)]
    hop_nodes <- unique(c(start_gene, hop_nodes))
    ranked_subset <- ranked[ranked$gene %in% hop_nodes, , drop = FALSE]
    ranked_subset <- ranked_subset[order(ranked_subset$score, decreasing = TRUE), , drop = FALSE]
    keep <- unique(c(start_gene, head(ranked_subset$gene, max_nodes)))
    return(list(nodes = keep[keep %in% graph_nodes], paths = list(), mode = graph_scope, start_gene = start_gene, end_gene = end_gene))
  }

  if (identical(graph_scope, "paths")) {
    end_gene <- end_gene %||% res$focal_gene
    if (!end_gene %in% igraph::V(candidate_graph)$name) {
      end_gene <- setdiff(base_nodes, start_gene)[1L] %||% start_gene
    }
    path_info <- find_semgraph_paths(candidate_graph, from = start_gene, to = end_gene, top_n = path_n, max_steps = path_max_steps)
    keep <- unique(c(start_gene, end_gene, path_info$nodes))
    if (length(keep) < 2L) {
      keep <- unique(c(start_gene, end_gene, base_nodes))
    }
    ranked_subset <- ranked[ranked$gene %in% keep, , drop = FALSE]
    ranked_subset <- ranked_subset[order(ranked_subset$score, decreasing = TRUE), , drop = FALSE]
    keep <- unique(c(start_gene, end_gene, head(ranked_subset$gene, max_nodes)))
    return(list(nodes = keep[keep %in% graph_nodes], paths = path_info$paths, mode = graph_scope, start_gene = start_gene, end_gene = end_gene))
  }

  out <- list(nodes = head(base_nodes, max_nodes), paths = list(), mode = graph_scope, start_gene = start_gene, end_gene = end_gene)
  if (isTRUE(selected_only) || isTRUE(restrict_to_focus)) {
    connector_seed_nodes <- unique(c(selected_nodes, focus_anchor_nodes))
    focus_nodes <- select_connector_nodes(candidate_graph, selected_nodes = connector_seed_nodes, max_hops = selected_hops)
    if (length(focus_nodes)) {
      ranked_subset <- ranked[ranked$gene %in% focus_nodes, , drop = FALSE]
      ranked_subset <- ranked_subset[order(ranked_subset$score, decreasing = TRUE), , drop = FALSE]
      out$nodes <- unique(c(connector_seed_nodes, head(ranked_subset$gene, max_nodes)))
    }
  }
  out
}

build_wgcna_module_graph <- function(res, module_view, edge_quantile = 0.9) {
  gm <- res$gene_modules
  mod_genes <- gm$gene_id[gm$module == module_view]
  if (length(mod_genes) < 3L) return(NULL)
  mod_genes <- head(gm$gene_id[order(ifelse(gm$module == module_view, gm$kWithin, -Inf), decreasing = TRUE)], 40L)
  mod_genes <- mod_genes[gm$module[match(mod_genes, gm$gene_id)] == module_view]
  adj <- res$adjacency[mod_genes, mod_genes, drop = FALSE]
  diag(adj) <- 0
  cutoff <- stats::quantile(adj[upper.tri(adj)], probs = edge_quantile, na.rm = TRUE)
  edge_idx <- which(adj >= cutoff, arr.ind = TRUE)
  edge_idx <- edge_idx[edge_idx[, 1L] < edge_idx[, 2L], , drop = FALSE]
  edges <- if (nrow(edge_idx)) {
    data.frame(
      from = rownames(adj)[edge_idx[, 1L]],
      to = colnames(adj)[edge_idx[, 2L]],
      weight = adj[edge_idx],
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(from = character(), to = character(), weight = numeric(), stringsAsFactors = FALSE)
  }
  verts <- data.frame(
    name = mod_genes,
    stringsAsFactors = FALSE
  )
  igraph::graph_from_data_frame(edges, directed = FALSE, vertices = verts)
}

split_context_values <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(x)) return(character())
  vals <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE), use.names = FALSE))
  vals <- vals[nzchar(vals)]
  unique(vals)
}

summarize_modes <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x) & x != "no interaction found"])
  if (!length(x)) return("no direct gene-target drug evidence")
  if (length(x) == 1L) return(x)
  paste(x, collapse = ", ")
}

make_downloadable_table <- function(df, page_length = 10L, rownames = FALSE) {
  DT::datatable(
    df,
    rownames = rownames,
    filter = "top",
    extensions = c("Buttons", "SearchBuilder"),
    options = list(
      pageLength = page_length,
      scrollX = TRUE,
      dom = "QBfrtip",
      buttons = c("copy", "csv"),
      searchBuilder = list(depthLimit = 4)
    ),
    class = "compact stripe hover"
  )
}

make_downloadable_plot <- function(p, filename = "rnaseq_workbench_plot") {
  safe_name <- gsub("[^A-Za-z0-9_]+", "_", filename)
  plotly::config(
    p,
    displaylogo = FALSE,
    responsive = TRUE,
    toImageButtonOptions = list(
      format = "png",
      filename = safe_name,
      height = 900,
      width = 1400,
      scale = 2
    )
  )
}

package_versions_df <- function() {
  data.frame(
    package = all_required_packages,
    version = vapply(
      all_required_packages,
      function(pkg) {
        if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
        as.character(utils::packageVersion(pkg))
      },
      character(1)
    ),
    stringsAsFactors = FALSE
  )
}

append_analysis_log <- function(tbl, step, params = list()) {
  row <- data.frame(
    timestamp = as.character(Sys.time()),
    step = step,
    params_json = jsonlite::toJSON(params, auto_unbox = TRUE, null = "null"),
    stringsAsFactors = FALSE
  )
  rbind(tbl, row)
}

write_text_file <- function(path, lines) {
  writeLines(enc2utf8(lines), path, useBytes = TRUE)
}

report_parse_params <- function(log_df, step_name) {
  if (is.null(log_df) || !nrow(log_df) || !step_name %in% log_df$step) return(list())
  idx <- max(which(log_df$step == step_name))
  raw <- log_df$params_json[idx]
  if (is.null(raw) || is.na(raw) || !nzchar(raw)) return(list())
  tryCatch(jsonlite::fromJSON(raw, simplifyVector = FALSE), error = function(...) list())
}

r_literal <- function(x) {
  if (is.null(x)) return("NULL")
  if (length(x) == 0L) return("character()")
  if (is.list(x) && !is.data.frame(x)) {
    inner <- vapply(x, r_literal, character(1))
    return(sprintf("list(%s)", paste(inner, collapse = ", ")))
  }
  if (is.character(x)) {
    vals <- vapply(x, function(v) encodeString(v, quote = "\""), character(1), USE.NAMES = FALSE)
    if (length(vals) == 1L) return(vals)
    return(sprintf("c(%s)", paste(vals, collapse = ", ")))
  }
  if (is.logical(x)) {
    vals <- ifelse(is.na(x), "NA", ifelse(x, "TRUE", "FALSE"))
    if (length(vals) == 1L) return(vals)
    return(sprintf("c(%s)", paste(vals, collapse = ", ")))
  }
  if (is.numeric(x)) {
    vals <- ifelse(is.na(x), "NA", format(x, scientific = FALSE, trim = TRUE))
    if (length(vals) == 1L) return(vals)
    return(sprintf("c(%s)", paste(vals, collapse = ", ")))
  }
  encodeString(as.character(x), quote = "\"")
}

commented_code_block <- function(lines) {
  paste(lines, collapse = "\n")
}

save_plot_widget <- function(plot_obj, file, title = NULL) {
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    write_text_file(file, c(
      "<html><body>",
      sprintf("<h3>%s</h3>", title %||% "Plot unavailable"),
      "<p>The htmlwidgets package was not available, so this interactive plot could not be saved.</p>",
      "</body></html>"
    ))
    return(invisible(file))
  }
  htmlwidgets::saveWidget(plot_obj, file = file, selfcontained = FALSE, title = title %||% "Plot")
  invisible(file)
}

zip_dir_contents <- function(zipfile, dir_path) {
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(dir_path)
  files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
  utils::zip(zipfile = zipfile, files = files)
}

build_semgraph_pathway_context <- function(res, de_df = NULL, mode = "de_connectors", alpha = 0.05, lfc_cutoff = 0) {
  display_nodes <- select_semgraph_display_nodes(res, de_df, mode = mode, alpha = alpha, lfc_cutoff = lfc_cutoff)
  display_nodes <- display_nodes[display_nodes %in% rownames(res$local_annotations)]
  ann <- res$local_annotations[match(display_nodes, rownames(res$local_annotations)), , drop = FALSE]
  if (nrow(ann) == 0L) {
    return(list(summary = data.frame(), edges = data.frame()))
  }
  roles <- ifelse(
    ann$hgnc_symbol == res$focal_gene, "focal",
    ifelse(ann$hgnc_symbol %in% res$upstream, "upstream",
    ifelse(ann$hgnc_symbol %in% res$downstream, "downstream", "other"))
  )
  de_supported <- rep(FALSE, nrow(ann))
  log2fc <- rep(NA_real_, nrow(ann))
  padj <- rep(NA_real_, nrow(ann))
  if (!is.null(de_df) && nrow(de_df) > 0L) {
    de_match <- de_df[match(ann$gene_id, de_df$gene_id), , drop = FALSE]
    de_supported <- deg_mask(de_match, alpha = alpha, lfc_cutoff = lfc_cutoff)
    log2fc <- de_match$log2FoldChange
    padj <- de_match$padj
  }
  pathway_rows <- lapply(seq_len(nrow(ann)), function(i) {
    pathways <- split_context_values(ann$reactome_context[i])
    pathways <- pathways[!grepl("^No Reactome", pathways, ignore.case = TRUE)]
    if (!length(pathways)) return(NULL)
    data.frame(
      pathway = pathways,
      gene = ann$hgnc_symbol[i],
      gene_label = ann$gene_label[i] %||% ann$hgnc_symbol[i],
      full_gene_name = ann$full_gene_name[i] %||% "Full gene name unavailable",
      role = roles[i],
      de_supported = de_supported[i],
      log2FoldChange = log2fc[i],
      padj = padj[i],
      is_druggable = ann$is_druggable[i],
      drug_mode = ann$drug_mode[i],
      drug_summary = ann$drug_summary[i],
      drug_source = ann$drug_source[i],
      stringsAsFactors = FALSE
    )
  })
  pathway_edges <- do.call(rbind, pathway_rows[!vapply(pathway_rows, is.null, logical(1))])
  if (is.null(pathway_edges) || !nrow(pathway_edges)) {
    return(list(summary = data.frame(), edges = data.frame()))
  }
  split_by_path <- split(pathway_edges, pathway_edges$pathway)
  pathway_summary <- do.call(rbind, lapply(names(split_by_path), function(pw) {
    cur <- split_by_path[[pw]]
    direct_gene_rows <- cur[cur$is_druggable & !grepl("Reactome", cur$drug_source, fixed = TRUE), , drop = FALSE]
    data.frame(
      pathway = pw,
      genes_in_pathway_view = nrow(cur),
      de_supported_genes = sum(cur$de_supported, na.rm = TRUE),
      contains_focal_gene = sum(cur$role == "focal", na.rm = TRUE) > 0,
      upstream_genes = sum(cur$role == "upstream", na.rm = TRUE),
      downstream_genes = sum(cur$role == "downstream", na.rm = TRUE),
      context_genes = sum(cur$role == "other", na.rm = TRUE),
      direct_druggable_genes = sum(direct_gene_rows$is_druggable, na.rm = TRUE),
      dominant_de_direction = if (all(is.na(cur$log2FoldChange))) "unknown" else if (mean(cur$log2FoldChange, na.rm = TRUE) > 0) "mostly upregulated" else "mostly downregulated",
      pathway_drug_effect_hint = summarize_modes(direct_gene_rows$drug_mode),
      supporting_genes = paste(unique(cur$gene), collapse = ", "),
      direct_drug_examples = if (nrow(direct_gene_rows)) paste(unique(direct_gene_rows$drug_summary), collapse = " | ") else "No direct member-gene drug examples found",
      drug_sources = paste(unique(cur$drug_source[cur$is_druggable]), collapse = ", "),
      interpretation = paste(
        "This pathway overlaps the displayed SEMgraph neighborhood.",
        if (sum(cur$de_supported, na.rm = TRUE) > 0) "Some member genes are supported by the current DE result." else "Current support is mainly contextual rather than DE-driven.",
        if (nrow(direct_gene_rows) > 0) "Direct member-gene drug evidence exists, so pathway modulation may be indirectly approachable through those genes." else "No direct member-gene drug evidence was found here, so treat this pathway mainly as mechanistic context."
      ),
      stringsAsFactors = FALSE
    )
  }))
  pathway_summary <- pathway_summary[order(pathway_summary$de_supported_genes, pathway_summary$direct_druggable_genes, decreasing = TRUE), , drop = FALSE]
  list(summary = pathway_summary, edges = pathway_edges)
}

build_formula <- function(terms, interaction = FALSE) {
  if (length(terms) == 0L) return("~ 1")
  rhs <- terms
  if (interaction && length(terms) >= 2L) rhs <- c(rhs, paste(terms, collapse = ":"))
  paste("~", paste(rhs, collapse = " + "))
}

group_colors <- function(groups) {
  lev <- unique(as.character(groups))
  pal <- c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b", "#e377c2", "#17becf")
  out <- rep(pal, length.out = length(lev)); names(out) <- lev; out
}

tick_text <- function(labels, groups, colors) {
  vapply(seq_along(labels), function(i) sprintf("<span style='color:%s;'>%s</span>", colors[[as.character(groups[[i]])]], labels[[i]]), character(1))
}

build_metadata_annotation_plot <- function(meta, sample_ids, metadata_cols = NULL, title = "Sample metadata") {
  if (is.null(meta) || !nrow(meta) || !length(sample_ids)) return(NULL)
  meta_sub <- meta[match(sample_ids, meta$run_id), , drop = FALSE]
  metadata_cols <- metadata_cols %||% setdiff(names(meta_sub), "run_id")
  metadata_cols <- metadata_cols[metadata_cols %in% names(meta_sub)]
  metadata_cols <- metadata_cols[metadata_cols != "run_id"]
  if (!length(metadata_cols)) return(NULL)

  ann <- meta_sub[, metadata_cols, drop = FALSE]
  ann[] <- lapply(ann, function(x) {
    x <- as.character(x)
    x[is.na(x) | !nzchar(x)] <- "missing"
    x
  })

  ann_values <- unlist(ann, use.names = FALSE)
  lev <- unique(ann_values)
  pal <- grDevices::hcl.colors(max(3L, length(lev)), "Set 3")
  cols <- stats::setNames(pal[seq_along(lev)], lev)
  value_to_num <- stats::setNames(seq_along(lev), lev)
  z <- t(vapply(ann, function(col) unname(value_to_num[col]), numeric(length(sample_ids))))
  rownames(z) <- metadata_cols
  colnames(z) <- sample_ids

  text_mat <- matrix("", nrow = nrow(z), ncol = ncol(z), dimnames = dimnames(z))
  for (i in seq_len(nrow(z))) {
    for (j in seq_len(ncol(z))) {
      text_mat[i, j] <- paste0("Sample: ", sample_ids[j], "<br>Metadata: ", metadata_cols[i], "<br>Value: ", ann[[i]][j])
    }
  }

  n_cols <- length(lev)
  colorscale <- unlist(lapply(seq_len(n_cols), function(i) {
    low <- if (n_cols == 1L) 0 else (i - 1) / n_cols
    high <- i / n_cols
    list(list(low, cols[[lev[i]]]), list(high, cols[[lev[i]]]))
  }), recursive = FALSE)

  plotly::plot_ly(
    x = sample_ids,
    y = rev(metadata_cols),
    z = z[nrow(z):1, , drop = FALSE],
    type = "heatmap",
    colorscale = colorscale,
    showscale = FALSE,
    text = text_mat[nrow(text_mat):1, , drop = FALSE],
    hoverinfo = "text"
  ) |>
    plotly::layout(
      title = list(text = title, x = 0),
      xaxis = list(title = "", showticklabels = FALSE),
      yaxis = list(title = "", automargin = TRUE)
    )
}

design_note <- function(meta, terms, interaction, samples, counts) {
  if (length(terms) == 0L) return("Add at least one metadata term for a contrastable model.")
  notes <- c()
  libs <- colSums(counts[, samples, drop = FALSE])
  if (max(libs) / max(1, min(libs)) > 3) notes <- c(notes, "Library sizes vary by more than three-fold across selected samples.")
  for (term in terms) {
    if (term %in% names(meta) && is.factor(meta[[term]])) {
      tab <- table(meta[match(samples, meta$run_id), term, drop = TRUE])
      if (any(tab < 2L)) notes <- c(notes, sprintf("Some `%s` groups have fewer than two replicates.", term))
    }
  }
  if (interaction && length(terms) < 2L) notes <- c(notes, "Interactions need at least two terms.")
  paste(c(
    paste(sprintf("`%s`", terms), collapse = " + "),
    if (interaction && length(terms) >= 2L) "with an interaction term" else NULL,
    if (length(notes)) paste(notes, collapse = " ") else "looks reasonable for a first-pass model."
  ), collapse = " ")
}

deg_mask <- function(df, alpha = 0.05, lfc_cutoff = 0) {
  !is.na(df$padj) &
    df$padj < alpha &
    !is.na(df$log2FoldChange) &
    abs(df$log2FoldChange) >= lfc_cutoff
}

run_de <- function(counts, meta, samples, formula_str, gene_map, var = NULL, num = NULL, den = NULL) {
  meta_sub <- meta[match(samples, meta$run_id), , drop = FALSE]
  rownames(meta_sub) <- meta_sub$run_id
  dds <- DESeq2::DESeqDataSetFromMatrix(round(counts[, samples, drop = FALSE]), meta_sub, stats::as.formula(formula_str))
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  if (!is.null(var) && nzchar(var)) {
    res <- DESeq2::results(dds, contrast = c(var, num, den))
    label <- sprintf("%s: %s vs %s", var, num, den)
  } else {
    rn <- DESeq2::resultsNames(dds)
    if (length(rn) < 2L) stop("No contrastable coefficient in this design.")
    res <- DESeq2::results(dds, name = rn[2L]); label <- rn[2L]
  }
  df <- as.data.frame(res); df$gene_id <- rownames(df); df$sig <- deg_mask(df, alpha = 0.05, lfc_cutoff = 0); df <- df[order(df$padj, na.last = TRUE), ]
  df <- attach_gene_labels(df, gene_map)
  vst_obj <- DESeq2::vst(dds, blind = TRUE)
  vst <- SummarizedExperiment::assay(vst_obj)
  vst_pca <- prcomp(t(vst), center = TRUE, scale. = TRUE)
  disp_df <- as.data.frame(SummarizedExperiment::mcols(dds))
  disp_df$gene_id <- rownames(dds)
  disp_df <- attach_gene_labels(disp_df, gene_map)
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)
  vst_dist <- as.matrix(dist(t(vst)))
  size_factors <- DESeq2::sizeFactors(dds)
  top_abs_lfc <- head(df$gene_id[order(abs(df$log2FoldChange), decreasing = TRUE, na.last = TRUE)], 20L)
  mean_sd_df <- data.frame(
    gene_id = rownames(vst),
    vst_mean = rowMeans(vst, na.rm = TRUE),
    vst_sd = apply(vst, 1L, stats::sd, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  mean_sd_df <- attach_gene_labels(mean_sd_df, gene_map)
  top <- head(df$gene_id[!is.na(df$padj)], 30L)
  list(
    df = df,
    label = label,
    sig = sum(df$sig, na.rm = TRUE),
    heat = if (length(top)) vst[top, samples, drop = FALSE] else NULL,
    names = DESeq2::resultsNames(dds),
    vst = vst,
    vst_pca = vst_pca,
    dispersions = disp_df,
    norm_counts = norm_counts,
    vst_dist = vst_dist,
    size_factors = size_factors,
    top_abs_lfc = top_abs_lfc,
    mean_sd = mean_sd_df,
    samples = samples,
    transformed_box = stack(as.data.frame(vst)),
    pvalue_fraction = mean(df$pvalue < 0.05, na.rm = TRUE),
    padj_fraction = mean(df$padj < 0.05, na.rm = TRUE),
    default_deg_fraction = mean(df$sig, na.rm = TRUE)
  )
}

ui <- shiny::fluidPage(
  shiny::tags$head(shiny::tags$style(shiny::HTML("
    body{background:#f6f7fb}.card{background:#fff;border:1px solid #d7dbe8;border-radius:10px;padding:14px;margin-bottom:14px}
    .help{display:inline-block;width:18px;height:18px;border-radius:50%;background:#4c566a;color:#fff;text-align:center;line-height:18px;font-size:12px;cursor:help;margin-left:6px}
    .formula{font-family:monospace;background:#f1f4fb;border:1px solid #d7dbe8;border-radius:8px;padding:10px}.legend-item{display:flex;gap:8px;align-items:center;margin-bottom:6px}
    .sw{width:12px;height:12px;border-radius:50%;display:inline-block}.muted{color:#58627a}
    .status-note{padding:10px 12px;border-radius:8px;background:#eef3fb;border:1px solid #d4deef;margin-bottom:12px}
    .top-input-note{margin-top:24px}
    .launcher-grid{background:#fbfcfe;border:1px solid #dfe6f3;border-radius:10px;padding:12px;margin-bottom:12px}
    .launcher-section{font-size:12px;font-weight:700;letter-spacing:.04em;text-transform:uppercase;color:#66728a;margin-bottom:8px}
    .launcher-hint{background:#f1f5fb;border:1px solid #d8e1ef;border-radius:8px;padding:10px 12px;color:#556178}
    .launcher-btn{margin-top:25px}
    .explorer-panel{background:#fbfcfe;border:1px solid #dfe6f3;border-radius:10px;padding:12px;margin-bottom:12px}
    .explorer-title{font-size:12px;font-weight:700;letter-spacing:.04em;text-transform:uppercase;color:#66728a;margin-bottom:6px}
    .explorer-note{font-size:12px;color:#5e6a81;margin-bottom:8px}
    .explorer-tight .form-group{margin-bottom:10px}
    .live-color-preview{width:18px;height:18px;border-radius:4px;border:1px solid #b7c1d6;background:#ffffff;margin-bottom:12px;flex:0 0 auto}
    .js-plotly-plot .modebar{opacity:1 !important}
  "))),
  shiny::tags$head(shiny::tags$script(shiny::HTML("
    window.openPlotlyGraphWindow = function(outputId, titleText) {
      var host = document.getElementById(outputId);
      if (!host) return;
      var plotDiv = (host.classList && host.classList.contains('js-plotly-plot')) ? host : host.querySelector('.js-plotly-plot');
      if (!plotDiv || !window.Plotly) return;
      var buildFigure = function() {
        try {
          if (window.Plotly.Plots && typeof window.Plotly.Plots.graphJson === 'function') {
            var graphObj = window.Plotly.Plots.graphJson(plotDiv, false, 'keepdata');
            if (graphObj && graphObj.data) {
              return {
                data: graphObj.data || [],
                layout: graphObj.layout || {}
              };
            }
          }
        } catch (e) {}
        return {
          data: JSON.parse(JSON.stringify(plotDiv.data || [])),
          layout: JSON.parse(JSON.stringify(plotDiv.layout || {}))
        };
      };
      var figure = buildFigure();
      if (!figure.data || !figure.data.length) return;
      var popup = window.open('', '_blank', 'width=1500,height=1000,resizable=yes,scrollbars=yes');
      if (!popup) return;
      var plotlyScript = null;
      var plotlyScripts = document.querySelectorAll('script[src]');
      for (var i = 0; i < plotlyScripts.length; i++) {
        var src = plotlyScripts[i].getAttribute('src') || '';
        if (src.toLowerCase().indexOf('plotly') !== -1) {
          plotlyScript = src;
          break;
        }
      }
      popup.document.open();
      popup.document.write('<!doctype html><html><head><meta charset=\"utf-8\"><title>' + (titleText || 'Graph explorer') + '</title><style>html,body{margin:0;padding:0;width:100%;height:100%;font-family:Arial,sans-serif;background:#f6f7fb;}#header{padding:12px 16px;background:#ffffff;border-bottom:1px solid #d7dbe8;}#plot{width:100vw;height:calc(100vh - 54px);}</style>' + (plotlyScript ? ('<script src=\"' + plotlyScript + '\"></scr' + 'ipt>') : '') + '</head><body><div id=\"header\">Use mouse wheel / trackpad to zoom, drag to pan, and hover to inspect nodes and edges.</div><div id=\"plot\"></div></body></html>');
      popup.document.close();
      popup.focus();
      var renderPopupPlot = function() {
        var target = popup.document.getElementById('plot');
        if (!target) return;
        var data = JSON.parse(JSON.stringify(figure.data || []));
        var layout = JSON.parse(JSON.stringify(figure.layout || {}));
        layout.width = null;
        layout.height = null;
        var config = {displaylogo: false, responsive: true, scrollZoom: true};
        var popupPlotly = popup.Plotly || window.Plotly;
        if (!popupPlotly) return;
        popupPlotly.newPlot(target, data, layout, config);
        popup.onresize = function() { popupPlotly.Plots.resize(target); };
      };
      if (popup.document.readyState === 'complete') {
        window.setTimeout(renderPopupPlot, 150);
      } else {
        popup.onload = renderPopupPlot;
      }
      window.setTimeout(renderPopupPlot, 300);
    };
    Shiny.addCustomMessageHandler('openExternalUrl', function(message) {
      if (!message || !message.url) return;
      window.open(message.url, '_blank', 'noopener=yes');
    });
    window.updateColorPreviewSwatches = function() {
      var previews = document.querySelectorAll('[data-color-input]');
      previews.forEach(function(preview) {
        var inputId = preview.getAttribute('data-color-input');
        var input = document.getElementById(inputId);
        if (!input) return;
        var val = (input.value || '').trim();
        preview.style.background = val || '#ffffff';
        preview.style.boxShadow = val ? 'inset 0 0 0 999px ' + val : 'none';
      });
    };
    window.bindColorPreviewSwatches = function() {
      var previews = document.querySelectorAll('[data-color-input]');
      previews.forEach(function(preview) {
        var inputId = preview.getAttribute('data-color-input');
        var input = document.getElementById(inputId);
        if (!input || input.dataset.colorPreviewBound === '1') return;
        var update = function() { window.updateColorPreviewSwatches(); };
        input.addEventListener('input', update);
        input.addEventListener('change', update);
        input.dataset.colorPreviewBound = '1';
      });
      window.updateColorPreviewSwatches();
    };
    document.addEventListener('DOMContentLoaded', function() {
      window.bindColorPreviewSwatches();
    });
    document.addEventListener('shiny:connected', function() {
      window.bindColorPreviewSwatches();
    });
    document.addEventListener('shiny:value', function() {
      window.bindColorPreviewSwatches();
    });
  "))),
  shiny::titlePanel("RNA Hypothesis Workbench: QC, DESeq2, WGCNA, SEMgraph, Enrichment, and Reporting"),
  shiny::div(
    class = "card",
    shiny::fluidRow(
      shiny::column(8,
        shiny::tags$h4("Dataset and metadata")
      ),
      shiny::column(4,
        shiny::tags$p(class = "muted", style = "text-align:right;margin-top:10px;", "Load data first, then move from QC to DESeq2, biological interpretation, causal prioritization, and finally reproducible reporting.")
      )
    ),
    shiny::div(
      class = "status-note",
      shiny::tags$strong("How to use this workbench: "),
      "start by loading counts and metadata, confirm sample quality in QC, define a sensible DESeq2 design, then use enrichment, WGCNA, and SEMgraph as progressively more interpretive layers. Earlier steps describe what changed in the samples; later steps help explain what those changes may mean biologically."
    ),
    shiny::fluidRow(
      shiny::column(6,
        shiny::div(class = "launcher-grid",
          shiny::div(class = "launcher-section", "Expression counts"),
          shiny::textInput("counts_path", "Counts file path", default_counts),
          shiny::fileInput("counts_upload", "Browse counts file", accept = c(".csv", ".tsv", ".txt"))
        )
      ),
      shiny::column(6,
        shiny::div(class = "launcher-grid",
          shiny::div(class = "launcher-section", "Sample metadata"),
          shiny::textInput("meta_path", "Metadata file path", default_meta),
          shiny::fileInput("meta_upload", "Browse metadata file", accept = c(".csv", ".tsv", ".txt"))
        )
      )
    ),
    shiny::fluidRow(
      shiny::column(2,
        shiny::numericInput("min_count", "Minimum count", 10, min = 0)
      ),
      shiny::column(3,
        shiny::numericInput("min_samples", "Minimum samples at threshold", 2, min = 1)
      ),
      shiny::column(3,
        shiny::selectInput("qc_group", "QC color grouping", choices = "treatment")
      ),
      shiny::column(2,
        shiny::div(class = "launcher-btn",
          shiny::actionButton("load_data", "Load dataset")
        )
      ),
      shiny::column(2,
        shiny::div(class = "launcher-hint",
          "Set the low-count filter here, then use the QC grouping to see the same biological question consistently across PCA, distances, and later diagnostics."
        )
      )
    )
  ),
  shiny::verbatimTextOutput("status"),
  shiny::tabsetPanel(
        shiny::tabPanel("QC",
          DT::dataTableOutput("sample_table"),
          plotly::plotlyOutput("lib_plot", height = "260px"),
          plotly::plotlyOutput("pca_plot", height = "420px"),
          shiny::fluidRow(
            shiny::column(9,
              shiny::div(shiny::tags$strong("Distances"), shiny::tags$span("?", class = "help", title = "Distances are computed on filtered log2(CPM + 1). Axis colors follow the selected metadata grouping.")),
              plotly::plotlyOutput("dist_plot", height = "460px")
            ),
            shiny::column(3, shiny::uiOutput("dist_legend"))
          )
        ),
        shiny::tabPanel("Design Builder",
          shiny::fluidRow(
            shiny::column(5,
              shiny::div(class = "card",
                shiny::checkboxGroupInput("de_samples", "Samples used in DESeq2", choices = NULL),
                shiny::selectizeInput("design_terms", "Design terms", choices = NULL, multiple = TRUE),
                shiny::checkboxInput("design_interaction", "Add interaction term", FALSE),
                shiny::textInput("design_name", "Design name", "Primary design"),
                shiny::actionButton("save_design", "Save design")
              )
            ),
            shiny::column(7,
              shiny::div(class = "card",
                shiny::tags$h4("Formula preview"),
                shiny::div(class = "formula", shiny::textOutput("formula_preview")),
                shiny::tags$p(class = "muted", shiny::textOutput("formula_help"))
              ),
              shiny::div(class = "card",
                shiny::tags$h4("Saved designs"),
                shiny::tags$span("?", class = "help", title = "A design formula tells DESeq2 which metadata columns should explain changes in expression. Save a simple baseline model first, then compare it to an adjusted model such as treatment + batch."),
                DT::dataTableOutput("design_table")
              )
            )
          )
        ),
        shiny::tabPanel("DESeq2 Results",
          shiny::fluidRow(
            shiny::column(4,
              shiny::div(class = "card",
                shiny::selectInput("active_design", "Saved design", choices = character()),
                shiny::selectInput("contrast_var", "Comparison variable", choices = character()),
                shiny::tags$div(style = "margin:8px 0 4px 0;font-weight:600;", "DEG thresholds"),
                shiny::numericInput("de_alpha", "Adjusted p-value cutoff (alpha)", value = 0.05, min = 0.001, max = 0.5, step = 0.005),
                shiny::numericInput("de_lfc_cutoff", "Absolute log2 fold-change cutoff", value = 0, min = 0, max = 10, step = 0.25),
                shiny::div(
                  shiny::tags$strong("Contrast direction"),
                  shiny::tags$span(
                    "?",
                    class = "help",
                    title = "Numerator is the group you want to measure relative to the reference. Denominator is the reference group. For example, numerator = treated and denominator = control gives log2 fold changes for treated versus control."
                  )
                ),
                shiny::selectInput("contrast_num", "Numerator", choices = character()),
                shiny::selectInput("contrast_den", "Denominator", choices = character()),
                shiny::actionButton("run_de", "Run DESeq2"),
                shiny::tags$p(class = "muted", "This is the central inferential step of the workflow. Choose a saved design, define the biological comparison carefully, and set the DEG thresholds that later enrichment and SEMgraph views will treat as sample-supported evidence.")
              ),
              shiny::div(class = "card", DT::dataTableOutput("run_history"))
            ),
            shiny::column(8,
              shiny::tabsetPanel(
                shiny::tabPanel("1. Fit Summary",
                  shiny::div(class = "card", shiny::uiOutput("de_summary"))
                ),
                shiny::tabPanel("2. Transform Checks",
                  shiny::div(class = "card",
                    shiny::tags$h4("Variance-stabilized transform"),
                    shiny::tags$span("?", class = "help", title = "DESeq2 fits the statistical model on raw counts, but transformed values such as VST are useful for diagnostics and visualization because they reduce the dependence of variance on mean expression."),
                    shiny::uiOutput("transform_summary")
                  ),
                  shiny::div(class = "card", plotly::plotlyOutput("vst_pca_plot", height = "340px")),
                  shiny::div(class = "card", plotly::plotlyOutput("vst_box_plot", height = "320px")),
                  shiny::div(class = "card", plotly::plotlyOutput("mean_sd_plot", height = "320px"))
                ),
                shiny::tabPanel("3. Statistical Diagnostics",
                  shiny::div(class = "card", shiny::uiOutput("diag_summary")),
                  shiny::div(class = "card", plotly::plotlyOutput("dispersion_plot", height = "320px")),
                  shiny::div(class = "card", shiny::uiOutput("dispersion_help")),
                  shiny::div(class = "card", plotly::plotlyOutput("size_factor_plot", height = "280px")),
                  shiny::div(class = "card", plotly::plotlyOutput("pvalue_hist", height = "300px")),
                  shiny::div(class = "card", plotly::plotlyOutput("padj_hist", height = "300px")),
                  shiny::div(class = "card", plotly::plotlyOutput("pvalue_ma_plot", height = "320px")),
                  shiny::div(class = "card", shiny::uiOutput("improvement_suggestions"))
                ),
                shiny::tabPanel("4. Differential Signals",
                  shiny::div(class = "card", plotly::plotlyOutput("ma_plot", height = "320px")),
                  shiny::div(class = "card", plotly::plotlyOutput("volcano_plot", height = "320px")),
                  shiny::div(class = "card",
                    shiny::fluidRow(
                      shiny::column(10, plotly::plotlyOutput("heat_plot", height = "460px")),
                      shiny::column(2, shiny::uiOutput("heat_legend"))
                    )
                  ),
                  shiny::div(class = "card", DT::dataTableOutput("de_table"))
                )
              )
            )
          )
        ),
        shiny::tabPanel("Pathway Enrichment",
          shiny::fluidRow(
            shiny::column(4,
              shiny::div(class = "card",
                shiny::tags$h4("Enrichment setup"),
                shiny::checkboxGroupInput("enrich_methods", "Methods", choices = c("GO ORA", "KEGG ORA", "GO GSEA", "KEGG GSEA"), selected = c("GO ORA", "GO GSEA")),
                shiny::numericInput("enrich_padj", "Adjusted p-value cutoff", value = 0.05, min = 0.001, max = 0.5, step = 0.01),
                shiny::numericInput("enrich_lfc", "Minimum absolute log2 fold change for ORA genes", value = 1, min = 0, step = 0.25),
                shiny::numericInput("enrich_min_gs", "Minimum pathway size", value = 10, min = 5, step = 1),
                shiny::numericInput("enrich_max_gs", "Maximum pathway size", value = 500, min = 20, step = 10),
                shiny::selectInput("enrich_adjust_method", "Multiple-testing adjustment", choices = c("BH", "BY", "bonferroni", "holm"), selected = "BH"),
                shiny::numericInput("enrich_top_n", "Top pathways to emphasize", value = 15, min = 5, max = 50, step = 1),
                shiny::checkboxInput("refresh_enrichment_cache", "Refresh cached enrichment databases/results", value = FALSE),
                shiny::actionButton("run_enrichment", "Run enrichment"),
                shiny::tags$p(class = "muted", "Use this after DESeq2 to ask what the changing genes are doing biologically. ORA focuses on genes that passed your current DEG cutoffs, while GSEA keeps the full ranked result and can detect broader coordinated pathway shifts. Cached database-heavy results are reused unless you ask to refresh them.")
              ),
              shiny::div(class = "card", shiny::uiOutput("enrich_help"))
            ),
            shiny::column(8,
              shiny::tabsetPanel(
                shiny::tabPanel("1. Summary",
                  shiny::div(class = "card", shiny::uiOutput("enrich_summary")),
                  shiny::div(class = "card", shiny::uiOutput("enrich_mechanism")),
                  shiny::div(class = "card", DT::dataTableOutput("consensus_table"))
                ),
                shiny::tabPanel("2. Diagnostics",
                  shiny::div(class = "card", plotly::plotlyOutput("rank_hist_plot", height = "300px")),
                  shiny::div(class = "card", plotly::plotlyOutput("method_balance_plot", height = "300px"))
                ),
                shiny::tabPanel("3. By Method",
                  shiny::div(class = "card",
                    shiny::selectInput("enrich_view_method", "Inspect one method", choices = character()),
                    shiny::tags$p(class = "muted", "This view shows the pathways returned by exactly one selected method.")
                  ),
                  shiny::div(class = "card", plotly::plotlyOutput("enrich_method_plot", height = "420px")),
                  shiny::div(class = "card", DT::dataTableOutput("enrich_method_table"))
                ),
                shiny::tabPanel("4. Consensus",
                  shiny::div(class = "card",
                    shiny::tags$p("Consensus does not merge p-values across methods. It highlights pathways that recur across multiple selected methods, then orders them by overlap count and significance.")
                  ),
                  shiny::div(class = "card", plotly::plotlyOutput("enrich_consensus_plot", height = "420px")),
                  shiny::div(class = "card", DT::dataTableOutput("enrich_table"))
                )
              )
            )
          )
        ),
        shiny::tabPanel("WGCNA Modules",
          shiny::fluidRow(
            shiny::column(4,
              shiny::div(class = "card",
                shiny::tags$h4("Network setup"),
                shiny::radioButtons("wgcna_power_mode", "Soft-threshold power", choices = c("Automatic" = "auto", "Manual" = "manual"), selected = "auto", inline = TRUE),
                shiny::numericInput("wgcna_manual_power", "Manual power", value = 6, min = 1, max = 30, step = 1),
                shiny::numericInput("wgcna_top_var_genes", "Most variable genes to include", value = 4000, min = 500, step = 250),
                shiny::numericInput("wgcna_min_module_size", "Minimum module size", value = 30, min = 10, step = 5),
                shiny::sliderInput("wgcna_deep_split", "Cluster splitting sensitivity", min = 0, max = 4, value = 2, step = 1),
                shiny::sliderInput("wgcna_merge_cut_height", "Module merge threshold", min = 0.05, max = 0.45, value = 0.2, step = 0.05),
                shiny::checkboxGroupInput("wgcna_traits", "Traits to relate to modules", choices = character()),
                shiny::actionButton("run_wgcna", "Run WGCNA"),
                shiny::tags$p(class = "muted", "Use WGCNA after you understand the DE result. This layer asks which genes move together across samples, even if not every gene is strongly differential. Lower merge threshold keeps modules more separate; higher values merge more similar modules. Increase split sensitivity to create more modules from the dendrogram.")
              ),
              shiny::div(class = "card", shiny::uiOutput("wgcna_help"))
            ),
            shiny::column(8,
              shiny::tabsetPanel(
                shiny::tabPanel("1. Power",
                  shiny::div(class = "card", shiny::uiOutput("wgcna_summary")),
                  shiny::div(class = "card", plotly::plotlyOutput("wgcna_power_plot", height = "320px"))
                ),
                shiny::tabPanel("2. Module Tree",
                  shiny::div(class = "card", shiny::uiOutput("wgcna_tree_help")),
                  shiny::div(class = "card", shiny::plotOutput("wgcna_dendrogram_plot", height = "520px"))
                ),
                shiny::tabPanel("3. Module Traits",
                  shiny::div(class = "card",
                    shiny::selectInput("wgcna_trait_view", "Focus trait/contrast", choices = character()),
                    shiny::tags$p(class = "muted", "Use this to inspect only the most relevant trait columns instead of every possible module-trait comparison.")
                  ),
                  shiny::div(class = "card", plotly::plotlyOutput("wgcna_trait_heatmap", height = "420px")),
                  shiny::div(class = "card", plotly::plotlyOutput("wgcna_trait_focus_plot", height = "320px"))
                ),
                shiny::tabPanel("4. Modules",
                  shiny::div(class = "card",
                    shiny::selectInput("wgcna_module_view", "Inspect module", choices = character()),
                    shiny::sliderInput("wgcna_edge_quantile", "Show only stronger module edges", min = 0.7, max = 0.99, value = 0.9, step = 0.01),
                    shiny::tags$p(class = "muted", "Raising the threshold shows only the strongest co-expression links. Lowering it reveals more of the internal module structure.")
                  ),
                  shiny::div(class = "card",
                    shiny::tags$h4("Module biological annotation"),
                    shiny::textInput("wgcna_module_annotation", "Interpretation label", value = ""),
                shiny::textAreaInput("wgcna_module_annotation_notes", "Notes", value = "", rows = 3, resize = "vertical"),
                shiny::actionButton("save_wgcna_annotation", "Save annotation"),
                shiny::tags$p(class = "muted", "Use this to assign a plain-language biological meaning to the selected module, such as 'ECM remodeling', 'cell cycle', or 'interferon response'.")
              ),
              shiny::div(class = "card",
                shiny::checkboxInput("refresh_module_cache", "Refresh cached module/program enrichment", value = FALSE),
                shiny::tags$p(class = "muted", "If unchecked, module and program downstream enrichment will reuse cached results when the same gene set and parameters were already analyzed.")
              ),
                  shiny::div(class = "card", plotly::plotlyOutput("wgcna_module_size_plot", height = "320px")),
                  shiny::div(class = "card", shiny::uiOutput("wgcna_module_help")),
                  shiny::div(class = "card", plotly::plotlyOutput("wgcna_connectivity_plot", height = "280px")),
                  shiny::div(class = "card", plotly::plotlyOutput("wgcna_membership_plot", height = "320px")),
                  shiny::div(class = "card",
                    shiny::downloadButton("download_wgcna_module_graphml", "Download module network as GraphML"),
                    plotly::plotlyOutput("wgcna_module_network_plot", height = "520px")
                  ),
                  shiny::div(class = "card",
                    shiny::tags$h4("Top hub genes in selected module"),
                    shiny::tags$p(class = "muted", "This table is focused. It shows the strongest candidate hub genes only for the currently selected module, ranked by within-module connectivity."),
                    DT::dataTableOutput("wgcna_hubgene_table")
                  ),
                  shiny::div(class = "card",
                    shiny::tags$h4("All genes and module assignments"),
                    shiny::tags$p(class = "muted", "This table is global. It lists every analyzed gene, its assigned module, and its within-module connectivity so you can compare across modules, not just within the selected one."),
                    DT::dataTableOutput("wgcna_module_table")
                  ),
                  shiny::div(class = "card",
                    shiny::tags$h4("Selected module downstream analysis"),
                    shiny::tags$p(class = "muted", "This section asks what the selected module may represent biologically. It summarizes the module, relates it back to traits, and runs pathway enrichment on the module genes."),
                    shiny::uiOutput("wgcna_module_downstream_help")
                  ),
                  shiny::div(class = "card", shiny::uiOutput("wgcna_module_summary")),
                  shiny::div(class = "card", plotly::plotlyOutput("wgcna_module_enrichment_plot", height = "380px")),
                  shiny::div(class = "card", DT::dataTableOutput("wgcna_module_enrichment_table"))
                ),
                shiny::tabPanel("5. Programs",
                  shiny::div(class = "card",
                    shiny::tags$h4("Suggested module-combination programs"),
                    shiny::selectInput("wgcna_program_suggestion", "Suggested combination", choices = character()),
                    shiny::textInput("wgcna_program_name", "Program label", value = ""),
                    shiny::actionButton("create_wgcna_program", "Create combined program"),
                    shiny::tags$p(class = "muted", "These suggestions do not rewrite WGCNA modules. They create a higher-level biological program when modules look like opposite arms of the same process.")
                  ),
                  shiny::div(class = "card", shiny::uiOutput("wgcna_program_help")),
                  shiny::div(class = "card", DT::dataTableOutput("wgcna_program_suggestion_table")),
                  shiny::div(class = "card",
                    shiny::selectInput("wgcna_program_view", "Inspect saved program", choices = character())
                  ),
                  shiny::div(class = "card", shiny::uiOutput("wgcna_program_summary")),
                  shiny::div(class = "card", plotly::plotlyOutput("wgcna_program_diag_plot", height = "360px")),
                  shiny::div(class = "card", plotly::plotlyOutput("wgcna_program_trait_plot", height = "340px")),
                  shiny::div(class = "card", plotly::plotlyOutput("wgcna_program_enrichment_plot", height = "380px")),
                  shiny::div(class = "card", DT::dataTableOutput("wgcna_program_table"))
                )
              )
            )
          )
        ),
        shiny::tabPanel("SEMgraph Causal",
          shiny::fluidRow(
            shiny::column(4,
              shiny::div(class = "card",
                shiny::tags$h4("Causal setup"),
                shiny::selectInput("semgraph_unit", "Analysis unit", choices = character()),
                shiny::selectInput("semgraph_group", "Binary group for perturbation testing", choices = stats::setNames("", "")),
                shiny::radioButtons("semgraph_display_mode", "Display focus", choices = c("Strict DE-only" = "de_only", "DE + minimal connectors" = "de_connectors", "Full local neighborhood" = "full_local"), selected = "de_connectors"),
                shiny::numericInput("semgraph_model_max_genes", "Maximum genes entering SEMgraph fit", value = 300, min = 100, max = 2000, step = 100),
                shiny::numericInput("semgraph_annotation_max_genes", "Maximum genes receiving external annotation", value = 180, min = 50, max = 1000, step = 25),
                shiny::radioButtons("semgraph_table_scope", "Table evidence scope", choices = c("DE-only" = "de_only", "DE + key regulators" = "de_plus_regulators", "Broad relevant evidence" = "evidence_context"), selected = "evidence_context"),
                shiny::numericInput("semgraph_table_max_rows", "Maximum table rows", value = 120, min = 20, max = 500, step = 20),
                shiny::checkboxInput("semgraph_use_reactome_live", "Use live Reactome lookup if cache is missing", value = FALSE),
                shiny::checkboxInput("semgraph_use_string_live", "Use live STRING lookup if cache is missing", value = FALSE),
                shiny::fluidRow(
                  shiny::column(6, shiny::actionButton("update_reactome_cache", "Update Reactome cache")),
                  shiny::column(6, shiny::actionButton("semgraph_download_string_local", "Download local STRING file"))
                ),
                shiny::checkboxInput("refresh_semgraph_cache", "Refresh cached pharmacology/context databases", value = FALSE),
                shiny::checkboxInput("refresh_semgraph_fit", "Refit SEMgraph model instead of reusing cached fit", value = FALSE),
                shiny::actionButton("run_semgraph", "Run causal analysis"),
                shiny::tags$p(class = "muted", "Use SEMgraph only after you have a believable DESeq2 result and a biologically interpretable module or program. This layer does not replace DE or WGCNA: it adds a directional hypothesis layer on top of them. Large units are prioritized before fitting so the model stays tractable, and external annotations are limited to the top-ranked genes. Use the buttons above to prepare or refresh local validation resources that affect the SEMgraph summary and downstream validation views.")
              ),
              shiny::div(class = "card", shiny::uiOutput("semgraph_help"))
            ),
            shiny::column(8,
              shiny::tabsetPanel(
                shiny::tabPanel("1. Summary",
                  shiny::div(class = "card", shiny::uiOutput("semgraph_summary"))
                ),
                shiny::tabPanel("2. Diagnostics",
                  shiny::div(class = "card", shiny::uiOutput("semgraph_diag_help")),
                  shiny::div(class = "card", plotly::plotlyOutput("semgraph_support_plot", height = "280px")),
                  shiny::div(class = "card", plotly::plotlyOutput("semgraph_degree_scatter", height = "320px")),
                  shiny::div(class = "card", plotly::plotlyOutput("semgraph_role_support_plot", height = "300px"))
                ),
                shiny::tabPanel("3. Causal Graph",
                  shiny::div(class = "card", shiny::uiOutput("semgraph_graph_help")),
                  shiny::div(class = "card",
                    shiny::fluidRow(
                      shiny::column(6,
                        shiny::div(class = "explorer-panel explorer-tight",
                          shiny::div(class = "explorer-title", "1. Choose the graph view"),
                          shiny::div(class = "explorer-note", "Start with an overview, then switch to a local neighborhood or directed paths when you want to inspect a smaller mechanism more closely."),
                          shiny::fluidRow(
                            shiny::column(6,
                              shiny::selectInput(
                                "semgraph_layout",
                                "Layout",
                                choices = c(
                                  "Spring layout (Fruchterman-Reingold)" = "fr",
                                  "Spring layout (Kamada-Kawai)" = "kk",
                                  "Graphopt layout" = "graphopt",
                                  "Large-graph layout" = "lgl",
                                  "Circle layout" = "circle"
                                ),
                                selected = "fr"
                              )
                            ),
                            shiny::column(6,
                              shiny::radioButtons(
                                "semgraph_graph_scope",
                                "Exploration mode",
                                choices = c(
                                  "Overview" = "overview",
                                  "Neighborhood by hops" = "neighborhood",
                                  "Directed paths" = "paths"
                                ),
                                selected = "overview"
                              )
                            )
                          )
                        )
                      ),
                      shiny::column(6,
                        shiny::div(class = "explorer-panel explorer-tight",
                          shiny::div(class = "explorer-title", "2. Open a larger explorer"),
                          shiny::div(class = "explorer-note", "Use the in-app larger explorer for quick expansion, or open a separate browser window when you want maximum space for zooming and panning."),
                          shiny::fluidRow(
                            shiny::column(6,
                              shiny::div(class = "top-input-note",
                                shiny::actionButton("semgraph_open_large", "Open larger graph explorer")
                              )
                            ),
                            shiny::column(6,
                              shiny::div(class = "top-input-note",
                                shiny::tags$button(
                                  type = "button",
                                  class = "btn btn-default action-button",
                                  onclick = "openPlotlyGraphWindow('semgraph_graph_plot', 'SEMgraph Explorer');",
                                  "Open in separate browser window"
                                )
                              )
                            )
                          )
                        )
                      )
                    ),
                    shiny::div(class = "explorer-panel explorer-tight",
                      shiny::div(class = "explorer-title", "3. Build a focused graph"),
                      shiny::div(class = "explorer-note", "Use this single control area to decide which part of the SEMgraph should be visible. Choose the exploration mode, set the graph size, restrict to STRING-supported edges if needed, then define the anchor genes and categories that should remain in view."),
                      shiny::uiOutput("semgraph_graph_scope_ui")
                    )
                  ),
                  shiny::div(class = "card",
                    shiny::downloadButton("download_semgraph_graphml", "Download displayed SEMgraph as GraphML"),
                    plotly::plotlyOutput("semgraph_graph_plot", height = "560px")
                  ),
                  shiny::div(class = "card",
                    shiny::tags$h4("Visible graph edges and STRING action support"),
                    shiny::tags$p(class = "muted", "This table follows the graph you are currently looking at. It shows only the edges in the visible SEMgraph view, together with STRING pair support, STRING action details when available, and whether STRING action direction agrees with the SEMgraph arrow."),
                    DT::dataTableOutput("semgraph_graph_edge_table")
                  ),
                  shiny::div(class = "card", DT::dataTableOutput("semgraph_graph_analysis_table"))
                ),
                shiny::tabPanel("4. Regulators",
                  shiny::div(class = "card", shiny::uiOutput("semgraph_regulator_help")),
                  shiny::div(class = "card", DT::dataTableOutput("semgraph_degree_table"))
                ),
                shiny::tabPanel("5. Effect Estimates",
                  shiny::div(class = "card", shiny::uiOutput("semgraph_ace_help")),
                  shiny::div(class = "card", plotly::plotlyOutput("semgraph_effect_proxy_plot", height = "340px")),
                  shiny::div(class = "card", DT::dataTableOutput("semgraph_ace_table"))
                ),
                shiny::tabPanel("6. Pathway Context",
                  shiny::div(class = "card", shiny::uiOutput("semgraph_pathway_help")),
                  shiny::div(class = "card", plotly::plotlyOutput("semgraph_pathway_summary_plot", height = "340px")),
                  shiny::div(class = "card", plotly::plotlyOutput("semgraph_pathway_gene_plot", height = "460px")),
                  shiny::div(class = "card", DT::dataTableOutput("semgraph_pathway_table"))
                ),
                shiny::tabPanel("7. Edge Validation",
                  shiny::div(class = "card", shiny::uiOutput("semgraph_validation_help")),
                  shiny::div(class = "card", plotly::plotlyOutput("semgraph_validation_plot", height = "320px")),
                  shiny::div(class = "card", DT::dataTableOutput("semgraph_edge_table"))
                )
              )
            )
          )
        ),
        shiny::tabPanel("Reports",
          shiny::fluidRow(
            shiny::column(4,
              shiny::div(class = "card",
                shiny::tags$h4("Reproducibility"),
                shiny::tags$p(class = "muted", "Use these exports after the analysis story is coherent. They are meant to help someone else understand what was run, why those settings mattered, and how to recreate or audit the same workflow outside the app."),
                shiny::downloadButton("download_repro_bundle", "Download reproducibility R Markdown bundle"),
                shiny::tags$br(),
                shiny::tags$br(),
                shiny::downloadButton("download_html_bundle", "Download HTML report bundle (.zip)")
              ),
              shiny::div(class = "card",
                shiny::tags$h4("Package versions"),
                shiny::tags$p(class = "muted", "These are the package versions used by this running app environment."),
                DT::dataTableOutput("package_versions_table")
              )
            ),
            shiny::column(8,
              shiny::div(class = "card",
                shiny::tags$h4("What the report bundles contain"),
                shiny::tags$ul(
                  shiny::tags$li("An analysis-action log describing what the user actually ran in the app."),
                  shiny::tags$li("A package-version summary for recreating the environment elsewhere."),
                  shiny::tags$li("Snapshot tables and graph files for the currently available results."),
                  shiny::tags$li("A reproducibility-focused R Markdown bundle with explanatory text and code chunks that reconstruct the saved report outputs from the exported snapshot."),
                  shiny::tags$li("A shareable HTML bundle with explanations, collapsible code sections, plots, and table previews with full CSV exports.")
                )
              ),
              shiny::div(class = "card",
                shiny::tags$h4("Analysis log"),
                DT::dataTableOutput("analysis_log_table")
              )
            )
          )
        )
  )
)

server <- function(input, output, session) {
  saved <- shiny::reactiveVal(data.frame(name = character(), formula = character(), terms = character(), stringsAsFactors = FALSE))
  history <- shiny::reactiveVal(data.frame(design = character(), formula = character(), contrast = character(), sig = integer(), stringsAsFactors = FALSE))
  wgcna_annotations <- shiny::reactiveVal(data.frame(module = character(), annotation = character(), notes = character(), stringsAsFactors = FALSE))
  wgcna_programs <- shiny::reactiveVal(data.frame(name = character(), modules = character(), rationale = character(), stringsAsFactors = FALSE))
  analysis_log <- shiny::reactiveVal(data.frame(timestamp = character(), step = character(), params_json = character(), stringsAsFactors = FALSE))

  ds <- shiny::eventReactive(input$load_data, {
    counts_path <- resolve_input_path(input$counts_upload, input$counts_path)
    meta_path <- resolve_input_path(input$meta_upload, input$meta_path)
    cur <- read_inputs(counts_path, meta_path, input$min_count, input$min_samples)
    ann <- setdiff(names(cur$meta), "run_id")
    shiny::updateSelectInput(session, "qc_group", choices = ann, selected = if ("treatment" %in% ann) "treatment" else ann[1L])
    shiny::updateCheckboxGroupInput(session, "de_samples", choices = colnames(cur$counts), selected = colnames(cur$counts))
    shiny::updateSelectizeInput(session, "design_terms", choices = ann, selected = if ("treatment" %in% ann) "treatment" else ann[1L], server = TRUE)
    cur
  }, ignoreNULL = FALSE)
  shiny::observeEvent(input$load_data, {
    analysis_log(append_analysis_log(analysis_log(), "load_data", list(
      counts_path = input$counts_path,
      metadata_path = input$meta_path,
      min_count = input$min_count,
      min_samples = input$min_samples
    )))
  }, ignoreInit = TRUE)

  formula_now <- shiny::reactive(build_formula(input$design_terms %||% character(), isTRUE(input$design_interaction)))
  output$formula_preview <- shiny::renderText(formula_now())
  output$formula_help <- shiny::renderText({
    cur <- ds()
    design_note(cur$meta, input$design_terms %||% character(), isTRUE(input$design_interaction), input$de_samples %||% colnames(cur$counts), cur$filtered)
  })
  output$status <- shiny::renderText({
    cur <- ds()
    sprintf("Loaded %s samples and %s genes. Filtered set retains %s genes.", ncol(cur$counts), nrow(cur$counts), nrow(cur$filtered))
  })
  output$package_versions_table <- DT::renderDataTable({
    make_downloadable_table(package_versions_df(), page_length = 15)
  }, server = FALSE)
  output$analysis_log_table <- DT::renderDataTable({
    make_downloadable_table(analysis_log(), page_length = 12)
  }, server = FALSE)
  output$sample_table <- DT::renderDataTable({
    cur <- ds()
    tab <- data.frame(run_id = cur$meta$run_id, library_size = colSums(cur$counts), detected_genes = colSums(cur$counts > 0), cur$meta[setdiff(names(cur$meta), "run_id")], check.names = FALSE)
    make_downloadable_table(tab, page_length = 8)
  }, server = FALSE)
  output$lib_plot <- plotly::renderPlotly({
    cur <- ds(); grp <- cur$meta[[input$qc_group]]
    pd <- data.frame(run_id = cur$meta$run_id, lib = colSums(cur$counts), grp = grp)
    plotly::layout(plotly::plot_ly(pd, x = ~run_id, y = ~lib, color = ~grp, type = "bar", text = ~paste(run_id, "<br>", grp, "<br>", lib), hoverinfo = "text"), title = "Library size", xaxis = list(title = ""), yaxis = list(title = "Counts"))
  })
  output$pca_plot <- plotly::renderPlotly({
    cur <- ds(); grp <- cur$meta[[input$qc_group]]; scores <- as.data.frame(cur$pca$x[, 1:2, drop = FALSE]); scores$run_id <- rownames(scores); scores$grp <- grp; vars <- (cur$pca$sdev ^ 2) / sum(cur$pca$sdev ^ 2)
    plotly::layout(plotly::plot_ly(scores, x = ~PC1, y = ~PC2, color = ~grp, type = "scatter", mode = "markers+text", text = ~run_id, textposition = "top center"), title = "PCA", xaxis = list(title = sprintf("PC1 (%.1f%%)", vars[1L] * 100)), yaxis = list(title = sprintf("PC2 (%.1f%%)", vars[2L] * 100)))
  })
  output$dist_plot <- plotly::renderPlotly({
    cur <- ds()
    ord <- if (ncol(cur$dist) > 1L) stats::hclust(stats::as.dist(cur$dist))$order else 1L
    ids <- colnames(cur$dist)[ord]
    grp <- cur$meta[[input$qc_group]]
    names(grp) <- cur$meta$run_id
    grp <- grp[ids]
    cols <- group_colors(grp)
    ticks <- tick_text(ids, grp, cols)
    z <- cur$dist[ids, ids, drop = FALSE]
    ann_plot <- build_metadata_annotation_plot(cur$meta, ids, title = "Sample metadata")
    main_plot <- plotly::plot_ly(
      x = ids,
      y = rev(ids),
      z = z[nrow(z):1, , drop = FALSE],
      type = "heatmap",
      colorscale = "Blues"
    ) |>
      plotly::layout(
        title = "Sample distance heatmap",
        xaxis = list(title = "", tickmode = "array", tickvals = ids, ticktext = ticks),
        yaxis = list(title = "", tickmode = "array", tickvals = rev(ids), ticktext = rev(ticks))
      )
    if (is.null(ann_plot)) {
      main_plot
    } else {
      plotly::subplot(ann_plot, main_plot, nrows = 2, shareX = TRUE, heights = c(0.18, 0.82), titleY = TRUE)
    }
  })
  output$dist_legend <- shiny::renderUI({
    cur <- ds(); cols <- group_colors(cur$meta[[input$qc_group]])
    shiny::tagList(shiny::tags$h4("Group colors"), lapply(names(cols), function(nm) shiny::tags$div(class = "legend-item", shiny::tags$span(class = "sw", style = paste("background:", cols[[nm]], ";")), shiny::tags$span(nm))))
  })

  shiny::observeEvent(input$save_design, {
    nm <- trimws(input$design_name); if (!nzchar(nm)) nm <- paste("Design", nrow(saved()) + 1L)
    tbl <- saved(); tbl <- tbl[tbl$name != nm, , drop = FALSE]
    tbl <- rbind(tbl, data.frame(name = nm, formula = formula_now(), terms = paste(input$design_terms %||% character(), collapse = " + "), stringsAsFactors = FALSE))
    saved(tbl); shiny::updateSelectInput(session, "active_design", choices = tbl$name, selected = nm)
    analysis_log(append_analysis_log(analysis_log(), "save_design", list(
      design_name = nm,
      formula = formula_now(),
      terms = input$design_terms %||% character(),
      interaction = isTRUE(input$design_interaction)
    )))
  })
  output$design_table <- DT::renderDataTable(make_downloadable_table(saved(), page_length = 5), server = FALSE)
  shiny::observe({
    tbl <- saved()
    selected <- if (!is.null(input$active_design) && input$active_design %in% tbl$name) input$active_design else if (nrow(tbl)) tbl$name[1L] else character()
    shiny::updateSelectInput(session, "active_design", choices = tbl$name, selected = selected)
  })
  shiny::observe({
    cur <- ds(); tbl <- saved(); if (!nrow(tbl)) return()
    row <- tbl[tbl$name == input$active_design, , drop = FALSE]; terms <- trimws(unlist(strsplit(row$terms, "\\+"))); terms <- terms[nzchar(terms)]
    meta_sub <- cur$meta[match(input$de_samples %||% colnames(cur$counts), cur$meta$run_id), , drop = FALSE]
    fac_terms <- terms[vapply(terms, function(t) t %in% names(meta_sub) && is.factor(meta_sub[[t]]) && length(unique(meta_sub[[t]])) >= 2L, logical(1))]
    if (!length(fac_terms)) {
      shiny::updateSelectInput(session, "contrast_var", choices = ""); shiny::updateSelectInput(session, "contrast_num", choices = ""); shiny::updateSelectInput(session, "contrast_den", choices = "")
    } else {
      var <- fac_terms[1L]; lv <- levels(meta_sub[[var]])
      shiny::updateSelectInput(session, "contrast_var", choices = fac_terms, selected = var)
      shiny::updateSelectInput(session, "contrast_num", choices = lv, selected = lv[min(2L, length(lv))])
      shiny::updateSelectInput(session, "contrast_den", choices = lv, selected = lv[1L])
    }
  })
  shiny::observeEvent(input$contrast_var, {
    cur <- ds(); if (!nzchar(input$contrast_var)) return()
    meta_sub <- cur$meta[match(input$de_samples %||% colnames(cur$counts), cur$meta$run_id), , drop = FALSE]
    lv <- levels(meta_sub[[input$contrast_var]])
    shiny::updateSelectInput(session, "contrast_num", choices = lv, selected = lv[min(2L, length(lv))])
    shiny::updateSelectInput(session, "contrast_den", choices = lv, selected = lv[1L])
  }, ignoreNULL = FALSE)

  de <- shiny::eventReactive(input$run_de, {
    shiny::withProgress(message = "Running DESeq2", detail = "Fitting model and extracting results. Please wait.", value = 0, {
      cur <- ds(); tbl <- saved(); row <- tbl[tbl$name == input$active_design, , drop = FALSE]
      shiny::incProgress(0.2, detail = "Preparing design and count data")
      res <- run_de(cur$filtered, cur$meta, input$de_samples %||% colnames(cur$counts), row$formula, cur$gene_map, if (nzchar(input$contrast_var)) input$contrast_var else NULL, input$contrast_num, input$contrast_den)
      shiny::incProgress(0.9, detail = "Finalizing diagnostics")
      current_sig <- sum(deg_mask(res$df, alpha = input$de_alpha %||% 0.05, lfc_cutoff = input$de_lfc_cutoff %||% 0), na.rm = TRUE)
      hist <- history(); hist <- rbind(hist, data.frame(design = row$name, formula = row$formula, contrast = res$label, sig = current_sig, stringsAsFactors = FALSE)); history(hist); res
    })
  })
  shiny::observeEvent(input$run_de, {
    analysis_log(append_analysis_log(analysis_log(), "run_deseq2", list(
      active_design = input$active_design,
      comparison_variable = input$contrast_var,
      numerator = input$contrast_num,
      denominator = input$contrast_den,
      selected_samples = input$de_samples %||% character(),
      alpha = input$de_alpha %||% 0.05,
      abs_log2fc_cutoff = input$de_lfc_cutoff %||% 0
    )))
  }, ignoreInit = TRUE)
  output$run_history <- DT::renderDataTable(make_downloadable_table(history(), page_length = 5), server = FALSE)
  output$de_summary <- shiny::renderUI({
    res <- de()
    alpha <- input$de_alpha %||% 0.05
    lfc_cut <- input$de_lfc_cutoff %||% 0
    is_deg <- deg_mask(res$df, alpha = alpha, lfc_cutoff = lfc_cut)
    sig_n <- sum(is_deg, na.rm = TRUE)
    padj_fraction <- mean(is_deg, na.rm = TRUE)
    shiny::tagList(
      shiny::tags$p(sprintf("Contrast: %s", res$label)),
      shiny::tags$p(sprintf("DEGs under current thresholds (padj < %.3f and |log2FC| >= %.2f): %s", alpha, lfc_cut, sig_n)),
      shiny::tags$p(sprintf("Fraction with raw p < 0.05: %.3f", res$pvalue_fraction)),
      shiny::tags$p(sprintf("Fraction meeting current DEG thresholds: %.3f", padj_fraction)),
      shiny::tags$p(sprintf("Coefficients: %s", paste(res$names, collapse = ", ")))
    )
  })
  output$transform_summary <- shiny::renderUI({
    res <- de()
    shiny::tagList(
      shiny::tags$p("These plots use VST-transformed expression values. They are diagnostic views only; the DESeq2 model itself was fitted on raw counts."),
      shiny::tags$p("Use this step to check whether major group structure remains interpretable after transformation and whether one or two samples dominate the spread. This is the last major quality gate before you start telling a biological story from the DE result."),
      shiny::tags$p("If the VST PCA looks almost identical to the earlier QC PCA, that is not automatically a problem. It often means the same biological or technical structure is genuinely dominating both views.")
    )
  })
  output$vst_pca_plot <- plotly::renderPlotly({
    cur <- ds()
    res <- de()
    scores <- as.data.frame(res$vst_pca$x[, 1:2, drop = FALSE])
    scores$run_id <- rownames(scores)
    color_var <- if (!is.null(input$contrast_var) && nzchar(input$contrast_var) && input$contrast_var %in% names(cur$meta)) input$contrast_var else input$qc_group
    scores$group <- cur$meta[[color_var]][match(scores$run_id, cur$meta$run_id)]
    vars <- (res$vst_pca$sdev ^ 2) / sum(res$vst_pca$sdev ^ 2)
    plotly::layout(
      plotly::plot_ly(scores, x = ~PC1, y = ~PC2, color = ~group, type = "scatter", mode = "markers+text", text = ~run_id, textposition = "top center"),
      title = "VST PCA",
      xaxis = list(title = sprintf("PC1 (%.1f%%)", vars[1L] * 100)),
      yaxis = list(title = sprintf("PC2 (%.1f%%)", vars[2L] * 100))
    )
  })
  output$vst_box_plot <- plotly::renderPlotly({
    res <- de()
    bx <- res$transformed_box
    names(bx) <- c("value", "run_id")
    plotly::layout(
      plotly::plot_ly(bx, x = ~run_id, y = ~value, type = "box", boxpoints = FALSE),
      title = "VST distribution by sample",
      xaxis = list(title = ""),
      yaxis = list(title = "VST expression")
    )
  })
  output$mean_sd_plot <- plotly::renderPlotly({
    msd <- de()$mean_sd
    plotly::layout(
      plotly::plot_ly(
        msd,
        x = ~vst_mean,
        y = ~vst_sd,
        type = "scatter",
        mode = "markers",
        text = ~paste("Gene:", gene_label, "<br>Full gene name:", full_gene_name, "<br>Ensembl:", gene_id, "<br>VST mean:", signif(vst_mean, 3), "<br>VST SD:", signif(vst_sd, 3)),
        hoverinfo = "text"
      ),
      title = "Mean-SD trend after VST",
      xaxis = list(title = "Mean VST expression"),
      yaxis = list(title = "SD across samples")
    )
  })
  output$diag_summary <- shiny::renderUI({
    res <- de()
    notes <- c(
      "A left-enriched raw p-value histogram can reflect real global biology, but it can also point to confounding, sample duplication, or an under-modeled design.",
      "If adjusted p-values remain extremely concentrated near zero, compare this design against a more adjusted saved formula, or temporarily exclude suspicious samples and rerun.",
      "Check the VST PCA, size factors, and dispersion trend together before interpreting a very large number of DE genes as purely biological."
    )
    shiny::tagList(
      shiny::tags$p(sprintf("Current result has %.1f%% of tested genes below raw p < 0.05 and %.1f%% below adjusted p < 0.05.", 100 * res$pvalue_fraction, 100 * res$padj_fraction)),
      shiny::tags$p("Read these diagnostics together rather than in isolation. The goal is not simply to reduce the number of significant genes, but to decide whether the current model is calibrated, interpretable, and biologically believable enough for downstream pathway, module, and causal analyses."),
      shiny::tags$ul(lapply(notes, shiny::tags$li))
    )
  })
  output$dispersion_plot <- plotly::renderPlotly({
    disp <- de()$dispersions
    keep <- !is.na(disp$baseMean) & !is.na(disp$dispGeneEst)
    disp <- disp[keep, , drop = FALSE]
    plotly::layout(
      plotly::plot_ly(
        disp,
        x = ~log10(baseMean + 1),
        y = ~log10(dispGeneEst + 1e-8),
        type = "scatter",
        mode = "markers",
        text = ~paste("Gene:", gene_label, "<br>Full gene name:", full_gene_name, "<br>Ensembl:", gene_id, "<br>baseMean:", signif(baseMean, 3), "<br>Gene-wise dispersion:", signif(dispGeneEst, 3), if ("dispFit" %in% names(disp)) paste0("<br>Fitted dispersion: ", signif(dispFit, 3)) else ""),
        hoverinfo = "text"
      ),
      title = "Dispersion trend",
      xaxis = list(title = "log10(baseMean + 1)"),
      yaxis = list(title = "log10(gene-wise dispersion)")
    )
  })
  output$dispersion_help <- shiny::renderUI({
    shiny::tagList(
      shiny::tags$p("How to read this plot: each point is one gene. Genes with lower average counts usually have higher uncertainty, so they often sit higher on the left."),
      shiny::tags$p("A healthy trend usually shows dispersion gradually decreasing as mean expression increases, without a large cloud of extreme outliers dominating the plot."),
      shiny::tags$p("If many genes stay very dispersed even at high expression, or if the pattern looks unusually ragged, revisit sample QC, replicate balance, and whether the design formula is missing an important source of variation."),
      shiny::tags$p("Why this matters for the whole tool: if this plot looks problematic, enrichment, WGCNA, and SEMgraph will all be harder to trust because they are being built on top of a questionable DESeq2 result.")
    )
  })
  output$size_factor_plot <- plotly::renderPlotly({
    cur <- ds()
    res <- de()
    color_var <- if (!is.null(input$contrast_var) && nzchar(input$contrast_var) && input$contrast_var %in% names(cur$meta)) input$contrast_var else input$qc_group
    pd <- data.frame(
      run_id = names(res$size_factors),
      size_factor = as.numeric(res$size_factors),
      group = cur$meta[[color_var]][match(names(res$size_factors), cur$meta$run_id)],
      stringsAsFactors = FALSE
    )
    plotly::layout(
      plotly::plot_ly(
        pd,
        x = ~run_id,
        y = ~size_factor,
        color = ~group,
        type = "bar",
        text = ~paste("Sample:", run_id, "<br>Group:", group, "<br>Size factor:", signif(size_factor, 3)),
        hoverinfo = "text"
      ),
      title = "DESeq2 size factors",
      xaxis = list(title = ""),
      yaxis = list(title = "Size factor")
    )
  })
  output$pvalue_hist <- plotly::renderPlotly({
    df <- de()$df
    pv <- df$pvalue[!is.na(df$pvalue)]
    plotly::layout(
      plotly::plot_ly(x = pv, type = "histogram", nbinsx = 30),
      title = "Raw p-value histogram",
      xaxis = list(title = "p-value"),
      yaxis = list(title = "Gene count")
    )
  })
  output$pvalue_ma_plot <- plotly::renderPlotly({
    df <- de()$df
    df$neg_log10_p <- -log10(pmax(df$pvalue, .Machine$double.xmin))
    plotly::layout(
      plotly::plot_ly(
        df,
        x = ~log10(baseMean + 1),
        y = ~neg_log10_p,
        color = ~sig,
        type = "scatter",
        mode = "markers",
        text = ~paste("Gene:", gene_label, "<br>Full gene name:", full_gene_name, "<br>Ensembl:", gene_id, "<br>baseMean:", signif(baseMean, 3), "<br>p-value:", signif(pvalue, 3), "<br>padj:", signif(padj, 3)),
        hoverinfo = "text"
      ),
      title = "Significance versus mean expression",
      xaxis = list(title = "log10(baseMean + 1)"),
      yaxis = list(title = "-log10 raw p-value")
    )
  })
  output$improvement_suggestions <- shiny::renderUI({
    cur <- ds()
    res <- de()
    suggestions <- character()

    if (res$padj_fraction > 0.2) {
      suggestions <- c(suggestions, "More than 20% of genes are significant after adjustment. Compare this design against a model that adds a technical term such as batch, donor, lane, or run metadata if available.")
    }
    if (max(res$size_factors, na.rm = TRUE) / max(min(res$size_factors, na.rm = TRUE), 1e-8) > 3) {
      suggestions <- c(suggestions, "Size factors vary strongly across samples. Check whether a few libraries dominate sequencing depth and whether excluding clear outliers changes the result.")
    }
    if (length(input$design_terms %||% character()) <= 1L) {
      suggestions <- c(suggestions, "The current design is simple. Save a second design that adjusts for another plausible source of variation and compare the p-value histograms and number of significant genes.")
    }
    if (length(input$de_samples %||% cur$meta$run_id) == nrow(cur$meta)) {
      suggestions <- c(suggestions, "Try a sensitivity rerun excluding any sample that looks separated in PCA or the distance heatmap. A better-calibrated model should show more stable clustering and a less extreme p-value concentration if an outlier was driving the signal.")
    }
    suggestions <- c(
      suggestions,
      "Improvement is supported when the saved alternative design still captures the main biological contrast, but diagnostic plots look more balanced and the results are easier to explain biologically."
    )

    shiny::tagList(
      shiny::tags$h4("What to try next"),
      shiny::tags$ul(lapply(unique(suggestions), shiny::tags$li))
    )
  })
  output$padj_hist <- plotly::renderPlotly({
    df <- de()$df
    pv <- df$padj[!is.na(df$padj)]
    plotly::layout(
      plotly::plot_ly(x = pv, type = "histogram", nbinsx = 30),
      title = "Adjusted p-value histogram",
      xaxis = list(title = "adjusted p-value"),
      yaxis = list(title = "Gene count")
    )
  })
  output$ma_plot <- plotly::renderPlotly({
    alpha <- input$de_alpha %||% 0.05
    lfc_cut <- input$de_lfc_cutoff %||% 0
    df <- de()$df
    df$state <- ifelse(deg_mask(df, alpha = alpha, lfc_cutoff = lfc_cut), sprintf("DEG: padj < %.3f and |log2FC| >= %.2f", alpha, lfc_cut), "Outside current DEG cutoff")
    df$m <- log10(df$baseMean + 1)
    make_downloadable_plot(plotly::layout(plotly::plot_ly(df, x = ~m, y = ~log2FoldChange, color = ~state, type = "scatter", mode = "markers", text = ~paste("Gene:", gene_label, "<br>Full gene name:", full_gene_name, "<br>Ensembl:", gene_id, "<br>baseMean:", signif(baseMean, 3), "<br>log2FC:", signif(log2FoldChange, 3), "<br>padj:", signif(padj, 3)), hoverinfo = "text"), title = "MA-style plot", xaxis = list(title = "log10(baseMean + 1)"), yaxis = list(title = "log2 fold change")), "deseq2_ma_plot")
  })
  output$volcano_plot <- plotly::renderPlotly({
    alpha <- input$de_alpha %||% 0.05
    lfc_cut <- input$de_lfc_cutoff %||% 0
    df <- de()$df
    df$state <- ifelse(deg_mask(df, alpha = alpha, lfc_cutoff = lfc_cut), sprintf("DEG: padj < %.3f and |log2FC| >= %.2f", alpha, lfc_cut), "Outside current DEG cutoff")
    df$nlp <- -log10(pmax(df$padj, .Machine$double.xmin))
    make_downloadable_plot(plotly::layout(plotly::plot_ly(df, x = ~log2FoldChange, y = ~nlp, color = ~state, type = "scatter", mode = "markers", text = ~paste("Gene:", gene_label, "<br>Full gene name:", full_gene_name, "<br>Ensembl:", gene_id, "<br>log2FC:", signif(log2FoldChange, 3), "<br>padj:", signif(padj, 3)), hoverinfo = "text"), title = "Volcano plot", xaxis = list(title = "log2 fold change"), yaxis = list(title = "-log10 adjusted p-value")), "deseq2_volcano_plot")
  })
  output$heat_plot <- plotly::renderPlotly({
    cur <- ds()
    res <- de()
    alpha <- input$de_alpha %||% 0.05
    lfc_cut <- input$de_lfc_cutoff %||% 0
    top_deg <- head(res$df$gene_id[deg_mask(res$df, alpha = alpha, lfc_cutoff = lfc_cut)], 30L)
    shiny::req(length(top_deg) > 0L)
    hm <- res$vst[top_deg, res$samples, drop = FALSE]
    col_ord <- if (ncol(hm) > 1L) stats::hclust(stats::dist(t(hm)))$order else 1L
    row_ord <- if (nrow(hm) > 1L) stats::hclust(stats::dist(hm))$order else 1L
    hm <- hm[row_ord, col_ord, drop = FALSE]
    sample_ids <- colnames(hm)
    color_var <- if (!is.null(input$contrast_var) && nzchar(input$contrast_var) && input$contrast_var %in% names(cur$meta)) input$contrast_var else input$qc_group
    sample_groups <- cur$meta[[color_var]]
    names(sample_groups) <- cur$meta$run_id
    sample_groups <- sample_groups[sample_ids]
    sample_colors <- group_colors(sample_groups)
    sample_ticktext <- tick_text(sample_ids, sample_groups, sample_colors)
    gene_labels <- cur$gene_map$gene_label[match(rownames(hm), cur$gene_map$gene_id)]
    ann_plot <- build_metadata_annotation_plot(cur$meta, sample_ids, title = "Sample metadata")
    main_plot <- plotly::plot_ly(
      x = sample_ids,
      y = rev(gene_labels),
      z = hm[nrow(hm):1, , drop = FALSE],
      type = "heatmap",
      colorscale = "RdBu"
    ) |>
      plotly::layout(
        title = "Top DE genes heatmap (VST)",
        xaxis = list(title = "Samples", tickmode = "array", tickvals = sample_ids, ticktext = sample_ticktext),
        yaxis = list(title = "Genes")
      )
    if (is.null(ann_plot)) {
      main_plot
    } else {
      plotly::subplot(ann_plot, main_plot, nrows = 2, shareX = TRUE, heights = c(0.18, 0.82), titleY = TRUE)
    }
  })
  output$heat_legend <- shiny::renderUI({
    cur <- ds()
    res <- de()
    alpha <- input$de_alpha %||% 0.05
    lfc_cut <- input$de_lfc_cutoff %||% 0
    shiny::req(any(deg_mask(res$df, alpha = alpha, lfc_cutoff = lfc_cut), na.rm = TRUE))
    color_var <- if (!is.null(input$contrast_var) && nzchar(input$contrast_var) && input$contrast_var %in% names(cur$meta)) input$contrast_var else input$qc_group
    cols <- group_colors(cur$meta[[color_var]])
    shiny::tagList(
      shiny::tags$div(class = "muted", style = "font-size:12px;margin-bottom:6px;", paste("Colors:", color_var)),
      lapply(names(cols), function(nm) {
        shiny::tags$div(
          class = "legend-item",
          style = "margin-bottom:4px;",
          shiny::tags$span(class = "sw", style = paste("background:", cols[[nm]], ";")),
          shiny::tags$span(style = "font-size:12px;", nm)
        )
      })
    )
  })
  output$de_table <- DT::renderDataTable({
    alpha <- input$de_alpha %||% 0.05
    lfc_cut <- input$de_lfc_cutoff %||% 0
    df <- de()$df[, c("gene_label", "gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), drop = FALSE]
    df$current_deg <- deg_mask(de()$df, alpha = alpha, lfc_cutoff = lfc_cut)
    names(df)[1:2] <- c("gene", "ensembl_id")
    names(df)[ncol(df)] <- "meets_current_deg_cutoff"
    make_downloadable_table(df, page_length = 10)
  }, server = FALSE)

  enrich <- shiny::eventReactive(input$run_enrichment, {
    shiny::withProgress(message = "Running pathway enrichment", detail = "Evaluating selected enrichment methods. Please wait.", value = 0, {
      cur <- ds()
      de_res <- de()
      shiny::incProgress(0.2, detail = "Preparing ranked and significant gene sets")
      run_enrichment(
        de_res = de_res,
        gene_map = cur$gene_map,
        methods = input$enrich_methods %||% character(),
        padj_cutoff = input$enrich_padj,
        lfc_cutoff = input$enrich_lfc,
        min_gs = input$enrich_min_gs,
        max_gs = input$enrich_max_gs,
        p_adjust_method = input$enrich_adjust_method,
        top_n = input$enrich_top_n,
        refresh = isTRUE(input$refresh_enrichment_cache)
      )
    })
  })
  shiny::observeEvent(input$run_enrichment, {
    analysis_log(append_analysis_log(analysis_log(), "run_enrichment", list(
      methods = input$enrich_methods %||% character(),
      padj_cutoff = input$enrich_padj,
      abs_log2fc_cutoff = input$enrich_lfc,
      min_pathway_size = input$enrich_min_gs,
      max_pathway_size = input$enrich_max_gs,
      adjust_method = input$enrich_adjust_method,
      top_n = input$enrich_top_n
    )))
  }, ignoreInit = TRUE)
  output$enrich_help <- shiny::renderUI({
    shiny::tagList(
      shiny::tags$p("How to choose a method: use ORA when you want to focus on the subset of genes that already passed your DE thresholds. Use GSEA when you want pathway-level trends from the full ranked gene list, including coordinated but smaller shifts."),
      shiny::tags$p("A practical novice workflow is to start with GO ORA and GO GSEA together, then add KEGG if you want a second pathway database for consensus."),
      shiny::tags$p("Role in the overall workflow: this step translates DESeq2 output into biological processes. It is often the clearest bridge between a gene list and later network-level interpretation in WGCNA and SEMgraph.")
    )
  })
  output$enrich_summary <- shiny::renderUI({
    res <- enrich()
    shiny::tagList(
      if (nrow(res$combined) == 0L) shiny::tags$div(class = "status-note", "No pathways passed the current enrichment settings. That can happen when the DE signal is weak, the selected thresholds are strict, or the tested pathway databases do not capture the biology well."),
      shiny::tags$p(sprintf("Significant genes entering ORA: %s", res$sig_gene_count)),
      shiny::tags$p(sprintf("Ranked genes entering GSEA: %s", res$ranked_gene_count)),
      shiny::tags$p(sprintf("Enrichment rows returned across methods: %s", nrow(res$combined))),
      shiny::tags$p("Consensus pathways are those that recur across more than one enrichment method in this run."),
      shiny::tags$p("Interpretation tip: use this summary to decide which biological themes are strong enough to carry forward into module interpretation, combined programs, and causal prioritization.")
    )
  })
  output$enrich_mechanism <- shiny::renderUI({
    methods <- input$enrich_methods %||% character()
    shiny::tagList(
      shiny::tags$h4("How selected methods are used"),
      shiny::tags$p(sprintf("You selected %s method(s): %s.", length(methods), if (length(methods)) paste(methods, collapse = ", ") else "none")),
      shiny::tags$ul(
        shiny::tags$li("Each selected method runs separately on the same DE result."),
        shiny::tags$li("The `By Method` tab shows one method at a time, with its own pathways and significance values."),
        shiny::tags$li("The `Consensus` tab does not recompute enrichment. It summarizes overlap across the already-run methods."),
        shiny::tags$li("A pathway appearing in several methods can be more robust, but the methods are still conceptually different and should not be treated as identical tests.")
      )
    )
  })
  output$consensus_table <- DT::renderDataTable({
    res <- enrich()
    make_downloadable_table(res$consensus, page_length = 8)
  }, server = FALSE)
  output$rank_hist_plot <- plotly::renderPlotly({
    ranked <- enrich()$ranked
    plotly::layout(
      plotly::plot_ly(x = ranked, type = "histogram", nbinsx = 40),
      title = "Ranked statistic distribution for GSEA",
      xaxis = list(title = "DESeq2 Wald statistic"),
      yaxis = list(title = "Gene count")
    )
  })
  output$method_balance_plot <- plotly::renderPlotly({
    res <- enrich()$combined
    shiny::req(nrow(res) > 0L)
    counts <- as.data.frame(table(res$method), stringsAsFactors = FALSE)
    names(counts) <- c("method", "pathways")
    plotly::layout(
      plotly::plot_ly(counts, x = ~method, y = ~pathways, type = "bar", text = ~pathways, textposition = "auto"),
      title = "Returned pathways by method",
      xaxis = list(title = ""),
      yaxis = list(title = "Pathway count")
    )
  })
  shiny::observe({
    method_names <- names(enrich()$method_tables)
    selected <- if (!is.null(input$enrich_view_method) && input$enrich_view_method %in% method_names) input$enrich_view_method else if (length(method_names)) method_names[1L] else character()
    shiny::updateSelectInput(session, "enrich_view_method", choices = method_names, selected = selected)
  })
  output$enrich_method_plot <- plotly::renderPlotly({
    method_tables <- enrich()$method_tables
    shiny::req(length(method_tables) > 0L)
    shiny::req(input$enrich_view_method %in% names(method_tables))
    top <- head(method_tables[[input$enrich_view_method]], input$enrich_top_n)
    top$pathway <- factor(top$pathway, levels = rev(unique(top$pathway)))
    top$score <- -log10(pmax(top$p.adjust, .Machine$double.xmin))
    plotly::layout(
      plotly::plot_ly(
        top,
        x = ~score,
        y = ~pathway,
        color = ~method,
        size = ~gene_ratio_value,
        type = "scatter",
        mode = "markers",
        text = ~paste("Pathway:", pathway, "<br>Method:", method, "<br>Adjusted p:", signif(p.adjust, 3), if ("GeneRatio" %in% names(top)) paste0("<br>Gene ratio: ", GeneRatio) else "", if ("NES" %in% names(top)) paste0("<br>NES: ", signif(NES, 3)) else ""),
        hoverinfo = "text"
      ),
      title = paste("Top pathways for", input$enrich_view_method),
      xaxis = list(title = "-log10 adjusted p-value"),
      yaxis = list(title = "")
    )
  })
  output$enrich_consensus_plot <- plotly::renderPlotly({
    res <- enrich()$consensus
    shiny::req(nrow(res) > 0L)
    top <- head(res, input$enrich_top_n)
    top$pathway <- factor(top$pathway, levels = rev(unique(top$pathway)))
    top$consensus_hits <- as.numeric(top$consensus_hits)
    top$p.adjust <- suppressWarnings(as.numeric(top$p.adjust))
    top$score <- -log10(pmax(top$p.adjust, .Machine$double.xmin))
    point_sizes <- pmax(10, top$score * 6)
    plotly::layout(
      plotly::plot_ly(
        top,
        x = ~consensus_hits,
        y = ~pathway,
        type = "scatter",
        mode = "markers",
        marker = list(
          size = point_sizes,
          color = "#2a6fdb",
          line = list(color = "#163a70", width = 1)
        ),
        text = ~paste("Pathway:", pathway, "<br>Methods supporting it:", consensus_hits, "<br>Best adjusted p-value:", signif(p.adjust, 3)),
        hoverinfo = "text"
      ),
      title = "Consensus pathways across selected methods",
      xaxis = list(title = "Number of supporting methods"),
      yaxis = list(title = "")
    )
  })
  output$enrich_method_table <- DT::renderDataTable({
    method_tables <- enrich()$method_tables
    shiny::req(length(method_tables) > 0L)
    shiny::req(input$enrich_view_method %in% names(method_tables))
    res <- method_tables[[input$enrich_view_method]]
    keep <- intersect(c("pathway", "Description", "GeneRatio", "NES", "pvalue", "p.adjust", "qvalue", "Count", "core_enrichment", "geneID"), names(res))
    make_downloadable_table(res[, keep, drop = FALSE], page_length = 12)
  }, server = FALSE)
  output$enrich_table <- DT::renderDataTable({
    res <- enrich()$combined
    shiny::req(nrow(res) > 0L)
    keep <- intersect(c("method", "pathway", "Description", "GeneRatio", "NES", "pvalue", "p.adjust", "qvalue", "Count", "core_enrichment", "geneID", "consensus_hits"), names(res))
    make_downloadable_table(res[, keep, drop = FALSE], page_length = 12)
  }, server = FALSE)

  shiny::observe({
    cur <- ds()
    trait_choices <- setdiff(names(cur$meta), "run_id")
    shiny::updateCheckboxGroupInput(session, "wgcna_traits", choices = trait_choices, selected = if ("treatment" %in% trait_choices) "treatment" else trait_choices)
  })
  wgcna_res <- shiny::eventReactive(input$run_wgcna, {
    shiny::withProgress(message = "Running WGCNA", detail = "Building network and detecting modules. Please wait.", value = 0, {
      cur <- ds()
      shiny::incProgress(0.15, detail = "Selecting genes and estimating soft-threshold power")
      run_wgcna_analysis(
        expr_matrix = cur$filtered,
        metadata = cur$meta,
        trait_cols = input$wgcna_traits %||% character(),
        power_mode = input$wgcna_power_mode,
        manual_power = input$wgcna_manual_power,
        min_module_size = input$wgcna_min_module_size,
        deep_split = input$wgcna_deep_split,
        merge_cut_height = input$wgcna_merge_cut_height,
        top_var_genes = input$wgcna_top_var_genes
      )
    })
  })
  shiny::observeEvent(input$run_wgcna, {
    analysis_log(append_analysis_log(analysis_log(), "run_wgcna", list(
      power_mode = input$wgcna_power_mode,
      manual_power = input$wgcna_manual_power,
      top_variable_genes = input$wgcna_top_var_genes,
      min_module_size = input$wgcna_min_module_size,
      deep_split = input$wgcna_deep_split,
      merge_cut_height = input$wgcna_merge_cut_height,
      traits = input$wgcna_traits %||% character()
    )))
  }, ignoreInit = TRUE)
  output$wgcna_help <- shiny::renderUI({
    shiny::tagList(
      shiny::tags$p("WGCNA groups genes into modules based on similar expression patterns across samples. A module can then be compared against traits such as treatment or batch."),
      shiny::tags$p("Beginner workflow: start with automatic power, leave the merge threshold near 0.2, inspect the module tree, then adjust only if modules look obviously over-split or over-merged."),
      shiny::tags$p("Role in the overall workflow: DESeq2 asks what changed between groups; WGCNA asks which genes move together across samples. Use WGCNA to organize the DE signal into coordinated systems rather than interpreting genes one by one.")
    )
  })
  output$wgcna_summary <- shiny::renderUI({
    res <- wgcna_res()
    shiny::tagList(
      if (length(unique(res$module_colors)) <= 2L) shiny::tags$div(class = "status-note", "Very few modules were detected. This may mean the data are dominated by one broad pattern, the sample size is limited, or the split/merge settings are too conservative for the current dataset."),
      shiny::tags$p(sprintf("Automatic suggested power: %s", res$auto_power)),
      shiny::tags$p(sprintf("Power used in this run: %s", res$selected_power)),
      shiny::tags$p(sprintf("Detected modules (including grey/unassigned): %s", length(unique(res$module_colors)))),
      shiny::tags$p(sprintf("Genes used for network construction: %s", ncol(res$dat_expr))),
      shiny::tags$p("Beginner tip: a good power is usually one that gives a high scale-free fit without making connectivity collapse too sharply. The automatic suggestion is a reasonable starting point, not a final truth."),
      shiny::tags$p("Interpretation tip: do not judge WGCNA only by the number of modules. A useful run is one that produces modules you can connect back to traits, DE patterns, and biological enrichment.")
    )
  })
  output$wgcna_power_plot <- plotly::renderPlotly({
    fit <- wgcna_res()$fit_tbl
    plotly::layout(
      plotly::plot_ly(
        fit,
        x = ~Power,
        y = ~SFT.R.sq,
        type = "scatter",
        mode = "lines+markers",
        text = ~paste("Power:", Power, "<br>Scale-free fit:", signif(SFT.R.sq, 3), "<br>Mean connectivity:", signif(mean.k., 3)),
        hoverinfo = "text"
      ),
      title = "Soft-threshold selection",
      xaxis = list(title = "Power"),
      yaxis = list(title = "Scale-free topology fit")
    )
  })
  output$wgcna_tree_help <- shiny::renderUI({
    shiny::tagList(
      shiny::tags$p("How to use this tree: branches group genes with similar expression. The color bars underneath show module assignments before and after merging."),
      shiny::tags$p("If two neighboring color blocks are well separated in the tree, keep them separate. If they are highly interwoven and the eigengene heatmap also shows them as very similar, increasing the merge threshold can be reasonable."),
      shiny::tags$p("If a large module looks like it contains several clearly distinct branches, increase split sensitivity or lower the merge threshold and rerun.")
    )
  })
  output$wgcna_dendrogram_plot <- shiny::renderPlot({
    res <- wgcna_res()
    WGCNA::plotDendroAndColors(
      res$dendrogram,
      cbind(Unmerged = res$unmerged_colors, Merged = res$module_colors),
      groupLabels = c("Before merge", "After merge"),
      dendroLabels = FALSE,
      hang = 0.03,
      addGuide = TRUE,
      guideHang = 0.05
    )
  })
  shiny::observe({
    res <- wgcna_res()
    if (!is.null(res$module_trait_cor)) {
      choices <- colnames(res$module_trait_cor)
      selected <- if (!is.null(input$wgcna_trait_view) && input$wgcna_trait_view %in% choices) input$wgcna_trait_view else choices[1L]
      shiny::updateSelectInput(session, "wgcna_trait_view", choices = choices, selected = selected)
    } else {
      shiny::updateSelectInput(session, "wgcna_trait_view", choices = character(), selected = character())
    }
  })
  output$wgcna_trait_heatmap <- plotly::renderPlotly({
    res <- wgcna_res()
    shiny::req(!is.null(res$module_trait_cor))
    txt <- matrix("", nrow = nrow(res$module_trait_cor), ncol = ncol(res$module_trait_cor))
    for (i in seq_len(nrow(res$module_trait_cor))) {
      for (j in seq_len(ncol(res$module_trait_cor))) {
        txt[i, j] <- paste0(signif(res$module_trait_cor[i, j], 2), "<br>p=", signif(res$module_trait_p[i, j], 2))
      }
    }
    plotly::layout(
      plotly::plot_ly(
        x = colnames(res$module_trait_cor),
        y = rownames(res$module_trait_cor),
        z = res$module_trait_cor,
        text = txt,
        type = "heatmap",
        colorscale = "RdBu",
        zmin = -1,
        zmax = 1,
        hoverinfo = "text"
      ),
      title = "Module-trait correlations",
      xaxis = list(title = "Traits / contrasts"),
      yaxis = list(title = "Modules")
    )
  })
  output$wgcna_trait_focus_plot <- plotly::renderPlotly({
    res <- wgcna_res()
    shiny::req(!is.null(res$module_trait_cor))
    shiny::req(input$wgcna_trait_view %in% colnames(res$module_trait_cor))
    vals <- data.frame(
      module = rownames(res$module_trait_cor),
      correlation = res$module_trait_cor[, input$wgcna_trait_view],
      pvalue = res$module_trait_p[, input$wgcna_trait_view],
      stringsAsFactors = FALSE
    )
    vals <- vals[order(abs(vals$correlation), decreasing = TRUE), , drop = FALSE]
    plotly::layout(
      plotly::plot_ly(
        vals,
        x = ~reorder(module, correlation),
        y = ~correlation,
        type = "bar",
        text = ~paste("Module:", module, "<br>Correlation:", signif(correlation, 3), "<br>p-value:", signif(pvalue, 3)),
        hoverinfo = "text"
      ),
      title = paste("Module association with", input$wgcna_trait_view),
      xaxis = list(title = "", tickangle = -45),
      yaxis = list(title = "Correlation")
    )
  })
  output$wgcna_module_size_plot <- plotly::renderPlotly({
    sizes <- as.data.frame(wgcna_res()$module_sizes, stringsAsFactors = FALSE)
    names(sizes) <- c("module", "genes")
    ann_tbl <- wgcna_annotations()
    sizes$annotation <- ann_tbl$annotation[match(sizes$module, ann_tbl$module)]
    sizes$display_module <- ifelse(is.na(sizes$annotation) | !nzchar(sizes$annotation), sizes$module, paste0(sizes$module, " - ", sizes$annotation))
    plotly::layout(
      plotly::plot_ly(sizes, x = ~display_module, y = ~genes, type = "bar", text = ~genes, textposition = "auto", hovertext = ~paste("Module:", module, ifelse(is.na(annotation) | !nzchar(annotation), "", paste0("<br>Annotation: ", annotation)), "<br>Genes:", genes), hoverinfo = "text"),
      title = "Module sizes",
      xaxis = list(title = "Module"),
      yaxis = list(title = "Gene count")
    )
  })
  output$wgcna_module_help <- shiny::renderUI({
    shiny::req(nzchar(input$wgcna_module_view))
    ann_tbl <- wgcna_annotations()
    row <- ann_tbl[ann_tbl$module == input$wgcna_module_view, , drop = FALSE]
    shiny::tagList(
      shiny::tags$h4(paste("How to interpret module", input$wgcna_module_view)),
      if (nrow(row) && nzchar(row$annotation[1L])) shiny::tags$p(sprintf("Current biological label: %s", row$annotation[1L])),
      if (nrow(row) && nzchar(row$notes[1L])) shiny::tags$p(sprintf("Notes: %s", row$notes[1L])),
      shiny::tags$ul(
        shiny::tags$li("Connectivity distribution: shows whether the module is broadly cohesive or driven by only a few highly connected genes."),
        shiny::tags$li("Eigengene membership: shows how strongly each gene follows the overall module pattern. Higher values usually indicate more representative module genes."),
        shiny::tags$li("Network graph: a simplified view of the strongest links only. Dense clusters suggest a tight module; thin fragmented structure suggests a weaker or more heterogeneous module."),
        shiny::tags$li("Hub genes: genes with especially high within-module connectivity. They are useful candidates for follow-up, but should still be checked biologically rather than assumed to be causal.")
      )
    )
  })
  shiny::observe({
    gm <- wgcna_res()$gene_modules
    modules <- sort(unique(gm$module))
    modules <- modules[modules != "grey"]
    selected <- if (!is.null(input$wgcna_module_view) && input$wgcna_module_view %in% modules) input$wgcna_module_view else if (length(modules)) modules[1L] else character()
    shiny::updateSelectInput(session, "wgcna_module_view", choices = modules, selected = selected)
  })
  shiny::observe({
    shiny::req(nzchar(input$wgcna_module_view))
    ann_tbl <- wgcna_annotations()
    row <- ann_tbl[ann_tbl$module == input$wgcna_module_view, , drop = FALSE]
    shiny::updateTextInput(session, "wgcna_module_annotation", value = if (nrow(row)) row$annotation[1L] else "")
    shiny::updateTextAreaInput(session, "wgcna_module_annotation_notes", value = if (nrow(row)) row$notes[1L] else "")
  })
  shiny::observeEvent(input$save_wgcna_annotation, {
    shiny::req(nzchar(input$wgcna_module_view))
    tbl <- wgcna_annotations()
    tbl <- tbl[tbl$module != input$wgcna_module_view, , drop = FALSE]
    tbl <- rbind(
      tbl,
      data.frame(
        module = input$wgcna_module_view,
        annotation = trimws(input$wgcna_module_annotation),
        notes = trimws(input$wgcna_module_annotation_notes),
        stringsAsFactors = FALSE
      )
    )
    wgcna_annotations(tbl)
  })
  output$wgcna_module_network_plot <- plotly::renderPlotly({
    cur <- ds()
    res <- wgcna_res()
    shiny::req(nzchar(input$wgcna_module_view))
    ann_tbl <- wgcna_annotations()
    row_ann <- ann_tbl[ann_tbl$module == input$wgcna_module_view, , drop = FALSE]
    module_title <- if (nrow(row_ann) && nzchar(row_ann$annotation[1L])) paste(input$wgcna_module_view, "-", row_ann$annotation[1L]) else input$wgcna_module_view
    gm <- res$gene_modules
    g_mod <- build_wgcna_module_graph(res, input$wgcna_module_view, edge_quantile = input$wgcna_edge_quantile)
    shiny::req(!is.null(g_mod))
    mod_genes <- igraph::V(g_mod)$name
    adj <- igraph::as_adjacency_matrix(g_mod, attr = "weight", sparse = FALSE)
    diag(adj) <- 0
    edge_tbl <- igraph::as_data_frame(g_mod, what = "edges")

    dist_mat <- 1 - adj
    diag(dist_mat) <- 0
    coords <- cmdscale(stats::as.dist(dist_mat), k = 2)
    coords <- as.data.frame(coords)
    coords$gene_id <- rownames(coords)
    coords <- attach_gene_labels(coords, cur$gene_map)
    coords$kWithin <- gm$kWithin[match(coords$gene_id, gm$gene_id)]

    edge_df <- if (nrow(edge_tbl) > 0L) {
      data.frame(
        x = coords$V1[match(edge_tbl$from, coords$gene_id)],
        y = coords$V2[match(edge_tbl$from, coords$gene_id)],
        xend = coords$V1[match(edge_tbl$to, coords$gene_id)],
        yend = coords$V2[match(edge_tbl$to, coords$gene_id)],
        weight = edge_tbl$weight,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame()
    }

    p <- plotly::plot_ly()
    if (nrow(edge_df) > 0L) {
      for (i in seq_len(nrow(edge_df))) {
        p <- plotly::add_segments(
          p,
          x = edge_df$x[i], y = edge_df$y[i],
          xend = edge_df$xend[i], yend = edge_df$yend[i],
          inherit = FALSE,
          line = list(color = "rgba(120,120,120,0.35)", width = 1)
        )
      }
    }
    p <- plotly::add_markers(
      p,
      data = coords,
      x = ~V1,
      y = ~V2,
      color = I(input$wgcna_module_view),
      text = ~paste("Gene:", gene_label, "<br>Full gene name:", full_gene_name, "<br>Ensembl:", gene_id, "<br>Within-module connectivity:", signif(kWithin, 3)),
      hoverinfo = "text",
      marker = list(
        size = {
          kw <- coords$kWithin
          kw[!is.finite(kw)] <- 0
          rng <- range(kw, na.rm = TRUE)
          if (!all(is.finite(rng)) || diff(rng) == 0) {
            rep(14, length(kw))
          } else {
            10 + 18 * (kw - rng[1]) / diff(rng)
          }
        },
        color = input$wgcna_module_view,
        line = list(color = "#1f1f1f", width = 0.5)
      )
    )
    plotly::layout(
      p,
      title = paste("Module network:", module_title),
      xaxis = list(title = "", showticklabels = FALSE, zeroline = FALSE),
      yaxis = list(title = "", showticklabels = FALSE, zeroline = FALSE)
    )
  })
  output$download_wgcna_module_graphml <- shiny::downloadHandler(
    filename = function() sprintf("wgcna_module_%s.graphml", input$wgcna_module_view %||% "network"),
    content = function(file) {
      g_mod <- build_wgcna_module_graph(wgcna_res(), input$wgcna_module_view, edge_quantile = input$wgcna_edge_quantile)
      if (is.null(g_mod)) stop("No module network available for export.")
      igraph::write_graph(g_mod, file, format = "graphml")
    }
  )
  output$wgcna_connectivity_plot <- plotly::renderPlotly({
    res <- wgcna_res()
    shiny::req(nzchar(input$wgcna_module_view))
    gm <- res$gene_modules
    mod <- gm[gm$module == input$wgcna_module_view, , drop = FALSE]
    plotly::layout(
      plotly::plot_ly(
        x = mod$kWithin,
        type = "histogram",
        nbinsx = 25
      ),
      title = paste("Connectivity distribution:", input$wgcna_module_view),
      xaxis = list(title = "Within-module connectivity"),
      yaxis = list(title = "Gene count")
    )
  })
  output$wgcna_membership_plot <- plotly::renderPlotly({
    cur <- ds()
    res <- wgcna_res()
    shiny::req(nzchar(input$wgcna_module_view))
    gm <- res$gene_modules
    mod <- gm[gm$module == input$wgcna_module_view, , drop = FALSE]
    shiny::req(nrow(mod) > 0L)
    me_col <- paste0("ME", input$wgcna_module_view)
    shiny::req(me_col %in% colnames(res$module_eigengenes))

    expr_sub <- res$dat_expr[, mod$gene_id, drop = FALSE]
    memberships <- apply(expr_sub, 2L, function(x) suppressWarnings(stats::cor(x, res$module_eigengenes[, me_col], use = "pairwise.complete.obs")))
    plot_df <- data.frame(
      gene_id = names(memberships),
      membership = as.numeric(memberships),
      kWithin = mod$kWithin[match(names(memberships), mod$gene_id)],
      stringsAsFactors = FALSE
    )
    plot_df <- attach_gene_labels(plot_df, cur$gene_map)

    plotly::layout(
      plotly::plot_ly(
        plot_df,
        x = ~membership,
        y = ~kWithin,
        type = "scatter",
        mode = "markers",
        text = ~paste("Gene:", gene_label, "<br>Full gene name:", full_gene_name, "<br>Ensembl:", gene_id, "<br>Module membership:", signif(membership, 3), "<br>Within-module connectivity:", signif(kWithin, 3)),
        hoverinfo = "text"
      ),
      title = paste("Membership vs connectivity:", input$wgcna_module_view),
      xaxis = list(title = "Correlation with module eigengene"),
      yaxis = list(title = "Within-module connectivity")
    )
  })
  output$wgcna_hubgene_table <- DT::renderDataTable({
    cur <- ds()
    res <- wgcna_res()
    shiny::req(nzchar(input$wgcna_module_view))
    gm <- res$gene_modules
    mod <- gm[gm$module == input$wgcna_module_view, , drop = FALSE]
    mod <- mod[order(mod$kWithin, decreasing = TRUE), , drop = FALSE]
    mod <- head(mod, 25L)
    mod <- attach_gene_labels(mod, cur$gene_map)
    keep <- mod[, c("gene_label", "gene_id", "module", "kWithin"), drop = FALSE]
    names(keep) <- c("gene", "ensembl_id", "module", "within_module_connectivity")
    make_downloadable_table(keep, page_length = 10)
  }, server = FALSE)
  output$wgcna_module_table <- DT::renderDataTable({
    cur <- ds()
    gm <- wgcna_res()$gene_modules
    gm <- attach_gene_labels(gm, cur$gene_map)
    ann_tbl <- wgcna_annotations()
    gm$annotation <- ann_tbl$annotation[match(gm$module, ann_tbl$module)]
    keep <- gm[, c("gene_label", "gene_id", "module", "annotation", "kWithin"), drop = FALSE]
    names(keep) <- c("gene", "ensembl_id", "module", "module_annotation", "within_module_connectivity")
    make_downloadable_table(keep, page_length = 12)
  }, server = FALSE)
  selected_module_enrichment <- shiny::reactive({
    cur <- ds()
    res <- wgcna_res()
    shiny::req(nzchar(input$wgcna_module_view))
    module_genes <- res$gene_modules$gene_id[res$gene_modules$module == input$wgcna_module_view]
    run_module_enrichment(
      module_genes = module_genes,
      gene_map = cur$gene_map,
      padj_cutoff = 0.05,
      min_gs = 10,
      max_gs = 500,
      p_adjust_method = "BH",
      top_n = 12,
      refresh = isTRUE(input$refresh_module_cache),
      cache_label = paste0("module_", input$wgcna_module_view)
    )
  })
  output$wgcna_module_downstream_help <- shiny::renderUI({
    shiny::tagList(
      shiny::tags$p("This is the bridge from network structure to biological meaning. It helps answer what the selected module may represent, how it relates to the phenotype, and whether its central genes and pathways tell a coherent story."),
      shiny::tags$ul(
        shiny::tags$li("If a module is strongly associated with a trait and also enriched for a pathway family, that gives you a more interpretable biological story."),
        shiny::tags$li("Hub genes tell you which genes are central inside the module. Enrichment tells you what the module is doing as a group."),
        shiny::tags$li("These results are most useful for colored modules with clear internal structure. Grey genes are usually unassigned background and are less informative.")
      )
    )
  })
  output$wgcna_module_summary <- shiny::renderUI({
    res <- wgcna_res()
    shiny::req(nzchar(input$wgcna_module_view))
    gm <- res$gene_modules
    mod <- gm[gm$module == input$wgcna_module_view, , drop = FALSE]
    enrich_res <- selected_module_enrichment()
    ann_tbl <- wgcna_annotations()
    row_ann <- ann_tbl[ann_tbl$module == input$wgcna_module_view, , drop = FALSE]
    trait_lines <- NULL
    if (!is.null(res$module_trait_cor)) {
      me_name <- paste0("ME", input$wgcna_module_view)
      if (me_name %in% rownames(res$module_trait_cor)) {
        trait_df <- data.frame(
          trait = colnames(res$module_trait_cor),
          correlation = res$module_trait_cor[me_name, ],
          pvalue = res$module_trait_p[me_name, ],
          stringsAsFactors = FALSE
        )
        trait_df <- trait_df[order(abs(trait_df$correlation), decreasing = TRUE), , drop = FALSE]
        trait_df <- head(trait_df, 3L)
        trait_lines <- lapply(seq_len(nrow(trait_df)), function(i) {
          shiny::tags$li(sprintf("%s: correlation %.3f, p = %.3g", trait_df$trait[i], trait_df$correlation[i], trait_df$pvalue[i]))
        })
      }
    }
    shiny::tagList(
      if (!isTRUE(enrich_res$has_signal)) shiny::tags$div(class = "status-note", "This selected module did not yield significant pathway enrichment under the current settings. That may reflect a small or broad module, limited pathway coverage, or simply weaker functional coherence."),
      shiny::tags$p(sprintf("Selected module: %s", input$wgcna_module_view)),
      if (nrow(row_ann) && nzchar(row_ann$annotation[1L])) shiny::tags$p(sprintf("Biological annotation: %s", row_ann$annotation[1L])),
      if (nrow(row_ann) && nzchar(row_ann$notes[1L])) shiny::tags$p(sprintf("Annotation notes: %s", row_ann$notes[1L])),
      shiny::tags$p(sprintf("Genes in module: %s", nrow(mod))),
      shiny::tags$p(sprintf("Genes mapped for enrichment: %s", enrich_res$gene_count)),
      if (!is.null(trait_lines)) shiny::tagList(
        shiny::tags$p("Strongest trait associations for this module:"),
        shiny::tags$ul(trait_lines)
      )
    )
  })
  output$wgcna_module_enrichment_plot <- plotly::renderPlotly({
    enrich_res <- selected_module_enrichment()$combined
    shiny::req(nrow(enrich_res) > 0L)
    top <- head(enrich_res, 12L)
    top$pathway <- factor(top$pathway, levels = rev(unique(top$pathway)))
    top$score <- -log10(pmax(as.numeric(top$p.adjust), .Machine$double.xmin))
    plotly::layout(
      plotly::plot_ly(
        top,
        x = ~score,
        y = ~pathway,
        color = ~method,
        size = ~gene_ratio_value,
        type = "scatter",
        mode = "markers",
        text = ~paste("Pathway:", pathway, "<br>Method:", method, "<br>Adjusted p:", signif(p.adjust, 3), if ("GeneRatio" %in% names(top)) paste0("<br>Gene ratio: ", GeneRatio) else ""),
        hoverinfo = "text"
      ),
      title = paste("Pathway enrichment for module", input$wgcna_module_view),
      xaxis = list(title = "-log10 adjusted p-value"),
      yaxis = list(title = "")
    )
  })
  output$wgcna_module_enrichment_table <- DT::renderDataTable({
    enrich_res <- selected_module_enrichment()$combined
    shiny::req(nrow(enrich_res) > 0L)
    keep <- intersect(c("method", "pathway", "Description", "GeneRatio", "pvalue", "p.adjust", "qvalue", "Count", "geneID", "consensus_hits"), names(enrich_res))
    make_downloadable_table(enrich_res[, keep, drop = FALSE], page_length = 10)
  }, server = FALSE)
  wgcna_program_suggestions <- shiny::reactive({
    de_df <- tryCatch(de()$df, error = function(...) NULL)
    suggest_module_programs(wgcna_res(), de_df = de_df, alpha = input$de_alpha %||% 0.05, lfc_cutoff = input$de_lfc_cutoff %||% 0)
  })
  shiny::observe({
    sug <- wgcna_program_suggestions()
    choices <- if (nrow(sug)) stats::setNames(sug$modules, paste0(sug$modules, " | score ", sprintf("%.2f", sug$suggestion_score))) else character()
    selected <- if (length(choices)) unname(choices[[1L]]) else character()
    shiny::updateSelectInput(session, "wgcna_program_suggestion", choices = choices, selected = selected)
    if (length(selected)) {
      shiny::updateTextInput(session, "wgcna_program_name", value = paste("Program:", selected))
    }
  })
  shiny::observeEvent(input$create_wgcna_program, {
    shiny::req(nzchar(input$wgcna_program_suggestion))
    sug <- wgcna_program_suggestions()
    row <- sug[sug$modules == input$wgcna_program_suggestion, , drop = FALSE]
    label <- trimws(input$wgcna_program_name)
    if (!nzchar(label)) label <- paste("Program:", input$wgcna_program_suggestion)
    tbl <- wgcna_programs()
    tbl <- tbl[tbl$name != label, , drop = FALSE]
    tbl <- rbind(
      tbl,
      data.frame(
        name = label,
        modules = input$wgcna_program_suggestion,
        rationale = if (nrow(row)) row$rationale[1L] else "User-created combined program from suggested modules.",
        stringsAsFactors = FALSE
      )
    )
    wgcna_programs(tbl)
    shiny::updateSelectInput(session, "wgcna_program_view", choices = tbl$name, selected = label)
    analysis_log(append_analysis_log(analysis_log(), "create_wgcna_program", list(
      program_name = label,
      suggested_modules = input$wgcna_program_suggestion
    )))
  })
  shiny::observe({
    tbl <- wgcna_programs()
    choices <- tbl$name
    selected <- if (!is.null(input$wgcna_program_view) && input$wgcna_program_view %in% choices) input$wgcna_program_view else if (length(choices)) choices[1L] else character()
    shiny::updateSelectInput(session, "wgcna_program_view", choices = choices, selected = selected)
  })
  selected_wgcna_program <- shiny::reactive({
    tbl <- wgcna_programs()
    shiny::req(nrow(tbl) > 0L, nzchar(input$wgcna_program_view))
    row <- tbl[tbl$name == input$wgcna_program_view, , drop = FALSE]
    shiny::req(nrow(row) > 0L)
    modules <- trimws(unlist(strsplit(row$modules[1L], "\\+", perl = TRUE)))
    modules <- modules[nzchar(modules)]
    de_df <- tryCatch(de()$df, error = function(...) NULL)
    prog <- build_module_program(wgcna_res(), ds()$gene_map, modules = modules, name = row$name[1L], de_df = de_df, trait_cols = input$wgcna_traits %||% character(), alpha = input$de_alpha %||% 0.05, lfc_cutoff = input$de_lfc_cutoff %||% 0)
    prog$enrichment <- run_module_enrichment(
      module_genes = unique(prog$gene_table$gene_id),
      gene_map = ds()$gene_map,
      padj_cutoff = 0.05,
      min_gs = 10,
      max_gs = 500,
      p_adjust_method = "BH",
      top_n = 15,
      refresh = isTRUE(input$refresh_module_cache),
      cache_label = paste0("program_", row$name[1L])
    )
    prog$rationale <- row$rationale[1L]
    prog
  })
  output$wgcna_program_help <- shiny::renderUI({
    shiny::tagList(
      shiny::tags$p("What a combined program is: a higher-level biological interpretation layer built from two WGCNA modules that appear to behave like opposite arms of the same process."),
      shiny::tags$p("What the app checks before suggesting one: strong inverse eigengene correlation, similar trait association pattern, and opposite DE tendency when DESeq2 results are available."),
      shiny::tags$p("What does not happen: the original WGCNA modules are not overwritten. You still keep the network result, and this new layer simply adds a more biologically coherent downstream interpretation when it is justified."),
      shiny::tags$p("Beginner guide: if a suggestion makes sense biologically, create the program and then inspect its aligned eigengene pattern, trait summary, and enrichment before using it in interpretation."),
      shiny::tags$p("Role in the overall workflow: programs are useful when one biological process seems split into opposite module arms. They create a more realistic unit for downstream enrichment and SEMgraph without erasing the original WGCNA result.")
    )
  })
  output$wgcna_program_suggestion_table <- DT::renderDataTable({
    sug <- wgcna_program_suggestions()
    shiny::req(nrow(sug) > 0L)
    make_downloadable_table(sug[, c("modules", "eigengene_correlation", "dominant_trait", "opposite_trait_sign", "opposite_de_direction", "suggestion_score", "rationale"), drop = FALSE], page_length = 8)
  }, server = FALSE)
  output$wgcna_program_summary <- shiny::renderUI({
    prog <- selected_wgcna_program()
    comp <- infer_module_de_direction(wgcna_res(), tryCatch(de()$df, error = function(...) NULL))
    comp <- comp[comp$module %in% prog$modules, , drop = FALSE]
    shiny::tagList(
      shiny::tags$h4(prog$name),
      shiny::tags$p(sprintf("Modules combined: %s", paste(prog$modules, collapse = " + "))),
      shiny::tags$p(sprintf("Why this was suggested: %s", prog$rationale)),
      shiny::tags$p("How to interpret this: the app flips inverse module eigengenes to the same biological direction before summarizing them together. This helps represent one process that may contain both upregulated and downregulated arms in the original WGCNA result."),
      shiny::tags$p(sprintf("Genes in this combined program: %s", nrow(prog$gene_table))),
      if (nrow(comp)) shiny::tags$ul(lapply(seq_len(nrow(comp)), function(i) shiny::tags$li(sprintf("Module %s shows mean DE direction %s (mean log2FC %.2f).", comp$module[i], comp$direction[i], comp$mean_log2fc[i]))))
    )
  })
  output$wgcna_program_diag_plot <- plotly::renderPlotly({
    prog <- selected_wgcna_program()
    df <- as.data.frame(prog$aligned_eigengenes, stringsAsFactors = FALSE)
    df$sample <- rownames(prog$aligned_eigengenes)
    long <- do.call(rbind, lapply(setdiff(names(df), "sample"), function(nm) {
      data.frame(sample = df$sample, component = sub("^ME", "", nm), value = df[[nm]], stringsAsFactors = FALSE)
    }))
    long <- rbind(long, data.frame(sample = names(prog$score), component = "combined_program", value = as.numeric(prog$score), stringsAsFactors = FALSE))
    plotly::layout(
      plotly::plot_ly(long, x = ~sample, y = ~value, color = ~component, type = "scatter", mode = "lines+markers", hovertext = ~paste("Sample:", sample, "<br>Component:", component, "<br>Aligned eigengene score:", signif(value, 3)), hoverinfo = "text"),
      title = "Aligned module eigengene view for the combined program",
      xaxis = list(title = "", tickangle = -45),
      yaxis = list(title = "Aligned eigengene / program score")
    )
  })
  output$wgcna_program_trait_plot <- plotly::renderPlotly({
    prog <- selected_wgcna_program()
    shiny::req(!is.null(prog$trait_cor), length(prog$trait_cor) > 0L)
    trait_names <- colnames(prog$trait_p)
    if (is.null(trait_names) || !length(trait_names)) trait_names <- names(as.numeric(prog$trait_cor))
    vals <- data.frame(trait = trait_names, correlation = as.numeric(prog$trait_cor), pvalue = as.numeric(prog$trait_p[1, ]), stringsAsFactors = FALSE)
    vals <- vals[order(abs(vals$correlation), decreasing = TRUE), , drop = FALSE]
    plotly::layout(
      plotly::plot_ly(vals, x = ~reorder(trait, correlation), y = ~correlation, type = "bar", text = ~paste("Trait:", trait, "<br>Correlation:", signif(correlation, 3), "<br>p-value:", signif(pvalue, 3)), hoverinfo = "text"),
      title = "Combined program trait associations",
      xaxis = list(title = "", tickangle = -45),
      yaxis = list(title = "Correlation")
    )
  })
  output$wgcna_program_enrichment_plot <- plotly::renderPlotly({
    prog <- selected_wgcna_program()
    enrich_res <- prog$enrichment$combined
    shiny::req(nrow(enrich_res) > 0L)
    top <- head(enrich_res, 12L)
    top$pathway <- factor(top$pathway, levels = rev(unique(top$pathway)))
    top$score <- -log10(pmax(as.numeric(top$p.adjust), .Machine$double.xmin))
    plotly::layout(
      plotly::plot_ly(top, x = ~score, y = ~pathway, color = ~method, size = ~gene_ratio_value, type = "scatter", mode = "markers", text = ~paste("Pathway:", pathway, "<br>Method:", method, "<br>Adjusted p:", signif(p.adjust, 3)), hoverinfo = "text"),
      title = paste("Pathway enrichment for", prog$name),
      xaxis = list(title = "-log10 adjusted p-value"),
      yaxis = list(title = "")
    )
  })
  output$wgcna_program_table <- DT::renderDataTable({
    prog <- selected_wgcna_program()
    gt <- prog$gene_table[, c("gene_label", "gene_id", "module", "module_orientation", "kWithin", "log2FoldChange", "padj"), drop = FALSE]
    names(gt) <- c("gene", "ensembl_id", "original_module", "program_arm", "within_module_connectivity", "log2FoldChange", "padj")
    make_downloadable_table(gt, page_length = 12)
  }, server = FALSE)

  shiny::observe({
    gm <- wgcna_res()$gene_modules
    mods <- sort(unique(gm$module))
    mods <- mods[mods != "grey"]
    program_tbl <- wgcna_programs()
    module_choices <- if (length(mods)) stats::setNames(paste0("module:", mods), paste0("Module: ", mods)) else character()
    program_choices <- if (nrow(program_tbl)) stats::setNames(paste0("program:", program_tbl$name), paste0("Program: ", program_tbl$name)) else character()
    choices <- c(module_choices, program_choices)
    selected <- if (!is.null(input$semgraph_unit) && input$semgraph_unit %in% unname(choices)) input$semgraph_unit else if (length(choices)) unname(choices[[1L]]) else character()
    shiny::updateSelectInput(session, "semgraph_unit", choices = choices, selected = selected)
    meta_cols <- names(ds()$meta)
    binary_cols <- meta_cols[vapply(meta_cols, function(col) !is.null(pick_binary_group(ds()$meta, col)), logical(1))]
    shiny::updateSelectInput(session, "semgraph_group", choices = c(stats::setNames("", ""), stats::setNames(binary_cols, binary_cols)), selected = if ("treatment" %in% binary_cols) "treatment" else "")
  })
  shiny::observeEvent(input$update_reactome_cache, {
    shiny::withProgress(message = "Updating Reactome cache", detail = "Downloading and parsing the local Reactome mapping file. Please wait.", value = 0, {
      progress_step <- function(value, detail) shiny::setProgress(value = min(max(value, 0), 1), detail = detail)
      tryCatch({
        load_reactome_mapping(force_refresh = TRUE, progress = progress_step)
        shiny::showNotification("Reactome cache updated successfully.", type = "message", duration = 5)
        analysis_log(append_analysis_log(analysis_log(), "update_reactome_cache", list(force_refresh = TRUE)))
      }, error = function(e) {
        shiny::showNotification(
          paste(
            "Reactome cache update failed:",
            conditionMessage(e),
            "SEMgraph can still run without refreshing Reactome, but Reactome pathway-context annotations may be limited until a local cache is available."
          ),
          type = "error",
          duration = 12
        )
      })
    })
  })
  shiny::observeEvent(input$semgraph_download_string_local, {
    res <- tryCatch(semgraph_res(), error = function(...) NULL)
    if (is.null(res) || igraph::vcount(res$graph) < 2L) {
      shiny::showNotification("Run SEMgraph first so the app knows which genes to use for the local STRING export.", type = "warning", duration = 8)
      return()
    }
    symbols <- sort(unique(igraph::V(res$graph)$name))
    destfile <- file.path(data_dir, "STRING_edges.tsv")
    shiny::withProgress(message = "Preparing local STRING validation file", detail = "Collecting current SEMgraph genes", value = 0, {
      progress_step <- function(value, detail) {
        shiny::setProgress(value = min(max(value, 0), 1), detail = detail)
      }
      progress_step(0.02, sprintf("Preparing STRING query for %d genes", length(symbols)))
      endpoint <- paste0(
        "https://string-db.org/api/tsv/network?identifiers=",
        utils::URLencode(paste(symbols, collapse = "\r"), reserved = TRUE),
        "&species=9606&required_score=0&network_type=functional"
      )
      tryCatch({
        download_url_with_progress(endpoint, destfile, progress = function(frac, detail) {
          progress_step(0.05 + 0.85 * frac, detail)
        })
        progress_step(0.94, "Parsing and validating downloaded STRING file")
        tbl <- utils::read.delim(destfile, stringsAsFactors = FALSE, check.names = FALSE)
        if (is.null(tbl) || !nrow(tbl)) stop("STRING returned an empty interaction table.")
        progress_step(1, sprintf("Saved local STRING file to %s", basename(destfile)))
        shiny::showNotification("Saved data/STRING_edges.tsv. Rerun SEMgraph with external refresh if you want to rebuild edge validation immediately.", type = "message", duration = 10)
      }, error = function(e) {
        if (file.exists(destfile) && file.info(destfile)$size == 0) file.remove(destfile)
        shiny::showNotification(paste("STRING download failed:", conditionMessage(e)), type = "error", duration = 12)
      })
    })
  })
  semgraph_res <- shiny::eventReactive(input$run_semgraph, {
    shiny::withProgress(message = "Running SEMgraph causal analysis", detail = "Orienting prior graph and estimating regulator structure. Please wait.", value = 0, {
      cur <- ds()
      wres <- wgcna_res()
      de_df <- tryCatch(de()$df, error = function(...) NULL)
      alpha <- input$de_alpha %||% 0.05
      lfc_cut <- input$de_lfc_cutoff %||% 0
      progress_step <- function(value, detail) {
        shiny::setProgress(value = min(max(value, 0), 1), detail = detail)
      }
      shiny::incProgress(0.02, detail = "Preparing selected analysis unit")
      run_semgraph_analysis(
        unit_id = input$semgraph_unit,
        wgcna_res = wres,
        gene_map = cur$gene_map,
        metadata = cur$meta,
        program_tbl = wgcna_programs(),
        de_df = de_df,
        trait_cols = input$wgcna_traits %||% character(),
        refresh_external = isTRUE(input$refresh_semgraph_cache),
        refresh_fit = isTRUE(input$refresh_semgraph_fit),
        max_model_genes = input$semgraph_model_max_genes,
        annotation_max_genes = input$semgraph_annotation_max_genes,
        reactome_live = isTRUE(input$semgraph_use_reactome_live),
        use_live_string = isTRUE(input$semgraph_use_string_live),
        group_col = input$semgraph_group,
        alpha = alpha,
        lfc_cutoff = lfc_cut,
        progress = progress_step
      )
    })
  })
  shiny::observeEvent(input$run_semgraph, {
    analysis_log(append_analysis_log(analysis_log(), "run_semgraph", list(
      analysis_unit = input$semgraph_unit,
      group = input$semgraph_group,
      display_mode = input$semgraph_display_mode,
      max_model_genes = input$semgraph_model_max_genes,
      max_annotation_genes = input$semgraph_annotation_max_genes,
      use_live_string = isTRUE(input$semgraph_use_string_live),
      graph_scope = input$semgraph_graph_scope,
      graph_max_nodes = input$semgraph_graph_max_nodes
    )))
  }, ignoreInit = TRUE)
  shiny::observeEvent(semgraph_res(), {
    res <- semgraph_res()
    shiny::req(!is.null(res))
    node_choices <- stats::setNames(res$degree_df$gene, res$degree_df$gene)
    default_end <- setdiff(res$degree_df$gene, res$focal_gene)[1L] %||% res$focal_gene
    shiny::updateSelectInput(session, "semgraph_start_gene", choices = node_choices, selected = res$focal_gene)
    shiny::updateSelectInput(session, "semgraph_end_gene", choices = node_choices, selected = default_end)
    shiny::updateSelectizeInput(session, "semgraph_selected_nodes", choices = node_choices, selected = character(), server = TRUE)
  })
  shiny::observeEvent(input$semgraph_clear_selected_nodes, {
    shiny::updateSelectizeInput(session, "semgraph_selected_nodes", selected = character(), server = TRUE)
  })
  semgraph_filter_state <- shiny::reactiveVal(NULL)
  semgraph_color_state <- shiny::reactiveVal(NULL)
  current_semgraph_focus_inputs <- function() {
    list(
      selected_nodes = input$semgraph_selected_nodes %||% character(),
      selected_hops = input$semgraph_selected_hops %||% 3L,
      selected_only = isTRUE(input$semgraph_selected_only),
      filter_tf = isTRUE(input$semgraph_filter_tf_only),
      filter_druggable = isTRUE(input$semgraph_filter_druggable_only),
      filter_hub_mode = if (isTRUE(input$semgraph_filter_hub_only)) "wgcna_hubs" else "off",
      focal_top_n = if (isTRUE(input$semgraph_filter_regulator_only)) (input$semgraph_filter_focal_top_n %||% 0L) else 0L,
      restrict_to_focus = isTRUE(input$semgraph_restrict_focus_only),
      edge_rule = input$semgraph_edge_validation_filter %||% "all",
      anchor_color = input$semgraph_anchor_color %||% "#0f9d8a",
      tf_color = input$semgraph_tf_color %||% "#7e57c2",
      druggable_color = input$semgraph_druggable_color %||% "#ffd24a",
      hub_color = input$semgraph_hub_color %||% "#ff7043",
      regulator_color = input$semgraph_regulator_color %||% "#d81b60"
    )
  }
  current_semgraph_filter_rule <- function() {
    filter_state <- semgraph_filter_state()
    if (is.null(filter_state)) "all" else (filter_state$edge_rule %||% "all")
  }
  current_semgraph_color_rule <- function() {
    color_state <- semgraph_color_state()
    if (is.null(color_state)) "all" else (color_state$edge_rule %||% "all")
  }
  shiny::observeEvent(input$semgraph_apply_filter, {
    semgraph_filter_state(current_semgraph_focus_inputs())
  }, ignoreInit = TRUE)
  shiny::observeEvent(input$semgraph_apply_color, {
    semgraph_color_state(current_semgraph_focus_inputs())
  }, ignoreInit = TRUE)
  shiny::observeEvent(input$semgraph_reset_focus_mode, {
    semgraph_filter_state(NULL)
    semgraph_color_state(NULL)
  }, ignoreInit = TRUE)
  semgraph_graph_selection <- shiny::reactive({
    res <- semgraph_res()
    shiny::req(!is.null(res))
    de_df <- tryCatch(de()$df, error = function(...) NULL)
    alpha <- input$de_alpha %||% 0.05
    lfc_cut <- input$de_lfc_cutoff %||% 0
    focus_criteria <- semgraph_filter_state() %||% list()
    filter_active <- !is.null(semgraph_filter_state())
    active_edge_rule <- if (filter_active) (focus_criteria$edge_rule %||% "all") else "all"
    sel <- select_semgraph_graph_nodes(
      res,
      de_df = de_df,
      display_mode = input$semgraph_display_mode,
      graph_scope = input$semgraph_graph_scope,
      max_nodes = input$semgraph_graph_max_nodes,
      start_gene = input$semgraph_start_gene,
      end_gene = input$semgraph_end_gene,
      hops = input$semgraph_neighbor_hops %||% 2L,
      path_n = input$semgraph_path_n %||% 5L,
      path_max_steps = input$semgraph_path_max_steps %||% 6L,
      selected_nodes = focus_criteria$selected_nodes %||% character(),
      selected_hops = focus_criteria$selected_hops %||% 3L,
      selected_only = isTRUE(focus_criteria$selected_only) && filter_active,
      filter_tf = isTRUE(focus_criteria$filter_tf) && filter_active,
      filter_druggable = isTRUE(focus_criteria$filter_druggable) && filter_active,
      filter_hub_mode = if (filter_active) (focus_criteria$filter_hub_mode %||% "off") else "off",
      focal_top_n = if (filter_active) (focus_criteria$focal_top_n %||% 0L) else 0L,
      restrict_to_focus = isTRUE(focus_criteria$restrict_to_focus) && filter_active,
      alpha = alpha,
      lfc_cutoff = lfc_cut
    )
    if (!identical(active_edge_rule, "all")) {
      full_edge_tbl <- filter_validated_edges(
        annotate_semgraph_edges(res$graph, res),
        mode = active_edge_rule,
        string_cutoff = input$semgraph_string_cutoff %||% 0.4
      )
      current_edge_tbl <- filter_validated_edges(
        annotate_semgraph_edges(igraph::induced_subgraph(res$graph, vids = intersect(sel$nodes, igraph::V(res$graph)$name)), res),
        mode = active_edge_rule,
        string_cutoff = input$semgraph_string_cutoff %||% 0.4
      )
      if (!nrow(current_edge_tbl) && nrow(full_edge_tbl)) {
        ranked <- rank_semgraph_nodes(res, de_df, alpha = alpha, lfc_cutoff = lfc_cut)
        full_nodes <- unique(c(full_edge_tbl$source, full_edge_tbl$target))
        ranked <- ranked[ranked$gene %in% full_nodes, , drop = FALSE]
        ranked <- ranked[order(ranked$de_supported, ranked$total_degree, decreasing = TRUE), , drop = FALSE]
        keep_nodes <- unique(c(res$focal_gene, head(ranked$gene, input$semgraph_graph_max_nodes %||% 60L)))
        sel$nodes <- keep_nodes[keep_nodes %in% igraph::V(res$graph)$name]
      }
    }
    sel
  })
  semgraph_color_targets <- shiny::reactive({
    res <- semgraph_res()
    shiny::req(!is.null(res))
    focus_criteria <- semgraph_color_state()
    if (is.null(focus_criteria)) {
      return(list(nodes = character(), apply_mode = "none", edge_rule = "all", categories = list()))
    }
    de_df <- tryCatch(de()$df, error = function(...) NULL)
    alpha <- input$de_alpha %||% 0.05
    lfc_cut <- input$de_lfc_cutoff %||% 0
    ranked <- rank_semgraph_nodes(res, de_df, alpha = alpha, lfc_cutoff = lfc_cut)
    ann <- res$local_annotations[match(ranked$gene, rownames(res$local_annotations)), , drop = FALSE]
    graph_nodes <- igraph::V(res$graph)$name
    selected_nodes <- unique((focus_criteria$selected_nodes %||% character())[(focus_criteria$selected_nodes %||% character()) %in% graph_nodes])
    tf_nodes <- if (isTRUE(focus_criteria$filter_tf)) unique(ann$hgnc_symbol[ifelse(is.na(ann$is_tf), FALSE, ann$is_tf)]) else character()
    tf_nodes <- tf_nodes[tf_nodes %in% graph_nodes]
    drug_nodes <- if (isTRUE(focus_criteria$filter_druggable)) unique(ann$hgnc_symbol[ifelse(is.na(ann$is_druggable), FALSE, ann$is_druggable)]) else character()
    drug_nodes <- drug_nodes[drug_nodes %in% graph_nodes]
    hub_nodes <- if (identical(focus_criteria$filter_hub_mode %||% "off", "wgcna_hubs")) unique(ann$hgnc_symbol[ifelse(is.na(ann$is_hub), FALSE, ann$is_hub)]) else character()
    hub_nodes <- hub_nodes[hub_nodes %in% graph_nodes]
    regulator_nodes <- character()
    reg_n <- focus_criteria$focal_top_n %||% 0L
    if (isTRUE(reg_n > 0L)) {
      regulator_score <- ranked$out_degree - ranked$in_degree + ranked$total_degree / max(1, max(ranked$total_degree, na.rm = TRUE))
      ord <- order(regulator_score, decreasing = TRUE, na.last = NA)
      regulator_nodes <- unique(head(ranked$gene[ord], reg_n))
      regulator_nodes <- regulator_nodes[regulator_nodes %in% graph_nodes]
    }
    categories <- list(
      anchor = list(nodes = selected_nodes, color = focus_criteria$anchor_color %||% "#0f9d8a"),
      tf = list(nodes = tf_nodes, color = focus_criteria$tf_color %||% "#7e57c2"),
      druggable = list(nodes = drug_nodes, color = focus_criteria$druggable_color %||% "#ffd24a"),
      hub = list(nodes = hub_nodes, color = focus_criteria$hub_color %||% "#ff7043"),
      regulator = list(nodes = regulator_nodes, color = focus_criteria$regulator_color %||% "#d81b60")
    )
    focus_nodes <- unique(unlist(lapply(categories, `[[`, "nodes"), use.names = FALSE))
    list(
      nodes = focus_nodes,
      apply_mode = "color",
      edge_rule = focus_criteria$edge_rule %||% "all",
      categories = categories
    )
  })
  output$semgraph_focus_mode_status <- shiny::renderUI({
    filter_state <- semgraph_filter_state()
    color_state <- semgraph_color_state()
    edge_label <- function(rule) switch(
      rule %||% "all",
      any_support = "Any STRING support",
      string_only = "STRING only above cutoff",
      "All SEMgraph edges"
    )
    parts <- character()
    if (!is.null(filter_state)) parts <- c(parts, paste("Filter active | edge rule:", edge_label(filter_state$edge_rule)))
    if (!is.null(color_state)) parts <- c(parts, paste("Color active | edge rule:", edge_label(color_state$edge_rule)))
    label <- if (length(parts)) paste(parts, collapse = " || ") else "Active mode: none"
    shiny::tags$div(class = "muted top-input-note", label)
  })
  visible_semgraph_edge_table <- shiny::reactive({
    res <- semgraph_res()
    shiny::req(!is.null(res))
    sel <- semgraph_graph_selection()
    sub_nodes <- intersect(sel$nodes, igraph::V(res$graph)$name)
    if (!length(sub_nodes)) return(data.frame())
    g_sub <- igraph::induced_subgraph(res$graph, vids = sub_nodes)
    edge_tbl <- annotate_semgraph_edges(g_sub, res)
    if (!nrow(edge_tbl)) return(edge_tbl)
    filter_rule <- current_semgraph_filter_rule()
    if (!identical(filter_rule, "all")) {
      edge_tbl <- filter_validated_edges(edge_tbl, mode = filter_rule, string_cutoff = input$semgraph_string_cutoff %||% 0.7)
    }
    color_rule <- current_semgraph_color_rule()
    cutoff_now <- input$semgraph_string_cutoff %||% 0.4
    edge_tbl$matches_current_color_rule <- if (identical(color_rule, "string_only")) {
      edge_tbl$string_supported & edge_tbl$string_score >= cutoff_now
    } else if (identical(color_rule, "any_support")) {
      edge_tbl$string_supported
    } else {
      FALSE
    }
    edge_tbl
  })
  output$semgraph_graph_scope_ui <- shiny::renderUI({
    res <- semgraph_res()
    shiny::req(!is.null(res))
    choices <- stats::setNames(res$degree_df$gene, res$degree_df$gene)
    selection_rows <- shiny::tagList(
      shiny::fluidRow(
        shiny::column(6, shiny::selectizeInput("semgraph_selected_nodes", "Anchor genes", choices = choices, selected = input$semgraph_selected_nodes %||% character(), multiple = TRUE)),
        shiny::column(2, shiny::numericInput("semgraph_selected_hops", "Hops", value = 3, min = 0, max = 10, step = 1)),
        shiny::column(2, shiny::checkboxInput("semgraph_selected_only", "Use anchors", value = FALSE)),
        shiny::column(2, color_text_input("semgraph_anchor_color", "Color", value = "#0f9d8a"))
      ),
      shiny::fluidRow(
        shiny::column(6, shiny::checkboxInput("semgraph_filter_tf_only", "TF genes", value = FALSE)),
        shiny::column(2, shiny::tags$div(class = "muted top-input-note", "All matched TFs")),
        shiny::column(2, shiny::tags$div(class = "muted top-input-note", "")),
        shiny::column(2, color_text_input("semgraph_tf_color", "Color", value = "#7e57c2"))
      ),
      shiny::fluidRow(
        shiny::column(6, shiny::checkboxInput("semgraph_filter_druggable_only", "Druggable genes", value = FALSE)),
        shiny::column(2, shiny::tags$div(class = "muted top-input-note", "All matched targets")),
        shiny::column(2, shiny::tags$div(class = "muted top-input-note", "")),
        shiny::column(2, color_text_input("semgraph_druggable_color", "Color", value = "#ffd24a"))
      ),
      shiny::fluidRow(
        shiny::column(6, shiny::checkboxInput("semgraph_filter_hub_only", "WGCNA hub genes", value = FALSE)),
        shiny::column(2, shiny::tags$div(class = "muted top-input-note", "Annotated hubs")),
        shiny::column(2, shiny::tags$div(class = "muted top-input-note", "")),
        shiny::column(2, color_text_input("semgraph_hub_color", "Color", value = "#ff7043"))
      ),
      shiny::fluidRow(
        shiny::column(6, shiny::checkboxInput("semgraph_filter_regulator_only", "Top SEMgraph regulators", value = FALSE)),
        shiny::column(2, shiny::numericInput("semgraph_filter_focal_top_n", "Top N", value = 5, min = 0, max = 50, step = 1)),
        shiny::column(2, shiny::tags$div(class = "muted top-input-note", "By out-in degree")),
        shiny::column(2, color_text_input("semgraph_regulator_color", "Color", value = "#d81b60"))
      )
    )
    common_focus <- shiny::fluidRow(
      shiny::column(4, shiny::sliderInput("semgraph_graph_max_nodes", "Graph node limit", min = 15, max = 120, value = 60, step = 5)),
      shiny::column(4, shiny::selectInput(
        "semgraph_edge_validation_filter",
        "STRING edge criterion",
        choices = c(
          "All SEMgraph edges" = "all",
          "Any STRING support" = "any_support",
          "STRING only above cutoff" = "string_only"
        ),
        selected = "all"
      )),
      shiny::column(4, shiny::sliderInput("semgraph_string_cutoff", "Minimum STRING confidence", min = 0, max = 1, value = 0.4, step = 0.05)),
      shiny::column(12, selection_rows),
      shiny::column(12, shiny::checkboxInput("semgraph_restrict_focus_only", "When filtering, keep only chosen categories and their connectors", value = FALSE)),
      shiny::column(12,
        shiny::div(class = "top-input-note",
          shiny::tags$label("Apply current criteria"),
          shiny::div(
            style = "display:flex; gap:8px; flex-wrap:wrap; align-items:center;",
            shiny::actionButton("semgraph_apply_filter", "Apply Filter"),
            shiny::actionButton("semgraph_apply_color", "Apply Color"),
            shiny::actionButton("semgraph_clear_selected_nodes", "Clear anchors"),
            shiny::actionButton("semgraph_reset_focus_mode", "Reset")
          ),
          shiny::uiOutput("semgraph_focus_mode_status")
        )
      )
    )
    if (identical(input$semgraph_graph_scope, "neighborhood")) {
      return(shiny::tagList(
        shiny::fluidRow(
          shiny::column(6, shiny::selectInput("semgraph_start_gene", "Seed gene", choices = choices, selected = input$semgraph_start_gene %||% res$focal_gene)),
          shiny::column(3, shiny::numericInput("semgraph_neighbor_hops", "Explore out to this many hops", value = 2, min = 1, max = 8, step = 1)),
          shiny::column(3, shiny::tags$p(class = "muted top-input-note", "Use this when you want to ask: what sits near this gene in the fitted causal graph?"))
        ),
        common_focus
      ))
    }
    if (identical(input$semgraph_graph_scope, "paths")) {
      default_end <- setdiff(res$degree_df$gene, res$focal_gene)[1L] %||% res$focal_gene
      return(shiny::tagList(
        shiny::fluidRow(
          shiny::column(4, shiny::selectInput("semgraph_start_gene", "Start gene", choices = choices, selected = input$semgraph_start_gene %||% res$focal_gene)),
          shiny::column(4, shiny::selectInput("semgraph_end_gene", "End gene", choices = choices, selected = input$semgraph_end_gene %||% default_end)),
          shiny::column(2, shiny::numericInput("semgraph_path_n", "Show top path count", value = 5, min = 1, max = 20, step = 1)),
          shiny::column(2, shiny::numericInput("semgraph_path_max_steps", "Allow path length up to", value = 6, min = 2, max = 12, step = 1)),
          shiny::column(12, shiny::tags$p(class = "muted", "Use this when you want to ask: how could the graph connect one candidate regulator to one downstream target?"))
        ),
        common_focus
      ))
    }
    shiny::tagList(
      shiny::fluidRow(
      shiny::column(12, shiny::tags$p(class = "muted", "Overview mode is best for the first pass. Use the focus-gene selector below if you want to collapse the graph to only selected genes plus their connector genes."))
      ),
      common_focus
    )
  })
  output$semgraph_help <- shiny::renderUI({
    shiny::tagList(
      shiny::tags$p("What this step tries to do: starting from one original WGCNA module or one saved combined program, it builds a graph-informed SEM model to suggest which genes may sit upstream or downstream of the main biological signal."),
      shiny::tags$p("Role in the overall workflow: this is the most interpretive layer in the tool. It should usually be read only after QC, DESeq2, enrichment, and WGCNA already suggest a believable biological story."),
      shiny::tags$p("This is exploratory causal analysis, not proof. Use it to prioritize regulators and pathways for follow-up, especially when they agree with module hubs and enrichment results."),
      shiny::tags$p("Analysis-unit guide: choose a single module when you want the cleanest network-specific view. Choose a saved combined program when two modules appear to be opposite arms of the same biology and you want SEMgraph to consider them together."),
      shiny::tags$p("Fit-size guide: SEMgraph becomes slow on very large gene sets. The app therefore prioritizes genes before fitting the causal model, favoring DE-supported genes and highly connected genes. This keeps the model tractable while still targeting the most informative part of the biology."),
      shiny::tags$p("Annotation-size guide: external drug/pathway annotation can also be slow, especially for Reactome. The app therefore annotates only the top-ranked fitted genes by default. Genes outside that annotation cap can still appear in SEMgraph, but may show fewer external context fields."),
      shiny::tags$p("Validation guide: this app can now validate SEMgraph edges against STRING functional associations. STRING supports biological relatedness, but not causal direction by itself."),
      shiny::tags$p("Validation-data guide: use the setup button to download a local STRING file for the current SEMgraph genes. The STRING downloader now reports downloaded size, speed, and ETA during transfer."),
      shiny::tags$p("Reactome-speed guide: the app now prefers a local downloaded Reactome mapping cache. Use `Update Reactome cache` once, then later SEMgraph runs can reuse that local file quickly. Live Reactome lookup remains optional only as a fallback."),
      shiny::tags$p("Progress guide: the SEMgraph progress bar now reports the main fit stages. If you see `SEMgraph node ordering in progress`, the model is still actively working on the slowest internal step rather than being idle."),
      shiny::tags$p("Graph-vs-table guide: the graph is constrained for readability and zooming, while the tables can be broader and are filtered by evidence rather than by what fits visually. This means the tables may contain relevant genes that are not currently drawn in the graph."),
      shiny::tags$p("Pharmacology layer: the app can combine DGIdb, a local HCDT export when available in `data/`, and Reactome drug-mediated pathway context. Each displayed pharmacology hit carries its source so you can see what came from where.")
    )
  })
  output$semgraph_summary <- shiny::renderUI({
    res <- semgraph_res()
    shiny::req(!is.null(res))
    de_df <- tryCatch(de()$df, error = function(...) NULL)
    alpha <- input$de_alpha %||% 0.05
    lfc_cut <- input$de_lfc_cutoff %||% 0
    source_status <- validation_source_status(igraph::V(res$graph)$name, allow_live_string = isTRUE(input$semgraph_use_string_live))
    display_nodes <- select_semgraph_display_nodes(res, de_df, mode = input$semgraph_display_mode, max_nodes = input$semgraph_graph_max_nodes, alpha = alpha, lfc_cutoff = lfc_cut)
    supported_n <- NA_integer_
    local_n <- length(display_nodes)
    display_ann <- res$local_annotations[match(display_nodes, rownames(res$local_annotations)), , drop = FALSE]
    druggable_n <- sum(display_ann$is_druggable, na.rm = TRUE)
    validated_edges <- if (!is.null(res$edge_validation)) filter_validated_edges(res$edge_validation, mode = current_semgraph_filter_rule(), string_cutoff = input$semgraph_string_cutoff %||% 0.7) else data.frame()
    raw_edge_validation <- res$edge_validation %||% data.frame()
    string_supported_n <- if (nrow(raw_edge_validation)) sum(raw_edge_validation$string_supported, na.rm = TRUE) else 0L
    string_above_cutoff_n <- if (nrow(raw_edge_validation)) sum(raw_edge_validation$string_supported & raw_edge_validation$string_score >= (input$semgraph_string_cutoff %||% 0.4), na.rm = TRUE) else 0L
    best_string_score <- if (nrow(raw_edge_validation)) max(raw_edge_validation$string_score, na.rm = TRUE) else 0
    if (!is.finite(best_string_score)) best_string_score <- 0
    if (!is.null(de_df) && nrow(de_df) > 0L) {
      local_gene_ids <- display_ann$gene_id
      de_match <- de_df[match(local_gene_ids, de_df$gene_id), , drop = FALSE]
      supported_n <- sum(deg_mask(de_match, alpha = alpha, lfc_cutoff = lfc_cut), na.rm = TRUE)
    }
    shiny::tagList(
      shiny::tags$div(class = "status-note", "SEMgraph fitting has completed for the selected module. If you can see this summary, the causal model run is finished."),
      if (is.finite(supported_n)) shiny::tags$div(
        class = "status-note",
        if (supported_n <= max(2L, floor(0.2 * local_n))) {
          sprintf("Only %s of %s genes in the displayed causal neighborhood meet the current DEG thresholds (padj < %.3f and |log2FC| >= %.2f). Treat most nodes here as prior/module context rather than transcriptomics-driven discoveries.", supported_n, local_n, alpha, lfc_cut)
        } else {
          sprintf("%s of %s genes in the displayed causal neighborhood meet the current DEG thresholds (padj < %.3f and |log2FC| >= %.2f). This gives the causal view a stronger data-supported component.", supported_n, local_n, alpha, lfc_cut)
        }
      ),
      if (is.finite(supported_n) && supported_n == 0L) shiny::tags$div(class = "status-note", "No displayed SEMgraph genes are supported by the current DE result. In this case, use the causal graph only as prior-knowledge/module context and avoid interpreting it as a transcriptomics-driven regulatory discovery."),
      if (druggable_n > 0L) shiny::tags$div(class = "status-note", sprintf("%s displayed genes have external drug-target or drug-pathway annotations. Those gold-ring annotations come from sources such as DGIdb, HCDT, and Reactome and should be interpreted as actionability context, not as proof from your RNA-seq samples.", druggable_n)),
      if (igraph::ecount(res$graph) < 3L) shiny::tags$div(class = "status-note", "The fitted causal graph is sparse, so upstream/downstream interpretation should be treated cautiously. This often means the current module has limited directional information under the present data and prior assumptions."),
      if (!source_status$string_local && !source_status$string_cache && !source_status$string_live_enabled) shiny::tags$div(class = "status-note", "STRING validation is currently unavailable: no local STRING export was found in `data/`, no cached STRING network is available, and live STRING lookup is off. Zero STRING-supported edges in this state do not diagnose the biology."),
      if ((source_status$string_local || source_status$string_cache || source_status$string_live_enabled) && string_supported_n > 0L && string_above_cutoff_n == 0L) shiny::tags$div(class = "status-note", sprintf("STRING found %s matched SEMgraph edges, but none reach the current STRING confidence cutoff of %.2f. The best matched STRING score is %.3f. Lower the STRING cutoff if you want a less strict validation view.", string_supported_n, input$semgraph_string_cutoff %||% 0.4, best_string_score)),
      if (nrow(validated_edges) == 0L && (source_status$string_local || source_status$string_cache || source_status$string_live_enabled)) shiny::tags$div(class = "status-note", "No SEMgraph edges passed the current STRING validation filter. That can still happen even when STRING is available, for example because the current confidence cutoff is strict, the selected causal edges are novel, or the current SEMgraph edges are not represented in STRING."),
      shiny::tags$p(sprintf("Analyzed unit: %s", res$analysis_unit)),
      shiny::tags$p(sprintf("Analysis-unit type: %s", ifelse(identical(res$analysis_unit_type, "program"), "combined program", "original module"))),
      if (identical(res$analysis_unit_type, "program")) shiny::tags$p(sprintf("Program source modules: %s", res$source_modules)),
      shiny::tags$p(sprintf("Genes in original analysis unit: %s", res$semgraph_original_gene_count)),
      shiny::tags$p(sprintf("Genes entered into SEMgraph fit after prioritization: %s", res$semgraph_input_gene_count)),
      shiny::tags$p(sprintf("Genes receiving external pharmacology/pathway annotation: up to %s", input$semgraph_annotation_max_genes)),
      shiny::tags$p(sprintf("Genes in causal graph: %s", igraph::vcount(res$graph))),
      shiny::tags$p(sprintf("Edges in causal graph: %s", igraph::ecount(res$graph))),
      shiny::tags$p(sprintf("Edges passing current external-validation filter: %s", nrow(validated_edges))),
      shiny::tags$p(sprintf("Matched STRING-supported edges before cutoff filtering: %s", string_supported_n)),
      shiny::tags$p(sprintf("Matched STRING-supported edges at current cutoff %.2f: %s", input$semgraph_string_cutoff %||% 0.4, string_above_cutoff_n)),
      shiny::tags$p(sprintf("STRING source status: %s", if (source_status$string_local) paste("local export", basename(source_status$string_local_path)) else if (source_status$string_cache) "cached STRING network" else if (source_status$string_live_enabled) "live STRING lookup enabled" else "not available")),
      shiny::tags$p(sprintf("Genes displayed after DE-focused filtering: %s", local_n)),
      shiny::tags$p(sprintf("Current graph node limit: %s", input$semgraph_graph_max_nodes)),
      shiny::tags$p(sprintf("Current table evidence scope: %s", switch(input$semgraph_table_scope, de_only = "DE-only", de_plus_regulators = "DE + key regulators", evidence_context = "Broad relevant evidence", input$semgraph_table_scope))),
      shiny::tags$p(sprintf("Displayed genes with external drug-target context: %s", druggable_n)),
      shiny::tags$p(sprintf("Current display mode: %s", switch(input$semgraph_display_mode, de_only = "Strict DE-only", de_connectors = "DE + minimal connectors", full_local = "Full local neighborhood", input$semgraph_display_mode))),
      shiny::tags$p(sprintf("Focal gene: %s", res$focal_gene)),
      shiny::tags$p("What 'focal gene' means: this is the current central regulator candidate in the fitted local causal graph, chosen here as the gene with the strongest overall directed connectivity. It is a prioritization anchor for exploration, not guaranteed proof of being the true master regulator."),
      shiny::tags$p(sprintf("Binary group used for perturbation: %s", ifelse(nzchar(res$group_col), res$group_col, "none"))),
      shiny::tags$p("Interpretation: genes with many outgoing links are candidate upstream regulators, while genes with many incoming links can behave more like downstream responders."),
      shiny::tags$p("P-value guide: DESeq2 adjusted p-values shown in this SEMgraph section come from the current differential expression contrast. WGCNA module-trait p-values come from correlations between module eigengenes and traits. SEMgraph does not reuse those p-values; its own effect-estimate outputs come from the fitted causal model."),
      shiny::tags$p("Drug-target guide: drug names and activating/inhibiting labels come from external resources. DGIdb and HCDT contribute drug-target evidence, while Reactome contributes drug-mediated pathway context. These annotations are not inferred from this dataset, but they can help you prioritize which data-supported regulators may be pharmacologically tractable.")
    )
  })
  output$semgraph_diag_help <- shiny::renderUI({
    shiny::tagList(
      shiny::tags$p("Use these diagnostics before trusting the causal picture. The goal is to check whether the displayed SEMgraph story is still anchored in your DESeq2-supported data rather than drifting into weak context-only interpretation."),
      shiny::tags$ul(
        shiny::tags$li("Evidence support: shows how many displayed genes are actually supported by the current DESeq2 result versus kept mostly as context."),
        shiny::tags$li("Degree scatter: helps identify regulator-like genes with high out-degree and responder-like genes with high in-degree."),
      shiny::tags$li("Role by DE support: checks whether the upstream/downstream/context interpretation is mostly data-supported or mainly prior-graph context. The single focal gene is intentionally not treated as a group here.")
      )
    )
  })
  output$semgraph_support_plot <- plotly::renderPlotly({
    res <- semgraph_res()
    shiny::req(!is.null(res))
    de_df <- tryCatch(de()$df, error = function(...) NULL)
    alpha <- input$de_alpha %||% 0.05
    lfc_cut <- input$de_lfc_cutoff %||% 0
    local_nodes <- select_semgraph_display_nodes(res, de_df, mode = input$semgraph_display_mode, max_nodes = input$semgraph_graph_max_nodes, alpha = alpha, lfc_cutoff = lfc_cut)
    ann <- res$local_annotations[match(local_nodes, rownames(res$local_annotations)), , drop = FALSE]
    support <- rep(FALSE, length(local_nodes))
    if (!is.null(de_df) && nrow(de_df) > 0L) {
      de_match <- de_df[match(ann$gene_id, de_df$gene_id), , drop = FALSE]
      support <- deg_mask(de_match, alpha = alpha, lfc_cutoff = lfc_cut)
    }
    plot_df <- data.frame(
      evidence = c("Supported by current DE result", "Mainly prior/module context"),
      n = c(sum(support, na.rm = TRUE), sum(!support, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
    plotly::layout(
      plotly::plot_ly(plot_df, x = ~evidence, y = ~n, type = "bar", text = ~n, textposition = "auto"),
      title = "Displayed-gene evidence support",
      xaxis = list(title = ""),
      yaxis = list(title = "Gene count")
    )
  })
  output$semgraph_degree_scatter <- plotly::renderPlotly({
    res <- semgraph_res()
    shiny::req(!is.null(res))
    de_df <- tryCatch(de()$df, error = function(...) NULL)
    alpha <- input$de_alpha %||% 0.05
    lfc_cut <- input$de_lfc_cutoff %||% 0
    df <- res$degree_df
    df <- df[df$gene %in% select_semgraph_display_nodes(res, de_df, mode = input$semgraph_display_mode, max_nodes = input$semgraph_graph_max_nodes, alpha = alpha, lfc_cutoff = lfc_cut), , drop = FALSE]
    role <- ifelse(df$gene == res$focal_gene, "focal", ifelse(df$gene %in% res$upstream, "upstream", ifelse(df$gene %in% res$downstream, "downstream", "other")))
    ann <- res$local_annotations[match(df$gene, rownames(res$local_annotations)), c("gene_label", "full_gene_name", "is_druggable", "drug_mode", "drug_summary", "drug_source"), drop = FALSE]
    plot_df <- data.frame(df, role = role, gene_label = ifelse(is.na(ann$gene_label), df$gene, ann$gene_label), full_gene_name = ifelse(is.na(ann$full_gene_name), "Full gene name unavailable", ann$full_gene_name), is_druggable = ifelse(is.na(ann$is_druggable), FALSE, ann$is_druggable), drug_mode = ifelse(is.na(ann$drug_mode), "no interaction found", ann$drug_mode), drug_summary = ifelse(is.na(ann$drug_summary), "No external drug-target interaction found", ann$drug_summary), drug_source = ifelse(is.na(ann$drug_source), "None", ann$drug_source), stringsAsFactors = FALSE)
    plotly::layout(
      plotly::plot_ly(
        plot_df,
        x = ~in_degree,
        y = ~out_degree,
        color = ~role,
        symbol = ~is_druggable,
        symbols = c("circle", "diamond"),
        type = "scatter",
        mode = "markers+text",
        text = ~gene_label,
        textposition = "top center",
        hovertext = ~paste("Gene:", gene_label, "<br>Full gene name:", full_gene_name, "<br>In-degree:", in_degree, "<br>Out-degree:", out_degree, "<br>Total degree:", total_degree, "<br>Druggable:", ifelse(is_druggable, "yes", "no"), "<br>Drug mode:", drug_mode, "<br>Source(s):", drug_source, "<br>Example drugs/context:", drug_summary),
        hoverinfo = "text"
      ),
      title = "Regulator-responder geometry",
      xaxis = list(title = "In-degree"),
      yaxis = list(title = "Out-degree")
    )
  })
  output$semgraph_role_support_plot <- plotly::renderPlotly({
    res <- semgraph_res()
    shiny::req(!is.null(res))
    de_df <- tryCatch(de()$df, error = function(...) NULL)
    alpha <- input$de_alpha %||% 0.05
    lfc_cut <- input$de_lfc_cutoff %||% 0
    local_nodes <- select_semgraph_display_nodes(res, de_df, mode = input$semgraph_display_mode, max_nodes = input$semgraph_graph_max_nodes, alpha = alpha, lfc_cutoff = lfc_cut)
    ann <- res$local_annotations[match(local_nodes, rownames(res$local_annotations)), , drop = FALSE]
    role <- ifelse(local_nodes == res$focal_gene, "focal", ifelse(local_nodes %in% res$upstream, "upstream", ifelse(local_nodes %in% res$downstream, "downstream", "other")))
    support <- rep(FALSE, length(local_nodes))
    if (!is.null(de_df) && nrow(de_df) > 0L) {
      de_match <- de_df[match(ann$gene_id, de_df$gene_id), , drop = FALSE]
      support <- deg_mask(de_match, alpha = alpha, lfc_cutoff = lfc_cut)
    }
    plot_df <- as.data.frame(table(role = role, support = ifelse(support, "DE-supported", "Context-only")), stringsAsFactors = FALSE)
    plot_df <- plot_df[plot_df$role != "focal", , drop = FALSE]
    names(plot_df)[3] <- "n"
    role_label <- c(upstream = "Upstream candidates", downstream = "Downstream candidates", other = "Context genes")
    plot_df$role_label <- role_label[plot_df$role]
    plotly::layout(
      plotly::plot_ly(plot_df, x = ~role_label, y = ~n, color = ~support, type = "bar"),
      title = "Displayed causal roles by data support",
      xaxis = list(title = ""),
      yaxis = list(title = "Gene count"),
      barmode = "stack"
    )
  })
  output$semgraph_graph_help <- shiny::renderUI({
    res <- semgraph_res()
    shiny::req(!is.null(res))
    shiny::tagList(
      shiny::tags$p("This graph is a focused causal neighborhood, not the full fitted graph. The full graph can be too large to render clearly, so the app shows a readable subset centered on the current exploration mode."),
      shiny::tags$p("If you chose a combined program, the graph is fitted on the union of its source-module genes. That can reveal regulators spanning both arms of the biology, but it can also be broader and more heterogeneous than a single-module SEMgraph view."),
      shiny::tags$p("Important interpretation rule: this SEMgraph view is not built from differential expression alone. It combines your module-derived relationships with graph-orientation assumptions and prior-knowledge context. Differential expression is shown separately as a data-support layer."),
      shiny::tags$p("Display rule: the graph and tables are filtered to prioritize genes supported by your current DESeq2 result. Additional non-DE genes are only kept when they are needed to connect the DE-supported genes into a readable local causal neighborhood."),
      shiny::tags$p("New scope rule: the graph and the tables are no longer forced to use the same node set. The graph uses a readability limit, while the tables use a separate evidence filter so they can stay more comprehensive."),
      shiny::tags$p("Display modes: `Strict DE-only` shows only DE-supported genes plus the focal regulator. `DE + minimal connectors` adds the fewest context genes needed to connect the DE-supported neighborhood. `Full local neighborhood` shows the broader SEMgraph local context."),
      shiny::tags$p("Edge-validation rule: `STRING only` keeps edges with a STRING confidence above the chosen cutoff. `Any STRING support` keeps all matched STRING edges regardless of confidence threshold."),
      shiny::tags$p("Filter-versus-color rule: use the same chosen criteria either to filter the currently visible graph or to recolor the currently visible graph. Filtering changes which nodes and edges remain visible. Coloring does not remove anything from that current view; it only changes the appearance of the matching nodes and edges."),
      shiny::tags$p("Graph exploration modes: `Overview` gives the best high-level summary. `Neighborhood by hops` starts from one chosen gene and expands stepwise through nearby edges. `Directed paths` highlights the top short directed routes between a chosen start gene and end gene."),
      shiny::tags$p("Focused-graph rule: use the anchor-gene selector plus the biological category rows to define one focus set, then press either `Apply Filter` or `Apply Color`. A common workflow is to filter first to a smaller mechanism and then recolor within that filtered view."),
      shiny::tags$p("Graph size control: increase the maximum graph nodes when you want a broader zoomable overview. If the graph becomes hard to read, keep the graph smaller and use the broader tables to inspect additional relevant genes."),
      shiny::tags$p("Visual legend: node fill color shows differential expression log2 fold change, with blue for downregulation, white for near-zero change, and red for upregulation when a DESeq2 result is available."),
      shiny::tags$p("Node opacity shows how much the current DESeq2 result supports a node. More opaque nodes have stronger evidence from your transcriptomics contrast; faded nodes are mainly context from the fitted causal/prior graph."),
      shiny::tags$p("Border color meaning: dark red = focal master-regulator candidate, blue = upstream candidate, green = downstream candidate, gray = nearby context gene."),
      shiny::tags$p("Shape/style meaning: star = focal master-regulator candidate, diamond = transcription factor from GO.db/org.Hs.eg.db, circle = other genes. A thicker border marks a WGCNA hub-like gene in the selected module."),
      shiny::tags$p("Edge color meaning: gray edges have no STRING support under the current view. Blue edges have STRING pair support, but the STRING score is below the current confidence threshold. Red edges have STRING pair support and also pass the current confidence threshold. Dark slate edges are path-navigation aids shown when you use the directed-path exploration mode."),
      shiny::tags$p("How to interpret red versus blue: both red and blue edges are supported by STRING. The difference is support strength, not direction. Blue means lower-confidence STRING support. Red means higher-confidence STRING support according to the threshold slider. SEMgraph still provides the direction hypothesis in both cases."),
      shiny::tags$p("Arrow color meaning: the arrowhead always belongs to the same edge line and uses the same color code as the edge. So a red arrowhead means a higher-confidence STRING-supported SEMgraph edge, a blue arrowhead means a lower-confidence STRING-supported SEMgraph edge, a gray arrowhead means no STRING support, and a dark slate arrowhead marks a path-navigation edge."),
      shiny::tags$p("How to change edge and arrow colors: adjust the STRING confidence slider to move edges between blue and red. Use `Apply Filter` if you want the graph to keep only edges matching the selected STRING rule. Use `Apply Color` if you want to keep the same graph structure and only strengthen the color of edges that match that rule while the other visible edges stay in place with their usual gray/blue/red meaning. The arrowheads update automatically with the same color logic."),
      shiny::tags$p("Visible-edge support table: just under the graph there is now a dedicated edge table for the currently displayed graph. Use it when you want to see which visible edges have STRING pair support, which have STRING action annotations, and whether any directed STRING action agrees with the SEMgraph arrow."),
      shiny::tags$p("Drug-target clue: druggable genes now use a brighter gold glow plus a double gold outer ring so they remain visible even in dense graphs. If you also apply category colors, those colors can tint the text labels without removing the gold druggability marker. Hover over the node to see which sources contributed the annotation, example drugs, and whether a source described the interaction as inhibiting, activating, mixed, or unclear."),
      shiny::tags$p("Arrow meaning: the arrow still points from putative upstream to putative downstream in the fitted SEMgraph DAG, even when the edge is shown because STRING supports the pair. STRING supports the pairwise biological plausibility; SEMgraph still provides the direction hypothesis."),
      shiny::tags$p("How direction was inferred: the app first builds a prior graph from the selected module's strongest gene-gene relationships, then converts that prior graph into a directed acyclic graph and fits SEMgraph with topological ordering. So the arrow direction is a model-based orientation under those assumptions, not direct experimental proof."),
      shiny::tags$p("Drug-target provenance: DGIdb and HCDT provide external drug-target knowledge, while Reactome provides drug-mediated pathway context. None of these sources prove drug action in your samples. Use them as prioritization layers on top of your DE, WGCNA, and SEMgraph results."),
      shiny::tags$p("Beginner interpretation: prioritize genes that combine several evidence layers at once: strong DE support, meaningful causal position, WGCNA hub behavior, and biological plausibility from TF annotation, pathway results, or drug-target tractability. Do not over-interpret faded non-DE nodes as transcriptomics discoveries by themselves."),
      shiny::tags$p("Why this is separate from the graph summary tab: the summary tells you what model was fit and how much of the display is supported by your data; the graph tab focuses on directional structure and visual interpretation.")
    )
  })
  build_semgraph_graph_plot <- function() {
    res <- semgraph_res()
    shiny::req(!is.null(res))
    de_df <- tryCatch(de()$df, error = function(...) NULL)
    alpha <- input$de_alpha %||% 0.05
    lfc_cut <- input$de_lfc_cutoff %||% 0
    filter_state <- semgraph_filter_state()
    color_targets <- semgraph_color_targets()
    focus_apply_mode <- color_targets$apply_mode %||% "filter"
    active_edge_rule <- if (!is.null(filter_state)) (filter_state$edge_rule %||% "all") else "all"
    color_edge_rule <- color_targets$edge_rule %||% "all"
    g <- res$graph
    selection <- semgraph_graph_selection()
    local_nodes <- selection$nodes
    local_nodes <- local_nodes[local_nodes %in% igraph::V(g)$name]
    shiny::req(length(local_nodes) >= 1L)
    g <- igraph::induced_subgraph(g, vids = local_nodes)
    lay <- compute_semgraph_layout(g, method = input$semgraph_layout)
    coords <- data.frame(x = lay[, 1], y = lay[, 2], gene = igraph::V(g)$name, stringsAsFactors = FALSE)
    coords$role <- ifelse(coords$gene == res$focal_gene, "focal", ifelse(coords$gene %in% res$upstream, "upstream", ifelse(coords$gene %in% res$downstream, "downstream", "other")))
    local_ann <- res$local_annotations[match(coords$gene, rownames(res$local_annotations)), , drop = FALSE]
    coords$is_hub <- ifelse(is.na(local_ann$is_hub), FALSE, local_ann$is_hub)
    coords$is_tf <- ifelse(is.na(local_ann$is_tf), FALSE, local_ann$is_tf)
    coords$is_druggable <- ifelse(is.na(local_ann$is_druggable), FALSE, local_ann$is_druggable)
    coords$drug_mode <- ifelse(is.na(local_ann$drug_mode), "no interaction found", local_ann$drug_mode)
    coords$drug_summary <- ifelse(is.na(local_ann$drug_summary), "No external drug-target interaction found", local_ann$drug_summary)
    coords$drug_count <- ifelse(is.na(local_ann$drug_count), 0L, local_ann$drug_count)
    coords$drug_source <- ifelse(is.na(local_ann$drug_source), "None", local_ann$drug_source)
    coords$drug_source_details <- ifelse(is.na(local_ann$drug_source_details), "No DGIdb, HCDT, or Reactome drug context found", local_ann$drug_source_details)
    coords$reactome_context <- ifelse(is.na(local_ann$reactome_context), "No Reactome drug-mediated pathway context found", local_ann$reactome_context)
    coords$gene_label <- ifelse(is.na(local_ann$gene_label), coords$gene, local_ann$gene_label)
    coords$full_gene_name <- ifelse(is.na(local_ann$full_gene_name), "Full gene name unavailable", local_ann$full_gene_name)
    coords$tf_source <- ifelse(is.na(local_ann$tf_source), "No TF annotation found", local_ann$tf_source)
    coords$symbol <- ifelse(coords$role == "focal", "star", ifelse(coords$is_tf, "diamond", "circle"))
    coords$border_color <- ifelse(
      coords$role == "focal", "#8e1b10",
      ifelse(coords$role == "upstream", "#1f77b4",
      ifelse(coords$role == "downstream", "#2ca02c", "#6c757d"))
    )
    coords$line_width <- ifelse(coords$is_hub, 3, 1)
    coords$line_color <- coords$border_color
    de_df <- tryCatch(de()$df, error = function(...) NULL)
    coords$log2fc <- 0
    coords$de_source <- "No DESeq2 result available"
    if (!is.null(de_df) && nrow(de_df) > 0L) {
      de_match <- de_df[match(local_ann$gene_id, de_df$gene_id), , drop = FALSE]
      coords$log2fc <- ifelse(is.na(de_match$log2FoldChange), 0, de_match$log2FoldChange)
      coords$de_source <- ifelse(is.na(de_match$log2FoldChange), "No matching DE value in current result", "Current DESeq2 contrast log2 fold change")
      coords$padj <- de_match$padj
      coords$de_supported <- deg_mask(de_match, alpha = alpha, lfc_cutoff = lfc_cut)
    }
    if (!"padj" %in% names(coords)) coords$padj <- NA_real_
    if (!"de_supported" %in% names(coords)) coords$de_supported <- FALSE
    fc_lim <- max(abs(coords$log2fc), na.rm = TRUE)
    if (!is.finite(fc_lim) || fc_lim == 0) {
      fc_lim <- 1
    }
    deg_df <- res$degree_df[match(coords$gene, res$degree_df$gene), , drop = FALSE]
    coords$node_size <- 12 + 18 * (deg_df$total_degree - min(deg_df$total_degree, na.rm = TRUE)) / max(1, diff(range(deg_df$total_degree, na.rm = TRUE)))
    coords$drug_ring_size <- coords$node_size + 8
    coords$is_focus_colored <- coords$gene %in% (color_targets$nodes %||% character())
    coords$focus_ring_size <- coords$node_size + 14
    coords$text_color <- "#404040"
    if (identical(focus_apply_mode, "color") && length(color_targets$categories)) {
      for (nm in names(color_targets$categories)) {
        cat_info <- color_targets$categories[[nm]]
        if (is.null(cat_info) || !length(cat_info$nodes)) next
        idx <- coords$gene %in% cat_info$nodes
        coords$text_color[idx] <- cat_info$color %||% "#404040"
      }
    }
    coords$opacity <- ifelse(coords$de_supported, 0.95, 0.28)
    coords$evidence_layer <- ifelse(coords$de_supported, "Supported by current DE result", "Mainly prior/module context in current contrast")
    hub_quantile <- switch(
      "off",
      NA_real_
    )
    coords$is_central_hub <- FALSE
    if (is.finite(hub_quantile)) {
      cutoff <- stats::quantile(deg_df$total_degree, probs = hub_quantile, na.rm = TRUE)
      coords$is_central_hub <- deg_df$total_degree >= cutoff
    }
    coords$is_highlighted <- FALSE
    edge_df <- annotate_semgraph_edges(g, res)
    if (!identical(active_edge_rule, "all")) {
      edge_df <- filter_validated_edges(edge_df, mode = active_edge_rule, string_cutoff = input$semgraph_string_cutoff %||% 0.7)
    }
    if (nrow(edge_df) > 0L) {
      keep_nodes <- unique(c(edge_df$source, edge_df$target))
      coords <- coords[coords$gene %in% keep_nodes, , drop = FALSE]
      edge_df$x <- coords$x[match(edge_df$source, coords$gene)]
      edge_df$y <- coords$y[match(edge_df$source, coords$gene)]
      edge_df$xend <- coords$x[match(edge_df$target, coords$gene)]
      edge_df$yend <- coords$y[match(edge_df$target, coords$gene)]
      if (identical(selection$mode, "paths") && length(selection$paths)) {
        path_keys <- unique(unlist(lapply(selection$paths, function(path) {
          if (length(path) < 2L) return(character())
          paste(utils::head(path, -1L), utils::tail(path, -1L), sep = "||")
        }), use.names = FALSE))
        edge_df$is_path_edge <- paste(edge_df$source, edge_df$target, sep = "||") %in% path_keys
      } else {
        edge_df$is_path_edge <- FALSE
      }
    } else if (!identical(active_edge_rule, "all")) {
      coords <- coords[0, , drop = FALSE]
    }
    standard_edge_color <- function(row) {
      string_above_cutoff <- isTRUE(row$string_supported) && isTRUE(row$string_score >= (input$semgraph_string_cutoff %||% 0.4))
      if (isTRUE(row$is_path_edge)) {
        "rgba(52,73,94,0.9)"
      } else if (string_above_cutoff) {
        "rgba(200,30,30,0.75)"
      } else if (isTRUE(row$string_supported)) {
        "rgba(45,95,160,0.55)"
      } else {
        "rgba(130,130,130,0.25)"
      }
    }
    standard_edge_width <- function(row) {
      string_above_cutoff <- isTRUE(row$string_supported) && isTRUE(row$string_score >= (input$semgraph_string_cutoff %||% 0.4))
      if (isTRUE(row$is_path_edge)) 2.8 else if (string_above_cutoff) 1.9 else if (isTRUE(row$string_supported)) 1.6 else 1
    }
    standard_arrow_color <- function(row) {
      string_above_cutoff <- isTRUE(row$string_supported) && isTRUE(row$string_score >= (input$semgraph_string_cutoff %||% 0.4))
      if (isTRUE(row$is_path_edge)) "rgba(52,73,94,0.95)" else if (string_above_cutoff) "rgba(200,30,30,0.95)" else if (isTRUE(row$string_supported)) "rgba(45,95,160,0.85)" else "rgba(100,100,100,0.65)"
    }
    standard_arrow_size <- function(row) {
      string_above_cutoff <- isTRUE(row$string_supported) && isTRUE(row$string_score >= (input$semgraph_string_cutoff %||% 0.4))
      if (isTRUE(row$is_path_edge)) 1.1 else if (string_above_cutoff) 1 else if (isTRUE(row$string_supported)) 0.95 else 0.8
    }
    standard_arrow_width <- function(row) {
      string_above_cutoff <- isTRUE(row$string_supported) && isTRUE(row$string_score >= (input$semgraph_string_cutoff %||% 0.4))
      if (isTRUE(row$is_path_edge)) 2 else if (string_above_cutoff) 1.8 else if (isTRUE(row$string_supported)) 1.5 else 1
    }
    p <- plotly::plot_ly(source = "semgraph_graph")
    if (nrow(edge_df) > 0L) {
      cutoff_now <- input$semgraph_string_cutoff %||% 0.4
      for (i in seq_len(nrow(edge_df))) {
        string_above_cutoff <- isTRUE(edge_df$string_supported[i]) && isTRUE(edge_df$string_score[i] >= cutoff_now)
        edge_matches_color_rule <- if (identical(color_edge_rule, "string_only")) {
          string_above_cutoff
        } else if (identical(color_edge_rule, "any_support")) {
          isTRUE(edge_df$string_supported[i])
        } else {
          FALSE
        }
        cur_row <- edge_df[i, , drop = FALSE]
        edge_col <- standard_edge_color(cur_row)
        edge_w <- standard_edge_width(cur_row)
        if (identical(focus_apply_mode, "color") && !identical(color_edge_rule, "all") && edge_matches_color_rule) {
          edge_col <- if (identical(color_edge_rule, "string_only")) {
            "rgba(200,30,30,0.95)"
          } else if (string_above_cutoff) {
            "rgba(200,30,30,0.95)"
          } else {
            "rgba(45,95,160,0.95)"
          }
          edge_w <- edge_w + 0.9
        }
        p <- plotly::add_segments(
          p,
          x = edge_df$x[i], y = edge_df$y[i],
          xend = edge_df$xend[i], yend = edge_df$yend[i],
          inherit = FALSE,
          line = list(
            color = edge_col,
            width = edge_w
          ),
          showlegend = FALSE
        )
      }
    }
    if (identical(focus_apply_mode, "color") && length(color_targets$categories)) {
      ring_layers <- list(
        regulator = 16,
        hub = 12,
        tf = 8,
        anchor = 4,
        druggable = 0
      )
      for (nm in names(ring_layers)) {
        cat_info <- color_targets$categories[[nm]]
        if (is.null(cat_info) || !length(cat_info$nodes)) next
        idx <- coords$gene %in% cat_info$nodes
        if (!any(idx)) next
        ring_size <- coords$node_size[idx] + ring_layers[[nm]]
        ring_width <- if (identical(nm, "druggable")) 4 else 3
        p <- plotly::add_markers(
          p,
          data = coords[idx, , drop = FALSE],
          x = ~x,
          y = ~y,
          inherit = FALSE,
          hoverinfo = "skip",
          marker = list(
            size = ring_size,
            color = "rgba(0,0,0,0)",
            line = list(color = cat_info$color %||% "#0f9d8a", width = ring_width),
            symbol = "circle-open"
          ),
          showlegend = FALSE
        )
      }
    }
    if (any(coords$is_druggable)) {
      p <- plotly::add_markers(
        p,
        data = coords[coords$is_druggable, , drop = FALSE],
        x = ~x,
        y = ~y,
        inherit = FALSE,
        hoverinfo = "skip",
        marker = list(
          size = coords$drug_ring_size[coords$is_druggable] + 14,
          color = "rgba(255,210,74,0.58)",
          line = list(color = "rgba(0,0,0,0)", width = 0),
          symbol = "circle"
        ),
        showlegend = FALSE
      )
    }
    if (any(coords$is_druggable) && !identical(focus_apply_mode, "color")) {
      p <- plotly::add_markers(
        p,
        data = coords[coords$is_druggable, , drop = FALSE],
        x = ~x,
        y = ~y,
        inherit = FALSE,
        hoverinfo = "skip",
        marker = list(
          size = coords$drug_ring_size[coords$is_druggable] + 7,
          color = "rgba(0,0,0,0)",
          line = list(color = "#a66b00", width = 6),
          symbol = "circle-open"
        ),
        showlegend = FALSE
      )
      p <- plotly::add_markers(
        p,
        data = coords[coords$is_druggable, , drop = FALSE],
        x = ~x,
        y = ~y,
        inherit = FALSE,
        hoverinfo = "skip",
        marker = list(
          size = coords$drug_ring_size[coords$is_druggable] + 1,
          color = "rgba(0,0,0,0)",
          line = list(color = "#ffd24a", width = 3),
          symbol = "circle-open"
        ),
        showlegend = FALSE
      )
    }
    p <- plotly::add_trace(
      p,
      data = coords,
      x = ~x,
      y = ~y,
      type = "scatter",
      mode = "markers",
      hovertext = ~paste(
        "Gene:", gene_label,
        "<br>Full gene name:", full_gene_name,
        "<br>Role:", role,
        "<br>log2 fold change:", signif(log2fc, 3),
        "<br>Adjusted p-value:", ifelse(is.na(padj), "NA", signif(padj, 3)),
        "<br>Data-support layer:", evidence_layer,
        "<br>DE source:", de_source,
        "<br>WGCNA hub:", ifelse(is_hub, "yes", "no"),
        "<br>Hub source: WGCNA within-module connectivity (top local decile)",
        "<br>Transcription factor:", ifelse(is_tf, "yes", "no"),
        "<br>TF source:", tf_source,
        "<br>Druggable target:", ifelse(is_druggable, "yes", "no"),
        "<br>Drug interaction mode:", drug_mode,
        "<br>Example drugs:", drug_summary,
        "<br>Drug-target source(s):", drug_source,
        "<br>Source detail:", drug_source_details,
        "<br>Reactome context:", reactome_context,
        "<br>SEMgraph central hub by selected cutoff:", ifelse(is_central_hub, "yes", "no"),
        "<br>Matches current focus-color criteria:", ifelse(is_focus_colored, "yes", "no"),
        "<br>Included in current focused graph:", "yes"
      ),
      hoverinfo = "text",
      marker = list(
        size = coords$node_size,
        color = coords$log2fc,
        colorscale = "RdBu",
        cmin = -fc_lim,
        cmax = fc_lim,
        opacity = coords$opacity,
        colorbar = list(title = "log2FC"),
        symbol = coords$symbol,
        line = list(color = coords$line_color, width = coords$line_width)
      ),
      showlegend = FALSE
    )
    if (nrow(coords) > 0L) {
      text_groups <- split(coords, coords$text_color)
      for (col_name in names(text_groups)) {
        grp <- text_groups[[col_name]]
        if (is.null(grp) || !nrow(grp)) next
        p <- plotly::add_trace(
          p,
          data = grp,
          x = ~x,
          y = ~y,
          type = "scatter",
          mode = "text",
          text = ~gene,
          textposition = "top center",
          hoverinfo = "skip",
          textfont = list(size = 10, color = col_name),
          showlegend = FALSE,
          inherit = FALSE
        )
      }
    }
    layout_revision <- paste(
      input$semgraph_layout %||% "fr",
      paste(sort(coords$gene), collapse = "||"),
      input$semgraph_graph_scope %||% "overview",
      sep = "::"
    )
    plotly::config(plotly::layout(
      p,
      title = paste(
        "SEMgraph causal explorer:",
        res$analysis_unit,
        if (identical(selection$mode, "neighborhood")) sprintf("| %s within %s hops", selection$start_gene, input$semgraph_neighbor_hops %||% 2L)
        else if (identical(selection$mode, "paths")) sprintf("| %s to %s", selection$start_gene, selection$end_gene)
        else ""
      ),
      xaxis = list(title = "", showticklabels = FALSE, zeroline = FALSE),
      yaxis = list(title = "", showticklabels = FALSE, zeroline = FALSE),
      showlegend = FALSE,
      dragmode = "pan",
      uirevision = layout_revision,
      annotations = if (nrow(edge_df) > 0L) {
        lapply(seq_len(nrow(edge_df)), function(i) {
          dx <- edge_df$xend[i] - edge_df$x[i]
          dy <- edge_df$yend[i] - edge_df$y[i]
          frac <- 0.18
          string_above_cutoff <- isTRUE(edge_df$string_supported[i]) && isTRUE(edge_df$string_score[i] >= (input$semgraph_string_cutoff %||% 0.4))
          edge_matches_color_rule <- if (identical(color_edge_rule, "string_only")) {
            string_above_cutoff
          } else if (identical(color_edge_rule, "any_support")) {
            isTRUE(edge_df$string_supported[i])
          } else {
            FALSE
          }
          list(
            x = edge_df$xend[i],
            y = edge_df$yend[i],
            ax = edge_df$xend[i] - dx * frac,
            ay = edge_df$yend[i] - dy * frac,
            xref = "x",
            yref = "y",
            axref = "x",
            ayref = "y",
            text = "",
            showarrow = TRUE,
            arrowhead = 2,
            arrowsize = {
              base_size <- standard_arrow_size(edge_df[i, , drop = FALSE])
              if (identical(focus_apply_mode, "color") && !identical(color_edge_rule, "all") && edge_matches_color_rule) base_size + 0.12 else base_size
            },
            arrowwidth = {
              base_width <- standard_arrow_width(edge_df[i, , drop = FALSE])
              if (identical(focus_apply_mode, "color") && !identical(color_edge_rule, "all") && edge_matches_color_rule) base_width + 0.4 else base_width
            },
            arrowcolor = {
              base_col <- standard_arrow_color(edge_df[i, , drop = FALSE])
              if (identical(focus_apply_mode, "color") && !identical(color_edge_rule, "all") && edge_matches_color_rule) {
                if (identical(color_edge_rule, "string_only")) {
                  "rgba(200,30,30,0.98)"
                } else if (string_above_cutoff) {
                  "rgba(200,30,30,0.98)"
                } else {
                  "rgba(45,95,160,0.98)"
                }
              } else {
                base_col
              }
            }
          )
        })
      } else {
        list()
      }
    ), displaylogo = FALSE, responsive = TRUE, scrollZoom = TRUE)
  }
  output$semgraph_graph_plot <- plotly::renderPlotly({
    build_semgraph_graph_plot()
  })
  output$semgraph_graph_plot_large <- plotly::renderPlotly({
    build_semgraph_graph_plot()
  })
  output$download_semgraph_graphml <- shiny::downloadHandler(
    filename = function() sprintf("semgraph_%s.graphml", gsub("[^A-Za-z0-9_]+", "_", semgraph_res()$analysis_unit %||% "subgraph")),
    content = function(file) {
      res <- semgraph_res()
      de_df <- tryCatch(de()$df, error = function(...) NULL)
      alpha <- input$de_alpha %||% 0.05
      lfc_cut <- input$de_lfc_cutoff %||% 0
      focus_criteria <- semgraph_filter_state() %||% list()
      focus_mode <- if (!is.null(semgraph_filter_state())) "filter" else "none"
      active_edge_rule <- if (identical(focus_mode, "filter")) (focus_criteria$edge_rule %||% "all") else "all"
      sel <- select_semgraph_graph_nodes(
        res,
        de_df = de_df,
        display_mode = input$semgraph_display_mode,
        graph_scope = input$semgraph_graph_scope,
        max_nodes = input$semgraph_graph_max_nodes,
        start_gene = input$semgraph_start_gene,
        end_gene = input$semgraph_end_gene,
        hops = input$semgraph_neighbor_hops %||% 2L,
        path_n = input$semgraph_path_n %||% 5L,
        path_max_steps = input$semgraph_path_max_steps %||% 6L,
        selected_nodes = if (identical(focus_mode, "filter")) (focus_criteria$selected_nodes %||% character()) else character(),
        selected_hops = if (identical(focus_mode, "filter")) (focus_criteria$selected_hops %||% 3L) else 3L,
        selected_only = isTRUE(focus_criteria$selected_only) && identical(focus_mode, "filter"),
        filter_tf = isTRUE(focus_criteria$filter_tf) && identical(focus_mode, "filter"),
        filter_druggable = isTRUE(focus_criteria$filter_druggable) && identical(focus_mode, "filter"),
        filter_hub_mode = if (identical(focus_mode, "filter")) (focus_criteria$filter_hub_mode %||% "off") else "off",
        focal_top_n = if (identical(focus_mode, "filter")) (focus_criteria$focal_top_n %||% 0L) else 0L,
        restrict_to_focus = isTRUE(focus_criteria$restrict_to_focus) && identical(focus_mode, "filter"),
        alpha = alpha,
        lfc_cutoff = lfc_cut
      )
      sub_nodes <- sel$nodes[sel$nodes %in% igraph::V(res$graph)$name]
      if (!length(sub_nodes)) stop("No SEMgraph subgraph available for export.")
      g_sub <- igraph::induced_subgraph(res$graph, vids = sub_nodes)
      edge_tbl <- annotate_semgraph_edges(g_sub, res)
      if (!identical(active_edge_rule, "all")) {
        edge_tbl <- filter_validated_edges(edge_tbl, mode = active_edge_rule, string_cutoff = input$semgraph_string_cutoff %||% 0.7)
      }
      keep_edges <- integer()
      if (nrow(edge_tbl)) {
        keep_edges <- unname(vapply(seq_len(nrow(edge_tbl)), function(i) {
          eid <- igraph::get.edge.ids(g_sub, vp = c(edge_tbl$source[i], edge_tbl$target[i]), directed = TRUE, error = FALSE)
          if (length(eid) && eid[1L] > 0) eid[1L] else NA_integer_
        }, integer(1)))
        keep_edges <- keep_edges[is.finite(keep_edges)]
      }
      if (length(keep_edges)) {
        g_sub <- igraph::subgraph.edges(g_sub, eids = keep_edges, delete.vertices = TRUE)
        edge_tbl <- annotate_semgraph_edges(g_sub, res)
        if (nrow(edge_tbl)) {
          edge_map <- paste(edge_tbl$source, edge_tbl$target, sep = "||")
          cur_map <- paste(igraph::as_data_frame(g_sub, what = "edges")$from, igraph::as_data_frame(g_sub, what = "edges")$to, sep = "||")
          for (col in setdiff(names(edge_tbl), c("source", "target"))) {
            igraph::edge_attr(g_sub, col) <- edge_tbl[[col]][match(cur_map, edge_map)]
          }
        }
      } else {
        g_sub <- igraph::delete_edges(g_sub, igraph::E(g_sub))
      }
      igraph::write_graph(g_sub, file, format = "graphml")
    }
  )
  shiny::observeEvent(input$semgraph_open_large, {
    shiny::showModal(
      shiny::modalDialog(
        title = "Larger SEMgraph explorer",
        size = "l",
        easyClose = TRUE,
        shiny::tags$p(class = "muted", "Use this larger view when the standard graph is too dense. Zoom, pan, and hover to inspect genes and edges more comfortably."),
        plotly::plotlyOutput("semgraph_graph_plot_large", height = "860px")
      )
    )
  })
  output$semgraph_graph_analysis_table <- DT::renderDataTable({
    res <- semgraph_res()
    shiny::req(!is.null(res))
    sel <- semgraph_graph_selection()
    de_df <- tryCatch(de()$df, error = function(...) NULL)
    ranked <- rank_semgraph_nodes(res, de_df, alpha = input$de_alpha %||% 0.05, lfc_cutoff = input$de_lfc_cutoff %||% 0)
    validation_tbl <- summarize_node_validation(
      filter_validated_edges(
        annotate_semgraph_edges(res$graph, res),
        mode = current_semgraph_filter_rule(),
        string_cutoff = input$semgraph_string_cutoff %||% 0.7
      )
    )
    display_nodes <- sel$nodes
    if (!identical(current_semgraph_filter_rule(), "all")) {
      g_sub <- igraph::induced_subgraph(res$graph, vids = intersect(sel$nodes, igraph::V(res$graph)$name))
      edge_tbl <- annotate_semgraph_edges(g_sub, res)
      edge_tbl <- filter_validated_edges(edge_tbl, mode = current_semgraph_filter_rule(), string_cutoff = input$semgraph_string_cutoff %||% 0.7)
      if (nrow(edge_tbl)) {
        display_nodes <- unique(c(edge_tbl$source, edge_tbl$target))
      }
    }
    ranked <- ranked[ranked$gene %in% display_nodes, , drop = FALSE]
    ranked$role <- ifelse(ranked$gene == res$focal_gene, "focal", ifelse(ranked$gene %in% res$upstream, "upstream", ifelse(ranked$gene %in% res$downstream, "downstream", "other")))
    ann <- res$local_annotations[match(ranked$gene, rownames(res$local_annotations)), c("gene_label", "gene_id", "is_hub", "is_tf", "is_druggable", "drug_summary"), drop = FALSE]
    ranked <- merge(ranked, validation_tbl, by = "gene", all.x = TRUE, sort = FALSE)
    ranked$gene_label <- ifelse(is.na(ann$gene_label), ranked$gene, ann$gene_label)
    ranked$is_hub <- ifelse(is.na(ann$is_hub), FALSE, ann$is_hub)
    ranked$is_tf <- ifelse(is.na(ann$is_tf), FALSE, ann$is_tf)
    ranked$is_druggable <- ifelse(is.na(ann$is_druggable), FALSE, ann$is_druggable)
    ranked$drug_summary <- ifelse(is.na(ann$drug_summary), "No external drug-target interaction found", ann$drug_summary)
    ranked$externally_supported_edges <- ifelse(is.na(ranked$externally_supported_edges), 0, ranked$externally_supported_edges)
    ranked$string_supported_edges <- ifelse(is.na(ranked$string_supported_edges), 0, ranked$string_supported_edges)
    ranked$best_string_score <- ifelse(is.na(ranked$best_string_score), 0, ranked$best_string_score)
    if (!is.null(de_df) && nrow(de_df) > 0L) {
      de_match <- de_df[match(ann$gene_id, de_df$gene_id), , drop = FALSE]
      ranked$log2FoldChange <- de_match$log2FoldChange
      ranked$padj <- de_match$padj
    } else {
      ranked$log2FoldChange <- NA_real_
      ranked$padj <- NA_real_
    }
    if (identical(sel$mode, "paths") && length(sel$paths)) {
      path_rows <- do.call(rbind, lapply(seq_along(sel$paths), function(i) {
        path <- sel$paths[[i]]
        data.frame(
          view = "path",
          item = paste0("Path ", i),
          details = paste(path, collapse = " -> "),
          length = length(path) - 1L,
          stringsAsFactors = FALSE
        )
      }))
      make_downloadable_table(path_rows, page_length = 8, rownames = FALSE)
    } else {
      ranked$role <- ifelse(ranked$role == "focal", "focal anchor", ifelse(ranked$role == "upstream", "upstream candidate", ifelse(ranked$role == "downstream", "downstream candidate", "context gene")))
      ranked <- ranked[order(ranked$de_supported, ranked$externally_supported_edges, ranked$total_degree, decreasing = TRUE), c("gene_label", "role", "de_supported", "log2FoldChange", "padj", "is_hub", "is_tf", "is_druggable", "externally_supported_edges", "string_supported_edges", "best_string_score", "distance_to_focal", "total_degree", "drug_summary"), drop = FALSE]
      make_downloadable_table(ranked, page_length = 10, rownames = FALSE)
    }
  }, server = FALSE)
  output$semgraph_graph_edge_table <- DT::renderDataTable({
    edge_tbl <- visible_semgraph_edge_table()
    if (is.null(edge_tbl) || !nrow(edge_tbl)) {
      return(make_downloadable_table(
        data.frame(message = "No edges are visible under the current graph filter.", stringsAsFactors = FALSE),
        page_length = 5,
        rownames = FALSE
      ))
    }
    keep <- edge_tbl[, c(
      "source", "target", "validation_tier", "string_supported", "string_score",
      "string_direction_available", "string_direction_matches", "string_action_mode",
      "string_action_effect", "matches_current_color_rule", "string_source", "string_action_source"
    ), drop = FALSE]
    names(keep) <- c(
      "source_gene", "target_gene", "validation_tier", "string_supported", "string_score",
      "string_direction_available", "string_direction_matches_semgraph", "string_action_mode",
      "string_action_effect", "matches_current_color_rule", "string_pair_source", "string_action_source"
    )
    make_downloadable_table(keep, page_length = 10, rownames = FALSE)
  }, server = FALSE)
  output$semgraph_effect_proxy_plot <- plotly::renderPlotly({
    res <- semgraph_res()
    shiny::req(!is.null(res))
    de_df <- tryCatch(de()$df, error = function(...) NULL)
    df <- res$degree_df
    df <- df[df$gene %in% select_semgraph_table_nodes(res, de_df, scope = input$semgraph_table_scope, max_rows = input$semgraph_table_max_rows, alpha = input$de_alpha %||% 0.05, lfc_cutoff = input$de_lfc_cutoff %||% 0), , drop = FALSE]
    shiny::req(nrow(df) > 0L)
    df$role <- ifelse(df$gene == res$focal_gene, "focal", ifelse(df$gene %in% res$upstream, "upstream", ifelse(df$gene %in% res$downstream, "downstream", "other")))
    df$net_flow <- df$out_degree - df$in_degree
    ann <- res$local_annotations[match(df$gene, rownames(res$local_annotations)), c("gene_label", "full_gene_name", "is_druggable", "drug_mode", "drug_summary", "drug_source"), drop = FALSE]
    df$gene_label <- ifelse(is.na(ann$gene_label), df$gene, ann$gene_label)
    df$full_gene_name <- ifelse(is.na(ann$full_gene_name), "Full gene name unavailable", ann$full_gene_name)
    df$is_druggable <- ifelse(is.na(ann$is_druggable), FALSE, ann$is_druggable)
    df$drug_mode <- ifelse(is.na(ann$drug_mode), "no interaction found", ann$drug_mode)
    df$drug_summary <- ifelse(is.na(ann$drug_summary), "No external drug-target interaction found", ann$drug_summary)
    df$drug_source <- ifelse(is.na(ann$drug_source), "None", ann$drug_source)
    df <- df[order(df$net_flow, decreasing = TRUE), , drop = FALSE]
    top_df <- head(df, 20L)
    plotly::layout(
      plotly::plot_ly(
        top_df,
        x = ~reorder(gene_label, net_flow),
        y = ~net_flow,
        color = ~role,
        type = "bar",
        text = ~paste("Gene:", gene_label, "<br>Full gene name:", full_gene_name, "<br>Out-degree:", out_degree, "<br>In-degree:", in_degree, "<br>Net flow:", net_flow, "<br>Druggable:", ifelse(is_druggable, "yes", "no"), "<br>Drug mode:", drug_mode, "<br>Source(s):", drug_source, "<br>Example drugs/context:", drug_summary),
        hoverinfo = "text"
      ),
      title = "Local regulator priority view",
      xaxis = list(title = "", tickangle = -45),
      yaxis = list(title = "Out-degree minus in-degree")
    )
  })
  output$semgraph_regulator_help <- shiny::renderUI({
    shiny::tagList(
      shiny::tags$p("How to read this table: this is the main ranked regulator view for the current SEMgraph unit. It combines position in the directed graph with DE support, STRING support, and druggability so you can decide which genes deserve follow-up first."),
      shiny::tags$ul(
        shiny::tags$li("Higher out-degree suggests a gene may influence more downstream targets."),
        shiny::tags$li("Higher in-degree suggests a gene may integrate more upstream signals."),
        shiny::tags$li("Total degree is a simple importance score within the fitted causal graph."),
        shiny::tags$li("The DESeq2 columns in this table, such as log2 fold change and adjusted p-value, come from the current differential expression result, not from SEMgraph itself."),
        shiny::tags$li("External-validation columns summarize how many currently retained SEMgraph edges around a gene are supported by STRING. These are support layers for plausibility, not proof."),
        shiny::tags$li("Druggability columns come from multiple external resources. DGIdb and HCDT provide drug-target evidence, while Reactome contributes drug-mediated pathway context."),
        shiny::tags$li("Combine SEMgraph position with WGCNA hub status, DE support, pathway enrichment, and druggability before prioritizing a regulator.")
      )
    )
  })
  output$semgraph_degree_table <- DT::renderDataTable({
    cur <- ds()
    res <- semgraph_res()
    shiny::req(!is.null(res))
    de_df <- tryCatch(de()$df, error = function(...) NULL)
    df <- res$degree_df
    validation_tbl <- summarize_node_validation(
      filter_validated_edges(
        annotate_semgraph_edges(res$graph, res),
        mode = current_semgraph_filter_rule(),
        string_cutoff = input$semgraph_string_cutoff %||% 0.7
      )
    )
    map_df <- unique(cur$gene_map[, c("hgnc_symbol", "gene_label")])
    df <- merge(df, map_df, by.x = "gene", by.y = "hgnc_symbol", all.x = TRUE, sort = FALSE)
    ann <- res$local_annotations[, c("hgnc_symbol", "is_hub", "is_tf", "tf_source", "is_druggable", "drug_mode", "drug_summary", "drug_count", "drug_source", "drug_source_details", "reactome_context"), drop = FALSE]
    ann <- ann[!duplicated(ann$hgnc_symbol), , drop = FALSE]
    df <- merge(df, ann, by.x = "gene", by.y = "hgnc_symbol", all.x = TRUE, sort = FALSE)
    de_df <- tryCatch(de()$df, error = function(...) NULL)
    if (!is.null(de_df) && nrow(de_df) > 0L) {
      map_back <- unique(cur$gene_map[, c("hgnc_symbol", "gene_id")])
      df <- merge(df, map_back, by.x = "gene", by.y = "hgnc_symbol", all.x = TRUE, sort = FALSE)
      de_match <- de_df[match(df$gene_id, de_df$gene_id), , drop = FALSE]
      df$log2FoldChange <- de_match$log2FoldChange
      df$padj <- de_match$padj
      df$de_supported <- deg_mask(de_match, alpha = input$de_alpha %||% 0.05, lfc_cutoff = input$de_lfc_cutoff %||% 0)
    } else {
      df$log2FoldChange <- NA_real_
      df$padj <- NA_real_
      df$de_supported <- FALSE
      df$gene_id <- NA_character_
    }
    df <- merge(df, validation_tbl, by = "gene", all.x = TRUE, sort = FALSE)
    df <- df[df$gene %in% select_semgraph_table_nodes(res, de_df, scope = input$semgraph_table_scope, max_rows = input$semgraph_table_max_rows, alpha = input$de_alpha %||% 0.05, lfc_cutoff = input$de_lfc_cutoff %||% 0), , drop = FALSE]
    df$display_gene <- ifelse(is.na(df$gene_label), df$gene, df$gene_label)
    df$is_hub <- ifelse(is.na(df$is_hub), FALSE, df$is_hub)
    df$is_tf <- ifelse(is.na(df$is_tf), FALSE, df$is_tf)
    df$tf_source <- ifelse(is.na(df$tf_source), "No TF annotation found", df$tf_source)
    df$is_druggable <- ifelse(is.na(df$is_druggable), FALSE, df$is_druggable)
    df$drug_mode <- ifelse(is.na(df$drug_mode), "no interaction found", df$drug_mode)
    df$drug_summary <- ifelse(is.na(df$drug_summary), "No external drug-target interaction found", df$drug_summary)
    df$drug_count <- ifelse(is.na(df$drug_count), 0L, df$drug_count)
    df$drug_source <- ifelse(is.na(df$drug_source), "None", df$drug_source)
    df$drug_source_details <- ifelse(is.na(df$drug_source_details), "No DGIdb, HCDT, or Reactome drug context found", df$drug_source_details)
    df$reactome_context <- ifelse(is.na(df$reactome_context), "No Reactome drug-mediated pathway context found", df$reactome_context)
    df$externally_supported_edges <- ifelse(is.na(df$externally_supported_edges), 0, df$externally_supported_edges)
    df$string_supported_edges <- ifelse(is.na(df$string_supported_edges), 0, df$string_supported_edges)
    df$best_string_score <- ifelse(is.na(df$best_string_score), 0, df$best_string_score)
    df <- df[order(df$externally_supported_edges, df$total_degree, decreasing = TRUE), c("display_gene", "de_supported", "log2FoldChange", "padj", "is_hub", "is_tf", "is_druggable", "externally_supported_edges", "string_supported_edges", "best_string_score", "drug_mode", "drug_count", "drug_summary", "tf_source", "drug_source", "drug_source_details", "reactome_context", "out_degree", "in_degree", "total_degree"), drop = FALSE]
    names(df)[1] <- "gene"
    make_downloadable_table(df, page_length = 10)
  }, server = FALSE)
  output$semgraph_ace_table <- DT::renderDataTable({
    res <- semgraph_res()
    shiny::req(!is.null(res))
    de_df <- tryCatch(de()$df, error = function(...) NULL)
    df <- res$ace
    validation_tbl <- summarize_node_validation(
      filter_validated_edges(
        annotate_semgraph_edges(res$graph, res),
        mode = current_semgraph_filter_rule(),
        string_cutoff = input$semgraph_string_cutoff %||% 0.7
      )
    )
    ann <- res$local_annotations[, c("hgnc_symbol", "is_druggable", "drug_mode", "drug_summary", "drug_count", "drug_source", "drug_source_details", "reactome_context"), drop = FALSE]
    ann <- ann[!duplicated(ann$hgnc_symbol), , drop = FALSE]
    if (nrow(df) == 0L) {
      fallback <- res$degree_df
      fallback <- fallback[fallback$gene %in% select_semgraph_table_nodes(res, de_df, scope = input$semgraph_table_scope, max_rows = input$semgraph_table_max_rows, alpha = input$de_alpha %||% 0.05, lfc_cutoff = input$de_lfc_cutoff %||% 0), , drop = FALSE]
      fallback$regulator_bias_score <- fallback$out_degree - fallback$in_degree
      fallback <- merge(fallback, ann, by.x = "gene", by.y = "hgnc_symbol", all.x = TRUE, sort = FALSE)
      fallback <- merge(fallback, validation_tbl, by = "gene", all.x = TRUE, sort = FALSE)
      fallback <- fallback[order(fallback$regulator_bias_score, decreasing = TRUE), , drop = FALSE]
      fallback$is_druggable <- ifelse(is.na(fallback$is_druggable), FALSE, fallback$is_druggable)
      fallback$drug_mode <- ifelse(is.na(fallback$drug_mode), "no interaction found", fallback$drug_mode)
      fallback$drug_summary <- ifelse(is.na(fallback$drug_summary), "No external drug-target interaction found", fallback$drug_summary)
      fallback$drug_count <- ifelse(is.na(fallback$drug_count), 0L, fallback$drug_count)
      fallback$drug_source <- ifelse(is.na(fallback$drug_source), "None", fallback$drug_source)
      fallback$drug_source_details <- ifelse(is.na(fallback$drug_source_details), "No DGIdb, HCDT, or Reactome drug context found", fallback$drug_source_details)
      fallback$reactome_context <- ifelse(is.na(fallback$reactome_context), "No Reactome drug-mediated pathway context found", fallback$reactome_context)
      fallback$externally_supported_edges <- ifelse(is.na(fallback$externally_supported_edges), 0, fallback$externally_supported_edges)
      fallback$string_supported_edges <- ifelse(is.na(fallback$string_supported_edges), 0, fallback$string_supported_edges)
      fallback$best_string_score <- ifelse(is.na(fallback$best_string_score), 0, fallback$best_string_score)
      return(make_downloadable_table(fallback, page_length = 10))
    }
    merge_candidates <- intersect(c("source", "target", "from", "to", "gene"), names(df))
    if (length(merge_candidates) > 0L) {
      gene_col <- merge_candidates[1L]
      df <- merge(df, ann, by.x = gene_col, by.y = "hgnc_symbol", all.x = TRUE, sort = FALSE)
      df$is_druggable <- ifelse(is.na(df$is_druggable), FALSE, df$is_druggable)
      df$drug_mode <- ifelse(is.na(df$drug_mode), "no interaction found", df$drug_mode)
      df$drug_summary <- ifelse(is.na(df$drug_summary), "No external drug-target interaction found", df$drug_summary)
      df$drug_count <- ifelse(is.na(df$drug_count), 0L, df$drug_count)
      df$drug_source <- ifelse(is.na(df$drug_source), "None", df$drug_source)
      df$drug_source_details <- ifelse(is.na(df$drug_source_details), "No DGIdb, HCDT, or Reactome drug context found", df$drug_source_details)
      df$reactome_context <- ifelse(is.na(df$reactome_context), "No Reactome drug-mediated pathway context found", df$reactome_context)
    }
    make_downloadable_table(df, page_length = 10)
  }, server = FALSE)
  output$semgraph_ace_help <- shiny::renderUI({
    res <- semgraph_res()
    shiny::req(!is.null(res))
    shiny::tagList(
      shiny::tags$p(sprintf("ACE means Average Causal Effect. In plain language, it tries to estimate how much changing one node would be expected to influence another node inside the fitted causal model.")),
      shiny::tags$p(res$ace_note),
      shiny::tags$p("Why this is a separate tab: the graph tab is for directional structure, while this tab is for effect-size style interpretation or, when ACE is unavailable, a simpler local regulator-priority summary."),
      shiny::tags$p("If an ACE table is shown, its numbers come from SEMgraph's fitted causal model. If ACE is skipped, the fallback table/plot instead summarize local graph directionality such as out-degree minus in-degree, which is not the same thing as a causal effect estimate."),
      shiny::tags$p("Druggability is not an ACE output. Any drug names, activating/inhibiting labels, or Reactome drug-mediated pathway notes shown elsewhere come from external resources and should be used as separate pharmacology-prioritization layers."),
      shiny::tags$p("Use this tab after the graph and regulator tabs. It is best for deciding which already plausible regulators may also have stronger local effect-like influence inside the fitted causal model.")
    )
  })
  semgraph_pathway_context <- shiny::reactive({
    res <- semgraph_res()
    shiny::req(!is.null(res))
    de_df <- tryCatch(de()$df, error = function(...) NULL)
    build_semgraph_pathway_context(res, de_df = de_df, mode = input$semgraph_display_mode, alpha = input$de_alpha %||% 0.05, lfc_cutoff = input$de_lfc_cutoff %||% 0)
  })
  output$semgraph_pathway_help <- shiny::renderUI({
    shiny::tagList(
      shiny::tags$p("This section connects the displayed SEMgraph neighborhood to Reactome pathway context without overcrowding the main causal graph."),
      shiny::tags$ul(
        shiny::tags$li("Pathway summary plot: each point is a Reactome pathway overlapping the currently displayed causal neighborhood."),
        shiny::tags$li("Gene-pathway map: shows which causal-graph genes sit inside each overlapping pathway and whether those genes are DE-supported, direct-druggable, or mainly contextual."),
        shiny::tags$li("Direct drug-target evidence means a member gene has external drug-target support from DGIdb or HCDT."),
        shiny::tags$li("Pathway-context evidence means Reactome links the gene to a drug-mediated pathway, but that does not prove the whole pathway is cleanly activated or inhibited in your samples."),
        shiny::tags$li("Use this tab to understand where the DE-related causal subgraph touches known pathways and whether there are plausible pharmacologic leverage points around those components.")
      )
    )
  })
  output$semgraph_pathway_summary_plot <- plotly::renderPlotly({
    ctx <- semgraph_pathway_context()
    shiny::req(nrow(ctx$summary) > 0L)
    top <- head(ctx$summary, 15L)
    top$pathway <- as.character(top$pathway)
    top$genes_in_pathway_view <- suppressWarnings(as.numeric(top$genes_in_pathway_view))
    top$de_supported_genes <- suppressWarnings(as.numeric(top$de_supported_genes))
    top$direct_druggable_genes <- suppressWarnings(as.numeric(top$direct_druggable_genes))
    top$pathway <- factor(top$pathway, levels = rev(top$pathway))
    marker_cols <- ifelse(top$direct_druggable_genes > 0, "#0b7285", "#d9ecf2")
    plotly::layout(
      plotly::plot_ly(
        top,
        x = ~de_supported_genes,
        y = ~pathway,
        size = ~pmax(genes_in_pathway_view, 1),
        type = "scatter",
        mode = "markers",
        marker = list(
          color = marker_cols,
          line = list(color = "#205072", width = 1)
        ),
        text = ~paste(
          "Pathway:", pathway,
          "<br>Genes in displayed graph:", genes_in_pathway_view,
          "<br>DE-supported genes:", de_supported_genes,
          "<br>Directly druggable member genes:", direct_druggable_genes,
          "<br>Dominant DE direction:", dominant_de_direction,
          "<br>Likely drug-effect hint:", pathway_drug_effect_hint,
          "<br>Supporting genes:", supporting_genes,
          "<br>Example direct drug evidence:", direct_drug_examples
        ),
        hoverinfo = "text"
      ),
      title = "Reactome pathway overlap around the displayed causal neighborhood",
      xaxis = list(title = "DE-supported genes in pathway overlap"),
      yaxis = list(title = "")
    )
  })
  output$semgraph_pathway_gene_plot <- plotly::renderPlotly({
    ctx <- semgraph_pathway_context()
    shiny::req(nrow(ctx$edges) > 0L)
    edges <- ctx$edges
    top_paths <- head(ctx$summary$pathway, 12L)
    edges <- edges[edges$pathway %in% top_paths, , drop = FALSE]
    shiny::req(nrow(edges) > 0L)
    edges$pathway <- factor(edges$pathway, levels = rev(top_paths))
    edges$gene_label <- ifelse(is.na(edges$gene_label) | !nzchar(edges$gene_label), edges$gene, edges$gene_label)
    edges$full_gene_name <- ifelse(is.na(edges$full_gene_name) | !nzchar(edges$full_gene_name), "Full gene name unavailable", edges$full_gene_name)
    edges$role <- factor(edges$role, levels = c("focal", "upstream", "downstream", "other"))
    edges$role_symbol <- c(focal = "star", upstream = "triangle-up", downstream = "triangle-down", other = "circle")[as.character(edges$role)]
    edges$marker_size <- ifelse(edges$de_supported, 16, 11)
    edges$fc <- ifelse(is.na(edges$log2FoldChange), 0, edges$log2FoldChange)
    fc_lim <- max(abs(edges$fc), na.rm = TRUE)
    if (!is.finite(fc_lim) || fc_lim == 0) fc_lim <- 1
    y_levels <- rev(unique(edges$gene_label[order(edges$fc, na.last = TRUE)]))
    edges$gene_axis <- factor(edges$gene_label, levels = y_levels)
    p <- plotly::plot_ly()
    if (any(edges$is_druggable)) {
      p <- plotly::add_markers(
        p,
        data = edges[edges$is_druggable, , drop = FALSE],
        x = ~pathway,
        y = ~gene_axis,
        inherit = FALSE,
        hoverinfo = "skip",
        marker = list(
          size = edges$marker_size[edges$is_druggable] + 7,
          color = "rgba(0,0,0,0)",
          line = list(color = "#d4a017", width = 3),
          symbol = "circle-open"
        ),
        showlegend = FALSE
      )
    }
    p <- plotly::add_markers(
      p,
      data = edges,
      x = ~pathway,
      y = ~gene_axis,
      type = "scatter",
      mode = "markers",
      text = ~paste(
        "Gene:", gene_label,
        "<br>Full gene name:", full_gene_name,
        "<br>Pathway:", pathway,
        "<br>Causal role:", role,
        "<br>DE-supported:", ifelse(de_supported, "yes", "no"),
        "<br>log2FC:", signif(fc, 3),
        "<br>Adjusted p-value:", ifelse(is.na(padj), "NA", signif(padj, 3)),
        "<br>Direct druggable member:", ifelse(is_druggable, "yes", "no"),
        "<br>Drug mode hint:", drug_mode,
        "<br>Drug source(s):", drug_source,
        "<br>Drug/context detail:", drug_summary
      ),
      hoverinfo = "text",
      marker = list(
        size = edges$marker_size,
        color = edges$fc,
        colorscale = "RdBu",
        cmin = -fc_lim,
        cmax = fc_lim,
        colorbar = list(title = "log2FC"),
        symbol = edges$role_symbol,
        line = list(color = ifelse(edges$de_supported, "#111111", "#999999"), width = ifelse(edges$de_supported, 2, 1))
      ),
      showlegend = FALSE
    )
    plotly::layout(
      p,
      title = "Where displayed causal genes sit across overlapping Reactome pathways",
      xaxis = list(title = "", tickangle = -35),
      yaxis = list(title = "")
    )
  })
  output$semgraph_pathway_table <- DT::renderDataTable({
    ctx <- semgraph_pathway_context()
    shiny::req(nrow(ctx$summary) > 0L)
    keep <- ctx$summary[, c("pathway", "genes_in_pathway_view", "de_supported_genes", "contains_focal_gene", "upstream_genes", "downstream_genes", "direct_druggable_genes", "dominant_de_direction", "pathway_drug_effect_hint", "supporting_genes", "direct_drug_examples", "drug_sources", "interpretation"), drop = FALSE]
    make_downloadable_table(keep, page_length = 10)
  }, server = FALSE)
  output$semgraph_validation_help <- shiny::renderUI({
    res <- semgraph_res()
    shiny::req(!is.null(res))
    source_status <- validation_source_status(igraph::V(res$graph)$name, allow_live_string = isTRUE(input$semgraph_use_string_live))
    shiny::tagList(
      shiny::tags$p("This validation layer helps judge whether a SEMgraph edge also looks plausible in external biology resources. It should be used to strengthen or weaken confidence in the SEMgraph hypothesis, not to replace the SEMgraph model itself."),
      shiny::tags$ul(
        shiny::tags$li("STRING support means the same source-target gene pair is also seen as a known functional association. It supports biological relatedness, while SEMgraph still provides the arrow direction."),
        shiny::tags$li("STRING action data, when available, can add biological meaning such as activation, inhibition, binding, reaction, or expression-related action. This is a stronger validation layer than pair support alone, but it is available for only a subset of STRING interactions."),
        shiny::tags$li("`STRING direction matches SEMgraph` means STRING action data indicated a direction and that direction agreed with the SEMgraph source-to-target arrow. This is stronger than ordinary STRING pair support."),
        shiny::tags$li("Use the graph edge filter to keep only STRING-supported edges when you want a stricter causal view whose directed SEMgraph arrows are also backed by known STRING connectivity.")
      ),
      shiny::tags$p(sprintf("Current STRING source: %s", if (source_status$string_local) basename(source_status$string_local_path) else if (source_status$string_cache) "cached STRING network" else if (source_status$string_live_enabled) "live STRING lookup enabled" else "none available"))
    )
  })
  output$semgraph_validation_plot <- plotly::renderPlotly({
    res <- semgraph_res()
    shiny::req(!is.null(res), !is.null(res$edge_validation), nrow(res$edge_validation) > 0L)
    df <- res$edge_validation
    df <- filter_validated_edges(df, mode = current_semgraph_filter_rule(), string_cutoff = input$semgraph_string_cutoff %||% 0.7)
    shiny::req(nrow(df) > 0L)
    plot_df <- as.data.frame(table(validation_tier = df$validation_tier), stringsAsFactors = FALSE)
    names(plot_df)[2] <- "n"
    plotly::layout(
      plotly::plot_ly(plot_df, x = ~validation_tier, y = ~n, type = "bar", text = ~n, textposition = "auto"),
      title = "External validation of SEMgraph edges",
      xaxis = list(title = ""),
      yaxis = list(title = "Edge count")
    )
  })
  output$semgraph_edge_table <- DT::renderDataTable({
    res <- semgraph_res()
    shiny::req(!is.null(res), !is.null(res$edge_validation))
    df <- res$edge_validation
    df <- filter_validated_edges(df, mode = current_semgraph_filter_rule(), string_cutoff = input$semgraph_string_cutoff %||% 0.7)
    if (!nrow(df)) {
      return(make_downloadable_table(data.frame(message = "No SEMgraph edges passed the current external-validation filter.", stringsAsFactors = FALSE), page_length = 5))
    }
    keep <- df[, c("source", "target", "validation_tier", "validation_score", "string_supported", "string_score", "string_direction_available", "string_direction_matches", "string_action_mode", "string_action_effect", "string_source", "string_action_source"), drop = FALSE]
    names(keep) <- c("source_gene", "target_gene", "validation_tier", "validation_score", "string_supported", "string_score", "string_direction_available", "string_direction_matches_semgraph", "string_action_mode", "string_action_effect", "string_pair_source", "string_action_source")
    make_downloadable_table(keep, page_length = 12)
  }, server = FALSE)

  current_report_snapshot <- shiny::reactive({
    cur <- tryCatch(ds(), error = function(...) NULL)
    de_res <- tryCatch(de(), error = function(...) NULL)
    enrich_res <- tryCatch(enrich(), error = function(...) NULL)
    wres <- tryCatch(wgcna_res(), error = function(...) NULL)
    sres <- tryCatch(semgraph_res(), error = function(...) NULL)

    out <- list(
      generated_at = as.character(Sys.time()),
      package_versions = package_versions_df(),
      analysis_log = analysis_log(),
      designs = saved(),
      inputs = list(
        counts_path = input$counts_path,
        metadata_path = input$meta_path,
        min_count = input$min_count,
        min_samples = input$min_samples,
        qc_group = input$qc_group,
        de_alpha = input$de_alpha %||% 0.05,
        de_lfc_cutoff = input$de_lfc_cutoff %||% 0
      )
    )

    if (!is.null(cur)) {
      vars <- (cur$pca$sdev ^ 2) / sum(cur$pca$sdev ^ 2)
      pca_scores <- as.data.frame(cur$pca$x[, 1:2, drop = FALSE])
      pca_scores$run_id <- rownames(pca_scores)
      pca_scores$group <- cur$meta[[input$qc_group]]
      out$qc <- list(
        sample_table = data.frame(run_id = cur$meta$run_id, library_size = colSums(cur$counts), detected_genes = colSums(cur$counts > 0), cur$meta[setdiff(names(cur$meta), "run_id")], check.names = FALSE),
        lib_df = data.frame(run_id = cur$meta$run_id, library_size = colSums(cur$counts), group = cur$meta[[input$qc_group]], stringsAsFactors = FALSE),
        pca_df = pca_scores,
        pca_var = vars[1:2],
        dist = cur$dist
      )
    }

    if (!is.null(de_res)) {
      current_deg <- deg_mask(de_res$df, alpha = input$de_alpha %||% 0.05, lfc_cutoff = input$de_lfc_cutoff %||% 0)
      out$de <- list(
        label = de_res$label,
        sig = de_res$sig,
        pvalue_fraction = de_res$pvalue_fraction,
        padj_fraction = de_res$padj_fraction,
        current_deg_count = sum(current_deg, na.rm = TRUE),
        current_deg_fraction = mean(current_deg, na.rm = TRUE),
        table = de_res$df,
        mean_sd = de_res$mean_sd,
        size_factors = data.frame(sample = names(de_res$size_factors), size_factor = as.numeric(de_res$size_factors), stringsAsFactors = FALSE)
      )
    }

    if (!is.null(enrich_res)) {
      out$enrichment <- list(
        settings = list(
          methods = input$enrich_methods %||% character(),
          padj = input$enrich_padj,
          abs_log2fc = input$enrich_lfc,
          min_gs = input$enrich_min_gs,
          max_gs = input$enrich_max_gs,
          adjust_method = input$enrich_adjust_method,
          top_n = input$enrich_top_n
        ),
        consensus = enrich_res$consensus,
        combined = enrich_res$combined
      )
    }

    if (!is.null(wres)) {
      out$wgcna <- list(
        module_sizes = data.frame(module = names(wres$module_sizes), genes = as.numeric(wres$module_sizes), stringsAsFactors = FALSE),
        selected_module = input$wgcna_module_view,
        gene_modules = wres$gene_modules
      )
      if (nzchar(input$wgcna_module_view)) {
        mod_graph <- build_wgcna_module_graph(wres, input$wgcna_module_view, edge_quantile = input$wgcna_edge_quantile)
        if (!is.null(mod_graph)) {
          adj <- igraph::as_adjacency_matrix(mod_graph, attr = "weight", sparse = FALSE)
          dist_mat <- 1 - adj
          diag(dist_mat) <- 0
          coords <- as.data.frame(cmdscale(stats::as.dist(dist_mat), k = 2))
          coords$gene_id <- rownames(coords)
          out$wgcna$module_graph <- list(
            nodes = coords,
            edges = igraph::as_data_frame(mod_graph, what = "edges")
          )
        }
      }
    }

    if (!is.null(sres)) {
      sel <- semgraph_graph_selection()
      sub_nodes <- sel$nodes[sel$nodes %in% igraph::V(sres$graph)$name]
      sg <- igraph::induced_subgraph(sres$graph, vids = sub_nodes)
      sg_edges <- annotate_semgraph_edges(sg, sres)
      sg_edges <- filter_validated_edges(sg_edges, mode = current_semgraph_filter_rule(), string_cutoff = input$semgraph_string_cutoff %||% 0.7)
      if (nrow(sg_edges)) {
        keep_nodes <- unique(c(sg_edges$source, sg_edges$target))
        sg <- igraph::induced_subgraph(sg, vids = keep_nodes)
      } else if (!identical(current_semgraph_filter_rule(), "all")) {
        sg <- igraph::delete_vertices(sg, igraph::V(sg))
      }
      if (igraph::vcount(sg) > 0L) {
        lay <- compute_semgraph_layout(sg, method = input$semgraph_layout)
        node_df <- data.frame(x = lay[, 1], y = lay[, 2], gene = igraph::V(sg)$name, stringsAsFactors = FALSE)
        node_df$role <- ifelse(node_df$gene == sres$focal_gene, "focal anchor", ifelse(node_df$gene %in% sres$upstream, "upstream candidate", ifelse(node_df$gene %in% sres$downstream, "downstream candidate", "context gene")))
      } else {
        node_df <- data.frame(x = numeric(), y = numeric(), gene = character(), role = character(), stringsAsFactors = FALSE)
      }
      sg_edges <- annotate_semgraph_edges(sg, sres)
      sg_edges <- filter_validated_edges(sg_edges, mode = current_semgraph_filter_rule(), string_cutoff = input$semgraph_string_cutoff %||% 0.7)
      out$semgraph <- list(
        analysis_unit = sres$analysis_unit,
        focal_gene = sres$focal_gene,
        degree_table = tryCatch({
          df <- sres$degree_df
          df[df$gene %in% select_semgraph_table_nodes(sres, tryCatch(de()$df, error = function(...) NULL), scope = input$semgraph_table_scope, max_rows = input$semgraph_table_max_rows, alpha = input$de_alpha %||% 0.05, lfc_cutoff = input$de_lfc_cutoff %||% 0), , drop = FALSE]
        }, error = function(...) data.frame()),
        edge_validation_table = tryCatch({
          filter_validated_edges(sres$edge_validation, mode = current_semgraph_filter_rule(), string_cutoff = input$semgraph_string_cutoff %||% 0.7)
        }, error = function(...) data.frame()),
        ace_table = sres$ace,
        pathway_summary = tryCatch(semgraph_pathway_context()$summary, error = function(...) data.frame()),
        graph = list(
          nodes = node_df,
          edges = sg_edges
        )
      )
    }

    out
  })

  build_report_bundle <- function(snapshot, out_dir, html_mode = FALSE) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    tables_dir <- file.path(out_dir, "tables")
    plots_dir <- file.path(out_dir, "plots")
    dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(snapshot, file.path(out_dir, "report_snapshot.rds"))
    utils::write.csv(snapshot$package_versions, file.path(out_dir, "package_versions.csv"), row.names = FALSE)
    utils::write.csv(snapshot$analysis_log, file.path(out_dir, "analysis_log.csv"), row.names = FALSE)

    load_params <- report_parse_params(snapshot$analysis_log, "load_data")
    de_params <- report_parse_params(snapshot$analysis_log, "run_deseq2")
    enrich_params <- report_parse_params(snapshot$analysis_log, "run_enrichment")
    wgcna_params <- report_parse_params(snapshot$analysis_log, "run_wgcna")
    semgraph_params <- report_parse_params(snapshot$analysis_log, "run_semgraph")

    active_design_name <- de_params$active_design %||% ""
    design_tbl <- snapshot$designs %||% data.frame()
    design_formula <- if (nrow(design_tbl) && nzchar(active_design_name) && active_design_name %in% design_tbl$name) {
      design_tbl$formula[match(active_design_name, design_tbl$name)]
    } else {
      "~ 1"
    }
    load_counts_path <- snapshot$inputs$counts_path %||% load_params$counts_path %||% ""
    load_meta_path <- snapshot$inputs$metadata_path %||% load_params$metadata_path %||% ""
    qc_code_lines <- c(
      "source('report_methods.R')",
      sprintf("counts_path <- %s", r_literal(load_counts_path)),
      sprintf("metadata_path <- %s", r_literal(load_meta_path)),
      sprintf("min_count <- %s", r_literal(snapshot$inputs$min_count %||% load_params$min_count %||% 10)),
      sprintf("min_samples <- %s", r_literal(snapshot$inputs$min_samples %||% load_params$min_samples %||% 2)),
      "cur <- read_inputs(counts_path, metadata_path, min_count, min_samples)",
      "qc_sample_table <- data.frame(",
      "  run_id = cur$meta$run_id,",
      "  library_size = colSums(cur$counts),",
      "  detected_genes = colSums(cur$counts > 0),",
      "  cur$meta[setdiff(names(cur$meta), 'run_id')],",
      "  check.names = FALSE",
      ")",
      sprintf("qc_group <- %s", r_literal(snapshot$inputs$qc_group %||% "treatment")),
      "qc_lib_df <- data.frame(run_id = cur$meta$run_id, library_size = colSums(cur$counts), group = cur$meta[[qc_group]], stringsAsFactors = FALSE)",
      "qc_vars <- (cur$pca$sdev ^ 2) / sum(cur$pca$sdev ^ 2)",
      "qc_scores <- as.data.frame(cur$pca$x[, 1:2, drop = FALSE])",
      "qc_scores$run_id <- rownames(qc_scores)",
      "qc_scores$group <- cur$meta[[qc_group]]"
    )
    de_code_lines <- c(
      "source('report_methods.R')",
      sprintf("counts_path <- %s", r_literal(load_counts_path)),
      sprintf("metadata_path <- %s", r_literal(load_meta_path)),
      sprintf("min_count <- %s", r_literal(snapshot$inputs$min_count %||% load_params$min_count %||% 10)),
      sprintf("min_samples <- %s", r_literal(snapshot$inputs$min_samples %||% load_params$min_samples %||% 2)),
      "cur <- read_inputs(counts_path, metadata_path, min_count, min_samples)",
      sprintf("selected_samples <- %s", r_literal(de_params$selected_samples %||% character())),
      sprintf("formula_str <- %s", r_literal(design_formula)),
      sprintf("contrast_var <- %s", r_literal(de_params$comparison_variable %||% NULL)),
      sprintf("contrast_num <- %s", r_literal(de_params$numerator %||% NULL)),
      sprintf("contrast_den <- %s", r_literal(de_params$denominator %||% NULL)),
      sprintf("alpha <- %s", r_literal(snapshot$inputs$de_alpha %||% de_params$alpha %||% 0.05)),
      sprintf("lfc_cutoff <- %s", r_literal(snapshot$inputs$de_lfc_cutoff %||% de_params$abs_log2fc_cutoff %||% 0)),
      "de_res <- run_de(cur$filtered, cur$meta, samples = selected_samples, formula_str = formula_str, gene_map = cur$gene_map, var = contrast_var, num = contrast_num, den = contrast_den)",
      "de_table <- de_res$df",
      "de_table$meets_current_deg_cutoff <- deg_mask(de_table, alpha = alpha, lfc_cutoff = lfc_cutoff)",
      "volcano_df <- de_table",
      "volcano_df$state <- ifelse(deg_mask(volcano_df, alpha = alpha, lfc_cutoff = lfc_cutoff) & volcano_df$log2FoldChange > 0, sprintf('Upregulated DEG: padj < %.3f and |log2FC| >= %.2f', alpha, lfc_cutoff), ifelse(deg_mask(volcano_df, alpha = alpha, lfc_cutoff = lfc_cutoff) & volcano_df$log2FoldChange < 0, sprintf('Downregulated DEG: padj < %.3f and |log2FC| >= %.2f', alpha, lfc_cutoff), 'Outside current DEG cutoff'))",
      "volcano_df$nlp <- -log10(pmax(volcano_df$padj, .Machine$double.xmin))",
      "volcano_colors <- c('#c0392b', '#1f78b4', '#9aa5b1')",
      "names(volcano_colors) <- c(sprintf('Upregulated DEG: padj < %.3f and |log2FC| >= %.2f', alpha, lfc_cutoff), sprintf('Downregulated DEG: padj < %.3f and |log2FC| >= %.2f', alpha, lfc_cutoff), 'Outside current DEG cutoff')",
      "plot_ly(volcano_df, x = ~log2FoldChange, y = ~nlp, color = ~state, colors = volcano_colors, type = 'scatter', mode = 'markers', text = ~paste('Gene:', gene_label, '<br>Full gene name:', full_gene_name, '<br>Ensembl:', gene_id, '<br>log2FC:', signif(log2FoldChange, 3), '<br>padj:', signif(padj, 3)), hoverinfo = 'text')"
    )
    enrich_code_lines <- c(
      "source('report_methods.R')",
      "# Run the DESeq2 block above first so `de_res` is available.",
      sprintf("enrich_methods <- %s", r_literal(enrich_params$methods %||% character())),
      sprintf("enrich_padj_cutoff <- %s", r_literal(enrich_params$padj_cutoff %||% snapshot$enrichment$settings$padj %||% 0.05)),
      sprintf("enrich_abs_lfc_cutoff <- %s", r_literal(enrich_params$abs_log2fc_cutoff %||% snapshot$enrichment$settings$abs_log2fc %||% 0)),
      sprintf("min_pathway_size <- %s", r_literal(enrich_params$min_pathway_size %||% snapshot$enrichment$settings$min_gs %||% 10)),
      sprintf("max_pathway_size <- %s", r_literal(enrich_params$max_pathway_size %||% snapshot$enrichment$settings$max_gs %||% 500)),
      sprintf("adjust_method <- %s", r_literal(enrich_params$adjust_method %||% snapshot$enrichment$settings$adjust_method %||% "BH")),
      "sig_tbl <- subset(de_res$df, deg_mask(de_res$df, alpha = enrich_padj_cutoff, lfc_cutoff = enrich_abs_lfc_cutoff) & !is.na(entrez_id))",
      "ora_genes <- unique(sig_tbl$entrez_id)",
      "rank_tbl <- de_res$df[!is.na(de_res$df$stat) & !is.na(de_res$df$entrez_id), c('entrez_id','stat'), drop = FALSE]",
      "rank_tbl <- rank_tbl[order(abs(rank_tbl$stat), decreasing = TRUE), , drop = FALSE]",
      "rank_tbl <- rank_tbl[!duplicated(rank_tbl$entrez_id), , drop = FALSE]",
      "ranked_stats <- rank_tbl$stat; names(ranked_stats) <- rank_tbl$entrez_id; ranked_stats <- sort(ranked_stats, decreasing = TRUE)",
      "if ('GO ORA' %in% enrich_methods) ego_ora <- clusterProfiler::enrichGO(gene = ora_genes, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = 'ENTREZID', ont = 'BP', pAdjustMethod = adjust_method, pvalueCutoff = enrich_padj_cutoff, qvalueCutoff = enrich_padj_cutoff, minGSSize = min_pathway_size, maxGSSize = max_pathway_size, readable = TRUE)",
      "if ('KEGG ORA' %in% enrich_methods) kegg_ora <- clusterProfiler::enrichKEGG(gene = ora_genes, organism = 'hsa', keyType = 'ncbi-geneid', pAdjustMethod = adjust_method, pvalueCutoff = enrich_padj_cutoff, qvalueCutoff = enrich_padj_cutoff, minGSSize = min_pathway_size, maxGSSize = max_pathway_size)",
      "if ('GO GSEA' %in% enrich_methods) ego_gsea <- clusterProfiler::gseGO(geneList = ranked_stats, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = 'ENTREZID', ont = 'BP', pAdjustMethod = adjust_method, minGSSize = min_pathway_size, maxGSSize = max_pathway_size, eps = 0)",
      "if ('KEGG GSEA' %in% enrich_methods) kegg_gsea <- clusterProfiler::gseKEGG(geneList = ranked_stats, organism = 'hsa', keyType = 'ncbi-geneid', pAdjustMethod = adjust_method, minGSSize = min_pathway_size, maxGSSize = max_pathway_size, eps = 0)"
    )
    wgcna_code_lines <- c(
      "source('report_methods.R')",
      "# Run the QC block above first so `cur` is available.",
      sprintf("top_variable_genes <- %s", r_literal(wgcna_params$top_variable_genes %||% 4000)),
      sprintf("power_mode <- %s", r_literal(wgcna_params$power_mode %||% "auto")),
      sprintf("manual_power <- %s", r_literal(wgcna_params$manual_power %||% 6)),
      sprintf("min_module_size <- %s", r_literal(wgcna_params$min_module_size %||% 30)),
      sprintf("deep_split <- %s", r_literal(wgcna_params$deep_split %||% 2)),
      sprintf("merge_cut_height <- %s", r_literal(wgcna_params$merge_cut_height %||% 0.2)),
      "wgcna_expr <- t(cur$log_cpm)",
      "gene_var <- apply(wgcna_expr, 2, stats::var, na.rm = TRUE)",
      "keep_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(top_variable_genes, length(gene_var)))]",
      "wgcna_expr <- wgcna_expr[, keep_genes, drop = FALSE]",
      "soft <- WGCNA::pickSoftThreshold(wgcna_expr, powerVector = c(1:10, 12, 14, 16, 18, 20), verbose = 0)",
      "soft_power <- if (identical(power_mode, 'manual')) manual_power else soft$powerEstimate %||% manual_power",
      "adjacency <- WGCNA::adjacency(wgcna_expr, power = soft_power, type = 'signed', corFnc = 'cor', corOptions = list(use = 'p'))",
      "tom <- WGCNA::TOMsimilarity(adjacency, TOMType = 'signed')",
      "diss_tom <- 1 - tom",
      "gene_tree <- stats::hclust(stats::as.dist(diss_tom), method = 'average')",
      "dynamic_modules <- dynamicTreeCut::cutreeDynamic(dendro = gene_tree, distM = diss_tom, deepSplit = deep_split, pamRespectsDendro = FALSE, minClusterSize = min_module_size)",
      "module_colors <- WGCNA::labels2colors(dynamic_modules)",
      "merged <- WGCNA::mergeCloseModules(wgcna_expr, module_colors, cutHeight = merge_cut_height, verbose = 0)",
      "module_colors <- merged$colors",
      sprintf("selected_module <- %s", r_literal(snapshot$wgcna$selected_module %||% "")),
      "module_gene_ids <- colnames(wgcna_expr)[module_colors == selected_module]",
      "module_map <- cur$gene_map[match(module_gene_ids, cur$gene_map$gene_id), , drop = FALSE]",
      "module_kwithin <- WGCNA::softConnectivity(wgcna_expr[, module_gene_ids, drop = FALSE], power = soft_power)",
      "module_plot_df <- data.frame(gene_id = module_gene_ids, kWithin = as.numeric(module_kwithin), stringsAsFactors = FALSE)",
      "module_plot_df <- attach_gene_labels(module_plot_df, cur$gene_map)"
    )
    semgraph_code_lines <- c(
      "source('report_methods.R')",
      "# Run the QC and WGCNA blocks above first so `cur` and `module_colors` are available.",
      sprintf("analysis_unit <- %s", r_literal(semgraph_params$analysis_unit %||% snapshot$semgraph$analysis_unit %||% "")),
      sprintf("display_mode <- %s", r_literal(semgraph_params$display_mode %||% "de_connectors")),
      sprintf("max_model_genes <- %s", r_literal(semgraph_params$max_model_genes %||% 300)),
      sprintf("max_annotation_genes <- %s", r_literal(semgraph_params$max_annotation_genes %||% 180)),
      sprintf("use_live_string <- %s", r_literal(semgraph_params$use_live_string %||% FALSE)),
      "# The app then prioritizes genes, builds a prior graph from the selected WGCNA unit, orients it to a DAG with SEMgraph, and annotates the resulting edges with STRING support and STRING action data when available.",
      "# See the exported semgraph graph/table files for the exact visible state that was rendered in the app."
    )

    if (!is.null(snapshot$qc)) {
      utils::write.csv(snapshot$qc$sample_table, file.path(tables_dir, "qc_sample_table.csv"), row.names = FALSE)
      utils::write.csv(snapshot$qc$lib_df, file.path(tables_dir, "qc_library_sizes.csv"), row.names = FALSE)
    }
    if (!is.null(snapshot$de)) {
      utils::write.csv(snapshot$de$table, file.path(tables_dir, "deseq2_results.csv"), row.names = FALSE)
      utils::write.csv(snapshot$de$mean_sd, file.path(tables_dir, "deseq2_mean_sd.csv"), row.names = FALSE)
    }
    if (!is.null(snapshot$enrichment)) {
      utils::write.csv(snapshot$enrichment$consensus, file.path(tables_dir, "enrichment_consensus.csv"), row.names = FALSE)
      utils::write.csv(snapshot$enrichment$combined, file.path(tables_dir, "enrichment_combined.csv"), row.names = FALSE)
    }
    if (!is.null(snapshot$wgcna)) {
      utils::write.csv(snapshot$wgcna$module_sizes, file.path(tables_dir, "wgcna_module_sizes.csv"), row.names = FALSE)
      utils::write.csv(snapshot$wgcna$gene_modules, file.path(tables_dir, "wgcna_gene_modules.csv"), row.names = FALSE)
    }
    if (!is.null(snapshot$semgraph)) {
      utils::write.csv(snapshot$semgraph$degree_table, file.path(tables_dir, "semgraph_degree_table.csv"), row.names = FALSE)
      utils::write.csv(snapshot$semgraph$edge_validation_table, file.path(tables_dir, "semgraph_edge_validation_table.csv"), row.names = FALSE)
      utils::write.csv(snapshot$semgraph$ace_table, file.path(tables_dir, "semgraph_effects_table.csv"), row.names = FALSE)
      utils::write.csv(snapshot$semgraph$pathway_summary, file.path(tables_dir, "semgraph_pathway_summary.csv"), row.names = FALSE)
    }

    report_methods_lines <- c(
      "`%||%` <- function(x, y) if (is.null(x)) y else x",
      "",
      "map_gene_symbols <- function(gene_ids) {",
      "  clean_ids <- sub('\\\\..*$', '', gene_ids)",
      "  symbol_map <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = unique(clean_ids), column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')",
      "  entrez_map <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = unique(clean_ids), column = 'ENTREZID', keytype = 'ENSEMBL', multiVals = 'first')",
      "  fullname_map <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = unique(clean_ids), column = 'GENENAME', keytype = 'ENSEMBL', multiVals = 'first')",
      "  symbols <- unname(symbol_map[clean_ids])",
      "  entrez <- unname(entrez_map[clean_ids])",
      "  full_names <- unname(fullname_map[clean_ids])",
      "  symbols[is.na(symbols) | !nzchar(symbols)] <- gene_ids[is.na(symbols) | !nzchar(symbols)]",
      "  full_names[is.na(full_names) | !nzchar(full_names)] <- 'Full gene name unavailable'",
      "  display <- ifelse(duplicated(symbols) | duplicated(symbols, fromLast = TRUE), paste0(symbols, ' (', gene_ids, ')'), symbols)",
      "  data.frame(gene_id = gene_ids, ensembl_id = clean_ids, entrez_id = entrez, hgnc_symbol = symbols, full_gene_name = full_names, gene_label = display, stringsAsFactors = FALSE)",
      "}",
      "",
      "attach_gene_labels <- function(df, gene_map) merge(df, gene_map, by = 'gene_id', all.x = TRUE, sort = FALSE)",
      "",
      "read_inputs <- function(counts_path, meta_path, min_count, min_samples) {",
      "  counts_df <- read.csv(counts_path, check.names = FALSE, stringsAsFactors = FALSE)",
      "  meta <- read.csv(meta_path, check.names = FALSE, stringsAsFactors = FALSE)",
      "  sample_key <- if ('run_id' %in% names(meta)) 'run_id' else if ('sample_id' %in% names(meta)) 'sample_id' else stop('Metadata must contain run_id or sample_id')",
      "  if (sample_key != 'run_id') names(meta)[names(meta) == sample_key] <- 'run_id'",
      "  counts <- as.matrix(counts_df[, -1L, drop = FALSE])",
      "  storage.mode(counts) <- 'numeric'",
      "  rownames(counts) <- counts_df[[1L]]",
      "  samples <- colnames(counts)",
      "  meta <- merge(data.frame(run_id = samples, ord = seq_along(samples), stringsAsFactors = FALSE), meta, by = 'run_id', all.x = TRUE, sort = FALSE)",
      "  meta <- meta[order(meta$ord), , drop = FALSE]",
      "  meta$ord <- NULL",
      "  for (nm in setdiff(names(meta), 'run_id')) {",
      "    meta[[nm]][is.na(meta[[nm]]) | meta[[nm]] == ''] <- 'Unannotated'",
      "    if (is.character(meta[[nm]])) meta[[nm]] <- factor(meta[[nm]])",
      "  }",
      "  keep <- rowSums(counts >= min_count) >= min_samples",
      "  filtered <- counts[keep, , drop = FALSE]",
      "  lib <- colSums(filtered)",
      "  log_cpm <- log2(sweep(filtered, 2L, lib / 1e6, '/') + 1)",
      "  pca <- prcomp(t(log_cpm), center = TRUE, scale. = TRUE)",
      "  gene_map <- map_gene_symbols(rownames(counts))",
      "  list(counts = counts, filtered = filtered, meta = meta, log_cpm = log_cpm, pca = pca, gene_map = gene_map)",
      "}",
      "",
      "deg_mask <- function(df, alpha = 0.05, lfc_cutoff = 0) {",
      "  !is.na(df$padj) & df$padj < alpha & !is.na(df$log2FoldChange) & abs(df$log2FoldChange) >= lfc_cutoff",
      "}",
      "",
      "run_de <- function(counts, meta, samples, formula_str, gene_map, var = NULL, num = NULL, den = NULL) {",
      "  meta_sub <- meta[match(samples, meta$run_id), , drop = FALSE]",
      "  rownames(meta_sub) <- meta_sub$run_id",
      "  dds <- DESeq2::DESeqDataSetFromMatrix(round(counts[, samples, drop = FALSE]), meta_sub, stats::as.formula(formula_str))",
      "  dds <- DESeq2::DESeq(dds, quiet = TRUE)",
      "  if (!is.null(var) && !is.na(var) && nzchar(var)) {",
      "    res <- DESeq2::results(dds, contrast = c(var, num, den))",
      "  } else {",
      "    rn <- DESeq2::resultsNames(dds); res <- DESeq2::results(dds, name = rn[2L])",
      "  }",
      "  df <- as.data.frame(res); df$gene_id <- rownames(df); df <- attach_gene_labels(df, gene_map)",
      "  vst_obj <- DESeq2::vst(dds, blind = TRUE)",
      "  list(df = df, vst = SummarizedExperiment::assay(vst_obj), size_factors = DESeq2::sizeFactors(dds))",
      "}"
    )
    write_text_file(file.path(out_dir, "report_methods.R"), report_methods_lines)
    report_methods_block <- commented_code_block(report_methods_lines)
    write_text_file(file.path(out_dir, "analysis_commands.R"), c(
      "# RNA-seq Workbench commands generated from the Shiny session",
      "",
      "# QC / input loading",
      commented_code_block(qc_code_lines),
      "",
      "# DESeq2",
      commented_code_block(de_code_lines),
      "",
      "# Pathway enrichment",
      commented_code_block(enrich_code_lines),
      "",
      "# WGCNA",
      commented_code_block(wgcna_code_lines),
      "",
      "# SEMgraph",
      commented_code_block(semgraph_code_lines)
    ))

    html_method_card <- function(title, logic_lines, arg_lines = NULL, usage_lines = NULL) {
      c(
        "<div class='method-box'>",
        sprintf("<h3>%s</h3>", htmltools::htmlEscape(title)),
        vapply(logic_lines, function(x) sprintf("<p>%s</p>", htmltools::htmlEscape(x)), character(1)),
        if (!is.null(arg_lines) && length(arg_lines)) c(
          "<h4>Why the chosen arguments matter</h4>",
          "<ul>",
          vapply(arg_lines, function(x) sprintf("<li>%s</li>", htmltools::htmlEscape(x)), character(1)),
          "</ul>"
        ) else NULL,
        if (!is.null(usage_lines) && length(usage_lines)) c(
          "<h4>How to use the elements in this section</h4>",
          "<ul>",
          vapply(usage_lines, function(x) sprintf("<li>%s</li>", htmltools::htmlEscape(x)), character(1)),
          "</ul>"
        ) else NULL,
        "</div>"
      )
    }

    if (html_mode) {
      live_cur <- tryCatch(ds(), error = function(...) NULL)
      live_de <- tryCatch(de(), error = function(...) NULL)
      live_wgcna <- tryCatch(wgcna_res(), error = function(...) NULL)
      if (!is.null(snapshot$qc)) {
        qc_df <- if (!is.null(live_cur)) {
          data.frame(run_id = live_cur$meta$run_id, library_size = colSums(live_cur$counts), group = live_cur$meta[[input$qc_group]], stringsAsFactors = FALSE)
        } else {
          snapshot$qc$lib_df
        }
        p <- make_downloadable_plot(plotly::layout(plotly::plot_ly(qc_df, x = ~run_id, y = ~library_size, color = ~group, type = "bar"), title = "QC library sizes"), "qc_library_sizes")
        save_plot_widget(p, file.path(plots_dir, "qc_library_sizes.html"), "QC library sizes")
      }
      if (!is.null(snapshot$de)) {
        ddf <- if (!is.null(live_de)) live_de$df else snapshot$de$table
        alpha_now <- snapshot$inputs$de_alpha %||% 0.05
        lfc_now <- snapshot$inputs$de_lfc_cutoff %||% 0
        ddf$state <- ifelse(
          deg_mask(ddf, alpha = alpha_now, lfc_cutoff = lfc_now) & ddf$log2FoldChange > 0,
          sprintf("Upregulated DEG: padj < %.3f and |log2FC| >= %.2f", alpha_now, lfc_now),
          ifelse(
            deg_mask(ddf, alpha = alpha_now, lfc_cutoff = lfc_now) & ddf$log2FoldChange < 0,
            sprintf("Downregulated DEG: padj < %.3f and |log2FC| >= %.2f", alpha_now, lfc_now),
            "Outside current DEG cutoff"
          )
        )
        ddf$nlp <- -log10(pmax(ddf$padj, .Machine$double.xmin))
        volcano_colors <- c(
          "Upregulated DEG: padj < 0.050 and |log2FC| >= 0.00" = "#c0392b",
          "Downregulated DEG: padj < 0.050 and |log2FC| >= 0.00" = "#1f78b4",
          "Outside current DEG cutoff" = "#9aa5b1"
        )
        names(volcano_colors)[1] <- sprintf("Upregulated DEG: padj < %.3f and |log2FC| >= %.2f", alpha_now, lfc_now)
        names(volcano_colors)[2] <- sprintf("Downregulated DEG: padj < %.3f and |log2FC| >= %.2f", alpha_now, lfc_now)
        p <- make_downloadable_plot(plotly::layout(
          plotly::plot_ly(
            ddf,
            x = ~log2FoldChange,
            y = ~nlp,
            color = ~state,
            colors = volcano_colors,
            type = "scatter",
            mode = "markers",
            text = ~paste(
              "Gene:", gene_label,
              "<br>Full gene name:", full_gene_name,
              "<br>Ensembl:", gene_id,
              "<br>log2FC:", signif(log2FoldChange, 3),
              "<br>padj:", signif(padj, 3)
            ),
            hoverinfo = "text"
          ),
          shapes = list(
            list(type = "line", x0 = -lfc_now, x1 = -lfc_now, y0 = 0, y1 = max(ddf$nlp, na.rm = TRUE), line = list(color = "#6b7280", dash = "dot")),
            list(type = "line", x0 = lfc_now, x1 = lfc_now, y0 = 0, y1 = max(ddf$nlp, na.rm = TRUE), line = list(color = "#6b7280", dash = "dot")),
            list(type = "line", x0 = min(ddf$log2FoldChange, na.rm = TRUE), x1 = max(ddf$log2FoldChange, na.rm = TRUE), y0 = -log10(alpha_now), y1 = -log10(alpha_now), line = list(color = "#6b7280", dash = "dot"))
          ),
          title = "Volcano plot",
          xaxis = list(title = "log2 fold change"),
          yaxis = list(title = "-log10 adjusted p-value")
        ), "deseq2_volcano")
        save_plot_widget(p, file.path(plots_dir, "deseq2_volcano.html"), "DESeq2 volcano-style view")
      }
      if (!is.null(snapshot$enrichment) && nrow(snapshot$enrichment$consensus) > 0L) {
        top <- head(snapshot$enrichment$consensus, 12L)
        p <- make_downloadable_plot(plotly::layout(plotly::plot_ly(top, x = ~consensus_hits, y = ~pathway, type = "scatter", mode = "markers"), title = "Enrichment consensus"), "enrichment_consensus")
        save_plot_widget(p, file.path(plots_dir, "enrichment_consensus.html"), "Enrichment consensus")
      }
      if (!is.null(snapshot$wgcna) && nzchar(snapshot$wgcna$selected_module %||% "") && !is.null(live_wgcna) && !is.null(live_cur)) {
        g_mod <- build_wgcna_module_graph(live_wgcna, snapshot$wgcna$selected_module, edge_quantile = input$wgcna_edge_quantile %||% 0.9)
        if (!is.null(g_mod)) {
          gm <- live_wgcna$gene_modules
          adj <- igraph::as_adjacency_matrix(g_mod, attr = "weight", sparse = FALSE)
          diag(adj) <- 0
          edge_tbl <- igraph::as_data_frame(g_mod, what = "edges")
          dist_mat <- 1 - adj
          diag(dist_mat) <- 0
          coords <- cmdscale(stats::as.dist(dist_mat), k = 2)
          coords <- as.data.frame(coords)
          coords$gene_id <- rownames(coords)
          coords <- attach_gene_labels(coords, live_cur$gene_map)
          coords$kWithin <- gm$kWithin[match(coords$gene_id, gm$gene_id)]
          edge_df <- if (nrow(edge_tbl) > 0L) {
            data.frame(
              x = coords$V1[match(edge_tbl$from, coords$gene_id)],
              y = coords$V2[match(edge_tbl$from, coords$gene_id)],
              xend = coords$V1[match(edge_tbl$to, coords$gene_id)],
              yend = coords$V2[match(edge_tbl$to, coords$gene_id)],
              stringsAsFactors = FALSE
            )
          } else data.frame()
          p <- plotly::plot_ly()
          if (nrow(edge_df) > 0L) {
            for (i in seq_len(nrow(edge_df))) {
              p <- plotly::add_segments(
                p,
                x = edge_df$x[i], y = edge_df$y[i],
                xend = edge_df$xend[i], yend = edge_df$yend[i],
                inherit = FALSE,
                line = list(color = "rgba(120,120,120,0.35)", width = 1)
              )
            }
          }
          p <- plotly::add_markers(
            p,
            data = coords,
            x = ~V1,
            y = ~V2,
            color = I(snapshot$wgcna$selected_module),
            text = ~paste("Gene:", gene_label, "<br>Full gene name:", full_gene_name, "<br>Ensembl:", gene_id, "<br>Within-module connectivity:", signif(kWithin, 3)),
            hoverinfo = "text",
            marker = list(
              size = {
                kw <- coords$kWithin
                kw[!is.finite(kw)] <- 0
                rng <- range(kw, na.rm = TRUE)
                if (!all(is.finite(rng)) || diff(rng) == 0) rep(14, length(kw)) else 10 + 18 * (kw - rng[1]) / diff(rng)
              },
              color = snapshot$wgcna$selected_module,
              line = list(color = "#1f1f1f", width = 0.5)
            )
          )
          p <- plotly::add_text(p, data = coords, x = ~V1, y = ~V2, text = ~gene_label, textposition = "top center", inherit = FALSE, showlegend = FALSE)
          save_plot_widget(make_downloadable_plot(plotly::layout(p, title = "Selected WGCNA module network", xaxis = list(title = "", showticklabels = FALSE, zeroline = FALSE), yaxis = list(title = "", showticklabels = FALSE, zeroline = FALSE)), "wgcna_module_network"), file.path(plots_dir, "wgcna_module_network.html"), "Selected WGCNA module network")
        }
      }
      if (!is.null(snapshot$semgraph) && !is.null(snapshot$semgraph$graph)) {
        save_plot_widget(build_semgraph_graph_plot(), file.path(plots_dir, "semgraph_subgraph.html"), "Displayed SEMgraph subgraph")
      }
    }

    html_code_panel <- function(id, title, code_lines, note = NULL) {
      escaped <- htmltools::htmlEscape(commented_code_block(code_lines))
      c(
        sprintf("<details class='code-dropdown'><summary>%s</summary>", htmltools::htmlEscape(title)),
        if (!is.null(note) && nzchar(note)) sprintf("<p class='code-note'>%s</p>", htmltools::htmlEscape(note)) else NULL,
        sprintf("<div class='code-toolbar'><button type='button' class='copy-btn' onclick=\"copyCodeBlock('%s', this)\">Copy code</button></div>", id),
        sprintf("<pre class='code-block'><code id='%s'>%s</code></pre>", id, escaped),
        "</details>"
      )
    }

    rmd_lines <- c(
      "---",
      "title: \"RNA Hypothesis Workbench Reproducibility Report\"",
      "output:",
      "  html_document:",
      "    toc: true",
      "    toc_float: true",
      "    code_folding: hide",
      "---",
      "",
      "```{r setup, include=FALSE}",
      "knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)",
      "library(jsonlite)",
      "library(plotly)",
      "snap <- readRDS('report_snapshot.rds')",
      "source('report_methods.R')",
      "logs <- snap$analysis_log",
      "parse_params <- function(x) {",
      "  if (is.null(x) || is.na(x) || !nzchar(x)) return(list())",
      "  jsonlite::fromJSON(x, simplifyVector = FALSE)",
      "}",
      "log_rows <- if (nrow(logs)) {",
      "  data.frame(",
      "    timestamp = logs$timestamp,",
      "    step = logs$step,",
      "    stringsAsFactors = FALSE",
      "  )",
      "} else data.frame()",
      "last_params <- function(step_name) {",
      "  idx <- which(logs$step == step_name)",
      "  if (!length(idx)) return(list())",
      "  parse_params(logs$params_json[max(idx)])",
      "}",
      "show_params <- function(x) {",
      "  if (!length(x)) return('No recorded parameters for this step.')",
      "  capture.output(str(x))",
      "}",
      "```",
      "",
      "# How to use this report",
      "",
      "This report is designed to let you inspect the exact app settings that were used, understand why each method exists, and rerun the saved workflow with explicit R expressions and arguments.",
      "",
      "Each section has three parts:",
      "",
      "1. What the method means and why it matters.",
      "2. Which arguments were chosen in the Shiny app.",
      "3. Executable R chunks together with copy-ready commands generated from the Shiny session.",
      "",
      "# Generated commands",
      "",
      "The bundle also includes `analysis_commands.R`, which contains copy-ready commands for another R session, and `report_methods.R`, which provides the helper functions those commands use.",
      "",
      "How to reuse the code in this report:",
      "",
      "1. Define the helper functions first or source `report_methods.R`.",
      "2. Run the generated workflow commands in `analysis_commands.R` from top to bottom in a fresh R session.",
      "3. Compare the recreated objects and plots with the exported tables and HTML widgets in this bundle.",
      "4. If you want to adapt the workflow, change the argument values shown in the generated commands rather than editing the report snapshot directly.",
      "",
      "# Core helper functions used by the generated commands",
      "",
      "These helper functions are required by the command blocks below. They are shown here directly so the report remains self-contained and the code is visible at first use.",
      "",
      "```{r, results='asis'}",
      "cat('```r\\n')",
      "cat(paste(readLines('report_methods.R', warn = FALSE), collapse = '\\n'))",
      "cat('\\n```\\n')",
      "```",
      "",
      "```{r, results='asis'}",
      "cat('```r\\n')",
      "cat(paste(readLines('analysis_commands.R', warn = FALSE), collapse = '\\n'))",
      "cat('\\n```\\n')",
      "```",
      "",
      "# Environment summary",
      "",
      "These package versions were used by the app instance that generated this bundle. Reproducing the workflow on another computer is easiest when the same major package versions are installed first.",
      "",
      "```{r}",
      "snap$package_versions",
      "```",
      "",
      "# Actions run in the app",
      "",
      "This is the chronological analysis log captured from the Shiny workflow. Use it as the audit trail showing which steps were actually run before you inspect the detailed code and plots below.",
      "",
      "```{r}",
      "log_rows",
      "```",
      "",
      "# Input settings",
      "",
      "These are the dataset-level settings used when the data were loaded.",
      "",
      "```{r}",
      "str(snap$inputs)",
      "```"
    )
    if (!is.null(snapshot$qc)) {
      rmd_lines <- c(rmd_lines,
        "",
        "# QC",
        "",
        "## Method meaning",
        "",
        "Quality control is the first pass for checking whether sample-level behavior looks sensible before model fitting. In this workflow, QC focuses on library size, detected genes, sample clustering, and PCA structure.",
        "",
        "## Chosen arguments",
        "",
        "```{r}",
        "show_params(last_params('load_data'))",
        "```",
        "",
        "Argument interpretation:",
        "",
        "- `min_count`: minimum count threshold used for gene retention.",
        "- `min_samples`: how many samples must pass that threshold for a gene to be kept.",
        "- `qc_group`: metadata column used to color the sample-level plots.",
        "",
        "How to use the elements below:",
        "",
        "- The table is the sample-level reference sheet for checking sequencing depth, detected genes, and metadata alignment.",
        "- The QC plot is meant to reveal obvious sample imbalance before any statistical modeling. Colors come from the chosen QC grouping, so groups should be read as metadata categories rather than significance labels.",
        "- The code blocks below can be copied into another R session to rebuild the same QC objects from the same inputs and thresholds.",
        "",
        "## Stepwise recreation",
        "",
        "### QC table",
        "",
        "```{r}",
        "head(snap$qc$sample_table)",
        "```",
        "",
        "### Library-size plot data",
        "",
        "```{r}",
        "head(snap$qc$lib_df)",
        "```",
        "",
        "### Recreate the library-size plot",
        "",
        "```{r}",
        "plot_ly(snap$qc$lib_df, x = ~run_id, y = ~library_size, color = ~group, type = 'bar')",
        "```",
        "",
        "### PCA score table used for plotting",
        "",
        "```{r}",
        "head(snap$qc$pca_df)",
        "```",
        "",
        "### Recreate the PCA plot",
        "",
        "```{r}",
        "plot_ly(snap$qc$pca_df, x = ~PC1, y = ~PC2, color = ~group, type = 'scatter', mode = 'markers+text', text = ~run_id, textposition = 'top center')",
        "```"
      )
    }
    if (!is.null(snapshot$de)) {
      rmd_lines <- c(rmd_lines,
        "",
        "# DESeq2",
        "",
        "## Method meaning",
        "",
        "DESeq2 models gene-wise count variation to estimate differential expression between conditions while accounting for the selected design formula. The key outputs here are `log2FoldChange` and `padj`.",
        "",
        "## Chosen arguments",
        "",
        "```{r}",
        "show_params(last_params('run_deseq2'))",
        "```",
        "",
        "Argument interpretation:",
        "",
        "- `active_design`: the saved design formula chosen by the user.",
        "- `comparison_variable`: the metadata column used for the contrast.",
        "- `numerator` and `denominator`: the direction of the fold change.",
        "- `selected_samples`: which samples entered the DESeq2 fit.",
        "",
        "How to use the elements below:",
        "",
        "- The DE table is the main ranked result. It should be used to inspect effect size, adjusted significance, and gene identity together rather than relying on p-values alone.",
        "- In the volcano-style plot, point color separates upregulated DEGs, downregulated DEGs, and genes outside the current DEG thresholds. The vertical and horizontal guide lines mark the exact log2FC and adjusted-significance cutoffs used in the app.",
        "- The code blocks below are copy-ready: run the helper functions first, then the DESeq2 code to recreate the design, contrast, DEG calls, and plot styling in another session.",
        "",
        "## Stepwise recreation",
        "",
        "### DE result table",
        "",
        "```{r}",
        "head(snap$de$table[, c('gene_label','log2FoldChange','padj')])",
        "```",
        "",
        "### Rebuild the volcano-style diagnostic table",
        "",
        "```{r}",
        "d <- snap$de$table",
        "d$score <- -log10(pmax(d$padj, .Machine$double.xmin))",
        "head(d[, c('gene_label','log2FoldChange','score')])",
        "```",
        "",
        "### Recreate the volcano-style plot",
        "",
        "```{r}",
        "plot_ly(d, x = ~log2FoldChange, y = ~score, type = 'scatter', mode = 'markers', text = ~gene_label, hoverinfo = 'text')",
        "```",
        "",
        "### Size factors used for normalization diagnostics",
        "",
        "```{r}",
        "snap$de$size_factors",
        "```"
      )
    }
    if (!is.null(snapshot$enrichment)) {
      rmd_lines <- c(rmd_lines,
        "",
        "# Pathway enrichment",
        "",
        "## Method meaning",
        "",
        "The enrichment stage summarizes the DE result in terms of pathways. ORA tests a selected DE gene subset, while GSEA evaluates the full ranked list. Consensus here means overlap across methods, not a recomputed merged p-value.",
        "",
        "## Chosen arguments",
        "",
        "```{r}",
        "show_params(last_params('run_enrichment'))",
        "```",
        "",
        "Argument interpretation:",
        "",
        "- `methods`: enrichment algorithms selected in the app.",
        "- `padj_cutoff`: threshold used to define significant pathway results.",
        "- `abs_log2fc_cutoff`: ORA gene-selection effect-size cutoff.",
        "- `min_pathway_size` and `max_pathway_size`: pathway size bounds.",
        "- `adjust_method`: multiple-testing correction method.",
        "",
        "How to use the elements below:",
        "",
        "- The consensus table is for interpretation across methods, not for re-running the statistical test itself.",
        "- The pathway plot uses point size and position to summarize pathway prominence and recurrence across methods. Treat it as a prioritization view, then use the CSV for the exact values.",
        "- The code blocks below can be reused to rebuild the same ORA/GSEA inputs and apply the same pathway arguments in another R session.",
        "",
        "## Stepwise recreation",
        "",
        "### Consensus pathway table",
        "",
        "```{r}",
        "head(snap$enrichment$consensus)",
        "```",
        "",
        "### Recreate the consensus plot",
        "",
        "```{r}",
        "top_cons <- head(snap$enrichment$consensus, 12)",
        "plot_ly(top_cons, x = ~consensus_hits, y = ~pathway, type = 'scatter', mode = 'markers', text = ~pathway, hoverinfo = 'text')",
        "```",
        "",
        "### Inspect the combined per-method table",
        "",
        "```{r}",
        "head(snap$enrichment$combined)",
        "```"
      )
    }
    if (!is.null(snapshot$wgcna)) {
      rmd_lines <- c(rmd_lines,
        "",
        "# WGCNA",
        "",
        "## Method meaning",
        "",
        "WGCNA identifies co-expression modules, then summarizes them through module eigengenes for trait association and downstream biology. In this app, WGCNA is also used as the starting graph prior for later SEMgraph analysis.",
        "",
        "## Chosen arguments",
        "",
        "```{r}",
        "show_params(last_params('run_wgcna'))",
        "```",
        "",
        "Argument interpretation:",
        "",
        "- `power_mode` and `manual_power`: how the soft-threshold was chosen.",
        "- `top_variable_genes`: how many genes entered network construction.",
        "- `min_module_size`, `deep_split`, `merge_cut_height`: module granularity controls.",
        "- `traits`: metadata variables used for module-trait correlation.",
        "",
        "How to use the elements below:",
        "",
        "- The module-size and module-membership outputs summarize the network structure and should be interpreted as coexpression organization, not as differential-expression testing.",
        "- In the exported module network plot, labels are shown as gene symbols and the visual styling mirrors the app view so hubs and stronger within-module edges remain recognizable.",
        "- The code blocks below can be copied into another R session to rebuild the variable-gene set, soft-thresholding, TOM-based clustering, and module assignments.",
        "",
        "## Stepwise recreation",
        "",
        "### Module-size table",
        "",
        "```{r}",
        "snap$wgcna$module_sizes",
        "```"
      )
      if (!is.null(snapshot$wgcna$module_graph)) {
        rmd_lines <- c(rmd_lines,
          "",
          "### Selected module graph snapshot",
          "",
          "```{r}",
          "head(snap$wgcna$module_graph$nodes)",
          "head(snap$wgcna$module_graph$edges)",
          "```",
          "",
          "### Recreate the selected WGCNA module graph",
          "",
          "```{r}",
          "nd <- snap$wgcna$module_graph$nodes; ed <- snap$wgcna$module_graph$edges",
          "p <- plot_ly()",
          "if (nrow(ed) > 0) {",
          "  for (i in seq_len(nrow(ed))) {",
          "    p <- add_segments(p,",
          "      x = nd$V1[match(ed$from[i], nd$gene_id)], y = nd$V2[match(ed$from[i], nd$gene_id)],",
          "      xend = nd$V1[match(ed$to[i], nd$gene_id)], yend = nd$V2[match(ed$to[i], nd$gene_id)],",
          "      inherit = FALSE, line = list(color = 'rgba(120,120,120,0.35)', width = 1))",
          "  }",
          "}",
          "add_trace(p, data = nd, x = ~V1, y = ~V2, type = 'scatter', mode = 'markers+text', text = ~gene_id, textposition = 'top center')",
          "```"
        )
      }
    }
    if (!is.null(snapshot$semgraph)) {
      rmd_lines <- c(rmd_lines,
        "",
        "# SEMgraph",
        "",
        "## Method meaning",
        "",
        "SEMgraph starts from a graph prior and orients/model-fits it to generate a directed hypothesis network. In this workflow, those directions are model-based causal candidates, not direct proof from the RNA-seq data alone.",
        "",
        "## Chosen arguments",
        "",
        "```{r}",
        "show_params(last_params('run_semgraph'))",
        "```",
        "",
        "Argument interpretation:",
        "",
        "- `analysis_unit`: original module or combined program used as the prior gene set.",
        "- `group`: optional binary variable used for local perturbation-style causal effect estimation.",
        "- `display_mode`: how strongly DE-supported genes are prioritized in the shown graph.",
        "- `max_model_genes`: fit-size cap applied before SEMgraph.",
        "- `max_annotation_genes`: cap for slower external annotation layers.",
        "",
        "How to use the elements below:",
        "",
        "- The regulator and edge tables are the precise reference for the currently exported SEMgraph view. Use them when you want the exact values behind the graph styling.",
        "- In the SEMgraph plot, node fill shows DE direction and magnitude, node role styling reflects causal position, gold halos mark druggability context, and edge colors reflect STRING support strength while SEMgraph still provides arrow direction.",
        "- The code blocks below can be reused to rebuild the selected SEMgraph setup and reproduce the same directed hypothesis layer, provided the helper functions and upstream objects are defined first.",
        "",
        "## Stepwise recreation",
        "",
        "### Degree table",
        "",
        "```{r}",
        "head(snap$semgraph$degree_table)",
        "```"
      )
      if (!is.null(snapshot$semgraph$graph)) {
        rmd_lines <- c(rmd_lines,
          "",
          "### Displayed SEMgraph subgraph snapshot",
          "",
          "```{r}",
          "head(snap$semgraph$graph$nodes)",
          "head(snap$semgraph$graph$edges)",
          "```",
          "",
          "### Recreate the displayed SEMgraph subgraph",
          "",
          "```{r}",
          "nd <- snap$semgraph$graph$nodes; ed <- snap$semgraph$graph$edges",
          "p <- plot_ly()",
          "if (nrow(ed) > 0) {",
          "  for (i in seq_len(nrow(ed))) {",
          "    p <- add_segments(p,",
          "      x = nd$x[match(ed$from[i], nd$gene)], y = nd$y[match(ed$from[i], nd$gene)],",
          "      xend = nd$x[match(ed$to[i], nd$gene)], yend = nd$y[match(ed$to[i], nd$gene)],",
          "      inherit = FALSE, line = list(color = 'rgba(120,120,120,0.4)', width = 1))",
          "  }",
          "}",
          "add_trace(p, data = nd, x = ~x, y = ~y, type = 'scatter', mode = 'markers+text', text = ~gene, textposition = 'top center')",
          "```"
        )
      }
      rmd_lines <- c(rmd_lines,
        "",
        "### Effect or regulator-priority table",
        "",
        "```{r}",
        "head(snap$semgraph$ace_table)",
        "```",
        "",
        "### Pathway-context summary",
        "",
        "```{r}",
        "head(snap$semgraph$pathway_summary)",
        "```"
      )
    }
    write_text_file(file.path(out_dir, "report.Rmd"), rmd_lines)
    helper_lines <- c(
      "library(jsonlite)",
      "library(plotly)",
      "",
      "message('For rerunning the workflow with explicit commands, see analysis_commands.R and report_methods.R in this bundle.')",
      "",
      "snap <- readRDS('report_snapshot.rds')",
      "logs <- snap$analysis_log",
      "parse_params <- function(x) {",
      "  if (is.null(x) || is.na(x) || !nzchar(x)) return(list())",
      "  jsonlite::fromJSON(x, simplifyVector = FALSE)",
      "}",
      "last_params <- function(step_name) {",
      "  idx <- which(logs$step == step_name)",
      "  if (!length(idx)) return(list())",
      "  parse_params(logs$params_json[max(idx)])",
      "}",
      "",
      "print('Package versions:')",
      "print(snap$package_versions)",
      "print('Analysis log:')",
      "print(logs)",
      "",
      "if (!is.null(snap$qc)) {",
      "  print('QC sample table:')",
      "  print(utils::head(snap$qc$sample_table))",
      "}",
      "if (!is.null(snap$de)) {",
      "  print('DESeq2 result table:')",
      "  print(utils::head(snap$de$table[, c('gene_label','log2FoldChange','padj')]))",
      "}",
      "if (!is.null(snap$enrichment)) {",
      "  print('Enrichment consensus:')",
      "  print(utils::head(snap$enrichment$consensus))",
      "}",
      "if (!is.null(snap$wgcna)) {",
      "  print('WGCNA module sizes:')",
      "  print(snap$wgcna$module_sizes)",
      "}",
      "if (!is.null(snap$semgraph)) {",
      "  print('SEMgraph degree table:')",
      "  print(utils::head(snap$semgraph$degree_table))",
      "}"
    )
    write_text_file(file.path(out_dir, "recreate_from_snapshot.R"), helper_lines)

    html_lines <- c(
      "<!doctype html><html><head><meta charset='utf-8'><title>RNA-seq Workbench HTML Report</title>",
      "<style>body{font-family:Arial,sans-serif;background:#f6f7fb;color:#1f2937;margin:0;padding:24px;} .card{background:#fff;border:1px solid #d7dbe8;border-radius:10px;padding:16px;margin-bottom:16px;} iframe{width:100%;height:700px;border:1px solid #d7dbe8;border-radius:8px;background:#fff;} table{border-collapse:collapse;width:100%;} th,td{border:1px solid #d7dbe8;padding:6px 8px;text-align:left;} h1,h2,h3{margin-top:0;} a{color:#0b7285;} .method-box{margin-top:14px;padding:14px;border:1px solid #dbe4f0;border-radius:8px;background:#f8fbff;} .method-box h4{margin-bottom:8px;} .code-dropdown{margin-top:12px;border:1px solid #d7dbe8;border-radius:8px;background:#fbfcfe;} .code-dropdown summary{cursor:pointer;padding:10px 12px;font-weight:700;color:#334155;} .code-dropdown[open] summary{border-bottom:1px solid #e2e8f0;} .code-note{margin:10px 12px 0 12px;color:#556178;} .code-toolbar{display:flex;justify-content:flex-end;padding:8px 12px 0 12px;} .copy-btn{background:#0b7285;color:#fff;border:none;border-radius:6px;padding:6px 10px;cursor:pointer;font-size:12px;} .copy-btn.copied{background:#2f9e44;} .code-block{margin:8px 12px 12px 12px;padding:12px;background:#0f172a;color:#e2e8f0;border-radius:8px;overflow:auto;} .code-block code{font-family:Consolas,Menlo,monospace;font-size:12px;white-space:pre;} .report-links a{display:inline-block;margin-right:12px;margin-top:6px;}</style>",
      "<script>function copyCodeBlock(id, btn){var el=document.getElementById(id);if(!el)return;var txt=el.innerText||el.textContent||'';navigator.clipboard.writeText(txt).then(function(){if(btn){btn.classList.add('copied');btn.textContent='Copied';setTimeout(function(){btn.classList.remove('copied');btn.textContent='Copy code';},1500);}});}</script>",
      "</head><body>",
      "<h1>RNA-seq Workbench HTML Report</h1>",
      sprintf("<p>Generated at %s.</p>", htmltools::htmlEscape(snapshot$generated_at)),
      "<div class='card'><h2>Reproducibility files</h2><p>This bundle includes both the saved outputs and the code needed to rebuild them. Use <code>report_methods.R</code> first for helper definitions, then <code>analysis_commands.R</code> for the session-specific workflow commands.</p><div class='report-links'><a href='analysis_commands.R'>Open analysis_commands.R</a><a href='report_methods.R'>Open report_methods.R</a></div></div>",
      "<div class='card'><h2>Core helper functions</h2><p>The generated report code depends on helper functions such as <code>read_inputs()</code>, <code>deg_mask()</code>, and <code>run_de()</code>. They are shown here directly so the report remains self-contained and the copied code is understandable at first use.</p>",
      html_code_panel("code_report_methods", "Show required helper functions", report_methods_lines, "Reuse guide: copy or source these helpers first in a new R session. The later workflow code boxes assume these functions already exist."),
      "</div>",
      "<div class='card'><h2>Generated workflow commands</h2><p>This is the main copy-ready command sequence produced from the Shiny session. It captures the actual argument values chosen in the app and is meant to be rerun after the helper functions are defined.</p>",
      html_code_panel("code_analysis_commands", "Show generated workflow commands", c(readLines(file.path(out_dir, "analysis_commands.R"), warn = FALSE)), "Reuse guide: run the helper functions first, then execute these commands from top to bottom in a fresh R session. Adjust argument values here if you want to adapt the workflow while keeping the same method logic."),
      "</div>",
      "<div class='card'><h2>Package versions</h2><p>This file captures the package versions needed to recreate the environment.</p><p><a href='package_versions.csv'>Download CSV</a></p></div>",
      "<div class='card'><h2>Analysis log</h2><p>This file records which steps were run in the app and with which saved arguments.</p><p><a href='analysis_log.csv'>Download CSV</a></p></div>"
    )
    if (file.exists(file.path(plots_dir, "qc_library_sizes.html"))) {
      html_lines <- c(html_lines,
        "<div class='card'><h2>QC plot</h2><p>This plot shows sample-level library sizes, colored by the selected QC grouping.</p><iframe src='plots/qc_library_sizes.html'></iframe>",
        html_method_card(
          "How to read this QC step",
          c(
            "Quality control is the first pass for checking whether the sequencing depth and sample-level structure look plausible before any downstream model fitting.",
            "In this workflow, low-count genes are removed first, then sample behavior is inspected through library-size summaries, PCA, and distance-based clustering."
          ),
          c(
            sprintf("Minimum count threshold: %s. This removes genes with very weak evidence of expression.", snapshot$inputs$min_count %||% load_params$min_count %||% 10),
            sprintf("Minimum samples at threshold: %s. A gene is retained only if it reaches the minimum count in at least this many samples.", snapshot$inputs$min_samples %||% load_params$min_samples %||% 2),
            sprintf("QC grouping variable: %s. This affects how samples are colored in the exported QC plots.", snapshot$inputs$qc_group %||% "treatment")
          ),
          c(
            "Use the plot to check whether some samples are obviously different in sequencing depth before trusting downstream statistics.",
            "Plot colors indicate the chosen QC grouping only; they are for interpretation, not significance.",
            "Use the CSV when you want exact numeric sample values, and use the code dropdown to rebuild the same QC objects elsewhere."
          )
        ),
        html_code_panel("code_qc", "Show copy-ready QC code", qc_code_lines, "This block reloads the counts and metadata, applies the chosen QC thresholds, and rebuilds the sample-level QC objects. Copy it into a fresh R session after defining the helper functions."),
        "<p><a href='tables/qc_sample_table.csv'>Open QC sample table CSV</a></p></div>"
      )
    }
    if (file.exists(file.path(plots_dir, "deseq2_volcano.html"))) {
      html_lines <- c(html_lines,
        "<div class='card'><h2>DESeq2 plot</h2><p>This plot summarizes differential expression by effect size and adjusted significance.</p><iframe src='plots/deseq2_volcano.html'></iframe>",
        html_method_card(
          "How to read this DESeq2 step",
          c(
            "DESeq2 models gene-wise count variation to estimate differential expression between the chosen conditions while accounting for the selected design formula.",
            "The volcano plot combines effect size and adjusted significance. Upregulated and downregulated DEGs are separated so direction is easier to interpret."
          ),
          c(
            sprintf("Design formula: %s. This determines which metadata terms DESeq2 adjusts for before testing the selected contrast.", design_formula),
            sprintf("Contrast variable: %s; numerator: %s; denominator: %s. These determine the direction of the reported log2 fold changes.", de_params$comparison_variable %||% "NA", de_params$numerator %||% "NA", de_params$denominator %||% "NA"),
            sprintf("Adjusted p-value cutoff: %s and absolute log2FC cutoff: %s. These two thresholds define which genes are currently treated as DEGs in the exported plot.", snapshot$inputs$de_alpha %||% de_params$alpha %||% 0.05, snapshot$inputs$de_lfc_cutoff %||% de_params$abs_log2fc_cutoff %||% 0)
          ),
          c(
            "Use the plot as a quick map of effect size versus adjusted significance, then confirm exact values in the CSV table.",
            "Point colors separate upregulated DEGs, downregulated DEGs, and genes outside the current thresholds. The guide lines mark the exact cutoffs used in the app.",
            "The code dropdown recreates both the DEG calls and the styled volcano-ready data, so it can be copied into another R session."
          )
        ),
        html_code_panel("code_de", "Show copy-ready DESeq2 code", de_code_lines, "This block rebuilds the active DESeq2 design, uses the chosen samples and contrast, and regenerates the DEG table plus the styled volcano input data. Copy it into a fresh R session after the helper functions are defined."),
        "<p><a href='tables/deseq2_results.csv'>Open DESeq2 results CSV</a></p></div>"
      )
    }
    if (file.exists(file.path(plots_dir, "enrichment_consensus.html"))) {
      html_lines <- c(html_lines,
        "<div class='card'><h2>Pathway enrichment plot</h2><p>This view shows the overlap-supported pathways returned by the selected enrichment methods.</p><iframe src='plots/enrichment_consensus.html'></iframe>",
        html_method_card(
          "How to read this pathway-enrichment step",
          c(
            "ORA tests whether a selected DEG subset is overrepresented in known pathways, while GSEA evaluates whether pathway members collectively drift toward one end of the ranked DE result.",
            "Consensus in this report means recurrence across methods, not a merged p-value."
          ),
          c(
            sprintf("Methods selected: %s.", paste(enrich_params$methods %||% character(), collapse = ", ")),
            sprintf("Pathway significance cutoff: %s; absolute log2FC cutoff for ORA gene selection: %s.", enrich_params$padj_cutoff %||% snapshot$enrichment$settings$padj %||% 0.05, enrich_params$abs_log2fc_cutoff %||% snapshot$enrichment$settings$abs_log2fc %||% 0),
            sprintf("Pathway size limits: min %s, max %s. This controls how broad or narrow the tested pathway sets can be.", enrich_params$min_pathway_size %||% snapshot$enrichment$settings$min_gs %||% 10, enrich_params$max_pathway_size %||% snapshot$enrichment$settings$max_gs %||% 500)
          ),
          c(
            "Use the plot to prioritize pathway themes, then use the CSV for the exact pathway statistics and member lists.",
            "Point size and ordering summarize pathway prominence in the exported view; they are not a substitute for reading the underlying adjusted p-values.",
            "The code dropdown rebuilds the same ORA and GSEA inputs with the same settings, so it is ready to reuse elsewhere."
          )
        ),
        html_code_panel("code_enrich", "Show copy-ready enrichment code", enrich_code_lines, "This block rebuilds the ORA/GSEA inputs from the current DE result and applies the same enrichment-method arguments chosen in the app. Copy it into another R session after the helper functions are defined."),
        "<p><a href='tables/enrichment_consensus.csv'>Open enrichment consensus CSV</a></p></div>"
      )
    }
    if (file.exists(file.path(plots_dir, "wgcna_module_network.html"))) {
      html_lines <- c(html_lines,
        "<div class='card'><h2>WGCNA module network</h2><p>This graph is the selected module-level coexpression network snapshot exported from the app.</p><iframe src='plots/wgcna_module_network.html'></iframe>",
        html_method_card(
          "How to read this WGCNA step",
          c(
            "WGCNA identifies co-expression modules, then summarizes each module as a coherent expression program whose central genes may serve as useful biological anchors.",
            "This network view focuses on one selected module so the strongest internal co-expression structure can be inspected without the rest of the network overwhelming the plot."
          ),
          c(
            sprintf("Most variable genes included: %s. This limits WGCNA to the genes with the strongest expression variation across samples.", wgcna_params$top_variable_genes %||% 4000),
            sprintf("Soft-threshold mode: %s; manual power fallback: %s. This shapes how strongly pairwise correlations are transformed into network connectivity.", wgcna_params$power_mode %||% "auto", wgcna_params$manual_power %||% 6),
            sprintf("Minimum module size: %s; deep split: %s; merge cut height: %s. These control how easily modules are split or merged.", wgcna_params$min_module_size %||% 30, wgcna_params$deep_split %||% 2, wgcna_params$merge_cut_height %||% 0.2),
            sprintf("Selected module for this export: %s. Gene labels in the plot are shown as HGNC-style labels rather than raw Ensembl identifiers.", snapshot$wgcna$selected_module %||% "NA")
          ),
          c(
            "Use this network plot to inspect module structure and hub-like genes, not to read statistical significance directly.",
            "The exported styling mirrors the app view so gene-symbol labels, stronger edges, and module-centered structure remain visually consistent.",
            "The code dropdown reconstructs the variable-gene set, TOM network, clustering, and module definitions used in the app."
          )
        ),
        html_code_panel("code_wgcna", "Show copy-ready WGCNA code", wgcna_code_lines, "This block reconstructs the variable-gene selection, soft-thresholding, TOM network, dendrogram cut, and module merge settings used in the app. Copy it into a fresh R session after the helper functions are defined."),
        "<p><a href='tables/wgcna_gene_modules.csv'>Open WGCNA gene module CSV</a></p></div>"
      )
    }
    if (file.exists(file.path(plots_dir, "semgraph_subgraph.html"))) {
      html_lines <- c(html_lines,
        "<div class='card'><h2>SEMgraph subgraph</h2><p>This graph is the displayed SEMgraph view exported from the app after all current filtering rules were applied, including any STRING edge-validation filter.</p><iframe src='plots/semgraph_subgraph.html'></iframe>",
        html_method_card(
          "How to read this SEMgraph step",
          c(
            "SEMgraph uses the selected module or combined program as a graph-informed prior, orients that graph under structural assumptions, and fits a causal-style model to prioritize upstream and downstream relationships.",
            "STRING support in this view validates biological relatedness of a pair, while SEMgraph still provides the direction hypothesis. When STRING action data are available, the report also notes whether a directed action agrees with the SEMgraph arrow."
          ),
          c(
            sprintf("Analysis unit: %s. This determines which WGCNA module or combined program was used as the SEMgraph prior set.", semgraph_params$analysis_unit %||% snapshot$semgraph$analysis_unit %||% "NA"),
            sprintf("Display mode: %s. This controls whether the graph view stays strictly DEG-focused or allows connector genes for readability.", semgraph_params$display_mode %||% "de_connectors"),
            sprintf("Maximum genes entering SEMgraph fit: %s; maximum genes receiving external annotation: %s. These settings cap computational cost while prioritizing the most relevant genes.", semgraph_params$max_model_genes %||% 300, semgraph_params$max_annotation_genes %||% 180)
          ),
          c(
            "Use the graph for directional hypothesis inspection, then use the regulator and edge CSV tables for exact numeric support.",
            "Node styling carries multiple meanings at once: DE support, causal role, and druggability context. Edge colors summarize STRING support strength, while arrow direction still comes from SEMgraph.",
            "The code dropdown captures the chosen SEMgraph setup and can be reused after the helper functions and upstream DE/WGCNA objects are already defined."
          )
        ),
        html_code_panel("code_semgraph", "Show copy-ready SEMgraph code", semgraph_code_lines, "This block records the chosen SEMgraph unit and key fit arguments. The exported graph and validation tables capture the exact visible SEMgraph state rendered in the app. Copy it into another R session after the helper functions and upstream objects are available."),
        "<p><a href='tables/semgraph_degree_table.csv'>Open SEMgraph degree table CSV</a></p><p><a href='tables/semgraph_edge_validation_table.csv'>Open SEMgraph edge validation CSV</a></p></div>"
      )
    }
    html_lines <- c(html_lines, "</body></html>")
    write_text_file(file.path(out_dir, "index.html"), html_lines)
  }

  output$download_repro_bundle <- shiny::downloadHandler(
    filename = function() sprintf("rnaseq_reproducibility_bundle_%s.zip", format(Sys.time(), "%Y%m%d_%H%M%S")),
    content = function(file) {
      bundle_dir <- file.path(tempdir(), paste0("repro_bundle_", format(Sys.time(), "%Y%m%d_%H%M%S")))
      build_report_bundle(current_report_snapshot(), bundle_dir, html_mode = FALSE)
      zip_dir_contents(file, bundle_dir)
    }
  )
  output$download_html_bundle <- shiny::downloadHandler(
    filename = function() sprintf("rnaseq_html_report_bundle_%s.zip", format(Sys.time(), "%Y%m%d_%H%M%S")),
    content = function(file) {
      bundle_dir <- file.path(tempdir(), paste0("html_bundle_", format(Sys.time(), "%Y%m%d_%H%M%S")))
      build_report_bundle(current_report_snapshot(), bundle_dir, html_mode = TRUE)
      zip_dir_contents(file, bundle_dir)
    }
  )
}

message("Starting RNA Hypothesis Workbench on http://127.0.0.1:3838")
shiny::runApp(shiny::shinyApp(ui = ui, server = server), host = "127.0.0.1", port = 3838, launch.browser = TRUE)
