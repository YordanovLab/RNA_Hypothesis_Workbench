#!/usr/bin/env Rscript

required_packages <- c("shiny", "plotly", "DT")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
  message("Missing required packages: ", paste(missing_packages, collapse = ", "))
  message("Install them with:")
  message("R -e \"install.packages(c(", paste(sprintf("'%s'", missing_packages), collapse = ", "), "))\"")
  quit(status = 1L)
}

data_dir <- file.path(getwd(), "data")
default_counts <- file.path(data_dir, "GSE267213_counts.csv")
default_metadata <- file.path(data_dir, "SraRunTable.csv")

read_counts_file <- function(path) {
  counts_df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(counts_df) < 2L) {
    stop("Counts file must contain a gene identifier column and at least one sample column.")
  }

  gene_col <- names(counts_df)[1L]
  gene_ids <- counts_df[[gene_col]]
  count_matrix <- as.matrix(counts_df[, -1L, drop = FALSE])
  storage.mode(count_matrix) <- "numeric"
  rownames(count_matrix) <- gene_ids

  if (anyNA(count_matrix)) {
    stop("Counts file contains missing or non-numeric values in sample columns.")
  }

  list(
    gene_column = gene_col,
    gene_ids = gene_ids,
    counts = count_matrix,
    samples = colnames(count_matrix)
  )
}

read_metadata_file <- function(path) {
  metadata <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"run_id" %in% names(metadata)) {
    stop("Metadata file must contain a 'run_id' column.")
  }
  metadata
}

align_metadata <- function(sample_ids, metadata) {
  aligned <- merge(
    data.frame(run_id = sample_ids, sample_order = seq_along(sample_ids), stringsAsFactors = FALSE),
    metadata,
    by = "run_id",
    all.x = TRUE,
    sort = FALSE
  )
  aligned <- aligned[order(aligned$sample_order), , drop = FALSE]
  aligned$sample_order <- NULL

  annotation_cols <- setdiff(names(aligned), "run_id")
  if (length(annotation_cols) == 0L) {
    aligned$annotation_status <- "Unannotated"
  }

  for (column in setdiff(names(aligned), "run_id")) {
    aligned[[column]][is.na(aligned[[column]]) | aligned[[column]] == ""] <- "Unannotated"
  }

  aligned
}

compute_lib_sizes <- function(counts) {
  colSums(counts)
}

compute_detected_genes <- function(counts) {
  colSums(counts > 0)
}

filter_counts <- function(counts, min_count = 10, min_samples = 2) {
  keep <- rowSums(counts >= min_count) >= min_samples
  counts[keep, , drop = FALSE]
}

compute_log_cpm <- function(counts, prior_count = 1) {
  lib_sizes <- compute_lib_sizes(counts)
  normalized <- sweep(counts, 2L, lib_sizes / 1e6, "/")
  log2(normalized + prior_count)
}

compute_pca <- function(log_cpm) {
  pca <- prcomp(t(log_cpm), center = TRUE, scale. = TRUE)
  variance <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2)
  list(scores = pca$x, variance = variance)
}

compute_sample_distance <- function(log_cpm) {
  as.matrix(dist(t(log_cpm)))
}

summarize_dataset <- function(counts, filtered_counts, metadata) {
  list(
    n_genes = nrow(counts),
    n_filtered_genes = nrow(filtered_counts),
    n_samples = ncol(counts),
    groups = if ("treatment" %in% names(metadata)) table(metadata$treatment) else NULL
  )
}

make_sample_qc <- function(counts, metadata) {
  lib_sizes <- compute_lib_sizes(counts)
  detected <- compute_detected_genes(counts)
  qc <- data.frame(
    run_id = colnames(counts),
    library_size = as.numeric(lib_sizes),
    detected_genes = as.integer(detected),
    stringsAsFactors = FALSE
  )
  merge(qc, metadata, by = "run_id", all.x = TRUE, sort = FALSE)
}

make_group_colors <- function(groups) {
  unique_groups <- unique(as.character(groups))
  palette <- c(
    "#1f77b4", "#d62728", "#2ca02c", "#ff7f0e",
    "#9467bd", "#8c564b", "#e377c2", "#17becf"
  )
  colors <- palette[seq_along(unique_groups)]
  names(colors) <- unique_groups
  colors
}

make_colored_ticktext <- function(labels, groups, colors) {
  vapply(
    seq_along(labels),
    function(i) {
      sprintf(
        "<span style='color:%s;'>%s</span>",
        colors[[as.character(groups[[i]])]],
        labels[[i]]
      )
    },
    character(1)
  )
}

ui <- shiny::fluidPage(
  shiny::tags$head(
    shiny::tags$style(shiny::HTML("
      .help-badge {
        display: inline-block;
        width: 18px;
        height: 18px;
        border-radius: 50%;
        background: #4c566a;
        color: #ffffff;
        text-align: center;
        font-size: 12px;
        line-height: 18px;
        cursor: help;
        margin-left: 6px;
      }
      .legend-wrap {
        padding-top: 8px;
      }
      .legend-item {
        display: flex;
        align-items: center;
        gap: 8px;
        margin-bottom: 6px;
      }
      .legend-swatch {
        width: 12px;
        height: 12px;
        border-radius: 50%;
        display: inline-block;
      }
    "))
  ),
  shiny::titlePanel("RNA-seq Exploratory Analysis Launcher"),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::textInput("counts_path", "Counts file", value = default_counts),
      shiny::textInput("metadata_path", "Metadata file", value = default_metadata),
      shiny::actionButton("load_data", "Load dataset"),
      shiny::tags$hr(),
      shiny::numericInput("min_count", "Minimum count per gene", value = 10, min = 0, step = 1),
      shiny::numericInput("min_samples", "Minimum samples meeting threshold", value = 2, min = 1, step = 1),
      shiny::selectInput("group_column", "Annotation column", choices = "treatment", selected = "treatment"),
      shiny::tags$p("Use the controls above to prepare normalized exploratory views before differential expression analysis.")
    ),
    shiny::mainPanel(
      shiny::verbatimTextOutput("status"),
      shiny::tabsetPanel(
        shiny::tabPanel(
          "Overview",
          shiny::uiOutput("dataset_summary"),
          DT::dataTableOutput("sample_table")
        ),
        shiny::tabPanel(
          "Sample QC",
          plotly::plotlyOutput("library_plot", height = "320px"),
          plotly::plotlyOutput("detected_plot", height = "320px")
        ),
        shiny::tabPanel(
          "PCA",
          plotly::plotlyOutput("pca_plot", height = "500px")
        ),
        shiny::tabPanel(
          "Distances",
          shiny::fluidRow(
            shiny::column(
              width = 9,
              shiny::div(
                shiny::tags$strong("Sample Distance Heatmap"),
                shiny::tags$span(
                  "?",
                  class = "help-badge",
                  title = "This heatmap shows pairwise Euclidean distances between samples using filtered log2(CPM + 1) values. Lower distances suggest more similar expression profiles. Axis label colors reflect the selected metadata grouping."
                )
              ),
              plotly::plotlyOutput("distance_heatmap", height = "500px")
            ),
            shiny::column(
              width = 3,
              shiny::div(class = "legend-wrap",
                shiny::tags$strong("Group colors"),
                shiny::uiOutput("distance_legend")
              )
            )
          )
        ),
        shiny::tabPanel(
          "Distributions",
          plotly::plotlyOutput("box_plot", height = "420px"),
          plotly::plotlyOutput("mean_variance_plot", height = "420px")
        ),
        shiny::tabPanel(
          "Gene Explorer",
          shiny::selectizeInput("gene_choice", "Gene ID", choices = NULL, selected = NULL, options = list(maxOptions = 5000)),
          plotly::plotlyOutput("gene_plot", height = "420px")
        ),
        shiny::tabPanel(
          "DE Prep Notes",
          shiny::uiOutput("prep_notes")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  dataset <- shiny::eventReactive(input$load_data, {
    counts_info <- read_counts_file(input$counts_path)
    metadata <- read_metadata_file(input$metadata_path)
    metadata <- align_metadata(counts_info$samples, metadata)

    annotation_choices <- setdiff(names(metadata), "run_id")
    if (length(annotation_choices) == 0L) {
      annotation_choices <- "run_id"
    }
    shiny::updateSelectInput(session, "group_column", choices = annotation_choices, selected = annotation_choices[1L])

    filtered_counts <- filter_counts(
      counts_info$counts,
      min_count = input$min_count,
      min_samples = input$min_samples
    )
    if (nrow(filtered_counts) == 0L) {
      stop("Filtering removed all genes. Lower the thresholds and reload the dataset.")
    }

    log_cpm <- compute_log_cpm(filtered_counts)
    pca <- compute_pca(log_cpm)
    distances <- compute_sample_distance(log_cpm)
    sample_qc <- make_sample_qc(counts_info$counts, metadata)
    summary <- summarize_dataset(counts_info$counts, filtered_counts, metadata)

    shiny::updateSelectizeInput(session, "gene_choice", choices = rownames(filtered_counts), selected = rownames(filtered_counts)[1L], server = TRUE)

    list(
      counts = counts_info$counts,
      filtered_counts = filtered_counts,
      metadata = metadata,
      log_cpm = log_cpm,
      pca = pca,
      distances = distances,
      sample_qc = sample_qc,
      summary = summary
    )
  }, ignoreNULL = FALSE)

  selected_group <- shiny::reactive({
    current <- dataset()
    if (is.null(current)) {
      return(NULL)
    }
    group_column <- input$group_column
    if (is.null(group_column) || !group_column %in% names(current$metadata)) {
      rep("Samples", nrow(current$metadata))
    } else {
      current$metadata[[group_column]]
    }
  })

  output$status <- shiny::renderText({
    current <- dataset()
    paste(
      "Loaded",
      ncol(current$counts), "samples and",
      nrow(current$counts), "genes.",
      "Filtered set retains", nrow(current$filtered_counts), "genes."
    )
  })

  output$dataset_summary <- shiny::renderUI({
    current <- dataset()
    group_lines <- NULL
    if (!is.null(current$summary$groups)) {
      group_lines <- lapply(
        seq_along(current$summary$groups),
        function(i) shiny::tags$li(sprintf("%s: %s samples", names(current$summary$groups)[i], current$summary$groups[[i]]))
      )
    }

    shiny::tagList(
      shiny::tags$h4("Dataset Summary"),
      shiny::tags$p(sprintf("Raw genes: %s", format(current$summary$n_genes, big.mark = ","))),
      shiny::tags$p(sprintf("Genes after abundance filter: %s", format(current$summary$n_filtered_genes, big.mark = ","))),
      shiny::tags$p(sprintf("Samples: %s", current$summary$n_samples)),
      if (!is.null(group_lines)) shiny::tags$ul(group_lines)
    )
  })

  output$sample_table <- DT::renderDataTable({
    DT::datatable(dataset()$sample_qc, options = list(pageLength = 8, scrollX = TRUE))
  })

  output$library_plot <- plotly::renderPlotly({
    current <- dataset()
    group_values <- selected_group()
    plot_data <- current$sample_qc
    plot_data$group <- group_values

    plotly::plot_ly(
      plot_data,
      x = ~run_id,
      y = ~library_size,
      color = ~group,
      type = "bar",
      text = ~paste("Sample:", run_id, "<br>Group:", group, "<br>Library size:", format(library_size, big.mark = ",")),
      hoverinfo = "text"
    ) |>
      plotly::layout(title = "Library Sizes", xaxis = list(title = ""), yaxis = list(title = "Total counts"))
  })

  output$detected_plot <- plotly::renderPlotly({
    current <- dataset()
    plot_data <- current$sample_qc
    plot_data$group <- selected_group()

    plotly::plot_ly(
      plot_data,
      x = ~run_id,
      y = ~detected_genes,
      color = ~group,
      type = "bar",
      text = ~paste("Sample:", run_id, "<br>Group:", group, "<br>Detected genes:", detected_genes),
      hoverinfo = "text"
    ) |>
      plotly::layout(title = "Detected Genes per Sample", xaxis = list(title = ""), yaxis = list(title = "Genes with count > 0"))
  })

  output$pca_plot <- plotly::renderPlotly({
    current <- dataset()
    scores <- as.data.frame(current$pca$scores[, 1:2, drop = FALSE])
    scores$run_id <- rownames(scores)
    scores$group <- selected_group()
    scores$library_size <- compute_lib_sizes(current$counts)
    variance <- current$pca$variance

    plotly::plot_ly(
      scores,
      x = ~PC1,
      y = ~PC2,
      color = ~group,
      type = "scatter",
      mode = "markers+text",
      text = ~run_id,
      textposition = "top center",
      hovertext = ~paste("Sample:", run_id, "<br>Group:", group, "<br>Library size:", format(library_size, big.mark = ",")),
      hoverinfo = "text"
    ) |>
      plotly::layout(
        title = "Principal Component Analysis",
        xaxis = list(title = sprintf("PC1 (%.1f%%)", variance[1L] * 100)),
        yaxis = list(title = sprintf("PC2 (%.1f%%)", variance[2L] * 100))
      )
  })

  output$distance_heatmap <- plotly::renderPlotly({
    current <- dataset()
    sample_ids <- colnames(current$distances)
    groups <- selected_group()
    names(groups) <- current$metadata$run_id
    ordered_groups <- groups[sample_ids]
    group_colors <- make_group_colors(ordered_groups)
    colored_ticktext <- make_colored_ticktext(sample_ids, ordered_groups, group_colors)

    plotly::plot_ly(
      x = sample_ids,
      y = rownames(current$distances),
      z = current$distances,
      type = "heatmap",
      colorscale = "Blues"
    ) |>
      plotly::layout(
        title = "Sample Distance Heatmap",
        xaxis = list(
          title = "",
          tickmode = "array",
          tickvals = sample_ids,
          ticktext = colored_ticktext
        ),
        yaxis = list(
          title = "",
          tickmode = "array",
          tickvals = rownames(current$distances),
          ticktext = colored_ticktext
        )
      )
  })

  output$distance_legend <- shiny::renderUI({
    current <- dataset()
    groups <- selected_group()
    shiny::req(groups)
    group_colors <- make_group_colors(groups)

    legend_items <- lapply(
      names(group_colors),
      function(group_name) {
        shiny::tags$div(
          class = "legend-item",
          shiny::tags$span(
            class = "legend-swatch",
            style = paste("background:", group_colors[[group_name]], ";")
          ),
          shiny::tags$span(group_name)
        )
      }
    )

    shiny::tagList(legend_items)
  })

  output$box_plot <- plotly::renderPlotly({
    current <- dataset()
    log_cpm_long <- stack(as.data.frame(current$log_cpm))
    names(log_cpm_long) <- c("log_cpm", "run_id")
    log_cpm_long$group <- selected_group()[match(log_cpm_long$run_id, dataset()$metadata$run_id)]

    plotly::plot_ly(
      log_cpm_long,
      x = ~run_id,
      y = ~log_cpm,
      color = ~group,
      type = "box",
      boxpoints = FALSE
    ) |>
      plotly::layout(title = "Log2 CPM Distribution by Sample", xaxis = list(title = ""), yaxis = list(title = "log2(CPM + 1)"))
  })

  output$mean_variance_plot <- plotly::renderPlotly({
    current <- dataset()
    gene_means <- rowMeans(current$log_cpm)
    gene_vars <- apply(current$log_cpm, 1L, stats::var)
    mv <- data.frame(
      gene_id = rownames(current$log_cpm),
      mean_log_cpm = gene_means,
      variance_log_cpm = gene_vars,
      stringsAsFactors = FALSE
    )

    plotly::plot_ly(
      mv,
      x = ~mean_log_cpm,
      y = ~variance_log_cpm,
      type = "scatter",
      mode = "markers",
      text = ~paste("Gene:", gene_id, "<br>Mean logCPM:", round(mean_log_cpm, 3), "<br>Variance:", round(variance_log_cpm, 3)),
      hoverinfo = "text"
    ) |>
      plotly::layout(title = "Mean-Variance Trend", xaxis = list(title = "Mean log2(CPM + 1)"), yaxis = list(title = "Variance"))
  })

  output$gene_plot <- plotly::renderPlotly({
    current <- dataset()
    shiny::req(input$gene_choice)
    gene_id <- input$gene_choice
    if (!gene_id %in% rownames(current$counts)) {
      return(NULL)
    }

    plot_data <- data.frame(
      run_id = colnames(current$counts),
      count = as.numeric(current$counts[gene_id, ]),
      log_cpm = as.numeric(log2((current$counts[gene_id, ] / compute_lib_sizes(current$counts)) * 1e6 + 1)),
      group = selected_group(),
      stringsAsFactors = FALSE
    )

    plotly::plot_ly(
      plot_data,
      x = ~group,
      y = ~count,
      color = ~group,
      type = "box",
      boxpoints = "all",
      jitter = 0.2,
      pointpos = 0,
      text = ~paste("Sample:", run_id, "<br>Count:", count, "<br>logCPM:", round(log_cpm, 3)),
      hoverinfo = "text"
    ) |>
      plotly::layout(title = paste("Gene Explorer:", gene_id), xaxis = list(title = "Group"), yaxis = list(title = "Raw counts"))
  })

  output$prep_notes <- shiny::renderUI({
    current <- dataset()
    group_column <- input$group_column
    groups <- if (!is.null(group_column) && group_column %in% names(current$metadata)) unique(current$metadata[[group_column]]) else "Samples"
    enough_groups <- length(groups) >= 2L

    shiny::tagList(
      shiny::tags$h4("Ready-for-DE Checklist"),
      shiny::tags$ul(
        shiny::tags$li(sprintf("Sample annotation column in use: %s", ifelse(is.null(group_column), "run_id", group_column))),
        shiny::tags$li(sprintf("Abundance filter retained %s of %s genes.", format(nrow(current$filtered_counts), big.mark = ","), format(nrow(current$counts), big.mark = ","))),
        shiny::tags$li("Review PCA and distance heatmap for outliers or mislabeled samples."),
        shiny::tags$li("Confirm library sizes and detected genes are balanced enough for the planned contrast."),
        shiny::tags$li("Use the selected annotation column as the design factor for DESeq2, edgeR, or limma-voom downstream.")
      ),
      shiny::tags$p(
        if (enough_groups) {
          sprintf("Current annotation provides %s distinct groups, which is suitable for setting up a differential expression design.", length(groups))
        } else {
          "Only one group is currently visible in the selected annotation column, so add or choose a grouping column before differential expression analysis."
        }
      )
    )
  })
}

message("Starting RNA-seq EDA launcher on http://127.0.0.1:3838")
message("If you are in WSL, open that URL in your browser if it does not launch automatically.")

shiny::runApp(
  shiny::shinyApp(ui = ui, server = server),
  host = "127.0.0.1",
  port = 3838,
  launch.browser = TRUE
)
