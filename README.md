# RNA_Hypothesis_Workbench
This is a compilation of R packages, integrated in a Shiny App-based browser-accessible interface for transcriptomic data analysis, combined with prior knowledge to narrow down potential mechanisms of substance effects on biological systems and based on them, create hypotheses, testable in the lab.

This workspace includes browser-based R launchers for exploratory analysis of the files in `data/`.
The workbench also creates a `cache/` directory to store reusable external/database-backed results.

## Start Here

To avoid clutter for new users, only focus on these files:

- `START_HERE.txt`
- `START_WORKBENCH_Windows.bat`
- `START_WORKBENCH_Linux.desktop`

Windows:

1. Open `START_HERE.txt` if you need setup help.
2. Double-click `START_WORKBENCH_Windows.bat`.

Linux desktop:

1. Open `START_HERE.txt` if you need setup help.
2. Double-click `START_WORKBENCH_Linux.desktop`.

WSL or Linux terminal:

```bash
bash launchers/start_workbench.sh
```

Both launchers ultimately run:

```text
Rscript deseq2_workbench.R
```

The browser should then open automatically. If it does not, visit [http://127.0.0.1:3838](http://127.0.0.1:3838).

## Run from WSL Ubuntu

```bash
cd /mnt/c/Users/USER/Documents/Codex_Try
Rscript rnaseq_eda_launcher.R
```

If the browser does not open automatically in WSL, visit [http://127.0.0.1:3838](http://127.0.0.1:3838).

To launch the DESeq2 workflow instead:

```bash
Rscript deseq2_workbench.R
```

## Prepare SEMgraph validation files

To create a local STRING file for SEMgraph edge validation, first export a gene list from the app or save a simple text file with one HGNC symbol per line, then run:

```bash
Rscript download_semgraph_validation_data.R --genes-file data/my_semgraph_genes.txt --out-dir data
```

This helper creates:

- `data/STRING_edges.tsv` from the live STRING network API for the supplied genes

For a CSV/TSV gene list, you can also specify the column:

```bash
Rscript download_semgraph_validation_data.R --genes-file data/my_semgraph_genes.csv --column gene_label --out-dir data
```

## Required R packages

The EDA launcher checks for:

- `shiny`
- `plotly`
- `DT`

If any are missing, the script prints an install command and exits.

The DESeq2 workbench also requires:

- `DESeq2`
- `clusterProfiler`
- `WGCNA`
- `SEMgraph`
- `AnnotationDbi`
- `org.Hs.eg.db`
- `GO.db`
- `jsonlite`

## Manual verification

- The app loads both files from `data/` without changing the paths.
- The DESeq2 workbench accepts either pasted file paths or files chosen through the browse picker.
- Gene-facing tables and plots prefer HGNC symbols, while still retaining Ensembl IDs where needed for traceability.
- The overview shows 8 samples and the available treatment annotations.
- Library size and detected-gene plots show one bar per sample with hover details.
- PCA separates samples interactively and tooltips show sample annotations.
- The distance heatmap is interactive and labels all samples.
- The gene explorer updates when a different gene ID is selected.
- The DESeq2 workbench lets you save multiple design formulas and rerun different contrasts.
- The DESeq2 results view lets you tune both the adjusted p-value cutoff and the absolute log2 fold-change cutoff used for DEG selection.
- The DESeq2 results area is split into fit summary, transform checks, statistical diagnostics, and signal views.
- VST PCA and VST distribution plots render after a DESeq2 run.
- The workbench includes a mean-SD plot after VST, a dispersion interpretation panel, DESeq2 size factor diagnostics, and a significance-versus-mean plot.
- Raw p-value and adjusted p-value histograms render after a DESeq2 run.
- A separate pathway enrichment stage supports ORA and GSEA, exposes beginner-friendly tuning parameters, and shows a simple consensus across methods.
- In pathway enrichment, selected methods run separately; the app now provides both a per-method inspection view and a separate consensus-overlap view for clarity.
- A separate WGCNA stage supports automatic soft-threshold suggestion, module split/merge tuning, module-trait heatmaps, and dendrogram-guided refinement for beginner users.
- The WGCNA stage now also includes per-module network-style visualization and hub-gene ranking to help interpret selected modules.
- The WGCNA stage also includes a diagnostic-guided `Programs` view that suggests when two inverse modules may represent opposite arms of the same biological process, then lets you inspect a combined interpretation layer without overwriting the original WGCNA modules.
- A SEMgraph causal-analysis stage helps explore likely upstream and downstream regulators for a selected module using graph-based structural equation modeling.
- The SEMgraph causal graph now layers in multi-source pharmacology annotations: gold-ring genes have external drug-target or drug-pathway support, and hover/table views show which source contributed the evidence.
- DGIdb is used as a live drug-gene interaction source with action-style labels when available.
- HCDT is supported as an optional local export placed in `data/` with a name such as `HCDT_drug_gene.tsv` or `HCDT_drug_gene.csv`.
- Reactome contributes drug-mediated pathway context, which is shown separately from direct drug-target evidence.
- For SEMgraph edge validation, STRING can be queried from a local export such as `data/STRING_edges.tsv` or `data/STRING_edges.csv`.
- The main SEMgraph pane now also includes side buttons to:
  - download `data/STRING_edges.tsv` for the current SEMgraph genes
- The in-app STRING download progress now reports downloaded size, total size when available, transfer speed, and ETA.
- If those validation files are missing and live STRING lookup is not enabled, the app now explicitly reports that zero validated edges reflect missing validation sources rather than evidence against the biology.
- The SEMgraph stage now also includes a separate `Pathway Context` view that highlights which displayed causal-graph genes overlap Reactome pathways, how those genes are positioned as focal/upstream/downstream/context nodes, and what direct druggability evidence exists around those pathway-overlapping components.
- SEMgraph can now be run on either one original module or one saved combined program, so users can compare a narrow module-centric causal view against a broader biology-centric combined-program view.
- SEMgraph now separates graph scope from table scope: the graph can be expanded to a larger zoomable node set for overview, while the tables use their own evidence-based cutoffs so they remain more comprehensive without being forced to match the exact displayed node subset.
- SEMgraph now also includes an external edge-validation layer based on STRING functional-association support. The causal graph can be filtered to show only STRING-supported edges.
- Relevant panes now include refresh checkboxes. If left unchecked, the app reuses cached KEGG/enrichment and pharmacology/context lookups when the same inputs were already analyzed, which makes repeat exploration much faster.
- SEMgraph also caches fitted causal models and uses a prioritized gene cap before fitting, so very large combined programs do not hang indefinitely during graph orientation. The UI reports how many genes were in the original unit and how many actually entered the SEMgraph fit.
- SEMgraph also separately caps the number of genes sent to slow external annotation sources such as Reactome and updates the progress bar during those annotation batches, so long pharmacology/pathway-context steps are easier to distinguish from a true stall.
- Live Reactome lookup is optional in SEMgraph. By default the app prefers cached Reactome context or skips live Reactome requests, because that stage is network-bound and often much slower than the rest of the analysis.
- Reactome can now be downloaded once into `cache/reactome/` using the in-app `Update Reactome cache` button. After that, SEMgraph uses the local Reactome mapping cache for much faster pathway-context annotation.
- A `Reports` tab now summarizes package versions and analysis actions, and can generate two export bundles:
  - a reproducibility-oriented R Markdown bundle with snapshot data and code chunks
  - a shareable HTML bundle zip with plots, tables, explanations, and CSV exports
- Tables now support CSV export plus richer filtering via top-of-column filters and SearchBuilder rules such as exact match, greater than, smaller than, or equals.
- Plotly plots expose their modebar so they can be downloaded as images from the app, and selected network views can also be exported as GraphML.
