#' @description Limma analysis with/without blocking, using
#' duplicateCorrelation, QC, heatmap, and optional GSEA.
#'
#' Runs limma on an abundance matrix. If blocks (repeated measures) are
#' detected, then function accounts for this via duplicateCorrelation().
#'
#' Output files:
#'  - limma_<numerator>_vs_<denominator>.csv
#'  - summary_limma.csv
#'  - heatmap_limma_<numerator>-vs-<denominator>.png
#'  - GSEA_<group>-enriched.csv (only if GSEA runs)
#'
#' @param mtx; numeric matrix; rows are features, columns are samples.
#'
#' @param meta; dataframe; must contain column "variable", which must contain
#' elements passed to numerator and denominator arguments. Rownames of meta
#' must match (not order specific) columns of mtx. If column "block" exists in
#' meta, then patient blocking is performed; otherwise limma without adjusments
#' is performed.
#'
#' @param numerator; character; name of the numerator condition. Must be in
#' meta$variable.
#'
#' @param denominator; character; name of the denominator condition. Must be in
#' meta$variable.
#'
#' @param gsea_gmt_files; character; path to a directory with GMT files, or
#' NULL to skip GSEA.
#'
#' @param output_dir character; directory where output files will be written.
#'
#' @returns list; list has two elements with the following names:
#'  - "limma.results", limma results dataframe with abundance matrix attached
#'  - "gsea", GSEA results.
#'
#'
#' Example run:
#'
#'   set.seed(42)
#'
#'   mtx <- matrix(
#'     data = round(rnorm(20 * 6, mean = 8, sd = 2), 2),
#'     nrow = 20,
#'     ncol = 6,
#'     dimnames = list(
#'       paste("Protein", 1:100),
#'       paste("Sample", 1:6)
#'     )
#'   )
#'
#'   meta <- data.frame(
#'     variable = c("A", "A", "B", "B", "A", "B"),
#'     block = c(
#'       "Patient1", "Patient1", "Patient2",
#'       "Patient2", "Patient3", "Patient3"
#'     ),
#'     row.names = paste("Sample", 1:6)
#'   )
#'
#'  limma_analysis(
#'    mtx = mtx,
#'    meta = meta,
#'    numerator = "A",
#'    denominator = "B",
#'    gsea_gmt_files = NULL,
#'    output_dir = "test_limma"
#'  )
#'
#' @import limma
#' @import fgsea
#' @import dplyr
#' @import magrittr
#' @importFrom stats setNames
#' @importFrom utils write.csv
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_count str_split_1
#' @importFrom glue glue
#' @importFrom paletteer paletteer_c
#' @importFrom pheatmap pheatmap
#' @export
limma_analysis <- function(
    mtx = NULL,
    meta = NULL,
    numerator = NULL,
    denominator = NULL,
    gsea_gmt_files = NULL,
    output_dir = NULL) {
  # Setup ---------------------------------------------------------------------
  # Need these libraries
  library(dplyr)
  library(magrittr)
  library(limma)
  # Short glue
  g <- function(x) glue::glue(x)

  # List of results, we append to this as function progresses.
  results <- list()

  # Make sure output directory exists
  dir.create(
    output_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )

  if (is.null(mtx)) stop("mtx is NULL.")
  if (is.null(meta)) stop("meta is NULL.")
  if (is.null(numerator) || is.null(denominator)) {
    stop("numerator and denominator must be provided.")
  }
  if (is.null(output_dir)) stop("output_dir is NULL.")

  # Ensure mtx is a numeric matrix
  if (is.data.frame(mtx)) mtx <- as.matrix(mtx)
  if (!is.matrix(mtx) || !is.numeric(mtx)) {
    stop("mtx must be a numeric matrix (rows = features, cols = samples).")
  }

  # Ensure meta has required columns and values
  if (!"variable" %in% names(meta)) {
    stop("meta must contain a column named 'variable'.")
  }
  if ("block" %in% names(meta)) {
    print("Paired/Repeated measures detected (blocking).")
    blocking <- TRUE
  } else {
    print("No paired/repeated measures detected (no blocking).")
    blocking <- FALSE
  }
  if (!all(c(numerator, denominator) %in% meta$variable)) {
    missing_levels <- setdiff(c(numerator, denominator), unique(meta$variable))
    stop(paste0(
      "meta$variable must contain both numerator and denominator. Missing: ",
      paste(missing_levels, collapse = ", ")
    ))
  }

  # Exact element match between rownames(meta) and colnames(mtx)
  if (is.null(rownames(meta))) stop("meta must have row names.")
  if (!setequal(rownames(meta), colnames(mtx))) {
    stop("Row names of meta and column names of mtx must match.")
  }

  # Define some functions, because I want to keep the variables used in these
  # functions local, and not clog up the global variable space.
  #
  #' @description sanity check your limma results.
  #' @param l data.frame, limma results from makeLimma() function.
  #' @param which.p string, which p value to use, 'P.Value' or 'adj.P.Val'.
  #' @param p.cut numeric, p value cutoff.
  #' @returns tibble with three columns comparing avg expression values for
  #' top 5 high/low significant logFC values.
  is_my_limma_good <- function(l = NULL,
                               which.p = "P.Value",
                               p.cut = 0.05) {
    log.fc.column <- grep(x = colnames(l), pattern = "logFC", value = TRUE)
    l <- l %>%
      filter(get(which.p) < p.cut) %>%
      arrange(log.fc.column %>% get() %>% desc())
    groups <- grep(x = colnames(l), pattern = "logFC", value = TRUE) %>%
      sub("logFC ", "", .) %>%
      gsub("[\\(,\\)]", "", .) %>%
      str_split_1(pattern = " VS ")
    head.table <- tibble(
      logfc = l %>%
        select(starts_with("logFC")) %>%
        head(5) %>%
        unlist(),
      top = lapply(1:5, function(i) {
        mean_up <- l[
          i,
          grepl(
            x = colnames(l),
            pattern = groups[1],
            fixed = TRUE
          )
        ] %>%
          unlist() %>%
          mean()
        return(mean_up)
      }) %>% unlist(),
      bottom = lapply(1:5, function(i) {
        mean_down <- l[
          i,
          grepl(
            x = colnames(l),
            pattern = groups[2],
            fixed = TRUE
          )
        ] %>%
          unlist() %>%
          mean()
        return(mean_down)
      }) %>% unlist()
    )
    tail.table <- tibble(
      logfc = l %>%
        select(starts_with("logFC")) %>%
        tail(5) %>%
        unlist(),
      top = lapply(
        (nrow(l) - 4):nrow(l), function(i) {
          mean_up <- l[
            i,
            grepl(
              x = colnames(l),
              pattern = groups[1],
              fixed = TRUE
            )
          ] %>%
            unlist() %>%
            mean()
          return(mean_up)
        }
      ) %>% unlist(),
      bottom = lapply(
        (nrow(l) - 4):nrow(l), function(i) {
          mean_down <- l[
            i,
            grepl(
              x = colnames(l),
              pattern = groups[2],
              fixed = TRUE
            )
          ] %>%
            unlist() %>%
            mean()
          return(mean_down)
        }
      ) %>% unlist()
    )
    full.table <- rbind(
      head.table,
      tail.table
    ) %>%
      `colnames<-`(., c(
        "logFC",
        paste0(groups[1], ".sample.mean"),
        paste0(groups[2], ".sample.mean")
      )) %>%
      as_tibble()
    return(full.table)
  }

  #' @description run GSEA analysis on gene_list against GO_file database
  #' @param gene_list numeric vector, values are rank and names are gene
  #' symbols
  #' @param GO_file list; list of paths to GO file(s)
  #' @param significant_only boolean; whether or not to keep only significant
  #' results; default TRUE
  #' @param pcut numeric; pvalue cutoff; default 0.05
  #' @param weight; numeric; how much weight to give the ith gene in the ranked
  #' list. If weight = 0 then all genes are treated equally (use only when you
  #' care about the membership/enrichment, not the strength of the ranking
  #' metric)
  #' @return results from GSEA as dataframe
  GSEA <- function(gene_list = NULL,
                   GO_files = list(),
                   significant_only = TRUE,
                   weight = 1,
                   score = "pos",
                   pcut = 0.05) {
    GO <- lapply(
      GO_files, function(G) {
        fgsea::gmtPathways(G) %>%
          return(.)
      }
    ) %>%
      unlist(recursive = FALSE)
    # GSEA results dataframe
    gsea <- fgsea::fgsea(
      pathways = GO,
      stats = gene_list,
      minSize = 4, # minimum gene set size
      maxSize = 400, # maximum gene set size
      scoreType = score,
      gseaParam = weight
    ) %>%
      as.data.frame() %>%
      mutate(leadingEdge = gsub("\\|", ",", leadingEdge)) %>%
      arrange(pval) %>%
      rename("p-value" = pval) %>%
      rename("Enrichment score" = ES)
    # to get full GSEA path size, build a dataframe of GO pathways and their
    # corresponding size, then merge with GSEA results dataframe above by
    # pathway
    gsea <- data.frame(
      pathway = names(GO),
      pathway_size = sapply(GO, length)
    ) %>%
      left_join(
        gsea, .,
        by = "pathway"
      ) %>%
      mutate(
        pathway = pathway %>%
          tolower() %>%
          gsub("_", " ", .) %>%
          sub("hallmark\\.", "", .) %>%
          sub("canonical\\.", "", .) %>%
          sub("hallmark ", "", .) %>%
          sub("reactome ", "", .) %>%
          sub("kegg ", "", .) %>%
          sub("wp ", "", .),
        k = stringr::str_count(leadingEdge, ",") + 1,
        K = pathway_size,
        `Percent enriched` = round(100 * (k / K), 1),
        `Proteins in pathway` = leadingEdge %>%
          sub("c\\(", "", .) %>%
          sub("\\)", "", .) %>%
          gsub('"', "", .)
      ) %>%
      select(-size, -pathway_size) %>%
      select(
        "pathway", "p-value", "padj", "log2err", "Enrichment score", "NES",
        "k", "K", "Percent enriched", "Proteins in pathway"
      ) %>%
      rename("Pathway" = "pathway") %>%
      arrange(desc(`Enrichment score`))
    if (significant_only) {
      gsea %>%
        filter(`p-value` <= pcut) %>%
        return(.)
    } else {
      return(gsea)
    }
  }

  # limma ---------------------------------------------------------------------
  # Ensures limma directionality.
  meta$variable <- factor(
    meta$variable,
    levels = c(numerator, denominator)
  )

  # Align metadata and matrix.
  meta <- meta[colnames(mtx), , drop = FALSE]
  stopifnot(colnames(mtx) == row.names(meta))

  # Create model matrix (with treatment as the only predictor variable).
  design <- model.matrix(~ 0 + variable, data = meta)
  colnames(design) <- c(numerator, denominator)

  if (isTRUE(blocking)) {
    print("Running limma and correcting for repeated measures.")

    # Model the correlation between samples within the same block
    # (repeated measures).
    corfit <- limma::duplicateCorrelation(mtx, design, block = meta$block)

    # Run limma, correcting for correlation between samples of the same block.
    fit <- lmFit(
      mtx,
      design,
      block = meta$block,
      correlation = corfit$consensus
    )
  } else {
    print("Running limma with no corrections for repeated measures.")
    # Run limma, no corrections.
    fit <- lmFit(mtx, design)
  }

  # Finish limma
  limma.results <- fit %>%
    contrasts.fit(
      contrasts = makeContrasts(
        contrasts = g("{numerator}-{denominator}"),
        levels = colnames(design)
      )
    ) %>%
    eBayes(
      trend = TRUE,
      robust = TRUE
    ) %>%
    topTable(number = Inf)

  # limma post processing -----------------------------------------------------
  print("Cleaning limma")

  # Crystal clear logFC column name
  colnames(limma.results) <- colnames(limma.results) %>%
    sub("logFC", g("logFC {numerator} VS {denominator}"), .)
  stopifnot(row.names(meta) == colnames(mtx))

  # Add contrast group to each sample column in abundance matrix
  colnames(mtx) <- paste(
    row.names(meta),
    meta[["variable"]],
    sep = "; "
  )

  # Make sure protein order is the same between abundance matrix and limma
  # results so we can cbind the abundance matrix to the limma results
  limma.results <- limma.results[limma.results %>%
    row.names() %>%
    order(), ]
  mtx <- mtx[mtx %>%
    row.names() %>%
    order(), ]
  stopifnot(row.names(limma.results) == row.names(mtx))

  # Attach numeric matrix to limma results
  limma.results <- cbind(limma.results, mtx) %>%
    arrange("P.Value") %>%
    mutate(Gene = sub(".*;", "", row.names(.)))

  # Which column is the logFC column
  logfc_col <- grep("^logFC", names(limma.results), value = TRUE)

  # Make limma look pretty
  limma.results <- limma.results %>%
    rownames_to_column("Accession") %>%
    mutate(
      Gene = sub(".*;", "", Accession),
      Accession = sub(";.*", "", Accession),
      rank = -log10(adj.P.Val) * .data[[logfc_col]]
    ) %>%
    arrange(desc(rank)) %>%
    select(Gene, Accession, rank, everything())

  # QC check on limma
  print("Is your limma good?")
  is_my_limma_good(
    l = limma.results,
    which.p = "P.Value",
    p.cut = 0.05
  ) %>% print()
  g("In the top half, everything should be higher in column 2, and column 1
    should be positive") %>%
    print()
  g("In the bottom half, everything should be higher in column 3, and column 1
    should be negative") %>%
    print()

  # Write limma flat file
  limma.results %>%
    rename("Rank (P.Value * logFC)" = "rank") %>%
    write.csv(
      file = file.path(
        output_dir, g("limma_{numerator}_vs_{denominator}.csv")
      ),
      row.names = FALSE
    )

  results[["limma.results"]] <- limma.results

  # Write limma summary -------------------------------------------------------
  horse <- tibble(
    g("{numerator} vs {denominator}"),
    limma.results %>% filter(adj.P.Val < 0.01) %>% nrow(),
    limma.results %>% filter(adj.P.Val < 0.05) %>% nrow(),
    limma.results %>% filter(adj.P.Val < 0.1) %>% nrow(),
    limma.results %>% filter(P.Value < 0.01) %>% nrow(),
    limma.results %>% filter(P.Value < 0.05) %>% nrow(),
    limma.results %>% filter(P.Value < 0.1) %>% nrow()
  )
  colnames(horse) <- c(
    "Comparison",
    "Features passing FDR < 0.01",
    "Features passing FDR < 0.05",
    "Features passing FDR < 0.1",
    "Features passing p-value < 0.01",
    "Features passing p-value < 0.05",
    "Features passing p-value < 0.1"
  )

  summary_path <- file.path(output_dir, "summary_limma.csv")

  # if file doesn't exist, write with header
  if (!file.exists(summary_path)) {
    write.csv(horse, summary_path, row.names = FALSE)
  } else {
    # append only the data, no header
    write.table(
      horse,
      summary_path,
      row.names = FALSE,
      col.names = FALSE,
      sep = ",",
      append = TRUE
    )
  }

  # Plot limma results --------------------------------------------------------
  print("Unsupervised limma analysis.")
  # Prepare matrix
  hm_mtx <- limma.results %>%
    filter(P.Value < 0.05) %>%
    mutate(
      "Accession;Gene" = paste(Accession, Gene, sep = ";")
    ) %>%
    column_to_rownames("Accession;Gene")
  hm_mtx <- hm_mtx[, grep(";", colnames(hm_mtx))]
  colnames(hm_mtx) <- colnames(hm_mtx) %>% sub(";.*", "", .)

  if (nrow(hm_mtx) > 2) { # Because you need >2 to cluster
    # Create colors for comparison
    annotation.palette <- list(
      variable = c("grey30", "grey60")
    )
    names(annotation.palette[["variable"]]) <- c(numerator, denominator)
    class(annotation.palette)

    # Create heatmap colors
    heatmap_palette <- paletteer::paletteer_c(
      palette = "grDevices::Blue-Red 2",
      n = 30
    ) %>%
      as.vector()

    # Create heatmap and save
    png(
      file = file.path(
        output_dir,
        g("heatmap_limma_{numerator}-vs-{denominator}.png")
      ),
      width = 4.5, height = 3.5, res = 300, unit = "in"
    )
    pheatmap::pheatmap(
      mat = hm_mtx,
      main = g("{numerator} vs {denominator}
        {nrow(hm_mtx)} proteins, p-value < 0.01"),
      clustering_method = "ward.D",
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      cellwidth = 2,
      cellheight = 150 / nrow(hm_mtx),
      border_color = "white",
      color = heatmap_palette,
      treeheight_row = 0,
      treeheight_col = 20,
      show_rownames = FALSE,
      show_colnames = FALSE,
      cutree_cols = 2,
      angle_col = 90,
      legend = TRUE,
      annotation_col = meta %>% select(variable),
      annotation_names_col = FALSE,
      annotation_colors = annotation.palette,
      breaks = seq(
        # center scale at 0
        -(hm_mtx %>% abs() %>% max(na.rm = TRUE)),
        hm_mtx %>% abs() %>% max(na.rm = TRUE),
        length.out = 30
      )
    )
    dev.off()
  }


  # GSEA Analysis -------------------------------------------------------------
  print("Performing GSEA.")
  if (!is.null(gsea_gmt_files)) {
    # This list is passed to GSEA search function. It contains the gene sets
    # that are searched against.
    gsea_files <- list(

      # You can add as many gene sets as you want in here. You can call them
      # whatever you like.
      hallmark = file.path(
        gsea_gmt_files,
        "h.all.v2023.2.Hs.symbols.gmt"
      ),
      canonical = file.path(
        gsea_gmt_files,
        "c2.cp.v2023.2.Hs.symbols.gmt"
      )
    )

    # Make sure GSEA files exist.
    lapply(seq_along(gsea_files), function(i) {
      if (!file.exists(gsea_files[[i]])) {
        g("ERROR: {names(gsea_files)[i]} file does not exist
          {gsea_files[[i]]}") %>%
          stop()
      }
      return(NULL)
    })

    # Create a list of genes and order it by logFC * -log10(pval).
    genes <- limma.results %>%
      filter(adj.P.Val < 0.05) %>%
      distinct(`Gene`, .keep_all = TRUE) %>%
      `rownames<-`(., .[["Gene"]]) %>%
      arrange(desc(`rank`)) %>%
      select(`rank`)
    ranked_genes <- unlist(genes)
    names(ranked_genes) <- sub(
      ";.*$", "", row.names(genes)
    )

    # GSEA analysis on the biggest positive logFC proteins (pos) and biggest
    # negative logFC proteins (neg).
    gsea_param <- c("pos", "neg")
    names(gsea_param) <- c(numerator, denominator)

    gsea_combined <- lapply(seq_along(gsea_param), function(i) {
      enriched <- GSEA(
        gene_list = ranked_genes,
        score = gsea_param[i],
        GO_files = gsea_files
      )

      # add labels for downstream clarity
      if (!is.null(enriched) && nrow(enriched) > 0) {
        enriched <- enriched %>%
          mutate(
            enriched.in = names(gsea_param)[i],
            contrast = g("{numerator} vs {denominator}")
          ) %>%
          select(
            contrast, enriched.in, everything()
          )
      }
      return(enriched)
    }) %>%
      bind_rows()

    results[["gsea"]] <- gsea_combined
    write.csv(
      x = gsea_combined,
      file = file.path(
        output_dir,
        g("GSEA_{numerator}-vs-{denominator}.csv")
      ),
      row.names = FALSE
    )
  }

  print("Done.")
  return(results)
}
