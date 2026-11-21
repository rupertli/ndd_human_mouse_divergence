required <- c("DBI", "RSQLite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
docs <- dirname(dirname(script_path))
workspace <- normalizePath(file.path(docs, "..", "..", ".."), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
db_path <- if (length(args) >= 1) args[[1]] else file.path(workspace, "brainspan-dev", "brainspan-dev", "extdata", "brainspan.db")
gene_path <- if (length(args) >= 2) args[[2]] else file.path(workspace, "scRNASeq", "NDD_human_mouse_trajectory_cerebellum_6_Phases", "ndd_gene_provenance.csv")
if (!file.exists(db_path) || !file.exists(gene_path)) stop("Database or NDD provenance file is missing; pass both paths as arguments.")
out <- file.path(docs, "data", "BrainSpan_NDD")
dir.create(file.path(out, "genes"), recursive = TRUE, showWarnings = FALSE)
write_csv <- function(data, path) write.csv(data, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")

export <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path, flags = RSQLite::SQLITE_RO)
  on.exit(DBI::dbDisconnect(con))
  provenance <- read.csv(gene_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("gene", "sources") %in% names(provenance))) stop("Expected gene and sources in NDD provenance.")
  if (anyNA(provenance$gene) || anyDuplicated(provenance$gene) || any(!nzchar(provenance$gene))) stop("Invalid NDD gene list.")
  provenance <- provenance[order(provenance$gene), ]
  genes <- provenance$gene
  samples <- DBI::dbReadTable(con, "phenotype_data")
  columns <- c("sample_id", "donor_id", "age", "gender", "Days", "age_group", "structure_acronym", "structure_name", "structure_group1", "structure_group2")
  if (!all(columns %in% names(samples)) || anyDuplicated(samples$sample_id) || any(!is.finite(samples$Days) | samples$Days <= 0)) stop("Invalid phenotype metadata.")
  samples <- samples[, columns]
  sample_ids <- as.character(samples$sample_id)
  for (table in c("expression_data", "activity_data")) {
    if (!setequal(setdiff(DBI::dbListFields(con, table), c("gene_id", "driver_id")), sample_ids)) stop("Sample columns do not match phenotype IDs in ", table)
  }
  drivers <- DBI::dbReadTable(con, "ms_tab")
  drivers <- drivers[drivers$originalID %in% genes, c("originalID", "driver_id", "funcType")]
  if (anyDuplicated(drivers$driver_id)) stop("Duplicate driver mapping.")
  names(drivers) <- c("gene", "driver_id", "function_type")
  drivers <- drivers[order(drivers$gene, drivers$driver_id), ]
  quoted_genes <- paste(DBI::dbQuoteString(con, genes), collapse = ",")
  networks <- DBI::dbGetQuery(con, paste0(
    "SELECT source, target, \"source.symbol\" AS source_gene, \"target.symbol\" AS target_gene, MI, pearson, spearman FROM network_data ",
    "WHERE \"source.symbol\" IN (", quoted_genes, ") AND \"target.symbol\" IN (", quoted_genes, ")"
  ))
  networks <- networks[networks$source %in% drivers$driver_id, ]
  if (any(!is.finite(networks$MI) | networks$MI < 0)) stop("Invalid network MI.")
  write_csv(provenance, file.path(out, "ndd_gene_provenance.csv"))
  write_csv(samples, file.path(out, "samples.csv"))
  write_csv(drivers, file.path(out, "drivers.csv"))
  manifest <- data.frame(gene = genes, sources = provenance$sources, file = "", expression = FALSE, activity_drivers = 0L, network_edges = 0L)
  profile_columns <- c("gene", "profile", "id", sample_ids)
  batch_size <- 100L
  batch_number <- 1L
  profile_batch <- list()
  network_batch <- list()
  write_batch <- function() {
    stem <- sprintf("batch-%05d", batch_number)
    write_csv(do.call(rbind, profile_batch), file.path(out, "genes", paste0(stem, "-profiles.csv")))
    write_csv(do.call(rbind, network_batch), file.path(out, "genes", paste0(stem, "-network.csv")))
  }
  for (i in seq_along(genes)) {
    gene <- genes[[i]]
    expr <- DBI::dbGetQuery(con, "SELECT * FROM expression_data WHERE gene_id = ?", params = list(gene))
    ids <- drivers$driver_id[drivers$gene == gene]
    if (length(ids)) {
      act <- DBI::dbGetQuery(con, paste0("SELECT * FROM activity_data WHERE driver_id IN (", paste(DBI::dbQuoteString(con, ids), collapse = ","), ")"))
    } else {
      act <- data.frame()
    }
    if (nrow(expr) > 1 || (nrow(act) && anyDuplicated(act$driver_id))) stop("Duplicate profile for ", gene)
    profile <- data.frame(matrix(nrow = 0, ncol = length(profile_columns)))
    names(profile) <- profile_columns
    if (nrow(expr)) profile <- rbind(profile, data.frame(gene = gene, profile = "Expression", id = expr$gene_id, expr[, sample_ids, drop = FALSE], check.names = FALSE))
    if (nrow(act)) profile <- rbind(profile, data.frame(gene = gene, profile = "Activity", id = act$driver_id, act[, sample_ids, drop = FALSE], check.names = FALSE))
    values <- as.matrix(profile[, sample_ids, drop = FALSE])
    if (any(!is.finite(values) & !is.na(values))) stop("Non-finite profile values for ", gene)
    edges <- networks[networks$source %in% ids, ]
    edges <- edges[order(edges$source, edges$target), ]
    manifest$expression[[i]] <- nrow(expr) > 0
    manifest$activity_drivers[[i]] <- nrow(act)
    manifest$network_edges[[i]] <- nrow(edges)
    if (nrow(profile) || nrow(edges)) {
      stem <- sprintf("batch-%05d", batch_number)
      manifest$file[[i]] <- stem
      profile_batch[[length(profile_batch) + 1L]] <- profile
      network_batch[[length(network_batch) + 1L]] <- edges
      if (length(profile_batch) == batch_size) {
        write_batch()
        batch_number <- batch_number + 1L
        profile_batch <- list()
        network_batch <- list()
      }
    }
  }
  if (length(profile_batch)) write_batch()
  expected <- unlist(lapply(unique(manifest$file[nzchar(manifest$file)]), function(stem) paste0(stem, c("-profiles.csv", "-network.csv"))))
  existing <- list.files(file.path(out, "genes"), pattern = "^(gene|batch)-[0-9]{5}-(profiles|network)\\.csv$")
  write_csv(manifest, file.path(out, "manifest.csv"))
  stale <- setdiff(existing, expected)
  if (length(stale)) {
    for (name in stale) {
      path <- file.path(out, "genes", name)
      if (unlink(path) != 0) stop("Could not remove stale generated shard: ", path)
    }
  }
  cat("Exported", length(genes), "NDD genes;", sum(manifest$expression), "expression profiles;",
      sum(manifest$activity_drivers), "activity profiles;", nrow(networks), "NDD-only edges;",
      nrow(samples), "samples;", length(expected), "batch CSV files.\n")
}
export()
