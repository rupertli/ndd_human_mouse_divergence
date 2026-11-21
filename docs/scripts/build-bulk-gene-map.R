required <- c("AnnotationDbi", "org.Hs.eg.db", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
data_dir <- file.path(dirname(dirname(script_path)), "data", "ASD_Bulk_analysis")
assignments <- read.csv(file.path(data_dir, "cluster_result", "expression_kmeans_all_clusters_genesv2.csv"))
genes <- sort(unique(trimws(assignments$Gene)))
mapped <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db, keys = genes, keytype = "SYMBOL",
  column = "ENTREZID", multiVals = "list"
)
mapped <- lapply(mapped, function(ids) sort(unique(as.character(stats::na.omit(ids)))))
unmapped <- names(mapped)[lengths(mapped) == 0]
if (length(unmapped)) stop("Unmapped cluster genes: ", paste(unmapped, collapse = ", "))
jsonlite::write_json(
  list(annotation_package = "org.Hs.eg.db",
       annotation_version = as.character(utils::packageVersion("org.Hs.eg.db")), genes = mapped),
  file.path(data_dir, "gene-entrez-map.json"), auto_unbox = FALSE, pretty = TRUE
)
cat("Mapped", length(mapped), "cluster genes without contacting external services.\n")
