suppressPackageStartupMessages({
  library(msigdbr)
})

options(timeout=3600)

root <- "/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA"
out_dir <- file.path(root, "02_metadata/msigdb_2026_1")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

go_file <- file.path(out_dir, "MSigDB_C5_GO_BP_HS.rds")
kegg_file <- file.path(out_dir, "MSigDB_C2_KEGG_LEGACY_HS.rds")

download_with_retry <- function(label, expression, output_file, max_attempts=3) {
  if (file.exists(output_file)) {
    x <- readRDS(output_file)

    if (is.data.frame(x) && nrow(x) > 0) {
      cat(
        "[EXISTS] ", label,
        " rows=", nrow(x),
        " file=", output_file, "\n",
        sep=""
      )
      return(x)
    }
  }

  last_error <- NULL

  for (attempt in seq_len(max_attempts)) {
    cat(
      "[DOWNLOAD] ", label,
      " attempt=", attempt,
      "/", max_attempts, "\n",
      sep=""
    )

    result <- tryCatch(
      eval(expression),
      error=function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )

    if (!is.null(result) && nrow(result) > 0) {
      saveRDS(result, output_file)

      cat(
        "[SAVED] ", label,
        " rows=", nrow(result),
        " pathways=", length(unique(result$gs_name)),
        " file=", output_file, "\n",
        sep=""
      )

      return(result)
    }

    cat(
      "[RETRY] ", label,
      " error=", last_error, "\n",
      sep=""
    )
  }

  stop(
    "Failed to download ",
    label,
    " after ",
    max_attempts,
    " attempts. Last error: ",
    last_error
  )
}

go_bp <- download_with_retry(
  label="GO_BP",
  expression=quote(
    msigdbr(
      db_species="HS",
      species="Homo sapiens",
      collection="C5",
      subcollection="GO:BP"
    )
  ),
  output_file=go_file
)

kegg <- download_with_retry(
  label="KEGG_LEGACY",
  expression=quote(
    msigdbr(
      db_species="HS",
      species="Homo sapiens",
      collection="C2",
      subcollection="CP:KEGG_LEGACY"
    )
  ),
  output_file=kegg_file
)

stopifnot(
  nrow(go_bp) > 0,
  nrow(kegg) > 0,
  length(unique(go_bp$gs_name)) > 1000,
  length(unique(kegg$gs_name)) > 100
)

manifest <- data.frame(
  database=c("GO_BP", "KEGG_LEGACY"),
  msigdbr_version=as.character(packageVersion("msigdbr")),
  database_version=c(
    paste(unique(go_bp$db_version), collapse=","),
    paste(unique(kegg$db_version), collapse=",")
  ),
  rows=c(nrow(go_bp), nrow(kegg)),
  pathways=c(
    length(unique(go_bp$gs_name)),
    length(unique(kegg$gs_name))
  ),
  file=c(go_file, kegg_file)
)

write.table(
  manifest,
  file.path(out_dir, "MSigDB_manifest.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

cat("===== MSIGDB MANIFEST =====\n")
print(manifest)
cat("STATUS=PASS\n")
