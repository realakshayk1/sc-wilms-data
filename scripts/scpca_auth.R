#!/usr/bin/env Rscript
# One-time ScPCA auth helper. Opens terms, prompts for email, stores SCPCA_AUTH_TOKEN.
args <- commandArgs(trailingOnly = TRUE)
email <- if (length(args)) args[1] else ""

main <- function() {
  if (!requireNamespace("ScPCAr", quietly = TRUE)) {
    stop("Install ScPCAr first: scripts\\rscript.bat scripts\\install_r_packages.R")
  }
  message("Opening ScPCA terms of use in your browser...")
  ScPCAr::view_terms()
  if (!nzchar(email)) email <<- readline("Your academic email: ")
  if (!nzchar(trimws(email))) stop("Email required.")
  token <- ScPCAr::get_auth(email = trimws(email), agree = TRUE)
  message("[ok] Token stored in SCPCA_AUTH_TOKEN (length ", nchar(token), ")")
  invisible(token)
}

if (sys.nframe() == 0) main()
