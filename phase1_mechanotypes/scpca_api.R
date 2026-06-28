# ScPCAr API wrappers — aligned with package docs as of 2026-06.
# https://alexslemonade.github.io/ScPCAr/reference/

ensure_scpcar <- function() {
  if (!requireNamespace("ScPCAr", quietly = TRUE)) {
    stop(
      "ScPCAr not installed. Run:\n",
      "  remotes::install_github('AlexsLemonade/ScPCAr')"
    )
  }
  invisible(TRUE)
}

ensure_scpca_auth <- function() {
  token <- get_scpca_token()
  if (nzchar(token)) return(invisible(token))
  stop(
    "No ScPCA auth token. In R, run once:\n",
    "  library(ScPCAr)\n",
    "  ScPCAr::view_terms()\n",
    "  ScPCAr::get_auth(email = 'you@upenn.edu', agree = TRUE)\n",
    "get_auth() stores SCPCA_AUTH_TOKEN automatically."
  )
}

has_create_dataset_api <- function() {
  ensure_scpcar()
  exists("create_dataset", where = asNamespace("ScPCAr"), inherits = FALSE)
}

#' Download ScPCA data using current API patterns.
#'
#' - Full project: download_project() (pre-built CCDL zip; fast for exploration)
#' - Single sample: create_dataset() + download_dataset(await_processing=TRUE)
#'   because download_sample() is deprecated (computed-files endpoint removed)
download_scpca_data <- function(
  project_id,
  destination,
  sample_id = "",
  format = c("sce", "anndata", "spatial"),
  merged = TRUE,
  email = NULL,
  timeout_minutes = 120
) {
  ensure_scpcar()
  ensure_scpca_auth()
  format <- match.arg(format)
  if (!dir.exists(destination)) dir.create(destination, recursive = TRUE, showWarnings = FALSE)

  if (nzchar(sample_id)) {
    return(download_scpca_samples(
      sample_ids = sample_id,
      destination = destination,
      format = format,
      email = email,
      timeout_minutes = timeout_minutes
    ))
  }

  message("[scpca] download_project format=", format,
          if (format != "spatial") paste0(" merged=", merged) else "")
  paths <- ScPCAr::download_project(
    project_id = project_id,
    destination = destination,
    format = format,
    merged = merged && format != "spatial"
  )
  invisible(paths)
}

download_scpca_samples <- function(
  sample_ids,
  destination,
  format = c("sce", "anndata", "spatial"),
  email = NULL,
  timeout_minutes = 120
) {
  format <- match.arg(format)

  if (!has_create_dataset_api()) {
    warning(
      "Installed ScPCAr lacks create_dataset(); using deprecated download_sample(). ",
      "Upgrade: remotes::install_github('AlexsLemonade/ScPCAr')"
    )
    paths <- ScPCAr::download_sample(
      sample_id = sample_ids[1],
      destination = destination,
      format = format
    )
    return(invisible(paths))
  }

  # create_dataset() accepts only sce/anndata; spatial libraries return Space Ranger output.
  ds_format <- if (format == "anndata") "anndata" else "sce"
  sample_ids <- unique(sample_ids)

  message("[scpca] create_dataset samples=", paste(sample_ids, collapse = ", "))
  ds_id <- ScPCAr::create_dataset(
    samples = sample_ids,
    format = ds_format,
    email = email
  )
  message("[scpca] dataset id=", ds_id, " — processing (often 1–30+ min)...")

  paths <- ScPCAr::download_dataset(
    ds_id,
    destination = destination,
    await_processing = TRUE,
    poll_interval = 0.5,
    timeout = timeout_minutes
  )
  invisible(paths)
}

filter_wilms_samples <- function(samples) {
  # Nucleus snRNA vs Visium spots for SCPCP000006 exploration.
  out <- list(all = samples)
  if ("seq_units" %in% colnames(samples)) {
    su <- samples$seq_units
    out$nucleus <- samples[grepl("nucleus", su, ignore.case = TRUE), , drop = FALSE]
    out$spot <- samples[grepl("spot", su, ignore.case = TRUE), , drop = FALSE]
  }
  if ("has_spatial_data" %in% colnames(samples)) {
    out$spatial <- samples[samples$has_spatial_data %in% TRUE, , drop = FALSE]
  }
  out
}
