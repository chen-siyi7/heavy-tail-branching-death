main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (identical(args, "--help")) {
    cat(
      "Usage: Rscript --vanilla run.R [--output-dir PATH]\n",
      "Default output parent: results/ beside this script.\n",
      "An explicit PATH is relative to the current working directory.\n",
      "The audit creates simulation_audit_output/ inside that parent.\n",
      "Existing audit directories are not overwritten.\n",
      sep = ""
    )
    return(invisible(NULL))
  }
  if (!(length(args) == 0L ||
        (length(args) == 2L && args[1] == "--output-dir" &&
         nzchar(args[2]) && !startsWith(args[2], "--")))) {
    stop("Invalid arguments. Run Rscript --vanilla run.R --help.",
         call. = FALSE)
  }

  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_arg) != 1L) {
    stop("Please run this entry point with Rscript.", call. = FALSE)
  }
  script_path <- normalizePath(sub("^--file=", "", script_arg),
                               winslash = "/", mustWork = TRUE)
  project_dir <- dirname(script_path)
  source_path <- file.path(project_dir, "simulation_audit.R")
  if (!file.exists(source_path)) {
    stop("simulation_audit.R must be beside run.R.", call. = FALSE)
  }

  output_parent <- if (length(args)) path.expand(args[2]) else {
    file.path(project_dir, "results")
  }
  if (!dir.exists(output_parent) &&
      !dir.create(output_parent, recursive = TRUE)) {
    stop("Could not create the output directory.", call. = FALSE)
  }
  output_parent <- normalizePath(output_parent, winslash = "/",
                                 mustWork = TRUE)
  output_dir <- file.path(output_parent, "simulation_audit_output")
  if (file.exists(output_dir)) {
    stop(paste0("Output already exists: ", output_dir,
                ". Choose a fresh --output-dir."), call. = FALSE)
  }

  original_dir <- getwd()
  on.exit(setwd(original_dir), add = TRUE)
  setwd(output_parent)

  # Explicitly retain the generator settings used in the manuscript audit.
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  started <- proc.time()[["elapsed"]]
  audit <- new.env(parent = globalenv())
  message("Running the complete audit. Results: ", output_dir)
  sys.source(source_path, envir = audit, chdir = FALSE)

  expected_counts <- c(csv = 15L, pdf = 6L, tex = 8L)
  actual_counts <- vapply(names(expected_counts), function(ext) {
    length(list.files(output_dir, pattern = paste0("\\.", ext, "$")))
  }, integer(1))
  stopifnot(identical(actual_counts, expected_counts))
  generated <- list.files(output_dir, full.names = TRUE)
  stopifnot(all(file.info(generated)$size > 0))

  writeLines(capture.output(sessionInfo()),
             file.path(output_dir, "session_info.txt"))
  writeLines(c(
    "Reproducibility run information",
    paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Master seed:", audit$master_seed),
    paste("RNG kinds:", paste(RNGkind(), collapse = ", ")),
    paste("Elapsed seconds:", round(proc.time()[["elapsed"]] - started, 2)),
    paste("simulation_audit.R MD5:", unname(tools::md5sum(source_path))),
    "No --publish option was passed; manuscript files were not copied."
  ), file.path(output_dir, "run_info.txt"))

  generated <- sort(list.files(output_dir, full.names = TRUE))
  writeLines(paste(unname(tools::md5sum(generated)), basename(generated)),
             file.path(output_dir, "checksums.md5"))
  message("Complete. Inspect audit_summary.txt and session_info.txt in ",
          output_dir)
}

main()
