# Resolve and set the project root as the working directory.
# Works when a script is sourced() or run with Rscript; otherwise
# assumes the current working directory is already the project root.

local({
  root <- NULL

  # Rscript --file=...
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
    if (file.exists(script_path)) {
      root <- normalizePath(file.path(dirname(script_path), ".."))
    }
  }

  # source("scripts/....R")
  if (is.null(root)) {
    ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
    if (!is.null(ofile) && nzchar(ofile)) {
      root <- normalizePath(file.path(dirname(ofile), ".."))
    }
  }

  # Fallback: walk up from getwd() looking for marker files
  if (is.null(root)) {
    d <- normalizePath(getwd())
    for (i in seq_len(6)) {
      if (file.exists(file.path(d, "Heliosperma_smRNAs.Rproj")) ||
          file.exists(file.path(d, "scripts", "00_set_project_root.R"))) {
        root <- d
        break
      }
      parent <- dirname(d)
      if (identical(parent, d)) break
      d <- parent
    }
  }

  if (!is.null(root)) {
    setwd(root)
  } else {
    message(
      "Could not auto-detect project root. ",
      "Please setwd() to the Heliosperma_smRNAs repository root before running scripts."
    )
  }
})
