# Validate that scripts are executed from the repository root.
# Usage: source("scripts/validate_paths.R")

repo_root_markers <- c("README.md", "code", "configs")
missing <- repo_root_markers[!file.exists(repo_root_markers)]
if (length(missing) > 0) {
  stop(
    "Please run scripts from the repository root. Missing: ",
    paste(missing, collapse = ", ")
  )
}
message("OK: repository root detected.")
