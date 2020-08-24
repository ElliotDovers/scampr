.onLoad <- function(lib, pkg) {
  cat("Loading compiled code...\n")
  library.dynam("scampr", pkg, lib)
}
