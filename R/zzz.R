.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to the 'VarReg' package to perform semi-parametric regression")
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  op.varreg <- list(
    varreg.path = "~/R-dev",
    varreg.install.args = "",
    varreg.name = "Kristy Robledo",
    varreg.desc.author = '"Kristy Robledo <robledo.kristy@gmail.com> [aut, cre]"',
    varreg.desc.license = "GPL (>= 2)",
    varreg.desc.suggests = NULL,
    varreg.desc = list()
  )
  toset <- !(names(op.varreg) %in% names(op))
  if(any(toset)) options(op.varreg[toset])

  invisible()
}