## "sfareg" is an S3 class, and texreg's extract() is an S4 generic. Declaring
## the S3 class to the S4 system is what lets it appear in a method signature;
## without it setMethod() warns "no definition for class sfareg" and the
## signature is matched only by accident.
methods::setOldClass("sfareg")

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "sfa version ", utils::packageVersion("sfa"), "\n",
    "Type citation('sfa') for citing this package in publications."
  )
}

## sandwich is in Suggests, so its generics cannot be imported. Register the
## methods at load time when it is present -- the standard pattern for making
## an optional package's generics work without depending on it.
.onLoad <- function(libname, pkgname) {
  if (requireNamespace("sandwich", quietly = TRUE)) {
    registerS3method("bread", "sfareg", bread.sfareg,
                     envir = asNamespace("sandwich"))
    registerS3method("estfun", "sfareg", estfun.sfareg,
                     envir = asNamespace("sandwich"))
  }
  ## texreg::extract is an S4 generic, so it needs setMethod() and not
  ## registerS3method() or NAMESPACE's S3method(). S3method() writes into the
  ## S3 methods table and S4 dispatch never reads that table, so the line
  ## NAMESPACE carried was INERT: texreg::screenreg(list(fit)) on an "sfareg"
  ## fit reported "Neither texreg nor broom supports models of class sfareg"
  ## for a full release while the registration looked present (gap A79).
  ##
  ## tryCatch so that nothing here can stop the package loading. An optional
  ## package's generic is not worth a load failure, and a texreg whose S4
  ## generic changes signature would otherwise take sfa down with it.
  ## setMethod() is passed the GENERIC ITSELF rather than the name "extract".
  ## Both alternatives fail, and were measured failing on an installed build
  ## before this form was chosen:
  ##   setMethod("extract", ..., where = asNamespace("sfa"))
  ##     -> "the environment 'sfa' is locked; cannot assign methods"
  ##   setMethod("extract", ...)                  # no `where`
  ##     -> "no existing definition for function 'extract'", because texreg is
  ##        loaded but not attached, so the generic is not on the search path
  if (requireNamespace("texreg", quietly = TRUE)) {
    tryCatch(
      methods::setMethod(texreg::extract, signature = "sfareg",
        definition = function(model, ...) extract.sfareg(model, ...)
      ),
      error = function(e) NULL
    )
  }
  invisible()
}
