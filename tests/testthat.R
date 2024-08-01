
pkg_name <- "bsafe"

library(testthat)
library(pkg_name, character.only = T)

is_CI <- isTRUE(as.logical(Sys.getenv("CI")))

if (is_CI) {
  utils::capture.output(
    test_package(pkg_name, reporter = "junit"),
    file = "./tests/test-out.xml"
  )
} else {
  test_check(pkg_name)
}
