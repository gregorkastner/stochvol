context("C++")
test_that("Catch unit tests pass", {
    skip_on_cran()
    expect_cpp_tests_pass("stochvol")
})
