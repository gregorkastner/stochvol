context("Validation")

test_that("assert_numeric passes good input", {
  expect_true(assert_numeric(4L, "x"))
  expect_true(assert_numeric(c(4, -1), "x"))
  expect_true(assert_numeric(c(4, -Inf), "x"))
})

test_that("assert_numeric failes bad input", {
  expect_error(assert_numeric(list(), "x"))
  expect_error(assert_numeric(c(), "x"))
  expect_error(assert_numeric("a", "x"))
  expect_error(assert_numeric(NULL, "x"))
  expect_error(assert_numeric(c(1, NA), "x"))
  expect_error(assert_numeric(c(1, NA_real_), "x"))
  expect_error(assert_numeric(c(1, NaN), "x"))
})

test_that("assert_logical passes good input", {
  expect_true(assert_logical(TRUE, "x"))
  expect_true(assert_logical(FALSE, "x"))
  expect_true(assert_logical(c(TRUE, FALSE), "x"))
})

test_that("assert_logical failes bad input", {
  expect_error(assert_logical(list(), "x"))
  expect_error(assert_logical(c(), "x"))
  expect_error(assert_logical("a", "x"))
  expect_error(assert_logical(NULL, "x"))
  expect_error(assert_logical(c(TRUE, NA), "x"))
})

test_that("assert_single passes good input", {
  expect_true(assert_single(4L, "y"))
  expect_true(assert_single(-10, "y"))
  expect_true(assert_single("a", "y"))
  expect_true(assert_single(list(list(a = "a", b = 4)), "y"))
  expect_true(assert_single(TRUE, "y"))
})

test_that("assert_single failes bad input", {
  expect_error(assert_single(c(), "y"))
  expect_error(assert_single(list(), "y"))
  expect_error(assert_single(list(a = list(a = "a"), b = list(4)), "y"))
  expect_error(assert_single(NULL, "y"))
  expect_error(assert_single(c(1, -1), "y"))
})

test_that("assert_length passes good input", {
  expect_true(assert_length(4L, 1, "y"))
  expect_true(assert_length(c(4, 1), 2, "y"))
  expect_true(assert_length(c("b", "ca"), 2, "y"))
})

test_that("assert_length failes bad input", {
  expect_error(assert_length(0, 0, "y"))
  expect_error(assert_length(0, 2, "y"))
  expect_error(assert_length(-1L, 8, "y"))
  expect_error(assert_length("a", 10, "y"))
  expect_error(assert_length(NULL, 1, "y"))
  expect_error(assert_length(c(1, -1), 5, "y"))
})

test_that("assert_infinite passes good input", {
  expect_true(assert_infinite(Inf, "y"))
  expect_true(assert_infinite(c(Inf, -Inf), "y"))
})

test_that("assert_infinite failes bad input", {
  expect_error(assert_infinite(c(0, Inf), "w"))
  expect_error(assert_infinite(-1L, "w"))
  expect_error(assert_infinite("a", "w"))
  expect_error(assert_infinite(NULL, "w"))
  expect_error(assert_infinite(c(1, -1), "w"))
})

test_that("assert_ge passes good input", {
  expect_true(assert_ge(0, 0, "w"))
  expect_true(assert_ge(4L, 3, "w"))
  expect_true(assert_ge(c(4, 0.1), -1, "w"))
  expect_true(assert_ge(c(4, 0.1), c(4, -1), "w"))
})

test_that("assert_ge failes bad input", {
  expect_error(assert_ge(0, 0.1, "w"))
  expect_error(assert_ge(-1L, 0.5, "w"))
  expect_error(assert_ge("a", 0, "w"))
  expect_error(assert_ge(NULL, NA, "w"))
  expect_error(assert_ge(c(1, -1), c(1.1, -2), "w"))
})

test_that("assert_gt passes good input", {
  expect_true(assert_gt(4L, 3, "w"))
  expect_true(assert_gt(c(4, 0.1), -1, "w"))
  expect_true(assert_gt(c(4, 0.1), c(3, -1), "w"))
})

test_that("assert_gt failes bad input", {
  expect_error(assert_gt(0, 0, "w"))
  expect_error(assert_gt(-1L, 0.5, "w"))
  expect_error(assert_gt("a", 0, "w"))
  expect_error(assert_gt(NULL, NA, "w"))
  expect_error(assert_gt(c(1, -1), c(1, -2), "w"))
})

test_that("assert_lt passes good input", {
  expect_true(assert_lt(4L, 5, "w"))
  expect_true(assert_lt(c(4, 0.1), 10, "w"))
  expect_true(assert_lt(c(-4, 0.1), c(3, 1), "w"))
})

test_that("assert_lt failes bad input", {
  expect_error(assert_lt(0, 0, "w"))
  expect_error(assert_lt(1L, -0.5, "w"))
  expect_error(assert_lt("a", 0, "w"))
  expect_error(assert_lt(NULL, NA, "w"))
  expect_error(assert_lt(c(1, -10), c(1, -2), "w"))
})

test_that("assert_positive passes good input", {
  expect_true(assert_positive(4L, "w"))
  expect_true(assert_positive(c(4, 0.1), "w"))
})

test_that("assert_positive failes bad input", {
  expect_error(assert_positive(0, "w"))
  expect_error(assert_positive(-1L, "w"))
  expect_error(assert_positive("a", "w"))
  expect_error(assert_positive(NULL, "w"))
  expect_error(assert_positive(c(1, -1), "w"))
})

test_that("assert_nonnegative passes good input", {
  expect_true(assert_nonnegative(0, "w"))
  expect_true(assert_nonnegative(4L, "w"))
  expect_true(assert_nonnegative(c(4, 0.1), "w"))
})

test_that("assert_nonnegative failes bad input", {
  expect_error(assert_nonnegative(-1L, "w"))
  expect_error(assert_nonnegative("a", "w"))
  expect_error(assert_nonnegative(NULL, "w"))
  expect_error(assert_nonnegative(c(1, -1), "w"))
})

test_that("assert_element passes good input", {
  expect_true(assert_element(c("a", "aa"), c("aa", "a", "aaa"), "xx", "yy"))
  expect_true(assert_element("a", "a", "xx", "yy"))
  expect_true(assert_element(1L, seq_len(5), "xx", "yy"))
})

test_that("assert_element failes bad input", {
  expect_error(assert_element(0, Inf, "xx", "yy"))
  expect_error(assert_element(-1L, seq_len(15), "xx", "yy"))
  expect_error(assert_element("1", seq_len(15), "xx", "yy"))
  expect_error(assert_element("a", c(), "xx", "yy"))
  expect_error(assert_element(c("a", "aa"), c("a", "aaa"), "xx", "yy"))
  expect_error(assert_element(c("a", "aa"), NULL, "xx", "yy"))
})

