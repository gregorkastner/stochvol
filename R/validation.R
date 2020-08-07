# Functions that assert some properties of the inputs

sv_stop <- function (what, name, should, and_not) {
  stop(name, " = ", prettify(what), " should ", should, "; not ", and_not)
}

prettify <- function (x) {
  out <- paste0("c(", paste(x, collapse = ", "), ")")
  if (nchar(out) > 40L) {
    out <- paste(substr(out, start = 1, stop = 40L), "...)")
  }
  out
}

assert_numeric <- function (x, name) {
  if (!(isTRUE(is.numeric(x)) &&
        isTRUE(all(!is.na(x))))) {
    sv_stop(what = x, name = name, should = "be numeric", and_not = typeof(x))
  }
  TRUE
}

assert_single <- function (x, name) {
  assert_length(x, 1L, name)
}

assert_length <- function (x, len, name) {
  if (!isTRUE(length(x) == len)) {
    sv_stop(what = x, name = name, should = paste("have length", len), and_not = length(x))
  }
  TRUE
}

assert_infinite <- function (x, name) {
  assert_numeric(x, name)

  if (!isTRUE(all(is.infinite(x)))) {
    sv_stop(what = x, name = name, should = "be infinite", and_not = x)
  }
  TRUE
}

assert_positive <- function (x, name) {
  assert_numeric(x, name)

  if (!isTRUE(all(x > 0))) {
    sv_stop(what = x, name = name, should = "be positive", and_not = x)
  }
  TRUE
}

assert_element <- function (x, v, name_x, name_v) {
  if (!(isTRUE(typeof(x) == typeof(v)) &&
        isTRUE(all(x %in% v)))) {
    sv_stop(what = x, name = name_x, should = paste("be a subset of", name_v, "=", prettify(v)), and_not = "otherwise")
  }
  TRUE
}

