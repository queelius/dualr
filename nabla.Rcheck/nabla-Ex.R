pkgname <- "nabla"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "nabla-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('nabla')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("D")
### * D

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: D
### Title: Composable total derivative operator
### Aliases: D

### ** Examples

f <- function(x) x[1]^2 * x[2]
D(f, c(3, 4))
D(f, c(3, 4), order = 2)

g <- function(x) list(x[1] * x[2], x[1]^2)
D(g, c(2, 3))

Df <- D(f)
DDf <- D(Df)
DDf(c(3, 4))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("D", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("beta")
### * beta

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: beta
### Title: Beta function for dual numbers
### Aliases: beta beta,ANY,ANY-method beta,numeric,numeric-method

### ** Examples

beta(2, 3)
a <- dual_variable(2)
value(beta(a, 3))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("beta", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("deriv")
### * deriv

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: deriv
### Title: Extract the derivative (tangent) part of a dual number
### Aliases: deriv deriv,dualr-method deriv,numeric-method

### ** Examples

deriv(dual(3, 1))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("deriv", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("deriv_n")
### * deriv_n

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: deriv_n
### Title: Extract the k-th derivative from a nested dual result
### Aliases: deriv_n

### ** Examples

x <- dual_variable_n(1, order = 3)
r <- exp(x)
deriv_n(r, 0)
deriv_n(r, 1)
deriv_n(r, 2)
deriv_n(r, 3)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("deriv_n", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("differentiate_n")
### * differentiate_n

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: differentiate_n
### Title: Compute a function value and all derivatives up to order n
### Aliases: differentiate_n

### ** Examples

differentiate_n(sin, pi/4, order = 4)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("differentiate_n", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dual-arithmetic")
### * dual-arithmetic

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dual-arithmetic
### Title: Arithmetic and comparison operators for dual numbers
### Aliases: dual-arithmetic Ops,dualr,dualr-method
###   Ops,dualr,numeric-method Ops,numeric,dualr-method
###   +,dualr,dualr-method +,dualr,numeric-method +,numeric,dualr-method
###   -,dualr,dualr-method -,dualr,numeric-method -,numeric,dualr-method
###   *,dualr,dualr-method *,dualr,numeric-method *,numeric,dualr-method
###   /,dualr,dualr-method /,dualr,numeric-method /,numeric,dualr-method
###   ^,dualr,dualr-method ^,dualr,numeric-method ^,numeric,dualr-method
###   +,dualr,missing-method -,dualr,missing-method !,dualr-method

### ** Examples

x <- dual_variable(3)
y <- dual_variable(4)

value(x + y)
deriv(x * x)
value(x^2)
deriv(x^2)

x < y
x == y




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dual-arithmetic", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dual-atan2")
### * dual-atan2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dual-atan2
### Title: Two-argument arctangent for dual numbers
### Aliases: dual-atan2 atan2,dualr,dualr-method atan2,dualr,numeric-method
###   atan2,numeric,dualr-method

### ** Examples

y <- dual_variable(1)
x <- dual_constant(1)
result <- atan2(y, x)
value(result)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dual-atan2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dual-coerce")
### * dual-coerce

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dual-coerce
### Title: Coerce dual to numeric
### Aliases: dual-coerce as.numeric,dualr-method

### ** Examples

x <- dual(3.14, 1)
as.numeric(x)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dual-coerce", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dual-combine")
### * dual-combine

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dual-combine
### Title: Combine dual numbers into a dual_vector
### Aliases: dual-combine c,dualr-method

### ** Examples

x <- dual_variable(1)
y <- dual_variable(2)
dv <- c(x, y)
length(dv)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dual-combine", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dual-is-numeric")
### * dual-is-numeric

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dual-is-numeric
### Title: Check if a dual number is numeric
### Aliases: dual-is-numeric is.numeric,dualr-method

### ** Examples

is.numeric(dual(1, 0))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dual-is-numeric", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dual-log")
### * dual-log

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dual-log
### Title: Logarithm with optional base for dual numbers
### Aliases: dual-log log,dualr-method

### ** Examples

x <- dual_variable(8)
value(log(x, base = 2))
deriv(log(x, base = 2))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dual-log", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dual-math")
### * dual-math

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dual-math
### Title: Math group generic for dual numbers
### Aliases: dual-math Math,dualr-method exp,dualr-method sqrt,dualr-method

### ** Examples

x <- dual_variable(pi / 4)
value(sin(x))
deriv(sin(x))

y <- dual_variable(2)
value(exp(y))
deriv(exp(y))
deriv(log(y))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dual-math", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dual-math2")
### * dual-math2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dual-math2
### Title: Math2 group generic for dual numbers
### Aliases: dual-math2 Math2,dualr-method

### ** Examples

x <- dual_variable(3.14159)
value(round(x, 2))
deriv(round(x, 2))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dual-math2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dual-maxmin")
### * dual-maxmin

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dual-maxmin
### Title: Piecewise max and min for dual numbers
### Aliases: dual-maxmin max,dualr-method min,dualr-method

### ** Examples

x <- dual_variable(3)
y <- dual_variable(5)
value(max(x, y))
value(min(x, y))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dual-maxmin", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dual-show")
### * dual-show

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dual-show
### Title: Display a dual number
### Aliases: dual-show show,dualr-method show,dual_vector-method

### ** Examples

x <- dual(3, 1)
x

dv <- dual_vector(dual(1, 0), dual(2, 1))
dv




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dual-show", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dual-summary")
### * dual-summary

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dual-summary
### Title: Summary group generic for dual numbers
### Aliases: dual-summary Summary,dualr-method

### ** Examples

x <- dual_variable(2)
y <- dual_variable(5)
value(sum(x, y))
value(prod(x, y))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dual-summary", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dual")
### * dual

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dual
### Title: Create a dual number
### Aliases: dual

### ** Examples

x <- dual(3, 1)
value(x)
deriv(x)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dual", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dual_constant")
### * dual_constant

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dual_constant
### Title: Create a dual constant (derivative seed = 0)
### Aliases: dual_constant

### ** Examples

k <- dual_constant(5)
deriv(k)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dual_constant", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dual_constant_n")
### * dual_constant_n

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dual_constant_n
### Title: Create a constant dual for n-th order differentiation
### Aliases: dual_constant_n

### ** Examples

k <- dual_constant_n(5, order = 3)
deriv_n(k, 1)
deriv_n(k, 2)
deriv_n(k, 3)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dual_constant_n", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dual_variable")
### * dual_variable

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dual_variable
### Title: Create a dual variable (derivative seed = 1)
### Aliases: dual_variable

### ** Examples

x <- dual_variable(2)
deriv(x^2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dual_variable", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dual_variable_n")
### * dual_variable_n

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dual_variable_n
### Title: Create a dual seeded for n-th order differentiation
### Aliases: dual_variable_n

### ** Examples

x <- dual_variable_n(2, order = 3)
r <- x^4
deriv_n(r, 3)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dual_variable_n", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dual_vector-access")
### * dual_vector-access

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dual_vector-access
### Title: Indexing and length for dual_vector
### Aliases: dual_vector-access [,dual_vector,numeric-method
###   length,dual_vector-method

### ** Examples

dv <- dual_vector(dual(10, 1), dual(20, 0), dual(30, 0))
value(dv[1])
length(dv)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dual_vector-access", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dual_vector")
### * dual_vector

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dual_vector
### Title: Create a vector of dual numbers
### Aliases: dual_vector

### ** Examples

dv <- dual_vector(dual(1, 0), dual(2, 1))
length(dv)
value(dv[1])



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dual_vector", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("erf")
### * erf

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: erf
### Title: Error function
### Aliases: erf erf,numeric-method erf,dualr-method

### ** Examples

erf(1)
x <- dual_variable(1)
value(erf(x))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("erf", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("erfc")
### * erfc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: erfc
### Title: Complementary error function
### Aliases: erfc erfc,numeric-method erfc,dualr-method

### ** Examples

erfc(1)
x <- dual_variable(0)
value(erfc(x))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("erfc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gradient")
### * gradient

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gradient
### Title: Compute the gradient of a scalar-valued function
### Aliases: gradient

### ** Examples

f <- function(x) -(x[1] - 3)^2 - (x[2] - 5)^2
gradient(f, c(1, 2))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gradient", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("hessian")
### * hessian

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: hessian
### Title: Compute the Hessian of a scalar-valued function
### Aliases: hessian

### ** Examples

f <- function(x) -(x[1] - 3)^2 - (x[2] - 5)^2
hessian(f, c(1, 2))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("hessian", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("is_dual")
### * is_dual

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: is_dual
### Title: Test whether an object is a dual number
### Aliases: is_dual

### ** Examples

is_dual(dual(1, 0))
is_dual(42)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("is_dual", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("jacobian")
### * jacobian

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: jacobian
### Title: Compute the Jacobian of a vector-valued function
### Aliases: jacobian

### ** Examples

f <- function(x) {
  a <- x[1]; b <- x[2]
  list(a * b, a^2, sin(b))
}
jacobian(f, c(2, pi/4))

g <- function(x) x[1]^2 + x[2]^2
jacobian(g, c(3, 4))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("jacobian", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lbeta")
### * lbeta

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lbeta
### Title: Log-beta function for dual numbers
### Aliases: lbeta lbeta,numeric,numeric-method lbeta,dualr,dualr-method
###   lbeta,dualr,numeric-method lbeta,numeric,dualr-method

### ** Examples

lbeta(2, 3)
a <- dual_variable(2)
value(lbeta(a, 3))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lbeta", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("psigamma")
### * psigamma

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: psigamma
### Title: Polygamma function for dual numbers
### Aliases: psigamma psigamma,numeric-method psigamma,dualr-method

### ** Examples

psigamma(1, deriv = 0)
x <- dual_variable(2)
value(psigamma(x, deriv = 1))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("psigamma", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("value")
### * value

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: value
### Title: Extract the value (primal) part of a dual number
### Aliases: value value,dualr-method value,numeric-method

### ** Examples

value(dual(3, 1))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("value", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
