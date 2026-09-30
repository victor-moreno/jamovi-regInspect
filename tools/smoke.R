# Smoke test, run by tools/install.sh right after installing the module into
# jamovi desktop and/or the docker container (the runner has put jamovi's
# libraries, jmv's included, and the installed module on .libPaths() and
# attached it). Moved here from the docker target of the previous install.sh.

# model coefficients against a plain lm()
data(mtcars)
results <- regInspect::regInspect(
    data = mtcars, dep = mpg, covs = c("wt", "hp"),
    blocks = list(list("wt", "hp"))
)
tbl <- results$models[[1]]$coef$asDF
oracle <- coef(lm(mpg ~ wt + hp, data = mtcars))
stopifnot(isTRUE(all.equal(tbl$est, unname(oracle), tolerance = 1e-6)))
cat(sprintf("   coefficients smoke test passed: wt = %.3f, hp = %.4f\n",
            tbl$est[2], tbl$est[3]))

# partial-residual points on the emm plot
options <- regInspect:::regInspectOptions$new(
    dep = "mpg", covs = "wt", blocks = list(list("wt")),
    emMeans = list(list("wt")), emmPoints = "partial"
)
analysis <- regInspect:::regInspectClass$new(options = options, data = mtcars)
analysis$run()

priv <- analysis$.__enclos_env__$private
pts <- priv$.emmPointData(1, c("wt", NA, NA), c(jmvcore::toB64("wt"), NA, NA),
                           list(x = jmvcore::toB64("wt")))
stopifnot(isTRUE(all.equal(sort(pts$.y), sort(mtcars$mpg), tolerance = 1e-6)))
cat("   partial-residual smoke test passed (single predictor == raw y)\n")
