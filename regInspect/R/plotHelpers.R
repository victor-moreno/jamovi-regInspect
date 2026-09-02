# Pure statistics/formatting helpers for regInspect's descriptive plots.
# Adapted from corrInspect's R/corrCompute.R (same author) -- trimmed to
# just what the Plots tab's annotation needs (Pearson r/p only; the
# descriptive scatters always use an OLS fit, so a rank correlation doesn't
# apply here the way it does in corrInspect's general-purpose scatter).
# No jmvcore, no R6, no self -- testable standalone against base R.

pearsonFit <- function(x, y) {
    ok <- stats::complete.cases(x, y)
    x <- x[ok]
    y <- y[ok]

    out <- list(r = NA_real_, n = length(x), p = NA_real_)
    if (length(x) < 3 || stats::sd(x) == 0 || stats::sd(y) == 0)
        return(out)

    fit <- tryCatch(stats::cor.test(x, y, method = 'pearson'), error = function(e) NULL)
    if (is.null(fit))
        return(out)

    out$r <- unname(fit$estimate)
    out$p <- fit$p.value
    out
}

# Compact p-value text for a plot annotation (not jamovi's own pvalue
# format, which isn't reachable from plain R).
pText <- function(p) {
    if (is.na(p))
        return('')
    if (p < .001)
        return('p<.001')
    sprintf('p=%.3f', p)
}

# ggtheme (jmvcore's renderFun 3rd-arg-equivalent) is theme + 2 discrete
# palette scales, packed as a list -- `p + ggtheme` is only safe when the
# plot declares no scale of its own. Pull just the theme back out for the
# annotation text colour and for panels (the scatterplot matrix) that need
# their own scale.
themeOnly <- function(ggtheme) {
    if (is.null(ggtheme))
        return(NULL)
    if (inherits(ggtheme, 'theme'))
        return(ggtheme)
    parts <- Filter(function(x) inherits(x, 'theme'), ggtheme)
    if (length(parts) == 0L) NULL else parts[[1L]]
}
