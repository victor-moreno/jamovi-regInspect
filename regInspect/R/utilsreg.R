#' @import jmvcore
NULL

#' Update reference levels for a set of variables
#'
#' @param data The data frame containing the variables
#' @param vars The names of the variables to update
#' @param refLevels A list of reference levels to use for each variable
#' @return A list of updated reference levels and a list of variables that
#'   had their reference levels changed
#' @keywords internal
getReferenceLevels = function(data, vars, refLevels) {
    updatedRefLevels <- list()
    changedVars <- c()

    # Create a named list from the refLevels input for easier access
    refLevelsList <- stats::setNames(
        lapply(refLevels, function(ref) ref$ref),
        sapply(refLevels, function(ref) ref$var)
    )

    for (var in vars) {
        factorLevels <- levels(data[[var]])
        refLevel <- refLevelsList[[var]]

        # If no refLevel is provided or the provided level is invalid, use the first level
        if (is.null(refLevel) || ! (refLevel %in% factorLevels)) {
            refLevel <- factorLevels[1]
            changedVars <- c(changedVars, var)
        }

        updatedRefLevels[[ length(updatedRefLevels) + 1 ]] <- list(
            var = var, ref = refLevel
        )
    }

    return(list(refLevels=updatedRefLevels, changedVars=changedVars))
}


#' Set a warning notice for reference level changes
#'
#' @param self The analysis object
#' @param changedVars The variables that had their reference levels changed
#' @keywords internal
setRefLevelWarning = function(self, changedVars) {
    message <- jmvcore::format(
        .("The specified reference level was not found for the following variable(s): {vars}. Defaulting to the first available level. To use a custom reference level, ensure the defined reference level is available in the data."),
        vars=listItems(self, changedVars)
    )

    setAnalysisNotice(
        self, message, name="refLevelWarning", type=jmvcore::NoticeType$STRONG_WARNING
    )
}

#' Build a natural-language list ('a', 'a and b', 'a, b, and c')
#'
#' @param self The analysis object
#' @param items The items to list
#' @param quote The quote character to wrap each item in
#' @keywords internal
listItems <- function(self, items, quote = '\'') {

    items <- vapply(items, function(x) paste0(quote, x, quote), '', USE.NAMES=FALSE)

    if (length(items) == 1) {

        l <- items

    } else if (length(items) == 2) {

        l <- jmvcore::format(.('{item1} and {item2}'), item1=items[1], item2=items[2])

    } else {

        l <- items[1]
        for (i in seq(2, length(items)-1))
            l <- jmvcore::format(.('{list}, {nextItem}'), list=l, nextItem=items[i])
        l <- jmvcore::format(.('{list}, and {lastItem}'), list=l, lastItem=items[length(items)])
    }

    return(l)
}

#' Set a global analysis notice
#'
#' @param self The analysis object
#' @param message The message to set
#' @param name The name of the notice
#' @param type The type of the notice
#' @keywords internal
setAnalysisNotice = function(self, message, name, type) {
    notice <- jmvcore::Notice$new(
        options = self$options,
        name = 'warningMessage',
        type = jmvcore::NoticeType$STRONG_WARNING
    )
    notice$setContent(message)
    self$results$insert(1, notice)
}
