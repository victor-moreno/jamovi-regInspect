#' @importFrom jmvcore .
regInspectClass <- R6::R6Class(
    "regInspectClass",
    inherit = regInspectBase,
    #### Active bindings ----
    active = list(
        dataProcessed = function() {
            if (is.null(private$.dataProcessed))
                private$.dataProcessed <- private$.cleanData()

            return(private$.dataProcessed)
        },
        weights = function() {
            if (is.null(private$.weights))
                private$.weights <- private$.computeWeights()

            return(private$.weights)
        },
        models = function() {
            if (is.null(private$.models))
                private$.models <- private$.computeModels()

            return(private$.models)
        },
        modelsScaled = function() {
            if (is.null(private$.modelsScaled))
                private$.modelsScaled <- private$.computeModels(scaled = TRUE)

            return(private$.modelsScaled)
        },
        nModels = function() {
            if (is.null(private$.nModels))
                private$.nModels <- length(self$options$blocks)

            return(private$.nModels)
        },
        residuals = function() {
            if (is.null(private$.residuals))
                private$.residuals <- private$.computeResiduals()

            return(private$.residuals)
        },
        fitted = function() {
            if (is.null(private$.fitted))
                private$.fitted <- private$.computeFitted()

            return(private$.fitted)
        },
        predicted = function() {
            if (is.null(private$.predicted))
                private$.predicted <- private$.computePredicted()

            return(private$.predicted)
        },
        cooks = function() {
            if (is.null(private$.cooks))
                private$.cooks <- private$.computeCooks()

            return(private$.cooks)
        },
        mahal = function() {
            if (is.null(private$.mahal))
                private$.mahal <- private$.computeMahal()

            return(private$.mahal)
        },
        anovaModelComparison = function() {
            if (is.null(private$.anovaModelComparison))
                private$.anovaModelComparison <- do.call(stats::anova, self$models)

            return(private$.anovaModelComparison)
        },
        anovaModelTerms = function() {
            if (is.null(private$.anovaModelTerms))
                private$.anovaModelTerms <- private$.computeAnovaModelTerms()

            return(private$.anovaModelTerms)
        },
        AIC = function() {
            if (is.null(private$.AIC))
                private$.AIC <- private$.computeAIC()

            return(private$.AIC)
        },
        BIC = function() {
            if (is.null(private$.BIC))
                private$.BIC <- private$.computeBIC()

            return(private$.BIC)
        },
        CICoefEst = function() {
            if (is.null(private$.CICoefEst))
                private$.CICoefEst <- private$.computeCICoefEst()

            return(private$.CICoefEst)
        },
        CICoefEstScaled = function() {
            if (is.null(private$.CICoefEstScaled))
                private$.CICoefEstScaled <- private$.computeCICoefEst(scaled = TRUE)

            return(private$.CICoefEstScaled)
        },
        VIF = function() {
            if (is.null(private$.VIF))
                private$.VIF <- private$.computeVIF()

            return(private$.VIF)
        },
        durbinWatson = function() {
            if (is.null(private$.durbinWatson))
                private$.durbinWatson <- private$.computeDurbinWatson()

            return(private$.durbinWatson)
        },
        emMeans = function() {
            if (is.null(private$.emMeans))
                private$.emMeans <- private$.computeEmMeans()

            return(private$.emMeans)
        },
        refLevels = function() {
            if (is.null(private$.refLevels)) {
                refLevels <- getReferenceLevels(
                    self$data, self$options$factors, self$options$refLevels
                )
                private$.refLevels <- refLevels$refLevels

                if (length(refLevels$changedVars) > 0)
                    setRefLevelWarning(self, refLevels$changedVars)
            }

            return(private$.refLevels)
        }
    ),
    private = list(
        #### Member variables ----
        .dataProcessed = NULL,
        .dataRowNums = NULL,
        .weights = NULL,
        .models = NULL,
        .modelsScaled = NULL,
        .isAliased = NULL,
        .nModels = NULL,
        .residuals = NULL,
        .fitted = NULL,
        .predicted = NULL,
        .cooks = NULL,
        .mahal = NULL,
        .anovaModelComparison = NULL,
        .anovaModelTerms = NULL,
        .AIC = NULL,
        .BIC = NULL,
        .CICoefEst = NULL,
        .CICoefEstScaled = NULL,
        .VIF = NULL,
        .durbinWatson = NULL,
        .modelTerms = NULL,
        .rowNamesModel = NULL,
        .emMeans = NULL,
        .emMeansForPlot = NULL,
        .refLevels = NULL,

        #### Init + run functions ----
        .init = function() {
            private$.initModelFitTable()
            private$.initModelCompTable()
            private$.initModelSpec()
            private$.initAnovaTables()
            private$.initCoefTable()
            private$.initMahalTable()
            private$.initCollinearityTable()
            private$.initResPlots()
            private$.initEmm()
            private$.initEmmTable()
            private$.initOutputs()
            private$.initPlots()
        },
        .run = function() {
            if (is.null(self$options$dep) || self$nModels < 1 || length(self$options$blocks[[1]]) == 0)
                return()

            if (any(unlist(private$.getIsAliased())))
                setSingularityWarning(self)

            private$.populateModelFitTable()
            private$.populateModelCompTable()
            private$.populateAnovaTables()
            private$.populateCoefTables()

            private$.populateCooksTable()
            private$.populateMahalTable()
            private$.populateDurbinWatsonTable()
            private$.populateCollinearityTable()
            private$.populateNormality()
            private$.populateEmmTables()

            private$.prepareQQPlot()
            private$.prepareResPlots()
            private$.prepareEmmPlots()

            private$.populateOutputs()
        },

        #### Compute results ----
        .computeModels = function(modelNo = NULL, scaled = FALSE) {
            if (scaled)
                data <- private$.scaleData(self$dataProcessed)
            else
                data <- self$dataProcessed

            formulas <- private$.formulas()

            if (is.numeric(modelNo))
                formulas <- formulas[modelNo]

            models <- list()
            for (i in seq_along(formulas))
                models[[i]] <- lm(formulas[[i]], data=data, weights=self$weights)

            return(models)
        },
        .computeWeights = function() {
            legacy_weights <- self$options$weights
            global_weights <- attr(self$data, "jmv-weights")

            if (is.null(legacy_weights) && is.null(global_weights))
                return()

            if (! is.null(legacy_weights)) {
                weights <- self$dataProcessed[[jmvcore::toB64(legacy_weights)]]
            } else {
                weights <- self$dataProcessed[[".WEIGHTS"]]
            }

            if (any(weights < 0)) {
                jmvcore::reject(
                    .("'{var}' contains negative values. Negative weights are not permitted."),
                    var=legacy_weights
                )
            }

            return(weights)
        },
        .computeAnovaModelTerms = function() {
            anovaTerms <- list()
            for (i in seq_along(self$models)) {
                suppressMessages({
                    anovaTerms[[i]] <- car::Anova(
                        self$models[[i]], type=3, singular.ok=TRUE
                    )
                })
            }

            return(anovaTerms)
        },
        .computeResiduals = function() {
            res <- list()
            for (i in seq_along(self$models))
                res[[i]] <- self$models[[i]]$residuals

            return(res)
        },
        .computeFitted = function() {
            fitted <- list()
            for (i in seq_along(self$models))
                fitted[[i]] <- self$models[[i]]$fitted.values

            return(fitted)
        },
        .computePredicted = function() {
            data <- private$.cleanData(naSkip=jmvcore::toB64(self$options$dep))
            fitted <- list()
            for (i in seq_along(self$models))
                fitted[[i]] <- predict(self$models[[i]], data)

            return(fitted)
        },
        .computeCooks = function() {
            cooks <- list()
            for (i in seq_along(self$models))
                cooks[[i]] <- stats::cooks.distance(self$models[[i]])

            return(cooks)
        },
        .computeCooksSummary = function() {
            cooksSummary <- list()
            for (i in seq_along(self$cooks)) {
                cooks <- self$cooks[[i]]
                cooksSummary[[i]] <- list(
                    'mean' = mean(cooks),
                    'median' = median(cooks),
                    'sd' = sd(cooks),
                    'min' = min(cooks),
                    'max' = max(cooks)
                )
            }

            return(cooksSummary)
        },
        .computeMahal = function() {
            mahal <- vector(mode = "list", length(self$nModels))
            data <- self$dataProcessed
            for (i in seq_along(self$models)) {
                terms <- attr(self$models[[i]]$terms, "term.labels")
                covs <- intersect(terms, jmvcore::toB64(self$options$covs))
                # skip if number of variables less than 2
                if (length(covs) >= 2) {
                    mahalChiSq <- stats::mahalanobis(
                        data[, covs], colMeans(data[, covs]), cov(data[, covs])
                    )
                    mahal[[i]] <- data.frame(
                        chisq = mahalChiSq,
                        p = pchisq(mahalChiSq, length(covs), lower.tail = FALSE),
                        row.names = rownames(data)
                    )
                }
            }

            return(mahal)
        },
        .computeAIC = function() {
            AIC <- list()
            for (i in seq_along(self$models))
                AIC[[i]] <- stats::AIC(self$models[[i]])

            return(AIC)
        },
        .computeBIC = function() {
            BIC <- list()
            for (i in seq_along(self$models))
                BIC[[i]] <- stats::BIC(self$models[[i]])

            return(BIC)
        },
        .computeCICoefEst = function(scaled = FALSE) {
            if (scaled)
                models <- self$modelsScaled
            else
                models <- self$models

            CICoefEst <- list()
            for (i in seq_along(models))
                CICoefEst[[i]] <- stats::confint(models[[i]], level = self$options$ciWidth / 100)

            return(CICoefEst)
        },
        .computeVIF = function() {
            modelTerms <- private$.getModelTerms()

            VIF <- list()
            for (i in seq_along(self$models)) {
                if ( ! private$.getIsAliased()[[i]] && length(modelTerms[[i]]) > 1 ) {
                    VIF[[i]] <- car::vif(self$models[[i]])
                } else {
                    VIF[[i]] <- c(VIF, list(NULL))
                }
            }

            return(VIF)
        },
        .computeDurbinWatson = function() {
            durbinWatson <- list()
            for (i in seq_along(self$models)) {
                durbinWatson[[i]] <- car::durbinWatsonTest(self$models[[i]])
            }

            return(durbinWatson)
        },
        .computeEmMeans = function(forPlot = FALSE) {
            covs <- self$options$covs
            termsAll <- private$.getModelTerms()
            emMeansTerms <- self$options$emMeans

            emMeans <- list()
            for (i in seq_along(termsAll)) {
                terms <- unique(unlist(termsAll[[i]]))
                model <- self$models[[i]]

                emmTable <- list()
                for (j in seq_along(emMeansTerms)) {
                    term <- emMeansTerms[[j]]

                    if ( ! is.null(term) && all(term %in% terms)) {
                        termB64 <- jmvcore::toB64(term)

                        FUN <- list()
                        for (k in seq_along(termB64)) {
                            if (term[k] %in% covs) {
                                if (forPlot && k == 1) {
                                    FUN[[termB64[k]]] <- function(x) pretty(x, 25)
                                } else {
                                    FUN[[termB64[k]]] <- function(x) c(mean(x)-sd(x), mean(x), mean(x)+sd(x))
                                }
                            }
                        }

                        formula <- formula(paste('~', jmvcore::composeTerm(termB64)))

                        if (self$options$emmWeights)
                            weights <- 'equal'
                        else
                            weights <- 'cells'

                        suppressMessages({
                            emmeans::emm_options(sep=",", parens="a^", cov.keep=1)

                            mm <- try(
                                emmeans::emmeans(
                                    model,
                                    formula,
                                    cov.reduce = FUN,
                                    options = list(level=self$options$ciWidthEmm / 100),
                                    weights = weights,
                                    data = self$dataProcessed,
                                    non.nuis = all.vars(formula),
                                ),
                                silent = TRUE
                            )

                            emmTable[[j]] <- try(as.data.frame(summary(mm)), silent = TRUE)
                        })
                    }
                }

                emMeans[[i]] <- emmTable
            }

            return(emMeans)
        },

        #### Init tables/plots functions ----
        .initModelFitTable = function() {
            table <- self$results$modelFit

            for (i in seq_along(self$options$blocks))
                table$addRow(rowKey=i, values=list(model = i))

            table$setNote(
                "n",
                jmvcore::format(
                    .("Models estimated using sample size of N={n}"), n="..."
                )
            )
        },
        .initModelCompTable = function() {
            table <- self$results$modelComp
            blocks <- self$options$blocks

            if (length(blocks) <= 1) {
                table$setVisible(visible = FALSE)
                return()
            }

            for (i in 1:(length(blocks)-1))
                table$addRow(rowKey=i, values=list(model1 = i, model2 = as.integer(i+1)))
        },
        .initModelSpec = function() {
            groups <- self$results$models

            for (i in seq_along(self$options$blocks)) {
                groups$addItem(key=i)
                group <- groups$get(key=i)
                group$setTitle(paste("Model",i))
            }
        },
        .initAnovaTables = function() {
            groups <- self$results$models
            termsAll <- private$.getModelTerms()

            for (i in seq_along(termsAll)) {

                table <- groups$get(key=i)$anova
                terms <- termsAll[[i]]

                for (j in seq_along(terms))
                    table$addRow(rowKey=paste0(terms[[j]]), values=list(term = jmvcore::stringifyTerm(terms[[j]])))

                table$addRow(rowKey='.RES', values=list(term = .('Residuals')))
                table$addFormat(col=1, rowKey='.RES', format=Cell.BEGIN_GROUP)

                table$setNote("ss", .("Type 3 sum of squares"))
            }
        },
        .initCoefTable = function() {
            groups <- self$results$models
            factors <- self$options$factors
            weights <- self$options$weights
            termsAll <- private$.getModelTerms()
            rowNamesModel <- private$.getRowNamesModel()

            for (i in seq_along(termsAll)) {
                table <- groups$get(key=i)$coef

                ciWidth <- self$options$ciWidth
                ciWidthTitle <- jmvcore::format(.('{ciWidth}% Confidence Interval'), ciWidth=ciWidth)
                table$getColumn('lower')$setSuperTitle(ciWidthTitle)
                table$getColumn('upper')$setSuperTitle(ciWidthTitle)

                ciWidthStdEst <- self$options$ciWidthStdEst
                ciWidthStdEstTitle <- jmvcore::format(.('{ciWidth}% Confidence Interval'), ciWidth=ciWidthStdEst)
                table$getColumn('stdEstLower')$setSuperTitle(ciWidthStdEstTitle)
                table$getColumn('stdEstUpper')$setSuperTitle(ciWidthStdEstTitle)

                coefTerms <- rowNamesModel[[i]]

                table$addRow(rowKey="`(Intercept)`", values=list(term = .("Intercept")))

                if (! is.null(weights)) {
                    private$.weightsName <- weights
                    table$setNote(
                        "weights",
                        jmvcore::format(.("Weighted by '{varName}'"), varName=weights)
                    )
                }

                if ( ! is.null(factors)) {
                    note <- ifelse(self$options$intercept == 'refLevel',
                                   .('Represents reference level'),
                                   .('Represents grand mean'))
                    table$addFootnote(rowKey="`(Intercept)`", 'term', note)
                }

                terms <- termsAll[[i]]

                iter <- 1
                for (j in seq_along(terms)) {
                    if (any(terms[[j]] %in% factors)) {
                        table$addRow(rowKey=terms[[j]],
                                     values=list(term = paste0(jmvcore::stringifyTerm(terms[[j]]), ':'),
                                                 est='', se='', t='', p='', lower='', upper='', stdEst='',
                                                 stdEstLower='', stdEstUpper=''))

                        contrastNames <- private$.contrastsCoefTable(terms[[j]])

                        for (k in seq_along(contrastNames)) {
                            iter <-iter + 1
                            rowKey <- jmvcore::composeTerm(coefTerms[[iter]])
                            table$addRow(rowKey=rowKey, values=list(term = contrastNames[[k]]))
                            table$addFormat(rowKey=rowKey, col=1, Cell.INDENTED)
                        }
                    } else {
                        iter <- iter + 1
                        rowKey <- jmvcore::composeTerm(jmvcore::toB64(terms[[j]]))
                        table$addRow(rowKey=rowKey, values=list(term = jmvcore::stringifyTerm(terms[[j]])))
                    }
                }
            }
        },
        .initMahalTable = function() {
            if ( ! self$options$mahal)
                return()

            groups <- self$results$models
            termsAll <- private$.getModelTerms()
            covs <- self$options$covs

            for (i in seq_along(termsAll)) {
                modelTerms <- termsAll[[i]]
                modelCovs <- intersect(modelTerms, covs)

                group <- groups[[i]]$dataSummary

                if (length(modelCovs) == 1) {
                    mahalNote <- jmvcore::Notice$new(
                        options = self$options,
                        name = 'warningMessage',
                        type = jmvcore::NoticeType$WARNING
                    )
                    mahalNote$setContent(
                        .("Mahalanobis distance can only be calculated for models with two or more covariates")
                    )

                    group$insert(1, mahalNote)
                    group$mahal$setVisible(FALSE)
                } else {
                    group$mahal$setVisible(TRUE)
                }
            }
        },
        .initCollinearityTable = function() {
            groups <- self$results$models
            termsAll <- private$.getModelTerms()

            for (i in seq_along(termsAll)) {

                table <- groups$get(key=i)$assump$collin
                terms <- termsAll[[i]]

                if (length(terms) < 1)
                    terms <- ''

                for (i in seq_along(terms))
                    table$addRow(rowKey=i, values=list(term = jmvcore::stringifyTerm(terms[i])))
            }
        },
        .initResPlots=function() {
            groups <- self$results$models
            termsAll <- private$.getModelTerms()
            covs <- self$options$covs

            for (i in seq_along(termsAll)) {
                modelTerms <- termsAll[[i]]
                if (length(modelTerms) < 1) {
                    terms <- ''
                } else {
                    terms <- c(.('Fitted'), self$options$dep)
                    for (term in modelTerms) {
                        if (length(term) == 1 && term %in% covs)
                            terms <- c(terms, term)
                    }
                }

                images <- groups$get(key=i)$assump$resPlots

                for (term in terms)
                    images$addItem(term)
            }
        },
        .initEmm = function() {
            groups <- self$results$models
            termsAll <- private$.getModelTerms()
            emMeans <- self$options$emMeans

            for (i in seq_along(termsAll)) {
                group <- groups$get(key=i)$emm
                terms <- unique(unlist(termsAll[[i]]))

                for (j in seq_along(emMeans)) {
                    emm <- emMeans[[j]]

                    if ( ! is.null(emm) && all(emm %in% terms)) {
                        group$addItem(key=j)
                        emmGroup <- group$get(key=j)
                        emmGroup$setTitle(jmvcore::stringifyTerm(emm))

                        image <- emmGroup$emmPlot
                        size <- private$.plotSize(emm)
                        image$setSize(size[1], size[2])
                    }
                }
            }
        },
        .initEmmTable = function() {
            groups <- self$results$models
            termsAll <- private$.getModelTerms()
            emMeans <- self$options$emMeans
            factors <- self$options$factors

            emMeansTableTitle <- .('Estimated Marginal Means - {term}')
            ciWidthTitle <- jmvcore::format(.('{ciWidth}% Confidence Interval'), ciWidth=self$options$ciWidthEmm)

            for (i in seq_along(termsAll)) {
                group <- groups$get(key=i)$emm
                terms <- unique(unlist(termsAll[[i]]))

                for (j in seq_along(emMeans)) {
                    emm <- emMeans[[j]]

                    if ( ! is.null(emm) && all(emm %in% terms)) {
                        emmGroup <- group$get(key=j)

                        table <- emmGroup$emmTable
                        table$setTitle(jmvcore::format(emMeansTableTitle, term=jmvcore::stringifyTerm(emm)))

                        nLevels <- numeric(length(emm))
                        for (k in rev(seq_along(emm))) {
                            if (emm[k] %in% factors) {
                                table$addColumn(name=emm[k], title=emm[k], type='text', combineBelow=TRUE)
                                nLevels[k] <- length(levels(self$data[[ emm[k] ]]))
                            } else {
                                table$addColumn(name=emm[k], title=emm[k], type='number', combineBelow=TRUE)
                                nLevels[k] <- 3
                            }
                        }

                        table$addColumn(name='emmean', title=.('Marginal Mean'), type='number')
                        table$addColumn(name='se', title=.('SE'), type='number')
                        table$addColumn(name='lower', title=.('Lower'), type='number', superTitle=ciWidthTitle, visible="(ciEmm)")
                        table$addColumn(name='upper', title=.('Upper'), type='number', superTitle=ciWidthTitle, visible="(ciEmm)")

                        nRows <- prod(nLevels)

                        for (k in 1:nRows) {
                            row <- list()
                            table$addRow(rowKey=k, row)
                        }
                    }
                }
            }
        },
        .initOutputs = function() {
            description = function(varType, modelNo=NULL) {
                return(
                    jmvcore::format(
                        .("{varType} of linear regression model{modelNo}"),
                        varType=varType,
                        modelNo=ifelse(is.null(modelNo), "", paste0(" ", modelNo))
                    )
                )
            }

            title = function(varType, modelNo) {
                return(jmvcore::format("{} {}", varType, modelNo))
            }

            if (self$nModels > 1) {
                keys <- seq_len(self$nModels)
                measureTypes <- rep('continuous', self$nModels)

                titles <- vapply(keys, function(key) title(.('Residuals'), key), '')
                descriptions <- vapply(keys, function(key) description(.('Residuals'), key), '')
                self$results$residsOV$set(keys, titles, descriptions, measureTypes)

                titles <- vapply(keys, function(key) title(.('Predicted values'), key), '')
                descriptions <- vapply(keys, function(key) description(.('Predicted values'), key), '')
                self$results$predictOV$set(keys, titles, descriptions, measureTypes)

                titles <- vapply(keys, function(key) title(.("Cook's distance"), key), '')
                descriptions <- vapply(keys, function(key) description(.("Cook's distance"), key), '')
                self$results$cooksOV$set(keys, titles, descriptions, measureTypes)

                titles <- c(
                    vapply(keys, function(key) title(.("Mahalanobis ChiSq"), key), ''),
                    vapply(keys, function(key) title(.("Mahalanobis p"),     key), '')
                )
                descriptions <- c(
                    vapply(keys, function(key) description(.("Mahalanobis distance (Chi sq.)"), key), ''),
                    vapply(keys, function(key) description(.("Mahalanobis distance (p)"),       key), '')
                )
                order <- as.vector(
                    vapply(seq(self$nModels), function(p) c(p, p + self$nModels), numeric(2))
                )
                self$results$mahalOV$set(
                    seq(2 * self$nModels), titles[order], descriptions[order], rep(measureTypes, 2)
                )
            } else {
                titles <- c(.("Mahalanobis ChiSq"), .("Mahalanobis p"))
                descriptions <- c(
                    description(.("Mahalanobis distance (Chi sq.)")),
                    description(.("Mahalanobis distance (p)"))
                )
                self$results$mahalOV$set(seq(2), titles, descriptions, rep('continuous', 2))
            }
        },

        #### Populate tables functions ----
        .populateModelFitTable = function() {
            table <- self$results$modelFit

            models <- self$models
            AIC <- self$AIC
            BIC <- self$BIC

            for (i in seq_along(models)) {
                row <- list()
                row[["aic"]] <- AIC[[i]]
                row[["bic"]] <- BIC[[i]]
                row[["r"]] <- sqrt(summary(models[[i]])$r.squared)
                row[["r2"]] <- summary(models[[i]])$r.squared
                row[["r2Adj"]] <- summary(models[[i]])$adj.r.squared
                row[["rmse"]] <- sqrt(mean(models[[i]]$residuals^2))

                F <- summary(models[[i]])$fstatistic

                if ( ! is.null(F)) {
                    row[["f"]] <- as.numeric(F[1])
                    row[["df1"]] <- as.numeric(F[2])
                    row[["df2"]] <- as.numeric(F[3])
                    row[["p"]] <- stats::pf(F[1], F[2], F[3], lower.tail=FALSE)
                } else {
                    row[["f"]] <- ""
                    row[["df1"]] <- ""
                    row[["df2"]] <- ""
                    row[["p"]] <- ""
                }

                table$setRow(rowNo=i, values = row)
            }

            table$setNote(
                "n",
                jmvcore::format(
                    "Models estimated using sample size of N={n}",
                    n=length(models[[1]]$fitted.values)
                )
            )
        },
        .populateModelCompTable = function() {
            table <- self$results$modelComp

            models <- self$models
            ANOVA <- self$anovaModelComparison
            r <- ANOVA[-1,]

            if (length(models) <= 1)
                return()

            for (i in 1:(length(models)-1)) {
                row <- list()
                row[["r2"]] <- abs(summary(models[[i]])$r.squared - summary(models[[i+1]])$r.squared)
                row[["f"]] <- (r[i,4] / r[i,3]) / (r[i,2] / r[i,1])
                row[["df1"]] <- r[i,3]
                row[["df2"]] <- r[i,1]
                row[["p"]] <- stats::pf(row[["f"]], row[["df1"]], row[["df2"]], lower.tail=FALSE)

                table$setRow(rowNo=i, values = row)
            }
        },
        .populateAnovaTables = function() {
            groups <- self$results$models
            termsAll <- private$.getModelTerms()
            anova <- self$anovaModelTerms

            for (i in seq_along(termsAll)) {

                table <- groups$get(key=i)$anova

                terms <- termsAll[[i]]
                termsB64 <- lapply(terms, jmvcore::toB64)
                r <- anova[[i]]
                rowTerms <- jmvcore::decomposeTerms(rownames(r))
                resIndex <- length(rowTerms)

                for (j in seq_along(terms)) {

                    term <- termsB64[[j]]

                    # check which rows have the same length + same terms
                    index <- which(length(term) == sapply(rowTerms, length) &
                                       sapply(rowTerms, function(x) all(term %in% x)))

                    ss <- r[index, 1]
                    df <- r[index, 2]

                    if (df == 0) {
                        ms <- NaN
                        F <- NaN
                        p <- NaN
                    } else {
                        ms <- ss / df
                        F  <- r[index, 3]
                        p  <- r[index, 4]
                    }

                    row <- list(ss=ss, df=df, ms=ms, F=F, p=p)

                    table$setRow(rowKey=paste0(terms[[j]]), values = row)
                }

                ss <- r[resIndex, 1]
                df <- r[resIndex, 2]
                ms <- ss / df

                row <- list(ss=ss, df=df, ms=ms, F='', p='')

                table$setRow(rowKey='.RES', values = row)

                if ( private$.getIsAliased()[[i]] )
                    table$setNote("alias", .(SINGULAR_WARNING))
            }
        },
        .populateCoefTables = function() {
            groups <- self$results$models
            termsAll <- private$.getRowNamesModel()
            models <- self$models
            modelsScaled <- self$modelsScaled

            for (i in seq_along(termsAll)) {
                table <- groups$get(key=i)$coef

                model <- summary(models[[i]])
                modelScaled <- summary(modelsScaled[[i]])

                CI <- self$CICoefEst[[i]]
                CIScaled <- self$CICoefEstScaled[[i]]
                coef<- model$coef
                coefScaled <- modelScaled$coef
                terms <- termsAll[[i]]
                rowTerms <- jmvcore::decomposeTerms(rownames(coef))

                for (j in seq_along(terms)) {
                    term <- terms[[j]]

                    # check which rows have the same length + same terms
                    index <- which(length(term) == sapply(rowTerms, length) &
                                       sapply(rowTerms, function(x) all(term %in% x)))

                    if (length(index) > 0) {
                        row <- list()
                        row[["est"]] <- coef[index, 1]
                        row[["se"]] <- coef[index, 2]
                        row[["t"]] <- coef[index, 3]
                        row[["p"]] <- coef[index, 4]
                        row[["lower"]] <- CI[index, 1]
                        row[["upper"]] <- CI[index, 2]

                        if (rowTerms[index] == "(Intercept)") {
                            row[["stdEst"]] <- ""
                            row[["stdEstLower"]] <- ""
                            row[["stdEstUpper"]] <- ""
                        } else {
                            row[["stdEst"]] <- coefScaled[index, 1]
                            row[["stdEstLower"]] <- CIScaled[index, 1]
                            row[["stdEstUpper"]] <- CIScaled[index, 2]
                        }
                    } else {
                        row <- list(
                            est = NaN,
                            se = NaN,
                            t = NaN,
                            p = NaN,
                            lower = NaN,
                            upper = NaN,
                            stdEst = NaN,
                            stdEstLower = NaN,
                            stdEstUpper = NaN
                        )
                    }

                    table$setRow(rowKey=jmvcore::composeTerm(term), values=row)

                    if ( private$.getIsAliased()[[i]] )
                        table$setNote("alias", .(SINGULAR_WARNING))
                }
            }
        },
        .populateCooksTable = function() {
            if (! self$options$cooks)
                return()

            groups <- self$results$models
            termsAll <- private$.getModelTerms()
            cooks <- private$.computeCooksSummary()

            for (i in seq_along(termsAll)) {
                table <- groups$get(key=i)$dataSummary$cooks
                table$setRow(rowNo=1, values=cooks[[i]])
            }
        },
        .populateMahalTable = function() {
            if ( ! self$options$mahal)
                return()

            mahal <- self$mahal
            mahalp <- as.numeric(self$options$mahalp)
            groups <- self$results$models
            termsAll <- private$.getModelTerms()

            for (i in seq_along(termsAll)) {
                table <- groups$get(key=i)$dataSummary$mahal
                if (i <= length(mahal) && ! is.null(mahal[[i]])) {
                    mahalRow <- rownames(mahal[[i]][mahal[[i]]$p <= mahalp, ])
                    mahalVal = list(
                        mean = mean(mahal[[i]]$chisq),
                        median = median(mahal[[i]]$chisq),
                        sd = sd(mahal[[i]]$chisq),
                        min = min(mahal[[i]]$chisq),
                        max = max(mahal[[i]]$chisq),
                        excRow = paste(mahalRow, collapse = ", ")
                    )
                    table$setRow(rowNo = 1, values = mahalVal)
                    if (length(mahalRow) == 0) {
                        table$setNote(
                            "excRow",
                            jmvcore::format(
                                .("There were no Mahalanobis distances with a threshold below p < {p}."),
                                p = mahalp)
                            )
                    }
                }
            }
        },
        .populateDurbinWatsonTable = function() {
            if (! self$options$durbin)
                return()

            groups <- self$results$models
            termsAll <- private$.getModelTerms()

            for (i in seq_along(termsAll)) {
                table <- groups$get(key=i)$assump$durbin
                dwTest <- self$durbinWatson[[i]]

                row <- list()
                row[["autoCor"]] <- as.numeric(dwTest[1])
                row[["dw"]] <- as.numeric(dwTest[2])
                row[["p"]] <- as.numeric(dwTest[3])

                table$setRow(rowNo=1, values=row)
            }
        },
        .populateCollinearityTable = function() {
            groups <- self$results$models
            termsAll <- private$.getModelTerms()

            for (i in seq_along(termsAll)) {
                table <- groups$get(key=i)$assump$collin
                terms <- lapply(termsAll[[i]], jmvcore::toB64)

                if (length(self$VIF) == 0)
                    VIF <- NULL
                else
                    VIF <- self$VIF[[i]]

                if (length(dim(VIF)) > 1) {
                    names <- rownames(VIF)
                    VIF <- VIF[,3]
                    names(VIF) <- names
                }

                rowTerms <- jmvcore::decomposeTerms(names(VIF))

                for (j in seq_along(terms)) {
                    row <- list()

                    if (private$.getIsAliased()[[i]]) {
                        row[["tol"]] <- NaN
                        row[["vif"]] <- NaN
                    } else if (length(terms) <= 1) {
                        row[["tol"]] <- 1
                        row[["vif"]] <- 1
                    } else {
                        # check which rows have the same length + same terms
                        index <- which(length(terms[[j]]) == sapply(rowTerms, length) &
                                           sapply(rowTerms, function(x) all(terms[[j]] %in% x)))

                        row[["tol"]] <- 1 / as.numeric(VIF[index])
                        row[["vif"]] <- as.numeric(VIF[index])
                    }

                    table$setRow(rowNo=j, values=row)
                }

                if ( private$.getIsAliased()[[i]] )
                    table$setNote("alias", .(SINGULAR_WARNING))
            }
        },
        .populateEmmTables = function() {
            if (! self$options$emmTables)
                return()

            groups <- self$results$models
            termsAll <- private$.getModelTerms()
            emMeans <- self$options$emMeans
            factors <- self$options$factors
            covs <- self$options$covs

            for (i in seq_along(termsAll)) {
                group <- groups$get(key=i)$emm
                terms <- unique(unlist(termsAll[[i]]))

                for (j in seq_along(emMeans)) {
                    emm <- emMeans[[j]]

                    if ( ! is.null(emm) && all(emm %in% terms)) {
                        emmGroup <- group$get(key=j)
                        table <- emmGroup$emmTable

                        emmTable <- self$emMeans[[i]][[j]]

                        covValues <- list()
                        for (k in seq_along(emm)) {
                            if (emm[k] %in% covs)
                                covValues[[ emm[k] ]] <- sort(unique(emmTable[, jmvcore::toB64(emm[k])]))
                        }

                        for (k in 1:nrow(emmTable)) {
                            row <- list()
                            sign <- list()

                            for (l in seq_along(emm)) {
                                value <- emmTable[k, jmvcore::toB64(emm[l])]

                                if (emm[l] %in% factors) {
                                    row[[emm[l]]] <- jmvcore::fromB64(value)
                                } else {
                                    row[[emm[l]]] <- value

                                    if (value == covValues[[ emm[l] ]][1])
                                        sign[[ emm[l] ]] <- '\u207B'
                                    else if (value == covValues[[ emm[l] ]][3])
                                        sign[[ emm[l] ]] <- '\u207A'
                                    else
                                        sign[[ emm[l] ]] <- '<sup>\u03BC</sup>'
                                }
                            }

                            row[['emmean']] <- emmTable[k, 'emmean']
                            row[['se']] <- emmTable[k, 'SE']
                            row[['lower']] <- emmTable[k, 'lower.CL']
                            row[['upper']] <- emmTable[k, 'upper.CL']

                            table$setRow(rowNo=k, values=row)

                            if (length(covValues) > 0) {
                                table$setNote("sub", .("\u207B mean - 1SD, <sup>\u03BC</sup> mean, \u207A mean + 1SD"))

                                for (l in seq_along(emm)) {
                                    if (emm[l] %in% covs)
                                        table$addSymbol(rowNo=k, emm[l], sign[[ emm[l] ]])
                                }
                            }
                        }
                    }
                }
            }
        },
        .populateNormality = function() {
            groups <- self$results$models
            termsAll <- private$.getModelTerms()

            for (i in seq_along(termsAll)) {
                residuals <- self$residuals[[i]]
                table <- groups$get(key=i)$assump$get('norm')

                res <- try(shapiro.test(residuals), silent=TRUE)
                if (jmvcore::isError(res)) {
                    values <- list(`s[sw]`=NaN, `p[sw]`='')
                } else {
                    values <- list(`s[sw]`=res$statistic, `p[sw]`=res$p.value)
                }
                table$setRow(rowNo=1, values)
            }
        },
        .populateOutputs = function() {
            if (self$options$residsOV && self$results$residsOV$isNotFilled()) {
                self$results$residsOV$setRowNums(private$.getDataRowNums())
                for (i in seq_along(self$residuals))
                    self$results$residsOV$setValues(index=i, self$residuals[[i]])
            }

            if (self$options$predictOV && self$results$predictOV$isNotFilled()) {
                self$results$predictOV$setRowNums(names(self$predicted[[1]]))
                for (i in seq_along(self$fitted))
                    self$results$predictOV$setValues(index=i, as.numeric(self$predicted[[i]]))
            }

            if (self$options$cooksOV && self$results$cooksOV$isNotFilled()) {
                self$results$cooksOV$setRowNums(private$.getDataRowNums())
                for (i in seq_along(self$cooks))
                    self$results$cooksOV$setValues(index=i, self$cooks[[i]])
            }

            if (self$options$mahalOV && self$results$mahalOV$isNotFilled()) {
                self$results$mahalOV$setRowNums(private$.getDataRowNums())
                for (i in seq_along(self$mahal)) {
                    if (is.null(self$mahal[[i]])) next
                    self$results$mahalOV$setValues(index=(i * 2 - 1), self$mahal[[i]]$chisq)
                    self$results$mahalOV$setValues(index=(i * 2 - 0), self$mahal[[i]]$p)
                }
            }
        },

        #### Plot functions ----
        .prepareQQPlot = function() {
            groups <- self$results$models
            termsAll <- private$.getModelTerms()

            for (i in seq_along(termsAll)) {
                image <- groups$get(key=i)$assump$get('qqPlot')
                image$setState(i)
            }
        },
        .qqPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state))
                return(FALSE)

            df <- as.data.frame(qqnorm(scale(self$residuals[[image$state]]), plot.it=FALSE))

            p <- ggplot(data=df, aes(x=x, y=y)) +
                      geom_abline(slope=1, intercept=0, colour=theme$color[1]) +
                      geom_point(aes(x=x,y=y), size=2, colour=theme$color[1]) +
                      xlab(.("Theoretical Quantiles")) +
                      ylab(.("Standardized Residuals")) +
                      ggtheme

            return(p)
        },
        .prepareResPlots = function() {
            groups <- self$results$models
            termsAll <- private$.getModelTerms()

            for (i in seq_along(termsAll)) {
                images <- groups$get(key=i)$assump$resPlots
                for (term in images$itemKeys) {
                    image <- images$get(key=term)
                    image$setState(list(modelNo=i, term=term))
                }
            }
        },
        .resPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state))
                return(FALSE)

            res <- self$residuals[[image$state$modelNo]]

            if (image$state$term == .('Fitted')) {
                x <- self$fitted[[image$state$modelNo]]
            } else {
                x <- self$dataProcessed[[jmvcore::toB64(image$state$term)]]
            }

            df <- data.frame(y=res, x=x)

            p <- ggplot(data=df, aes(y=y, x=x)) +
                      geom_point(aes(x=x,y=y), colour=theme$color[1]) +
                      xlab(image$state$term) +
                      ylab(.("Residuals")) +
                      ggtheme

            return(p)
        },
        .prepareEmmPlots = function() {
            if (! self$options$emmPlots)
                return()

            covs <- self$options$covs
            dep <- self$options$dep

            groups <- self$results$models
            termsAll <- private$.getModelTerms()
            emMeans <- self$options$emMeans
            emMeansTables <- private$.getEmMeansForPlot()

            for (i in seq_along(termsAll)) {
                group <- groups$get(key=i)$emm
                terms <- unique(unlist(termsAll[[i]]))

                for (j in seq_along(emMeans)) {
                    term <- emMeans[[j]]

                    if ( ! is.null(term) && all(term %in% terms)) {
                        image <- group$get(key=j)$emmPlot

                        termB64 <- jmvcore::toB64(term)
                        cont <- FALSE
                        if (term[1] %in% covs)
                            cont <- TRUE

                        d <- emMeansTables[[i]][[j]]

                        if (self$options$piEmm) {
                            sigma <- summary(self$models[[i]])$sigma
                            resDf <- if ('df' %in% names(d))
                                d[['df']]
                            else
                                stats::df.residual(self$models[[i]])
                            alpha <- 1 - self$options$ciWidthEmm / 100
                            tCrit <- stats::qt(1 - alpha / 2, resDf)
                            halfWidth <- tCrit * sqrt(d[['SE']]^2 + sigma^2)
                            d[['lower.PI']] <- d[['emmean']] - halfWidth
                            d[['upper.PI']] <- d[['emmean']] + halfWidth
                        }

                        for (k in 1:3) {
                            if ( ! is.na(termB64[k])) {
                                if (term[k] %in% covs) {
                                    if (k > 1) {
                                        d[[ termB64[k] ]] <- factor(d[[ termB64[k] ]])
                                        levels(d[[ termB64[k] ]]) <- c('-1SD', 'Mean', '+1SD')
                                    }
                                } else {
                                    d[[ termB64[k] ]] <- factor(
                                        jmvcore::fromB64(d[[ termB64[k] ]]),
                                        jmvcore::fromB64(levels(d[[ termB64[k] ]]))
                                    )
                                }
                            }
                        }

                        names <- list(
                            'x'=termB64[1], 'y'='emmean', 'lines'=termB64[2], 'plots'=termB64[3],
                            'lower'='lower.CL', 'upper'='upper.CL'
                        )
                        names <- lapply(names, function(x) if (is.na(x)) NULL else x)

                        labels <- list('x'=term[1], 'y'=dep, 'lines'=term[2], 'plots'=term[3])
                        labels <- lapply(labels, function(x) if (is.na(x)) NULL else x)

                        image$setState(list(
                            data=d, names=names, labels=labels, cont=cont,
                            modelNo=i, term=term, termB64=termB64,
                            otherVars=setdiff(terms, term)
                        ))
                    }
                }
            }
        },
        .emmPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state))
                return(FALSE)

            data <- image$state$data
            names <- image$state$names
            labels <- image$state$labels
            cont <- image$state$cont
            modelNo <- image$state$modelNo
            term <- image$state$term
            termB64 <- image$state$termB64
            otherVars <- image$state$otherVars
            covs <- self$options$covs

            dodge <- position_dodge(0.4)

            if (theme$bw) {
                lty <- names$lines
                color <- NULL
            } else {
                lty <- NULL
                color <- names$lines
            }

            captions <- character(0)

            ptData <- NULL
            if (self$options$emmPoints != 'none') {
                groupIsCov <- ( ! is.na(term[2]) && term[2] %in% covs)
                facetIsCov <- ( ! is.na(term[3]) && term[3] %in% covs)

                if (groupIsCov || facetIsCov) {
                    suppressedVar <- if (groupIsCov) term[2] else term[3]
                    captions <- c(captions, jmvcore::format(
                        .("Observed data not shown: '{var}' is a covariate here, binned to -1SD/Mean/+1SD -- no observation has a real position on it."),
                        var=suppressedVar))
                } else {
                    ptData <- private$.emmPointData(modelNo, term, termB64, names)
                }
            }

            if (self$options$piEmm && length(otherVars) > 0) {
                captions <- c(captions, jmvcore::format(
                    .("Prediction interval for a new observation at this reference point ({vars} held at their mean, or averaged across levels)."),
                    vars=listItems(self, otherVars)))
            }

            p <- ggplot(data=data,
                        aes_string(x=names$x, y=names$y, group=names$lines),
                        inherit.aes = FALSE)

            if (cont) {
                if (self$options$piEmm && is.null(names$plots) && is.null(names$lines))
                    p <- p +
                        geom_ribbon(
                            aes_string(x=names$x, ymin='lower.PI', ymax='upper.PI'),
                            show.legend=FALSE, alpha=.2, fill='#F48FB1'
                        )

                if (self$options$ciEmm && is.null(names$plots) && is.null(names$lines))
                    p <- p +
                        geom_ribbon(
                            aes_string(x=names$x, ymin=names$lower, ymax=names$upper),
                            show.legend=TRUE, alpha=.3
                        )

                if ( ! is.null(ptData))
                    p <- p + geom_point(
                        data=ptData, mapping=aes_string(x=names$x, y='.y', color=color, fill=color),
                        inherit.aes=FALSE, alpha=.4
                    )

                p <- p + geom_line(aes_string(color=color, fill=color, linetype=lty))
            } else {
                if (self$options$piEmm)
                    p <- p +
                        geom_errorbar(
                            aes_string(
                                x=names$x, ymin='lower.PI', ymax='upper.PI',
                                color=names$lines, fill=names$lines
                            ),
                            width=.05, size=.5, alpha=.4, position=dodge, show.legend=FALSE
                        )

                if ( ! is.null(ptData))
                    p <- p + geom_point(
                        data=ptData, mapping=aes_string(x=names$x, y='.y', color=names$lines, fill=names$lines),
                        inherit.aes=FALSE, alpha=.4,
                        position=position_jitterdodge(jitter.width=.1, dodge.width=.4)
                    )

                p <- p +
                    geom_point(
                        aes_string(color=names$lines, fill=names$lines),
                        position = dodge
                    )

                if (self$options$ciEmm)
                    p <- p +
                        geom_errorbar(
                            aes_string(
                                x=names$x, ymin=names$lower, ymax=names$upper,
                                color=names$lines, fill=names$lines
                            ),
                            width=.1, size=.8, position=dodge
                        )
            }

            if ( ! is.null(names$plots)) {
                formula <- as.formula(paste(". ~", names$plots))
                p <- p + facet_grid(formula)
            }

            p <- p +
                labs(x=labels$x, y=labels$y, fill=labels$lines, color=labels$lines, linetype=labels$lines) +
                ggtheme + theme(panel.spacing = unit(2, "lines"))

            if (length(captions) > 0)
                p <- p + labs(caption=paste(captions, collapse='\n'))

            return(p)
        },

        # observed points to overlay on an emm plot for `term` (up to 3
        # variable names, NA-padded): x = the observation's own value for
        # term[1]; y = partial residuals (this term's own fitted
        # contribution + the observation's residual, jmvcore/car's usual
        # component+residual definition) or the raw dependent variable,
        # per `emmPoints`. Colour/facet columns reuse the B64 names already
        # used by the emm summary table `data`, so the same aes_string()
        # mappings apply to both layers. Returns NULL if this term has no
        # own contribution to the model (e.g. only appears in a higher-order
        # interaction with a variable outside `term`).
        .emmPointData = function(modelNo, term, termB64, names) {
            dp <- self$dataProcessed
            depB64 <- jmvcore::toB64(self$options$dep)
            covs <- self$options$covs

            x64 <- termB64[1]
            if (term[1] %in% covs) {
                x <- dp[[x64]]
            } else {
                col <- dp[[x64]]
                x <- factor(jmvcore::fromB64(as.character(col)), levels=jmvcore::fromB64(levels(col)))
            }

            if (self$options$emmPoints == 'raw') {
                y <- dp[[depB64]]
            } else {
                fit <- private$.partialFitForTerm(modelNo, termB64[ ! is.na(termB64)])
                if (is.null(fit))
                    return(NULL)
                y <- fit + self$residuals[[modelNo]]
            }

            df <- data.frame(x=x, .y=y)
            names(df)[1] <- names$x

            if ( ! is.na(termB64[2])) {
                col <- dp[[ termB64[2] ]]
                df[[ names$lines ]] <- factor(jmvcore::fromB64(as.character(col)), levels=jmvcore::fromB64(levels(col)))
            }
            if ( ! is.na(termB64[3])) {
                col <- dp[[ termB64[3] ]]
                df[[ names$plots ]] <- factor(jmvcore::fromB64(as.character(col)), levels=jmvcore::fromB64(levels(col)))
            }

            return(df)
        },

        # partial fit for a set of variables (B64 names, no NAs): the sum of
        # every model term whose own variables are entirely contained in
        # `termB64`, evaluated at each observation's own data (car/termplot's
        # "component" -- add residuals for the "component + residual", i.e.
        # partial residual). With a single-predictor model this equals the
        # model's full fitted value, so partial residuals reduce to the raw
        # dependent variable exactly.
        .partialFitForTerm = function(modelNo, termB64) {
            model <- self$models[[modelNo]]
            tt <- stats::terms(model)
            labels <- attr(tt, 'term.labels')
            if (length(labels) == 0)
                return(NULL)

            factorsMat <- attr(tt, 'factors')
            varNames <- rownames(factorsMat)

            covered <- vapply(labels, function(lbl) {
                vars <- varNames[factorsMat[, lbl] > 0]
                all(vars %in% termB64)
            }, logical(1))

            if ( ! any(covered))
                return(NULL)

            # predict(..., type='terms') mean-centres every term column
            # (attr(,'constant') is what was subtracted out, summed across
            # ALL terms); summing only the covered subset and adding that
            # constant back gives a prediction anchored on the response's
            # own scale, using each observation's own value for the terms in
            # `termB64` and every other term's mean contribution -- i.e.
            # exactly the reference-grid definition this needs, matching a
            # continuous covariate held at its mean exactly, and a factor at
            # its frequency-weighted average (a mild mismatch from emmeans'
            # default equal-weighting, only when 'termB64' excludes an
            # unbalanced factor -- flagged in task_plan.md, not chased for v1).
            ptTerms <- stats::predict(model, type='terms')
            as.numeric(rowSums(ptTerms[, labels[covered], drop=FALSE]) + attr(ptTerms, 'constant'))
        },

        #### Descriptive plots ----
        # Model-free views of the analysed variables (dep vs. each covariate,
        # dep by each factor, a scatterplot matrix), independent of the model
        # blocks -- so they work even before a block has been built. No
        # image$state: each render function reads self$data/self$options
        # directly (requiresData: true in the .r.yaml), the same
        # state-free pattern corrInspect uses for its own scatterplots.
        .initPlots = function() {
            scatterArray <- self$results$plots$scatter
            for (cov in self$options$covs)
                scatterArray$addItem(key=cov)

            boxArray <- self$results$plots$box
            for (factor in self$options$factors)
                boxArray$addItem(key=factor)

            vars <- c(self$options$covs, self$options$dep)
            if (length(vars) >= 2) {
                side <- max(450, 150 * length(vars))
                self$results$plots$matrix$setSize(side, side)
            }
        },
        .descScatterPlot = function(image, ggtheme, theme, ...) {
            if ( ! isTRUE(self$options$scatDep) || is.null(self$options$dep))
                return(FALSE)

            xvar <- image$key
            if (length(xvar) == 0)
                return(FALSE)

            private$.drawDescScatter(xvar, self$options$dep, ggtheme, theme)
        },
        .descBoxPlot = function(image, ggtheme, theme, ...) {
            if ( ! isTRUE(self$options$boxDep) || is.null(self$options$dep))
                return(FALSE)

            fvar <- image$key
            if (length(fvar) == 0)
                return(FALSE)

            private$.drawDescBox(fvar, self$options$dep, ggtheme, theme)
        },
        .descMatrixPlot = function(image, ggtheme, theme, ...) {
            if ( ! isTRUE(self$options$scatMatrix) || is.null(self$options$dep))
                return(FALSE)

            # dep last: the lower-triangle layout puts the higher-indexed
            # variable of each pair on y, so this keeps the dependent
            # variable on y whenever it appears, instead of on x against
            # every covariate.
            vars <- c(self$options$covs, self$options$dep)
            if (length(vars) < 2)
                return(FALSE)

            private$.drawDescMatrix(vars, ggtheme, theme)
        },

        # an annotated scatter (fitted line, r/p/equation, CI/PI bands),
        # adapted from corrInspect's .drawAnnotatedScatter with an optional
        # colour-by-group factor and, with groupLines, one fitted line/band
        # per group instead of a single overall one. The r/p/equation
        # annotation always summarises the OVERALL relationship (ignoring
        # grouping), since a single label can't stand in for several groups.
        .drawDescScatter = function(xvar, yvar, ggtheme, theme) {
            data <- self$data
            x <- jmvcore::toNumeric(data[[xvar]])
            y <- jmvcore::toNumeric(data[[yvar]])
            groupVar <- self$options$groupVar
            g <- if ( ! is.null(groupVar)) data[[groupVar]] else NULL

            ok <- stats::complete.cases(x, y)
            if ( ! is.null(g)) ok <- ok & ! is.na(g)
            x <- x[ok]
            y <- y[ok]
            if ( ! is.null(g)) g <- droplevels(g[ok])

            if (length(x) < 3 || stats::sd(x) == 0 || stats::sd(y) == 0)
                return(FALSE)

            df <- data.frame(x=x, y=y)
            if ( ! is.null(g)) df$g <- g

            textColour <- ggplot2::calc_element('text', themeOnly(ggtheme))$colour
            pointColour <- '#5DADE2'
            lineColour <- 'black'
            ciLevel <- self$options$ciWidth / 100
            perGroup <- ! is.null(g) && isTRUE(self$options$groupLines)

            p <- if (is.null(g))
                ggplot2::ggplot(df, ggplot2::aes(x=x, y=y))
            else
                ggplot2::ggplot(df, ggplot2::aes(x=x, y=y, colour=g))

            if (isTRUE(self$options$descLine)) {
                fitSets <- if (perGroup) split(df, df$g) else list(all=df)
                for (sub in fitSets) {
                    if (nrow(sub) < 3 || stats::sd(sub$x) == 0)
                        next

                    fit <- stats::lm(y ~ x, data=sub)
                    grid <- data.frame(x=seq(min(sub$x), max(sub$x), length.out=100))
                    bandAlpha <- if (perGroup) .2 else .35

                    if (isTRUE(self$options$descPIBand)) {
                        pi <- stats::predict(fit, newdata=grid, interval='prediction', level=ciLevel)
                        ribbon <- cbind(grid, as.data.frame(pi[, c('lwr', 'upr'), drop=FALSE]))
                        p <- p + ggplot2::geom_ribbon(
                            data=ribbon, mapping=ggplot2::aes(x=x, ymin=lwr, ymax=upr),
                            inherit.aes=FALSE, fill='#F48FB1', alpha=bandAlpha)
                    }
                    if (isTRUE(self$options$descCIBand)) {
                        ci <- stats::predict(fit, newdata=grid, interval='confidence', level=ciLevel)
                        ribbon <- cbind(grid, as.data.frame(ci[, c('lwr', 'upr'), drop=FALSE]))
                        p <- p + ggplot2::geom_ribbon(
                            data=ribbon, mapping=ggplot2::aes(x=x, ymin=lwr, ymax=upr),
                            inherit.aes=FALSE, fill=pointColour, alpha=bandAlpha)
                    }
                }
            }

            p <- p + ggplot2::geom_point(alpha=.7)

            if (isTRUE(self$options$descLine)) {
                if (perGroup)
                    p <- p + ggplot2::geom_smooth(method='lm', formula=y ~ x, se=FALSE)
                else
                    p <- p + ggplot2::geom_smooth(method='lm', formula=y ~ x, se=FALSE, colour=lineColour)
            }

            if (isTRUE(self$options$descAnnot)) {
                fit <- pearsonFit(x, y)
                coefv <- stats::coef(stats::lm(y ~ x, data=df))
                lines <- c(
                    sprintf('r = %.3f', fit$r), pText(fit$p),
                    sprintf('%s = %.3g + %.3g·%s', yvar, coefv[1], coefv[2], xvar))
                label <- paste(lines, collapse='\n')
                xRange <- range(x)
                xPad <- diff(xRange) * .02
                p <- p + ggplot2::annotate('text', x=xRange[1] + xPad, y=Inf, hjust=0, vjust=1.3,
                                            label=label, colour=textColour, size=3.6)
            }

            yRange <- range(y)
            ySpan <- diff(yRange)
            topPad <- if (isTRUE(self$options$descAnnot)) .15 else .05
            p <- p + ggplot2::coord_cartesian(ylim=c(yRange[1] - ySpan * .05, yRange[2] + ySpan * topPad))

            p <- p +
                ggplot2::labs(x=xvar, y=yvar, colour=groupVar) +
                ggtheme +
                ggplot2::theme(
                    axis.text=ggplot2::element_text(colour='black'),
                    axis.title=ggplot2::element_text(colour='black'),
                    axis.ticks=ggplot2::element_line(colour='black'))

            print(p)
            TRUE
        },

        # the dependent variable by one factor: box or violin, with optional
        # jittered points.
        .drawDescBox = function(fvar, yvar, ggtheme, theme) {
            data <- self$data
            y <- jmvcore::toNumeric(data[[yvar]])
            f <- data[[fvar]]

            ok <- stats::complete.cases(y) & ! is.na(f)
            y <- y[ok]
            f <- droplevels(f[ok])

            if (length(y) < 2 || nlevels(f) < 1)
                return(FALSE)

            df <- data.frame(f=f, y=y)
            pointColour <- '#5DADE2'

            p <- ggplot2::ggplot(df, ggplot2::aes(x=f, y=y))

            if (identical(self$options$boxType, 'violin'))
                p <- p + ggplot2::geom_violin(fill=pointColour, alpha=.3, colour='black')
            else
                p <- p + ggplot2::geom_boxplot(fill=pointColour, alpha=.3, colour='black', outlier.shape=NA)

            if (isTRUE(self$options$boxPoints))
                p <- p + ggplot2::geom_jitter(width=.15, height=0, alpha=.5, colour=pointColour)

            p <- p +
                ggplot2::labs(x=fvar, y=yvar) +
                ggtheme +
                ggplot2::theme(
                    axis.text=ggplot2::element_text(colour='black'),
                    axis.title=ggplot2::element_text(colour='black'),
                    axis.ticks=ggplot2::element_line(colour='black'))

            print(p)
            TRUE
        },

        # scatterplot matrix for the dependent variable + covariates,
        # adapted from corrInspect's .drawScatterMatrix (lower triangle =
        # mini scatter + lm line + r, diagonal = variable name, upper
        # triangle blank), with points optionally coloured by groupVar.
        .drawDescMatrix = function(vars, ggtheme, theme) {
            data <- self$data
            n <- length(vars)
            textColour <- ggplot2::calc_element('text', themeOnly(ggtheme))$colour
            pointColour <- '#5DADE2'
            lineColour <- 'black'
            groupVar <- self$options$groupVar
            g <- if ( ! is.null(groupVar)) data[[groupVar]] else NULL

            idx <- function(row, col) (row - 1) * n + col
            panels <- vector('list', n * n)

            for (row in seq_len(n)) {
                for (col in seq_len(n)) {
                    if (row == col) {
                        panels[[idx(row, col)]] <- grid::textGrob(
                            vars[row], gp=grid::gpar(col=textColour, fontsize=10))
                    } else if (row > col) {
                        x <- jmvcore::toNumeric(data[[ vars[col] ]])
                        y <- jmvcore::toNumeric(data[[ vars[row] ]])
                        ok <- stats::complete.cases(x, y)
                        if ( ! is.null(g)) ok <- ok & ! is.na(g)
                        x <- x[ok]
                        y <- y[ok]
                        gp <- if ( ! is.null(g)) g[ok] else NULL

                        if (length(x) < 3 || stats::sd(x) == 0 || stats::sd(y) == 0) {
                            panels[[idx(row, col)]] <- grid::nullGrob()
                        } else {
                            r <- stats::cor(x, y)
                            df <- data.frame(x=x, y=y)
                            if ( ! is.null(gp)) df$g <- gp
                            xPanelRange <- range(x)
                            xPanelPad <- diff(xPanelRange) * 0.05

                            pnl <- if (is.null(gp))
                                ggplot2::ggplot(df, ggplot2::aes(x=x, y=y)) +
                                    ggplot2::geom_point(colour=pointColour, alpha=0.6, size=0.8)
                            else
                                ggplot2::ggplot(df, ggplot2::aes(x=x, y=y, colour=g)) +
                                    ggplot2::geom_point(alpha=0.6, size=0.8)

                            pnl <- pnl +
                                ggplot2::geom_smooth(method='lm', formula=y ~ x, se=FALSE,
                                                      colour=lineColour, linewidth=0.5) +
                                ggplot2::annotate('text', x=xPanelRange[1] + xPanelPad, y=Inf,
                                                   hjust=0, vjust=1.3,
                                                   label=sprintf('r = %.2f', r),
                                                   colour=textColour, size=2.6) +
                                # ggtheme first for its discrete colour/fill
                                # scales (so a grouped panel matches the
                                # descriptive scatters' group colours);
                                # theme_minimal() only replaces the theme
                                # part, the scales from ggtheme stick.
                                ggtheme +
                                ggplot2::theme_minimal(base_size=6) +
                                ggplot2::theme(
                                    legend.position='none',
                                    axis.text=ggplot2::element_text(colour='black', size=5),
                                    axis.title=ggplot2::element_blank(),
                                    axis.ticks=ggplot2::element_line(colour='black'),
                                    panel.grid=ggplot2::element_blank(),
                                    panel.border=ggplot2::element_rect(colour='grey85', fill=NA))

                            panels[[idx(row, col)]] <- ggplot2::ggplotGrob(pnl)
                        }
                    } else {
                        panels[[idx(row, col)]] <- grid::nullGrob()
                    }
                }
            }

            gridExtra::grid.arrange(grobs=panels, nrow=n, ncol=n)
            TRUE
        },

        #### Helper functions ----
        .getModelTerms = function() {
            if (is.null(private$.modelTerms)) {
                blocks <- self$options$blocks

                terms <- list()
                if (is.null(blocks)) {
                    terms[[1]] <- c(self$options$covs, self$options$factors)
                } else {
                    for (i in seq_along(blocks)) {
                        terms[[i]] <- unlist(blocks[1:i], recursive = FALSE)
                    }
                }
                private$.modelTerms <- terms
            }

            return(private$.modelTerms)
        },
        .getRowNamesModel = function() {
            if (is.null(private$.rowNamesModel)) {
                factors <- self$options$factors
                termsAll <- private$.getModelTerms()

                rowNamesAll <- list()
                for (i in seq_along(termsAll)) {
                    rowNames <- list()
                    rowNames[[1]] <- "(Intercept)"

                    terms <- termsAll[[i]]
                    for (term in terms) {
                        if (any(term %in% factors))
                            rowNames <- c(rowNames, private$.contrastModel(term))
                        else
                            rowNames[[length(rowNames) + 1]] <- jmvcore::toB64(term)
                    }
                    rowNamesAll[[i]] <- rowNames
                }
                private$.rowNamesModel <- rowNamesAll
            }

            return(private$.rowNamesModel)
        },
        .getEmMeansForPlot = function() {
            if (is.null(private$.emMeansForPlot))
                private$.emMeansForPlot <- private$.computeEmMeans(forPlot=TRUE)

            return(private$.emMeansForPlot)
        },
        .getDataRowNums = function() {
            if (is.null(private$.dataRowNums))
                private$.dataRowNums <- rownames(self$dataProcessed)

            return(private$.dataRowNums)
        },
        .getIsAliased = function() {
            if (is.null(private$.isAliased)) {
                aliased <- list()
                for (i in seq_along(self$models))
                    aliased[[i]] <- any(summary(self$models[[i]])$aliased)

                private$.isAliased <- aliased
            }

            return(private$.isAliased)
        },
        .contrastModel = function(terms) {
            # Returns the names of the factor contrasts as they are
            # defined by the lm function

            factors <- self$options$factors

            nLevels <- list()
            for (factor in factors)
                nLevels[[factor]] <- length(levels(self$data[[factor]]))

            contrast <- list()
            for (term in terms) {
                if (term %in% factors) {
                    contrast[[term]] <- paste0(jmvcore::toB64(term), 1:(nLevels[[term]] - 1))
                } else {
                    contrast[[term]] <- jmvcore::toB64(term)
                }
            }

            contrastGrid <- expand.grid(contrast)
            contrastTerms <- list()
            for (i in 1:nrow(contrastGrid))
                contrastTerms[[i]] <- as.character(unlist(contrastGrid[i,]))

            return(contrastTerms)
        },
        .contrastsCoefTable = function(terms) {
            # Returns the names of the factor contrasts as they are
            # displayed in the coef table

            factors <- self$options$factors
            refLevels <- self$refLevels
            refVars <- sapply(refLevels, function(x) x$var)

            levels <- list()
            for (factor in factors)
                levels[[factor]] <- levels(self$data[[factor]])

            contrast <- list()
            for (term in terms) {
                if (term %in% factors) {
                    ref <- refLevels[[which(term == refVars)]][['ref']]
                    refNo <- which(ref == levels[[term]])

                    contrLevels <- levels[[term]][-refNo]
                    refLevel <- levels[[term]][refNo]

                    c <- paste(contrLevels, refLevel, sep = ' \u2013 ')
                    if (length(terms) > 1)
                        c <- paste0('(', c, ')')

                    contrast[[term]] <- c
                } else {
                    contrast[[term]] <- term
                }
            }
            contrastGrid <- expand.grid(contrast)
            contrastNames <- apply(contrastGrid, 1, jmvcore::stringifyTerm)

            return(contrastNames)
        },
        .formulas = function() {
            dep <- self$options$dep
            depB64 <- jmvcore::toB64(dep)
            terms <- private$.getModelTerms()

            formulas <- list();
            for (i in seq_along(terms)) {
                termsB64 <- lapply(terms[[i]], jmvcore::toB64)
                composedTerms <- jmvcore::composeTerms(termsB64)
                formulas[[i]] <- as.formula(paste(depB64, paste0(composedTerms, collapse ="+"), sep="~"))
            }

            return(formulas)
        },
        .cleanData = function(naOmit=TRUE, naSkip=NULL) {
            dep <- self$options$dep
            covs <- self$options$covs
            factors <- self$options$factors
            weights <- self$options$weights
            refLevels <- self$refLevels

            dataRaw <- self$data

            data <- list()

            refVars <- sapply(refLevels, function(x) x$var)

            for (factor in factors) {
                if (length(levels(dataRaw[[factor]])) <= 1)
                    stop(jmvcore::format(.("Factor '{factor}' needs to have at least 2 levels"), factor=factor))

                ref <- refLevels[[which(factor == refVars)]][['ref']]
                column <- factor(
                    dataRaw[[factor]],
                    ordered = FALSE,
                    levels = levels(dataRaw[[factor]])
                )
                levels(column) <- jmvcore::toB64(levels(column))
                column <- relevel(column, ref = jmvcore::toB64(ref))

                data[[jmvcore::toB64(factor)]] <- column
                stats::contrasts(data[[jmvcore::toB64(factor)]]) <- private$.createContrasts(levels(column))
            }

            for (cov in c(dep, covs, weights))
                data[[jmvcore::toB64(cov)]] <- jmvcore::toNumeric(dataRaw[[cov]])

            global_weights <- attr(dataRaw, "jmv-weights")
            if (is.null(weights) && ! is.null(global_weights))
                data[[".WEIGHTS"]] = jmvcore::toNumeric(global_weights)

            attr(data, 'row.names') <- rownames(self$data)
            attr(data, 'class') <- 'data.frame'

            if (naOmit) {
                data <- tibble::rownames_to_column(data)
                data <- tidyr::drop_na(data, -tidyselect::all_of(naSkip))
                data <- tibble::column_to_rownames(data)
            }

            return(data)
        },
        .plotSize = function(emm) {
            data <- self$data
            factors <- self$options$factors

            levels <- list()
            for (i in seq_along(emm)) {
                column <- data[[ emm[i] ]]

                if (emm[i] %in% factors) {
                    levels[[ emm[i] ]] <- levels(column)
                } else {
                    if (i == 1)
                        levels[[ emm[i] ]] <- ''
                    else
                        levels[[ emm[i] ]] <- c('-1SD', 'Mean', '+1SD')
                }
            }

            nLevels <- as.numeric(sapply(levels, length))
            nLevels <- ifelse(is.na(nLevels[1:3]), 1, nLevels[1:3])
            nCharLevels <- as.numeric(sapply(lapply(levels, nchar), max))
            nCharLevels <- ifelse(is.na(nCharLevels[1:3]), 0, nCharLevels[1:3])
            nCharNames <- as.numeric(nchar(names(levels)))
            nCharNames <- ifelse(is.na(nCharNames[1:3]), 0, nCharNames[1:3])

            xAxis <- 30 + 20
            yAxis <- 30 + 20

            if (emm[1] %in% factors) {
                width <- max(350, 25 * nLevels[1] * nLevels[2] * nLevels[3])
                height <- 300 + ifelse(nLevels[3] > 1, 20, 0)
            } else {
                width <- max(350, 300 * nLevels[3])
                height <- 300 + ifelse(nLevels[3] > 1, 20, 0)
            }

            legend <- max(25 + 21 + 3.5 + 8.3 * nCharLevels[2] + 28, 25 + 10 * nCharNames[2] + 28)
            width <- yAxis + width + ifelse(nLevels[2] > 1, legend, 0)
            height <- xAxis + height

            return(c(width, height))
        },
        .createContrasts=function(levels) {
            if (self$options$intercept == 'refLevel') {
                contrast <- contr.treatment(levels)
                dimnames(contrast) <- NULL
            } else {
                nLevels <- length(levels)
                dummy <- contr.treatment(levels)
                dimnames(dummy) <- NULL
                coding <- matrix(rep(1/nLevels, prod(dim(dummy))), ncol=nLevels-1)
                contrast <- (dummy - coding)
            }

            return(contrast)
        },
        .scaleData = function(data) {
            for (col in names(data)) {
                if ( ! is.factor(data[[col]]) && length(unique(data[[col]])) > 1 )
                    data[[col]] <- scale(data[[col]])
            }

            return(data)
        })
)
