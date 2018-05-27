# Bootstrap Magnitude Based Predictions (BMBP)
# Copyright 2018 Mladen Jovanovic
#
# To reference plese use:
#
#      Jovanovic, M. (2018, May, 25). Bootstrap Magnitude Based Predictions. 
#          Complementary Training. Retrieved from http://www.complementarytraining.com/bmbp
#
#

require(dplyr)
require(ggplot2)
require(ggthemes)
require(gridExtra)
require(boot)
require(irr)
require(ggridges)
require(reshape2)

# GGplot theme for all the graphs
theme_mladen <- function(...) {
    theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank())
}

# Get bootstraped reliability estimates using 2 repeated measures 
bootReliability <- function(df,
                            var1,
                            var2,
                            xlabel = var1,
                            ylabel = var2,
                            SWCcoef = 0.2,
                            SWC = NA, # Use this if you want to use specific SWC (not est)
                            logTrans = FALSE, # Log transformation - next version
                            measurementError = 0, # Random error to add in boot samples
                            paired = TRUE, # Will be used in next version (now it is paired)
                            confidence = 0.95,
                            samples = 2000,
                            plot = TRUE,
                            BAplotShow = "both", # What variable to show 
                            iter = TRUE,
                            na.rm = TRUE,
                            errorNoise = 10^-6) { # The random noise for boot not to complain
    # Plot    
    df$diff <- df[[var2]] - df[[var1]]
    sd.diff <- sd(df$diff, na.rm = na.rm)
    
    # What is on x axis of BA plot?
    if(BAplotShow == "both") {
        df$avg <- (df[[var1]] + df[[var2]]) / 2
        BAxlabel <- paste("(", xlabel, " + ", ylabel, ") / 2", sep = "")
    } else {
        df$avg <- df[[var1]]
        BAxlabel <- xlabel
    }
    
    t.conf <- qt(1-((1-confidence)/2), df = length(na.omit(df$diff))-1) 
    
    # For SWC use pooled between Variable 1 and Variable 2
    N <- nrow(df)
    
    if(is.na(SWC)) { # If SWC not provided estimate it
        SWC_plot <- sqrt(((var(df[[var1]], na.rm = na.rm)*(N-1)) + (var(df[[var2]], na.rm = na.rm)*(N-1))) /(2*(N-1))) * SWCcoef
    } else {
        SWC_plot <- SWC
    }

    # Create graphs
    graphDF <- data.frame(x = df[[var1]], y = df[[var2]], diff = df$diff, avg = df$avg)
    plot1 <- ggplot(graphDF, aes(x = x, y = y)) +
        theme_mladen() +
        geom_point(alpha = 0.3) +  
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey") +
        geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
        xlab(xlabel) + 
        ylab(ylabel) + 
        theme(legend.position="none")

    plot2 <- ggplot(graphDF, aes(x = avg, y = diff)) + 
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = -SWC_plot, ymax = SWC_plot, fill = "grey", alpha = 0.3) +
        geom_point(alpha = 0.3) +
        theme_mladen() +
        geom_hline(yintercept = mean(df$diff, na.rm = na.rm), color = "grey") +
        geom_hline(yintercept = mean(df$diff, na.rm = na.rm) + t.conf*sd.diff, color = "grey", linetype = "dashed") +
        geom_hline(yintercept = mean(df$diff, na.rm = na.rm) - t.conf*sd.diff, color = "grey", linetype = "dashed") +
        geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) +
        ylab(paste(ylabel, " - ", xlabel, sep = "")) + 
        xlab(BAxlabel) +
        theme(legend.position="right")
    
    BAgraph <- arrangeGrob(plot1, plot2, ncol=2)
    
    # Bootstrap
    bs <- function(data, indices) {
        if(iter) cat(".")
        
        # Data for reliability
        df <- data[indices,] # allows boot to select sample 
        
        # Add random measurement error
        mError1 <- rnorm(n = nrow(df), sd = measurementError, mean = 0)
        mError2 <- rnorm(n = nrow(df), sd = measurementError, mean = 0)
        df[[var1]] <- df[[var1]] + mError1
        df[[var2]] <- df[[var2]] + mError2
        
        # Get the difference 
        df$diff <- df[[var2]] - df[[var1]]
        
        # Get the thresholds for LOA
        t.conf <- qt(1-((1-confidence)/2), df = length(na.omit(df$diff))-1) 
        
        ## Estimations
        
        if(is.na(SWC)) { # If SWC not provided estimate it
            SWC <- sqrt(((var(df[[var1]], na.rm = na.rm)*(N-1)) + (var(df[[var2]], na.rm = na.rm)*(N-1))) /(2*(N-1))) * SWCcoef
        } else { # If provided add some random noise so boot works
            SWC <- SWC + rnorm(n = 1, mean = 0, sd = errorNoise)     
        }
        
        sd.diff <- sd(df$diff, na.rm = na.rm)
        meanDiff <- mean(df$diff, na.rm = na.rm)
        medianDiff <- median(df$diff, na.rm = na.rm) + rnorm(n = 1, mean = 0, sd = errorNoise) # Small 'noise' for boot
        MADdiff <- mad(df$diff, na.rm = na.rm) + rnorm(n = 1, mean = 0, sd = errorNoise) # Small 'noise' for boot
        TE <- sd.diff / sqrt(2)
        LOA <- t.conf*sd.diff
        CORR <- cor(df[[var1]], df[[var2]], use = "complete.obs") # Correlation
        ICC <-  icc(data.frame(df[[var1]], df[[var2]]), 
                    "twoway", "agreement")$value # ICC
        # Get the hit rate 
        hit.rate <- ifelse(df$diff <= SWC & df$diff >= -SWC, 1, 0)
        hit.rate <- sum(hit.rate) / length(hit.rate)
        hit.rate <- hit.rate + rnorm(n = 1, mean = 0, sd = errorNoise)
        # Vector for saving results
        results <- numeric(13) 
        
        results[1] <- meanDiff
        results[2] <- sd.diff
        results[3] <- medianDiff
        results[4] <- MADdiff
        results[5] <- LOA
        results[6] <- TE
        results[7] <- CORR
        results[8] <- ICC
        results[9] <- SWC
        results[10] <- SWC / TE
        results[11] <- SWC / LOA
        results[12] <- SWC / MADdiff
        results[13] <- hit.rate
        
        return(results)
    }
    
    results <- boot(data = df, statistic = bs,
                    R = samples)
    
    meanDiff <- boot.ci(results, type="perc", index = 1, conf = confidence)
    SDdiff <- boot.ci(results, type="perc", index = 2, conf = confidence)
    medianDiff <- boot.ci(results, type="perc", index = 3, conf = confidence)
    MADdiff <- boot.ci(results, type="perc", index = 4, conf = confidence)
    LOA <- boot.ci(results, type="perc", index = 5, conf = confidence)
    TE <- boot.ci(results, type="perc", index = 6, conf = confidence)
    CORR <- boot.ci(results, type="perc", index = 7, conf = confidence)
    ICC <- boot.ci(results, type="perc", index = 8, conf = confidence)
    SWC <- boot.ci(results, type="perc", index = 9, conf = confidence)
    SWCtoTE <- boot.ci(results, type="perc", index = 10, conf = confidence)
    SWCtoLOA <- boot.ci(results, type="perc", index = 11, conf = confidence)
    SWCtoMAD <- boot.ci(results, type="perc", index = 12, conf = confidence)
    HitRate <- boot.ci(results, type="perc", index = 13, conf = confidence)
    
    estimates <- data.frame(meanDiff = c(meanDiff$t0, meanDiff$percent[4], meanDiff$percent[5]),
                            SDdiff = c(SDdiff$t0, SDdiff$percent[4], SDdiff$percent[5]),
                            medianDiff = c(medianDiff$t0, medianDiff$percent[4], medianDiff$percent[5]),
                            MADdiff = c(MADdiff$t0, MADdiff$percent[4], MADdiff$percent[5]),
                            LOA = c(LOA$t0, LOA$percent[4], LOA$percent[5]),
                            TE = c(TE$t0, TE$percent[4], TE$percent[5]),
                            CORR = c(CORR$t0, CORR$percent[4], CORR$percent[5]),
                            ICC = c(ICC$t0, ICC$percent[4], ICC$percent[5]),
                            SWC = c(SWC$t0, SWC$percent[4], SWC$percent[5]),
                            SWCtoTE = c(SWCtoTE$t0, SWCtoTE$percent[4], SWCtoTE$percent[5]),
                            SWCtoLOA = c(SWCtoLOA$t0, SWCtoLOA$percent[4], SWCtoLOA$percent[5]),
                            SWCtoMAD = c(SWCtoMAD$t0, SWCtoMAD$percent[4], SWCtoMAD$percent[5]),
                            HitRate = c(HitRate$t0, HitRate$percent[4], HitRate$percent[5]))
    
    estimates <- t(estimates) # Transpose
    row.names(estimates) <- NULL
    metrics <- c("meanDiff",
                 "SDdiff",
                 "medianDiff",
                 "MADdiff",
                 "LOA",
                 "TE",
                 "CORR",
                 "ICC",
                 "SWC",
                 "SWCtoTE",
                 "SWCtoLOA",
                 "SWCtoMAD",
                 "HitRate")
    estimates <- data.frame(Estimate = factor(metrics, levels = metrics),
                            estimates)
    colnames(estimates) <- c("Estimate", "Value", "Lower", "Upper")
    
    # Create graphs
    signalVSnoiseGraph <- ggplot(filter(estimates, Estimate %in% c("SWCtoTE", "SWCtoLOA", "SWCtoMAD", "HitRate")),
                                 aes(x = Value, y = Estimate)) +
        theme_mladen() +
        geom_errorbarh(aes(xmax = Upper, xmin = Lower), color = "black", height = 0.3) + 
        geom_point(shape=21, size=1.5, fill="white") +
        ylab("") + xlab("")
    
    if (plot == TRUE) {        
        plot(BAgraph)
        plot(signalVSnoiseGraph)
    }
    
    return(list(dfRows = nrow(df),
                Var1N = sum(!is.na(df[[var1]])),
                Var2N = sum(!is.na(df[[var2]])),
                DiffN = sum(!is.na(df$diff)),
                
                bootSamples = samples,
                measurementError = measurementError,
                paired = paired,
                confidence = confidence,
                SWCcoef = SWCcoef,
                logTrans = logTrans,
                errorNoise = errorNoise,
                estimates = estimates,
                graphs = list(BAgraph = BAgraph,
                              signalVSnoiseGraph = signalVSnoiseGraph)))
}


# Get bootstraped reliability estimates
bootReliabilityMultiple <- function(df,
                                    SWCcoef = 0.2,
                                    logTrans = FALSE, # Log transformation - next version
                                    measurementError = 0, # Random error to add in boot samples
                                    paired = TRUE, # Will be used in next version (now it is paired)
                                    type = "pairwise", # Use c("pairwise", "sequential")
                                    confidence = 0.95,
                                    samples = 2000,
                                    plot = TRUE,
                                    iter = TRUE,
                                    na.rm = TRUE,
                                    errorNoise = 10^-6) {
    
    
    # Define functions
    # Function for pairwise difference
    pairwiseDifference <- function(df, type = "pairwise") {
        x <- as.matrix(df)
        
        # create column combination pairs
        if (type == "pairwise") {
            prs <- cbind(rep(1:ncol(x), each = ncol(x)), 1:ncol(x))
            col.diffs <- prs[prs[, 1] < prs[, 2], , drop = FALSE]
        } else { # "sequential"
            col.diffs <- matrix(nrow = ncol(df) - 1, ncol = 2)
            col.diffs[, 1] <- seq(1, ncol(df) - 1)
            col.diffs[, 2] <- seq(2, ncol(df))
        }
        # do pairwise differences 
        result <- x[, col.diffs[, 2]] - x[, col.diffs[, 1], drop = FALSE]
        
        # set colnames
        if(is.null(colnames(x)))
            colnames(x) <- 1:ncol(x)
        
        colnames(result) <- paste(colnames(x)[col.diffs[, 2]], ".VS.", 
                                  colnames(x)[col.diffs[, 1]], sep = "")
        return(data.frame(result))
    }
    
    # Function for pairwise correlation
    pairwiseCorrelation <- function(df, type = "pairwise") {
        x <- as.matrix(df)
        
        if (type == "pairwise") {
            prs <- cbind(rep(1:ncol(x), each = ncol(x)), 1:ncol(x))
            col.cor <- prs[prs[, 1] < prs[, 2], , drop = FALSE]
        } else { # "sequential"
            col.cor <- matrix(nrow = ncol(df) - 1, ncol = 2)
            col.cor[, 1] <- seq(1, ncol(df) - 1)
            col.cor[, 2] <- seq(2, ncol(df))
        }
        
        pair.cor <- matrix(nrow = 1, ncol = nrow(col.cor))
        
        # do pairwise correlation 
        for( i in seq(1, nrow(col.cor))) {
            pair.cor[1, i] <- cor(x[,col.cor[i, 1]], x[,col.cor[i, 2]], use = "complete.obs")
        }
        
        # set colnames
        if(is.null(colnames(x)))
            colnames(x) <- 1:ncol(x)
        
        colnames(pair.cor) <- paste(colnames(x)[col.cor[, 2]], ".VS.", 
                                    colnames(x)[col.cor[, 1]], sep = "")
        return(data.frame(pair.cor))
    }
    
    # Function for pairwise pooled SD
    pairwiseSD <- function(df, type = "pairwise", na.rm = TRUE) {
        x <- as.matrix(df)
        
        # create column combination pairs
        if (type == "pairwise") {
            prs <- cbind(rep(1:ncol(x), each = ncol(x)), 1:ncol(x))
            col.comp <- prs[prs[, 1] < prs[, 2], , drop = FALSE]
        } else { # "sequential"
            col.comp <- matrix(nrow = ncol(df) - 1, ncol = 2)
            col.comp[, 1] <- seq(1, ncol(df) - 1)
            col.comp[, 2] <- seq(2, ncol(df))
        }
        
        pair.sd <- matrix(nrow = 1, ncol = nrow(col.comp))
        
        # do pairwise SD 
        for( i in seq(1, nrow(col.comp))) {
            var1 <- x[,col.comp[i, 1]]
            N1 <- sum(!is.na(var1)) - 1
            
            var2 <- x[,col.comp[i, 2]]
            N2 <- sum(!is.na(var2)) -1 
            
            pair.sd[1, i] <- sqrt(((var(var1, na.rm = na.rm) * N1) + (var(var2, na.rm = na.rm) * N2)) / (N1+N2)) 
        }
        
        # set colnames
        if(is.null(colnames(x)))
            colnames(x) <- 1:ncol(x)
        
        colnames(pair.sd) <- paste(colnames(x)[col.comp[, 2]], ".VS.", 
                                   colnames(x)[col.comp[, 1]], sep = "")
        return(data.frame(pair.sd))
    }
    
    # Find pairwise hit rate
    pairwiseHitRate <- function(pairDiff, SWC, na.rm = TRUE) {
        HitRate <- numeric(ncol(pairDiff))
        
        for (i in seq(1, ncol(pairDiff))) {
            HitRate[i] <- sum(ifelse(as.numeric(pairDiff[,i]) <= SWC[[i]] & as.numeric(pairDiff[,i]) >= -SWC[[i]], 1, 0), na.rm = na.rm) / sum(!is.na(pairDiff[,i]))
        }
        return(HitRate)   
    }
    
    
    # Get summary data
    dfN <- nrow(df)
    varName <- colnames(df)
    varN <- length(varName)
    observationsN <- df %>% summarize_all(function(x)sum(!is.na(x)))
    
    differences <- pairwiseDifference(df, type = type)
    diffName <- colnames(differences)
    diffPairs <- length(diffName)
    diffN <- differences %>% summarize_all(function(x)sum(!is.na(x)))
    
    ## Boot functions
    bs <- function(data, indices) {
        if(iter) cat(".")
        
        # Data for reliability
        df <- data[indices,] # allows boot to select sample 
        
        # Add measurement error
        df <- df %>%
            mutate_all(function(x)x + rnorm(n = nrow(df), sd = measurementError, mean = 0))
        
        # get pairwise differences
        differences <- pairwiseDifference(df, type = type)
        
        # Convert both to long
        dfLong <- melt(df, id.vars = NULL)
        differencesLong <- melt(differences, id.vars = NULL)
        
        # Summary stats for variables
        varSummary <- dfLong %>%
            group_by(variable) %>%
            summarize(mean = mean(value, na.rm = na.rm),
                      SD = sd(value, na.rm = na.rm),
                      median = median(value, na.rm = na.rm) + rnorm(n = 1, mean = 0, sd = errorNoise),
                      MAD = mad(value, na.rm = na.rm) + rnorm(n = 1, mean = 0, sd = errorNoise),
                      SWC = sd(value, na.rm = na.rm) * SWCcoef)
        
        # Save SWC for later
        varSD <- varSummary$SD
        
        # rotate summary
        tmp <- data.frame(varSummary)
        rownames(tmp) <- tmp$variable
        tmp$variable <- NULL
        varSummary <- t(tmp)
        
        # Summary stats for variables
        diffSummary <- differencesLong %>%
            group_by(variable) %>%
            summarize(meanDiff = mean(value, na.rm = na.rm),
                      SDdiff = sd(value, na.rm = na.rm),
                      medianDiff = median(value, na.rm = na.rm) + rnorm(n = 1, mean = 0, sd = errorNoise),
                      MADdiff = mad(value, na.rm = na.rm) + rnorm(n = 1, mean = 0, sd = errorNoise),
                      LOA = SDdiff * qt(1-((1-confidence)/2), df = sum(!is.na(value))-1), 
                      TE = SDdiff / sqrt(2))
        
        # save meanDiff for later
        meanDiff <- diffSummary$meanDiff
        
        # Rotate the summary
        tmp <- data.frame(diffSummary)
        rownames(tmp) <- tmp$variable
        tmp$variable <- NULL
        diffSummary <- t(tmp)
        
        # Get pairwise SWC
        pairSWC <- pairwiseSD(df, type = type, na.rm = na.rm) * SWCcoef
        rownames(pairSWC) <- "SWC"
        
        # Get pairwise correlation
        pairCor <- pairwiseCorrelation(df, type = type)
        rownames(pairCor) <- "CORR"
        
        
        # Combine with pairSWC and pairCor
        diffSummary <- rbind(diffSummary, pairSWC, pairCor)
        
        # Convert back to data frame
        diffSummary <- t(diffSummary)
        
        # Convert back to data.frame
        diffSummary <- data.frame(pair = rownames(diffSummary),
                                  diffSummary)
        rownames(diffSummary) <- NULL
        
        # Create new variables
        diffSummary <- diffSummary %>%
            mutate(SWCtoTE = SWC / TE,
                   SWCtoLOA = SWC / LOA,
                   SWCtoMAD = SWC / MADdiff)
        
        # HitRate
        HitRate <- pairwiseHitRate(differences, pairSWC, na.rm = na.rm)
        HitRate <- HitRate + rnorm(n = length(HitRate), mean = 0, sd = errorNoise)
        
        ## Combine everything in summary df
        diffSummary$HitRate <- HitRate
        # Save TE for later
        pairTE <- diffSummary$TE
        # combine to data.frame
        tmp <- data.frame(diffSummary)
        rownames(tmp) <- tmp$pair
        tmp$pair <- NULL
        diffSummary <- t(tmp)
        
        # Create summary SWC
        SWC <- sqrt(sum(varSD^2 * (as.numeric(observationsN)-1)) / (sum(as.numeric(observationsN)-1))) * SWCcoef
        
        # Create summary TE
        TE <- sqrt(sum(pairTE^2 * (as.numeric(diffN)-1)) / (sum(as.numeric(diffN)-1)))
        
        # Get mean absolute mean difference [bias]
        avgAbsMeanDiff <- mean(abs(meanDiff), na.rm = na.rm)
        
        # Combine to numeric
        varSummary <- as.numeric(as.matrix(varSummary))
        diffSummary <- as.numeric(as.matrix(diffSummary))
        
        # ICC
        ICC <- icc(df,"twoway", "agreement")$value
        
        # Avg HitRate
        avgHitRate <- mean(HitRate) 
        
        # Save everything in one vector for boot
        results <- c(avgAbsMeanDiff, SWC, TE, SWC / TE, ICC,  avgHitRate, varSummary, diffSummary)
        
        return(results)
    }
    
    # Bootrstap
    results <- boot(data = df, statistic = bs,
                    R = samples)
    
    
    # Get estimates and confidene intervals
    estimates <- matrix(ncol = 3, nrow = ncol(results$t))
    colnames(estimates) <- c("Value", "Lower", "Upper")
    
    for (i in seq(1, ncol(results$t))) {
        bootRes <- boot.ci(results, type="perc", index = i, conf = confidence)
        estimates[i, 1] <- bootRes$t0
        estimates[i, 2] <- bootRes$percent[4]
        estimates[i, 3] <- bootRes$percent[5]
    }
    
    # Overall Results
    overallSummary <- data.frame(Estimate = factor(c("avgAbsMeanDiff", "SWC", "TE", "SWCtoTE", "ICC", "avgHitRate"),
                                                   levels = c("avgAbsMeanDiff", "SWC", "TE", "SWCtoTE", "ICC", "avgHitRate")),
                                 estimates[1:6,])
    # Variables summaries
    indexStart <- 7
    indexLength <- varN * 5
    indexStop <- indexStart + indexLength 
    
    variableSummary <- expand.grid(Estimate = c("mean", "SD", "median", "MAD", "SWC"),
                                   Variable = varName)
    variableSummary$Value <- NA
    variableSummary$Lower <- NA
    variableSummary$Upper <- NA
    
    for (i in seq(1, indexLength) ) {
        variableSummary$Value[i] <- estimates[indexStart + i - 1, 1]  
        variableSummary$Lower[i] <- estimates[indexStart + i - 1, 2] 
        variableSummary$Upper[i] <- estimates[indexStart + i - 1, 3] 
    }
    
    
    # Pairwise summaries
    indexStart <- indexStop
    indexLength <- diffPairs * 12
    indexStop <- indexStart + indexLength
    
    pairwiseSummary <- expand.grid(Estimate = c("meanDiff", "SDdiff", "medianDiff", "MADdiff",
                                                "LOA", "TE", "SWC", "CORR", "SWCtoTE",
                                                "SWCtoLOA", "SWCtoMAD", "HitRate"),
                                   Variable = diffName)
    pairwiseSummary$Value <- NA
    pairwiseSummary$Lower <- NA
    pairwiseSummary$Upper <- NA
    
    for (i in seq(1, indexLength) ) {
        pairwiseSummary$Value[i] <- estimates[indexStart + i - 1, 1]  
        pairwiseSummary$Lower[i] <- estimates[indexStart + i - 1, 2] 
        pairwiseSummary$Upper[i] <- estimates[indexStart + i - 1, 3] 
    }
    
    # Plot
    overallGraph <- ggplot(overallSummary,
                           aes(x = Value, y = 3)) +
        theme_mladen() +
        geom_errorbarh(aes(xmax = Upper, xmin = Lower), color = "black", height = 0.3) + 
        geom_point(shape=21, size=2, fill="white") +
        facet_wrap(~Estimate, scales = "free_x") + 
        ylab("") + xlab("") + ylim(1, 5) +
        theme(axis.title.y=  element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    
    variableGraph <- ggplot(variableSummary, aes(x = Value, y = reorder(Variable, desc(Variable)))) +
        theme_mladen() +
        geom_errorbarh(aes(xmax = Upper, xmin = Lower), color = "black", height = 0.3) + 
        geom_point(shape=21, size=2, fill="white") +
        facet_wrap(~Estimate, scales = "free_x") + 
        ylab("") + xlab("")
    
    pairwiseGraph <- ggplot(pairwiseSummary, aes(x = Value, y = reorder(Variable, desc(Variable)))) +
        theme_mladen() +
        geom_errorbarh(aes(xmax = Upper, xmin = Lower), color = "black", height = 0.3) + 
        geom_point(shape=21, size=2, fill="white") +
        facet_wrap(~Estimate, scales = "free_x") + 
        ylab("") + xlab("")
    
    if(plot) {
        plot(overallGraph)
        plot(variableGraph)
        plot(pairwiseGraph)
        
    }
    
    
    return(list(dfRows = dfN,
                numberOfVariables = varN,
                variableNames = varName,
                variableObservations = observationsN,
                
                pairwiseComparisons = diffPairs,
                pairwiseNames = diffName,
                pairwiseObservations = diffN,
                pairwiseDF = differences,
                
                bootSamples = samples,
                measurementError = measurementError,
                paired = paired,
                type = type,
                confidence = confidence,
                SWCcoef = SWCcoef,
                logTrans = logTrans,
                errorNoise = errorNoise,
                overallSummary = overallSummary,
                variableSummary = variableSummary,
                pairwiseSummary = pairwiseSummary,
                graphs = list(overallGraph = overallGraph,
                              variableGraph = variableGraph,
                              pairwiseGraph = pairwiseGraph)))
    
}

# Get bootstraped estimates of effects
bootChange <- function(df,
                       pre,
                       post,
                       xlabel = pre,
                       ylabel = post,
                       logTrans = FALSE, # Log transformation - next version
                       SWCcoef = 0.2,
                       SWC = NA, # Use this if you want to use specific SWC (not est)
                       confidence = 0.95,
                       samples = 2000,
                       TE = 0, # Typical Error
                       paired = TRUE, # Will be used in next version (now it is paired)
                       iter = TRUE,
                       na.rm = TRUE,
                       plot = TRUE,
                       errorNoise = 10^-6) {
    
    
    # Bootstrap
    bs <- function(data, indices) {
        if(iter) cat(".")
        
        # Get the data
        df <- data[indices,] # allows boot to select sample 
        
        # Add Typical Error around data pre and post scores
        TypicalError1 <- rnorm(n = nrow(df), sd = TE, mean = 0)
        TypicalError2 <- rnorm(n = nrow(df), sd = TE, mean = 0)
        
        df[[pre]] <- df[[pre]] + TypicalError1
        df[[post]] <- df[[post]] + TypicalError2
        
        # Get the difference 
        df$diff <- df[[post]] - df[[pre]]
        
        ## Calculate the estimates
        
        if(is.na(SWC)) { # If SWC not provided estimate it
            # Estimate SWC from Pre
            SWC <- sd(df[[pre]], na.rm = na.rm) * SWCcoef
        } else { # If provided add some random noise so boot works
            SWC <- SWC + rnorm(n = 1, mean = 0, sd = errorNoise)     
        }
        
        N <- nrow(df)
        meanDiff <- mean(df$diff, na.rm = na.rm)
        # Get the MBIs for mean change
        sd.diff <- sd(df$diff, na.rm = na.rm)
        
        medianDiff <- median(df$diff, na.rm = na.rm) + rnorm(n = 1, mean = 0, sd = errorNoise) # Small 'noise' for boot
        MADdiff <- mad(df$diff, na.rm = na.rm) + rnorm(n = 1, mean = 0, sd = errorNoise) # Small 'noise' for boot
        
        # Estimate Cohen's D using pooled SD
        pooledSD <- sqrt(((var(df[[pre]], na.rm = na.rm)*(N-1)) + (var(df[[post]], na.rm = na.rm)*(N-1))) /(2*(N-1)))
        CohenD <- (mean(df[[post]], na.rm = na.rm) - mean(df[[pre]], na.rm = na.rm)) / pooledSD 
        
        # 'Novel metric' - robust efect
        robustEffect <- medianDiff / MADdiff 
        
        # Magnitude based predictions
        MBP.Harmful.Count <- sum(ifelse(df$diff < -SWC, 1, 0), na.rm = na.rm) / N
        MBP.Beneficial.Count <- sum(ifelse(df$diff > SWC, 1, 0), na.rm = na.rm) / N
        MBP.Trivial.Count <- sum(ifelse(df$diff >= -SWC & df$diff <= SWC, 1, 0), na.rm = na.rm) / N
        
        # Get the mean harmful and beneficial effect
        MBP.Harmful.Mean <- mean(df$diff[df$diff < -SWC], na.rm = na.rm)
        MBP.Beneficial.Mean <- mean(df$diff[df$diff > SWC], na.rm = na.rm)
        
        # Get the median harmful and beneficial effect
        MBP.Harmful.Median <- median(df$diff[df$diff < -SWC], na.rm = na.rm)
        MBP.Beneficial.Median <- median(df$diff[df$diff > SWC], na.rm = na.rm)
        
        #### The "TRICKY PART" - In the case of no 'data' in harmful/beneficil
        #### the boot will complain due NAs or NaNs
        #### The solution would be to either remove NAs from the boot sample
        #### Or to return threshold value (plus some random noise for boot)
        if(is.na(MBP.Harmful.Mean)) MBP.Harmful.Mean <- -SWC + rnorm(n = 1, mean = 0, sd = errorNoise)
        if(is.na(MBP.Harmful.Median)) MBP.Harmful.Median <- -SWC + rnorm(n = 1, mean = 0, sd = errorNoise)
        if(is.na(MBP.Beneficial.Mean)) MBP.Beneficial.Mean <- SWC + rnorm(n = 1, mean = 0, sd = errorNoise)
        if(is.na(MBP.Beneficial.Median)) MBP.Beneficial.Median <- SWC + rnorm(n = 1, mean = 0, sd = errorNoise)
        #####
        #####
        #####
        
        
        # Vector for saving results
        results <- numeric(10)     
        
        results[1] <- meanDiff
        results[2] <- sd.diff
        results[3] <- medianDiff
        results[4] <- MADdiff
        results[5] <- SWC
        results[6] <- CohenD
        results[7] <- robustEffect
        results[8] <- MBP.Harmful.Count
        results[9] <- MBP.Trivial.Count
        results[10] <- MBP.Beneficial.Count
        results[11] <- MBP.Harmful.Mean
        results[12] <- MBP.Beneficial.Mean 
        results[13] <- MBP.Harmful.Median
        results[14] <- MBP.Beneficial.Median
        
        return(results)
    }
    
    results <- boot(data = df, statistic = bs,
                    R = samples)
    
    meanDiff <- boot.ci(results, type="perc", index = 1, conf = confidence)
    
    # Get MBI for mean
    SWC <- boot.ci(results, type="perc", index = 5, conf = confidence)
    SWCpre <- SWC$t0
    
    meanDiff.Harmful <- sum(ifelse(results$t[,1] < -SWCpre, 1, 0)) / samples
    meanDiff.Trivial <- sum(ifelse(results$t[,1] >= -SWCpre & results$t[,1] <= SWCpre, 1, 0)) / samples
    meanDiff.Beneficial <- sum(ifelse(results$t[,1] > SWCpre, 1, 0)) / samples
    
    SDdiff <- boot.ci(results, type="perc", index = 2, conf = confidence)
    
    # Get MBI for median
    medianDiff <- boot.ci(results, type="perc", index = 3, conf = confidence)
    medianDiff.Harmful <- sum(ifelse(results$t[,3] < -SWCpre, 1, 0)) / samples
    medianDiff.Trivial <- sum(ifelse(results$t[,3] >= -SWCpre & results$t[,3] <= SWCpre, 1, 0)) / samples
    medianDiff.Beneficial <- sum(ifelse(results$t[,3] > SWCpre, 1, 0)) / samples
    
    MADdiff <- boot.ci(results, type="perc", index = 4, conf = confidence)

    
    # Get MBI for Cohen's D
    CohenD <- boot.ci(results, type="perc", index = 6, conf = confidence)
    CohenD.Harmful <- sum(ifelse(results$t[,6] < -SWCcoef, 1, 0)) / samples
    CohenD.Trivial <- sum(ifelse(results$t[,6] >= -SWCcoef & results$t[,6] <= SWCcoef, 1, 0)) / samples
    CohenD.Beneficial <- sum(ifelse(results$t[,6] > SWCcoef, 1, 0)) / samples
    
    robustEffect <- boot.ci(results, type="perc", index = 7, conf = confidence)
    
    # MBP
    MBP.Harmful.Count <- boot.ci(results, type="perc", index = 8, conf = confidence)
    MBP.Trivial.Count <- boot.ci(results, type="perc", index = 9, conf = confidence)
    MBP.Beneficial.Count <- boot.ci(results, type="perc", index = 10, conf = confidence)
    MBP.Harmful.Mean <- boot.ci(results, type="perc", index = 11, conf = confidence)
    MBP.Beneficial.Mean <- boot.ci(results, type="perc", index = 12, conf = confidence) 
    MBP.Harmful.Median <- boot.ci(results, type="perc", index = 13, conf = confidence)
    MBP.Beneficial.Median <- boot.ci(results, type="perc", index = 14, conf = confidence)  
    
    changeEstimate <- data.frame(meanDiff = c(meanDiff$t0, meanDiff$percent[4], meanDiff$percent[5]),
                                 meanDiff.Harmful = c(meanDiff.Harmful, NA, NA),
                                 meanDiff.Trivial = c(meanDiff.Trivial, NA, NA),
                                 meanDiff.Beneficial = c(meanDiff.Beneficial, NA, NA),
                                 
                                 SDdiff = c(SDdiff$t0, SDdiff$percent[4], SDdiff$percent[5]),
                                 
                                 medianDiff = c(medianDiff$t0, medianDiff$percent[4], medianDiff$percent[5]),
                                 medianDiff.Harmful = c(medianDiff.Harmful, NA, NA),
                                 medianDiff.Trivial = c(medianDiff.Trivial, NA, NA),
                                 medianDiff.Beneficial = c(medianDiff.Beneficial, NA, NA),
                                 
                                 MADdiff = c(MADdiff$t0, MADdiff$percent[4], MADdiff$percent[5]),
                                 SWC = c(SWC$t0, SWC$percent[4], SWC$percent[5]),
                                 
                                 CohenD = c(CohenD$t0, CohenD$percent[4], CohenD$percent[5]),
                                 CohenD.Harmful = c(CohenD.Harmful, NA, NA),
                                 CohenD.Trivial = c(CohenD.Trivial, NA, NA),
                                 CohenD.Beneficial = c(CohenD.Beneficial, NA, NA),
                                 
                                 robustEffect = c(robustEffect$t0, robustEffect$percent[4], robustEffect$percent[5]),
                                 
                                 MBP.Harmful.Count = c(MBP.Harmful.Count$t0, MBP.Harmful.Count$percent[4], MBP.Harmful.Count$percent[5]),
                                 MBP.Trivial.Count = c(MBP.Trivial.Count$t0, MBP.Trivial.Count$percent[4], MBP.Trivial.Count$percent[5]),
                                 MBP.Beneficial.Count = c(MBP.Beneficial.Count$t0, MBP.Beneficial.Count$percent[4], MBP.Beneficial.Count$percent[5]),
                                 MBP.Harmful.Mean = c(MBP.Harmful.Mean$t0, MBP.Harmful.Mean$percent[4], MBP.Harmful.Mean$percent[5]),
                                 MBP.Beneficial.Mean = c(MBP.Beneficial.Mean$t0, MBP.Beneficial.Mean$percent[4], MBP.Beneficial.Mean$percent[5]), 
                                 MBP.Harmful.Median = c(MBP.Harmful.Median$t0, MBP.Harmful.Median$percent[4], MBP.Harmful.Median$percent[5]),
                                 MBP.Beneficial.Median = c(MBP.Beneficial.Median$t0, MBP.Beneficial.Median$percent[4], MBP.Beneficial.Median$percent[5]))
    changeEstimate <- t(changeEstimate)
    changeEstimate <- data.frame(Estimate = factor(row.names(changeEstimate),
                                                   levels = row.names(changeEstimate)),
                                 changeEstimate)
    colnames(changeEstimate) <- c("Estimate", "Value", "Lower", "Upper")
    row.names(changeEstimate) <- NULL
    
    # Create graphs
    plotDF <- changeEstimate %>%
        filter(Estimate %in% c("MBP.Harmful.Count",
                               "MBP.Trivial.Count",
                               "MBP.Beneficial.Count"))
    plotDF$Magnitude <- factor(c("Harmful", "Trivial", "Beneficial"),
                               levels = c("Harmful", "Trivial", "Beneficial"))
    plotDF$Effect <- changeEstimate[13:15, 2]
    
    gg <- ggplot(plotDF, aes(y = Magnitude, x = Value)) +
        theme_bw() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank()) +
        geom_errorbarh(aes(xmax = Upper, xmin = Lower), color = "black", height = 0.05) + 
        ylab("") + 
        xlab("Probability") +
        geom_point(shape=21, size=3, fill="white")+
        geom_point(aes(x = Effect), color = "grey", shape = "+", size = 7)
    
    effectsGraph <- gg
    if (plot == TRUE) plot(gg)
    
    plotDF <- changeEstimate %>%
        filter(Estimate %in% c("meanDiff",
                               "medianDiff",
                               "MBP.Harmful.Mean",
                               "MBP.Beneficial.Mean",
                               "MBP.Harmful.Median",
                               "MBP.Beneficial.Median"))
    
    plotDF$Change <- factor(c("Change Mean",
                              "Change Median",
                              "Harmful Mean",
                              "Beneficial Mean",
                              "Harmful Median",
                              "Beneficial Median"),
                            levels = c(
                                "Harmful Median",
                                "Harmful Mean",
                                "Beneficial Median",
                                "Beneficial Mean",
                                "Change Median",
                                "Change Mean"))
    
    gg <- ggplot(plotDF, aes(y = Change, x = Value)) +
        theme_bw() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank()) +
        geom_vline(xintercept = 0, color = "grey") +
        geom_vline(xintercept = -SWCpre, color = "grey", linetype = "dashed") +
        geom_vline(xintercept = SWCpre, color = "grey", linetype = "dashed") +
        geom_errorbarh(aes(xmax = Upper, xmin = Lower), color = "black", height = 0.1) + 
        ylab("") + 
        xlab("Change") +
        geom_point(shape=21, size=3, fill="white")
    changeGraph <- gg
    if (plot == TRUE) plot(gg)
    
    df$diff <- df[[post]] - df[[pre]]
    
    return(list(dfRows = nrow(df),
                PreN = sum(!is.na(df[[pre]])),
                Post2N = sum(!is.na(df[[post]])),
                DiffN = sum(!is.na(df$diff)),
                bootSamples = samples,
                TE = TE,
                paired = paired,
                confidence = confidence,
                SWCcoef = SWCcoef,
                logTrans = logTrans,
                errorNoise = errorNoise,
                changeEstimate = changeEstimate,
                graphs = list(effectsGraph = effectsGraph,
                              changeGraph = changeGraph)))
}

