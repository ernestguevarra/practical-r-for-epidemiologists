################################################################################
#
# prfe.r
#
# This file contains the sample R language functions developed in the
# handbook 'Practical R for Epidemiologists'.
#


################################################################################
#
# Two-by-Two table (display table, RR, sample OR, MLE OR)
#
tab2by2 <- function(exposure, outcome)
  {
  tab <- table(exposure, outcome)
  a <- tab[1,1]; b <- tab[1,2]; c <- tab[2,1]; d <- tab[2,2]
  rr <- (a / (a + b)) / (c / (c + d))
  se.log.rr <- sqrt((b / a) / (a + b) + (d / c) / (c + d))
  lci.rr <- exp(log(rr) - 1.96 * se.log.rr)
  uci.rr <- exp(log(rr) + 1.96 * se.log.rr)
  or <- (a / b) / (c / d)
  se.log.or <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
  lci.or <- exp(log(or) - 1.96 * se.log.or)
  uci.or <- exp(log(or) + 1.96 * se.log.or)
  ft <- fisher.test(tab)
  cat("\n")
  print(tab)
  cat("\nRelative Risk              :", rr, 
      "\n95% CI                     :", lci.rr, uci.rr, "\n")
  cat("\nSample Odds Ratio          :", or,
      "\n95% CI                     :", lci.or, uci.or, "\n")
  cat("\nConditional MLE Odds Ratio :", ft$estimate, 
      "\n95% CI                     :", ft$conf.int, "\n\n")
  }


################################################################################
#
# Calculate and display odds ratios and 95% CIs from logistic regression
#
lreg.or <- 	function(model, digits = 2)
  {
  lreg.coeffs <- coef(summary(model))
  OR  <- exp(lreg.coeffs[ ,1])
  LCI <- exp(lreg.coeffs[ ,1] - 1.96 * lreg.coeffs[ ,2])
  UCI <- exp(lreg.coeffs[ ,1] + 1.96 * lreg.coeffs[ ,2])
  p   <- lreg.coeffs[ ,4]
  lreg.or <- round(cbind(OR, LCI, UCI, p), digits = digits)        
  lreg.or
  }


################################################################################
#
# Sample OR from a two-by-two table
#
or <- function(x) 
  {
  (x[1,1] / x[1,2]) / (x[2,1] / x[2,2])
  }


################################################################################
#
# Sample RR from a two-by-two table
#
rr <- function(x)
  {
  (x[1,1] / (x[1,1] + x[1,2])) / (x[2,1] / (x[2,1] + x[2,2]))
  }


################################################################################
#
# Woolf's test for the homogeneity of odds ratios from a 2-by-2-by-k table
#
woolf.test <- function(x) 
  {
  x <- x + 0.5
  k <- dim(x)[3]
  or <- apply(x, 3, function(x) {(x[1, 1] / x[1, 2]) / (x[2, 1] / x[2, 2])})
  w <- apply(x, 3, function(x) {1 / sum(1 / x)})
  chi.sq <- sum(w * (log(or) - weighted.mean(log(or), w))^2)
  p <- pchisq(chi.sq, df = k - 1, lower.tail = FALSE)
  cat("\nWoolf's X2 :", chi.sq, "\np-value    :", p, "\n")
  }


################################################################################
#
# Function to test the null hypothesis that the variance to mean ratio of a
# vector of numbers is equal to one
#
v2m.test <- function(data)
  {
  nsu <- length(data); obs <- sum(data)
  obs.nsu.mean <- obs / nsu; obs.nsu.var <- var(data)
  var.mean.ratio <- obs.nsu.var / obs.nsu.mean
  chi2 <- sum((data - obs.nsu.mean)^2) / obs.nsu.mean
  df <- nsu - 1; p <- 1 - pchisq(chi2, df)
  names(chi2) <- "Chi-square"
  names(df) <- "df"
  names(var.mean.ratio) <- "Variance : mean ratio"
  v2m <- list(method = "Variance to mean test",
              data.name = deparse(substitute(data)),
              statistic = chi2,
              parameter = df,
              p.value = p, 
              estimate = var.mean.ratio
             )
  class(v2m) <- "htest"
  return(v2m)
  }


################################################################################
#
# Example of defining functions and a new class to deal with 2-by-2 tables
#
rr22 <- function(exposure, outcome)
  {
  tab <- table(exposure, outcome)
  a <- tab[1,1]; b <- tab[1,2]; c <- tab[2,1]; d <- tab[2,2]
  rr <- (a / (a + b)) / (c / (c + d))
  se.log.rr <- sqrt((b / a) / (a + b) + (d / c) / (c + d))
  lci.rr <- exp(log(rr) - 1.96 * se.log.rr)
  uci.rr <- exp(log(rr) + 1.96 * se.log.rr)
  rr22.output <- list(estimate = rr, conf.int = c(lci.rr, uci.rr))
  class(rr22.output) <- "rr22"
  return(rr22.output)
  }
#
# Print method for objects of class "rr22"
#
print.rr22 <- function(x)
  {
  cat("RR     : ", x$estimate, "\n", 
      "95% CI : ", x$conf.int[1], "; ", x$conf.int[2], "\n",
      sep = "")
  }


################################################################################
#
# Bland & Altman plot
#
ba.plot <- function(a, b, title = "Bland and Altman Plot") 
  {
  a.txt <- deparse(substitute(a))
  b.txt <- deparse(substitute(b))
  x.lab <- paste("Mean of", a.txt, "and", b.txt)
  y.lab <- paste(a.txt, "-", b.txt)
  mean.two <- (a + b) / 2
  diff.two <- a - b
  plot(mean.two, diff.two, xlab = x.lab, ylab = y.lab, main = title, frame.plot = FALSE)
  mean.diff <- mean(diff.two)
  sd.diff <- sd(diff.two)
  upper <- mean.diff + 1.96 * sd.diff
  lower <- mean.diff - 1.96 * sd.diff
  lines(x = range(mean.two), y = c(mean.diff, mean.diff), lty = 3)     
  lines(x = range(mean.two), y = c(upper, upper), lty = 3)
  lines(x = range(mean.two), y = c(lower, lower), lty = 3)
  mean.text  <- round(mean.diff, digits = 1)
  upper.text <- round(upper , digits = 1)
  lower.text <- round(lower, digits = 1)
  text(max(mean.two), mean.diff, mean.text, adj = c(1, 1)) 
  text(max(mean.two), upper, upper.text, adj = c(1, 1))
  text(max(mean.two), lower, lower.text, adj = c(1, 1))
  ba <- list(mean = mean.diff, limits = c(lower, upper))
  class(ba) <- "ba"
  return(ba)
  }
#
# Print method for objects of class "ba"
#
print.ba <- function(x)
  {
  cat("Mean difference     : ", x$mean, "\n", 
      "Limits of agreement : ", x$limits[1], "; ", x$limits[2], "\n",
      sep = "")
  }


################################################################################
#
# Plot two (y) variables on a common (x) axis 
#
plot2var <- function(y1, y2, x.ticks, x.lab = deparse(substitute(x.ticks)), y1.lab = deparse(substitute(y1)), y2.lab = deparse(substitute(y2)), main = paste(y1.lab, "&", y2.lab))
  {
  old.par.mar <- par("mar")
  old.par.lab <- par("lab")
  par(mar = c(5, 5, 4, 5))
  if(!missing(x.ticks))
    {
    par(lab = c(length(x.ticks), 5, 7))
    }
  plot(y1, type = "l", lty = 1, axes = FALSE, xlab = "", ylab = "", main = main)
  axis(side = 2)
  mtext(text = y1.lab, side = 2, line = 3)
  par(new = TRUE)
  plot(y2, type = "l", lty = 2, axes = FALSE, ylab = "", xlab = x.lab)
  axis(side = 4)
  mtext(text = y2.lab, side = 4, line = 3)
  if(!missing(x.ticks))
    {
    axis(side = 1, labels = as.character(x.ticks), at = 1:length(x.ticks))
    } else {axis(side = 1)}
  par(mar = old.par.mar)
  par(lab = old.par.lab)
  }


################################################################################
#
# Pyramid plot 
#
pyramid.plot <- function(x, g, main = paste("Pyramid plot of", deparse(substitute(x)), "by", deparse(substitute(g))), xlab = paste(deparse(substitute(g)), "(", levels(g)[1], "/", levels(g)[2],")"), ylab = deparse(substitute(x)))
  {
  tab <- table(x, g); tab[ ,1] <- -tab[ ,1]
  barplot(tab, horiz = TRUE, beside = TRUE, space = c(0, -nrow(tab)), names.arg = c(dimnames(tab)$x, dimnames(tab)$x), xlim = c(min(tab) * 1.2, max(tab) * 1.2), col = "white", main = main, xlab = xlab, ylab = ylab, axes = FALSE)
  axis(side = 1, labels = abs(axTicks(side = 1)), at = (axTicks(side = 1)))
  }


################################################################################
#
# Pareto plot 
#
pareto <- function(x, xlab = deparse(substitute(x)), ylab = "Count", main = paste("Pareto Chart of", deparse(substitute(x))))
  {
  barplot(rev(sort(table(x))), xlab = xlab, ylab = ylab, main = main, col = "white")
  }

################################################################################
#
# XY plot with confidence intervals 
#
plot.ci <- function(x, y, y.lci, y.uci, ylim = c(min(y.lci), max(y.uci)), xlab = deparse(substitute(x)), ylab = deparse(substitute(y)), main = paste(ylab, "by", xlab), type = "l", lty = 1, col = "black", axes = TRUE, pch = 21, bg = "white")
  {
  plot(x, y, type = type, ylim = ylim, xlab = xlab, ylab = ylab, main = main, lty = lty, col = col, axes = axes, frame.plot = FALSE)
  arrows(x, y.lci, x, y.uci, code=3, angle=90, length=0.1, lty = lty, col = col)
  points(x, y, pch = pch, bg = bg, col = col)
  }


################################################################################
#
# Bar plot with confidence intervals 
#
barplot.ci <- function(y, bar.names, lci, uci, ylim = c(min(0, lci), max(0, uci)), xlab = deparse(substitute(bar.names)), ylab = deparse(substitute(y)), main = paste(ylab, "by", xlab))
  {
  ylim <- ylim * 1.1
  bp <- barplot(y, names.arg = bar.names, ylim = ylim, xlab = xlab, ylab = ylab, main = main, col = "white")
  arrows(bp, lci, bp, uci, code=3, angle=90, length=0.1)
  points(bp, y, pch = 21, bg = "white", col = "black")
  }
