library(deSolve);


epsdevice <- function(epsname, width = 8, height = 6)
{
  postscript(epsname, width = width, height = height, paper = "special", onefile = FALSE, horizontal = FALSE);
  par(cex = 1.5, cex.main = 1.3);
}


commaSep <- function(x, fmt)
{
  return(paste(sprintf(fmt, x), collapse=", "));
}


binaryAtoi <- function(s)
{
  return(sapply(strsplit(s, ""), function(d) { return(sum(as.integer(d) * 2L^((length(d) - 1L):0L))); }));
}


tippingCubic <- function(x, a, b, c)
{
  return(a * x - b * x * x * x + c);
}


tippingCubicDerivative1 <- function(x, a, b)
{
  return(as.numeric(a - 3 * b * x * x));
}


tippingCubicExtrema <- function(a, b, c)
{
  e <- sqrt(a / 3 * b);
  return(c(-e, e));
}


tippingCubicCcrit <- 2 / 3 / sqrt(3);


cardanoCubicRoot <- function(p, q)
{
  ## compute roots of x^3 + p * x + q
  ## based on https://en.wikipedia.org/wiki/Cubic_equation
  ## further stuff at https://math.vanderbilt.edu/schectex/courses/cubic/, https://mathworld.wolfram.com/CubicFormula.html
  ## d <- 4 * p * p * p + 27 * q * q;
  u <- (-1 + sqrt(-3 + 0i)) / 2.0;
  r2 <- q * q / 4.0 + p * p * p / 27.0 + 0i;
  ## message(sprintf("r2 = %f+%fi", Re(r2), Im(r2)));
  c1 <- (-(q + 0i) / 2 + sqrt(r2))^(1 / 3);
  c2  <- c1 * u;
  c3  <- c2 * u;
  return(c(c1 - p / (3.0 * c1), c2 - p / (3.0 * c2), c3 - p / (3.0 * c3)));
}


tippingCubicRoot <- function(a, b, c)
{
  return(cardanoCubicRoot(-b, 0, a, c));
}


cubicExtremaXY <- function(a, b, c)
{
  xExt <- tippingCubicExtrema(a, b, c);
  yExt <- tippingCubic(xExt, a, b, c);
  return(list(xExt=xExt, yExt=yExt));
}


findStable <- function(xExt, yExt, r)
{
  xStable <- numeric();
  if (yExt[1] < 0)
  {
    xStable <- min(r);
  }
  if (yExt[2] > 0)
  {
    xStable <- c(xStable, max(r));
  }
  return(xStable);
}


wrapPolyroot <- function(x)
{
  return(polyroot(x));
}


tippingCubicStableFixedPoints <- function(a, b, c)
{
  xyExt <- cubicExtremaXY(a, b, c);
  xExt  <- xyExt$xExt;
  yExt  <- xyExt$yExt;
  ## print(sprintf("a = %f, b = %f, c = %f", a, b, c));
  ## r <- Re(polyroot(c(c, a, 0, -b)));
  r <- Re(wrapPolyroot(c(c, a, 0, -b)));
  xStable  <- findStable(xExt, yExt, r);
  return(xStable);
}


plotTippingCubic <- function(a, b, c, ...)
{
  
  x <- -1000:1000 / 500;
  # FIXME: code repetition to tippingPointDiffEq
  y <- tippingCubic(x, a, b, c);
  xExt  <- tippingCubicExtrema(a, b, c);
  yExt  <- tippingCubic(xExt, a, b, c);
  xStable <- tippingCubicStableFixedPoints(a, b, c);
  yStable <- tippingCubic(xStable, a, b, c);
  plot(x, y, type="l", ...);
  lines(c(-3, 3), c(0, 0), col="blue");
  points(xExt, yExt, col="red");
  points(xStable, yStable, col="green");
}


makeToyCoupledTippingParams <- function()
{
  toyTippingParams <- list(a=c(1, 1, 1), b=c(1, 1, 1), c=c(0, 0.1, -0.1), d=matrix(c(0, 0.1, 0.2, 0, 0, -0.1, 0, 0, 0), nrow=3, ncol=3));
  names(toyTippingParams$a) <- names(toyTippingParams$b) <- names(toyTippingParams$c) <- rownames(toyTippingParams$d)  <- colnames(toyTippingParams$d) <- sprintf("x%02d", 1:length(toyTippingParams$a));
  return(toyTippingParams);
}


makeToyHoppingCoupledTippingParams <- function()
{
  n <- 3L;
  d <- matrix(0, nrow=n, ncol=n);
  d[2, 1] <- 1.0;
  d[3, 2] <- 1.5;
  params <- list(a=rep(1.0, n), b=rep(1.0, n), c=c(0, -2, 2.5), d=d);
  names(params$a) <- names(params$b) <- names(params$c) <- rownames(params$d)  <- colnames(params$d) <- sprintf("x%02d", 1:n);
  return(params);
}


hoppingDemo <- function(ctp, agentImpact, nSteps, dtime)
{
  ## ctp <- makeToyHoppingCoupledTippingParams();
  n <- coupledTippingDim(ctp);
  fpList <- coupledTippingAllStableFixedPoints(ctp);
  actuation <- rep(0, n);
  actuation[1] <- 1;
  s <- agentImpactedCoupledTippingTimeSeries(ctp, actuation, agentImpact, fpList[[1]], nSteps, nSteps, dtime);
  plotODESeries(s, ylim=c(-2, 2));
  return(invisible(s));
}


makeHoppingCoupledTippingParams <- function(n)
{
  d <- matrix(0, nrow=n, ncol=n);
  stop("unfinished");
}


makeRandomCoupledTippingParams <- function(n, cMin, cMax, dMin, dMax)
{
  cRange <- cMax - cMin;
  dRange <- dMax - dMin;
  d <- matrix(0, nrow=n, ncol=n);
  d[lower.tri(d)] <- runif(n * (n - 1) / 2) * dRange + dMin;
  params <- list(a=rep(1.0, n), b=rep(1.0, n), c=runif(n) * cRange + cMin, d=d);
  names(params$a) <- names(params$b) <- names(params$c) <- rownames(params$d)  <- colnames(params$d) <- sprintf("x%02d", 1:n);
  return(params);
}


makeRandomChainCoupledTippingParams <- function(n, cMin, cMax, dMin, dMax)
{
  cRange <- cMax - cMin;
  dRange <- dMax - dMin;
  d <- matrix(0, nrow=n, ncol=n);
  r <- runif(n - 1) * dRange + dMin;
  for (i in 1L:(n - 1))
  {
    d[i + 1, i] <- r[i];
  }
  params <- list(a=rep(1.0, n), b=rep(1.0, n), c=runif(n) * cRange + cMin, d=d);
  names(params$a) <- names(params$b) <- names(params$c) <- rownames(params$d)  <- colnames(params$d) <- sprintf("x%02d", 1:n);
  return(params);
}


coupledTippingDim <- function(coupledTippingParams)
{
  return(length(coupledTippingParams$a));
}


coupledTippingDimnames <- function(coupledTippingParams)
{
  return(names(coupledTippingParams$a));
}


makeCoupledTippingInitialState <- function(y, coupledTippingParams)
{
  if (length(y) != coupledTippingDim(coupledTippingParams))
  {
    stop("incompatible dimensions");
  }
  names(y) <- names(coupledTippingParams$a);
  return(y);
}


coupledTippingImpact <- function(y, d)
{
  dyCoupling <- rep(0, length(y));
  for (i in seq(along=y))
  {
    dyCoupling[i] <- sum(d[i, ] * y);
  }
  return(dyCoupling);
}


coupledTippingDiffEq <- function(t, y, params)
{
  ## following notation in Klose2019_interactingtippingelements, with y rather than x:
  ## dy[i] = a[i] * y[i] - b[i] * y[i]^3 + c[i] + sum(d[j][i] * y[j])
  ## the condition j != i in the sum is not enforced, users must ensure that d[i][i] = 0
  dy  <-  tippingCubic(y, params$a, params$b, params$c) + coupledTippingImpact(y, params$d);
  return(list(dy));
}


upstreamImpact <- function(coupledTippingParams, coupledTippingState, i)
{
  if (i == 1L)
  {
    return(0);
  }
  j <- 1L:(i - 1L);
  return(sum(coupledTippingParams$d[i, j] * coupledTippingState[j]));
}


coupledTippingStableFixedPoint <- function(coupledTippingParams, initialState)
{
  n <- coupledTippingDim(coupledTippingParams);
  coupledTippingState <- numeric();
  for (i in 1L:n)
  {
    u <- upstreamImpact(coupledTippingParams, coupledTippingState, i);
    ## message(sprintf("i = %d, u = %f, coupledTippingState = %s", i, u, paste(sprintf("%f", coupledTippingState), collapse=", ")));
    tcsfp <- tippingCubicStableFixedPoints(coupledTippingParams$a[i], coupledTippingParams$b[i], coupledTippingParams$c[i] + u);
    if (length(tcsfp) == 1L)
    {
      ## message(sprintf("i=%d: single state: %f", i, tcsfp));
      coupledTippingState <- c(coupledTippingState, tcsfp);
    }
    else
    {
      y <- tippingCubic(initialState[i], coupledTippingParams$a[i], coupledTippingParams$b[i], coupledTippingParams$c[i] + u);
      ## message(sprintf("i=%d: multiple states: %s, initialState[%d]=%f, u=%f, y=%f", i, paste(sprintf("%f", tcsfp), collapse=", "), i, initialState[i], u, y));
      if (initialState[i] <= tcsfp[1L])
      {
        coupledTippingState <- c(coupledTippingState, tcsfp[1L]);
      }
      else if (initialState[i] >= tcsfp[2L])
      {
        coupledTippingState <- c(coupledTippingState, tcsfp[2L]);
      }
      else
      {
        if (y <= 0)
        {
          coupledTippingState <- c(coupledTippingState, tcsfp[1L]);
        }
        else
        {
          coupledTippingState <- c(coupledTippingState, tcsfp[2L]);
        }
      }
    }
  }
  return(coupledTippingState);
}


coupledTippingAllStableFixedPoints <- function(coupledTippingParams)
{
  if (any(coupledTippingParams$d[upper.tri(coupledTippingParams$d, diag=TRUE)] != 0))
  {
    stop("d is not a lower triangular matrix");
  }
  n <- coupledTippingDim(coupledTippingParams);
  fixedPointPrefixList <- as.list(tippingCubicStableFixedPoints(coupledTippingParams$a[1L], coupledTippingParams$b[1L], coupledTippingParams$c[1L]));
  if (n > 1L)
  {
    for (i in 2L:n)
    {
      newPrefixList <- list();
      for (fixedPointPrefix in fixedPointPrefixList)
      {
        upstreamImpact <- sum(coupledTippingParams$d[i, 1L:(i - 1L)] * fixedPointPrefix);
        for (fixedPoint in tippingCubicStableFixedPoints(coupledTippingParams$a[i], coupledTippingParams$b[i], coupledTippingParams$c[i] + upstreamImpact))
        {
          newPrefixList[[length(newPrefixList) + 1L]] <- c(fixedPointPrefix, fixedPoint);
        }
      }
      fixedPointPrefixList <- newPrefixList;
    }
  }
  return(fixedPointPrefixList);
}


plotCoupledTippingCubics <- function(coupledTippingParams, coupledTippingState, ...)
{
  cpDim <- coupledTippingDim(coupledTippingParams);
  cpDimnames <- coupledTippingDimnames(coupledTippingParams);
  for (i in 1L:cpDim)
  {
    u <- upstreamImpact(coupledTippingParams, coupledTippingState, i);
    plotTippingCubic(coupledTippingParams$a[i], coupledTippingParams$b[i], coupledTippingParams$c[i] + u, main=cpDimnames[i], ...);
    readline(sprintf("%s -- hit return", cpDimnames[i]));
  }
}


coupledTippingTimeSeries <- function(coupledTippingParams, nSteps, dTime, y0=NULL)
{
  tpDim <- coupledTippingDim(coupledTippingParams);
  if (is.null(y0))
  {
    y0 <- rnorm(tpDim);
  }
  y0 <- makeCoupledTippingInitialState(y0, coupledTippingParams);
  ## y0 <- c(x = -1e-3);
  y <- rk4(y0, 0:(nSteps - 1) * dTime, coupledTippingDiffEq, coupledTippingParams);
  return(y);
}


discretiseCoupledTippingState <- function(y)
{
  s <- sign(y);
  s[s < 0] <- 0;
  s <- as.integer(s);
  names(s) <- names(y);
  return(s);
}


discretiseCoupledTippingStateStr <- function(y)
{
  s <- discretiseCoupledTippingState(y);
  return(paste(as.character(s), collapse=""));
}


discretiseCoupledTippingStateInt <- function(y)
{
  s <- discretiseCoupledTippingState(y);
  ## print(s);
  i <- as.integer(0);
  for (j in seq(along=s))
  {
    if (s[length(y) - j + 1] != 0)
    {
      i <- i + 2L^(j - 1L);
      ## print(i);
    }
  }
  return(i);
}


coupledTippingStateIntSet <- function(coupledTippingParams)
{
  return(as.integer(lapply(coupledTippingAllStableFixedPoints(coupledTippingParams), discretiseCoupledTippingStateInt)));
}


coupledTippingDot <- function(coupledTippingParams, dotFname)
{
  n <- coupledTippingDim(coupledTippingParams);
  d <- coupledTippingDimnames(coupledTippingParams);
  l <- "digraph coupledTipping {";
  for (i in 1L:n)
  {
    l <- c(l, sprintf("%s;", d[i]));
    if (i > 1L)
    {
      for (j in 1:(i - 1L))
      {
        l <- c(l, sprintf("%s -> %s;", d[j], d[i]));
      }
    }
  }
  l <- c(l, "}");
  return(l);
}


plotODESeries <- function(odeResult, ...)
{
  odeComponentList <- colnames(odeResult);
  odeComponentList <- odeComponentList[!(odeComponentList %in% "time")];
  opar <- par(no.readonly=TRUE);
  par(mfrow=c(length(odeComponentList), 1));
  for (odeComponent in odeComponentList)
  {
    plot(odeResult[, "time"], odeResult[, odeComponent], type = "l", main=odeComponent, xlab="time", ylab=odeComponent, ...);
    ## readline(sprintf("%s -- hit return", odeComponent));
  }
  par(opar);
}


nextActuation <- function(actuation)
{
  i <- length(actuation);
  while (i > 0L)
  {
    if (actuation[i] == 0)
    {
      actuation[i] <- -1;
      return(actuation);
    }
    else if (actuation[i] == -1)
    {
      actuation[i] <- 1;
      return(actuation);
    }
    else if (actuation[i] == 1)
    {
      actuation[i] <- 0;
      i <- i - 1L;
    }
  }
  return(NULL);
}


agentImpactedCoupledTippingParams <- function(coupledTippingParams, actuation, agentImpact)
{
  agentImpactedCtp <- coupledTippingParams;
  agentImpactedCtp$c <- agentImpactedCtp$c + actuation * agentImpact;
  return(agentImpactedCtp);
}


agentImpactedCoupledTippingTimeSeries <- function(coupledTippingParams, actuation, agentImpact, initialState, nStepsImpact, nStepsPostImpact, dTime)
{
  agentImpactedCtp <- agentImpactedCoupledTippingParams(coupledTippingParams, actuation, agentImpact);
  s <- coupledTippingTimeSeries(agentImpactedCtp, nStepsImpact, dTime, initialState);
  s <- rbind(s, coupledTippingTimeSeries(coupledTippingParams, nStepsPostImpact, dTime, s[nrow(s), 2:ncol(s)]));
  s[, "time"] <- 0:(nrow(s) - 1) * dTime;
  return(s);
}


transientEmpowerment <- function(coupledTippingParams, agentImpact, initialState, nSteps, dTime)
{
  n <- coupledTippingDim(coupledTippingParams);
  finalStateSet <- integer();
  cpDim <- coupledTippingDim(coupledTippingParams);
  actuation <- rep(0, cpDim);
  while (!is.null(actuation))
  {
    s <- agentImpactedCoupledTippingTimeSeries(coupledTippingParams, actuation, agentImpact, initialState, nSteps, nSteps, dTime) ;
    finalState <- s[nrow(s), 2:ncol(s)];
    ## print(finalState);
    finalStateInt <- discretiseCoupledTippingStateInt(finalState);
    finalStateSet <- union(finalStateSet, finalStateInt);
    actuation <- nextActuation(actuation);
    ## print(finalStateSet);
    ## print(actuation);
    ## print(coupledTippingParams$c);
    ## print(agentCpParams$c);
    ## plotODESeries(s, ylim=c(-2, 2));
    ## readline(sprintf("finalState: %d, hit return", finalStateInt));
  }
  return(log2(length(finalStateSet)));
}


agentReachableStateSetAnalytic <- function(coupledTippingParams, agentImpact, initialState)
{
  n <- coupledTippingDim(coupledTippingParams);
  finalStateSet <- integer();
  cpDim <- coupledTippingDim(coupledTippingParams);
  actuation <- rep(0, cpDim);
  while (!is.null(actuation))
  {
    agentImpactedCtp <- agentImpactedCoupledTippingParams(coupledTippingParams, actuation, agentImpact);
    agentImpactedFinalState = coupledTippingStableFixedPoint(agentImpactedCtp, initialState);
    relaxedFinalState = coupledTippingStableFixedPoint(coupledTippingParams, agentImpactedFinalState);
    ## message(sprintf("agentImpactedFinalState: %s, relaxedFinalState: %s", commaSep(agentImpactedFinalState, "%f"), commaSep(relaxedFinalState, "%f")));
    finalState <- discretiseCoupledTippingStateInt(relaxedFinalState);
    finalStateSet <- union(finalStateSet, finalState);
    ## print(finalStateSet);
    actuation <- nextActuation(actuation);
    ## print(finalStateSet);
    ## print(actuation);
    ## print(coupledTippingParams$c);
    ## print(agentCpParams$c);
    ## plotODESeries(s, ylim=c(-2, 2));
    ## readline(sprintf("finalState: %d, hit return", finalState));
  }
  return(finalStateSet);
}


stateTransitionMatrix <- function(coupledTippingParams, agentImpact)
{
  fpList <- coupledTippingAllStableFixedPoints(coupledTippingParams);
  numStates <- length(fpList);
  fpListInt <- sapply(fpList, discretiseCoupledTippingStateInt);
  m <- matrix(FALSE, nrow=numStates, ncol=numStates);
  rownames(m) <- colnames(m) <- sapply(fpList, discretiseCoupledTippingStateStr);
  for (i in seq(along=fpList))
  {
    m[i, ] <- fpListInt %in% agentReachableStateSetAnalytic(coupledTippingParams, agentImpact, fpList[[i]]);
  }
  return(m);
}


stateTransitionMatrixEmpowerment <- function(sm)
{
  return(data.frame(initialState=binaryAtoi(rownames(sm)), transientEmpowerment=log2(rowSums(sm)), sustainableEmpowerment=log2(sapply(1L:nrow(sm), function(i) { return(sum(sm[i, ] & sm[, i])); })), stringsAsFactors=FALSE));
}


stateTransitionMatrixSweep <- function(coupledTippingParams, agentImpactList)
{
  l <- list();
  for (agentImpact in agentImpactList)
  {
    l[[length(l) + 1]] <- stateTransitionMatrix(coupledTippingParams, agentImpact);
  }
  attr(l, "coupledTippingParams") <- coupledTippingParams;
  attr(l, "agentImpactList") <- agentImpactList;
  return(l);
}


stateTransitionMatrixEmpowermentSweep <- function(smList)
{
  agentImpactList <- attr(smList, "agentImpactList");
  d <- NULL;
  for (i in seq(along=agentImpactList))
  {
    a <- stateTransitionMatrixEmpowerment(smList[[i]]);
    a[["agentImpact"]] <- rep(agentImpactList[i], nrow(a));
    if (is.null(d))
    {
      d <- a;
    }
    else
    {
      d <- rbind(d, a);
    }
  }
  return(d);
}


transientEmpowermentAnalytic <- function(coupledTippingParams, agentImpact, initialState)
{
  return(log2(length(agentReachableStateSetAnalytic(coupledTippingParams, agentImpact, initialState))));
}


allTransientEmpowerment <- function(coupledTippingParams, agentImpact)
{
  fpList <- coupledTippingAllStableFixedPoints(coupledTippingParams);
  fpIntList <- integer();
  empList <- numeric();
  for (fp in fpList)
  {
    fpIntList <- c(fpIntList, discretiseCoupledTippingStateInt(fp));
    empList <- c(empList, transientEmpowermentAnalytic(coupledTippingParams, agentImpact, fp));
  }
  return(data.frame(initialState=fpIntList, transientEmpowerment=empList, stringsAsFactors=FALSE));
}


allTransientEmpowermentSweep <- function(coupledTippingParams, agentImpactList)
{
  d <- NULL;
  for (agentImpact in agentImpactList)
  {
    a <- allTransientEmpowerment(coupledTippingParams, agentImpact);
    a[["agentImpact"]] <- rep(agentImpact, nrow(a));
    if (is.null(d))
    {
      d <- a;
    }
    else
    {
      d <- rbind(d, a);
    }
  }
  return(d);
}


coupledTippingDemo <- function(coupledTippingParams, e, nSteps, dTime, fnamePattern=NULL)
{
  fpList <- coupledTippingAllStableFixedPoints(coupledTippingParams);
  fp <- fpList[[1]];
  for (i in 1:coupledTippingDim(coupledTippingParams))
  {
    ctpMod <- coupledTippingParams;
    ctpMod$c[i]  <- ctpMod$c[i] + e;
    s <- coupledTippingTimeSeries(coupledTippingParams, nSteps, dTime, fp);
    s <- rbind(s, coupledTippingTimeSeries(ctpMod, nSteps, dTime, s[nrow(s), 2:ncol(s)]));
    s <- rbind(s, coupledTippingTimeSeries(coupledTippingParams, nSteps, dTime, s[nrow(s), 2:ncol(s)]));
    s[, "time"] <- 0:(nrow(s) - 1) * dTime;
    if (!is.null(fnamePattern))
    {
      fname <- sprintf(fnamePattern, coupledTippingDimnames(coupledTippingParams)[i]);
      epsdevice(fname, width=6, height=10);
      par(mar=c(3, 5, 2, 3));
    }
    plotODESeries(s, ylim=c(-2.5, 2.5));
    if (!is.null(fnamePattern))
    {
      dev.off();
    }
    else
    {
      readline(sprintf("%s -- hit return", coupledTippingDimnames(coupledTippingParams)[i]));
    }
  }
}


coupledTippingDemoFigs <- function()
{
  set.seed(3);
  p10 <- makeRandomCoupledTippingParams(10, -0.1, 0.1, -0.2, 0.2);
  coupledTippingDemo(p10, 0.5, 1000, 0.01, fnamePattern="coupledtippingdemo_%s.eps");
}


agentImpactTESweep <- function(coupledTippingParams, agentImpactList, initialState, nSteps, dTime)
{
  return(data.frame(agentImpact=agentImpactList, transientEmpowerment=sapply(agentImpactList, function(agentImpact) { return(transientEmpowerment(coupledTippingParams, agentImpact, initialState, nSteps, dTime)); })));
}


agentImpactTESweepAnalytic <- function(coupledTippingParams, agentImpactList, initialState)
{
  return(data.frame(agentImpact=agentImpactList, transientEmpowerment=sapply(agentImpactList, function(agentImpact) { return(transientEmpowermentAnalytic(coupledTippingParams, agentImpact, initialState)); })));
}


tippingAnalysis <- function(coupledTippingParams, nSteps, dTime)
{
  fpList <- coupledTippingAllStableFixedPoints(coupledTippingParams);
  agentImpactTe <- agentImpactTESweep(coupledTippingParams, 0:20 / 20, fpList[[1]], nSteps, dTime);
  return(invisible(list(coupledTippingParams=coupledTippingParams, fpList=fpList, agentImpactTe=agentImpactTe)));
}


tippingAnalysisAnalytic <- function(coupledTippingParams)
{
  fpList <- coupledTippingAllStableFixedPoints(coupledTippingParams);
  ## agentImpactTe <- agentImpactTESweepAnalytic(coupledTippingParams, 0:20 / 20, fpList[[1]]);
  ## agentImpactTe <- allTransientEmpowermentSweep(coupledTippingParams, 0:20 / 20);
  stateTransitionMatrixList <- stateTransitionMatrixSweep(coupledTippingParams, 0:20 / 20);
  agentImpactTe <- stateTransitionMatrixEmpowermentSweep(stateTransitionMatrixList);
  return(invisible(list(coupledTippingParams=coupledTippingParams, fpList=fpList, stateTransitionMatrixList=stateTransitionMatrixList, agentImpactTe=agentImpactTe)));
}


toyAnalysis <- function()
{
  toyCtp <- makeToyCoupledTippingParams();
  toyFpList <- coupledTippingAllStableFixedPoints(toyCtp);
  toyAgentImpact <- agentImpactTESweep(toyCtp, 0:20 / 20, toyFpList[[1]], 500, 0.1);
  barplot(toyAgentImpact$transientEmpowerment, names.arg=sprintf("%3.1f", toyAgentImpact$agentImpact));
  return(invisible(list(toyCtp=toyCtp, toyFpList=toyFpList, toyAgentImpact=toyAgentImpact)));
  ## tippingAnalysis(toyCtp);
}


fixedPointScan <- function(nList, cList, dList)
{
  nCol <- integer();
  cCol <- numeric();
  dCol <- numeric();
  numStatesCol <- integer();
  for (n in nList)
  {
    for (c in cList)
    {
      for (d in dList)
      {
        p <- makeRandomCoupledTippingParams(n, -c, c, -d, d);
        nCol <- c(nCol, n);
        cCol <- c(cCol, c);
        dCol <- c(dCol, d);
        numStatesCol <- c(numStatesCol, length(coupledTippingAllStableFixedPoints(p)));
      }
    }
  }
  return(data.frame(n=nCol, c=cCol, d=dCol, numStates=numStatesCol));
}


empowermentBarplot <- function(e, numNodes, numStates=NULL, ...)
{
  barplot(e, ylim=c(0, numNodes), ...);
  if (!is.null(numStates))
  {
    eMax <- log2(numStates);
    lines(c(-1, length(e) * 1.2 + 1L), c(eMax, eMax), col="blue");
  }
}


plotIeeeAlife2020Analysis <- function(d)
{
  opar <- par(no.readonly=TRUE);
  par(mfrow=c(3, 2));
  empowermentBarplot(d8$adSmall$agentImpactTe$transientEmpowerment, 8, length(d8$adSmall$fpList), main="standard");
  empowermentBarplot(d8$adChainSmall$agentImpactTe$transientEmpowerment, 8, length(d8$adChainSmall$fpList), main="chain");
  empowermentBarplot(d8$adMedium$agentImpactTe$transientEmpowerment, 8, length(d8$adMedium$fpList), main="standard");
  empowermentBarplot(d8$adChainMedium$agentImpactTe$transientEmpowerment, 8, length(d8$adChainMedium$fpList), main="chain");
  empowermentBarplot(d8$adLarge$agentImpactTe$transientEmpowerment, 8, length(d8$adLarge$fpList), main="standard");
  empowermentBarplot(d8$adChainLarge$agentImpactTe$transientEmpowerment, 8, length(d8$adChainLarge$fpList), main="chain");
  par(opar);
}


ieeealife2020Analysis <- function(n)
{
  d <- list();
  d$n <- n;
  ## d$nSteps <- 300L;
  ## d$dTime <- 0.1;
  d$dSmall <- 0.2
  ## d$dMedium <- 0.4;
  d$dLarge <- 0.8;
  set.seed(6);
  d$ctpSmall <- makeRandomCoupledTippingParams(n, -0.1, 0.1, -d$dSmall, d$dSmall);
  d$adSmall <- tippingAnalysisAnalytic(d$ctpSmall);
  ## set.seed(6);
  ## d$ctpMedium <- makeRandomCoupledTippingParams(n, -0.1, 0.1, -d$dMedium, d$dMedium);
  ## d$adMedium <- tippingAnalysisAnalytic(d$ctpMedium);
  set.seed(6);
  d$ctpLarge <- makeRandomCoupledTippingParams(n, -0.1, 0.1, -d$dLarge, d$dLarge);
  d$adLarge <- tippingAnalysisAnalytic(d$ctpLarge);
  set.seed(6);
  d$ctpChainSmall <- makeRandomChainCoupledTippingParams(n, -0.1, 0.1, -d$dSmall, d$dSmall);
  d$adChainSmall <- tippingAnalysisAnalytic(d$ctpChainSmall);
  ## set.seed(6);
  ## d$ctpChainMedium <- makeRandomChainCoupledTippingParams(n, -0.1, 0.1, -d$dMedium, d$dMedium);
  ## d$adChainMedium <- tippingAnalysisAnalytic(d$ctpChainMedium);
  set.seed(6);
  d$ctpChainLarge <- makeRandomChainCoupledTippingParams(n, -0.1, 0.1, -d$dLarge, d$dLarge);
  d$adChainLarge <- tippingAnalysisAnalytic(d$ctpChainLarge);
  return(invisible(d));
}


alife2020Figs <- function(n, d=NULL)
{
  if (is.null(d))
  {
    d <- list();
    d$n <- n;
    d$nSteps <- 300L;
    d$dTime <- 0.1;
    d$dSmall <- 0.2
    d$dMedium <- 0.4;
    d$dLarge <- 0.8;
    set.seed(6);
    d$pdSmall <- makeRandomCoupledTippingParams(n, -0.1, 0.1, -d$dSmall, d$dSmall);
    d$adSmall <- tippingAnalysis(d$pdSmall, d$nSteps, d$dTime);
    set.seed(6);
    d$pdMedium <- makeRandomCoupledTippingParams(n, -0.1, 0.1, -d$dMedium, d$dMedium);
    d$adMedium <- tippingAnalysis(d$pdMedium, d$nSteps, d$dTime);
    set.seed(6);
    d$pdLarge <- makeRandomCoupledTippingParams(n, -0.1, 0.1, -d$dLarge, d$dLarge);
    d$adLarge <- tippingAnalysis(d$pdLarge, d$nSteps, d$dTime);
  }
  else if (d$n != n)
  {
    stop("incompatible n");
  }
  epsdevice(sprintf("empprof%02ds.eps", d$n));
  barplot(d$adSmall$agentImpactTe$transientEmpowerment, names.arg=sprintf("%3.1f", d$adSmall$agentImpactTe$agentImpact), las=2, sub=sprintf("d_ji in [%3.1f, %3.1f]", -d$dSmall, d$dSmall), xlab="E", ylab="empowerment", ylim=c(0, d$n));
  dev.off();
  epsdevice(sprintf("empprof%02dm.eps", d$n));
  barplot(d$adMedium$agentImpactTe$transientEmpowerment, names.arg=sprintf("%3.1f", d$adMedium$agentImpactTe$agentImpact), las=2, sub=sprintf("d_ji in [%3.1f, %3.1f]", -d$dMedium, d$dMedium), xlab="E", ylab="empowerment", ylim=c(0, d$n));
  dev.off();
  epsdevice(sprintf("empprof%02dl.eps", d$n));
  barplot(d$adLarge$agentImpactTe$transientEmpowerment, names.arg=sprintf("%3.1f", d$adLarge$agentImpactTe$agentImpact), las=2, sub=sprintf("d_ji in [%3.1f, %3.1f]", -d$dLarge, d$dLarge), xlab="E", ylab="empowerment", ylim=c(0, d$n));
  dev.off();
  return(invisible(d));
}


alife2020FigsAlt <- function(n, d=NULL)
{
  if (is.null(d))
  {
    d <- list();
    d$n <- n;
    d$nSteps <- 1000L;
    d$dTime <- 0.1;
    d$dSmall <- 0.2
    d$dMedium <- 0.4;
    d$dLarge <- 0.8;
    set.seed(6);
    d$pdSmall <- makeRandomCoupledTippingParams(n, -0.1, 0.1, -d$dSmall, d$dSmall);
    d$adSmall <- tippingAnalysis(d$pdSmall, d$nSteps, d$dTime);
    set.seed(6);
    d$pdMedium <- makeRandomCoupledTippingParams(n, -0.1, 0.1, -d$dMedium, d$dMedium);
    d$adMedium <- tippingAnalysis(d$pdMedium, d$nSteps, d$dTime);
    set.seed(6);
    d$pdLarge <- makeRandomCoupledTippingParams(n, -0.1, 0.1, -d$dLarge, d$dLarge);
    d$adLarge <- tippingAnalysis(d$pdLarge, d$nSteps, d$dTime);
  }
  else if (d$n != n)
  {
    stop("incompatible n");
  }
  epsdevice(sprintf("alt_empprof%02ds.eps", d$n));
  barplot(d$adSmall$agentImpactTe$transientEmpowerment, names.arg=sprintf("%3.1f", d$adSmall$agentImpactTe$agentImpact), las=2, sub=sprintf("d_ji in [%3.1f, %3.1f]", -d$dSmall, d$dSmall), xlab="E", ylab="empowerment", ylim=c(0, d$n));
  dev.off();
  epsdevice(sprintf("alt_empprof%02dm.eps", d$n));
  barplot(d$adMedium$agentImpactTe$transientEmpowerment, names.arg=sprintf("%3.1f", d$adMedium$agentImpactTe$agentImpact), las=2, sub=sprintf("d_ji in [%3.1f, %3.1f]", -d$dMedium, d$dMedium), xlab="E", ylab="empowerment", ylim=c(0, d$n));
  dev.off();
  epsdevice(sprintf("alt_empprof%02dl.eps", d$n));
  barplot(d$adLarge$agentImpactTe$transientEmpowerment, names.arg=sprintf("%3.1f", d$adLarge$agentImpactTe$agentImpact), las=2, sub=sprintf("d_ji in [%3.1f, %3.1f]", -d$dLarge, d$dLarge), xlab="E", ylab="empowerment", ylim=c(0, d$n));
  dev.off();
  return(invisible(d));
}


alife2020FigsAnalytic <- function(n, d=NULL)
{
  if (is.null(d))
  {
    d <- list();
    d$n <- n;
    d$dSmall <- 0.2
    d$dMedium <- 0.4;
    d$dLarge <- 0.8;
    set.seed(6);
    d$pdSmall <- makeRandomCoupledTippingParams(n, -0.1, 0.1, -d$dSmall, d$dSmall);
    d$adSmall <- tippingAnalysisAnalytic(d$pdSmall);
    set.seed(6);
    d$pdMedium <- makeRandomCoupledTippingParams(n, -0.1, 0.1, -d$dMedium, d$dMedium);
    d$adMedium <- tippingAnalysisAnalytic(d$pdMedium);
    set.seed(6);
    d$pdLarge <- makeRandomCoupledTippingParams(n, -0.1, 0.1, -d$dLarge, d$dLarge);
    d$adLarge <- tippingAnalysisAnalytic(d$pdLarge);
  }
  else if (d$n != n)
  {
    stop("incompatible n");
  }
  epsdevice(sprintf("aempprof%02ds.eps", d$n));
  barplot(d$adSmall$agentImpactTe$transientEmpowerment, names.arg=sprintf("%3.1f", d$adSmall$agentImpactTe$agentImpact), las=2, sub=sprintf("d_ji in [%3.1f, %3.1f]", -d$dSmall, d$dSmall), xlab="E", ylab="empowerment", ylim=c(0, d$n));
  dev.off();
  epsdevice(sprintf("aempprof%02dm.eps", d$n));
  barplot(d$adMedium$agentImpactTe$transientEmpowerment, names.arg=sprintf("%3.1f", d$adMedium$agentImpactTe$agentImpact), las=2, sub=sprintf("d_ji in [%3.1f, %3.1f]", -d$dMedium, d$dMedium), xlab="E", ylab="empowerment", ylim=c(0, d$n));
  dev.off();
  epsdevice(sprintf("aempprof%02dl.eps", d$n));
  barplot(d$adLarge$agentImpactTe$transientEmpowerment, names.arg=sprintf("%3.1f", d$adLarge$agentImpactTe$agentImpact), las=2, sub=sprintf("d_ji in [%3.1f, %3.1f]", -d$dLarge, d$dLarge), xlab="E", ylab="empowerment", ylim=c(0, d$n));
  dev.off();
  return(invisible(d));
}


allTransientEmpowermentDemo <- function(n, cRange, dRange, agentImpactList)
{
  set.seed(3L);
  coupledTippingParams <- makeRandomCoupledTippingParams(n, -cRange, cRange, -dRange, dRange);
  atp <- allTransientEmpowermentSweep(coupledTippingParams, agentImpactList);
  return(invisible(list(coupledTippingParams=coupledTippingParams, atp=atp)));
}


statePerturbationAnalysis <- function(coupledTippingParams, initialState, standardDeviation, nSteps, dTime)
{
  perturbedState <- initialState + rnorm(length(initialState), sd=standardDeviation);
  ctpTs <- coupledTippingTimeSeries(coupledTippingParams, nSteps, dTime, perturbedState);
  ## plotODESeries(ctpTs, ...);
  sse <- sum(apply(ctpTs[, 2:ncol(ctpTs)], 2L, function(x) { xDiff <- diff(x); return(sum(xDiff * xDiff)); }));
  return(invisible(list(initialState=initialState, ctpTs=ctpTs, sse=sse)));
}


statePerturbationAllSse <- function(coupledTippingParams, standardDeviation, numRepeats, nSteps, dTime)
{
  fixedPointList <- coupledTippingAllStableFixedPoints(coupledTippingParams);
  l <- list();
  for (initialState in fixedPointList)
  {
    sseList <- numeric();
    for (r in 1:numRepeats)
    {
      sse <- statePerturbationAnalysis(coupledTippingParams, initialState, standardDeviation, nSteps, dTime)$sse;
      sseList <- c(sseList, sse);
    }
    l[[length(l) + 1L]] <- sseList;
  }
  return(l);
}


coupledTippingStateDeriv1SquaredSum <- function(coupledTippingParams, state)
{
  ss <- 0.0;
  n <- coupledTippingDim(coupledTippingParams);
  for (i in 1L:n)
  {
    d1 <- tippingCubicDerivative1(state[i], coupledTippingParams$a[i], coupledTippingParams$b[i]);
    ss <- ss + d1 * d1;
  }
  return(ss);
}
