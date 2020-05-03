library(deSolve);


epsdevice <- function(epsname, width = 8, height = 6)
{
  postscript(epsname, width = width, height = height, paper = "special", onefile = FALSE, horizontal = FALSE);
  par(cex = 1.5, cex.main = 1.3);
}


tippingCubic <- function(x, a, b, c)
{
  return(a * x - b * x * x * x + c);
}


tippingCubicExtrema <- function(a, b, c)
{
  return(c(-sqrt(a / 3 * b), sqrt(a / 3 * b)));
}


tippingCubicStableFixedPoints <- function(a, b, c)
{
  xExt  <- tippingCubicExtrema(a, b, c);
  yExt  <- tippingCubic(xExt, a, b, c);
  xRoot <- sort(Re(polyroot(c(c, a, 0, -b))));
  xStable <- numeric();
  if (yExt[1] < 0)
  {
    xStable <- c(xStable, xRoot[1]);
  }
  if (yExt[2] > 0)
  {
    xStable <- c(xStable, xRoot[3]);
  }
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


coupledTippingStableFixedPoints <- function(coupledTippingParams)
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


upstreamImpact <- function(coupledTippingParams, coupledTippingState, i)
{
  if (i == 1L)
  {
    return(0);
  }
  j <- 1L:(i - 1L);
  return(sum(coupledTippingParams$d[i, j] * coupledTippingState[j]));
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


coupledTippingStateIntSet <- function(coupledTippingParams)
{
  return(as.integer(lapply(coupledTippingStableFixedPoints(coupledTippingParams), discretiseCoupledTippingStateInt)));
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


agentImpactedCoupledTippingTimeSeries <- function(coupledTippingParams, actuation, agentImpact, initialState, nStepsImpact, nStepsPostImpact, dTime)
{
  agentImpactedCtp <- coupledTippingParams;
  agentImpactedCtp$c <- agentImpactedCtp$c + actuation * agentImpact;
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
    finalState <- discretiseCoupledTippingStateInt(s[nrow(s), 2:ncol(s)]);
    finalStateSet <- union(finalStateSet, finalState);
    actuation <- nextActuation(actuation);
    ## print(finalStateSet);
    ## print(actuation);
    ## print(coupledTippingParams$c);
    ## print(agentCpParams$c);
    ## plotODESeries(s, ylim=c(-2, 2));
    ## readline(sprintf("finalState: %d, hit return", finalState));
  }
  return(log2(length(finalStateSet)));
}


allTransientEmpowerment <- function(coupledTippingParams, agentImpact, nSteps, dTime)
{
  ## FIXME: function assumes that initialState is a valid fixed point
  fpList <- coupledTippingStableFixedPoints(coupledTippingParams);
  fpIntList <- integer();
  empList <- numeric();
  for (fp in fpList)
  {
    fpIntList <- c(fpIntList, discretiseCoupledTippingStateInt(fp));
    empList <- c(empList, transientEmpowerment(coupledTippingParams, agentImpact, fp, nSteps, dTime));
  }
  return(data.frame(initialState=fpIntList, transientEmpowerment=empList, stringsAsFactors=FALSE));
}


coupledTippingDemo <- function(coupledTippingParams, e, nSteps, dTime, fnamePattern=NULL)
{
  fpList <- coupledTippingStableFixedPoints(coupledTippingParams);
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


tippingAnalysis <- function(coupledTippingParams)
{
  fpList <- coupledTippingStableFixedPoints(coupledTippingParams);
  agentImpactTe <- agentImpactTESweep(coupledTippingParams, 0:20 / 20, fpList[[1]], 300, 0.1);
  return(invisible(list(coupledTippingParams=coupledTippingParams, fpList=fpList, agentImpactTe=agentImpactTe)));
}


toyAnalysis <- function()
{
  toyCtp <- makeToyCoupledTippingParams();
  toyFpList <- coupledTippingStableFixedPoints(toyCtp);
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
        numStatesCol <- c(numStatesCol, length(coupledTippingStableFixedPoints(p)));
      }
    }
  }
  return(data.frame(n=nCol, c=cCol, d=dCol, numStates=numStatesCol));
}


alife2020Figs <- function(n, d=NULL)
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
    d$adSmall <- tippingAnalysis(d$pdSmall);
    set.seed(6);
    d$pdMedium <- makeRandomCoupledTippingParams(n, -0.1, 0.1, -d$dMedium, d$dMedium);
    d$adMedium <- tippingAnalysis(d$pdMedium);
    set.seed(6);
    d$pdLarge <- makeRandomCoupledTippingParams(n, -0.1, 0.1, -d$dLarge, d$dLarge);
    d$adLarge <- tippingAnalysis(d$pdLarge);
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

if (!exists("d5"))
{
  d5 <- alife2020Figs(5);
} else {
  d5 <- alife2020Figs(5, d5);
}

