source("tipping.R");
if (!("d8" %in% ls()))
{
  d8 <- ieeealife2020Analysis(8, 0:50 / 50);
}
allIeeeAlife2020Plots(d8);

