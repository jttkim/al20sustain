source("tipping.R");


if (!exists("d5"))
{
  d5 <- alife2020Figs(5);
} else {
  d5 <- alife2020Figs(5, d5);
}
