## Wavelet linear regression (WLR) for wavelet selection and soil prediction


```{r, eval = FALSE}
## install and library the pacakge
install.packages("devtools")

library(devtools)
install_github("yongzesong/WS")

library("WS")

## Example 1
data("soilcarbon")
wn <- seq(350, 2500, 10)
wlr1 <- wlr(y = y, x = x, wn = wn, ns = 6) # ~ 9.5s
wlr1
plot(wlr1)
```


&nbsp; 

To cite "WS" R package in publications, please use:

Song, Y., Shen, Z., Wu, P., and  R. A. V. Rossel. Wavelet geographically weighted regression for spectroscopic modelling of soil properties. Sci Rep 11, 17503 (2021). https://doi.org/10.1038/s41598-021-96772-z
