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
wlr1 <- wlr(y = y, x = x, wn = wn, ns = 6)
wlr1
plot(wlr1)
```


&nbsp; 
