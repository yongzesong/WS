#' Wavelet linear regression (WLR)
#'
#' @description Function for wavelet linear regression.
#'
#' @usage wlr(y, x, wn, ns = 6)
#' \method{print}{wlr}(x, ...)
#' \method{plot}{wlr}(x, ...)
#' \method{predict}{wlr}(object, newdata, ...)
#'
#' @aliases wlr print.wlr plot.wlr predict.wlr
#'
#' @param y A vector of response variable
#' @param x A matrix of independent variables. Row: samples; Column: spectra
#' @param wn A vector of wavelength of x
#' @param ns Number of scale of multiresolution analysis. Default: ns = 6.
#' @param object A list of results
#' @param newdata A matrix of new independent variables for prediction
#' @param ... Ignore
#'
#' @import RcppEigen
#' @importFrom inline cxxfunction
#' @importFrom stats cov2cor cor
#' @importFrom methods signature
#' @importFrom waveslim mra
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_point scale_x_continuous scale_y_continuous xlab ylab theme_classic geom_bar scale_fill_brewer geom_line geom_vline geom_text facet_grid vars
#'
#' @examples
#' data("soilcarbon")
#' wn <- seq(350, 2500, 10)
#' # test
#' wlr1 <- wlr(y = y, x = x[,1], wn = wn, ns = 6)
#' # all data
#' \dontrun{wlr1 <- wlr(y = y, x = x, wn = wn, ns = 6)} # about 10s
#'
#' @export
#'

wlr <- function(y, x, wn, ns = 6){

  if (class(x)[1] == "numeric"){
    message("x should be a matrix")
    result <- list()
  } else {


    ##===========================================
    ## prepare functions
    ##===========================================

    ## fast linear regression
    lltLSCpp <- '
const MapMatd         X(as<MapMatd>(XX));
const MapVecd         y(as<MapVecd>(yy));
const int             n(X.rows()), p(X.cols());
const LLT<MatrixXd> llt(AtA(X));
const VectorXd  betahat(llt.solve(X.adjoint() * y));
const VectorXd   fitted(X * betahat);
const VectorXd    resid(y - fitted);
const int            df(n - p);
const double        ssq(resid.squaredNorm() / double(df));
const MatrixXd     vcov(ssq * llt.solve(MatrixXd::Identity(p, p)));
return     List::create(Named("coefficients")   = betahat,
                        Named("fitted.values")  = fitted,
                        Named("residuals")      = resid,
                        Named("s")              = sqrt(ssq),
                        Named("df.residual")    = df,
                        Named("rank")           = p,
                        Named("vcov")           = vcov);
'
    incl <- '
using   Eigen::LLT;
using   Eigen::Lower;
using   Eigen::Map;
using   Eigen::MatrixXd;
using   Eigen::MatrixXi;
using   Eigen::Upper;
using   Eigen::VectorXd;
typedef Map<MatrixXd>  MapMatd;
typedef Map<MatrixXi>  MapMati;
typedef Map<VectorXd>  MapVecd;
inline  MatrixXd AtA(const MatrixXd& A) {
  int    n(A.cols());
  return   MatrixXd(n,n).setZero().selfadjointView<Lower>()
                        .rankUpdate(A.adjoint());
}
inline  MatrixXd AAt(const MatrixXd& A) {
  int    n(A.cols());
  return   MatrixXd(n,n).setZero().selfadjointView<Lower>()
                        .rankUpdate(A);
}
'


lltLS <- cxxfunction(signature(XX = "matrix", yy = "numeric"), lltLSCpp, "RcppEigen", incl)

## vif
lltLS.vif <- function(mod) {
  v <- mod$vcov[-1,-1]
  n.terms <- nrow(v)
  assign <- 1:n.terms
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) *
      det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- result[, 1]
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  result
}

## select x variables based on multicollinearity < 10
# x: row: a group of observations
xselectmcl <- function(ys, xs, ctr.vif = 10){
  nys <- length(ys)
  xs <- t(xs)
  nxs <- ncol(xs)

  m <- c()
  for (i in 2:nxs){
    k <- 1:i
    if (length(m) > 0) {
      k <- k[-m]
    }

    x1 <- xs[, k]

    xx <- cbind(rep(1, nys), x1)
    f1 <- lltLS(XX = xx, y = ys)
    f1.vif <- lltLS.vif(f1)

    if (max(f1.vif) > ctr.vif){ ## remove multicollinearity
      m <- c(m, i)
    }
  }
  xk <- c(1:nxs)[-m]
  return(xk)
}

##===========================================
## prepare wavelet analysis
##===========================================

# select maximum wavelets in the analysis
L <- length(x[1,])
L <- floor(L/2^ns)*2^ns

# all combinations of wavelets: 'cow'
lambda <- c(as.character(2^c(1:ns)), "sv") # multiresolution

cow.lambda <- rep(lambda, each = L)
cow.wavenumber <- rep(1:L, ns + 1)
cow.raw.wn <- wn[cow.wavenumber]

# prepare folds of cross validation
n.total <- length(y)
k0 <- sample(n.total)
folds <- cut(seq(1,n.total), breaks = 10, labels = FALSE)[k0]

# statistical functions
MAE <- function(o, p) mean(abs(o - p))
RMSE <- function(o, p) sqrt(mean((o - p)^2))
R2 <- function(o, p) 1 - sum((o - p)^2)/sum((o - mean(o))^2)

##===========================================
## prepare x
##===========================================

# all spectra
xi <- list()
for (i in 1:n.total){
  spec.v1.la8 <- mra(x[i, ][1:L], "la8", ns, "dwt", boundary = "reflection")
  xi[[i]] <- spec.v1.la8
}

xi.lambda <- list()
for (i in 1:(ns + 1)){
  xi.lambda[[i]] <- do.call(cbind, lapply(xi, function(x) x[[i]]))
}
names(xi.lambda) <- lambda # important: x

xall <- do.call(rbind, xi.lambda)

cow.cor <- apply(xall, 1, FUN = function(u) cor(y, u))
result.cor <- data.frame(lambda = cow.lambda, wavenumber = cow.wavenumber, raw.wn = cow.raw.wn, cor = cow.cor)

cow.id <- 1:(L * (ns + 1))
la3 <- cow.id[rev(order(abs(cow.cor)))]
result.cor.rank <- la3

xnew <- xall[la3,]

##===========================================
## select x by vif
##===========================================

xk <- xselectmcl(ys = y, xs = xnew, ctr.vif = 10) # 15s: select x by vif
xselect <- t(xnew[xk, ])

result.selectbymcl <- xk

##===========================================
## select x by cv
##===========================================

## 10-fold cv to select best variables
n.xs <- length(xk)

xs.m0 <- matrix(rep(1:n.xs, each = n.xs), n.xs, n.xs)
xs.m1 <- upper.tri(xs.m0)
xs.m0[xs.m1] <- NA

cv.train.r2 <- c(); cv.test.r2 <- c(); cv.train.rmse <- c(); cv.test.rmse <- c()
for (i in 1:n.xs){
  sk <- xs.m0[i,]
  sk <- sk[which(!is.na(sk))]
  si <- xselect[, sk]
  sx <- cbind(rep(1, n.total), si)

  s.train.r2 <- c(); s.test.r2 <- c(); s.train.rmse <- c(); s.test.rmse <- c()

  for (kf in 1:10){
    testIndexes <- which(folds == kf, arr.ind = TRUE)

    y.train <- y[-testIndexes]; y.test <- y[testIndexes]
    x.train <- sx[-testIndexes, ]; x.test <- sx[testIndexes, ]
    n.train <- length(y.train); n.test <- length(y.test)

    f1 <- lltLS(XX = x.train, yy = y.train) ## fastest lm
    coef1 <- f1$coefficients

    y.train.pred <- x.train %*% coef1
    y.test.pred <- x.test %*% coef1

    s.train.r2[kf] <- R2(y.train, y.train.pred)
    s.test.r2[kf] <- R2(y.test, y.test.pred)
    s.train.rmse[kf] <- RMSE(y.train, y.train.pred)
    s.test.rmse[kf] <- RMSE(y.test, y.test.pred)
  }
  cv.train.r2[i] <- mean(s.train.r2); cv.test.r2[i] <- mean(s.test.r2);
  cv.train.rmse[i] <- mean(s.train.rmse); cv.test.rmse[i] <- mean(s.test.rmse);
}

cv.track <- data.frame(id = la3[xk], lambda = cow.lambda[la3[xk]],
                       wavenumber = cow.wavenumber[la3[xk]],
                       raw.wn = cow.raw.wn[la3[xk]],
                       cv.train.r2, cv.train.rmse, cv.test.r2, cv.test.rmse)

cv.k <- which(cv.test.r2 == max(cv.test.r2))

xk2 <- xk[1:cv.k]
result.selectbycv <- xk2

xselect2 <- t(xnew[xk2, ])
best.wavelet <- cv.track[1:cv.k, 1:4]

sx2 <- cbind(rep(1, n.total), xselect2)
f2 <- lltLS(XX = sx2, yy = y) ## fastest lm
y.pred2 <- sx2 %*% f2$coefficients
f2.vif <- lltLS.vif(f2)

xselect2 <- data.frame(xselect2)
names(xselect2) <- paste(best.wavelet$lambda, best.wavelet$raw.wn, sep = "_")
fitness <- list("R2" = R2(y, y.pred2), "MAE" = MAE(y, y.pred2), "RMSE" = RMSE(y, y.pred2))

result <- list("fun" = f2, "fitness" = fitness, "vif" = f2.vif,
               "best.x" = xselect2, "best.wavelets" = best.wavelet,
               "spec.cor" = result.cor, "spec.cor.rank" = result.cor.rank,
               "x.selectbyvif" = result.selectbymcl, "x.selectbycv" = result.selectbycv,
               "cv.track" = cv.track, "y" = y, "x" = x, "wn" = wn, "ns" = ns)

  }
## define class
class(result) <- "wlr"
result
}

print.wlr <- function(x, ...){
  cobw <- round(x$fun$coefficients, digits = 4)
  names(cobw) <- c("(Intercept)", names(x$best.x))
  cat("Coefficients of best wavelets:\n")
  print(cobw)
  invisible(x)
}

plot.wlr <- function(x, ...){
  ns <- x$ns
  lambda <- c(as.character(2^c(1:ns)), "sv")
  #######################################
  ## frequency of wavelets by R
  #######################################
  a <- x$spec.cor
  a <- a[x$spec.cor.rank,]
  a$abs.cor <- abs(a$cor)
  a$id <- 1:nrow(a)
  a$lambda <- factor(a$lambda, levels = lambda)
  names(a)[1] <- "Scale"
  plot1 <- ggplot(a, aes(x = id, y = abs.cor, color = Scale)) + geom_point() +
    scale_x_continuous() + scale_y_continuous() +
    xlab("Rank of wavelets") + ylab("|R|") +
    theme_classic()

  id <- 1:floor(nrow(a)/100)
  a1 <- seq(1, nrow(a), 100)
  a2 <- seq(100, nrow(a), 100)
  a1 <- a1[1:length(a2)]
  pctrn <- matrix(0, length(id), ns + 1)
  for (i in id){
    di <- a$Scale[a1[i]:a2[i]]
    pctrn[i,] <- table(di)
  }
  pctrn <- data.frame(pctrn)
  names(pctrn) <- lambda
  pctrn$id <- id * 100
  pctrn <- reshape2::melt(pctrn, id = "id")
  plot2 <- ggplot(data=pctrn, aes(x=id, y=value, fill=variable)) +
    geom_bar(stat="identity") +
    scale_fill_brewer(palette="Spectral", direction = -1, name = "Scale") +
    xlab("Rank of wavelets") +
    ylab("Frequency") +
    theme_classic()

  #######################################
  ## cross validation for selecting optimal combinations
  #######################################

  cv.track <- x$cv.track
  cvt.r2 <- data.frame(id = 1:nrow(cv.track), cv.track[, c(5,7)])
  names(cvt.r2) = c("id", "Training", "Testing")
  plot.cvt.r2 <- melt(cvt.r2, id = "id")
  plot3 <- ggplot(plot.cvt.r2, aes(x = id, y = value, color = variable)) +
    geom_point() + geom_line() +
    scale_x_continuous(limits = c(0, nrow(cv.track)), breaks = seq(0, nrow(cv.track), 4)) +
    scale_y_continuous() +
    geom_vline(xintercept = cvt.r2$id[which(cvt.r2$Testing == max(cvt.r2$Testing))]) +
    xlab("Number of wavelets") +
    ylab("Cross-validation R-squared") +
    theme_classic()

  #######################################
  ## distributions of optimal wavelets
  #######################################
  spec.v1 <- colMeans(x$x)
  FunDiv2J <- function(x, nscale){ ## use the dividing function before MRA
    L <- length(x)
    L <- floor(L/2^nscale)*2^nscale
    x[1:L]
  }
  spec.v1 <- FunDiv2J(spec.v1, nscale = ns)
  spec.v1.la8 <- mra(spec.v1, "la8", ns, "dwt", boundary = "reflection") ## debug: use reflection

  ## must shift LA(8) coefficients
  spec.v1.plot <- data.frame(cbind(wavelength = x$wn[1:length(spec.v1)], wavelet = spec.v1, do.call(cbind, spec.v1.la8))) #
  svp <- melt(spec.v1.plot, id = "wavelength")


  svp.b <- svp[-which(svp$variable == "wavelet"),]
  names(svp.b)[2] <- "lambda"
  svp.b$lambda <- as.character(svp.b$lambda)
  k <- match(svp.b$lambda, unique(svp.b$lambda))
  svp.b$lambda <- lambda[k]

  a1 <- x$cv.track
  a2 <- x$best.wavelets
  a1$value <- 0
  a2$value <- 0

  for (i in 1:nrow(a1)){
    k <- which(a1$lambda[i] == svp.b$lambda & a1$raw.wn[i] == svp.b$wavelength)
    a1$value[i] <- svp.b$value[k]
  }
  for (i in 1:nrow(a2)){
    k <- which(a2$lambda[i] == svp.b$lambda & a2$raw.wn[i] == svp.b$wavelength)
    a2$value[i] <- svp.b$value[k]
  }

  svp.b$lambda <- factor(svp.b$lambda, levels = lambda)
  a1$lambda <- factor(a1$lambda, levels = lambda)
  a2$lambda <- factor(a2$lambda, levels = lambda)

  plot4 <- ggplot(svp.b, aes(x = wavelength, y = value)) +
    geom_line() +
    geom_point(data = a1, aes(x = raw.wn, y = value), color = 'blue') +
    geom_text(data = a1, aes(x = raw.wn, y = value, label = raw.wn), color = 'blue') +
    facet_grid(rows = vars(lambda), scales = "free_y") +
    theme_classic()

  plot5 <- ggplot(svp.b, aes(x = wavelength, y = value)) +
    geom_line() +
    geom_point(data = a2, aes(x = raw.wn, y = value), color = '#EB1103') +
    geom_text(data = a2, aes(x = raw.wn, y = value, label = raw.wn), color = '#EB1103') +
    facet_grid(rows = vars(lambda), scales = "free_y") +
    theme_classic()

  cat("Ranked wavelets by absolute values of correlations ...\n")
  print(plot1)
  cat("Frequency of wavelet scales by the rank of wavelets ...\n")
  print(plot2)
  cat("Wavelets selected by multicollinearity analysis ...\n")
  print(plot4)
  cat("Ten-fold cross validation for selecting wavelets ...\n")
  print(plot3)
  cat("Optimal wavelets ...\n")
  print(plot5)
}

predict.wlr <- function(object, newdata, ...){
  ns <- object$ns
  n.total <- nrow(newdata)
  # select maximum wavelets in the analysis
  L <- length(newdata[1,])
  L <- floor(L/2^ns)*2^ns

  # all combinations of wavelets: 'cow'
  lambda <- c(as.character(2^c(1:ns)), "sv") # multiresolution

  xi <- list()
  for (i in 1:n.total){
    spec.v1.la8 <- mra(newdata[i, ][1:L], "la8", ns, "dwt", boundary = "reflection") ## x
    xi[[i]] <- spec.v1.la8
  }

  xi.lambda <- list()
  for (i in 1:(ns + 1)){
    xi.lambda[[i]] <- do.call(cbind, lapply(xi, function(x) x[[i]]))
  }
  names(xi.lambda) <- lambda # important: x

  xall <- do.call(rbind, xi.lambda)
  x1 <- t(xall[object$best.wavelets$id,])

  xx <- cbind(rep(1, n.total), x1)
  coef1 <- object$fun$coefficients
  y.pred <- as.vector(xx %*% coef1)

  x1 <- data.frame(x1)
  names(x1) <- paste(object$best.wavelets$lambda, object$best.wavelets$raw.wn, sep = "_")

  result <- list("y.pred" = y.pred, "x" = x1)
  return(result)
}
