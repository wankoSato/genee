#' Regularized Regression.
#'
#' @param betas A vector of GWAS effect size estimators (OLS betas).
#' @param ld A matrix of pairwise correlations statistics (i.e., r).
#' @param alpha A tuning parameter determining the type of regularization. More specifically, alpha = 0: Ridge Regression; 0 < alpha < 1: Elastic Net; alpha = 1: LASSO.
#' @param nfolds A number indicates how many folds to do for cross validation for regularize regression. The default is 10. If it's fed with 0, then it will perform no cross validation and choose lambda with smallest penalty.
#' @return A vector of regularized effect size estimators.
#' @export
#' @examples
#' x1 = c(0, 1, 1)
#' x2 = c(0, 1, 2)
#' x3 = c(1, 1, 2)
#' x = cbind(x1 ,x2, x3)
#' ld = cor(x)
#' betas = c(0.1, 0.2, 0.1)
#' genee_regularize_regression(betas, ld, 0.5)
genee_regularize_regression<-function(betas, ld, alpha, nfolds){
  # fit a regularized model
  # alpha is the tuning parameter determining which type of the shrinkage is used
  # alpha = 0: Ridge Regression
  # 0 < alpha < 1: Elastic Net
  # alpha = 1: LASSO

  # Technically, glmnet assumes alpha = 1 for LASSO. However, the LASSO solution can unstable when setting alpha = 1. Thus, we recommend using 0.99 or 0.98 instead.
  if(nfolds == 0){
    fit = glmnet(x = ld, y = betas, alpha = alpha, intercept = FALSE)

    # use effect size with smallest penalty
    beta_regularized = as.numeric(fit$beta[,ncol(fit$beta)])
  }else{
    fit = cv.glmnet(x = ld, y = betas, alpha = alpha, intercept = FALSE, nfolds=5)
    # use effect size with smallest error in cross valiation
    beta_regularized = as.numeric(fit$glmnet.fit$beta[,which(fit$lambda==fit$lambda.min)])
  }

  # return results
  return(beta_regularized)
}
