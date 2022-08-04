

# Non-Monotone Unified Estimate in R
# Summary -----------------------------------------------------------------
# The nmuer() function fits the unified estimate to any combination of
# glm() and coxph() model specifications.
#
# The required arguments are:
# - data: a data frame containing one row for each individual
# - main_model: either a glm() or coxph() function that contains a formula
#               for the model and other optional arguments (eg glm family)
# - working_models: a list() of glm() and/or coxph() functions that have
#                   formulas. They may have other arguments like family
# *Note: data arguments are not required for the main_model or
#   working_models, only for the main nmuer() function.
#
# The optional arguments are:
# - main_weight_model: a function that returns weights for the main_model
# - working_weight_models: a list of functions that return weights for
#   the working models.
# - num_boots: the number of bootstraps to estimate variance. the default
#              is 500
# - DEBUG: Setting to TRUE will  print extra information. The code can
#          be edited to use browser() to step through the function.
#
# By default the nmuer() function assumes MCAR and uses constant weights.
# The other option is to use nmuer_weights_logistic() with a
# "missing_formula" argument that gives the logistic formula.
# See Examples for explanation of how to use nmuer_weights_logistic().
# Other weight options are easy to add in the future.
#
#
# If one working_weight_models is supplied, then all must be supplied.
# The same number of working_weight_models as the number of working_models.
# You can use the function nmuer_weights_constant() to use constant weights.
# See example.
#
# The data must not have any variables named "R", "R0", "main_weights",
# or "working_weights".
#
# The output is:
# - beta_bar, the unified estimate of parameter
# - beta_bar_var, the estimated covariance matrix of beta_bar.
# - variance estimate, the type of variance estimation used
# - num_boots, the number of bootstraps used for variance estimate
# - unadjusted_betahat, the original beta_hat without unified adjustment
# -------------------------------------------------------------------------


# TODO: -------------------------------------------------------------------
# Clean up warnings in nmuer_weights_logistic()
# Double check tied event times for Coxph, I think it works right now
# See if we can use "data masking" to avoid name conflicts
# More checks that data is supplied with appropriate missing/observed
# More checks that user supplies working models correctly
# Use coxph.control?
# -------------------------------------------------------------------------


# Code Details ------------------------------------------------------------
#' Non-Monotone Unified Estimate
#'
#' @description  This function fits the unified estimate for either a Cox
#' or a GLM model. It uses bootstrapping to estimate variance.
#'
#' @details   TODO
#'
#' @export
#'
#' @importFrom stats  na.fail update.formula coef complete.cases cov2cor
#' fitted.values glm binomial pchisq printCoefmat
#'
#' @param main_model (Required) This is the main analysis model. This should
#' be either a glm() or coxph() function which includes a formula but does
#' not include a data argument
#' @param working_models (Required) A list containing the working models,
#' similar to the main_model specification.
#' @param data (Required) A dataframe with 1 row for each individual.
#' @param main_weight_model (Optional) A function that takes at minimum a
#' Data argument and a observation_vector argument, and outputs a vector of
#' weights. These weights are used in the main_model.
#' @param working_weight_models (Optional) A list of functions that each
#' take DATA and observation_vector arguments and output a vector of
#' weights.
#' @param variance_estimate (Optional) There are two different ways of
#' estimating the variance. Only "bootstrap" is currently allowed. 
#' "bootstrap" uses bootstrap resampling on rows of the data to estimate
#' the variance. "sandwich" will use a robust sandwich estimate of the form
#' A^-1 * B * A^T-1, where A is the matrix of partial
#' derivatives of the score with respect to the parameters and B is the
#' matrix of crossproducts of the score vectors. 
#' @param num_boots (Optional) The number of bootstrap replications to use
#' when using the bootstrap method of estimating variance. 200 was
#' recommended as the minimum by Efron, pg 52.
#' @param DEBUG (Optional) Prints extra model fitting and call information
#' @examples 
#' if(require(mice)) {
#'   nmuer_glm_mod <- nmuer(
#'     main_model = glm(hm ~ age*sex),
#'     working_models = list(glm(hr ~ age*sex)),
#'     data = mice::selfreport)
#'   summary(nmuer_glm_mod)
#' }
nmuer <- function(
    main_model,
    working_models,
    data,
    main_weight_model = NULL,
    working_weight_models = NULL,
    variance_estimate = c("bootstrap","sandwich"),
    num_boots = 500,
    DEBUG = FALSE
) {
  # Start debugger --------------------------------------------------------
  # Uncomment to enable stepping through the code when DEBUG is TRUE
  # if(DEBUG == TRUE) {
  #   browser()
  # }
  # Read and verify input -------------------------------------------------
  cal <- match.call()
  variance_estimate <- match.arg(variance_estimate)
  if(variance_estimate == "sandwich") {
    stop("Sandwich variance not implemented yet")
  }

  if(any(c("R0","R", "main_weights","working_weights") %in%
         names(data))) {
    stop("Conflict between nmuer internal variables and
         variables defined in data. Please rename any variables named
         R, R0, main_weights, or working_weights")
  }

  # Standardize model specifications --------------------------------------
  # Need to capture the code with enexpr before proceeding, otherwise R
  # will try to run the code and encounter errors.
  main_model <- rlang::enexpr(main_model)
  working_models <- rlang::enexpr(working_models)
  main_weight_model <- rlang::enexpr(main_weight_model)
  working_weight_models <- rlang::enexpr(working_weight_models)

  # If main_weight_model or working_weight_models are NULL, then use
  # the constant weight functions. IE, assume MCAR.
  if(is.null(main_weight_model)) {
    main_weight_model <- rlang::expr(nmuer_weights_constant())
  }
  if(is.null(working_weight_models)) {
    working_weight_models <- list()
    for(i in 1:length(working_models)) {
      working_weight_models[[i]] <- rlang::expr(nmuer_weights_constant())
    }
  }

  if(DEBUG == TRUE) {
    print("Before call_standardise()")
    print("main model:")
    print(main_model)
    print("working models:")
    print(working_models)
    print("main_weight_model:")
    print(main_weight_model)
    print("working_weight_models:")
    print(working_weight_models)
  }
  if(length(working_models) != length(working_weight_models)) {
    stop("working_models and working_weight_models must be the same length")
  }

  main_model <- rlang::call_standardise(main_model)
  main_weight_model <- rlang::call_standardise(main_weight_model)
  for(i in 2:length(working_models)) {
    working_models[[i]] <- rlang::call_standardise(working_models[[i]])
    working_weight_models[[i]] <-
      rlang::call_standardise(working_weight_models[[i]])
  }
  if(DEBUG == TRUE) {
    print("After call_standa")
    print("main model:")
    print(main_model)
    print("working models:")
    print(working_models)
    print("main_weight_model:")
    print(main_weight_model)
    print("working_weight_models:")
    print(working_weight_models)
  }
  # Check if models are adequately specified ------------------------------
  nmuer_check_specification(main_model)
  for(i in 2:length(working_models)) {
    nmuer_check_specification(working_models[[i]])
  }

  # Check if there is missing data ----------------------------------------
  # There should be some missing data in the main model, and the
  # working models should not have exactly the same observed data
  # as the main model
  R0 <- stats::complete.cases(data[,all.vars(main_model$formula)])
  if(all(R0)) {
    stop("the main_model has no missing data")
  }
  R <- matrix(data = NA, nrow = nrow(data), ncol = length(working_models)-1)
  for(i in 1:(length(working_models)-1)) {
    R[,i] <- stats::complete.cases(data[,all.vars(working_models[[i+1]])])
    if(all(R0 == R[,i])) {
      stop(paste0("main_model missingness exactly matches some ",
                  "working_model missingness"))
    }
  }

  # Fit models ------------------------------------------------------------
  main_fit <-
    nmuer_fitcoef(data = data,
                  main_model = main_model,
                  working_models = working_models,
                  main_weight_model = main_weight_model,
                  working_weight_models = working_weight_models,
                  return_fit = TRUE,
                  DEBUG = DEBUG)
  all_params <- main_fit[c("main_params",
                           "working_params_cc",
                           "working_params_full")]
  num_params <- c(length(all_params$main_params),
                  length(all_params$working_params_cc),
                  length(all_params$working_params_full))
  num_working_params <- main_fit$num_working_params
  param_vec <- c(all_params$main_params,
                 all_params$working_params_cc,
                 all_params$working_params_full)
  # if(DEBUG == TRUE) {
  #   browser()
  # }

  # Fit Variance ----------------------------------------------------------
  # Note: If any bootstrap parameters are NA, it tries to bootstrap again
  maxit <- 100
  if(variance_estimate == "bootstrap") {
    boot_params <- matrix(data = 0,
                          nrow = num_boots,
                          ncol = sum(num_params))
    for(i in 1:num_boots) {
      continue <- TRUE
      curit <- 1
      while(continue) {
        if(curit > maxit) {
          stop("Bootstrapping failed because of NA parameters")
        }
        boot_data <- data[sample(1:nrow(data),replace=TRUE),]
        boot_params[i,] <- unlist(nmuer_fitcoef(
          data = boot_data,
          main_model = main_model,
          working_models = working_models,
          main_weight_model = main_weight_model,
          working_weight_models = working_weight_models,
          return_fit = FALSE))
        curit <- curit + 1
        if(!(any(is.na(boot_params[i,])))) {
          continue <- FALSE
        }
      }
    }
    mean_boot_est <- colMeans(boot_params)
    boot_var <- matrix(data = 0,
                       nrow = sum(num_params),
                       ncol = sum(num_params))
    for(i in 1:num_boots) {
      boot_var <- boot_var +
        tcrossprod(boot_params[i,] - mean_boot_est)
    }
    sigma <- boot_var/(num_boots - 1)
  } else {
    stop("only bootstrapping variance is currently supported")
  }

  theta_cor_matrix <- stats::cov2cor(sigma)
  if(DEBUG == TRUE) {
    print("the estimated theta covariance matrix is:")
    print(sigma)
    print("the estimated theta correlation matrix is:")
    print(theta_cor_matrix)
  }
  # calculate the correlations between beta_hat and gamma_hat and gamma_bar
  beta_gammahat_cors <- list()
  beta_gammabar_cors <- list()

  for(i in seq_along(num_working_params)) {
    beta_gammahat_cors[[i]] <-
      theta_cor_matrix[1:num_params[1],
                       (num_params[1]+1):
                         (num_params[1]+sum(num_working_params[0:i]))]
    beta_gammabar_cors[[i]] <-
      theta_cor_matrix[1:num_params[1],
                       (num_params[1]+1+num_params[2]):
                         (num_params[1]+num_params[2] +
                            sum(num_working_params[0:i]))]
  }



  # Calculate Variance ----------------------------------------------------
  # Multiply the covariance matrix of theta by a permutation matrix to
  # get the covariance matrix of (beta, gamma_hat - gamma_bar), called
  # omega
  permutation_matrix <-
    matrix(data = 0,
           nrow = (num_params[1] + num_params[2]),
           ncol = sum(num_params))
  permutation_matrix[1:num_params[1],
                     1:num_params[1]] <-
    diag(1, num_params[1])
  permutation_matrix[(num_params[1] + 1):(num_params[1] + num_params[2]),
                     (num_params[1] + 1):(num_params[1] + num_params[2])] <-
    diag(1, num_params[2])
  permutation_matrix[(num_params[1] + 1):(num_params[1] + num_params[2]),
                     (num_params[1] + num_params[2] + 1):sum(num_params)] <-
    diag(-1, num_params[2])
  omega <- permutation_matrix %*% sigma %*% t(permutation_matrix)
  if(DEBUG == TRUE) {
    print("the estimated omega covariance matrix is:")
    print(omega)
    print("the estimated omega correlation matrix is:")
    print(stats::cov2cor(omega))
  }


  # Calculate beta_bar ----------------------------------------------------
  beta_hat <- all_params$main_params
  gamma_hat <- all_params$working_params_cc
  gamma_bar <- all_params$working_params_full
  omega_11 <- omega[
    1:num_params[[1]],
    1:num_params[[1]]]
  omega_12 <- omega[
    1:num_params[[1]],
    (num_params[[1]] + 1):(num_params[[1]] + num_params[[2]])]
  omega_22 <- omega[
    (num_params[[1]] + 1):(num_params[[1]] + num_params[[2]]),
    (num_params[[1]] + 1):(num_params[[1]] + num_params[[2]])]
  beta_bar <- beta_hat - (omega_12 %*% solve(omega_22) %*%
                            (gamma_hat - gamma_bar))


  # create output ---------------------------------------------------------
  final_fit <- main_fit$main_model_fit
  final_fit$coefficients <- beta_bar



  # Calculate beta_bar variance -------------------------------------------
  beta_var_reduction <- omega_12 %*% solve(omega_22) %*% t(omega_12)
  beta_bar_var <- omega_11 - beta_var_reduction
  # Below code gives exactly the same variance estimate, as expected.
  # lambda_transform_matrix <-
  #   matrix(0,
  #          nrow = num_params[[1]],
  #          ncol = num_params[[1]] + num_params[[2]])
  # lambda_transform_matrix[1:num_params[[1]],
  #                         1:num_params[[1]]] <-
  #   diag(1,
  #        nrow = num_params[[1]],
  #        ncol = num_params[[1]])
  # lambda_transform_matrix[
  #   1:num_params[[1]],
  #   (num_params[[1]]+1):(num_params[[1]] + num_params[[2]])] <-
  #   -(omega_12 %*% solve(omega_22))
  # beta_bar_var <-
  #   lambda_transform_matrix %*% omega %*% t(lambda_transform_matrix)

  # Add names to output ---------------------------------------------------
  beta_bar_names <- names(all_params$main_params)
  beta_bar <- as.vector(beta_bar)
  names(beta_bar) <- beta_bar_names
  dimnames(beta_bar_var) <- list(beta_bar_names,beta_bar_names)

  # Return Output ---------------------------------------------------------
  out <- list(main_model_type = main_model[[1]],
              final_fit = final_fit,
              beta_bar = beta_bar,
              beta_bar_var = beta_bar_var,
              variance_estimate = variance_estimate,
              num_boots = num_boots,
              unadjusted_betahat = beta_hat,
              call = cal,
              num_working_models = length(working_models)-1,
              beta_gammahat_cors = beta_gammahat_cors,
              beta_gammabar_cors = beta_gammabar_cors,
              unadjusted_betahat_var = omega_11,
              beta_var_reduction = beta_var_reduction)
  class(out) <- "nmuer"
  return(out)
}

#' Weighting functions for unified estimate
#' 
#' @description 
#' Functions to supply inverse probability weights for unified
#' estimates
#'
#' @details
#' When missingness is MAR, inverse probability weights (ipw) can
#' be used to rebalance the estimating equations to allow for unbiased
#' estimation of parameters. These functions are designed to allow 
#' nmuer() to automatically estimate the weights. Users can also write
#' their own weight functions. The two functions supplied here allow
#' for constant weights and weights estimated by a logistic model.
#' 
#' Each weight function MUST take an argument called DATA and an argument
#' called observation_vector. The observation_vector and DATA should be 
#' left NULL, the nmuer() code automatically calculates and includes them.
#' Other arguments may be defined depending on the method of estimating 
#' weights.
#' 
#' Each weight function must output a vector of weights the same
#' length as the number of rows in DATA. These weights are then
#' supplied to the model fitting functions glm() and coxph().
#' 
#' nmuer_weights_constant() gives all observations the same weight.
#' This is appropriate for MCAR data. The user does not need to
#' enter any additional information to use this method. This is the
#' default behaviour if no weight functions are supplied, but explicitly
#' calling nmuer_weights_constat() is required if you want to combine
#' some constant weights with some non-constant weights.
#' 
#' nmuer_weights_logistic() uses a GLM model with a binomial family and
#' a logistic link function to estimate the probability a row has complete
#' data, then takes the inverse of that probability as the weight. A
#' missing_formula argument is required. The right-hand side of the 
#' formula will be used to model observation probability. The left-hand
#' side of the formula will be automatically replaced by nmuer() to use
#' the appropriate observation vector. If the working_model has no
#' missing data, then nmuer_weights_logistic() will return a vector of 1s
#' with a warning.
#'
#' @export
#' 
#' @param DATA A dataframe on which the weighting function will
#' be fit to calculate the weights. This is automatically supplied
#' by nmuer()
#' @param observation_vector A logical vector indicating
#' if the ith row has all variables in the analysis or working model
#' observed. This is automatically supplied by nmuer().
#' @param missing_formula A formula, the right-hand side defines
#' the variables to be used in the logistic model. This must be supplied
#' by the user.
nmuer_weights_constant <- function(observation_vector = NULL,
                                   DATA = NULL) {
  if(is.null(observation_vector) | is.null(DATA)) {
    stop("Missing observation_vector or DATA")
  }
  return(rep(1,nrow(DATA)))
}

#' @rdname nmuer_weights_constant
#' @export
nmuer_weights_logistic <- function(missing_formula,
                                   observation_vector = NULL,
                                   DATA = NULL) {
  if(is.null(observation_vector) | is.null(DATA)) {
    stop("Missing observation_vector or DATA")
  }
  if(all(observation_vector)) {
    warning(paste0("No missing values in nmuer_weights_logistic ",
                   "detected. Using constant weights."))
    return(rep(1,nrow(DATA)))
  }
  if("R0" %in% names(DATA)) {
    stop("R0 already exists in DATA")
  }
  DATA$R0 <- observation_vector
  environment(missing_formula) <- environment()
  missing_formula <- stats::update.formula(old = missing_formula,
                                           new = R0 ~.)
  if(any(!stats::complete.cases(DATA[,all.vars(missing_formula)]))) {
    stop("A variable in missing_formula has missing values")
  }
  mod <- stats::glm(formula = missing_formula,
                    family = stats::binomial(link = "logit"),
                    data = DATA,
                    na.action = stats::na.fail)
  weights <- 1/stats::fitted.values(mod)
  if(any(weights > 1000)) {
    warning("Weights > 1000 detected in nmuer_weights_logistic")
  }
  return(weights)
}

# nmuer_check_specification checks if the supplied model has enough
# information to be fit.
# For now I only allow glm() and coxph().
# They must have a formula.
# In the future other models will probably be allowed, and they
# may have other options.
nmuer_check_specification <- function(model) {
  if(!(as.character(model[[1]]) %in% c("glm","coxph"))) {
    stop("only glm and coxph models are currently supported")
  }
  if(is.null(model$formula)) {
    stop("all analysis and working models must have a formula")
  }
  # if(as.character(model[[1]]) == "geeglm") {
  #   stop("geeglm models are not currently supported")
  #   if(!(id %in% names(model))) {
  #     stop("geeglm models must have an id specified")
  #   }
  # }
}

# nmuer_fitcoef function takes the validated model specifications and data
# and fits the models. It calculates weights (if they were specified) and
# fits the models on the appropriate subsets of the data.
# It returns only the coefficients from the analysis
# and working models. It does not calculate variance.
nmuer_fitcoef <- function(
  data,
  main_model,
  working_models,
  main_weight_model,
  working_weight_models,
  return_fit = FALSE,
  DEBUG = FALSE
) {
  if(DEBUG == TRUE) {
    browser()
  }
  # Calculate weights for the main model and complete case working
  # models. This happens even if missingness is MCAR.
  main_vars <- all.vars(main_model$formula)
  R0 <- stats::complete.cases(data[,main_vars])
  main_weight_model$observation_vector <- R0
  main_weight_model$DATA <- data
  main_weights <- eval(main_weight_model,envir = data)

  # Calculate weights for available case working models
  # Note: i starts at 2, because of strangeness in rlang package.
  # Note: if there is no missingness in working_models, then
  # nmuer_weights_logistic() will return all 1s.
  working_vars <- list()
  R <- matrix(data = FALSE, nrow = nrow(data), ncol = (length(working_models)-1))
  working_weights <-
    matrix(data = 0, nrow = nrow(data), ncol = (length(working_models)-1))
  for(i in 2:length(working_models)) {
    working_vars[[i-1]] <- all.vars(working_models[[i]]$formula)
    R[,(i-1)] <- stats::complete.cases(data[,working_vars[[i-1]]])
    working_weight_models[[i]]$observation_vector <- R[,(i-1)]
    working_weight_models[[i]]$DATA <- data
    working_weights[,(i-1)] <- eval(working_weight_models[[i]],
                                    envir = data)
  }
  working_models_cc <- working_models #uses only analysis complete cases
  working_models_full <- working_models #uses any available cases

  # Modify function calls with the appropriate weights, data, and subset.
  # main_model$weights <- main_weights
  # main_model$data <- data
  # main_model$subset <- R0
  # main_model$na.action <- stats::na.fail
  main_model$weights <- quote(main_weights)
  main_model$data <- quote(data)
  main_model$subset <- quote(R0)
  main_model$na.action <- quote(stats::na.fail)


  for(i in 2:length(working_models)) {
    working_models_cc[[i]]$weights <- main_weights
    working_models_cc[[i]]$data <- data
    working_models_cc[[i]]$subset <- R0
    working_models_cc[[i]]$na.action <- stats::na.fail
    working_models_full[[i]]$weights <- working_weights[,(i-1)]
    working_models_full[[i]]$data <- data
    working_models_full[[i]]$subset <-  R[,(i-1)]
    working_models_full[[i]]$na.action <- stats::na.fail
  }

  # If returning the model to use as final_fit, include the x matrix
  main_model$x <- return_fit

  # Fit the models
  main_model_fit <- eval(main_model, envir = data)
  working_model_fits_cc <-
    vector(mode = "list", length = length(working_models))
  working_model_fits_full <-
    vector(mode = "list", length = length(working_models))
  for(i in 2:length(working_models)) {
    working_model_fits_cc[[i]] <-
      eval(working_models_cc[[i]], envir = data)
    working_model_fits_full[[i]] <-
      eval(working_models_full[[i]], envir = data)
  }

  # Extract parameters and return results
  main_params <- stats::coef(main_model_fit)
  working_params_cc <- c()
  working_params_full <- c()
  num_working_params <- c()
  for(i in 2:length(working_models)) {
    working_params_cc <- c(working_params_cc,
                           stats::coef(working_model_fits_cc[[i]]))
    working_params_full <- c(working_params_full,
                             stats::coef(working_model_fits_full[[i]]))
    num_working_params <- c(num_working_params,
                            length(stats::coef(working_model_fits_cc[[i]])))
  }
  out <- list(main_params = main_params,
              working_params_cc = working_params_cc,
              working_params_full = working_params_full)
  if(return_fit == TRUE) {
    out$main_model_fit <- main_model_fit
    out$num_working_params <- num_working_params
  }
  return(out)
}



# S3 Methods --------------------------------------------------------------
#' @export
coef.nmuer <- function(object, ...) {
  object$beta_bar
}
#' @export
vcov.nmuer <- function(object, ...) {
  object$beta_bar_var
}
#' @export
summary.nmuer <- function(object, 
                          digits = max(1, getOption("digits") - 2), ...) {
  print(object, digits = digits, ...)
}
#' @export
print.nmuer <- function(x, digits = max(1, getOption("digits") - 2), ...) {
  cat("A non-monotone unified estimate with a ", 
      deparse(x$main_model_type), " analysis model and ", 
      x$num_working_models,  " working model",
      ifelse(x$num_working_models == 1, ".\n", "s.\n"),
      sep = "")
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if(x$main_model_type == "coxph") {
    cat("coefficients:\n")
    coef <- coef(x)
    se <- sqrt(diag(x$beta_bar_var))
    result <- cbind(coef, se, exp(coef), coef/se,
                    pchisq((coef/se)^2,
                           1, lower.tail = FALSE))
    dimnames(result) <- list(names(coef(x)),
                             c("coef","se(coef)","exp(coef)", "z", "p"))
    stats::printCoefmat(x = result, digits = digits, signif.stars = TRUE,
                        P.values = TRUE, has.Pvalue = TRUE)
  } else if(x$main_model_type == "glm") {
    cat("coefficients:\n")
    coef <- coef(x)
    se <- sqrt(diag(x$beta_bar_var))
    result <- cbind(coef, se, coef/se,
                    pchisq((coef/se)^2,
                           1, lower.tail = FALSE))
    dimnames(result) <- list(names(coef(x)),
                             c("coef","se(coef)", "z", "p"))
    stats::printCoefmat(x = result, digits = digits, signif.stars = TRUE,
                        P.values = TRUE, has.Pvalue = TRUE)
  } else {
    print("unsupported main_model type")
  }
  invisible(x)
}




