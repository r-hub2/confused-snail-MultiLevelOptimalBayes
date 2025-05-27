#' Multi-Level Optimal Bayes Function (MLOB)
#'
#'
#' Implements a regularized Bayesian approach that optimizes
#' the estimation of between-group coefficients by minimizing 
#' Mean Squared Error (MSE), balancing both variance and bias.
#' This method provides more reliable estimates 
#' in scenarios with limited data, offering a robust solution for 
#' accurate parameter estimation in multilevel models. The package is designed for researchers 
#' in psychology, education, and related fields who face challenges in 
#' estimating between-group effects in two-level latent variable models, particularly
#' in scenarios with small sample sizes and low intraclass correlation coefficients. 
#'
#' @param formula an object of class "\link{formula}" (or one that can be coerced to that class): a symbolic description of the model to be fitted. Formula specifies the model (e.g., \code{Y ~ X + C...}), where Y is the dependent variable, X is the context variable, which is the focus of most applications of the model  (always included), and C includes all additional covariates.
#' @param data a data frame (or object converted by \link{as.data.frame} to a data frame) containing the variables referenced in the formula. All variables used in the model, including the dependent variable, context variable, covariates, and grouping variable must be present in this data frame.
#' @param group a name of the variable that defines the affiliation of an individual (row) to the specific group.
#' @param balancing.limit a number that represents the threshold of the maximum relative part of the dataset that can be deleted to balance the data. Defaults to \code{0.2}
#' @param conf.level a numeric value representing the confidence level used to calculate confidence intervals for the estimators. Defaults to \code{0.95}, corresponding to a \code{95\%} confidence level.
#' @param jackknife logical variable. If \code{TRUE}, the jackknife re-sampling method will be applied to  calculate the standard error of the between-group and its confidence interval coefficient. Defaults to \code{FALSE}.
#' @param punish.coeff a multiplier that punishes the balancing procedure when deleting the whole group. If punish.coeff is equal to \code{1}, no additional punishment is applied for deleting the group. Higher values intensify the penalty. Defaults to \code{2}.
#' @param ... additional arguments passed to the function.
#'
#' @details
#' This function also verifies whether the data is balanced (i.e., whether each group contains the same number of individuals). If the data is unbalanced, the balancing procedure
#' comes into effect, and identifies the optimal number of items and groups to delete based on the punishment coefficient. If the amount of data deleted is more than defined by threshold
#' (balancing.limit) then results should be interpreted with caution.
#' 
#' The \code{summary()} function produces output similar to:
#' 
#' \preformatted{
#' Summary of Coefficients:
#'                     Estimate Std. Error Lower CI (99%) Upper CI (99%)   Z value   Pr(>|z|) Significance
#' beta_b             0.4279681  0.7544766     -1.5154349       2.371371 0.5672384 0.57055223
#' gamma_Petal.Length 0.4679522  0.2582579     -0.1972762       1.133181 1.8119567 0.06999289            .
#'
#' For comparison, summary of coefficients from unoptimized analysis (ML):
#'                    Estimate   Std. Error Lower CI (99%) Upper CI (99%)      Z value   Pr(>|z|) Significance
#' beta_b             0.6027440 5.424780e+15  -1.397331e+16   1.397331e+16 1.111094e-16 1.00000000
#' gamma_Petal.Length 0.4679522 2.582579e-01  -1.972762e-01   1.133181e+00 1.811957e+00 0.06999289            .
#'
#' Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#' }
#'
#' @return A list containing the results of the regularized Bayesian estimation, which includes the model formula,dependent and context variables,and other relevant details from the analysis.
#
#' @author 
#' Valerii Dashuk \email{vadashuk@gmail.com},
#' Binayak Timilsina \email{binayak.timilsina001@gmail.com},
#' Martin Hecht, and
#' Steffen Zitzmann
#' 
#' @references
#' Dashuk, V., Hecht, M., Luedtke, O., Robitzsch, A., & Zitzmann, S. (2024). \doi{10.13140/RG.2.2.18148.39048}
#' 
#' Dashuk, V., Hecht, M., Luedtke, O., Robitzsch, A., & Zitzmann, S. (2024). \doi{10.13140/RG.2.2.34350.01604}
#' 
#' Luedtke, O., Marsh, H. W., Robitzsch, A., Trautwein, U., Asparouhov, T., & Muthen, B. (2008). \doi{10.1037/a0012869}
#'
#' @examples
#' 
#' # Example 1: usage with the iris dataset
#'
#' result_iris <- mlob(
#' Sepal.Length ~ Sepal.Width + Petal.Length, 
#' data = iris, group = 'Species',
#' conf.level = 0.01,
#' jackknife = FALSE)
#' 
#' # View summary statistics (similar to summary of a linear model);
#' 
#' summary(result_iris)
#' 
#' # Example 2: usage with slightly unbalanced ChickWeight dataset
#' 
#' result_ChickWeight <- mlob(
#' weight ~ Diet, 
#' data = ChickWeight, 
#' group = 'Time', 
#' punish.coeff = 1.5, 
#' jackknife = FALSE)
#' 
#' # View the results
#' 
#' print(result_ChickWeight)
#' 
#' # View summary statistics
#' 
#' summary(result_ChickWeight)
#' 
#' # Example 3: usage with highly unbalanced mtcars dataset (adjusted balancing.limit)
#' 
#' result_mtcars <- mlob(
#' mpg ~ hp + wt + am + hp:wt + hp:am, 
#' data = mtcars, group = 'cyl', 
#' balancing.limit = 0.35)
#' 
#' # View summary statistics
#' 
#' summary(result_mtcars)
#' 
#' @export
mlob <- function(formula, data, group, balancing.limit=0.2, conf.level = 0.05, jackknife = FALSE, punish.coeff = 2, ...) {

  # Save the name of dataframe for displaying the output
  data_name  <- deparse(substitute(data))
  
  # Ensure data is a data frame
  if (!is.data.frame(data)) {
    warning("The 'data' argument is not a data frame. Converting to data frame.")
    data <- as.data.frame(data)
  }

  # Check if all columns in the data are numeric
  #if (!all(sapply(data, is.numeric))) {
   # stop("All columns in the 'data' must be numeric.")
  #}

  # Ensure the 'conf.level' is numeric and between 0 and 1
  if (!is.numeric(conf.level) || length(conf.level) != 1 || conf.level <= 0 || conf.level >= 1) {
    stop("The 'conf.level' argument must be a numeric scalar with value between 0 and 1 (exclusive).")
  }

  # Ensure 'jackknife' is a logical (TRUE or FALSE)
  if (!is.logical(jackknife) || length(jackknife) != 1) {
    stop("The 'jackknife' argument must be a scalar Boolean variable, i.e. TRUE or FALSE.")

  }

  #Ensure balancing.limit is a numeric single value between 0 and 1
  if (!is.numeric(balancing.limit) || length(balancing.limit) != 1 || balancing.limit < 0 || balancing.limit > 1) {
    stop("The 'balancing.limit' argument must be a numeric scalar with value between 0 and 1 (inclusive).")
  }

  #Ensure the punish.coeff is numeric single value greater than 0.
  if (!is.numeric(punish.coeff) || length(punish.coeff) != 1 || punish.coeff <= 0) {
    stop("The 'punish.coeff' should be a numeric value greater than 0.")
  }

  
  if (is.factor(data[[group]])) {
    data[[group]] = as.numeric(data[[group]])
  }
  
  # Parse the formula (Y ~ X + C...)
  # to incorporate model.matrix in mlob function,we can modify formula parsing part to automatically handle interaction terms,factors,
  # we parse covariates from the Matrix A not from the data

  all_vars <- all.vars(formula)
  response_var <- all_vars[1]

  # Create the model matrix using the formula
  A <- stats::model.matrix(formula, data = data)
  A <- as.data.frame(A)
  

  # Check if intercept exists, rename it to 'y' or create a new column response_var if not present
  if ("(Intercept)" %in% colnames(A)) {
    colnames(A)[colnames(A) == "(Intercept)"] <- response_var
    A[[response_var]] = data[[response_var]]
  } else {
    A[[response_var]] = data[[response_var]] # Creates a new column response_var
    # Move 'response_var' to the first position
    A <- A[, c(response_var, setdiff(colnames(A), response_var))]
  }


  # Get all predictors from the formula (i.e., X + control variables)
  # All columns from the model matrix

  predictors <- colnames(A)
  

  # Exclude the intercept (if included in the model matrix)
  predictors <- predictors[predictors != response_var]

  # Extract response variable from data(dependent variable)
  response <- data[[response_var]]


  # Ensure the response variable exists in the data
  if (is.null(response)) {
    stop(paste("The response variable", response_var, "is not found in the data."))
  }

  # Now we have `model_matrix` and `response`

  # Extract the group variable from the data using the specified group name
  group_var <- data[[group]]

  # Add the group variable to the model matrix
  A$group <- group_var

  # Check group sizes
  group_counts <- table(A$group)
  # cat("Groups and their sizes:")
  # print(group_counts)

  # Number of groups
  group_num = length(group_counts)


  # Check for balancing
  if (length(unique(A$group)) < 2) {
    stop("Not enough groups to balance the data.")
  }

  # Initialize variables for deletion logic
  tab <- as.numeric(group_counts)  # Convert table to numeric vector
  s <- rep(0, max(tab) - min(tab) + 1)  # Prepare a vector for items to delete
  s2 <- rep(0, max(tab) - min(tab) + 1)
  s3 <- rep(0, max(tab) - min(tab) + 1)

  # Function to detect unbalanced data
  check_balance <- function(tab) {
    # Check if all elements in `tab` are the same
    length(unique(tab)) == 1
  }


  # Balancing procedure for unbalanced data

  # Calculate items to delete for each group size
  for (i in min(tab):max(tab)) {
    delta <- tab - i
    delta2 <- tab - i
    delta3 <- rep(0,length(delta))
    for (j in 1:length(tab)) {
      if (delta[j] < 0) { # this means we need to delete the whole group
        delta[j] = punish.coeff*tab[j]  # with punishment
        delta2[j] = tab[j] # without punishment
        delta3[j] = 1
      }
    }
    s[i - min(tab) + 1] <- sum(delta) # Total items to delete for this group size with punishment
    s2[i - min(tab) + 1] <- sum(delta2) # without punishment
    s3[i - min(tab) + 1] <- sum(delta3) # total number of group we delete
  }

  # Create a summary data frame for group sizes and items to delete
  S <- data.frame(
    Group_size = c(min(tab):max(tab)),
    Items_to_delete_with_punishment = s,
    Items_to_delete_without_punishment = s2,
    Groups_to_delete= s3
  )

  # View the data frame in R
  # View(S)

  # Display the minimum items to delete, if the data is unbalanced
  # Find the minimum entries to delete with punishment, and for this case show the real number of entries to delete (without punishment)
  if (check_balance(tab)==FALSE){

    min_items_to_delete_without_punishment <- S$Items_to_delete_without_punishment[which.min(S$Items_to_delete_with_punishment)]
    # cat("\nOptimal number of entries to delete for balance:", min_items_to_delete_without_punishment, "\n")

    min_groups_to_delete <- S$Groups_to_delete [which.min(S$Items_to_delete_with_punishment)]
    # cat("\nOptimal number of groups to delete for balance:", min_groups_to_delete, "\n")
  }

  # If data is unbalanced issue a warning with number of entries, groups to delete
  if (check_balance(tab)==FALSE) {
    percentage_items_deleted <- 100*min_items_to_delete_without_punishment / nrow(A)
    
    warning(sprintf(" Your data is unbalanced. Balancing procedure was used and %.1f%% of data was deleted.\n   Deleted entries: %d\n   Deleted groups: %d", percentage_items_deleted, min_items_to_delete_without_punishment, min_groups_to_delete))


    if (balancing.limit != 0.2)
    {
      warning(" You changed the balancing limit to fine-tune the balancing.\n   Increasing the balancing limit might result in the loss of data connections.\n")
    }
  }

  # Check if the ratio of minimum items to delete is within the balancing.limit
  if (check_balance(tab)==FALSE) {
    if (min_items_to_delete_without_punishment / nrow(A) > balancing.limit || min_groups_to_delete / group_num > balancing.limit ) {
      
      balancing.percentage <-balancing.limit*100 # to display balancing limit in percents
      
      stop(sprintf("The share of data that should be deleted is %.1f%% and exceeds the balancing.limit, %s%%.\n\n  Data may not be balanced. If you want to run, adjust the balancing.limit. \n\n", percentage_items_deleted, format(balancing.percentage, trim = TRUE)))

      return(A)  # Optionally return the unbalanced matrix
    }

    # Balance the model matrix A based on optimal size (number of groups and elements in the group)
    target_group_size <- S$Group_size[which.min(S$Items_to_delete_with_punishment)]

    for (group_name in names(group_counts)) {
      group_indices <- which(A$group == group_name)

      if (length(group_indices) > target_group_size) {
        # If the group size is larger than the target, remove excess elements
        to_remove <- length(group_indices) - target_group_size
        group_indices_to_remove <- sample(group_indices, to_remove)
        A <- A[-group_indices_to_remove, ]

      } else if (length(group_indices) < target_group_size) {
        # If the group size is smaller than the target, remove the entire group
        group_indices_to_remove <- group_indices
        A <- A[-group_indices_to_remove, ]
      }
    }
    
      #  delete groups that are unused for the case when groups are given by sequential factor type of the data
      
      if (is.factor(A$group)) {
        A$group <- droplevels(A$group)
      }
    
    
      # Print final counts after balancing
      balanced_group_counts <- table(A$group)
      
      
  #    cat("\nGroups and their sizes after balancing:\n")
  #    print(balanced_group_counts)

    # The balanced model matrix and group sizes for further use
    #return(A)
    group_num <- length(balanced_group_counts)
    group_size <- unique(balanced_group_counts)
    if (length(group_size) > 1){
      warning("Data was not balanced correctly.")
    }
    
  } else {
    group_num <- length(group_counts)
    group_size <- unique(group_counts)
  }
  

  if(FALSE){
   #Use na.omit to remove rows with NA values in relevant columns
    relevant_vars <- c(response_var, predictor_var, control_vars)
    complete_data <- na.omit(data[relevant_vars]) # complete_data <- data[complete.cases(data[relevant_vars]), ]
  
    #Recalculate number of rows
    n_rows <- nrow(complete_data)
  
    #Ensure we have enough data points after removing NAs
    if (n_rows < 1) {
      stop("No complete cases available after removing missing values.")
    }
  
  }

  # changing the model matrix into a data frame called data
  data <-as.data.frame(A)

  predictor_var <- names(data)[2]
  
  # check that there are control variables in the data
  if (ncol(data) > 3) {
    control_vars <- names(data)[3:(length(names(data)) - 1)] # do not include last column - it defines groups
  }
  # Extract relevant data columns
  y <- data[[response_var]]
  x <- data[[predictor_var]]
  # if control variables exist
  if (exists("control_vars")) {
    C <- data[control_vars]
    C <- as.matrix(C)
  }
  # Ensure that all variables (response, predictor, and controls) are numeric
  # Check response variable
  if (!is.numeric(y)) {
    stop(paste("The response variable", response_var, "must be numeric."))
  }

  # Check predictor variable
  if (!is.numeric(x)) {
    stop(paste("The predictor variable", predictor_var, "must be numeric."))
  }

  # Check all control variables (if any exist)
  if (exists("control_vars")) {
    nonnum_values <- sapply(C, function(x) !is.numeric(x))

    if (any(nonnum_values)) {
      stop(paste("The following control variables contain non-numeric values:",
                 paste(colnames(C[nonnum_values]), collapse = ", ")))
    }
  }

  if (sum(is.na(y))>0) {
    stop(paste("The response variable", response_var," has null (nan) values."))
  }

  # Check predictor variable
  if (sum(is.na(x))>0) {
    stop(paste("The predictor variable", predictor_var, "has null (nan) values."))
  }
  exists("control_vars")
  # Checking for missing values in control variables
  if (exists("control_vars")) {
    missing_values <- sapply(C, anyNA)

    if (any(missing_values)) {
      stop(paste("The following control variables contain missing values:",
                 paste(colnames(C[missing_values]), collapse = ", ")))
    }
  }

  # If all checks pass, proceed with the function...
  # message("All checks passed. Function is ready to proceed.")


  # Calculate group size (n) and total observations (kn)
  n <- group_size  # Size of each group
  k <- group_num # Number of groups
  kn <- k * n  # Total observations
  
  # Prepare data_CV list
  data_CV <- list(
    y = y,
    x = x,
    k = k,
    n = n,
    kn = kn
  )
  
  # Check if C exists
  if (exists("control_vars")) {
    data_CV$C <- C  # Add C to the list
    data_CV$kc <- length(control_vars)  # Add number of control variables
  }

  
  ML <- estimate_ML_CV(data_CV) # run ML estimator and get a preliminary estimation of b_b for estimate_Bay_CV
  #check if there is any variation between group
  if (ML$tau_x2==0 || ML$tau_yx==0){
    stop(sprintf("The provided data does not include between-group variation for the given model.\n  A two-level model is no longer necessary because the grouping structure does not add information.\n  Consider using a single-level model, treating all observations as independent."))
  }
  #check if there is any variation within group
  if (ML$sigma_x2==0 || ML$sigma_yx==0){
    stop(sprintf("The provided data does not include within-group variation for the given model.\n  All observations within a group are identical (or perfectly correlated).\n  Consider using a simpler model."))
  }
  
  
  data_CV$b_b = ML$beta_b_ML_CV # dummy real value of beta_b
  
  # Call estimate_Bay_CV function with data_CV
  Bay <- estimate_Bay_CV(data_CV)
  
  # recalculate SE of Bayesian estimator with jackknife if a
  if (jackknife == TRUE){
    Bay_jackknife    <- estimate_Bay_CV_SE_jackknife_individual(data_CV)
    Bay$SE_beta_Bay  <- Bay_jackknife$SE_beta_Bay_ML_jackknife_individual
    # Bay$SE_beta_ML   <- Bay_jackknife$SE_beta_ML_jackknife_individual # open in case jackknife for ML needed
    # If there were any gamma‐covariates, copy over their SEs
    if (!is.null(Bay_jackknife$SE_gamma_jackknife_individual)) {
      Bay$SE_gamma   <- Bay_jackknife$SE_gamma_jackknife_individual
    }
    
  }

  # Generate the result output

  # Number of control variables (kc)
  if (exists("control_vars")){
    kc <- data_CV$kc
  } else {
    kc<-0
  }
  
  
  # If control variables are present
  if (kc>0) {
    # Create the list of estimated values dynamically
    Coefficients <- data.frame(
      beta_b = Bay$beta_b_Bay,
      t(sapply(1:kc, function(i) Bay$gamma[i]))   # Dynamic number of gamma columns
    )
    
    colnames(Coefficients) <- c("beta_b", paste0("gamma_", control_vars))  # Adjust gamma column names
  } else {
    Coefficients <- data.frame(
      beta_b = Bay$beta_b_Bay
    )
    colnames(Coefficients) <- c("beta_b")  # Adjust column names
  }
  
  
  if (kc>0) {
    Standard_Error <- data.frame(
      beta_b = Bay$SE_beta_Bay,
      t(sapply(1:kc, function(i) Bay$SE_gamma[i]))
    )
  
    colnames(Standard_Error) <- c("beta_b", paste0("gamma_", control_vars))
  } else {
    Standard_Error <- data.frame(
      beta_b = Bay$SE_beta_Bay
    )
    
    colnames(Standard_Error) <- c("beta_b")
  }
    
  
  if (kc>0) {
    Confidence_Interval <- data.frame(
      Lower = c(Bay$beta_b_Bay - stats::qnorm(1-conf.level/2) * Bay$SE_beta_Bay, Bay$gamma - stats::qnorm(1-conf.level/2) * Bay$SE_gamma),
      Upper = c(Bay$beta_b_Bay + stats::qnorm(1-conf.level/2) * Bay$SE_beta_Bay, Bay$gamma + stats::qnorm(1-conf.level/2) * Bay$SE_gamma)
    )
    
    rownames(Confidence_Interval) <- c("beta_b", paste0("gamma_", control_vars))
  } else {
    Confidence_Interval <- data.frame(
      Lower = c(Bay$beta_b_Bay - stats::qnorm(1-conf.level/2) * Bay$SE_beta_Bay),
      Upper = c(Bay$beta_b_Bay + stats::qnorm(1-conf.level/2) * Bay$SE_beta_Bay)
    )
    
    rownames(Confidence_Interval) <- c("beta_b")
  }

  Confidence_level <- paste0((1 - conf.level) * 100, "%")
  
  if (kc>0) {
    Z_value <- data.frame(
      beta_b = Bay$beta_b_Bay / Bay$SE_beta_Bay,
      t(sapply(1:kc, function(i) Bay$gamma[i] / Bay$SE_gamma[i]))
    )
  
    colnames(Z_value) <- c("beta_b", paste0("gamma_", control_vars))
  } else {
    Z_value <- data.frame(
      beta_b = Bay$beta_b_Bay / Bay$SE_beta_Bay
    )
    
    colnames(Z_value) <- c("beta_b")
  }
  
  if (kc>0) {
    p_value <- data.frame(
      beta_b = 2 * (1 - stats::pnorm(abs(Bay$beta_b_Bay / Bay$SE_beta_Bay))),
      t(sapply(1:kc, function(i) 2 * (1 - stats::pnorm(abs(Bay$gamma[i] / Bay$SE_gamma[i])))))
    )
  
    colnames(p_value) <- c("beta_b", paste0("gamma_", control_vars))
  } else {
    p_value <- data.frame(
      beta_b = 2 * (1 - stats::pnorm(abs(Bay$beta_b_Bay / Bay$SE_beta_Bay)))
    )
    
    colnames(p_value) <- c("beta_b")
  }

  # Create the dynamic call_info string
  call_info <- paste0("mlob(", deparse(formula), ", data = ", data_name, ", group = ", group)

  # Add balancing.limit if it is not default value
  if (!missing(balancing.limit) && balancing.limit != 0.2) {
    call_info <- paste0(call_info, ", balancing.limit = ", balancing.limit)
  }

  # Conditionally add `conf.level` if it's not the default value
  if (!missing(conf.level) && conf.level != 0.05) {
    call_info <- paste0(call_info, ", conf.level = ", conf.level)
  }

  # Conditionally add `jackknife` if it's not TRUE
  if (!missing(jackknife) && !jackknife) {
    call_info <- paste0(call_info, ", jackknife = ", jackknife)
  }

  call_info <- paste0(call_info, ")")

  result <- list(
    Coefficients = Coefficients,
    Standard_Error = Standard_Error,
    Confidence_Interval = Confidence_Interval,
    Confidence_level = Confidence_level,
    Z_value = Z_value,
    p_value = p_value,
    call_info = call_info

  )
  
  class(result) <- "mlob_result"  # Assign custom class
  
  # return(result)

  # For unoptimized estimator  estimate_ML
  
  if (kc>0) {
    Coefficients_ML <- data.frame(
      beta_b = Bay$beta_b_ML,  # Using beta_b_ML here
      t(sapply(1:kc, function(i) Bay$gamma[i]))   # Dynamic number of gamma columns
    )
    colnames(Coefficients_ML) <- c("beta_b_ML", paste0("gamma_", control_vars))
  
    Standard_Error_ML <- data.frame(
      beta_b = Bay$SE_beta_ML,
      t(sapply(1:kc, function(i) Bay$SE_gamma[i]))
    )
    colnames(Standard_Error_ML) <- c("beta_b_ML", paste0("gamma_", control_vars))
  
    Confidence_Interval_ML <- data.frame(
      Lower = c(Bay$beta_b_ML - stats::qnorm(1 - conf.level / 2) * Bay$SE_beta_ML, Bay$gamma - stats::qnorm(1 - conf.level / 2) * Bay$SE_gamma),
      Upper = c(Bay$beta_b_ML + stats::qnorm(1 - conf.level / 2) * Bay$SE_beta_ML, Bay$gamma + stats::qnorm(1 - conf.level / 2) * Bay$SE_gamma)
    )
    rownames(Confidence_Interval_ML) <- c("beta_b_ML", paste0("gamma_", control_vars))
  
    Confidence_level_ML <- paste0((1 - conf.level) * 100, "%")
  
    Z_value_ML <- data.frame(
      beta_b = Bay$beta_b_ML / Bay$SE_beta_ML,
      t(sapply(1:kc, function(i) Bay$gamma[i] / Bay$SE_gamma[i]))
    )
    colnames(Z_value_ML) <- c("beta_b_ML", paste0("gamma_", control_vars))
  
    p_value_ML <- data.frame(
      beta_b = 2 * (1 - stats::pnorm(abs(Bay$beta_b_ML / Bay$SE_beta_ML))),
      t(sapply(1:kc, function(i) 2 * (1 - stats::pnorm(abs(Bay$gamma[i] / Bay$SE_gamma[i])))))
    )
    colnames(p_value_ML) <- c("beta_b", paste0("gamma_", control_vars))
    
  } else {
    
    Coefficients_ML <- data.frame(
      beta_b = Bay$beta_b_ML  # Using beta_b_ML here
    )
    
    colnames(Coefficients_ML) <- c("beta_b_ML")
    
    Standard_Error_ML <- data.frame(
      beta_b = Bay$SE_beta_ML
    )
    
    colnames(Standard_Error_ML) <- c("beta_b_ML")
    
    Confidence_Interval_ML <- data.frame(
      Lower = c(Bay$beta_b_ML - stats::qnorm(1 - conf.level / 2) * Bay$SE_beta_ML),
      Upper = c(Bay$beta_b_ML + stats::qnorm(1 - conf.level / 2) * Bay$SE_beta_ML)
    )
    
    rownames(Confidence_Interval_ML) <- c("beta_b_ML")
    
    Confidence_level_ML <- paste0((1 - conf.level) * 100, "%")
    
    Z_value_ML <- data.frame(
      beta_b = Bay$beta_b_ML / Bay$SE_beta_ML
    )
    
    colnames(Z_value_ML) <- c("beta_b_ML")
    
    p_value_ML <- data.frame(
      beta_b = 2 * (1 - stats::pnorm(abs(Bay$beta_b_ML / Bay$SE_beta_ML)))
    )
    
    colnames(p_value_ML) <- c("beta_b")
  }

  result$Coefficients_ML = Coefficients_ML
  result$Standard_Error_ML = Standard_Error_ML
  result$Confidence_Interval_ML = Confidence_Interval_ML
  result$Confidence_level_ML = Confidence_level_ML
  result$Z_value_ML = Z_value_ML
  result$p_value_ML = p_value_ML

  return(result)

}


#' @exportS3Method
print.mlob_result <- function(x, ...) {
  # Custom print method for mlob_result
  # Extract call information
  cat("Call:\n", x$call_info, "\n\n")

  # Print each section without row names
  cat("Coefficients\n")
  print(x$Coefficients, row.names = FALSE)

  cat("\nStandard_Error\n")
  print(x$Standard_Error, row.names = FALSE)

  cat("\nConfidence_Interval (", x$Confidence_level, ")\n", sep = "")
  print(x$Confidence_Interval)

  cat("\nZ value\n")
  print(x$Z_value, row.names = FALSE)

  cat("\np value\n")
  print(x$p_value, row.names = FALSE)
}


#' @exportS3Method
summary.mlob_result <- function(object, ...) {
  # Custom print method for summary_mlob_result
  # Function to assign stars based on p-values
  signif_stars <- function(pval) {
    if (pval < 0.001) {
      return("***")
    } else if (pval < 0.01) {
      return("**")
    } else if (pval < 0.05) {
      return("*")
    } else if (pval < 0.1) {
      return(".")
    } else {
      return(" ")
    }
  }


  # Prepare the summary table with significance stars
  summary_table <- data.frame(
    Estimate = as.numeric(object$Coefficients),
    `Std. Error` = as.numeric(object$Standard_Error),
    `Lower CI` = as.numeric(object$Confidence_Interval$Lower),
    `Upper CI` = as.numeric(object$Confidence_Interval$Upper),
    `Z value` = as.numeric(object$Z_value),
    `Pr(>|z|)` = as.numeric(object$p_value),
    Significance = sapply(object$p_value, signif_stars)
  )

  # If control variables are present in the model
  if (ncol(object$Coefficients)>1){
  # Extract control variable names from the Coefficients (removing the 'beta_b' column)
  control_vars <- colnames(object$Coefficients)[-1]  # Get control variable names (gamma_Age, gamma_Education, etc.)

  # Update the row names for the summary table using the actual control variable names
  rownames(summary_table) <- c("beta_b", control_vars)
  } else {
    # Update the row names for the summary table using the actual control variable names
    rownames(summary_table) <- c("beta_b")
  }
  
  # Redact column names for consistency
  colnames(summary_table) <- c("Estimate", "Std. Error", paste0("Lower CI (", object$Confidence_level, ")"),
                               paste0("Upper CI (", object$Confidence_level, ")"), "Z value", "Pr(>|z|)", "Significance")


  # Print the summary header
  cat("Call:\n", object$call_info, "\n\n")
  cat("Summary of Coefficients:\n")

  print(summary_table, row.names = TRUE)


  invisible(summary_table)  # Return the table invisibly so it doesn't clutter the output


  # Prepare the summary table For unoptimized estimator ML with significance stars
  summary_table_ML<- data.frame(
    Estimate = as.numeric(object$Coefficients_ML),
    `Std. Error` = as.numeric(object$Standard_Error_ML),
    `Lower CI` = as.numeric(object$Confidence_Interval_ML$Lower),
    `Upper CI` = as.numeric(object$Confidence_Interval_ML$Upper),
    `Z value` = as.numeric(object$Z_value_ML),
    `Pr(>|z|)` = as.numeric(object$p_value_ML),
    Significance = sapply(object$p_value_ML, signif_stars)
  )

  ## Extract control variable names from the Coefficients (removing the 'beta_b' column)
  # control_vars <- colnames(object$Coefficients)[-1]  # Get control variable names (gamma_Age, gamma_Education, etc.)
  
  # If control variables are present in the model
  if (ncol(object$Coefficients)>1){
    # Define the number of rows based on the number of coefficients (beta_b + gamma terms)
    rownames(summary_table_ML) <- c("beta_b", control_vars)
  } else {
    rownames(summary_table_ML) <- c("beta_b")
  }
  
  # Redact column names for consistency
  colnames(summary_table_ML) <- c("Estimate", "Std. Error", paste0("Lower CI (", object$Confidence_level, ")"),
                               paste0("Upper CI (", object$Confidence_level, ")"), "Z value", "Pr(>|z|)", "Significance")


  # Print the ML part of summary

  cat("\n\nFor comparison, summary of coefficients from unoptimized analysis (ML):\n")

  print(summary_table_ML, row.names = TRUE)

  # Add significance codes
  cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  
  
  # Extract first SE of ML and of regularized Bayes (or NA if missing)
  se_ml  <- if (!is.null(object$Standard_Error_ML) &&
                length(object$Standard_Error_ML) >= 1) object$Standard_Error_ML[[1]] else NA_real_
  se_reg <- if (!is.null(object$Standard_Error) &&
                length(object$Standard_Error) >= 1) object$Standard_Error[[1]] else NA_real_
  
  # Compute ratio (will be NaN/Inf if se_reg is zero or NA)
  se_ratio <- se_ml / se_reg
  
  # Display the note header
  cat("\nNote:\n")
  
  # Three‐way branching without any stop()
  if (!is.finite(se_ratio)) {
    
    # Case 1: missing, infinite, zero, etc.
    cat(
      "  The standard errors cannot be compared because one or both are missing, zero, NA, NaN, or infinite.\n",
      sep = ""
    )
    
  } else if (se_ratio <= 1.1) {
    
    # Case 2: roughly equal
    cat(
      "  The standard error from unoptimized ML estimation is approximately the same as the standard error obtained through our optimization procedure,\n",
      "  meaning that both approaches yield similarly accurate estimates.\n",
      sep = ""
    )
    
  } else {
    
    # Case 3: ML SE substantially larger
    cat(
      sprintf(
        "  The standard error from unoptimized ML estimation is about %.1f%% larger than the standard error obtained through our optimization procedure,\n",
        (se_ratio - 1) * 100
      ),
      "  meaning that the optimized estimates are more accurate.\n",
      "  Concerning the estimates themselves, the unoptimized ML estimates may differ greatly from the optimized estimates and should not be reported.\n",
      "  As the optimized estimates are always at least as accurate as the unoptimized ML estimates,\n",
      "  please use them and their corresponding standard errors (first table of output) for interpretation and reporting.\n",
      "  For more information, see Dashuk et al. (2025).\n",
      sep = ""
    )
  }
}


estimate_ML_CV <- function(data) {
   
  
  data$j <- data$k
  # Reshape the data.x into an n x j matrix
  x <- matrix(data$x, nrow = data$n, ncol = data$k)

  av_x <- mean(x)
  av_x_j <- colMeans(x)

  SSA <- data$n * sum(av_x_j^2) - data$n * data$k * (av_x^2)
  SSD <- sum(data$x * data$x) - data$n * sum(av_x_j * av_x_j)


  MSA <- SSA / (data$j - 1)
  MSD <- SSD / ((data$n - 1) * data$j)

  tau_x2 <- (MSA - MSD) / data$n
  sigma_x2 <- MSD

  # Reshape the data.y into an n x j matrix
  y <- matrix(data$y, nrow = data$n, ncol = data$j)

  av_y <- mean(y)
  av_y_j <- colMeans(y)

  SPA <- data$n * sum(av_y_j * av_x_j) - data$n * data$j * (av_x * av_y)
  SPD <- sum(data$x *data$y) - data$n * sum(av_y_j * av_x_j)

  MPA <- SPA / (data$j - 1)
  MPD <- SPD / ((data$n - 1) * data$j)

  tau_yx <- (MPA - MPD) / data$n
  sigma_yx <- MPD

  # ML estimation of beta_b
  beta_b_ML <- tau_yx / tau_x2

  # Multiple control variables part
  
  if (!is.null(data$kc) && data$kc != 0) {
    
    beta_b_ML_tilde <- NULL
    data$C = matrix(data$C)
    
    phi <- matrix(0, nrow = 3, ncol = data$kc)
    data$C = matrix(data$C, nrow = data$k*data$n, ncol = data$kc)
    for (i in 1:data$kc) {

      if (data$kc>1){
        # Reshape data$C[,i] into a matrix with n rows and J columns
        C_i_C <- matrix(data$C[, i], nrow = data$n, ncol = data$k)
      }else{
        C_i_C <- matrix(data$C, nrow = data$n, ncol = data$k)
      }

      av_C <- mean(C_i_C)
      av_C_j <- colMeans(C_i_C)

      SPA_C <- data$n * sum(av_C_j * av_x_j) - data$n * data$j* (av_x * av_C)
      SPD_C <- sum(x * data$C[, i]) - data$n * sum(av_C_j * av_x_j)

      MPA_C <- SPA_C / (data$j - 1)
      MPD_C <- SPD_C / ((data$n - 1) * data$j)
      
      
      tau_C_x <- (MPA_C - MPD_C) / data$n
      sigma_C_x <- MPD_C
      
      
      phi[2, i] <- tau_C_x / tau_x2
      phi[3, i] <- sigma_C_x / sigma_x2
    }

    C <- data$C
    x_b <- rep(av_x_j, each = data$n)
    #x_w <- as.vector(t(x)) - x_b
    x_w <- as.vector(data$x) - x_b
    C_tilde <- C - cbind(x_b, x_w) %*% (phi[2:3, ])

    C_tilde_mat <- cbind(1, C_tilde)

    gamma_all <- solve(t(C_tilde_mat) %*% C_tilde_mat) %*% t(C_tilde_mat) %*% data$y
    gamma <- gamma_all[-1]

    y_tilde <- data$y - C %*% gamma
    y_tilde <- matrix(y_tilde, nrow = data$n, ncol = data$j)

    av_y_tilde <- mean(y_tilde)
    av_y_tilde_j <- colMeans(y_tilde)

    SPA_tilde <- data$n * sum(av_y_tilde_j * av_x_j) - data$n * data$j * (av_x * av_y_tilde)
    #SPD_tilde <- sum(x * as.vector(t(y_tilde))) - data$n * sum(av_y_tilde_j * av_x_j)
    SPD_tilde <- sum(as.vector(data$x) * as.vector(y_tilde)) - data$n * sum(av_y_tilde_j * av_x_j)

    MPA_tilde <- SPA_tilde / (data$j - 1)
    MPD_tilde <- SPD_tilde / ((data$n - 1) * data$j)

    tau_y_tilde_x <- (MPA_tilde - MPD_tilde) / data$n
    sigma_y_tilde_x <- MPD_tilde

    # ML estimation of beta_b for the model with control variables
    beta_b_ML_tilde <- tau_y_tilde_x / tau_x2
    
    return(list(beta_b_ML_CV = beta_b_ML_tilde, beta_b_ML_no_CV = beta_b_ML, gamma = gamma, tau_x2 = tau_x2, sigma_x2 = sigma_x2, tau_yx = tau_yx, sigma_yx = sigma_yx))
  
  } else {
    
    return(list(beta_b_ML_CV = beta_b_ML, beta_b_ML_no_CV = beta_b_ML, tau_x2 = tau_x2, sigma_x2 = sigma_x2, tau_yx = tau_yx, sigma_yx = sigma_yx))
  }

  
}


estimate_Bay_CV <- 
  function(data) {

  J <- data$k
  n <- data$n
  

  a <- rep(0, n * J + J + 1)

  for (i in 1:(n * J)) {
    a[i] <- -1 / (n * (n - 1) * J)
  }

  for (i in (n * J + 1):(n * J + J)) {
    a[i] <- (n * J - 1) / ((n - 1) * (J - 1) * J)
  }

  a[n * J + J + 1] <- -J / (J - 1)

  A <- diag(a)

  # B1 = [A, zeros(n*J+J+1); zeros(n*J+J+1,2*(n*J+J+1))];

  # B2 = [zeros(n*J+J+1), A/2; A/2, zeros(n*J+J+1)];

  # calculate covariance matrices:

  #cov(x):

  x <- matrix(data$x, n, J)

  av_x <- mean(x)
  av_x_j <- colMeans(x)

  SSA <- n * sum(av_x_j^2) - n * J * (av_x^2)
  SSD <- sum(data$x^2) - n * sum(av_x_j^2)

  MSA <- SSA / (J - 1)
  MSD <- SSD / ((n - 1) * J)

  tau_x2 <- (MSA - MSD) / n
  sigma_x2 <- MSD

  #Multiple control variables part
  # checks the existence of control variables and their number is not 0

  #Find 3 phi's for every of kc control variables from equation
  # C = phi1 + phi2*X_b + phi3*X_w + eps. Our aim is phi2 and phi3

  if ("kc" %in% names(data) && !is.null(data$kc) && data$kc != 0) {
    phi <- matrix(0, 3, data$kc)

    C <- list()

    for (i in 1:data$kc) {
        C[[i]] <- list()
  
        # 2 cases for number of control variables
        if (data$kc>1){
          # Reshape data$C[,i] into a matrix with n rows and J columns
          C[[i]]$C <- matrix(data$C[, i], nrow = n, ncol = J)
        } else{
          C[[i]]$C <- matrix(data$C, nrow = n, ncol = J)
        }
  
        # Compute the average of the entire C matrix
        C[[i]]$av_C <- mean(C[[i]]$C)
  
        # Compute the average of each column in the C matrix
        C[[i]]$av_C_j <- colMeans(C[[i]]$C)
  
        # Compute SPA_C and SPD_C
        C[[i]]$SPA_C <- n * sum(C[[i]]$av_C_j * av_x_j) - n * J * (av_x * C[[i]]$av_C)
        # 2 cases for number of control variables
        if (data$kc>1){
          C[[i]]$SPD_C <- sum(data$x * data$C[, i]) - n * sum(C[[i]]$av_C_j * av_x_j)
        }else{
          C[[i]]$SPD_C <- sum(data$x * data$C) - n * sum(C[[i]]$av_C_j * av_x_j)
        }
        # Compute MPA_C and MPD_C
        C[[i]]$MPA_C <- C[[i]]$SPA_C / (J - 1)
        C[[i]]$MPD_C <- C[[i]]$SPD_C / ((n - 1) * J)
  
        # Compute tau_C_x and sigma_C_x
        C[[i]]$tau_C_x <- (C[[i]]$MPA_C - C[[i]]$MPD_C) / n
        C[[i]]$sigma_C_x <- C[[i]]$MPD_C
  
        # Initialize phi as a vector of zeros
        C[[i]]$phi <- numeric(3)
  
        # Compute phi values
        C[[i]]$phi[2] <- C[[i]]$tau_C_x / tau_x2
        C[[i]]$phi[3] <- C[[i]]$sigma_C_x / sigma_x2
  
        # Store the phi vector in the phi matrix
        phi[, i] <- C[[i]]$phi
      }

    C <- data$C
    x_b <- matrix(rep(av_x_j, each = n), nrow = n * J, ncol = 1)
    #x_b <- matrix(rep(av_x_j, each = n), nrow = n, byrow = TRUE)
    x_w <- as.vector(matrix(data$x, nrow = n * J, ncol = 1)) - x_b
    C_tilde <- C - cbind(x_b, x_w) %*% phi[2:nrow(phi), ]
    C_tilde_mat <- cbind(rep(1, data$kn), C_tilde)

    #Estimate gamma from equation y = gamma1 + gamma*C_tilde + eps. Our
    #interest is vestor of parameters gamma without intersection gamma1
    gamma_all <- solve(t(C_tilde_mat) %*% C_tilde_mat) %*% t(C_tilde_mat) %*% data$y

    SE_gamma <- numeric(data$kc)

    for (i in 1:data$kc) {
      var_res_gamma <- sum((data$y - C_tilde_mat[, i + 1] * gamma_all[i + 1])^2) / (data$kn - 2)
      SE_gamma[i] <- sqrt(var_res_gamma / ((data$kn - 1) * stats::var(C_tilde_mat[, i + 1])))
    }

    gamma <- gamma_all[-1]
    # Find y_tilde - the difference between real and estimated
    # y (without intersection)
    if (data$kc>1){
      y_tilde <- data$y - C %*% gamma
    }else{
      y_tilde <- data$y - C * gamma
    }
    data$y <- y_tilde
  }
  #cov(x,y)
  y <- matrix(data$y, n, J)

  av_y <- mean(y)
  av_y_j <- colMeans(y)

  SPA <- n * sum(av_y_j * av_x_j) - n * J * (av_x * av_y)
  SPD <- sum(data$x * data$y) - n * sum(av_y_j * av_x_j)

  MPA <- SPA / (J - 1)
  MPD <- SPD / ((n - 1) * J)

  tau_yx <- (MPA - MPD) / n
  sigma_yx <- MPD

  #cov(y):
  SOA <- n * sum(av_y_j^2) - n * J * (av_y^2)
  SOD <- sum(data$y^2) - n * sum(av_y_j^2)

  MOA <- SOA / (J - 1)
  MOD <- SOD / ((n - 1) * J)

  tau_y2 <- (MOA - MOD) / n
  sigma_y2 <- MOD
  
  # truncate the estimated variances to avoid negatives
  
  #tau_x2 <- max(0, tau_x2)
  #tau_y2 <- max(0, tau_y2)
  
  #sigma_x2 <- max(0, sigma_x2)
  #sigma_y2 <- max(0, sigma_y2)
  
  # dummy to remember tau_x2 it was truncated -> no between-group variability
  if (tau_x2==0){
    truncated = 1
  }

  #create Sigma_x matrix of covariances of x_ij, mean(x_j), mean(x)

  # create Sigma blockvise with 9 blocks

  # block 1: cov(x_ij,x_ij):
  a1 <- diag(sigma_y2 * rep(1, n))
  #a1 <- diag(sigma_x2, n)

  a2 <- matrix(tau_x2, n, n)

  a <- a1 + a2
  block1 <- kronecker(diag(J), a)
  # Block 2: cov(y_ij, mean(y_j))

  block2 <- kronecker(diag(J), (tau_x2 + sigma_x2 / n) * matrix(1, n, 1))
  #block2 <- kronecker(diag(J), (tau_y2 + sigma_y2 / n) * rep(1, n))

  block3 <- matrix((tau_x2 / J + sigma_x2 / (n * J)), n * J, 1)
  #block3 <- (tau_y2 / J + sigma_y2 / (n * J)) * rep(1, n * J)
  block4 <- t(block2)

  block5 <- (tau_x2 + sigma_x2 / n) * diag(J)

  block6 <- (tau_x2 / J + sigma_x2 / (n * J)) * matrix(1, J, 1)
  #block6 <- (tau_y2 / J + sigma_y2 / (n * J)) * rep(1, J
  block7 <- t(block3)

  block8 <- t(block6)

  block9 <- tau_x2 / J + sigma_x2 / (n * J)

  Sigma_x <- rbind(cbind(block1, block2, block3), cbind(block4, block5, block6), cbind(block7, block8, block9))

  a1 <- diag(sigma_yx, n)
  a2 <- matrix(tau_yx, n, n)

  a <- a1 + a2
  block1 <- kronecker(diag(J), a)

  block2 <- kronecker(diag(J), (tau_yx + sigma_yx / n) * matrix(1, n, 1))

  block3 <- matrix((tau_yx / J + sigma_yx / (n * J)), n * J, 1)

  block4 <- t(block2)

  block5 <- (tau_yx + sigma_yx / n) * diag(J)

  block6 <- (tau_yx / J + sigma_yx / (n * J)) * matrix(1, J, 1)

  block7 <- t(block3)

  block8 <- t(block6)

  block9 <- tau_yx / J + sigma_yx / (n * J)

  Sigma_yx <- rbind(cbind(block1, block2, block3), cbind(block4, block5, block6), cbind(block7, block8, block9))

  a1 <- diag(sigma_y2, n)
  a2 <- matrix(tau_y2, n, n)

  a <- a1 + a2
  block1 <- kronecker(diag(J), a)

  block2 <- kronecker(diag(J), (tau_y2 + sigma_y2 / n) * matrix(1, n, 1))

  block3 <- matrix((tau_y2 / J + sigma_y2 / (n * J)), n * J, 1)

  block4 <- t(block2)

  block5 <- (tau_y2 + sigma_y2 / n) * diag(J)

  block6 <- (tau_y2 / J + sigma_y2 / (n * J)) * matrix(1, J, 1)

  block7 <- t(block3)

  block8 <- t(block6)

  block9 <- tau_y2 / J + sigma_y2 / (n * J)

  Sigma_y <- rbind(cbind(block1, block2, block3), cbind(block4, block5, block6), cbind(block7, block8, block9))

  #spectral decomposition of Sigma_x

  #compute eigenvalues matrix D_x of Sigma_x

  # Initialize D_x
  D_x <- rep(0, n * J + J + 1)

  # Set values for specific indices
  D_x[(J + 2):(n * J + 1)] <- sigma_x2
  D_x[(n * J + 2):(n * J + J)] <- (n + 1) * (tau_x2 + sigma_x2 / n)
  D_x[n * J + J + 1] <- (n * J + J + 1) * (tau_x2 + sigma_x2 / n) / J

  # Create diagonal matrix D_x
  D_x <- diag(D_x)
 
  

  # Initialize V_yx matrix
  V_yx <- matrix(0, n * J + J + 1, n * J + J + 1)

  # lambda_i = 0 ((J+1) pieces)
  for (i in 1:J) {
    V_yx[((n * (i - 1) + 1):(n * i)), i] <- -1 / sqrt(n * (n + 1))
    V_yx[(n * J + i), i] <- sqrt(n / (n + 1))
  }

  V_yx[1:(n * J + J), J + 1] <- -1 / sqrt((n * J + J) * (n * J + J + 1))
  V_yx[(n * J + J + 1), J + 1] <- sqrt((n * J + J) / (n * J + J + 1))

  # lambda_i = sigma_yx ((n-1)J pieces)
  AU <- matrix(0, n, n - 1)
  for (i in 1:(n - 1)) {
    AU[1:i, i] <- -1 / sqrt(i * (i + 1))
    AU[(i + 1), i] <- i / sqrt(i * (i + 1))
  }

  for (j in 1:J) {
    V_yx[(n * (j - 1) + 1):(n * j), (J + 2 + (n - 1) * (j - 1)):(J + 1 + (n - 1) * j)] <- AU
  }

  # lambda_i = (n + 1)(tau_yx + sigma_yx / n) ((J - 1) pieces)
  for (i in 1:(J - 1)) {
    V_yx[1:(i * n), (n * J + 1 + i)] <- -1 / sqrt(i * (i + 1) * (n + 1))
    V_yx[(n * J + 1):(n * J + i), (n * J + 1 + i)] <- -1 / sqrt(i * (i + 1) * (n + 1))
    V_yx[(i * n + 1):((i + 1) * n), (n * J + 1 + i)] <- i / sqrt(i * (i + 1) * (n + 1))
    V_yx[(n * J + i + 1), (n * J + 1 + i)] <- i / sqrt(i * (i + 1) * (n + 1))
  }

  # lambda_nJ + J + 1 = (n * J + J + 1)(tau_yx + sigma_yx / n) / J (1 piece)
  V_yx[1:(n * J + J + 1), (n * J + J + 1)] <- 1 / sqrt(n * J + J + 1)

  # Spectral decomposition of Sigma_y

  # Initialize D_y
  D_y <- rep(0, n * J + J + 1)

  # Set values for specific indices
  D_y[(J + 2):(n * J + 1)] <- sigma_y2
  D_y[(n * J + 2):(n * J + J)] <- (n + 1) * (tau_y2 + sigma_y2 / n)
  D_y[n * J + J + 1] <- (n * J + J + 1) * (tau_y2 + sigma_y2 / n) / J

  # Create diagonal matrix D_y
  D_y <- diag(D_y)

  # Initialize V_y matrix
  V_y <- matrix(0, n * J + J + 1, n * J + J + 1)

  # lambda_i = 0 ((J + 1) pieces)
  for (i in 1:J) {
    V_y[((n * (i - 1) + 1):(n * i)), i] <- -1 / sqrt(n * (n + 1))
    V_y[(n * J + i), i] <- sqrt(n / (n + 1))
  }

  V_y[1:(n * J + J), J + 1] <- -1 / sqrt((n * J + J) * (n * J + J + 1))
  V_y[(n * J + J + 1), J + 1] <- sqrt((n * J + J) / (n * J + J + 1))

  # lambda_i = sigma_yx ((n - 1)J pieces)
  for (i in 1:(n - 1)) {
    AU[1:i, i] <- -1 / sqrt(i * (i + 1))
    AU[(i + 1), i] <- i / sqrt(i * (i + 1))
  }

  for (j in 1:J) {
    V_y[(n * (j - 1) + 1):(n * j), (J + 2 + (n - 1) * (j - 1)):(J + 1 + (n - 1) * j)] <- AU
  }
  V_x <- matrix(0, nrow = n*sum(J)+sum(J)+1, ncol = n*sum(J)+sum(J)+1)

  # lambda_i = (n + 1)(tau_yx + sigma_yx / n) ((J - 1) pieces)
  for (i in 1:(J - 1)) {
    V_x[1:(i * n), (n * J + 1 + i)] <- -1 / sqrt(i * (i + 1) * (n + 1))
    V_x[(n * J + 1):(n * J + i), (n * J + 1 + i)] <- -1 / sqrt(i * (i + 1) * (n + 1))
    V_x[(i * n + 1):((i + 1) * n), (n * J + 1 + i)] <- i / sqrt(i * (i + 1) * (n + 1))
    V_x[(n * J + i + 1), (n * J + 1 + i)] <- i / sqrt(i * (i + 1) * (n + 1))
  }

  # lambda_nJ + J + 1 = (n * J + J + 1)(tau_yx + sigma_yx / n) / J (1 piece)
  V_x[1:(n * J + J + 1), (n * J + J + 1)] <- 1 / sqrt(n * J + J + 1)
  # check the correctness of spectral decomposition:

  # diff_x = sum(sum((Sigma_x - V_x*D_x*V_x').^2));

  # disp(['difference between covariance matrix for x and its decomposition is ', num2str(diff_x)]);


  # Spectral decomposition of Sigma_yx

  # compute eigenvalues matrix D_yx of Sigma_yx
  D_yx <- matrix(0, n * sum(J) + sum(J) + 1, 1)

  for (i in ((J+2 ): ((n*J)+1))){
    D_yx[i,] <- sigma_yx
  }

  for (i in (((n*J)+2) : ((n*J)+J))){
    D_yx[i,] <- (n+1)*(tau_yx + sigma_yx/n)
  }

  D_yx[(n*J)+J+1,] <- ((n*J)+J+1)*(tau_yx + sigma_yx/n)/J

  D_yx <- diag(D_yx[,1])

  # compute eigenvectors of Sigma_yx

  # create matrix of eigenvectors V_yx for covariance matrix Sigma_yx

  V_yx <- matrix(0, nrow = n*sum(J)+sum(J)+1, ncol = n*sum(J)+sum(J)+1)


  # lambda_i = 0 ((J+1) pieces)

  for (i in 1:(length(seq_len(J)))){
    V_yx[((n*(i-1))+1):(n*i),i] <- -1/sqrt(n*(n+1))
    V_yx[(n*J)+i,i] <- sqrt(n/(n+1))
  }

  V_yx[1:((n*J)+J),J+1] <- -1/sqrt((n*J+J)*(n*J+J+1))
  V_yx[(n*J)+J+1,J+1] <- sqrt((n*J+J)/(n*J+J+1))

  # lambda_i = sigma_yx ((n-1)J pieces)
  # auxiliary matrix AU

  AU <- matrix(0, nrow = n, ncol = n - 1)

  for (i in 1:(n-1)){
    AU[1:i,i] <- -1/sqrt(i*(i+1))
    AU[i+1,i] <- i/sqrt(i*(i+1))
  }

  for (j in 1:length(seq_len(J))){
    V_y[(n*(j-1)+1):(n*j),(J+2+((n-1)*(j-1))):(J+1+((n-1)*j))] <- AU
  }
  # lambda_i = (n+1)(tau_yx + sigma_yx/n) ((J-1) pieces)

  for (i in 1:(J-1)){
    V_y[1:(i*n),(n*J)+1+i] <- -1/sqrt(i*(i+1)*(n+1))
    V_y[((n*J)+1):((n*J)+i),(n*J)+1+i] <- -1/sqrt(i*(i+1)*(n+1))
    V_y[(i*n+1):((i+1)*n),(n*J)+1+i] <- i/sqrt(i*(i+1)*(n+1))
    V_y[(n*J)+i+1,(n*J)+1+i] <- i/sqrt(i*(i+1)*(n+1))
  }

  # lambda_nJ+J+1 = (nJ+J+1)(tau_yx + sigma_yx/n)/J (1 piece)

  V_yx[1:((n*J)+J+1),(n*J)+J+1] <- 1/sqrt(n*J+J+1)

  # check the correctness of spectral decomposition:

  # diff_yx = sum(sum((Sigma_yx - V_yx*D_yx*V_yx').^2));

  # disp(['difference between covariance matrix for x and y and its decomposition is ', num2str(diff_yx)]);


  # Spectral decomposition of Sigma_y

  # compute eigenvalues matrix D_y of Sigma_y

  D_y <- matrix(0, nrow = n*sum(J)+sum(J)+1, ncol = 1)

  for (i in ((J+2) : ((n*J)+1))){
    D_y[i,] <- sigma_y2
  }

  for (i in ((n*J+2) : ((n*J)+J))){
    D_y[i,] <- (n+1)*(tau_y2 + sigma_y2/n)
  }

  D_y[(n*J)+J+1,] <- ((n*J)+J+1)*(tau_y2 + sigma_y2/n)/J

  D_y <- diag(D_y[,1])

  # compute eigenvectors of Sigma_y

  # create matrix of eigenvectors V_y for covariance matrix Sigma_y

  V_y <- matrix(0, nrow = n*sum(J)+sum(J)+1, ncol = n*sum(J)+sum(J)+1)

  # lambda_i = 0 ((J+1) pieces)

  for (i in 1:(length(seq_len(J)))){
    V_y[((n*(i-1))+1):(n*i),i] <- -1/sqrt(n*(n+1))
    V_y[(n*J)+i,i] <- sqrt(n/(n+1))
  }

  V_y[1:((n*J)+J),J+1] <- -1/sqrt((n*J+J)*(n*J+J+1))
  V_y[(n*J)+J+1,J+1] <- sqrt((n*J+J)/(n*J+J+1))

  # lambda_i = sigma_yx ((n-1)J pieces)

  # auxiliary matrix AU

  AU <- matrix(0, nrow = n, ncol = n - 1)
  for (i in 1:(n-1)){
    AU[1:i,i] <- -1/sqrt(i*(i+1))
    AU[i+1,i] <- i/sqrt(i*(i+1))
  }

  for (j in 1:J){
    V_y[(n*(j-1)+1):(n*j),(J+2+((n-1)*(j-1))):(J+1+((n-1)*j))] <- AU
  }

  # lambda_i = (n+1)(tau_yx + sigma_yx/n) ((J-1) pieces)

  for (i in 1:(J-1)){
    V_y[1:(i*n),(n*J)+1+i] <- -1/sqrt(i*(i+1)*(n+1))
    V_y[((n*J)+1):((n*J)+i),(n*J)+1+i] <- -1/sqrt(i*(i+1)*(n+1))
    V_y[(i*n+1):((i+1)*n),(n*J)+1+i] <- i/sqrt(i*(i+1)*(n+1))
    V_y[(n*J)+i+1,(n*J)+1+i] <- i/sqrt(i*(i+1)*(n+1))
  }

  # lambda_nJ+J+1 = (nJ+J+1)(tau_yx + sigma_yx/n)/J (1 piece)

  V_y[1:((n*J)+J+1),(n*J)+J+1] <- 1/sqrt(n*J+J+1)

  # check the correctness of spectral decomposition:

  # diff_y = sum(sum((Sigma_y - V_y*D_y*V_y').^2));

  # disp(['difference between covariance matrix for y and its decomposition is ', num2str(diff_y)]);


  # calculate the eigenvalues of transformed sum Z_x'*A*Z_x =
  # W_x'*(Sigma_x^0.5)*A*(Sigma_x^0.5)*W_x = {use spectral decomposition and
  # substitute sqrt(D) with S=sqrt(D), because D is diagonal} =
  # = W_x*V_x*S_x*V_x'*A*V_x*S_x*V_x'*W_x = W_x*C_x*W_x
  # C_x = V_x*S_x*V_x'*A*V_x*S_x*V_x' with internal part
  # L1 = S_x*V_x'*A*V_x*S_x  - diagonal -> gives eigenvalues of C_x, which we use
  # to calculate MSE:

  S_x <- sqrt(abs(D_x))
  t_v_x <-t(V_x)
  L1 <- S_x %*% t_v_x %*% A %*% V_x %*% S_x

  L1 <- diag(L1)

  L1 <- Re(L1)# could be complex insignificant tails, because of computer precision
  
  for (i in 1:(n*J+J+1)){ # delete computer precision tails
    if (abs(L1[i])<0.0000001){
      L1[i] <- 0
    }
  }

  # calculate the eigenvalues of transformed sum Z_yx'*A*Z_yx =
  # W_x'*(Sigma_x^0.5)*A*(Sigma_x^0.5)*W_y = {use spectral decomposition and
  # substitute sqrt(D) with S=sqrt(D), because D is diagonal} =
  # {W_x, W_y ~N(0,I)} = W_x'*V_x*S_x*V_x'*A*V_x*S_y*V_y'*W_y =
  # {H_x = V_x'*W_x ~ N(0,I), H_y = V_y'*W_y ~ N(0,I)} =
  # = H_x'*S_x*V_x'*A*V_x*S_y*H_y = H_x'*Q*H_y = {H = [H_x, H_y]} = H'*Q1*H1
  # = H' * EH^(-0.5) * EH^(0.5) * Q1 * EH^(0.5) * EH^(-0.5) *H =
  # = {EH^(-0.5) * H = H1 ~ N(0,I)} = H1' *EH^(0.5) * Q1 * EH^(0.5) * H1 =
  # {use spectral decomposition and substitute: EH = V_H*D_H*V_H' ->
  # EH^(0.5) = V_H*S_H*V_H' with with S_H=sqrt(D_H), because D_H is diagonal}
  # H1' * V_H * S_H * V_H' * Q1 * V_H * S_H * V_H' * H1 =
  # {H1_tilde = V_H' * H1 ~ N(0,I)} =
  # = H1_tilde' * S_H * V_H' * Q1 * V_H * S_H * H1_tilde, with internal part
  # L2 = S_H * V_H' * Q1 * V_H * S_H  - diagonal -> gives eigenvalues of L2,
  # which we use to calculate MSE:

  S_y <- sqrt(D_y)

  Q <- S_x %*% t_v_x %*% A %*% V_y %*% S_y ## wrong

  #
  upper_left <- matrix(0, nrow(Q), ncol(Q))
  lower_right <- matrix(0, nrow(Q), ncol(Q))

  upper_right <- Q / 2
  lower_left <- Q / 2

  Q1 <- rbind(cbind(upper_left, upper_right), cbind(lower_left, lower_right))

  #


  #
  upper_left <- diag(nrow(D_yx))
  lower_right <- diag(nrow(D_yx))

  pinv_S_x_S_y_D_yx <- pracma::pinv(S_x) %*% pracma::pinv(S_y) %*% D_yx


  upper_right <- pinv_S_x_S_y_D_yx
  lower_left <- pinv_S_x_S_y_D_yx

  EH <- rbind(cbind(upper_left, upper_right), cbind(lower_left, lower_right))
  #


  D1 <- EH[(1:(n*J+J+1)), ((n*J+J+2):(2*(n*J+J+1)))]
  D_H <- matrix(0, nrow = nrow(EH), ncol = ncol(EH))
  for (i in 1:(n*J+J+1)){
    D_H[2*i-1,2*i-1] <- 1 + D1[i,i]
    D_H[2*i,2*i] <- 1 - D1[i,i]
  }
  S_H <- ifelse(D_H>0,sqrt(D_H),0)
  rm(D1)

  ##--


  V_H <- matrix(0, nrow = nrow(D_H), ncol = ncol(D_H))

  for (i in 1:(n*J+J+1)){
    V_H[i,(2*i)-1] <- 1/sqrt(2)
    V_H[i, 2*i] <- 1/sqrt(2)
    V_H[(n*J)+J+1 + i, (2*i)-1] <- 1/sqrt(2)
    V_H[(n*J)+J+1 + i, 2*i] <- -1/sqrt(2)
  }

  L2 <- S_H %*% t(V_H) %*% Q1 %*% V_H %*% S_H

  L2 <- diag(L2)
  L2 <- Re(L2)# could be complex insignificant tails, because of computer precision

  for (i in 1:(2*(n*J+J+1))){# delete computer precision tails
    if (abs(L2[i])<0.0000001){
      L2[i] <- 0
    }
  }


  # Use L1 and L2 and ML estimator of beta_b as coefficients for computing optimal MSE

  beta_b_ML <- tau_yx / tau_x2

  K_sum1 <- sum(L1)^2/(2*(t(L1)%*%L1))# unused

  T_sum1 <- (t(L1)%*%L1)/sum(L1)# unused

  K_sum2 <- sum(L2)^2/(2*(t(L2)%*%L2))

  T_sum2 <- (t(L2)%*%L2)/sum(L2)

  beta_b <- data$b2# real value of between parameter ##check if b_b is same as b2

  # this is the case when between-group effect is not zero
  if (is.finite(K_sum2) && is.finite(T_sum2)){
    
  # MSE of ML estimated parameter beta_b_ML

  MSE_add <- (((K_sum2 * T_sum2)/ (T_sum1*(K_sum1-1))) - beta_b)^2 +
    ((K_sum2 * T_sum2^2 * (K_sum1 + K_sum2 - 1)) / (T_sum1^2 * (K_sum1-1)^2 * (K_sum1-2)))

  # MSE_add = (K_sum1^2 * T_sum1/ (K_sum2 * T_sum2*(K_sum1-1)) - beta_b)^2 +...
  #    K_sum1^4 * T_sum1^2 * (K_sum1+K_sum2-1)/ (K_sum2^3 * T_sum2^2 * (K_sum1-1)^2 * (K_sum1-2)); # original


  # optimize MSE with beta_b_ML(old_version )


  #Tau02 <- kronecker(    matrix(1,  nrow = length(seq(0, 1, by = 0.01)), ncol = 1)   , t( matrix(seq(0.05, 10, by = 0.05), nrow = 1) )  )# restrict search interval of tau02 to 10

  #W <-kronecker(t(matrix(seq(0, 1, by = 0.01), nrow=1)), matrix(1,nrow = length(seq(0.5, by = 0.05)), ncol=1))
  ## recheck???
  
  # optimize MSE with beta_b_ML within 5*sigma search region from the distribution of tau_x2
  
  d_search <- 5
  radius   <- sqrt(K_sum1 * T_sum1^2)
  
  # Attempt the usual endpoints
  Tau02_min <- max(0.01, tau_x2 - d_search * radius)
  Tau02_max <- max(0.01, tau_x2 + d_search * radius)
  
  # Fallback to a fixed window if anything’s off
  if (!all(is.finite(c(Tau02_min, Tau02_max))) || Tau02_min >= Tau02_max) {
    Tau02_min <- 0
    Tau02_max <- 10
  }
  
  # Build the sequence
  Tau02_seq <- seq(
    from       = Tau02_min,
    to         = Tau02_max,
    length.out = 101
  )
  
  W_seq <- seq(0, 1, by = 0.01)
  
  grid <- expand.grid(Tau02 = Tau02_seq, W = W_seq)
  
  Tau02 = grid$Tau02
  
  W = grid$W

  MSE_ML <- matrix(0, nrow = length(W), ncol = 1)



  # use grid search to find w and tau02 that give smallest MSE
  for (i in 1:length(W)){
    MSE_ML[i] <- ((K_sum2 * T_sum2^2 * (K_sum2+1) * ((1-W[i])*Tau02[i] + W[i]*sum(L1)))/
                    ((((1-W[i])*Tau02[i] + W[i]*sum(L1))^2 - 2*W[i]^2*(t(L1)%*%L1)) *
                       (((1-W[i])*Tau02[i] + W[i]*sum(L1))^2 - 4*W[i]^2*(t(L1)%*%L1))) -
                    (2*beta_b_ML *K_sum2 * T_sum2 * ((1-W[i])*Tau02[i] + W[i]*sum(L1)))/
                    ((((1-W[i])*Tau02[i] + W[i]*sum(L1))^2) - (2*W[i]^2*(t(L1)%*%L1)))+ beta_b_ML)

  }

  # optimize MSE with beta_b

  MSE <- matrix(0, nrow = length(W), ncol = 1)


  # use grid search to find w and tau02 that give smallest MSE
  for (i in 1:length(W)){
    MSE[i] <- ((K_sum2 * T_sum2^2 * (K_sum2+1) * ((1-W[i])*Tau02[i] + W[i]*sum(L1)))/
                 ((((1-W[i])*Tau02[i] + W[i]*sum(L1))^2 - 2*W[i]^2*(t(L1)%*%L1)) *
                    (((1-W[i])*Tau02[i] + W[i]*sum(L1))^2 - 4*W[i]^2*(t(L1)%*%L1))) -
                 (2*data$b_b *K_sum2 * T_sum2 * ((1-W[i])*Tau02[i] + W[i]*sum(L1)))/
                 ((((1-W[i])*Tau02[i] + W[i]*sum(L1))^2) - (2*W[i]^2*(t(L1)%*%L1)))+data$b_b)



  }


  # use MSE to find beta_b_Bay with smallest MSE:

  #
  ind_MSEbetabBay <- which.min(MSE)
  MSE_beta_b_Bay <- MSE[ind_MSEbetabBay]
  #


  w_opt <- W[ind_MSEbetabBay]
  tau02_opt <- Tau02[ind_MSEbetabBay]

  beta_b_Bay <- tau_yx/((1-w_opt)*tau02_opt + w_opt*tau_x2)


  # use MSE_ML to find beta_b_Bay_ML with smallest MSE:

  #
  ind_MSEbetabBayML <- which.min(MSE_ML)
  MSE_beta_b_Bay_ML <- MSE_ML[ind_MSEbetabBayML]
  #

  w_opt_ML <- W[ind_MSEbetabBayML]
  tau02_opt_ML <- Tau02[ind_MSEbetabBayML]

  beta_b_Bay_ML <- tau_yx/((1-w_opt_ML)*tau02_opt_ML + w_opt_ML*tau_x2)



  # Compute standard errors of beta_b_Bay, beta_b_Bay_ML, beta_b_ML using their distributions

  if (w_opt != 0) {
    K_B_Bay <- (w_opt * T_sum1 * K_sum1 + (1 - w_opt) * tau02_opt)^2 / (w_opt^2 * T_sum1^2 * K_sum1)

    T_B_Bay <- (w_opt^2 * T_sum1^2 * K_sum1) / (w_opt * T_sum1 * K_sum1 + (1 - w_opt) * tau02_opt)
  } else {
    K_B_Bay <- 4  # to avoid negative values in SE_beta_Bay

    T_B_Bay <- 0.25
  }

  if (w_opt_ML != 0) {
    K_B_Bay_ML <- (w_opt_ML * T_sum1 * K_sum1 + (1 - w_opt_ML) * tau02_opt_ML)^2 / (w_opt_ML^2 * T_sum1^2 * K_sum1)

    T_B_Bay_ML <- (w_opt_ML^2 * T_sum1^2 * K_sum1) / (w_opt_ML * T_sum1 * K_sum1 + (1 - w_opt_ML) * tau02_opt_ML)
  } else {
    K_B_Bay_ML <- 4  # to avoid negative values in SE_beta_Bay_ML

    T_B_Bay_ML <- 0.25
  }



  SE_beta_ML <- abs((T_sum2 / (T_sum1 * (K_sum1 - 1))) * sqrt(abs((K_sum2 * (K_sum1 + K_sum2 - 1)) / (K_sum1 - 2))))

  SE_beta_Bay <- abs((T_sum2 / (T_B_Bay * (K_B_Bay - 1))) * sqrt(abs((K_sum2 * (K_B_Bay + K_sum2 - 1)) / (K_B_Bay - 2))))

  SE_beta_Bay_ML <- abs((T_sum2 / (T_B_Bay_ML * (K_B_Bay_ML - 1))) * sqrt(abs((K_sum2 * (K_B_Bay_ML + K_sum2 - 1)) / (K_B_Bay_ML - 2))))

  }
  
  # this is the case when L2 is zero, therefore tau_yx is also defined as zero ->
  #-> no between-group effect, and K_sum2 and T_sum2 are NaN
  else
  {
    MSE <- 1
    MSE_ML <- 1
    W <- 0
    Tau02 <- 1
    beta_b_Bay <- 0
    beta_b_Bay_ML <- 0
    MSE_beta_b_Bay <- 1
    MSE_beta_b_Bay_ML <- 1
    MSE_add <- 1
    SE_beta_ML <- sqrt(abs((tau_y2 + sigma_y2 / n) / (J * tau_x2)))
    SE_beta_Bay <- 1
    SE_beta_Bay_ML <- 1
    warning('Singularities in the between-group estimation')
  }  
  
  
  #Bay <- 0.1
  Bay <- list()
  Bay$A <- A
  Bay$Sigma_x <- Sigma_x
  Bay$Sigma_yx <- Sigma_yx
  Bay$Sigma_y <- Sigma_y
  Bay$tau_x2 <- tau_x2
  Bay$sigma_x2 <- sigma_x2
  Bay$tau_yx <- tau_yx
  Bay$sigma_yx <- sigma_yx
  Bay$tau_y2 <- tau_y2
  Bay$sigma_y2 <- sigma_y2
  Bay$D_x <- D_x
  Bay$V_x <- V_x
  Bay$D_yx <- D_yx
  Bay$V_yx <- V_yx
  Bay$D_y <- D_y
  Bay$V_y <- V_y
  Bay$L1 <- L1
  Bay$L2 <- L2
  Bay$MSE <- MSE
  Bay$MSE_ML <- MSE_ML
  Bay$W <- W
  Bay$Tau02 <- Tau02
  Bay$beta_b_ML <- beta_b_ML
  Bay$beta_b_Bay <- beta_b_Bay
  Bay$beta_b_Bay_ML <- beta_b_Bay_ML
  
  if (exists("gamma") && !is.function(gamma)){
    Bay$gamma <- gamma
  }
  
  Bay$MSE_beta_b_Bay <- MSE_beta_b_Bay
  Bay$MSE_beta_b_Bay_ML <- MSE_beta_b_Bay_ML
  Bay$MSE_add <- MSE_add
  Bay$SE_beta_ML = SE_beta_ML
  Bay$SE_beta_Bay = SE_beta_Bay
  Bay$SE_beta_Bay_ML = SE_beta_Bay_ML
  
  if (exists("SE_gamma")){
    Bay$SE_gamma = SE_gamma
  }

  #Bay_list <- list( 'A'=A,  'Sigma_x'=Sigma_x, 'Sigma_y'=Sigma_y, 'Sigma_yx'=Sigma_yx, 'tau_x2'=tau_x2,
  #                  'sigma_x2'=sigma_x2, 'tau_y2'=tau_y2, 'sigma_y2'=sigma_y2,  'tau_yx'=tau_yx,
  #                  'sigma_yx'=sigma_yx,'D_x'=D_x,  'D_y'=D_y, 'D_yx'=D_yx,'V_x'=V_x, 'V_y'=V_y,
  #                  'V_yx'=V_yx,'L1'=L1, 'L2'=L2,  'MSE'=MSE,  'MSE_ML'=MSE_ML,'W'=W,'Tau02'=Tau02,
  #                  'beta_b_ML'=beta_b_ML, 'beta_b_Bay'=beta_b_Bay,'beta_b_Bay_ML'=beta_b_Bay_ML,
  #                  'gamma' = gamma, 'MSE_beta_b_Bay'=MSE_beta_b_Bay,'MSE_beta_b_Bay_ML'=MSE_beta_b_Bay_ML,
  #                  'MSE_add'= MSE_add, 'SE_beta_ML' = SE_beta_ML, 'SE_beta_Bay' = SE_beta_Bay,
  #                  'SE_beta_Bay_ML' = SE_beta_Bay_ML, 'SE_gamma' = SE_gamma)

  #print(Bay_list)
  
  Bay_list = Bay
}


estimate_Bay_CV_SE_jackknife_individual <- function(data) {
 

  # Original estimator calculation (to get initial results)
  original_result <- estimate_Bay_CV(data)

  # Initialize parameters
  n <- data$n
  J <- data$k
  r <- min(n,50) # set up the number of replications for jackknife, r from 1 to n^J

  # Initialize vectors to store jackknife estimates for individual deletion
  jackknife_beta_b_Bay_individual <- numeric(r)
  jackknife_beta_b_Bay_ML_individual <- numeric(r)
  jackknife_beta_b_ML_individual <- numeric(r)
  
  # only if gamma‐covariates present
  has_gamma <- ("kc" %in% names(data) && data$kc > 0)
  if (has_gamma) {
    p_gamma <- data$kc
    jackknife_gamma_individual <- matrix(NA, nrow = r, ncol = p_gamma)
  }

  # Jackknife by deleting one individual from each group simultaneously
  for (i in seq_len(r)) {
    # Generate J random indices (one for each group) from 1 to n
    random_indices <- sample(1:n, J, replace = TRUE)

    # Initialize a list to store indices to keep for each group
    indices_to_keep <- vector("list", J)

    for (j in 1:J) {
      indices_to_keep[[j]] <- (n * (j - 1) + 1):(n * j)
      indices_to_keep[[j]] <- indices_to_keep[[j]][-random_indices[j]]
    }

    # Create new indices after deleting one individual from each group
    indices <- unlist(indices_to_keep)

    # Create a new dataset excluding the selected individuals
    data_jackknife <- data
    data_jackknife$x <- data$x[indices]
    data_jackknife$y <- data$y[indices]
    
    if ("kc" %in% names(data)){
      if (data$kc>1){
        data_jackknife$C <- data$C[indices, , drop = FALSE]
      } else {
        data_jackknife$C <- data$C[indices, drop = FALSE]
      }
    }
    data_jackknife$n <- n - 1
    data_jackknife$kn <- (n-1)*J

    # Recalculate estimators without the selected individuals
    result_jackknife <- estimate_Bay_CV(data_jackknife)

    # Store jackknife estimates
    jackknife_beta_b_Bay_individual[i] <- result_jackknife$beta_b_Bay
    jackknife_beta_b_Bay_ML_individual[i] <- result_jackknife$beta_b_Bay_ML
    jackknife_beta_b_ML_individual[i] <- result_jackknife$beta_b_ML
    
    if (has_gamma) {
      # assume res_j$gamma is a length‐kc numeric vector
      jackknife_gamma_individual[i, ] <- result_jackknife$gamma
    }
    
  }

  # Calculate jackknife means for individual deletion
  jackknife_mean_beta_b_Bay_individual <- mean(jackknife_beta_b_Bay_individual)
  jackknife_mean_beta_b_Bay_ML_individual <- mean(jackknife_beta_b_Bay_ML_individual)
  jackknife_mean_beta_b_ML_individual <- mean(jackknife_beta_b_ML_individual)
  if (has_gamma) {
    jackknife_mean_gamma_individual <- colMeans(jackknife_gamma_individual)
  }

  # Calculate jackknife standard errors for individual deletion
  
  SE <- function(vals, m) {
    sqrt((n - 1) / n * sum((vals - m)^2))
  }
  
  SE_beta_Bay_jackknife_individual       <- SE(jackknife_beta_b_Bay_individual, jackknife_mean_beta_b_Bay_individual)
  SE_beta_Bay_ML_jackknife_individual    <- SE(jackknife_beta_b_Bay_ML_individual, jackknife_mean_beta_b_Bay_ML_individual)
  SE_beta_ML_jackknife_individual        <- SE(jackknife_beta_b_ML_individual, jackknife_mean_beta_b_ML_individual)
  
  if (has_gamma) {
    
    SE_gamma_jackknife_individual <- sapply(
      seq_len(p_gamma),
      function(j) {
        SE(jackknife_gamma_individual[, j], m = jackknife_mean_gamma_individual[j])
      }
    )
    
    names(SE_gamma_jackknife_individual) <- paste0("SE_gamma_", seq_len(p_gamma))
  }
  
  # Output updated results with jackknife standard errors
  
  # 7. Return everything
  out <- list(
    beta_b_Bay                          = original_result$beta_b_Bay,
    beta_b_Bay_ML                       = original_result$beta_b_Bay_ML,
    beta_b_ML                           = original_result$beta_b_ML,
    SE_beta_Bay_jackknife_individual    = SE_beta_Bay_jackknife_individual,
    SE_beta_Bay_ML_jackknife_individual = SE_beta_Bay_ML_jackknife_individual,
    SE_beta_ML_jackknife_individual     = SE_beta_ML_jackknife_individual
  )
  
  if (has_gamma) {
    out$gamma                           <- original_result$gamma
    out$SE_gamma_jackknife_individual   <- SE_gamma_jackknife_individual
  }
  
  out
  
}

# Example usage:
# You would use your actual `estimate_Bay_CV` function here, which returns a complete `Bay` list.
# result <- mlob(y ~ X + C1 + C2, data = data, group = "Group", conf.level = 0.01, jackknife = FALSE)
# print(result)
# summary(result)

