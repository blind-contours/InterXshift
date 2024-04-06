#' @title Calculates the Inverse Variance Pooled Estimate Including Null Folds
#' @description Does a weighted combination estimate for folds with estimates
#' and null folds finding 0 using inverse variance
#' @param results_df Table of results
#' @param n_folds Total number of folds specified
#' @export

calculatePooledEstimate <- function(results_df, n_folds, delta = NULL) {
  # Identifying the row with the current pooled TMLE
  pooledRow <- results_df[results_df$Fold == "Pooled TMLE", ]
  calculatePooledEstimate <- function(results_df, n_folds, delta = NULL) {
    # Identifying the row with the current pooled TMLE
    pooledRow <- results_df[results_df$Fold == "Pooled TMLE", ]

    # Calculate n_0 and n_1
    n_1 <- nrow(results_df[results_df$Fold != "Pooled TMLE", ])
    n_0 <- n_folds - n_1

    # Variance of the current pooled estimates
    var_pooled <- pooledRow$Variance

    # Estimating the variance of the "null" estimate
    var_null <- var_pooled * n_1 / n_0

    # Standard deviation for the null distribution
    sd_null <- sqrt(var_null)

    # Sample from the null distribution for each null fold
    null_samples <- rnorm(n_0, mean = 0, sd = sd_null)

    # Weights
    w_1 <- 1 / var_pooled
    w_0 <- 1 / var_null

    # Calculate new pooled estimate and its variance
    pooled_estimate <- pooledRow$Psi

    # Calculate the weighted average using the sampled null values
    new_pooled_psi <- (w_1 * pooled_estimate + sum(w_0 * null_samples)) / (w_1 + sum(w_0))

    # Recalculate the variance
    var_new_pooled <- 1 / (w_1 + sum(w_0))

    # Calculate standard error and confidence intervals
    se_new_pooled <- sqrt(var_new_pooled)
    lower_ci <- new_pooled_psi - 1.96 * se_new_pooled
    upper_ci <- new_pooled_psi + 1.96 * se_new_pooled

    p_value_pooled <- 2 * stats::pnorm(abs(new_pooled_psi / se_new_pooled), lower.tail = F)

    # Add new row to the table
    new_row <- data.table(
      Condition = pooledRow$Condition,
      Psi = new_pooled_psi,
      Variance = var_new_pooled,
      SE = se_new_pooled,
      `Lower CI` = lower_ci,
      `Upper CI` = upper_ci,
      `P-value` = p_value_pooled,
      Fold = "Inverse Variance Pooled",
      Type = pooledRow$Type,
      Variables = pooledRow$Variables,
      N = pooledRow$N,
      Delta = pooledRow$Delta
    )

    bind_rows(results_df, new_row)
  }
}




#' @title Calculates the Individual Shift Parameter
#' @description Simply subtract the psi estimates of a shift compared to no
#' shift. Use the delta method to do the same thing for the eif to derive
#' variance estimates for this new parameter.
#' @param tmle_fit TMLE results for the individual shift
#' @param exposure The exposure identified
#' @param fold_k Fold exposure is being shifted in
#' @export
calc_final_ind_shift_param <- function(tmle_fit, exposure, fold_k) {
  condition <- exposure
  psi_param <-  tmle_fit$psi - tmle_fit$noshift_psi
  variance_est <- var(tmle_fit$eif - tmle_fit$noshift_eif) /
    length(tmle_fit$eif)
  se_est <- sqrt(variance_est)
  CI <- calc_CIs(psi_param, se_est)

  Lower_CI <- CI[1]
  Upper_CI <- CI[2]
  p.value <- 2 * stats::pnorm(abs(psi_param / se_est), lower.tail = F)

  n <- length(tmle_fit$eif)

  results <- data.table::data.table(
    condition, psi_param, variance_est, se_est,
    Lower_CI, Upper_CI, p.value, fold_k,
    "Indiv Shift", exposure, n
  )

  names(results) <- c(
    "Condition", "Psi", "Variance", "SE", "Lower CI",
    "Upper CI", "P-value", "Fold", "Type", "Variables", "N"
  )

  return(results)
}


#' @title Calculates the Joint Shift Parameter
#' @description Estimates the shift parameter for a joint shift
#' @param joint_shift_fold_results Results of the joint shift
#' @param rank ranking of the interaction found
#' @param exposures Exposures shifted
#' @param fold_k Fold the joint shift is identified
#' @param deltas_updated The new delta, could be updated if Hn has positivity
#' violations
#' @export

calc_final_joint_shift_param <- function(joint_shift_fold_results,
                                         rank,
                                         fold_k,
                                         deltas_updated,
                                         exposures) {
  results <- lapply(joint_shift_fold_results, calc_joint_results)
  results_table <- do.call(rbind, results)

  intxn_results <- calc_intxn_results(results_table, joint_shift_fold_results)

  joint_intxn_results <- rbind(
    as.data.frame(results_table),
    t(as.data.frame(unlist(intxn_results)))
  )

  joint_intxn_results <- as.data.frame(cbind(
    rep(rank, 4),
    joint_intxn_results,
    rep(fold_k, 4),
    length(joint_shift_fold_results[[3]]$eif),
    deltas_updated[[1]],
    deltas_updated[[2]],
    c(exposures[[1]], exposures[[2]], paste(exposures[[3]], collapse = "-"), "Interaction")
  ))

  names(joint_intxn_results) <- c(
    "Rank", "Psi", "Variance",
    "SE", "Lower CI", "Upper CI",
    "P-value", "Fold", "N",
    "Delta Exposure 1",
    "Delta Exposure 2",
    "Type"
  )
  rownames(joint_intxn_results) <- NULL
  return(joint_intxn_results)
}

#' @title Calculates the Mediation Shift Parameters
#' @description Estimates the shift parameter for a natural direct and indirect
#' effect
#' @param tmle_fit_a_shift TMLE results for a shift in A alone
#' @param tmle_fit_a_z_shift TMLE results for a shift in A and Z
#' @param exposure The exposure found
#' @param mediator The mediator found
#' @param fold_k Fold the exposure and mediator were found
#' @param y Outcome
#' @param delta The shift amount


#' @export

calc_mediation_param <- function(tmle_fit_a_shift,
                                 tmle_fit_a_z_shift,
                                 exposure,
                                 mediator,
                                 fold_k,
                                 y,
                                 delta) {
  condition <- paste(exposure, mediator, sep = "&")

  nde_est <- tmle_fit_a_shift$psi - mean(y)
  nie_est <- tmle_fit_a_z_shift$psi - tmle_fit_a_shift$psi

  no_shift_eif <- tmle_fit_a_z_shift$noshift_eif

  nde_eif <- tmle_fit_a_shift$eif - no_shift_eif
  nie_eif <- tmle_fit_a_z_shift$eif - tmle_fit_a_shift$eif

  var_nde <- var(nde_eif) /
    length(nde_eif)
  se_nde <- sqrt(var_nde)
  ci_nde <- calc_CIs(nde_est, se_nde)
  lower_ci_nde <- ci_nde[1]
  upper_ci_nde <- ci_nde[2]
  p_value_nde <- 2 * stats::pnorm(abs(nde_est / se_nde), lower.tail = F)

  var_nie <- var(nie_eif) /
    length(nie_eif)
  se_nie <- sqrt(var_nie)
  ci_nie <- calc_CIs(nie_est, se_nie)
  lower_ci_nie <- ci_nie[1]
  upper_ci_nie <- ci_nie[2]
  p_value_nie <- 2 * stats::pnorm(abs(nie_est / se_nie), lower.tail = F)

  n <- length(nie_eif)

  nde_results <- data.table::data.table(
    "NDE", nde_est, var_nde, se_nde,
    lower_ci_nde, upper_ci_nde, p_value_nde, fold_k,
    "Mediation", n
  )

  nie_results <- data.table::data.table(
    "NIE", nie_est, var_nie, se_nie,
    lower_ci_nie, upper_ci_nie, p_value_nie, fold_k,
    "Mediation", n
  )

  names(nde_results) <- c(
    "Parameter", "Estimate", "Variance", "SE", "Lower CI",
    "Upper CI", "P-value", "Fold", "Type", "N"
  )


  names(nie_results) <- c(
    "Parameter", "Estimate", "Variance", "SE", "Lower CI",
    "Upper CI", "P-value", "Fold", "Type", "N"
  )

  nde_results$Delta <- delta
  nie_results$Delta <- delta

  mediation_results_table <- rbind(nde_results, nie_results)

  return(mediation_results_table)
}
