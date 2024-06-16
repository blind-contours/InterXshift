#' @title InterXshift
#'
#' @description Under a fixed shift to exposures identify using g-computation the joint shift of
#' pairwise exposures in a mixed exposure compared to the additive individual shifts. Positive values indicate
#' synergy and negative antagonism, get the top synergy and antagonism results and use CV-TMLE to efficiently
#' estimate the interaction target parameter.
#'
#' @param w A \code{matrix}, \code{data.frame}, or similar containing a set of
#' baseline covariates. These variables are measured before exposures.
#' @param a \code{matrix}, \code{data.frame}, or similar containing individual or
#' multiple exposures.
#' @param z \code{matrix}, \code{data.frame}, or similar containing individual or
#' multiple mediators (optional).
#' @param y \code{numeric} vector of observed outcomes.
#' @param deltas A \code{numeric} value indicating the shift in exposures to
#' define the target parameter, with respect to the scale of the exposures (A). If adaptive_delta
#' is true, these values will be reduced.
#' @param var_sets A list specifying variable sets for deterministic InterXshift usage.
#' Example: var_sets <- c("A_1", "A_1-Z_2") where the analyst provides variable sets
#' for exposures, exposure-mediator, or exposure-covariate relationships.
#' @param estimator The type of estimator to fit: \code{"tmle"} for targeted
#' maximum likelihood estimation, or \code{"onestep"} for a one-step estimator.
#' @param fluctuation Method used in the targeting step for TML estimation: "standard" or "weighted".
#' This determines where to place the auxiliary covariate in the logistic tilting regression.
#' @param pi_learner Learners for fitting Super Learner ensembles to densities via \pkg{sl3}.
#' @param mu_learner Learners for fitting Super Learner ensembles to the outcome model via \pkg{sl3}.
#' @param g_learner Learners for fitting Super Learner ensembles to the g-mechanism
#' g(A|W) (a probability estimator, not a density estimator) for mediation via \pkg{sl3}.
#' @param e_learner Learners for fitting Super Learner ensembles to the e-mechanism
#' g(A|Z,W) (a probability estimator, not a density estimator) for mediation via \pkg{sl3}.
#' @param zeta_learner Learners for fitting Super Learner ensembles to the outcome model via \pkg{sl3}..
#' @param n_folds Number of folds to use in cross-validation, default is 2.
#' @param outcome_type Data type of the outcome, default is "continuous".
#' @param parallel Whether to parallelize across cores (default: TRUE).
#' @param parallel_type Type of parallelization to use if parallel is TRUE:
#' "multi_session" (default), "multicore", or "sequential".
#' @param num_cores Number of CPU cores to use in parallelization (default: 2).
#' @param seed \code{numeric} seed value to be passed to all functions.
#' @param hn_trunc_thresh Truncation level for the clever covariate (default: 10).
#' @param adaptive_delta If TRUE, reduces the user-specified delta until
#' the Hn calculated for a shift does not have any observation greater
#' than hn_trunc_thresh (default: FALSE).
#'
#' @return An S3 object of class \code{InterXshift} containing the results of the
#' procedure to compute a TML or one-step estimate of the counterfactual mean
#' under a modified treatment policy that shifts a continuous-valued exposure
#' by a scalar amount \code{delta}. These exposures are data-adaptively
#' identified using the CV-TMLE procedure.
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%
#' @importFrom stats as.formula glm p.adjust plogis predict qlogis qnorm qunif rnorm runif
#' @importFrom rlang :=
#' @importFrom stringr str_count
#' @import furrr
#' @importFrom purrr map
#' @importFrom data.table rbindlist

InterXshift <- function(w,
                      a,
                      y,
                      deltas,
                      estimator = "tmle",
                      fluctuation = "standard",
                      var_sets = NULL,
                      pi_learner = NULL,
                      mu_learner = NULL,
                      g_learner = NULL,
                      n_folds = 2,
                      outcome_type = "continuous",
                      parallel = TRUE,
                      parallel_type = "multi_session",
                      num_cores = 2,
                      seed = seed,
                      hn_trunc_thresh = 10,
                      adaptive_delta = FALSE,
                      discover_only = FALSE,
                      top_n = 2,
                      density_classification = TRUE) {
  # check arguments and set up some objects for programmatic convenience
  call <- match.call(expand.dots = TRUE)
  estimator <- match.arg(estimator)
  fluctuation <- match.arg(fluctuation)

  # coerce W to matrix and, if no names in W, assign them generically
  if (!is.data.frame(w)) w <- as.data.frame(w)
  w_names <- colnames(w)
  if (is.null(w_names)) {
    w_names <- paste0("w", seq_len(ncol(w)))
    colnames(w) <- w_names
  }


  # coerce W to matrix and, if no names in W, assign them generically
  a <- data.frame(a)
  a_names <- colnames(a)

  if (is.null(a_names)) {
    a_names <- paste0("a", seq_len(ncol(a)))
    colnames(a) <- a_names
  }

  if (is.null(pi_learner)) {
    sls <- create_sls()
    pi_learner <- sls$pi_learner
  }

  if (is.null(mu_learner)) {
    sls <- create_sls()
    mu_learner <- sls$mu_learner
  }


  if (is.null(g_learner)) {
    sls <- create_sls()
    g_learner <- sls$g_learner
  }

  if (parallel == TRUE) {
    if (parallel_type == "multi_session") {
      future::plan(future::multisession,
                   workers = num_cores,
                   gc = TRUE
      )
    } else {
      future::plan(future::multicore,
                   workers = num_cores,
                   gc = TRUE
      )
    }
  } else {
    future::plan(future::sequential,
                 gc = TRUE
    )
  }

  data_internal <- data.table::data.table(w, a, y)
  `%notin%` <- Negate(`%in%`)

  if (outcome_type == "binary") {
    ## create the CV folds
    data_internal$folds <- create_cv_folds(n_folds, data_internal$y)
  } else {
    data_internal$folds <- create_cv_folds(n_folds, data_internal$y)
  }

  fold_basis_results <- furrr::future_map(unique(data_internal$folds),
                                          function(fold_k) {
                                            at <- data_internal[data_internal$folds != fold_k, ]
                                            av <- data_internal[data_internal$folds == fold_k, ]

                                            intxn_results <- find_synergy_antagonism(
                                              data = at,
                                              deltas = deltas,
                                              a_names = a_names,
                                              w_names = w_names,
                                              outcome = "y",
                                              outcome_type = outcome_type,
                                              mu_learner = mu_learner,
                                              seed = seed,
                                              top_n = top_n
                                            )

                                          },
                                          .options = furrr::furrr_options(seed = seed, packages = "InterXshift")
  )

  pos_rank_fold_results <- list()
  neg_rank_fold_results <- list()
  synergy_rank_fold_results <- list()
  antagonism_rank_fold_results <- list()

  fold_InterXshift_results <- furrr::future_map(
    unique(data_internal$folds), function(fold_k) {

      fold_intxn_results <- fold_basis_results[[fold_k]]

      fold_positive_effects <- fold_intxn_results$top_positive_effects
      fold_negative_effects <- fold_intxn_results$top_negative_effects
      fold_synergy_effects <- fold_intxn_results$top_synergistic_interactions
      fold_antagonism_effects <- fold_intxn_results$top_antagonistic_interactions

      ## Calculate top positive marginal results:

      for (i in 1:nrow(fold_positive_effects)) {

        exposure <- fold_positive_effects$Variable[i]

        at <- data_internal[data_internal$folds != fold_k, ]
        av <- data_internal[data_internal$folds == fold_k, ]

        delta <- deltas[[exposure]]

        rank <- fold_positive_effects$Rank[i]

        lower_bound <- min(min(av[[exposure]]), min(at[[exposure]]))
        upper_bound <- max(max(av[[exposure]]), max(at[[exposure]]))

        if (density_classification) {
          ind_gn_exp_estim <- estimate_density_ratio(at = at, av = av, delta =  delta, var = exposure, covars = c(exposure, w_names), classifier = mu_learner)
        } else {


        ind_gn_exp_estim <- indiv_stoch_shift_est_g_exp(
          exposure = exposure,
          delta = delta,
          g_learner = pi_learner,
          covars = w_names,
          av = av,
          at = at,
          adaptive_delta = adaptive_delta,
          hn_trunc_thresh = hn_trunc_thresh,
          lower_bound = lower_bound,
          upper_bound = upper_bound,
          outcome_type = "continuous"
        )

        delta <- ind_gn_exp_estim$delta
        }



        covars <- c(a_names, w_names)

        ind_qn_estim <- indiv_stoch_shift_est_Q(
          exposure = exposure,
          delta = delta,
          mu_learner = mu_learner,
          covars = covars,
          av = av,
          at = at,
          lower_bound = lower_bound,
          upper_bound = upper_bound,
          outcome_type = outcome_type
        )

        Hn <- ind_gn_exp_estim$Hn_av

        tmle_fit <- tmle_exposhift(
          data_internal = av,
          delta = delta,
          Qn_scaled = ind_qn_estim$q_av,
          Qn_unscaled = scale_to_original(ind_qn_estim$q_av, min_orig = min(av$y), max_orig = max(av$y)),
          Hn = Hn,
          fluctuation = fluctuation,
          y = av$y
        )

        tmle_fit$call <- call

        pos_rank_shift_in_fold <- calc_final_ind_shift_param(
          tmle_fit,
          exposure,
          fold_k
        )

        pos_rank_shift_in_fold$Delta <- delta

        pos_rank_fold_results[[
          paste("Rank", rank, ":", exposure)
        ]] <- list(
          "data" = av,
          "Qn_scaled" = ind_qn_estim$q_av,
          "Hn" = Hn,
          "k_fold_result" = pos_rank_shift_in_fold,
          "Delta" = delta
        )

      }

      for (i in 1:nrow(fold_negative_effects)) {

        exposure <- fold_negative_effects$Variable[i]

        at <- data_internal[data_internal$folds != fold_k, ]
        av <- data_internal[data_internal$folds == fold_k, ]

        delta <- deltas[[exposure]]

        rank <- fold_positive_effects$Rank[i]


        lower_bound <- min(min(av[[exposure]]), min(at[[exposure]]))
        upper_bound <- max(max(av[[exposure]]), max(at[[exposure]]))


        if (density_classification) {
          ind_gn_exp_estim <- estimate_density_ratio(at = at, av = av, delta =  delta, var = exposure, covars = c(exposure, w_names), classifier = mu_learner)
        } else {


        ind_gn_exp_estim <- indiv_stoch_shift_est_g_exp(
          exposure = exposure,
          delta = delta,
          g_learner = pi_learner,
          covars = w_names,
          av = av,
          at = at,
          adaptive_delta = adaptive_delta,
          hn_trunc_thresh = hn_trunc_thresh,
          lower_bound = lower_bound,
          upper_bound = upper_bound,
          outcome_type = "continuous"
        )

        delta <- ind_gn_exp_estim$delta
        }

        covars <- c(a_names, w_names)

        ind_qn_estim <- indiv_stoch_shift_est_Q(
          exposure = exposure,
          delta = delta,
          mu_learner = mu_learner,
          covars = covars,
          av = av,
          at = at,
          lower_bound = lower_bound,
          upper_bound = upper_bound,
          outcome_type = outcome_type
        )

        Hn <- ind_gn_exp_estim$Hn_av

        tmle_fit <- tmle_exposhift(
          data_internal = av,
          delta = delta,
          Qn_scaled = ind_qn_estim$q_av,
          Qn_unscaled = scale_to_original(ind_qn_estim$q_av, min_orig = min(av$y), max_orig = max(av$y)),
          Hn = Hn,
          fluctuation = fluctuation,
          y = av$y
        )

        tmle_fit$call <- call

        neg_rank_shift_in_fold <- calc_final_ind_shift_param(
          tmle_fit,
          exposure,
          fold_k
        )

        neg_rank_shift_in_fold$Delta <- delta

        neg_rank_fold_results[[
          paste("Rank", rank, ":", exposure)
        ]] <- list(
          "data" = av,
          "Qn_scaled" = ind_qn_estim$q_av,
          "Hn" = Hn,
          "k_fold_result" = neg_rank_shift_in_fold,
          "Delta" = delta
        )

      }

      for (i in 1:nrow(fold_synergy_effects)) {

        at <- data_internal[data_internal$folds != fold_k, ]
        av <- data_internal[data_internal$folds == fold_k, ]

        exposure_1 <- fold_synergy_effects$Variable1[[i]]
        exposure_2 <- fold_synergy_effects$Variable2[[i]]

        exposures <- as.list(c(exposure_1, exposure_2))
        delta <- deltas[unlist(exposures)]
        exposures[[3]] <- unlist(exposures)

        rank <- fold_positive_effects$Rank[i]

        covars <- c(w_names)

        if (density_classification) {
          joint_gn_exp_estims <- list()
          for (i in 1:length(exposures)) {
            exposure_i <- exposures[[i]]
            delta_i <- deltas[exposure_i]
            results <- estimate_density_ratio(at = at, av = av, delta =  delta_i, var = exposure_i, covars = c(exposure_i, w_names), classifier = mu_learner)
            joint_gn_exp_estims[["Hn_results"]][[i]] <- results$Hn_av
          }

          deltas_updated <- delta
          deltas_updated[[3]] <- c(deltas_updated[[1]], deltas_updated[[2]])
        } else {

        joint_gn_exp_estims <- joint_stoch_shift_est_g_exp(
          exposures,
          deltas,
          g_learner = pi_learner,
          covars = covars,
          av = av,
          at = at,
          adaptive_delta = adaptive_delta,
          hn_trunc_thresh = hn_trunc_thresh,
          outcome_type = "continuous"
        )

        joint_gn_exp_estims$gn_results[[3]] <- mapply(
          `*`,
          joint_gn_exp_estims$gn_results[[1]],
          joint_gn_exp_estims$gn_results[[3]]
        )


        deltas_updated <- joint_gn_exp_estims$delta_results
        deltas_updated[[3]] <- c(deltas_updated[[1]], deltas_updated[[2]])
        }

        covars <- c(a_names, w_names)

        joint_qn_estims <- joint_stoch_shift_est_Q(
          exposures,
          deltas = deltas_updated,
          mu_learner = mu_learner,
          covars,
          av,
          at,
          outcome_type = outcome_type
        )

        intxn_results_list <- list()
        qn_estim_scaled_list <- list()
        joint_hn_estims <- list()

        for (i in 1:length(joint_qn_estims)) {
          hn_estim <- joint_gn_exp_estims$Hn_results[[i]]

          qn_estim <- joint_qn_estims[[i]]
          delta <- deltas_updated[[i]]

          tmle_fit <- tmle_exposhift(
            data_internal = av,
            delta = delta,
            Qn_scaled = qn_estim,
            Hn = hn_estim,
            fluctuation = fluctuation,
            y = av$y
          )

          intxn_results_list[[i]] <- tmle_fit
          joint_hn_estims[[i]] <- hn_estim
        }

        syngery_in_fold <- calc_final_joint_shift_param(
          joint_shift_fold_results = intxn_results_list,
          rank = rank,
          fold_k = fold_k,
          deltas_updated = deltas_updated,
          exposures = exposures
        )

        synergy_rank_fold_results[[
          paste("Rank", rank, ":", paste(exposure_1, exposure_2, sep = "-"))
        ]] <- list(
          "data" = av,
          "Qn_scaled" = joint_qn_estims,
          "Hn" = joint_hn_estims,
          "k_fold_result" = syngery_in_fold,
          "deltas" = deltas_updated
        )

      }

      for (i in 1:nrow(fold_antagonism_effects)) {

        at <- data_internal[data_internal$folds != fold_k, ]
        av <- data_internal[data_internal$folds == fold_k, ]

        exposure_1 <- fold_antagonism_effects$Variable1[[i]]
        exposure_2 <- fold_antagonism_effects$Variable2[[i]]

        exposures <- as.list(c(exposure_1, exposure_2))
        delta <- deltas[unlist(exposures)]
        exposures[[3]] <- unlist(exposures)

        rank <- fold_antagonism_effects$Rank[i]

        covars <- c(w_names)


        if (density_classification) {
          joint_gn_exp_estims <- list()
          for (i in 1:length(exposures)) {
            exposure_i <- exposures[[i]]
            delta_i <- deltas[exposure_i]
            results <- estimate_density_ratio(at = at, av = av, delta =  delta_i, var = exposure_i, covars = c(exposure_i, w_names), classifier = mu_learner)
            joint_gn_exp_estims[["Hn_results"]][[i]] <- results$Hn_av
          }

          deltas_updated <- delta
          deltas_updated[[3]] <- c(deltas_updated[[1]], deltas_updated[[2]])
        } else {

        joint_gn_exp_estims <- joint_stoch_shift_est_g_exp(
          exposures,
          deltas,
          g_learner = pi_learner,
          covars = covars,
          av = av,
          at = at,
          adaptive_delta = adaptive_delta,
          hn_trunc_thresh = hn_trunc_thresh,
          outcome_type = "continuous"
        )

        joint_gn_exp_estims$gn_results[[3]] <- mapply(
          `*`,
          joint_gn_exp_estims$gn_results[[1]],
          joint_gn_exp_estims$gn_results[[3]]
        )


        deltas_updated <- joint_gn_exp_estims$delta_results
        deltas_updated[[3]] <- c(deltas_updated[[1]], deltas_updated[[2]])

        }

        covars <- c(a_names, w_names)

        joint_qn_estims <- joint_stoch_shift_est_Q(
          exposures,
          deltas = deltas_updated,
          mu_learner = mu_learner,
          covars,
          av,
          at,
          outcome_type = outcome_type
        )

        intxn_results_list <- list()
        qn_estim_scaled_list <- list()
        joint_hn_estims <- list()

        for (i in 1:length(joint_qn_estims)) {
          hn_estim <- joint_gn_exp_estims$Hn_results[[i]]

          qn_estim <- joint_qn_estims[[i]]
          delta <- deltas_updated[[i]]

          tmle_fit <- tmle_exposhift(
            data_internal = av,
            delta = delta,
            Qn_scaled = qn_estim,
            Hn = hn_estim,
            fluctuation = fluctuation,
            y = av$y
          )

          intxn_results_list[[i]] <- tmle_fit
          joint_hn_estims[[i]] <- hn_estim
        }

        antagonism_in_fold <- calc_final_joint_shift_param(
          joint_shift_fold_results = intxn_results_list,
          rank = rank,
          fold_k = fold_k,
          deltas_updated = deltas_updated,
          exposures = exposures
        )

        antagonism_rank_fold_results[[
          paste("Rank", rank, ":", paste(exposure_1, exposure_2, sep = "-"))
        ]] <- list(
          "data" = av,
          "Qn_scaled" = joint_qn_estims,
          "Hn" = joint_hn_estims,
          "k_fold_result" = antagonism_in_fold,
          "deltas" = deltas_updated
        )

      }

      results_list <- list(
        pos_rank_fold_results,
        neg_rank_fold_results,
        synergy_rank_fold_results,
        antagonism_rank_fold_results
      )

      names(results_list) <- c(
        "top_positive_effects",
        "top_negative_effects",
        "top_synergy_effects",
        "top_antagonism_effects"
      )

      results_list
    },
    .options = furrr::furrr_options(seed = seed, packages = "InterXshift")
  )

  top_positive_results <- purrr::map(fold_InterXshift_results, c("top_positive_effects"))
  top_negative_results <- purrr::map(fold_InterXshift_results, c("top_negative_effects"))
  top_synergy_results <- purrr::map(fold_InterXshift_results, c("top_synergy_effects"))
  top_antagonism_results <- purrr::map(fold_InterXshift_results, c("top_antagonism_effects"))

  top_positive_results <- unlist(top_positive_results, recursive = FALSE)
  top_negative_results <- unlist(top_negative_results, recursive = FALSE)
  top_synergy_results <- unlist(top_synergy_results, recursive = FALSE)
  top_antagonism_results <- unlist(top_antagonism_results, recursive = FALSE)


  ## calculate pooled results based on rank
  pooled_pos_shift_results_rank <- calc_pooled_indiv_shifts(
    indiv_shift_results = top_positive_results,
    estimator = estimator,
    fluctuation = fluctuation,
    n_folds = n_folds,
    rank = TRUE
  )

  pooled_neg_shift_results_rank <- calc_pooled_indiv_shifts(
    indiv_shift_results = top_negative_results,
    estimator = estimator,
    fluctuation = fluctuation,
    n_folds = n_folds,
    rank = TRUE
  )

  pooled_synergy_shift_results_rank <- calc_pooled_intxn_shifts(
    intxn_shift_results = top_synergy_results,
    estimator = estimator,
    a_names = a_names,
    w_names = w_names,
    z_names = z_names,
    fluctuation = fluctuation,
    n_folds = n_folds,
    rank = TRUE
  )

  pooled_antagonism_shift_results_rank <- calc_pooled_intxn_shifts(
    intxn_shift_results = top_antagonism_results,
    estimator = estimator,
    a_names = a_names,
    w_names = w_names,
    z_names = z_names,
    fluctuation = fluctuation,
    n_folds = n_folds,
    rank = TRUE
  )


  ## calculate pooled results based on variable set
  pooled_pos_shift_results_varset <- calc_pooled_indiv_shifts(
    indiv_shift_results = top_positive_results,
    estimator = estimator,
    fluctuation = fluctuation,
    n_folds = n_folds,
    rank = FALSE
  )

  pooled_neg_shift_results_varset <- calc_pooled_indiv_shifts(
    indiv_shift_results = top_negative_results,
    estimator = estimator,
    fluctuation = fluctuation,
    n_folds = n_folds,
    rank = FALSE
  )

  pooled_synergy_shift_results_varset <- calc_pooled_intxn_shifts(
    intxn_shift_results = top_synergy_results,
    estimator = estimator,
    a_names = a_names,
    w_names = w_names,
    z_names = z_names,
    fluctuation = fluctuation,
    n_folds = n_folds,
    rank = FALSE
  )

  pooled_antagonism_shift_results_varset <- calc_pooled_intxn_shifts(
    intxn_shift_results = top_antagonism_results,
    estimator = estimator,
    a_names = a_names,
    w_names = w_names,
    z_names = z_names,
    fluctuation = fluctuation,
    n_folds = n_folds,
    rank = FALSE
  )



  results_list <- list(
    "Pos Shift Results by Rank" = pooled_pos_shift_results_rank,
    "Neg Shift Results by Rank" = pooled_neg_shift_results_rank,
    "Pos Shift Results by Exposure Set" = pooled_pos_shift_results_varset,
    "Neg Shift Results by Exposure Set" = pooled_neg_shift_results_varset,
    "K Fold Synergy Results" = pooled_synergy_shift_results_rank$k_fold_results,
    "Pooled Synergy Results by Rank" = pooled_synergy_shift_results_rank$pooled_results,
    "K Fold Antagonism Results" = pooled_antagonism_shift_results_rank$k_fold_results,
    "Pooled Antagonism Results by Rank" = pooled_antagonism_shift_results_rank$pooled_results,
    "Pooled Synergy Results by Exposure Set" = pooled_synergy_shift_results_varset$pooled_results,
    "Pooled Antagonism Results by Exposure Set" = pooled_antagonism_shift_results_varset$pooled_results
  )

  return(results_list)
}
