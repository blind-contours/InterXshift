fit_estimators <- function(w,
                           a,
                           y,
                           seed,
                           top_rank_effects,
                           m14_effect_truth,
                           m14_intxn_truth,
                           true_em_effects,
                           deltas,
                           cv_folds,
                           var_sets){

  sim_results <- InterXshift(
    w = w,
    a = a,
    y = y,
    delta = deltas,
    estimator = "tmle",
    fluctuation = "standard",
    n_folds = cv_folds,
    num_cores = 5,
    outcome_type = "continuous",
    parallel = TRUE,
    seed = seed,
    density_classification = TRUE,
    top_n = 1
  )


  top_rank_pos <- sim_results$`Pos Shift Results by Rank`$`Rank 1`
  joint_shift_results <- sim_results$`Pooled Synergy Results by Rank`$`Rank 1`

  indiv_biases <- vector()
  indiv_coverage <- vector()
  indiv_estimates <- vector()


  exposure_results_pooled <- top_rank_pos[top_rank_pos$Fold == "Pooled TMLE", ]
  pooled_exposure_est <- exposure_results_pooled$Psi
  indiv_bias <- top_rank_effects - pooled_exposure_est
  indiv_coverage <- ifelse(
      (exposure_results_pooled$`Lower CI` <= top_rank_effects &
         top_rank_effects <= exposure_results_pooled$`Upper CI`), 1, 0
    )

  indiv_estimates <- pooled_exposure_est



  indiv_results <- list(
    "indiv_est" = indiv_estimates,
    "indiv_bias" = indiv_bias,
    "indiv_cov" = indiv_coverage
  )

  pooled_tmle_joint <- joint_shift_results[joint_shift_results$Type == "Joint", ]
  pooled_tmle_intxn <- joint_shift_results[joint_shift_results$Type == "Interaction", ]

  intxn_bias <- pooled_tmle_intxn$Psi - m14_intxn_truth
  joint_bias <- pooled_tmle_joint$Psi - m14_effect_truth

  intxn_cov <- ifelse(
    (pooled_tmle_intxn$`Lower CI` <= m14_intxn_truth &
       m14_intxn_truth <= pooled_tmle_intxn$`Upper CI`), 1, 0
  )

  joint_cov <- ifelse(
    (pooled_tmle_joint$`Lower CI` <= m14_effect_truth &
       m14_effect_truth <= pooled_tmle_joint$`Upper CI`), 1, 0
  )

  joint_results <- list(
    "joint_bias" = joint_bias,
    "joint_est" = pooled_tmle_joint$Psi,
    "joint_cov" = joint_cov
  )

  intxn_results <- list(
    "intxn_bias" = intxn_bias,
    "intxn_est" = pooled_tmle_intxn$Psi,
    "intxn_cov" = intxn_cov
  )

  # bundle estimators in list
  estimates <- c(indiv_results, joint_results, intxn_results)

  return(estimates)
}


fit_estimators_mediation <- function(w,
                                     a,
                                     z,
                                     y,
                                     seed,
                                     nde_effects,
                                     nie_effects,
                                     ate_effects,
                                     deltas,
                                     cv_folds,
                                     num_cores,
                                     var_sets,
                                     exposure_quantized,
                                     mediator_quantized,
                                     n_mc_sample,
                                     density_type,
                                     integration_method = "MC",
                                     n_bins,
                                     use_multinomial){

  sim_results <- InterXshift(
    w = w,
    a = a,
    z = z,
    y = y,
    deltas = deltas,
    estimator = "tmle",
    fluctuation = "standard",
    n_folds = cv_folds,
    outcome_type = "continuous",
    quantile_thresh = 0,
    verbose = TRUE,
    parallel = TRUE,
    seed = seed,
    var_sets = var_sets,
    num_cores = num_cores,
    exposure_quantized = exposure_quantized,
    mediator_quantized = mediator_quantized,
    density_type = density_type,
    n_mc_sample = n_mc_sample,
    integration_method = integration_method,
    n_bins = n_bins,
    use_multinomial = use_multinomial
  )

  med_shift_results <- sim_results$`Mediation Shift Results`$az
  med_shift_results_pooled <- med_shift_results[med_shift_results$id == "Pooled", ]

  pseudo_reg_nde <- med_shift_results_pooled[med_shift_results_pooled$Parameter == "NDE-Pseudo-Reg"]
  pseudo_reg_nie <- med_shift_results_pooled[med_shift_results_pooled$Parameter == "NIE-Pseudo-Reg"]

  integrated_nde <- med_shift_results_pooled[med_shift_results_pooled$Parameter == "NDE-Integrated"]
  integrated_nie <- med_shift_results_pooled[med_shift_results_pooled$Parameter == "NIE-Integrated"]

  pooled_total <- med_shift_results_pooled[med_shift_results_pooled$Parameter == "Total-Pooled-TMLE"]

  pseudo_reg_nde_bias <- pseudo_reg_nde$Psi - nde_effects
  pseudo_reg_nie_bias <- pseudo_reg_nie$Psi - nie_effects

  integrated_nde_bias <- integrated_nde$Psi - nde_effects
  integrated_nie_bias <- integrated_nie$Psi - nie_effects

  ate_bias <- pooled_total$Psi - ate_effects

  ate_cov <- ifelse(
    (pooled_total$`Lower CI` <= ate_effects &
       ate_effects <= pooled_total$`Upper CI`), 1, 0
  )

  pseudo_reg_nde_cov <- ifelse(
    (pseudo_reg_nde$`Lower CI` <= nde_effects &
       nde_effects <= pseudo_reg_nde$`Upper CI`), 1, 0
  )

  pseudo_reg_nie_cov <- ifelse(
    (pseudo_reg_nie$`Lower CI` <= nie_effects &
       nie_effects <= pseudo_reg_nie$`Upper CI`), 1, 0
  )

  integrated_nde_cov <- ifelse(
    (integrated_nde$`Lower CI` <= nde_effects &
       nde_effects <= integrated_nde$`Upper CI`), 1, 0
  )

  integrated_nie_cov <- ifelse(
    (integrated_nie$`Lower CI` <= nie_effects &
       nie_effects <= integrated_nie$`Upper CI`), 1, 0
  )



  bias_results <- list(
    "pseudo_reg_nde_bias" = pseudo_reg_nde_bias,
    "pseudo_reg_nie_bias" = pseudo_reg_nie_bias,
    "integrated_nde_bias" = integrated_nde_bias,
    "integrated_nie_bias" = integrated_nie_bias,
    "ate_bias" = ate_bias
  )

  coverage_results <- list(
    "pseudo_reg_nde_cov" = pseudo_reg_nde_cov,
    "pseudo_reg_nie_cov" = pseudo_reg_nie_cov,
    "integrated_nde_cov" = integrated_nde_cov,
    "integrated_nie_cov" = integrated_nie_cov,
    "ate_cov" = ate_cov
  )

  estimate_results <- list(
    "pseudo_reg_nde_est" = pseudo_reg_nde$Psi,
    "pseudo_reg_nie_est" = pseudo_reg_nie$Psi,
    "integrated_nde_est" = integrated_nde$Psi,
    "integrated_nie_est" = integrated_nie$Psi,
    "ate_est" = pooled_total$Psi
  )

  # bundle estimators in list
  estimates <- c(bias_results, coverage_results, estimate_results)

  return(estimates)
}

fit_estimators_mediation_complicated <- function(w,
                                     a,
                                     z,
                                     y,
                                     seed,
                                     nde_a1,
                                     nie_a1,
                                     ate_a1,
                                     nde_a2,
                                     nie_a2,
                                     ate_a2,
                                     deltas,
                                     cv_folds,
                                     num_cores,
                                     var_sets,
                                     exposure_quantized,
                                     mediator_quantized,
                                     n_mc_sample,
                                     density_type,
                                     integration_method = "MC",
                                     n_bins,
                                     use_multinomial){

  sim_results <- InterXshift(
    w = w,
    a = a,
    z = z,
    y = y,
    deltas = deltas,
    estimator = "tmle",
    fluctuation = "standard",
    n_folds = cv_folds,
    outcome_type = "continuous",
    quantile_thresh = 0,
    verbose = TRUE,
    parallel = TRUE,
    seed = seed,
    var_sets = "a_1_quant-z_1_quant",
    num_cores = num_cores,
    quantize_exposure = exposure_quantized,
    mediator_quantized = mediator_quantized,
    density_type = density_type,
    n_mc_sample = n_mc_sample,
    integration_method = integration_method,
    n_bins = n_bins,
    use_multinomial = use_multinomial
  )

  med_shift_results <- sim_results$`Mediation Shift Results`$az
  med_shift_results_pooled <- med_shift_results[med_shift_results$id == "Pooled", ]

  pseudo_reg_nde <- med_shift_results_pooled[med_shift_results_pooled$Parameter == "NDE-Pseudo-Reg"]
  pseudo_reg_nie <- med_shift_results_pooled[med_shift_results_pooled$Parameter == "NIE-Pseudo-Reg"]

  integrated_nde <- med_shift_results_pooled[med_shift_results_pooled$Parameter == "NDE-Integrated"]
  integrated_nie <- med_shift_results_pooled[med_shift_results_pooled$Parameter == "NIE-Integrated"]

  pooled_total <- med_shift_results_pooled[med_shift_results_pooled$Parameter == "Total-Pooled-TMLE"]

  pseudo_reg_nde_bias <- pseudo_reg_nde$Psi - nde_effects
  pseudo_reg_nie_bias <- pseudo_reg_nie$Psi - nie_effects

  integrated_nde_bias <- integrated_nde$Psi - nde_effects
  integrated_nie_bias <- integrated_nie$Psi - nie_effects

  ate_bias <- pooled_total$Psi - ate_effects

  ate_cov <- ifelse(
    (pooled_total$`Lower CI` <= ate_effects &
       ate_effects <= pooled_total$`Upper CI`), 1, 0
  )

  pseudo_reg_nde_cov <- ifelse(
    (pseudo_reg_nde$`Lower CI` <= nde_effects &
       nde_effects <= pseudo_reg_nde$`Upper CI`), 1, 0
  )

  pseudo_reg_nie_cov <- ifelse(
    (pseudo_reg_nie$`Lower CI` <= nie_effects &
       nie_effects <= pseudo_reg_nie$`Upper CI`), 1, 0
  )

  integrated_nde_cov <- ifelse(
    (integrated_nde$`Lower CI` <= nde_effects &
       nde_effects <= integrated_nde$`Upper CI`), 1, 0
  )

  integrated_nie_cov <- ifelse(
    (integrated_nie$`Lower CI` <= nie_effects &
       nie_effects <= integrated_nie$`Upper CI`), 1, 0
  )



  bias_results <- list(
    "pseudo_reg_nde_bias" = pseudo_reg_nde_bias,
    "pseudo_reg_nie_bias" = pseudo_reg_nie_bias,
    "integrated_nde_bias" = integrated_nde_bias,
    "integrated_nie_bias" = integrated_nie_bias,
    "ate_bias" = ate_bias
  )

  coverage_results <- list(
    "pseudo_reg_nde_cov" = pseudo_reg_nde_cov,
    "pseudo_reg_nie_cov" = pseudo_reg_nie_cov,
    "integrated_nde_cov" = integrated_nde_cov,
    "integrated_nie_cov" = integrated_nie_cov,
    "ate_cov" = ate_cov
  )

  estimate_results <- list(
    "pseudo_reg_nde_est" = pseudo_reg_nde$Psi,
    "pseudo_reg_nie_est" = pseudo_reg_nie$Psi,
    "integrated_nde_est" = integrated_nde$Psi,
    "integrated_nie_est" = integrated_nie$Psi,
    "ate_est" = pooled_total$Psi
  )

  # bundle estimators in list
  estimates <- c(bias_results, coverage_results, estimate_results)

  return(estimates)
}


