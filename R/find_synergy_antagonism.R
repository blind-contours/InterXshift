#' @title Data-adaptive Discovery of Interactions Based on Joint vs. Individual Shift Interventions
#'
#' @description The `find_synergy_antagonism` function provides a g-computation approach to finding interactions.
#' This implementation fits an SL and then predicts two way joint shifts in exposures and compares this to individual shifts,
#' ranks the interactions based on the difference and then this becomes the interactions we want to estimate
#' using CV-TMLE.
#'
#' @param data A \code{data.frame} containing all the variables needed for the analysis, including
#' baseline covariates, exposures, and the outcome.
#' @param deltas A named \code{list} or \code{vector} specifying the shift in exposures to define the target parameter.
#' Each element should correspond to an exposure variable specified in \code{a_names}, detailing the amount
#' by which that exposure is to be shifted.
#' @param a_names A \code{character} vector specifying the names of the exposure variables within \code{data}.
#' @param w_names A \code{character} vector specifying the names of the covariate variables within \code{data}.
#' @param outcome The name of the outcome variable in \code{data}.
#' @param outcome_type A \code{character} string indicating the type of the outcome variable; either "continuous",
#' "binary", or "count".
#' @param mu_learner A list of \code{\link[sl3]{Lrnr_sl}} learners specifying the ensemble machine learning models to be used
#' for outcome prediction within the Super Learner framework.
#' @param top_n An \code{integer} specifying the number of top positive and negative effects to return.
#' @param seed An \code{integer} value to set the seed for reproducibility.
#'
#' @return A list containing the top effects and interactions identified and estimated by the function.
#' It includes elements for top positive and negative individual effects as well as top synergistic and antagonistic interactions.
#'
#' @examples
#' \dontrun{
#' data <- data.frame(matrix(rnorm(100 * 10), ncol = 10))
#' names(data) <- c(paste0("X", 1:8), "exposure", "outcome")
#' deltas <- list(exposure = 0.1)
#' a_names <- "exposure"
#' w_names <- paste0("X", 1:8)
#' outcome <- "outcome"
#' outcome_type <- "continuous"
#' mu_learner <- list(sl3::Lrnr_mean$new(), sl3::Lrnr_glm$new())
#' top_n <- 3
#' seed <- 123
#'
#' results <- find_synergy_antagonism(data, deltas, a_names, w_names, outcome, outcome_type, mu_learner, top_n, seed)
#' print(results)
#' }
#'
#' @importFrom sl3 make_sl3_Task Lrnr_sl
#' @export

find_synergy_antagonism <- function(data, deltas, a_names, w_names, outcome, outcome_type, mu_learner, top_n = 3, seed) {
  future::plan(future::sequential, gc = TRUE)
  set.seed(seed)

  # Prepare Super Learner Task with specified data, exposures, and outcome
  task <- sl3::make_sl3_Task(
    data = data,
    covariates = c(a_names, w_names), # w_names assumed to be predefined
    outcome = outcome,
    outcome_type = outcome_type,
    folds = 5
  )

  # Train Super Learner
  sl <- sl3::Lrnr_sl$new(
    learners = mu_learner,
    metalearner = sl3::Lrnr_nnls$new()
  )

  sl_fit <- sl$train(task)

  # Initialize data frames for storing individual and interaction effects
  individual_effects_df <- data.frame(Variable = character(), Effect = numeric(), stringsAsFactors = FALSE)
  interaction_effects_df <- data.frame(Variable1 = character(), Variable2 = character(), Effect = numeric(), stringsAsFactors = FALSE)

  # Calculate individual effects
  for (var in a_names) {
    shifted_data <- as.data.frame(data)
    shifted_data[[var]] <- shifted_data[[var]] + deltas[[var]]
    predictions <- sl_fit$predict(sl3::make_sl3_Task(data = shifted_data, covariates = c(a_names, w_names), outcome = outcome))
    effect <- mean(predictions - data[[outcome]])
    individual_effects_df <- rbind(individual_effects_df, data.frame(Variable = var, Effect = effect))
  }

  # Rank individual effects
  individual_effects_df <- individual_effects_df[order(-individual_effects_df$Effect), ]
  top_positive_effects <- head(individual_effects_df, top_n)
  top_negative_effects <- tail(individual_effects_df, top_n)

  top_positive_effects$Rank <- seq(top_n)
  top_negative_effects$Rank <- rev(seq(top_n))

  # Calculate interaction effects
  for (indices in combn(seq_along(a_names), 2, simplify = FALSE)) {
    shifted_data <- as.data.frame(data)
    vars <- a_names[indices]
    shifted_data[vars] <- shifted_data[vars] + deltas[indices]
    joint_predictions <- sl_fit$predict(sl3::make_sl3_Task(data = shifted_data, covariates = c(a_names, w_names), outcome = outcome))
    individual_sum <- sum(individual_effects_df$Effect[individual_effects_df$Variable %in% vars])
    interaction_effect <- mean((joint_predictions - data[[outcome]]) - individual_sum)
    interaction_effects_df <- rbind(interaction_effects_df, data.frame(Variable1 = vars[1], Variable2 = vars[2], Effect = interaction_effect))
  }

  # Rank interaction effects
  interaction_effects_df <- interaction_effects_df[order(-interaction_effects_df$Effect), ]
  top_synergistic_interactions <- head(interaction_effects_df, top_n)
  top_antagonistic_interactions <- tail(interaction_effects_df, top_n)

  top_synergistic_interactions$Rank <- seq(top_n)
  top_antagonistic_interactions$Rank <- rev(seq(top_n))

  # Return the results
  return(list(
    top_positive_effects = top_positive_effects,
    top_negative_effects = top_negative_effects,
    top_synergistic_interactions = top_synergistic_interactions,
    top_antagonistic_interactions = top_antagonistic_interactions
  ))
}
