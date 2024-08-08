AXES_LABEL_SIZE <- 16

#' @title Plot comparison MAP vs. robust MAP
#'
#' @description  Compare robust MAP prior to MAP prior
#'
#' @param rob_comp A dataframe containing values for the MAP and robust MAP prior
#' @param saf_topic Selected safety topic to analyze/the adverse event of interest
#' @param select_btrt selected backgroundtreatment
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#' @export
robust_map_prior_plot <- function(rob_comp, saf_topic, select_btrt, select_analysis) {
  # Plot densities
  p <- ggplot2::ggplot(rob_comp, ggplot2::aes(x = Probability, y = Value, color = Density, linetype = Density)) +
    ggplot2::geom_line(size = 1.65) +
    ggplot2::scale_linetype_manual(values = c("dotted", "dashed")) +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_continuous(labels = function(x) formatC(x, digits = 1, format = "f")) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = AXES_LABEL_SIZE)) +
    ggplot2::theme(text = ggplot2::element_text(size = AXES_LABEL_SIZE))
  if (select_analysis == "Incidence proportion") {
    p + ggplot2::labs(
      x = paste0("Proportion of Patients with ", saf_topic, " with treatment ", select_btrt),
      y = "Probability Density", color = "Density", linetype = "Density"
    )
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    p + ggplot2::labs(
      x = paste0("Log Scale of Incidence Rates of Patients with ", saf_topic, " with treatment ", select_btrt),
      y = "Probability Density", color = "Density", linetype = "Density"
    )
  }
}


#' @title Parametric Mix Display
#'
#' @description  Display parametric mixture density in the MAP Prior Tab
#'
#' @param param_approx map prior; best fitting mixture model from parametric_approx()
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#' @param saf_topic Selected safety topic to analyze/the adverse event of interest
#' @param select_btrt selected background treatment
#'
#' @export
param_mix_density_display <- function(param_approx, select_analysis, saf_topic, select_btrt) {
  # MAP prior
  mixture_mat <- data.frame(param_approx[, 1:ncol(param_approx)])
  str_vec <- vector(length = ncol(mixture_mat))

  # MAP Prior density function
  if (select_analysis == "Incidence proportion") {
    x <- seq(0.009, 1, length = 10000)
    for (j in 1:ncol(mixture_mat)) {
      if (j != ncol(mixture_mat)) {
        str_vec[j] <- paste0(mixture_mat[1, j], "*", "dbeta(x, shape1 = ", mixture_mat[2, j], ", shape2 = ", mixture_mat[3, j], ")", " + ")
      } else {
        str_vec[ncol(mixture_mat)] <- paste0(mixture_mat[1, ncol(mixture_mat)], "*", "dbeta(x, shape1 = ", mixture_mat[2, ncol(mixture_mat)], ", shape2 = ", mixture_mat[3, ncol(mixture_mat)], ")")
      }
    }
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    a <- RBesT::qmix(param_approx, 0.009)
    b <- RBesT::qmix(param_approx, 0.991)
    x <- seq(a, b, length = 10000)
    rm(a, b)
    for (j in 1:ncol(mixture_mat)) {
      if (j != ncol(mixture_mat)) {
        str_vec[j] <- paste0(mixture_mat[1, j], "*", "dnorm(x, mean = ", mixture_mat[2, j], ", sd = ", mixture_mat[3, j], ")", " + ")
      } else {
        str_vec[ncol(mixture_mat)] <- paste0(mixture_mat[1, ncol(mixture_mat)], "*", "dnorm(x, mean = ", mixture_mat[2, ncol(mixture_mat)], ", sd = ", mixture_mat[3, ncol(mixture_mat)], ")")
      }
    }
  }

  Prior <- eval(parse(text = paste(str_vec, collapse = "")))

  # Create dataframe for ggplot
  df <- data.frame(
    Probability = x,
    Density = rep("MAP Prior", each = 500),
    Value = Prior
  )
  # Plot density
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Probability, y = Value)) +
    ggplot2::geom_line(size = 1.65) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = AXES_LABEL_SIZE)) +
    ggplot2::theme(text = ggplot2::element_text(size = AXES_LABEL_SIZE))
  if (select_analysis == "Incidence proportion") {
    p + ggplot2::labs(x = paste0("Proportion of Patients with ", saf_topic, " with treatment ", select_btrt), y = "Probability Density") +
      ggplot2::scale_x_continuous(labels = function(x) paste0(formatC(x * 100, digits = 1, format = "f"), "%"))
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    p + ggplot2::labs(x = "Log Scale of Incidence Rates", y = "Probability Density")
  }
}


#' @title Forest plot display
#'
#' @description  Display Forestplot in the MAP Prior Tab
#'
#' @param map_object The gMAP element from RBesT
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#' @param saf_topic Selected safety topic to analyze/the adverse event of interest
#' @param select_btrt Selected background treatment
#'
#' @export
forest_plot_display <- function(map_object, select_analysis, saf_topic, select_btrt) {
  if (is.null(map_object)) {
    return(NULL)
  } else if (select_analysis == "Incidence proportion") {
    p <- plot(map_object)$forest_model + bayesplot::legend_move("right")
    p + ggplot2::ylab(paste0("Proportion of Patients with ", saf_topic, " with treatment ", select_btrt)) +
      ggplot2::scale_y_continuous(labels = function(x) paste0(formatC(x * 100, digits = 1, format = "f"), "%")) +
      ggplot2::theme(axis.text = ggplot2::element_text(size = AXES_LABEL_SIZE)) +
      ggplot2::theme(text = ggplot2::element_text(size = AXES_LABEL_SIZE))
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    p <- plot(map_object)$forest_model + bayesplot::legend_move("right")
    p + ggplot2::ylab(paste0("Hazard for ", saf_topic, " with treatment ", select_btrt)) +
      ggplot2::theme(axis.text = ggplot2::element_text(size = AXES_LABEL_SIZE)) +
      ggplot2::theme(text = ggplot2::element_text(size = AXES_LABEL_SIZE))
  }
}


#' @title Graphic assessment of Prior Data conflict
#'
#' @description  Prior data conflict assessment - compare prior, likelihood, and posterior in the NTA Tab
#'
#' @param new_trial_analysis Coordinates from the nta tab
#' @param saf_topic Selected safety topic to analyze/the adverse event of interest
#' @param select_btrt Selected backgroundtreatments
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#'
#' @export
nta_data_conflict_assassment_plot <- function(new_trial_analysis, saf_topic, select_btrt, select_analysis) {
  # Plot densities
  p <- ggplot2::ggplot(new_trial_analysis, ggplot2::aes(x = Probability, y = Value, color = Density, linetype = Density)) +
    ggplot2::geom_line(size = 1.65) +
    ggplot2::scale_linetype_manual(values = c("dotted", "dashed", "dotdash")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = AXES_LABEL_SIZE)) +
    ggplot2::theme(text = ggplot2::element_text(size = AXES_LABEL_SIZE))
  if (select_analysis == "Incidence proportion") {
    p + ggplot2::labs(x = paste0("Proportion of Patients with ", saf_topic, " with treatment ", select_btrt), y = "Probability Density") +
      ggplot2::scale_x_continuous(labels = function(x) paste0(formatC(x, digits = 1, format = "f"), "%")) +
      ggplot2::scale_x_continuous(labels = function(x) paste0(formatC(x * 100, digits = 1, format = "f"), "%"))
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    p + ggplot2::scale_x_continuous(labels = function(x) paste0(formatC(x, digits = 1, format = "f"), " ")) +
      ggplot2::labs(
        x = paste0("Log Scale of Incidence Rates of Patients with ", saf_topic, " with treatment ", select_btrt),
        y = "Probability Density"
      )
  }
}

#' @title Plot density
#'
#' @description  Plot density
#'
#' @param stat_inf_dist Data frame with statistical information of the distribution
#' @param ae_prop Vector of quantiles
#' @param saf_topic Selected safety topic to analyze/the adverse event of interest
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#' @param select_btrt selected background treatment
#'
#' @export
decision_making_density_plot <- function(select_analysis, stat_inf_dist, ae_prop, saf_topic, select_btrt) {
  if (select_analysis == "Incidence proportion") {
    ggplot2::ggplot(stat_inf_dist, ggplot2::aes(x = Probability, y = Value)) +
      ggplot2::geom_line(size = 1.65) +
      ggplot2::xlim(0, 1) +
      ggplot2::ylim(0, max(stat_inf_dist %>% dplyr::select(Value))) +
      ggplot2::geom_area(mapping = ggplot2::aes(x = ifelse(Probability >= ae_prop[1] & Probability <= ae_prop[2], Probability, 0)), fill = "salmon") +
      ggplot2::labs(x = paste0("Proportion of Patients with ", saf_topic, " with treatment ", select_btrt), y = "Probability Density") +
      ggplot2::theme_minimal() +
      ggplot2::scale_x_continuous(labels = function(x) paste0(formatC(x * 100, digits = 1, format = "f"), "%")) +
      ggplot2::theme(axis.text = ggplot2::element_text(size = AXES_LABEL_SIZE)) +
      ggplot2::theme(text = ggplot2::element_text(size = AXES_LABEL_SIZE))
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    ae_prop <- ae_prop * 100 # for text display they are devided by 100 in the shiny package
    stat_inf_dist_shade <- stat_inf_dist %>% filter(Probability >= ae_prop[1] &
      Probability <= ae_prop[2])


    ggplot2::ggplot(stat_inf_dist, ggplot2::aes(x = Probability, y = Value)) +
      ggplot2::geom_line(size = 1.65) +
      ggplot2::xlim(min(stat_inf_dist %>% dplyr::select(Probability)), max(stat_inf_dist %>% dplyr::select(Probability))) +
      ggplot2::ylim(0, max(stat_inf_dist %>% dplyr::select(Value))) +
      ggplot2::geom_area(data = stat_inf_dist_shade, aes(x = Probability, y = Value), fill = "salmon") +
      ggplot2::labs(x = paste0("log(incidence rate) for Patients with ", saf_topic, " with treatment ", select_btrt), y = "Probability Density") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text = ggplot2::element_text(size = AXES_LABEL_SIZE)) +
      ggplot2::theme(text = ggplot2::element_text(size = AXES_LABEL_SIZE))
  }
}


#' @title App plots
#'
#' @description Only includes plots generated in the app
#'
#' @param select_analysis event proportions or rates
#' @param map_object MAP Prior object
#' @param saf_topic Selected safety topic to analyze/the adverse event of interest
#' @param select_btrt selected treatment arm
#' @param new_trial_analysis nta
#' @param stat_inf_dist Data frame with statistical information of the distribution
#' @param mix mixture distribution object
#' @export
app_plots <- function(select_analysis, map_object, saf_topic, select_btrt, new_trial_analysis, stat_inf_dist, mix) {
  grDevices::pdf(file)

  if (select_analysis == "Incidence proportion") {
    p <- plot(map_object$forest_model + bayesplot::legend_move("right"))
    print(p + ggplot2::ylab(paste0("Proportion of Patients with ", saf_topic)))
  }
  print(ggplot2::ggplot(new_trial_analysis, ggplot2::aes(x = Probability, y = Value, color = Density, linetype = Density)) +
    ggplot2::geom_line(size = 1.65) +
    ggplot2::labs(x = paste0("Proportion of Patients with ", saf_topic, " with treatment ", select_btrt), y = "Probability Density", color = "Density", linetype = "Density") +
    ggplot2::scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
    ggplot2::theme_minimal())

  print(ggplot2::ggplot(stat_inf_dist, ggplot2::aes(x = Probability, y = Value)) +
    ggplot2::geom_line(size = 1.65) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, max(stat_inf_dist %>% dplyr::select(Value))) +
    ggplot2::geom_area(mapping = ggplot2::aes(x = ifelse(Probability >= RBesT::qmix(mix, 0.10, lower.tail = TRUE), Probability, 0)), fill = "salmon") +
    ggplot2::labs(x = paste0("Proportion of Patients with ", saf_topic, " with treatment ", select_btrt), y = "Probability Density") +
    ggplot2::theme_minimal())

  print(ggplot2::ggplot(stat_inf_dist, ggplot2::aes(x = Probability, y = Value)) +
    ggplot2::geom_line(size = 1.65) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, max(stat_inf_dist %>% dplyr::select(Value))) +
    ggplot2::geom_area(mapping = ggplot2::aes(x = ifelse(Probability >= RBesT::qmix(mix, 0.05, lower.tail = TRUE), Probability, 0)), fill = "salmon") +
    ggplot2::labs(x = paste0("Proportion of Patients with ", saf_topic, " with treatment ", select_btrt), y = "Probability Density") +
    ggplot2::theme_minimal())

  print(ggplot2::ggplot(stat_inf_dist, ggplot2::aes(x = Probability, y = Value)) +
    ggplot2::geom_line(size = 1.65) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, max(stat_inf_dist %>% dplyr::select(Value))) +
    ggplot2::geom_area(mapping = ggplot2::aes(x = ifelse(Probability >= RBesT::qmix(mix, 0.01, lower.tail = TRUE), Probability, 0)), fill = "salmon") +
    ggplot2::labs(x = paste0("Proportion of Patients with ", saf_topic, " with treatment ", select_btrt), y = "Probability Density") +
    ggplot2::theme_minimal())

  grDevices::dev.off()
}
