#' @title Posterior Summary Statistics for Safety Topics of Interest
#'
#' @description Creates a table of posterior summary statistics for each
#' safety topic of interest. Creates the arrays fot thr pdf.
#'
#' @param input_data dataframe from data_table_prep() including the new trial
#' @param saf_topic the safety topic
#' @param seed a seed
#' @param cb_list_ctrl expects List with control Arm indicators
#' @param cb_list_trt expects List with treatment Arm indicators
#' @param explore for exploring prior data conflict, expected by posterior_dist()
#'
#' @return a table of posterior summary statistics for each event type for rates and proportions
#' @export
ae_summary_table <-
  function(input_data,
           cb_list_ctrl,
           cb_list_trt,
           saf_topic,
           explore = FALSE,
           seed = NULL) {
    if (length(cb_list_ctrl) != length(cb_list_trt)) {
      stop("Same amount of compared groups necessary.")
    }

    # setting up overall variables

    n_group <- length(cb_list_ctrl)
    tau <- tau_adjust(
      select_analysis = "Incidence proportion",
      hist_borrow = "Large"
    )
    # boundaries for credible intervals
    lb <- 0.025
    ub <- 0.975

    # sample size for posterior comparison
    sample_size <- 1000

    if (is.null(seed)) {
      seed <- as.numeric(Sys.time())
    }

    # Initiating the arrays for the download tables

    array_inci <-
      array(
        data = NA,
        dim = c(length(saf_topic), 7, 2, n_group),
        dimnames = list(
          saf_topic,
          c("r", "n", "%", "post_mean", "post_median", "cri_l", "cri_u"),
          c("Arm A", "Arm B"),
          group = 1:n_group
        )
      )

    array_inci_comp <- array(
      data = NA,
      dim = c(length(saf_topic), 6, 2, n_group),
      dimnames = list(
        saf_topic,
        c("naive", "post_mean", "post_median", "cri_l", "cri_u", "A>B%"),
        c("Risk Diff", "Risk Ratio"),
        group = 1:n_group
      )
    )

    array_er <-
      array(
        data = NA,
        dim = c(length(saf_topic), 5, 2, n_group),
        dimnames = list(
          saf_topic,
          c("naive", "post_mean", "post_median", "cri_l", "cri_u"),
          c("Arm A", "Arm B"),
          group = 1:n_group
        )
      )

    array_er_comp <-
      array(
        data = NA,
        dim = c(length(saf_topic), 6, 2, n_group),
        dimnames = list(
          saf_topic,
          c("naive", "post_mean", "post_median", "cri_l", "cri_u", "A>B%"),
          c("Risk Diff", "Risk Ratio"),
          group = 1:n_group
        )
      )

    ID_ctr_list <- list(grp1 = NULL, grp2 = NULL, grp3 = NULL, grp4 = NULL)
    ID_trt_list <- list(grp1 = NULL, grp2 = NULL, grp3 = NULL, grp4 = NULL)
    for (group in 1:n_group) {
      trt_arm <- input_data %>% dplyr::filter(ARM == cb_list_trt[[group]])
      ID_trt_list[[group]] <- unique(trt_arm$STUDYID)
      ctr_arm <- input_data %>% dplyr::filter(ARM == cb_list_ctrl[[group]])
      ID_ctr_list[[group]] <- unique(ctr_arm$STUDYID)
      for (topic in saf_topic) {
        trt_trials <- data_table_prep(
          input_data = input_data,
          select_analysis = "Incidence proportion",
          saf_topic = topic,
          select_btrt = cb_list_trt[[group]]
        )

        ctr_trials <- data_table_prep(
          input_data = input_data,
          select_analysis = "Incidence proportion",
          saf_topic = topic,
          select_btrt = cb_list_ctrl[[group]]
        )

        array_inci[topic, "r", "Arm A", group] <- sum(trt_trials$N_WITH_AE)
        array_inci[topic, "r", "Arm B", group] <- sum(ctr_trials$N_WITH_AE)
        array_inci[topic, "n", "Arm A", group] <- sum(trt_trials$N)
        array_inci[topic, "n", "Arm B", group] <- sum(ctr_trials$N)
        array_inci[topic, "%", "Arm A", group] <- array_inci[topic, "r", "Arm A", group] / array_inci[topic, "n", "Arm A", group]
        array_inci[topic, "%", "Arm B", group] <- array_inci[topic, "r", "Arm B", group] / array_inci[topic, "n", "Arm B", group]

        ### Binomial
        # weakly informative prior Be(1,1)
        weak_inf_bin <- RBesT::mixbeta(c(1, 1, 1))
        robust_weight <- 0.2
        robust_mean_bin <- 0.5

        hist_trt_trials <- trt_trials %>% dplyr::filter(HIST == 1)
        hist_ctr_trials <- ctr_trials %>% dplyr::filter(HIST == 1)

        ### Arm A #Binomial

        if (nrow(hist_trt_trials) > 0) {
          trt_map_prior <- map_prior_func(
            input_data      = hist_trt_trials,
            select_analysis = "Incidence proportion",
            tau_dist        = "HalfNormal",
            adj_tau         = tau,
            seed            = seed
          )
          trt_approx <- parametric_approx("Incidence proportion", trt_map_prior)
        } else {
          trt_approx <- weak_inf_bin
        }


        trt_rob_prior <- RBesT::robustify(trt_approx, weight = robust_weight, mean = robust_mean_bin)

        trt_current_trial <- trt_trials %>% dplyr::filter(HIST == 0)
        trt_new_n <- sum(trt_current_trial$N)
        trt_new_r <- sum(trt_current_trial$N_WITH_AE)

        trt_post <-
          posterior_dist(
            select_analysis = "Incidence proportion",
            robust_map_prior = trt_rob_prior,
            new_v1 = trt_new_n,
            new_v2 = trt_new_r,
            seed = seed
          )

        array_inci[topic, "cri_l", "Arm A", group] <- RBesT::qmix(trt_post, lb)
        array_inci[topic, "cri_u", "Arm A", group] <- RBesT::qmix(trt_post, ub)
        array_inci[topic, "post_median", "Arm A", group] <- RBesT::qmix(trt_post, 0.5)
        array_inci[topic, "post_mean", "Arm A", group] <- sum(trt_post[1, ] * (trt_post[2, ] / (trt_post[2, ] + trt_post[3, ])))

        ## Arm B # Binomial

        if (nrow(hist_ctr_trials) > 0) {
          ctr_map_prior <- map_prior_func(
            input_data = hist_ctr_trials,
            select_analysis = "Incidence proportion",
            tau_dist = "HalfNormal",
            adj_tau = tau,
            seed = seed
          )
          ctr_approx <- parametric_approx("Incidence proportion", ctr_map_prior)
        } else {
          # ctr_map_prior <- weak_inf_bin
          ctr_approx <- weak_inf_bin
        }

        ctr_rob_prior <-
          RBesT::robustify(ctr_approx, weight = robust_weight, mean = robust_mean_bin)

        ctr_current_trial <- ctr_trials %>% dplyr::filter(HIST == 0)
        ctr_new_n <- sum(ctr_current_trial$N)
        ctr_new_r <- sum(ctr_current_trial$N_WITH_AE)

        ctr_post <-
          posterior_dist(
            select_analysis = "Incidence proportion",
            robust_map_prior = ctr_rob_prior,
            explore = FALSE,
            new_v1 = ctr_new_n,
            new_v2 = ctr_new_r,
            seed = seed
          )
        array_inci[topic, "cri_l", "Arm B", group] <- RBesT::qmix(ctr_post, lb)
        array_inci[topic, "cri_u", "Arm B", group] <- RBesT::qmix(ctr_post, ub)
        array_inci[topic, "post_median", "Arm B", group] <- RBesT::qmix(ctr_post, 0.5)
        array_inci[topic, "post_mean", "Arm B", group] <- sum(ctr_post[1, ] * (ctr_post[2, ] / (ctr_post[2, ] + ctr_post[3, ])))

        # propotion posterior samples
        trt_pr_post_samp <- RBesT::rmix(trt_post, sample_size)
        ctr_pr_post_samp <- RBesT::rmix(ctr_post, sample_size)

        pr_diffe_samp <- trt_pr_post_samp - ctr_pr_post_samp
        pr_ratio_samp <- trt_pr_post_samp / ctr_pr_post_samp

        array_inci_comp[topic, "naive", "Risk Diff", group] <- as.numeric(array_inci[topic, "%", "Arm A", group]) - as.numeric(array_inci[topic, "%", "Arm B", group])
        array_inci_comp[topic, "post_mean", "Risk Diff", group] <- mean(pr_diffe_samp)
        array_inci_comp[topic, "post_median", "Risk Diff", group] <- median(pr_diffe_samp)
        array_inci_comp[topic, "cri_l", "Risk Diff", group] <- quantile(pr_diffe_samp, lb)
        array_inci_comp[topic, "cri_u", "Risk Diff", group] <- quantile(pr_diffe_samp, ub)
        array_inci_comp[topic, "A>B%", "Risk Diff", group] <- sum(pr_diffe_samp > 0) / sample_size

        array_inci_comp[topic, "naive", "Risk Ratio", group] <- as.numeric(array_inci[topic, "%", "Arm A", group]) / as.numeric(array_inci[topic, "%", "Arm B", group])
        array_inci_comp[topic, "post_mean", "Risk Ratio", group] <- mean(pr_ratio_samp)
        array_inci_comp[topic, "post_median", "Risk Ratio", group] <- median(pr_ratio_samp)
        array_inci_comp[topic, "cri_l", "Risk Ratio", group] <- quantile(pr_ratio_samp, lb)
        array_inci_comp[topic, "cri_u", "Risk Ratio", group] <- quantile(pr_ratio_samp, ub)
        array_inci_comp[topic, "A>B%", "Risk Ratio", group] <- sum(pr_diffe_samp > 0) / sample_size

        ### Event Rates

        trt_trials <- data_table_prep(
          input_data = input_data,
          select_analysis = "Exposure-adjusted AE rate",
          saf_topic = topic,
          select_btrt = cb_list_trt[[group]]
        )

        ctr_trials <- data_table_prep(
          input_data = input_data,
          select_analysis = "Exposure-adjusted AE rate",
          saf_topic = topic,
          select_btrt = cb_list_ctrl[[group]]
        )

        tau_er <- tau_adjust(
          select_analysis = "Exposure-adjusted AE rate",
          hist_borrow = "Large"
        )

        hist_trt_trials <- trt_trials %>% dplyr::filter(HIST == 1)
        trt_current_trial <- trt_trials %>% dplyr::filter(HIST == 0)

        # For a reasonable posterioe historical information should be available,
        # as well as new data.
        trt_post_check <- nrow(trt_current_trial) > 0 && nrow(hist_trt_trials) > 0

        if (trt_post_check == TRUE) {
          ### Arm A #Exponential

          er_trt_prior <-
            RBesT::gMAP(N_WITH_AE ~ 1 + offset(log(TOT_EXP)) | STUDYID,
              data = hist_trt_trials,
              tau.dist = "LogNormal",
              tau.prior = tau_er, ## assuming moderate heterogeniety
              beta.prior = 1,
              family = "poisson"
            )

          er_trt_prior_fit <- parametric_approx("Exposure-adjusted AE rate", er_trt_prior)
          robust_mean_er <- summary(er_trt_prior_fit)["mean"]

          # found error
          trt_rob_prior <- RBesT::robustify(er_trt_prior_fit, weight = robust_weight, mean = robust_mean_er)

          trt_new_nwae <- sum(trt_current_trial$N_WITH_AE)
          trt_new_texp <- sum(trt_current_trial$TOT_EXP)

          trt_post <-
            posterior_dist(
              select_analysis = "Exposure-adjusted AE rate",
              robust_map_prior = trt_rob_prior,
              new_v1 = trt_new_nwae,
              new_v2 = trt_new_texp,
              seed = seed
            )

          array_er[topic, "naive", "Arm A", group] <- sum(trt_trials$N_WITH_AE) / sum(trt_trials$TOT_EXP)
          desc_stats <- summary(trt_post)
          # stats_mat_post_a  <- data.frame(desc_stats$theta.pred)
          array_er[topic, "post_mean", "Arm A", group] <- desc_stats["mean"]
          array_er[topic, "post_median", "Arm A", group] <- RBesT::qmix(trt_post, 0.5)
          array_er[topic, "cri_l", "Arm A", group] <- RBesT::qmix(trt_post, lb)
          array_er[topic, "cri_u", "Arm A", group] <- RBesT::qmix(trt_post, ub)
        } else {
          array_er[topic, "naive", "Arm A", group] <- NA
          array_er[topic, "post_mean", "Arm A", group] <- NA
          array_er[topic, "post_median", "Arm A", group] <- NA
          array_er[topic, "cri_l", "Arm A", group] <- NA
          array_er[topic, "cri_u", "Arm A", group] <- NA
        }

        ## Arm B # Exponential
        hist_ctr_trials <- ctr_trials %>% dplyr::filter(HIST == 1)
        ctr_current_trial <- ctr_trials %>% dplyr::filter(HIST == 0)
        ctr_post_check <- nrow(ctr_current_trial) > 0 && nrow(hist_ctr_trials) > 0

        if (ctr_post_check == TRUE) {
          er_ctr_prior <-
            RBesT::gMAP(N_WITH_AE ~ 1 + offset(log(TOT_EXP)) | STUDYID,
              data = ctr_trials,
              tau.dist = "LogNormal",
              tau.prior = tau_er, ## assuming moderate heterogeniety
              beta.prior = 1,
              family = "poisson"
            )

          er_ctr_prior_fit <- parametric_approx("Exposure-adjusted AE rate", er_ctr_prior)
          robust_mean_er <- summary(er_ctr_prior_fit)["mean"]

          ctr_rob_prior <- RBesT::robustify(er_ctr_prior_fit, weight = robust_weight, mean = robust_mean_er)

          ctr_new_nwae <- sum(ctr_current_trial$N_WITH_AE)
          ctr_new_texp <- sum(ctr_current_trial$TOT_EXP)

          ctr_post <-
            posterior_dist(
              select_analysis = "Exposure-adjusted AE rate",
              robust_map_prior = ctr_rob_prior,
              new_v1 = ctr_new_nwae,
              new_v2 = ctr_new_texp,
              seed = seed
            )

          array_er[topic, "naive", "Arm B", group] <- sum(ctr_trials$N_WITH_AE) / sum(ctr_trials$TOT_EXP)
          desc_stats <- summary(ctr_post)
          array_er[topic, "post_mean", "Arm B", group] <- desc_stats["mean"]
          array_er[topic, "post_median", "Arm B", group] <- RBesT::qmix(ctr_post, 0.5)
          array_er[topic, "cri_l", "Arm B", group] <- RBesT::qmix(ctr_post, lb)
          array_er[topic, "cri_u", "Arm B", group] <- RBesT::qmix(ctr_post, ub)
        } else {
          array_er[topic, "naive", "Arm B", group] <- NA
          array_er[topic, "post_mean", "Arm B", group] <- NA
          array_er[topic, "post_median", "Arm B", group] <- NA
          array_er[topic, "cri_l", "Arm B", group] <- NA
          array_er[topic, "cri_u", "Arm B", group] <- NA
        }


        comp_check <- ctr_post_check == TRUE && trt_post_check == TRUE

        if (comp_check == TRUE) {
          # propotion posterior samples
          trt_pr_post_samp <- RBesT::rmix(trt_post, sample_size)
          ctr_pr_post_samp <- RBesT::rmix(ctr_post, sample_size)

          er_diffe_samp <- trt_pr_post_samp - ctr_pr_post_samp
          er_ratio_samp <- trt_pr_post_samp / ctr_pr_post_samp

          array_er_comp[topic, "naive", "Risk Diff", group] <- as.numeric(array_er[topic, "naive", "Arm A", group]) - as.numeric(array_er[topic, "naive", "Arm B", group])
          array_er_comp[topic, "post_mean", "Risk Diff", group] <- mean(er_diffe_samp)
          array_er_comp[topic, "post_median", "Risk Diff", group] <- median(er_diffe_samp)
          array_er_comp[topic, "cri_l", "Risk Diff", group] <- quantile(er_diffe_samp, lb)
          array_er_comp[topic, "cri_u", "Risk Diff", group] <- quantile(er_diffe_samp, ub)
          array_er_comp[topic, "A>B%", "Risk Diff", group] <- sum(er_diffe_samp > 0) / sample_size

          array_er_comp[topic, "naive", "Risk Ratio", group] <- as.numeric(array_er[topic, "naive", "Arm A", group]) / as.numeric(array_er[topic, "naive", "Arm B", group])
          array_er_comp[topic, "post_mean", "Risk Ratio", group] <- mean(er_ratio_samp)
          array_er_comp[topic, "post_median", "Risk Ratio", group] <- median(er_ratio_samp)
          array_er_comp[topic, "cri_l", "Risk Ratio", group] <- quantile(er_ratio_samp, lb)
          array_er_comp[topic, "cri_u", "Risk Ratio", group] <- quantile(er_ratio_samp, ub)
          array_er_comp[topic, "A>B%", "Risk Ratio", group] <- sum(er_diffe_samp > 0) / sample_size
        } else {
          array_er_comp[topic, "naive", "Risk Diff", group] <- NA
          array_er_comp[topic, "post_mean", "Risk Diff", group] <- NA
          array_er_comp[topic, "post_median", "Risk Diff", group] <- NA
          array_er_comp[topic, "cri_l", "Risk Diff", group] <- NA
          array_er_comp[topic, "cri_u", "Risk Diff", group] <- NA
          array_er_comp[topic, "A>B%", "Risk Diff", group] <- NA

          array_er_comp[topic, "naive", "Risk Ratio", group] <- NA
          array_er_comp[topic, "post_mean", "Risk Ratio", group] <- NA
          array_er_comp[topic, "post_median", "Risk Ratio", group] <- NA
          array_er_comp[topic, "cri_l", "Risk Ratio", group] <- NA
          array_er_comp[topic, "cri_u", "Risk Ratio", group] <- NA
          array_er_comp[topic, "A>B%", "Risk Ratio", group] <- NA
        }
      }
    }

    return(list(
      InciProp = array_inci, InciPropComp = array_inci_comp,
      EventRates = array_er, EvenRatesComp = array_er_comp,
      ID_ctr_list = ID_ctr_list, ID_trt_list = ID_trt_list,
      cb_list_ctrl = cb_list_ctrl, cb_list_trt = cb_list_trt
    ))
  }

# Display summary stats of robust MAP prior and MAP prior
#' Title
#'
#' @param map_object
#' @param param_approx
#' @param ess_method
#' @param robust_map_object
#' @param rob_ess_method
#' @param download
#' @export
summary_stats_robust_map_prior_display <- function(map_object, select_analysis, param_approx, ess_method, robust_map_object, rob_ess_method, download = FALSE) {
  if (select_analysis == "Incidence proportion") {
    # Summary statistics for MAP prior
    desc_stats <- summary(map_object)
    stats_mat <- data.frame(desc_stats$theta.pred)

    Prior_mean <- paste0(formatC(stats_mat[1, 1] * 100, digits = 2, format = "f"), "%")
    Prior_sd <- paste0(formatC(stats_mat[1, 2] * 100, digits = 2, format = "f"), "%")
    Prior_median <- paste0(formatC(stats_mat[1, 4] * 100, digits = 2, format = "f"), "%")
    LB <- as.numeric(format(unlist(stats_mat$X2.5)))
    LB_perc <- formatC(LB * 100, digits = 2, format = "f")
    UB <- as.numeric(format(unlist(stats_mat$X97.5)))
    UB_perc <- formatC(UB * 100, digits = 2, format = "f")
    Prior_95ci <- paste0("[", LB_perc, "%, ", UB_perc, "%]")
    Prior_ess <- round(RBesT::ess(param_approx, method = ess_method))
    # Summary statistics for robust MAP prior
    rob_mixture_mat <- data.frame(robust_map_object[, 1:ncol(robust_map_object)])
    rob_vec1 <- vector(length = ncol(rob_mixture_mat)) # Robust MAP prior mean
    rob_vec2 <- vector(length = ncol(rob_mixture_mat)) # Robust MAP prior sd
    rob_vec3 <- vector(length = ncol(rob_mixture_mat)) # Robust MAP prior median
    # Robust MAP prior mean
    for (j in 1:ncol(rob_mixture_mat)) {
      rob_vec1[j] <- rob_mixture_mat[1, j] * (rob_mixture_mat[2, j] / (rob_mixture_mat[2, j] + rob_mixture_mat[3, j]))
    }
    # Robust MAP prior variance
    for (j in 1:ncol(rob_mixture_mat)) {
      rob_vec2[j] <- rob_mixture_mat[1, j] * ((rob_mixture_mat[2, j] * rob_mixture_mat[3, j]) / (((rob_mixture_mat[2, j] + rob_mixture_mat[3, j])^2) * (rob_mixture_mat[2, j] + rob_mixture_mat[3, j] + 1)))
    }
    # Robust MAP prior median
    for (j in 1:ncol(rob_mixture_mat)) {
      rob_vec3[j] <- rob_mixture_mat[1, j] * ((rob_mixture_mat[2, j] - (1 / 3)) / (rob_mixture_mat[2, j] + rob_mixture_mat[3, j] - (2 / 3)))
    }
    rob_Prior_mean <- paste0(formatC(sum(rob_vec1) * 100, digits = 2, format = "f"), "%")
    rob_Prior_sd <- paste0(formatC(sqrt(sum(rob_vec2)) * 100, digits = 2, format = "f"), "%")
    rob_Prior_median <- paste0(formatC(sum(rob_vec3) * 100, digits = 2, format = "f"), "%")
    rob_Prior_95ci <- paste0("[", formatC(RBesT::qmix(robust_map_object, 0.025) * 100, digits = 2, format = "f"), "%, ", formatC(RBesT::qmix(robust_map_object, 0.975) * 100, digits = 2, format = "f"), "%]")
    rob_Prior_ess <- round(RBesT::ess(robust_map_object, method = rob_ess_method))
    sum_stats <- data.frame(
      Mean = c(Prior_mean, rob_Prior_mean),
      SD = c(Prior_sd, rob_Prior_sd),
      Median = c(Prior_median, rob_Prior_median),
      conf_int = c(Prior_95ci, rob_Prior_95ci),
      ESS = c(Prior_ess, rob_Prior_ess)
    )
    if (download == FALSE) {
      rownames(sum_stats) <- c("MAP Prior: log(hazard)", "Robust MAP Prior")
      sum_stats
      sum_stats %>%
        dplyr::rename("95% CrI" = conf_int)
    } else {
      sum_stats
    }
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    # Rbest gives the model alreday out for exp(theta_pred), but we need
    # theta_pred for further analysis
    theta_fit <- as.data.frame(map_object$fit)
    theta_pred <- theta_fit$theta_pred
    theta_resp <- theta_fit$theta_resp_pred

    stats_mat_prop <- array(data = NA, dim = c(4, 5), dimnames = list(
      c(
        "MAP Prior: log(hazard)", "robust MAP Prior: log(hazard)",
        "MAP Prior: hazard", "robust MAP Prior: hazard"
      ),
      c("Mean", "SD", "Median", "95% CrI", "ESS")
    ))
    stats_mat_prop <- as.data.frame(stats_mat_prop)

    stats_mat_prop["MAP Prior: log(hazard)", "Mean"] <- round(mean(theta_pred), 2)
    stats_mat_prop["MAP Prior: log(hazard)", "SD"] <- round(sd(theta_pred), 2)
    stats_mat_prop["MAP Prior: log(hazard)", "Median"] <- round(median(theta_pred), 2)
    stats_mat_prop["MAP Prior: log(hazard)", "95% CrI"] <- paste0(
      "[",
      round(quantile(theta_pred, 0.25), 2), " , ",
      round(quantile(theta_pred, 0.75), 2), "]"
    )
    stats_mat_prop["MAP Prior: log(hazard)", "ESS"] <- round(RBesT::ess(param_approx, method = ess_method, sigma = 1))


    stats_mat_prop["MAP Prior: hazard", "Mean"] <- round(mean(theta_resp), 2)
    stats_mat_prop["MAP Prior: hazard", "SD"] <- round(sd(theta_resp), 2)
    stats_mat_prop["MAP Prior: hazard", "Median"] <- round(median(theta_resp), 2)
    stats_mat_prop["MAP Prior: hazard", "95% CrI"] <- paste0(
      "[",
      round(quantile(theta_resp, 0.25), 2), " , ",
      round(quantile(theta_resp, 0.75), 2), "]"
    )
    stats_mat_prop["MAP Prior: hazard", "ESS"] <- c("Not applicable.")


    rob_Prior_mean <- round(mean(RBesT::rmix(robust_map_object, 1000)), digits = 2)
    rob_Prior_sd <- round(sd(RBesT::rmix(robust_map_object, 1000)), digits = 2)
    rob_Prior_median <- round(RBesT::qmix(robust_map_object, 0.5), digits = 2)
    rob_Prior_95ci <- paste0(
      "[", formatC(RBesT::qmix(robust_map_object, 0.025),
        digits = 2, format = "f"
      ), " , ",
      formatC(RBesT::qmix(robust_map_object, 0.975),
        digits = 2, format = "f"
      ), "]"
    )
    rob_Prior_ess <- round(RBesT::ess(robust_map_object, method = ess_method, sigma = 1))


    stats_mat_prop["robust MAP Prior: log(hazard)", "Mean"] <- rob_Prior_mean
    stats_mat_prop["robust MAP Prior: log(hazard)", "SD"] <- rob_Prior_sd
    stats_mat_prop["robust MAP Prior: log(hazard)", "Median"] <- rob_Prior_median
    stats_mat_prop["robust MAP Prior: log(hazard)", "95% CrI"] <- rob_Prior_95ci
    stats_mat_prop["robust MAP Prior: log(hazard)", "ESS"] <- round(RBesT::ess(robust_map_object, method = ess_method, sigma = 1))


    stats_mat_prop["robust MAP Prior: hazard", "Mean"] <- round(exp(rob_Prior_mean), 2)
    stats_mat_prop["robust MAP Prior: hazard", "SD"] <- round(exp(rob_Prior_sd), 2)
    stats_mat_prop["robust MAP Prior: hazard", "Median"] <- round(exp(rob_Prior_median), 2)
    stats_mat_prop["robust MAP Prior: hazard", "95% CrI"] <- paste0(
      "[", formatC(exp(RBesT::qmix(robust_map_object, 0.025)),
        digits = 2, format = "f"
      ), " , ",
      formatC(exp(RBesT::qmix(robust_map_object, 0.975)),
        digits = 2, format = "f"
      ), "]"
    )
    stats_mat_prop["robust MAP Prior: hazard", "ESS"] <- c("Not applicable.")
    #

    # sum_stats <- data.frame(
    #   Mean = c(Prior_mean, rob_Prior_mean),
    #   ExpMean = exp(c(Prior_mean, rob_Prior_mean)),
    #   SD = c(Prior_sd, rob_Prior_sd),
    #   ExpSD = exp(c(Prior_sd, rob_Prior_sd)),
    #   Median = c(Prior_median, rob_Prior_median),
    #   ExpMedian = exp(c(Prior_median, rob_Prior_median)),
    #   Conf_int = c(Prior_95ci, rob_Prior_95ci),
    #   ExpConf_int = exp(c(Prior_95ci, rob_Prior_95ci)),
    #   ESS = c(Prior_ess, rob_Prior_ess)
    # )

    sum_stats <- stats_mat_prop

    if (download == FALSE) {
      sum_stats
      sum_stats
      # dplyr::rename("95% CrI" = Conf_int) %>%
      # dplyr::rename("95% exp(CrI)" = ExpConf_int) %>%
      # dplyr::rename("exp(SD)" = ExpSD) %>%
      # dplyr::rename("exp(Median)" = ExpMedian) %>%
    } else {
      sum_stats
    }
  }
}
# Display model summary output in the MAP Prior Tab
#' Title
#'
#' @param map_object
#' @param select_analysis
#' @param param_approx
#' @param ess_method
#' @export
model_summary_display <- function(map_object, select_analysis, param_approx, ess_method) {
  if (is.null(map_object)) {
    return(NULL)
  } else if (select_analysis == "Incidence proportion") {
    desc_stats <- summary(map_object)
    stats_mat_prop <- data.frame(desc_stats$theta.pred)
    rownames(stats_mat_prop) <- "MAP Prior"
    LB <- as.numeric(format(unlist(stats_mat_prop$X2.5)))
    UB <- as.numeric(format(unlist(stats_mat_prop$X97.5)))
    stats_mat_prop$mean <- paste0(formatC(stats_mat_prop$mean * 100, digits = 2, format = "f"), "%")
    stats_mat_prop$sd <- paste0(formatC(stats_mat_prop$sd * 100, digits = 2, format = "f"), "%")
    stats_mat_prop$X50. <- paste0(formatC(stats_mat_prop$X50. * 100, digits = 2, format = "f"), "%")
    LB_perc <- formatC(LB * 100, digits = 2, format = "f")
    UB_perc <- formatC(UB * 100, digits = 2, format = "f")
    stats_mat_prop$conf_int <- paste0("[", LB_perc, "%, ", UB_perc, "%]")
    stats_mat_prop$ESS <- round(RBesT::ess(param_approx, method = ess_method))
    stats_mat_prop %>%
      dplyr::rename(Mean = mean, SD = sd, Median = X50., "95% CrI" = conf_int) %>%
      dplyr::select(Mean, SD, Median, "95% CrI", ESS)
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    # Rbest gives the model alreday out for exp(theta_pred), but we need
    # theta_pred for further analysis
    theta_fit <- as.data.frame(map_object$fit)
    theta_pred <- theta_fit$theta_pred
    theta_resp <- theta_fit$theta_resp_pred

    stats_mat_prop <- array(data = NA, dim = c(2, 5), dimnames = list(
      c("MAP Prior: log(hazard)", "MAP Prior: hazard"),
      c("Mean", "SD", "Median", "95% CrI", "ESS")
    ))
    stats_mat_prop <- as.data.frame(stats_mat_prop)

    stats_mat_prop["MAP Prior: log(hazard)", "Mean"] <- round(mean(theta_pred), 2)
    stats_mat_prop["MAP Prior: log(hazard)", "SD"] <- round(sd(theta_pred), 2)
    stats_mat_prop["MAP Prior: log(hazard)", "Median"] <- round(median(theta_pred), 2)
    stats_mat_prop["MAP Prior: log(hazard)", "95% CrI"] <- paste0(
      "[",
      round(quantile(theta_pred, 0.25), 2), " , ",
      round(quantile(theta_pred, 0.75), 2), "]"
    )
    stats_mat_prop["MAP Prior: log(hazard)", "ESS"] <- round(RBesT::ess(param_approx, method = ess_method, sigma = 1))


    stats_mat_prop["MAP Prior: hazard", "Mean"] <- round(mean(theta_resp), 2)
    stats_mat_prop["MAP Prior: hazard", "SD"] <- round(sd(theta_resp), 2)
    stats_mat_prop["MAP Prior: hazard", "Median"] <- round(median(theta_resp), 2)
    stats_mat_prop["MAP Prior: hazard", "95% CrI"] <- paste0(
      "[",
      round(quantile(theta_resp, 0.25), 2), " , ",
      round(quantile(theta_resp, 0.75), 2), "]"
    )
    stats_mat_prop["MAP Prior: hazard", "ESS"] <- c("Not applicable.")



    stats_mat_prop %>%
      dplyr::select(Mean, SD, Median, "95% CrI", ESS)
  }
}
# Display input data in the historical data tab
#' Title
#'
#' @param data
#' @param select_analysis
#' @param saf_topic
#' @export
input_data_display <- function(data, select_analysis, saf_topic) {
  tab <- data
  rownames(tab) <- NULL
  if (select_analysis == "Incidence proportion") {
    colnames(tab) <- c("STUDY ID", "Number of Patients in Arm", paste0("Number of Patients with ", saf_topic, " in Arm"), "Historical")
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    colnames(tab) <- c("STUDY ID", "Number of Patients in Arm", paste0("Number of Patients with ", saf_topic, " in Arm"), "Total Exposure Time", "Historical")
  }
  tab
}
# Summary statistics for prior, likelihood, and posterior
#' Title
#' @param select_analysis
#' @param robust_map_object
#' @param ess_method
#' @param current_trial_data
#' @param select_analysis
#' @param post_dist
#' @param download
#' @export
summary_stat_all_display <- function(select_analysis, robust_map_object, ess_method, current_trial_data, post_dist, download = FALSE) {
  # assign overall variables
  if (select_analysis == "Incidence proportion") {
    new_n <- current_trial_data[["new_v1"]]
    new_r <- current_trial_data[["new_v2"]]
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    new_n_with_ae <- current_trial_data[["new_v1"]]
    new_tot_exp <- current_trial_data[["new_v2"]]
  }

  if (select_analysis == "Incidence proportion") {
    # Summary statistics for robust MAP prior
    rob_mixture_mat <- data.frame(robust_map_object[, 1:ncol(robust_map_object)])
    rob_vec1 <- vector(length = ncol(rob_mixture_mat)) # Robust MAP prior mean
    rob_vec2 <- vector(length = ncol(rob_mixture_mat)) # Robust MAP prior sd
    rob_vec3 <- vector(length = ncol(rob_mixture_mat)) # Robust MAP prior median
    # Robust MAP prior mean
    for (j in 1:ncol(rob_mixture_mat)) {
      rob_vec1[j] <- rob_mixture_mat[1, j] * (rob_mixture_mat[2, j] / (rob_mixture_mat[2, j] + rob_mixture_mat[3, j]))
    }
    # Robust MAP prior variance
    for (j in 1:ncol(rob_mixture_mat)) {
      rob_vec2[j] <- rob_mixture_mat[1, j] * ((rob_mixture_mat[2, j] * rob_mixture_mat[3, j]) / (((rob_mixture_mat[2, j] + rob_mixture_mat[3, j])^2) * (rob_mixture_mat[2, j] + rob_mixture_mat[3, j] + 1)))
    }
    # Robust MAP prior median
    for (j in 1:ncol(rob_mixture_mat)) {
      rob_vec3[j] <- rob_mixture_mat[1, j] * ((rob_mixture_mat[2, j] - (1 / 3)) / (rob_mixture_mat[2, j] + rob_mixture_mat[3, j] - (2 / 3)))
    }

    # Compute statistics
    rob_Prior_mean <- paste0(formatC(sum(rob_vec1) * 100, digits = 1, format = "f"), "%")
    rob_Prior_sd <- paste0(formatC(sqrt(sum(rob_vec2)) * 100, digits = 1, format = "f"), "%")
    rob_Prior_median <- paste0(formatC(sum(rob_vec3) * 100, digits = 1, format = "f"), "%")
    rob_Prior_95ci <- paste0("[", formatC(RBesT::qmix(robust_map_object, 0.025) * 100, digits = 1, format = "f"), "%, ", formatC(RBesT::qmix(robust_map_object, 0.975) * 100, digits = 1, format = "f"), "%]")
    rob_Prior_ess <- round(RBesT::ess(robust_map_object, method = ess_method))
    # Summary statistics for likelihood
    n_new <- current_trial_data[["new_v1"]]
    r_new <- current_trial_data[["new_v2"]]
    alpha <- r_new + 1
    beta <- n_new - r_new + 1
    lik <- RBesT::mixbeta(inf = c(1, alpha, beta))
    Lik_mean <- paste0(formatC(alpha / (alpha + beta) * 100, digits = 1, format = "f"), "%")
    Lik_sd <- paste0(formatC(sqrt((alpha * beta) / (((alpha + beta)^2) * (alpha + beta + 1))) * 100, digits = 1, format = "f"), "%")
    Lik_median <- paste0(formatC((alpha - (1 / 3)) / (alpha + beta - (2 / 3)) * 100, digits = 1, format = "f"), "%")
    Lik_95ci <- paste0("[", formatC(RBesT::qmix(lik, 0.025) * 100, digits = 1, format = "f"), "%, ", formatC(RBesT::qmix(lik, 0.975) * 100, digits = 1, format = "f"), "%]")
    Lik_ess <- NA # Likelihood does not have an effective sample size
    # Summary statistics for posterior
    post <- post_dist
    post_mixture_mat <- data.frame(post[, 1:ncol(post)])
    post_vec1 <- vector(length = ncol(post_mixture_mat)) # Posterior mean
    post_vec2 <- vector(length = ncol(post_mixture_mat)) # Posterior sd
    post_vec3 <- vector(length = ncol(post_mixture_mat)) # Posterior median
    # Posterior mean
    for (j in 1:ncol(post_mixture_mat)) {
      post_vec1[j] <- post_mixture_mat[1, j] * (post_mixture_mat[2, j] / (post_mixture_mat[2, j] + post_mixture_mat[3, j]))
    }
    # Posterior variance
    for (j in 1:ncol(post_mixture_mat)) {
      post_vec2[j] <- post_mixture_mat[1, j] * ((post_mixture_mat[2, j] * post_mixture_mat[3, j]) / (((post_mixture_mat[2, j] + post_mixture_mat[3, j])^2) * (post_mixture_mat[2, j] + post_mixture_mat[3, j] + 1)))
    }
    # Posterior median
    for (j in 1:ncol(post_mixture_mat)) {
      post_vec3[j] <- post_mixture_mat[1, j] * ((post_mixture_mat[2, j] - (1 / 3)) / (post_mixture_mat[2, j] + post_mixture_mat[3, j] - (2 / 3)))
    }
    # Compute statistics
    Posterior_mean <- paste0(formatC(sum(post_vec1) * 100, digits = 1, format = "f"), "%")
    Posterior_sd <- paste0(formatC(sqrt(sum(post_vec2)) * 100, digits = 1, format = "f"), "%")
    Posterior_median <- paste0(formatC(sum(post_vec3) * 100, digits = 1, format = "f"), "%")
    Posterior_95ci <- paste0("[", formatC(RBesT::qmix(post, 0.025) * 100, digits = 1, format = "f"), "%, ", formatC(RBesT::qmix(post, 0.975) * 100, digits = 1, format = "f"), "%]")
    Posterior_ess <- NA

    # Create table
    sum_stats <- data.frame(
      Mean = c(rob_Prior_mean, Lik_mean, Posterior_mean),
      SD = c(rob_Prior_sd, Lik_sd, Posterior_sd),
      Median = c(rob_Prior_median, Lik_median, Posterior_median),
      conf_int = c(rob_Prior_95ci, Lik_95ci, Posterior_95ci),
      ESS = c(rob_Prior_ess, Lik_ess, Posterior_ess)
    )

    rownames(sum_stats) <- c("Robust MAP Prior", "Likelihood", "Posterior")

    if (download == FALSE) {
      # Likelihood does not have an effective sample size - blank cell

      sum_stats %>%
        dplyr::rename("95% CrI" = conf_int)
    } else {
      sum_stats
    }
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    # Summary statistics for robust MAP prior


    rob_Prior_mean <- round(mean(RBesT::rmix(robust_map_object, 1000)), digits = 2)
    rob_Prior_sd <- round(sd(RBesT::rmix(robust_map_object, 1000)), digits = 2)
    rob_Prior_median <- round(RBesT::qmix(robust_map_object, 0.5), digits = 2)
    rob_Prior_95ci <- paste0(
      "[", formatC(RBesT::qmix(robust_map_object, 0.025),
        digits = 2, format = "f"
      ), " , ",
      formatC(RBesT::qmix(robust_map_object, 0.975),
        digits = 2, format = "f"
      ), "]"
    )
    rob_Prior_ess <- round(RBesT::ess(robust_map_object, method = ess_method, sigma = 1))


    # Summary statistics for likelihood
    new_n_with_ae <- current_trial_data[["new_v1"]]
    new_tot_exp <- current_trial_data[["new_v2"]]

    theta <- new_n_with_ae / new_tot_exp
    like_sample <- rnorm(1000, mean = log(theta), sd = 1 / new_n_with_ae)
    Likelihood <- RBesT::automixfit(like_sample, Nc = seq(3, 3), type = "norm")
    desc_stats <- summary(Likelihood)

    Lik_mean <- formatC(desc_stats["mean"], digits = 2, format = "f")
    Lik_sd <- formatC(desc_stats["sd"], digits = 2, format = "f")
    Lik_median <- formatC(desc_stats["50.0%"], digits = 2, format = "f")
    Lik_95ci <- paste0("[", formatC(RBesT::qmix(Likelihood, 0.025), digits = 2, format = "f"), " , ", formatC(RBesT::qmix(Likelihood, 0.975), digits = 2, format = "f"), "]")
    Lik_ess <- NA

    post_er <-
      posterior_dist(
        select_analysis = "Exposure-adjusted AE rate",
        robust_map_prior = robust_map_object,
        new_v1 = new_n_with_ae,
        new_v2 = new_tot_exp,
        seed = seed
      )

    desc_stats <- summary(post_er)

    Posterior_mean <- formatC(desc_stats["mean"], digits = 2, format = "f")
    Posterior_sd <- formatC(desc_stats["sd"], digits = 2, format = "f")
    Posterior_median <- formatC(desc_stats["50.0%"], digits = 2, format = "f")
    Posterior_95ci <- paste0("[", formatC(RBesT::qmix(post_er, 0.025), digits = 2, format = "f"), " , ", formatC(RBesT::qmix(post_er, 0.975), digits = 2, format = "f"), "]")

    Posterior_ess <- round(RBesT::ess(post_er, method = ess_method, sigma = 1))

    # Create table
    sum_stats <- data.frame(
      Mean = c(rob_Prior_mean, Lik_mean, Posterior_mean),
      SD = c(rob_Prior_sd, Lik_sd, Posterior_sd),
      Median = c(rob_Prior_median, Lik_median, Posterior_median),
      conf_int = c(rob_Prior_95ci, Lik_95ci, Posterior_95ci),
      ESS = c(rob_Prior_ess, Lik_ess, Posterior_ess)
    )

    rownames(sum_stats) <- c("Robust MAP Prior", "Likelihood", "Posterior")

    if (download == FALSE) {
      sum_stats %>%
        dplyr::rename("95% CrI" = conf_int)
    } else {
      sum_stats
    }
  }
}

# Table of preset statistical inference statements
#' Title
#' @param mix
#' @param saf_topic
#' @param select_analysis
#' @export
preset_stat_table <- function(mix, saf_topic, select_analysis) {
  certainty90 <- round(100 * RBesT::qmix(mix, 0.10, lower.tail = TRUE))
  certainty95 <- round(100 * RBesT::qmix(mix, 0.05, lower.tail = TRUE))
  certainty99 <- round(100 * RBesT::qmix(mix, 0.01, lower.tail = TRUE))

  postfix <- ""
  denominator <- 1

  if (select_analysis == "Incidence proportion") {
    postfix <- "%."
    denominator <- 1
    midfix <- " proportion of patients with "
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    postfix <- "."
    denominator <- 100
    midfix <- " incidence rate of patients with "
  }

  inf_mat <- as.data.frame(
    matrix(
      c(
        paste0("We are at least 90% certain that the", midfix, saf_topic, " is greater than ", certainty90 / denominator, postfix),
        paste0("We are at least 95% certain that the", midfix, saf_topic, " is greater than ", certainty95 / denominator, postfix),
        paste0("We are at least 99% certain that the", midfix, saf_topic, " is greater than ", certainty99 / denominator, postfix)
      ),
      nrow = 3, ncol = 1
    )
  )

  colnames(inf_mat) <- ""

  inf_mat
}
