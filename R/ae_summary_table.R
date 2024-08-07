#' @title Posterior Summary Statistics for Safety Topics of Interest
#'
#' @description Creates a kable table of posterior summary statistics for each
#' safety topic of interest. Creates the arrays for the pdf-download file.
#'
#' @param input_data dataframe from data_table_prep() including the new trial
#' @param saf_topic Selected safety topic to analyze/the adverse event of interest
#' @param seed a seed
#' @param cb_list_ctrl expects List with control Arm indicators
#' @param cb_list_trt expects List with treatment Arm indicators
#'
#' @return kable table of posterior summary statistics for each event type for rates and proportions
#' @export
ae_summary_table <-
  function(input_data,
           cb_list_ctrl,
           cb_list_trt,
           saf_topic,
           seed = NA) {
    if (length(cb_list_ctrl) != length(cb_list_trt)) {
      stop("Same amount of compared groups necessary.")
    }

    # Setup -------------------------------------------------------------------

    # Number of comparisons
    n_group <- length(cb_list_ctrl)

    # setting up overall variables
    tau <- tau_adjust(
      select_analysis = "Incidence proportion",
      hist_borrow = "Large"
    )

    # Weight for robustification
    robust_weight <- 0.2

    if (is.null(seed)) {
      seed <- as.numeric(Sys.time())
    }

    # Credible Interval boundaries
    lb <- 0.025
    ub <- 0.975

    # For later sampling of risk ratio and risk difference
    sample_size <- 10000

    # Warning texts for missing data
    warn_txt <- data.frame(
      Issue = character(), Topic = character(),
      Analysis = character(), Message = character(),
      Group = character(), stringsAsFactors = FALSE
    )

    # Initiating the arrays for the download tables ---------------------------

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

    # Filtering for respective group ------------------------------------------

    ID_ctr_list <- list(grp1 = NULL, grp2 = NULL, grp3 = NULL, grp4 = NULL)
    ID_trt_list <- list(grp1 = NULL, grp2 = NULL, grp3 = NULL, grp4 = NULL)
    for (group in 1:n_group) {
      # for the log
      print(paste0("Comparison ", group, " started at ", Sys.time()))

      trt_arm <- input_data %>% dplyr::filter(ARM == cb_list_trt[[group]])
      ID_trt_list[[group]] <- unique(trt_arm$STUDYID)
      ctr_arm <- input_data %>% dplyr::filter(ARM == cb_list_ctrl[[group]])
      ID_ctr_list[[group]] <- unique(ctr_arm$STUDYID)

      # Proportions Prep -----------------------------------------------------------

      for (topic in saf_topic) {
        # for the log
        print(paste0(
          "Variable analysis ", topic,
          " in comparison ", group, " started at ", Sys.time()
        ))

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

        hist_trt_trials <- trt_trials %>% dplyr::filter(HIST == 1)
        hist_ctr_trials <- ctr_trials %>% dplyr::filter(HIST == 1)


        trt_current_trial <- trt_trials %>% dplyr::filter(HIST == 0)
        ctr_current_trial <- ctr_trials %>% dplyr::filter(HIST == 0)

        data_check_bin <- data_available(
          hist_trt = hist_trt_trials,
          hist_ctr = hist_ctr_trials,
          trt_current = trt_current_trial,
          ctr_current = ctr_current_trial
        )
        dc_bin <- sum(data_check_bin)

        warn_txt <- add_row(
          df = warn_txt, analysis = "Incidence proportion",
          data_check = data_check_bin, group = group, topic = topic
        )

        # If there is no prop information, there won't be rate info, jump to the next topic
        if (dc_bin == 0) {
          next
        }

        # Naive estimation --------------------------------------------------------

        array_inci <- inci_naiv(
          data = trt_trials, array_inci = array_inci, arm = "Arm A",
          topic = topic, group = group
        )

        array_inci <- inci_naiv(
          data = ctr_trials, array_inci = array_inci, arm = "Arm B",
          topic = topic, group = group
        )

        ### Binomial
        # weakly informative prior Be(1,1)
        weak_inf_bin <- RBesT::mixbeta(c(1, 1, 1))
        robust_mean_bin <- 0.5

        # Arm A Binomial ----------------------------------------------------------

        # for the log
        print(paste0(
          "ARM A: Binomial analysis for ", topic,
          " in comparison ", group, " started at ", Sys.time()
        ))

        if (data_check_bin["hist_trt"] == TRUE) {
          trt_map_prior <- map_prior_func(
            input_data      = hist_trt_trials,
            select_analysis = "Incidence proportion",
            tau_dist        = "HalfNormal",
            adj_tau         = tau,
            seed            = seed
          )
          trt_approx <- parametric_approx("Incidence proportion", trt_map_prior)
        } else {
          # If there is no historical information
          trt_approx <- weak_inf_bin
        }

        if (data_check_bin["hist_trt"] == TRUE) {
          trt_rob_prior <- RBesT::robustify(trt_approx, weight = robust_weight, mean = robust_mean_bin)
        } else {
          trt_rob_prior <- weak_inf_bin
        }

        if (data_check_bin["trt_current"] == TRUE) {
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
        } else {
          trt_post <- trt_rob_prior
        }

        array_inci <- array_rmix(
          rmix_obj = trt_post, array = array_inci, arm = "Arm A",
          topic = topic, group = group, lb = lb, ub = ub
        )

        # Arm B Proportions -------------------------------------------------------

        # for the log
        print(paste0(
          "ARM B: Binomial analysis for ", topic,
          " in comparison ", group, " started at ", Sys.time()
        ))

        if (data_check_bin["hist_ctr"] == TRUE) {
          ctr_map_prior <- map_prior_func(
            input_data = hist_ctr_trials,
            select_analysis = "Incidence proportion",
            tau_dist = "HalfNormal",
            adj_tau = tau,
            seed = seed
          )
          ctr_approx <- parametric_approx("Incidence proportion", ctr_map_prior)
        } else {
          ctr_approx <- weak_inf_bin
        }


        if (data_check_bin["hist_trt"] == TRUE) {
          ctr_rob_prior <- RBesT::robustify(ctr_approx, weight = robust_weight, mean = robust_mean_bin)
        } else {
          ctr_rob_prior <- weak_inf_bin
        }


        if (data_check_bin["ctr_current"] == TRUE) {
          ctr_new_n <- sum(ctr_current_trial$N)
          ctr_new_r <- sum(ctr_current_trial$N_WITH_AE)

          ctr_post <-
            posterior_dist(
              select_analysis = "Incidence proportion",
              robust_map_prior = ctr_rob_prior,
              new_v1 = ctr_new_n,
              new_v2 = ctr_new_r,
              seed = seed
            )
        } else {
          ctr_post <- ctr_rob_prior
        }

        array_inci <- array_rmix(
          rmix_obj = ctr_post, array = array_inci, arm = "Arm B",
          topic = topic, group = group, lb = lb, ub = ub
        )


        # proportion posterior samples ---------------------------------------------

        set.seed(seed)
        trt_pr_post_samp <- RBesT::rmix(trt_post, sample_size)
        ctr_pr_post_samp <- RBesT::rmix(ctr_post, sample_size)

        pr_diffe_samp <- trt_pr_post_samp - ctr_pr_post_samp
        pr_ratio_samp <- trt_pr_post_samp / ctr_pr_post_samp

        array_inci_comp <- array_comp(
          ana = "inci",
          array_comp = array_inci_comp, comp = "Risk Diff", crilb = lb, criub = ub,
          pr_sample = pr_diffe_samp, array_ana = array_inci
        )

        array_inci_comp <- array_comp(
          ana = "inci",
          array_comp = array_inci_comp, comp = "Risk Ratio", crilb = lb, criub = ub,
          pr_sample = pr_ratio_samp, array_ana = array_inci
        )

        # Exposure adjusted event rates prep -------------------------------------------

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
        hist_ctr_trials <- ctr_trials %>% dplyr::filter(HIST == 1)

        trt_current_trial <- trt_trials %>% dplyr::filter(HIST == 0)
        ctr_current_trial <- ctr_trials %>% dplyr::filter(HIST == 0)

        data_check_rate <- data_available(
          hist_trt = hist_trt_trials,
          hist_ctr = hist_ctr_trials,
          trt_current = trt_current_trial,
          ctr_current = ctr_current_trial
        )
        dc_rate <- sum(data_check_rate)

        warn_txt <- add_row(
          df = warn_txt, analysis = "Exposure-adjusted AE rate",
          data_check = data_check_rate, group = group, topic = topic
        )

        # If there is no rate information, there won't be rate info, jump to the next topic
        if (dc_rate == 0) {
          array_er[topic, , , group] <- "No info available"
          next
        }


        # Naive estimation --------------------------------------------------------

        array_er[topic, "naive", "Arm A", group] <- sum(trt_trials$N_WITH_AE) / sum(trt_trials$TOT_EXP)
        array_er[topic, "naive", "Arm B", group] <- sum(ctr_trials$N_WITH_AE) / sum(ctr_trials$TOT_EXP)

        # For a reasonable posterior historical information should be available,
        # as well as new data.
        trt_post_check <- nrow(trt_current_trial) > 0 && nrow(hist_trt_trials) > 0


        # Rates Arm A -------------------------------------------------------------

        # for the log
        print(paste0(
          "ARM A: Rate analysis for ", topic,
          " in comparison ", group, " started at ", Sys.time()
        ))

        if (data_check_rate["hist_trt"] == TRUE) {
          er_trt_prior <- map_prior_func(
            input_data = hist_trt_trials,
            select_analysis = "Exposure-adjusted AE rate",
            tau_dist = "HalfNormal",
            adj_tau = tau,
            seed = seed
          )
          er_trt_prior_fit <- parametric_approx("Exposure-adjusted AE rate", er_trt_prior)
        } else {
          # If there is no historical information
          er_trt_mean <- log(array_er[topic, "naive", "Arm A", group])
          er_trt_prior_fit <- RBesT::mixnorm(
            inf1 = c(1 / 3, er_trt_mean, 1),
            inf2 = c(1 / 3, er_trt_mean, 1),
            inf3 = c(1 / 3, er_trt_mean, 1), sigma = 1
          )
        }

        if (data_check_rate["hist_trt"] == TRUE) {
          robust_mean_er <- summary(er_trt_prior_fit)["mean"]
          trt_rob_prior <- RBesT::robustify(er_trt_prior_fit, weight = robust_weight, mean = robust_mean_er, sigma = 1)
        } else {
          trt_rob_prior <- RBesT::robustify(er_trt_prior_fit, weight = er_trt_mean, mean = robust_mean_er, sigma = 1)
        }


        if (data_check_rate["trt_current"] == TRUE) {
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
        } else {
          trt_post <- trt_rob_prior
        }

        array_er <- array_rmix(
          rmix_obj = trt_post, array = array_er, arm = "Arm A",
          topic = topic, group = group, lb = lb, ub = ub
        )


        # Rates Arm B -------------------------------------------------------------

        # for the log
        print(paste0(
          "ARM B: Rate analysis for ", topic,
          " in comparison ", group, " started at ", Sys.time()
        ))

        if (data_check_rate["hist_ctr"] == TRUE) {
          er_ctr_prior <- map_prior_func(
            input_data = hist_ctr_trials,
            select_analysis = "Exposure-adjusted AE rate",
            tau_dist = "HalfNormal",
            adj_tau = tau,
            seed = seed
          )
          er_ctr_prior_fit <- parametric_approx("Exposure-adjusted AE rate", er_ctr_prior)
        } else {
          # If there is no historical information
          er_ctr_mean <- log(array_er[topic, "naive", "Arm B", group])
          er_ctr_prior_fit <- RBesT::mixnorm(
            inf1 = c(1 / 3, er_ctr_mean, 1),
            inf2 = c(1 / 3, er_ctr_mean, 1),
            inf3 = c(1 / 3, er_ctr_mean, 1), sigma = 1
          )
        }

        if (data_check_rate["hist_ctr"] == TRUE) {
          robust_mean_er <- summary(er_ctr_prior_fit)["mean"]
          ctr_rob_prior <- RBesT::robustify(er_ctr_prior_fit, weight = robust_weight, mean = robust_mean_er, sigma = 1)
        } else {
          ctr_rob_prior <- RBesT::robustify(er_ctr_prior_fit, weight = er_ctr_mean, mean = robust_mean_er, sigma = 1)
        }

        if (data_check_rate["ctr_current"] == TRUE) {
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
        } else {
          ctr_post <- ctr_rob_prior
        }

        array_er <- array_rmix(
          rmix_obj = ctr_post, array = array_er, arm = "Arm B",
          topic = topic, group = group, lb = lb, ub = ub
        )


        # rates posterior samples ---------------------------------------------------------

        set.seed(seed)
        trt_pr_post_samp <- RBesT::rmix(trt_post, sample_size)
        ctr_pr_post_samp <- RBesT::rmix(ctr_post, sample_size)

        pr_diffe_samp <- trt_pr_post_samp - ctr_pr_post_samp
        pr_ratio_samp <- trt_pr_post_samp / ctr_pr_post_samp

        array_er_comp <- array_comp(
          ana = "rate",
          array_comp = array_er_comp, comp = "Risk Diff", crilb = lb, criub = ub,
          pr_sample = pr_diffe_samp, array_ana = array_er
        )
        array_er_comp <- array_comp(
          ana = "rate",
          array_comp = array_er_comp, comp = "Risk Ratio", crilb = lb, criub = ub,
          pr_sample = pr_ratio_samp, array_ana = array_er
        )
      }
    }

    return(list(
      InciProp = array_inci, InciPropComp = array_inci_comp,
      EventRates = array_er, EvenRatesComp = array_er_comp,
      ID_ctr_list = ID_ctr_list, ID_trt_list = ID_trt_list,
      cb_list_ctrl = cb_list_ctrl, cb_list_trt = cb_list_trt,
      warn_txt
    ))
  }
