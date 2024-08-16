# template file - please delete after adding real tests
test_that(
  "incidence rates work as expected",
  {
    hist_borrow_choices <- c(
      "Small", "Moderate", "Substantial", "Large",
      "Very Large"
    )
    path <- "scenarios/"
    csv_choices <- list.files(path, pattern = "*.csv")
    csv_choices <- csv_choices[-c(2, 4, 12)]
    testing_list_rates <- readRDS("thresholds/testing_list_rates.rds")

    for (i in c(1:c(length(testing_list_rates)))[-c(3, 5, 13)]) {
      tp_rates <- testing_list_rates[[i]]$parameters
      thresholds_rates <- testing_list_rates[[i]]$treshholds


      data <- NULL
      if (grepl("any", tp_rates$csv)) {
        data <- read.csv(paste0(path, sample(csv_choices, 1)))
      } else {
        data <- read.csv(paste0(path, tp_rates$csv))
      }
      print(tp_rates$saf_topic)
      if (i > 1) {
        data <- data %>%
          dplyr::filter(HIST == 1)
        data <- data_table_prep(
          input_data = data,
          select_analysis = tp_rates$analysis,
          saf_topic = tp_rates$saf_topic,
          select_btrt = tp_rates$group,
          bool_pooled = TRUE # might change later
        )
      }

      adj_tau <- ifelse(tp_rates$heterog == "any",
        tau_adjust(select_analysis = tp_rates$analysis, hist_borrow = sample(hist_borrow_choices, 1)),
        tau_adjust(select_analysis = tp_rates$analysis, hist_borrow = tp_rates$heterog)
      )

      current_trial_data <- list(
        new_v1 = tp_rates$nta_event,
        new_v2 = tp_rates$nta_time
      )

      map_object <- map_prior_func(
        input_data = data,
        select_analysis = tp_rates$analysis,
        tau_dist = tp_rates$tau,
        adj_tau = adj_tau, # scenario specific
        seed = tp_rates$seed,
        testing = FALSE
      )
      param_approx <- parametric_approx(
        select_analysis = tp_rates$analysis,
        map_prior = map_object
      )

      robust_map <- robust_map(
        select_analysis = tp_rates$analysis,
        param_approx = param_approx,
        input_data = data,
        robust_weight = tp_rates$rob_weight,
        robust_mean = tp_rates$rob_mean,
        adj_tau = adj_tau,
        seed = tp_rates$seed,
        testing = FALSE
      )

      post_dist <- posterior_dist(
        select_analysis = tp_rates$analysis,
        input_data = data,
        robust_map_prior = robust_map,
        new_v1 = tp_rates$nta_event,
        new_v2 = tp_rates$nta_time,
        seed = tp_rates$seed
      )

      result_map <- model_summary_display(
        map_object = map_object,
        select_analysis = tp_rates$analysis,
        param_approx = param_approx,
        ess_method = tp_rates$ESS,
        numerical = TRUE
      )

      rerun_counter <- 0
      while (sum(data$N_WITH_AE) < result_map$ESS[1] && rerun_counter < 2) {
        rerun_counter <- rerun_counter + 1

        map_object <- map_prior_func(
          input_data = data,
          select_analysis = tp_rates$analysis,
          tau_dist = tp_rates$tau,
          adj_tau = adj_tau, # scenario specific
          seed = tp_rates$seed + rerun_counter,
          testing = FALSE
        )
        param_approx <- parametric_approx(
          select_analysis = tp_rates$analysis,
          map_prior = map_object
        )

        robust_map <- robust_map(
          select_analysis = tp_rates$analysis,
          param_approx = param_approx,
          input_data = data,
          robust_weight = tp_rates$rob_weight,
          robust_mean = tp_rates$rob_mean,
          adj_tau = adj_tau,
          seed = tp_rates$seed + rerun_counter,
          testing = FALSE
        )

        post_dist <- posterior_dist(
          select_analysis = tp_rates$analysis,
          input_data = data,
          robust_map_prior = robust_map,
          new_v1 = tp_rates$nta_event,
          new_v2 = tp_rates$nta_time,
          seed = tp_rates$seed + rerun_counter
        )

        result_map <- model_summary_display(
          map_object = map_object,
          select_analysis = tp_rates$analysis,
          param_approx = param_approx,
          ess_method = tp_rates$ESS,
          numerical = TRUE
        )
      }

      result_rob <- summary_stats_robust_map_prior_display(
        map_object = map_object,
        select_analysis = tp_rates$analysis,
        param_approx = param_approx,
        ess_method = tp_rates$ESS,
        robust_map_object = robust_map,
        rob_ess_method = tp_rates$ESS, # twice the same?
        numerical = TRUE,
        seed = tp_rates$seed
      )

      result_nta <- summary_stat_all_display(
        select_analysis = tp_rates$analysis,
        robust_map_object = robust_map,
        ess_method = tp_rates$ESS,
        current_trial_data = current_trial_data,
        post_dist = post_dist,
        numerical = TRUE,
        seed = tp_rates$seed
      )

      # transformations ---------------------------------------------------------

      # map

      log_map_string <- "log MAP Prior"
      log_rob_string <- "log Robustified"
      log_rob_string_comp <- "robust MAP Prior: log(hazard)"
      log_likeli_string <- "log Likelihood"
      log_post_string <- "log Posterior"
      exp_map_string <- "exp MAP Prior"
      exp_rob_string <- "exp Robustified"
      exp_rob_string_comp <- "robust MAP Prior: hazard"
      exp_likeli_string <- "exp Likelihood"
      exp_post_string <- "exp Posterior"

      res_string_without_ess <- c("mean", "sd", "median", "crilb", "criub")
      th_string_without_ess_lb <- c("mean_lb", "sd_lb", "median_lb", "cri_lb_lb", "cri_ub_lb")
      th_string_without_ess_ub <- c("mean_ub", "sd_ub", "median_ub", "cri_lb_ub", "cri_ub_ub")

      th_rates_lb <- thresholds_rates[, c("mean_lb", "sd_lb", "median_lb", "cri_lb_lb", "cri_ub_lb", "ess_lb")]
      th_rates_ub <- thresholds_rates[, c("mean_ub", "sd_ub", "median_ub", "cri_lb_ub", "cri_ub_ub", "ess_ub")]
      test_that(paste0("MAP Prior row is higher or equal than lower bound thresholds in ",tp_rates$saf_topic), {
        expect_true(
          all(
            result_map[, res_string_without_ess] >= (th_rates_lb[c(log_map_string, exp_map_string), th_string_without_ess_lb] - 8e-02),
            na.rm = TRUE
          )
        )
      })

      test_that(paste0("MAP Prior ESS is higher or equal than lower bound thresholds in ", tp_rates$saf_topic), {
        expect_true(
          all(
            result_map[, c("ESS")] >= (thresholds_rates[log_map_string, c("ess_min")] * 0.9),
            na.rm = TRUE
          )
        )
      })
      test_that(paste0("NA values are in the same positions in both data frames in ", tp_rates$saf_topic), {
        expect_equal(which(is.na(result_map)), which(is.na(th_rates_lb[c(log_map_string, exp_map_string), ])))
      })
      test_that(paste0("MAP Prior row is lower or equal than upper bound thresholds in ", tp_rates$saf_topic), {
        expect_true(
          all(
            result_map[, res_string_without_ess] <= (th_rates_ub[c(log_map_string, exp_map_string), th_string_without_ess_ub] + 8e-02), # gt/lt tolerance
            na.rm = TRUE
          )
        )
      })
      test_that(paste0("MAP Prior ESS is lower or equal than upper bound thresholds in ", tp_rates$saf_topic), {
        expect_true(
          all(
            result_map[, c("ESS")] <= (th_rates_ub[log_map_string, c("ess_ub")] * 1.1),
            na.rm = TRUE
          )
        )
      })

      test_that(paste0("NA values are in the same positions in both data frames in ", tp_rates$saf_topic), {
        expect_equal(which(is.na(result_map)), which(is.na(th_rates_ub[c(log_map_string, exp_map_string), ])))
      })
      test_that(paste0("MAP Prior Row is identical in both MAP Prior as well as Robust MAP Prior Table in ", tp_rates$saf_topic), {
        expect_true(
          all(
            result_map == result_rob[c(log_map_string, exp_map_string), ],
            na.rm = TRUE
          )
        )
      })
      test_that(paste0("NA values are in the same positions in both data frames in ", tp_rates$saf_topic), {
        expect_equal(which(is.na(result_nta)), which(is.na(th_rates_lb[c(
          log_likeli_string, log_rob_string, log_post_string, exp_rob_string,
          exp_post_string
        ), ])))
      })
      test_that(paste0("NTA Table is higher or equal than lower bound thresholds in ", tp_rates$saf_topic), {
        expect_true(
          all(
            result_nta[, res_string_without_ess] >= (th_rates_lb[c(
              log_likeli_string, log_rob_string, log_post_string, exp_rob_string,
              exp_post_string
            ), th_string_without_ess_lb] - 8e-02),
            na.rm = TRUE
          )
        )
      })
      test_that(paste0("NTA Table ESS is higher or equal than lower bound thresholds in ", tp_rates$saf_topic), {
        expect_true(
          all(
            result_nta[, c("ESS")] >= (thresholds_rates[c(
              log_likeli_string, log_rob_string, log_post_string, exp_rob_string,
              exp_post_string
            ), c("ess_min")] * 0.9),
            na.rm = TRUE
          )
        )
      })
      test_that(paste0("NA values are in the same positions in both data frames in", tp_rates$saf_topic), {
        expect_equal(which(is.na(result_nta)), which(is.na(th_rates_lb[c(
          log_likeli_string, log_rob_string, log_post_string, exp_rob_string,
          exp_post_string
        ), ])))
      })
      # FALSE - replace with difference
      test_that(paste0("NTA Table is lower or equal than upper bound thresholds in ", tp_rates$saf_topic), {
        expect_true(
          all(
            result_nta[, res_string_without_ess] <= (th_rates_ub[c(
              log_likeli_string, log_rob_string, log_post_string, exp_rob_string,
              exp_post_string
            ), th_string_without_ess_ub] + 8e-02),
            na.rm = TRUE
          )
        )
      })

      test_that(paste0("NTA Table ESS is lower or equal than upper bound thresholds in ", tp_rates$saf_topic), {
        expect_true(
          all(
            result_nta[, c("ESS")] <= (th_rates_ub[c(
              log_likeli_string, log_rob_string, log_post_string, exp_rob_string,
              exp_post_string
            ), c("ess_ub")] * 1.1),
            na.rm = TRUE
          )
        )
      })

      test_that(paste0("NA values are in the same positions in both data frames in ", tp_rates$saf_topic), {
        expect_equal(which(is.na(result_nta)), which(is.na(th_rates_ub[c(
          log_likeli_string, log_rob_string, log_post_string, exp_rob_string,
          exp_post_string
        ), ])))
      })
      test_that(paste0("Robust MAP Prior Row is identical in both Robust MAP Prior as well as NTA Table in ", tp_rates$saf_topic), {
        expect_true(
          all(
            result_nta[c(log_rob_string_comp, exp_rob_string_comp), ] == result_rob[c(log_rob_string_comp, exp_rob_string_comp), ],
            na.rm = TRUE
          )
        )
      })
      test_that(paste0("NA values are in the same positions in both data frames in", tp_rates$saf_topic), {
        expect_equal(
          which(is.na(result_nta[c(log_rob_string_comp, exp_rob_string_comp), ])),
          which(is.na(result_rob[c(log_rob_string_comp, exp_rob_string_comp), ]))
        )
      })
    }
  }
)
