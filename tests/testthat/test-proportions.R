# template file - please delete after adding real tests
test_that(
  "proportions work as expected",
  {
    hist_borrow_choices <- c(
      "Small", "Moderate", "Substantial", "Large",
      "Very Large"
    )
    path <- "scenarios/"
    csv_choices <- list.files(path, pattern = "*.csv")
    csv_choices <- csv_choices[-c(4, 12, 13)]
    robust_mean <- 0.5
    testing_list_prop <- readRDS("thresholds/testing_list_props.rds")

    for (i in c(1:length(testing_list_prop))[-c(5, 13, 14)]) {
      tp_prop <- testing_list_prop[[i]]$parameters
      thresholds_prop <- testing_list_prop[[i]]$treshholds

      data <- NULL
      if (grepl("any", tp_prop$csv)) {
        data <- read.csv(paste0(path, sample(csv_choices, 1)))
      } else {
        data <- read.csv(paste0(path, tp_prop$csv))
      }
      print(tp_prop$saf_topic)
      if (i > 1) {
        data <- data %>%
          dplyr::filter(HIST == 1)
        data <- data_table_prep(
          input_data = data,
          select_analysis = tp_prop$analysis,
          saf_topic = tp_prop$saf_topic,
          select_btrt = tp_prop$group,
          bool_pooled = TRUE # might change later
        )
      }

      # sum(data$n_pat) >= ess
      adj_tau <- ifelse(tp_prop$heterog == "any",
        tau_adjust(select_analysis = tp_prop$analysis, hist_borrow = sample(hist_borrow_choices, 1)),
        tau_adjust(select_analysis = tp_prop$analysis, hist_borrow = tp_prop$heterog)
      )

      current_trial_data <- list(
        new_v1 = tp_prop$nta_npat,
        new_v2 = tp_prop$nta_event
      )

      map_object <- map_prior_func(
        input_data = data,
        select_analysis = tp_prop$analysis,
        tau_dist = tp_prop$tau,
        adj_tau = adj_tau, # scenario specific
        seed = tp_prop$seed,
        testing = FALSE
      )
      param_approx <- parametric_approx(
        select_analysis = tp_prop$analysis,
        map_prior = map_object
      )

      robust_map <- robust_map(
        select_analysis = tp_prop$analysis,
        param_approx = param_approx,
        input_data = data,
        robust_weight = tp_prop$rob_weight,
        robust_mean = robust_mean, # needs to be adjusted later
        adj_tau = adj_tau,
        seed = tp_prop$seed,
        testing = FALSE
      )

      post_dist <- posterior_dist(
        select_analysis = tp_prop$analysis,
        input_data = data,
        robust_map_prior = robust_map,
        new_v1 = tp_prop$nta_npat,
        new_v2 = tp_prop$nta_event,
        seed = tp_prop$seed
      )

      result_map <- model_summary_display(
        map_object = map_object,
        select_analysis = tp_prop$analysis,
        param_approx = param_approx,
        ess_method = tp_prop$ESS,
        numerical = TRUE
      )

      rerun_counter <- 0
      while (sum(data$N) < result_map$ESS && rerun_counter < 2) {
        rerun_counter <- rerun_counter + 1

        map_object <- map_prior_func(
          input_data = data,
          select_analysis = tp_prop$analysis,
          tau_dist = tp_prop$tau,
          adj_tau = adj_tau, # scenario specific
          seed = tp_prop$seed + Sys.time(),
          testing = FALSE
        )

        param_approx <- parametric_approx(
          select_analysis = tp_prop$analysis,
          map_prior = map_object
        )

        robust_map <- robust_map(
          select_analysis = tp_prop$analysis,
          param_approx = param_approx,
          input_data = data,
          robust_weight = tp_prop$rob_weight,
          robust_mean = robust_mean, # needs to be adjusted later
          adj_tau = adj_tau,
          seed = tp_prop$seed + Sys.time(),
          testing = FALSE
        )
        post_dist <- posterior_dist(
          select_analysis = tp_prop$analysis,
          input_data = data,
          robust_map_prior = robust_map,
          new_v1 = tp_prop$nta_npat,
          new_v2 = tp_prop$nta_event,
          seed = tp_prop$seed + Sys.time()
        )

        result_map <- model_summary_display(
          map_object = map_object,
          select_analysis = tp_prop$analysis,
          param_approx = param_approx,
          ess_method = tp_prop$ESS,
          numerical = TRUE
        )
      }

      result_rob <- summary_stats_robust_map_prior_display(
        map_object = map_object,
        select_analysis = tp_prop$analysis,
        param_approx = param_approx,
        ess_method = tp_prop$ESS,
        robust_map_object = robust_map,
        rob_ess_method = tp_prop$ESS, # twice the same?
        numerical = TRUE,
        seed = tp_prop$seed
      )

      result_nta <- summary_stat_all_display(
        select_analysis = tp_prop$analysis,
        robust_map_object = robust_map,
        ess_method = tp_prop$ESS,
        current_trial_data = current_trial_data,
        post_dist = post_dist,
        numerical = TRUE,
        seed = tp_prop$seed
      )



      #
      #       # transformations ---------------------------------------------------------
      #
      #       # map
      #
      map_string <- "MAP Prior"
      rob_string <- "Robustified"
      rob_string_rob <- "Robustified MAP Prior"
      rob_string_nta <- "Robust MAP Prior"
      likeli_string <- "Likelihood"
      post_string <- "Posterior"

      res_string_without_ess <- c("mean", "sd", "median", "crilb", "criub")
      th_string_without_ess_lb <- c("mean_lb", "sd_lb", "median_lb", "cri_lb_lb", "cri_ub_lb")
      th_string_without_ess_ub <- c("mean_ub", "sd_ub", "median_ub", "cri_lb_ub", "cri_ub_ub")


      th_prop_lb <- thresholds_prop[, c("mean_lb", "sd_lb", "median_lb", "cri_lb_lb", "cri_ub_lb", "ess_lb")]
      th_prop_ub <- thresholds_prop[, c("mean_ub", "sd_ub", "median_ub", "cri_lb_ub", "cri_ub_ub", "ess_ub")]
      test_that(paste0("MAP Prior row is higher or equal than lower bound thresholds in ", tp_prop$saf_topic), {
        expect_true(
          all(
            result_map[, res_string_without_ess] >= (th_prop_lb[map_string, th_string_without_ess_lb] - 3e-02)
          )
        )
      })

      test_that(paste0("MAP Prior ESS is higher or equal than lower bound thresholds in ", tp_prop$saf_topic), {
        expect_true(
          all(
            result_map[, c("ESS")] >= (thresholds_prop[map_string, c("ess_min")] * 0.9)
          )
        )
      })
      test_that(paste0("MAP Prior row is lower or equal than upper bound thresholds in ", tp_prop$saf_topic), {
        expect_true(
          all(
            result_map[, res_string_without_ess] <= (th_prop_ub[map_string, th_string_without_ess_ub] + 3e-02)
          )
        )
      })
      test_that(paste0("MAP Prior ESS is lower or equal than upper bound thresholds in ", tp_prop$saf_topic), {
        expect_true(
          all(
            result_map[, c("ESS")] <= (th_prop_ub[map_string, c("ess_ub")] * 1.1)
          )
        )
      })
      test_that(paste0("RMAP Prior Row is identical in both MAP Prior as well as Robust MAP Prior Table in ", tp_prop$saf_topic), {
        expect_true(
          all(
            result_map == result_rob[map_string, ]
          )
        )
      })

      # test_that("RMAP Prior Row is identical in both MAP Prior as well as Robust MAP Prior Table", {
      #   expect_equal(
      #     result_map, result_rob[map_string, ],
      #     tolerance = 1e-03 # tolerance_var
      #   )
      # })

      test_that(paste0("NTA Table is higher or equal than lower bound thresholds in ", tp_prop$saf_topic), {
        expect_true(
          all(
            result_nta[, res_string_without_ess] >= (th_prop_lb[c(rob_string, likeli_string, post_string), th_string_without_ess_lb] - 3e-02),
            na.rm = TRUE
          )
        )
      })
      test_that(paste0("NTA Table ESS is higher or equal than lower bound thresholds in ", tp_prop$saf_topic), {
        expect_true(
          all(
            result_nta[, c("ESS")] >= (thresholds_prop[c(rob_string, likeli_string, post_string), c("ess_min")] * 0.9),
            na.rm = TRUE
          )
        )
      })
      test_that(paste0("NA values are in the same positions in both data frames in ", tp_prop$saf_topic), {
        expect_equal(which(is.na(result_nta)), which(is.na(th_prop_lb[c(rob_string, likeli_string, post_string), ])))
      })
      # FALSE - replace with difference
      test_that(paste0("NTA Table is lower or equal than upper bound thresholds in ", tp_prop$saf_topic), {
        expect_true(
          all(
            result_nta[, res_string_without_ess] <= (th_prop_ub[c(rob_string, likeli_string, post_string), th_string_without_ess_ub] + 3e-02),
            na.rm = TRUE
          )
        )
      })
      test_that(paste0("NTA Table ESS is lower or equal than upper bound thresholds in ", tp_prop$saf_topic), {
        expect_true(
          all(
            result_nta[, c("ESS")] <= (th_prop_ub[c(rob_string, likeli_string, post_string), c("ess_ub")] * 1.1),
            na.rm = TRUE
          )
        )
      })

      test_that(paste0("NA values are in the same positions in both data frames in ", tp_prop$saf_topic), {
        expect_equal(which(is.na(result_nta)), which(is.na(th_prop_ub[c(rob_string, likeli_string, post_string), ])))
      })
      test_that(paste0("Robust MAP Prior Row is identical in both Robust MAP Prior as well as NTA Table in ", tp_prop$saf_topic), {
        expect_true(
          all(
            result_nta[rob_string_nta, ] == result_rob[rob_string_rob, ]
          )
        )
      })
    }
  }
)
