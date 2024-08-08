test_that("display.R works as expected", {
  data <- read.csv("scenarios/Scen07.csv")
  pa_prop <- parametric_approx(
    select_analysis = "Incidence proportion",
    map_prior = readRDS("img/prop/map_object_prop3.rds")
  )

  pa_rates <- parametric_approx(
    select_analysis = "Exposure-adjusted AE rate",
    map_prior = readRDS("img/rates/map_object_rate3.rds")
  )

  test_mp_func_display_prop <- map_prior_function_display(
    param_approx = pa_prop,
    select_analysis = "Incidence proportion"
  )

  test_mp_func_display_rates <- map_prior_function_display(
    param_approx = pa_rates,
    select_analysis = "Exposure-adjusted AE rate"
  )

  expect_no_error(test_mp_func_display_prop)
  expect_no_error(test_mp_func_display_rates)

  expect_type(test_mp_func_display_prop, "list")
  expect_type(test_mp_func_display_rates, "list")

  adj_tau_prop <- tau_adjust(select_analysis = "Incidence proportion", hist_borrow = "Very Large")
  adj_tau_rates <- tau_adjust(select_analysis = "Exposure-adjusted AE rate", hist_borrow = "Very Large")

  robust_map_prop <- robust_map(
    select_analysis = "Incidence proportion",
    param_approx = pa_prop,
    input_data = data,
    robust_weight = 0.2,
    robust_mean = 0.5,
    adj_tau = adj_tau_prop,
    seed = 1701416989,
    testing = FALSE
  )

  robust_map_rates <- robust_map(
    select_analysis = "Exposure-adjusted AE rate",
    param_approx = pa_rates,
    input_data = data,
    robust_weight = 0.2,
    robust_mean = 0.19,
    adj_tau = adj_tau_rates,
    seed = 1701416989,
    testing = FALSE
  )

  rmp_mix_dens_display_prop <- robust_map_prior_mix_dens_display(
    robust_map_object = robust_map_prop,
    select_analysis = "Incidence proportion"
  )
  rmp_mix_dens_display_rates <- robust_map_prior_mix_dens_display(
    robust_map_object = robust_map_rates,
    select_analysis = "Exposure-adjusted AE rate"
  )

  expect_no_error(rmp_mix_dens_display_prop)
  expect_no_error(rmp_mix_dens_display_rates)

  expect_type(rmp_mix_dens_display_prop, "list")
  expect_type(rmp_mix_dens_display_rates, "list")

  current_trial_data_prop <- list(
    new_v1 = 200,
    new_v2 = 35
  )

  current_trial_data_rates <- list(
    new_v1 = 35,
    new_v2 = 200
  )

  post_dist_prop <- posterior_dist(
    select_analysis = "Incidence proportion",
    input_data = data,
    robust_map_prior = robust_map_prop,
    new_v1 = 200,
    new_v2 = 35,
    seed = 1701416989
  )

  post_dist_rates <- posterior_dist(
    select_analysis = "Exposure-adjusted AE rate",
    input_data = data,
    robust_map_prior = robust_map_rates,
    new_v1 = 35,
    new_v2 = 200,
    seed = 1701416989
  )


  # MP ----------------------------------------------------------------------

  dist <- "MAP Prior"

  mix_prop <- mix_distribution_all(
    select_analysis = "Incidence proportion",
    current_trial_data = current_trial_data_prop,
    select_dist = dist,
    param_approx = pa_prop,
    robust_map_object = robust_map_prop,
    post_dist = post_dist_prop,
    seed = 1701416989
  )

  mix_rates <- mix_distribution_all(
    select_analysis = "Exposure-adjusted AE rate",
    current_trial_data = current_trial_data_rates,
    select_dist = dist,
    param_approx = pa_rates,
    robust_map_object = robust_map_rates,
    post_dist = post_dist_rates,
    seed = 1701416989
  )

  auc_prop <- area_under_the_curve(
    ae_prop = c(30, 100) / 100,
    mix = mix_prop,
    saf_topic = "Scen07"
  )

  auc_rates <- area_under_the_curve(
    ae_prop = c(30, 100) / 100,
    mix = mix_rates,
    saf_topic = "Scen07"
  )

  expect_no_error(auc_prop)
  expect_no_error(auc_rates)

  expect_type(auc_prop, "character")
  expect_type(auc_rates, "character")


  # RMP ---------------------------------------------------------------------

  dist <- "Robust MAP Prior"

  mix_prop <- mix_distribution_all(
    select_analysis = "Incidence proportion",
    current_trial_data = current_trial_data_prop,
    select_dist = dist,
    param_approx = pa_prop,
    robust_map_object = robust_map_prop,
    post_dist = post_dist_prop,
    seed = 1701416989
  )

  mix_rates <- mix_distribution_all(
    select_analysis = "Exposure-adjusted AE rate",
    current_trial_data = current_trial_data_rates,
    select_dist = dist,
    param_approx = pa_rates,
    robust_map_object = robust_map_rates,
    post_dist = post_dist_rates,
    seed = 1701416989
  )

  auc_prop <- area_under_the_curve(
    ae_prop = c(30, 100) / 100,
    mix = mix_prop,
    saf_topic = "Scen07"
  )

  auc_rates <- area_under_the_curve(
    ae_prop = c(30, 100) / 100,
    mix = mix_rates,
    saf_topic = "Scen07"
  )

  expect_no_error(auc_prop)
  expect_no_error(auc_rates)

  expect_type(auc_prop, "character")
  expect_type(auc_rates, "character")


  # Likelihood --------------------------------------------------------------

  dist <- "Likelihood"

  mix_prop <- mix_distribution_all(
    select_analysis = "Incidence proportion",
    current_trial_data = current_trial_data_prop,
    select_dist = dist,
    param_approx = pa_prop,
    robust_map_object = robust_map_prop,
    post_dist = post_dist_prop,
    seed = 1701416989
  )

  mix_rates <- mix_distribution_all(
    select_analysis = "Exposure-adjusted AE rate",
    current_trial_data = current_trial_data_rates,
    select_dist = dist,
    param_approx = pa_rates,
    robust_map_object = robust_map_rates,
    post_dist = post_dist_rates,
    seed = 1701416989
  )

  auc_prop <- area_under_the_curve(
    ae_prop = c(30, 100) / 100,
    mix = mix_prop,
    saf_topic = "Scen07"
  )

  auc_rates <- area_under_the_curve(
    ae_prop = c(30, 100) / 100,
    mix = mix_rates,
    saf_topic = "Scen07"
  )

  expect_no_error(auc_prop)
  expect_no_error(auc_rates)

  expect_type(auc_prop, "character")
  expect_type(auc_rates, "character")


  # Posterior ---------------------------------------------------------------

  dist <- "Posterior"

  mix_prop <- mix_distribution_all(
    select_analysis = "Incidence proportion",
    current_trial_data = current_trial_data_prop,
    select_dist = dist,
    param_approx = pa_prop,
    robust_map_object = robust_map_prop,
    post_dist = post_dist_prop,
    seed = 1701416989
  )

  mix_rates <- mix_distribution_all(
    select_analysis = "Exposure-adjusted AE rate",
    current_trial_data = current_trial_data_rates,
    select_dist = dist,
    param_approx = pa_rates,
    robust_map_object = robust_map_rates,
    post_dist = post_dist_rates,
    seed = 1701416989
  )

  auc_prop <- area_under_the_curve(
    ae_prop = c(30, 100) / 100,
    mix = mix_prop,
    saf_topic = "Scen07"
  )

  auc_rates <- area_under_the_curve(
    ae_prop = c(30, 100) / 100,
    mix = mix_rates,
    saf_topic = "Scen07"
  )

  expect_no_error(auc_prop)
  expect_no_error(auc_rates)

  expect_type(auc_prop, "character")
  expect_type(auc_rates, "character")
})
