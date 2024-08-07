library(dplyr)
library(ggplot2)
# proportions -------------------------------------------------------------



hist_borrow_choices <- c(
  "Small", "Moderate", "Substantial", "Large",
  "Very Large"
)
path <- "scenarios/"
csv_choices <- list.files(path, pattern = "*.csv")
robust_mean <- 0.5
testing_list_prop <- readRDS("thresholds/testing_list_props.rds")

tp_prop <- testing_list_prop[[4]]$parameters
thresholds_prop <- testing_list_prop[[4]]$treshholds

data <- read.csv(paste0(path, tp_prop$csv))

data <- data_table_prep(
  input_data = data,
  select_analysis = tp_prop$analysis,
  saf_topic = tp_prop$saf_topic,
  select_btrt = tp_prop$group,
  bool_pooled = TRUE # might change later
)

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

# rerun_counter <- 0
# while (sum(data$N) < result_map$ESS && rerun_counter < 2) {
#   rerun_counter <- rerun_counter + 1
#
#   map_object <- map_prior_func(
#     input_data = data,
#     select_analysis = tp_prop$analysis,
#     tau_dist = tp_prop$tau,
#     adj_tau = adj_tau, # scenario specific
#     seed = tp_prop$seed + Sys.time(),
#     testing = FALSE
#   )
#
#   param_approx <- parametric_approx(
#     select_analysis = tp_prop$analysis,
#     map_prior = map_object
#   )
#
#   robust_map <- robust_map(
#     select_analysis = tp_prop$analysis,
#     param_approx = param_approx,
#     input_data = data,
#     robust_weight = tp_prop$rob_weight,
#     robust_mean = robust_mean, # needs to be adjusted later
#     adj_tau = adj_tau,
#     seed = tp_prop$seed + Sys.time(),
#     testing = FALSE
#   )
#   post_dist <- posterior_dist(
#     select_analysis = tp_prop$analysis,
#     input_data = data,
#     robust_map_prior = robust_map,
#     new_v1 = tp_prop$nta_npat,
#     new_v2 = tp_prop$nta_event,
#     seed = tp_prop$seed + Sys.time()
#   )
#
#   result_map <- model_summary_display(
#     map_object = map_object,
#     select_analysis = tp_prop$analysis,
#     param_approx = param_approx,
#     ess_method = tp_prop$ESS,
#     numerical = TRUE
#   )
# }


rob_comp <- robust_compare(
  select_analysis = tp_prop$analysis,
  robust_map_prior = robust_map,
  param_approx = param_approx
)

new_trial_compare <- new_trial_compare(
  select_analysis = tp_prop$analysis,
  robust_map_prior = robust_map,
  new_v1 = current_trial_data$new_v1,
  new_v2 = current_trial_data$new_v2,
  post_dist = post_dist
)

forest_plot <- forest_plot_display(
  map_object = map_object,
  select_analysis = tp_prop$analysis,
  saf_topic = tp_prop$saf_topic,
  select_btrt = tp_prop$group
)

rob_map_plot <- robust_map_prior_plot(
  rob_comp = rob_comp,
  saf_topic = tp_prop$saf_topic,
  select_btrt = tp_prop$group,
  select_analysis = tp_prop$analysis
)

param_mix_density <- param_mix_density_display(
  param_approx = param_approx,
  select_analysis = tp_prop$analysis,
  saf_topic = tp_prop$saf_topic,
  select_btrt = tp_prop$group
)

nta_plot <- nta_data_conflict_assassment_plot(
  new_trial_analysis = new_trial_compare,
  saf_topic = tp_prop$saf_topic,
  select_btrt = tp_prop$group,
  select_analysis = tp_prop$analysis
)

choices <- c("Likelihood", "MAP Prior", "Robust MAP Prior", "Posterior")
stat_inf_dist <- list()
decision_making_plot <- list()

for (i in 1:length(choices)) {
  stat_inf_dist[[i]] <- sampling_all_plot(
    select_analysis = tp_prop$analysis,
    select_dist = choices[i],
    param_approx = param_approx,
    new_trial_analysis = new_trial_compare
  )

  decision_making_plot[[i]] <- decision_making_density_plot(
    select_analysis = tp_prop$analysis,
    stat_inf_dist = stat_inf_dist[[i]],
    ae_prop = c(30, 100) / 100,
    saf_topic = tp_prop$saf_topic,
    select_btrt = tp_prop$group
  )
}



# plots -------------------------------------------------------------------


# reference images --------------------------------------------------------

# png(filename = paste0("img/prop/forest_plot_", tp_prop$saf_topic, ".png"))
# forest_plot
# dev.off()
#
# png(filename = paste0("img/prop/rob_map_plot_", tp_prop$saf_topic, ".png"))
# rob_map_plot
# dev.off()
#
# png(filename = paste0("img/prop/param_mix_density_", tp_prop$saf_topic, ".png"))
# param_mix_density
# dev.off()
#
# png(filename = paste0("img/prop/nta_plot_", tp_prop$saf_topic, ".png"))
# nta_plot
# dev.off()
#
# png(filename = paste0("img/prop/decision_making_plot_1_", tp_prop$saf_topic, ".png"))
# decision_making_plot[[1]]
# dev.off()
#
# png(filename = paste0("img/prop/decision_making_plot_2_", tp_prop$saf_topic, ".png"))
# decision_making_plot[[2]]
# dev.off()
#
# png(filename = paste0("img/prop/decision_making_plot_3_", tp_prop$saf_topic, ".png"))
# decision_making_plot[[3]]
# dev.off()
#
# png(filename = paste0("img/prop/decision_making_plot_4_", tp_prop$saf_topic, ".png"))
# decision_making_plot[[4]]
# dev.off()

# actual images -----------------------------------------------------------
#
png(filename = paste0("test_img/prop/forest_plot_", tp_prop$saf_topic, ".png"))
forest_plot
dev.off()

png(filename = paste0("test_img/prop/rob_map_plot_", tp_prop$saf_topic, ".png"))
rob_map_plot
dev.off()

png(filename = paste0("test_img/prop/param_mix_density_", tp_prop$saf_topic, ".png"))
param_mix_density
dev.off()

png(filename = paste0("test_img/prop/nta_plot_", tp_prop$saf_topic, ".png"))
nta_plot
dev.off()

png(filename = paste0("test_img/prop/decision_making_plot_1_", tp_prop$saf_topic, ".png"))
decision_making_plot[[1]]
dev.off()

png(filename = paste0("test_img/prop/decision_making_plot_2_", tp_prop$saf_topic, ".png"))
decision_making_plot[[2]]
dev.off()

png(filename = paste0("test_img/prop/decision_making_plot_3_", tp_prop$saf_topic, ".png"))
decision_making_plot[[3]]
dev.off()

png(filename = paste0("test_img/prop/decision_making_plot_4_", tp_prop$saf_topic, ".png"))
decision_making_plot[[4]]
dev.off()



# rates -------------------------------------------------------------------



hist_borrow_choices <- c(
  "Small", "Moderate", "Substantial", "Large",
  "Very Large"
)
path <- "scenarios/"
csv_choices <- list.files(path, pattern = "*.csv")
testing_list_rates <- readRDS("thresholds/testing_list_rates.rds")

tp_rates <- testing_list_rates[[4]]$parameters
thresholds_rates <- testing_list_rates[[4]]$treshholds

data <- read.csv(paste0(path, tp_rates$csv))

data <- data_table_prep(
  input_data = data,
  select_analysis = tp_rates$analysis,
  saf_topic = tp_rates$saf_topic,
  select_btrt = tp_rates$group,
  bool_pooled = TRUE # might change later
)

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

# rerun_counter <- 0
# while(sum(data$N_WITH_AE) < result_map$ESS[1] && rerun_counter < 2) {
#   rerun_counter <- rerun_counter + 1
#
#   map_object <- map_prior_func(
#     input_data = data,
#     select_analysis = tp_rates$analysis,
#     tau_dist = tp_rates$tau,
#     adj_tau = adj_tau, # scenario specific
#     seed = tp_rates$seed + Sys.time(),
#     testing = FALSE
#   )
#   param_approx <- parametric_approx(
#     select_analysis = tp_rates$analysis,
#     map_prior = map_object
#   )
#
#   robust_map <- robust_map(
#     select_analysis = tp_rates$analysis,
#     param_approx = param_approx,
#     input_data = data,
#     robust_weight = tp_rates$rob_weight,
#     robust_mean = tp_rates$rob_mean,
#     adj_tau = adj_tau,
#     seed = tp_rates$seed + Sys.time(),
#     testing = FALSE
#   )
#
#   post_dist <- posterior_dist(
#     select_analysis = tp_rates$analysis,
#     input_data = data,
#     robust_map_prior = robust_map,
#     new_v1 = tp_rates$nta_event,
#     new_v2 = tp_rates$nta_time,
#     seed = tp_rates$seed + Sys.time()
#   )
#
#   result_map <- model_summary_display(
#     map_object = map_object,
#     select_analysis = tp_rates$analysis,
#     param_approx = param_approx,
#     ess_method = tp_rates$ESS,
#     numerical = TRUE
#   )
# }

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


rob_comp <- robust_compare(
  select_analysis = tp_rates$analysis,
  robust_map_prior = robust_map,
  param_approx = param_approx
)

new_trial_compare <- new_trial_compare(
  select_analysis = tp_rates$analysis,
  robust_map_prior = robust_map,
  new_v1 = current_trial_data$new_v1,
  new_v2 = current_trial_data$new_v2,
  post_dist = post_dist
)

forest_plot <- forest_plot_display(
  map_object = map_object,
  select_analysis = tp_prop$analysis,
  saf_topic = tp_prop$saf_topic,
  select_btrt = tp_prop$group
)

rob_map_plot <- robust_map_prior_plot(
  rob_comp = rob_comp,
  saf_topic = tp_rates$saf_topic,
  select_btrt = tp_rates$group,
  select_analysis = tp_rates$analysis
)

param_mix_density <- param_mix_density_display(
  param_approx = param_approx,
  select_analysis = tp_rates$analysis,
  saf_topic = tp_rates$saf_topic,
  select_btrt = tp_rates$group
)

nta_plot <- nta_data_conflict_assassment_plot(
  new_trial_analysis = new_trial_compare,
  saf_topic = tp_rates$saf_topic,
  select_btrt = tp_rates$group,
  select_analysis = tp_rates$analysis
)

choices <- c("Likelihood", "MAP Prior", "Robust MAP Prior", "Posterior")
stat_inf_dist <- list()
decision_making_plot <- list()

for (i in 1:length(choices)) {
  stat_inf_dist[[i]] <- sampling_all_plot(
    select_analysis = tp_rates$analysis,
    select_dist = choices[i],
    param_approx = param_approx,
    new_trial_analysis = new_trial_compare
  )

  decision_making_plot[[i]] <- decision_making_density_plot(
    select_analysis = tp_rates$analysis,
    stat_inf_dist = stat_inf_dist[[i]],
    ae_prop = c(30, 100) / 100,
    saf_topic = tp_rates$saf_topic,
    select_btrt = tp_rates$group
  )
}



# plots -------------------------------------------------------------------


# reference images --------------------------------------------------------

# png(filename = paste0("img/rates/forest_plot_", tp_rates$saf_topic, ".png"))
# forest_plot
# dev.off()
#
# png(filename = paste0("img/rates/rob_map_plot_", tp_rates$saf_topic, ".png"))
# rob_map_plot
# dev.off()
#
# png(filename = paste0("img/rates/param_mix_density_", tp_rates$saf_topic, ".png"))
# param_mix_density
# dev.off()
#
# png(filename = paste0("img/rates/nta_plot_", tp_rates$saf_topic, ".png"))
# nta_plot
# dev.off()
#
# png(filename = paste0("img/rates/decision_making_plot_1_", tp_rates$saf_topic, ".png"))
# decision_making_plot[[1]]
# dev.off()
#
# png(filename = paste0("img/rates/decision_making_plot_2_", tp_rates$saf_topic, ".png"))
# decision_making_plot[[2]]
# dev.off()
#
# png(filename = paste0("img/rates/decision_making_plot_3_", tp_rates$saf_topic, ".png"))
# decision_making_plot[[3]]
# dev.off()
#
# png(filename = paste0("img/rates/decision_making_plot_4_", tp_rates$saf_topic, ".png"))
# decision_making_plot[[4]]
# dev.off()

# actual images -----------------------------------------------------------
#
png(filename = paste0("test_img/rates/forest_plot_", tp_rates$saf_topic, ".png"))
forest_plot
dev.off()

png(filename = paste0("test_img/rates/rob_map_plot_", tp_rates$saf_topic, ".png"))
rob_map_plot
dev.off()

png(filename = paste0("test_img/rates/param_mix_density_", tp_rates$saf_topic, ".png"))
param_mix_density
dev.off()

png(filename = paste0("test_img/rates/nta_plot_", tp_rates$saf_topic, ".png"))
nta_plot
dev.off()

png(filename = paste0("test_img/rates/decision_making_plot_1_", tp_rates$saf_topic, ".png"))
decision_making_plot[[1]]
dev.off()

png(filename = paste0("test_img/rates/decision_making_plot_2_", tp_rates$saf_topic, ".png"))
decision_making_plot[[2]]
dev.off()

png(filename = paste0("test_img/rates/decision_making_plot_3_", tp_rates$saf_topic, ".png"))
decision_making_plot[[3]]
dev.off()

png(filename = paste0("test_img/rates/decision_making_plot_4_", tp_rates$saf_topic, ".png"))
decision_making_plot[[4]]
dev.off()
