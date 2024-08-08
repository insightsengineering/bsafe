test_that("test-ae_summary_table_helpers.R works as expected", {

  data <- read.csv("scenarios/Scen07.csv")
  cb_list_ctrl <- "g1"
  cb_list_trt <- "g2"
  saf_topic <- "Scen07"
  seed <- 8888

  trt_trials <- data_table_prep(
    input_data = data,
    select_analysis = "Incidence proportion",
    saf_topic = saf_topic,
    select_btrt = "g2"
  )

  ctr_trials <- data_table_prep(
    input_data = data,
    select_analysis = "Incidence proportion",
    saf_topic = saf_topic,
    select_btrt = "g1"
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

  expect_no_error(data_check_bin)

  dc_bin <- sum(data_check_bin)

  expect_type(dc_bin, "integer")

  warn_txt <- data.frame(
    Issue = character(), Topic = character(),
    Analysis = character(), Message = character(),
    Group = character(), stringsAsFactors = FALSE
  )

  warn_txt <- add_row(
    df = warn_txt,
    analysis = "Incidence proportion",
    data_check = data_check_bin,
    group = 1,
    topic = saf_topic
  )
  expect_no_error(warn_txt)

})
