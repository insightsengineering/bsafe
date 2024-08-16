test_that("ae summary table function works as expected", {
  data <- read.csv("scenarios/Scen07.csv")
  cb_list_ctrl <- "g1"
  cb_list_trt <- "g2"
  saf_topic <- "Scen07"
  seed <- 8888

  ae_result <- ae_summary_table(
    input_data = data,
    cb_list_ctrl = cb_list_ctrl,
    cb_list_trt = cb_list_trt,
    saf_topic = saf_topic,
    seed = seed
  )

  expect_no_error(ae_result)
  expect_type(ae_result, "list")
})
