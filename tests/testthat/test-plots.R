test_that("forest plot work as expected in proportions", {

# props -------------------------------------------------------------------


  # forest plot
  expect_true(
    compare_file_binary(
      paste0("img/prop/forest_plot_", tp_prop$saf_topic, ".png"),
      paste0("test_img/prop/forest_plot_", tp_prop$saf_topic, ".png")
    )
  )
})

test_that("rob map plot work as expected in proportions", {
  expect_true(
    compare_file_binary(
      paste0("img/prop/rob_map_plot_", tp_prop$saf_topic, ".png"),
      paste0("test_img/prop/rob_map_plot_", tp_prop$saf_topic, ".png")
    )
  )
})

test_that("param_mix plot work as expected in proportions", {
  expect_true(
    compare_file_binary(
      paste0("img/prop/param_mix_density_", tp_prop$saf_topic, ".png"),
      paste0("test_img/prop/param_mix_density_", tp_prop$saf_topic, ".png")
    )
  )
})
expect_true(
  compare_file_binary(
    paste0("img/prop/nta_plot_", tp_prop$saf_topic, ".png"),
    paste0("test_img/prop/nta_plot_", tp_prop$saf_topic, ".png")
  )
)
test_that("decision making plots work as expected in proportions", {
  expect_true(
    compare_file_binary(
      paste0("img/prop/decision_making_plot_1_", tp_prop$saf_topic, ".png"),
      paste0("test_img/prop/decision_making_plot_1_", tp_prop$saf_topic, ".png")
    )
  )

  expect_true(
    compare_file_binary(
      paste0("img/prop/decision_making_plot_2_", tp_prop$saf_topic, ".png"),
      paste0("test_img/prop/decision_making_plot_2_", tp_prop$saf_topic, ".png")
    )
  )

  expect_true(
    compare_file_binary(
      paste0("img/prop/decision_making_plot_3_", tp_prop$saf_topic, ".png"),
      paste0("test_img/prop/decision_making_plot_3_", tp_prop$saf_topic, ".png")
    )
  )

  expect_true(
    compare_file_binary(
      paste0("img/prop/decision_making_plot_4_", tp_prop$saf_topic, ".png"),
      paste0("test_img/prop/decision_making_plot_4_", tp_prop$saf_topic, ".png")
    )
  )


# rates -------------------------------------------------------------------

  expect_true(
    compare_file_binary(
      paste0("img/rates/forest_plot_", tp_prop$saf_topic, ".png"),
      paste0("test_img/rates/forest_plot_", tp_prop$saf_topic, ".png")
    )
  )
})

test_that("rob map plot work as expected in proportions", {
  expect_true(
    compare_file_binary(
      paste0("img/rates/rob_map_plot_", tp_prop$saf_topic, ".png"),
      paste0("test_img/rates/rob_map_plot_", tp_prop$saf_topic, ".png")
    )
  )
})

test_that("param_mix plot work as expected in proportions", {
  expect_true(
    compare_file_binary(
      paste0("img/rates/param_mix_density_", tp_prop$saf_topic, ".png"),
      paste0("test_img/rates/param_mix_density_", tp_prop$saf_topic, ".png")
    )
  )
})
expect_true(
  compare_file_binary(
    paste0("img/rates/nta_plot_", tp_prop$saf_topic, ".png"),
    paste0("test_img/rates/nta_plot_", tp_prop$saf_topic, ".png")
  )
)
test_that("decision making plots work as expected in proportions", {
  expect_true(
    compare_file_binary(
      paste0("img/rates/decision_making_plot_1_", tp_prop$saf_topic, ".png"),
      paste0("test_img/rates/decision_making_plot_1_", tp_prop$saf_topic, ".png")
    )
  )

  expect_true(
    compare_file_binary(
      paste0("img/rates/decision_making_plot_2_", tp_prop$saf_topic, ".png"),
      paste0("test_img/rates/decision_making_plot_2_", tp_prop$saf_topic, ".png")
    )
  )

  expect_true(
    compare_file_binary(
      paste0("img/rates/decision_making_plot_3_", tp_prop$saf_topic, ".png"),
      paste0("test_img/rates/decision_making_plot_3_", tp_prop$saf_topic, ".png")
    )
  )

  expect_true(
    compare_file_binary(
      paste0("img/rates/decision_making_plot_4_", tp_prop$saf_topic, ".png"),
      paste0("test_img/rates/decision_making_plot_4_", tp_prop$saf_topic, ".png")
    )
  )
})
