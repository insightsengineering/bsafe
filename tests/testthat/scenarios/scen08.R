# Scen8
scen8 <- SimTestData(
  SimStudy_nPat = c(g1 = 200, g2 = 200),
  SimStudy_hz = c(g1 = 0.1, g2 = 0.2),
  SimStudy_dropout = c(rate = 0.3, time = 12),
  SimStudy_accr = 6,
  SimStudy_accr_method = "Uniform",
  SimStudy_surv_method = "Exponential",
  SimStudy_intensity = NA,
  SimStudy_accr_timepoint = NA,
  SimStudy_time_cutoff = 18,
  SimStudy_NObsEvt = 93,
  SimStudy_censor_type = 2,
  nStudy = 6,
  tau = 0.15,
  prior_data_conflict = FALSE,
  pdc_hz = NA,
  SAF_TOPIC = "Scen08",
  seed = 1701652217
)


write.csv(scen8, file = "Scen08.csv")
