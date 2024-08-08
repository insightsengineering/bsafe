# - Heavy Prior Data conflict
#
# Quantify the above conditions
# Scen3
# st3 <- Sys.time()
# seed_scen3 <- as.numeric(st3)
# 1701621871

st3data <- SimTestData(
  SimStudy_nPat = c(g1 = 200, g2 = 200),
  SimStudy_hz = c(g1 = 0.1, g2 = 0.2),
  SimStudy_dropout = c(rate = 0.05, time = 12),
  SimStudy_accr = 6,
  SimStudy_accr_method = "Uniform",
  SimStudy_surv_method = "Exponential",
  SimStudy_intensity = NA,
  SimStudy_accr_timepoint = NA,
  SimStudy_time_cutoff = 18,
  SimStudy_NObsEvt = 93,
  SimStudy_censor_type = 2,
  nStudy = 6,
  tau = 0.02,
  prior_data_conflict = FALSE,
  pdc_hz = NA,
  SAF_TOPIC = "Scen03",
  seed = 1701621871
)

write.csv(st3data, file = "Scen03.csv")
