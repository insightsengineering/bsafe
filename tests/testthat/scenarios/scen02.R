# - Heavy Prior Data conflict
#
# Quantify the above conditions
# Scen2
# st2 <- Sys.time()
# seed_scen2 <- as.numeric(st2)
# 1701611344

st2data <- SimTestData(
  SimStudy_nPat = c(g1 = 200, g2 = 200),
  SimStudy_hz = c(g1 = 0.1, g2 = 0.3),
  SimStudy_dropout = c(rate = 0, time = 12),
  SimStudy_accr = 1,
  SimStudy_accr_method = "Uniform",
  SimStudy_surv_method = "Exponential",
  SimStudy_intensity = NA,
  SimStudy_accr_timepoint = NA,
  SimStudy_time_cutoff = 18,
  SimStudy_NObsEvt = 0.9,
  SimStudy_censor_type = 2,
  nStudy = 10,
  tau = 0.01,
  prior_data_conflict = TRUE,
  pdc_hz = c(g1 = 0.4, g2 = 0.05),
  SAF_TOPIC = "Scen02",
  seed = 1701611344
)

write.csv(st2data, file="Scen02.csv")

