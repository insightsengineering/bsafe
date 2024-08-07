# Huge censoring
# Huge Noise
# Little events observed, Power = 0.9, not regarded
# Heterogeneous Data
# Huge data conflict

# Quantify the above conditions
# Scen12
scen12 <- SimTestData(
  SimStudy_nPat = c(g1 = 200, g2 = 200),
  SimStudy_hz = c(g1 = 0.1, g2 = 0.2),
  SimStudy_dropout = c(rate = 0.05, time = 12),
  SimStudy_accr = 6,
  SimStudy_accr_method = "Uniform",
  SimStudy_surv_method = "Exponential",
  SimStudy_intensity = NA,
  SimStudy_accr_timepoint = NA,
  SimStudy_time_cutoff = NA,
  SimStudy_NObsEvt = 400,
  SimStudy_censor_type = 2,
  nStudy = 6,
  tau = 0.15,
  prior_data_conflict = TRUE,
  pdc_hz = c(g1 = 0.05, g2 = 0.1),
  SAF_TOPIC = "Scen12",
  seed = 1701878308
)


write.csv(scen12, file = "Scen12.csv")
