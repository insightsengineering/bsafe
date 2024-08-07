# Low censoring
# Small Noise
# Majority Events observed (for beta 0.9 we would need 93, but censor after 200)
# Quite Homogenous DataQuantify the above conditions

# Scen9
scen9 <- SimTestData(
  SimStudy_nPat = c(g1 = 300, g2 = 300),
  SimStudy_hz = c(g1 = 0.1, g2 = 0.2),
  SimStudy_dropout = c(rate = 0, time = 12),
  SimStudy_accr = 6,
  SimStudy_accr_method = "Uniform",
  SimStudy_surv_method = "Exponential",
  SimStudy_intensity = NA,
  SimStudy_accr_timepoint = NA,
  SimStudy_time_cutoff = 24,
  SimStudy_NObsEvt = 0.999,
  SimStudy_censor_type = 1,
  nStudy = 8,
  tau = 0.01,
  prior_data_conflict = FALSE,
  pdc_hz = NA,
  SAF_TOPIC = "Scen09",
  seed = 1701655293
)


write.csv(scen9, file = "Scen09.csv")
