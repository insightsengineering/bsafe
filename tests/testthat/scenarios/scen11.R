# Censoring: dropout = 0.05
# Small noise: tau = 0.02
# Power = 90%, calculate total # of events observed, N = 93, not regarded nevertheless studyy duration for 2 Years
# Homogeneous Historical Data
# No planned Prior Data conflict

# Quantify the above conditions
# Scen11
scen11 <- SimTestData(
  SimStudy_nPat = c(g1 = 200, g2 = 200),
  SimStudy_hz = c(g1 = 0.1, g2 = 0.2),
  SimStudy_dropout = c(rate = 0.05, time = 12),
  SimStudy_accr = 6,
  SimStudy_accr_method = "Uniform",
  SimStudy_surv_method = "Exponential",
  SimStudy_intensity = NA,
  SimStudy_accr_timepoint = NA,
  SimStudy_time_cutoff = 24,
  SimStudy_NObsEvt = 93,
  SimStudy_censor_type = 1,
  nStudy = 6,
  tau = 0.02,
  prior_data_conflict = FALSE,
  pdc_hz = NA,
  SAF_TOPIC = "Scen11",
  seed = 1701876972
)


write.csv(scen11, file = 'Scen11.csv')
