# - Worst case scenario
#
# Quantify the above conditions
# Scen4
# st4 <- Sys.time()
# seed_scen4 <- as.numeric(st4)
# 1701626683

st4data <- SimTestData(
  SimStudy_nPat = c(g1 = 50, g2 = 100),
  SimStudy_hz = c(g1 = 0.1, g2 = 0.2),
  SimStudy_dropout = c(rate = 0.2, time = 12),
  SimStudy_accr = 6,
  SimStudy_accr_method = "Uniform",
  SimStudy_surv_method = "Exponential",
  SimStudy_intensity = NA,
  SimStudy_accr_timepoint = NA,
  SimStudy_time_cutoff = 18,
  SimStudy_NObsEvt = 112,
  SimStudy_censor_type = 2,
  nStudy = 3,
  tau = 0.15,
  prior_data_conflict = TRUE,
  pdc_hz = c(g1 = 0.05, g2 = 0.1),
  SAF_TOPIC = "Scen04",
  seed = 1701626683
)

write.csv(st4data, file="Scen04.csv")
