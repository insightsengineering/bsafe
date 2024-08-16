# Scen1
# st1 <- Sys.time()
# as.numeric(st1)
seed_scen1 <- 1699874539

# Scen1
st1data <- SimTestData(
  SimStudy_nPat = c(g1 = 300, g2 = 300),
  SimStudy_hz = c(g1 = 0.1, g2 = 0.2),
  SimStudy_dropout = c(rate = 0, time = 12),
  SimStudy_accr = 1,
  SimStudy_accr_method = "Uniform",
  SimStudy_surv_method = "Exponential",
  SimStudy_intensity = NA,
  SimStudy_accr_timepoint = NA,
  SimStudy_time_cutoff = 18,
  SimStudy_NObsEvt = 0.99,
  SimStudy_censor_type = 2,
  nStudy = 10,
  tau = 0,
  prior_data_conflict = FALSE,
  pdc_hz = NA,
  SAF_TOPIC = "Scen01",
  seed = 1699874539
)

write.csv(st1data, file = "Scen01.csv")
