# - Medium Heterogeneity scenario
#
# Quantify the above conditions
# Scen5
# st6 <- Sys.time()
# as.numeric(st6)
# seed_scen6 <- 1701629671

st6data <- SimTestData(
  SimStudy_nPat = c(g1 = 200, g2 = 200),
  SimStudy_hz = c(g1 = 0.1, g2 = 0.2),
  SimStudy_dropout = c(rate = 0.3, time = 12),
  SimStudy_accr = 6,
  SimStudy_accr_method = "Uniform",
  SimStudy_surv_method = "Exponential",
  SimStudy_intensity = NA,
  SimStudy_accr_timepoint = NA,
  SimStudy_time_cutoff = 18,
  SimStudy_NObsEvt = 95,
  SimStudy_censor_type = 2,
  nStudy = 6,
  tau = 0.02,
  prior_data_conflict = FALSE,
  pdc_hz = NA,
  SAF_TOPIC = "Scen06",
  seed = 1701628373
)

write.csv(st6data, file="Scen06.csv")
