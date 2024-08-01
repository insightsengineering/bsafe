# - Medium Heterogeneity scenario
#
# Quantify the above conditions
# Scen5
# st5 <- Sys.time()
# as.numeric(st5)
seed_scen5 <- 1701628373

st5data <- SimTestData(
  SimStudy_nPat = c(g1 = 200, g2 = 200),
  SimStudy_hz = c(g1 = 0.1, g2 = 0.2),
  SimStudy_dropout = c(rate = 0.05, time = 12),
  SimStudy_accr = 1,
  SimStudy_accr_method = "Uniform",
  SimStudy_surv_method = "Exponential",
  SimStudy_intensity = NA,
  SimStudy_accr_timepoint = NA,
  SimStudy_time_cutoff = 18,
  SimStudy_NObsEvt = 93,
  SimStudy_censor_type = 2,
  nStudy = 6,
  tau = 0.05,
  prior_data_conflict = FALSE,
  pdc_hz = NA,
  SAF_TOPIC = "Scen05",
  seed = 1701628373
)

write.csv(st5data, file="Scen05.csv")
