# - Various Treatment Length
#
# Quantify the above conditions
# Scen13
# st14 <- Sys.time()
seed_scen14 <- as.numeric(st14)
# 1718356066

st14data <- SimTestData(
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
  SimStudy_censor_type = 1,
  nStudy = 6,
  tau = 0.02,
  prior_data_conflict = FALSE,
  diff_trt_length = TRUE,
  pdc_hz = NA,
  SAF_TOPIC = "Scen13",
  seed = 1718356066
)

write.csv(st14data, file="Scen13.csv")

