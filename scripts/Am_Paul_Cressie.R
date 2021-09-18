PaC2011_pred <- read.csv("./data/Am_Paul_Cressie_2011_pred.csv",
                         head = TRUE)
PaC2011_pred$Scheme <- ifelse(PaC2011_pred$p_mu < 0.6, "1", "2")
PaC2011_pred$variable <- "p_mu"
PaC2011_pred$p_mu <- PaC2011_pred$p_mu * 10^5
PaC2011_RMSPE <- read.csv("./data/Am_Paul_Cressie_2011_RMSPE.csv",
                          head = TRUE)
PaC2011_RMSPE$Scheme <- ifelse(PaC2011_RMSPE$RMSPE_mu < 2, "1", "2")
PaC2011_RMSPE$variable <- "RMSPE_mu"
PaC2011_RMSPE$RMSPE_mu <- PaC2011_RMSPE$RMSPE_mu * 10^4
PaC2011_pred <- dplyr::rename(PaC2011_pred, "value" = "p_mu")
PaC2011_RMSPE <- dplyr::rename(PaC2011_RMSPE, "value" = "RMSPE_mu")
PaC2011 <- rbind(PaC2011_pred, PaC2011_RMSPE)
PaC2011$Framework <- "Paul and Cressie (2011): log-normal block-kriging"

write.csv(PaC2011,
          file = "./data/Am_Paul_Cressie_2011_block_results.csv",
          row.names = FALSE)
