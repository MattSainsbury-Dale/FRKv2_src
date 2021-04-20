Am_data <- readRDS("./intermediates/Am_data.rds")
GZ_df   <- readRDS("./intermediates/Am_GZ.rds")


# ---- Figure 1 from Paul and Cressie (2011) ----

## Data plot
p_locations <- p_basic +
  geom_point(shape = 4, size = 2) +
  geom_point(data = GZ_df, shape = 8, size = 5)

# Paul and Cressie (2011) let the diameter scale with Americium.
p_data_bubble <- p_basic +
  geom_point(aes(size = Americium), shape = 1)


# ---- Figure 2 from Paul and Cressie (2011) ----

## Data plot
p_hist <- ggplot(data = Am_data) +
  geom_histogram(aes(Americium), bins = 10, color="black", fill="white") +
  theme_bw()

p_hist_log <- ggplot(data = Am_data) +
  geom_histogram(aes(log(Americium)), bins = 10, color="black", fill="white") +
  theme_bw()

# ---- Figure 3 from Paul and Cressie (2011) ----

# Now, create covariates
d_cutoff <- 30.48  

## Sanity check: should be 6 observations within d_cutoff
sum(Am_data$d <= d_cutoff)

## Sanity check: check the prediction close to ground zero
subset(Am_data, Easting == round(GZ_df$Easting, 1), Northing == round(GZ_df$Northing, 1))
Am_data$d[which.min(Am_data$d)]
Am_data$Am[which.min(Am_data$d)]

## Fit the linear model
r.lm <- lm(log(Americium) ~ -1 + x1 + x2 + x3 + x4, Am_data)
Am_data$resids <- residuals(r.lm)

## plot
p_trend <- ggplot(Am_data, aes(x = d, y = log(Americium))) +
  geom_point(aes(shape = (d <= d_cutoff)), size = 4) +
  scale_shape(solid = FALSE) +
  geom_segment(x = min(Am_data$d), y = coef(r.lm)[1] + coef(r.lm)[2] * min(Am_data$d),
               xend = 25, yend = coef(r.lm)[1] + coef(r.lm)[2] * 25) +
  geom_segment(x = 30.48, y = coef(r.lm)[3] + coef(r.lm)[4] * 30.48,
               xend = max(Am_data$d), yend = coef(r.lm)[3] + coef(r.lm)[4] * max(Am_data$d)) +
  labs(shape = "Within 30.48m from GZ",
       x = "Distance from Ground zero (m)",
       y = "log(Am) concetrations") +
  theme_bw()

p_residualsEast <- ggplot(Am_data, aes(x = Easting / 10^5, y = resids)) +
  geom_point(aes(shape = (d <= d_cutoff)), size = 4) +
  scale_shape(solid = FALSE) +
  lab1 + theme_bw() + theme(legend.position = "none")

p_residualsNorth <- ggplot(Am_data, aes(x = Northing / 10^5, y = resids)) +
  geom_point(aes(shape = (d <= d_cutoff)), size = 4) +
  scale_shape(solid = FALSE) + theme_bw() + theme(legend.position = "none")

p_histResiduals <- ggplot(data = Am_data) +
  geom_histogram(aes(resids), bins = 10, color="black", fill="white") +
  theme_bw()



