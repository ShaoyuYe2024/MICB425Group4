library(tidyverse)
library(ggplot2)
library(dbplyr)

# Load data
saanich_data <- read.csv("Saanich_Data.csv")

# Convert Depth from kilometers to meters
saanich_data <- saanich_data %>%
  mutate(Depth = Depth * 1000)

# 1. Summarize the data to calculate mean and standard deviation for each depth and chemical
summary_data <- saanich_data %>%
  group_by(Depth) %>%
  summarise(
    mean_WS_H2S = mean(WS_H2S, na.rm = TRUE),
    sd_WS_H2S = sd(WS_H2S, na.rm = TRUE),
    mean_WS_NO3 = mean(WS_NO3, na.rm = TRUE),
    sd_WS_NO3 = sd(WS_NO3, na.rm = TRUE),
    mean_WS_O2 = mean(WS_O2, na.rm = TRUE),
    sd_WS_O2 = sd(WS_O2, na.rm = TRUE)
  )

# 2. Pivot longer to get standard deviations and means
mean_data <- summary_data %>%
  pivot_longer(
    cols = starts_with("mean_WS_"),
    names_to = "Chemical",
    names_prefix = "mean_WS_",
    values_to = "Concentration_uM"
  )

sd_data <- summary_data %>%
  pivot_longer(
    cols = starts_with("sd_WS_"),
    names_to = "Chemical",
    names_prefix = "sd_WS_",
    values_to = "sd_Concentration_uM"
  )

# 3. Combine mean and standard deviation data
combined_data <- left_join(mean_data, sd_data, by = c("Depth", "Chemical"))

# 4. Filter out NaN values
combined_data <- combined_data %>%
  filter(!is.nan(Concentration_uM))

# 5. Order data by Chemical and Depth
combined_data <- combined_data %>%
  arrange(Chemical, Depth)

# 6. Plot all chemicals (O2, NO3, H2S)
ggplot(combined_data, aes(x = Concentration_uM, y = Depth, colour = Chemical)) +
  geom_point() +
  geom_line(aes(group = Chemical), orientation = "y") +
  facet_grid(. ~ Chemical, scales = "free_x") +
  scale_y_reverse() +
  geom_errorbar(aes(xmin = Concentration_uM - sd_Concentration_uM,
                    xmax = Concentration_uM + sd_Concentration_uM)) +
  labs(y = "Depth (m)", x = "Concentration (µM)") +
  theme_bw()

# Plot with GAM smoothing
ggplot(combined_data, aes(x = Concentration_uM, y = Depth, colour = Chemical)) +
  geom_point() +
  geom_line(aes(group = Chemical), orientation = "y") +
  facet_grid(. ~ Chemical, scales = "free_x") +
  scale_y_reverse() +
  geom_errorbar(aes(xmin = Concentration_uM - sd_Concentration_uM,
                    xmax = Concentration_uM + sd_Concentration_uM)) +
  geom_smooth(method = "gam", se = TRUE) +
  labs(y = "Depth (m)", x = "Concentration (µM)") +
  theme_bw()
