
library(tidyverse)
library(ggplot2)
library(dbplyr)

saanich_data <- read.csv("Saanich_Data.csv")

#1Summarize the data to calculate mean and standard deviation for each depth and chemical using group_by() and summarise() :
  
  summary_data <- saanich_data %>%
    group_by(Depth) %>%
    summarise(
      mean_WS_H2S = mean(WS_H2S, na.rm = TRUE),
      sd_WS_H2S = sd(WS_H2S, na.rm = TRUE),
      mean_WS_NO3 = mean(WS_NO3, na.rm = TRUE),
      sd_WS_NO3 = sd(WS_NO3, na.rm = TRUE) ,
      mean_WS_O2 = mean(WS_O2, na.rm = TRUE),
      sd_WS_O2 = sd(WS_O2, na.rm = TRUE),
    )
  
#2 pivot longer to get standard deviations and means

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

#3 Combine mean and standard deviation data
combined_data <- left_join(mean_data, sd_data, by = c("Depth", "Chemical"))

#4 Recode Location values
combined_data <- combined_data %>%
  select(-sd_WS_H2S, -sd_WS_NO3, -sd_WS_O2, 
         -mean_WS_H2S, -mean_WS_NO3, -mean_WS_O2) %>%
  filter(Concentration_uM != "NaN")

#5 Order your data by Depth

combined_data <- combined_data %>% arrange(Chemical, Depth)

#6 Plot
  
  ggplot(combined_data, aes(x = Concentration_uM, y = Depth, colour = Chemical)) +
  geom_point() +
  geom_line(aes(group = Chemical), orientation = "y") + # This line connects your dots with a line using a Y-axis order
  facet_grid(. ~ Chemical, scales = "free_x") +
  scale_y_reverse() + 
  geom_errorbar(aes(xmin = Concentration_uM - sd_Concentration_uM, 
                    xmax = Concentration_uM + sd_Concentration_uM)) + #This is how you add your error bars
  theme_bw()

 # smooth the data using a generalized additive model (GAM):

ggplot(combined_data, aes(x = Concentration_uM, y = Depth, colour = Chemical)) +
  geom_point() +
  geom_line(aes(group = Chemical), orientation = "y") +  # Ensure it's connecting by Depth
facet_grid(. ~ Chemical, scales = "free_x") +
  scale_y_reverse() + 
  geom_errorbar(aes(xmin = Concentration_uM - sd_Concentration_uM, 
                    xmax = Concentration_uM + sd_Concentration_uM)) +
  geom_smooth(method = "gam", se = TRUE) +
  theme_bw()

no3_data <- combined_data %>%
  filter(Chemical == "NO3")

ggplot(no3_data, 
       aes(x = Concentration_uM, y = Depth, colour = Chemical)) +
  geom_point() +
  geom_line(aes(group = Chemical), orientation = "y") +
  facet_grid(. ~ Chemical, scales = "free_x") +
  scale_y_reverse() + 
  geom_errorbar(aes(xmin = Concentration_uM - sd_Concentration_uM, 
                    xmax = Concentration_uM + sd_Concentration_uM)) +
  theme_bw()

ggplot(no3_data, 
       aes(x = Concentration_uM, y = Depth, colour = Chemical)) +
  geom_point() +
  geom_line(aes(group = Chemical), orientation = "y") +
  facet_grid(. ~ Chemical, scales = "free_x") +
  scale_y_reverse() + 
  geom_errorbar(aes(xmin = Concentration_uM - sd_Concentration_uM, 
                    xmax = Concentration_uM + sd_Concentration_uM)) +
  geom_smooth(method = "gam", se = TRUE) +
  theme_bw()
