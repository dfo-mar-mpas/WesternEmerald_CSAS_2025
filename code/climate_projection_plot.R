#make the figure of climate projections for WEBMR based on the ensemble extraction from CMIP 6 (Lewis et al. 2023)
#extracts generate from CMIP outputs in MAR_thermal_emerg project (Lewis et al) and code/cmip_extract.R

#load libraries
library(tidyverse)
library(zoo)
library(scales)



#load the outputs from cmip_extract

load("output/climate_extracts/climate_extracts_ensemble_2-6.RData")
extract_26 <- cmip_extract
rm(cmip_extract)

load("output/climate_extracts/climate_extracts_ensemble_8-5.RData")
extract_85 <- cmip_extract
rm(cmip_extract)

cmip_extract <- rbind(extract_26,extract_85)

mean_temp <- cmip_extract %>%
  group_by(climate_proj, year, month) %>%
  summarise(meant = weighted.mean(temp, cell_area, na.rm = TRUE), .groups = "drop") %>%  # monthly, area-weighted
  group_by(climate_proj, year) %>%                                                        # keep climate_proj through to the annual summary
  summarise(n     = n(),
            sd    = sd(meant, na.rm = TRUE),
            meant = mean(meant, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(se = sd / sqrt(n)) %>%                                                            # matches "± s.e." in the caption
  group_by(climate_proj) %>%
  arrange(year) %>%
  mutate(smooth5 = rollmean(meant, k = 5, fill = NA, align = "center")) %>%                # 5-yr smoothed line
  ungroup() %>%
  data.frame()


#create the plot

p1 <- ggplot(mean_temp, aes(x = year)) +
  # shaded background marking the baseline period (drawn per panel via facet)
  annotate("rect", xmin = 2015, xmax = 2025, ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.4) +
  # ± s.e. ribbon around the raw annual mean (jagged, not smoothed)
  geom_ribbon(aes(ymin = meant - se, ymax = meant + se, fill = climate_proj,col=climate_proj),
              alpha = 0.25, show.legend = FALSE) +
  # raw annual means as points
  geom_point(aes(y = meant, color = climate_proj,fill=climate_proj), alpha=0.5,shape=21, size = 1.3, show.legend = FALSE) +
  # 5-year smoothed line
  geom_line(aes(y = smooth5, color = climate_proj), linewidth = 1, show.legend = FALSE) +
  facet_wrap(~climate_proj, labeller = labeller(climate_proj = c("2.6" = "RCP 2.6", "8.5" = "RCP 8.5"))) +
  scale_color_manual(values = c("2.6" = "#3B7DD8", "8.5" = "#D8483B")) +
  scale_fill_manual(values = c("2.6" = "#3B7DD8", "8.5" = "#D8483B")) +
  labs(x = NULL, y = "Temperature (°C) ± se") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank())

p1_french <- ggplot(mean_temp, aes(x = year)) +
  # shaded background marking the baseline period (drawn per panel via facet)
  annotate("rect", xmin = 2015, xmax = 2025, ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.4) +
  # ± s.e. ribbon around the raw annual mean (jagged, not smoothed)
  geom_ribbon(aes(ymin = meant - se, ymax = meant + se, fill = climate_proj,col=climate_proj),
              alpha = 0.25, show.legend = FALSE) +
  # raw annual means as points
  geom_point(aes(y = meant, color = climate_proj,fill=climate_proj), alpha=0.5,shape=21, size = 1.3, show.legend = FALSE) +
  # 5-year smoothed line
  geom_line(aes(y = smooth5, color = climate_proj), linewidth = 1, show.legend = FALSE) +
  facet_wrap(~climate_proj, labeller = labeller(climate_proj = c("2.6" = "RCP 2.6", "8.5" = "RCP 8.5"))) +
  scale_color_manual(values = c("2.6" = "#3B7DD8", "8.5" = "#D8483B")) +
  scale_fill_manual(values = c("2.6" = "#3B7DD8", "8.5" = "#D8483B")) +
  labs(x = NULL, y = "Température (°C) ± e.s.") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank())

#save plots
ggsave("output/Figure12_eng.jpg",p1,units="in",dpi=300)
ggsave("output/Figure12_fr.jpg",p1_french,units="in",dpi=300)
