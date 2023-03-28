library(crqa)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)

source("R/lorenz_harmonic.R")
source("R/utils.R")

# None of the time series need to be further embedded,
# so delay and embed are set to 1.
delay <- 1
embed <- 1

# Set radius for time series 1
radius_1 <- 0.2
# Radius for time series 2 should be fixed to give same RR, but here we
# just do a conceptual example, so we set it to the same value as ts 1.
radius_2 <- 0.2

data_length <- 1000

model_data <- lorenz_harmonic(n = data_length, skip = 1000, coupling = 0.4)


ts_1 <- extract_lorenz(model_data)
ts_2 <- extract_oscillator(model_data)

rqa_ts_1 <- mdrqa(ts_1, delay, embed, radius_1)
rp_ts_1 <- as.matrix(rqa_ts_1$RP)

rqa_ts_2 <- mdrqa(ts_2, delay, embed, radius_2)
rp_ts_2 <- as.matrix(rqa_ts_2$RP)  

# Compute the JRP as element wise product (Hadamard product)
jrp_ts_1_ts_2 <- rp_ts_1 * rp_ts_2


system_1_data <- model_data %>% 
  select(time, x, y, z) %>% 
  mutate(i = row_number()) %>% 
  pivot_longer(cols = c("x", "y", "z"), names_to = "variable", values_to = "value")

ggplot(system_1_data, aes(x = i, y = value)) +
  geom_line() +
  ylab("") +
  theme_minimal() +
  facet_wrap(~ variable, nrow = 3, strip.position = "right") +
  theme(strip.text.y = element_text(angle = 0))
ggsave("Figures/Conceptual/ts_1.pdf", width = 4, height = 4)

system_2_data <- model_data %>% 
  select(time, u, v) %>% 
  mutate(i = row_number()) %>% 
  pivot_longer(cols = c("u", "v"), names_to = "variable", values_to = "value")

ggplot(system_2_data, aes(x = i, y = value)) +
  geom_line() +
  ylab("") +
  theme_minimal() +
  facet_wrap(~ variable, nrow = 3, strip.position = "right") +
  theme(strip.text.y = element_text(angle = 0))
ggsave("Figures/Conceptual/ts_2.pdf", width = 4, height = 4)

#
# 3D plotly plot which matches ggplot aesthetically
#
scene <-  list(camera = list(eye = list(x = -1.6, y = 1.6, z = 0.4)))
ps_1 <- plot_ly(model_data, x=~x, y=~y, z=~z) %>%
  # add_markers(size=1) %>% 
  add_trace(x=~x, x=~y, z=~z,
            type="scatter3d", mode="lines",
            line = list(width=8, color="black")) %>% 
  layout(scene = scene)
print(ps_1)
# Save to PDF. This needs to be installed in a python environment
# kaleido(ps_1, "ps_1.pdf")

# For now: Click "Download plot as png" in the viewer and crop the image.

ggplot(model_data, aes(x = u, y = v)) +
  geom_path(linewidth = 1) + 
  theme_minimal()
ggsave("Figures/Conceptual/ps_2.pdf", width = 4, height = 4)


plot_rp(rp_ts_1) + 
  xlim(c(0, data_length)) +
  ylim(c(0, data_length)) +
  xlab("i") + 
  ylab("j") + 
  theme_minimal() +
  theme(legend.position = "none")
ggsave("Figures/Conceptual/rp_1.pdf", width = 4, height = 4)

plot_rp(rp_ts_2) + 
  xlim(c(0, data_length)) +
  ylim(c(0, data_length)) +
  xlab("i") + 
  ylab("j") + 
  theme_minimal() +
  theme(legend.position = "none")
ggsave("Figures/Conceptual/rp_2.pdf", width = 4, height = 4)

plot_rp(jrp_ts_1_ts_2) + 
  xlim(c(0, data_length)) +
  ylim(c(0, data_length)) +
  xlab("i") + 
  ylab("j") + 
  theme_minimal() +
  theme(legend.position = "none")
ggsave("Figures/Conceptual/jrp.pdf", width = 4, height = 4)

