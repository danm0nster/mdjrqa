# load libraries
library(crqa)
library(tidyr)
library(dplyr)
library(ggplot2)

source("R/lorenz_harmonic.R")
source("R/utils.R")

# None of the time series need to be further embedded,
# so delay and embed are set to 1.
delay <- 1
embed <- 1

# Set radius for time series 1
radius_1 <- 0.2


rp_1_list <- list()
rp_2_list <- list()
rp_j_list <- list()
lagged_rr_plot <- list()
joint_rr <- data.frame()

# Loop over values of the coupling constant, cc.
couplings <- seq(0, 0.38, 0.02)

for (cc in couplings) {
  # Run the model to produce data
  model_data <- lorenz_harmonic(n = 1000, skip = 100, coupling = cc)
  
  # Extract the two time series from model_data
  ts_1 <- extract_lorenz(model_data)
  ts_2 <- extract_oscillator(model_data)
  
  rqa_ts_1 <- mdrqa(ts_1, delay, embed, radius_1)
  
  ## Find the radius that gives the same RR for ts_2 as for ts_1
  
  # This function computed the RR for a time series as a function of
  # radius, r. Note that delay and embed come from the calling environment
  # so they must be defined in this scope.
  rr <- function(time_series, r) {
    rqa <- mdrqa(time_series, delay, embed, r)
    return(rqa$RR)
  }
  
  # This is the function whose root must be found to get the radius for ts_2
  f <- function(r) rr(ts_2, r) - rqa_ts_1$RR
  
  root_2 <- uniroot(f, c(0.01, 2), maxiter = 10)
  radius_2 <- root_2$root
  rqa_ts_2 <- mdrqa(ts_2, delay, embed, radius_2)
  
  rr_1 <- rqa_ts_1$RR
  rr_2 <- rqa_ts_2$RR
  
  # converts sparse-matrix recurrence plots to standard matrix
  rp_ts_1 <- as.matrix(rqa_ts_1$RP)
  rp_ts_2 <- as.matrix(rqa_ts_2$RP)  
  
  # Compute the JRP as element wise product (Hadamard product)
  jrp_ts_1_ts_2 <- rp_ts_1 * rp_ts_2
  
  rr_j <- recurrence_rate(jrp_ts_1_ts_2)
  joint_rr <- rbind(joint_rr,
                    data.frame(c = cc, rr_1 = rr_1, rr_2 = rr_2, rr_j = rr_j))
  
  rp_1_list[[length(rp_1_list) + 1]] <- plot_rp(rp_ts_1)
  rp_2_list[[length(rp_2_list) + 1]] <- plot_rp(rp_ts_2)
  rp_j_list[[length(rp_j_list) + 1]] <- plot_rp(jrp_ts_1_ts_2)
  
  # plot_rp(rp_ts_1, plot_title = "MdRP for Lorenz system")
  # plot_rp(rp_ts_2, plot_title = "MdRP for perturbed harmonic oscillator")
  # plot_rp(jrp_ts_1_ts_2, plot_title = "MdJRP for coupled systems")
  
  
  # plot lagged profile
  lagged_RR <- lagged_joint_rr(rp_ts_1, rp_ts_2)
  second_largest_value <- sort(lagged_RR$RR, decreasing = TRUE)[2]
  lagged_rr_plot[[length(lagged_rr_plot) + 1]] <- ggplot(lagged_RR,
                                                         aes(x = w, y = RR)) +
    geom_line() +
    geom_point() +
    labs(title = paste("c = ", cc, "RR = ", rr_j)) +
    # annotate("text", x = -9, y = second_largest_value,
    #          hjust = 0,
    #          label = paste("c = ", cc, "RR = ", rr_j)) +
    theme_classic()

} # End loop over coupling



# Normalize by dividing by value at first point (c = 0)
joint_rr_norm <- joint_rr %>% 
  mutate(rr_1 = rr_1 / rr_1[1], rr_2 = rr_2 / rr_2[1], rr_j = rr_j / rr_j[1])

joint_rr_norm_long <- joint_rr_norm %>%
  pivot_longer(cols = c(rr_1, rr_2, rr_j),
               names_to = "RR_type",
               names_prefix = "rr_",
               values_to = "rr")

ggplot(joint_rr_norm_long,
       aes(x = c, y = rr, shape = RR_type)) +
  # geom_hline(yintercept = 1, colour = "grey") +
  geom_segment(aes(x = 0, xend = max(joint_rr$c), y = 1, yend = 1),
               colour = "black", linetype = "solid") +
  geom_point(alpha = 0.8) +
  scale_shape_manual(values = c("1" = 16, "2" = 16, "j" = 15)) +
  geom_line(data = joint_rr_norm_long %>% filter(RR_type == "j")) +
  xlab("coupling strength (c)") +
  ylab("Relative RR") +
  scale_y_continuous(breaks = c(1.0, 1.2, 1.4, 1.6)) +
  # labs(title = "RR relative to c = 0") +
  annotate("text", x =  max(joint_rr$c), y = 1.05, hjust = 1,
           label = "Individual sub-system recurrence") +
  annotate("text", x =  max(joint_rr$c), y = max(joint_rr_norm$rr_j) * 1.01,
           hjust = 1, vjust = 0,
           label = "Joint recurrence") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.line.x = element_line(arrow=arrow(length = unit(3, 'mm'))), 
axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm')))
)
ggsave(filename = "Figures/joint_and_single_rr_relative.pdf",
       width = 1.5 * 4, height = 1.5 * 3)
