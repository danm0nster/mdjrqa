# wipe work space
rm(list = ls())

# load libraries
library(crqa)
library(ggplot2)

# loop through computation and save results
beta_1 <- seq(0,10, by = 0.5)
no_samples <- 100
sample_size <- 100
target <- c(9.75,10.25) # target recurrence rate
delay_1 <- 1
embed_1 <- 1
radius_1 <- 0.22
delay_2 <- 1
embed_2 <- 1
radius_2 <- 0.22
jrqa_RR <- NULL
rqa_RR <- NULL
index_beta <- NULL
index_simu <- NULL
k <- 0

for(j in 1:length(beta_1)) {
  
  for(i in 1:no_samples) {
    
    k <- k+1
    
    # create data
    epsilon_4 <- rnorm(sample_size)
    x_1 <- rnorm(sample_size)
    y_1 <- beta_1[j] * x_1 + epsilon_4 + rnorm(sample_size)
    y_2 <- beta_1[j] * x_1 + epsilon_4 + rnorm(sample_size)
    Y <- as.matrix(cbind(y_1,y_2))
    
    
    # print(k) # print progress
    
    
    
    # set rqa parameters for (unidimensional) time series x_1
    
    
    # run rqa
    RR <- 0
    while(RR < target[1] | RR > target[2]) {
      rqa_ts_1 <-  crqa(x_1, x_1,
                        delay_1, embed_1, radius_1,
                        rescale     = 0,           # 0 (none) | 1 (mean distance) | 2 (max distance) | 3 (min distance) | 4 (euclidean distance)
                        normalize   = 2,           # 0 (none) | 1 (unit interval) | 2 (z-score)
                        mindiagline = 2,           # min recurrent points for diagonal lines
                        minvertline = 2,           # min recurrent points for vertical lines
                        tw          = 1,           # Theiler window: 0 (incl LOI) | 1 (excl ROI)
                        whiteline   = TRUE,        # calculate empty vertical lines or not
                        recpt       = FALSE,       # calculate measures directly from RP or not
                        side        = "both",      # base measures on "upper" or "lower" triangle, or "both"
                        method      = "rqa",       # "rqa" | "crqa" | "mdrqa"
                        metric      = "euclidean", # distance metric: "euclidean" | "maximum" | "minkowski" | ...
                        datatype    = "continuous")
      # print(c(1,j,rqa_ts_1$RR)) # print progress
      RR <- rqa_ts_1$RR
      if(RR < target[1]) {
        radius_1 <- radius_1 +0.004
      } else {
        radius_1 <- radius_1 -0.004
      }
    }
    
    # run rqa
    RR <- 0
    while(RR < target[1] | RR > target[2]) {
      rqa_ts_2 <-  crqa(Y, Y,
                        delay_2, embed_2, radius_2,
                        rescale     = 0,           # 0 (none) | 1 (mean distance) | 2 (max distance) | 3 (min distance) | 4 (euclidean distance)
                        normalize   = 2,           # 0 (none) | 1 (unit interval) | 2 (z-score)
                        mindiagline = 2,           # min recurrent points for diagonal lines
                        minvertline = 2,           # min recurrent points for vertical lines
                        tw          = 1,           # Theiler window: 0 (incl LOI) | 1 (excl ROI)
                        whiteline   = TRUE,        # calculate empty vertical lines or not
                        recpt       = FALSE,       # calculate measures directly from RP or not
                        side        = "both",      # base measures on "upper" or "lower" triangle, or "both"
                        method      = "mdrqa",      # "rqa" | "crqa" | "mdrqa"
                        metric      = "euclidean", # distance metric: "euclidean" | "maximum" | "minkowski" | ...
                        datatype    = "continuous")
      # print(c(2,j,rqa_ts_1$RR)) # print progress
      RR <- rqa_ts_2$RR
      if(RR < target[1]) {
        radius_2 <- radius_2 +0.004
      } else {
        radius_2 <- radius_2 -0.004
      }
    }  
    
    
    
    
    
    # converts sparse-matrix recurrence plots to standard matrix
    rp_ts_1 <- as.matrix(rqa_ts_1$RP)
    rp_ts_2 <- as.matrix(rqa_ts_2$RP)  
    
    # get joint recurrence plot
    jrp_ts_1_ts_2 <- rp_ts_1 * rp_ts_2
    
    # get joint RR
    jrqa_RR[k] <- 100*sum(jrp_ts_1_ts_2)/(sample_size^2)
    
    # get average individual RR
    rqa_RR[k] <- (rqa_ts_1$RR + rqa_ts_2$RR) / 2
    
    # get beta value
    index_beta[k] <- beta_1[j]
    
    # get simulation trial
    index_simu[k] <- k
    
  }
}

# collect data
df <- data.frame(index_beta, index_simu, jrqa_RR, rqa_RR)


# Take mean of RR and jRR
# First get the mean of the RR and JRR for beta == 0
RR_0 <- df %>% filter(index_beta == 0) %>% pull(rqa_RR) %>% mean()
JRR_0 <- df %>% filter(index_beta == 0) %>% pull(jrqa_RR) %>% mean()

df_mean <- df %>% 
  group_by(index_beta) %>% 
  summarise(jrqa_RR = mean(jrqa_RR) / JRR_0,
            rqa_RR = mean(rqa_RR) / RR_0) %>% 
  pivot_longer(cols = c("jrqa_RR", "rqa_RR"),
               names_to = "RR_type",
               values_to = "rr")

ggplot(df_mean,
       aes(x = index_beta, y = rr, shape = RR_type)) +
  geom_segment(aes(x = 0, xend = max(df$index_beta), y = 1, yend = 1),
               colour = "black", linetype = "solid") +
  geom_point(alpha = 0.8, size = 2) +
  scale_shape_manual(values = c("rqa_RR" = 16,  "jrqa_RR" = 15)) +
  geom_line(data = df_mean %>% filter(RR_type == "jrqa_RR")) +
  xlab(expression(paste(beta [1]))) +
  ylab("Relative RR") +
  scale_y_continuous(breaks = 1:6) +
  annotate("text", x =  max(df$index_beta), y = 1.2, hjust = 1,
           label = "Individual sub-system recurrence") +
  annotate("text", x =  max(df$index_beta), y = max(df_mean$rr) * 1.01,
           hjust = 1, vjust = 0,
           label = "Joint recurrence") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.line.x = element_line(arrow=arrow(length = unit(3, 'mm'))), 
        axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm')))
  )

ggsave(filename = "Figures/linear_rr_relative.pdf",
       width = 1.5 * 4, height = 1.5 * 3)


# plots of aggregated data
rqa_RR_mean <- aggregate(df$rqa_RR, list(df$index_beta), FUN=mean)[,2]
jrqa_RR_mean <- aggregate(df$jrqa_RR, list(df$index_beta), FUN=mean)[,2]
index_beta_group <- aggregate(df$rqa_RR, list(df$index_beta), FUN=mean)[,1]

pdf("Figures/linear_rr.pdf")
par(mfrow = c(2,1))
plot(index_beta_group,rqa_RR_mean,
     type = "o",
     ylim=c(9.75,10.25),
     xlab = expression(paste(beta [1])),
     ylab = "RR",
     main = "Individual sub-system recurrence")
plot(index_beta_group,jrqa_RR_mean,
     type = "o",
     ylim=c(0,7),
     xlab = expression(paste(beta [1])),
     ylab = "RR",
     main = "Joint recurrence")
dev.off()

