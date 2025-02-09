library(meta)
library(metafor)
library(ggplot2)
library(dplyr)
library(forcats)
library(grid)
library(ggrepel)

hrs <- read.csv("hrs_by_cohort.csv")
hrs$logHR <- log(hrs$HR_OS)
# Calculate standard error of the log HR
hrs$SE <- (log(hrs$HR_OS_95high) - log(hrs$HR_OS_95low)) / (2 * 1.96)

meta_analysis <- metagen(TE = logHR, 
                         seTE = SE, 
                         studlab = hrs$Study, 
                         sm = "HR", 
                         common = TRUE, 
                         random = FALSE, 
                         method.random.ci = TRUE, 
                         method.tau = "REML", 
                         data = hrs)

forest(meta_analysis)


# Effect sizes
effects <- exp(meta_analysis$TE)
# Confidence intervals
ci_lower <- exp(meta_analysis$lower)
ci_upper <- exp(meta_analysis$upper)

# Study labels
studies <- meta_analysis$studlab

# Weights
weights <- (meta_analysis$w.common / sum(meta_analysis$w.common)) * 100  # Use random-effects weights

# Summary effect (if needed)
summary_effect <- exp(meta_analysis$TE.common)
summary_ci_lower <- exp(meta_analysis$lower.common)
summary_ci_upper <- exp(meta_analysis$upper.common)

# Prediction interval (optional)
pred_lower <- exp(meta_analysis$lower.predict)
pred_upper <- exp(meta_analysis$upper.predict)



forest_data <- data.frame(
  Study = studies,
  Effect = effects,
  CI_Lower = ci_lower,
  CI_Upper = ci_upper,
  Weight = weights  # Summary does not have a weight
)

summary_data <- data.frame(
  Study = c("Summary"),
  Effect = summary_effect,
  CI_Lower = summary_ci_lower,
  CI_Upper = summary_ci_upper,
  pred_lower = pred_lower,
  pred_upper = pred_upper
)


forest_data <- forest_data[order(-forest_data$Effect) , ]

custom_labels <- c(
  "Tlemsani 2021" = "Tlemsani et al, 2021",
  "Ul Haq 2024" = "Wei et al, 2023",
  "Summary" = expression(bold("Summary (" * italic("I")^"2"~"= 8.3%)"))
)


p <- ggplot(forest_data, aes(x = Effect, y = fct_reorder(Study, Effect))) +
  
  # Summary effect value
  geom_vline(xintercept = summary_data$Effect, color="#636566", linetype=5) +
  
  geom_linerange(aes(xmin=CI_Lower, xmax=CI_Upper), color="#87898a") +
  geom_point(shape=22, size=3, color="#0c94b3", stroke=0.8, fill="#FFF") +
  labs(x = "Hazard ratio", y = "") +
  
  # the weights are labelled
  geom_text(aes(label = paste0(sprintf(Weight, fmt = '%#.1f'), "%"), 
                x = max(CI_Upper) + 0.1), # Offset text for space
            hjust = 0, vjust = 0.5, size = 3, color = "black") +
  
  annotate("text", 
           x = max(forest_data$CI_Upper) + 0.1, y = length(forest_data$Study) + 1, 
           label = "Weights", size = 3, color = "black", hjust = 0, fontface = "bold") +
  
  # Summary point
  geom_linerange(summary_data, mapping=aes(xmin=CI_Lower, xmax=CI_Upper), color="#87898a") +
  geom_point(summary_data, mapping=aes(x=Effect, y=Study), shape=18, size=8, color="#0c94b3") +
  
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(), # minor x axis grid line
        plot.margin = margin(5, 20, 5, 5), # top, right, bottom, left
        axis.text.y = element_text(size = 10, hjust = 0),
        axis.line.x = element_line(linewidth=0.5)) +
  coord_fixed(ratio = 0.2)  + # Adjust the aspect ratio (increasing the ratio makes the plot wider relative to its height) +
  scale_y_discrete(labels = custom_labels) # Apply custom labels

# Disable clipping
ggplot_gtable <- ggplotGrob(p)
ggplot_gtable$layout$clip[ggplot_gtable$layout$name == "panel"] <- "off"

pdf("Prettier figures/HR_by_cohort.pdf", width = 7, height = 3)
# Draw the plot
grid.draw(ggplot_gtable)
dev.off()







# Extract data from meta_analysis
funnel_data <- data.frame(
  Effect = meta_analysis$TE,         # Effect sizes
  SE = meta_analysis$seTE,           # Standard errors
  Study = meta_analysis$studlab
)

pooled_effect <- meta_analysis$TE.common

# Compute confidence limits
funnel_data <- funnel_data %>%
  mutate(
    Lower95 = 0 - qnorm(0.975) * (1 / SE),  # 95% lower bound
    Upper95 = 0 + qnorm(0.975) * (1 / SE)   # 95% upper bound
  )

# Calculate the confidence funnel bounds
funnel_data <- funnel_data %>%
  mutate(
    Lower95 = meta_analysis$TE.common - qnorm(0.975) * SE, # 95% lower confidence limit
    Upper95 = meta_analysis$TE.common + qnorm(0.975) * SE  # 95% upper confidence limit
  )



# Base plot
p <- ggplot(funnel_data, aes(x = Effect, y = 1 / SE)) + # Use 1/SE for vertical axis
  geom_point(size = 2, color = "black", alpha = 0.7) +  # Add points for studies
  
  # Add lines for symmetry
  geom_vline(xintercept = pooled_effect, linetype = "dashed", color = "grey50") +
  
  # Labels and theme
  labs(
    x = "Effect Size",
    y = "Precision (1 / SE)") +
  
  # true less straight confidence intervals
  geom_line(data = funnel_data, aes(x = Lower95, y = 1 / SE), linetype = "dotted", color = "red") +
  geom_line(data = funnel_data, aes(x = Upper95, y = 1 / SE), linetype = "dotted", color = "red") +
  
  geom_text_repel(
    aes(label = Study),          # Use Study labels
    size = 3,                    # Text size
    color = "black",
    box.padding = 0.3,           # Space around text
    point.padding = 0.3,         # Space between text and points
    max.overlaps = Inf           # Allow all labels to display
  ) +
  
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linetype = "dotted"),
    panel.grid.major.x = element_line(linetype = "dotted"))

p
ggsave("Prettier figures/HR_by_cohort_funnel.pdf")
