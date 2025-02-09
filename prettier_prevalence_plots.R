library(meta)
library(metafor)
library(ggplot2)
library(dplyr)
library(forcats)
library(grid)
library(ggrepel)

blah <- read.csv("prevalence.by.cohort.csv")

meta_result <- metaprop(
  event = blah[, 5],
  n = blah[, 4],
  studlab = blah[, 2],   # Study name
  sm = "PLOGIT",         # Logit transformation for proportions
  method = "Inverse",    # Inverse variance method
  common = F,     # Perform fixed-effects meta-analysis
  random = TRUE,    # Perform random-effects meta-analysis
  method.random.ci = TRUE,           # Use Hartung-Knapp adjustment
  prediction = TRUE      # Include prediction interval
)


# Effect sizes
effects <- plogis(meta_result$TE)
# Confidence intervals
ci_lower <- plogis(meta_result$lower)
ci_upper <- plogis(meta_result$upper)

# Study labels
studies <- meta_result$studlab

# Weights
weights <- (meta_result$w.fixed / sum(meta_result$w.fixed)) * 100  # Use random-effects weights

# Summary effect (if needed)
summary_effect <- plogis(meta_result$TE.random)
summary_ci_lower <- plogis(meta_result$lower.random)
summary_ci_upper <- plogis(meta_result$upper.random)

# Prediction interval (optional)
pred_lower <- plogis(meta_result$lower.predict)
pred_upper <- plogis(meta_result$upper.predict)

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
  "Tlemsani 2021 (George cohort)" = "Tlemsani et al (George cohort), 2021",
  "Tlemsani 2021" = "Tlemsani et al, 2021",
  "Li 2020" = "Li et al, 2020",
  "Mukherjee 2022" = "Mukherjee et al, 2022",
  "Ul Haq 2024" = "Ul Haq et al, 2024",
  "Liu 2022" = "Liu et al, 2022",
  "Peng 2022" = "Peng et al, 2022",
  "Tian 2020" = "Tian et al, 2020",
  "Wu 2023" = "Wu et al, 2023",
  "Wei 2023" = "Wei et al, 2023",
  "Summary" = expression(bold("Summary (" * italic("I")^"2"~"= 94%)"))
)


p <- ggplot(forest_data, aes(x = Effect, y = fct_reorder(Study, Effect))) +
  
  # Summary effect value
  geom_vline(xintercept = summary_data$Effect, color="#636566", linetype=5) +
  
  geom_linerange(aes(xmin=CI_Lower, xmax=CI_Upper), color="#87898a") +
  geom_point(shape=22, size=3, color="#0c94b3", stroke=0.8, fill="#FFF") +
  labs(x = "Proportion of germline variants in cohort", y = "") +
  
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
  scale_x_continuous(trans='log10', labels = scales::percent, minor_breaks = NULL) +
  coord_fixed(ratio = 0.2) +  # Adjust the aspect ratio (increasing the ratio makes the plot wider relative to its height) +
  scale_y_discrete(labels = custom_labels) # Apply custom labels
  


# Disable clipping
ggplot_gtable <- ggplotGrob(p)
ggplot_gtable$layout$clip[ggplot_gtable$layout$name == "panel"] <- "off"

pdf("Prettier figures/germline_prevalence_by_cohort.pdf", width = 9, height = 5)
# Draw the plot
grid.draw(ggplot_gtable)
dev.off()

"#89023E"
"#EA638C"


"#519E8A"
"#476A6F"

"#907F9F"
"#8D5A97"

"#B392AC"
"#735D78"

"#DE9151"
"#BC5D2E"




# Extract data from meta_result
funnel_data <- data.frame(
  Effect = meta_result$TE,         # Effect sizes
  SE = meta_result$seTE,           # Standard errors
  Study = meta_result$studlab
)

pooled_effect <- meta_result$TE.random

# Compute confidence limits
funnel_data <- funnel_data %>%
  mutate(
    Lower95 = 0 - qnorm(0.975) * (1 / SE),  # 95% lower bound
    Upper95 = 0 + qnorm(0.975) * (1 / SE)   # 95% upper bound
  )

# Calculate the confidence funnel bounds
funnel_data <- funnel_data %>%
  mutate(
    Lower95 = meta_result$TE.random - qnorm(0.975) * SE, # 95% lower confidence limit
    Upper95 = meta_result$TE.random + qnorm(0.975) * SE  # 95% upper confidence limit
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
ggsave("Prettier figures/germline_prevalence_by_cohort_funnel.pdf")



