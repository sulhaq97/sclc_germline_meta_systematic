library(meta)
library(metafor)


blah <- read.csv("prevalence.by.cohort.csv")

meta_result <- metaprop(
  event = blah[, 5],
  n = blah[, 4],
  studlab = blah[, 2],   # Study name
  sm = "PLOGIT",         # Logit transformation for proportions
  method = "Inverse",    # Inverse variance method
  comb.fixed = F,     # Perform fixed-effects meta-analysis
  comb.random = TRUE,    # Perform random-effects meta-analysis
  hakn = TRUE,           # Use Hartung-Knapp adjustment
  prediction = TRUE      # Include prediction interval
)

pdf("germline_variants_prevalence_by_cohort.pdf", width=10, height=5)
forest(meta_result,
       xlim = c(0, 1),   # Define x-axis limits
       smlab = "Proportion of Events",   # Label for summary measure
       digits = 2,        # Number of decimal places in the plot
       col.diamond = "blue",  # Color for the summary diamond
       col.study = "black",   # Color for study results
       col.predict = "red")   # Color for prediction interval
dev.off()

pdf("funnel_germline_variants_prevalence_by_cohort.pdf", width=5, height=5)
funnel(meta_result, studlab = F, backtransf = F)
dev.off()

metabias(meta_result, method.bias = "linreg", k.min = 2)






genes <- c("ATM", "BRCA1", "BRCA2", "BRIP1", "CHEK1", "CHEK2", "GJB2", 
           "MUTYH", "PARK2", "POLQ", "SLFN11", "TP53")


bygene <- read.csv("prevalence.by.gene.csv")

for(eacho in genes){
  blah <- subset(bygene, bygene$Gene %in% eacho)
  
  meta_result <- metaprop(
    event = blah[, 5],  # Events in that population
    n = blah[, 3],       # Total in cohort
    studlab = blah[, 2],  # Study name
    sm = "PLOGIT",                 # Logit transformation for proportions
    method = "Inverse",            # Inverse variance method
    comb.fixed = TRUE,             # Perform fixed-effects meta-analysis
    comb.random = F,            # Perform random-effects meta-analysis
    hakn = TRUE,                   # Use Hartung-Knapp adjustment
    prediction = TRUE               # Include prediction interval
  )
  
  customname <- paste0("forest plots by gene/forest_plot_", eacho, ".pdf")
  pdf(customname, width = 10, height = 4)
  
  customlabel <- paste0(eacho)
  
  # Generate the forest plot
  forest(meta_result,
         xlim = c(0, 0.2),   # Define x-axis limits
         smlab = customlabel,
         digits = 2,
         col.diamond = "blue",
         col.study = "black",
         col.predict = "red")
  dev.off()
  
  customname2 <- paste0("funnel plots by gene/funnel_plot_", eacho, ".pdf")
  
  pdf(customname2, width=10, height=10)
  funnel(meta_result, studlab = T)
  dev.off()
  
  print(eacho)
  print(metabias(meta_result, method.bias = "linreg", k.min = 2))
  
}


hrs <- read.csv("hrs_by_cohort.csv")
hrs$logHR <- log(hrs$HR_OS)
# Calculate standard error of the log HR
hrs$SE <- (log(hrs$HR_OS_95high) - log(hrs$HR_OS_95low)) / (2 * 1.96)

meta_analysis <- metagen(TE = logHR, seTE = SE, studlab = hrs$Study, 
                         sm = "HR", comb.fixed = TRUE, comb.random = FALSE, 
                         hakn = TRUE, method.tau = "REML", data = hrs)

# Create a forest plot (regular HR is automatically shown)
pdf("forest_plot_hrs_by_cohort.pdf", width = 10 , height = 5)
forest(meta_analysis)
dev.off()

pdf("funnel_hrs_by_cohort.pdf", width=5, height=5)
funnel(meta_analysis, studlab = F)
dev.off()



library(dplyr)
# Combined BRCA1 and BRCA2 prevalence
blah <- subset(bygene, bygene$Gene %in% c("BRCA1", "BRCA2"))
blah$Gene <- factor(blah$Gene)
levels(blah$Gene) <- c("BRCA1/2", "BRCA1/2")

blah2 <- blah %>%
  group_by(Study) %>%  # Group by the "Study" column
  summarize(
    Total_SCLC = first(total_sclc),  # Retain the value from the total_sclc column,
    Total_Events = sum(Event, na.rm = TRUE))  # Summarize the "Event" column

meta_result <- metaprop(
  event = blah2$Total_Events,  # Events in that population
  n = blah2$Total_SCLC,       # Total in cohort
  studlab = blah2$Study,  # Study name
  sm = "PLOGIT",                 # Logit transformation for proportions
  method = "Inverse",            # Inverse variance method
  comb.fixed = TRUE,             # Perform fixed-effects meta-analysis
  comb.random = TRUE,            # Perform random-effects meta-analysis
  hakn = TRUE,                   # Use Hartung-Knapp adjustment
  prediction = TRUE               # Include prediction interval
)

# Generate the forest plot
forest(meta_result,
       xlim = c(0, 1),   # Define x-axis limits
       smlab = "BRCA1/2",
       digits = 2,
       col.diamond = "blue",
       col.study = "black",
       col.predict = "red")

funnel(meta_result, studlab = T)

