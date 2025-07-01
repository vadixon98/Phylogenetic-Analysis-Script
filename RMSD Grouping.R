# 1. Enter your data
df <- data.frame(
  range = c("77–87","88–112","113–115","116–135","187–252","325–347","495–516","575–595"),
  type  = c("Intracellular domain",
            "Transmembrane helix 1",
            "Extracellular domain",
            "Transmembrane helix 2",
            "Extracellular domain",
            "Transmembrane helix",
            "Transmembrane helix",
            "Intracellular Domain"),
  rmsd  = c(0.901, 0.479, 0.084, 0.343, 0.454, 0.424, 0.369, 0.537)
)

# 2. Collapse into three groups
library(dplyr)
df <- df %>%
  mutate(group = case_when(
    grepl("Transmembrane", type)    ~ "Transmembrane",
    grepl("Extracellular", type)    ~ "Extracellular",
    grepl("Intracellular", type, ignore.case=TRUE) ~ "Intracellular"
  ))

# 3. Draw boxplots by group
library(ggplot2)
ggplot(df, aes(x = group, y = rmsd)) +
  geom_boxplot() +
  xlab("") + ylab("RMSD (Å)") +
  ggtitle("RMSD Distribution by Domain Class") +
  theme_minimal()
