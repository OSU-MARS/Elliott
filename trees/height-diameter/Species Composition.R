library(dplyr)
library(ggplot2)
library(readr)
library(cowplot)

# Read your dataset
df <- read_csv("Elliott_timber_cruise_2015-16.csv")  # Replace with your file path

# Filter NA values
filtered_df <- df %>%
  filter(!is.na(totalHt), !is.na(DBH))

# Define species full names
species_lookup <- c(
  "PSME" = "Douglas-fir",
  "ALRU2" = "Red alder",
  "TSHE" = "Western hemlock",
  "ACMA3" = "Bigleaf maple",
  "UMCA" = "California bay",
  "THPL" = "Western redcedar",
  "PISI" = "Sitka spruce",
  "RHPU" = "Cascara buckthorn",
  "ARME" = "Pacific madrone",
  "Other" = "Other"
)

# Summarize counts with full name labels
species_counts <- filtered_df %>%
  mutate(species_category = ifelse(species %in% names(species_lookup), species, "Other")) %>%
  count(species_category, name = "count") %>%
  mutate(
    species_full = species_lookup[species_category],
    species_full = factor(species_full, levels = species_lookup[unique(species_category)])
  )

# MAIN BARPLOT
main_bar <- ggplot(species_counts, aes(x = species_full, y = count)) +
  geom_col(fill = "darkgrey") +
  geom_text(aes(label = count), vjust = -0.5, size = 4)+
  labs(title = "Overall Species Composition", x = "Species", y = "Number of Trees") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=16))

# ZOOMED-IN BARPLOT FOR RARE SPECIES
zoom_species <- species_counts %>%
  filter(species_category %in% c("ARME", "RHPU", "PISI"))

zoom_bar <- ggplot(zoom_species, aes(x = species_full, y = count)) +
  geom_col(fill = "darkorange") +
  geom_text(aes(label = count), vjust = -0.5, size = 3.5)+
  labs(title = "Species of Interest", x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(size=12,angle = 30, hjust = 1),
    plot.title = element_text(size = 12, face = "bold")
  )

# COMBINE PLOTS WITH SMALLER INSET
final_plot <- ggdraw() +
  draw_plot(main_bar, 0, 0, 1, 1) +
  draw_plot(zoom_bar, x = 0.68, y = 0.68, width = 0.28, height = 0.28)

# Show it
print(final_plot)

# Save final plot to PNG
ggsave("species_composition_inset.png", plot = final_plot, width = 10, height = 7, dpi = 300)

