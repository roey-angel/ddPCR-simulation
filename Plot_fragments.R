# Plot genome fragment distribution
library(ggplot2) # ggplot aes geom_point geom_density2d scale_color_discrete labs theme_bw theme element_text guides guide_legend ggsave
library(readr) # read_tsv
library(stringr) # str_replace %>%
library(dplyr) # %>% pull filter
library(ggExtra) # ggMarginal
library(ragg) # scaling option in ggsave

# wd <- "/proj/Methods_dev/ddPCR/SIPSim_example/Fragments/"
wd <- "./Fragments/"
out_dir <- "./Results/"
Fragments <- read_tsv(paste0(wd, "Frags_table.txt"))
Taxa <- read_tsv(paste0(out_dir, "Random_assembly_ACC.txt"), col_names = FALSE)
Fragments$Taxon <- Taxa$X2[match(str_replace(Fragments$taxon_name,
                                             "(^.*?_.*?)_.*", "\\1"),
                                 Taxa$X1)]

# retain no more than 100 random taxa
if (length(unique(Fragments$taxon_name)) > 100) {
  Fragments %>%
    pull(taxon_name) %>%
    unique() %>%
    sample(100) -> subset
  Fragments <- filter(Fragments, taxon_name %in% subset)
}

p <- ggplot(Fragments, aes(fragGC, fragLength, color = Taxon, group = taxon_name)) +
  geom_point(alpha = 0) +
  geom_density2d(alpha = 1/5) +
  scale_color_discrete('Taxon') +
  labs(x = 'Fragment G+C (%)', y = 'Fragment length (bp)') +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.text = element_text(margin = margin(r = 1, unit = 'cm')))  + # space out the legend
  guides(colour = guide_legend(override.aes = list(alpha = 1/2),
         nrow = length(unique(Fragments$Taxon))/3)) +
  # guides(colour = guide_legend(override.aes = list(alpha = 1/2))) +
  theme(legend.position = "bottom") 
p <- ggMarginal(p, groupColour = TRUE, groupFill = TRUE, alpha = 1/5)

# delete a file created (prob.) by ggMarginal
file.remove("Rplots.pdf")

pngfile <- fs::path(out_dir,  "Frag_dist.png")
ggsave(
  pngfile, 
  p, 
  device = agg_png, 
  width = 22, 
  height = 22, 
  units = "cm",
  res = 300,
  scaling = 0.9
)