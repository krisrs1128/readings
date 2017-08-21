#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## We've run some simulations to better understand nonparametric anova (adonis)
## and the situations where its permutation procedure may or may not be valid.
## Here, we'll study some examples in public phyloseq tutorials about how to use
## this method, and investigate whether they may be artificially inflating
## significance. Our examples come from,
##
## http://joey711.github.io/phyloseq-demo/Restroom-Biogeography.html
##
## author: sankaran.kris@gmail.com
## date: 08/16/2017

###############################################################################
## generic setup
###############################################################################
library("phyloseq")
library("tidyverse")
library("ape")
library("vegan")
scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2")
scale_fill_discrete <- function(...)
  scale_fill_brewer(..., palette="Set2")

theme_set(theme_bw())
min_theme <- theme_update(
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 8),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  legend.key = element_blank()
)

###############################################################################
## The restroom biogeography example
## http://joey711.github.io/phyloseq-demo/Restroom-Biogeography.html
###############################################################################

## first get the data (code directly from link)
zipftp = "ftp://ftp.microbio.me/pub/restroom-data/study_1335_split_library_seqs_and_mapping.zip"
zipfile = tempfile("RestroomBiogeography")
download.file(zipftp, zipfile)
zipfile = "study_1335_split_library_seqs_and_mapping.zip"
import_dir <- tempdir()
unzip(zipfile, exdir = import_dir)
biomfile = paste0(import_dir, "/study_1335_closed_reference_otu_table.biom")
biom = import_biom(biomfile) ## <-- this line is breaking...aborting now
sdfile = paste0(import_dir, "/study_1335_mapping_file.txt")
sample_metadata = import_qiime_sample_data(sdfile)
restroom = merge_phyloseq(biom, sample_metadata)

###############################################################################
## A tutorial from genoweb toulouse
## http://genoweb.toulouse.inra.fr/~formation/15_FROGS/5-June2016/FROGS_phyloseq_23062016.pdf
###############################################################################

## download the data
download_base <- function(url_base, file) {
  tmp <- tempdir()
  download.file(
    file.path(url_base, file),
    file.path(tmp, file)
  )
  file.path(tmp, file)
}

url_base <- "http://genoweb.toulouse.inra.fr/~formation/15_FROGS/5-June2016/phyloseq_June_2016/data/chaillou/"
biom_path <- download_base(url_base, "chaillou.biom")
sample_data_path <- download_base(url_base, "sample_data.tsv")
tree_path <- download_base(url_base, "tree.nwk")

## import it and load into a phyloseq object
food <- import_biom(biom_path)
sample_data(food) <- read.csv(
  sample_data_path,
  sep = "\t",
  row.names = 1
)
phy_tree(food) <- read.tree(tree_path)

## perform adonis
dist.uf <- distance(food, method = "unifrac")
metadata <- as(sample_data(food), "data.frame")
adonis(dist.uf ~ EnvType, data = metadata, perm = 9999)
adonis(dist.uf ~ EnvType + Description, data = metadata, perm = 9999)
adonis(dist.uf ~ Description + EnvType, data = metadata, perm = 9999)

x <- food %>%
  get_taxa() %>%
  asinh()
x_order <- hclust(dist(x))$order
sample_order <- hclust(dist(t(x)))$order
x_df <- x %>%
  data.frame() %>%
  rownames_to_column("species") %>%
  as_data_frame() %>%
  gather(sample, value, -species)
x_df$species <- factor(
  x_df$species,
  levels = taxa_names(food)[x_order]
)
x_df$samples <- factor(
  x_df$sample,
  levels = sample_names(food)[sample_order]
)

x_df <- metadata %>%
  rownames_to_column("sample") %>%
  right_join(x_df)

ggplot(x_df) +
  geom_tile(
    aes(x = species, y = sample, alpha = value, fill = EnvType)
  ) +
  scale_alpha(range = c(0.0, 1)) +
  theme(
    axis.text = element_blank(),
    legend.position =  "bottom"
  )

###############################################################################
## An example with the enterotypes data
## https://nioo.knaw.nl/sites/default/files/downloads/Lab%202_%20phyloseq,%20distances,%20Ordination,%20Adonis%20(46).pdf
###############################################################################
data(enterotype)
ent = prune_taxa(!(taxa_names(enterotype) %in% "-1"),enterotype)
ent = prune_taxa(taxa_sums(ent)>0.2,ent)
ent <-transform_sample_counts(ent,function(x)x/sum(x))
ent <-subset_samples(ent,!is.na(Enterotype))
df = data.frame(sample_data(ent))

OTU <- t(get_taxa(ent))
adonis(OTU~Enterotype*SeqTech,data=df)

x <- get_taxa(ent) %>%
  sqrt()
x_order <- hclust(dist(x))$order
sample_order <- hclust(dist(t(x)))$order
x_df <- x %>%
  data.frame() %>%
  rownames_to_column("species") %>%
  as_data_frame() %>%
  gather(sample, value, -species)
x_df$species <- factor(
  x_df$species,
  levels = taxa_names(ent)[x_order]
)
x_df$samples <- factor(
  x_df$sample,
  levels = sample_names(ent)[sample_order]
)

x_df <- sample_data(ent) %>%
  as("data.frame") %>%
  rownames_to_column("sample") %>%
  right_join(x_df)

ggplot(x_df) +
  geom_tile(
    aes(x = species, y = sample, alpha = value, fill = Enterotype)
  ) +
  scale_alpha(range = c(0.0, 1)) +
  theme(
    axis.text = element_blank(),
    legend.position =  "bottom"
  )

ggplot(x_df) +
  geom_tile(
    aes(x = species, y = sample, alpha = value, fill = SeqTech)
  ) +
  scale_alpha(range = c(0.0, 1)) +
  theme(
    axis.text = element_blank(),
    legend.position =  "bottom"
  )
