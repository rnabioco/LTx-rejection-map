library(cowplot)
library(RColorBrewer)
library(tidyverse)
library(Seurat) #v3
library(tidyverse)
library(readxl)
library(here)
library(scbp)
library(presto)
library(qs)
theme_set(theme_cowplot())

palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

tableu_classic_palatte <-
  c("#1f77b4",
    "#aec7e8",
    "#ff7f0e",
    "#ffbb78",
    "#2ca02c",
    "#98df8a",
    "#d62728",
    "#ff9896",
    "#9467bd",
    "#c5b0d5",
    "#8c564b",
    "#c49c94",
    "#e377c2",
    "#f7b6d2",
    "#7f7f7f",
    "#c7c7c7",
    "#bcbd22",
    "#dbdb8d",
    "#17becf",
    "#9edae5")

discrete_palette_default <- c(tableu_classic_palatte,
                              brewer.pal(8, "Dark2"),
                              palette_OkabeIto)

