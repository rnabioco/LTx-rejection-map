library(cowplot)
library(RColorBrewer)
library(tidyverse)
library(Seurat) #v3
library(tidyverse)
library(readxl)
library(here)
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

#' Plot cells in reduced dimensionality 2D space
#'
#' @description Cells can be colored by gene or feature in meta.data dataframe
#'
#' @param seurat_obj object of class Seurat
#' @param feature feature to plot, either gene name or column in seurat_obj@meta.data
#' @param plot_dat supplemental data.frame containing feature to plot.
#' Must have a column named cell that contains matching colnames in seurat_obj@data
#' @param pt_size size of points produced by geom_point
#' @param pt_alpha alpha value for points plotted by geom_point
#' @param label_text if TRUE display feature labels on plot
#' @param label_size size of label text
#' @param label_color color of label text
#' @param .cols vector of colors to use for plot.
#' @param cell_filter character vector of cell names to include in plot
#' @param palette_type color palette type to use (either viridis, brewer, or cloupe)
#' defaults to using cellranger loupe-like colors
#' @param col_pal palette name to use if palette_type is brewer
#' @param max_y maximum feature value to set scale to. Defaults to max of the feature
#' @param legend_title string to supply for title for the legend
#' @param embedding dimensionality reduction to extract from seurat_obj. Can be any
#' dr method present in seurat_obj@dr (e.g. umap, pca, tsne). defaults to tsne
#' @param show_negative By default the legend value for continuous features will be clipped at zero.
#' If false, then the minumum value for the plotted feature will be used.
plot_feature <- function(seurat_obj,
                         feature = NULL,
                         plot_dat = NULL,
                         pt_size = 0.001,
                         pt_alpha = 1,
                         label_text = FALSE,
                         label_size = 6,
                         label_color = "grey",
                         .cols = NULL,
                         cell_filter = NULL,
                         palette_type = "cloupe",
                         col_pal = "Reds",
                         max_y = NULL,
                         legend_title = NULL,
                         embedding = "tsne",
                         show_negative = FALSE){

  mdata <- seurat_obj@meta.data %>% tibble::rownames_to_column("cell")

  if(!embedding %in% names(seurat_obj@reductions)){
    stop(paste0(embedding, " not found in seurat object"))
  }

  embed_dat <- seurat_obj@reductions[[embedding]]@cell.embeddings %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell")

  embed_cols <- colnames(embed_dat)
  xcol <- embed_cols[2]
  ycol <- embed_cols[3]

  embed_dat <- left_join(mdata, embed_dat, by = "cell")

  if (!is.null(cell_filter)){
    embed_dat <- dplyr::filter(embed_dat,
                               cell %in% cell_filter)
  }

  meta_data_col <- feature %in% colnames(embed_dat)

  if (!is.null(feature) & !meta_data_col) {
    feature_dat <- FetchData(seurat_obj, feature) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell")
    embed_dat <- left_join(embed_dat, feature_dat, by = "cell")
  }

  if (!is.null(plot_dat)){
    embed_dat <- left_join(embed_dat, plot_dat, by = "cell")
  }

  color_aes_str <- feature

  color_aes_str_q <- quo(color_aes_str)
  embed_dat <- embed_dat %>% arrange_at(.vars = color_aes_str)

  p <- ggplot(embed_dat,
              aes_string(xcol, ycol)) +
    geom_point(aes_string(color = paste0("`", color_aes_str, "`")),
               size = pt_size,
               alpha = pt_alpha)

  ## discrete or continuous data?
  if (typeof(embed_dat[[feature]]) %in% c(
    "character",
    "logical"
  ) | is.factor(embed_dat[[feature]])) {
    discrete <- TRUE
  } else {
    discrete <- FALSE
  }

  ## increase legend size
  if (discrete) {
    p <- p + guides(colour = guide_legend(override.aes = list(size = 4))) +
      theme(legend.title = element_blank())
  }

  if (label_text) {
    if(discrete) {
      embed_mean_dat <- embed_dat %>%
        group_by_at(vars(one_of(feature))) %>%
        summarize(med_dim_1 = median(get(xcol)),
                  med_dim_2 = median(get(ycol)))

      p <- p +
        geom_text(data = embed_mean_dat,
                  aes_string(x = "med_dim_1",
                             y = "med_dim_2",
                             label = feature),
                  size = label_size,
                  color = label_color)
    } else {
      warning("label_text not compatible with continuous features")
    }
  }

  ## handle legend limit
  if (is.null(max_y) & !discrete) {
    min_value <- ifelse(show_negative, min(embed_dat[[color_aes_str]]), 0L)
    max_y <- c(min_value, max(embed_dat[[color_aes_str]]))
  } else if (discrete & is.null(max_y)){
    max_y <- c(NA, NA)
  }

  # loupe-like colors
  cols <- rev(brewer.pal(11, "RdGy")[c(1:5, 7)])

  #handle legend name
  if(is.null(legend_title)) legend_title <- color_aes_str

  ## handle zero expression
  if (!all(is.na(max_y)) && all(max_y == c(0, 0))){
    p <- p + scale_color_gradient(low = cols[1], high = cols[1], name = legend_title)
    return(p)
  }

  ## handle colors
  if (is.null(.cols) && !discrete){
    if (palette_type == "viridis") {
      p <- p + scale_color_viridis(discrete = F,
                                   direction = -1,
                                   option = col_pal,
                                   limits = max_y, name = legend_title)
    } else if (palette_type == "brewer") {
      p <- p + scale_color_distiller(limits = max_y,
                                     palette = col_pal,
                                     direction = 1, name = legend_title)
    } else if (palette_type == "cloupe") {
      p <- p + scale_color_gradientn(limits = max_y,
                                     colors = cols, name = legend_title)
    }
  } else if (!is.null(.cols) && !discrete){
    p <- p + scale_color_gradientn(limits = max_y,
                                   colors = .cols, name = legend_title)
  } else {

    if(!is.null(.cols)) {
      # use colors provided
      p <- p + scale_color_manual(
        values = .cols,
        name = legend_title
      )
    } else {
      p <- p + scale_color_manual(
        values = discrete_palette_default,
        name = legend_title
      )
    }
  }
  p
}


plot_umap <- function(...){
  plot_feature(..., embedding = "umap")
}

plot_tsne <- function(...){
  plot_feature(..., embedding = "tsne")
}

plot_pca <- function(...){
  plot_feature(..., embedding = "pca")
}

plot_violin <- function(df, .x, .y,
                        .fill = NULL,
                        .size = 0.50,
                        .width = 1,
                        .scale = "width",
                        .alpha = 1,
                        cols = ggplot2::scale_fill_viridis_d(),
                        single_col = NULL,
                        jitter = F,
                        rotate_x_text = TRUE,
                        arrange_by_fill = TRUE){

  if (arrange_by_fill && !is.null(.fill)){
    tmp <- sym(.fill)
    df <- arrange(df, !!tmp)
    df[[.x]] <- factor(df[[.x]], levels = unique(df[[.x]]))
  }

  p <- ggplot(df, aes_string(x = .x, y = .y))

  if (jitter){
    p <- p  + geom_jitter(size = 0.1, alpha = 0.2, color = "black")
  }

  if (!is.null(single_col)){
    p <- p +
      geom_violin(size = .size,
                  scale = .scale,
                  fill = single_col,
                  alpha = .alpha)
  } else {
    p <- p +
      geom_violin(aes_string(fill = .fill),
                  size = .size,
                  scale = .scale,
                  alpha = .alpha) +
      cols
  }

  if(rotate_x_text){
    p <- p + theme(axis.text.x = element_text(angle = 90,
                                              hjust = 1,
                                              vjust = 0.5))
  }
  p <- p + theme(legend.title = element_blank())
  p
}



ilovehue_pal <- c(
  "#e03d6e",
  "#e27c8b",
  "#a64753",
  "#da453e",
  "#db8364",
  "#a54423",
  "#dc652e",
  "#de8c31",
  "#d2a46c",
  "#8f672b",
  "#cea339",
  "#b2b939",
  "#717822",
  "#627037",
  "#a3b46c",
  "#7ba338",
  "#67c042",
  "#3d8829",
  "#35773e",
  "#55c267",
  "#5ca76a",
  "#277257",
  "#5fcea4",
  "#399d82",
  "#40c2d1",
  "#5099cf",
  "#7490df",
  "#615ea5",
  "#716bdf",
  "#c291d6",
  "#984db6",
  "#d558c2",
  "#e17fc0",
  "#995580",
  "#bd3c80"
)
get_distinct_cols <- function(vec, seed = 42) {

  seq_col_pals <- c("Blues", "Greens", "Oranges", "Purples", "Reds", "Greys")
  #seq_cols <- map(seq_col_pals, ~brewer.pal(9, .x) %>% .[1:9] %>% rev(.))

  vec <- sort(vec)
  n_needed <- rle(as.character(vec))$lengths
  n_groups <- length(levels(factor(vec)))

  if(n_groups > 6){
    stop("not enough palettes for ", n_groups, " groups", call. = FALSE)
  }

  seq_col_pals <- seq_col_pals[order(n_needed, decreasing = T)]

  vals <- list()
  for (i in 1:n_groups){
    n <- n_needed[i]
    cols <- suppressWarnings(brewer.pal(n, seq_col_pals[i]))
    if (n < 3){
      cols <- cols[n:3]
    }
    vals[[i]] <- cols
  }
  unlist(vals)
}



set_xlsx_class <- function(df, col, xlsx_class){
  for(i in seq_along(col)){
    class(df[[col[i]]]) <- xlsx_class
  }
  df
}


#' Extract out reduced dimensions and cell metadata to tibble
#'
#' @param obj Seurat Object
#' @param embedding dr slot to extract (defaults to all embeddings (2D))
#'
get_metadata <- function(obj, embedding = NULL) {

  mdata <- as_tibble(obj@meta.data, rownames = "cell")

  if (!is.null(embedding)) {
    if (!embedding %in% names(obj@reductions)) {
      stop(paste0(embedding, " not found in seurat object"), call. = FALSE)
    }

    embed_dat <- obj@reductions[[embedding]]@cell.embeddings %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell")

  } else {
    embed_dat <- map(names(obj@reductions),
                     ~obj@reductions[[.x]]@cell.embeddings[, 1:2]) %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell")

  }

  embed_dat <- left_join(mdata,
                         embed_dat,
                         by = "cell")
  embed_dat
}


#' plot feature across multiple panels split by group
plot_features_split <- function(sobj, feature, group = "orig.ident",
                                embedding = "umap", cols = NULL, add_title = FALSE,
                                plot_grid = FALSE,
                                ...) {

  # get max coordinates
  dim_reduc <- sobj@reductions[[embedding]]@cell.embeddings[, 1:2]
  x_lims <- c(min(dim_reduc[, 1]), max(dim_reduc[, 1]))
  y_lims <- c(min(dim_reduc[, 2]), max(dim_reduc[, 2]))

  groups <- sort(unique(sobj@meta.data[[group]]))

  if(!is.null(cols)) {
    cols <- cols[1:length(groups)]
    plts <- map2(groups, cols, function(x, y) {
      cells <- rownames(sobj@meta.data)[sobj@meta.data[[group]] == x]
      plot_feature(sobj,
                   feature = feature,
                   embedding = embedding,
                   cell_filter = cells,
                   .cols = y,
                   ...) +
        coord_cartesian(xlim = x_lims, y = y_lims)
    })
  } else {
    plts <- map(groups, function(x) {
      cells <- rownames(sobj@meta.data)[sobj@meta.data[[group]] == x]
      plot_feature(sobj,
                   feature = feature,
                   embedding = embedding,
                   cell_filter = cells,
                   ...) +
        coord_cartesian(xlim = x_lims, y = y_lims)
    })
  }

  if(add_title){
    plts <- map2(plts, groups, ~.x + labs(title = .y))
  }

  if(plot_grid){
    n_grps <- length(unique(sobj[[group]]))
    nr <- max(1, ceiling(n_grps / 2))
    nc <- ceiling(n_grps / nr)

    plts <- cowplot::plot_grid(plotlist = plts,
                       nrow = nr,
                       ncol =  nc)
  }
  plts
}


#'  get clonotype information per cell
#'
#' @param vdj_dir path to clonotypes file
#' @param prefix prefix to add to cell_id
#' @export
get_clonotypes <- function(vdj_dir, prefix = "") {
  ctypes <- read_csv(file.path(vdj_dir, "filtered_contig_annotations.csv"))

  ctypes <- mutate(ctypes,
                   cell_id = str_remove(barcode, "-1$"),
                   cell_id = str_c(prefix, cell_id)) %>%
    select(cell_id,
           raw_clonotype_id) %>%
    unique()

  consensus_ctypes <- read_csv(file.path(vdj_dir,
                                         "clonotypes.csv"))

  ctypes <- left_join(ctypes,
                      consensus_ctypes,
                      by = c("raw_clonotype_id" = "clonotype_id"))

  ctypes
}

#'  add clonotype information per cell
#'
#' @param sobj Seurat object
#' @param vdj_dirs path to clonotypes file
#' @param prefixes prefixes to add to cell_id
#' @export
add_clonotypes <- function(sobj, vdj_dirs, prefixes = NULL) {

  if(is.null(prefixes)){
    prefixes <- rep("", length(vdj_dirs))
  } else {
    if(length(prefixes) != length(vdj_dirs)){
      stop("prefixes must be the same length as the # of vdj_dirs")
    }
  }
  out <- map2_dfr(vdj_dirs, prefixes, get_clonotypes)

  cells <- tibble(cell_id = Cells(sobj))

  res <- left_join(cells, out, by = "cell_id") %>%
    as.data.frame() %>%
    tibble::column_to_rownames("cell_id")

  sobj <- AddMetaData(sobj, metadata = res)

  sobj

}




#' Export Seurat object for UCSC cell browser
#' optionally uses Readr instead of base R for writing to disk
#'
#' @param object Seurat object
#' @param dir path to directory where to save exported files. These are:
#' exprMatrix.tsv, tsne.coords.tsv, meta.tsv, markers.tsv and a default cellbrowser.conf
#' @param dataset.name name of the dataset. Defaults to Seurat project name
#' @param reductions vector of reduction names to export
#' @param markers.file path to file with marker genes
#' @param cluster.field name of the metadata field containing cell cluster
#' @param cb.dir path to directory where to create UCSC cellbrowser static
#' website content root, e.g. an index.html, .json files, etc. These files
#' can be copied to any webserver. If this is specified, the cellbrowser
#' package has to be accessible from R via reticulate.
#' @param port on which port to run UCSC cellbrowser webserver after export
#' @param skip.expr.matrix whether to skip exporting expression matrix
#' @param skip.metadata whether to skip exporting metadata
#' @param skip.reductions whether to skip exporting reductions
#' @param ... specifies the metadata fields to export. To supply field with
#' human readable name, pass name as \code{field="name"} parameter.
#'
#' @return This function exports Seurat object as a set of tsv files
#' to \code{dir} directory, copying the \code{markers.file} if it is
#' passed. It also creates the default \code{cellbrowser.conf} in the
#' directory. This directory could be read by \code{cbBuild} to
#' create a static website viewer for the dataset. If \code{cb.dir}
#' parameter is passed, the function runs \code{cbBuild} (if it is
#' installed) to create this static website in \code{cb.dir} directory.
#' If \code{port} parameter is passed, it also runs the webserver for
#' that directory and opens a browser.
#'
#' @author Maximilian Haeussler, Nikolay Markov
#'
#' @importFrom utils browseURL
#' @importFrom reticulate py_module_available import
#' @importFrom tools file_ext
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ExportToCellbrowser(object = pbmc_small, dataset.name = "PBMC", dir = "out")
#' }
#'
ExportToCellbrowserFast <- function(
  object,
  dir,
  dataset.name = Project(object = object),
  reductions = "tsne",
  markers.file = NULL,
  cluster.field = "Cluster",
  cb.dir = NULL,
  port = NULL,
  skip.expr.matrix = FALSE,
  skip.metadata = FALSE,
  skip.reductions = FALSE,
  ...
) {
  vars <- c(...)

  use_readr <- Seurat:::PackageCheck("readr", error = FALSE)

  if (is.null(x = vars)) {
    vars <- c("nCount_RNA", "nFeature_RNA")
    if (length(x = levels(x = Seurat::Idents(object = object))) > 1) {
      vars <- c(vars, cluster.field)
      names(x = vars) <- c("", "", "ident")
    }
  }
  names(x = vars) <- names(x = vars) %||% vars
  names(x = vars) <- sapply(
    X = 1:length(x = vars),
    FUN = function(i) {
      return(ifelse(
        test = nchar(x = names(x = vars)[i]) > 0,
        yes = names(x = vars[i]),
        no = vars[i]
      ))
    }
  )
  if (!is.null(x = port) && is.null(x = cb.dir)) {
    stop("cb.dir parameter is needed when port is set")
  }
  if (!dir.exists(paths = dir)) {
    dir.create(path = dir)
  }
  if (!dir.exists(paths = dir)) {
    stop("Output directory ", dir, " cannot be created or is a file")
  }
  if (dataset.name == "SeuratProject") {
    warning("Using default project name means that you may overwrite project with the same name in the cellbrowser html output folder")
  }
  order <- colnames(x = object)
  enum.fields <- c()
  # Export expression matrix:
  if (!skip.expr.matrix) {
    # Relatively memory inefficient - maybe better to convert to sparse-row and write in a loop, row-by-row?
    df <- as.data.frame(x = as.matrix(x = Seurat::GetAssayData(object = object)))
    df <- data.frame(gene = rownames(x = object), df, check.names = FALSE)
    gzPath <- file.path(dir, "exprMatrix.tsv.gz")

    message("Writing expression matrix to ", gzPath)

    if(!use_readr){
      z <- gzfile(gzPath, "w")
      write.table(x = df, sep = "\t", file = z, quote = FALSE, row.names = FALSE)
      close(con = z)
    } else {
      readr::write_tsv(x = df, path = gzPath)
    }

  }
  # Export cell embeddings
  embeddings.conf <- c()
  for (reduction in reductions) {
    if (!skip.reductions) {
      df <- Seurat::Embeddings(object = object, reduction = reduction)
      if (ncol(x = df) > 2) {
        warning(
          'Embedding ',
          reduction,
          ' has more than 2 coordinates, taking only the first 2'
        )
        df <- df[, 1:2]
      }
      colnames(x = df) <- c("x", "y")
      df <- data.frame(cellId = rownames(x = df), df)
      fname <- file.path(dir, paste0(reduction, '.coords.tsv'))
      message("Writing embeddings to ", fname)
      if(!use_readr){
       write.table(
        x = df[order, ],
        sep = "\t",
        file = fname,
        quote = FALSE,
        row.names = FALSE
       )
      } else {
        readr::write_tsv(x = df[order, ], path = fname)
      }

    }
    conf <- sprintf(
      '{"file": "%s.coords.tsv", "shortLabel": "Seurat %1$s"}',
      reduction
    )
    embeddings.conf <- c(embeddings.conf, conf)
  }
  # Export metadata
  df <- data.frame(row.names = rownames(x = object[[]]))
  df <- Seurat::FetchData(object = object, vars = names(x = vars))
  colnames(x = df) <- vars
  enum.fields <- Filter(
    f = function(name) {!is.numeric(x = df[[name]])},
    x = vars
  )
  if (!skip.metadata) {
    fname <- file.path(dir, "meta.tsv")
    message("Writing meta data to ", fname)
    df <- data.frame(Cell = rownames(x = df), df, check.names = FALSE)

    if(!use_readr){
      write.table(
        x = df[order, ],
        sep = "\t",
        file = fname,
        quote = FALSE,
        row.names = FALSE
      )
    } else {
      readr::write_tsv(x = df[order, ],
                       path = fname)
    }

  }
  # Export markers
  markers.string <- ''
  if (!is.null(x = markers.file)) {
    ext <- tools::file_ext(x = markers.file)
    fname <- paste0("markers.", ext)
    file.copy(from = markers.file, to = file.path(dir, fname))
    markers.string <- sprintf(
      'markers = [{"file": "%s", "shortLabel": "Seurat Cluster Markers"}]',
      fname
    )
  }
  config <- c(
    'name="%s"',
    'shortLabel="%1$s"',
    'exprMatrix="exprMatrix.tsv.gz"',
    '#tags = ["10x", "smartseq2"]',
    'meta="meta.tsv"',
    '# possible values: "gencode-human", "gencode-mouse", "symbol" or "auto"',
    'geneIdType="auto"',
    'clusterField="%s"',
    'labelField="%2$s"',
    'enumFields=%s',
    '%s',
    'coords=%s'
  )
  config <- paste(config, collapse = '\n')
  enum.string <- paste0(
    "[",
    paste(paste0('"', enum.fields, '"'), collapse = ", "),
    "]"
  )
  coords.string <- paste0(
    "[",
    paste(embeddings.conf, collapse = ",\n"),
    "]"
  )
  config <- sprintf(
    config,
    dataset.name,
    cluster.field,
    enum.string,
    markers.string,
    coords.string
  )
  fname <- file.path(dir, "cellbrowser.conf")
  if (file.exists(fname)) {
    message(
      "`cellbrowser.conf` already exists in target directory, refusing to ",
      "overwrite it"
    )
  } else {
    cat(config, file = fname)
  }
  message("Prepared cellbrowser directory ", dir)
  if (!is.null(x = cb.dir)) {
    if (!reticulate::py_module_available(module = "cellbrowser")) {
      stop(
        "The Python package `cellbrowser` is required to prepare and run ",
        "Cellbrowser. Please install it ",
        "on the Unix command line with `sudo pip install cellbrowser` (if root) ",
        "or `pip install cellbrowser --user` (as a non-root user). ",
        "To adapt the Python that is used, you can either set the env. variable RETICULATE_PYTHON ",
        "or do `require(reticulate) and use one of these functions: use_python(), use_virtualenv(), use_condaenv(). ",
        "See https://rstudio.github.io/reticulate/articles/versions.html; ",
        "at the moment, R's reticulate is using this Python: ",
        reticulate::import(module = 'sys')$executable,
        ". "
      )
    }
    if (!is.null(x = port)) {
      port <- as.integer(x = port)
    }
    message("Converting cellbrowser directory to html/json files")
    cb <- reticulate::import(module = "cellbrowser")
    cb$cellbrowser$build(dir, cb.dir)
    if (!is.null(port)) {
      message("Starting http server")
      cb$cellbrowser$stop()
      cb$cellbrowser$serve(cb.dir, port)
      Sys.sleep(time = 0.4)
      utils::browseURL(url = paste0("http://localhost:", port))
    }
  }
}


