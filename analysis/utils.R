require(dplyr)
require(ggplot2)
require(tibble)
require(stringr)
require(tidyr)

my.plot_96_profile_absolute <-
  function(mut_matrix,
           colors = NA,
           condensed = FALSE,
           facet_grid_scales = NULL)
  {
    freq <- full_context <- substitution <- context <- NULL
    if (is.na(colors))
    {
      colors = c("#2EBAED",
                 "#000000",
                 "#DE1C14",
                 "#D4D2D2",
                 "#ADCC54",
                 "#F0D0CE")
    }
    if (length(colors) != 6)
    {
      stop("Provide colors vector with length 6", call. = FALSE)
    }
    tb <-
      mut_matrix %>% as.data.frame() %>% tibble::rownames_to_column("full_context") %>%
      dplyr::mutate(
        substitution = stringr::str_replace(full_context,
                                            "\\w\\[(.*)\\]\\w", "\\1"),
        context = gsub("\\[([A-Z])>[A-Z]\\]", "\\1", full_context)
      ) %>% dplyr::select(-full_context) %>%
      tidyr::pivot_longer(c(-substitution,-context),
                          names_to = "sample",
                          values_to = "freq") %>% dplyr::mutate(sample = factor(sample,
                                                                                levels = unique(sample)))
    if (condensed == TRUE)
    {
      width = 1
      spacing = 0
    }
    else
    {
      width = 0.6
      spacing = 0.5
    }
    facet_grid_scales = ifelse(is.null(facet_grid_scales), "free", facet_grid_scales)
    plot = ggplot(data = tb,
                  aes(
                    x = context,
                    y = freq,
                    fill = substitution,
                    width = width
                  )) +
      geom_bar(stat = "identity",
               colour = "black",
               size = 0.2) +
      scale_fill_manual(values = colors) + facet_grid(sample ~ substitution, scales = facet_grid_scales) +
      ylab("Nr of SNVs") + guides(fill = "none") + theme_bw() +
      theme(
        axis.title.y = element_text(size = 12, vjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(
          size = 5,
          angle = 90,
          vjust = 0.5
        ),
        strip.text.x = element_text(size = 9),
        strip.text.y = element_text(size = 9),
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(spacing, "lines")
      )
    return(plot)
  }

my.plot_contribution_absolute <-
  function(contribution,
           mutation_counts,
           manual_color = NA,
           linewidth = 1)
  {
    m_contribution = melt(t(t(contribution) / colSums(contribution) * mutation_counts))
    colnames(m_contribution) = c("Signature", "Sample", "Contribution")
    plot = ggplot(m_contribution,
                  aes(
                    x = factor(Sample),
                    y = Contribution,
                    fill = factor(Signature),
                    order = Sample
                  )) +
      geom_bar(stat = "identity", color = "black", linewidth = linewidth) +
      labs(x = "", y = "Absolute contribution \n (# of sSNVs / cell)") +
      theme_bw() + theme(panel.grid.minor.x = element_blank(),
                         panel.grid.major.x = element_blank()) + theme(panel.grid.minor.y = element_blank(),
                                                                       panel.grid.major.y = element_blank())
    if (is.na(manual_color))
    {
      plot = plot + xlim(levels(factor(m_contribution$Sample))) + scale_fill_discrete(name = "Signature")
    }
    else
    {
      plot = plot + xlim(levels(factor(m_contribution$Sample))) + scale_fill_manual(values = manual_color, name = "Signature")
    }
    return(plot)
  }

my.plot_contribution_absolute_indel <-
  function(contribution,
           mutation_counts,
           manual_color = NA)
  {
    m_contribution = melt(t(t(contribution) / colSums(contribution) * mutation_counts))
    colnames(m_contribution) = c("Signature", "Sample", "Contribution")
    plot = ggplot(m_contribution,
                  aes(
                    x = factor(Sample),
                    y = Contribution,
                    fill = factor(Signature),
                    order = Sample
                  )) +
      geom_bar(stat = "identity", color = "black") + labs(x = "",
                                                          y = "Absolute contribution \n (# of sIndels / cell)") +
      theme_bw() + theme(panel.grid.minor.x = element_blank(),
                         panel.grid.major.x = element_blank()) + theme(panel.grid.minor.y = element_blank(),
                                                                       panel.grid.major.y = element_blank())
    if (is.na(manual_color))
    {
      plot = plot + xlim(levels(factor(m_contribution$Sample))) + scale_fill_discrete(name = "Signature")
    }
    else
    {
      plot = plot + xlim(levels(factor(m_contribution$Sample))) + scale_fill_manual(values = manual_color, name = "Signature")
    }
    return(plot)
  }

my.DoMultiBarHeatmap <- function (object,
                                  color_set,
                                  features = NULL,
                                  cells = NULL,
                                  group.by = "ident",
                                  additional.group.by = NULL,
                                  group.bar = TRUE,
                                  disp.min = -2.5,
                                  disp.max = NULL,
                                  slot = "scale.data",
                                  assay = NULL,
                                  label = TRUE,
                                  size = 5.5,
                                  hjust = 0,
                                  angle = 45,
                                  raster = TRUE,
                                  draw.lines = TRUE,
                                  lines.width = NULL,
                                  group.bar.height = 0.02,
                                  combine = TRUE)
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data",
                                   yes = 2.5,
                                   no = 6)
  possible.features <- rownames(x = GetAssayData(object = object,
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ",
           slot,
           " slot for the ",
           assay,
           " assay.")
    }
    warning(
      "The following features were omitted as they were not found in the ",
      slot,
      " slot for the ",
      assay,
      " assay: ",
      paste(bad.features,
            collapse = ", ")
    )
  }
  data <-
    as.data.frame(x = as.matrix(x = t(x = GetAssayData(
      object = object,
      slot = slot
    )[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object,
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <-
    object[[c(group.by, additional.group.by)]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    group.use <-
      groups.use[, c(i, additional.group.by), drop = FALSE]
    
    for (colname in colnames(group.use)) {
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }
    }
    
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) *
                                                0.0025)
      placeholder.cells <-
        sapply(
          X = 1:(length(x = levels(x = group.use[[i]])) *
                   lines.width),
          FUN = function(x) {
            return(Seurat:::RandomName(length = 20))
          }
        )
      placeholder.groups <-
        data.frame(foo = rep(x = levels(x = group.use[[i]]), times = lines.width))
      placeholder.groups[additional.group.by] = NA
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells
      
      group.levels <- levels(x = group.use[[i]])
      
      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells
      
      group.use <- rbind(group.use, placeholder.groups)
      
      na.data.group <-
        matrix(
          data = NA,
          nrow = length(x = placeholder.cells),
          ncol = ncol(x = data.group),
          dimnames = list(placeholder.cells,
                          colnames(x = data.group))
        )
      data.group <- rbind(data.group, na.data.group)
    }
    
    #group.use = group.use[order(group.use[[i]]), , drop=F]
    group.use <-
      group.use[with(group.use, eval(parse(text = paste(
        'order(', paste(c(i, additional.group.by), collapse = ', '), ')', sep =
          ''
      )))), , drop = F]
    
    plot <-
      Seurat:::SingleRasterMap(
        data = data.group,
        raster = raster,
        disp.min = disp.min,
        disp.max = disp.max,
        feature.order = features,
        cell.order = rownames(x = group.use),
        group.by = group.use[[i]]
      )
    
    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in colnames(group.use2)) {
        if (colname == group.by) {
          colid = colname
        } else {
          colid = colname
        }
        #Create labels
        levels(x = group.use2[[colname]]) <-
          unique(group.use2[[colname]])
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <-
            c(levels(x = group.use2[[colname]]), na.group)
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <-
            c(scales::hue_pal()(length(x = levels(x = group.use2[[colname]])) - 1), "#FFFFFF")
          
          if (colname == group.by) {
            cols[[colname]] <-
              c(scales::hue_pal()(length(x = levels(x = group.use2[[colname]])) - 1), "#FFFFFF")
          } else{
            cols[[colname]] <-
              c(color_set[1], color_set[2], color_set[1], color_set[2])
          }
          
        } else {
          if (colname == group.by) {
            cols[[colname]] <-
              c(scales::hue_pal()(length(x = levels(x = group.use2[[colname]])) - 1), "#FFFFFF")
          } else{
            cols[[colname]] <- c(color_set[1], color_set[2], "#FFFFFF")
          }
        }
        names(x = cols[[colname]]) <-
          levels(x = group.use2[[colname]])
        
        
        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <-
          max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <-
          c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)
        
        plot <- suppressMessages(
          plot +
            annotation_raster(
              raster = t(x = cols[[colname]][group.use2[[colname]]]),
              xmin = -Inf,
              xmax = Inf,
              ymin = y.pos,
              ymax = y.max
            ) +
            annotation_custom(
              grob = grid::textGrob(
                label = colid,
                hjust = 0,
                gp = grid::gpar(cex = 0.75)
              ),
              ymin = mean(c(y.pos, y.max)),
              ymax = mean(c(y.pos, y.max)),
              xmin = Inf,
              xmax = Inf
            ) +
            coord_cartesian(ylim = c(0, y.max), clip = "off")
        )
        
        #temp <- as.data.frame(cols[[colname]][levels(group.use[[colname]])])
        #colnames(temp) <- 'color'
        #temp$x <- temp$y <- 1
        #temp[['name']] <- as.factor(rownames(temp))
        
        #temp <- ggplot(temp, aes(x=x, y=y, fill=name)) + geom_point(shape=21, size=5) + labs(fill=colname) + theme(legend.position = "bottom")
        #legend <- get_legend(temp)
        #multiplot(plot, legend, heights=3,1)
        
        if ((colname != group.by) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <-
            pbuild$layout$panel_params[[1]]$x$break_positions()
          group.use$x <- x.divs
          label.unique <-
            paste(group.use[[colname]], group.use[[group.by]], sep = "+=$")
          label.x.pos <-
            tapply(X = group.use$x,
                   INDEX = label.unique,
                   FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos),
                                    label.x.pos)
          label.x.pos <-
            label.x.pos[!grepl("NA", label.x.pos$group), ]
          label.x.pos$group <-
            label.x.pos$group %>% lapply(function(x)
              gsub("\\+\\=\\$.*", "", x))
          plot <- plot + geom_text(
            stat = "identity",
            data = label.x.pos,
            aes_string(label = "group",
                       x = "label.x.pos"),
            y = y.max + y.max *
              0.03 * 0.5,
            angle = angle,
            hjust = hjust,
            size = size
          )
          
          plot <-
            suppressMessages(plot + coord_cartesian(
              ylim = c(0,
                       y.max + y.max * 0.002 * max(nchar(
                         x = unique(group.use[[colname]])
                       ), na.rm = TRUE) *
                         size),
              clip = "off"
            ))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}


my.plot_indel_contexts <-
  function (counts,
            same_y = FALSE,
            extra_labels = FALSE,
            condensed = FALSE)
  {
    count <- muttype <- muttype_sub <- muttype_total <- sample <- NULL
    counts <-
      counts %>% as.data.frame() %>% tibble::rownames_to_column("muttype_total") %>%
      tidyr::separate(muttype_total, c("muttype", "muttype_sub"),
                      sep = "_(?=[0-9])") %>% dplyr::mutate(muttype = factor(muttype,
                                                                             levels = unique(muttype))) %>% tidyr::gather(key = "sample",
                                                                                                                          value = "count", -muttype, -muttype_sub) %>% dplyr::mutate(sample = factor(sample,
                                                                                                                                                                                                     levels = unique(sample)))
    nr_muts <-
      counts %>% dplyr::group_by(sample) %>% dplyr::summarise(nr_muts = round(sum(count)))
    facet_labs_y <- stringr::str_c(nr_muts$sample)
    names(facet_labs_y) <- nr_muts$sample
    facet_labs_x <- c("1: C", "1: T", "1: C", "1: T", 2, 3,
                      4, "5+", 2, 3, 4, "5+", 2, 3, 4, "5+")
    names(facet_labs_x) <- levels(counts$muttype)
    if (same_y) {
      facet_scale <- "free_x"
    }
    else {
      facet_scale <- "free"
    }
    colors <- c(
      "#FDBE6F",
      "#FF8001",
      "#B0DD8B",
      "#36A12E",
      "#FDCAB5",
      "#FC8A6A",
      "#F14432",
      "#BC141A",
      "#D0E1F2",
      "#94C4DF",
      "#4A98C9",
      "#1764AB",
      "#E2E2EF",
      "#B6B6D8",
      "#8683BD",
      "#61409B"
    )
    if (extra_labels) {
      title <-
        stringr::str_c("Deletion           ",
                       "Insertion          ",
                       "Deletion                                   ",
                       "Insertion                                  ",
                       "Deletion (MH)")
      x_lab <-
        stringr::str_c("Homopolymer length                            ",
                       "Number of repeat units                                                                               ",
                       "Microhomology length")
    }
    else {
      title <- x_lab <- ""
    }
    if (condensed == TRUE) {
      width <- 1
      spacing <- 0
    }
    else {
      width <- 0.6
      spacing <- 0.5
    }
    fig <-
      ggplot(counts,
             aes(
               x = muttype_sub,
               y = count,
               fill = muttype,
               width = width
             )) + geom_bar(stat = "identity",
                           color = "black",
                           size = 0.2) + facet_grid(
                             sample ~
                               muttype,
                             scales = facet_scale,
                             space = "free_x",
                             labeller = labeller(muttype = facet_labs_x,
                                                 sample = facet_labs_y)
                           ) + scale_fill_manual(values = colors) +
      theme_bw() + labs(
        fill = "Mutation type",
        title = title,
        y = "Nr of indels",
        x = x_lab
      ) + theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing.x = unit(spacing,
                               "lines"),
        panel.spacing.y = unit(spacing, "lines")
      )
    return(fig)
  }

my.plot_dbs_contexts <- function (counts,
          same_y = FALSE,
          condensed = FALSE)
{
  count <- REF <- ALT <- muttype_total <- sample <- NULL
  counts <-
    counts %>% as.data.frame() %>% tibble::rownames_to_column("muttype_total") %>%
    tidyr::separate(muttype_total, c("REF", "ALT"), sep = "_") %>%
    dplyr::mutate(REF = factor(REF, levels = BiocGenerics::unique(REF)))
  bases <- c("A", "C", "G", "T")
  bases1 <- bases
  bases_combi <- tidyr::crossing(bases, bases1)
  counts$ALT <-
    factor(counts$ALT,
           levels = stringr::str_c(bases_combi$bases,
                                   bases_combi$bases1))
  counts <-
    counts %>% tidyr::gather(key = "sample", value = "count",-REF,-ALT) %>% dplyr::mutate(sample = factor(sample,
                                                                                                          levels = unique(sample)))
  nr_muts <-
    counts %>% dplyr::group_by(sample) %>% dplyr::summarise(nr_muts = round(sum(count)))
  if (same_y) {
    facet_scale <- "free_x"
  }
  else {
    facet_scale <- "free"
  }
  facet_labs_y <- stringr::str_c(nr_muts$sample)
  names(facet_labs_y) <- nr_muts$sample
  facet_labs_x <- stringr::str_c(levels(counts$REF), ">NN")
  names(facet_labs_x) <- levels(counts$REF)
  if (condensed == TRUE) {
    width <- 1
    spacing <- 0
  }
  else {
    width <- 0.6
    spacing <- 0.5
  }
  colors <- c(
    "#A6CEE3",
    "#1F78B4",
    "#B2DF8A",
    "#33A02C",
    "#FB9A99",
    "#E31A1C",
    "#FDBF6F",
    "#FF7F00",
    "#CAB2D6",
    "#6A3D9A"
  )
  fig <- ggplot(counts, aes(
    x = ALT,
    y = count,
    fill = REF,
    width = width
  )) + geom_bar(stat = "identity",
                color = "black",
                size = 0.2) + facet_grid(
                  sample ~
                    REF,
                  scales = facet_scale,
                  space = "free_x",
                  labeller = labeller(REF = facet_labs_x,
                                      sample = facet_labs_y)
                ) + scale_fill_manual(guide = "none",
                                      values = colors) + labs(
                                        fill = "Mutation type",
                                        title = "",
                                        y = "Nr of DBSs",
                                        x = ""
                                      ) + theme_bw() + theme(
                                        panel.grid.major.x = element_blank(),
                                        panel.grid.minor.y = element_blank(),
                                        panel.spacing.x = unit(spacing,
                                                               "lines"),
                                        axis.text.x = element_text(
                                          angle = 90,
                                          hjust = 1,
                                          vjust = 0.5
                                        )
                                      )
  return(fig)
}

fit_signatures_until_stable <- function(
    mut_matrix,
    signatures,
    max_delta = 0.004,
    min_signatures = 1
) {
  # ------------------------------------------------------------
  # Input validation
  # ------------------------------------------------------------
  mut_matrix <- as.matrix(mut_matrix)
  signatures <- as.matrix(signatures)
  
  storage.mode(mut_matrix) <- "numeric"
  storage.mode(signatures) <- "numeric"
  
  if (nrow(mut_matrix) == 0 || ncol(mut_matrix) == 0) {
    stop("mut_matrix must contain at least one row and one sample.")
  }
  
  if (nrow(signatures) == 0 || ncol(signatures) == 0) {
    stop("signatures must contain at least one row and one signature.")
  }
  
  if (!is.numeric(max_delta) ||
      length(max_delta) != 1 ||
      is.na(max_delta) ||
      max_delta < 0) {
    stop("max_delta must be a single non-negative number.")
  }
  
  if (!is.numeric(min_signatures) ||
      length(min_signatures) != 1 ||
      is.na(min_signatures) ||
      min_signatures < 1 ||
      min_signatures != as.integer(min_signatures)) {
    stop("min_signatures must be a positive integer.")
  }
  
  min_signatures <- as.integer(min_signatures)
  
  if (min_signatures > ncol(signatures)) {
    stop(
      "min_signatures cannot be greater than the number of input signatures."
    )
  }
  
  if (any(!is.finite(mut_matrix))) {
    stop("mut_matrix contains NA, NaN, or infinite values.")
  }
  
  if (any(!is.finite(signatures))) {
    stop("signatures contains NA, NaN, or infinite values.")
  }
  
  if (any(mut_matrix < 0)) {
    stop("mut_matrix contains negative values.")
  }
  
  if (any(signatures < 0)) {
    stop("signatures contains negative values.")
  }
  
  # ------------------------------------------------------------
  # Match mutation channels
  # ------------------------------------------------------------
  if (!is.null(rownames(mut_matrix)) &&
      !is.null(rownames(signatures))) {
    
    if (!setequal(rownames(mut_matrix), rownames(signatures))) {
      stop(
        paste0(
          "The row names of mut_matrix and signatures do not match. ",
          "Both matrices must contain the same mutation channels."
        )
      )
    }
    
    signatures <- signatures[rownames(mut_matrix), , drop = FALSE]
    
  } else if (nrow(mut_matrix) != nrow(signatures)) {
    stop(
      paste0(
        "mut_matrix and signatures must have the same number of rows ",
        "when mutation-channel row names are unavailable."
      )
    )
  }
  
  if (is.null(colnames(mut_matrix))) {
    colnames(mut_matrix) <- paste0(
      "sample_",
      seq_len(ncol(mut_matrix))
    )
  }
  
  if (anyDuplicated(colnames(mut_matrix))) {
    stop("Sample names in mut_matrix must be unique.")
  }
  
  if (is.null(colnames(signatures))) {
    colnames(signatures) <- paste0(
      "signature_",
      seq_len(ncol(signatures))
    )
  }
  
  if (anyDuplicated(colnames(signatures))) {
    stop("Signature names must be unique.")
  }
  
  # ------------------------------------------------------------
  # Check that fit_to_signatures() is available
  # ------------------------------------------------------------
  if (!exists("fit_to_signatures", mode = "function")) {
    stop(
      paste0(
        "fit_to_signatures() was not found. ",
        "Load the package that provides it, for example ",
        "library(MutationalPatterns)."
      )
    )
  }
  
  # ------------------------------------------------------------
  # Helper functions
  # ------------------------------------------------------------
  cosine_sim <- function(x, y) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    
    if (length(x) != length(y)) {
      stop(
        "Observed and reconstructed mutation profiles have different lengths."
      )
    }
    
    if (any(!is.finite(x)) || any(!is.finite(y))) {
      return(NA_real_)
    }
    
    denominator <- sqrt(sum(x^2)) * sqrt(sum(y^2))
    
    if (!is.finite(denominator) || denominator == 0) {
      return(0)
    }
    
    sum(x * y) / denominator
  }
  
  get_contribution_matrix <- function(
    fit_result,
    expected_signatures = NULL
  ) {
    contribution <- fit_result$contribution
    
    if (is.null(contribution)) {
      stop(
        "fit_to_signatures() returned no contribution matrix."
      )
    }
    
    if (is.null(dim(contribution))) {
      contribution_names <- names(contribution)
      
      contribution <- matrix(
        as.numeric(contribution),
        ncol = 1
      )
      
      if (!is.null(contribution_names)) {
        rownames(contribution) <- contribution_names
      }
    } else {
      contribution <- as.matrix(contribution)
    }
    
    storage.mode(contribution) <- "numeric"
    
    if (ncol(contribution) != 1) {
      stop(
        paste0(
          "Expected one contribution column because one sample is ",
          "fitted at a time, but received ",
          ncol(contribution),
          " columns."
        )
      )
    }
    
    if (is.null(rownames(contribution))) {
      if (!is.null(expected_signatures) &&
          nrow(contribution) == length(expected_signatures)) {
        rownames(contribution) <- expected_signatures
      } else {
        stop(
          "The contribution matrix has no signature row names."
        )
      }
    }
    
    contribution
  }
  
  get_reconstructed_vector <- function(
    fit_result,
    expected_length
  ) {
    reconstructed <- fit_result$reconstructed
    
    if (is.null(reconstructed)) {
      stop(
        "fit_to_signatures() returned no reconstructed profile."
      )
    }
    
    reconstructed <- as.matrix(reconstructed)
    
    if (nrow(reconstructed) == expected_length) {
      if (ncol(reconstructed) != 1) {
        stop(
          paste0(
            "Expected one reconstructed profile because one sample ",
            "is fitted at a time."
          )
        )
      }
      
      reconstructed_vector <- reconstructed[, 1]
      
    } else if (ncol(reconstructed) == expected_length &&
               nrow(reconstructed) == 1) {
      reconstructed_vector <- reconstructed[1, ]
      
    } else {
      stop(
        paste0(
          "The reconstructed profile has unexpected dimensions: ",
          paste(dim(reconstructed), collapse = " x "),
          "."
        )
      )
    }
    
    reconstructed_vector <- as.numeric(reconstructed_vector)
    
    if (length(reconstructed_vector) != expected_length) {
      stop(
        "The reconstructed profile has an unexpected length."
      )
    }
    
    reconstructed_vector
  }
  
  # ------------------------------------------------------------
  # Fit each sample independently
  # ------------------------------------------------------------
  all_results <- vector(
    mode = "list",
    length = ncol(mut_matrix)
  )
  
  names(all_results) <- colnames(mut_matrix)
  
  for (i in seq_len(ncol(mut_matrix))) {
    sample_name <- colnames(mut_matrix)[i]
    mut_sample <- mut_matrix[, i, drop = FALSE]
    
    current_sigs <- signatures
    
    removed_signatures <- character(0)
    removed_for_table <- "None"
    cosine_similarities <- numeric(0)
    
    # Initial fit using all signatures
    fit_res <- fit_to_signatures(
      mut_sample,
      current_sigs
    )
    
    reconstructed_vector <- get_reconstructed_vector(
      fit_res,
      expected_length = nrow(mut_matrix)
    )
    
    current_similarity <- cosine_sim(
      mut_sample[, 1],
      reconstructed_vector
    )
    
    cosine_similarities <- c(
      cosine_similarities,
      current_similarity
    )
    
    repeat {
      if (ncol(current_sigs) <= min_signatures) {
        break
      }
      
      contribution <- get_contribution_matrix(
        fit_res,
        expected_signatures = colnames(current_sigs)
      )
      
      contribution_vector <- contribution[, 1]
      names(contribution_vector) <- rownames(contribution)
      
      # Align contribution values with the current signature matrix
      contribution_vector <- contribution_vector[
        colnames(current_sigs)
      ]
      
      if (length(contribution_vector) == 0) {
        break
      }
      
      if (all(!is.finite(contribution_vector))) {
        break
      }
      
      contribution_for_selection <- contribution_vector
      contribution_for_selection[
        !is.finite(contribution_for_selection)
      ] <- Inf
      
      if (all(is.infinite(contribution_for_selection))) {
        break
      }
      
      weakest_sig <- names(contribution_for_selection)[
        which.min(contribution_for_selection)
      ]
      
      weakest_idx <- match(
        weakest_sig,
        colnames(current_sigs)
      )
      
      if (is.na(weakest_idx)) {
        stop(
          paste0(
            "Could not match signature '",
            weakest_sig,
            "' between the contribution matrix and ",
            "the signature matrix."
          )
        )
      }
      
      # Test the fit after removing the weakest signature
      test_sigs <- current_sigs[
        ,
        -weakest_idx,
        drop = FALSE
      ]
      
      fit_test <- fit_to_signatures(
        mut_sample,
        test_sigs
      )
      
      reconstructed_test <- get_reconstructed_vector(
        fit_test,
        expected_length = nrow(mut_matrix)
      )
      
      new_similarity <- cosine_sim(
        mut_sample[, 1],
        reconstructed_test
      )
      
      delta <- current_similarity - new_similarity
      
      if (!is.finite(delta)) {
        delta <- Inf
      }
      
      # Remove the signature when the loss in cosine similarity
      # is no greater than max_delta.
      if (delta <= max_delta) {
        current_sigs <- test_sigs
        fit_res <- fit_test
        current_similarity <- new_similarity
        
        removed_signatures <- c(
          removed_signatures,
          weakest_sig
        )
        
        removed_for_table <- c(
          removed_for_table,
          weakest_sig
        )
        
        cosine_similarities <- c(
          cosine_similarities,
          new_similarity
        )
      } else {
        break
      }
    }
    
    # fit_res already corresponds to current_sigs, so there is no need
    # to perform the same fit again.
    final_fit <- fit_res
    
    decay_table <- data.frame(
      step = seq_along(cosine_similarities),
      removed_signature = removed_for_table,
      cosine_similarity = cosine_similarities,
      sample = sample_name,
      stringsAsFactors = FALSE
    )
    
    all_results[[i]] <- list(
      sample = sample_name,
      selected_signatures = colnames(current_sigs),
      removed_signatures = removed_signatures,
      decay_table = decay_table,
      fit_res = final_fit
    )
  }
  
  # ------------------------------------------------------------
  # Combine contribution matrices without dplyr/full_join
  # ------------------------------------------------------------
  contribution_list <- lapply(
    all_results,
    function(result) {
      contribution <- get_contribution_matrix(
        result$fit_res,
        expected_signatures = result$selected_signatures
      )
      
      colnames(contribution) <- result$sample
      contribution
    }
  )
  
  all_signature_names <- unique(
    unlist(
      lapply(contribution_list, rownames),
      use.names = FALSE
    )
  )
  
  # Keep the original signature order where possible
  all_signature_names <- c(
    colnames(signatures)[
      colnames(signatures) %in% all_signature_names
    ],
    setdiff(
      all_signature_names,
      colnames(signatures)
    )
  )
  
  contribution <- matrix(
    0,
    nrow = length(all_signature_names),
    ncol = length(all_results),
    dimnames = list(
      all_signature_names,
      names(all_results)
    )
  )
  
  for (j in seq_along(contribution_list)) {
    sample_contribution <- contribution_list[[j]]
    
    contribution[
      rownames(sample_contribution),
      j
    ] <- as.numeric(sample_contribution[, 1])
  }
  
  # ------------------------------------------------------------
  # Combine reconstructed mutation profiles
  # ------------------------------------------------------------
  reconstructed <- vapply(
    all_results,
    function(result) {
      get_reconstructed_vector(
        result$fit_res,
        expected_length = nrow(mut_matrix)
      )
    },
    numeric(nrow(mut_matrix))
  )
  
  if (is.null(dim(reconstructed))) {
    reconstructed <- matrix(
      reconstructed,
      ncol = 1
    )
  }
  
  rownames(reconstructed) <- rownames(mut_matrix)
  colnames(reconstructed) <- names(all_results)
  
  # ------------------------------------------------------------
  # Combine decay tables
  # ------------------------------------------------------------
  decay_table <- do.call(
    rbind,
    lapply(
      all_results,
      function(result) result$decay_table
    )
  )
  
  rownames(decay_table) <- NULL
  
  # ------------------------------------------------------------
  # Return results
  # ------------------------------------------------------------
  list(
    contribution = contribution,
    reconstructed = reconstructed,
    decay_table = decay_table,
    per_sample = all_results
  )
}