# install.packages("R6")

### IMPORT LIBRARIES
library(R6)
library(biomformat)
library(tidyverse)
library(microeco)
library(magrittr)
library(ape)
library(seqinr)

# FUNTION FOR MODIFY THE ANGLE OF THE XTEXT
ggplot_xtext_anglesize <- function(xtext_angle, xtext_size, text_color = "black"){
  if(xtext_angle == 0){
    theme(axis.text.x = element_text(colour = text_color, size = xtext_size))
  }else{
    theme(axis.text.x = element_text(angle = xtext_angle, colour = text_color, vjust = 1, hjust = 1, size = xtext_size))
  }
}

# AUXILIAR FUNTION FOR CREATE A NEW AESTHETIC
new_aesthetic <- function (x, env = globalenv()){
  if (rlang::is_quosure(x)) {
    if (!rlang::quo_is_symbolic(x)) {
      x <- rlang::quo_get_expr(x)
    }
    return(x)
  }
  if (rlang::is_symbolic(x)) {
    x <- rlang::new_quosure(x, env = env)
    return(x)
  }
  x
}

# FUNTION FOR CREATE A MICROECO AESTHETIC
aes_meco <- function(x, y, ...){
  mapping <- list(...)
  if (!missing(x)) 
    mapping["x"] <- list(x)
  if (!missing(y)) 
    mapping["y"] <- list(y)
  caller_env <- parent.frame()
  mapping <- lapply(mapping, function(x) {
    if (is.character(x)) {
      x <- rlang::parse_expr(x)
    }
    new_aesthetic(x, env = caller_env)
  })
  structure(mapping, class = "uneval")
}

# FUNCTION FOR EXPAND THE COLOR PALETTES
expand_colors <- function(color_values, output_length){
  if(output_length <= length(color_values)){
    total_colors <- color_values[1:output_length]
  }else{
    message("Input colors are not enough to use. Add more colors automatically via color interpolation ...")
    ceiling_cycle_times <- ceiling(output_length/length(color_values))
    total_cycle_times <- lapply(seq_along(color_values), function(x){
      if((ceiling_cycle_times - 1) * length(color_values) + x <= output_length){
        ceiling_cycle_times
      }else{
        ceiling_cycle_times - 1
      }
    }) %>% unlist
    total_color_list <- lapply(seq_along(color_values), function(x){
      colorRampPalette(c(color_values[x], "white"))(total_cycle_times[x] + 1)
    })
    total_colors <- lapply(seq_len(ceiling_cycle_times), function(x){
      unlist(lapply(total_color_list, function(y){
        if((x + 1) <= length(y)){
          y[x]
        }
      }))
    }) %>% unlist
  }
  total_colors
}

#' @title
#' Create \code{trans_abund_2} object for taxonomic abundance visualization for italic nested legend.
#'
#' @description
#' This class is a wrapper for the taxonomic abundance transformations and visualization (bar plot).
#' The converted data style is the long-format for \code{ggplot2} plot.
#'
#' @export
trans_abund_2 <- R6Class(classname = "trans_abund",
                       public = list(
                         #' @param dataset default NULL; the object of \code{\link{microtable}} class.
                         #' @param taxrank default "Phylum"; taxonomic level, i.e. a column name in \code{tax_table} of the input object.
                         #'   The function extracts the abundance from the \code{taxa_abund} list according to the names in the list. 
                         #'   If the \code{taxa_abund} list is NULL, the function can automatically calculate the relative abundance to generate \code{taxa_abund} list.
                         #' @param show default 0; the relative abundance threshold for filtering the taxa with low abundance.
                         #' @param ntaxa default 10; how many taxa are selected to show. Taxa are ordered by abundance from high to low. 
                         #'   This parameter does not conflict with the parameter \code{show}. Both can be used. \code{ntaxa = NULL} means the parameter will be invalid.
                         #' @param groupmean default NULL; calculate mean abundance for each group. Select a column name in \code{microtable$sample_table}.
                         #' @param group_morestats default FALSE; only available when \code{groupmean} parameter is provided; 
                         #'   Whether output more statistics for each group, including min, max, median and quantile;
                         #'   Thereinto, quantile25 and quantile75 denote 25\% and 75\% quantiles, respectively.
                         #' @param delete_taxonomy_lineage default TRUE; whether delete the taxonomy lineage in front of the target level.
                         #' @param delete_taxonomy_prefix default TRUE; whether delete the prefix of taxonomy, such as "g__".
                         #' @param prefix default NULL; character string; available when \code{delete_taxonomy_prefix = T}; 
                         #'   default NULL represents using the "letter+__", e.g. "k__" for Phylum level;
                         #'   Please provide the customized prefix when it is not standard, otherwise the program can not correctly recognize it.
                         #' @param use_percentage default TRUE; show the abundance percentage.
                         #' @param input_taxaname default NULL; character vector; input taxa names to select some taxa.
                         #' @param high_level default NULL; a taxonomic rank, such as "Phylum", used to add the taxonomic information of higher level.
                         #'   It is necessary for the legend with nested taxonomic levels in the bar plot.
                         #' @param high_level_fix_nsub default NULL; an integer, used to fix the number of selected abundant taxa in each taxon from higher taxonomic level.
                         #'   If the total number under one taxon of higher level is less than the high_level_fix_nsub, the total number will be used.
                         #'   When \code{high_level_fix_nsub} is provided, the taxa number of higher level is calculated as: \code{ceiling(ntaxa/high_level_fix_nsub)}.
                         #'   Note that \code{ntaxa} means either the parameter \code{ntaxa} or the taxonomic number obtained by filtering according to the \code{show} parameter.
                         #' @return \code{data_abund} stored in the object. The column 'all_mean_abund' represents mean relative abundance across all the samples.
                         #'   So the values in one taxon are all same across all the samples.
                         #'   If the sum of column 'Abundance' in one sample is larger than 1, the 'Abundance', 'SD' and 'SE' has been multiplied by 100.
                         #' @examples
                         #' \donttest{
                         #' data(dataset)
                         #' t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10)
                         #' }
                         initialize = function(
    dataset = NULL, 
    taxrank = "Phylum", 
    show = 0, 
    ntaxa = 10, 
    groupmean = NULL,
    group_morestats = FALSE,
    delete_taxonomy_lineage = TRUE,
    delete_taxonomy_prefix = TRUE,
    prefix = NULL,
    use_percentage = TRUE, 
    input_taxaname = NULL,
    high_level = NULL,
    high_level_fix_nsub = NULL
                         ){

                           if(! taxrank %in% names(dataset$taxa_abund)){
                             stop("The input taxrank: ", taxrank, " is not found! Please check it!")
                           }
                           sample_table <- dataset$sample_table
                           if("Sample" %in% colnames(sample_table)){
                             colnames(sample_table)[colnames(sample_table) == "Sample"] <- "Sample_replace"
                           }
                           abund_data <- dataset$taxa_abund[[taxrank]] %>% 
                             rownames_to_column(var = "Taxonomy") %>% 
                             reshape2::melt(id.vars = "Taxonomy") %>% 
                             `colnames<-`(c("Taxonomy", "Sample", "Abundance"))
                           check_nd <- grepl("__$", abund_data$Taxonomy)
                           if(any(check_nd)){
                             abund_data$Taxonomy[check_nd] %<>% paste0(., "unidentified")
                           }
                           if(delete_taxonomy_lineage | delete_taxonomy_prefix){
                             if(delete_taxonomy_prefix){
                               if(is.null(prefix)){
                                 prefix <- ".__"
                               }
                               abund_data$Taxonomy %<>% gsub(prefix, "", .)
                             }
                             if(delete_taxonomy_lineage){
                               abund_data$Taxonomy %<>% gsub(".*\\|", "", .)
                             }
                           }
                           abund_data %<>% dplyr::group_by(!!! syms(c("Taxonomy", "Sample"))) %>% 
                             dplyr::summarise(Abundance = sum(Abundance)) %>%
                             as.data.frame(stringsAsFactors = FALSE)
                           abund_data$Taxonomy %<>% as.character
                           mean_abund <- tapply(abund_data$Abundance, abund_data$Taxonomy, FUN = mean)
                           # add mean abundance for all samples
                           all_mean_abund <- data.frame(Taxonomy = names(mean_abund), all_mean_abund = mean_abund)
                           rownames(all_mean_abund) <- NULL
                           abund_data %<>% {suppressWarnings(dplyr::left_join(., rownames_to_column(sample_table), by = c("Sample" = "rowname")))}
                           if(!is.null(groupmean)){
                             message(groupmean, " column is used to calculate mean abundance ...")
                             abund_data <- microeco:::summarySE_inter(abund_data, measurevar = "Abundance", groupvars = c("Taxonomy", groupmean), more = group_morestats)
                             colnames(abund_data)[colnames(abund_data) == "Mean"] <- "Abundance"
                             colnames(abund_data)[colnames(abund_data) == groupmean] <- "Sample"
                             if(is.factor(sample_table[, groupmean])){
                               abund_data$Sample %<>% factor(., levels = levels(sample_table[, groupmean]))
                             }
                           }
                           abund_data <- dplyr::left_join(abund_data, all_mean_abund, by = c("Taxonomy" = "Taxonomy"))
                           if(!is.null(high_level)){
                             if(length(high_level) > 1){
                               warning("Input high_level has multiple elements! Only select the first one!")
                               high_level <- high_level[1]
                             }
                             message("Add higher taxonomic level into the table ...")
                             if(! high_level %in% colnames(dataset$tax_table)){
                               stop("Provided high_level must be a colname of input dataset$tax_table!")
                             }else{
                               extract_tax_table <- dataset$tax_table[, c(high_level, taxrank)] %>% unique
                               if(!delete_taxonomy_lineage){
                                 stop("The delete_taxonomy_lineage should be TRUE when high_level is provided!")
                               }
                               if(delete_taxonomy_prefix){
                                 extract_tax_table[, taxrank] %<>% gsub(prefix, "", .)
                               }
                               abund_data <- dplyr::left_join(abund_data, extract_tax_table, by = c("Taxonomy" = taxrank))
                             }
                           }
                           use_taxanames <- as.character(rev(names(sort(mean_abund))))
                           if(!is.null(ntaxa)){
                             ntaxa_theshold <- ntaxa_use <- ntaxa
                           }else{
                             ntaxa_theshold <- ntaxa_use <- sum(mean_abund > show)
                           }
                           if(ntaxa_use > sum(mean_abund > show)){
                             ntaxa_use <- sum(mean_abund > show)
                           }
                           use_taxanames %<>% .[!grepl("unidentified|unculture|Incertae.sedis", .)]
                           if(is.null(input_taxaname)){
                             if(is.null(high_level_fix_nsub)){
                               if(length(use_taxanames) > ntaxa_use){
                                 use_taxanames %<>% .[1:ntaxa_use]
                               }
                             }else{
                               high_level_n <- ceiling(ntaxa_use/high_level_fix_nsub)
                               high_level_ordered_taxa <- abund_data[match(names(mean_abund), abund_data$Taxonomy), high_level] %>% unique %>% .[1:high_level_n]
                               use_taxanames <- lapply(high_level_ordered_taxa, function(x){
                                 tmp <- abund_data[abund_data[, high_level] == x, ]
                                 tmp <- names(mean_abund) %>% .[. %in% tmp$Taxonomy]
                                 tmp[1:ifelse(length(tmp) < high_level_fix_nsub, length(tmp), high_level_fix_nsub)]
                               }) %>% unlist
                             }
                           }else{
                             if(!any(input_taxaname %in% use_taxanames)){
                               stop("The input_taxaname does not match to taxa names! Please check the input!")
                             }else{
                               use_taxanames <- input_taxaname[input_taxaname %in% use_taxanames]
                             }
                           }
                           if(!is.null(high_level)){
                             # sort the taxa in high levels according to abundance sum
                             tmp <- abund_data[abund_data$Taxonomy %in% use_taxanames, c(high_level, "Abundance")]
                             tmp <- tapply(tmp$Abundance, tmp[, high_level], FUN = sum)
                             data_taxanames_highlevel <- as.character(names(sort(tmp, decreasing = TRUE)))
                             self$data_taxanames_highlevel <- data_taxanames_highlevel
                           }
                           if(ntaxa_theshold < sum(mean_abund > show) | show == 0){
                             if(use_percentage == T){
                               abund_data$Abundance %<>% {. * 100}
                               if("SE" %in% colnames(abund_data)) abund_data$SE %<>% {. * 100}
                               if("SD" %in% colnames(abund_data)) abund_data$SD %<>% {. * 100}
                               ylabname <- "Relative abundance (%)"
                             }else{
                               ylabname <- "Relative abundance"
                             }
                           }else{
                             ylabname <- paste0("Relative abundance (", taxrank, " > ", show*100, "%)")
                           }
                           self$use_percentage <- use_percentage
                           self$ylabname <- ylabname
                           self$taxrank <- taxrank
                           self$data_abund <- abund_data
                           self$data_taxanames <- use_taxanames
                           self$high_level <- high_level
                           message('The transformed abundance data is stored in object$data_abund ...')
                         },
    #' @description
    #' Bar plot.
    #'
    #' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); colors palette for the bars.
    #' @param bar_full default TRUE; Whether the bar shows all the features (including 'Others'). 
    #'    Default \code{TRUE} means total abundance are summed to 1 or 100 (percentage). \code{FALSE} means 'Others' will not be shown.
    #' @param others_color default "grey90"; the color for "Others" taxa.
    #' @param facet default NULL; a character vector for the facet; group column name of \code{sample_table}, such as, \code{"Group"};
    #'    If multiple facets are needed, please provide ordered names, such as \code{c("Group", "Type")}.
    #'    The latter should have a finer scale than the former one;
    #'    Please adjust the facet orders in the plot by assigning factors in \code{sample_table} before creating \code{trans_abund} object or 
    #'    assigning factors in the \code{data_abund} table of \code{trans_abund} object.
    #'    When multiple facets are used, please first install package \code{ggh4x} using the command \code{install.packages("ggh4x")}.
    #' @param order_x default NULL; vector; used to order the sample names in x axis; must be the samples vector, such as \code{c("S1", "S3", "S2")}.
    #' @param x_axis_name NULL; a character string; a column name of sample_table in dataset; used to show the sample names in x axis.
    #' @param barwidth default NULL; bar width, see \code{width} in \code{geom_bar}.
    #' @param use_alluvium default FALSE; whether add alluvium plot. If \code{TRUE}, please first install \code{ggalluvial} package.
    #' @param clustering default FALSE; whether order samples by the clustering.
    #' @param clustering_plot default FALSE; whether add clustering plot.
    #'     If \code{clustering_plot = TRUE}, \code{clustering} will be also TRUE in any case for the clustering.
    #' @param cluster_plot_width default 0.2, the dendrogram plot width; available when \code{clustering_plot = TRUE}.
    #' @param facet_color default "grey95"; facet background color.
    #' @param strip_text default 11; facet text size.
    #' @param legend_text_italic default FALSE; whether use italic in legend.
    #' @param xtext_angle default 0; number ranging from 0 to 90; used to adjust x axis text angle to reduce text overlap; 
    #' @param xtext_size default 10; x axis text size.
    #' @param xtext_keep default TRUE; whether retain x text.
    #' @param xtitle_keep default TRUE; whether retain x title.
    #' @param ytitle_size default 17; y axis title size.
    #' @param coord_flip default FALSE; whether flip cartesian coordinates so that horizontal becomes vertical, and vertical becomes horizontal.
    #' @param ggnested default FALSE; whether use nested legend. Need \code{ggnested} package to be installed (https://github.com/gmteunisse/ggnested).
    #'   To make it available, please assign \code{high_level} parameter when creating the object.
    #' @param high_level_add_other default FALSE; whether add 'Others' (all the unknown taxa) in each taxon of higher taxonomic level.
    #'   Only available when \code{ggnested = TRUE}.
    #' @param ... Capture unknown parameters.
    #' @return ggplot2 object. 
    #' @examples
    #' \donttest{
    #' t1$plot_bar(facet = "Group", xtext_keep = FALSE)
    #' }
    plot_bar = function(
    color_values = RColorBrewer::brewer.pal(8, "Dark2"),
    bar_full = TRUE,
    others_color = "grey90",
    facet = NULL,
    order_x = NULL,
    x_axis_name = NULL,
    barwidth = NULL,
    use_alluvium = FALSE,
    clustering = FALSE,
    clustering_plot = FALSE,
    cluster_plot_width = 0.2,
    facet_color = "grey95",
    strip_text = 11,
    legend_text_italic = FALSE,
    xtext_angle = 0,
    xtext_size = 10,
    xtext_keep = TRUE,
    xtitle_keep = TRUE,
    ytitle_size = 17,
    coord_flip = FALSE,
    ggnested = FALSE,
    high_level_add_other = FALSE,
    ...
    ){
      all_parameters <- c(as.list(environment()), list(...))
      if("bar_type" %in% names(all_parameters)){
        warning("Parameter bar_type is deprecated! Please use bar_full instead of it!")
        if(all_parameters["bar_type"] == "full"){
          bar_full <- TRUE
        }else{
          bar_full <- FALSE
        }
      }
      plot_data <- self$data_abund
      # try to filter useless columns
      plot_data %<>% .[, ! colnames(.) %in% c("N", "SD", "SE", "Median", "Min", "Max", "quantile25", "quantile75", "all_mean_abund")]
      use_taxanames <- self$data_taxanames
      if(ggnested){
        if(is.null(self$high_level)){
          stop("The high_level is necessary when ggnested = TRUE! Please assign high_level parameter when creating the object!")
        }
        if(high_level_add_other){
          plot_data$Taxonomy[!plot_data$Taxonomy %in% use_taxanames] <- "Others"
          use_taxanames %<>% c(., "Others")
          new_data <- plot_data %>% dplyr::group_by(!!! syms(c(self$high_level, "Taxonomy", "Sample"))) %>% 
            dplyr::summarise(Abundance = sum(Abundance)) %>%
            as.data.frame(stringsAsFactors = FALSE)
          plot_data_merge <- plot_data[, ! colnames(plot_data) %in% c(self$high_level, "Taxonomy", "Abundance"), drop = FALSE] %>% unique
          plot_data <- dplyr::left_join(new_data, plot_data_merge, by = c("Sample" = "Sample"))
        }
        bar_full <- FALSE
      }
      if(bar_full){
        # make sure that taxonomy info are all in selected use_taxanames in case of special data
        if(!all(plot_data$Taxonomy %in% use_taxanames)){
          plot_data$Taxonomy[!plot_data$Taxonomy %in% use_taxanames] <- "Others"
          new_data <- plot_data %>% dplyr::group_by(!!! syms(c("Taxonomy", "Sample"))) %>% 
            dplyr::summarise(Abundance = sum(Abundance)) %>%
            as.data.frame(stringsAsFactors = FALSE)
          plot_data_merge <- plot_data[, ! colnames(plot_data) %in% c("Taxonomy", "Abundance"), drop = FALSE] %>% unique
          plot_data <- dplyr::left_join(new_data, plot_data_merge, by = c("Sample" = "Sample"))
          plot_data$Taxonomy %<>% factor(., levels = rev(c(use_taxanames, "Others")))
        }else{
          plot_data$Taxonomy %<>% factor(., levels = rev(use_taxanames))
        }
      }else{
        if(ggnested){
          plot_data %<>% .[.[, self$high_level] %in% self$data_taxanames_highlevel, ]
          plot_data[, self$high_level] %<>% factor(., levels = self$data_taxanames_highlevel)
        }
        plot_data %<>% {.[.$Taxonomy %in% use_taxanames, ]}
        # two legend ordering types depending on ggnested
        if(ggnested){
          plot_data$Taxonomy %<>% factor(., levels = use_taxanames)
        }else{
          plot_data$Taxonomy %<>% factor(., levels = rev(use_taxanames))
        }
      }
      # order x axis samples
      plot_data <- private$adjust_axis_facet(
        plot_data = plot_data, 
        x_axis_name = x_axis_name, 
        order_x = order_x
      )
      # arrange plot_data--Abundance according to the Taxonomy-group column factor-levels
      plot_data <- plot_data[unlist(lapply(levels(plot_data$Taxonomy), function(x) which(plot_data$Taxonomy == x))),]
      if(!ggnested){
        if(any(grepl("Others", as.character(plot_data$Taxonomy)))){
          bar_colors_use <- expand_colors(color_values, length(unique(plot_data$Taxonomy)) - 1)
          bar_colors_use <- c(bar_colors_use, others_color)
        }else{
          bar_colors_use <- expand_colors(color_values, length(unique(plot_data$Taxonomy)))
        }
      }else{
        # high_level determine the colors
        bar_colors_use <- expand_colors(color_values, length(unique(plot_data[, self$high_level])))
      }
      if(clustering | clustering_plot){
        data_clustering <- reshape2::dcast(plot_data, Sample ~ Taxonomy, value.var = "Abundance", fun.aggregate = sum) %>% 
          `row.names<-`(.[,1]) %>% .[, -1]
        tmp_hclust <- hclust(dist(data_clustering)) 
        order_x_clustering <- tmp_hclust %>% {.$labels[.$order]} %>% as.character
        plot_data$Sample %<>% factor(., levels = order_x_clustering)
      }
      if(use_alluvium){
        p <- ggplot(plot_data, aes(
          x = Sample, y = Abundance, 
          fill = Taxonomy, color = Taxonomy, 
          weight = Abundance, 
          alluvium = Taxonomy, stratum = Taxonomy
        )) +
          ggalluvial::geom_flow(alpha = .4, width = 3/15) +
          ggalluvial::geom_stratum(width = .2) +
          scale_color_manual(values = rev(bar_colors_use))
      }else{
        if(ggnested){
          p <- ggnested::ggnested(plot_data, aes_meco(x = "Sample", y = "Abundance", main_group = self$high_level, sub_group = "Taxonomy"), main_palette = bar_colors_use)  + 
          ggnested::theme_nested(legend.text = element_text(face = "italic"))
          # MODIFICATION TO OBTAIN LEGEND IN ITALIC
        }else{
          p <- ggplot(plot_data, aes_meco(x = "Sample", y = "Abundance", fill = "Taxonomy"))
        }
        if(bar_full){
          if(self$use_percentage == T){
            p <- p + geom_bar(stat = "identity", position = "stack", show.legend = T, width = barwidth)
          }else{
            p <- p + geom_bar(stat = "identity", position = "fill", show.legend = T, width = barwidth)
          }
        }else{
          p <- p + geom_bar(stat = "identity", position = "stack", show.legend = T, width = barwidth)
        }
      }
      if(!ggnested){
        p <- p + scale_fill_manual(values = rev(bar_colors_use)) + 
          ggnested::theme_nested(legend.text = element_text(face = "italic"))
      }
      p <- p + xlab("") + ylab(self$ylabname)
      if(!is.null(facet)){
        if(coord_flip){
          facet_formula <- reformulate(".", paste0(facet, collapse = " + "))
        }else{
          facet_formula <- reformulate(facet, ".")
        }
        if(length(facet) == 1){
          p <- p + facet_grid(facet_formula, scales = "free", space = "free")
        }else{
          p <- p + ggh4x::facet_nested(facet_formula, nest_line = element_line(linetype = 2), scales = "free", space = "free")
        }
        p <- p + theme(strip.background = element_rect(fill = facet_color, color = facet_color), strip.text = element_text(size=strip_text))
        p <- p + scale_y_continuous(expand = c(0, 0.01))
      }else{
        if(bar_full & self$use_percentage == FALSE){
          p <- p + scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
        }else{
          p <- p + scale_y_continuous(expand = c(0, 0))
        }
      }
      p <- p + theme(panel.grid = element_blank(), panel.border = element_blank()) + 
        theme(axis.line.y = element_line(color = "grey60", linetype = "solid", lineend = "square"))
      if(legend_text_italic == T) {
        p <- p + theme(legend.text = element_text(face = 'italic'))
      }
      if(clustering_plot){
        if(! coord_flip){
          message("Rotate the axis automatically to add the clustering plot ...")
          coord_flip <- TRUE
        }
      }
      p <- p + private$ggplot_xtext_type(xtext_angle = xtext_angle, xtext_size = xtext_size, xtext_keep = xtext_keep, coord_flip = coord_flip)
      p <- p + theme(axis.title.y = element_text(size = ytitle_size))
      if(xtitle_keep == F){
        p <- p + theme(axis.title.x = element_blank())
      }
      p <- p + guides(fill = guide_legend(title = self$taxrank))
      if(use_alluvium | ggnested){
        p <- p + guides(color = guide_legend(title = self$taxrank))
      }
      if(coord_flip){
        p <- p + coord_flip()
      }
      if(clustering_plot){
        left_plot <- ggtree::ggtree(tmp_hclust, hang = 0)
        p %<>% aplot::insert_left(left_plot, width = cluster_plot_width)
      }
      p
    }
  ),
    private = list(
      adjust_axis_facet = function(plot_data, x_axis_name, order_x){
        # order x axis samples and facet
        if(!is.null(x_axis_name)){
          colnames(plot_data)[colnames(plot_data) == "Sample"] <- "Sample_rownames_before"
          if(! x_axis_name %in% colnames(plot_data)){
            stop(paste("No", x_axis_name, "found in the column names of sample_table!"))
          }else{
            colnames(plot_data)[colnames(plot_data) == x_axis_name] <- "Sample"
          }
        }
        if(!is.null(order_x)){
          if(length(order_x) == 1){
            stop("This may be wrong. Only one sample used to order the samples!")
          }else{
            plot_data$Sample %<>% factor(., levels = order_x)
          }
        }
        plot_data
      },
      ggplot_xtext_type = function(xtext_angle, xtext_size, xtext_keep = TRUE, coord_flip = FALSE){
        if(coord_flip){
          if(xtext_keep){
            theme(axis.text.y = element_text(colour = "black", size = xtext_size))
          }else{
            theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
          }
        }else{
          if(xtext_keep){
            ggplot_xtext_anglesize(xtext_angle, xtext_size)
          }else{
            theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
          }
        }
      },
      blank_theme = 
        theme_minimal() +
        theme(
          axis.title = element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank(),
          axis.ticks = element_blank(),
          legend.position="right",
          plot.title=element_text(size=14, face="bold")
        )
    ),
    lock_objects = FALSE,
    lock_class = FALSE
)


####
## IMPORT DATA
mt <- readRDS("../Data/metagenome.Rds")

# Calculate metrics
mt$cal_abund()
mt$cal_alphadiv()
mt$cal_betadiv()

# CREATE PLOT 
t1 <- trans_abund_2$new(mt, taxrank = "Species", ntaxa = 20,  high_level = "Genus", delete_taxonomy_prefix = T)
g_nested <- t1$plot_bar(ggnested = T, xtext_angle = 90, facet = c("Group"),
                        color_values = RColorBrewer::brewer.pal(12, "Set3"))

# SAVE FIGURE
ggsave(filename = "../Results/Figure 1b.tiff", plot = g_nested, dpi = 300,  height = 10, width = 20)













