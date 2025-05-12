#' scale_override for facet_custom
#'
#' @name scale_override
#' @param which  which grid in the facet
#' @param scale  x or y scale
#' @return a structure
scale_override <- function(which, scale) {
    if (!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
        stop("which must be an integer of length 1")
    }

    if (is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
        stop("scale must be an x or y position scale")
    }

    structure(list(which = which, scale = scale), class = "scale_override")
}

CustomFacetWrap <- ggproto(
    "CustomFacetWrap", FacetWrap,
    init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
        # make the initial x, y scales list
        scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)

        if(is.null(params$scale_overrides)) return(scales)

        max_scale_x <- length(scales$x)
        max_scale_y <- length(scales$y)

        # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
        for (scale_override in params$scale_overrides) {
            which <- scale_override$which
            scale <- scale_override$scale

            if ("x" %in% scale$aesthetics) {
                if (!is.null(scales$x)) {
                    if (which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
                    scales$x[[which]] <- scale$clone()
                }
            } else if ("y" %in% scale$aesthetics) {
                if (!is.null(scales$y)) {
                    if (which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
                    scales$y[[which]] <- scale$clone()
                }
            } else {
                stop("Invalid scale")
            }
        }

        # return scales
        scales
    }
)

#' facet_wrap_custom
#'
#' @name facet_wrap_custom
#' @param scale_overrides  scale_overrides
facet_wrap_custom <- function(..., scale_overrides = NULL) {
    # take advantage of the sanitizing that happens in facet_wrap
    facet_super <- facet_wrap(...)

    # sanitize scale overrides
    if (inherits(scale_overrides, "scale_override")) {
        scale_overrides <- list(scale_overrides)
    } else if (!is.list(scale_overrides) || !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
        stop("scale_overrides must be a scale_override object or a list of scale_override objects")
    }

    facet_super$params$scale_overrides <- scale_overrides
    ggproto(NULL, CustomFacetWrap, shrink = facet_super$shrink, params = facet_super$params)
}

#' Create barplot filler samples
#'
#' Adds filler rows to a data frame to ensure consistent bar counts per facet in a barplot.
#'
#' @name add_filler_samples
#'
#' @param df A data frame containing the data to be plotted.
#' @param bars_per_facet An integer specifying the number of bars per facet.
#'
#' @return A list with:
#' \describe{
#'   \item{df:}{A data frame including the original and any filler samples.}
#'   \item{filler_names:}{A named character vector with blank labels for filler samples.}
#' }

add_filler_samples <- function(df, bars_per_facet) {

    # How many facets in plot
    facet_groups <- ceiling(nrow(df) / bars_per_facet)

    # Number of missing samples
    n_missing <- bars_per_facet * facet_groups - nrow(df)

    # Create df with filler samples (so that that all bars have equal width)
    if (n_missing > 0) {

      # DF with only filler samples
      col_names <- colnames(df)
      filler_samples <- as.data.frame(matrix(NA, nrow = n_missing, ncol = length(col_names)))
      colnames(filler_samples) <- col_names
      filler_samples$samples <- paste("filler", 1:n_missing, sep = "")

      # Combine DFs
      df <- rbind(df, filler_samples)
      rownames(df) <- NULL

      # Blank x-label names for filler samples
      filler_names <- filler_samples$samples
      filler_names <- setNames(rep("", length(filler_names)), filler_names)
    } else {
      filler_names <- setNames(character(0), character(0))
    }

    # Assign samples to facet groups (rows)
    df$facet_group <- rep(1:facet_groups, each = bars_per_facet)[1:nrow(df)]

    # Return list
    list(df = df, filler_names = filler_names)
}
