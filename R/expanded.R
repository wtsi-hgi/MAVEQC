#' transparent color function
#'
#' @name t_col
#' @param col  color name
#' @param rate alpha rate
#' @return transparent color
t_col <- function(col, rate) {
    newcol <- rgb(col2rgb(col)["red",], col2rgb(col)["green",], col2rgb(col)["blue",], as.integer(rate*255), maxColorValue = 255)
    return(newcol)
}

#' not in function
#'
#' @name %nin%
#' @param x X
#' @param y Y
#' @return True or False
`%nin%` <- function(x, y) !(x %in% y)

#' reverse complement
#'
#' @name revcomp
#' @param seq sequence
#' @return string
revcomp <- function(seq) {
    seq <- toupper(seq)
    splits <- strsplit(seq, "")[[1]]
    reversed <- rev(splits)
    seq_rev <- paste(reversed, collapse = "")
    seq_rev_comp <- chartr("ATCG", "TAGC", seq_rev)
    return(seq_rev_comp)
}

#' trim adaptor sequences
#'
#' @name trim_adaptor
#' @param seq    sequence
#' @param adapt5 5 prime adaptor sequence
#' @param adapt3 3 prime adaptor sequence
#' @return string
trim_adaptor <- function(seq, adapt5, adapt3) {
    adapt5_pos <- regexpr(adapt5, seq, fixed = TRUE)[1]
    adapt3_pos <- regexpr(adapt3, seq, fixed = TRUE)[1]

    is_revcomp <- FALSE
    # ? could adaptor revcomp ?
    if (adapt5_pos < 0 & adapt3_pos < 0) {
        adapt5_revcomp <- revcomp(adapt5)
        adapt3_revcomp <- revcomp(adapt3)

        adapt5_revcomp_pos <- regexpr(adapt5_revcomp, seq, fixed = TRUE)[1]
        adapt3_revcomp_pos <- regexpr(adapt3_revcomp, seq, fixed = TRUE)[1]

        if (adapt3_revcomp_pos < 0 & adapt5_revcomp_pos < 0) {
            return(seq)
        } else {
            is_revcomp <- TRUE
        }
    }

    if (is_revcomp == FALSE) {
        if (adapt3_pos > adapt5_pos) {
            if (adapt5_pos > 0 & adapt3_pos > 0) {
                return(substr(seq, adapt5_pos + nchar(adapt5), adapt3_pos - 1))
            } else if (adapt5_pos > 0 & adapt3_pos < 0) {
                return(substr(seq, adapt5_pos + nchar(adapt5), nchar(seq)))
            } else if (adapt5_pos < 0 & adapt3_pos > 0) {
                return(substr(seq, 1, adapt3_pos - 1))
            }
        } else {
            stop(paste0("====> Error: 3 prime adaptor found before 5 prime adaptor in the sequence: ", seq))
        }
    } else {
        if (adapt5_revcomp_pos > adapt3_revcomp_pos) {
            if (adapt3_revcomp_pos > 0 & adapt5_revcomp_pos > 0) {
                return(substr(seq, adapt3_revcomp_pos + nchar(adapt3_revcomp), adapt5_revcomp_pos - 1))
            } else if (adapt3_revcomp_pos > 0 & adapt5_revcomp_pos < 0) {
                return(substr(seq, adapt3_revcomp_pos + nchar(adapt3_revcomp), nchar(seq)))
            } else if (adapt3_revcomp_pos < 0 & adapt5_revcomp_pos > 0) {
                return(substr(seq, 1, adapt3_revcomp_pos - 1))
            }
        } else {
            stop(paste0("====> Error: 5 prime adaptor (RC) found before 3 prime adaptor (RC) in the sequence: ", seq))
        }
    }
}

#' reverse complement
#'
#' @name cbind_fill
#' @return matrix
cbind_fill <- function(...) {
    nm <- list(...)
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow))
    do.call(cbind, lapply(nm, function (x) rbind(x, matrix(, n - nrow(x), ncol(x)))))
}

#' calculate gini coefficiency for a sample
#'
#' @param x a vector
#' @return a value
cal_gini <- function(x, corr = FALSE, na.rm = TRUE) {
    if (!na.rm && any(is.na(x))) return(NA_real_)
    x <- as.numeric(na.omit(x))
    n <- length(x)
    x <- sort(x)
    G <- sum(x * 1L:n)
    G <- 2 * G/sum(x) - (n + 1L)

    if (corr) {
        return(G / (n - 1L))
    } else {
        return(G / n)
    }
}

#' merge a list of vector values into a data frame
#' values may have different names in the list
#'
#' @param objects a list of vector values
#' @return a data frame
merge_list_to_df <- function(list_vals) {
    dt_out <- data.table()

    for (i in 1:length(list_vals)) {
        dt_val <- as.data.table(list_vals[[i]])
        dt_val$seq <- names(list_vals[[i]])

        if (nrow(dt_out) == 0) {
            dt_out <- dt_val[, c(2, 1)]
            colnames(dt_out) <- c("seq", names(list_vals)[i])
        } else {
            cols <- colnames(dt_out)
            dt_out <- merge(dt_out, dt_val, by = "seq", all = TRUE)
            colnames(dt_out) <- c(cols, names(list_vals)[i])
        }
    }

    df_out <- as.data.frame(dt_out)
    rownames(df_out) <- df_out$seq
    df_out <- subset(df_out, select = -seq)

    return(df_out)
}

#' color blind friendly
#'
#' @param col_id a character to select colors
#' @return a vector of colors
select_colorblind <- function(col_id) {
    col8 <- c("#D55E00", "#56B4E9", "#F0E442",
              "#009E73", "#E69F00", "#0072B2",
              "#CC79A7", "#000000")

    col12 <- c("#88CCEE", "#CC6677", "#DDCC77",
               "#117733", "#332288", "#AA4499",
               "#44AA99", "#999933", "#882255",
               "#661100", "#6699CC", "#888888")

    col21 <- c("#560133", "#EF0096", "#000000",
               "#65019F", "#DA00FD", "#FF92FD",
               "#F60239", "#FF6E3A", "#FFDC3D",
               "#005745", "#00AF8E", "#00EBC1",
               "#00489E", "#0079FA", "#00E5F8",
               "#005A01", "#009503", "#AFFF2A",
               "#00F407", "#9900E6", "#009FFA")

    if (col_id == "col8") {
        return(col8)
    } else if (col_id == "col12") {
        return(col12)
    } else if (col_id == "col21") {
        return(col21)
    } else {
        stop(paste0("====> Error: wrong col_id"))
    }
}

#' fetch objects from list by the names or indexes
#'
#' @export
#' @param objects a list of objects
#' @param tags    a vector of names or indexes
#' @return a list of objects
select_objects <- function(objects, tags) {
    if (length(objects) == 0) {
        stop(paste0("====> Error: no object found in the list"))
    }

    if (length(tags) == 0) {
        stop(paste0("====> Error: please provide tags to fetch objects"))
    } else {
        if (class(tags) %in% c("character", "numeric")) {
            if (class(tags) == "numeric") {
                tags <- as.integer(tags)
            }
        } else {
            stop(paste0("====> Error: wrong tag type, must be integer or character"))
        }
    }

    list_select <- list()
    if (class(tags) == "character") {
        for (i in 1:length(tags)) {
            for (j in 1:length(objects)) {
                if (tags[i] == objects[[j]]@sample) {
                    list_select <- append(list_select, objects[[j]])
                    break
                }
            }
        }
    } else {
        for (i in 1:length(tags)) {
            list_select <- append(list_select, objects[[tags[i]]])
        }
    }

    return(list_select)
}
