#' initialize function
setGeneric("qcplot_samqc_readlens", function(object, ...) {
  standardGeneric("qcplot_samqc_readlens")
})

#' create the read length plot
#'
#' @export
#' @param object   sampleQC object
#' @param len_bins the bins of length distribution
#' @param plot_dir the output plot directory
setMethod(
    "qcplot_samqc_readlens",
    signature = "sampleQC",
    definition = function(object,
                          len_bins = seq(0, 300, 50),
                          plot_dir = NULL) {
        read_lens <- data.table()
        for (i in 1:length(object@lengths)) {
            tmp_lens <- object@lengths[[i]][, "length", drop = FALSE]
            tmp_lens$sample <- names(object@lengths)[i]
            tmp_lens <- as.data.table(tmp_lens)

            if (nrow(read_lens) == 0) {
                read_lens <- tmp_lens
            } else {
                read_lens <- rbind(read_lens, tmp_lens)
            }
        }

        sample_names <- vector()
        for (s in object@samples) {
            sample_names <- append(sample_names, s@sample)
        }

        read_lens <- as.data.frame(read_lens)
        read_lens$sample <- factor(read_lens$sample, levels = sample_names)

        p1 <- ggplot(read_lens, aes(x = length)) +
                geom_histogram(aes(y = after_stat(width * density)), breaks = len_bins, color = "black", fill = "grey") +
                geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), linetype = "dashed", color = "yellowgreen", linewidth = 0.4) +
                scale_y_continuous(labels = scales::percent) +
                coord_trans(y = "sqrt") +
                labs(x = "Length Distribution", y = "Composition Percentage", title = "Sample QC read lengths") +
                theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                theme(axis.text = element_text(size = 8, face = "bold")) +
                facet_wrap(~sample, scales = "free", dir = "v", ncol = 4)

        p2 <- ggplot(read_lens, aes(x = length)) +
                geom_histogram(aes(y = after_stat(width * density)), breaks = len_bins, color = "black", fill = "grey") +
                geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), linetype = "dashed", color = "yellowgreen", linewidth = 0.4) +
                scale_y_continuous(labels = scales::percent) +
                labs(x = "Length Distribution", y = "Composition Percentage", title = "Sample QC read lengths") +
                theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                theme(axis.text = element_text(size = 8, face = "bold")) +
                facet_wrap(~sample, scales = "free", dir = "v", ncol = 4)

        pheight <- 400 * as.integer((length(sample_names) / 3))

        if (is.null(plot_dir)) {
            ggplotly(p2)
        } else {
            png(paste0(plot_dir, "/", "sample_qc_read_length.png"), width = 1200, height = pheight, res = 200)
            print(p1)
            dev.off()
        }
    }
)

#' initialize function
setGeneric("qcplot_samqc_clusters", function(object, ...) {
  standardGeneric("qcplot_samqc_clusters")
})

#' create the sequence counts and clusters plot
#'
#' @export
#' @param object    sampleQC object
#' @param qc_type   qc type for plot
#' @param plot_dir  the output plot directory
setMethod(
    "qcplot_samqc_clusters",
    signature = "sampleQC",
    definition = function(object,
                          qc_type = c("plasmid", "screen"),
                          plot_dir = NULL) {
        qc_type <- match.arg(qc_type)

        if (qc_type == "screen") {
            seq_clusters <- object@seq_clusters[[1]]
            seq_clusters_1 <- seq_clusters[seq_clusters$cluster == 1, ]
            seq_clusters_2 <- seq_clusters[seq_clusters$cluster == 2, ]
            seq_clusters_new <- rbind(seq_clusters_1, seq_clusters_2)

            select_colors <- select_colorblind("col8")[1:2]
            fill_colors <- sapply(select_colors, function(x) t_col(x, 0.5), USE.NAMES = FALSE)

            p1 <- ggplot(seq_clusters_new, aes(x = 1:dim(seq_clusters_new)[1], y = count_log2, color = factor(cluster))) +
                    geom_point(shape = 21, size = 1, aes(fill = factor(cluster), color = factor(cluster))) +
                    scale_fill_manual(values = fill_colors) +
                    scale_color_manual(values = select_colors) +
                    labs(x = "sequence index", y = "log2(count+1)", title = "Sample QC clusters") +
                    theme(legend.position = "none", panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 8, face = "bold"))

            p2 <- ggplot(seq_clusters_new, aes(x = count_log2, color = factor(cluster))) +
                    geom_density(aes(fill = factor(cluster), color = factor(cluster))) +
                    scale_fill_manual(values = c(t_col("tomato", 0.5), t_col("royalblue", 0.5))) +
                    scale_color_manual(values = c("tomato", "royalblue")) +
                    labs(x = "log2(count+1)", y = "frequency", title = "Sample QC clusters") +
                    theme(legend.position = "none", panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 12, face = "bold"))
        } else {
            seq_clusters <- data.table()
            for (i in 1:length(object@seq_clusters)) {
                tmp_cluster <- object@seq_clusters[[i]][, c("count_log2", "cluster")]
                tmp_cluster$samples <- names(object@seq_clusters)[i]
                tmp_cluster <- as.data.table(tmp_cluster)
                tmp_cluster$samples <- object@samples[[i]]@sample
                tmp_cluster[cluster == 1, group := "low-count cluster"]
                tmp_cluster[cluster == 2, group := "high-count cluster"]

                if (nrow(seq_clusters) == 0) {
                    seq_clusters <- tmp_cluster
                } else {
                    seq_clusters <- rbind(seq_clusters, tmp_cluster)
                }
            }

            p1 <- ggplot(seq_clusters, aes(x = count_log2, color = samples)) +
                    geom_density() +
                    labs(x = "log2(count+1)", y = "frequency", title = "Sample QC clusters") +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 12, face = "bold")) +
                    facet_wrap(~group, scales = "free", dir = "v")

            p2 <- p1
        }

        if (is.null(plot_dir)) {
            ggplotly(p2)
        } else {
            png(paste0(plot_dir, "/", "sample_qc_seq_clusters.png"), width = 1200, height = 1200, res = 200)
            print(p1)
            dev.off()
        }
    }
)

#' initialize function
setGeneric("qcplot_samqc_total", function(object, ...) {
  standardGeneric("qcplot_samqc_total")
})

#' create the stats plot
#'
#' @export
#' @param object   sampleQC object
#' @param plot_dir the output plot directory
setMethod(
    "qcplot_samqc_total",
    signature = "sampleQC",
    definition = function(object,
                          plot_dir = NULL) {
        df_total <- object@stats[, c("excluded_reads", "accepted_reads")]
        df_total$samples <- rownames(df_total)
        dt_total <- reshape2::melt(as.data.table(df_total), id.vars = "samples", variable.name = "types", value.name = "counts")

        dt_total$samples <- factor(dt_total$samples, levels = mixedsort(levels(factor(dt_total$samples))))

        select_colors <- select_colorblind("col8")[1:2]
        fill_colors <- sapply(select_colors, function(x) t_col(x, 0.5), USE.NAMES = FALSE)

        p1 <- ggplot(dt_total,  aes(x = samples, y = counts, fill = types)) +
                geom_bar(stat = "identity") +
                scale_fill_manual(values = fill_colors) +
                scale_color_manual(values = select_colors) +
                labs(x = "samples", y = "counts", title = "Sample QC Stats") +
                scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
                theme(legend.position = "right", legend.title = element_blank()) +
                theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                theme(axis.text = element_text(size = 8, face = "bold")) +
                theme(axis.text.x = element_text(angle = 90))

        pwidth <- 150 * nrow(df_total)

        if (is.null(plot_dir)) {
            ggplotly(p1)
        } else {
            png(paste0(plot_dir, "/", "sample_qc_stats_total.png"), width = pwidth, height = 1200, res = 200)
            print(p1)
            dev.off()
        }
    }
)

#' initialize function
setGeneric("qcplot_samqc_accepted", function(object, ...) {
  standardGeneric("qcplot_samqc_accepted")
})

#' create the stats plot
#'
#' @export
#' @param object   sampleQC object
#' @param plot_dir the output plot directory
setMethod(
    "qcplot_samqc_accepted",
    signature = "sampleQC",
    definition = function(object,
                          plot_dir = NULL) {
        df_accepted <- object@stats[, c("per_unmapped_reads", "per_ref_reads", "per_pam_reads", "per_library_reads")]
        colnames(df_accepted) <- c("unmapped_reads", "ref_reads", "pam_reads", "library_reads")
        df_accepted <- round(df_accepted * 100, 1)
        df_accepted$samples <- rownames(df_accepted)
        dt_filtered <- reshape2::melt(as.data.table(df_accepted), id.vars = "samples", variable.name = "types", value.name = "percent")

        dt_filtered$samples <- factor(dt_filtered$samples, levels = mixedsort(levels(factor(dt_filtered$samples))))

        df_cov <- object@stats[, c("total_reads", "library_reads", "library_cov")]
        colnames(df_cov) <- c("num_total_reads", "num_library_reads", "library_cov")
        df_cov$samples <- rownames(df_cov)
        df_cov$type <- "coverage"

        select_colors <- select_colorblind("col8")[1:4]
        fill_colors <- sapply(select_colors, function(x) t_col(x, 0.5), USE.NAMES = FALSE)

        y_scale <- max(df_cov$library_cov) * 2

        p1 <- ggplot(dt_filtered,  aes(x = samples, y = percent, fill = types)) +
                geom_bar(stat = "identity", position = "fill") +
                geom_line(data = df_cov, aes(x = samples, y = library_cov / y_scale, group = 1), linetype = "dashed", color = "red", inherit.aes = FALSE) +
                geom_point(data = df_cov, aes(x = samples, y = library_cov / y_scale, color = type), shape = 18, size = 3, inherit.aes = FALSE) +
                scale_y_continuous(labels = scales::percent, sec.axis = sec_axis(~. * y_scale, name = "library coverage")) +
                scale_fill_manual(values = fill_colors) +
                scale_color_manual(values = "red") +
                labs(x = "samples", y = "percent", title = "Sample QC Stats") +
                theme(legend.position = "right", legend.title = element_blank()) +
                theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                theme(axis.text = element_text(size = 8, face = "bold")) +
                theme(axis.text.x = element_text(angle = 90)) +
                geom_text(aes(label = percent), position = position_fill(vjust = 0.5), size = 3)

        pwidth <- 150 * nrow(df_accepted)

        if (is.null(plot_dir)) {
            dt_filtered$types <- factor(dt_filtered$types, levels = rev(levels(dt_filtered$types)))

            ay <- list(overlaying = "y",
                       side = "right",
                       title = "Library Coverage")

            mk <- list(size = 12,
                       symbol = "diamond",
                       color = "red")

            plot_ly(data = dt_filtered, x = ~samples, y = ~percent, color = ~types, type = "bar", colors = rev(fill_colors)) %>%
                layout(barmode = "stack") %>%
                add_markers(data = df_cov, x = ~samples, y = ~library_cov, inherit = FALSE, yaxis = "y2", marker = mk, name = "library") %>%
                layout(yaxis2 = ay)
        } else {
            png(paste0(plot_dir, "/", "sample_qc_stats_accepted.png"), width = pwidth, height = 1200, res = 200)
            print(p1)
            dev.off()
        }

        # bubble plot, may be useful, leave it here

        # p2 <- ggplot(df_cov,  aes(x = total_reads, y = library_reads, color = samples)) +
        #         geom_point(alpha = 0.7, aes(size = library_cov)) +
        #         geom_text(size = 2, color = "black", aes(label = library_cov)) +
        #         geom_text(size = 2, color = "black", vjust = -1, aes(label = samples)) +
        #         labs(x = "total reads", y = "library reads", title = "Sample QC Stats") +
        #         scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
        #         scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
        #         scale_size_continuous(range = c(6, 12)) +
        #         theme(legend.position = "right") +
        #         theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
        #         theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
        #         theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
        #         theme(axis.text = element_text(size = 12, face = "bold"))

        # png(paste0(plot_dir, "/", "sample_qc_stats_cov.png"), width = 1200, height = 1200, res = 200)
        # print(p2)
        # dev.off()
    }
)

#' initialize function
setGeneric("qcplot_samqc_gini", function(object, ...) {
  standardGeneric("qcplot_samqc_gini")
})

#' create the gini plot
#'
#' @export
#' @param object   sampleQC object
#' @param plot_dir the output plot directory
setMethod(
    "qcplot_samqc_gini",
    signature = "sampleQC",
    definition = function(object,
                          plot_dir = NULL) {
        sample_names <- character()
        all_gini <- character()
        for (s in object@samples) {
            sample_names <- append(sample_names, s@sample)
            all_gini <- append(all_gini, s@allstats_qc$gini_coeff)
        }
        names(all_gini) <- sample_names

        lib_gini <- object@stats$gini_coeff_before_qc
        names(lib_gini) <- rownames(object@stats)
        qc_gini <- object@stats$gini_coeff_after_qc
        names(qc_gini) <- rownames(object@stats)

        num_samples <- length(sample_names)
        df_gini <- data.frame(matrix(NA, num_samples * 3, 3))
        colnames(df_gini) <- c("gini", "sample", "type")
        df_gini$gini <- c(all_gini, lib_gini, qc_gini)
        df_gini$sample <- c(names(all_gini), names(lib_gini), names(qc_gini))
        df_gini$type <- c(rep("independent", num_samples), rep("dependent", num_samples), rep("after_qc", num_samples))

        df_gini$gini <- as.numeric(df_gini$gini)
        df_gini$sample <- factor(df_gini$sample, levels = sample_names)
        df_gini$type <- factor(df_gini$type, levels = c("independent", "dependent", "after_qc"))

        gg_colors_fill <- c(t_col("tomato", 0.5), t_col("royalblue", 0.5), t_col("yellowgreen", 0.5))
        gg_colors <- c(c("tomato", "royalblue", "yellowgreen"))
        p1 <- ggplot(df_gini,  aes(x = sample, y = gini, fill = type)) +
                geom_bar(position = "dodge", stat = "identity") +
                scale_fill_manual(values = gg_colors_fill) +
                scale_color_manual(values = gg_colors) +
                labs(x = "samples", y = "score", title = "Sample QC Gini Efficiency") +
                theme(legend.position = "right", legend.title = element_blank()) +
                theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                theme(axis.text = element_text(size = 12, face = "bold")) +
                theme(axis.text.x = element_text(angle = 90)) +
                scale_y_continuous(limits = c(0, 1))

        pwidth <- 150 * num_samples

        if (is.null(plot_dir)) {
            ggplotly(p1)
        } else {
            png(paste0(plot_dir, "/", "sample_qc_gini.png"), width = pwidth, height = 1200, res = 200)
            print(p1)
            dev.off()
        }
    }
)

#' initialize function
setGeneric("qcplot_samqc_pos_cov", function(object, ...) {
  standardGeneric("qcplot_samqc_pos_cov")
})

#' create the position plot
#'
#' @export
#' @param object   sampleQC object
#' @param qc_type  plot type, screen or plasmid
#' @param plot_dir the output plot directory
setMethod(
    "qcplot_samqc_pos_cov",
    signature = "sampleQC",
    definition = function(object,
                          qc_type = c("plasmid", "screen"),
                          plot_dir = NULL) {
        if (is.null(plot_dir)) {
            stop(paste0("====> Error: plot_dir is not provided, no output directory."))
        }

        qc_type <- match.arg(qc_type)

        sample_names <- vector()
        libcounts_pos <- data.table()
        for (s in object@samples) {
            sample_names <- append(sample_names, s@sample)

            tmp_counts <- object@library_counts_pos[[s@sample]][, c("seq", "position", "count")]
            tmp_counts <- as.data.table(tmp_counts)
            tmp_counts[, sample := s@sample]

            if (nrow(libcounts_pos) == 0) {
                libcounts_pos <- tmp_counts
            } else {
                libcounts_pos <- rbind(libcounts_pos, tmp_counts)
            }
        }
        libcounts_pos[, log2p1 := log2(count+1)]

        if (qc_type == "plasmid") {
            libcounts_dependent_pos <- data.table()

            for (s in object@samples) {
                tmp_counts <- s@libcounts[, c("name", "count")]
                colnames(tmp_counts) <- c("oligo_name", "count")
                tmp_counts <- as.data.table(tmp_counts)
                tmp_counts[, sample := s@sample]

                tmp_meta <- s@valiant_meta[, c("oligo_name", "mut_position")]
                tmp_meta <- as.data.table(tmp_meta)

                tmp_counts[tmp_meta, position := i.mut_position, on = .(oligo_name)]
                setorder(tmp_counts, cols = "position")

                if (nrow(libcounts_dependent_pos) == 0) {
                libcounts_dependent_pos <- tmp_counts
                } else {
                    libcounts_dependent_pos <- rbind(libcounts_dependent_pos, tmp_counts)
                }
            }

            libcounts_pos <- libcounts_dependent_pos
            libcounts_pos[, log2p1 := log2(count+1)]
        }

        libcounts_pos$sample <- factor(libcounts_pos$sample, levels = mixedsort(levels(factor(libcounts_pos$sample))))

        p1 <- ggplot(libcounts_pos, aes(x = position, y = log2p1)) +
                geom_point(shape = 16, size = 0.5, color = "tomato", alpha = 0.8) +
                geom_smooth(linewidth = 0.5, se = FALSE, color = "royalblue", fill = "lightblue") +
                geom_hline(yintercept = log2(object@cutoffs$seq_low_count+1), linetype = "dashed", color = "springgreen4", linewidth = 0.4) +
                labs(x = "Genomic Coordinate", y = "log2(count+1)", title = "Sample QC position coverage") +
                ylim(0, as.integer(max(libcounts_pos$log2p1)) + 1) +
                theme(legend.position = "none", panel.grid.major = element_blank()) +
                theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                theme(axis.text = element_text(size = 6, face = "bold")) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                facet_wrap(~sample)

        pheight <- 400 * as.integer((length(sample_names) / 3))

        if (is.null(plot_dir)) {
            stop(paste0("====> Error: plot_dir is not provided, no output directory."))
        } else {
            png(paste0(plot_dir, "/", "sample_qc_position_cov.dots.png"), width = 2400, height = pheight, res = 200)
            print(p1)
            dev.off()
        }
    }
)

#' initialize function
setGeneric("qcplot_samqc_pos_anno", function(object, ...) {
  standardGeneric("qcplot_samqc_pos_anno")
})

#' create the position plot
#'
#' @export
#' @param object    sampleQC object
#' @param samples   a vector of sample names
#' @param type      plot type, lof or all
#' @param plot_dir  the output plot directory
setMethod(
    "qcplot_samqc_pos_anno",
    signature = "sampleQC",
    definition = function(object,
                          samples = NULL,
                          type = "lof",
                          plot_dir = NULL) {
        if (is.null(plot_dir)) {
            stop(paste0("====> Error: plot_dir is not provided, no output directory."))
        }

        if (is.null(samples)) {
            stop(paste0("====> Error: please provide samples, a vector."))
        }

        if (type %nin% c("lof", "all")) {
            stop(paste0("====> Error: wrong type, please use lof or all."))
        }

        libcounts_pos <- as.data.frame(object@library_counts_pos_anno)
        libcounts_pos <- libcounts_pos[, c(samples, "position", "consequence")]

        if (type == "lof") {
            libcounts_pos$consequence <- ifelse(libcounts_pos$consequence == "LOF", "LOF", "Others")

            # be careful, df / vec is by row, not column
            libcounts_pos[, samples] <- t(t(libcounts_pos[, samples]) / object@stats[samples, ]$accepted_reads * 100)

            df_libcounts_pos <- reshape2::melt(libcounts_pos, id.vars = c("consequence", "position"), variable.name = "samples", value.name = "counts")
            df_libcounts_pos$samples <- factor(df_libcounts_pos$samples, levels = samples)

            tmp_cutoff <- object@cutoffs$low_abundance_per * 100

            p1 <- ggplot(df_libcounts_pos, aes(x = position, y = counts)) +
                    geom_point(shape = 19, size = 0.5, aes(color = factor(consequence))) +
                    geom_hline(yintercept = tmp_cutoff, linetype = "dashed", color = "springgreen4", linewidth = 0.4) +
                    scale_color_manual(values = c(t_col("red", 1), t_col("royalblue", 0.2)), labels = c("LOF", "Others")) +
                    labs(x = "Genomic Coordinate", y = "Percentage", title = "Sample QC position percentage", color = "Type") +
                    coord_trans(y = "log2") +
                    scale_y_continuous(breaks = c(0.001, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1)) +
                    theme(legend.position = "right", panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 6, face = "bold")) +
                    facet_wrap(~samples, dir = "v")

            pheight <- 300 * length(samples)

            if (is.null(plot_dir)) {
                stop(paste0("====> Error: plot_dir is not provided, no output directory."))
            } else {
                png(paste0(plot_dir, "/", "sample_qc_position_anno.lof_dots.png"), width = 1200, height = pheight, res = 200)
                print(p1)
                dev.off()
            }
        } else {
            libcounts_pos[, samples] <- t(t(libcounts_pos[, samples]) / object@stats[samples, ]$accepted_reads * 100)

            df_libcounts_pos <- reshape2::melt(libcounts_pos, id.vars = c("consequence", "position"), variable.name = "samples", value.name = "counts")
            df_libcounts_pos$samples <- factor(df_libcounts_pos$samples, levels = samples)

            df_libcounts_pos[df_libcounts_pos == 0] <- NA

            num_colors <- length(unique(libcounts_pos$consequence))
            index_colors <- sample(seq(1, length(select_colorblind("col21"))), num_colors)
            select_colors <- select_colorblind("col21")[index_colors]

            freq_cons <- table(libcounts_pos$consequence)
            names(select_colors) <- names(freq_cons)

            freq_cons <- sort(freq_cons, decreasing = TRUE)
            freq_cons <- names(freq_cons)
            rate_cons <- seq(0.2, 0.1 + length(freq_cons)/10, 0.1)
            names(rate_cons) <- freq_cons

            for (i in 1:(length(select_colors) - 1)) {
                select_colors[i] <- t_col(select_colors[i], rate_cons[names(select_colors[i])])
            }
            select_colors <- as.vector(select_colors)

            tmp_cutoff <- object@cutoffs$low_abundance_per * 100

            p1 <- ggplot(df_libcounts_pos, aes(x = position, y = counts)) +
                    geom_point(shape = 19, size = 0.5, aes(color = factor(consequence))) +
                    geom_hline(yintercept = tmp_cutoff, linetype = "dashed", color = "springgreen4", linewidth = 0.4) +
                    scale_color_manual(values = select_colors) +
                    labs(x = "Genomic Coordinate", y = "Percentage", title = "Sample QC position percentage", color = "Type") +
                    coord_trans(y = "log2") +
                    scale_y_continuous(breaks = c(0.001, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1)) +
                    theme(legend.position = "right", panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 6, face = "bold")) +
                    facet_wrap(~samples, dir = "v")

            pheight <- 300 * length(samples)

            if (is.null(plot_dir)) {
                stop(paste0("====> Error: plot_dir is not provided, no output directory."))
            } else {
                png(paste0(plot_dir, "/", "sample_qc_position_anno.all_dots.png"), width = 1200, height = pheight, res = 200)
                print(p1)
                dev.off()
            }
        }
    }
)

#' initialize function
setGeneric("qcplot_expqc_sample_corr", function(object, ...) {
  standardGeneric("qcplot_expqc_sample_corr")
})

#' create the heatmap of samples
#'
#' @export
#' @param object   experimentQC object
#' @param plot_dir the output plot directory
setMethod(
    "qcplot_expqc_sample_corr",
    signature = "experimentQC",
    definition = function(object,
                          plot_dir = NULL) {
        sample_dist <- as.matrix(dist(t(object@deseq_rlog)))
        sample_rlog <- as.matrix(object@deseq_rlog)

        min_rlog <- round(min(sample_rlog))
        max_rlog <- round(max(sample_rlog))

        # heatmap, leave it temporarily

        # pwidth <- 100 * ncol(sample_rlog)
        # png(paste0(plot_dir, "/", "sample_qc_distance_samples.heatmap.png"), width = pwidth, height = 1200, res = 200)
        # lmat <- rbind(c(4, 3), c(2, 1))
        # lhei <- c(3, 8)
        # lwid <- c(3, 8)

        # heatmap.2(sample_rlog,
        #           distfun = function(x) dist(x, method = "euclidean"),
        #           hclustfun = function(x) hclust(x, method = "ward.D2"),
        #           col = colorpanel(100, "royalblue", "ivory", "tomato"),
        #           na.color = "grey",
        #           breaks = seq(min_rlog, max_rlog, length.out = 101),
        #           density.info = "none", trace = "none", dendrogram = "both",
        #           Rowv = TRUE, Colv = TRUE, labRow = FALSE,
        #           cexCol = 0.6, cexRow = 0.6,
        #           key.xlab = "deseq rlog", key.title = "", key.par = list(cex.lab = 1),
        #           margins = c(8, 1),
        #           colsep = 1:ncol(sample_dist),
        #           sepwidth=c(0.01, 0.01),
        #           lmat = lmat, lhei = lhei, lwid = lwid)
        # dev.off()

        sample_corr <- cor(scale(sample_rlog))
        min_corr <- floor(min(sample_corr) * 10) / 10

        p <- ggcorrplot(sample_corr,
                        method = "square",
	                    hc.method = "ward.D2",
                        hc.order = TRUE,
  		                lab = TRUE,
  		                lab_col = "black",
  		                lab_size = 3,
		                p.mat = cor_pmat(sample_corr),
		                sig.level = 0.05,
                        tl.col = "black",
                        tl.cex = 12)
        p1 <- p + scale_fill_gradient2(limit = c(min_corr, 1),
                                       low = "royalblue",
                                       high =  "tomato",
                                       mid = "ivory",
                                       midpoint = (1 + min_corr) / 2,
                                       name = "correlation")

        if (is.null(plot_dir)) {
            ggplotly(p1)
        } else {
            png(paste0(plot_dir, "/", "sample_qc_samples_corr.png"), width = 1200, height = 1200, res = 200)
            corrplot(sample_corr,
                     method = "color",
                     order = "hclust",
                     col = colorpanel(100, "royalblue", "ivory", "tomato"),
                     col.lim = c(min_corr, 1),
                     is.corr = FALSE,
                     addrect = 3,
                     rect.col = "black",
                     rect.lwd = 1.5,
                     addgrid.col = "white",
                     tl.col = "black",
                     tl.cex = 0.75,
                     addCoef.col = "black",
                     number.cex = 0.75)
            dev.off()
        }
    }
)

#' initialize function
setGeneric("qcplot_expqc_sample_pca", function(object, ...) {
  standardGeneric("qcplot_expqc_sample_pca")
})

#' create the pca of samples
#'
#' @export
#' @param object     experimentQC object
#' @param ntop       the number of top variances
#' @param plot_dir   the output plot directory
setMethod(
    "qcplot_expqc_sample_pca",
    signature = "experimentQC",
    definition = function(object,
                          ntop = 500,
                          plot_dir = NULL) {
        if (is.null(plot_dir)) {
            stop(paste0("====> Error: plot_dir is not provided, no output directory."))
        }

        pca_input <- as.matrix(object@deseq_rlog)
        rv <- rowVars(pca_input)
        select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]

        pca <- prcomp(t(pca_input[select, ]), center = TRUE, scale = TRUE)
        percentVar <- pca$sdev^2 / sum(pca$sdev^2)
        percentVar <- round(percentVar, digits = 3) * 100

        pc1_set <- c((min(pca$x[, 1]) - sd(pca$x[, 1])), (max(pca$x[, 1]) + sd(pca$x[, 1])))
        pc2_set <- c((min(pca$x[, 2]) - sd(pca$x[, 2])), (max(pca$x[, 2]) + sd(pca$x[, 2])))
        pc3_set <- c((min(pca$x[, 3]) - sd(pca$x[, 3])), (max(pca$x[, 3]) + sd(pca$x[, 3])))

        ds_coldata <- object@coldata
        # mark conditions
        default_colors <- c("tomato", "royalblue", "yellowgreen", "orange", "pink", "purple", "coral", "cyan")
        select_colors <- default_colors[1:length(levels(ds_coldata$condition))]
        names(select_colors) <- levels(ds_coldata$condition)

        pca_colors <- 1:nrow(ds_coldata)
        for (i in 1:nrow(ds_coldata)) {
            pca_colors[i] <- select_colors[ds_coldata[i, ]$condition]
        }

        pca_bgs <- sapply(pca_colors, function(x) t_col(x, 0.5))

        # mark replicates
        default_pchs <- c(21, 22, 23, 24, 25)
        select_pchs <- default_pchs[1:length(levels(ds_coldata$replicate))]
        names(select_pchs) <- levels(ds_coldata$replicate)

        pca_pchs <- 1:nrow(ds_coldata)
        for (i in 1:nrow(ds_coldata)) {
            pca_pchs[i] <- select_pchs[ds_coldata[i, ]$replicate]
        }

        png(paste0(plot_dir, "/", "sample_qc_pca_samples.png"), width = 1200, height = 1200, res = 200)
        par(mfrow = c(2, 2), mar = c(4, 4, 4, 1))
        plot(pca$x[, 1], pca$x[, 2], xlab = "PC1", ylab = "PC2", pch = pca_pchs, col = pca_colors, bg = pca_bgs, lwd = 1, cex = 2, xlim = pc1_set, ylim = pc2_set, main = "PC1 vs PC2")
        plot(pca$x[, 2], pca$x[, 3], xlab = "PC2", ylab = "PC3", pch = pca_pchs, col = pca_colors, bg = pca_bgs, lwd = 1, cex = 2, xlim = pc2_set, ylim = pc3_set, main = "PC2 vs PC3")
        plot(pca$x[, 1], pca$x[, 3], xlab = "PC1", ylab = "PC3", pch = pca_pchs, col = pca_colors, bg = pca_bgs, lwd = 1, cex = 2, xlim = pc1_set, ylim = pc3_set, main = "PC1 vs PC3")
        b <- barplot(percentVar, col = t_col("royalblue", 0.5), border = "royalblue", ylim = c(0, 105))
        text(b, percentVar + 5, paste0(percentVar, "%"), cex = 0.6)
        legend("topright", legend = levels(ds_coldata$replicate), pch = select_pchs, cex = 1, bty = "n")
        legend("top", legend = levels(ds_coldata$condition), pch = 19, col = select_colors, cex = 1, bty = "n")
        dev.off()
    }
)

#' initialize function
setGeneric("qcplot_expqc_deseq_fc", function(object, ...) {
  standardGeneric("qcplot_expqc_deseq_fc")
})

#' create fold change and consequence plot
#'
#' @export
#' @param object   experimentQC object
#' @param cons     a vector of consequences showed in the figure
#' @param plot_dir the output plot directory
setMethod(
    "qcplot_expqc_deseq_fc",
    signature = "experimentQC",
    definition = function(object,
                          cons = c("Synonymous_Variant", "LOF", "Missense_Variant"),
                          plot_dir = NULL) {
        if (length(plot_dir) == 0) {
            stop(paste0("====> Error: plot_dir is not provided, no output directory."))
        }

        comparisions <- names(object@deseq_res_anno)
        for (i in 1:length(object@deseq_res_anno)) {
            res <- as.data.frame(object@deseq_res_anno[[i]])
            rownames(res) <- res$seq

            res_cons <- res[res$consequence %in% cons, ]

            p1 <- ggplot(res_cons, aes(x = consequence, y = log2FoldChange)) +
                    geom_violinhalf(trim = FALSE, scale = "width", fill = t_col("yellowgreen", 0.5), color = "yellowgreen", position = position_nudge(x = .2, y = 0)) +
                    geom_jitter(width = 0.15, size = 0.75, aes(color = factor(stat))) +
                    scale_color_manual(values = c(t_col("grey", 0.3), t_col("tomato", 0.8), t_col("royalblue", 0.8))) +
                    ylim(-4, 1) +
                    coord_flip() +
                    labs(x = "log2FoldChange", title = comparisions[i], color = "Type") +
                    theme(legend.position = "right", panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 12, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 8, face = "bold"))

            pheight <- 200 * length(cons)

            if (is.null(plot_dir)) {
                stop(paste0("====> Error: plot_dir is not provided, no output directory."))
            } else {
                png(paste0(plot_dir, "/", "sample_qc_deseq_fc.", comparisions[i], ".violin.png"), width = 1500, height = pheight, res = 200)
                print(p1)
                dev.off()
            }

            # volcano plot has pvaj, may be useful

            # res_cons_volcano <- res_cons
            # res_cons_volcano$padj <- -log10(res_cons_volcano$padj)
            # p2 <- ggplot(res_cons_volcano, aes(x = log2FoldChange, y = padj)) +
            #         geom_point(shape = 19, size = 0.5, aes(color = factor(stat))) +
            #         scale_color_manual(values = c(t_col("grey", 0.3), t_col("tomato", 0.8), t_col("royalblue", 0.8))) +
            #         labs(x = "log2FoldChange", y = "-log10(padj)", title = comparisions[i]) +
            #         theme(legend.position = "right", panel.grid.major = element_blank()) +
            #         theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
            #         theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
            #         theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
            #         theme(axis.text = element_text(size = 12, face = "bold")) +
            #         xlim(-2, 2) +
            #         facet_wrap(~consequence, dir = "v")

            # pheight <- 300 * length(cons)
            # png(paste0(plot_dir, "/", "sample_qc_deseq_fc.", comparisions[i], ".volcano.png"), width = 1200, height = pheight, res = 200)
            # print(p2)
            # dev.off()
        }
    }
)


#' initialize function
setGeneric("qcplot_expqc_deseq_fc_pos", function(object, ...) {
  standardGeneric("qcplot_expqc_deseq_fc_pos")
})

#' create fold change and consequence plot
#'
#' @export
#' @param object        experimentQC object
#' @param consequences  a vector of all the consequences in the vep annotation file
#' @param plot_dir      the output plot directory
setMethod(
    "qcplot_expqc_deseq_fc_pos",
    signature = "experimentQC",
    definition = function(object,
                          consequences = c("Synonymous_Variant",
                                           "LOF",
                                           "Missense_Variant",
                                           "Intronic_Variant",
                                           "Inframe_Deletion",
                                           "Splice_Variant",
                                           "Splice_Polypyrimidine_Tract_Variant",
                                           "Others"),
                          plot_dir = NULL) {
        if (length(plot_dir) == 0) {
            stop(paste0("====> Error: plot_dir is not provided, no output directory."))
        }

        comparisions <- names(object@deseq_res_anno)
        colors <- select_colorblind("col21")[1:length(consequences)]
        select_colors <- sapply(colors, function(x) t_col(x, 0.3), USE.NAMES = FALSE)
        fill_colors <- sapply(colors, function(x) t_col(x, 0.8), USE.NAMES = FALSE)

        for (i in 1:length(object@deseq_res_anno)) {
            dt_res <- object@deseq_res_anno[[i]]
            dt_res$consequence <- factor(dt_res$consequence, levels = consequences)

            p1 <- ggplot(dt_res, aes(x = position, y = log2FoldChange)) +
                    geom_point(aes(size = factor(stat), shape = factor(stat), fill = factor(consequence), color = factor(consequence))) +
                    scale_size_manual(values = c(0.5, 2, 2)) +
                    scale_shape_manual(values = c(16, 24, 25)) +
                    scale_color_manual(values = select_colors) +
                    scale_fill_manual(values = fill_colors) +
                    labs(x = "Genomic Coordinate", y = "Log2 Fold Change", title = comparisions[i]) +
                    theme(legend.position = "right", legend.title = element_blank(), panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title = element_text(size = 12, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 12, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 8, face = "bold")) +
                    ylim(-6, 2) + guides(fill = guide_legend(override.aes = list(shape = 21)))

            if (is.null(plot_dir)) {
                stop(paste0("====> Error: plot_dir is not provided, no output directory."))
            } else {
                png(paste0(plot_dir, "/", "sample_qc_deseq_fc.", comparisions[i], ".position.png"), width = 1500, height = 1000, res = 200)
                print(p1)
                dev.off()
            }
        }
    }
)