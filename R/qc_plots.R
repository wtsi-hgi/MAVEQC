#' initialize function
setGeneric("qcplot_readlens", function(object, ...) {
  standardGeneric("qcplot_readlens")
})

#' create the read length plot
#'
#' @export
#' @param object  sampleQC object
#' @param plotdir the output plot directory
setMethod(
    "qcplot_readlens",
    signature = "sampleQC",
    definition = function(object,
                          plotdir) {
        if (length(plotdir) == 0) {
            stop(paste0("====> Error: plotdir is not provided, no output directory."))
        }

        read_lens <- data.table()
        for (i in 1:length(object@lengths)) {
            tmp_lens <- object@lengths[[i]][, "length", drop = FALSE]
            tmp_lens$samples <- names(object@lengths)[i]
            tmp_lens <- as.data.table(tmp_lens)

            if (nrow(read_lens) == 0) {
                read_lens <- tmp_lens
            } else {
                read_lens <- rbind(read_lens, tmp_lens)
            }
        }

        read_lens <- as.data.frame(read_lens)
        read_lens$samples <- factor(read_lens$samples, levels = colnames(object@effective_counts))

        p1 <- ggplot(read_lens, aes(x = factor(samples), y = length)) +
                geom_violinhalf(alpha = 0.3, scale = "width", color = "royalblue", fill = t_col("royalblue", 0.5)) +
                coord_flip() +
                labs(x = "read length", y = "frequency", title = "Sample QC read lengths") +
                theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                theme(axis.title = element_text(size = 16,face = "bold", family = "Arial")) +
                theme(plot.title = element_text(size = 16,face = "bold.italic", family = "Arial")) +
                theme(axis.text = element_text(size = 10, face = "bold"))

        pheight <- 150 * length(object@lengths)
        png(paste0(plotdir, "/", "sample_qc_read_length.violin.png"), width = 1200, height = pheight, res = 200)
        print(p1)
        dev.off()
    }
)

#' initialize function
setGeneric("qcplot_clusters", function(object, ...) {
  standardGeneric("qcplot_clusters")
})

#' create the sequence counts and clusters plot
#'
#' @export
#' @param object  sampleQC object
#' @param qctype  qc type for plot
#' @param plotdir the output plot directory
setMethod(
    "qcplot_clusters",
    signature = "sampleQC",
    definition = function(object,
                          qctype,
                          plotdir) {
        if (length(qctype) == 0) {
            stop(paste0("====> Error: qctype is not provided, plasmid or screen."))
        } else {
            if (qctype %nin% c("plasmid", "screen")) {
                stop(paste0("====> Error: wrong qctype, plasmid or screen."))
            } else {
                if (qctype == "screen") {
                    if (is.null(object@seq_clusters[["ref"]])) {
                        stop(paste0("====> Error: selected wrong qctype, no ref in object seq_clusters."))
                    }
                }
            }
        }

        if (length(plotdir) == 0) {
            stop(paste0("====> Error: plotdir is not provided, no output directory."))
        }

        if (qctype == "screen") {
            seq_clusters <- object@seq_clusters[["ref"]]
            seq_clusters_1 <- seq_clusters[seq_clusters$cluster==1, ]
            seq_clusters_2 <- seq_clusters[seq_clusters$cluster==2, ]
            seq_clusters_new <- rbind(seq_clusters_1, seq_clusters_2)

            p1 <- ggplot(seq_clusters_new, aes(x = 1:dim(seq_clusters_new)[1], y = count_log2, color = factor(cluster))) +
                    geom_point(shape = 21, size = 1, aes(fill = factor(cluster), color = factor(cluster))) +
                    scale_fill_manual(values = c(t_col("tomato", 0.5), t_col("royalblue", 0.5))) +
                    scale_color_manual(values = c("tomato", "royalblue")) +
                    labs(x = "sequence index", y = "log2(count+1)", title = "Sample QC clusters") +
                    theme(legend.position = "none", panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 12, face = "bold"))

            png(paste0(plotdir, "/", "sample_qc_seq_clusters.point.png"), width = 1200, height = 1200, res = 200)
            print(p1)
            dev.off()

            p2 <- ggplot(seq_clusters_new, aes(x = count_log2, color = factor(cluster))) +
                    geom_histogram(aes(fill = factor(cluster), color = factor(cluster))) +
                    scale_fill_manual(values = c(t_col("tomato", 0.5), t_col("royalblue", 0.5))) +
                    scale_color_manual(values = c("tomato", "royalblue")) +
                    labs(x = "log2(count+1)", y = "frequency", title = "Sample QC clusters") +
                    theme(legend.position = "none", panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 12, face = "bold"))

            #png(paste0(plotdir, "/", "sample_qc_seq_clusters.hist.png"), width = 1200, height = 1200, res = 200)
            #print(p2)
            #dev.off()

            p3 <- ggplot(seq_clusters_new, aes(x = count_log2, color = factor(cluster))) +
                    geom_density(aes(fill = factor(cluster), color = factor(cluster))) +
                    scale_fill_manual(values = c(t_col("tomato", 0.5), t_col("royalblue", 0.5))) +
                    scale_color_manual(values = c("tomato", "royalblue")) +
                    labs(x = "log2(count+1)", y = "frequency", title = "Sample QC clusters") +
                    theme(legend.position = "none", panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 12, face = "bold"))

            png(paste0(plotdir, "/", "sample_qc_seq_clusters.density.png"), width = 1200, height = 1200, res = 200)
            print(p3)
            dev.off()
        } else {
            seq_clusters <- data.table()
            for (i in 1:length(object@seq_clusters)) {
                tmp_cluster <- object@seq_clusters[[i]][, c("count_log2", "cluster")]
                tmp_cluster$samples <- names(object@seq_clusters)[i]
                tmp_cluster <- as.data.table(tmp_cluster)

                if (nrow(seq_clusters) == 0) {
                    seq_clusters <- tmp_cluster
                } else {
                    seq_clusters <- rbind(clusters, tmp_cluster)
                }
            }

            p1 <- ggplot(seq_clusters, aes(x = count_log2, color = samples, fill = samples)) +
                    geom_density(alpha = 0.2) +
                    labs(x = "log2(count+1)", y = "frequency", title = "Sample QC clusters") +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 12, face = "bold")) +
                    facet_wrap(~cluster, scales = "free")

            png(paste0(plotdir, "/", "sample_qc_seq_clusters.density.png"), width = 1200, height = 1200, res = 200)
            print(p1)
            dev.off()
        }
    }
)

#' initialize function
setGeneric("qcplot_stats", function(object, ...) {
  standardGeneric("qcplot_stats")
})

#' create the stats plot
#'
#' @export
#' @param object  sampleQC object
#' @param plotdir the output plot directory
setMethod(
    "qcplot_stats",
    signature = "sampleQC",
    definition = function(object,
                          plotdir) {
        if (length(plotdir) == 0) {
            stop(paste0("====> Error: plotdir is not provided, no output directory."))
        }

        df_total <- object@stats[, c("failed_reads", "filtered_reads")]
        df_total$samples <- rownames(df_total)
        dt_total <- reshape2::melt(as.data.table(df_total), id.vars = "samples", variable.name = "types", value.name = "counts")

        dt_total$samples <- factor(dt_total$samples, levels = df_total$samples)

        p1 <- ggplot(dt_total,  aes(x = samples, y = counts, fill = types)) +
                geom_bar(stat = "identity") +
                scale_fill_manual(values = c(t_col("tomato", 0.5), t_col("royalblue", 0.5))) +
                scale_color_manual(values = c("tomato", "royalblue")) +
                labs(x = "samples", y = "counts", title = "Sample QC Stats") +
                scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
                theme(legend.position = "right", legend.title = element_blank()) +
                theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                theme(axis.text = element_text(size = 12, face = "bold")) +
                theme(axis.text.x = element_text(angle = 90))

        pwidth <- 150 * nrow(df_total)
        png(paste0(plotdir, "/", "sample_qc_stats_total.png"), width = pwidth, height = 1200, res = 200)
        print(p1)
        dev.off()

        df_filtered <- object@stats[, c("per_unmapped_reads", "per_ref_reads", "per_pam_reads", "per_effective_reads")]
        colnames(df_filtered) <- c("unmapped_reads", "ref_reads", "pam_reads", "effective_reads")
        df_filtered <- round(df_filtered*100, 1)
        df_filtered$samples <- rownames(df_filtered)
        dt_filtered <- reshape2::melt(as.data.table(df_filtered), id.vars = "samples", variable.name = "types", value.name = "percent")

        dt_filtered$samples <- factor(dt_filtered$samples, levels = df_filtered$samples)

        gg_colors_fill <- c(t_col("tomato", 0.5), t_col("grey", 0.5), t_col("yellowgreen", 0.5), t_col("royalblue", 0.5))
        gg_colors <- c(c("tomato", "grey", "yellowgreen", "royalblue"))
        p2 <- ggplot(dt_filtered,  aes(x = samples, y = percent, fill = types)) +
                geom_bar(stat = "identity", position = "fill") +
                scale_fill_manual(values = gg_colors_fill) +
                scale_color_manual(values = gg_colors) +
                labs(x = "samples", y = "percent", title = "Sample QC Stats") +
                scale_y_continuous(labels = scales::percent) +
                theme(legend.position = "right", legend.title = element_blank()) +
                theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                theme(axis.text = element_text(size = 12, face = "bold")) +
                theme(axis.text.x = element_text(angle = 90)) +
                geom_text(aes(label = paste0(percent, "%")), position = position_fill(vjust = 0.5), size = 3)

        pwidth <- 150 * nrow(df_filtered)
        png(paste0(plotdir, "/", "sample_qc_stats_filtered.png"), width = pwidth, height = 1200, res = 200)
        print(p2)
        dev.off()

        df_cov <- object@stats[, c("total_reads", "effective_reads", "effective_cov")]
        df_cov$samples <- rownames(df_cov)

        p3 <- ggplot(df_cov,  aes(x = total_reads, y = effective_reads, size = effective_cov, color = samples)) +
                geom_point(alpha = 0.7) +
                labs(x = "total reads", y = "effective reads", title = "Sample QC Stats") +
                scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
                scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
                theme(legend.position = "right") +
                theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                theme(axis.text = element_text(size = 12, face = "bold"))

        png(paste0(plotdir, "/", "sample_qc_stats_cov.png"), width = 1200, height = 1200, res = 200)
        print(p3)
        dev.off()
    }
)

#' initialize function
setGeneric("qcplot_gini", function(object, ...) {
  standardGeneric("qcplot_gini")
})

#' create the gini plot
#'
#' @export
#' @param object  sampleQC object
#' @param plotdir the output plot directory
setMethod(
    "qcplot_gini",
    signature = "sampleQC",
    definition = function(object,
                          plotdir) {
        if (length(plotdir) == 0) {
            stop(paste0("====> Error: plotdir is not provided, no output directory."))
        }

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
        png(paste0(plotdir, "/", "sample_qc_gini.png"), width = pwidth, height = 1200, res = 200)
        print(p1)
        dev.off()
    }
)

#' initialize function
setGeneric("qcplot_position", function(object, ...) {
  standardGeneric("qcplot_position")
})

#' create the position plot
#'
#' @export
#' @param object  sampleQC object
#' @param plotdir the output plot directory
setMethod(
    "qcplot_position",
    signature = "sampleQC",
    definition = function(object,
                          plotdir) {
        if (length(plotdir) == 0) {
            stop(paste0("====> Error: plotdir is not provided, no output directory."))
        }

        # normalisation?
        effcounts_pos <- object@effective_counts_pos
        effcounts_pos <- apply(effcounts_pos, 2, function(x) x / (sum(x, na.rm = TRUE) / 1000000))
        effcounts_pos_log <- log2(effcounts_pos + 1)

        sample_names <- character()
        for (s in object@samples) {
            sample_names <- append(sample_names, s@sample)
        }

        pheight <- 100 * length(sample_names)
        png(paste0(plotdir, "/", "sample_qc_position_cov.heatmap.png"), width = 2400, height = pheight, res = 200)
        lmat <- rbind(c(3, 4), c(2, 1))
        lhei <- c(2, 8)
        lwid <- c(2, 8)

        heatmap.2(t(as.matrix(effcounts_pos_log)),
                  distfun=function(x) dist(x, method = "euclidean"),
                  hclustfun=function(x) hclust(x, method = "ward.D2"),
                  col = colorpanel(100, "royalblue", "ivory", "tomato"),
                  na.color = "grey",
                  breaks = seq(0, 10, length.out = 101),
                  density.info = "none", trace = "none", dendrogram = "none",
                  Rowv = FALSE, Colv = FALSE,
                  labCol = FALSE, cexRow = 1,
                  key = FALSE,
                  #key.xlab = "Log2(count+1)", key.title = "", key.par = list(cex.lab = 1.2),
                  margins = c(4, 8), rowsep = 1:length(sample_names),
                  lmat = lmat, lhei = lhei, lwid = lwid)
        dev.off()

        rownames(effcounts_pos_log) <- 1:nrow(effcounts_pos_log)
        dt_effcounts_pos_log <- reshape2::melt(effcounts_pos_log)
        colnames(dt_effcounts_pos_log) <- c("index", "samples", "log_counts")

        dt_effcounts_pos_log$samples <- factor(dt_effcounts_pos_log$samples, levels = sample_names)

        p1 <- ggplot(dt_effcounts_pos_log, aes(x = index, y = log_counts)) +
                geom_point(shape = 16, size = 0.5, color = "tomato", alpha = 0.8) +
                labs(x = "sequence position index", y = "log2(count+1)", title = "Sample QC position coverage") +
                theme(legend.position = "none", panel.grid.major = element_blank()) +
                theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                theme(axis.text = element_text(size = 6, face = "bold")) +
                facet_wrap(~samples, dir = "v", ncol = 2)

        pheight <- 300 * length(sample_names)
        png(paste0(plotdir, "/", "sample_qc_position_cov.dots.png"), width = 2400, height = pheight, res = 200)
        print(p1)
        dev.off()
    }
)

#' initialize function
setGeneric("qcplot_position_anno", function(object, ...) {
  standardGeneric("qcplot_position_anno")
})

#' create the position plot
#'
#' @export
#' @param object    sampleQC object
#' @param samples   a vector of sample names
#' @param type      plot type, lof or all
#' @param major_cut the qc cutoff of major variants
#' @param plotdir   the output plot directory
setMethod(
    "qcplot_position_anno",
    signature = "sampleQC",
    definition = function(object,
                          samples,
                          type,
                          major_cut = 0.005,
                          plotdir) {
        if (length(plotdir) == 0) {
            stop(paste0("====> Error: plotdir is not provided, no output directory."))
        }

        if (type %nin% c("lof", "all")) {
            stop(paste0("====> Error: wrong type, please use lof or all."))
        }

        effcounts_pos <- object@effective_counts_pos_anno
        effcounts_pos <- effcounts_pos[, c(samples, "consequence")]

        if (type == "lof") {
            effcounts_pos$consequence <- ifelse(effcounts_pos$consequence == "lof", "lof", "others")
            effcounts_pos[, samples] <- effcounts_pos[, samples] / object@stats[samples, ]$filtered_reads * 100

            dt_effcounts_pos <- reshape2::melt(effcounts_pos, id.vars = "consequence", variable.name = "samples", value.name = "counts")
            dt_effcounts_pos$index <- 1:nrow(effcounts_pos)
            dt_effcounts_pos$samples <- factor(dt_effcounts_pos$samples, levels = samples)

            dt_effcounts_pos[dt_effcounts_pos == 0] <- NA

            p1 <- ggplot(dt_effcounts_pos, aes(x = index, y = counts)) +
                    geom_point(shape = 19, size = 0.5, aes(fill = factor(consequence), color = factor(consequence))) +
                    geom_hline(yintercept = major_cut, linetype = "dashed", color = "springgreen4", size = 0.4) +
                    scale_color_manual(values = c(t_col("red", 1), t_col("royalblue", 0.2))) +
                    labs(x = "sequence position", y = "percentage", title = "Sample QC position percentage") +
                    coord_trans(y = "log2") +
                    scale_y_continuous(breaks = c(0.001, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1)) +
                    theme(legend.position = "right", panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 6, face = "bold")) +
                    facet_wrap(~samples, dir = "v")

            pheight <- 300 * length(samples)
            png(paste0(plotdir, "/", "sample_qc_position_anno.lof_dots.png"), width = 2400, height = pheight, res = 200)
            print(p1)
            dev.off()
        } else {
            effcounts_pos[, samples] <- effcounts_pos[, samples] / object@stats[samples, ]$filtered_reads * 100

            dt_effcounts_pos <- reshape2::melt(effcounts_pos, id.vars = "consequence", variable.name = "samples", value.name = "counts")
            dt_effcounts_pos$index <- 1:nrow(effcounts_pos)
            dt_effcounts_pos$samples <- factor(dt_effcounts_pos$samples, levels = samples)

            dt_effcounts_pos[dt_effcounts_pos == 0] <- NA

            default_colors <- c("tomato", "royalblue", "yellowgreen", "orange",
                                "pink", "purple", "coral", "cyan",
                                "violet", "springgreen", "skyblue", "lightgrey")
            select_colors <- default_colors[1:length(unique(effcounts_pos$consequence))]

            freq_cons <- table(effcounts_pos$consequence)
            names(select_colors) <- names(freq_cons)

            freq_cons <- sort(freq_cons, decreasing = TRUE)
            freq_cons <- names(freq_cons)
            rate_cons <- seq(0.2, 0.1 + length(freq_cons)/10, 0.1)
            names(rate_cons) <- freq_cons

            for (i in 1:(length(select_colors) - 1)) {
                select_colors[i] <- t_col(select_colors[i], rate_cons[names(select_colors[i])])
            }
            select_colors <- as.vector(select_colors)

            p1 <- ggplot(dt_effcounts_pos, aes(x = index, y = counts)) +
                    geom_point(shape = 19, size = 0.5, aes(fill = factor(consequence), color = factor(consequence))) +
                    geom_hline(yintercept = major_cut, linetype = "dashed", color = "springgreen4", linewidth = 0.4) +
                    scale_color_manual(values = select_colors) +
                    labs(x = "sequence position", y = "percentage", title = "Sample QC position percentage") +
                    coord_trans(y = "log2") +
                    scale_y_continuous(breaks = c(0.001, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1)) +
                    theme(legend.position = "right", panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 6, face = "bold")) +
                    facet_wrap(~samples, dir = "v")

            pheight <- 300 * length(samples)
            png(paste0(plotdir, "/", "sample_qc_position_anno.all_dots.png"), width = 2400, height = pheight, res = 200)
            print(p1)
            dev.off()
        }
    }
)

#' initialize function
setGeneric("qcplot_dist_samples", function(object, ...) {
  standardGeneric("qcplot_dist_samples")
})

#' create the heatmap of samples
#'
#' @export
#' @param object  sampleQC object
#' @param plotdir the output plot directory
setMethod(
    "qcplot_dist_samples",
    signature = "sampleQC",
    definition = function(object,
                          plotdir) {
        if (length(plotdir) == 0) {
            stop(paste0("====> Error: plotdir is not provided, no output directory."))
        }

        sample_dist <- as.matrix(dist(t(object@deseq_rlog)))
        sample_rlog <- as.matrix(object@deseq_rlog)

        min_rlog <- round(min(sample_rlog))
        max_rlog <- round(max(sample_rlog))

        pwidth <- 100 * ncol(sample_rlog)
        png(paste0(plotdir, "/", "sample_qc_distance_samples.heatmap.png"), width = pwidth, height = 1200, res = 200)
        lmat <- rbind(c(4, 3), c(2, 1))
        lhei <- c(3, 8)
        lwid <- c(3, 8)

        heatmap.2(sample_rlog,
                  distfun = function(x) dist(x, method = "euclidean"),
                  hclustfun = function(x) hclust(x, method = "ward.D2"),
                  col = colorpanel(100, "royalblue", "ivory", "tomato"),
                  na.color = "grey",
                  breaks = seq(min_rlog, max_rlog, length.out = 101),
                  density.info = "none", trace = "none", dendrogram = "both",
                  Rowv = TRUE, Colv = TRUE, labRow = FALSE,
                  cexCol = 0.6, cexRow = 0.6,
                  key.xlab = "deseq rlog", key.title = "", key.par = list(cex.lab = 1),
                  margins = c(8, 1),
                  colsep = 1:ncol(sample_dist),
                  sepwidth=c(0.01, 0.01),
                  lmat = lmat, lhei = lhei, lwid = lwid)
        dev.off()

        sample_corr <- cor(scale(sample_rlog))
        min_corr <- floor(min(sample_corr) * 10) / 10

        png(paste0(plotdir, "/", "sample_qc_distance_samples.corr.png"), width = 1200, height = 1200, res = 200)
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
)

#' initialize function
setGeneric("qcplot_pca_samples", function(object, ...) {
  standardGeneric("qcplot_pca_samples")
})

#' create the pca of samples
#'
#' @export
#' @param object     sampleQC object
#' @param ds_coldata deseq coldata
#' @param ntop       the number of top variances
#' @param plotdir    the output plot directory
setMethod(
    "qcplot_pca_samples",
    signature = "sampleQC",
    definition = function(object,
                          ds_coldata,
                          ntop = 500,
                          plotdir) {
        if (length(plotdir) == 0) {
            stop(paste0("====> Error: plotdir is not provided, no output directory."))
        }

        if ("condition" %nin% colnames(ds_coldata)) {
            stop(paste0("====> Error: coldata must have condition values!"))
        } else {
            ds_coldata <- as.data.frame(ds_coldata)

            ds_coldata$condition <- factor(ds_coldata$condition)
            ds_coldata$condition <- factor(ds_coldata$condition, levels = mixsort(levels(ds_coldata$condition)))

            ds_coldata$replicate <- factor(ds_coldata$replicate)
            ds_coldata$replicate <- factor(ds_coldata$replicate, levels = mixsort(levels(ds_coldata$replicate)))
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

        png(paste0(plotdir, "/", "sample_qc_pca_samples.png"), width = 1200, height = 1200, res = 200)
        par(mfrow = c(2, 2), mar = c(4, 4, 4, 1))
        plot(pca$x[, 1], pca$x[, 2], xlab = "PC1", ylab = "PC2", pch = pca_pchs, col = pca_colors, bg = pca_bgs, lwd = 1, cex = 2, xlim = pc1_set, ylim = pc2_set, main = "PC1 vs PC2")
        plot(pca$x[, 2], pca$x[, 3], xlab = "PC2", ylab = "PC3", pch = pca_pchs, col = pca_colors, bg = pca_bgs, lwd = 1, cex = 2, xlim = pc2_set, ylim = pc3_set, main = "PC2 vs PC3")
        plot(pca$x[, 1], pca$x[, 3], xlab = "PC1", ylab = "PC3", pch = pca_pchs, col = pca_colors, bg = pca_bgs, lwd = 1, cex = 2, xlim = pc1_set, ylim = pc3_set, main = "PC1 vs PC3")
        b <- barplot(percentVar, col = t_col("royalblue", 0.5), border = "royalblue", ylim = c(0, 105))
        text(b, percentVar + 5, paste0(percentVar, "%"), cex = 0.6)
        legend("topright", legend = levels(ds_coldata$replicate), pch = select_pchs, cex = 1, bty = "n")
        legend("right", legend = levels(ds_coldata$condition), pch = 19, col = select_colors, cex = 1, bty = "n")
        dev.off()
    }
)

#' initialize function
setGeneric("qcplot_deseq_fc", function(object, ...) {
  standardGeneric("qcplot_deseq_fc")
})

#' create fold change and consequence plot
#'
#' @export
#' @param object  sampleQC object
#' @param cons    a vector of consequences showed in the figure
#' @param pcut    the padj cutoff
#' @param dcut    the depleted cutoff
#' @param ecut    the enriched cutoff
#' @param plotdir the output plot directory
setMethod(
    "qcplot_deseq_fc",
    signature = "sampleQC",
    definition = function(object,
                          cons = c("synonymous", "lof", "missense"),
                          pcut = 0.05,
                          dcut = 0,
                          ecut = 0,
                          plotdir) {
        if (length(plotdir) == 0) {
            stop(paste0("====> Error: plotdir is not provided, no output directory."))
        }

        comparisions <- names(object@deseq_res)
        for (i in 1:length(object@deseq_res)) {
            res <- object@deseq_res[[i]]$shrunken[, c("log2FoldChange", "padj")]
            res <- as.data.frame(res)

            res$consequence <- object@effective_counts_anno[rownames(res), ]$consequence
            res$stat <- "no impact"
            res[(res$padj < pcut) & (res$log2FoldChange > ecut), ]$stat <- "enriched"
            res[(res$padj < pcut) & (res$log2FoldChange < dcut), ]$stat <- "depleted"

            res_cons <- res[res$consequence %in% cons, ]
            res_cons$stat <- factor(res_cons$stat, levels = c("no impact", "enriched", "depleted"))

            p1 <- ggplot(res_cons, aes(x = consequence, y = log2FoldChange)) +
                    geom_violin(trim = FALSE, scale = "width", fill = 'white', color = "yellowgreen") +
                    geom_quasirandom(size = 1, varwidth = TRUE, method = "tukeyDense", aes(color = factor(stat))) +
                    scale_color_manual(values = c(t_col("grey", 0.3), t_col("tomato", 0.8), t_col("royalblue", 0.8))) +
                    labs(y = "log2FoldChange", title = comparisions[i]) +
                    theme(legend.position = "right", panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 8, face = "bold")) +
                    ylim(-4, 1)

            pwidth <- 300 * length(cons)
            png(paste0(plotdir, "/", "sample_qc_deseq_fc.", comparisions[i], ".beeswarm.png"), width = pwidth, height = 1200, res = 200)
            print(p1)
            dev.off()

            p2 <- ggplot(res_cons, aes(x = consequence, y = log2FoldChange)) +
                    geom_violinhalf(trim = FALSE, scale = "width", fill = t_col("yellowgreen", 0.5), color = "yellowgreen", position = position_nudge(x = .2, y = 0)) +
                    geom_jitter(width = 0.15, size = 0.75, aes(color = factor(stat))) +
                    scale_color_manual(values = c(t_col("grey", 0.3), t_col("tomato", 0.8), t_col("royalblue", 0.8))) +
                    labs(y = "log2FoldChange", title = comparisions[1]) +
                    theme(legend.position = "right", panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 8, face = "bold")) +
                    ylim(-4, 1)

            pwidth <- 300 * length(cons)
            png(paste0(plotdir, "/", "sample_qc_deseq_fc.", comparisions[i], ".violin.png"), width = pwidth, height = 1200, res = 200)
            print(p2)
            dev.off()

            res_cons_volcano <- res_cons
            res_cons_volcano$padj <- -log10(res_cons_volcano$padj)
            p3 <- ggplot(res_cons_volcano, aes(x = log2FoldChange, y = padj)) +
                    geom_point(shape = 19, size = 0.5, aes(color = factor(stat))) +
                    scale_color_manual(values = c(t_col("grey", 0.3), t_col("tomato", 0.8), t_col("royalblue", 0.8))) +
                    labs(x = "log2FoldChange", y = "-log10(padj)", title = comparisions[i]) +
                    theme(legend.position = "right", panel.grid.major = element_blank()) +
                    theme(panel.background = element_rect(fill = "ivory", colour = "white")) +
                    theme(axis.title = element_text(size = 16, face = "bold", family = "Arial")) +
                    theme(plot.title = element_text(size = 16, face = "bold.italic", family = "Arial")) +
                    theme(axis.text = element_text(size = 12, face = "bold")) +
                    xlim(-2, 2) +
                    facet_wrap(~consequence, dir = "v")

            pheight <- 300 * length(cons)
            png(paste0(plotdir, "/", "sample_qc_deseq_fc.", comparisions[i], ".volcano.png"), width = 1200, height = pheight, res = 200)
            print(p3)
            dev.off()
        }
    }
)