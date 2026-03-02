require(statmod)
require(scProportionTest)
require(miloR)
require(speckle)
require(dplyr)
require(ggplot2)
require(sta)
# Now wilcoxon rank-sum
source("pdf_and_png.R")

## Below no longer used - use standard packages
patRatios2pvals <- function(patRatios, groupA, groupB) {
    pvals <- vector()

    for (type in levels(patRatios$cellID)) {
        cellType <- subset(patRatios, subset = cellID == type)
        FP <- subset(cellType, subset = stim == groupA)$ratio
        NP <- subset(cellType, subset = stim == groupB)$ratio
        out <- wilcox.test(NP, FP)
        pvals[type] <- out[[3]]
    }
    return(pvals)
}

patRatios2pvals_welch <- function(patRatios, groupA, groupB) {
    pvals <- vector()

    for (type in levels(patRatios$cellID)) {
        cellType <- subset(patRatios, subset = cellID == type)
        FP <- subset(cellType, subset = stim == groupA)$ratio
        NP <- subset(cellType, subset = stim == groupB)$ratio
        out <- t.test(NP, FP)
        pvals[type] <- out[[3]]
    }
    return(pvals)
}

patRatios2pvals_diff <- function(patRatios, groupA, groupB) {
    pvals <- vector()

    for (type in levels(patRatios$cellID)) {
        cellType <- subset(patRatios, subset = cellID == type)
        FP <- subset(cellType, subset = stim == groupA)$ratio.diff
        NP <- subset(cellType, subset = stim == groupB)$ratio.diff
        out <- t.test(NP, FP)
        pvals[type] <- out[[3]]
    }
    return(pvals)
}

patRatios2pvals_paired_blfu <- function(patRatios) {
    pvals <- vector()
    for (type in levels(patRatios$cellID)) {
        cellType <- subset(patRatios, subset = cellID == type)
        out <- t.test(x = cellType$ratio.bl, y = cellType$ratio.fu, paired = TRUE)
        pvals[type] <- out[[3]]
    }
    return(pvals)
}


do_scPropTest <- function(merged, groupA, groupB, group_var, labelSlot = "seurat_clusters") {
    prop_test <- sc_utils(merged)

    prop_test <- permutation_test(
        prop_test,
        cluster_identity = labelSlot,
        sample_1 = groupA, sample_2 = groupB,
        sample_identity = group_var
    )

    permPlot <- permutation_plot(prop_test)
}


# prop_out - dataframe from propeller. Base proportions under '^BaseProp'. P.Value and FDR for Significance. Pairwise has log ratio
do_propeller <- function(merged, group_var, label_slot = "seurat_clusters", sample_slot = "VisitID", output = ".", SKIP_SAVE_PLOTS = T) {
    if (is.factor(merged@meta.data[[label_slot]])) {
        merged@meta.data[[label_slot]] <- droplevels(merged@meta.data[[label_slot]]) # need to figure out how to handle this generally....
    }
    prop_out <- propeller(
        clusters = merged@meta.data[, label_slot], sample = merged@meta.data[, sample_slot],
        group = merged@meta.data[, group_var]
    )


    props <- getTransformedProps(
        clusters = merged@meta.data[, label_slot], sample = merged@meta.data[, sample_slot],
        transform = "logit"
    )

    transformed_props <- props$Proportions |> as.data.frame()
    per_patientMD <- merged@meta.data %>% distinct(.data[[sample_slot]], .keep_all = T)
    per_pat_info <- data.frame(sample = per_patientMD[[sample_slot]], stim = per_patientMD[[group_var]])
    transformed_props <- left_join(transformed_props, per_pat_info, by = "sample")

    summary_df <- transformed_props %>%
        dplyr::group_by(stim, clusters) %>%
        dplyr::summarise(
            mean = mean(Freq), std_error = plotrix::std.error(Freq),
            n = n(),
            .groups = "drop_last"
        )

    summary_df <- summary_df %>%
        mutate(error = std_error)

    plotList <- lapply(transformed_props$clusters, function(i) {
        summary_df_sub <- summary_df %>% subset(subset = (clusters %in% i), drop = T)

        patRatios_sub <- subset(transformed_props,
            subset = (clusters %in% i),
            drop = T
        )
        limitsub <- aes(ymax = mean + error, ymin = mean - error)
        gtmp <- ggplot(summary_df_sub, aes(x = stim, y = mean, fill = stim))
        # g <- g + geom_boxplot(notch = TRUE)
        gtmp <- gtmp + geom_bar(stat = "identity", position = position_dodge())
        gtmp <- gtmp + geom_errorbar(limitsub, width = 0.25, position = position_dodge(width = 0.9))
        gtmp <- gtmp + geom_jitter(width = 0.25, data = patRatios_sub, shape = 1, aes(x = stim, y = Freq))
        gtmp <- gtmp + facet_grid(~clusters) + scale_color_prism("floral") +
            scale_fill_prism("floral") +
            guides(y = "prism_offset_minor") +
            theme_prism(base_size = 16) +
            theme(legend.position = "none")
        if (!SKIP_SAVE_PLOTS) {
            pdf_and_png(gtmp, output, paste0(i, "_only"), pdfWidth = 4, pdfHeight = 5)
        }
        return(gtmp)
    })

    names(plotList) <- transformed_props$clusters
    return(list("prop_out" = prop_out, "plotList" = plotList))
}

do_propeller_ttest_mmrf <- function(merged, group_var, batch_var = "siteXbatch", label_slot = "seurat_clusters", sample_slot = "VisitID", output = ".", SKIP_SAVE_PLOTS = F) {
    if (is.factor(merged@meta.data[[label_slot]])) {
        merged@meta.data[[label_slot]] <- droplevels(merged@meta.data[[label_slot]]) # need to figure out how to handle this generally....
    }

    design <- model.matrix(~ 0 + !!group_var + !!batch_var)
    numGroups <- length(unique(merged@meta.data[, group_var]))

    props <- getTransformedProps(
        clusters = merged@meta.data[, label_slot], sample = merged@meta.data[, sample_slot],
        transform = "logit"
    )
    if (numGroups > 2) {
        prop_out <- propeller.anova(prop.list = props, design = design, coef = c(1:length(numGroups)), robust = T, trend = F, sort = T)
    } else {
        contr <- makeContrasts(!!colnames(design1) - !!colnames(design2), levels = design)
        prop_out <- propeller.ttest(props, design, contrasts = contr, robust = T, trend = F, sort = T)
    }

    transformed_props <- props$Proportions |> as.data.frame()
    per_patientMD <- merged@meta.data %>% distinct(.data[[sample_slot]], .keep_all = T)
    per_pat_info <- data.frame(sample = per_patientMD[[sample_slot]], stim = per_patientMD[[group_var]])
    transformed_props <- left_join(transformed_props, per_pat_info, by = "sample")

    summary_df <- transformed_props %>%
        dplyr::group_by(stim, clusters) %>%
        dplyr::summarise(
            mean = mean(Freq), std_error = plotrix::std.error(Freq),
            n = n(),
            .groups = "drop_last"
        )

    summary_df <- summary_df %>%
        mutate(error = std_error)

    plotList <- lapply(transformed_props$clusters, function(i) {
        summary_df_sub <- summary_df %>% subset(subset = (clusters %in% i), drop = T)

        patRatios_sub <- subset(transformed_props,
            subset = (clusters %in% i),
            drop = T
        )
        limitsub <- aes(ymax = mean + error, ymin = mean - error)
        gtmp <- ggplot(summary_df_sub, aes(x = stim, y = mean, fill = stim))
        # g <- g + geom_boxplot(notch = TRUE)
        gtmp <- gtmp + geom_bar(stat = "identity", position = position_dodge())
        gtmp <- gtmp + geom_errorbar(limitsub, width = 0.25, position = position_dodge(width = 0.9))
        gtmp <- gtmp + geom_jitter(width = 0.25, data = patRatios_sub, shape = 1, aes(x = stim, y = Freq))
        gtmp <- gtmp + facet_grid(~clusters) + scale_color_prism("floral") +
            scale_fill_prism("floral") +
            guides(y = "prism_offset_minor") +
            theme_prism(base_size = 16) +
            theme(legend.position = "none")
        if (!SKIP_SAVE_PLOTS) {
            pdf_and_png(gtmp, output, paste0(i, "_only"), pdfWidth = 4, pdfHeight = 5)
        }
        return(gtmp)
    })

    return(list("prop_out" = prop_out, "plotList" = plotList))
}

library(gtools) # Figure out why geom_boxplot only draws median
barChartBuilder <- function(
    merged, group_var, groupA = "", groupB = "", string_mod = "",
    patientSlot = "MMRFID", mmrfid_slot = "VisitID", labelSlot = "seurat_clusters",
    subset_bar_chart_list = NULL, width = 40, height = 5, output = ".", SKIP_SAVE_PLOTS = T) {
    # noMSSM <- subset(noMSSM, subset = (cellID %ni%
    # c('Erythrocyte','Erythroblast')))
    if (is.factor(merged@meta.data[[labelSlot]])) {
        merged@meta.data[[labelSlot]] <- droplevels(merged@meta.data[[labelSlot]]) # need to figure out how to handle this generally....
    }
    patTable <- table(merged@meta.data[[mmrfid_slot]], merged@meta.data[[labelSlot]])
    append_str <- labelSlot
    patTable_orig <- patTable
    patTable_orig
    write.csv(patTable, paste0(output, paste0(
        "/patRatios", append_str, "_",
        string_mod, "_allsites.csv"
    )))

    patSums <- rowSums(patTable)
    patRatTable <- patTable / patSums

    patRatios <- as.data.frame(patTable / patSums)
    colnames(patRatios) <- c("Pat", "cellID", "ratio")
    # Pull from metadata slot instead
    per_patientMD <- merged@meta.data %>% distinct(.data[[mmrfid_slot]], .keep_all = T)

    per_pat_info <- data.frame(Pat = per_patientMD[[mmrfid_slot]], stim = per_patientMD[[group_var]])

    patRatios <- left_join(patRatios, per_pat_info, by = "Pat")
    # patRatios <- patRatios %>% mutate(Resp_Group =
    # sample_info[as.vector(Pat), 'Response_Group'])

    if (!is.null(subset_bar_chart_list)) {
        patRatios <- subset(patRatios,
            subset = (cellID %in% subset_bar_chart_list),
            drop = T
        )
        patRatios <- droplevels(patRatios)
    }

    summary_df <- patRatios %>%
        dplyr::group_by(stim, cellID) %>%
        dplyr::summarise(
            mean = mean(ratio), std_error = plotrix::std.error(ratio),
            n = n(),
            .groups = "drop_last"
        )

    summary_df <- summary_df %>%
        mutate(error = std_error)

    head(summary_df)
    library(reshape2)
    # df2 <- melt(summary_df, id.vars = c('cellID', 'stim'))


    # print(head(summary_df))
    # print(head(patRatios))
    # print(patRatios$cellID)
    # print(summary_df$cellID)

    limits <- aes(ymax = mean + error, ymin = mean - error)
    g <- ggplot(summary_df, aes(x = stim, y = mean, fill = stim))
    # g <- g + geom_boxplot(colour = "black")
    g <- g + geom_bar(stat = "identity", position = position_dodge())
    g <- g + geom_errorbar(limits, width = 0.25, position = position_dodge(width = 0.9))
    g <- g + geom_jitter(width = 0.25, data = patRatios, shape = 1, aes(x = stim, y = ratio))
    g <- g + facet_grid(~cellID) + scale_color_prism("floral") +
        scale_fill_prism("floral") +
        guides(y = "prism_offset_minor") +
        theme_prism(base_size = 16) +
        theme(legend.position = "none")
    g

    if (!SKIP_SAVE_PLOTS) {
        pdf_and_png(g, output, paste0("BarChart_Full_", string_mod), pdfWidth = width, pdfHeight = height)
    }
    # future.apply::future_
    cells_to_loop <- unique(patRatios$cellID)
    plotList <- lapply(cells_to_loop, function(i) {
        patRatios_sub <- subset(patRatios,
            subset = (cellID %in% i),
            drop = T
        )
        summary_df_sub <- summary_df %>% subset(subset = (cellID %in% i), drop = T)

        limitsub <- aes(ymax = mean + error, ymin = mean - error)
        gtmp <- ggplot(summary_df_sub, aes(x = stim, y = mean, fill = stim))
        # gtmp <- gtmp + geom_boxplot()
        gtmp <- gtmp + geom_bar(stat = "identity", position = position_dodge())
        gtmp <- gtmp + geom_errorbar(limitsub, width = 0.25, position = position_dodge(width = 0.9))
        gtmp <- gtmp + geom_jitter(width = 0.25, data = patRatios_sub, shape = 1, aes(x = stim, y = ratio))
        gtmp <- gtmp + facet_grid(~cellID) + scale_color_prism("floral") +
            scale_fill_prism("floral") +
            guides(y = "prism_offset_minor") +
            theme_prism(base_size = 16) +
            theme(legend.position = "none")
        if (!SKIP_SAVE_PLOTS) {
            pdf_and_png(gtmp, output, paste0(i, "_only"), pdfWidth = 4, pdfHeight = 5)
        }
        return(gtmp)
    })
    names(plotList) <- cells_to_loop



    # tryCatch(exp = {
    if (groupA != "" && groupB != "") {
        pvalDF <- patRatios2pvals(patRatios, groupA, groupB)
        write.csv(pvals, paste0(
            output, "/welch_pval_table", append_str, "_",
            group_var, string_mod, ".csv"
        ))
    } else {
        print("Doing all group_var combos")
        numCombos <- length(unique(patRatios$stim))
        comboSize <- 2
        allComparisons <- gtools::combinations(numCombos, comboSize, unique(patRatios$stim))
        pvalList <- list()
        pvalList_welch <- list()
        for (idx in 1:nrow(allComparisons)) {
            group1 <- allComparisons[idx, 1]
            group2 <- allComparisons[idx, 2]
            pvalList[[idx]] <- c(patRatios2pvals(patRatios, group1, group2))
            pvalList_welch[[idx]] <- c(patRatios2pvals_welch(patRatios, group1, group2))
        }
        colNameList <- paste0(allComparisons[, 1], "-vs-", allComparisons[
            ,
            2
        ])
        pvalDF <- do.call(data.frame, pvalList)
        colnames(pvalDF) <- colNameList
        write.csv(pvalDF, paste0(
            output, "/wilcox_pval_table_mix", append_str,
            "_", group_var, string_mod, ".csv"
        ))
        pvalDF_welch <- do.call(data.frame, pvalList_welch)
        colnames(pvalDF_welch) <- colNameList
        write.csv(pvalDF_welch, paste0(
            output, "/welch_pval_table_mix", append_str,
            "_", group_var, string_mod, ".csv"
        ))
    }
    # }, error = function(e) {
    #     print(e)
    # }, warning = function(e) {
    #     print(e)
    # })

    returnList <- list("PatRatios" = patRatios, "FullPlot" = g, "PlotList" = plotList, "pvals" = pvalDF, "pvals_welch" = pvalDF_welch, "summary_df" = summary_df)
    return(returnList)
}

barChartBuilder_paired <- function(
    merged, group_var, time_var, time1 = "Baseline", time2 = "Followup", string_mod = "",
    patientSlot = "MMRFID", mmrfid_slot = "VisitID", labelSlot = "seurat_clusters",
    width = 40, height = 5, output = ".") {
    # noMSSM <- subset(noMSSM, subset = (cellID %ni%
    # c('Erythrocyte','Erythroblast')))
    merged@meta.data[[labelSlot]] <- droplevels(merged@meta.data[[labelSlot]])
    patTable <- table(merged@meta.data[[mmrfid_slot]], merged@meta.data[[labelSlot]])
    append_str <- labelSlot
    patTable_orig <- patTable
    patTable_orig
    write.csv(patTable, paste0(output, paste0(
        "/patRatios", append_str, "_",
        string_mod, "_allsites.csv"
    )))

    patSums <- rowSums(patTable)
    patRatTable <- patTable / patSums

    patRatios <- as.data.frame(patTable / patSums)
    colnames(patRatios) <- c("Sample", "cellID", "ratio")
    # Pull from metadata slot instead

    per_patientMD <- merged@meta.data %>% distinct(.data[[mmrfid_slot]], .keep_all = T)

    per_pat_info <- data.frame(Sample = per_patientMD[[mmrfid_slot]], Pat = per_patientMD[[patientSlot]], stim = per_patientMD[[group_var]], visit_type = per_patientMD[[time_var]])
    patRatios <- left_join(patRatios, per_pat_info, by = "Sample")
    patRatios$hybrid_stim <- paste0(patRatios$stim, ".", patRatios$visit_type)

    patRatios_bl <- patRatios %>% dplyr::filter(visit_type == time1)
    patRatios_fu <- patRatios %>% dplyr::filter(visit_type == time2)

    patRatios_bl <- dplyr::select(patRatios_bl, -c("Sample"))
    patRatios_bl$pat.cellID <- paste0(patRatios_bl$Pat, ".", patRatios_bl$cellID)
    patRatios_fu <- dplyr::select(patRatios_fu, -c("Sample"))
    patRatios_fu$pat.cellID <- paste0(patRatios_fu$Pat, ".", patRatios_fu$cellID)

    rownames(patRatios_bl) <- patRatios_bl$pat.cellID
    rownames(patRatios_fu) <- patRatios_fu$pat.cellID

    patRatios_paired <- patRatios_bl
    colnames(patRatios_paired)[names(patRatios_paired) == "ratio"] <- "ratio.bl"
    patRatios_paired$ratio.fu <- patRatios_fu[rownames(patRatios_paired), "ratio"]

    patRatios_paired$ratio.diff <- patRatios_paired$ratio.fu - patRatios_paired$ratio.bl
    # patRatios - grouped by Sample, ratios in ratio, hybrid_stim: RP.<visit_type>. visit_type = baseline, followup
    # patRatios_paired - grouped by Pat, ratios in ratio.bl and ratio.fu, stim = RP or NP

    summary_df <- patRatios %>%
        dplyr::group_by(hybrid_stim, cellID) %>%
        dplyr::summarise(
            mean = mean(ratio), std_error = plotrix::std.error(ratio),
            n = n(),
            .groups = "drop_last"
        )

    summary_df <- summary_df %>%
        mutate(error = std_error)

    head(summary_df)
    library(reshape2)

    limits <- aes(ymax = mean + error, ymin = mean - error)
    g <- ggplot(summary_df, aes(x = hybrid_stim, y = mean, fill = hybrid_stim))
    # g <- g + geom_boxplot(notch = TRUE)
    g <- g + geom_bar(stat = "identity", position = position_dodge())
    g <- g + geom_errorbar(limits, width = 0.25, position = position_dodge(width = 0.9))
    g <- g + geom_point(data = patRatios, shape = 1, aes(x = hybrid_stim, y = ratio))
    g <- g + facet_grid(~cellID)

    pdf_and_png(g, output, paste0("Barchart_blfu", string_mod), pdfWidth = width, pdfHeight = height)

    summary_df_paired <- patRatios_paired %>%
        dplyr::group_by(stim, cellID) %>%
        dplyr::summarise(
            mean = mean(ratio.diff), error = plotrix::std.error(ratio.diff),
            n = n(),
            .groups = "drop_last"
        )

    limits_paired <- aes(ymax = mean + error, ymin = mean - error)
    g2 <- ggplot(summary_df_paired, aes(x = stim, y = mean, fill = stim))
    # g <- g + geom_boxplot(notch = TRUE)
    g2 <- g2 + geom_bar(stat = "identity", position = position_dodge())
    g2 <- g2 + geom_errorbar(limits_paired, width = 0.25, position = position_dodge(width = 0.9))
    g2 <- g2 + geom_point(data = patRatios_paired, shape = 1, aes(x = stim, y = ratio.diff))
    g2 <- g2 + facet_grid(~cellID)

    pdf_and_png(g2, output, paste0("Barchart_diff", string_mod), pdfWidth = width, pdfHeight = height)

    g3 <- ggplot(patRatios, aes(x = visit_type, y = ratio))
    g3 <- g3 + geom_boxplot(notch = FALSE)
    g3 <- g3 + geom_point(shape = 1, aes(x = visit_type, y = ratio, group = Pat, fill = Pat, color = Pat))
    g3 <- g3 + geom_line(aes(x = visit_type, y = ratio, group = Pat, color = Pat)) # data=patRatios, aes(x=Pat, y=ratio))
    g3 <- g3 + facet_grid(stim ~ cellID)

    pdf_and_png(g3, output, paste0("Barchart_box_paired", string_mod), pdfWidth = width, pdfHeight = height * 3)


    allPlots <- lapply(unique(patRatios$cellID), function(i) {
        patRatios_sub <- subset(patRatios,
            subset = (cellID %in% i),
            drop = T
        )

        patRatios_paired_sub <- subset(patRatios_paired, subset = (cellID %in% i), drop = T)

        summary_df_sub <- summary_df %>% subset(subset = (cellID %in% i), drop = T)
        summary_df_paired_sub <- summary_df_paired %>% subset(subset = (cellID %in% i), drop = T)

        limitsub <- aes(ymax = mean + error, ymin = mean - error)
        gtmp <- ggplot(summary_df_sub, aes(x = hybrid_stim, y = mean, fill = hybrid_stim))
        # g <- g + geom_boxplot(notch = TRUE)
        gtmp <- gtmp + geom_bar(stat = "identity", position = position_dodge())
        gtmp <- gtmp + geom_errorbar(limitsub, width = 0.25, position = position_dodge(width = 0.9))
        gtmp <- gtmp + geom_point(data = patRatios_sub, shape = 1, aes(x = hybrid_stim, y = ratio))
        gtmp <- gtmp + facet_grid(~cellID)


        # pdf_and_png(gtmp, output, paste0(i, "_only"), pdfWidth = 4, pdfHeight = 5)

        limits_paired <- aes(ymax = mean + error, ymin = mean - error)
        gtmp2 <- ggplot(summary_df_paired_sub, aes(x = stim, y = mean, fill = stim))
        # g <- g + geom_boxplot(notch = TRUE)
        gtmp2 <- gtmp2 + geom_bar(stat = "identity", position = position_dodge())
        gtmp2 <- gtmp2 + geom_errorbar(limits_paired, width = 0.25, position = position_dodge(width = 0.9))
        gtmp2 <- gtmp2 + geom_point(data = patRatios_paired_sub, shape = 1, aes(x = stim, y = ratio.diff))
        gtmp2 <- gtmp2 + facet_grid(~cellID)
        # pdf_and_png(gtmp2, output, paste0(i, "_only_diff"), pdfWidth = 4, pdfHeight = 5)

        gtmp3 <- ggplot(patRatios_sub, aes(x = visit_type, y = ratio))
        gtmp3 <- gtmp3 + geom_boxplot(notch = FALSE)
        gtmp3 <- gtmp3 + geom_point(shape = 1, aes(x = visit_type, y = ratio, group = Pat, fill = Pat, color = Pat))
        gtmp3 <- gtmp3 + geom_line(aes(x = visit_type, y = ratio, group = Pat, color = Pat)) # data=patRatios, aes(x=Pat, y=ratio))
        gtmp3 <- gtmp3 + facet_grid(~stim)
        # pdf_and_png(gtmp3, output, paste0(i, "_only_box"), pdfWidth = 7, pdfHeight = 5)

        return(list("g1" = gtmp, "g2" = gtmp2, "g3" = gtmp3))
    })

    plotList <- lapply(allPlots, function(t) t$g1)
    plotList_diff <- lapply(allPlots, function(t) t$g3)
    plotList_box <- lapply(allPlots, function(t) t$g3)

    rm(allPlots)


    # tryCatch(exp = {
    print("Doing all group_var combos")
    numCombos <- length(unique(patRatios$stim))
    comboSize <- 2
    pval_blfu_paired <- list()
    for (idx in unique(patRatios$stim)) {
        pval_blfu_paired[[idx]] <- patRatios2pvals_paired_blfu(patRatios_paired %>% dplyr::filter(stim == idx))
    }
    pval_blfu_paired <- do.call(data.frame, pval_blfu_paired)

    write.csv(pval_blfu_paired, paste0(
        output, "/welch_pval_table_blfu", append_str,
        "_", group_var, string_mod, ".csv"
    ))

    allComparisons <- gtools::combinations(numCombos, comboSize, unique(patRatios$stim))
    pvalList <- list()
    for (idx in 1:nrow(allComparisons)) {
        group1 <- allComparisons[idx, 1]
        group2 <- allComparisons[idx, 2]


        pvalList[[idx]] <- patRatios2pvals_diff(patRatios_paired, groupA = group1, groupB = group2)
    }
    colNameList <- paste0(allComparisons[, 1], "-vs-", allComparisons[
        ,
        2
    ])
    pvalDF <- do.call(data.frame, pvalList)
    colnames(pvalDF) <- colNameList
    # print(rownames(pvalDF))
    write.csv(pvalDF, paste0(
        output, "/welch_pval_table_mix", append_str,
        "_", group_var, string_mod, ".csv"
    ))
    # }, error = function(e) {
    #     print(e)
    # }, warning = function(e) {
    #     print(e)
    # })

    returnList <- list("PatRatios" = patRatios, "FullPlot" = g, "FullPlot_diff" = g2, "FullPlot_box" = g3, "PlotList" = plotList, "PlotList_diff" = plotList_diff, "PlotList_box" = plotList_box, "pvals.paired" = pval_blfu_paired, "pvals.diff" = pvalDF)
    return(returnList)
}


doMilo <- function(
    seurat_object, string_append, group_var, samples_var, output,
    REBUILD_OBJECT = F, dims = 30, k = 75, prop = 0.05) {
    if (!REBUILD_OBJECT && file.exists(file.path(output, paste0(
        "MILO", string_append,
        "traj.rds"
    )))) {
        print("Reloading")
        traj_milo <- readRDS(file.path(output, paste0("MILO", string_append, "traj.rds")))
        da_results <- readRDS(file.path(output, paste0("MILO", string_append, "results.rds")))
    } else {
        traj_milo <- Milo(as.SingleCellExperiment(seurat_object))
        traj_milo <- buildGraph(traj_milo, k = k, d = dims)
        traj_milo <- makeNhoods(traj_milo, prop = prop, k = k, d = dims, refined = TRUE)

        p <- plotNhoodSizeHist(traj_milo)
        pdf(file.path(output, paste0("MILO_", string_append, "nhoodsize.pdf")),
            width = 7, height = 7
        )
        print(p)
        dev.off()

        traj_milo <- countCells(traj_milo,
            meta.data = seurat_object@meta.data,
            samples = samples_var
        )
        # Check the below later
        seurat_object@meta.data$Comparison <- seurat_object@meta.data[[group_var]]
        traj_design <- seurat_object@meta.data[, c(samples_var, "Comparison")]
        traj_design <- distinct(traj_design)
        rownames(traj_design) <- traj_design[[samples_var]]
        ## Reorder rownames to match columns of nhoodCounts(milo)
        traj_design <- traj_design[colnames(nhoodCounts(traj_milo)), , drop = FALSE]

        traj_design


        traj_milo <- calcNhoodDistance(traj_milo, d = dims)
        rownames(traj_design) <- traj_design[[samples_var]]
        da_results <- testNhoods(traj_milo, design = ~Comparison, design.df = traj_design)

        traj_milo <- buildNhoodGraph(traj_milo)

        da_results %>%
            arrange(-SpatialFDR) %>%
            head()
        saveRDS(traj_milo, file.path(output, paste0("MILO", string_append, "traj.rds")))
        saveRDS(da_results, file.path(output, paste0("MILO", string_append, "results.rds")))
    }
    p <- scater::plotUMAP(traj_milo, colour_by = "seurat_clusters") + plotNhoodGraphDA(traj_milo,
        da_results,
        alpha = 0.05
    ) + plot_layout(guides = "collect")

    pdf(file.path(output, paste0("MILO", string_append, "_TEST_CLINICAL_GROUP.pdf", k)),
        width = 12
    )
    print(p)
    dev.off()
}
