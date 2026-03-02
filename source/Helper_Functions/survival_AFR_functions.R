library(ggpubr)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Comprehensive survival analysis function
perform_gene_survival_analysis <- function(sce_object = NULL, gene_name = NULL, expression_values = NULL, meta_data = NULL, covariates = NULL, cutpoint_method = "median") {
    # Check input parameters
    if (is.null(sce_object) && (is.null(expression_values) || is.null(meta_data))) {
        stop("Must provide either 'sce_object' and 'gene_name', or 'expression_values' and 'meta_data'")
    }

    if (!is.null(sce_object) && (!is.null(expression_values) || !is.null(meta_data))) {
        warning("Both SCE and vector inputs provided. Using SCE object.")
    }

    # Extract expression data and metadata based on input type
    if (!is.null(sce_object)) {
        # SCE object approach
        if (is.null(gene_name)) {
            stop("gene_name must be provided when using sce_object")
        }
        expr_data <- assay(sce_object, "normalized")[gene_name, ]
        meta_data <- as.data.frame(colData(sce_object))
        analysis_name <- gene_name
    } else {
        # Vector + metadata approach
        expr_data <- expression_values
        analysis_name <- deparse(substitute(expression_values))

        # Ensure same number of samples
        if (length(expr_data) != nrow(meta_data)) {
            stop("Length of expression_values must match number of rows in meta_data")
        }
    }

    # Add expression data to metadata
    meta_data$gene_expr <- expr_data
    meta_data$AFR_group <- ifelse(meta_data$ANCESTRY.AFR < 0.5, "AFR Low", "AFR High")
    meta_data$AFR_group <- factor(meta_data$AFR_group, levels = c("AFR Low", "AFR High"))
    # Remove samples with missing survival data
    meta_data <- meta_data[!is.na(meta_data$censpfs) & !is.na(meta_data$ttcpfs), ]

    meta_data$censpfs <- as.numeric(meta_data$censpfs)

    meta_data <- meta_data |> dplyr::filter(!is.na(censpfs) & !is.na(ttcpfs))


    if (nrow(meta_data) == 0) {
        stop("No samples with complete survival data")
    }

    # Initialize results list
    results <- list()
    results$analysis_name <- analysis_name
    results$cutpoint_method <- cutpoint_method
    results$n_samples <- nrow(meta_data)

    # Determine cutpoints based on method
    if (cutpoint_method == "median") {
        # Calculate medians for different groups
        cutpoints <- list(
            all = median(meta_data$gene_expr, na.rm = TRUE),
            afr_low = median(meta_data$gene_expr[meta_data$AFR_group == "AFR Low"], na.rm = TRUE),
            afr_high = median(meta_data$gene_expr[meta_data$AFR_group == "AFR High"], na.rm = TRUE)
        )
        results$cutpoints <- cutpoints
    } else if (cutpoint_method == "optimal") {
        # Function to test cutpoints between 25th and 75th percentile
        find_optimal_cutpoint <- function(data, expr_col, event_col, time_col, covars = NULL) {
            expr_vals <- data[[expr_col]]
            q25 <- quantile(expr_vals, 0.25, na.rm = TRUE)
            q75 <- quantile(expr_vals, 0.75, na.rm = TRUE)

            # Test cutpoints in 1 percentile intervals
            cutpoint_range <- seq(q25, q75, length.out = 51) # 1% intervals

            cutpoint_results <- data.frame(
                cutpoint = cutpoint_range,
                hr = NA,
                pvalue = NA,
                ci_lower = NA,
                ci_upper = NA
            )

            for (i in seq_along(cutpoint_range)) {
                cutpoint <- cutpoint_range[i]
                data$expr_group <- ifelse(data[[expr_col]] > cutpoint, "High", "Low")
                data$expr_group <- factor(data$expr_group, levels = c("Low", "High"))

                # Build Cox model formula
                if (is.null(covars)) {
                    formula_str <- paste("Surv(", time_col, ",", event_col, ") ~ expr_group")
                } else {
                    covar_str <- paste(covars, collapse = " + ")
                    formula_str <- paste("Surv(", time_col, ",", event_col, ") ~ expr_group +", covar_str)
                }

                tryCatch(
                    {
                        cox_model <- coxph(as.formula(formula_str), data = data)
                        cox_summary <- summary(cox_model)

                        cutpoint_results$hr[i] <- cox_summary$coefficients["expr_groupHigh", "exp(coef)"]
                        cutpoint_results$pvalue[i] <- cox_summary$coefficients["expr_groupHigh", "Pr(>|z|)"]
                        cutpoint_results$ci_lower[i] <- cox_summary$conf.int["expr_groupHigh", "lower .95"]
                        cutpoint_results$ci_upper[i] <- cox_summary$conf.int["expr_groupHigh", "upper .95"]
                    },
                    error = function(e) {
                        # Skip if model fails
                    }
                )
            }

            # Find optimal cutpoint (smallest p-value)
            best_idx <- which.min(cutpoint_results$pvalue)
            optimal_cutpoint <- cutpoint_results$cutpoint[best_idx]

            return(list(
                optimal_cutpoint = optimal_cutpoint,
                cutpoint_results = cutpoint_results
            ))
        }

        # Find optimal cutpoints for different groups
        opt_all <- find_optimal_cutpoint(meta_data, "gene_expr", "censpfs", "ttcpfs", covariates)

        # Check if AFR groups have sufficient data
        afr_low_data <- meta_data[meta_data$AFR_group == "AFR Low", ]
        afr_high_data <- meta_data[meta_data$AFR_group == "AFR High", ]

        if (nrow(afr_low_data) < 10) {
            warning("Insufficient AFR Low samples for optimal cutpoint analysis")
            opt_afr_low <- list(
                optimal_cutpoint = median(afr_low_data$gene_expr, na.rm = TRUE),
                cutpoint_results = data.frame(
                    cutpoint = numeric(0), hr = numeric(0),
                    pvalue = numeric(0), ci_lower = numeric(0), ci_upper = numeric(0)
                )
            )
        } else {
            opt_afr_low <- find_optimal_cutpoint(afr_low_data, "gene_expr", "censpfs", "ttcpfs", covariates)
        }

        if (nrow(afr_high_data) < 10) {
            warning("Insufficient AFR High samples for optimal cutpoint analysis")
            opt_afr_high <- list(
                optimal_cutpoint = median(afr_high_data$gene_expr, na.rm = TRUE),
                cutpoint_results = data.frame(
                    cutpoint = numeric(0), hr = numeric(0),
                    pvalue = numeric(0), ci_lower = numeric(0), ci_upper = numeric(0)
                )
            )
        } else {
            opt_afr_high <- find_optimal_cutpoint(afr_high_data, "gene_expr", "censpfs", "ttcpfs", covariates)
        }

        cutpoints <- list(
            all = opt_all$optimal_cutpoint,
            afr_low = opt_afr_low$optimal_cutpoint,
            afr_high = opt_afr_high$optimal_cutpoint
        )

        results$cutpoints <- cutpoints
        results$cutpoint_analysis <- list(
            all = opt_all$cutpoint_results,
            afr_low = opt_afr_low$cutpoint_results,
            afr_high = opt_afr_high$cutpoint_results
        )
    } else if (cutpoint_method == "categorical") {
        cutpoints <- list(
            all = NA,
            afr_low = NA,
            afr_high = NA
        )
        results$cutpoints <- cutpoints
    } else {
        stop("Invalid cutpoint_method. Choose 'median', 'optimal', or 'categorical'.")
    }

    if (cutpoint_method %in% c("median", "optimal")) {
        # Create box plot
        p_boxplot <- ggplot(meta_data, aes(x = AFR_group, y = gene_expr, fill = AFR_group)) +
            geom_boxplot(alpha = 0.7) +
            geom_jitter(width = 0.2, alpha = 0.5) +
            scale_fill_manual(values = c("AFR Low" = fill_AFR_low, "AFR High" = fill_AFR_high)) +
            scale_color_manual(values = c("AFR Low" = color_AFR_low, "AFR High" = color_AFR_high)) +
            labs(
                title = paste("Expression:", analysis_name),
                x = "AFR Ancestry Group",
                y = "Log Normalized Expression"
            ) +
            theme_classic() +
            theme(legend.position = "none")
    } else if (cutpoint_method %in% c("categorical")) {
        # Create confusion matrix plot showing frequency of categories by AFR group
        confusion_df <- as.data.frame(table(meta_data$AFR_group, meta_data$gene_expr))
        colnames(confusion_df) <- c("AFR_group", "gene_expr", "Freq")

        # Calculate percentage within each AFR group
        confusion_df <- confusion_df %>%
            group_by(AFR_group) %>%
            mutate(
                Total_in_Group = sum(Freq),
                Percent = round((Freq / Total_in_Group) * 100, 1)
            ) %>%
            ungroup()

        p_boxplot <- ggplot(confusion_df, aes(x = AFR_group, y = gene_expr, fill = Percent)) +
            geom_tile(color = "white", size = 0.5) +
            geom_text(aes(label = paste0(Freq, "\n(", Percent, "%)")),
                color = "black", size = 4, fontface = "bold"
            ) +
            scale_fill_gradient(low = "#f7f7f7", high = "#2166ac", name = "Percent") +
            labs(
                title = paste("Category Frequency by AFR Group:", analysis_name),
                x = "AFR Ancestry Group",
                y = "Gene Expression Category"
            ) +
            theme_classic() +
            theme(
                axis.text.x = element_text(angle = 0, hjust = 0.5),
                axis.text.y = element_text(angle = 0, hjust = 1),
                legend.position = "right"
            )
    }

    # Add cutpoint lines based on method (only for continuous methods - not relevant for categorical)
    if (cutpoint_method == "median") {
        p_boxplot <- p_boxplot +
            geom_hline(yintercept = cutpoints$all, linetype = "dashed", color = "black", size = 1) +
            annotate("text", x = 1.5, y = cutpoints$all, label = "Overall Median", vjust = -0.5)
    } else if (cutpoint_method == "optimal") {
        p_boxplot <- p_boxplot +
            geom_hline(yintercept = cutpoints$all, linetype = "dashed", color = "black", size = 1) +
            geom_hline(yintercept = cutpoints$afr_low, linetype = "dashed", color = fill_AFR_low, size = 1) +
            geom_hline(yintercept = cutpoints$afr_high, linetype = "dashed", color = fill_AFR_high, size = 1) +
            annotate("text", x = 1.5, y = cutpoints$all, label = "All Patients", vjust = -0.5) +
            annotate("text", x = 1, y = cutpoints$afr_low, label = "AFR Low", vjust = -0.5, color = fill_AFR_low) +
            annotate("text", x = 2, y = cutpoints$afr_high, label = "AFR High", vjust = -0.5, color = fill_AFR_high)
    }


    results$boxplot <- p_boxplot

    if (cutpoint_method %in% c("median", "optimal")) {
        # Create KM curves
        # 1. KM by gene expression groups
        meta_data$expr_group_all <- ifelse(meta_data$gene_expr > cutpoints$all, "High", "Low")
        meta_data$expr_group_all <- factor(meta_data$expr_group_all, levels = c("Low", "High"))


        surv_fit_expr <- surv_fit(Surv(ttcpfs, censpfs) ~ expr_group_all, data = meta_data)


        p_km_expr <- ggsurvplot(surv_fit_expr,
            conf.int = FALSE,
            linetype = c(1, 2), # 1 = solid, 2 = dashed
            palette = c("black", "red"), # Order matches factor levels: Low, High
            title = paste("KM by", analysis_name, "Expression"),
            xlab = "Time (Days)",
            ylab = "Progression-Free Survival",
            legend.title = "Expression",
            legend.labs = c("Low", "High"),
            ggtheme = theme_classic()
        )$plot

        # 2. KM by AFR ancestry and expression interaction
        meta_data$combined_group <- paste(meta_data$AFR_group, meta_data$expr_group_all, sep = ", ")
        meta_data$combined_group <- factor(meta_data$combined_group,
            levels = c(
                "AFR Low, Low", "AFR Low, High",
                "AFR High, Low", "AFR High, High"
            )
        )

        surv_fit_combined <- surv_fit(Surv(ttcpfs, censpfs) ~ combined_group, data = meta_data)

        p_km_afr <- ggsurvplot(surv_fit_combined,
            conf.int = FALSE,
            palette = c(
                fill_AFR_low,
                color_AFR_low,
                fill_AFR_high,
                color_AFR_high
            ),
            linetype = c(1, 2, 1, 2), # solid for Low expr, dashed for High expr
            title = paste("KM by AFR Ancestry and", analysis_name, "Expression"),
            xlab = "Time (Days)",
            ylab = "Progression-Free Survival",
            legend.title = "Group",
            legend.labs = c(
                "AFR Low\nExpr Low", "AFR Low\nExpr High",
                "AFR High\nExpr Low", "AFR High\nExpr High"
            ),
            ggtheme = theme_classic()
        )$plot

        results$km_expression <- p_km_expr
        results$km_ancestry <- p_km_afr

        # 3. KM plots using AFR Low cutpoint
        meta_data$expr_group_afr_low <- ifelse(meta_data$gene_expr > cutpoints$afr_low, "High", "Low")
        meta_data$expr_group_afr_low <- factor(meta_data$expr_group_afr_low, levels = c("Low", "High"))

        surv_fit_expr_afr_low <- surv_fit(Surv(ttcpfs, censpfs) ~ expr_group_afr_low, data = meta_data)

        p_km_expr_afr_low <- ggsurvplot(surv_fit_expr_afr_low,
            conf.int = FALSE,
            linetype = c(1, 2), # 1 = solid, 2 = dashed
            palette = c("black", "red"), # Order matches factor levels: Low, High
            title = paste("KM by", analysis_name, "Expression (AFR Low Cutpoint)"),
            xlab = "Time (Days)",
            ylab = "Progression-Free Survival",
            legend.title = "Expression",
            legend.labs = c("Low", "High"),
            ggtheme = theme_classic()
        )$plot

        # KM by AFR ancestry and expression interaction (AFR Low cutpoint)
        meta_data$combined_group_afr_low <- paste(meta_data$AFR_group, meta_data$expr_group_afr_low, sep = ", ")
        meta_data$combined_group_afr_low <- factor(meta_data$combined_group_afr_low,
            levels = c(
                "AFR Low, Low", "AFR Low, High",
                "AFR High, Low", "AFR High, High"
            )
        )

        surv_fit_combined_afr_low <- surv_fit(Surv(ttcpfs, censpfs) ~ combined_group_afr_low, data = meta_data)

        p_km_afr_combined_afr_low <- ggsurvplot(surv_fit_combined_afr_low,
            conf.int = FALSE,
            palette = c(
                fill_AFR_low,
                color_AFR_low,
                fill_AFR_high,
                color_AFR_high
            ),
            linetype = c(1, 2, 1, 2), # solid for Low expr, dashed for High expr
            title = paste("KM by AFR Ancestry and", analysis_name, "Expression (AFR Low Cutpoint)"),
            xlab = "Time (Days)",
            ylab = "Progression-Free Survival",
            legend.title = "Group",
            legend.labs = c(
                "AFR Low\nExpr Low",
                "AFR Low\nExpr High",
                "AFR High\nExpr Low",
                "AFR High\nExpr High"
            ),
            ggtheme = theme_classic()
        )$plot

        # 4. KM plots using AFR High cutpoint
        meta_data$expr_group_afr_high <- ifelse(meta_data$gene_expr > cutpoints$afr_high, "High", "Low")
        meta_data$expr_group_afr_high <- factor(meta_data$expr_group_afr_high, levels = c("Low", "High"))

        surv_fit_expr_afr_high <- surv_fit(Surv(ttcpfs, censpfs) ~ expr_group_afr_high, data = meta_data)

        p_km_expr_afr_high <- ggsurvplot(surv_fit_expr_afr_high,
            conf.int = FALSE,
            linetype = c(1, 2), # 1 = solid, 2 = dashed
            palette = c("black", "red"), # Order matches factor levels: Low, High
            title = paste("KM by", analysis_name, "Expression (AFR High Cutpoint)"),
            xlab = "Time (Days)",
            ylab = "Progression-Free Survival",
            legend.title = "Expression",
            legend.labs = c("Low", "High"),
            ggtheme = theme_classic()
        )$plot

        # KM by AFR ancestry and expression interaction (AFR High cutpoint)
        meta_data$combined_group_afr_high <- paste(meta_data$AFR_group, meta_data$expr_group_afr_high, sep = ", ")
        meta_data$combined_group_afr_high <- factor(meta_data$combined_group_afr_high,
            levels = c(
                "AFR Low, Low", "AFR Low, High",
                "AFR High, Low", "AFR High, High"
            )
        )

        surv_fit_combined_afr_high <- surv_fit(Surv(ttcpfs, censpfs) ~ combined_group_afr_high, data = meta_data)

        p_km_afr_combined_afr_high <- ggsurvplot(surv_fit_combined_afr_high,
            conf.int = FALSE,
            palette = c(
                fill_AFR_low,
                color_AFR_low,
                fill_AFR_high,
                color_AFR_high
            ),
            linetype = c(1, 2, 1, 2), # solid for Low expr, dashed for High expr
            title = paste("KM by AFR Ancestry and", analysis_name, "Expression (AFR High Cutpoint)"),
            xlab = "Time (Days)",
            ylab = "Progression-Free Survival",
            legend.title = "Group",
            legend.labs = c(
                "AFR Low\nExpr Low", "AFR Low\nExpr High",
                "AFR High\nExpr Low", "AFR High\nExpr High"
            ),
            ggtheme = theme_classic()
        )$plot

        # Store additional plots in results
        results$km_expression_afr_low <- p_km_expr_afr_low
        results$km_ancestry_afr_low <- p_km_afr_combined_afr_low
        results$km_expression_afr_high <- p_km_expr_afr_high
        results$km_ancestry_afr_high <- p_km_afr_combined_afr_high

        # Create HR vs cutpoint plot if optimal method
        if (cutpoint_method == "optimal") {
            # Combine cutpoint results for plotting, only include groups with valid results
            hr_plot_data <- data.frame()

            if (nrow(results$cutpoint_analysis$all) > 0) {
                hr_plot_data <- rbind(
                    hr_plot_data,
                    data.frame(results$cutpoint_analysis$all, group = "All Patients")
                )
            }

            if (nrow(results$cutpoint_analysis$afr_low) > 0) {
                hr_plot_data <- rbind(
                    hr_plot_data,
                    data.frame(results$cutpoint_analysis$afr_low, group = "AFR Low")
                )
            }

            if (nrow(results$cutpoint_analysis$afr_high) > 0) {
                hr_plot_data <- rbind(
                    hr_plot_data,
                    data.frame(results$cutpoint_analysis$afr_high, group = "AFR High")
                )
            }

            if (nrow(hr_plot_data) > 0) {
                p_hr_cutpoint <- ggplot(hr_plot_data, aes(x = cutpoint, y = hr, color = group)) +
                    geom_line(size = 1) +
                    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = group), alpha = 0.2) +
                    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
                    scale_color_manual(values = c("All Patients" = "black", "AFR Low" = fill_AFR_low, "AFR High" = fill_AFR_high)) +
                    scale_fill_manual(values = c("All Patients" = "black", "AFR Low" = fill_AFR_low, "AFR High" = fill_AFR_high)) +
                    scale_y_log10() +
                    labs(
                        title = paste("Hazard Ratio vs Cutpoint:", analysis_name),
                        x = "Expression Cutpoint",
                        y = "Hazard Ratio (95% CI) [Log Scale]",
                        color = "Group",
                        fill = "Group"
                    ) +
                    theme_classic() +
                    theme(legend.position = "bottom")

                # Add vertical lines for optimal cutpoints
                p_hr_cutpoint <- p_hr_cutpoint +
                    geom_vline(xintercept = cutpoints$all, linetype = "dashed", color = "black")

                # Collect cutpoint annotation data
                annotation_data <- data.frame()

                # Get y-axis range for positioning annotations at the top
                y_max <- max(hr_plot_data$ci_upper, na.rm = TRUE) * 2 # Position at top with some padding

                # Calculate x-axis range for padding
                x_min <- min(hr_plot_data$cutpoint, na.rm = TRUE)
                x_max <- max(hr_plot_data$cutpoint, na.rm = TRUE)
                x_range <- x_max - x_min
                x_padding <- x_range * 0.15 # 15% padding on each side

                # Collect all available cutpoints for collision analysis
                available_cutpoints <- list()
                available_cutpoints[["all"]] <- cutpoints$all
                if (!is.na(cutpoints$afr_low)) {
                    available_cutpoints[["afr_low"]] <- cutpoints$afr_low
                }
                if (!is.na(cutpoints$afr_high)) {
                    available_cutpoints[["afr_high"]] <- cutpoints$afr_high
                }

                # Determine y-positions based on collision analysis
                y_positions <- list()

                if (length(available_cutpoints) >= 2) {
                    # Check if cutpoints are within 10% x-range window
                    cutpoint_values <- unlist(available_cutpoints)
                    cutpoint_names <- names(available_cutpoints)

                    # Find which cutpoints are within 10% of each other
                    collision_groups <- list()
                    used_indices <- c()

                    for (i in seq_along(cutpoint_values)) {
                        if (i %in% used_indices) next

                        current_group <- i
                        for (j in seq_along(cutpoint_values)) {
                            if (j <= i || j %in% used_indices) next
                            if (abs(cutpoint_values[i] - cutpoint_values[j]) < (x_range * 0.1)) {
                                current_group <- c(current_group, j)
                            }
                        }
                        collision_groups[[length(collision_groups) + 1]] <- current_group
                        used_indices <- c(used_indices, current_group)
                    }

                    # Assign y-positions based on collision groups
                    for (group in collision_groups) {
                        if (length(group) == 1) {
                            # No collision, use y_max
                            y_positions[[cutpoint_names[group]]] <- y_max
                        } else {
                            # Sort by cutpoint value (ascending)
                            group_sorted <- group[order(cutpoint_values[group])]

                            if (length(group_sorted) == 2) {
                                # Two cutpoints colliding: lower at y_max, higher at 0.8*y_max
                                y_positions[[cutpoint_names[group_sorted[1]]]] <- y_max
                                y_positions[[cutpoint_names[group_sorted[2]]]] <- 10^(log10(y_max) * 0.8)
                            } else if (length(group_sorted) == 3) {
                                # Three cutpoints colliding: lowest at y_max, middle at 0.8*y_max, highest at 0.6*y_max
                                y_positions[[cutpoint_names[group_sorted[1]]]] <- y_max
                                y_positions[[cutpoint_names[group_sorted[2]]]] <- 10^(log10(y_max) * 0.8)
                                y_positions[[cutpoint_names[group_sorted[3]]]] <- 10^(log10(y_max) * 0.6)
                            }
                        }
                    }
                } else {
                    # Only one cutpoint, place at y_max
                    y_positions[["all"]] <- y_max
                }

                # Add annotation for All Patients cutpoint
                all_cutpoint_idx <- which.min(abs(results$cutpoint_analysis$all$cutpoint - cutpoints$all))
                if (length(all_cutpoint_idx) > 0) {
                    all_hr <- results$cutpoint_analysis$all$hr[all_cutpoint_idx]
                    all_pval <- results$cutpoint_analysis$all$pvalue[all_cutpoint_idx]
                    all_label <- paste0("All: HR=", round(all_hr, 2), "\np=", format(all_pval, digits = 3, scientific = TRUE))

                    annotation_data <- rbind(annotation_data, data.frame(
                        x = cutpoints$all,
                        y = y_positions[["all"]],
                        label = all_label,
                        color = "black",
                        stringsAsFactors = FALSE
                    ))
                }

                if (!is.na(cutpoints$afr_low)) {
                    p_hr_cutpoint <- p_hr_cutpoint +
                        geom_vline(xintercept = cutpoints$afr_low, linetype = "dashed", color = fill_AFR_low)

                    # Add annotation for AFR Low cutpoint
                    if (nrow(results$cutpoint_analysis$afr_low) > 0) {
                        afr_low_cutpoint_idx <- which.min(abs(results$cutpoint_analysis$afr_low$cutpoint - cutpoints$afr_low))
                        if (length(afr_low_cutpoint_idx) > 0) {
                            afr_low_hr <- results$cutpoint_analysis$afr_low$hr[afr_low_cutpoint_idx]
                            afr_low_pval <- results$cutpoint_analysis$afr_low$pvalue[afr_low_cutpoint_idx]
                            afr_low_label <- paste0("AFR Low: HR=", round(afr_low_hr, 2), "\np=", format(afr_low_pval, digits = 3, scientific = TRUE))

                            annotation_data <- rbind(annotation_data, data.frame(
                                x = cutpoints$afr_low,
                                y = y_positions[["afr_low"]],
                                label = afr_low_label,
                                color = fill_AFR_low,
                                stringsAsFactors = FALSE
                            ))
                        }
                    }
                }

                if (!is.na(cutpoints$afr_high)) {
                    p_hr_cutpoint <- p_hr_cutpoint +
                        geom_vline(xintercept = cutpoints$afr_high, linetype = "dashed", color = fill_AFR_high)

                    # Add annotation for AFR High cutpoint
                    if (nrow(results$cutpoint_analysis$afr_high) > 0) {
                        afr_high_cutpoint_idx <- which.min(abs(results$cutpoint_analysis$afr_high$cutpoint - cutpoints$afr_high))
                        if (length(afr_high_cutpoint_idx) > 0) {
                            afr_high_hr <- results$cutpoint_analysis$afr_high$hr[afr_high_cutpoint_idx]
                            afr_high_pval <- results$cutpoint_analysis$afr_high$pvalue[afr_high_cutpoint_idx]
                            afr_high_label <- paste0("AFR High: HR=", round(afr_high_hr, 2), "\np=", format(afr_high_pval, digits = 3, scientific = TRUE))

                            annotation_data <- rbind(annotation_data, data.frame(
                                x = cutpoints$afr_high,
                                y = y_positions[["afr_high"]],
                                label = afr_high_label,
                                color = fill_AFR_high,
                                stringsAsFactors = FALSE
                            ))
                        }
                    }
                }

                # Add annotations to the plot
                if (nrow(annotation_data) > 0) {
                    p_hr_cutpoint <- p_hr_cutpoint +
                        geom_text(
                            data = annotation_data,
                            aes(x = x, y = y, label = label),
                            inherit.aes = FALSE,
                            color = annotation_data$color,
                            hjust = 0, vjust = 1, # Set vjust = 1 for top anchoring
                            size = 3,
                            fontface = "bold"
                        ) +
                        # Add axis limits with padding for labels
                        coord_cartesian(
                            xlim = c(x_min - x_padding, x_max + x_padding),
                            ylim = c(
                                min(hr_plot_data$ci_lower, na.rm = TRUE) * 0.8,
                                y_max * 1.2
                            ), # Extra padding at top for annotations
                            expand = FALSE
                        )
                }

                results$hr_cutpoint_plot <- p_hr_cutpoint
            } else {
                warning("No valid cutpoint analysis results for HR plot")
                results$hr_cutpoint_plot <- NULL
            }

            # Create forest plots for different scenarios
            create_forest_plot <- function(data, cutpoint_value, cutpoint_name, subset_name, covars = NULL) {
                tryCatch(
                    {
                        # Create expression groups based on cutpoint
                        data$expr_group <- ifelse(data$gene_expr > cutpoint_value, "High", "Low")
                        data$expr_group <- factor(data$expr_group, levels = c("Low", "High"))

                        # Build Cox model formula
                        if (is.null(covars)) {
                            formula_str <- "Surv(ttcpfs, censpfs) ~ expr_group"
                        } else {
                            covar_str <- paste(covars, collapse = " + ")
                            formula_str <- paste("Surv(ttcpfs, censpfs) ~ expr_group +", covar_str)
                        }

                        # Fit Cox model
                        cox_model <- coxph(as.formula(formula_str), data = data)
                        cox_summary <- summary(cox_model)

                        # Extract all coefficients
                        coef_names <- rownames(cox_summary$coefficients)
                        hrs <- cox_summary$coefficients[, "exp(coef)"]
                        ci_lowers <- cox_summary$conf.int[, "lower .95"]
                        ci_uppers <- cox_summary$conf.int[, "upper .95"]
                        pvals <- cox_summary$coefficients[, "Pr(>|z|)"]

                        # Create forest plot data for all variables
                        forest_data <- data.frame(
                            Variable = coef_names,
                            HR = hrs,
                            CI_lower = ci_lowers,
                            CI_upper = ci_uppers,
                            P_value = pvals,
                            stringsAsFactors = FALSE
                        )

                        # Add row numbers for plotting order
                        forest_data$row_num <- rev(seq_len(nrow(forest_data)))

                        # Create forest plot
                        p_forest <- ggplot(forest_data, aes(x = HR, y = row_num)) +
                            geom_point(size = 3, color = "black") +
                            geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, color = "black") +
                            geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
                            geom_text(aes(label = paste0("HR: ", round(HR, 2), ", p=", format(P_value, digits = 3, scientific = TRUE))),
                                x = max(forest_data$CI_upper) * 1.1, hjust = 0, size = 3
                            ) +
                            scale_x_log10() +
                            scale_y_continuous(breaks = forest_data$row_num, labels = forest_data$Variable) +
                            coord_cartesian(xlim = c(min(forest_data$CI_lower) * 0.8, max(forest_data$CI_upper) * 1.8), expand = FALSE) +
                            labs(
                                title = paste("Forest Plot:", analysis_name),
                                subtitle = paste("Cutpoint:", cutpoint_name, "| Applied to:", subset_name, "(n =", nrow(data), ")"),
                                x = "Hazard Ratio (95% CI) [Log Scale]",
                                y = "Variables"
                            ) +
                            theme_classic() +
                            theme(
                                axis.text.y = element_text(size = 9, hjust = 1),
                                plot.subtitle = element_text(size = 9),
                                plot.margin = margin(5.5, 60, 5.5, 5.5, "pt")
                            )

                        return(p_forest)
                    },
                    error = function(e) {
                        warning(paste("Failed to create forest plot for", subset_name, "with", cutpoint_name, "cutpoint:", e$message))
                        return(NULL)
                    }
                )
            }

            # Create forest plot for continuous z-score model
            create_zscore_forest_plot <- function(data, subset_name, covars = NULL) {
                tryCatch(
                    {
                        # Calculate z-score across all patients (done once)
                        data$gene_expr_zscore <- scale(data$gene_expr)[, 1]

                        # Build Cox model formula with z-score
                        if (is.null(covars)) {
                            formula_str <- "Surv(ttcpfs, censpfs) ~ gene_expr_zscore"
                        } else {
                            covar_str <- paste(covars, collapse = " + ")
                            formula_str <- paste("Surv(ttcpfs, censpfs) ~ gene_expr_zscore +", covar_str)
                        }

                        # Fit Cox model
                        cox_model <- coxph(as.formula(formula_str), data = data)
                        cox_summary <- summary(cox_model)

                        # Extract all coefficients
                        coef_names <- rownames(cox_summary$coefficients)
                        hrs <- cox_summary$coefficients[, "exp(coef)"]
                        ci_lowers <- cox_summary$conf.int[, "lower .95"]
                        ci_uppers <- cox_summary$conf.int[, "upper .95"]
                        pvals <- cox_summary$coefficients[, "Pr(>|z|)"]

                        # Create forest plot data for all variables
                        forest_data <- data.frame(
                            Variable = coef_names,
                            HR = hrs,
                            CI_lower = ci_lowers,
                            CI_upper = ci_uppers,
                            P_value = pvals,
                            stringsAsFactors = FALSE
                        )

                        # Add row numbers for plotting order
                        forest_data$row_num <- rev(seq_len(nrow(forest_data)))

                        # Create forest plot
                        p_forest <- ggplot(forest_data, aes(x = HR, y = row_num)) +
                            geom_point(size = 3, color = "blue") +
                            geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, color = "blue") +
                            geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
                            geom_text(aes(label = paste0("HR: ", round(HR, 2), ", p=", format(P_value, digits = 3, scientific = TRUE))),
                                x = max(forest_data$CI_upper) * 1.1, hjust = 0, size = 3
                            ) +
                            scale_x_log10() +
                            scale_y_continuous(breaks = forest_data$row_num, labels = forest_data$Variable) +
                            coord_cartesian(xlim = c(min(forest_data$CI_lower) * 0.8, max(forest_data$CI_upper) * 1.8), expand = FALSE) +
                            labs(
                                title = paste("Forest Plot (Z-Score):", analysis_name),
                                subtitle = paste("Continuous Expression Z-Score | Applied to:", subset_name, "(n =", nrow(data), ")"),
                                x = "Hazard Ratio (95% CI) [Log Scale]",
                                y = "Variables"
                            ) +
                            theme_classic() +
                            theme(
                                axis.text.y = element_text(size = 9, hjust = 1),
                                plot.subtitle = element_text(size = 9),
                                plot.margin = margin(5.5, 60, 5.5, 5.5, "pt")
                            )

                        return(p_forest)
                    },
                    error = function(e) {
                        warning(paste("Failed to create z-score forest plot for", subset_name, ":", e$message))
                        return(NULL)
                    }
                )
            }

            # Initialize forest plot results
            results$forest_plots <- list()

            # Prepare data subsets
            all_data <- meta_data
            afr_low_data <- meta_data[meta_data$AFR_group == "AFR Low", ]
            afr_high_data <- meta_data[meta_data$AFR_group == "AFR High", ]

            # 1) Forest plots using cutpoint computed across all patients
            if (!is.na(cutpoints$all)) {
                results$forest_plots$all_cutpoint_all_patients <- create_forest_plot(
                    all_data, cutpoints$all, "All Patients Cutpoint", "All Patients", covariates
                )
                results$forest_plots$all_cutpoint_afr_low <- create_forest_plot(
                    afr_low_data, cutpoints$all, "All Patients Cutpoint", "AFR Low", covariates
                )
                results$forest_plots$all_cutpoint_afr_high <- create_forest_plot(
                    afr_high_data, cutpoints$all, "All Patients Cutpoint", "AFR High", covariates
                )
            }

            # 2) Forest plots using cutpoint computed across AFR High patients
            if (!is.na(cutpoints$afr_high)) {
                results$forest_plots$afr_high_cutpoint_all_patients <- create_forest_plot(
                    all_data, cutpoints$afr_high, "AFR High Cutpoint", "All Patients", covariates
                )
                results$forest_plots$afr_high_cutpoint_afr_low <- create_forest_plot(
                    afr_low_data, cutpoints$afr_high, "AFR High Cutpoint", "AFR Low", covariates
                )
                results$forest_plots$afr_high_cutpoint_afr_high <- create_forest_plot(
                    afr_high_data, cutpoints$afr_high, "AFR High Cutpoint", "AFR High", covariates
                )
            }

            # 3) Forest plots using cutpoint computed across AFR Low patients
            if (!is.na(cutpoints$afr_low)) {
                results$forest_plots$afr_low_cutpoint_all_patients <- create_forest_plot(
                    all_data, cutpoints$afr_low, "AFR Low Cutpoint", "All Patients", covariates
                )
                results$forest_plots$afr_low_cutpoint_afr_low <- create_forest_plot(
                    afr_low_data, cutpoints$afr_low, "AFR Low Cutpoint", "AFR Low", covariates
                )
                results$forest_plots$afr_low_cutpoint_afr_high <- create_forest_plot(
                    afr_high_data, cutpoints$afr_low, "AFR Low Cutpoint", "AFR High", covariates
                )
            }

            # 4) Forest plots using continuous z-score expression
            results$forest_plots$zscore_all_patients <- create_zscore_forest_plot(
                all_data, "All Patients", covariates
            )
            results$forest_plots$zscore_afr_low <- create_zscore_forest_plot(
                afr_low_data, "AFR Low", covariates
            )
            results$forest_plots$zscore_afr_high <- create_zscore_forest_plot(
                afr_high_data, "AFR High", covariates
            )
        }
    } else if (cutpoint_method == "categorical") {
        # NYI - Need to consider how to best handle categorical covariates generically
    }
    return(results)
}
