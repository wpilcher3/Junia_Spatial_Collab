# Create volcano plot for Cox regression results
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_cox_volcano(results_df, title="Cox Regression Volcano Plot", p_val_column='p_value',
                     p_threshold=0.05, log_hr_threshold=0.2, top_n_labels=50):
    """
    Create a volcano plot for Cox regression results

    Parameters:
    - results_df: DataFrame with Cox regression results
    - title: Plot title
    - p_threshold: P-value threshold for horizontal line
    - log_hr_threshold: Log hazard ratio threshold for vertical lines (default ±0.2)
    - top_n_labels: Number of top genes to label
    """
    import matplotlib.pyplot as plt
    import numpy as np

    # Try to import adjustText for label repulsion
    try:
        from adjustText import adjust_text
        use_adjust_text = True
    except ImportError:
        print("Note: Install adjustText package for better label positioning: pip install adjustText")
        use_adjust_text = False

    # Filter for successful results
    plot_df = results_df[results_df['status'] == 'success'].copy()

    if len(plot_df) == 0:
        print("No successful results to plot")
        return

    # Calculate log hazard ratio and -log10(p-value)
    plot_df['log_hr'] = np.log(plot_df['hazard_ratio'])
    plot_df['neg_log10_p'] = -np.log10(plot_df[p_val_column])

    # Define significance criteria
    plot_df['significant_p'] = plot_df[p_val_column] < p_threshold
    plot_df['significant_hr'] = np.abs(plot_df['log_hr']) >= log_hr_threshold
    plot_df['significant'] = plot_df['significant_p'] & plot_df['significant_hr']

    # Create color coding based on significance
    colors = []
    for _, row in plot_df.iterrows():
        if not row['significant_p']:
            colors.append('lightgrey')  # Not significant p-value
        elif not row['significant_hr']:
            colors.append('grey')       # Significant p but small effect size
        elif row['log_hr'] < -log_hr_threshold:
            colors.append('blue')       # Significant protective (HR < 1)
        elif row['log_hr'] > log_hr_threshold:
            colors.append('red')        # Significant risk (HR > 1)
        else:
            colors.append('grey')       # Should not happen with current logic

    # Create the plot
    plt.figure(figsize=(12, 8))

    # Scatter plot
    scatter = plt.scatter(plot_df['log_hr'], plot_df['neg_log10_p'],
                          c=colors, alpha=0.6, s=30)

    # Add horizontal line for p-value threshold
    p_line_y = -np.log10(p_threshold)
    plt.axhline(y=p_line_y, color='gray', linestyle='--', alpha=0.7,
                label=f'p = {p_threshold}')

    # Add vertical lines for effect size thresholds
    plt.axvline(x=0, color='gray', linestyle='--', alpha=0.7, label='HR = 1')
    plt.axvline(x=log_hr_threshold, color='darkgray', linestyle=':', alpha=0.7,
                label=f'HR = {np.exp(log_hr_threshold):.2f}')
    plt.axvline(x=-log_hr_threshold, color='darkgray', linestyle=':', alpha=0.7,
                label=f'HR = {np.exp(-log_hr_threshold):.2f}')

    # Label top significant genes by p-value
    significant_genes = plot_df[plot_df['significant']]
    if len(significant_genes) > 0:
        top_genes = significant_genes.nsmallest(
            min(top_n_labels, len(significant_genes)), p_val_column)

        # Create text annotations
        texts = []
        x = []
        y = []
        for _, row in top_genes.iterrows():
            text = plt.text(row['log_hr'], row['neg_log10_p'], row['gene'],
                            fontsize=8, alpha=0.8,
                            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
            # plt.annotate(row['gene'],
            #                   (row['log_hr'], row['neg_log10_p']),
            #                   xytext=(row['log_hr'], row['neg_log10_p']), textcoords='offset points',
            #                   fontsize=8, alpha=0.8,
            #                   bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
            texts.append(text)
            x.append(row['log_hr'])
            y.append(row['neg_log10_p'])

    else:
        print(
            f"No genes meet significance criteria (p < {p_threshold} and |log HR| >= {log_hr_threshold})")

    # Customize plot
    plt.xlabel('Log Hazard Ratio', fontsize=12)
    plt.ylabel('-Log10(P-value)', fontsize=12)
    plt.title(title, fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='blue',
              label=f'Significant Protective (HR < {np.exp(-log_hr_threshold):.2f})'),
        Patch(facecolor='red',
              label=f'Significant Risk (HR > {np.exp(log_hr_threshold):.2f})'),
        Patch(facecolor='grey', label='Significant p, small effect'),
        Patch(facecolor='lightgrey', label='Not significant'),
        plt.Line2D([0], [0], color='gray', linestyle='--',
                   label=f'p = {p_threshold}'),
        plt.Line2D([0], [0], color='darkgray', linestyle=':',
                   label=f'|log HR| = {log_hr_threshold}')
    ]
    plt.legend(handles=legend_elements, loc='lower right', fontsize=9)

    # Add summary statistics as text
    n_significant_p = (plot_df[p_val_column] < p_threshold).sum()
    n_significant_full = plot_df['significant'].sum()
    n_protective = ((plot_df['log_hr'] < -log_hr_threshold)
                    & (plot_df[p_val_column] < p_threshold)).sum()
    n_risk = ((plot_df['log_hr'] > log_hr_threshold) &
              (plot_df[p_val_column] < p_threshold)).sum()

    stats_text = f'Significant p-value (< {p_threshold}): {n_significant_p}\\n'
    stats_text += f'Significant p + effect (|logHR| ≥ {log_hr_threshold}): {n_significant_full}\\n'
    stats_text += f'Protective: {n_protective}, Risk: {n_risk}'

    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes,
             verticalalignment='top', fontsize=10,
             bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.8))

    plt.tight_layout()

    # Adjust text positions to avoid overlap if adjustText is available
    if use_adjust_text and len(texts) > 0:
        try:
            adjust_text(texts,
                        arrowprops=dict(arrowstyle='-', color='black', alpha=0.5), force_points=0.2, force_text=0.2)
        except:
            print("Note: adjustText failed, using default positioning")
    plt.show()

    # Print summary
    print(f"\\nVolcano Plot Summary:")
    print(f"- Total genes plotted: {len(plot_df)}")
    print(f"- Significant p-value (< {p_threshold}): {n_significant_p}")
    print(
        f"- Significant p + effect (|logHR| ≥ {log_hr_threshold}): {n_significant_full}")
    print(
        f"- Protective effects (HR < {np.exp(-log_hr_threshold):.2f}): {n_protective}")
    print(f"- Risk effects (HR > {np.exp(log_hr_threshold):.2f}): {n_risk}")

    significant_count = len(
        significant_genes) if 'significant_genes' in locals() else 0
    labeled_count = min(
        top_n_labels, significant_count) if significant_count > 0 else 0
    print(f"- Labeled top {labeled_count} significant genes by p-value")


def fix_dtype_from_imported_seurat_object(adata, verbose=True):
    import pandas as pd
    cont_cov_col = ['d_dx_amm_age', 'd_dx_amm_bmi', 'ANCESTRY.AFR',
                    'ANCESTRY.PEL', 'ANCESTRY.EAS', 'ttcpfs', 'ttcos']

    adata.obs[cont_cov_col].apply(pd.to_numeric, errors='coerce')
    adata.obs[['d_dx_amm_age', 'd_dx_amm_bmi', 'ANCESTRY.AFR', 'ANCESTRY.PEL', 'ANCESTRY.EAS', 'ttcpfs', 'ttcos']] = pd.to_numeric(
        adata.obs[['d_dx_amm_age', 'd_dx_amm_bmi', 'ANCESTRY.AFR', 'ANCESTRY.PEL', 'ANCESTRY.EAS']], errors="coerce")


# Function to perform genome-wide Cox regression analysis
# Generalize covariates section...
def genome_wide_cox_analysis(adata, layer='scaled', covariates=['Study_Site'],
                             min_patients=20, max_genes=None, split=None, verbose=True):
    """
    Perform Cox regression across all genes using specified layer

    Parameters:
    - adata: AnnData object with survival data and gene expression
    - layer: Layer to use for gene expression ('scaled', 'normalized', etc.)
    - covariates: List of covariates to adjust for (['Study_Site'] or ['d_dx_amm_age', 'Study_Site'])
    - min_patients: Minimum number of patients required for analysis
    - max_genes: Maximum number of genes to analyze (for testing, set None for all)
    - split: How to split gene expression ('median', 'cutpoint', or None for continuous)
    - verbose: Whether to print progress

    Returns:
    - results_df: DataFrame with Cox regression results for each gene
    """
    from lifelines import CoxPHFitter
    import warnings
    import pandas as pd

    warnings.filterwarnings('ignore')

    # Prepare base survival dataframe
    survival_df = pd.DataFrame({
        'public_id': adata.obs['public_id'].values,
        'ttcpfs': pd.to_numeric(adata.obs['ttcpfs'], errors='coerce'),
        'censpfs': pd.to_numeric(adata.obs['censpfs'], errors='coerce')
    })

    survival_df = survival_df.merge(adata.obs[['public_id'] + covariates], on="public_id",
                                    how="left", suffixes=("", "_DUPLICATE"))

    # Remove rows with missing survival data
    survival_df = survival_df.dropna(subset=['ttcpfs', 'censpfs'])
    survival_df = survival_df[survival_df['ttcpfs'] > 0]

    # Get gene expression data from specified layer
    if layer and layer not in adata.layers:
        raise ValueError(
            f"Layer '{layer}' not found. Available layers: {list(adata.layers.keys())}")

    if layer:
        expression_data = adata.layers[layer]
        if hasattr(expression_data, 'toarray'):
            expression_data = expression_data.toarray()
    else:
        expression_data = adata.X
        if hasattr(expression_data, 'toarray'):
            expression_data = expression_data.toarray()

    # Filter to patients with survival data
    patient_indices = [i for i, pid in enumerate(
        adata.obs['public_id']) if pid in survival_df['public_id'].values]
    expression_data = expression_data[patient_indices, :]

    print(f"Genome-wide Cox regression analysis:")
    print(f"- Using layer: '{layer}'")
    print(f"- Covariates: {covariates}")
    print(f"- Expression split: {split if split else 'continuous'}")
    print(f"- {len(survival_df)} patients with survival data")
    print(f"- {adata.n_vars} genes to analyze")

    if len(survival_df) < min_patients:
        raise ValueError(
            f"Insufficient patients ({len(survival_df)}) for analysis (minimum: {min_patients})")

    # Prepare results storage
    results = []

    # Determine number of genes to analyze
    n_genes = min(max_genes, adata.n_vars) if max_genes else adata.n_vars
    gene_names = adata.var_names[:n_genes]

    if verbose:
        print(f"\\nAnalyzing {n_genes} genes...")

    # Initialize Cox fitter
    cph = CoxPHFitter()

    for i, gene in enumerate(gene_names):
        if verbose and (i + 1) % 1000 == 0:
            print(f"  Processed {i + 1}/{n_genes} genes...")

        try:
            # Create analysis dataframe for this gene
            gene_df = survival_df.copy()
            raw_expression = expression_data[:, i]

            # Apply splitting if specified
            if split is None:
                # Use continuous expression
                gene_df['gene_expression'] = raw_expression
            elif split == 'median':
                # Split by median
                median_expr = np.median(raw_expression)
                gene_df['gene_expression'] = (
                    raw_expression >= median_expr).astype(int)
            elif split == 'cutpoint':
                # Split by optimal cutpoint using log-rank statistic
                try:
                    from lifelines.utils import concordance_index
                    # Find optimal cutpoint by testing percentiles
                    best_p = np.inf
                    best_cutpoint = np.median(raw_expression)

                    # Test percentiles from 25th to 75th
                    percentiles = np.arange(25, 76, 5)
                    for pct in percentiles:
                        cutpoint = np.percentile(raw_expression, pct)
                        binary_expr = (raw_expression >= cutpoint).astype(int)

                        # Quick log-rank test
                        if len(np.unique(binary_expr)) == 2:  # Ensure we have both groups
                            try:
                                from lifelines.statistics import logrank_test
                                group1_mask = binary_expr == 0
                                group2_mask = binary_expr == 1

                                if group1_mask.sum() > 5 and group2_mask.sum() > 5:
                                    lr_result = logrank_test(
                                        gene_df.loc[group1_mask, 'ttcpfs'],
                                        gene_df.loc[group2_mask, 'ttcpfs'],
                                        gene_df.loc[group1_mask, 'censpfs'],
                                        gene_df.loc[group2_mask, 'censpfs']
                                    )
                                    if lr_result.p_value < best_p:
                                        best_p = lr_result.p_value
                                        best_cutpoint = cutpoint
                            except:
                                continue

                    gene_df['gene_expression'] = (
                        raw_expression >= best_cutpoint).astype(int)
                except:
                    # Fallback to median if cutpoint optimization fails
                    median_expr = np.median(raw_expression)
                    gene_df['gene_expression'] = (
                        raw_expression >= median_expr).astype(int)
            else:
                raise ValueError(
                    f"Invalid split method: {split}. Use 'median', 'cutpoint', or None.")

            # Create dummy variables for categorical covariates
            analysis_df = gene_df.copy()

            for cov in covariates:
                dtype = analysis_df[cov].dtypes
                if dtype in ('category', 'string', 'boolean'):
                    dummies = pd.get_dummies(
                        analysis_df[cov], prefix=cov, drop_first=True
                    )
                    analysis_df = pd.concat([analysis_df, dummies], axis=1)

            # Define model columns
            model_cols = ['ttcpfs', 'censpfs', 'gene_expression']

            # Add dummy covariates if categorical, else just add the continuous covariate
            for cov in covariates:
                dtype = analysis_df[cov].dtypes
                if dtype in ('category', 'string', 'boolean'):
                    col_list = [
                        col for col in analysis_df.columns if col.startswith(f'{cov}_')]
                    model_cols.extend(col_list)
                else:
                    model_cols.append(cov)

            # Select final dataset
            final_df = analysis_df[model_cols].dropna()

            if len(final_df) < min_patients:
                results.append({
                    'gene': gene,
                    'n_patients': len(final_df),
                    'n_events': np.nan,
                    'hazard_ratio': np.nan,
                    'hr_lower_ci': np.nan,
                    'hr_upper_ci': np.nan,
                    'p_value': np.nan,
                    'log10_p': np.nan,
                    'status': 'insufficient_patients'
                })
                continue

            # Fit Cox model
            cph.fit(final_df, duration_col='ttcpfs', event_col='censpfs')

            # Extract results for gene expression
            gene_coef = cph.params_['gene_expression']
            gene_hr = np.exp(gene_coef)
            gene_p = cph.summary.loc['gene_expression', 'p']

            # Get confidence intervals
            ci_lower = np.exp(
                cph.confidence_intervals_.loc['gene_expression', '95% lower-bound'])
            ci_upper = np.exp(
                cph.confidence_intervals_.loc['gene_expression', '95% upper-bound'])

            results.append({
                'gene': gene,
                'n_patients': len(final_df),
                'n_events': int(final_df['censpfs'].sum()),
                'hazard_ratio': gene_hr,
                'hr_lower_ci': ci_lower,
                'hr_upper_ci': ci_upper,
                'p_value': gene_p,
                'log10_p': -np.log10(gene_p) if gene_p > 0 else np.nan,
                'status': 'success'
            })

        except Exception as e:
            results.append({
                'gene': gene,
                'n_patients': np.nan,
                'n_events': np.nan,
                'hazard_ratio': np.nan,
                'hr_lower_ci': np.nan,
                'hr_upper_ci': np.nan,
                'p_value': np.nan,
                'log10_p': np.nan,
                'status': f'error: {str(e)[:50]}'
            })

    # Convert to DataFrame
    results_df = pd.DataFrame(results)

    # Add multiple testing correction
    valid_p = results_df['p_value'].dropna()
    if len(valid_p) > 0:
        try:
            # Benjamini-Hochberg FDR correction
            from statsmodels.stats.multitest import multipletests
            _, fdr_p, _, _ = multipletests(
                valid_p, alpha=0.05, method='fdr_bh')

            # Add FDR-corrected p-values
            results_df.loc[results_df['p_value'].notna(),
                           'fdr_p_value'] = fdr_p
            results_df.loc[results_df['fdr_p_value'].notna(), 'log10_fdr_p'] = -np.log10(
                results_df.loc[results_df['fdr_p_value'].notna(), 'fdr_p_value'])
        except:
            print("Warning: Could not calculate FDR correction")

    if verbose:
        successful = (results_df['status'] == 'success').sum()
        print(f"\\nAnalysis complete:")
        print(f"- {successful}/{len(results_df)} genes analyzed successfully")
        if successful > 0:
            sig_raw = (results_df['p_value'] < 0.05).sum()
            print(f"- {sig_raw} genes with p < 0.05")
            if 'fdr_p_value' in results_df.columns:
                sig_fdr = (results_df['fdr_p_value'] < 0.05).sum()
                print(f"- {sig_fdr} genes with FDR < 0.05")

    return results_df


# Function to perform genome-wide Cox regression analysis with interaction term
def genome_wide_cox_analysis_with_interaction(adata, interaction_covariate='ANCESTRY.AFR', layer='scaled',
                                              covariates=['Study_Site'],
                                              min_patients=20, max_genes=None, split=None,
                                              verbose=True):
    """
    Perform Cox regression across all genes with interaction term between gene expression and a covariate

    Parameters:
    - adata: AnnData object with survival data and gene expression
    - interaction_covariate: Name of the covariate to create interaction term with gene expression
    - layer: Layer to use for gene expression ('scaled', 'normalized', etc.)
    - covariates: List of covariates to adjust for (should include interaction_covariate)
    - min_patients: Minimum number of patients required for analysis
    - max_genes: Maximum number of genes to analyze (for testing, set None for all)
    - split: How to split gene expression ('median', 'cutpoint', or None for continuous)
    - verbose: Whether to print progress

    Returns:
    - results_df: DataFrame with Cox regression results including interaction terms for each gene
    """
    from lifelines import CoxPHFitter
    import warnings
    import pandas as pd

    warnings.filterwarnings('ignore')

    # Ensure interaction covariate is in the covariates list
    if interaction_covariate not in covariates:
        covariates = covariates + [interaction_covariate]

    # Prepare base survival dataframe
    survival_df = pd.DataFrame({
        'public_id': adata.obs['public_id'].values,
        'ttcpfs': pd.to_numeric(adata.obs['ttcpfs'], errors='coerce'),
        'censpfs': pd.to_numeric(adata.obs['censpfs'], errors='coerce')
    })

    survival_df = survival_df.merge(adata.obs[['public_id'] + covariates], on="public_id",
                                    how="left", suffixes=("", "_DUPLICATE"))

    # Remove rows with missing survival data
    survival_df = survival_df.dropna(subset=['ttcpfs', 'censpfs'])
    survival_df = survival_df[survival_df['ttcpfs'] > 0]

    # Get gene expression data from specified layer
    if layer and layer not in adata.layers:
        raise ValueError(
            f"Layer '{layer}' not found. Available layers: {list(adata.layers.keys())}")

    if layer:
        expression_data = adata.layers[layer]
        if hasattr(expression_data, 'toarray'):
            expression_data = expression_data.toarray()
    else:
        expression_data = adata.X
        if hasattr(expression_data, 'toarray'):
            expression_data = expression_data.toarray()

    # Filter to patients with survival data
    patient_indices = [i for i, pid in enumerate(
        adata.obs['public_id']) if pid in survival_df['public_id'].values]
    expression_data = expression_data[patient_indices, :]

    print(f"Genome-wide Cox regression analysis with interaction:")
    print(f"- Using layer: '{layer}'")
    print(f"- Interaction covariate: {interaction_covariate}")
    print(f"- Covariates: {covariates}")
    print(f"- Expression split: {split if split else 'continuous'}")
    print(f"- {len(survival_df)} patients with survival data")
    print(f"- {adata.n_vars} genes to analyze")

    if len(survival_df) < min_patients:
        raise ValueError(
            f"Insufficient patients ({len(survival_df)}) for analysis (minimum: {min_patients})")

    # Prepare results storage
    results = []

    # Determine number of genes to analyze
    n_genes = min(max_genes, adata.n_vars) if max_genes else adata.n_vars
    gene_names = adata.var_names[:n_genes]

    if verbose:
        print(f"\\nAnalyzing {n_genes} genes...")

    # Initialize Cox fitter
    cph = CoxPHFitter()

    for i, gene in enumerate(gene_names):
        if verbose and (i + 1) % 1000 == 0:
            print(f"  Processed {i + 1}/{n_genes} genes...")

        try:
            # Create analysis dataframe for this gene
            gene_df = survival_df.copy()
            raw_expression = expression_data[:, i]

            # Apply splitting if specified
            if split is None:
                # Use continuous expression
                gene_df['gene_expression'] = raw_expression
            elif split == 'median':
                # Split by median
                median_expr = np.median(raw_expression)
                gene_df['gene_expression'] = (
                    raw_expression >= median_expr).astype(int)
            elif split == 'cutpoint':
                # Split by optimal cutpoint using log-rank statistic
                try:
                    from lifelines.utils import concordance_index
                    # Find optimal cutpoint by testing percentiles
                    best_p = np.inf
                    best_cutpoint = np.median(raw_expression)

                    # Test percentiles from 25th to 75th
                    percentiles = np.arange(25, 76, 5)
                    for pct in percentiles:
                        cutpoint = np.percentile(raw_expression, pct)
                        binary_expr = (raw_expression >= cutpoint).astype(int)

                        # Quick log-rank test
                        if len(np.unique(binary_expr)) == 2:  # Ensure we have both groups
                            try:
                                from lifelines.statistics import logrank_test
                                group1_mask = binary_expr == 0
                                group2_mask = binary_expr == 1

                                if group1_mask.sum() > 5 and group2_mask.sum() > 5:
                                    lr_result = logrank_test(
                                        gene_df.loc[group1_mask, 'ttcpfs'],
                                        gene_df.loc[group2_mask, 'ttcpfs'],
                                        gene_df.loc[group1_mask, 'censpfs'],
                                        gene_df.loc[group2_mask, 'censpfs']
                                    )
                                    if lr_result.p_value < best_p:
                                        best_p = lr_result.p_value
                                        best_cutpoint = cutpoint
                            except:
                                continue

                    gene_df['gene_expression'] = (
                        raw_expression >= best_cutpoint).astype(int)
                except:
                    # Fallback to median if cutpoint optimization fails
                    median_expr = np.median(raw_expression)
                    gene_df['gene_expression'] = (
                        raw_expression >= median_expr).astype(int)
            else:
                raise ValueError(
                    f"Invalid split method: {split}. Use 'median', 'cutpoint', or None.")

            # Create dummy variables for categorical covariates (excluding interaction covariate initially)
            analysis_df = gene_df.copy()

            # Track which covariates are categorical for proper handling
            categorical_covariates = []

            for cov in covariates:
                dtype = analysis_df[cov].dtypes
                if dtype in ('category', 'string', 'boolean'):
                    categorical_covariates.append(cov)
                    # Only create dummies for non-interaction covariates
                    if cov != interaction_covariate:
                        dummies = pd.get_dummies(
                            analysis_df[cov], prefix=cov, drop_first=True
                        )
                        analysis_df = pd.concat([analysis_df, dummies], axis=1)

            # Handle interaction covariate and create interaction term
            if interaction_covariate in categorical_covariates:
                # If categorical, create dummies for interaction covariate
                interaction_dummies = pd.get_dummies(
                    analysis_df[interaction_covariate],
                    prefix=interaction_covariate,
                    drop_first=True
                )
                analysis_df = pd.concat(
                    [analysis_df, interaction_dummies], axis=1)

                # Create interaction terms for each dummy variable
                for dummy_col in interaction_dummies.columns:
                    analysis_df[f'gene_x_{dummy_col}'] = analysis_df['gene_expression'] * \
                        analysis_df[dummy_col]
            else:
                # If continuous, create single interaction term
                analysis_df['gene_x_interaction'] = analysis_df['gene_expression'] * \
                    analysis_df[interaction_covariate]

            # Define model columns
            model_cols = ['ttcpfs', 'censpfs', 'gene_expression']

            # Add covariates (dummy or continuous)
            for cov in covariates:
                dtype = analysis_df[cov].dtypes
                if dtype in ('category', 'string', 'boolean'):
                    col_list = [
                        col for col in analysis_df.columns if col.startswith(f'{cov}_')]
                    model_cols.extend(col_list)
                else:
                    model_cols.append(cov)

            # Add interaction terms
            interaction_cols = [
                col for col in analysis_df.columns if col.startswith('gene_x_')]
            model_cols.extend(interaction_cols)

            # Select final dataset
            final_df = analysis_df[model_cols].dropna()

            if len(final_df) < min_patients:
                results.append({
                    'gene': gene,
                    'n_patients': len(final_df),
                    'n_events': np.nan,
                    'gene_hr': np.nan,
                    'gene_hr_lower_ci': np.nan,
                    'gene_hr_upper_ci': np.nan,
                    'gene_p': np.nan,
                    'interaction_hr': np.nan,
                    'interaction_hr_lower_ci': np.nan,
                    'interaction_hr_upper_ci': np.nan,
                    'interaction_p': np.nan,
                    'covariate_hr': np.nan,
                    'covariate_hr_lower_ci': np.nan,
                    'covariate_hr_upper_ci': np.nan,
                    'covariate_p': np.nan,
                    'status': 'insufficient_patients'
                })
                continue

            # Fit Cox model
            cph.fit(final_df, duration_col='ttcpfs', event_col='censpfs')

            # Extract results for gene expression
            gene_coef = cph.params_['gene_expression']
            gene_hr = np.exp(gene_coef)
            gene_p = cph.summary.loc['gene_expression', 'p']
            gene_ci_lower = np.exp(
                cph.confidence_intervals_.loc['gene_expression', '95% lower-bound'])
            gene_ci_upper = np.exp(
                cph.confidence_intervals_.loc['gene_expression', '95% upper-bound'])

            # Extract results for interaction term(s)
            # For categorical interaction covariate, we'll take the first dummy's interaction
            if len(interaction_cols) > 0:
                # Use first interaction term
                interaction_col = interaction_cols[0]
                interaction_coef = cph.params_[interaction_col]
                interaction_hr = np.exp(interaction_coef)
                interaction_p = cph.summary.loc[interaction_col, 'p']
                interaction_ci_lower = np.exp(
                    cph.confidence_intervals_.loc[interaction_col, '95% lower-bound'])
                interaction_ci_upper = np.exp(
                    cph.confidence_intervals_.loc[interaction_col, '95% upper-bound'])
            else:
                interaction_hr = np.nan
                interaction_p = np.nan
                interaction_ci_lower = np.nan
                interaction_ci_upper = np.nan

            # Extract results for the interaction covariate itself
            if interaction_covariate in categorical_covariates:
                # Get the first dummy variable for the interaction covariate
                cov_cols = [col for col in final_df.columns if col.startswith(
                    f'{interaction_covariate}_')]
                if len(cov_cols) > 0:
                    cov_col = cov_cols[0]
                    covariate_coef = cph.params_[cov_col]
                    covariate_hr = np.exp(covariate_coef)
                    covariate_p = cph.summary.loc[cov_col, 'p']
                    covariate_ci_lower = np.exp(
                        cph.confidence_intervals_.loc[cov_col, '95% lower-bound'])
                    covariate_ci_upper = np.exp(
                        cph.confidence_intervals_.loc[cov_col, '95% upper-bound'])
                else:
                    covariate_hr = np.nan
                    covariate_p = np.nan
                    covariate_ci_lower = np.nan
                    covariate_ci_upper = np.nan
            else:
                # Continuous covariate
                covariate_coef = cph.params_[interaction_covariate]
                covariate_hr = np.exp(covariate_coef)
                covariate_p = cph.summary.loc[interaction_covariate, 'p']
                covariate_ci_lower = np.exp(
                    cph.confidence_intervals_.loc[interaction_covariate, '95% lower-bound'])
                covariate_ci_upper = np.exp(
                    cph.confidence_intervals_.loc[interaction_covariate, '95% upper-bound'])

            results.append({
                'gene': gene,
                'n_patients': len(final_df),
                'n_events': int(final_df['censpfs'].sum()),
                'gene_hr': gene_hr,
                'gene_hr_lower_ci': gene_ci_lower,
                'gene_hr_upper_ci': gene_ci_upper,
                'gene_p': gene_p,
                'gene_log10_p': -np.log10(gene_p) if gene_p > 0 else np.nan,
                'interaction_hr': interaction_hr,
                'interaction_hr_lower_ci': interaction_ci_lower,
                'interaction_hr_upper_ci': interaction_ci_upper,
                'interaction_p': interaction_p,
                'interaction_log10_p': -np.log10(interaction_p) if interaction_p > 0 and not np.isnan(interaction_p) else np.nan,
                'covariate_hr': covariate_hr,
                'covariate_hr_lower_ci': covariate_ci_lower,
                'covariate_hr_upper_ci': covariate_ci_upper,
                'covariate_p': covariate_p,
                'covariate_log10_p': -np.log10(covariate_p) if covariate_p > 0 and not np.isnan(covariate_p) else np.nan,
                'status': 'success'
            })

        except Exception as e:
            results.append({
                'gene': gene,
                'n_patients': np.nan,
                'n_events': np.nan,
                'gene_hr': np.nan,
                'gene_hr_lower_ci': np.nan,
                'gene_hr_upper_ci': np.nan,
                'gene_p': np.nan,
                'gene_log10_p': np.nan,
                'interaction_hr': np.nan,
                'interaction_hr_lower_ci': np.nan,
                'interaction_hr_upper_ci': np.nan,
                'interaction_p': np.nan,
                'interaction_log10_p': np.nan,
                'covariate_hr': np.nan,
                'covariate_hr_lower_ci': np.nan,
                'covariate_hr_upper_ci': np.nan,
                'covariate_p': np.nan,
                'covariate_log10_p': np.nan,
                'status': f'error: {str(e)[:50]}'
            })

    # Convert to DataFrame
    results_df = pd.DataFrame(results)

    # Add multiple testing correction for each p-value column
    for p_col in ['gene_p', 'interaction_p', 'covariate_p']:
        valid_p = results_df[p_col].dropna()
        if len(valid_p) > 0:
            try:
                # Benjamini-Hochberg FDR correction
                from statsmodels.stats.multitest import multipletests
                _, fdr_p, _, _ = multipletests(
                    valid_p, alpha=0.05, method='fdr_bh')

                # Add FDR-corrected p-values
                fdr_col = p_col.replace('_p', '_fdr_p')
                results_df.loc[results_df[p_col].notna(), fdr_col] = fdr_p
                results_df.loc[results_df[fdr_col].notna(), fdr_col.replace('_p', '_log10_p')] = -np.log10(
                    results_df.loc[results_df[fdr_col].notna(), fdr_col])
            except:
                print(
                    f"Warning: Could not calculate FDR correction for {p_col}")

    if verbose:
        successful = (results_df['status'] == 'success').sum()
        print(f"\\nAnalysis complete:")
        print(f"- {successful}/{len(results_df)} genes analyzed successfully")
        if successful > 0:
            print(f"\\nGene expression main effect:")
            sig_raw = (results_df['gene_p'] < 0.05).sum()
            print(f"  - {sig_raw} genes with p < 0.05")
            if 'gene_fdr_p' in results_df.columns:
                sig_fdr = (results_df['gene_fdr_p'] < 0.05).sum()
                print(f"  - {sig_fdr} genes with FDR < 0.05")

            print(f"\\nInteraction term:")
            sig_raw_int = (results_df['interaction_p'] < 0.05).sum()
            print(f"  - {sig_raw_int} genes with interaction p < 0.05")
            if 'interaction_fdr_p' in results_df.columns:
                sig_fdr_int = (results_df['interaction_fdr_p'] < 0.05).sum()
                print(f"  - {sig_fdr_int} genes with interaction FDR < 0.05")

            print(f"\\nCovariate main effect:")
            sig_raw_cov = (results_df['covariate_p'] < 0.05).sum()
            print(f"  - {sig_raw_cov} genes with covariate p < 0.05")
            if 'covariate_fdr_p' in results_df.columns:
                sig_fdr_cov = (results_df['covariate_fdr_p'] < 0.05).sum()
                print(f"  - {sig_fdr_cov} genes with covariate FDR < 0.05")

    return results_df


# # Create volcano plot for Study_Site adjustment results
# print("=== Volcano Plot: Cox Regression Results (Study_Site Adjustment) ===")
# plot_cox_volcano(test_results_site,
#                 title="Cox Regression: Gene Expression vs Survival (Study_Site Adjusted)",
#                 p_threshold=0.05,
#                 top_n_labels=50)

# Enhanced scatter plot with labels for genes significant in both analyses


def create_survival_vs_deg_scatterplot_with_labels(cox_results, deg_results, cox_pval_column='p_value', deg_pval_column='padj',
                                                   cox_p_threshold=0.05, deg_p_threshold=0.05,
                                                   figsize=(14, 10), label_both_sig=True, compartment_name="NK Cells"):
    """
    Create scatter plot comparing Cox regression hazard ratios with differential expression log2FoldChange
    with labels for genes significant in both analyses

    Parameters:
    - cox_results: DataFrame with Cox regression results (must have 'gene', 'hazard_ratio', 'p_value' columns)
    - deg_results: DataFrame with DEG results (must have log2FoldChange, padj as columns, genes as index)
    - cox_p_threshold: P-value threshold for Cox significance
    - deg_p_threshold: P-value threshold for DEG significance
    - label_both_sig: Whether to add gene labels for genes significant in both analyses
    """

    # Prepare the data for merging
    cox_df = cox_results[cox_results['status'] == 'success'].copy()
    cox_df['log_hazard_ratio'] = np.log2(cox_df['hazard_ratio'])

    # Create DEG dataframe with gene names as a column
    deg_df = deg_results.copy()
    deg_df['gene'] = deg_df.index

    # Merge the datasets on gene names
    merged_df = pd.merge(cox_df[['gene', 'log_hazard_ratio', cox_pval_column, 'hazard_ratio']],
                         deg_df[['gene', 'log2FoldChange', deg_pval_column]],
                         on='gene')

    print(
        f"Merged dataset: {len(merged_df)} genes with both survival and DEG data")

    merged_df['DEG_directionality'] = 'No Difference'

    # Create significance categories
    merged_df['significance_category'] = 'Not significant'

    # Cox significant only (column names are unique, no suffixes added)
    cox_sig = merged_df[cox_pval_column] < cox_p_threshold
    deg_sig = merged_df[deg_pval_column] < deg_p_threshold

    merged_df.loc[cox_sig & ~deg_sig, 'significance_category'] = 'Cox only'
    merged_df.loc[~cox_sig & deg_sig, 'significance_category'] = 'DEG only'
    merged_df.loc[cox_sig & deg_sig,
                  'significance_category'] = 'Both significant'

    merged_df.loc[merged_df['log2FoldChange']
                  < 0, 'DEG_directionality'] = 'AFR Down'
    merged_df.loc[merged_df['log2FoldChange']
                  > 0, 'DEG_directionality'] = 'AFR Up'

    # Count categories
    category_counts = merged_df['significance_category'].value_counts()
    print("\\nSignificance categories:")
    for cat, count in category_counts.items():
        print(f"  {cat}: {count} genes")

    # Create the plot
    fig, ax = plt.subplots(figsize=figsize)

    # Define colors for categories
    colors = {
        'Not significant': '#CCCCCC',
        'Cox only': '#FF6B6B',
        'DEG only': '#4ECDC4',
        'Both significant': '#45B7D1'
    }
    colors_label_face = {
        'No Difference': '#CCCCCC',
        'AFR Up': "#FAA6A6",
        'AFR Down': "#95E6E1"
    }

    # Plot each category
    for category in colors.keys():
        if category in merged_df['significance_category'].values:
            subset = merged_df[merged_df['significance_category'] == category]
            ax.scatter(subset['log2FoldChange'], subset['log_hazard_ratio'],
                       c=colors[category], label=f'{category} (n={len(subset)})',
                       alpha=0.7, s=50, edgecolors='white', linewidth=0.5)

    # Add gene labels for genes significant in both analyses
    if label_both_sig:
        both_sig = merged_df[merged_df['significance_category']
                             == 'Both significant']
        if len(both_sig) > 0:
            # Use adjustText for automatic label positioning to avoid overlap
            try:
                from adjustText import adjust_text
                texts = []
                for _, row in both_sig.iterrows():
                    texts.append(ax.text(
                        row['log2FoldChange'], row['log_hazard_ratio'], row['gene'],
                        fontsize=9, fontweight='bold',
                        bbox=dict(boxstyle="round,pad=0.3",
                                  facecolor=colors_label_face[row['DEG_directionality']],
                                  alpha=0.8,
                                  edgecolor='black',
                                  linewidth=0.5)))

                # Adjust text positions to avoid overlaps
                adjust_text(texts, ax=ax,
                            expand_points=(1.2, 1.2),
                            expand_text=(1.2, 1.2),
                            arrowprops=dict(arrowstyle='->', color='gray', alpha=0.7, lw=0.5))

            except ImportError:
                # Fallback to simple annotation without adjustText
                print("Note: adjustText not available, using simple annotations")
                for _, row in both_sig.iterrows():
                    ax.annotate(row['gene'],
                                (row['log2FoldChange'],
                                 row['log_hazard_ratio']),
                                xytext=(5, 5), textcoords='offset points',
                                fontsize=9, fontweight='bold',
                                bbox=dict(boxstyle="round,pad=0.3",
                                          facecolor="white",
                                          alpha=0.8,
                                          edgecolor='black',
                                          linewidth=0.5))

    # Add reference lines
    ax.axhline(y=0, color='black', linestyle='--', alpha=0.5, linewidth=1)
    ax.axvline(x=0, color='black', linestyle='--', alpha=0.5, linewidth=1)

    # Labels and title
    ax.set_xlabel(
        'log2FoldChange (DEG Analysis: AFR High vs Low)', fontsize=12)
    ax.set_ylabel('log2(Hazard Ratio) (Cox Survival Analysis)', fontsize=12)
    ax.set_title(f'Survival Association vs Differential Expression\\n{compartment_name} - Genes Significant in Both Analyses Labeled',
                 fontsize=14, pad=20)

    # Add correlation
    correlation = merged_df['log2FoldChange'].corr(
        merged_df['log_hazard_ratio'])
    ax.text(0.05, 0.95, f'Correlation: r = {correlation:.3f}',
            transform=ax.transAxes, fontsize=12,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    # Legend
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    # Grid
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

    # Print detailed information about labeled genes
    both_sig = merged_df[merged_df['significance_category']
                         == 'Both significant']
    if len(both_sig) > 0:
        print(
            f"\\n=== GENES SIGNIFICANT IN BOTH ANALYSES (n={len(both_sig)}) ===")
        both_sig_detailed = both_sig.copy()
        both_sig_detailed['combined_score'] = np.abs(
            both_sig_detailed['log2FoldChange']) + np.abs(both_sig_detailed['log_hazard_ratio'])
        both_sig_detailed = both_sig_detailed.sort_values(
            'combined_score', ascending=False)

        print("\\nRanked by combined effect size (|log2FC| + |logHR|):")
        for i, (_, row) in enumerate(both_sig_detailed.iterrows(), 1):
            deg_direction = "↑" if row['log2FoldChange'] > 0 else "↓"
            surv_direction = "↑risk" if row['log_hazard_ratio'] > 0 else "↓risk"
            concordant = "✓" if ((row['log2FoldChange'] > 0) == (
                row['log_hazard_ratio'] > 0)) else "✗"

            print(f"{i:2d}. {row['gene']:<12} | DEG: {deg_direction} {row['log2FoldChange']:+.3f} (p={row[deg_pval_column]:.3e}) | "
                  f"Survival: {surv_direction} {row['log_hazard_ratio']:+.3f} (p={row[cox_pval_column]:.3e}) | "
                  f"Concordant: {concordant}")

    return merged_df


def ulm_pathway_cox_analysis(adata, ulm_obsm_key='score_ulm', covariates=['Study_Site'],
                              min_patients=20, verbose=True):
    """
    Perform Cox regression across all ULM pathway scores

    Parameters:
    - adata: AnnData object with survival data and ULM scores in obsm
    - ulm_obsm_key: Key in adata.obsm containing ULM scores (default: 'score_ulm')
    - covariates: List of covariates to adjust for. Continuous covariates will be z-scored.
    - min_patients: Minimum number of patients required for analysis
    - verbose: Whether to print progress

    Returns:
    - results_df: DataFrame with Cox regression results for each pathway
    """
    from lifelines import CoxPHFitter
    from scipy.stats import zscore
    import warnings
    import pandas as pd
    import numpy as np

    warnings.filterwarnings('ignore')

    # Get ULM scores
    if ulm_obsm_key not in adata.obsm:
        raise ValueError(f"ULM score key '{ulm_obsm_key}' not found in adata.obsm")
    
    ulm_scores = pd.DataFrame(adata.obsm[ulm_obsm_key], 
                             index=adata.obs_names,
                             columns=adata.obsm[ulm_obsm_key].columns if hasattr(adata.obsm[ulm_obsm_key], 'columns') else None)
    
    # Prepare base survival dataframe
    survival_df = pd.DataFrame({
        'public_id': adata.obs['public_id'].values,
        'ttcpfs': pd.to_numeric(adata.obs['ttcpfs'], errors='coerce'),
        'censpfs': pd.to_numeric(adata.obs['censpfs'], errors='coerce')
    })

    # Add covariates
    for cov in covariates:
        if cov in adata.obs.columns:
            survival_df[cov] = adata.obs[cov].values
        else:
            print(f"Warning: Covariate '{cov}' not found in adata.obs")

    # Remove rows with missing survival data
    survival_df = survival_df.dropna(subset=['ttcpfs', 'censpfs'])
    survival_df = survival_df[survival_df['ttcpfs'] > 0]

    if len(survival_df) < min_patients:
        raise ValueError(
            f"Insufficient patients ({len(survival_df)}) for analysis (minimum: {min_patients})")

    print(f"ULM Pathway Cox regression analysis:")
    print(f"- Covariates: {covariates}")
    print(f"- {len(survival_df)} patients with survival data")
    print(f"- {ulm_scores.shape[1]} pathways to analyze")

    # Prepare results storage
    results = []

    # Initialize Cox fitter
    cph = CoxPHFitter()

    for i, pathway in enumerate(ulm_scores.columns):
        if verbose and (i + 1) % 10 == 0:
            print(f"  Processed {i + 1}/{ulm_scores.shape[1]} pathways...")

        try:
            # Create analysis dataframe for this pathway
            pathway_df = survival_df.copy()
            
            # Get ULM scores for patients with survival data
            pathway_scores = ulm_scores.loc[survival_df['public_id'].values, pathway].values
            
            # Z-score the pathway scores
            pathway_df['pathway_score_z'] = zscore(pathway_scores, nan_policy='omit')
            
            # Process covariates: z-score continuous, dummy encode categorical
            analysis_df = pathway_df.copy()
            model_cols = ['ttcpfs', 'censpfs', 'pathway_score_z']

            for cov in covariates:
                if cov not in analysis_df.columns:
                    continue
                    
                dtype = analysis_df[cov].dtypes
                if dtype in ('category', 'string', 'boolean', 'object'):
                    # Categorical: create dummy variables
                    dummies = pd.get_dummies(analysis_df[cov], prefix=cov, drop_first=True)
                    analysis_df = pd.concat([analysis_df, dummies], axis=1)
                    col_list = [col for col in analysis_df.columns if col.startswith(f'{cov}_')]
                    model_cols.extend(col_list)
                else:
                    # Continuous: z-score and prefix with 'scaled_'
                    scaled_col_name = f'scaled_{cov}'
                    analysis_df[scaled_col_name] = zscore(analysis_df[cov].values, nan_policy='omit')
                    model_cols.append(scaled_col_name)

            # Select final dataset
            final_df = analysis_df[model_cols].dropna()

            if len(final_df) < min_patients:
                results.append({
                    'pathway': pathway,
                    'n_patients': len(final_df),
                    'n_events': np.nan,
                    'hazard_ratio': np.nan,
                    'hr_lower_ci': np.nan,
                    'hr_upper_ci': np.nan,
                    'p_value': np.nan,
                    'log10_p': np.nan,
                    'status': 'insufficient_patients'
                })
                continue

            # Fit Cox model
            cph.fit(final_df, duration_col='ttcpfs', event_col='censpfs')

            # Extract results for pathway score
            pathway_coef = cph.params_['pathway_score_z']
            pathway_hr = np.exp(pathway_coef)
            pathway_p = cph.summary.loc['pathway_score_z', 'p']

            # Get confidence intervals
            ci_lower = np.exp(cph.confidence_intervals_.loc['pathway_score_z', '95% lower-bound'])
            ci_upper = np.exp(cph.confidence_intervals_.loc['pathway_score_z', '95% upper-bound'])

            results.append({
                'pathway': pathway,
                'n_patients': len(final_df),
                'n_events': int(final_df['censpfs'].sum()),
                'hazard_ratio': pathway_hr,
                'hr_lower_ci': ci_lower,
                'hr_upper_ci': ci_upper,
                'p_value': pathway_p,
                'log10_p': -np.log10(pathway_p) if pathway_p > 0 else np.nan,
                'status': 'success'
            })

        except Exception as e:
            results.append({
                'pathway': pathway,
                'n_patients': np.nan,
                'n_events': np.nan,
                'hazard_ratio': np.nan,
                'hr_lower_ci': np.nan,
                'hr_upper_ci': np.nan,
                'p_value': np.nan,
                'log10_p': np.nan,
                'status': f'error: {str(e)[:50]}'
            })

    # Convert to DataFrame
    results_df = pd.DataFrame(results)

    # Add multiple testing correction
    valid_p = results_df['p_value'].dropna()
    if len(valid_p) > 0:
        try:
            from statsmodels.stats.multitest import multipletests
            _, fdr_p, _, _ = multipletests(valid_p, alpha=0.05, method='fdr_bh')
            results_df.loc[results_df['p_value'].notna(), 'fdr_p_value'] = fdr_p
            results_df.loc[results_df['fdr_p_value'].notna(), 'log10_fdr_p'] = -np.log10(
                results_df.loc[results_df['fdr_p_value'].notna(), 'fdr_p_value'])
        except:
            print("Warning: Could not calculate FDR correction")

    if verbose:
        successful = (results_df['status'] == 'success').sum()
        print(f"\nAnalysis complete:")
        print(f"- {successful}/{len(results_df)} pathways analyzed successfully")
        if successful > 0:
            sig_raw = (results_df['p_value'] < 0.05).sum()
            print(f"- {sig_raw} pathways with p < 0.05")
            if 'fdr_p_value' in results_df.columns:
                sig_fdr = (results_df['fdr_p_value'] < 0.05).sum()
                print(f"- {sig_fdr} pathways with FDR < 0.05")

    return results_df


def plot_pathway_forest_plot(results_df, title="Pathway Cox Regression Forest Plot", 
                            p_threshold=0.05, top_n=None, figsize=(10, 12)):
    """
    Create a forest plot for pathway Cox regression results
    
    Parameters:
    - results_df: DataFrame with Cox regression results (must have 'pathway', 'hazard_ratio', 'hr_lower_ci', 'hr_upper_ci', 'p_value')
    - title: Plot title
    - p_threshold: P-value threshold for highlighting significant results
    - top_n: If specified, only plot top N pathways by p-value
    - figsize: Figure size tuple
    """
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Filter for successful results
    plot_df = results_df[results_df['status'] == 'success'].copy()
    
    if len(plot_df) == 0:
        print("No successful results to plot")
        return
    
    # Sort by p-value
    plot_df = plot_df.sort_values('p_value')
    
    # Limit to top N if specified
    if top_n is not None and len(plot_df) > top_n:
        plot_df = plot_df.head(top_n)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Assign y positions
    y_pos = np.arange(len(plot_df))
    
    # Color by significance
    colors = ['red' if p < p_threshold else 'gray' for p in plot_df['p_value']]
    
    # Plot horizontal bars for confidence intervals
    for i, (idx, row) in enumerate(plot_df.iterrows()):
        ax.plot([row['hr_lower_ci'], row['hr_upper_ci']], [i, i], 
               color=colors[i], linewidth=2, alpha=0.7)
    
    # Plot points for hazard ratios
    ax.scatter(plot_df['hazard_ratio'], y_pos, 
              c=colors, s=100, zorder=3, alpha=0.9, edgecolors='black', linewidth=1)
    
    # Add vertical line at HR = 1
    ax.axvline(x=1, color='black', linestyle='--', linewidth=1, alpha=0.5)
    
    # Customize axes
    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df['pathway'], fontsize=9)
    ax.set_xlabel('Hazard Ratio (95% CI)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax.grid(axis='x', alpha=0.3)
    
    # Add p-values as text
    for i, (idx, row) in enumerate(plot_df.iterrows()):
        p_text = f"p={row['p_value']:.3e}" if row['p_value'] < 0.001 else f"p={row['p_value']:.3f}"
        ax.text(ax.get_xlim()[1] * 0.95, i, p_text, 
               verticalalignment='center', fontsize=8, 
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='red', label=f'p < {p_threshold}'),
        Patch(facecolor='gray', label=f'p ≥ {p_threshold}')
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=10)
    
    plt.tight_layout()
    plt.show()
    
    print(f"\nForest Plot Summary:")
    print(f"- Total pathways plotted: {len(plot_df)}")
    print(f"- Significant pathways (p < {p_threshold}): {(plot_df['p_value'] < p_threshold).sum()}")


def plot_pathway_km_curves(adata, pathway_name, ulm_obsm_key='score_ulm', 
                          cutpoint_method='median', covariates=None,
                          figsize=(10, 6), title=None):
    """
    Plot Kaplan-Meier curves for a pathway with specified cutpoint method
    
    Parameters:
    - adata: AnnData object with survival data and ULM scores
    - pathway_name: Name of the pathway to plot
    - ulm_obsm_key: Key in adata.obsm containing ULM scores
    - cutpoint_method: 'median', 'quartiles', or 'optimal'
    - covariates: List of covariates for CoxPH statistics (optional)
    - figsize: Figure size tuple
    - title: Plot title (auto-generated if None)
    """
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test
    from scipy.stats import zscore
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    
    # Get ULM scores
    ulm_scores = pd.DataFrame(adata.obsm[ulm_obsm_key], 
                             index=adata.obs_names)
    
    if pathway_name not in ulm_scores.columns:
        raise ValueError(f"Pathway '{pathway_name}' not found in ULM scores")
    
    # Prepare survival dataframe
    survival_df = pd.DataFrame({
        'public_id': adata.obs['public_id'].values,
        'ttcpfs': pd.to_numeric(adata.obs['ttcpfs'], errors='coerce'),
        'censpfs': pd.to_numeric(adata.obs['censpfs'], errors='coerce'),
        'pathway_score': ulm_scores[pathway_name].values
    })
    
    # Add covariates if specified
    if covariates:
        for cov in covariates:
            if cov in adata.obs.columns:
                survival_df[cov] = adata.obs[cov].values
    
    # Remove missing data
    survival_df = survival_df.dropna(subset=['ttcpfs', 'censpfs', 'pathway_score'])
    survival_df = survival_df[survival_df['ttcpfs'] > 0]
    
    # Define groups based on cutpoint method
    if cutpoint_method == 'median':
        median_score = survival_df['pathway_score'].median()
        survival_df['group'] = (survival_df['pathway_score'] >= median_score).astype(int)
        survival_df['group_label'] = survival_df['group'].map({0: 'Low', 1: 'High'})
        groups = [0, 1]
        group_labels = ['Low', 'High']
        
    elif cutpoint_method == 'quartiles':
        quartiles = survival_df['pathway_score'].quantile([0.25, 0.5, 0.75])
        survival_df['group'] = pd.cut(survival_df['pathway_score'], 
                                      bins=[-np.inf, quartiles.iloc[0], quartiles.iloc[1], 
                                           quartiles.iloc[2], np.inf],
                                      labels=['Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)'])
        groups = survival_df['group'].unique()
        group_labels = groups.tolist()
        
    elif cutpoint_method == 'optimal':
        # Find optimal cutpoint between 25th and 75th percentile
        p25 = survival_df['pathway_score'].quantile(0.25)
        p75 = survival_df['pathway_score'].quantile(0.75)
        
        # Test cutpoints
        test_cutpoints = np.linspace(p25, p75, 50)
        best_p = np.inf
        best_cutpoint = survival_df['pathway_score'].median()
        
        for cutpoint in test_cutpoints:
            binary_group = (survival_df['pathway_score'] >= cutpoint).astype(int)
            if len(np.unique(binary_group)) == 2:
                try:
                    group0_mask = binary_group == 0
                    group1_mask = binary_group == 1
                    
                    if group0_mask.sum() > 5 and group1_mask.sum() > 5:
                        lr_result = logrank_test(
                            survival_df.loc[group0_mask, 'ttcpfs'],
                            survival_df.loc[group1_mask, 'ttcpfs'],
                            survival_df.loc[group0_mask, 'censpfs'],
                            survival_df.loc[group1_mask, 'censpfs']
                        )
                        if lr_result.p_value < best_p:
                            best_p = lr_result.p_value
                            best_cutpoint = cutpoint
                except:
                    continue
        
        survival_df['group'] = (survival_df['pathway_score'] >= best_cutpoint).astype(int)
        survival_df['group_label'] = survival_df['group'].map({0: 'Low', 1: 'High'})
        groups = [0, 1]
        group_labels = ['Low', 'High']
        print(f"Optimal cutpoint: {best_cutpoint:.3f} (p={best_p:.3e})")
    
    else:
        raise ValueError(f"Invalid cutpoint_method: {cutpoint_method}")
    
    # Plot KM curves
    fig, ax = plt.subplots(figsize=figsize)
    kmf = KaplanMeierFitter()
    
    colors = ['blue', 'red', 'green', 'purple']
    for i, group in enumerate(groups):
        if cutpoint_method == 'quartiles':
            mask = survival_df['group'] == group
            label = str(group)
        else:
            mask = survival_df['group'] == group
            label = group_labels[i]
        
        kmf.fit(survival_df.loc[mask, 'ttcpfs'],
               event_observed=survival_df.loc[mask, 'censpfs'],
               label=label)
        kmf.plot_survival_function(ax=ax, ci_show=True, color=colors[i % len(colors)])
    
    # Calculate statistics
    if len(groups) == 2:
        # Log-rank test
        group0_mask = survival_df['group'] == groups[0]
        group1_mask = survival_df['group'] == groups[1]
        
        lr_result = logrank_test(
            survival_df.loc[group0_mask, 'ttcpfs'],
            survival_df.loc[group1_mask, 'ttcpfs'],
            survival_df.loc[group0_mask, 'censpfs'],
            survival_df.loc[group1_mask, 'censpfs']
        )
        
        # CoxPH test if covariates specified
        if covariates:
            try:
                cox_df = survival_df.copy()
                cox_df['group_binary'] = cox_df['group']
                
                # Process covariates
                model_cols = ['ttcpfs', 'censpfs', 'group_binary']
                for cov in covariates:
                    if cov not in cox_df.columns:
                        continue
                    dtype = cox_df[cov].dtypes
                    if dtype in ('category', 'string', 'boolean', 'object'):
                        dummies = pd.get_dummies(cox_df[cov], prefix=cov, drop_first=True)
                        cox_df = pd.concat([cox_df, dummies], axis=1)
                        model_cols.extend([col for col in cox_df.columns if col.startswith(f'{cov}_')])
                    else:
                        scaled_col = f'scaled_{cov}'
                        cox_df[scaled_col] = zscore(cox_df[cov].values, nan_policy='omit')
                        model_cols.append(scaled_col)
                
                final_cox_df = cox_df[model_cols].dropna()
                cph = CoxPHFitter()
                cph.fit(final_cox_df, duration_col='ttcpfs', event_col='censpfs')
                
                hr = np.exp(cph.params_['group_binary'])
                hr_ci_lower = np.exp(cph.confidence_intervals_.loc['group_binary', '95% lower-bound'])
                hr_ci_upper = np.exp(cph.confidence_intervals_.loc['group_binary', '95% upper-bound'])
                cox_p = cph.summary.loc['group_binary', 'p']
                
                stats_text = f'Log-rank p = {lr_result.p_value:.3e}\n'
                stats_text += f'Cox HR = {hr:.2f} (95% CI: {hr_ci_lower:.2f}-{hr_ci_upper:.2f})\n'
                stats_text += f'Cox p = {cox_p:.3e}'
                if covariates:
                    stats_text += f'\nAdjusted for: {", ".join(covariates)}'
                
            except Exception as e:
                stats_text = f'Log-rank p = {lr_result.p_value:.3e}\n(Cox model failed: {str(e)[:50]})'
        else:
            stats_text = f'Log-rank p = {lr_result.p_value:.3e}'
    else:
        stats_text = f'{len(groups)} groups'
    
    # Customize plot
    ax.set_xlabel('Time (months)', fontsize=12)
    ax.set_ylabel('Progression-Free Survival', fontsize=12)
    
    if title is None:
        title = f'{pathway_name}: {cutpoint_method.capitalize()} Split'
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Add statistics text box
    ax.text(0.98, 0.05, stats_text, transform=ax.transAxes,
           verticalalignment='bottom', horizontalalignment='right',
           fontsize=10, bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.show()
    
    return survival_df


def plot_optimal_cutpoint_hr_curve(adata, pathway_name, ulm_obsm_key='score_ulm',
                                   covariates=None, figsize=(10, 5)):
    """
    Plot HR vs cutpoint for optimal cutpoint analysis
    
    Parameters:
    - adata: AnnData object with survival data and ULM scores
    - pathway_name: Name of the pathway to analyze
    - ulm_obsm_key: Key in adata.obsm containing ULM scores
    - covariates: List of covariates for adjustment
    - figsize: Figure size tuple
    
    Returns:
    - results_df: DataFrame with cutpoint, HR, CI, and p-value for each tested cutpoint
    - optimal_cutpoint: The optimal cutpoint value
    """
    from lifelines import CoxPHFitter
    from scipy.stats import zscore
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    
    # Get ULM scores
    ulm_scores = pd.DataFrame(adata.obsm[ulm_obsm_key], 
                             index=adata.obs_names)
    
    # Prepare survival dataframe
    survival_df = pd.DataFrame({
        'public_id': adata.obs['public_id'].values,
        'ttcpfs': pd.to_numeric(adata.obs['ttcpfs'], errors='coerce'),
        'censpfs': pd.to_numeric(adata.obs['censpfs'], errors='coerce'),
        'pathway_score': ulm_scores[pathway_name].values
    })
    
    # Add covariates if specified
    if covariates:
        for cov in covariates:
            if cov in adata.obs.columns:
                survival_df[cov] = adata.obs[cov].values
    
    # Remove missing data
    survival_df = survival_df.dropna(subset=['ttcpfs', 'censpfs', 'pathway_score'])
    survival_df = survival_df[survival_df['ttcpfs'] > 0]
    
    # Define cutpoint range (25th to 75th percentile)
    p25 = survival_df['pathway_score'].quantile(0.25)
    p75 = survival_df['pathway_score'].quantile(0.75)
    test_cutpoints = np.linspace(p25, p75, 50)
    
    # Store results
    cutpoint_results = []
    
    for cutpoint in test_cutpoints:
        try:
            cox_df = survival_df.copy()
            cox_df['group_binary'] = (cox_df['pathway_score'] >= cutpoint).astype(int)
            
            # Skip if groups are too small
            if cox_df['group_binary'].sum() < 10 or (1 - cox_df['group_binary']).sum() < 10:
                continue
            
            # Build Cox model
            model_cols = ['ttcpfs', 'censpfs', 'group_binary']
            
            if covariates:
                for cov in covariates:
                    if cov not in cox_df.columns:
                        continue
                    dtype = cox_df[cov].dtypes
                    if dtype in ('category', 'string', 'boolean', 'object'):
                        dummies = pd.get_dummies(cox_df[cov], prefix=cov, drop_first=True)
                        cox_df = pd.concat([cox_df, dummies], axis=1)
                        model_cols.extend([col for col in cox_df.columns if col.startswith(f'{cov}_')])
                    else:
                        scaled_col = f'scaled_{cov}'
                        cox_df[scaled_col] = zscore(cox_df[cov].values, nan_policy='omit')
                        model_cols.append(scaled_col)
            
            final_cox_df = cox_df[model_cols].dropna()
            
            # Fit Cox model
            cph = CoxPHFitter()
            cph.fit(final_cox_df, duration_col='ttcpfs', event_col='censpfs')
            
            hr = np.exp(cph.params_['group_binary'])
            hr_ci_lower = np.exp(cph.confidence_intervals_.loc['group_binary', '95% lower-bound'])
            hr_ci_upper = np.exp(cph.confidence_intervals_.loc['group_binary', '95% upper-bound'])
            cox_p = cph.summary.loc['group_binary', 'p']
            
            cutpoint_results.append({
                'cutpoint': cutpoint,
                'hr': hr,
                'hr_ci_lower': hr_ci_lower,
                'hr_ci_upper': hr_ci_upper,
                'p_value': cox_p
            })
        except:
            continue
    
    results_df = pd.DataFrame(cutpoint_results)
    
    if len(results_df) == 0:
        print("Failed to find valid cutpoints")
        return None, None
    
    # Find optimal cutpoint (smallest p-value)
    optimal_idx = results_df['p_value'].idxmin()
    optimal_cutpoint = results_df.loc[optimal_idx, 'cutpoint']
    optimal_hr = results_df.loc[optimal_idx, 'hr']
    optimal_p = results_df.loc[optimal_idx, 'p_value']
    
    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot HR with confidence intervals
    ax.plot(results_df['cutpoint'], results_df['hr'], color='blue', linewidth=2, label='Hazard Ratio')
    ax.fill_between(results_df['cutpoint'], 
                    results_df['hr_ci_lower'], 
                    results_df['hr_ci_upper'],
                    alpha=0.3, color='blue', label='95% CI')
    
    # Add horizontal line at HR = 1
    ax.axhline(y=1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    
    # Add vertical line at optimal cutpoint
    ax.axvline(x=optimal_cutpoint, color='red', linestyle='--', linewidth=2, 
              label=f'Optimal cutpoint: {optimal_cutpoint:.3f}')
    
    # Customize
    ax.set_xlabel('Cutpoint Value', fontsize=12)
    ax.set_ylabel('Hazard Ratio (95% CI)', fontsize=12)
    ax.set_title(f'{pathway_name}: HR vs Cutpoint', fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.grid(alpha=0.3)
    
    # Add text box with optimal cutpoint info
    stats_text = f'Optimal Cutpoint: {optimal_cutpoint:.3f}\n'
    stats_text += f'HR = {optimal_hr:.2f}\n'
    stats_text += f'p = {optimal_p:.3e}'
    if covariates:
        stats_text += f'\nAdjusted for: {", ".join(covariates)}'
    
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
           verticalalignment='top', fontsize=10,
           bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.show()
    
    print(f"\nOptimal cutpoint analysis:")
    print(f"- Tested {len(results_df)} cutpoints between {p25:.3f} and {p75:.3f}")
    print(f"- Optimal cutpoint: {optimal_cutpoint:.3f}")
    print(f"- HR at optimal cutpoint: {optimal_hr:.2f} (p={optimal_p:.3e})")
    
    return results_df, optimal_cutpoint
