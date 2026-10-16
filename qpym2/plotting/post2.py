
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import arviz as az
import xarray as xr
from scipy import stats
from typing import List, Optional, Union, Dict, Tuple

from qpym2.io import fit as fit_io
from qpym2.fit.params import model_cfg

# =============================================================================
# Helper: Draw Prior
# =============================================================================

def _draw_prior(ax: plt.Axes, prior_tuple: Tuple, data_range: Tuple[float, float]):
    """
    Draws the analytic prior PDF on the given axes.
    prior_tuple format: (type_str, params_dict) e.g. ('normal', {'mu':0, 'sigma':1})
    """
    if not isinstance(prior_tuple, (list, tuple)) or len(prior_tuple) != 2:
        return

    ptype, params = prior_tuple
    x = np.linspace(data_range[0], data_range[1], 200)
    y = np.zeros_like(x)
    label = f'Prior ({ptype})'

    try:
        if ptype == 'normal':
            mu, sigma = params['mu'], params['sigma']
            y = stats.norm.pdf(x, loc=mu, scale=sigma)
        elif ptype == 'halfnormal':
            # PyMC HalfNormal is defined for x > 0.
            # Scipy halfnorm is folded at loc.
            sigma = params['sigma']
            y = stats.halfnorm.pdf(x, scale=sigma)
        elif ptype == 'flat':
            lower, upper = params['lower'], params['upper']
            mask = (x >= lower) & (x <= upper)
            y[mask] = 1.0 / (upper - lower)
            label += f'\n[{lower}, {upper}]'
        elif ptype == 'truncated_normal':
             # specialized handling if needed
             mu, sigma= params['mu'], params['sigma']
             lower, upper = params.get('lower', -np.inf), params.get('upper', np.inf)
             a, b = (lower - mu) / sigma, (upper - mu) / sigma
             y = stats.truncnorm.pdf(x, a, b, loc=mu, scale=sigma)
        else:
            return # unsupported type

        # Scale if necessary? No, assumes density plot for posterior (stat='density')
        ax.plot(x, y, 'r--', lw=2, label=label)

    except Exception as e:
        print(f"Failed to draw prior {ptype}: {e}")

# =============================================================================
# 1. Log Likelihood Validation
# =============================================================================

def plot_loglike_validation(
    df_fit: pd.DataFrame, 
    df_toys: pd.DataFrame, 
    model_name: Optional[str] = None,
    titles: bool = True,
    ax: Optional[plt.Axes] = None
) -> plt.Axes:
    """
    Compares the data fit log-likelihood against the distribution of toy log-likelihoods.
    
    Args:
        df_fit: DataFrame containing the data fit result (single row or filtered). 
                Must have 'logp' or 'log_likelihood' column.
        df_toys: DataFrame containing toy fit results.
        model_name: Label for the current model (optional).
        ax: Matplotlib axes.
        
    Returns:
        plt.Axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
        
    # extract logp
    logp_col = 'logp' if 'logp' in df_toys.columns else 'log_likelihood'
    if logp_col not in df_toys.columns:
        # fallback for 'logp_mode' etc
        cols = [c for c in df_toys.columns if 'logp' in c]
        if cols: logp_col = cols[0]
        else: raise ValueError(f"Log likelihood column not found in toys df. Columns: {df_toys.columns}")

    # Plot toy distribution
    sns.histplot(df_toys[logp_col], ax=ax, stat='density', kde=True, label='Toys', color='C0', alpha=0.6)
    
    # Plot data value
    # Assuming df_fit has one row or we take the mean/mode
    val_col = logp_col if logp_col in df_fit.columns else [c for c in df_fit.columns if 'logp' in c][0]
    data_logp = df_fit[val_col].iloc[0]

    # Calculate statistics
    toys_vals = df_toys[logp_col].dropna()
    p16, p84 = np.percentile(toys_vals, [16, 84])
    
    # Determine color based on interval (68% inner interval)
    if (data_logp < p16) or (data_logp > p84):
        line_color = 'red'
    else:
        line_color = 'green'
    
    ax.axvline(data_logp, color=line_color, linestyle='--', linewidth=2, label=f'Data: {data_logp:.1f}')
    
    # Calculate p-value or z-score
    z_score = (data_logp - toys_vals.mean()) / toys_vals.std()
    
    if titles:
        title = "Log-Likelihood Validation"
        if model_name:
            title += f"\n{model_name}"
        title += f" (Z-score: {z_score:.2f})"
        
        ax.set_title(title)
        ax.set_xlabel("Log Likelihood")
    ax.legend()
    
    return ax


# =============================================================================
# 2. Sensitivity Comparison
# =============================================================================

def plot_sensitivity_comparison(
    df_sens: pd.DataFrame,
    groupby: str = 'co_cut',
    metric: str = 'ndbd_p90',
    ax: Optional[plt.Axes] = None
) -> plt.Axes:
    """
    Compares sensitivity stats (e.g. 90% CL limit) across different model configurations.
    
    Args:
        df_sens: DataFrame with sensitivity results. Expected to have 'model_name' 
                 or config columns (co_cut, ndbd_cutoff) and the metric column.
        groupby: Column to group/hue by (e.g., 'co_cut').
        metric: Column name for the sensitivity limit (default 'ndbd_p90' or similar).
        ax: Matplotlib axes.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
        
    # Check if metric exists
    if metric not in df_sens.columns:
        # fuzzy search
        cands = [c for c in df_sens.columns if 'p90' in c]
        if cands: metric = cands[0]
        else: raise ValueError(f"Metric {metric} not found in dataframe.")

    # Determine x-axis variables. Usually usually we compare across models.
    # If index is meaningful, use it. Otherwise create a label.
    plot_df = df_sens.copy()
    
    # Create a cleaner label if needed
    if 'label' not in plot_df.columns:
        # Try to construct from config cols
        cfg_cols = [c for c in ['ndbd_cutoff', 'co_cut', 'co_shift'] if c in plot_df.columns]
        if cfg_cols:
            plot_df['label'] = plot_df[cfg_cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
        else:
            plot_df['label'] = plot_df.index
            
    # Boxplot or Violinplot if multiple entries per model (e.g. toys), 
    # or Pointplot/Barplot if summary dataframe.
    # Assuming df_sens summary here (one row per model).
    
    # Check if we have multiple rows per label
    is_summary = None #plot_df.groupby('label')[metric].count().max() == 1
    
    if is_summary:
        sns.scatterplot(data=plot_df, x='label', y=metric, hue=groupby if groupby in plot_df else None, s=100, ax=ax)
        # Optional: Add error bars if available
    else:
        sns.boxplot(data=plot_df, x='label', y=metric, hue=groupby if groupby in plot_df else None, ax=ax)
        
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_ylabel("Sensitivity (90% C.L.)")
    ax.set_title("Sensitivity Comparison")
    
    return ax


# =============================================================================
# 3. Co-Shift Validity (Fit Validity vs Co-Shift)
# =============================================================================

def plot_coshift_validity(
    df_fits: pd.DataFrame,
    df_toys: Optional[pd.DataFrame] = None,
    metric: str = 'logp',
    ax: Optional[plt.Axes] = None
) -> plt.Axes:
    """
    Plots the fit metric (e.g. LogLikelihood) vs Cobalt Shift.
    Optionally overlays toy distribution ranges (validity check).
    
    Args:
        df_fits: DataFrame with data fit results (must have 'co_shift' and metric).
        df_toys: (Optional) DataFrame with toy results to show expected bands.
        metric: 'logp' or other fit statistic.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    if 'co_shift' not in df_fits.columns:
        raise ValueError("DataFrame must have 'co_shift' column")
        
    # Clean metric name
    if metric not in df_fits.columns:
         cands = [c for c in df_fits.columns if metric in c]
         if cands: metric = cands[0]

    # Plot Data Line
    sns.lineplot(data=df_fits, x='co_shift', y=metric, marker='o', label='Data Fit', ax=ax, color='black')
    
    # Overlay Toys if provided (as a band)
    if df_toys is not None:
        if 'co_shift' not in df_toys.columns:
             # Try to merge/map if missing
             pass
        else:
            # Calculate aggregations per shift
            toys_agg = df_toys.groupby('co_shift')[metric].agg(['mean', 'std']).reset_index()
            # Plot bands
            ax.fill_between(
                toys_agg['co_shift'], 
                toys_agg['mean'] - toys_agg['std'], 
                toys_agg['mean'] + toys_agg['std'], 
                color='C0', alpha=0.3, label='Toy Expected (±1σ)'
            )
            ax.plot(toys_agg['co_shift'], toys_agg['mean'], '--', color='C0')

    ax.set_xlabel('Cobalt Shift (keV)')
    ax.set_ylabel(metric)
    ax.set_title('Fit Validity vs. Cobalt Shift')
    ax.legend()
    
    return ax


# =============================================================================
# 4. Prior vs Posterior
# =============================================================================

def plot_prior_posterior_file(
    fit_path: str,
    params: List[str],
    ax: Optional[plt.Axes] = None
) -> plt.Axes:
    """
    Plots Prior vs Posterior from a fit result file.
    
    Args:
        fit_path: Path to the fit directory (containing fit.nc / fit.h5).
        params: List of parameter names to plot.
    """
    # 1. Load Fit
    # Assuming qpym2.io.fit.read_fit returns (model_config, trace, ...) or similar
    # Adjust based on actua library signature
    try:
        # Standard QPyM2 read
        # returns: (model, trace, data?) - check signature in your env
        comps, trace, _ = fit_io.read_fit(fit_path)
    except Exception as e:
        print(f"Failed to load fit from {fit_path}: {e}")
        return ax

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    # 2. Iterate and Plot
    # Only supporting 1st param for 1D, or first 2 for 2D if requested separately
    # This function creates 1D overlay
    
    param = params[0] # Focus on primary
    
    # Get Posterior Samples
    post_samples = trace.posterior[param].values.flatten()
    
    # Plot Posterior
    sns.histplot(post_samples, kde=True, stat='density', label='Posterior', ax=ax, color='C1', alpha=0.5)
    
    # Plot Prior 
    if param in comps.index:
        prior_info = comps.loc[param].get('prior')
        if prior_info is not None:
             # Expand x range slightly for plotting context
             x_min, x_max = post_samples.min(), post_samples.max()
             padding = (x_max - x_min) * 0.2
             data_range = (x_min - padding, x_max + padding)
             
             _draw_prior(ax, prior_info, data_range)
        
    ax.set_title(f"Prior vs Posterior: {param}")
    ax.legend()
    
    return ax


def plot_correlation_2d(
    fit_path: str,
    x_param: str,
    y_param: str,
    ax: Optional[plt.Axes] = None
) -> plt.Axes:
    """
    Plots 2D correlation/corner plot for two parameters.
    """
    _, trace, _ = fit_io.read_fit(fit_path)
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 7))
        
    az.plot_pair(
        trace,
        var_names=[x_param, y_param],
        kind='kde',
        marginals=True, # Shows 1D dists on sides
        ax=ax
    )
    
    return ax
