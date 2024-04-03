##################################################
# PLOT UTILITIES
# Author: Suk Yee Yong
##################################################

from collections import namedtuple
from mapie.metrics import regression_coverage_score, regression_mean_width_score
import matplotlib.pyplot as plt
import numpy as np


def cb_qualitative_ptol():
    """Colorblind friendly muted qualitative color scheme from: https://personal.sron.nl/~pault/"""
    cset = namedtuple('qualitative_ptol', 'rose indigo sand green cyan wine teal olive purple pale_grey black')
    return cset('#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE',
                '#882255', '#44AA99', '#999933', '#AA4499', '#DDDDDD',
                '#000000')


# Set matplotlib default color
from cycler import cycler
import matplotlib as mpl
mpl.rcParams['axes.prop_cycle'] = cycler(color=cb_qualitative_ptol())


def histridgeline_lifespanpi(df_data, y_labelname, df_datareport=None, use_baggedpred=True, strategy=None, xlabel='', ylabel='', save_plotname=False):
    """Histogram ridgeline plot of lifespan predictions and prediction intervals for each species"""
    
    def add_plusminuspred(df, ax):
        """Add +/- prediction interval"""
        pi_upper = df['pred_upper']-df['pred']
        pi_lower = df['pred']-df['pred_lower']
        ax.text(0.98, 0.9, f"+{pi_upper:.2f}\n-{pi_lower:.2f}", ha='right', va='top', transform=ax.transAxes, fontsize=6)
    
    ncol = 2
    nrows = len(df_data)
    nrow = int(np.ceil(nrows/ncol))
    fracheight_subplot = 0.4
    # Plot results sorted by _y_labelname
    fig, axs = plt.subplots(nrow, ncol, figsize=(8, nrow*fracheight_subplot+2), sharex='col', gridspec_kw={'hspace': 0, 'wspace': 0.5})
    axs = axs.T.ravel()
    # Sort by known values else use predicted
    for i, (_, od) in enumerate(df_data.sort_values('mean_'+y_labelname if 'mean_'+y_labelname in df_data.columns else 'bagged_detrans_predicted_'+y_labelname, ascending=False).iterrows()):
        if df_datareport is not None:
            datareport = df_datareport.loc[od['organism_name']]
            datareport_lifespan = datareport[y_labelname].values
            lifespan_mean = datareport['mean_'+y_labelname].values[0]
            lifespan_min, lifespan_max = min(datareport_lifespan), max(datareport_lifespan)
            # If covers all reported data
            # if (od['pred_lower'] <= lifespan_min) and (lifespan_max <= od['pred_upper']):
            # If covers reported mean
            if (od['pred_lower'] <= lifespan_mean) and (lifespan_mean <= od['pred_upper']):
                add_plusminuspred(od, axs[i])
            # Reported data
            axs[i].hist(datareport_lifespan, bins=np.arange(np.floor(lifespan_min)-0.5, np.ceil(lifespan_max)+1, 1), facecolor='grey', edgecolor='none', alpha=0.6, label='reported')
            axs[i].axvline(x=lifespan_mean, ls='--', lw=2.2, color='grey', label='mean')
        else:
            add_plusminuspred(od, axs[i])
        # Prediction from elastic net
        axs[i].axvline(x=od[('bagged_' if use_baggedpred else '') + 'detrans_predicted_'+y_labelname], ls='--', lw=2, label='glmnet enet prediction')
        # Prediction and prediction interval
        axs[i].axvline(x=od['pred'], lw=1.5, color=cb_qualitative_ptol().indigo, label=(f"{strategy} " if strategy is not None else '')+'prediction')
        axs[i].axvspan(od['pred_lower'], od['pred_upper'], edgecolor='none', facecolor=cb_qualitative_ptol().indigo, alpha=0.2, label=(f"{strategy} " if strategy is not None else '')+'prediction interval')
        # Labels
        axs[i].set_ylabel(od['organism_name'].replace(' ', '\n'), ha='right', va='center', rotation=0, fontsize=8)
        label_bottom = True if (i==nrows-1) or (i==nrow-1) else False
        axs[i].tick_params(left=False, labelleft=False, bottom=label_bottom, labelbottom=label_bottom, length=2)
        axs[i].set_ylim(None, axs[i].get_ylim()[-1]+0.5) # Increase height of y-axis
        axs[i].grid(axis='x', color='grey', ls=':', lw=0.5, alpha=0.5)
    if nrows % 2 != 0: axs[-1].set_axis_off() # Remove empty cell
    fig.text(-0.05, 0.5, 'Species' if not xlabel else xlabel, ha='center', va='bottom', transform=fig.transFigure, rotation='vertical', fontsize=10)
    fig.text(-0.6, -0.8, 'Lifespan (years)' if not ylabel else ylabel, ha='center', va='top', transform=axs[-1].transAxes, fontsize=10)
    legend_loc = (1,2.8) # Center legend
    # Add metrics
    if df_data.get('mean_'+y_labelname) is not None:
        axs[0].annotate(f"PICP={regression_coverage_score(df_data['mean_'+y_labelname], df_data['pred_lower'], df_data['pred_upper']):.3f} \n MPIW={regression_mean_width_score(df_data['pred_lower'], df_data['pred_upper']):.3f}", xy=(2.45, 2), xycoords='axes fraction', va='center', ha='right', bbox=dict(boxstyle='square', fc='none', edgecolor=(0,0,0,0.2)))
        legend_loc = (0.6,2.8) # Move legend to left
    axs[0].legend(loc='upper center', ncol=3, bbox_to_anchor=legend_loc, fontsize=10)
    fig.tight_layout()
    if isinstance(save_plotname, str): fig.savefig(save_plotname, bbox_inches='tight')


