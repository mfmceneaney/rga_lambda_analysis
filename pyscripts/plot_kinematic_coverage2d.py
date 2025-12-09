# Basic imports
import numpy as np
import os
import sys
import uproot as ur
import matplotlib.pyplot as plt
from matplotlib import colors
import argparse

# Import saga modules
from saga.aggregate import get_binscheme_cuts_and_ids
from saga.data import load_yaml, load_th1, load_csv
from saga.plot import (
    set_default_plt_settings,
    plot_th2,
    get_lims_coords,
    plot_lines,
    get_bin_centers,
    plot_bin_ids,
    plot_watermark,
)

parser = argparse.ArgumentParser(description='Script to plot 2d kinematic coverage in x vs. Q2 and x vs. W from `getBinKinematicsTH2Ds` jobs on RGA data for the `Lambda -> proton pion` channel')
parser.add_argument('--watermark', action="store_true", help='Plot watermark on plots')
parser.add_argument('--rgs', default=["dt_rga"], help='Run group', nargs="+", choices=['dt_rga','mc_rga','mc_rga_sss'])
args = parser.parse_args()

# Set base directory from environment
RGA_LAMBDA_ANALYSIS_HOME = os.environ['RGA_LAMBDA_ANALYSIS_HOME']

# Set channels and beam suffixes to loop
chs = ['ppim']
ch_labels = {'ppim':'p\\pi^{-}'}
rgs = args.rgs
rg_labels = {'dt_rga':'Data RGA','mc_rga':'MC RGA','mc_rga_SSS':'MC RGA SS'}

# Loop run groups, channels, and beam suffixes
for rg in rgs:
    for ch in chs:

            # Set bin variable labels and limits inside loop since they depend on channel and energy
            binvar_labels = {
                'x':'$x$',
                'Q2':'$Q^{2}$ (GeV)$^{2}$',
                'W':'$W$ (GeV)',
            }
            binvar_lims = {
                'x':[0.0,1.0],
                'Q2':[1.0,11.0],
                'W':[2.0,4.0],
            }

            # Set bin scheme variable pairs to loop
            binvars_pairs = [['x','Q2'],['x','W']]

            for binvars in binvars_pairs:

                # Setup, modify these as needed for your specific binning scheme
                sep = '_'
                hist_path = os.path.abspath(
                    os.path.join(
                        RGA_LAMBDA_ANALYSIS_HOME,
                        f'jobs/saga/test_getBinKinematicsTH2Ds__{ch}/',
                        f'out_{rg}_fullbin_binscheme_kinematics.root'
                    )
                )
                hist_name       = 'h2_bin0_'+sep.join(binvars)
                outpath         = f'kinematic_coverage2d_{binvars[0]}_{binvars[1]}_{rg}_{ch}.pdf'
                var_keys        = binvars #NOTE: This should only be set in the case of a 2D grid scheme.
                start_idx       = 0
                id_key          = 'bin_id'

                # Load TH2 histogram with uproot
                h2 = load_th1(hist_path,name=hist_name)

                # Set plt settings
                set_default_plt_settings()

                # Open the figure
                f, ax = plt.subplots(figsize=(16,10))

                # Plot the 2D distribution
                plot_th2(h2, ax, norm=colors.LogNorm())
                ax.set_xlabel(binvar_labels[binvars[0]],usetex=True)
                ax.set_ylabel(binvar_labels[binvars[1]],usetex=True)
                ax.set_title(f'{rg_labels[rg]} ${ch_labels[ch]}$',usetex=True)

                # Plot watermark
                if args.watermark:
                    watermark = "CLAS12 Preliminary"
                    plot_watermark(
                        watermark,
                        size=75,
                        rotation=25.0,
                        color="gray",
                        alpha=0.5,
                    )

                # Save the figure
                f.savefig(outpath)

                # Close figures
                plt.close()
