# Basic imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import argparse

# Import saga modules
import saga.plot as sagap

parser = argparse.ArgumentParser(description='Script to plot 2D kinematic distributions from `getBinKinematicsTH2Ds` jobs on RGA data for the `Lambda -> proton pion` channel')
parser.add_argument('--watermark', action="store_true", help='Plot watermark on plots')
parser.add_argument('--rgs', default=["dt_rga"], help='Run group', nargs="+", choices=['dt_rga','mc_rga','mc_rga_sss'])
parser.add_argument('--binvars', default=["x"], help='Bin variables', nargs="+", choices=['Q2','W','y','x','xF_ppim','z_ppim'])
parser.add_argument('--hist_density', action="store_true", help='Plot normalized histograms')
parser.add_argument('--grid_shape', default=(5,1), help='Grid shape', nargs=2, type=int)
args = parser.parse_args()

# Set base directory from environment
RGA_LAMBDA_ANALYSIS_HOME = os.environ['RGA_LAMBDA_ANALYSIS_HOME']

# Set channels and run groups to loop
chs = ['ppim']
ch_labels = {'ppim':'p\\pi^{-}'}
rgs = args.rgs #['dt_rga','mc_rga','mc_rga_sss']
rg_labels = {'dt_rga':'Data RGA','mc_rga':'MC RGA CD','mc_rga_sss':'MC RGA SS'}
binvars = args.binvars # binvars = ['Q2','W','y','x','xF_ppim','z_ppim']
grid_shape = args.grid_shape #(5,1)

# Loop run groups, channels
for rg in rgs:
    for ch in chs:

        # Set binvar titles
        binvar_titles = {
            'x'            : 'x',
            'Q2'           : 'Q^{2}',
            f'W'           : 'W',
            f'y'           : 'y',
            f'xF_{ch}'     : 'x_{F '+ch_labels[ch]+'}',
            f'z_{ch}'      : 'z_{'+ch_labels[ch]+'}',
        }
        binvar_unit_titles = {
            'x'            : '',
            'Q2'           : ' (GeV$^{2}$)',
            'W'            : ' (GeV)',
            'y'            : '',
            f'xF_{ch}'     : '',
            f'z_{ch}'      : '',
        }

        # Loop 1D bin variables
        for binvar in binvars:

            # Setup, modify these as needed for your specific binning scheme
            csv_path = os.path.abspath(
                os.path.join(
                    RGA_LAMBDA_ANALYSIS_HOME,
                    f'jobs/saga/test_getBinKinematics__{ch}/out_{rg}_1d_bins_{binvar}_kinematics.csv'
                )
            )
            hist_path = os.path.abspath(
                os.path.join(
                    RGA_LAMBDA_ANALYSIS_HOME,
                    f'jobs/saga/test_getBinKinematicsTH2Ds__{ch}/out_{rg}_1d_bins_{binvar}_kinematics.root'
                )
            )
            print(f"INFO: csv_path: {csv_path}")
            print(f"INFO: hist_path: {csv_path}")

            # Set kinematic variable pairs to plot
            kinvars = binvars.copy()
            kinvars.remove(binvar)

            kinvar_pairs = [['x','Q2'],['x','W']]

            # Loop kinematic variable pairs
            for kinvars in kinvar_pairs:

                kinvar_x, kinvar_y = kinvars

                # Grab kinematic variables and set related info
                xlabels = {
                    'x'             : '$x$',
                    'Q2'            : '$Q^{2}$ (GeV$^{2}$)',
                    'W'             : '$W$ (GeV)',
                    'y'             : '$y$',
                    f'xF_{ch}'    : '$x_{F '+ch_labels[ch]+'}$',
                    f'z_{ch}'       : '$z_{'+ch_labels[ch]+'}$',
                }                    
                xlims = {
                    'x'             : [0.0,1.0],
                    'Q2'            : [1.0,10.0],
                    'W'             : [2.0,4.0],
                    'y'             : [0.0,1.0],
                    f'xF_{ch}'      : [0.0,1.0],
                    f'z_{ch}'       : [0.0,1.0],
                }
                hist_colors = {
                    'x'             : ['tab:purple'],
                    'Q2'            : ['tab:cyan'],
                    'W'             : ['tab:olive'],
                    'y'             : ['tab:blue'],
                    f'xF_{ch}'      : ['tab:red'],
                    f'z_{ch}'       : ['tab:orange'],
                }

                # Only add bin variable values above plots for now
                cols = [
                    binvar,
                ]
                col_titles = {
                    binvar:binvar_titles[binvar],
                }
                col_unit_titles = {
                    binvar:binvar_unit_titles[binvar],
                }

                # Load kinematics CSV
                df = pd.read_csv(csv_path)
                bin_ids = df['bin'].unique().tolist()

                # Set graph and plot_results arrays
                graph_array = [{} for j in range(len(bin_ids))]
                plot_results_kwargs_array = [
                        {
                            'hist_keys':[f'h2_bin{bin_id}_'+kinvar_x+'_'+kinvar_y],
                            'title':sagap.get_bin_kinematics_title(
                                    bin_id,df,
                                    cols=cols,
                                    col_titles=col_titles,
                                    col_unit_titles=col_unit_titles
                                ),
                            'xlims':xlims[kinvar_x],
                            'xlabel':xlabels[kinvar_x],
                            'ylims':ylims[kinvar_y],
                            'ylabel':ylabels[kinvar_y],
                        }
                        for bin_id in bin_ids
                ]

                def reshape_grid(grid_array, grid_shape):
                    # Reshape the grid array to match the specified shape
                    reshaped_grid = []
                    for i in range(grid_shape[0]):
                        el = grid_array[i * grid_shape[1]:min(len(bin_ids),(i + 1) * grid_shape[1])]
                        if len(el) < grid_shape[1]:
                            el.extend([None for _ in range(grid_shape[1] - len(el))])
                        reshaped_grid.append(el)
                    return reshaped_grid

                # grid_shape = (5,1)
                graph_array = reshape_grid(graph_array, grid_shape)
                plot_results_kwargs_array = reshape_grid(plot_results_kwargs_array, grid_shape)
                print("INFO: graph_array = ",graph_array)
                print("INFO: plot_results_kwargs_array = ", plot_results_kwargs_array)

                # Set base kwargs
                plot_results_kwargs_base = {
                    'show_injected_asymmetries':False,
                    'hist_clone_axis':False,
                    'hist_paths':[hist_path],
                    'hist_labels':[f'{rg_labels[rg]}'],
                    'watermark':'CLAS12 Preliminary' if args.watermark else '',
                    'hist_density':hist_density,
                    'axlinewidth':0,
                    'hist_dim':2,
                    'legend_loc':None #NOTE: Do not plot a legend if you are using 2d hists.
                }

                # Set additional kwargs
                figsize = (16*grid_shape[1],10*grid_shape[0])
                outpath = f'plot_TH2Ds__{rg}__{ch}__1d_{kinvar_x}_{kinvar_y}__1d_{binvar}.pdf'
                use_default_plt_settings = True
                use_grid_titles = False
                use_grid_xlabels = True

                # Plot an array of graphs
                sagap.plot_results_array(
                    graph_array,
                    plot_results_kwargs_array,
                    plot_results_kwargs_base = plot_results_kwargs_base,
                    figsize = figsize,
                    outpath = outpath,
                    use_default_plt_settings = use_default_plt_settings,
                    use_grid_titles = use_grid_titles,
                    use_grid_xlabels = use_grid_xlabels, #NOTE: Since you plot different x-axis variables in each row make sure that the labels are not dropped for rows below the top row.
                )

                # Close the figures to conserve memory
                plt.close('all')
