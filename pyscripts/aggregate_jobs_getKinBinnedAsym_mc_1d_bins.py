# Basic imports
import numpy as np
import matplotlib.pyplot as plt
import os

# Import saga modules
import saga.aggregate as sagas
from saga.data import load_yaml, load_csv, save_bin_mig_mat_to_csv
from saga.plot import set_default_plt_settings, plot_results

# Set base directory from environment
RGA_LAMBDA_ANALYSIS = os.environ['RGA_LAMBDA_ANALYSIS']

# Set up chaining for batched data (specifically `old_dat_path`)
nbatch = 1
nbatches = {"nbatches":[nbatch]}
ibatches = {"ibatch":[i for i in range(nbatch)]}
chain_keys = ["nbatches", "ibatch"]
chain_configs = dict(
    nbatches,
    **ibatches,
) if nbatch > 1 else {}

# Set base directories to aggregate
run_groups = ['dt_rga']
channels   = ['ppim']
methods    = ["Asym", "HB"]
base_dirs  = [
    os.path.abspath(
        os.path.join(
                RGA_LAMBDA_ANALYSIS,
                f'jobs/saga/test_getKinBinned{method}__{rg}__{ch}__1D/'
            )
        ) for rg in run_groups for ch in channels for method in methods
]

# Set list of channels for each base directory
chs = [ch for rg in run_groups for ch in channels for method in methods]

# Set channel label for each base directory
ch_sgasym_labels = {
    'ppim':'$D^{\Lambda}_{LL\'}$',
}
ch_sgasym_labels = [ch_sgasym_labels[ch] for rg in run_groups for ch in channels]

# Set x-axis labels for kinematic variables in all channels
xlabel_map = {
    'Q2':'$Q^{2}$ (GeV$^{2}$)', 'W':'$W$ (GeV)', 'y':'$y$', 'x':'$x$', 
    'z_ppim':'$z_{p\pi^{-}}$', 'xF_ppim':'$x_{F p\pi^{-}}$',
    'mass_ppim':'$M_{p\pi^{-}}$ (GeV)',
}

# Set up list of run groups
rgs = [ rg for rg in run_groups for ch in channels]

# Loop base directories
for rg, base_dir, ch_sgasym_label in zip(rgs,base_dirs,ch_sgasym_labels):

    # Now loop signal asymmetries
    for ch_sgasym_label_idx, result_name in enumerate(ch_sgasym_label):

        # Setup input paths
        submit_path  = os.path.join(base_dir,"submit.sh")
        yaml_path    = os.path.join(base_dir,"args.yaml")
        out_path     = os.path.join(base_dir,"jobs.txt")
        #NOTE: Set the bin migration path below since this is binscheme dependent

        # Set aggregate keys
        aggregate_keys = ["inject_seed"]

        # Load the binschemes from the path specified in the job yaml assuming there is only one given path and it is an absolute path
        binschemes_paths_name = "binschemes_paths"
        yaml_path = load_yaml(yaml_path)[binschemes_paths_name][0]
        binschemes = load_yaml(yaml_path)

        # # Arguments for sagas.get_config_list() #NOTE: Set this above
        # result_name = "a0" #NOTE: This also gets recycled as the asymmetry name

        # Arguments for sagas.get_out_dirs_list()
        sep='_'
        ext='.pdf'

        # Arguments for sagas.get_out_file_name()
        out_file_name_ext = '.csv'
        bin_mig_base_name="bin_mig_mat_"

        # Arguments for sagas.apply_bin_mig()
        use_bin_mig = False
        id_gen_key='binid_gen'
        id_rec_key='binid_rec'
        mig_key='mig'
        results_keys = [result_name] #NOTE: You can apply bin migration to multiple dataframe entries in one go.

        # Arguments for sagas.get_graph_data()
        id_key = 'bin_id'

        # Arguments for sagas.get_graph_data()
        count_key  = 'count'
        asym_key   = result_name #NOTE: This is set from above
        err_ext    = '_err'

        # Arguments for sagas.plot_results()
        plot_results_kwargs_base = {
            'ylims':[-0.2,0.2],
            'sgasyms':[0.0], #NOTE: This will be set below for each configuration
            'sgasym_idx':ch_sgasym_label_idx,
            'sgasym_labels':[ch_sgasym_label[el] for el in ch_sgasym_label],
            'sg_colors':['blue','red','green','tab:pink', 'tab:purple', 'tab:gray', 'tab:orange', 'tab:cyan'],
            'bgasyms':[],
            'bgasym_labels':[],
            'bg_colors':[],
            'show_injected_asymmetries':True,
            'hist_paths':[],
            'hist_colors':[],
            'hist_keys':[],
            'hist_labels':[],
            'watermark':'CLAS12 Preliminary',
            'hist_clone_axis':False,
            'old_dat_path':None
        }

        # Additional useful parameters for plotting
        figsize = (16,10)
        #NOTE: Set outpath within the loop for unique naming
        use_default_plt_settings = True

        # If you want to rescale your results using results from other base directories set the following arguments
        rescale = False
        if rescale:
            plot_results_kwargs_base = dict(
                plot_results_kwargs_base,
                **{
                    'old_dat_path':'',
                    'new_sim_path':'',
                    'old_sim_path':'',
                    'count_key':'count',
                    'yerr_key':'',
                    'xs_ratio':1.0,
                    'lumi_ratio':0.0,
                },
            )

        #---------- Set configurations ----------#
        # Setup configuration dictionary
        # Create job submission structure
        asyms = [-0.1,-0.01,0.0,0.01,0.1]
        asymfitvars = {"asymfitvars":["costheta1","costheta2","costhetaT","costhetaTy"]}
        sgasyms = {"sgasyms":[[a1] for a1 in asyms]}
        bgasyms = {"bgasyms":[[a1] for a1 in asyms]}
        seeds   = {"inject_seed":[2**i for i in range(16)]}

        # Set job file paths and configs
        configs = dict(
            asymfitvars,
            **sgasyms,
            **bgasyms,
            **seeds
        )

        #TODO: ALLOW OPTIONS FOR DIFFERENT CONFIGS

        # Get list of configurations
        config_list = sagas.get_config_list(configs,aggregate_keys=aggregate_keys)

        # Get aggregated list of directories
        out_dirs_list = sagas.get_out_dirs_list(
                                        configs,
                                        base_dir,
                                        aggregate_keys=aggregate_keys
                                    )

        #---------- Loop bin schemes ----------#
        for binscheme_idx, binscheme_name in enumerate(binschemes.keys()):

            # Get the bin scheme
            binscheme = binschemes[binscheme_name]
            proj_var  = list(binscheme.keys())[0] #NOTE: Assume projection variable is the only variable in the bin scheme
            nbins = len(binscheme[proj_var])-1

            # Arguments for sagas.get_graph_data()
            xvar_keys = [proj_var]

            # Set some bin scheme dependent plotting parameters
            binlims = binscheme[proj_var]
            plot_results_kwargs_base['xlims'] = [binlims[0],binlims[-1]]
            plot_results_kwargs_base['xlabel'] = xlabel_map[binscheme_name]
            #plot_results_kwargs_base['binlims'] = binlims
            #plot_results_kwargs_base['hist_paths'] = hist_paths_map[binscheme_name]

            # Get the bin migration path
            bin_mig_path = sagas.get_out_file_name(
                base_dir=base_dir,
                base_name=bin_mig_base_name,
                binscheme_name=binscheme_name,
                ext=out_file_name_ext
            )

            # Load bin migration matrix and invert
            bin_mig_df, bin_mig_mat, inv_bin_mig_mat = None, None, None
            if use_bin_mig:
                bin_mig_df = load_csv(bin_mig_path)
                bin_mig_mat = sagas.get_bin_mig_mat(
                    bin_mig_df,
                    id_gen_key=id_gen_key,
                    id_rec_key=id_rec_key,
                    mig_key=mig_key,
                )
                save_bin_mig_mat_to_csv(
                    bin_mig_mat,
                    base_dir='./',
                    basename=binscheme_name,
                    delimiter=",",
                    header=None,
                    fmt=None,
                    comments='',
                )
                inv_bin_mig_mat = np.linalg.inv(bin_mig_mat)

            #---------- Loop configurations ----------#
            # Loop each aggregate list
            for config_idx in range(len(config_list)):

                # Set the config you are interested in
                config = config_list[config_idx]
                out_dirs = out_dirs_list[config_idx]

                # Set the output path basename for this config
                config_out_path = sagas.get_config_out_path(
                        base_dir,
                        aggregate_keys,
                        binscheme_name+sep+rg+sep+result_name,
                        config,
                        sep=sep,
                        ext=ext,
                    )
                config_out_path = os.path.join(base_dir,config_out_path)

                # Get the name of the CSV file for the binning scheme you are interested in
                out_file_names = [sagas.get_out_file_name(
                        base_dir=outdir,
                        base_name='out_',
                        binscheme_name=binscheme_name,
                        ext=out_file_name_ext
                    ) for outdir in out_dirs]

                # Load pandas dataframes from the files
                dfs = [load_csv(out_file_name,config=config,chain_configs=chain_configs) for out_file_name in out_file_names]

                # Apply bin migration correction
                if use_bin_mig:
                    for df in dfs:
                        sagas.apply_bin_mig(df,inv_bin_mig_mat,results_keys=results_keys) #NOTE: THIS MODIFIES THE DATAFRAMES IN PLACE

                # Get an aggregate graph
                proj_ids = [i for i in range(nbins)]#NOTE: Assume bin scheme indices are simple
                sgasym_idx = plot_results_kwargs_base['sgasym_idx'] #NOTE: Assume this is in the kwargs base dictionary
                sgasym = config['sgasyms'][sgasym_idx] if 'sgasyms' in config else 0.0
                plot_results_kwargs_base['sgasyms'] = config['sgasyms']#NOTE: TODO:
                aggregate_graph = sagas.get_aggregate_graph(
                    [
                        sagas.get_graph_data(
                                    df,
                                    proj_ids,
                                    id_key=id_key,
                                    count_key=count_key,
                                    xvar_keys=xvar_keys,
                                    asym_key=asym_key,
                                    err_ext=err_ext
                        ) for df in dfs
                    ],
                    xvar_keys=xvar_keys,
                    sgasym=sgasym
                )

                # Use default plotting settings
                if use_default_plt_settings: set_default_plt_settings()

                # Create figure and axes
                f, ax = plt.subplots(figsize=figsize)

                # Set additional arguments for saga.plot.plot_results()
                plot_results_kwargs_base['sgasyms'] = config['sgasyms']
                plot_results_kwargs_base['outpath'] = config_out_path

                # Plot the graph
                plot_results(ax,**aggregate_graph,**plot_results_kwargs_base)

                # Save the graph
                f.savefig(config_out_path)

                # Close the graph
                plt.close()
