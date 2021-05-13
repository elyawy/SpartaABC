# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 13:02:55 2020

@author: gillo
"""

import click
import os
import logging

# set to environment variable 'DEBUG' to 1 for full list of runtime parameters.
DEBUG = True if os.environ.get("DEBUG") == "0" else False


pipeline_path = os.path.dirname(os.path.abspath(__file__))
pipeline_path = pipeline_path.replace('\\','/')
if pipeline_path[-1]!='/':
	pipeline_path = pipeline_path + '/'

os.chdir(pipeline_path)

#%%


def pipeline(skip_config, res_dir, clean_run,msa_filename,tree_filename,pipeline_path='/bioseq/spartaabc/pipeline/',
			 minIR=0,maxIR=0.05,op_sys='linux',verbose=0, filter_p=(0.9,15),
			 b_num_top=100, num_alignments=200, submodel_params="amino"):
	
	res_dir = os.path.join(res_dir, '')
	if os.path.isdir(res_dir + "logs/"):
		print("using existing log directory")
	else:
		os.mkdir(res_dir + "logs/")
		print("creating log directory")
	log_dir = res_dir + "logs/"
	log_id = pipeline_path.split(r"/")[-2]
	import sys
	if verbose!=2:
		if not sys.warnoptions:
			import warnings
			warnings.simplefilter("ignore")
	logging.basicConfig(filename=log_dir+log_id+'.log',level=logging.INFO,
					format='%(asctime)s %(name)s %(levelname)-8s %(message)s',
					datefmt='%Y-%m-%d %H:%M:%S')  
	logger = logging.getLogger(__name__)

	
	import tree_cleaner as treec
	

	import summarize_results as sumres
	
	treec.fix_tree(res_dir, tree_filename)

	if skip_config["sparta"]:
		import run_sparta_abc_single_folder_pipeline as runs

		runs.create_sims_from_data(data_name='', ow_flag=False,
							verbose=verbose, res_dir=res_dir,
							data_dir=log_dir,
							msa_filename=msa_filename,
							tree_filename=tree_filename,
							minIR=minIR,maxIR=maxIR, num_alignments=num_alignments,
							cwd=pipeline_path, 
							op_sys=op_sys) # Running spartaABC C++ code to get simulations (through another python script)
		if clean_run:
			os.remove(f'{res_dir}_eq.conf')
			os.remove(f'{res_dir}_dif.conf')
	else:
		logger.info("Skipping Sparta.")
		#check if sparta params file exists
		for model in ["eq", "dif"]:
			if os.path.isfile(res_dir + f'SpartaABC_data_name_id{model}.posterior_params'):
				print("retrieved existing params files.")
				logger.info("retrieved existing params files.")
			else:
				print("Could not find SpartaABC params file.")
				logger.error("Could not find SpartaABC params file.\nPlease provide the params files or run the full pipeline")
				return

	if skip_config["correct_bias"]:
		import msa_bias_corrector as corrector

		corrector.apply_correction(skip_config=skip_config,res_path=res_dir, clean_run=clean_run,
								pipeline_path=pipeline_path,
								tree_file=tree_filename,filter_p=filter_p, submodel_params=submodel_params)
	else:
		logger.info("Skipping msa bias correction.")
		#check if sparta params file exists
		for model in ["eq", "dif"]:
			if os.path.isfile(res_dir + f'SpartaABC_msa_corrected_id{model}.posterior_params'):
				print("retrieved corrected param files.")
				logger.info("retrieved corrected param files.")
			else:
				print("Could not find corrected params file.")
				logger.error("Could not find corrected params file.\nPlease provide the param files or run without the --skip-bc option")
				return	

	if skip_config["inference"]:
		import infer_abc_params_single_folder_pipeline as sinf

		sinf.calc_stats(csv_out_path=res_dir,lib='msa_corrected',path='' ,
				verbose = verbose , models_list=['ideq','iddif'],
				b_num_top=b_num_top) # Inferring parameters from the C++ simulations
	else:
		logger.info("Skipping inference step.")
		#check if sparta params file exists
		if os.path.isfile(res_dir + 'msa_corrected_res.csv'):
			print("retrieved inference results.")
			logger.info("retrieved inference results.")
		else:
			print("Could not find inference results.")
			logger.error("Could not find inference results.\nPlease provide the param files or run without the --skip-i option")
			return	
		
	sumres.get_stats_v2(results_file_path=res_dir,file_name='msa_corrected_res.csv',
						minIR=minIR, maxIR=maxIR, minAI=0, maxAI=2, msa_path=res_dir+msa_filename,
						clean_run=clean_run,verbose=verbose)
	
	return



@click.command()
@click.option('--path', help='Path of the phylogenetic tree and MSA (output will be saved to this path too).', required=True)
@click.option('--msaf', help='MSA filename', required=True)
@click.option('--trf', help='Phylogenetic tree filename (Newick without bootstrap)', required=True)
@click.option('--ver', default=0 ,help='Verbosity (0/1/2) Default: 1')
@click.option('--minr', default=0.0 ,help='Minimal INDEL rate. Default: 0')
@click.option('--maxr', default=0.05 ,help='Maximal INDEL rate. Default: 0.05')
@click.option('--bn', default=100 ,help='epsilon num to use for ABC. Default: 100')
@click.option('--numalign', default=200 ,help='Number of alignments for MSA bias correction')
@click.option('--skip-s', default='True', is_flag=True ,help='Skip SpartaABC simulation', hidden=DEBUG)
@click.option('--skip-m', default='True', is_flag=True ,help='Skip Mafft alignment', hidden=DEBUG)
@click.option('--skip-i', default='True', is_flag=True ,help='Skip inference step', hidden=DEBUG)
@click.option('--skip-bc', default='True', is_flag=True ,help='Skip MSA bias correction', hidden=DEBUG)
@click.option('--filterp', default=(0.9,15) ,help='MSA bias correction filtering parameter. Default: 0.9 15', hidden=DEBUG)
@click.option('--nonclean-run', default='True', is_flag=True ,help='Do not clean files at runtime', hidden=DEBUG)
@click.option('--mode', type=click.Choice(['amino', 'nuc']), required=True, help='Specify type of alignment, proteins or DNA.')
@click.option('--submodel', type=click.Choice(['JC', 'GTR']), default='JC', help='Specify substitution model.')
@click.option('--freq', default=(0.25,0.25,0.25,0.25), help='Specify rate parameters.')
@click.option('--rates', default=(1.0,1.0,1.0,1.0,1.0), help='Specify rate parameters.')
@click.option('--inv-prop', default=0.25, help='Specify invariable sites proportion.')
@click.option('--gamma-shape', default=0.50, help='Specify shape parameter for the gamma distribution.')
@click.option('--gamma-cats', default=10, help='Specify number of categories to use in the discrete gamma approximation.')
def pipeline_click(path,msaf,trf,ver,minr,maxr, bn,
				   numalign,
				   skip_s, skip_m, skip_i, skip_bc,
				   filterp, nonclean_run, 
				   mode, submodel, freq, rates, inv_prop, gamma_shape, gamma_cats):
	
	skip_config = {
		"sparta": skip_s and skip_m and skip_bc and skip_i,
		"mafft": skip_m and skip_i,
		"inference": skip_i ,
		"correct_bias": skip_bc and skip_i
	}

	submodel_params_ = {
		"mode": mode,
		"submodel": submodel,
		"freq": freq,
		"rates": rates,
		"inv_prop": inv_prop,
		"gamma_shape": gamma_shape,
		"gamma_cats": gamma_cats
	}

	verbose = ver
	res_dir = path #results path
	clean_run = nonclean_run
	op_sys= 'linux'
	msa_filename = msaf
	tree_filename = trf
	minIR=minr
	maxIR=maxr
	verbose = ver
	b_num_top=bn
	num_alignments = numalign
	filter_p = filterp

	pipeline(skip_config=skip_config,
			 pipeline_path=pipeline_path,
			 res_dir=res_dir, clean_run=clean_run,
			 msa_filename=msa_filename,
			 tree_filename=tree_filename,
			 minIR=minIR,maxIR=maxIR,
			 op_sys=op_sys,verbose=verbose,
			 b_num_top=b_num_top,
			 num_alignments=num_alignments,
			 filter_p=filter_p, submodel_params=submodel_params_)

	
#%% main
		
if __name__ == '__main__':
	pipeline_click()

