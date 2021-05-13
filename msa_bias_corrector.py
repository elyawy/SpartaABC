# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 13:19:36 2021

@author: gillo
"""
#imports
import subprocess
import os
import logging
logger = logging.getLogger(__name__)
import re
import pandas as pd
import numpy as np
from sklearn import linear_model 
from sklearn import model_selection
from scipy.stats import pearsonr




#%%

#TODO - add log

def parse_alignments_file(real_alignments_path):
	"""
	Parses aligments file to list of fasta.
	Returns this list and the maximal sequence
	"""
	with open(real_alignments_path,'r') as f:
		lines = f.read()
	align_list = lines.split('\n\n')[:-1] # do not change on server
	max_sim_seq_len = len(max(lines.split('\n'),key=len))
	return align_list, max_sim_seq_len

def prepare_indelible_control_file(res_path,tree_filename,indelible_out_file_name,
								   num_msa,max_sim_seq_len,pipeline_path,
								   indelible_template_file_name,
								   submodel_params):
	"""
	prepare indelible control file
	"""
	with open(pipeline_path+indelible_template_file_name,'r') as f:
		indelible_temp = f.read()
	with open(res_path+tree_filename,'r') as f:
		tree = f.read().rstrip()


	if submodel_params["mode"] == "amino":
		mode_to_replace = '[TYPE] AMINOACID 2'
		mode_str = '[TYPE] AMINOACID 2'
		submodel_to_replace = '[submodel] WAG'
		submodel_str = '[submodel] WAG'
		freq_to_replace = "[statefreq]   0.25  0.25  0.25  0.25"
		freqs_str = ""
		rates_to_replace = "[rates]         0.25 0.50 10"
		rates_str = ""
	else:
		if submodel_params["submodel"] == "GTR":
			sub_rates_str = f"{submodel_params['rates'][0]:.9f} {submodel_params['rates'][1]:.9f} {submodel_params['rates'][2]:.9f} {submodel_params['rates'][3]:.9f} {submodel_params['rates'][4]:.9f}"
			freqs_str = f"{submodel_params['freq'][0]:.6f} {submodel_params['freq'][1]:.6f} {submodel_params['freq'][2]:.6f} {submodel_params['freq'][3]:.6f}"
			rates_str = f"{submodel_params['inv_prop']} {submodel_params['gamma_shape']} {submodel_params['gamma_cats']}"

			
			mode_to_replace = '[TYPE] AMINOACID 2'
			mode_str = '[TYPE] NUCLEOTIDE 2'

			submodel_to_replace = '[submodel] WAG'
			submodel_str = f'[submodel] {submodel_params["submodel"]} {sub_rates_str}'

			freq_to_replace = "[statefreq]   0.25  0.25  0.25  0.25"
			freqs_str = f"[statefreq] {freqs_str}"

			rates_to_replace = "[rates]         0.25 0.50 10"
			rates_str = f"[rates] {rates_str}"
		else:
			mode_to_replace = '[TYPE] AMINOACID 2'
			mode_str = '[TYPE] NUCLEOTIDE 2'

			submodel_to_replace = '[submodel] WAG'
			submodel_str = f'[submodel] {submodel_params["submodel"]}'
			
			freq_to_replace = "[statefreq]   0.25  0.25  0.25  0.25"
			freqs_str = ""
			rates_to_replace = "[rates]         0.25 0.50 10"
			rates_str = ""



	tree_to_replace = '(A:0.1,B:0.1);'
	num_sim_str_to_replace = 'partitionname1 2000 outputname1'
	num_sim_str = f'partitionname1 {num_msa} outputname1'
	seq_len_to_replace = 'treename1 modelname 3000'
	seq_len_str = f'treename1 modelname {max_sim_seq_len}'
	
	with open(res_path+indelible_out_file_name,'w') as f:
		out_str = indelible_temp.replace(tree_to_replace,tree)
		out_str = out_str.replace(num_sim_str_to_replace,num_sim_str)
		out_str = out_str.replace(seq_len_to_replace, seq_len_str)
		out_str = out_str.replace(mode_to_replace, mode_str)
		out_str = out_str.replace(submodel_to_replace, submodel_str)
		out_str = out_str.replace(freq_to_replace, freqs_str)
		out_str = out_str.replace(rates_to_replace, rates_str)

		f.write(out_str)
		
def run_indelible(res_path,logger=None):
	"""
	runs indelible.
	Requires control.txt at res_path and indelible command
	"""
	os.chdir(res_path)
	cmd = "indelible"
	if logger!=None:
		logger.info(f'Starting indelible. Executed command is:\n{cmd}')
	subprocess.run(cmd, shell=True)
	indelible_msa_list = parse_indelible_output(res_path)
	# clean indelible files
	os.remove(f'{res_path}outputname1_TRUE.phy')
	os.remove(f'{res_path}outputname1.fas')
	os.remove(f'{res_path}trees.txt')
	os.remove(f'{res_path}LOG.txt')

	return indelible_msa_list

def parse_indelible_output(res_path):
	"""
	reads the output of indelible and parse it to list of msas
	"""
	with open(res_path+'outputname1.fas','r') as f:
		indelible_subs = f.read()
	indelible_msa_list = re.split('\n *\n', indelible_subs)
	return indelible_msa_list[:-1]

def prepare_sparta_conf_sumstat(res_path, pipeline_path, sum_stat_file_name='tmp_sum_stat.csv',msa_filename='realigned_msa_tmp.fasta',conf_filename_out='sum_stat.conf',conf_file_template='sparta_conf_template.conf'):
	"""
	prepare a configuration file for summary stats only
	of input msa
	"""
	with open(f'{pipeline_path}{conf_file_template}','r') as f:
		conf_str = f.read()

	conf_str = conf_str.replace('_indelibleTemplateControlFile control_indelible_template.txt',
								f'_indelibleTemplateControlFile {pipeline_path}control_indelible_template.txt')
	conf_str = conf_str.replace('_only_real_stats 0', '_only_real_stats 1')
	conf_str = conf_str.replace('SpartaABC_data_name_model_name.posterior_params',sum_stat_file_name)
	tmp_ind1 = conf_str.find('_inputRealMSAFile')
	tmp_ind2 = tmp_ind1 + conf_str[tmp_ind1:].find('\n')
	conf_str = conf_str.replace(conf_str[tmp_ind1:tmp_ind2],f'_inputRealMSAFile results/{msa_filename}')
	conf_str = conf_str.replace('results/', res_path)
	
	with open(f'{res_path}{conf_filename_out}','w') as f:
		conf_str = f.write(conf_str)

def process_raw_msa(raw_msa):
	split_msa = raw_msa.strip().split('\n')
	return split_msa[1::2],split_msa[::2]

def restructure_msa(msa_list, organism_list):
	new_msa = []
	for i in range(len(organism_list)):
		new_msa.append(organism_list[i])
		new_msa.append(msa_list[i])

	return "\n".join(new_msa)

def add_subs_to_sim_msa(raw_sim_msa,indelible_msa):
	"""
	add AA from indelible output file to the simulated alignment
	in all non-gapped locations
	"""
	sim_msa_list, organism_list = process_raw_msa(raw_sim_msa)
	indelible_msa_list = process_raw_msa(indelible_msa)[0]

	sub_sim_msa_list = []
	unaligned_sub_sim_msa_list = []
	for i,alignment in enumerate(sim_msa_list):
		merged_alignment = ""
		for j,c in enumerate(alignment):
			if c=='-':
				merged_alignment += "-"
			else:
				merged_alignment += indelible_msa_list[i][j]
		sub_sim_msa_list.append(merged_alignment)
		unaligned_sub_sim_msa_list.append(merged_alignment.replace('-',''))
	
	sub_sim_msa_list = restructure_msa(sub_sim_msa_list, organism_list)

	unaligned_sub_sim_msa = restructure_msa(unaligned_sub_sim_msa_list, organism_list)
	return unaligned_sub_sim_msa, sub_sim_msa_list

def restructure_mafft_output(mafft_output):
	restructured_output = ""
	restructured_temp = mafft_output.split(">")[1:]
	for item in restructured_temp:
		organism, sequence = item.split("\n", 1)
		sequence = sequence.replace("\n",'')
		restructured_output += f">{organism}\n{sequence}\n"

	return restructured_output

def reconstruct_msa(res_path, unaligned_msa, output_name,  align_mode,logger=None):

	tmp_file = "temp_unaligned.fasta"
	
	with open(f'{res_path}{tmp_file}', 'w') as f:
		f.write(unaligned_msa)

	cmd = f'mafft --auto --{align_mode} {res_path+tmp_file}'
	
	if logger!=None:
		logger.info(f'Starting MAFFT. Executed command is:\n{cmd}')
	results = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout.decode('utf-8')
	
	results = restructure_mafft_output(results)
	# TODO: remove file safely
	os.remove(res_path+tmp_file)
	return results
	
def run_sparta_sum_stat(input_msa, pipeline_path, conf_file_path):
	with open(conf_file_path, 'r') as f:
		msa_path = next(filter(lambda x: "_inputRealMSAFile" in x, f.readlines()))
	msa_path = msa_path.split(" ")[1].rstrip("\n")
	with open(msa_path,'w') as f:
		f.write(input_msa)
	#TODO in linux - remove exe
	cmd = f'{pipeline_path}SpartaABC {conf_file_path}'
	subprocess.run(cmd, shell=True)
	os.remove(msa_path)
	
def load_sim_res_file(sim_res_file_path):

	with open(sim_res_file_path) as f:
	  file_rows_num = sum(1 for line in f)
	df_real = pd.read_csv(sim_res_file_path, delimiter='\t',skiprows=[i for i in range(1,4)],nrows=(file_rows_num-11))
	df_meta = pd.read_csv(sim_res_file_path, delimiter='\t', nrows=2)
	return df_real, df_meta   
	
def lasso_reg(X,y):
	# X,y? should be normalized
	reg = linear_model.Lasso(normalize=False,fit_intercept=True)
	parameters = {'alpha':np.logspace(-7,4,20)}
	clf = model_selection.GridSearchCV(estimator = reg,
							   param_grid = parameters, cv=3,
							   scoring = 'neg_mean_squared_error')
	clf.fit(X,y)
	cv_rmse = np.min(np.sqrt(-clf.cv_results_['mean_test_score']))
	
	res = clf.predict(X)
	out_stat = (pearsonr(res,y),cv_rmse)
	return clf,out_stat
	
def correct_mafft_bias(res_path, sim_res_file_path, df_mafft, num_msa,model_type, filter_p, alignment_flag):

	df_real, df_meta = load_sim_res_file(sim_res_file_path)
	sstat_cols = list(df_mafft.columns)[6:]
	params_cols = list(df_mafft.columns)[1:6]
	train_cols =  params_cols+sstat_cols # sstat_cols
	if alignment_flag:
		X = df_real.iloc[int(num_msa/2):num_msa][train_cols].values
	else:
		X = df_real.iloc[:num_msa][train_cols].values
	X_full = df_real[train_cols].values
	Y = df_mafft[sstat_cols].values
	X_train = X
	X_train_mean = X_train.mean(axis=0)
	X_train_std = X_train.std(axis=0) 
	Y_train = Y
	epsilon = 1E-4
	X_train_reg = (X-X_train_mean+epsilon)/(X_train_std+epsilon)
	X_full_reg = (X_full-X_train_mean+epsilon)/(X_train_std+epsilon)
	X_full_reg[np.isnan(X_full_reg)] = 0
	X_full_reg[X_full_reg == 10000000] = 0

	
	df_trans = df_real.copy()
	for ind,p in enumerate(sstat_cols):
		if alignment_flag:
			y = Y_train[int(num_msa/2):num_msa,ind]
		else:
			y = Y_train[:,ind]
		clf,out_stat = lasso_reg(X_train_reg,y)
		df_trans[p] = clf.predict(X_full_reg)
	

	df_trans_subset = df_trans.iloc[:num_msa,:]
	min_num_sumstat = filter_p[1]
	correction_th = filter_p[0]

	# create correlation filter
	msa_correct_qual_dict = {}
	for col in sstat_cols:
		msa_correct_qual_dict[col] = pearsonr(df_mafft[col],df_trans_subset[col])[0]
	
	logger.info(msa_correct_qual_dict)

	sumstat_to_use = [x for x in msa_correct_qual_dict if msa_correct_qual_dict[x]>=correction_th]
	with open(res_path + f"used_features_{model_type}.txt", 'w') as f:
		f.write("\n".join(sumstat_to_use))


	if len(sumstat_to_use) < min_num_sumstat:
		sumstat_to_use = sorted(msa_correct_qual_dict,key=msa_correct_qual_dict.get,reverse=True)[:min_num_sumstat]

	sumstat_to_drop = [x for x in sstat_cols if x not in sumstat_to_use]
	sstat_cols_inds = [i for i,x in enumerate(df_trans.columns) if x in sstat_cols]
	sstat_cols_dict = dict(zip(sstat_cols,sstat_cols_inds))


	# Calculating new weights
	n_weights = 10000
	# TODO: figure out why std_dev is 0 for some inputs
	std_dev = df_trans.iloc[-n_weights:].std().apply(lambda x: x if x > 0.001 else 10**5) # hack

	weights = 1/(std_dev) # check - should be 1/sigma. check in cpp if indeed this the way (or squared)
	# Check - not sure if correct
	for i in sumstat_to_drop:
		weights.at[i] = 0


	df_trans['DISTANCE'] = np.sum(((df_meta.iloc[0][sumstat_to_use].values.reshape(1,-1).T-df_trans[sumstat_to_use].values.T)*weights[sumstat_to_use].values.reshape(-1,1))**2,axis=0)**0.5
	df_trans_string = df_trans.to_csv(index=False, header=False, sep='\t', float_format='%.6f')



	df_meta2 = df_meta.copy()
	df_meta2.loc[1,sstat_cols] = weights
	df_meta2_string = df_meta2.to_csv(index=False, header=False, sep='\t', float_format='%.6f')

	df_head_string = df_meta2.head(0).to_csv(index=False, header=True, sep='\t', float_format='%.6f')

	with open(sim_res_file_path,'r') as f:
		string_tmp = f.readlines()

	out_str = df_head_string+"\n"+df_meta2_string+df_trans_string+"".join(string_tmp[-7:])
	out_str = out_str.replace('\r','')
	del string_tmp
	file_name_out = f'{res_path}SpartaABC_msa_corrected_id{model_type}.posterior_params'
	with open(file_name_out,'w') as f:
		f.write(out_str)
		


def continuous_write(interation, file_path, to_write):
	with open(file_path, ('w' if interation==0 else 'a')) as f:
			if interation > 0:
				f.write("\n\n")
			f.write(to_write)

def remove_large_files(res_path,to_remove):
	for i in to_remove:
		os.remove(f'{res_path}{i}')

def msa_bias_correction(skip_config, clean_run, res_path,
				 real_alignments_filename,tree_filename,
				 pipeline_path,indelible_template_file_name,
				 model_type,filter_p,submodel_params,
				 indelible_out_file_name='control.txt'):

	align_list, max_sim_seq_len = parse_alignments_file(res_path+real_alignments_filename)
	logger.info(f'Maximal sequence length for model {model_type}: {max_sim_seq_len}')
	
	num_msa = len(align_list)
	logger.info(f'Number of simulated MSAs  for model {model_type}: {max_sim_seq_len}')
	
	df_mafft = None
	if skip_config["mafft"]:

		prepare_indelible_control_file(res_path,tree_filename,indelible_out_file_name,
									   num_msa,max_sim_seq_len,pipeline_path,
									   indelible_template_file_name, submodel_params)
		indelible_msa_full_list = run_indelible(res_path)

		with open(f"sparta_aligned_{model_type}.fasta",'w') as f:
			f.write("\n\n".join(indelible_msa_full_list))
		
		logger.info(f'Number of indelible MSAs for model {model_type}: {len(indelible_msa_full_list)}')

		prepare_sparta_conf_sumstat(res_path, pipeline_path)

		realigned_msa_tmp_filename = 'realigned_msa_tmp.fasta'
		# use indelible simul results to replace sparta alignment res.
		print("Running MAFFT...")
		for i in range(num_msa):
			raw_sim_msa = align_list[i]
			indelible_msa = indelible_msa_full_list[i]

			unaligned_sub_sim_msa, indelible_sparta_msa = add_subs_to_sim_msa(raw_sim_msa,indelible_msa)
			# write all unaligned msas to file.
			continuous_write(interation=i,
							file_path=f'{res_path}{f"all_unaligned_sims_{model_type}.txt"}',
							to_write=unaligned_sub_sim_msa)
			continuous_write(interation=i,
							file_path=f'{res_path}{f"indelible_sparta_{model_type}.txt"}',
							to_write=indelible_sparta_msa)

			# run mafft on unaligned sequences.
			realigned_msa = reconstruct_msa(res_path=res_path, 
											unaligned_msa=unaligned_sub_sim_msa,
											output_name=realigned_msa_tmp_filename,
											align_mode=submodel_params["mode"],
											logger=None)
			# write all realigned msas to file.
			continuous_write(interation=i,
							file_path=f'{res_path}{f"all_realigned_sims_{model_type}.txt"}',
							to_write=realigned_msa)
			
			run_sparta_sum_stat(input_msa=realigned_msa,
								pipeline_path=pipeline_path,
								conf_file_path=res_path+'sum_stat.conf')
			
			df_tmp = pd.read_csv(f'{res_path}tmp_sum_stat.csv',delimiter='\t')
			df_mafft = df_tmp if i==0 else pd.concat([df_mafft,df_tmp],ignore_index=True) # TODO: save as separate file.
			os.remove(f'{res_path}tmp_sum_stat.csv')
		os.remove(f'{res_path}sum_stat.conf')
		os.remove(f'{res_path}control.txt')
		print("Done.")
		df_mafft.to_csv(f"mafft_sum_stats_{model_type}.csv", sep="\t", index=False)
		logger.info(f'Done with MAFFT')
	else:
		logger.info("Skipping Mafft.")
		if df_mafft is None:
			try:
				df_mafft = pd.read_csv(res_path+f"mafft_sum_stats_{model_type}.csv", sep="\t")
			except Exception:
				logging.error("No mafft results found, please provide mafft_sum_stats file")
				return
	sim_res_file_path = f'{res_path}SpartaABC_data_name_id{model_type}.posterior_params'
	correct_mafft_bias(res_path,sim_res_file_path,df_mafft, num_msa,model_type,filter_p, alignment_flag=False)
	if clean_run:
		remove_large_files(res_path,to_remove=[
			f"all_realigned_sims_{model_type}.txt",
			f"all_unaligned_sims_{model_type}.txt",
			f"alignments_{model_type}.fasta",
			f"sparta_aligned_{model_type}.fasta",
			f"indelible_sparta_{model_type}.txt",
			f"mafft_sum_stats_{model_type}.csv",
			f"SpartaABC_data_name_id{model_type}.posterior_params",
			f"used_features_{model_type}.txt"
		])
	logger.info(f'Corrected MSA bias for model {model_type}')
	
#%%
def apply_correction(skip_config, clean_run,res_path, pipeline_path, tree_file,filter_p, submodel_params="amino"):

	res_path = res_path
	# The following files should be on the res_path 

	tree_filename = tree_file

	pipeline_path = pipeline_path
	# The following files should be on the pipeline_path
	indelible_template_file_name = 'control_indelible_template.txt'

	for model_type in ['dif',"eq"]:
		real_alignments_filename = f'alignments_{model_type}.fasta'
		msa_bias_correction(skip_config, clean_run,res_path,
						real_alignments_filename,tree_filename,
						pipeline_path,indelible_template_file_name,
						model_type,filter_p,submodel_params,indelible_out_file_name='control.txt')




