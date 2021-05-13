# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:43:59 2020

@author: gillo
"""

import os
import subprocess
import shutil
import logging
logger = logging.getLogger(__name__)

from Bio import AlignIO


def phylip2fasta(data_path,dry_run_flag=False, verbose=1):
	if dry_run_flag:
		if verbose:
			print('dry run: entered phylip2fast')
			logger.info(data_path+' dry run: entered phylip2fast')
		return

	AlignIO.convert(os.path.join(data_path,'ref_msa.aa.phy'), 'phylip-relaxed', os.path.join(data_path,'ref_msa.aa.fasta'),'fasta')
	logger.info(data_path + " Done convert phylip to fasta")
	if verbose:
		print("Done convert phylip to fasta")
	return

def run_sparta_abc(exe_path,conf_path):
	args = exe_path + ' ' + conf_path
	tmp = subprocess.call(args,shell=True)
	return tmp



def create_sims(data_name,verbose=1,res_dir='results/',
				data_dir='data',msa_filename='ref_msa.aa.fasta',
				tree_filename='RAxML_tree.tree',
				minIR = 0, maxIR = 0.05, num_alignments=200,
				cwd='/groups/pupko/gilloe/spartaABC/code/abc_nn/',
				op_sys='linux'):
	logger.info('Inside create_sims function.')

	out_dir = res_dir + data_name

	model_types = ["eq", "dif"]#["eq", "a_dif", "r_dif", "dif"]

	logger.info('Writing conf files.')
	for model in model_types:
		with open(cwd+'sparta_conf_template.conf', "rt") as fin:
			with open(out_dir+f'_{model}.conf', "wt") as fout:
				for line in fin:
					if line.startswith('_outputGoodParamsFile'):
						line_out = f'_outputGoodParamsFile {out_dir}SpartaABC_data_name_id{model}.posterior_params\n'
					elif line.startswith('_outputAlignmnetsFile'):
						line_out = f'_outputAlignmnetsFile {out_dir}alignments_{model}.fasta\n'
					elif line.startswith('_alignments_output'):
						line_out = f'_alignments_output {num_alignments}\n'
					else:
						line_out = line.replace('results/', res_dir).replace('data/', data_dir).replace('data_name/', data_name).replace('model_name','ideq').replace('ref_msa.aa.fasta',msa_filename).replace('RAxML_tree.tree',tree_filename).replace('_minIRVal 0',f'_minIRVal {round(minIR,2)}').replace('_maxIRVal 0.05',f'_maxIRVal {round(maxIR,2)}').replace('_modelType eq', f'_modelType {model}')
					fout.write(line_out)
		logger.info(f"Wrote {out_dir+f'_{model}.conf'}. Based on {cwd+'sparta_conf_template.conf'}")
			
	stat = "didn't run"
	if op_sys=='linux':
		logger.info(f"linux op system")
		for model in model_types:
			stat = run_sparta_abc(exe_path=cwd+'SpartaABC', conf_path=res_dir+data_name+data_name+f'_{model}.conf')
			logger.info(f"ran {cwd+f'SpartaABC'}, conf_path={res_dir+data_name+data_name+f'_{model}.conf'}, stat={stat}")
	else: #windows
		logger.info(f"windows op system")
		stat_eq = run_sparta_abc(exe_path=cwd+'SpartaABC', conf_path=res_dir+data_name+data_name+'_ideq.conf')
		logger.info(f"ran {cwd+'SpartaABC'}, conf_path={res_dir+data_name+data_name+'_ideq.conf'}")
		stat_dif = run_sparta_abc(exe_path=cwd+'SpartaABC', conf_path=res_dir+data_name+data_name+'_iddif.conf')
		logger.info(f"ran {cwd+'SpartaABC'}, conf_path={res_dir+data_name+data_name+'_iddif.conf'}")
	logger.info(f'{data_name} simulations are done - finish stats {stat}.')
	if verbose:
		print(f'{data_name} simulations are done - finish stats {stat}.')
	return stat




def create_sims_from_data(data_name,ow_flag=False,verbose=1,
						  res_dir='results',data_dir='data',
						  msa_filename='ref_msa.aa.fasta',
						  tree_filename='RAxML_tree.tree',
						  minIR = 0, maxIR = 0.05, num_alignments=200,
						  cwd='/groups/pupko/gilloe/spartaABC/code/abc_nn/',
						  op_sys='linux'):

	if not res_dir.endswith('/'):
		res_dir += '/'
	if not data_dir.endswith('/'):
		data_dir += '/'

	logger.info(data_name +' started.')
	stat_eq = stat_dif = "didn't run"
	if verbose:
		print(data_name +' started.')
	try:
		create_sims(data_name=data_name,verbose=verbose,
					res_dir=res_dir,data_dir=data_dir,
					msa_filename=msa_filename,
					tree_filename=tree_filename,
					minIR=minIR,maxIR=maxIR, num_alignments=num_alignments,
					cwd=cwd,op_sys=op_sys)
	except:
		if ow_flag:
			shutil.rmtree(os.getcwd() + '/'+res_dir + data_name)
			create_sims(data_name=data_name,verbose=verbose,
						res_dir=res_dir,data_dir=data_dir,
						msa_filename=msa_filename,
						tree_filename=tree_filename,
						minIR=minIR,maxIR=maxIR, num_alignments=num_alignments,
						cwd=cwd,op_sys=op_sys)
		else:
			logger.info(data_name +' skipped.')
			if verbose:
				print(data_name +' skipped.')