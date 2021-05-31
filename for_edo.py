from pipeline_click import pipeline
import random


num_of_samples = int(1e5)
num_of_samples_per_file = int(1e5)

print(num_of_samples)

for sample in range(num_of_samples):
	branch_length = random.uniform(0.015,0.5)
	seq_length = random.randint(100,300)
	tree = f'(A:{branch_length:.4f},B:{branch_length:.4f});'

	with open("/home/elyawy/development/Msc/Thesis/Working_dir/SpartaABC/edo/temp_tree.tree",'w') as f:
		f.write(tree)

	with open("/home/elyawy/development/Msc/Thesis/Working_dir/SpartaABC/edo/temp_msa.fasta",'w') as f:
		f.write(">A\n")
		f.write("T"*seq_length + "\n")
		f.write(">B\n")
		f.write("T"*seq_length + "\n")


	skip_config = {
		"sparta": True,
		"mafft": True,
		"inference": True ,
		"correct_bias": True
	}

	submodel_params_ = {
		"mode": "nuc",
		"submodel": "GTR",
		"freq": (0.369764, 0.165546, 0.306709, 0.157981),
		"rates": (0.443757853, 0.084329474, 0.115502265, 0.107429571, 0.000270340),
		"inv_prop": 0.0,
		"gamma_shape": 99.852225,
		"gamma_cats": 4
	}

	res_dir = "/home/elyawy/development/Msc/Thesis/Working_dir/SpartaABC/edo" #results path
	clean_run = False
	op_sys= 'linux'
	msa_filename = "temp_msa.fasta"
	tree_filename = "temp_tree.tree"
	minIR=0.0
	maxIR=0.05
	verbose = 0
	b_num_top=100
	num_alignments = 1
	filter_p = (0.9,15)
	num_simulations = 1
	num_burnin = 1
	




	pipeline(skip_config=skip_config,
			 pipeline_path="/home/elyawy/development/Msc/Thesis/Working_dir/SpartaABC/",
			 res_dir=res_dir, 
			 clean_run=clean_run,
			 msa_filename=msa_filename,
			 tree_filename=tree_filename,
			 minIR=minIR,
			 maxIR=maxIR,
			 op_sys=op_sys,
			 verbose=verbose,
			 b_num_top=b_num_top,
			 num_alignments=num_alignments,
			 filter_p=filter_p, 
			 submodel_params=submodel_params_,
			 num_simulations=num_simulations,
			 num_burnin=num_burnin
			 )

	if sample % num_of_samples_per_file == 0:
		file = f'all_data{sample}'
	with open(file, 'a') as f:
		f.write(tree + "\n")
		with open("indelible_sparta_eq.txt",'r') as f2:
			lines = f2.readlines()
			f.write(lines[1])
			f.write(lines[3] + "\n")