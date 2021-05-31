from collections import OrderedDict


def get_indelible_config():
	indelible_config = OrderedDict()

	indelible_config["[TYPE]"] = 'AMINOACID 2'
	indelible_config["[MODEL]"] = 'modelname'
	indelible_config["[submodel]"] = 'WAG'
	indelible_config["[indelmodel]"] = 'POW 1.7 500'
	indelible_config["[indelrate]"] = '0.0'
	indelible_config["[rates]"] = ' 0.25 0.50 10'
	indelible_config["[statefreq]"] = ' 0.25  0.25  0.25  0.25'
	indelible_config["[TREE]"] = 'treename (A:0.1,B:0.1);'
	indelible_config["[PARTITIONS]"] = 'partitionname\n[treename modelname 3000]'
	indelible_config["[EVOLVE]"] = "partitionname 100 outputname" + " " # note: indlible requires space at last command.

	return indelible_config

def get_sparta_config():
	sparta_config = OrderedDict()

	# sparta_config["_indelibleTemplateControlFile"] = "control_indelible_template.txt"
	# sparta_config["_dawgTemplateControlFile"] = ""
	# sparta_config["_dawgSimulator"] = "0"
	sparta_config["_inputRealMSAFile"] = "results/ref_msa.aa.fasta"
	sparta_config["_inputTreeFileName"] = "results/RAxML_tree.tree"
	sparta_config["_outputGoodParamsFile"] = "results/SpartaABC_data_name_model_name.posterior_params"
	sparta_config["_numberOfSamplesToKeep"] = "100"
	sparta_config["_alignmentMode"] = "0"
	sparta_config["_similarity_mode"] = "0"
	sparta_config["_modelType"] = "eq"
	sparta_config["_minRLVal"] = "50"
	sparta_config["_maxRLVal"] = "500.0"
	sparta_config["_minIRVal"] = "0"
	sparta_config["_maxIRVal"] = "0.05"
	sparta_config["_distanceCutOff"] = "1.0"

	sparta_config["_wAvgUniqueGapSize"] = "-1.0"
	sparta_config["_wMSAMin"] = "-1.0"
	sparta_config["_wNumGapsLenTwo"] = "-1.0"
	sparta_config["_wAvgGapSize"] = "-1.0"
	sparta_config["_wTotNumGaps"] = "-1.0"
	sparta_config["_wNumGapsLenAtLeastFour"] = "-1.0"
	sparta_config["_wNumGapsLenOne"] = "-1.0"
	sparta_config["_wMSAMax"] = "-1.0"
	sparta_config["_wMSALen"] = "-1.0"
	sparta_config["_wTotNumUniqueGaps"] = "-1.0"
	sparta_config["_wNumGapsLenThree"] = "-1.0"
	sparta_config["_wNumGapsLenOneIn1Pos"] = "-1.0"
	sparta_config["_wNumGapsLenOneIn2Pos"] = "-1.0"
	sparta_config["_wNumGapsLenOneInNMinus1Pos"] = "-1.0"
	sparta_config["_wNumGapsLenTwoIn1Pos"] = "-1.0"
	sparta_config["_wNumGapsLenTwoIn2Pos"] = "-1.0"
	sparta_config["_wNumGapsLenTwoInNMinus1Pos"] = "-1.0"
	sparta_config["_wNumGapsLenThreeIn1Pos"] = "-1.0"
	sparta_config["_wNumGapsLenThreeIn2Pos"] = "-1.0"
	sparta_config["_wNumGapsLenThreeInNMinus1Pos"] = "-1.0"
	sparta_config["_wNumGapsLenAtLeastFourIn1Pos"] = "-1.0"
	sparta_config["_wNumGapsLenAtLeastFourIn2Pos"] = "-1.0"
	sparta_config["_wNumGapsLenAtLeastFourInNMinus1Pos"] = "-1.0"
	sparta_config["_wNumberOfMSA_position_with_0_gaps"] = "-1.0"
	sparta_config["_wNumberOfMSA_position_with_1_gaps"] = "-1.0"
	sparta_config["_wNumberOfMSA_position_with_2_gaps"] = "-1.0"
	sparta_config["_wNumberOfMSA_position_with_n_minus_1_gaps"] = "-1.0"

	sparta_config["_only_real_stats"] = "0"
	sparta_config["_alignments_output"] = "10"
	sparta_config["_outputAlignmnetsFile"] = "results/alignments.fasta"

	sparta_config["_numSimulations"] = "100000"
	sparta_config["_numBurnIn"] = "10000"



	return sparta_config

class PipelineConfig:
	
	def __init__(self):
		# init in first block
		self.skip_config = {}
		...

pipeline_config = PipelineConfig()
