import re
import os

def scientific_to_decimal(expression):
	return f"{float(expression):.10f}"

def remove_bootstrap(expression):
	if expression[0] == ")":
		return "):"
	else:
		return ""

def fix_tree(file_path, file_name):
	scinot_regex = re.compile(r'[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)')
	bootstrap_regex = re.compile(r'\)(\d+(?:\.\d+)?):|:(\d+(?:\.\d+)?)\[(\d+(?:\.\d+)?)\]')
	with open(file_path+file_name, 'r') as f:
		tree = f.read()

	new_tree = scinot_regex.sub(lambda x: scientific_to_decimal(x.group()), tree)
	new_tree = bootstrap_regex.sub(lambda x: remove_bootstrap(x.group()), new_tree)

	with open(file_path+file_name, 'w') as f:
		f.write(new_tree)
