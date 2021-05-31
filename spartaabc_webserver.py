import sys
import os
import logging
import shutil
import Bio.SeqUtils
import pandas as pd
from time import sleep
import CONSTANTS as CONSTS  # from /bioseq/sincopa/

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')


def verify_fasta_format(fasta_path):
    logger.info('Validating FASTA format')
    Bio.SeqUtils.IUPACData.ambiguous_dna_letters += 'U-'
    legal_chars = set(Bio.SeqUtils.IUPACData.ambiguous_dna_letters.lower() +
                      Bio.SeqUtils.IUPACData.ambiguous_dna_letters +
                      Bio.SeqUtils.IUPACData.protein_letters.lower() +
                      Bio.SeqUtils.IUPACData.protein_letters)
    with open(fasta_path) as f:
        line_number = 0
        try:
            line = f.readline()
            line_number += 1
            if not line.startswith('>'):
                return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. First line in MSA starts with "{line[0]}" instead of ">".'
            previous_line_was_header = True
            putative_end_of_file = False
            curated_content = f'>{line[1:]}'.replace("|", "_")
            for line in f:
                line_number += 1
                line = line.strip()
                if not line:
                    if not putative_end_of_file: # ignore trailing empty lines
                        putative_end_of_file = line_number
                    continue
                if putative_end_of_file:  # non empty line after empty line
                    return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {putative_end_of_file} in MSA is empty.'
                if line.startswith('>'):
                    if previous_line_was_header:
                        return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. MSA contains an empty record. Both lines {line_number-1} and {line_number} start with ">".'
                    else:
                        previous_line_was_header = True
                        curated_content += f'>{line[1:]}\n'.replace("|", "_")
                        continue
                else:  # not a header
                    previous_line_was_header = False
                    for c in line:
                        if c not in legal_chars:
                            return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {line_number} in MSA contains illegal DNA character "{c}".'
                    curated_content += f'{line}\n'
        except UnicodeDecodeError as e:
            logger.info(e.args)
            line_number += 1  # the line that was failed to be read
            return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {line_number} in MSA contains one (or more) non <a href="https://en.wikipedia.org/wiki/ASCII" target="_blank">ascii</a> character(s).'
    # override the old file with the curated content
    with open(fasta_path, 'w') as f:
        f.write(curated_content)



