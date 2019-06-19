###Stavros Giannoukakos### 

#Version of the program
__version__ = "0.1.0"

import argparse
import subprocess
import shutil, time, glob, sys, os, re

raw_data_dir =  "/shared/projects/silvia_ont_umc/all_rna_sample2/20190604_1354_MN29521_FAK43904_3bcfbe09/fastq_pass"
refGRCh38 = "~/playground/progs/reference_files/reference_genome/GRCh38_primAssembly/GRCh38_primary_assembly_genome.fa"
refAnnot = "~/playground/progs/reference_files/gene_annotation/gencode.v29.primary_assembly.annotation.gtf"

usage = "ont_analysis [options]"
epilog = " -- June 2019 | Stavros Giannoukakos -- "
description = "DESCRIPTION"

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, usage=usage, description=description, epilog=epilog)
# Number of threads/CPUs to be used
parser.add_argument('-th', '--threads', dest='threads', default=20, metavar='', 
                	help="Number of threads to be used in the analysis")
# Display the version of the pipeline 
parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(__version__))
# Get the options and return them
args = parser.parse_args()

current_dir = os.getcwd()

# Main folder hosting the analysis
analysis_dir = os.path.join(current_dir, "analysis")
alignments_dir = os.path.join(analysis_dir, "alignments")
reports_dir = os.path.join(analysis_dir, "reports")
# Secondary subfolders
preprocessing_reports = os.path.join(reports_dir, "preprocessing_reports")
alignment_reports = os.path.join(reports_dir, "alignment_reports")

def preprocessing_samples():
	raw_data = [f for f in glob.glob(os.path.join(raw_data_dir, "*.fastq"))]

	return

def alignment_against_ref_genome():

	return

def main():
	
	quality_control()

	# alignment_against_ref_genome()

if __name__ == "__main__": main()