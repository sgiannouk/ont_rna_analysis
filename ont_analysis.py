###Stavros Giannoukakos### 

#Version of the program
__version__ = "0.1.5"

import argparse
import subprocess
import numpy as np
from pathlib import Path
from datetime import datetime
import shutil, time, glob, sys, os, re

ont_data =  "/home/stavros/playground/ont_basecalling/guppy_v3_basecalling"

refGenomeGRCh38 = "/home/stavros/references/reference_genome/GRCh38_GencodeV31_primAssembly/GRCh38.primary_assembly.genome.fa"
refTranscGRCh38 = "/home/stavros/references/reference_transcriptome/GRCh38_gencode.v31.transcripts_trna.fa"
refAnnot = "/home/stavros/references/reference_annotation/GRCh38_gencode.v31.primAssembly_pseudo_trna.annotation.gtf"
reference_annotation_bed = "/home/stavros/references/reference_annotation/hg38_gencode.v31.allComprehensive_pseudo.annotation.bed"

usage = "ont_analysis [options]"
epilog = " -- June 2019 | Stavros Giannoukakos -- "
description = "DESCRIPTION"

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, usage=usage, description=description, epilog=epilog)
# Number of threads/CPUs to be used
parser.add_argument('-t', '--threads', dest='threads', default=str(30), metavar='', 
                	help="Number of threads to be used in the analysis")
# Display the version of the pipeline 
parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(__version__))
# Get the options and return them
args = parser.parse_args()

current_dir = os.getcwd()
startTime = datetime.now()

# Main folder hosting the analysis
analysis_dir = os.path.join(current_dir, "batch3_analysis_TEST")
# analysis_dir = os.path.join(current_dir, "analysis_batch3")
prepr_dir = os.path.join(analysis_dir, "preprocessed_data")
alignments_dir = os.path.join(analysis_dir, "alignments")
reports_dir = os.path.join(analysis_dir, "reports")
postanalysis_dir = os.path.join(analysis_dir, "post-analysis")
methylation_dir = os.path.join(analysis_dir, "methylation_analysis")



def quality_control(i, seq_summary_file, sample_id, raw_data_dir):
	# Producing preliminary QC reports 
	if not os.path.exists(reports_dir): os.makedirs(reports_dir)
	print("\t{0} QUALITY CONTROL OF THE INPUT SAMPLES".format(datetime.now().strftime("%d.%m.%Y %H:%M")))

	if i == 1: print("{0}  pycoQC - Quality Control of the raw data: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	pycoQC = " ".join([
	"pycoQC",  # Call pycoQC
	"--summary_file", seq_summary_file,  # Input of <sequencing_summary> generated by Albacore1.0.0 / Guppy 2.1.3+
	"--html_outfile", os.path.join(reports_dir, "{0}_pycoQC-report.html".format(sample_id)),  # Create the report in the reports directory
	"--report_title", "\"{0} quality control report\"".format(sample_id),  # A title to be used in the html report
	"2>>", os.path.join(reports_dir, "pycoQC-prelim_report.txt")])
	subprocess.run(pycoQC, shell=True)

	if i == 1: print("{0}  nanoPlot - Quality Control of the raw data: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	nanoPlot = " ".join([
	"NanoPlot",  # Call pycoQC
	"--threads", args.threads,  # Number of threads to be used by the script
	"--summary", seq_summary_file,  # Input of <sequencing_summary> generated by Albacore1.0.0 / Guppy 2.1.3+
	"--color saddlebrown",  # Color for the plots
	"--colormap PuBuGn",  # Colormap for the heatmap
	"--prefix", os.path.join(reports_dir, "{0}_".format(sample_id)),  # Create the report in the reports directory
	"--format png",  # Output format of the plots
	"--dpi 900", # Set the dpi for saving images in high resolution
	"2>>", os.path.join(reports_dir, "nanoPlot-prelim_report.txt")])
	subprocess.run(nanoPlot, shell=True)
	# preprocessing_raw_data(raw_data_dir)
	return

def preprocessing_raw_data(raw_data_dir):
	if not os.path.exists(prepr_dir): os.makedirs(prepr_dir)
	raw_data = [sum_file for sum_file in glob.glob(os.path.join(raw_data_dir, "fastq_pass_file/*.fastq.gz"))][0]
	print("\t{0} PREPROCESSING THE INPUT SAMPLES".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	file_name = os.path.basename(file.split(".")[0]) 
	porechop = " ".join([
	"porechop",  # Call porechop to trim adapters
	"--threads", args.threads,  # Number of threads to be used by the script
	"--format fastq.gz",  # Output as compressed fastq files
	"--input", file,  # FASTQ of input reads
	"--output", os.path.join(prepr_dir, "{0}.tr.fastq.gz".format(file_name)),  # Filename for FASTQ of trimmed reads
	"--verbosity 3",  # Level of progress information
	"2>>", os.path.join(reports_dir, "porechop-report.txt")])
	subprocess.run(porechop, shell=True)
	return

def alignment_against_ref(i, sample_id, raw_data_dir):
	if not os.path.exists(alignments_dir): os.makedirs(alignments_dir)
	if i == 1: print("\n\t{0} ALIGNING AGAINSTE THE REF. GENOME AND TRANSCRIPTOME".format(datetime.now().strftime("%d.%m.%Y %H:%M")))

	dataset = [sum_file for sum_file in glob.glob(os.path.join(raw_data_dir, "pass/*.fastq.gz"))][0]

	### ALIGN THE RAW READS AGAINST THE REFERENCE GENOME
	print("{0}  minimap2 - Mapping {1} against the reference genome: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), sample_id))
	minimap2_genome = " ".join([
	"~/playground/progs/minimap2-2.17_x64-linux/minimap2",  # Call minimap2 (v2.17-r941)
	"-t", args.threads,  # Number of threds to use
	"-ax splice",   # Long-read spliced alignment mode and output in SAM format (-a)
	"-k 14",  # k-mer size
	"-uf",  # Find canonical splicing sites GT-AG - f: transcript strand
	"--secondary=no",  #
	# "-o", os.path.join(alignments_dir, "{0}.genome.paf".format(sample_id)),
	refGenomeGRCh38,  # Inputting the reference genome
	dataset,  # Input .fastq.gz file
	"|", "samtools sort",  # Calling 'samtools sort' to sort the output alignment file
	"--threads", args.threads,  # Number of threads to be used by 'samtools sort'
	"--output-fmt BAM",  # Output in BAM format
	"-o", os.path.join(alignments_dir, "{0}.genome.bam".format(sample_id)),  # Sorted output  BAM file
	"2>>", os.path.join(reports_dir, "minimap2_genome-report.txt")])  # Directory where all FastQC and Cutadapt reports reside
	subprocess.run(minimap2_genome, shell=True)

	### EXTRACTING THE UNMAPPED READS (GENOME)
	subprocess.run("samtools view -f4 --threads {2} {0}/{1}.genome.bam | samtools bam2fq --threads {2} - | gzip > {0}/{1}.genome.unmapped.fastq.gz".format(alignments_dir, sample_id, args.threads), shell=True)

	## ALIGN THE UNMAPPED READS AGAINST THE REFERENCE TRANSCRIPTOME
	# print(">> minimap2 - Mapping the analigend reads against the reference transcriptome: in progress ..")
	print("\n{0}  minimap2 - Mapping {1} against the reference transcriptome: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), sample_id))	
	minimap2_transcriptome = " ".join([
	"~/playground/progs/minimap2-2.17_x64-linux/minimap2",  # Call minimap2 (v2.17-r941)
	"-t", args.threads,  # Number of threds to use
	"-ax map-ont",
	"-k 14",  # k-mer size
	"--secondary=no",  # 
	# "-o", os.path.join(alignments_dir, "{0}.transcriptome.paf".format(sample_id)),
	refTranscGRCh38,  # Inputting the reference genome
	dataset,
	# "{0}/{1}.genome_unmapped.bam".format( alignments_dir, sample_id),  # Input .fastq.gz file
	"|", "samtools sort",  # Calling 'samtools sort' to sort the output alignment file
	"--threads", args.threads,  # Number of threads to be used by 'samtools sort'
	"--output-fmt BAM",  # Output in BAM format
	"-o", os.path.join(alignments_dir, "{0}.transcriptome.bam".format(sample_id)), "-",  # Sorted output  BAM file
	"2>>", os.path.join(reports_dir, "minimap2_transcriptome-report.txt")])  # Directory where all FastQC and Cutadapt reports reside
	subprocess.run(minimap2_transcriptome, shell=True)

	### EXTRACTING THE UNMAPPED READS (TRANSCRIPTOME)
	subprocess.run("samtools view -f4 --threads {2} {0}/{1}.transcriptome.bam | samtools bam2fq --threads {2} - | gzip > {0}/{1}.transcriptome.unmapped.fastq.gz".format(alignments_dir, sample_id, args.threads), shell=True)
	return

def mapping_qc():
	print("\n\t{0} GENERATING ALIGNMENT STATS".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	if not os.path.exists(postanalysis_dir): os.makedirs(postanalysis_dir)
	summary_files = [file_path for file_path in Path(ont_data).glob('**/sequencing_summary.txt')]

	aligned_data = [sum_file for sum_file in glob.glob(os.path.join(alignments_dir, "*ome.bam"))]
	genome_alignments = [sum_file for sum_file in glob.glob(os.path.join(alignments_dir, "*.genome.bam"))]

	# print("{0}  RSeQC, Picard, AlignQC & WUB - Generating post-alignment stats for {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), file_name))
	for i, file in enumerate(aligned_data, 1):
		file_name = os.path.basename(file).split(".")[0]
		aligned_to = ["transcriptome", "genome"][file.endswith(".genome.bam")]
		ref_file = [refTranscGRCh38, refGenomeGRCh38][file.endswith(".genome.bam")]
		ref_summary_file = [str(files) for files in summary_files if file_name in str(files)]

		samtools_index = " ".join([
		"samtools index",  # Indexing the concat_samples.bam file
		"-@", args.threads,  # Number of threads to be used
		file])  # Input BAM file
		subprocess.run(samtools_index, shell=True)

		### EXPORTING ALIGNMENT STATS
		print("{0}  PycoQC - Generating post-alignment stats of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), file_name))
		pycoQC = " ".join([
		"pycoQC",  # Call pycoQC
		"--summary_file", ref_summary_file[0], 
		"--bam_file", file,  # Input of bam files from Minimap2
		"--html_outfile", os.path.join(postanalysis_dir, "{0}.{1}.pycoQC-report.html".format(file_name, aligned_to)),  # Create the report in the reports directory
		"--report_title", "\"Post-alignment quality control report\"",  # A title to be used in the html report
		"2>>", os.path.join(postanalysis_dir, "pycoQC-report.txt")])
		subprocess.run(pycoQC, shell=True)
	
		# BAM stats
		print("{0}  BAM stats - Generating post-alignment stats of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), file_name))
		bam_stat = ' '.join([
		"bam_stat.py",
		"-i", file,  # Input BAM file
		"> {0}/{1}.{2}.bamstat.txt".format(postanalysis_dir, file_name, aligned_to),  # Output file
		"2>>", os.path.join(postanalysis_dir, "{0}_{1}_bamstats-report.txt".format(file_name, aligned_to))])
		subprocess.run(bam_stat, shell=True)

		# Picard CollectAlignmentSummaryMetrics
		print("{0}  Picard - Collecting alignment summary metrics of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), file_name))
		CollectAlignmentSummaryMetrics = ' '.join([
		"picard CollectAlignmentSummaryMetrics",  # Call picard CollectAlignmentSummaryMetrics
		"INPUT= {0}".format(file),  # Input BAM file
		"OUTPUT= {0}/{1}.{2}.alignment_metrics.txt".format(postanalysis_dir, file_name, aligned_to),  # Output
		"REFERENCE_SEQUENCE= {0}".format(ref_file),  # Reference sequence file
		"2>>", os.path.join(postanalysis_dir, "{0}_{1}_CollectAlignmentSummaryMetrics-report.txt".format(file_name, aligned_to))])
		subprocess.run(CollectAlignmentSummaryMetrics, shell=True) 

		# AlignQC 
		print("{0}  AlignQC - Generating post-alignment stats of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), file_name))
		alignqc = ' '.join([
		"alignqc analyze",  # Calling AlignQC analyze
		"--threads", str(int(int(args.threads)/2)),
		"--genome", ref_file,  # Reference in .fasta
		"--gtf", ''.join([refAnnot, ".gz"]),  # Input reference annotation in gzipped form
		"--output", "{0}/{1}.{2}.alignqc.xhtml".format(postanalysis_dir, file_name, aligned_to),  # Output pdf file
		file,
		"2>>", os.path.join(postanalysis_dir, "{0}_{1}_alignqc-report.txt".format(file_name, aligned_to))])
		subprocess.run(alignqc, shell=True)

		# Wub Alignment based QC plots
		# print("{0}  WUB - Alignment based QC plots of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), file_name))
		alignment_qc = ' '.join([
		"bam_alignment_qc.py",
		"-f", ref_file,  # Input reference file
		"-x -Q",  # Do not plot per-reference information/ Keep qiet
		"-r", "{0}/{1}.{2}.bam_alignment_qc.pdf".format(postanalysis_dir, file_name, aligned_to),  # Output pdf file
		"-p", "{0}/{1}.{2}.bam_alignment_qc.pk".format(postanalysis_dir, file_name, aligned_to),  # Output pk file
		file,
		"2>>", os.path.join(postanalysis_dir, "{0}_{1}_alignment_qc-report.txt".format(file_name, aligned_to))])
		# subprocess.run(alignment_qc, shell=True)

		### GENOME ALIGNMENTS STATS 
		if file.endswith(".genome.bam"):
			# BAM read distribution
			print("{0}  RSeQC - Generating read distribution stats of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), file_name))
			read_distribution = ' '.join([
			"read_distribution.py", 
			"-i", file,  # Input BAM file
			"-r", reference_annotation_bed,
			"> {0}/{1}.{2}.fragSize".format(postanalysis_dir, file_name, aligned_to),  # Output file
			"2>>", os.path.join(postanalysis_dir, "{0}_{1}_read_distribution-report.txt".format(file_name, aligned_to))])
			subprocess.run(read_distribution, shell=True)

			# Check the strandness of the reads
			print("{0}  RSeQC - Generating read strandness stats of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), file_name))
			strandness = ' '.join([
			"infer_experiment.py",  # Call samtools flagstat
			"-i", file,  # Input BAM file
			"-r", reference_annotation_bed,
			"> {0}/{1}.{2}.strandness.txt".format(postanalysis_dir, file_name, aligned_to),  # Output file
			"2>>", os.path.join(postanalysis_dir, "{0}_{1}_strandness-report.txt".format(file_name, aligned_to))])
			subprocess.run(strandness, shell=True)

			# Gene body coverage
			print("{0}  RSeQC - Generating gene body coverage stats of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), file_name))
			gene_coverage = ' '.join([
			"geneBody_coverage.py",  # Call samtools flagstat
			"-r", reference_annotation_bed,
			"-i", file,  # Input BAM file
			"-o", "{0}/{1}.{2}.gene_coverage".format(postanalysis_dir, file_name, aligned_to),  # Output file
			"2>>", os.path.join(postanalysis_dir, "{0}_{1}_mapping_position-report.txt".format(file_name, aligned_to))])
			subprocess.run(gene_coverage, shell=True)

			# Check duplicate reads
			print("{0}  Picard - Extracting read duplication stats of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), file_name))
			duplicate_reads = ' '.join([
			"picard MarkDuplicates",  # Call samtools flagstat
			"INPUT= {0}".format(file),  # Input BAM file
			"OUTPUT= {0}/{1}.genome.dedup.bam".format(alignments_dir, file_name),
			"METRICS_FILE= {0}/{1}.{2}.mark_duplicates.txt".format(postanalysis_dir, file_name, aligned_to),  # Output file
			"2>>", os.path.join(postanalysis_dir, "{0}_{1}_duplicate_reads-report.txt".format(file_name, aligned_to))])
			subprocess.run(duplicate_reads, shell=True)
			os.system('rm {0}/{1}.genome.dedup.bam'.format(alignments_dir, file_name))

			# Number of reads mapped to each chromosome
			print("{0}  RSeQC - Generating mapping stats of {1}: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M"), file_name))
			mapping_pos = ' '.join([
			"samtools idxstats",  # Call samtools flagstat
			"--threads", args.threads,
			file,  # Input BAM file
			"> {0}/{1}.{2}.samtools_idxstats.txt".format(postanalysis_dir, file_name, aligned_to),  # Output file
			"2>>", os.path.join(postanalysis_dir, "{0}_{1}_mapping_position-report.txt".format(file_name, aligned_to))])
			subprocess.run(mapping_pos, shell=True)
	return

class expression_matrix:

	def __init__(self, threads):
		if not os.path.exists(postanalysis_dir): os.makedirs(postanalysis_dir)
		
		# self.generate_perGene_expression_matrix(threads)
		# self.generate_perTranscript_expression_matrix(threads)
		# self.generate_expression_matrices_Salmon()
		# self.novel_transcripts_detection(args.threads)
		# self.clean_expression_dir()
		return

	def generate_perGene_expression_matrix(self, threads):

		genome_alignments = [sum_file for sum_file in glob.glob(os.path.join(alignments_dir, "*.genome.bam"))]
		print("\n\t{0} GENERATING THE PER-GENE EXPRESSION MATRIX".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		print("{0}  FeatureCounts perGene - Counting reads from the genome aligned data: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		featureCounts_gn = " ".join([
		"featureCounts",  # Call featureCounts
		"-T", threads,  # Number of threads to be used by the script
		"-a", refAnnot,  # Annotation file in GTF/GFF format
		"-L",  # Count long reads such as Nanopore and PacBio reads
		"-o", os.path.join(postanalysis_dir, "genome_alignments_perGene_sum.tab"),
		' '.join(genome_alignments),  # Input bam file
		"2>>", os.path.join(postanalysis_dir, "featureCounts_genome_summary_perGene-report.txt")]) 
		subprocess.run(featureCounts_gn, shell=True)
		
		print("{0}  Exporting the per gene (genomic) expression matrix: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		subprocess.run("cut -f1,7- {0}/genome_alignments_perGene_sum.tab | sed 1d > {0}/featureCounts_expression_perGene_matrix.txt".format(postanalysis_dir), shell=True)

		annot = {}
		with open(refAnnot) as ref_in:
			for i, line in enumerate(ref_in):
				if not line.startswith("#"):
					if line.split("\t")[2].strip() == 'transcript' or line.split("\t")[2].strip().endswith('tRNA'):
						gene_id = line.split("\t")[-1].split(";")[0].split(" ")[-1].strip("\"")
						gene_type = [line.split("\t")[-1].split(";")[1].split()[1].strip("\"").strip() ,line.split("\t")[-1].split(";")[2].split()[1].strip("\"").strip()]\
									[line.split("\t")[-1].split(";")[2].split()[0].strip()=="gene_type"]
						annot[gene_id] = gene_type

		data = {}
		header = []
		# nc = ["miRNA", "piRNA", "rRNA", "siRNA", "snRNA", "snoRNA", "tRNA", "vaultRNA"]
		with open("{0}/featureCounts_expression_perGene_matrix.txt".format(postanalysis_dir)) as mat_in:
			for i, line in enumerate(mat_in, 1):
				if i == 1:
					header = line.split()
				else:
					gene = line.strip().split()[0]
					values = line.strip().split()[1:]
					if not sum(map(int, values)) == 0:
						# Grouping all pseudogenes
						if annot[gene].endswith("pseudogene"):
							data[(gene, "pseudogene")] = values
						# Grouping all immunoglobin genes
						elif annot[gene].startswith("IG_"):
							data[(gene, "Immunoglobulin_gene")] = values
						# Grouping all T-cell receptor genes
						elif annot[gene].startswith("TR_"):
							data[(gene, "Tcell_receptor_gene")] = values
						# Reanming TEC group
						elif annot[gene].startswith("TEC"):
							data[(gene, "To_be_Experimentally_Confirmed")] = values
						# Rest
						else:
							data[(gene, annot[gene])] = values

		# Remove paths from sample names
		header = [elm.split("/")[-1].split(".")[0] for elm in header]
		header.insert(1, "gene_type")  # Inserting gene_type in header
		header[0] = "gene_id"  # Replacing Geneid with gene_id in header

		# Writing output to file 'expression_matrix.csv'
		with open("{0}/perGene_expression_matrix.csv".format(postanalysis_dir), "a") as fout:
			fout.write("{0}\n".format(','.join(header)))
			for key, values in data.items():
				fout.write("{0},{1}\n".format(','.join(key), ','.join(values)))

		print("{0}  Visualising the RNA categories found in the (genomic) expression matrix: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		gene_type_sum = " ".join([
		"Rscript",  # Call Rscript
		"gene_type_summary.R",  # Calling the script
		"{0}/perGene_expression_matrix.csv".format(postanalysis_dir),  # Input matrix
		postanalysis_dir,  # Output dir
		"2>>", os.path.join(postanalysis_dir, "R_gene_type_sum-report.txt")]) 
		subprocess.run(gene_type_sum, shell=True)
		return

	def generate_perTranscript_expression_matrix(self, threads):
		
		genome_alignments = [sum_file for sum_file in glob.glob(os.path.join(alignments_dir, "*.genome.bam"))]
		print("\n\t{0} GENERATING THE PER-TRANSCRIPT EXPRESSION MATRIX".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		print("{0}  FeatureCounts per Transcript - Counting reads from the genome aligned data: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		featureCounts_gntr = " ".join([
		"featureCounts",  # Call featureCounts
		"-T", threads,  # Number of threads to be used by the script
		"-g", "\'transcript_id\'",  # 
		"-a", refAnnot,  # Annotation file in GTF/GFF format
		"-L",  # Count long reads such as Nanopore and PacBio reads
		"-o", os.path.join(postanalysis_dir, "genome_alignments_perTranscript_sum.tab"),
		' '.join(genome_alignments),  # Input bam file
		"2>>", os.path.join(postanalysis_dir, "featureCounts_genome_summary_perTranscript-report.txt")]) 
		subprocess.run(featureCounts_gntr, shell=True)
		
		print("{0}  Exporting the per transcript (genomic) expression matrix: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		subprocess.run("cut -f1,7- {0}/genome_alignments_perTranscript_sum.tab | sed 1d > {0}/featureCounts_expression_perTranscript_matrix.txt".format(postanalysis_dir), shell=True)

		annot = {}
		with open(refAnnot) as ref_in:
			for i, line in enumerate(ref_in):
				if not line.startswith("#"):
					if line.split("\t")[2].strip() == 'transcript' or line.split("\t")[2].strip().endswith('tRNA'):
						transcript_id = line.split("\t")[-1].split(";")[1].split(" ")[-1].strip("\"")
						gene_id = line.split("\t")[-1].split(";")[0].split(" ")[-1].strip("\"")
						gene_name = line.split("\t")[-1].split(";")[3].split(" ")[-1].strip("\"")
						gene_type = [line.split("\t")[-1].split(";")[1].split()[1].strip("\"").strip() ,line.split("\t")[-1].split(";")[2].split()[1].strip("\"").strip()]\
									[line.split("\t")[-1].split(";")[2].split()[0].strip()=="gene_type"]
						annot[transcript_id] = [gene_id, gene_name, gene_type]

		data = {}
		header = []
		with open("{0}/featureCounts_expression_perTranscript_matrix.txt".format(postanalysis_dir)) as mat_in:
			for i, line in enumerate(mat_in, 1):
				if i == 1:
					header = line.split()
				else:
					transcript = line.strip().split()[0]
					values = line.strip().split()[1:]
					if not sum(map(int, values)) == 0:
						data[(transcript, ' '.join(annot[transcript]))] = values
		
		# Remove paths from sample names
		header = [elm.split("/")[-1].split(".")[0] for elm in header]
		header.insert(1, "gene_id")  # Inserting gene_id in header
		header.insert(2, "gene_name")  # Inserting gene_name in header
		header.insert(3, "gene_type")  # Inserting gene_type in header
		header[0] = "transcript_id"  # Replacing Geneid with transcript_id in header

		# Writing output to file 'expression_matrix.csv'
		with open("{0}/perTranscript_expression_matrix.csv".format(postanalysis_dir), "a") as fout:
			fout.write("{0}\n".format(','.join(header)))
			for key, values in data.items():
				fout.write("{0},{1}\n".format(','.join(key), ','.join(values)))
		return

	def generate_expression_matrices_Salmon(self):
		
		genome_alignments = [sum_file for sum_file in glob.glob(os.path.join(alignments_dir, "*.transcriptome.bam"))]
		print("{0}  Salmon quant (alignment-based) - Counting reads from the transcriptome aligned data: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
		salmon_quant = " ".join([
		"salmon quant",  # Call salmon quant
		"--threads 6",  # 6 cores to be used
		"--libType", "A",  # Format string describing the library type
		"--targets", refTranscGRCh38,  #
		"--geneMap", refAnnot,  # Annotation file in GTF format
		"--alignments", ' '.join(genome_alignments),  # Input bam file
		"--output", os.path.join(postanalysis_dir, "salmon_analysis_tr"),
		"2>>", os.path.join(postanalysis_dir, "salmonQuant_transcriptome_summary-report.txt")]) 
		subprocess.run(salmon_quant, shell=True)
		return	

	def novel_transcripts_detection(self):
		""" TALON takes transcripts from one or more long read datasets (SAM format) 
		and assigns them transcript and gene identifiers based on a database-bound
		annotation. Novel events are assigned new identifiers """

		return

	def clean_expression_dir(self):
		os.system('mv {0}/*report.txt {1}'.format(postanalysis_dir, reports_dir))
		os.system('mv {0}/*tab.summary {1}'.format(postanalysis_dir, reports_dir))
		os.system('mv {0}/*sum.tab {1}'.format(postanalysis_dir, reports_dir))
		return

class special_analysis:

	def __init__(self):
		self.structural_variation()
		self.methylation_detection()
		return

	def methylation_detection(self):
		if not os.path.exists(methylation_dir): os.makedirs(methylation_dir)
		single_fast5_data = [os.path.dirname(sum_file) for sum_file in glob.glob("/shared/projects/silvia_rna_ont_umc/basecalling/sample*_raw_data/*/")]

		# print(single_fast5_data)
		for dirs in single_fast5_data:
			if dirs.endswith("/0"):

				sample_id = dirs.split("/")[5].split("_")[0]
				# 1. Re-squiggling the raw reads
				tombo_resquiggle = " ".join([
				"tombo resquiggle",
				dirs,
				refGenomeGRCh38,
				"--processes", args.threads,  # Number of threads to be used
				"--num-most-common-errors 5",
				# "2>>", os.path.join(methylation_dir, "tombo_resquiggle-report.txt")
				])
				subprocess.run(tombo_resquiggle, shell=True)
			
			# 2. Calling tombo to do the methylation analysis
			tombo_methyl = " ".join([
			"tombo detect_modifications de_novo",  # Indexing the concat_samples.bam file
			"--fast5-basedirs", dirs,  # Directory containing fast5 files
			"--statistics-file-basename", os.path.join(methylation_dir,"{0}.tombo.stats".format(sample_id)),
			"--rna",  # Explicitly select canonical RNA mode
			"--processes", args.threads,  # Number of threads to be used
			# "2>>", os.path.join(methylation_dir, "tombo_methylation-report.txt")
			])
			subprocess.run(tombo_methyl, shell=True)

			# 3. Output reference sequence around most significantly modified sites
			tombo_sign = " ".join([
			"tombo text_output signif_sequence_context",  # Indexing the concat_samples.bam file
			"--fast5-basedirs", dirs,  # Directory containing fast5 files
			"--statistics-filename", os.path.join(methylation_dir,"{0}.tombo.stats".format(sample_id)),
			"--sequences-filename",  os.path.join(methylation_dir,"{0}.tombo_significant_regions.fasta".format(sample_id)),
			# "2>>", os.path.join(methylation_dir, "tombo_significants-report.txt")
			])
			subprocess.run(tombo_sign, shell=True)

			# # 4. Use line Meme to estimate modified motifs
			# tombo_meme = " ".join([
			# "meme",
	  #  		"-rna",
	  #  		"-mod zoops", 
			# "-oc", os.path.join(methylation_dir,"{0}_de_novo_meme".format(sample_id)),
	  #  		os.path.join(methylation_dir,"{0}.tombo_significant_regions.fasta".format(sample_id)),
	  #  		# "2>>", os.path.join(methylation_dir, "tombo_meme-report.txt")
			# ])
			# subprocess.run(tombo_meme, shell=True)

			# # 5. This plot will identify the sites in the reference
			# tombo_plot = " ".join([
			# "tombo plot motif_with_stats",
			# " --fast5-basedirs", dirs,
	  #  		"--motif CCWGG",
	  #  		"--genome-fasta", refGenomeGRCh38,
	  #  		"--statistics-filename", os.path.join(methylation_dir,"{0}.tombo.stats".format(sample_id)),
	  #  		"--pdf-filename", os.path.join(methylation_dir,"{0}.tombo.pdf".format(sample_id))
	  #  		# "2>>", os.path.join(methylation_dir, "tombo_plot-report.txt")
			# ])
			# subprocess.run(tombo_plot, shell=True)
		return

	def structural_variation(self):

		return

def summary(num_of_files):
	print("{0}  multiQC - Summing all QC reports: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	multiQC = " ".join([
	"multiqc",  # Call MultiQC
	"--quiet",  # Print only log warnings
	"--outdir", postanalysis_dir,  # Create report in the FastQC reports directory
	"--filename", "post-alignment_summarised_report",  # Name of the output report 
	postanalysis_dir,  # Directory where all FastQC and Cutadapt reports reside
	"2>>", os.path.join(postanalysis_dir, "post-alignment_multiQC-report.txt")])  # Output multiQC report
	subprocess.run(multiQC, shell=True)

	pickle_files = [pk_file for pk_file in glob.glob(os.path.join(postanalysis_dir, "*_qc.pk"))]
	# Wub Compare alignment QC statistics of multiple samples
	bam_multi_qc = ' '.join([
	"bam_multi_qc.py",
	"-r", "{0}/comparison_qc.pdf".format(postanalysis_dir),  # Output pdf file
	' '.join(pickle_files),
	"2>>", os.path.join(postanalysis_dir, "bam_multi_qc-report.txt")])
	# subprocess.run(bam_multi_qc, shell=True)

	## Cleaning up reports folder
	figures_pre_dir = os.path.join(reports_dir, "figures")
	if not os.path.exists(figures_pre_dir): os.makedirs(figures_pre_dir)

	os.system('mv {0}/*.png {1}'.format(reports_dir, figures_pre_dir))  # Moving png files to "figures" directory
	## REMOVING UNNECESSARY FILES & REPORTS (reports)
	for path, subdir, folder in os.walk(reports_dir):
		for name in folder:
			file = os.path.join(path, name)
			if os.stat(file).st_size == 0 or\
			(name.startswith("minimap2_") and os.stat(file).st_size == (int(num_of_files) * 64)) or\
			(name.startswith("pycoQC-") and os.stat(file).st_size == 1481):
				os.remove(file)
	
	## REMOVING UNNECESSARY FILES & REPORTS (postanalysis)
	for path, subdir, folder in os.walk(postanalysis_dir):
		for name in folder:
			file = os.path.join(path, name)
			if os.stat(file).st_size == 0 or\
			(name.endswith("multiQC-report.txt") and os.stat(file).st_size == 583) or\
			(name.endswith("bamstats-report.txt") and os.stat(file).st_size == 24) or\
	  		(name.endswith("ribution-report.txt") and os.stat(file).st_size == (255 or 256)) or\
	  		(name.endswith("Metrics-report.txt") and os.stat(file).st_size == 3) or\
	  		 name.endswith(".r") or name.endswith("DupRate.xls") or name.endswith("duplicate_reads-report.txt"):
				os.remove(file)

	## Cleaning up postanalysis folder
	qc_reports = os.path.join(postanalysis_dir, "qc_reports")
	if not os.path.exists(qc_reports): os.makedirs(qc_reports)

	# Moving files to "qc_reports" directory
	# os.system('mv {0}/*_qc.pk {1}'.format(postanalysis_dir, qc_reports))
	os.system('mv {0}/*log* {1}'.format(current_dir, qc_reports))
	os.system('mv {0}/*.fragSize {1}'.format(postanalysis_dir, qc_reports))
	os.system('mv {0}/*sum.* {1}'.format(postanalysis_dir, qc_reports))
	os.system('mv {0}/*curves.pdf {1}'.format(postanalysis_dir, qc_reports))
	os.system('mv {0}/*genome*.txt {1}'.format(postanalysis_dir, qc_reports))
	os.system('mv {0}/*transcriptome*.txt {1}'.format(postanalysis_dir, qc_reports))
	os.system('mv {0}/*summarised_report_data {1}'.format(postanalysis_dir, qc_reports))
	return


def main():
	
	# chosen_samples = ("NonTransf_1",  "NonTransf_2",  "NonTransf_3", "Tumour_1",  "Tumour_2",  "Tumour_3")
	chosen_samples = ("Tumour_1", "Tumour_2")

	summary_files = [str(file_path) for file_path in Path(ont_data).glob('**/sequencing_summary.txt')]
	num_of_samples = len(summary_files)

	for i, sum_file in enumerate([s for s in summary_files if os.path.dirname(s).endswith(chosen_samples)], 1):
		raw_data_dir = os.path.dirname(str(sum_file))
		sample_id = os.path.basename(raw_data_dir)

		# quality_control(i, sum_file, sample_id, raw_data_dir)

		# alignment_against_ref(i, sample_id, raw_data_dir)

	# mapping_qc()

	expression_matrix(args.threads)

	# summary(num_of_samples)

	# 

	# 

	

	print('\t--- The pipeline finisded after {0} ---'.format(datetime.now() - startTime))
if __name__ == "__main__": main()