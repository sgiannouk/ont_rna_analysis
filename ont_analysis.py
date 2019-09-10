###Stavros Giannoukakos### 

#Version of the program
__version__ = "0.1.4"

import argparse
import subprocess
from pathlib import Path
import shutil, time, glob, sys, os, re

ont_data =  "/shared/projects/silvia_ont_umc/"

refGenomeGRCh38 = "/home/stavros/references/reference_genome/GRCh38_GencodeV31_primAssembly/GRCh38.primary_assembly.genome.fa"
refTranscGRCh38 = "/home/stavros/references/reference_transcriptome/ensembl_cdna_ncrna/GRCh38_cdna_ncrna.fasta"
refAnnot = "/home/stavros/references/reference_annotation/GRCh38_gencode.v31.primAssembly_psudo_trna.annotation.gtf"
reference_annotation_bed = "/home/stavros/references/reference_annotation/hg38_gencode.v31.allComprehensive_pseudo.annotation.bed"

usage = "ont_analysis [options]"
epilog = " -- June 2019 | Stavros Giannoukakos -- "
description = "DESCRIPTION"

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, usage=usage, description=description, epilog=epilog)
# Number of threads/CPUs to be used
parser.add_argument('-th', '--threads', dest='threads', default=str(30), metavar='', 
                	help="Number of threads to be used in the analysis")
# Display the version of the pipeline 
parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(__version__))
# Get the options and return them
args = parser.parse_args()

current_dir = os.getcwd()

# Main folder hosting the analysis
analysis_dir = os.path.join(current_dir, "ont_data_analysis")
prepr_dir = os.path.join(analysis_dir, "preprocessed_data")
alignments_dir = os.path.join(analysis_dir, "alignments")
reports_dir = os.path.join(analysis_dir, "pre-analysis")
postanalysis_dir = os.path.join(analysis_dir, "post-analysis")



def quality_control(i, seq_summary_file, sample_id, raw_data_dir):
	if not os.path.exists(reports_dir): os.makedirs(reports_dir)

	if i == 1: print(">>> pycoQC - Quality Control of the raw data: in progress ..")
	pycoQC = " ".join([
	"pycoQC",  # Call pycoQC
	"--summary_file", seq_summary_file,  # Input of <sequencing_summary> generated by Albacore1.0.0 / Guppy 2.1.3+
	"--html_outfile", os.path.join(reports_dir, "{0}_pycoQC-report.html".format(sample_id)),  # Create the report in the reports directory
	"--report_title", "\"{0} quality control report\"".format(sample_id),  # A title to be used in the html report
	"2>>", os.path.join(reports_dir, "pycoQC-prelim_report.txt".format(sample_id))])
	subprocess.run(pycoQC, shell=True)

	if i == 1: print(">>> nanoPlot - Quality Control of the raw data: in progress ..")
	nanoPlot = " ".join([
	"NanoPlot",  # Call pycoQC
	"--threads", args.threads,  # Number of threads to be used by the script
	"--summary", seq_summary_file,  # Input of <sequencing_summary> generated by Albacore1.0.0 / Guppy 2.1.3+
	"--color saddlebrown",  # Color for the plots
	"--colormap PuBuGn",  # Colormap for the heatmap
	"--prefix", os.path.join(reports_dir, "{0}_".format(sample_id)),  # Create the report in the reports directory
	"--format png",  # Output format of the plots
	"--dpi 900", # Set the dpi for saving images in high resolution
	"2>>", os.path.join(reports_dir, "nanoPlot-prelim_report.txt".format(sample_id))])
	subprocess.run(nanoPlot, shell=True)
	# preprocessing_raw_data(raw_data_dir)
	return

def preprocessing_raw_data(raw_data_dir):
	if not os.path.exists(prepr_dir): os.makedirs(prepr_dir)
	raw_data = [sum_file for sum_file in glob.glob(os.path.join(raw_data_dir, "fastq_pass_file/*.fastq.gz"))][0]

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

	dataset = [sum_file for sum_file in glob.glob(os.path.join(raw_data_dir, "fastq_pass_file/*fastq.gz"))][0]
	file_name = os.path.basename(dataset.split(".")[0])

	### ALIGN THE RAW READS AGAINST THE REFERENCE GENOME
	print(">>> minimap2 - Mapping {0} against the reference genome: in progress ..".format(file_name))
	minimap2_genome = " ".join([
	"~/playground/progs/minimap2-2.17_x64-linux/minimap2",  # Call minimap2 (v2.17-r941)
	"-t", args.threads,  # Number of threds to use
	"-ax splice",   # Long-read spliced alignment mode and output in SAM format (-a)
	"-k 14",  # k-mer size
	"-uf",  # Find canonical splicing sites GT-AG - f: transcript strand
	"--secondary=no",  #
	# "-o", os.path.join(alignments_dir, "{0}.genome.paf".format(file_name)),
	refGenomeGRCh38,  # Inputting the reference genome
	dataset,  # Input .fastq.gz file
	"|", "samtools sort",  # Calling 'samtools sort' to sort the output alignment file
	"--threads", args.threads,  # Number of threads to be used by 'samtools sort'
	"--output-fmt BAM",  # Output in BAM format
	"-o", os.path.join(alignments_dir, "{0}.{1}.genome.bam".format(sample_id, file_name.split("_")[0])),  # Sorted output  BAM file
	"2>>", os.path.join(reports_dir, "minimap2_genome-report.txt")])  # Directory where all FastQC and Cutadapt reports reside
	subprocess.run(minimap2_genome, shell=True)

	# ### EXTRACTING THE UNMAPPED READS
	# subprocess.run("samtools view -f4 {0}/{1}.genome.bam > {0}/{1}.genome_unmapped.bam".format(alignments_dir, file_name), shell=True)

	## ALIGN THE UNMAPPED READS AGAINST THE REFERENCE TRANSCRIPTOME
	# print(">> minimap2 - Mapping the analigend reads against the reference transcriptome: in progress ..")
	print(">>> minimap2 - Mapping {0} against the reference transcriptome: in progress ..".format(file_name))	
	minimap2_transcriptome = " ".join([
	"~/playground/progs/minimap2-2.17_x64-linux/minimap2",  # Call minimap2 (v2.17-r941)
	"-t", args.threads,  # Number of threds to use
	"-ax map-ont",
	"-k 14",  # k-mer size
	"--secondary=no",  # 
	# "-o", os.path.join(alignments_dir, "{0}.transcriptome.paf".format(file_name)),
	refTranscGRCh38,  # Inputting the reference genome
	dataset,
	# "{0}/{1}.genome_unmapped.bam".format( alignments_dir, file_name),  # Input .fastq.gz file
	"|", "samtools sort",  # Calling 'samtools sort' to sort the output alignment file
	"--threads", args.threads,  # Number of threads to be used by 'samtools sort'
	"--output-fmt BAM",  # Output in BAM format
	"-o", os.path.join(alignments_dir, "{0}.{1}.transcriptome.bam".format(sample_id, file_name.split("_")[0])), "-",  # Sorted output  BAM file
	"2>>", os.path.join(reports_dir, "minimap2_transcriptome-report.txt")])  # Directory where all FastQC and Cutadapt reports reside
	subprocess.run(minimap2_transcriptome, shell=True)
	return

def mapping_qc():

	if not os.path.exists(postanalysis_dir): os.makedirs(postanalysis_dir)
	aligned_data = [sum_file for sum_file in glob.glob(os.path.join(alignments_dir, "*ome.bam"))]
	genome_alignments = [sum_file for sum_file in glob.glob(os.path.join(alignments_dir, "*.genome.bam"))]

	### EXPORTING ALIGNMENT STATS
	print(">>> RSeQC, Picard & WUB - Generating post-alignment stats: in progress ..")
	for i, file in enumerate(aligned_data, 1):
		file_name = os.path.basename(file).split(".")[0]
		aligned_to = ["transcriptome", "genome"][file.endswith(".genome.bam")]
		ref_file = [refTranscGRCh38, refGenomeGRCh38][file.endswith(".genome.bam")]
		
		samtools_index = " ".join([
		"samtools index",  # Indexing the concat_samples.bam file
		"-@", args.threads,  # Number of threads to be used
		file])  # Input BAM file
		subprocess.run(samtools_index, shell=True)

		# BAM stats
		bam_stat = ' '.join([
		"bam_stat.py",
		"-i", file,  # Input BAM file
		"> {0}/{1}.{2}.bamstat.txt".format(postanalysis_dir, file_name, aligned_to),  # Output file
		"2>>", os.path.join(postanalysis_dir, "{0}_{1}_bamstats-report.txt".format(file_name, aligned_to))])
		subprocess.run(bam_stat, shell=True)

		# Picard CollectAlignmentSummaryMetrics
		CollectAlignmentSummaryMetrics = ' '.join([
		"picard CollectAlignmentSummaryMetrics",  # Call picard CollectAlignmentSummaryMetrics
		"INPUT= {0}".format(file),  # Input BAM file
		"OUTPUT= {0}/{1}.{2}.alignment_metrics.txt".format(postanalysis_dir, file_name, aligned_to),  # Output
		"REFERENCE_SEQUENCE= {0}".format(ref_file),  # Reference sequence file
		"2>>", os.path.join(postanalysis_dir, "{0}_{1}_CollectAlignmentSummaryMetrics-report.txt".format(file_name, aligned_to))])
		subprocess.run(CollectAlignmentSummaryMetrics, shell=True) 

		# Wub Alignment based QC plots
		alignment_qc = ' '.join([
		"bam_alignment_qc.py",
		"-f", ref_file,  # Input reference file
		"-x",  # Do not plot per-reference information
		"-r", "{0}/{1}.{2}.bam_alignment_qc.pdf".format(postanalysis_dir, file_name, aligned_to),  # Output pdf file
		"-p", "{0}/{1}.{2}.bam_alignment_qc.pk".format(postanalysis_dir, file_name, aligned_to),  # Output pk file
		file,
		"> {0}/{1}.{2}.alignment_qc.txt".format(postanalysis_dir, file_name, aligned_to),  # Output file
		"2>>", os.path.join(postanalysis_dir, "{0}_{1}_alignment_qc-report.txt".format(file_name, aligned_to))])
		subprocess.run(alignment_qc, shell=True)


		### GENOME ALIGNMENTS STATS 
		if file.endswith(".genome.bam"):
			# BAM read distribution
			if i==1: print(">>> RSeQC - Generating read distribution stats: in progress ..")
			read_distribution = ' '.join([
			"read_distribution.py", 
			"-i", file,  # Input BAM file
			"-r", reference_annotation_bed,
			"> {0}/{1}.{2}.fragSize".format(postanalysis_dir, file_name, aligned_to),  # Output file
			"2>>", os.path.join(postanalysis_dir, "{0}_{1}_read_distribution-report.txt".format(file_name, aligned_to))])
			subprocess.run(read_distribution, shell=True)

			# Check the strandness of the reads
			if i==1: print(">>> RSeQC - Generating read strandness stats: in progress ..")
			strandness = ' '.join([
			"infer_experiment.py",  # Call samtools flagstat
			"-i", file,  # Input BAM file
			"-r", reference_annotation_bed,
			"> {0}/{1}.{2}.strandness.txt".format(postanalysis_dir, file_name, aligned_to),  # Output file
			"2>>", os.path.join(postanalysis_dir, "{0}_{1}_strandness-report.txt".format(file_name, aligned_to))])
			subprocess.run(strandness, shell=True)

			# Gene body coverage
			if i==1: print(">>> RSeQC - Generating gene body coverage stats: in progress ..")
			gene_coverage = ' '.join([
			"geneBody_coverage.py",  # Call samtools flagstat
			"-r", reference_annotation_bed,
			"-i", file,  # Input BAM file
			"-o", "{0}/{1}.{2}.gene_coverage".format(postanalysis_dir, file_name, aligned_to),  # Output file
			"2>>", os.path.join(postanalysis_dir, "{0}_{1}_mapping_position-report.txt".format(file_name, aligned_to))])
			subprocess.run(gene_coverage, shell=True)

			# Check duplicate reads
			if i==1: print(">>> Picard - Extracting read duplication stats: in progress ..")
			duplicate_reads = ' '.join([
			"picard MarkDuplicates",  # Call samtools flagstat
			"INPUT= {0}".format(file),  # Input BAM file
			"OUTPUT= {0}/{1}.genome.dedup.bam".format(alignments_dir, file_name),
			"METRICS_FILE= {0}/{1}.{2}.mark_duplicates.txt".format(postanalysis_dir, file_name, aligned_to),  # Output file
			"2>>", os.path.join(postanalysis_dir, "{0}_{1}_duplicate_reads-report.txt".format(file_name, aligned_to))])
			subprocess.run(duplicate_reads, shell=True)
			os.system('rm {0}/{1}.genome.dedup.bam'.format(alignments_dir, file_name))

			# Number of reads mapped to each chromosome
			if i==1: print(">>> RSeQC - Generating mapping stats: in progress ..")
			mapping_pos = ' '.join([
			"samtools idxstats",  # Call samtools flagstat
			"--threads", args.threads,
			file,  # Input BAM file
			"> {0}/{1}.{2}.samtools_idxstats.txt".format(postanalysis_dir, file_name, aligned_to),  # Output file
			"2>>", os.path.join(postanalysis_dir, "{0}_{1}_mapping_position-report.txt".format(file_name, aligned_to))])
			subprocess.run(mapping_pos, shell=True)
	return

def generate_expression_matrices():
	
	genome_alignments = [sum_file for sum_file in glob.glob(os.path.join(alignments_dir, "*.genome.bam"))]
	print("FeatureCounts - Counting reads from the genome aligned data: in progress ..")
	featureCounts_gn = " ".join([
	"featureCounts",  # Call featureCounts
	"-T", args.threads,  # Number of threads to be used by the script
	"-a", refAnnot,  # Annotation file in GTF/GFF format
	"-L",  # Count long reads such as Nanopore and PacBio reads
	"-o", os.path.join(postanalysis_dir, "genome_alignments_sum.tab"),
	' '.join(genome_alignments),  # Input bam file
	"2>>", os.path.join(postanalysis_dir, "featureCounts_genome_summary-report.txt")]) 
	subprocess.run(featureCounts_gn, shell=True)
	
	print("Exporting the (genomic) expression matrix: in progress ..")
	subprocess.run("cut -f1,7- {0}/genome_alignments_sum.tab | sed 1d > {0}/featureCounts_expression_matrix.txt".format(postanalysis_dir), shell=True)

	# transcriptome_alignments = [sum_file for sum_file in glob.glob(os.path.join(alignments_dir, "*.transcriptome.bam"))]
	# print("Cufflinks - Counting reads from the transcriptome aligned data: in progress ..")
	# for file in transcriptome_alignments:
	# 	file_name = os.path.basename(files.split(".")[0])
	# 	cufflinks = " ".join([
	# 	"cufflinks",  # Call pycoQC
	# 	"--num-threads", args.threads,  # Number of threads to be used by the script
	# 	"--output-dir", os.path.join(postanalysis_dir, file_name),  # Output directory
	# 	"--GTF", refAnnot,  # Annotation file in GTF format
	# 	file,
	# 	# "2>>", os.path.join(reports_dir, "cufflinks_transcriptome_summary-report.txt")
	# 	])
	# 	# subprocess.run(cufflinks, shell=True)
	
	# print("Exporting the (transcriptomic) expression matrix: in progress ..")
	# subprocess.run("ls {0}/*.gtf > {0}/cufflinks_outputs.txt".format(postanalysis_dir), shell=True)
	# subprocess.run("cuffmerge --num-threads {0} --ref-gtf {1} --ref-sequence {2} -o {3} {3}/cufflinks_outputs.txt"\
	# 							   .format(args.threads, refAnnot, refGenomeGRCh38, postanalysis_dir), shell=True)
	return

def summary(num_of_files):
	print(">>> multiQC - Summing all QC reports: in progress ..")
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
	"-r", "{0}/{1}.{2}.bam_alignment_qc.pdf".format(postanalysis_dir, file_name, aligned_to),  # Output pdf file
	' '.join(pickle_files),
	"> {0}/{1}.{2}.bam_multi_qc.txt".format(postanalysis_dir, file_name, aligned_to),  # Output file
	"2>>", os.path.join(postanalysis_dir, "{0}_{1}_bam_multi_qc-report.txt".format(file_name, aligned_to))])
	subprocess.run(bam_multi_qc, shell=True)

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
	os.system('mv {0}/*alignment_qc.pk {1}'.format(postanalysis_dir, qc_reports))
	os.system('mv {0}/*log* {1}'.format(current_dir, qc_reports))
	os.system('mv {0}/*.fragSize {1}'.format(postanalysis_dir, qc_reports))
	os.system('mv {0}/*sum.* {1}'.format(postanalysis_dir, qc_reports))
	os.system('mv {0}/*curves.pdf {1}'.format(postanalysis_dir, qc_reports))
	os.system('mv {0}/*genome*.txt {1}'.format(postanalysis_dir, qc_reports))
	os.system('mv {0}/*transcriptome*.txt {1}'.format(postanalysis_dir, qc_reports))
	os.system('mv {0}/*summarised_report_data {1}'.format(postanalysis_dir, qc_reports))
	return

def main():
	
	summary_files = [file_path for file_path in Path(ont_data).glob('**/*_sequencing_summary.txt')]

	for i, sum_file in enumerate(summary_files, 1):
		sample_id = sum_file.parts[4].split("_")[-1]
		raw_data_dir = str(sum_file.parents[1])
		sum_file = str(sum_file)

		quality_control(i, sum_file, sample_id, raw_data_dir)

		alignment_against_ref(i, sample_id, raw_data_dir)

	mapping_qc()

	generate_expression_matrices()

	summary(len(summary_files))

if __name__ == "__main__": main()