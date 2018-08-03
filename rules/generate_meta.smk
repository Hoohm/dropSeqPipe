import math
import platform
"""Generate all the meta data files"""

#Which rules will be run on the host computer and not sent to nodes
localrules: create_dict, reduce_gtf, create_refFlat, create_intervals

rule create_dict:
	input:
		reference_file
	output:
		'{}.dict'.format(reference_prefix)
	threads:1
	params:
		picard="$CONDA_PREFIX/share/picard-2.14.1-0/picard.jar",
		temp_directory=config['LOCAL']['temp-directory']
	conda: '../envs/picard.yaml'
	shell:
		"""java -jar -Djava.io.tmpdir={params.temp_directory} {params.picard} CreateSequenceDictionary\
		REFERENCE={input}\
		OUTPUT={output}
		"""

rule reduce_gtf:
	input:
		reference_dict='{}.dict'.format(reference_prefix),
		annotation=annotation_file
	params:
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	output:
		'{}.reduced.gtf'.format(annotation_prefix)
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && ReduceGTF -m {params.memory}\
		GTF={input.annotation}\
		OUTPUT={output}\
		SEQUENCE_DICTIONARY={input.reference_dict}\
		IGNORE_FUNC_TYPE='null'\
		ENHANCE_GTF='false'
		"""

rule create_refFlat:
	input:
		annotation=annotation_file,
		reference_dict='{}.dict'.format(reference_prefix)
	params:
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	output:
		'{}.refFlat'.format(annotation_prefix)
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && ConvertToRefFlat -m {params.memory}\
		ANNOTATIONS_FILE={input.annotation}\
		OUTPUT={output}\
		SEQUENCE_DICTIONARY={input.reference_dict}
		"""

rule create_intervals:
	input:
		annotation_reduced='{}.reduced.gtf'.format(annotation_prefix),
		reference_dict='{}.dict'.format(reference_prefix)
	params:
		memory=config['LOCAL']['memory'],
		reference_directory=config['META']['reference-directory'],
		reference_prefix=re.split(".fasta|.fa",config['META']['reference-file'])[0],
		temp_directory=config['LOCAL']['temp-directory']
		
	output:
		'{}.rRNA.intervals'.format(reference_prefix)
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && CreateIntervalsFiles -m {params.memory}\
		REDUCED_GTF={input.annotation_reduced}\
		SEQUENCE_DICTIONARY={input.reference_dict}\
		O={params.reference_directory}\
		PREFIX={params.reference_prefix}
		"""
def get_sjdbOverhang(wildcards):
	return(int(wildcards.read_length)-1)

def get_genomeChrBinNbits(file):
	if (platform.system() == 'Darwin'):
		genomeLength = shell("wc -c {} | cut -d' ' -f2".format(file), iterable=True)
	else:
		genomeLength = shell("wc -c {} | cut -d' ' -f1".format(file), iterable=True)
	genomeLength = int(next(genomeLength))
	referenceNumber = shell('grep "^>" {} | wc -l'.format(file), iterable=True)
	referenceNumber = int(next(referenceNumber))
	return(min([18,int(math.log2(genomeLength/referenceNumber))]))

rule get_genomeChrBinNbits:
	input:
		reference_file=reference_file
	params:
		samples_file='samples.csv',
		reference_directory=config['META']['reference-directory']
	output:
		'{params.reference_directory}/index_params.txt'
	run:
		"""
		from math import log2
		from platform import system
		if (system() == 'Darwin'):
			genomeLength = shell("wc -c {} | cut -d' ' -f2".format(snakemake.reference_file), iterable=True)
		else:
			genomeLength = shell("wc -c {} | cut -d' ' -f1".format(snakemake.reference_file), iterable=True)
		genomeLength = int(next(genomeLength))
		referenceNumber = shell('grep "^>" {} | wc -l'.format(snakemake.reference_file), iterable=True)
		referenceNumber = int(next(referenceNumber))
		value = min([18,int(log2(genomeLength/referenceNumber))])
		"""

rule create_star_index:
	input:
		reference_file=reference_file,
		annotation_file=annotation_file
	params:
		sjdbOverhang=lambda wildcards: get_sjdbOverhang(wildcards),
		genomeDir='{star_index_prefix}_{read_length}',
		genomeChrBinNbits=get_genomeChrBinNbits(reference_file)
	output:
		'{star_index_prefix}_{read_length}/SA'
	threads: 12
	conda: '../envs/star.yaml'
	shell:
		"""mkdir -p {params.genomeDir}; STAR\
		--runThreadN {threads}\
		--runMode genomeGenerate\
		--genomeDir {params.genomeDir}\
		--genomeFastaFiles {input.reference_file}\
		--sjdbGTFfile {input.annotation_file}\
		--limitGenomeGenerateRAM 30000000000\
		--sjdbOverhang {params.sjdbOverhang}\
		--genomeChrBinNbits {params.genomeChrBinNbits}
		"""