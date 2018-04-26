import math
"""Generate all the meta data files"""

#Which rules will be run on the host computer and not sent to nodes
localrules: create_dict, reduce_gtf, create_refFlat, create_intervals, create_star_index

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
		dropseq_wrapper=config['LOCAL']['dropseq-wrapper'],
		memory=config['LOCAL']['memory']
	output:
		'{}.reduced.gtf'.format(annotation_prefix)
	shell:
		"""{params.dropseq_wrapper} -m {params.memory} -p ReduceGTF\
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
		dropseq_wrapper=config['LOCAL']['dropseq-wrapper'],
		memory=config['LOCAL']['memory']
	output:
		'{}.refFlat'.format(annotation_prefix)
	shell:
		"""{params.dropseq_wrapper} -m {params.memory} -p ConvertToRefFlat\
		ANNOTATIONS_FILE={input.annotation}\
		OUTPUT={output}\
		SEQUENCE_DICTIONARY={input.reference_dict}
		"""

rule create_intervals:
	input:
		annotation_reduced='{}.reduced.gtf'.format(annotation_prefix),
		reference_dict='{}.dict'.format(reference_prefix)
	params:
		dropseq_wrapper=config['LOCAL']['dropseq-wrapper'],
		memory=config['LOCAL']['memory'],
		reference_directory=config['META']['reference-directory'],
		reference_prefix=re.split(".fasta|.fa",config['META']['reference-file'])[0]
		
	output:
		'{}.rRNA.intervals'.format(reference_prefix)
	shell:
		"""{params.dropseq_wrapper} -m {params.memory} -p CreateIntervalsFiles\
		REDUCED_GTF={input.annotation_reduced}\
		SEQUENCE_DICTIONARY={input.reference_dict}\
		O={params.reference_directory}\
		PREFIX={params.reference_prefix}
		"""
def get_sjdbOverhang(wildcards):
	return(int(wildcards.read_length)-1)

def get_genomeChrBinNbits(file):
	genomeLength = shell("wc -c {} | cut -d' ' -f1".format(file), iterable=True)
	genomeLength = int(next(genomeLength))
	referenceNumber = shell('grep "^>" {} | wc -l'.format(file), iterable=True)
	referenceNumber = int(next(referenceNumber))
	return(min([18,int(math.log2(genomeLength/referenceNumber))]))


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
	threads: 4
	conda: '../envs/star.yaml'
	shell:
		"""mkdir -p {params.genomeDir}; STAR\
		--runThreadN 4\
		--runMode genomeGenerate\
		--genomeDir {params.genomeDir}\
		--genomeFastaFiles {input.reference_file}\
		--sjdbGTFfile {input.annotation_file}\
		--limitGenomeGenerateRAM 30000000000\
		--sjdbOverhang {params.sjdbOverhang}\
		--genomeChrBinNbits {params.genomeChrBinNbits}
		"""