"""Generate all the meta data files"""


rule create_dict:
	input: reference_file
	output: '{reference_prefix}.dict'
	threads:1
	params:
		picard = "$CONDA_PREFIX/share/picard-2.14.1-0/picard.jar",
		TMPDIR = config['LOCAL']['TMPDIR']
	shell:
		"""java -jar -Djava.io.tmpdir={params.TMPDIR} {params.picard} CreateSequenceDictionary\
		REFERENCE={input}\
		OUTPUT={output}
		"""

rule reduce_gtf:
	input:
		reference_dict=expand('{reference_prefix}.dict', reference_prefix=reference_prefix),
		annotation=annotation_file
	params:
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY']
	output: '{annotation_prefix}.reduced.gtf'
	threads:1
	shell:
		"""{params.DROPSEQ_wrapper} -m {params.memory} -p ReduceGTF\
		GTF={input.annotation}\
		OUTPUT={output}\
		SEQUENCE_DICTIONARY={input.reference_dict}\
		IGNORE_FUNC_TYPE='null'\
		ENHANCE_GTF='false'
		"""

rule create_refFlat:
	input:
		annotation=annotation_file,
		reference_dict=expand('{reference_prefix}.dict', reference_prefix=reference_prefix)
	params:
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY']
	output: '{annotation_prefix}.refFlat'
	threads: 1
	shell:
		"""{params.DROPSEQ_wrapper} -m {params.memory} -p ConvertToRefFlat\
		ANNOTATIONS_FILE={input.annotation}\
		OUTPUT={output}\
		SEQUENCE_DICTIONARY={input.reference_dict}
		"""

rule create_intervals:
	input:
		annotation_reduced = expand('{annotation_prefix}.reduced.gtf', annotation_prefix=annotation_prefix),
		reference_dict = expand('{reference_prefix}.dict', reference_prefix=reference_prefix)
	params:
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY'],
		reference_folder = config['META']['reference_folder'],
		reference_prefix = config['META']['reference_file'].split('.fasta')[0]
	output:
		'{reference_prefix}.rRNA.intervals'
	threads: 1
	shell:
		"""{params.DROPSEQ_wrapper} -m {params.memory} -p CreateIntervalsFiles\
		REDUCED_GTF={input.annotation_reduced}\
		SEQUENCE_DICTIONARY={input.reference_dict}\
		O={params.reference_folder}\
		PREFIX={params.reference_prefix}
		"""
def get_sjdbOverhang(wildcards):
	return(int(wildcards.read_length)-1)

rule create_star_index:
	input:
		reference_file = reference_file,
		annotation_file = annotation_file
	params:
		sjdbOverhang = lambda wildcards: get_sjdbOverhang(wildcards),
		genomeDir = '{star_index_prefix}_{read_length}'
	output:
		'{star_index_prefix}_{read_length}/SA'
	threads: 4
	shell:
		"""STAR\
		--runThreadN 4\
		--runMode genomeGenerate\
		--genomeDir {params.genomeDir}\
		--genomeFastaFiles {input.reference_file}\
		--sjdbGTFfile {input.annotation_file}\
		--limitGenomeGenerateRAM 30000000000\
		--sjdbOverhang {params.sjdbOverhang}
		"""