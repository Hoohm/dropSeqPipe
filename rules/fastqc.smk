"""Get fastqc reports"""

#Which rules will be run on the host computer and not sent to nodes
localrules:
    multiqc_fastqc_reads,
    multiqc_fastqc_barcodes

rule fastqc_barcodes:
    """Create fastqc report"""
    input: 
        get_R1_files,
        'fastqc_adapter.tsv',
    output:
        html='{results_dir}/logs/fastqc/{sample}_R1_fastqc.html',
        zip='{results_dir}/logs/fastqc/{sample}_R1_fastqc.zip'
    params: '--extract -a fastqc_adapter.tsv'
    wrapper:
        '0.27.1/bio/fastqc'

rule fastqc_reads:
    """Create fastqc report"""
    input: 
        get_R2_files,
        'fastqc_adapter.tsv',
    output:
        html='{results_dir}/logs/fastqc/{sample}_R2_fastqc.html',
        zip='{results_dir}/logs/fastqc/{sample}_R2_fastqc.zip'
    params: '--extract -a fastqc_adapter.tsv'
    wrapper:
        '0.27.1/bio/fastqc'


rule multiqc_fastqc_barcodes:
    input:
        expand('{results_dir}/logs/fastqc/{sample}_R1_fastqc.html', sample=samples.index, results_dir=results_dir)
    output:
        html='{results_dir}/reports/fastqc_barcodes.html'
    params: '-m fastqc --ignore *_R2*'
    wrapper:
        '0.27.1/bio/multiqc'

rule multiqc_fastqc_reads:
    input: 
        expand('{results_dir}/logs/fastqc/{sample}_R2_fastqc.html', sample=samples.index, results_dir=results_dir)
    output:
        html='{results_dir}/reports/fastqc_reads.html'
    params: '-m fastqc --ignore *_R1*'
    wrapper:
        '0.27.1/bio/multiqc'

rule fasta_fastq_adapter:
    input:
        adapters=config['FILTER']['cutadapt']['adapters-file']
    output:
        file="fastqc_adapter.tsv"
    run:
        import sys
        from Bio import SeqIO
        with open( output.file, "w" ) as myoutput:
            for seq_record in SeqIO.parse(input.adapters, "fasta"):
                myline = (str(seq_record.id)) + "\t" + str(seq_record.seq[0:13]) + "\n"
                myoutput.write(myline)
