import re
import glob
import gzip
from collections import defaultdict

multi_lane_pattern = re.compile("../data\/(.*)_(L[0-9]{3})_(R[1-2])_001.fastq.gz")


def get_input_files(wildcards):
    samples = [f for f in glob.glob("../{results_dir}/samples/*.fastq.gz") if re.match(multi_lane_pattern,f)]
    return(samples)

lanes = sorted(list(set([re.findall(multi_lane_pattern,f)[0][1] for f in glob.glob("../{results_dir}/samples/*.fastq.gz") if re.match(multi_lane_pattern,f)])))
samples = [re.findall(multi_lane_pattern,f)[0][0] for f in glob.glob("../{results_dir}/samples/*.fastq.gz") if re.match(multi_lane_pattern,f)]




rule all:
    input:
        expand('{results_dir}/samples/{sample}_R1.fastq.gz',sample=samples),
        expand('{results_dir}/samples/{sample}_R2.fastq.gz',sample=samples)


rule generate_samples:
    input:
        get_input_files
    output:
        'samples.csv'
    run:
        samples = defaultdict(lambda: {'sample_lanes':[],'read_length':0})
        with open(output[0],'w') as sample_file:
            sample_file.write("samples,expected_cells,read_length,batch\n")
            for file in input:
                if('R2' in file):
                    with gzip.open(file) as fastq_file:
                        next(fastq_file)
                        read_length = len(next(fastq_file).strip())
                        re_results = re.findall(multi_lane_pattern,file)
                        samples[re_results[0][0]]['sample_lanes'].append(re_results[0][1])
                        samples[re_results[0][0]]['read_length']=read_length
            for sample_name in samples:
                sample_file.write("{},,{},\n".format(sample_name, read_length))

rule concat_lanes:
    input:
        R1=expand('{results_dir}/samples/{{sample}}_{lane}_R1_001.fastq.gz', lane=lanes),
        R2=expand('{results_dir}/samples/{{sample}}_{lane}_R2_001.fastq.gz', lane=lanes),
        lanes='samples.csv'
    output:
        R1='{results_dir}/samples/{sample}_R1.fastq.gz',
        R2='{results_dir}/samples/{sample}_R2.fastq.gz'
    shell:
        """cat {input.R1} > {output.R1};cat {input.R2} > {output.R2}"""