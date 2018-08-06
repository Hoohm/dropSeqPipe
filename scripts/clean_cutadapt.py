from collections import defaultdict
import re

def fill_results(snakemake,pair, adapter_results):
    adapter_pattern = re.compile(pattern="=== Adapter (.*) ===\n")
    with open(snakemake.input[pair], 'r') as logfile:
        line = logfile.readline()
        while(line):
            adapter_matched = re.findall(pattern=adapter_pattern,string=line)
            if(adapter_matched):
                logfile.readline()
                line = logfile.readline().rstrip('.\n')
                line_list = line.split(';')
                adapter_results[adapter_matched[0]]['Pair'] = pair
                adapter_results[adapter_matched[0]]['Sequence'] = line_list[0].split(':')[1].strip()
                if(adapter_matched[0] in adapter_results):
                    adapter_results[adapter_matched[0]]['Times'] += int(line_list[3].split(':')[1].split(' ')[1].strip())
                else:    
                    adapter_results[adapter_matched[0]]['Times'] += int(line_list[3].split(':')[1].split(' ')[1].strip())
            line = logfile.readline()
    return(adapter_results)

adapter_results_R1 = defaultdict(lambda :{'Pair':'','Sequence':'','Times':0})
adapter_results_R2 = defaultdict(lambda :{'Pair':'','Sequence':'','Times':0})


adapter_results_R1 = fill_results(snakemake, 'R1', adapter_results_R1)
adapter_results_R2 = fill_results(snakemake, 'R2', adapter_results_R2)

with open(snakemake.output[0],'w') as outfile:
    outfile.write('Adapter,Sequence,Pair,Count\n')
    for adapter in adapter_results_R1:
        outfile.write("{},{},{},{}\n".format(adapter, adapter_results_R1[adapter]['Sequence'],adapter_results_R1[adapter]['Pair'],adapter_results_R1[adapter]['Times']))
    for adapter in adapter_results_R2:
        outfile.write("{},{},{},{}\n".format(adapter, adapter_results_R2[adapter]['Sequence'],adapter_results_R2[adapter]['Pair'],adapter_results_R2[adapter]['Times']))