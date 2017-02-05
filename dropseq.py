import subprocess
import sys
import os
import json

try:
	sys.argv[1]
except:
	print('You need to provide a folder on which to run on. Exiting')
	exit()
wrkdir = sys.argv[1]
configfile = os.path.join(wrkdir, 'config.json')
local_config_file = 'local.json'
FOLDERS = ['summary', 'logs', 'plots']

if (not os.path.isfile(configfile)):
	print('No configfile. Exiting')
	exit()
for folder in FOLDERS:
	joined = os.path.join(wrkdir, folder)
	if(not os.path.isdir(joined)):
		os.mkdir(joined)

with open(local_config_file) as config_json:
    json_data = json.load(config_json)

#Optimized Run for dropseq
first = 'snakemake -s snakefiles/Dropseq_pre_align.snake --cores {} -pT -d {} --configfile local.json'.format(json_data['CORES'], sys.argv[1])
second = 'snakemake -s snakefiles/Star_align.snake --cores {} -pT -d {} --configfile local.json'.format(json_data['CORES'], sys.argv[1])
third = 'snakemake -s snakefiles/Dropseq_post_align.snake --cores {} -pT -d {} --configfile local.json'.format(json_data['CORES'], sys.argv[1])
knee_plot = 'Rscript Rscripts/knee_plot.R {}'.format(sys.argv[1])
print('Running pre_processing')
subprocess.call(first, shell=True)
print('Running Alignement')
subprocess.call(second, shell=True)
print('Running post_processing')
subprocess.call(third, shell=True)
print('Plotting knee_plots')
subprocess.call(knee_plot, shell=True)
print('Pipeline finished')
