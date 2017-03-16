#! /usr/binPython3
# coding=utf-8

"""Main."""

import argparse
import sys
import os
import subprocess
import yaml


def get_args():
    """Get args."""
    parser = argparse.ArgumentParser(description='dropSeqPip')
    parser.add_argument('-c', '--config_file_path',
                        action='store',
                        help='path to local config file',
                        dest='config_file_path',
                        required=True)
    parser.add_argument('-f', '--folder',
                        help='Folder containing the samples',
                        dest='folder_path',
                        required=True)
    parser.add_argument('-m', '--mode',
                        help='Which mode to run.',
                        choices=['pre-process', 'knee-plot', 'species-plot', 'extract-expression'],
                        action='store',
                        nargs='+',
                        required=True)
    parser.add_argument('--rerun',
                        action='store_true',
                        help='forces a rerun of the selected modes',
                        default=False)
    args = parser.parse_args()
    return args


def main():
    """Main."""
    scripts_dir = os.path.dirname(__file__)
    args = get_args()
    configfile = os.path.join(args.folder_path, 'config.yaml')
    sub_folders = ['summary', 'logs', 'plots']
    package_dir = os.path.dirname(__file__)
    for folder in sub_folders:
        joined = os.path.join(args.folder_path, folder)
        if(not os.path.isdir(joined)):
            os.mkdir(joined)
    if(args.rerun):
        rerun = '--forceall'
    else:
        rerun = ''
    #Load config files
    with open(args.config_file_path) as config_yaml:
        yaml_data = yaml.load(config_yaml)
    with open(os.path.join(args.folder_path, 'config.yaml')) as samples_config:
        samples_yaml = yaml.load(samples_config)
    #Select step and run
    step_list = []
    if("pre-process" in args.mode):
        print("Mode is {}.".format(args.mode))
        pre_align = ('snakemake -s {}/Snakefiles/pre_align.snake --cores {} -pT -d {} --configfile {} {}'.format(scripts_dir, yaml_data['CORES'], args.folder_path, args.config_file_path,rerun),'Pre-processing before alignement')
        star_align = ('snakemake -s {}/Snakefiles/star_align.snake --cores {} -pT -d {} --configfile {} {}'.format(scripts_dir, yaml_data['CORES'], args.folder_path, args.config_file_path, rerun),'Star alignement')
        post_align = ('snakemake -s {}/Snakefiles/post_align.snake --cores {} -pT -d {} --configfile {} {}'.format(scripts_dir, yaml_data['CORES'], args.folder_path, args.config_file_path, rerun),'Post alignement')

        print('Running pre-processing')
        subprocess.call(pre_align, shell=True)
        print('Running Alignement')
        subprocess.call(star_align, shell=True)
        print('Running post-alignement')
        subprocess.call(post_align, shell=True)
    if("knee-plot" in args.mode):
        knee_plot = 'Rscript {}/Rscripts/knee_plot.R {}'.format(package_dir, args.folder_path)
        print('Plotting knee plots')
        subprocess.call(knee_plot, shell=True)
    if("species-plot" in args.mode):
        if(len(samples_yaml['SPECIES']) == 2):
            extract_species = ('snakemake -s {}/Snakefiles/extract_species.snake --cores {} -pT -d {} --configfile {} {}'.format(scripts_dir, yaml_data['CORES'], args.folder_path, args.config_file_path, rerun), 'Extracting species')
            print('Extracting species')
            subprocess.call(extract_species, shell=True)
            species_plot = 'Rscript {}/Rscripts/species_plot.R {}'.format(package_dir, args.folder_path)
            print('Plotting species plots')
            subprocess.call(species_plot, shell=True)
        else:
            print('You cannot run this with a number of species different than 2.\nPlease change the config file')
    if("extract-expression" in args.mode):
        if(len(samples_yaml['SPECIES']) == 2):
            extract_expression = ('snakemake -s {}/Snakefiles/extract_expression.snake --cores {} -pT -d {} --configfile {} {}'.format(scripts_dir, yaml_data['CORES'], args.folder_path, args.config_file_path, rerun), 'Extracting species')  
            print('Extracting expression')
            subprocess.call(extract_expression, shell=True)
        if(len(samples_yaml['SPECIES']) == 1):
            extract_expression_single = ('snakemake -s {}/Snakefiles/extract_expression_single.snake --cores {} -pT -d {} --configfile {} {}'.format(scripts_dir, yaml_data['CORES'], args.folder_path, args.config_file_path, rerun), 'Extracting species')  
            print('Extracting expression')
            subprocess.call(extract_expression_single, shell=True)
    print('Pipeline finished')
if __name__ == "__main__":
    sys.exit(main())
