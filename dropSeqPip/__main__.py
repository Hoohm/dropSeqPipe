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
                        required=True)
    args = parser.parse_args()
    return args


def main():
    """Main."""
    scripts_dir = os.path.dirname(__file__)
    args = get_args()
    configfile = os.path.join(args.folder_path, 'config.yaml')
    sub_folders = ['summary', 'logs', 'plots']
    package_dir = os.path.dirname(__file__)

    #if (not os.path.isfile(configfile)):
        #print('No configfile. Exiting')
        #sys.exit()
    #Create folder if not present
    for folder in sub_folders:
        joined = os.path.join(args.folder_path, folder)
        if(not os.path.isdir(joined)):
            os.mkdir(joined)
    #Load config files
    with open(args.config_file_path) as config_yaml:
        yaml_data = yaml.load(config_yaml)
    with open(os.path.join(args.folder_path, 'config.yaml')) as samples_config:
        samples_yaml = yaml.load(samples_config)
    #Select step and run
    step_list = []
    if(args.mode == "pre-process"):
        print("Mode is {}.".format(args.mode))
        pre_align = ('snakemake -s {}/Snakefiles/pre_align.snake --cores {} -pT -d {} --configfile {}'.format(scripts_dir, yaml_data['CORES'], args.folder_path, args.config_file_path),'Pre-processing before alignement')
        star_align = ('snakemake -s {}/Snakefiles/star_align.snake --cores {} -pT -d {} --configfile {}'.format(scripts_dir, yaml_data['CORES'], args.folder_path, args.config_file_path),'Star alignement')
        post_align = ('snakemake -s {}/Snakefiles/post_align.snake --cores {} -pT -d {} --configfile {}'.format(scripts_dir, yaml_data['CORES'], args.folder_path, args.config_file_path),'Post alignement')

        print('Running pre-processing')
        subprocess.call(pre_align, shell=True)
        print('Running Alignement')
        subprocess.call(star_align, shell=True)
        print('Running post-alignement')
        subprocess.call(post_align, shell=True)
    if(args.mode == "knee-plot"):
        print(package_dir)
        knee_plot = 'Rscript {}/Rscripts/knee_plot.R {}'.format(package_dir, args.folder_path)
        print('Plotting knee plots')
        subprocess.call(knee_plot, shell=True)
    if(args.mode == "species-plot"):
        if(len(samples_yaml['SPECIES']) == 2):
            extract_species = ('snakemake -s {}/Snakefiles/extract_species.snake --cores {} -pT -d {} --configfile {}'.format(scripts_dir, yaml_data['CORES'], args.folder_path, args.config_file_path), 'Extracting species')
            print('Extracting species')
            subprocess.call(extract_species, shell=True)
            species_plot = 'Rscript {}/Rscripts/species_plot.R {}'.format(package_dir, args.folder_path)
            print('Plotting species plots')
            subprocess.call(species_plot, shell=True)
        else:
            print('You cannot run this with a number of species different than 2.\nPlease change the config file')
    if(args.mode == "extract-expression"):
        if(len(samples_yaml['SPECIES']) == 2):
            extract_expression = ('snakemake -s {}/Snakefiles/extract_expression.snake --cores {} -pT -d {} --configfile {}'.format(scripts_dir, yaml_data['CORES'], args.folder_path, args.config_file_path), 'Extracting species')  
            print('Extracting expression')
            subprocess.call(extract_expression, shell=True)
    print('Pipeline finished')
if __name__ == "__main__":
    sys.exit(main())
