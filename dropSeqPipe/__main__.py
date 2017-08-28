#! /usr/binPython3
# coding=utf-8

"""Main."""
import argparse
import sys
import os
import yaml
from snakemake.shell import shell


def get_args():
    """Get args."""
    parser = argparse.ArgumentParser(description='dropSeqPipe')
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
                        choices=[
                            'pre-process',
                            'generate-plots',
                            'species-plot',
                            'extract-expression',
                            'fastqc',
                            'generate-meta',
                            'test'],
                        action='store',
                        nargs='+',
                        required=True)
    parser.add_argument('--rerun',
                        action='store_true',
                        help='forces a rerun of the selected modes',
                        default=False)
    parser.add_argument('--notemp',
                        action='store_true',
                        help='Keeps the temp files',
                        default=False)
    parser.add_argument('-cc','--cluster_config',
                        action='store',
                        dest='cluster_config_file',
                        help='Runs the pipeline with cluster jobs instead of local commands')
    parser.add_argument('--dryrun',
                        action='store_true',
                        help='Runs only as dryrun. For testing purposes',
                        default=False)

    args = parser.parse_args()
    return args


def main():
    """Main."""
    scripts_dir = os.path.dirname(__file__)
    args = get_args()
    complementory_args = ''
    if(args.rerun):
        complementory_args += '--forceall '
    if(args.notemp):
        complementory_args += '--notemp '
    if(args.dryrun):
        complementory_args += '-n '
    if(args.cluster_config_file):
        complementory_args += '--cluster-config {0} '.format(args.cluster_config_file)
        complementory_args += '--cluster "qsub -v PATH={{cluster.PATH}} -l mem={{cluster.mem}} -e {{cluster.error}} -o {{cluster.output}}" '
    # Load config files
    with open(args.config_file_path) as config_yaml:
        yaml_data = yaml.load(config_yaml)
    if("generate-meta" in args.mode):
        shell('snakemake -s {0}/Snakefiles/generate_meta.snake --cores {1} -pT -d {2} --configfile {3} {4}'.format(
            scripts_dir,
            yaml_data['CORES'],
            args.folder_path,
            args.config_file_path,
            complementory_args))
    else:
        try:
            with open(os.path.join(args.folder_path, 'config.yaml')) as samples_config:
                samples_yaml = yaml.load(samples_config)
        except:
            sys.exit('Samples configuration file not found or not formatted properly. Exiting')
        sub_folders = ['summary', 'logs', 'plots']
        if(args.cluster_config_file):
            sub_folders.append('logs/cluster')
        package_dir = os.path.dirname(__file__)
        for folder in sub_folders:
            joined = os.path.join(args.folder_path, folder)
            if(not os.path.isdir(joined)):
                os.mkdir(joined)
    
    # Select step and run
    # First step is fastqc
    if("fastqc" in args.mode):
        print("Mode is fastqc.")
        fastqc = 'snakemake -s {}/Snakefiles/{}/fastqc.snake --cores {} -pT -d {} --configfile {} {}'.format(
            scripts_dir,
            samples_yaml['GLOBAL']['data_type'],
            yaml_data['CORES'],
            args.folder_path,
            args.config_file_path,
            complementory_args)
        fastqc_summary = 'Rscript {}/Rscripts/fastqc.R {}'.format(
            package_dir,
            args.folder_path)
        print('Running fastqc')
        shell(fastqc)
        shell(fastqc_summary)
    # Second step is pre-processing which also does aligning and post processing
    if("pre-process" in args.mode):
        print("Mode is pre-processing.")
        pre_align = 'snakemake -s {}/Snakefiles/{}/pre_align.snake --cores {} -pT -d {} --configfile {} {}'.format(
            scripts_dir,
            samples_yaml['GLOBAL']['data_type'],
            yaml_data['CORES'],
            args.folder_path,
            args.config_file_path,
            complementory_args)
        star_align = 'snakemake -s {}/Snakefiles/{}/star_align.snake --cores {} -pT -d {} --configfile {} {}'.format(
            scripts_dir,
            samples_yaml['GLOBAL']['data_type'],
            yaml_data['CORES'],
            args.folder_path,
            args.config_file_path,
            complementory_args)
        star_summary = 'Rscript {}/Rscripts/STAR_log_plot.R {}'.format(
            package_dir,
            args.folder_path)
        post_align = 'snakemake -s {0}/Snakefiles/{1}/post_align.snake --cores {2} -pT -d {3} --configfile {4} {5}'.format(
            scripts_dir,
            samples_yaml['GLOBAL']['data_type'],
            yaml_data['CORES'],
            args.folder_path,
            args.config_file_path,
            complementory_args)
        print('Running pre-processing')
        shell(pre_align)
        print('Running Alignement')
        try:
            shell(star_align)
        except:
            pass
        print('Plotting STAR logs')
        shell(star_summary)
        print('Running post-alignement')
        shell(post_align)
        if(samples_yaml['GLOBAL']['data_type'] == 'bulk'):
            merge_expression = 'python3 {}/Python/merge_results.py {}'.format(
                package_dir,
                args.folder_path)
            shell(merge_expression)
    if("test" in args.mode):
        test = 'snakemake -s {0}/Snakefiles/{1}/star_align.snake --cores {2} -pT -d {3} --configfile {4} --cluster-config {3}cluster.yaml {5}'.format(
            scripts_dir,
            samples_yaml['GLOBAL']['data_type'],
            yaml_data['CORES'],
            args.folder_path,
            args.config_file_path,
            complementory_args)
        print(test)
        shell(test)
    if("generate-plots" in args.mode):
        print('Mode is generate-plots')
        if(samples_yaml['GLOBAL']['data_type'] == 'singleCell'):
            knee_plot = 'Rscript {}/Rscripts/{}/knee_plot.R {}'.format(
                package_dir,
                samples_yaml['GLOBAL']['data_type'],
                args.folder_path)
            base_summary = 'Rscript {}/Rscripts/{}/rna_metrics.R {}'.format(
                package_dir,
                samples_yaml['GLOBAL']['data_type'],
                args.folder_path)
            print('Plotting knee plots')
            shell(knee_plot)
            print('Plotting base stats')
            shell(base_summary)
        multiqc = 'multiqc -o {0} {0}/logs {0}/summary --force'.format(args.folder_path)
        print('Generating multiqc report')
        shell(multiqc)
    if("species-plot" in args.mode):
        if(len(samples_yaml['SPECIES']) == 2):
            print('Mode is species-plots')
            extract_species = ('snakemake -s {}/Snakefiles/singleCell/extract_species.snake --cores {} -pT -d {} --configfile {} {}'.format(
                scripts_dir,
                yaml_data['CORES'],
                args.folder_path,
                args.config_file_path,
                complementory_args))
            print('Extracting species')
            shell(extract_species)
            species_plot = 'Rscript {}/Rscripts/{}/species_plot.R {}'.format(
                package_dir,
                samples_yaml['GLOBAL']['data_type'],
                args.folder_path)
            print('Plotting species plots')
            shell(species_plot)
        else:
            sys.exit('You cannot run this with a number of species different than 2.\nPlease change the config file')
    if("extract-expression" in args.mode):
        if(len(samples_yaml['SPECIES']) == 2):
            print('Mode is generate-extract-expression')
            extract_expression = 'snakemake -s {}/Snakefiles/singleCell/extract_expression.snake --cores {} -pT -d {} --configfile {} {}'.format(
                scripts_dir,
                yaml_data['CORES'],
                args.folder_path,
                args.config_file_path,
                complementory_args)
            print('Extracting expression')
            shell(extract_expression)
        if(len(samples_yaml['SPECIES']) == 1):
            extract_expression_single = 'snakemake -s {}/Snakefiles/singleCell/extract_expression_single.snake --cores {} -pT -d {} --configfile {} {}'.format(
                scripts_dir,
                yaml_data['CORES'],
                args.folder_path,
                args.config_file_path,
                complementory_args)
            print('Extracting expression')
            shell(extract_expression_single)
            merge_expression = 'Rscript {}/Rscripts/{}/merge_counts_single.R {}'.format(
                package_dir,
                samples_yaml['GLOBAL']['data_type'],
                args.folder_path)
            shell(merge_expression)
    print('Pipeline finished')
if __name__ == "__main__":
    sys.exit(main())
