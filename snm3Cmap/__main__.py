from snm3Cmap import __version__

import argparse
import inspect
import subprocess
import sys
import logging
import os


log = logging.getLogger()

#__version__ = 1.0

DESCRIPTION = """
Pipeline for mapping snm3Cseq data

"""

EPILOG = ''

class NiceFormatter(logging.Formatter):
    """
    From Cutadapt https://github.com/marcelm/cutadapt
    Do not prefix "INFO:" to info-level log messages (but do it for all other
    levels).
    Based on http://stackoverflow.com/a/9218261/715090 .
    """

    def format(self, record):
        if record.levelno != logging.INFO:
            record.msg = '{}: {}'.format(record.levelname, record.msg)
        return super().format(record)

def setup_logging(stdout=False, quiet=False, debug=False):
    """
    From Cutadapt https://github.com/marcelm/cutadapt
    Attach handler to the global logger object
    """
    # Due to backwards compatibility, logging output is sent to standard output
    # instead of standard error if the -o option is used.
    stream_handler = logging.StreamHandler(sys.stdout if stdout else sys.stderr)
    stream_handler.setFormatter(NiceFormatter())
    # debug overrides quiet
    if debug:
        level = logging.DEBUG
    elif quiet:
        level = logging.ERROR
    else:
        level = logging.INFO
    stream_handler.setLevel(level)
    log.setLevel(level)
    log.addHandler(stream_handler)

def prepare_demultiplex_register_subparser(subparser):
    parser = subparser.add_parser('prepare-demultiplex',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Setup demultiplexing")
    parser_req = parser.add_argument_group("required arguments")
    parser_req.add_argument('--fastq-directory', type=str, default=None, required=True,
                            help='Directory where raw fastq files are stored ending with /')
    parser_req.add_argument('--plate-info', type=str, default=None, required=True,
                            help='Path to plate ids in text file')
    parser_req.add_argument('--output-directory', type=str, default=None, required=True,
                        help='Directory for output fastq files')
    parser_req.add_argument('--barcodes', type=str, default=None, required=True,
                        help='Path to barcodes fasta')

    parser_req = parser.add_argument_group("optional arguments")
    parser_req.add_argument('--jobs', type=int, default=2, required=False,
                        help='Number of concurrent jobs to run (for Snakemake)')

def demultiplex_register_subparser(subparser):
    parser = subparser.add_parser('demultiplex',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Demultiplex a well plate fastq")

    # Required arguments
    parser_req = parser.add_argument_group("required arguments")
    
    parser_req.add_argument('--plate', type=str, default=None, required=True,
                            help='Plate identifier')
    
    parser_req.add_argument('--R1', type=str, default=None, required=True,
                            help='Path to R1 fastq file to be demultiplexed')

    parser_req.add_argument('--R2', type=str, default=None, required=True,
                            help='Path to R2 fastq file to be demultiplexed')
    
    parser_req.add_argument('--barcodes', type=str, default=None, required=True,
                            help='Path to barcodes fasta')
    
    parser_req.add_argument('--out-dir', type=str, default=None, required=True,
                            help='Directory where demultiplexed fastq will be written')

    parser_req = parser.add_argument_group("optional arguments")
    parser_req.add_argument('--threads', type=int, default=8, required=False,
                            help='Number of threads to be used for writing fastq files')


def prepare_mapping_register_subparser(subparser):
    parser = subparser.add_parser('prepare-mapping',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Setup mapping")

    # Required arguments
    parser_req = parser.add_argument_group("required arguments")
    
    parser_req.add_argument('--plate-info', type=str, default=None, required=True,
                            help='Path to plate ids in text file')
                            
    parser_req.add_argument('--demultiplex-directory', type=str, default=None, required=True,
                            help='Path to demultiplexed plates directory')
                            
    parser_req.add_argument('--barcodes', type=str, default=None, required=True,
                            help='Path to barcodes fasta')

    parser_req.add_argument('--reference-genome', type=str, default=None, required=True,
                            help='Path to indexed reference genome')

    parser_req = parser.add_argument_group("optional arguments")
    parser_req.add_argument('--jobs', type=int, default=2, required=False,
                        help='Number of concurrent jobs to run (for Snakemake)')

def contamination_filter_register_subparser(subparser):
    parser = subparser.add_parser('contamination-filter',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Filter out reads with high CH methylation")

    # Required arguments
    parser_req = parser.add_argument_group("required arguments")
    
    parser_req.add_argument('--bam', type=str, default=None, required=True,
                            help='Path to input bam')
                            
    parser_req.add_argument('--mapq-threshold', type=int, default=30, required=True,
                            help="MAPQ threshold for considering a read's CH methylation")
                            
    parser_req.add_argument('--out-prefix', type=str, default=None, required=True,
                            help='Path including name prefix for output filtered bam')

def call_contacts_register_subparser(subparser):
    parser = subparser.add_parser('call-contacts',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Call contacts from BAM file and trim chimeric alignments to remove same-mate overlaps")

    # Required arguments
    parser_req = parser.add_argument_group("required arguments")
    
    parser_req.add_argument('--bam', type=str, default=None, required=True,
                            help='Path to input bam')
                            
    parser_req.add_argument('--out-prefix', type=str, default=None, required=True,
                            help='Path including name prefix for output files')


def mask_overlaps_register_subparser(subparser):
    parser = subparser.add_parser('mask-overlaps',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Masks overlapping bases between mates of the same read pair")

    # Required arguments
    parser_req = parser.add_argument_group("required arguments")
    
    parser_req.add_argument('--bam', type=str, default=None, required=True,
                            help='Path to input bam')
                            
    parser_req.add_argument('--out-prefix', type=str, default=None, required=True,
                            help='Path including name prefix for output bam file')

def bam_to_allc_register_subparser(subparser):
    parser = subparser.add_parser('bam-to-allc',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Convert a biscuit-derived BAM file to ALLC format")

    # Required arguments
    parser_req = parser.add_argument_group("required arguments")
    
    parser_req.add_argument('--bam-path', type=str, default=None, required=True,
                            help='Path to input bam file')

    parser_req.add_argument('--reference-fasta', type=str, default=None, required=True,
                            help='Path to reference fasta file')
    
    parser_req.add_argument('--output-path', type=str, default=None, required=True,
                            help='Path to output ALLC')

    # Optional arguments
    parser_req = parser.add_argument_group("optional arguments")

    parser_req.add_argument('--num-upstr-bases', type=int, default=0, required=False,
                            help='Number of upstream bases for context')

    parser_req.add_argument('--num-downstr-bases', type=int, default=2, required=False,
                            help='Number of downstream bases for context')
    
    parser_req.add_argument('--min-mapq', type=int, default=30, required=False,
                            help='Minimum MAPQ score for including aligned reads')

    parser_req.add_argument('--min-base-quality', type=int, default=20, required=False,
                            help='Minimum base quality for including aligned nucleotides')
    
    parser_req.add_argument('--compress-level', type=int, default=5, required=False,
                            help='Compression level')
    
    parser_req.add_argument('--save-count-df', action="store_true", required=False,
                            help='If set, save context count summary file')

def aggregate_qc_stats_register_subparser(subparser):
    parser = subparser.add_parser('aggregate-qc-stats',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Aggregate QC stats generated during trimming, filtering, mapping, and contact/methylation calling")

    # Required arguments
    parser_req = parser.add_argument_group("required arguments")
    
    parser_req.add_argument('--cell', type=str, default=None, required=True,
                            help='Name of cell')
                            
    parser_req.add_argument('--out-prefix', type=str, default=None, required=True,
                            help='Path including name prefix for output stats file')
    

def main():
    parser = argparse.ArgumentParser(description=DESCRIPTION,
                                     epilog=EPILOG,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     )
    subparsers = parser.add_subparsers(
        title="functions",
        dest="command",
        metavar="",
        required=True
    )

    # add subparsers
    current_module = sys.modules[__name__]
    # get all functions in parser
    for name, register_subparser_func in inspect.getmembers(current_module, inspect.isfunction):
        if 'register_subparser' in name:
            register_subparser_func(subparsers)

    # initiate
    args = None
    if len(sys.argv) > 1:
        # print out version
        if sys.argv[1] in ['-v', '--version']:
            print(__version__)
            exit()
        else:
            args = parser.parse_args()
    else:
        # print out help
        parser.parse_args(["-h"])
        exit()

    # set up logging
    if not logging.root.handlers:
        setup_logging(stdout=True,
                      quiet=False)
    # execute command
    args_vars = vars(args)
    for k, v in args_vars.items():
        log.info(f'{k}\t{v}')

    cur_command = args_vars.pop('command').lower().replace('_', '-')
    # Do real import here:
    if cur_command in ['prepare-demultiplex']:
        from .demultiplex import prepare_demultiplex as func
    elif cur_command in ['demultiplex']:
        from .demultiplex import demultiplexer as func
    elif cur_command in ['prepare-mapping']:
        from .mapping import prepare_mapping as func
    elif cur_command in ['contamination-filter']:
        from .mapping import ContaminationFilter as func
    elif cur_command in ['call-contacts']:
        from .mapping import ContactGenerator as func
    elif cur_command in ['mask-overlaps']:
        from .mapping import OverlapMask as func
    elif cur_command in ['bam-to-allc']:
        from .mapping import bam_to_allc as func
    elif cur_command in ['aggregate-qc-stats']:
        from .mapping import parse_stats as func
    else:
        log.debug(f'{cur_command} is not an valid sub-command')
        parser.parse_args(["-h"])
        return

    # run the command
    log.info(f"Executing {cur_command}...")
    func(**args_vars)
    log.info(f"{cur_command} finished.")
    return
    

if __name__ == '__main__':
    main()