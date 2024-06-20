from snm3Cmap import __version__

import argparse
import inspect
import subprocess
import sys
import logging
import os


log = logging.getLogger()

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
    
    parser_req.add_argument('--config', type=str, default=None, required=True,
                        help='Path to config YML file')

    parser_opt = parser.add_argument_group("optional arguments")
    
    parser_opt.add_argument('--jobs', type=int, default=2,
                        help='If set, Snakemake is run with this many concurrent processes')
    
    parser_opt.add_argument('--nolock', action="store_true", 
                            help='If set, Snakemake is run with --nolock (working directory will not be locked)')
    
    parser_opt.add_argument('--rerun-incomplete', action="store_true",
                            help="""If set, Snakemake is run with --rerun-incomplete 
                                    (re-run all jobs the output of which is recognized as incomplete)""")

def prepare_mapping_register_subparser(subparser):
    parser = subparser.add_parser('prepare-mapping',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Setup mapping")

    # Required arguments
    parser_req = parser.add_argument_group("required arguments")

    parser_req.add_argument('--config', type=str, default=None, required=True,
                        help='Path to config YML file')

    parser_opt = parser.add_argument_group("optional arguments")
    
    parser_opt.add_argument('--jobs', type=int, default=2,
                        help='If set, Snakemake is run with this many concurrent processes')
    
    parser_opt.add_argument('--nolock', action="store_true", 
                            help='If set, Snakemake is run with --nolock (working directory will not be locked)')
    
    parser_opt.add_argument('--rerun-incomplete', action="store_true",
                            help="""If set, Snakemake is run with --rerun-incomplete 
                                    (re-run all jobs the output of which is recognized as incomplete)""")

def contamination_filter_register_subparser(subparser):
    parser = subparser.add_parser('contamination-filter',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Filter out reads with high CH methylation")

    # Required arguments
    parser_req = parser.add_argument_group("required arguments")
    
    parser_req.add_argument('--bam', type=str, default=None, required=True,
                            help='Path to input bam')
                            
    parser_req.add_argument('--min-mapq', type=int, default=30, required=True,
                            help="MAPQ threshold for considering a read's CH methylation")
                            
    parser_req.add_argument('--out-prefix', type=str, default=None, required=True,
                            help='Path including name prefix for output filtered bam')

def call_contacts_register_subparser(subparser):
    parser = subparser.add_parser('call-contacts',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Call contacts from BAM file and trim split alignments to remove within-mate overlaps")

    # Required arguments
    parser_req = parser.add_argument_group("required arguments")
    
    parser_req.add_argument('--bam', type=str, default=None, required=True,
                            help='Path to input bam')
                            
    parser_req.add_argument('--out-prefix', type=str, default=None, required=True,
                            help='Path including name prefix for output files')

    parser_req.add_argument('--chrom-sizes', type=str, default=None, required=True,
                            help='Path to chromosome sizes file')

    parser_req.add_argument('--restriction-sites', type=str, nargs="+", default=[], required=True,
                            help='Paths to restriction sites files')

    parser_opt = parser.add_argument_group("optional arguments")
    
    parser_opt.add_argument('--min-mapq', type=int, default=30,
                        help='Minimum MAPQ to consider alignment (Pairtools parameter)')

    parser_opt.add_argument('--max-molecule-size', type=int, default=750,
                        help="""The maximal size of a Hi-C molecule; used to rescue single ligations 
                                (from molecules with three alignments) and to rescue complex ligations.
                                Used for walks-policy mask, not walks-policy all (Pairtools parameter)""")
    
    parser_opt.add_argument('--max-inter-align-gap', type=int, default=20,
                      help="""Read segments that are not covered by any alignment and longer than the 
                              specified value are treated as “null” alignments. These null alignments 
                              convert otherwise linear alignments into walks, and affect how they get reported 
                              as a Hi-C pair (Pairtools parameter)""")

    parser_opt.add_argument('--trim-reporting', type=str, default="minimal", choices=['minimal', 'full'],
                            help="""If set, output BAM files have 4 extra tags for split alignments. ZU and ZD report the 
                                    number of basepairs trimmed off of the 5' and 3' end of the alignment, respectively. 
                                    ZL reports if an alignment was determined to have a cut site at the 5' end (U), 3' 
                                    end (D), both (B), or neither (N)""")
    
    parser_opt.add_argument('--min-intra-dist', type=int, default=1000,
                            help='Minimum distance for intrachromosomal contacts')

    parser_opt.add_argument('--read-type', type=str, default="bisulfite", choices=['bisulfite', 'wgs'],
                            help='Indicates that reads were bisulfite converted or not bisulfite converted (wgs)')

    parser_opt.add_argument('--max-cut-site-split-algn-dist', type=int, default = 20,
                            help="""Max allowed distance (bp) from nearest cut site to split alignment to be considered ligation event""")

    parser_opt.add_argument('--max-cut-site-whole-algn-dist', type=int, default = 20,
                            help="""Max allowed distance (bp) from nearest cut site to whole alignment to be considered ligation event""")

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

    parser_opt = parser.add_argument_group("optional arguments")

    parser_opt.add_argument('--min-mapq', type=int, default=30,
                            help='Minimum MAPQ to consider alignment')

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

    parser_opt = parser.add_argument_group("optional arguments")

    parser_opt.add_argument('--num-upstr-bases', type=int, default=0,
                            help='Number of upstream bases for context')

    parser_opt.add_argument('--num-downstr-bases', type=int, default=2,
                            help='Number of downstream bases for context')
    
    parser_opt.add_argument('--min-mapq', type=int, default=30,
                            help='Minimum MAPQ score for including aligned reads')

    parser_opt.add_argument('--min-base-quality', type=int, default=20,
                            help='Minimum base quality for including aligned nucleotides')
    
    parser_opt.add_argument('--compress-level', type=int, default=5,
                            help='Compression level')
    
    parser_opt.add_argument('--save-count-df', action="store_true",
                            help='If set, save context count summary file')

def pairtools_stats_register_subparser(subparser):
    parser = subparser.add_parser('pairtools-stats',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="""
                                        Compute stats for pairs files for contacts and non-ligation artefacts 
                                        """
                                 )
    # Required arguments
    parser_req = parser.add_argument_group("required arguments")

    parser_req.add_argument('--out-prefix', type=str, default=None, required=True,
                        help='Path including name prefix for output stats file')

    parser_req.add_argument('--contacts', type=str, default=None, required=True,
                            help='Path to contacts pairs file')
    
    parser_req.add_argument('--artefacts', type=str, default=None, required=True,
                            help='Path to non-ligation artefacts pairs file')

    parser_req.add_argument('--contacts-stats', type=str, default=None, required=True,
                            help='Path to pairtools dedup stats for contacts pairs file')

    parser_req.add_argument('--artefacts-stats', type=str, default=None, required=True,
                            help='Path to pairtools dedup stats for non-ligation-artefacts pairs file')
    
    parser_req.add_argument('--filterbycov-stats', type=str, default=None, required=True,
                            help='Path to pairtools filterbycov stats for contacts pairs file')
                            
    
def aggregate_qc_stats_register_subparser(subparser):
    parser = subparser.add_parser('aggregate-qc-stats',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="""
                                        Aggregate QC stats generated during trimming, filtering, mapping, and contact/methylation calling, 
                                        as well as generate a .short file for contacts.
                                        """
                                 )
    # Required arguments
    parser_req = parser.add_argument_group("required arguments")
    
    parser_req.add_argument('--cell', type=str, default=None, required=True,
                            help='Name of cell')
                            
    parser_req.add_argument('--out-prefix', type=str, default=None, required=True,
                            help='Path including name prefix for output stats file')

    parser_opt = parser.add_argument_group("optional arguments")

    parser_opt.add_argument('--min-mapq', type=int, default=30,
                            help='Minimum MAPQ score for including aligned reads')

    parser_opt.add_argument('--min-base-quality', type=int, default=20,
                            help='Minimum base quality for including aligned nucleotides')
    

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
        from .mapping import PrepareMapping as func
    elif cur_command in ['contamination-filter']:
        from .mapping import ContaminationFilter as func
    elif cur_command in ['call-contacts']:
        from .mapping import ContactGenerator as func
    elif cur_command in ['pairtools-stats']:
        from .mapping import pairtools_stats as func
    elif cur_command in ['mask-overlaps']:
        from .mapping import OverlapMask as func
    elif cur_command in ['bam-to-allc']:
        from .mapping import bam_to_allc as func
    elif cur_command in ['aggregate-qc-stats']:
        from .mapping import aggregate_qc_stats as func
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