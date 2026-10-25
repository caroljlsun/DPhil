import argparse
from loguru import logger
import os
import tempfile
from typing import Tuple, Literal
import pandas as pd
from capcruncher.api.storage import (
    CoolerBinner,
    create_cooler_cc,
    merge_coolers,
)
import cooler
import warnings

from modded_interactions_store import bins, fragments

# ... function definitions ...

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process genomic data.')
    subparsers = parser.add_subparsers(dest='command')

    # Subparser for the 'fragments' command
    parser_fragments = subparsers.add_parser('fragments', help='Store restriction fragment interaction combinations.')
    parser_fragments.add_argument('counts', help='Path to restriction fragment interactions counts .tsv file.')
    parser_fragments.add_argument('fragment_map', help='Path to restriction fragment .bed file.')
    parser_fragments.add_argument('output', help='Output file path for cooler hdf5 file.')
    parser_fragments.add_argument('viewpoint_path', help='Path to viewpoints bed file.')
    parser_fragments.add_argument('--viewpoint_name', default='', help='Name of viewpoint.')
    parser_fragments.add_argument('--genome', default='', help='Name of genome used for alignment.')
    parser_fragments.add_argument('--suffix', default='', help='Suffix to append to filename.')

    # Subparser for the 'bins' command
    parser_bins = subparsers.add_parser('bins', help='Process data in a cooler file and bin it.')
    parser_bins.add_argument('cooler_path', help='Path to cooler file.')
    parser_bins.add_argument('output', help='Output file path for cooler hdf5 file.')
    parser_bins.add_argument('binsizes', nargs='+', type=int, help='List of bin sizes.')
    parser_bins.add_argument('--normalise', action='store_true', help='Normalise the data.')
    parser_bins.add_argument('--scale_factor', type=int, default=1e6, help='Scale factor for normalisation.')
    parser_bins.add_argument('--overlap_fraction', type=float, default=1e-9, help='Overlap fraction for binning.')
    parser_bins.add_argument('--conversion_tables', help='Path to conversion tables file.')
    parser_bins.add_argument('--n_cores', type=int, default=1, help='Number of cores to use.')
    parser_bins.add_argument('--assay', choices=["capture", "tri", "tiled"], default="capture", help='Assay type.')

    args = parser.parse_args()

    if args.command == 'fragments':
        fragments(args.counts, args.fragment_map, args.output, args.viewpoint_path, args.viewpoint_name, args.genome, args.suffix)
    elif args.command == 'bins':
        bins(args.cooler_path, args.output, tuple(args.binsizes), args.normalise, args.scale_factor, args.overlap_fraction, args.conversion_tables, args.n_cores, args.assay)



# run fragments
#python modded_interactions_store.py fragments counts.tsv fragment_map.bed output.hdf5 viewpoints.bed --viewpoint_name viewpoint1 --genome hg19 --suffix _suffix

# run bins
#python modded_interactions_store.py bins cooler_file.hdf5 output.hdf5 10000 20000 --normalise --scale_factor 1000000 --overlap_fraction 0.000000001 --conversion_tables conversion_tables.txt --n_cores 4 --assay capture
#python test_int_store.py bins /nfs_home/users/czjs/pilot_capture_analysis/capcruncher_output/interim/pileups/CopyOfcounts_by_restriction_fragment/DAAP6preadipocytes_P3.hdf5 DAAP6preadipocytes_P3_binned_10000.hdf5 10000 20000 --normalise --scale_factor 1000000 --overlap_fraction 0.000000001 --assay capture  