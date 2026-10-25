import argparse
from pileup import CoolerBedGraph

def main():
    parser = argparse.ArgumentParser(description='Generate a bedgraph file from a cooler file.')
    parser.add_argument('uri', help='Path to cooler group in hdf5 file.')
    parser.add_argument('--sparse', action='store_true', help='Only output non-zero bins.')
    parser.add_argument('--only_cis', action='store_true', help='Only output cis interactions.')
    parser.add_argument('--region_to_limit', help='Region to limit.')
    parser.add_argument('--normalisation', choices=['raw', 'n_cis', 'region'], default='raw', help='Normalisation method.')
    parser.add_argument('--scale_factor', type=float, default=1e6, help='Scaling factor for normalisation.')
    parser.add_argument('--region', help='Region for normalisation.')
    parser.add_argument('--output', default='output.bedgraph', help='Output bedgraph file.')

    args = parser.parse_args()

    cbg = CoolerBedGraph(args.uri, args.sparse, args.only_cis, args.region_to_limit)
    bedgraph = cbg.extract_bedgraph(args.normalisation)
    cbg._normalise_bedgraph(bedgraph, args.scale_factor, args.normalisation, args.region)
    bedgraph.to_csv(args.output, sep='\t', header=False, index=False)


if __name__ == "__main__":
    main()

# DAAP6matureadipocytes_M1_binned_10000.hdf5
# python test_bed.py DAAP6matureadipocytes_M1_binned_10000.hdf5:: --output=DAAP6matureadipocytes_M1_binned_10000.bedgraph 