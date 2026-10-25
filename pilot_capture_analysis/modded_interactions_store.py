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



def fragments(
    counts: os.PathLike,
    fragment_map: os.PathLike,
    output: os.PathLike,
    viewpoint_path: os.PathLike,
    viewpoint_name: str = "",
    genome: str = "",
    suffix: str = "",
):
    """
    Stores restriction fragment interaction combinations at the restriction fragment level.

    Parses reporter restriction fragment interaction counts produced by
    "capcruncher reporters count" and gerates a cooler formatted group in an HDF5 File.
    See `https://cooler.readthedocs.io/en/latest/` for further details.


    \f
    Args:
     counts (os.PathLike): Path to restriction fragment interactions counts .tsv file.
     fragment_map (os.PathLike): Path to restriction fragment .bed file, generated with genome-digest command.
     output (os.PathLike): Output file path for cooler hdf5 file.
     viewpoint_name (str): Name of viewpoint.
     viewpoint_path (os.PathLike): Path to viewpoints bed file.
     genome (str, optional): Name of genome used for alignment e.g. hg19. Defaults to "".
     suffix (str, optional): Suffix to append to filename. Defaults to "".
    """
    # Load restriction fragments
    df_restriction_fragment_map = pd.read_csv(
        fragment_map,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name"],
    )

    # Load counts
    if counts.endswith(".hdf5"):

        with pd.HDFStore(counts) as store:

            if not viewpoint_name:
                viewpoints = {k.split("/")[1] for k in store.keys()}
            else:
                viewpoints = {
                    viewpoint_name,
                }

            for viewpoint in viewpoints:
                df_counts = store[viewpoint]

                create_cooler_cc(
                    output,
                    bins=df_restriction_fragment_map,
                    pixels=df_counts,
                    viewpoint_name=viewpoint,
                    viewpoint_path=viewpoint_path,
                    assembly=genome,
                    suffix=suffix,
                )

    else:
        df_counts = pd.read_csv(counts, sep="\t")
        # Create cooler file at restriction fragment resolution
        create_cooler_cc(
            output,
            bins=df_restriction_fragment_map,
            pixels=df_counts,
            viewpoint_name=viewpoint_name,
            viewpoint_path=viewpoint_path,
            assembly=genome,
            suffix=suffix,
        )



def _bin_cooler(clr_in: os.PathLike, clr_out: os.PathLike, binsize: int, **kwargs):
    clr_binner = CoolerBinner(
        cooler_group=clr_in,
        binsize=binsize,
        **kwargs,
    )
    clr_binner.to_cooler(clr_out)
    return clr_out

def bins(
    cooler_path: os.PathLike,
    output: os.PathLike,
    binsizes: Tuple = None,
    normalise: bool = True,
    scale_factor: int = 1e6,
    overlap_fraction: float = 1e-9,
    conversion_tables: os.PathLike = None,
    n_cores: int = 1,
    assay: Literal["capture", "tri", "tiled"] = "capture",
    **kwargs,
):
    clr_groups = cooler.api.list_coolers(cooler_path)

    assert clr_groups, "No cooler groups found in file"
    assert binsizes, "No binsizes provided"

    clr_tempfiles = []

    for binsize in binsizes:
        for clr_group in clr_groups:

            logger.info(f"Processing {clr_group}")
            clr_in = cooler.Cooler(f"{cooler_path}::{clr_group}")
            clr_out = tempfile.NamedTemporaryFile().name

            default_kwargs = dict(
                method="midpoint",
                minimum_overlap=0.51,
                n_cis_interaction_correction=True,
                n_rf_per_bin_correction=True,
                scale_factor=1_000_000,
                assay=assay,
            )

            clr_tempfiles.append(
                _bin_cooler(clr_in, clr_out, binsize, **default_kwargs)
            )

    merge_coolers(clr_tempfiles, output)