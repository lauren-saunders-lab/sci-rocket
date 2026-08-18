__version__ = "1.0"

# Import modules.
import argparse
import gzip
import logging
import os
import pickle
import re
import sys
from collections import defaultdict

import pandas as pd
import pysam

from frozendict import frozendict
from preprocess import retrieve_barcodes
from sciClasses import sciRecord


def init_barcode_dict(barcodes: pd.DataFrame, samples: pd.DataFrame, experiment_name: str):
    """
    Generate the barcode dictionaries.
    These are used for fast lookup of the barcodes.

    Parameters:
    \»
        barcodes (pd.DataFrame): Barcode sheet.
        samples (pd.DataFrame): Sample sheet.
        experiment_name (str): Experiment name.

    Returns:
        dict_barcodes (dict): Dictionary of barcodes (frozendict).
    """

    # Determine length of sequences in the ligation barcodes.
    barcodes.loc[:, "length"] = barcodes["sequence"].str.len()

    # Initialize the barcode dictionary.
    dict_barcodes = {}

    # Get the samples used in this sequencing run.
    samples_exp = samples.query("experiment_name == @experiment_name")

    # Generate k;v of ligation barcodes.
    dict_barcodes["ligation"] = {}
    dict_barcodes["ligation"]["ligation_10nt"] = dict(zip(barcodes.query("type == 'ligation' & length == 10")["sequence"], barcodes.query("type == 'ligation' & length == 10")["barcode"]))
    dict_barcodes["ligation"]["ligation_9nt"] = dict(zip(barcodes.query("type == 'ligation' & length == 9")["sequence"], barcodes.query("type == 'ligation' & length == 9")["barcode"]))

    # Get the possible list of unique RT barcodes contained in the sequencing run.
    rt_barcodes_sequencing = retrieve_barcodes(None, samples_exp["rt"].unique(), barcodes, "rt")
    dict_barcodes["rt"] = dict(zip(barcodes.query("type == 'rt' & barcode in @rt_barcodes_sequencing")["sequence"], barcodes.query("type == 'rt' & barcode in @rt_barcodes_sequencing")["barcode"]))

    # Retrieve the p5 barcodes used in this sequencing run.
    indexes_p5 = retrieve_barcodes(None, samples_exp["p5"].unique(), barcodes, "p5")
    barcodes_p5 = dict(zip(barcodes.query("type == 'p5' & barcode in @indexes_p5")["sequence"], barcodes.query("type == 'p5' & barcode in @indexes_p5")["barcode"]))

    # Reverse complement the p5 barcodes.
    dict_barcodes["p5"] = {k[::-1].translate(str.maketrans("ATCG", "TAGC")): v for k, v in barcodes_p5.items()}

    # Retrieve the p7 barcodes used in this sequencing run.
    indexes_p7 = retrieve_barcodes(None, samples_exp["p7"].unique(), barcodes, "p7")
    dict_barcodes["p7"] = dict(zip(barcodes.query("type == 'p7' & barcode in @indexes_p7")["sequence"], barcodes.query("type == 'p7' & barcode in @indexes_p7")["barcode"]))

    # Convert the barcode dictionaries to frozendict.
    dict_barcodes["p5"] = frozendict(dict_barcodes["p5"])
    dict_barcodes["p7"] = frozendict(dict_barcodes["p7"])
    dict_barcodes["ligation"]["ligation_10nt"] = frozendict(dict_barcodes["ligation"]["ligation_10nt"])
    dict_barcodes["ligation"]["ligation_9nt"] = frozendict(dict_barcodes["ligation"]["ligation_9nt"])
    dict_barcodes["rt"] = frozendict(dict_barcodes["rt"])

    # Return the barcode dictionary.
    return dict_barcodes

def generate_sample_dict(samples: pd.DataFrame, barcodes: pd.DataFrame):
    """
    Generates a dictionary of the samples and their respective barcodes.
    The dictionary is used for fast lookup of the barcodes.

    Parameters:
        samples (pd.DataFrame): Sample sheet.
        barcodes (pd.DataFrame): Barcode sheet.

    Returns:
        dict_samples (dict): Dictionary of samples and barcodes (frozendict).
    """

    # Expand p5 and p7 barcodes (A01:G12) to individual barcodes (A01, A02, A03, ..., G12).
    samples = samples.to_dict(orient="index")    
    sample_dict = {}

    for sample in samples:
        sample_p5 = retrieve_barcodes(None, list([samples[sample]["p5"]]), barcodes, "p5")
        sample_p7 = retrieve_barcodes(None, list([samples[sample]["p7"]]), barcodes, "p7")
        sample_rt = retrieve_barcodes(None, list([samples[sample]["rt"]]), barcodes, "rt")

        for p5 in sample_p5:
            for p7 in sample_p7:
                for rt in sample_rt:
                    sample_dict["{}_{}_{}".format(p5, p7, rt)] = samples[sample]["sample_name"]

    return frozendict(sample_dict)


def init_qc(experiment_name: str, dict_barcodes: dict, samples: pd.DataFrame, dict_hashing: str = None):
    """
    Generates the QC dictionary.

    Parameters:
        experiment_name (str): Experiment name.
        dict_barcodes (dict): Dictionary of barcodes (frozendict).
        samples (pd.DataFrame): Sample sheet.
        dict_hashing (dict): Dictionary of hashing sheets.

    Returns:
        qc (dict): Dictionary of QC metrics.
    """

    qc = {}
    qc["version"] = __version__
    qc["experiment_name"] = experiment_name
    qc["n_pairs"] = 0  # Total number of initial read-pairs.
    qc["n_pairs_success"] = 0  # Total number of read-pairs with correct RT, p5, p7 and ligation barcodes.
    qc["n_pairs_failure"] = 0  # Total number of discarded read-pairs due to various reason.
    qc["max_r2_read_length"] = 0  # Maximum observed R2 read length.

    qc["n_corrected_p5"] = 0  # Total number of read-pairs with 1bp mismatch in p5.
    qc["n_corrected_p7"] = 0  # Total number of read-pairs with 1bp mismatch in p7.
    qc["n_corrected_ligation"] = 0  # Total number of read-pairs with 1bp mismatch in ligation.
    qc["n_corrected_rt"] = 0  # Total number of read-pairs with 1bp mismatch in RT.
    qc["n_corrected_hashing"] = 0  # Total number of read-pairs with 1bp mismatch in hashing barcode.

    # For all (succesfull) reads, keep track of the number of times each index/barcode.
    qc["p5_index_counts"] = {k: 0 for k in dict_barcodes["p5"].values()}
    qc["p7_index_counts"] = {k: 0 for k in dict_barcodes["p7"].values()}
    qc["rt_barcode_counts"] = {}  # Keep track of the number of times each RT barcode is seen per plate.
    
    # Combine the ligation barcodes (9nt and 10nt) into one dictionary.
    qc["ligation_barcode_counts"] = {**{k: 0 for k in dict_barcodes["ligation"]["ligation_9nt"].values()}, **{k: 0 for k in dict_barcodes["ligation"]["ligation_10nt"].values()}}

    # Keep track of uncorrectable barcode(s) schemes within each R1.
    # Structure: ["p5", "p7", "ligation", "rt"]
    # Example: qc["uncorrectables"][True, False, True, False] to indicate that the p7 and RT barcodes are uncorrectable.
    qc["uncorrectables_sankey"] = {}

    # Add uncorrectables barcode combinations to the dict.
    for p5 in [True, False]:
        for p7 in [True, False]:
            for ligation in [True, False]:
                for rt in [True, False]:
                    qc["uncorrectables_sankey"][p5, p7, ligation, rt] = 0

    qc["n_uncorrectable_p5"] = 0  # Total number of read-pairs with >1bp mismatch in p5.
    qc["n_uncorrectable_p7"] = 0  # Total number of read-pairs with >1bp mismatch in p7.
    qc["n_uncorrectable_ligation"] = 0  # Total number of read-pairs with >1bp mismatch in ligation.
    qc["n_uncorrectable_rt"] = 0  # Total number of read-pairs with >1bp mismatch in RT.

    # Keep track of the sequence of p5, p7, ligation and RT barcodes that are not found in the respective lookup tables.
    # Structure: qc["uncorrectable_p5"]["ATCG"] = 10
    qc["uncorrectable_p5"] = defaultdict(int)
    qc["uncorrectable_p7"] = defaultdict(int)
    qc["uncorrectable_ligation"] = defaultdict(int)
    qc["uncorrectable_rt"] = defaultdict(int)

    # Keep track of the number of succesfull read-pairs per sample.
    qc["sample_success"] = {k: {"n_pairs_success": 0} for k in samples["sample_name"].unique()}

    # Hash-specific qc metrics and regex (per sample).
    if dict_hashing:
        qc["hashing"] = {}
        qc["n_hashing"] = 0  # Total number of read-pairs with hashing barcode.
        for hashing_sample in dict_hashing:
            # Total number of read-pairs with correct(ed) hashing barcode.
            # For all (succesfull) reads, keep track of the number of times each hashing/UMI combination is seen per cellular barcode.
            # The counts contain:
            # - number of times a hashing barcode is seen per cell,
            # - distinct UMI's per cell,
            # - barcode label representation (p7_p5_ligation_rt).
            # qc["hashing"][hashing_sample]["counts"][hashing_name][cellular_sequence] = {"umi": set(), "count": 0, "cell_barcode_label": str}
            qc["hashing"][hashing_sample] = {k: {"n_correct": 0, "n_corrected": 0, "n_correct_upstream": 0, "counts": {}} for k in dict_hashing[hashing_sample]["sheet"].values()}

    # Return the QC dictionary.
    return qc


def update_qc(qc:dict, x:sciRecord):
    """
    Update the QC dictionary.

    Parameters:
        qc (dict): Dictionary of QC metrics.
        x (sciRecord): sciRecord object.

    Returns:
        qc (dict): Updated dictionary of QC metrics.    
    """

    # Update the total number of read-pairs processed.
    qc["n_pairs"] += 1
    qc["max_r2_read_length"] = max(qc["max_r2_read_length"], len(x.read2.sequence))

    # Update the number of read-pairs with correct barcodes.
    # Also count the hashing reads for this.
    if x.sample_name != None:
        # Count as a successful read-pair.
        qc["n_pairs_success"] += 1

        # Keep track of the number of succesfull read-pairs per sample.
        qc["sample_success"][x.sample_name]["n_pairs_success"] += 1

        # Count the number of corrected barcodes.
        qc["n_corrected_p7"] += 1 if x.p7_status == "Corrected" else 0
        qc["n_corrected_p5"] += 1 if x.p5_status == "Corrected" else 0
        qc["n_corrected_ligation"] += 1 if x.ligation_status == "Corrected" else 0
        qc["n_corrected_rt"] += 1 if x.rt_status == "Corrected" else 0
        qc["n_corrected_hashing"] += 1 if x.hashing_status == "Corrected" else 0

        # For all (succesfull) reads, keep track of the number of times each index/barcode.
        qc["p5_index_counts"][x.p5_name] += 1
        qc["p7_index_counts"][x.p7_name] += 1
        qc["ligation_barcode_counts"][x.ligation_name] += 1

        # Count the occurence of the RT barcodes (per plate).
        plate_well, plate_index = x.rt_name.split("-")

        if plate_well not in qc["rt_barcode_counts"]:
            qc["rt_barcode_counts"][plate_well] = {}
            qc["rt_barcode_counts"][plate_well][plate_index] = 1
        else:
            if plate_index not in qc["rt_barcode_counts"][plate_well]:
                qc["rt_barcode_counts"][plate_well][plate_index] = 1
            else:
                qc["rt_barcode_counts"][plate_well][plate_index] += 1

        # Update hashing metrics (if applicable and if found).
        if x.hashing_name:
            qc["n_hashing"] += 1
            cell_barcode_label = f"{x.p7_name}_{x.p5_name}_{x.ligation_name}_{x.rt_name}"
            if x.cellular_sequence not in qc["hashing"][x.sample_name][x.hashing_name]["counts"]:
                qc["hashing"][x.sample_name][x.hashing_name]["counts"][x.cellular_sequence] = {"umi": set([x.umi_sequence]), "count": 1, "cell_barcode_label": cell_barcode_label}
            else:
                qc["hashing"][x.sample_name][x.hashing_name]["counts"][x.cellular_sequence]["umi"].add(x.umi_sequence)
                qc["hashing"][x.sample_name][x.hashing_name]["counts"][x.cellular_sequence]["count"] += 1

            # Update the number of correct/corrected/upstream hashing barcodes.
            qc["hashing"][x.sample_name][x.hashing_name]["n_correct"] += 1 if x.hashing_status == "Correct" else 0
            qc["hashing"][x.sample_name][x.hashing_name]["n_corrected"] += 1 if x.hashing_status == "Corrected" else 0
            qc["hashing"][x.sample_name][x.hashing_name]["n_correct_upstream"] += 1 if x.hashing_status == "Corrected upstream" else 0

    else:
        # Count as a failed read-pair.
        qc["n_pairs_failure"] += 1

        # Keep track which barcode fails.
        sankey_index = [x.p5_name, x.p7_name, x.ligation_name, x.rt_name]
        sankey_index = [True if x != None else False for x in sankey_index]
        qc["uncorrectables_sankey"][tuple(sankey_index)] += 1

        # Count the number of uncorrectable barcodes.
        qc["n_uncorrectable_p7"] += 1 if x.p7_status == None else 0
        qc["n_uncorrectable_p5"] += 1 if x.p5_status == None else 0
        qc["n_uncorrectable_ligation"] += 1 if x.ligation_status == None else 0
        qc["n_uncorrectable_rt"] += 1 if x.rt_status == None else 0

        # Keep track of the uncorrectable barcode sequences.
        qc["uncorrectable_p7"][x.p7_sequence] += 1 if x.p7_status == None else 0
        qc["uncorrectable_p5"][x.p5_sequence] += 1 if x.p5_status == None else 0
        qc["uncorrectable_ligation"][x.ligation_sequence] += 1 if x.ligation_status == None else 0
        qc["uncorrectable_rt"][x.rt_sequence] += 1 if x.rt_status == None else 0

    # Per 1M read-pairs, only keep the top 50 most frequent uncorrectable barcode sequences.
    # This reduces the size of the pickle file and we only display the top 15 in the QC report.
    if qc["n_pairs"] % 1000000 == 0:
        for key in ["p5", "p7", "ligation", "rt"]:
            qc[f"uncorrectable_{key}"] = defaultdict(int, dict(sorted(qc[f"uncorrectable_{key}"].items(), key=lambda kv: kv[1], reverse=True)[:50]))

    # Return the updated QC dictionary.
    return qc

def open_file_handlers(samples: pd.DataFrame, experiment_name: str, path_r1: str, path_r2: str, path_out: str, log: logging.Logger):
    """
    Open file handlers for all experiment and sample-specific output files.

    Parameters:
        samples (pd.DataFrame): Sample sheet.
        experiment_name (str): Experiment name.
        path_r1 (str): Path to R1 fastq file.
        path_r2 (str): Path to R2 fastq file.
        path_out (str): Path to output directory.
        log (logging.Logger): Logger.
    Returns:
        dict_fh (dict): Dictionary of file handlers.
    """
    
    dict_fh = {}
    
    # Input R1 / R2.
    try:
        dict_fh["r1"] = pysam.FastxFile(path_r1)
        dict_fh["r2"] = pysam.FastxFile(path_r2)
    except OSError:
        log.error("Could not find experiment-based R1 and R2 .fq.gz files, please check the paths:\n(R1) %s\n(R2) %s", path_r1, path_r2)
        sys.exit(1)

    # Discarded R1 / R2.
    path_r1_discarded = os.path.join(path_out, experiment_name + "_R1_discarded.fastq.gz")
    path_r2_discarded = os.path.join(path_out, experiment_name + "_R2_discarded.fastq.gz")
    path_log_discarded = os.path.join(path_out, "log_" + experiment_name + "_discarded_reads.tsv.gz")

    try:
        dict_fh["r1_discarded"] = gzip.open(path_r1_discarded, "wt", compresslevel=2)
        dict_fh["r2_discarded"] = gzip.open(path_r2_discarded, "wt", compresslevel=2)
        dict_fh["discarded_log"] = gzip.open(path_log_discarded, "wt", compresslevel=2)

    except OSError:
        log.error(
            "Could not generate experiment-based discard files: Please check the paths:\n(R1 discarded) %s\n(R2 discarded) %s\n(log discarded) %s",
            path_r1_discarded,
            path_r2_discarded,
            path_log_discarded
        )
        sys.exit(1)

    # Sample-specific R1 / R2 output.
    output_sample_names = sorted(samples["sample_name"].unique().tolist())
    for sample in output_sample_names:
        # Initialize the sample dictionary.
        dict_fh[sample] = {}

        # Generate the output paths.
        path_r1_out = os.path.join(path_out, sample + "_R1.fastq.gz")
        path_r2_out = os.path.join(path_out, sample + "_R2.fastq.gz")

        try:
            # Open file handlers for output R1 and R2 files.
            dict_fh[sample]["r1"] = gzip.open(path_r1_out, "wt", compresslevel=2)
            dict_fh[sample]["r2"] = gzip.open(path_r2_out, "wt", compresslevel=2)

        except OSError:
            log.error("Could not generate sample-based demultiplexed files, please check the paths:\n(R1) %s\n(R2) %s", path_r1_out, path_r2_out)
            sys.exit(1)

    # Return the dictionary of file handlers.
    return dict_fh
    

def retrieve_hashing_sheets(samples: pd.DataFrame):
    """
    Import the hashing barcodes for samples that have a hashing sheet attached.

    Parameters:
        samples (pd.DataFrame): Sample sheet.

    Returns:
        hashing (dict): Dictionary of hashing sheets:
            - hashing[sample_name]["sheet"] = frozendict of hashing barcodes.
            - hashing[sample_name]["regex"] = pre-compiled regex of hashing barcodes.
    """

    # Initialize the hashing dictionary.
    hashing = {}

    # For each sample that has a hashing sheet attached, import the sheet and compile the regex for matching upstream hashing barcodes.
    if "hashing" in samples.columns:
        for hashing_sample in samples[~samples["hashing"].isna()].sample_name.unique():
            path_hash = samples.query("sample_name == @hashing_sample")["hashing"].iloc[0]
            x = pd.read_csv(path_hash, sep="\t", header=0)

            # Generate a dictionary with the hashing barcodes as keys and the sample names as values.
            hashing[hashing_sample] = {}
            hashing[hashing_sample]["sheet"] = frozendict(dict(zip(x["barcode"], x["hash_name"])))

            # Generate hashing regex for matching upstream.
            hash_match = "|".join(hashing[hashing_sample]["sheet"].keys())
            hash_regex = re.compile(hash_match)
            hashing[hashing_sample]["regex"] = hash_regex

    return hashing


def sciseq_sample_demultiplexing(
    log: logging.Logger,
    experiment_name: str,
    samples: pd.DataFrame,
    barcodes: pd.DataFrame,
    path_r1: str,
    path_r2: str,
    path_out: str,
    ignore_p7: bool = False,
):
    """
    Performs demultiplexing of the raw fastq files based on the PCR indexes (p5, p7) and RT barcode to produce sample-specific R1 and R2 files.
    The ligation barcode can be either 9nt or 10nt long and this can affect the location of the UMI and RT barcodes.

    Read-pairs are discarded if:
        - The R1 read is empty.
        - One of the barcodes (p5, p7, ligation and/or RT) is not found within R1.

    When a experiment/sample has a hashing sheet attached, the hashing barcode is retrieved from the R2 sequence and added to the read-name of R2.
    Additional metrics are generated for hashing experiments to keep track of the hashing barcodes per cellular barcode.

    Parameters:
        log (logging.Logger): Logger.
        experiment_name (str): Experiment name.
        samples (pd.DataFrame): Sample sheet.
        barcodes (pd.DataFrame): Barcode sheet.
        path_r1 (str): Path to R1 fastq file.
        path_r2 (str): Path to R2 fastq file.
        path_out (str): Path to output directory.
        ignore_p7 (bool): If True, skip matching the sequenced i7 index and assign the single
            p7 well configured for this experiment directly. Requires exactly one p7 well to be
            configured across all samples in the experiment.
    Returns:
        None
    """

    log.info("Starting sample-based demultiplexing of %s:\n(R1) %s\n(R2) %s", experiment_name, path_r1, path_r2)

    # Subset samples used for matching and QC in this experiment.
    samples_exp = samples.query("experiment_name == @experiment_name")

    # Open the IO handlers.
    dict_fh = open_file_handlers(samples, experiment_name, path_r1, path_r2, path_out, log)

    # Generate the barcode dictionaries.
    dict_barcodes = init_barcode_dict(barcodes, samples_exp, experiment_name)

    # If ignoring p7, resolve the single configured p7 well up front (the sequenced i7 index is not used).
    fixed_p7_sequence = None
    fixed_p7_name = None
    if ignore_p7:
        if len(dict_barcodes["p7"]) != 1:
            log.error(
                "settings.ignore_p7 is enabled but experiment '%s' has %d distinct p7 well(s) configured; "
                "ignore_p7 requires exactly one shared p7 well across all samples in the experiment.",
                experiment_name,
                len(dict_barcodes["p7"]),
            )
            sys.exit(1)
        fixed_p7_sequence, fixed_p7_name = next(iter(dict_barcodes["p7"].items()))
        log.info("settings.ignore_p7 is enabled: assigning all reads to p7 well '%s' (%s) without matching the sequenced i7 index.", fixed_p7_name, fixed_p7_sequence)

    # Generate the sample dictionary.
    dict_samples = generate_sample_dict(samples_exp, barcodes)

    # Import hashing information (if applicable)
    dict_hashing = retrieve_hashing_sheets(samples_exp)

    if dict_hashing:
        log.info("Hashing sample(s) detected in this experiment, hashing subroutines are enabled for these samples. After counting, hash-reads are discarded.")

    # Initialize the QC dictionary.
    qc = init_qc(experiment_name, dict_barcodes, samples_exp, dict_hashing)

    # Iterate over the read-pairs and search for the barcodes within R1.
    # If any barcode is not found, try to rescue a respective barcode sequence with 1bp mismatch.
    # If even one barcode is not found (p5, p7, ligation or RT), discard the read-pair and report status.
    for read1, read2 in zip(dict_fh["r1"], dict_fh["r2"]):

        # Create a sciRecord object using the read-pair.
        x = sciRecord(read1, read2)

        # region Retrieve the sci-seq barcodes from R1 -------------------------------------------------------------------------------------------------
        
        # Retrieve the sci-seq barcodes from R1.
        x.determine_p5(dict_barcodes["p5"])
        if ignore_p7:
            x.set_p7_ignored(fixed_p7_sequence, fixed_p7_name)
        else:
            x.determine_p7(dict_barcodes["p7"])
        x.determine_ligation(dict_barcodes["ligation"])
        x.determine_rt(dict_barcodes["rt"])
        x.determine_umi()

        # Determine sample and set cellular barcode.
        x.determine_sample(dict_samples)

        # Determine the hash barcode in R2 (if applicable).
        x.determine_hash(dict_hashing)

        # endregion --------------------------------------------------------------------------------------------------------------------------------
        
        # region Writing output files -----------------------------------------------------------------------------------------------------------

        # If any barcode (or corresponding sample) is not found, discard the read-pair.
        if x.sample_name == None:

            # Write to discarded fq.gz files.
            dict_fh["r1_discarded"].write(str(read1) + "\n")
            dict_fh["r2_discarded"].write(str(read2) + "\n")

            # Write to log (.tsv.gz) to keep track of missing barcodes.
            dict_fh["discarded_log"].write(
                "\t".join(
                    (
                        str(x.read1.name),
                        str(x.p5_name or x.p5_sequence or "?"),
                        str(x.p7_name or x.p7_sequence or "?"),
                        str(x.ligation_name or x.ligation_sequence or "?"),
                        str(x.rt_name or x.rt_sequence or "?"),
                        x.umi_sequence,
                        str(x.sample_name or "?")
                    )
                )
                + "\n"
            )

        # If all barcodes are found (except hashing barcode), write the read-pair to the correct sample file.
        if x.sample_name != None and x.hashing_name == None:
            # Convert the read-pair to the correct format.
            x.convert_R1andR2()

            # Write R1 and R2 to sample-specific files.
            dict_fh[x.sample_name]["r1"].write(str(x.read1) + "\n")
            dict_fh[x.sample_name]["r2"].write(str(x.read2) + "\n")

        # endregion --------------------------------------------------------------------------------------------------------------------------------

        # region QC and Logging --------------------------------------------------------------------------------------------------------------------

        # Update the QC dictionary.
        qc = update_qc(qc, x)

        if qc["n_pairs"] % 1000000 == 0:
            if dict_hashing:
                log.info("Done: %d read-pairs processed (%d discarded, %d hashing reads)", qc["n_pairs"], qc["n_pairs_failure"], qc["n_hashing"])
            else:
                log.info("Done: %d read-pairs processed (%d discarded,)", qc["n_pairs"], qc["n_pairs_failure"])

        # endregion --------------------------------------------------------------------------------------------------------------------------------

    # region Write final files ---------------------------------------------------------------------------------------------------------------------

    # STARSolo requires the barcodes in separate whitelist files.
    with open(os.path.join(path_out, experiment_name + "_whitelist_p5.txt"), "w") as fh:
        for sequence, barcode in dict_barcodes["p5"].items():
            fh.write(sequence + "\n")

    with open(os.path.join(path_out, experiment_name + "_whitelist_p7.txt"), "w") as fh:
        for sequence, barcode in dict_barcodes["p7"].items():
            fh.write(sequence + "\n")

    # Combine the ligation barcodes (9nt and 10nt) into one dictionary.
    with open(os.path.join(path_out, experiment_name + "_whitelist_ligation.txt"), "w") as fh:
        for sequence, barcode in dict_barcodes["ligation"]["ligation_10nt"].items():
            fh.write(sequence + "\n")
        for sequence, barcode in dict_barcodes["ligation"]["ligation_9nt"].items():
            fh.write(sequence + "G" + "\n")

    with open(os.path.join(path_out, experiment_name + "_whitelist_rt.txt"), "w") as fh:
        for sequence, barcode in dict_barcodes["rt"].items():
            fh.write(sequence + "\n")

    # Save QC to pickle file.
    pickle.dump(qc, open(os.path.join(path_out, experiment_name) + "_qc.pickle", "wb"))

    # endregion --------------------------------------------------------------------------------------------------------------------------------

    if dict_hashing:
        log.info("Done: %d read-pairs processed (%d discarded, %d hashing reads)", qc["n_pairs"], qc["n_pairs_failure"], qc["n_hashing"])
    else:
        log.info("Done: %d read-pairs processed (%d discarded,)", qc["n_pairs"], qc["n_pairs_failure"])
    log.info(sciRecord._sciRecord__find_closest_match.cache_info())
        
    # Close the file handlers.
    for fh in dict_fh.values():
        if isinstance(fh, dict):
            for sub_fh in fh.values():
                sub_fh.close()
        else:
            fh.close()


def init_logger():
    """
    Initializes the logger.

    Parameters:
        None

    Returns:
        log (logging.Logger): Logger object.
    """
    log = logging.getLogger(__name__)
    time_format = "%Y-%m-%d %H:%M:%S"
    formatter = logging.Formatter(fmt='%(asctime)s - (sci-rocket) %(levelname)s: %(message)s', datefmt=time_format)

    # Add a stream handler.
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    log.addHandler(stream_handler)

    # Set the verbosity level.
    log.setLevel(logging.INFO)

    return log


def main(arguments):
    description = """
    Performs demultiplexing on R1/R2.fastq(.gz) files generated by sci-RNA-seq v3 protocol, based on:
        - p7, p5, RT and ligation barcodes.
        - Collects hashing metrics (if applicable).
            - Reads used for hashing are removed.
    
    The R1 sequence is modified to a fixed length sequence (48nt) which includes all (corrected) barcodes: p5(10nt), p7(10nt), ligation(10nt), RT(10nt) and UMI (8nt) for downstream processing.
    The read names for R2 are modified to include the barcodes and UMI.
    
    It requires that the p5 and p7 barcodes are present in the read headers of the .fastq files (bcl2fastq):
    @<read name> 1:N:0:ACGGNNGGCC+NTCATGGNGC
                      |----p7---|+|----p5----|: p5 is reverse-complemented.

    The R1 sequence should adhere to the following scheme:
    First 9 or 10nt:  Ligation barcode
    Next 6nt:    Primer
    Next 8nt:    UMI
    Last 10nt:   RT Barcode (sample-specific)

    Anatomy of R1 (ligation of 10nt):
    |ACTTGATTGT| |CAGAGC| |TTTGGTAT| |CCTACCAGTT|
    |-LIGATION-| |Primer| |---UMI--| |----RT----|

    Anatomy of R1 (ligation of 9nt):
    |CTCGTTGAT| |CAGAGC| |TTTGGTAT| |CCTACCAGTT| |T|
    |-LIGATION| |Primer| |---UMI--| |----RT----| |.| <- Extra base.

    Anatomy of R2 (hashing read):
    - Has hashing barcode (10nt) within the first 10nt (0 or 1 hamming distance) or elsewhere in R2 (no hamming).
    - Has poly-A sequence (8xA).
    """

    # Setup argument parser.
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    parser.add_argument("--r1", required=True, type=str, help="(.fq) Input fastq (R1).")
    parser.add_argument("--r2", required=True, type=str, help="(.fq) Input fastq (R2).")
    parser.add_argument("--experiment_name", required=True, type=str, help="(str) Experiment name.")
    parser.add_argument("--samples", required=True, type=str, help="(str) Path to sample-sheet.")
    parser.add_argument("--barcodes", required=True, type=str, help="(str) Path to barcodes file.")
    parser.add_argument("--out", required=True, type=str, help="(str) Path to output directory.")
    parser.add_argument(
        "--ignore_p7",
        action="store_true",
        default=False,
        help="(flag) Skip matching the sequenced i7 index and assign the single p7 well configured for the experiment directly. Requires exactly one p7 well across all samples in the experiment.",
    )
    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="Display help and exit.")
    parser.add_argument("-v", "--version", action="version", version=__version__, help="Display version and exit.")

    # Parse arguments.
    args = parser.parse_args()

    # Initialize logging.
    log = init_logger()

    # Open sample-sheet.
    samples = pd.read_csv(args.samples, sep="\t", dtype=str, comment="#")

    # Open barcode-sheet.
    barcodes = pd.read_csv(args.barcodes, sep="\t", dtype=str)

    # Generate output directory if not exists.
    if not os.path.exists(args.out):
        os.makedirs(args.out)

    # Run the program.
    sciseq_sample_demultiplexing(
        log=log,
        experiment_name=args.experiment_name,
        samples=samples,
        barcodes=barcodes,
        path_r1=args.r1,
        path_r2=args.r2,
        path_out=args.out,
        ignore_p7=args.ignore_p7,
    )


if __name__ == "__main__":
    main(sys.argv[1:])
    sys.exit()

# args = parser.parse_args(
#     [
#         "--r1",
#         "/home/j103t/odomLab/manuscript_scirocket/sx42b/raw_reads/Undetermined_S0_R1_001.fastq.gz",
#         "--r2",
#         "/home/j103t/odomLab/manuscript_scirocket/sx42b/raw_reads/Undetermined_S0_R2_001.fastq.gz",
#         "--experiment_name",
#         "sx42b",
#         "--samples",
#         "/omics/groups/OE0606/internal/jvanriet/git/manuscript_scirocket/config/testing_samplesheet.tsv",
#         "--barcodes",
#         "/omics/groups/OE0606/internal/jvanriet/git/sci-rocket/workflow/examples/example_barcodes.tsv",
#         "--out",
#         "/home/j103t/test/",
#     ]
# )
