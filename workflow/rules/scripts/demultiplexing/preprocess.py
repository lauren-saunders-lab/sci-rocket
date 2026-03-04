import re
import logging
import pandas as pd
import numpy as np

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

def retrieve_barcodes(log: logging.Logger, requested_barcodes: pd.Series, available_barcodes: pd.DataFrame, barcode_type: str):
    """
    Determines the sanity of the requested barcodes and returns all barcodes as determined by this format:
        - A single barcode (e.g., A01)
        - A range of barcodes (e.g., A01:A12)
        - Multiple barcodes (e.g., A01,A02,A03)
        - Multiple ranges of barcodes (e.g., A01:A12,B01:B12)

    Parameters:
        log (logging.Logger): Logger object.
        requested_barcodes (pd.Series): Series of requested barcodes.
        available_barcodes (pd.DataFrame): DataFrame of available barcodes.
        barcode_type (str): Barcode type (e.g., p5, p7, rt).

    Returns:
        list or boolean: List of all requested barcodes if the sanity check passed, False otherwise.
    """

    # If IGNORE is specified, return it as-is (resolved later in demux_rocket.py).
    if list(requested_barcodes) == ["IGNORE"]:
        return ["IGNORE"]

    # Get all the available barcodes of the requested type.
    available_barcodes = available_barcodes.query("type == @barcode_type")["barcode"].unique()

    # Split the requested barcodes in case of multiple strips (, separated; e.g., A01:A12,B01:B12).
    requested_barcodes = [i.split(",") for i in requested_barcodes]
    requested_barcodes = list(set([j for i in requested_barcodes for j in i]))

    # Generate the list of all requested barcodes based on supplied ranges.
    experiment_barcodes = []

    for barcode_range in requested_barcodes:
        # Split the requested barcodes in case of ranges (e.g., A01:A12).
        barcode_range = barcode_range.split(":")

        if len(barcode_range) == 1:
            experiment_barcodes += [barcode_range[0]]
        else:
            # Per barcode range, check if the start index is smaller than the end index (if applicable).
            if barcode_range[0] > barcode_range[1]:
                log.error("Sanity check (Sample sheet) - {} start index is larger than the end index: {} > {}".format(barcode_type.upper(), barcode_range[0], barcode_range[1]))
                return False
            else:
                # Determine the various parts of the barcode (e.g., A01:A12 -> A, 01, A, 12) and for RT also the P0x- number.
                if barcode_type == 'rt':
                    plate_start = barcode_range[0].split("-")[0]
                    plate_end = barcode_range[1].split("-")[0]
                    letter_start = barcode_range[0].split("-")[1][0]
                    letter_end = barcode_range[1].split("-")[1][0]
                    number_start = int(barcode_range[0].split("-")[1][1:])
                    number_end = int(barcode_range[1].split("-")[1][1:])
                else:
                    letter_start = barcode_range[0][0]
                    letter_end = barcode_range[1][0]
                    number_start = int(barcode_range[0][1:])
                    number_end = int(barcode_range[1][1:])

                # Sanity checks. -----------------------------------------------------------------------
                if barcode_type == 'p7' and letter_start != letter_end:
                    log.error("Sanity check (Sample sheet) - {} start and end index are not on the same row (should be {}: {} and {}".format(barcode_type.upper(), letter_start, barcode_range[0], barcode_range[1]))
                    return False

                if barcode_type == 'p5' and number_start != number_end:
                    log.error("Sanity check (Sample sheet) - {} start and end index are not on the same column (should be {}: {} and {}".format(barcode_type.upper(), "0{}".format(number_start) if number_start < 10 else number_start, barcode_range[0], barcode_range[1]))
                    return False
                
                if barcode_type == 'rt' and plate_start != plate_end:
                    log.error("Sanity check (Sample sheet) - {} start and end index are not on the same plate (should be {}: {} and {}".format(barcode_type.upper(), plate_start, barcode_range[0], barcode_range[1]))
                    return False
                
                # ---------------------------------------------------------------------------------------

                # Generate the barcodes.
                if barcode_type == 'p5':
                    experiment_barcodes += ["{}{}".format(letter, "0{}".format(number_start) if number_start < 10 else number_start) for letter in list(map(chr, range(ord(letter_start), ord(letter_end) + 1)))]
                elif barcode_type == 'p7':
                    experiment_barcodes += ["{}{}".format(letter_start, "0{}".format(number) if number < 10 else number) for number in range(number_start, number_end + 1)]
                elif barcode_type == 'rt':
                    # Generate all letter and number combinations for the start and end letter/number.
                    experiment_barcodes += ["{}-{}{}".format(plate_start, letter, "0{}".format(number) if number < 10 else number) for letter in list(map(chr, range(ord(letter_start), ord(letter_end) + 1))) for number in range(number_start, number_end + 1)]

    return experiment_barcodes


def sanity_samples(log, samples, barcodes, config):
    """
    Checks the provided sample sheet for sanity and inconsistensies with the supplied config.

    Args:
        log (logging.Logger): Logger object.
        samples (pandas.DataFrame): Imported sample-sheet.
        barcodes (pandas.DataFrame): Imported barcodes.
        config (dict): Imported config.

    Returns:
        bool: True if the sample sheet is valid, False otherwise.
    """

    # Check if the sample sheet is empty.
    if samples.empty:
        log.error("Sanity check (Sample sheet) - Sample sheet is empty!")
        return False

    # Check if the sample sheet contains the required columns.
    required_columns = set(["experiment_name", "p5", "p7", "rt", "sample_name", "species", "n_expected_cells"])
    if not required_columns.issubset(samples.columns):
        log.error("Sanity check (Sample sheet) - Missing required column(s): {}".format(", ".join(required_columns.difference(samples.columns))))
        return False
    
    # Check if path, path_bcl or path_fastq exist in the sample sheet.
    if "path_reads" not in samples.columns and "path_bcl" not in samples.columns and "path_fastq" not in samples.columns:
        log.error("Sanity check (Sample sheet) - Missing 'path' (or 'path_bcl' or 'path_fastq') column.")
        return False

    # Check if the sample sheet contains species which are not defined in the config.
    if not set(samples["species"].unique()).issubset(config["species"].keys()):
        log.error("Sanity check (Sample sheet) - Contains species which are not defined in the config: {}".format(", ".join(set(samples["species"].unique()).difference(config["species"].keys()))))
        return False

    # Check if different species are assigned to the same sample name.
    if samples.groupby("sample_name")["species"].apply(lambda x: len(x.unique()) > 1).any():
        # Check which samples have different species assigned.
        samples_with_different_species = samples.groupby("sample_name")["species"].apply(lambda x: ", ".join(x.unique())).reset_index()

        # Print the samples with different species assigned.
        for _, row in samples_with_different_species.iterrows():
            if len(row["species"].split(", ")) > 1:
                log.error("Sanity check (Sample sheet) - Different species assigned to {}: {}".format(row["sample_name"], row["species"]))

        return False

    # region Sanity of p5/p7 indexes ---------------------------------------------------------------

    indexes_p5 = retrieve_barcodes(log, samples["p5"].unique(), barcodes, "p5")
    indexes_p7 = retrieve_barcodes(log, samples["p7"].unique(), barcodes, "p7")
    indexes_rt = retrieve_barcodes(log, samples["rt"].unique(), barcodes, "rt")

    if not indexes_p5 or not indexes_p7 or not indexes_rt:
        return False

    # Validate IGNORE usage for p5/p7 per experiment.
    for experiment_name, samples_exp in samples.groupby("experiment_name"):
        for barcode_type in ["p5", "p7"]:
            values = samples_exp[barcode_type].unique()
            has_ignore = "IGNORE" in values
            if has_ignore and len(values) > 1:
                log.error(
                    "Sanity check (Sample sheet) - {} is set to IGNORE for some samples but not all in "
                    "experiment {}. IGNORE must be used for all samples in that experiment or none."
                    .format(barcode_type.upper(), experiment_name)
                )
                return False

        if samples_exp["p5"].unique().tolist() == ["IGNORE"] and samples_exp["p7"].unique().tolist() == ["IGNORE"]:
            log.error("Sanity check (Sample sheet) - Both P5 and P7 cannot be set to IGNORE in experiment {}."
                      .format(experiment_name))
            return False

    # endregion -------------------------------------------------------------------------------------

    # region Sanity of hash heets -------------------------------------------------------------------

    # For samples which have a designated hashing sheet, check sanity of hashing sheet.
    if "hashing" in samples.columns:
        for hash_sheet in samples[~samples['hashing'].isna()].hashing.unique():
            # Open the hashing sheets (.tsv) and check if hash_name and barcode columns are present.
            x = pd.read_csv(hash_sheet, sep="\t", header=0)
            required_columns = set(["hash_name", "barcode"])

            if not required_columns.issubset(x.columns):
                log.error("Sanity check (Sample sheet) - Hashing sheet {} is missing required column(s): {}".format(hash_sheet, ", ".join(required_columns.difference(x.columns))))
                return False

            # Check if the hashing sheet contains duplicate barcodes names.
            if x["barcode"].duplicated().any():
                duplicated_barcodes = x["barcode"].duplicated()
                log.error("Sanity check (Sample sheet) - Hashing sheet {} contains duplicate barcodes: {}".format(hash_sheet, ", ".join(x.loc[duplicated_barcodes, "barcode"])))
                return False

    # endregion ---------------------------------------------------------------------------------------

    # Otherwise, the sample sheet is valid.
    return True


def sanity_barcodes(log, barcodes):
    """
    Checks the provided barcodes file for sanity.

    Args:
        log (logging.Logger): Logger object.
        barcodes (pandas.DataFrame): Imported barcodes.

    Returns:
        bool: True if the barcodes file is valid, False otherwise.
    """

    # Check if the barcodes file is empty.
    if barcodes.empty:
        log.error("Sanity check (Barcodes) - Barcodes file is empty!")
        return False

    # Check if the barcodes file contains the required columns.
    required_columns = set(["type", "barcode", "sequence"])
    if not set(required_columns).issubset(barcodes.columns):
        log.error("Sanity check (Barcodes) - Barcodes file is missing required column(s): {}".format(", ".join(required_columns.difference(barcodes.columns))))
        return False

    # Check if all barcode types are defined.
    required_columns = set(["rt", "p5", "p7", "ligation"])
    if not set(required_columns).issubset(barcodes["type"].unique()):
        log.error("Sanity check (Barcodes) - Barcodes file is missing barcode type(s): {}".format(", ".join(required_columns.difference(barcodes["type"].unique()))))
        return False

    # Check if the barcodes file contains duplicate barcodes per type.
    if barcodes.groupby("type")["barcode"].apply(lambda x: x.duplicated().any()).any():
        # Print the duplicate barcodes per type with barcode (type)
        duplicated_barcodes = barcodes.groupby("type")["barcode"].apply(lambda x: x[x.duplicated()]).unique()
        for barcode in duplicated_barcodes:
            log.error("Sanity check (Barcodes) - Barcode {} is duplicated in type {}".format(barcode, barcodes.query("barcode == @barcode")["type"].unique()[0]))

        return False

    # Check if barcode sequences contain only valid nucleotides (ATCG).
    if not barcodes["sequence"].str.contains("^[ATCG]+$").all():
        # Print the invalid barcodes
        invalid_barcodes = barcodes[~barcodes["sequence"].str.contains("^[ATCG]+$")]
        log.error("Sanity check (Barcodes) - Barcodes file contains invalid nucleotides (only ATCG allowed):\n{}".format(invalid_barcodes))

        return False

    # If all checks pass, the barcodes file is valid.
    return True


def check_sanity(samples, barcodes, config):
    """
    Checks for the sanity of the sample sheet and barcodes file.

    Args:
        samples (pandas.DataFrame): Imported sample-sheet.
        barcodes (pandas.DataFrame): Imported barcodes.
        config (dict): Imported config.

    Returns:
        bool: True if the sample sheet and barcodes file are valid, False otherwise.
    """

    # Initialize logging.
    log = init_logger()

    if not sanity_barcodes(log, barcodes):
        raise ValueError("Barcode file is not valid.")

    if not sanity_samples(log, samples, barcodes, config):
        raise ValueError("Sample metadata is not valid.")

    # Close Logger.
    logging.shutdown()


def get_wildcard_constraints(samples: pd.DataFrame, cols: list[str]) -> dict[str,str]:
    """
    Generate wildcard constraints based on the values in each specified column of the samples DataFrame.

    Parameters:
        samples (pd.DataFrame): DataFrame of samples.
        cols (list[str]): List of columns to generate wildcard constraints for.

    Returns:
        dict: Dictionary of wildcard constraints.
    """

    return {col: "|".join([re.escape(str(x)) for x in samples[col]]) for col in cols}


def get_samples(config: dict) -> pd.DataFrame:
    """
    Load the sample sheet and barcodes file, perform sanity checks, return unique samples.

    Parameters:
        config (dict): Imported config.
    Returns:
        pd.DataFrame: DataFrame of unique samples.
    """

    # Load sample sheet.
    samples = pd.read_csv(config["path_samples"], sep="\t", dtype=str, comment="#")

    # Load barcodes file.
    barcodes = pd.read_csv(config["path_barcodes"], sep="\t", dtype=str)

    # Perform sanity checks.
    check_sanity(samples, barcodes, config)

    # If 'path_reads' column is not provided, use 'path_bcl' or 'path_fastq' column as 'path_reads'.
    if "path_reads" not in samples.columns:
        if "path_fastq" in samples.columns and "path_bcl" in samples.columns:
            print("Both BCL and FASTQ paths are specified in the samples file - using FASTQ.")
        for col in ["path_bcl", "path_fastq"]:
            if col in samples.columns:
                samples["path_reads"] = samples[col]
                samples.drop(col, axis=1, inplace=True)

    # Split any path_reads entries which contain multiple semicolon-separated paths into multiple rows.
    samples = (
        samples
        .assign(path_reads=samples["path_reads"].str.split(";"))
        .explode("path_reads", ignore_index=True)
        .assign(path_reads=lambda d: d["path_reads"].str.strip().str.rstrip("/"))
        .query("path_reads != ''")
    )

    # Drop duplicate rows based on experiment_name, sample_name, species and path_reads.
    samples.drop_duplicates(subset=["path_reads", "experiment_name", "sample_name", "species"], inplace=True)

    # Generate the sequencing_name based on the last folder name of the path
    sequencing_name = samples["path_reads"].str.split("/").str[-1]
    # Append a counter to sequencing_name in case of duplicates with differing paths within the same experiment.
    # This avoids overwriting data for samples with different paths but same sequencing_name within the same experiment.
    path_id = (
        samples.groupby(["experiment_name", sequencing_name])["path_reads"]
        .transform(lambda s: pd.factorize(s)[0])
    )
    samples["sequencing_name"] = sequencing_name + path_id.where(path_id.eq(0), "_" + path_id.astype(str)).replace(0, "")

    return samples


def get_haplotyping_samples(samples_unique: pd.DataFrame) -> pd.DataFrame:
    """
    Returns haplotyping samples from provided unique samples DataFrame.

    Parameters:
        pd.DataFrame: DataFrame of unique samples.
    Returns:
        pd.DataFrame: DataFrame of haplotyping samples.
    """
    if "strain1" in samples_unique.columns and "strain2" in samples_unique.columns:
        return samples_unique[
            samples_unique[["strain1", "strain2"]].notna().all(axis=1) &
            (samples_unique["species"] == "mouse")
            ]
    else:
        return pd.DataFrame(columns=["experiment_name", "sample_name", "strain1", "strain2"])
