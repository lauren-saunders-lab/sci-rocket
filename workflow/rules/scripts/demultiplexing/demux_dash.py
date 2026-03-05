import argparse
import json
import os
import pandas as pd
import pickle
import sys
import pathlib

HASH_RATIO_PASSING_THRESHOLD = 3.0
HASH_UMI_FILTER_THRESHOLD = 5
HASH_COUNT_BIN_LABELS = [f"{start}-{start + 9}" for start in range(1, 100, 10)] + ["100+"]
HASH_COUNT_BIN_BREAKS = list(range(0, 101, 10)) + [float("inf")]


def parse_float(value):
    return float(value.strip())


def parse_int(value):
    return int(float(value.strip()))


def parse_summary_metrics(path_summary):
    """
    Extract sample-wise STARsolo summary statistics from Summary.csv.
    """

    stats = {}
    with open(path_summary, "r") as handle:
        for line in handle:
            split_line = line.split(",", 1)
            if len(split_line) != 2:
                continue

            key, value = split_line[0], split_line[1].strip()
            if key == "Number of Reads":
                stats["total_reads"] = parse_int(value)
            elif key == "Sequencing Saturation":
                stats["sequencing_saturation"] = parse_float(value)
            elif key == "Reads Mapped to Genome: Unique+Multiple":
                stats["perc_mapped_reads_genome"] = parse_float(value)
            elif key == "Reads Mapped to Genome: Unique":
                stats["perc_unique_reads_genome_unique"] = parse_float(value)
            elif key == "Estimated Number of Cells":
                stats["estimated_cells"] = parse_int(value)
            elif key == "Mean Reads per Cell":
                stats["mean_reads_per_cell"] = parse_int(value)
            elif key == "Mean UMI per Cell":
                stats["mean_umi_per_cell"] = parse_int(value)
            elif (
                key.startswith("Reads Mapped to ")
                and "Genome" not in key
                and ": Unique+Multiple " in key
            ):
                stats["perc_mapped_reads_gene"] = parse_float(value)
            elif (
                key.startswith("Reads Mapped to ")
                and "Genome" not in key
                and ": Unique " in key
            ):
                stats["perc_unique_reads_gene_unique"] = parse_float(value)
            elif (
                key.startswith("Mean ")
                and key.endswith(" per Cell")
                and key not in {"Mean Reads per Cell", "Mean UMI per Cell"}
            ):
                stats["mean_genes_per_cell"] = parse_int(value)
    return stats


def resolve_solo_feature_paths(path_star, sample):
    """
    Resolve STARSolo summary and count files for one sample with feature fallback.
    """

    solo_dirs = sorted(pathlib.Path(path_star).glob(f"{sample}_*_Solo.out"))
    if not solo_dirs:
        raise FileNotFoundError(
            f"No STARSolo output directory found for sample '{sample}' in {path_star}"
        )

    solo_dir = solo_dirs[0]
    candidate_features = []
    for feature in ["GeneFull_Ex50pAS", "GeneFull_ExonOverIntron", "GeneFull", "Gene"]:
        if feature and feature not in candidate_features:
            candidate_features.append(feature)

    for summary in sorted(solo_dir.glob("*/Summary.csv")):
        feature = summary.parent.name
        if feature not in candidate_features:
            candidate_features.append(feature)

    for feature in candidate_features:
        feature_dir = solo_dir / feature
        summary = feature_dir / "Summary.csv"
        cellreads = feature_dir / "CellReads.stats"
        filtered_barcodes = feature_dir / "filtered" / "barcodes.tsv"
        if summary.exists() and cellreads.exists() and filtered_barcodes.exists():
            return {
                "feature": feature,
                "summary": summary,
                "cellreads": cellreads,
                "filtered_barcodes": filtered_barcodes,
            }

    raise FileNotFoundError(
        f"Could not resolve STARSolo Summary.csv, CellReads.stats and filtered/barcodes.tsv for sample '{sample}' in {solo_dir}"
    )


def get_benchmarks(path_benchmarks) -> list[dict]:
    frames = []
    for path in sorted(pathlib.Path(path_benchmarks).glob("*.txt")) + sorted(pathlib.Path(path_benchmarks).parent.glob("*.txt")):
        df = pd.read_csv(path, sep="\t")
        df.insert(0, "job", path.stem)
        frames.append(df)
    if not frames:
        return []
    df = pd.concat(frames, ignore_index=True, sort=False)
    return df.to_dict('records')


def calculate_hashing_cell_dominance_rows(df_hashing: pd.DataFrame) -> list[dict]:
    """
    Build cell-level dominance rows from hashing metrics.
    """

    if df_hashing.empty:
        return []

    rows = []
    for (sample_name, cell_barcode), group in df_hashing.groupby(["sample_name", "cell_barcode"], sort=False):
        ranked = group.sort_values(by=["count", "hashing_name"], ascending=[False, True]).reset_index(drop=True)
        hash_count_top = int(ranked.loc[0, "count"])
        hash_count_total = int(ranked["count"].sum())
        hash_umi_total = ranked.loc[0, "n_umi"] if "n_umi" in ranked.columns else hash_count_top
        hash_umi_total = int(hash_umi_total) if pd.notna(hash_umi_total) else 0
        top_hash = str(ranked.loc[0, "hashing_name"])

        hash_count_second = None
        if len(ranked) > 1:
            hash_count_second = int(ranked.loc[1, "count"])

        # Use a pseudocount of 1 when no second hash exists for this cell.
        if hash_count_second is None:
            hash_enrichment_ratio = float(hash_count_top)
        elif hash_count_second == 0:
            hash_enrichment_ratio = float("inf")
        else:
            hash_enrichment_ratio = float(hash_count_top / hash_count_second)

        rows.append(
            {
                "sample_name": str(sample_name),
                "cell_barcode": str(cell_barcode),
                "top_hash": top_hash,
                "hash_count_top": hash_count_top,
                "hash_count_second": hash_count_second,
                "hash_count_total": hash_count_total,
                "hash_umi_total": hash_umi_total,
                "hash_enrichment_ratio": hash_enrichment_ratio,
            }
        )

    return rows


def _get_hashing_cell_dominance_df(df_hashing: pd.DataFrame) -> pd.DataFrame:
    rows = calculate_hashing_cell_dominance_rows(df_hashing)
    columns = [
        "sample_name",
        "cell_barcode",
        "top_hash",
        "hash_count_top",
        "hash_count_second",
        "hash_count_total",
        "hash_umi_total",
        "hash_enrichment_ratio",
    ]
    return pd.DataFrame(rows, columns=columns)


def _summarize_hashing_cells(df_cells: pd.DataFrame, sample_names: list[str] | None = None) -> list[dict]:
    summary = []
    grouped = {sample_name: group for sample_name, group in df_cells.groupby("sample_name", sort=False)}
    ordered_sample_names = sample_names if sample_names is not None else list(grouped.keys())

    for sample_name in ordered_sample_names:
        if sample_name in grouped:
            group = grouped[sample_name]
        else:
            group = df_cells.iloc[0:0]

        ratios = group["hash_enrichment_ratio"].dropna()
        cells_passing = int((group["hash_enrichment_ratio"] >= HASH_RATIO_PASSING_THRESHOLD).fillna(False).sum())
        total_cells = int(len(group))

        mean_hash_count = float(group["hash_count_total"].mean()) if total_cells else None
        hash_count_total = int(group["hash_count_total"].sum()) if total_cells else 0
        hash_umi_total = int(group["hash_umi_total"].sum()) if total_cells else 0
        median_hash_ratio = float(ratios.median()) if len(ratios) else None
        mean_hash_ratio = float(ratios.mean()) if len(ratios) else None
        fraction_passing = float(cells_passing / total_cells) if total_cells else None

        summary.append(
            {
                "sample_name": str(sample_name),
                "mean_hash_count": mean_hash_count,
                "hash_count_total": hash_count_total,
                "hash_umi_total": hash_umi_total,
                "median_hash_ratio": median_hash_ratio,
                "mean_hash_ratio": mean_hash_ratio,
                "total_cells": total_cells,
                "cells_passing": cells_passing,
                "fraction_passing": fraction_passing,
            }
        )

    return summary


def _filter_hashing_cells(df_cells: pd.DataFrame, min_hash_umi_total: int) -> pd.DataFrame:
    if df_cells.empty:
        return df_cells.copy()
    return df_cells[df_cells["hash_umi_total"] >= min_hash_umi_total].copy()


def calculate_hashing_summary(
    df_hashing: pd.DataFrame,
    sample_names: list[str] | None = None,
) -> list[dict]:
    """
    Summarize cell-level hashing dominance statistics per sample.
    """
    df_cells = _get_hashing_cell_dominance_df(df_hashing)
    return _summarize_hashing_cells(df_cells, sample_names=sample_names)


def calculate_hashing_summary_filtered(
    df_hashing: pd.DataFrame,
    sample_names: list[str] | None = None,
    min_hash_umi_total: int = HASH_UMI_FILTER_THRESHOLD,
) -> list[dict]:
    """
    Summarize cell-level hashing dominance statistics per sample after
    filtering cells on hash_umi_total.
    """
    df_cells = _get_hashing_cell_dominance_df(df_hashing)
    df_cells_filt = _filter_hashing_cells(df_cells, min_hash_umi_total=min_hash_umi_total)
    return _summarize_hashing_cells(df_cells_filt, sample_names=sample_names)


def calculate_hashing_bin_summary(
    df_hashing: pd.DataFrame,
    sample_names: list[str] | None = None,
    min_hash_umi_total: int = HASH_UMI_FILTER_THRESHOLD,
) -> list[dict]:
    """
    Build sample/count-bin rows with passing fraction and n_cells.
    """
    df_cells = _get_hashing_cell_dominance_df(df_hashing)
    df_cells_filt = _filter_hashing_cells(df_cells, min_hash_umi_total=min_hash_umi_total)

    if sample_names is None:
        sample_names = list(dict.fromkeys(df_cells_filt["sample_name"].tolist()))

    if df_cells_filt.empty:
        by_sample_and_bin = {}
    else:
        df_cells_filt["count_bin"] = pd.cut(
            df_cells_filt["hash_count_total"],
            bins=HASH_COUNT_BIN_BREAKS,
            labels=HASH_COUNT_BIN_LABELS,
            include_lowest=True,
        )

        by_sample_and_bin = {}
        for (sample_name, count_bin), group in df_cells_filt.groupby(["sample_name", "count_bin"], sort=False, observed=True):
            n_cells = int(len(group))
            fraction_passing = float((group["hash_enrichment_ratio"] >= HASH_RATIO_PASSING_THRESHOLD).mean()) if n_cells else None
            by_sample_and_bin[(str(sample_name), str(count_bin))] = {
                "fraction_passing": fraction_passing,
                "n_cells": n_cells,
            }

    rows = []
    for sample_name in sample_names:
        sample_name = str(sample_name)
        for count_bin in HASH_COUNT_BIN_LABELS:
            cell = by_sample_and_bin.get((sample_name, count_bin))
            if cell is None:
                rows.append(
                    {
                        "sample_name": sample_name,
                        "count_bin": count_bin,
                        "fraction_passing": None,
                        "n_cells": 0,
                        "label": "·",
                    }
                )
                continue

            fraction_passing = cell["fraction_passing"]
            n_cells = cell["n_cells"]
            rows.append(
                {
                    "sample_name": sample_name,
                    "count_bin": count_bin,
                    "fraction_passing": fraction_passing,
                    "n_cells": n_cells,
                    "label": f"{fraction_passing * 100:.0f}% (n={n_cells})" if fraction_passing is not None else "·",
                }
            )

    return rows


def write_cell_hashing_table(qc, out):
    """
    Write the following cell-based metrics:
        - Total no. of hash reads.
        - Total no. of distinct UMI  per hashing barcode.

    Parameters:
        qc (dict): Dictionary containing the QC, incl. hashing metrics.
        out (str): Path to output hashing metrics.

    Returns:
        (dict): Dictionary of the hashing metrics for the dashboard hashing table.
        (list[dict]): Summary statistics derived from all cells.
        (list[dict]): Summary statistics with low hash_umi_total cells removed.
        (list[dict]): Per-bin sample summaries on filtered cells.
    """

    # Create a list of dictionaries.
    array_hashing = []
    
    def format_cell_barcode(cellular_barcode: str) -> str:
        # Cellular barcodes are stored internally as a 40nt concatenation:
        # p7(10nt) + p5(10nt) + ligation(10nt) + rt(10nt).
        # Export them as underscore-delimited components for readability.
        if len(cellular_barcode) == 40:
            return "_".join(
                (
                    cellular_barcode[0:10],
                    cellular_barcode[10:20],
                    cellular_barcode[20:30],
                    cellular_barcode[30:40],
                )
            )

        return cellular_barcode

    for sample_name in qc["hashing"]:
        for hashing_name in qc["hashing"][sample_name]:
            for cellular_barcode in qc["hashing"][sample_name][hashing_name]["counts"]:
                cell_metrics = qc["hashing"][sample_name][hashing_name]["counts"][cellular_barcode]
                array_hashing.append(
                    {
                        "experiment_name": qc["experiment_name"],
                        "sample_name": sample_name,
                        "hashing_name": hashing_name,
                        "cell_barcode": format_cell_barcode(cellular_barcode),
                        "cell_barcode_label": cell_metrics.get("cell_barcode_label", ""),
                        "count": cell_metrics["count"],
                        "n_umi": len(cell_metrics["umi"])
                    }
                )

    # Convert to pandas dataframe.
    df_hashing = pd.DataFrame(columns=["experiment_name", "sample_name", "hashing_name", "cell_barcode", "cell_barcode_label", "count", "n_umi"], data=array_hashing)

    # Order on total hash count (and by cell_barcode for equal counts).
    df_hashing = df_hashing.sort_values(by=["count", "cell_barcode"], ascending=False)

    # Generate folder.
    if not os.path.exists(os.path.dirname(out)):
        os.makedirs(os.path.dirname(out))

    # Write to file.
    df_hashing.to_csv(out, sep="\t", index=False, header=True, encoding="utf-8", mode="w")

    # Build sample-level summary statistics used in an extra dashboard table.
    hashing_summary = calculate_hashing_summary(df_hashing, sample_names=list(qc["hashing"].keys()))
    hashing_summary_filt = calculate_hashing_summary_filtered(df_hashing, sample_names=list(qc["hashing"].keys()))
    hashing_summary_bins = calculate_hashing_bin_summary(df_hashing, sample_names=list(qc["hashing"].keys()))

    # Transform into a dictionary for sci-dashboard.
    # Initialize the dictionary of sample_name and underlying hashing_name also a dictionary.
    dict_hashing = {sample_name : {hashing_name : {} for hashing_name in qc["hashing"][sample_name]} for sample_name in qc["hashing"]}

    for sample_name in qc["hashing"]:
        for hashing_name in qc["hashing"][sample_name]:
            dict_hashing[sample_name][hashing_name]["n_correct"] = qc["hashing"][sample_name][hashing_name]["n_correct"]
            dict_hashing[sample_name][hashing_name]["n_corrected"] = qc["hashing"][sample_name][hashing_name]["n_corrected"]
            dict_hashing[sample_name][hashing_name]["n_correct_upstream"] = qc["hashing"][sample_name][hashing_name]["n_correct_upstream"]

    return dict_hashing, hashing_summary, hashing_summary_filt, hashing_summary_bins


def combine_logs(path_pickle, path_star, path_hashing, path_benchmarks):
    """
    Combine the demuxxing logs with the STAR logs for the sci-dash.

    Parameters:
        path_pickle (str): Path to the pickled dictionaries.
        path_star (str): Path to the STAR output folder.
        path_hashing (str): Path to store the hashing metrics.
        path_benchmarks (str): Path to workflow benchmarks.

    Returns:
        qc (dict): Dictionary containing the demuxxing statistics (in JSON format).
    """

    # region Import demuxxing statistics. ---------------------------------------------------------------------------------
    qc_json = {}

    # Load the pickled dictionary
    with open(path_pickle, "rb") as handle:
        qc = pickle.load(handle)

    # endregion

    # region Calculate summary statistics. --------------------------------------------------------------------------------
    if qc != None:

        # Take over all the keys from the pickle which are not nested.
        for key in qc:
            if not isinstance(qc[key], dict):
                qc_json[key] = qc[key]

        # Determine the top recurrent uncorrectable barcodes.
        # Sort by descending value, and sort equal values by key to ensure the ordering is deterministic.
        qc_json["top_uncorrectables"] = {}
        top_n = 15
        # p5 barcodes are stored RC'd (as they appear in the fastq read header).
        # RC them back to match the orientation of the user-provided input barcodes.
        _rc = str.maketrans("ATCG", "TAGC")
        for key in ["p5", "p7", "ligation", "rt"]:
            barcodes = dict(qc[f"uncorrectable_{key}"])
            if key == "p5":
                barcodes = {seq[::-1].translate(_rc): count for seq, count in barcodes.items()}
            sorted_by_key = sorted(barcodes.items(), key=lambda item: (-item[1], item[0]))
            sorted_by_key_and_value = sorted(sorted_by_key, key=lambda x: x[1], reverse=True)[:top_n]
            qc_json["top_uncorrectables"][key] = [{"barcode": k, "frequency": v} for k, v in sorted_by_key_and_value]

        # Name the key/value pairs in the plate counts. Split up the key in row and col, and sort by descending frequency.
        def transform_plate_counts(dct):
            return [{"row": key[0], "col": key[1:].lstrip("0"), "frequency": dct[key]} for key in sorted(dct, key=lambda k: (-dct[k], k))]

        qc_json["p5_index_counts"] = transform_plate_counts(qc["p5_index_counts"])
        qc_json["p7_index_counts"] = transform_plate_counts(qc["p7_index_counts"])

        qc_json["ignore_p5"] = bool(qc["ignore_p5"])
        qc_json["ignore_p7"] = bool(qc["ignore_p7"])

        # Transform the rt_barcode_counts dictionary to a list of dictionaries.
        qc_json["rt_barcode_counts"] = {}

        for key in qc["rt_barcode_counts"]:
            qc_json["rt_barcode_counts"][key] = transform_plate_counts(qc["rt_barcode_counts"][key])

        # Transform the ligation_barcode_counts dictionary to a list of dictionaries.
        qc_json["ligation_barcode_counts"] = [{"barcode": barcode, "frequency": qc["ligation_barcode_counts"][barcode]} for barcode in qc["ligation_barcode_counts"]]

        # Remove all unnecessary data. ------------------------------------------------------------------------------------

        # Remove the ligation_barcode_counts with zero counts.
        qc_json["ligation_barcode_counts"] = [barcode for barcode in qc_json["ligation_barcode_counts"] if barcode["frequency"] > 0]

        # Convert the uncorrectables_sankey.
        qc_json["uncorrectables_sankey"] = [{"source": str(key), "value": qc["uncorrectables_sankey"][key]} for key in qc["uncorrectables_sankey"]]

    # endregion ----------------------------------------------------------------------------------------------------------

    # region Import STAR statistics. -------------------------------------------------------------------------------------

    # Per sample, load STARSolo summary and per-cell read stats.
    qc_json["sample_success"] = qc["sample_success"]
    for sample in qc_json["sample_success"]:
        solo_paths = resolve_solo_feature_paths(path_star, sample)
        summary_stats = parse_summary_metrics(solo_paths["summary"])
        qc_json["sample_success"][sample].update(summary_stats)

        # Load the CellReads.stats file and extract several sample-wise statistics.
        df_cellreads = pd.read_csv(solo_paths["cellreads"], sep="\t", header=0, index_col=0)

        # Only keep the filtered cells.
        filtered_index = pd.read_csv(
            solo_paths["filtered_barcodes"], sep="\t", header=None, index_col=0
        ).index
        df_cellreads = df_cellreads[df_cellreads.index.isin(filtered_index)]

        # Summarize all cells.
        df_cellreads_summed = df_cellreads.sum(axis=0)

        # Add the sample-wise statistics to the dictionary.
        qc_json["sample_success"][sample]["total_exonic_reads"] = int(df_cellreads_summed.exonic)
        qc_json["sample_success"][sample]["total_intronic_reads"] = int(df_cellreads_summed.intronic)
        qc_json["sample_success"][sample]["total_intergenic_reads"] = int(df_cellreads_summed.genomeU + df_cellreads_summed.genomeM - df_cellreads_summed.exonic - df_cellreads_summed.intronic)
        qc_json["sample_success"][sample]["total_mitochondrial_reads"] = int(df_cellreads_summed.mito)
        qc_json["sample_success"][sample]["total_exonicAS_reads"] = int(df_cellreads_summed.exonicAS)
        qc_json["sample_success"][sample]["total_intronicAS_reads"] = int(df_cellreads_summed.intronicAS)
        
    # endregion ----------------------------------------------------------------------------------------------------------

    # region Write / import hashing statistics. --------------------------------------------------------------------------
    if "hashing" in qc:
        dict_hashing, hashing_summary, hashing_summary_filt, hashing_summary_bins = write_cell_hashing_table(qc, path_hashing)
        qc_json["hashing"] = dict_hashing
        qc_json["hashing_summary"] = hashing_summary
        qc_json["hashing_summary_filt"] = hashing_summary_filt
        qc_json["hashing_summary_bins"] = hashing_summary_bins
        qc_json["hashing_summary_bin_labels"] = HASH_COUNT_BIN_LABELS
    else:
        qc_json["hashing"] = {}
        qc_json["hashing_summary"] = []
        qc_json["hashing_summary_filt"] = []
        qc_json["hashing_summary_bins"] = []
        qc_json["hashing_summary_bin_labels"] = HASH_COUNT_BIN_LABELS
    # endregion ----------------------------------------------------------------------------------------------------------

    qc_json["benchmarks"] = get_benchmarks(path_benchmarks)

    # Return JSON-structured dict.
    return qc_json


def main(arguments):
    # Setup argument parser.
    parser = argparse.ArgumentParser(description="Combine the scattered demultiplexing files into a JSON structure for the sci-dash.", add_help=False)
    parser.add_argument("--path_pickle", required=True, type=str, help="(str) Path to demultiplexing qc pickle.")
    parser.add_argument("--path_star", required=True, type=str, help="(str) Path to the star alignment folder.")
    parser.add_argument("--path_out", required=True, type=str, help="(str) Path to store JSON structure.")
    parser.add_argument("--path_hashing", required=True, type=str, help="(str) Path to store hashing metrics (if applicable).")
    parser.add_argument("--path_benchmarks", required=True, type=str, help="(str) Path to workflow benchmarks.")

    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help="Display help and exit.")

    # Parse arguments.
    args = parser.parse_args()

    # Combine the demuxxing logs with the STAR logs for the sci-dash.
    qc_json = combine_logs(
        args.path_pickle,
        args.path_star,
        args.path_hashing,
        args.path_benchmarks,
    )

    # Write the JSON structure to file.
    if not os.path.exists(os.path.dirname(args.path_out)):
        os.makedirs(os.path.dirname(args.path_out))

    with open(args.path_out, "w") as handle:
        handle.write("var data = ")
        json.dump(qc_json, handle, indent=4)

    # Close the file
    handle.close()


if __name__ == "__main__":
    main(sys.argv[1:])
    sys.exit()
