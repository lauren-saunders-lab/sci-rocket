import pytest
import pathlib
import workflow.rules.scripts.demultiplexing.preprocess as preprocess
import pandas as pd


def write_tsv(path, header, rows):
    pathlib.Path(path).write_text("\t".join(header) + "\n" + "\n".join("\t".join(map(str, r)) for r in rows) + "\n")


@pytest.fixture
def test_config(tmp_path) -> dict:
    # Minimal samplesheet
    samples_path = tmp_path / "samples.tsv"
    write_tsv(samples_path,
              header=["path_reads", "experiment_name", "p5", "p7", "rt", "sample_name", "species", "n_expected_cells"],
              rows=[["run", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "sample1", "mouse", "10000"]])

    # Minimal barcodes file
    barcodes_path = tmp_path / "barcodes.tsv"
    write_tsv(barcodes_path, header=["type", "barcode", "sequence"],
              rows=[["ligation", "LIG1", "ACTTGATTGT"],
                    ["p5", "A01", "CGTTCTATCA"],
                    ["p7", "A01", "CTAAGCCTTG"],
                    ["rt", "P01-B02", "GCCGCAACGA"],
                    ])

    return {"path_samples": str(samples_path),
            "path_barcodes": str(barcodes_path),
            "species": {"mouse": {"genome": "", "genome_gtf": "", "star_index": ""}}}


def test_get_samples_deduplicated_sequencing_names(test_config: dict):
    write_tsv(test_config["path_samples"],
              header=["path_reads", "experiment_name", "p5", "p7", "rt", "sample_name", "species", "n_expected_cells"],
              rows=[
                  ["/path/r1", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "sample1", "mouse", "10000"],
                  ["/other/r1", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "sample2", "mouse", "10000"],
                  ["/path/r1", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "sample3", "mouse", "10000"],
                  ["/x/r1", "exp2", "A02:H02", "H01:H12", "P01-A01:P01-D12", "sample4", "mouse", "10000"],
                  ["/other/r2", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "sample5", "mouse", "10000"],
                  ["/path/r1", "exp2", "A02:H02", "H01:H12", "P01-A01:P01-D12", "sample6", "mouse", "10000"],
                  ["/other/r1", "exp", "A01:B01", "H06:H12", "P01-A01:P01-D12", "sample7", "mouse", "10000"],
              ], )
    samples = preprocess.get_samples(test_config)
    assert samples["sequencing_name"][0] == "r1"  # exp/r1: first r1 in exp, so ok as is
    assert samples["sequencing_name"][
               1] == "r1_1"  # exp/r1_1: second r1 in exp & has different path (/other/r1) to first one, so renamed
    assert samples["sequencing_name"][2] == "r1"  # exp/r1: repeat of first r1 path in exp, so ok as is
    assert samples["sequencing_name"][3] == "r1"  # exp2/r1: different experiment, so ok as is
    assert samples["sequencing_name"][4] == "r2"  # exp/r2: first r2 in exp, so ok as is
    assert samples["sequencing_name"][5] == "r1_1"  # exp2/r1_1: second r1 in exp2 with different path, so renamed
    assert samples["sequencing_name"][6] == "r1_1"  # exp1/r1_1: same path & exp as sample 2


def test_get_samples_path_reads_with_multiple_paths(test_config: dict):
    write_tsv(
        test_config["path_samples"],
        header=["path_reads", "experiment_name", "p5", "p7", "rt", "sample_name", "species", "n_expected_cells"],
        rows=[
            # single path, no semicolon -> single row
            ["path/r1", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s1", "mouse", "10000"],
            # single path, trailing semicolon -> single row
            ["/path/r2;", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s2", "mouse", "10000"],
            # two paths separated by a semicolon with trailing slashes -> two rows without trailing slashes
            ["/path/a/;/path/b/", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s3", "mouse", "10000"],
            # four paths separated by semicolons -> four rows
            ["/path/a1;/path/a2;/path/a3;/path/a4", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s3", "mouse",
             "10000"],
            # duplicate paths in same cell -> single row with single path
            ["c://dup;c://dup", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s4", "mouse", "10000"],
            # whitespace around tokens -> single rows with trimmed paths
            [" /ws1 ;  /ws2  ", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s5", "mouse", "10000"],
            # empty tokens from double semicolons + trailing -> two rows with trimmed paths
            [";/x;;/y;", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s6", "mouse", "10000"],
        ],
    )

    expected_samples = pd.DataFrame([
        ["path/r1", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s1", "mouse", "10000", "r1"],
        ["/path/r2", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s2", "mouse", "10000", "r2"],
        ["/path/a", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s3", "mouse", "10000", "a"],
        ["/path/b", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s3", "mouse", "10000", "b"],
        ["/path/a1", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s3", "mouse", "10000", "a1"],
        ["/path/a2", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s3", "mouse", "10000", "a2"],
        ["/path/a3", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s3", "mouse", "10000", "a3"],
        ["/path/a4", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s3", "mouse", "10000", "a4"],
        ["c://dup", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s4", "mouse", "10000", "dup"],
        ["/ws1", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s5", "mouse", "10000", "ws1"],
        ["/ws2", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s5", "mouse", "10000", "ws2"],
        ["/x", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s6", "mouse", "10000", "x"],
        ["/y", "exp", "A02:H02", "H01:H12", "P01-A01:P01-D12", "s6", "mouse", "10000", "y"],
    ], columns=["path_reads", "experiment_name", "p5", "p7", "rt", "sample_name", "species", "n_expected_cells", "sequencing_name"])

    samples = preprocess.get_samples(test_config)
    pd.testing.assert_frame_equal(samples.reset_index(drop=True), expected_samples, check_dtype=False)


def test_get_samples_ignores_comment_lines(test_config: dict):
    pathlib.Path(test_config["path_samples"]).write_text(
        "\n".join(
            [
                "# comment line should be ignored",
                "path_reads\texperiment_name\tp5\tp7\trt\tsample_name\tspecies\tn_expected_cells",
                "/path/r1\texp\tA02:H02\tH01:H12\tP01-A01:P01-D12\tsample1\tmouse\t10000",
                "# second comment should also be ignored",
                "/path/r2\texp\tA02:H02\tH01:H12\tP01-A01:P01-D12\tsample2\tmouse\t10000",
            ]
        )
        + "\n"
    )

    samples = preprocess.get_samples(test_config)

    assert len(samples) == 2
    assert set(samples["sample_name"]) == {"sample1", "sample2"}

# --- IGNORE barcode tests ---


class TestRetrieveBarcodesIgnore:
    """Tests for IGNORE handling in retrieve_barcodes()."""

    def test_ignore_returns_ignore_list(self):
        """retrieve_barcodes returns ["IGNORE"] when given "IGNORE"."""
        log = preprocess.init_logger()
        barcodes = pd.DataFrame({"type": ["p5"], "barcode": ["A01"], "sequence": ["CGTTCTATCA"]})
        result = preprocess.retrieve_barcodes(log, pd.Series(["IGNORE"]), barcodes, "p5")
        assert result == ["IGNORE"]


class TestSanitySamplesIgnore:
    """Tests for IGNORE validation in sanity_samples()."""

    @pytest.fixture
    def log(self):
        return preprocess.init_logger()

    @pytest.fixture
    def barcodes(self):
        return pd.DataFrame({
            "type": ["p5", "p5", "p7", "p7", "ligation", "rt"],
            "barcode": ["A01", "A02", "A01", "A02", "LIG1", "P01-A01"],
            "sequence": ["CGTTCTATCA", "AATTCCGGAA", "CTAAGCCTTG", "GGAATTCCAA", "ACTTGATTGT", "GCCGCAACGA"],
        })

    @pytest.fixture
    def config(self):
        return {"species": {"mouse": {"genome": "", "genome_gtf": "", "star_index": ""}}}

    def test_ignore_p5_all_samples_passes(self, log, barcodes, config):
        """sanity_samples passes when p5 is IGNORE for all samples."""
        samples = pd.DataFrame({
            "experiment_name": ["exp", "exp"],
            "p5": ["IGNORE", "IGNORE"],
            "p7": ["A01", "A02"],
            "rt": ["P01-A01", "P01-A01"],
            "sample_name": ["s1", "s2"],
            "species": ["mouse", "mouse"],
            "n_expected_cells": ["1000", "1000"],
            "path_reads": ["/a", "/b"],
        })
        assert preprocess.sanity_samples(log, samples, barcodes, config) is True

    def test_ignore_p5_mixed_fails(self, log, barcodes, config):
        """sanity_samples fails when IGNORE is mixed with real p5 values."""
        samples = pd.DataFrame({
            "experiment_name": ["exp", "exp"],
            "p5": ["IGNORE", "A01"],
            "p7": ["A01", "A02"],
            "rt": ["P01-A01", "P01-A01"],
            "sample_name": ["s1", "s2"],
            "species": ["mouse", "mouse"],
            "n_expected_cells": ["1000", "1000"],
            "path_reads": ["/a", "/b"],
        })
        assert preprocess.sanity_samples(log, samples, barcodes, config) is False

    def test_ignore_p5_mixed_between_experiments_passes(self, log, barcodes, config):
        """sanity_samples passes when IGNORE and real p5 values are separated by experiment."""
        samples = pd.DataFrame({
            "experiment_name": ["exp_ignore", "exp_ignore", "exp_real", "exp_real"],
            "p5": ["IGNORE", "IGNORE", "A01", "A02"],
            "p7": ["A01", "A02", "A01", "A02"],
            "rt": ["P01-A01", "P01-A01", "P01-A01", "P01-A01"],
            "sample_name": ["s1", "s2", "s3", "s4"],
            "species": ["mouse", "mouse", "mouse", "mouse"],
            "n_expected_cells": ["1000", "1000", "1000", "1000"],
            "path_reads": ["/a", "/b", "/c", "/d"],
        })
        assert preprocess.sanity_samples(log, samples, barcodes, config) is True

    def test_ignore_both_p5_and_p7_in_single_experiment_fails(self, log, barcodes, config):
        """sanity_samples fails when both P5 and P7 are IGNORE in one experiment."""
        samples = pd.DataFrame({
            "experiment_name": ["exp_bad", "exp_bad", "exp_ok"],
            "p5": ["IGNORE", "IGNORE", "A01"],
            "p7": ["IGNORE", "IGNORE", "A01"],
            "rt": ["P01-A01", "P01-A01", "P01-A01"],
            "sample_name": ["s1", "s2", "s3"],
            "species": ["mouse", "mouse", "mouse"],
            "n_expected_cells": ["1000", "1000", "1000"],
            "path_reads": ["/a", "/b", "/c"],
        })
        assert preprocess.sanity_samples(log, samples, barcodes, config) is False

    def test_ignore_both_p5_and_p7_fails(self, log, barcodes, config):
        """sanity_samples fails when both P5 and P7 are IGNORE."""
        samples = pd.DataFrame({
            "experiment_name": ["exp", "exp"],
            "p5": ["IGNORE", "IGNORE"],
            "p7": ["IGNORE", "IGNORE"],
            "rt": ["P01-A01", "P01-A01"],
            "sample_name": ["s1", "s2"],
            "species": ["mouse", "mouse"],
            "n_expected_cells": ["1000", "1000"],
            "path_reads": ["/a", "/b"],
        })
        assert preprocess.sanity_samples(log, samples, barcodes, config) is False
