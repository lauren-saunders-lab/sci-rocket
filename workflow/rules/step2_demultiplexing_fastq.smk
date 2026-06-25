#   * Generate sample-specific .fastq.gz files by demultiplexing the sci-seq barcodes. *
#
#   1. split_reads:                         Split R1 and R2 files into smaller files which will be handled in parallel.
#   3. demultiplex_fastq_split:             Demultiplex each split R1 and R2 file to generate sample-specific fastq.gz files (in parallel).
#   4. gather_demultiplexed_sequencing:     Combine the discarded and whitelist files (from multiple parallel jobs).
#   5. gather_demultiplexed_samples:        Combine the sample-specific fastq.gz files (from multiple parallel jobs).
#############


def get_fastq_split_parts():
    n_parts = int(config["settings"]["scatter_fastq_split"])
    if n_parts < 1:
        raise ValueError("settings.scatter_fastq_split must be >= 1")
    return n_parts


def get_fastq_split_scatteritems():
    n_parts = get_fastq_split_parts()
    width = max(3, len(str(n_parts)))
    return [f"{i:0{width}d}" for i in range(1, n_parts + 1)]


def get_sequencing_runs(experiment_name):
    """Return a list of sequencing runs for an experiment."""
    return samples_unique.query("experiment_name == @experiment_name").sequencing_name.unique().tolist()


# ---- Split R1 and R2 files into smaller files which will be handled in parallel. ----
rule split_reads:
    input:
        fastqs=lambda w: expand(
            out("{experiment_name}/raw_reads/{sequencing_name}/Undetermined_S0_{read}_001.fastq.gz"),
            experiment_name=w.experiment_name,
            sequencing_name=get_sequencing_runs(w.experiment_name),
            read=w.read,
        ),
    output:
        temp(
            expand(
                out("{{experiment_name}}/raw_reads_split/{{read}}_{scatteritem}.fastq.gz"),
                scatteritem=get_fastq_split_scatteritems(),
            )
        ),
    wildcard_constraints:
        read="R[12]",
    threads: 5
    resources:
        mem_mb=1024 * 10,
    log:
        out("logs/step2_demultiplexing_reads/split_{experiment_name}_{read}.log"),
    benchmark:
        out("benchmarks/{experiment_name}/split_{read}.txt")
    params:
        n_parts=get_fastq_split_parts(),
        out_dir=out("{experiment_name}/raw_reads_split"),
    conda:
        "../envs/sci-rocket.yaml",
    shell:
        r"""
        exec > "{log}" 2>&1
        set -euo pipefail

        # Stream per-run gzipped fastqs into seqkit split2 via a named fifo.
        # cat-ing gzip files yields a valid multi-member gzip stream, and the
        # .fastq.gz suffix on the fifo path makes seqkit emit gzipped chunks.
        LOCAL_SPLIT="{resources.tmpdir}/split_{wildcards.experiment_name}_{wildcards.read}"
        mkdir -p "$LOCAL_SPLIT"
        mkdir -p "{params.out_dir}"

        fifo="$(mktemp -u --suffix=.fastq.gz)"
        mkfifo "$fifo"
        trap 'rm -f "$fifo"; rm -rf "$LOCAL_SPLIT"' EXIT

        cat {input.fastqs:q} > "$fifo" &
        cat_pid=$!

        seqkit split2 \
          -p {params.n_parts} \
          -j {threads} \
          -O "$LOCAL_SPLIT" \
          --by-part-prefix "{wildcards.read}_" \
          "$fifo"

        wait "$cat_pid"

        mv "$LOCAL_SPLIT"/*.fastq.gz "{params.out_dir}/"
        """

# ---- Helper for explicit per-sample scatter outputs in demultiplex_fastq_split. ----
def get_sample_names():
    return sorted(samples_unique["sample_name"].unique().tolist())


# ---- Demultiplex each split R1 and R2 file and emit sample-specific fastq.gz files. ----
rule demultiplex_fastq_split:
    input:
        R1=out("{experiment_name}/raw_reads_split/R1_{scatteritem}.fastq.gz"),
        R2=out("{experiment_name}/raw_reads_split/R2_{scatteritem}.fastq.gz"),
    output:
        discard_R1=temp(out("{experiment_name}/demux_reads_scatter/{scatteritem}/{experiment_name}_R1_discarded.fastq.gz")),
        discard_R2=temp(out("{experiment_name}/demux_reads_scatter/{scatteritem}/{experiment_name}_R2_discarded.fastq.gz")),
        discard_log=temp(out("{experiment_name}/demux_reads_scatter/{scatteritem}/log_{experiment_name}_discarded_reads.tsv.gz")),
        qc=temp(out("{experiment_name}/demux_reads_scatter/{scatteritem}/{experiment_name}_qc.pickle")),
        whitelist_p7=temp(out("{experiment_name}/demux_reads_scatter/{scatteritem}/{experiment_name}_whitelist_p7.txt")),
        whitelist_p5=temp(out("{experiment_name}/demux_reads_scatter/{scatteritem}/{experiment_name}_whitelist_p5.txt")),
        whitelist_ligation=temp(out("{experiment_name}/demux_reads_scatter/{scatteritem}/{experiment_name}_whitelist_ligation.txt")),
        whitelist_rt=temp(out("{experiment_name}/demux_reads_scatter/{scatteritem}/{experiment_name}_whitelist_rt.txt")),
        sample_R1=temp(expand(out("{{experiment_name}}/demux_reads_scatter/{{scatteritem}}/{sample_name}_R1.fastq.gz"), sample_name=get_sample_names())),
        sample_R2=temp(expand(out("{{experiment_name}}/demux_reads_scatter/{{scatteritem}}/{sample_name}_R2.fastq.gz"), sample_name=get_sample_names())),
    log:
        out("logs/step2_demultiplexing_reads/demultiplex_fastq_split_{experiment_name}_{scatteritem}.log"),
    threads: 1
    resources:
        mem_mb=1024 * 5,
    benchmark:
        out("benchmarks/{experiment_name}/demultiplex_fastq_split_{scatteritem}.txt")
    params:
        path_samples=config["path_samples"],
        path_barcodes=config["path_barcodes"],
        path_out=out("{experiment_name}/demux_reads_scatter/{scatteritem}"),
    conda:
        "../envs/sci-rocket.yaml",
    shell:
        """
        exec > "{log}" 2>&1
        set -euo pipefail
        $CONDA_PREFIX/bin/python {workflow.basedir}/scripts/demultiplexing/demux_rocket.py \
        --experiment_name {wildcards.experiment_name} \
        --samples {params.path_samples} \
        --barcodes {params.path_barcodes} \
        --r1 {input.R1} --r2 {input.R2} \
        --out {params.path_out}
        """


rule gather_demultiplexed_sequencing:
    input:
        discard_R1=expand(
            out("{{experiment_name}}/demux_reads_scatter/{scatteritem}/{{experiment_name}}_R1_discarded.fastq.gz"),
            scatteritem=get_fastq_split_scatteritems(),
        ),
        discard_R2=expand(
            out("{{experiment_name}}/demux_reads_scatter/{scatteritem}/{{experiment_name}}_R2_discarded.fastq.gz"),
            scatteritem=get_fastq_split_scatteritems(),
        ),
        discard_log=expand(
            out("{{experiment_name}}/demux_reads_scatter/{scatteritem}/log_{{experiment_name}}_discarded_reads.tsv.gz"),
            scatteritem=get_fastq_split_scatteritems(),
        ),
        qc_pickles=expand(
            out("{{experiment_name}}/demux_reads_scatter/{scatteritem}/{{experiment_name}}_qc.pickle"),
            scatteritem=get_fastq_split_scatteritems(),
        ),
        whitelist_p7=expand(
            out("{{experiment_name}}/demux_reads_scatter/{scatteritem}/{{experiment_name}}_whitelist_p7.txt"),
            scatteritem=get_fastq_split_scatteritems(),
        ),
        whitelist_p5=expand(
            out("{{experiment_name}}/demux_reads_scatter/{scatteritem}/{{experiment_name}}_whitelist_p5.txt"),
            scatteritem=get_fastq_split_scatteritems(),
        ),
        whitelist_ligation=expand(
            out("{{experiment_name}}/demux_reads_scatter/{scatteritem}/{{experiment_name}}_whitelist_ligation.txt"),
            scatteritem=get_fastq_split_scatteritems(),
        ),
        whitelist_rt=expand(
            out("{{experiment_name}}/demux_reads_scatter/{scatteritem}/{{experiment_name}}_whitelist_rt.txt"),
            scatteritem=get_fastq_split_scatteritems(),
        ),
    output:
        R1_discarded=out("{experiment_name}/demux_reads/{experiment_name}_R1_discarded.fastq.gz"),
        R2_discarded=out("{experiment_name}/demux_reads/{experiment_name}_R2_discarded.fastq.gz"),
        discarded_log=out("{experiment_name}/demux_reads/log_{experiment_name}_discarded_reads.tsv.gz"),
        qc=out("{experiment_name}/demux_reads/{experiment_name}_qc.pickle"),
        whitelist_p7=out("{experiment_name}/demux_reads/{experiment_name}_whitelist_p7.txt"),
        whitelist_p5=out("{experiment_name}/demux_reads/{experiment_name}_whitelist_p5.txt"),
        whitelist_ligation=out("{experiment_name}/demux_reads/{experiment_name}_whitelist_ligation.txt"),
        whitelist_rt=out("{experiment_name}/demux_reads/{experiment_name}_whitelist_rt.txt"),
    log:
        out("logs/step2_demultiplexing_reads/gather_demultiplexed_sequencing_{experiment_name}.log"),
    threads: 1
    resources:
        mem_mb=1024 * 10,
    benchmark:
        out("benchmarks/{experiment_name}/gather_demultiplexed_sequencing.txt")
    params:
        path_demux_scatter=out("{experiment_name}/demux_reads_scatter/"),
    conda:
        "../envs/sci-rocket.yaml",
    shell:
        """
        exec > "{log}" 2>&1
        set -euo pipefail
        # Combine pickles.
        $CONDA_PREFIX/bin/python {workflow.basedir}/scripts/demultiplexing/demux_gather.py --path_demux_scatter {params.path_demux_scatter} --path_out {output.qc}

        # Combine the sequencing-specific R1/R2 discarded reads and logs.
        cat {input.discard_R1} > {output.R1_discarded}
        cat {input.discard_R2} > {output.R2_discarded}
        printf "read_name\tp5\tp7\tligation\trt\tumi\tsample_name\n" | gzip -c > {output.discarded_log}
        cat {input.discard_log} >> {output.discarded_log}

        # Move the whitelist files.
        cat {input.whitelist_p7} | sort -u > {output.whitelist_p7}
        cat {input.whitelist_p5} | sort -u > {output.whitelist_p5}
        cat {input.whitelist_ligation} | sort -u > {output.whitelist_ligation}
        cat {input.whitelist_rt} | sort -u > {output.whitelist_rt}
        """

rule gather_demultiplexed_samples:
    input:
        R1_scatter=expand(
            out("{{experiment_name}}/demux_reads_scatter/{scatteritem}/{{sample_name}}_R1.fastq.gz"),
            scatteritem=get_fastq_split_scatteritems(),
        ),
        R2_scatter=expand(
            out("{{experiment_name}}/demux_reads_scatter/{scatteritem}/{{sample_name}}_R2.fastq.gz"),
            scatteritem=get_fastq_split_scatteritems(),
        ),
    output:
        R1=out("{experiment_name}/demux_reads/{sample_name}_R1.fastq.gz"),
        R2=out("{experiment_name}/demux_reads/{sample_name}_R2.fastq.gz"),
    log:
        out("logs/step2_demultiplexing_reads/gather_demultiplexed_samples_{experiment_name}_{sample_name}.log"),
    threads: 1
    resources:
        mem_mb=1024 * 10,
    benchmark:
        out("benchmarks/{experiment_name}/gather_demultiplexed_samples_{sample_name}.txt")
    shell:
        """
        exec > "{log}" 2>&1
        set -euo pipefail
        # Combine files.
        cat {input.R1_scatter} > {output.R1}
        cat {input.R2_scatter} > {output.R2}
        """
