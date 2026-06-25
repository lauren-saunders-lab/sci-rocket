#   * Perform pre-processing and alignment of sci-seq reads. *
#
#   1. trim_fastp:                      Trimming adapters and low-quality reads.
#   2. generate_index_STAR:             Generating (or symlinking) STAR indexes.
#   3. starSolo_align:                  Aligning reads with STARsolo and adjusting the cellular barcode scheming to resemble sci-seq scheme.
#   4. sambamba_index:                  Indexing BAM files.
#############

import math
import pickle


rule trim_fastp:
    input:
        R1=out("{experiment_name}/demux_reads/{sample_name}_R1.fastq.gz"),
        R2=out("{experiment_name}/demux_reads/{sample_name}_R2.fastq.gz"),
    output:
        R1=temp(out("{experiment_name}/fastp/{sample_name}_R1.fastq.gz")),
        R2=temp(out("{experiment_name}/fastp/{sample_name}_R2.fastq.gz")),
        html=out("{experiment_name}/fastp/{sample_name}.html"),
        json=out("{experiment_name}/fastp/{sample_name}.json"),
    log:
        out("logs/step3_alignment/fastp_{experiment_name}_{sample_name}.log"),
    threads: 10
    resources:
        mem_mb=1024 * 4,
    benchmark:
        out("benchmarks/{experiment_name}/trim_fastp_{sample_name}.txt")
    params:
        extra=config["settings"]["fastp"],
    conda:
        "../envs/sci-rocket.yaml",
    shell:
        """
        exec > "{log}" 2>&1
        set -euo pipefail
        fastp {params.extra} --html {output.html} --json {output.json} --thread {threads} --in1 {input.R1} --in2 {input.R2} --out1 {output.R1} --out2 {output.R2}
        """


# Retrieve the expected no. of cells for a given (demultiplexed) sample.
def get_expected_cells(wildcards):
    x = samples_unique[samples_unique["sample_name"] == wildcards.sample_name]
    return x["n_expected_cells"].values[0]

def get_star_solo_features():
    configured = config["settings"].get("star_solo_features", "").strip()
    return configured if configured else "GeneFull_Ex50pAS"

def get_species_qc_pickles_for_index(species):
    if config["species"][species]["star_index"]:
        return []
    experiments = sorted(samples_unique.query("species == @species")["experiment_name"].unique())
    return [out(f"{experiment_name}/demux_reads/{experiment_name}_qc.pickle") for experiment_name in experiments]

def get_star_align_mem_mb(wildcards, input, attempt=1):
    # default memory requirement of 60 GB, sufficient for most samples with up to hundreds of millions of reads.
    min_mem_mb = 60 * 1024
    # for very large numbers (billions) of reads, the solo counting step in STARsolo dominates the memory usage.
    # with the default STARsolo settings in the pipeline, each mapped read requires ~17 bytes.
    # assuming 1 mapped read per input read and adding some buffer we require 20 bytes per read:
    bytes_per_read = 20
    # in case the step fails, we will increase the memory by a factor of 1.5 for each retry attempt:
    retry_multiplier = 1.5
    attempt = max(1, int(attempt))
    configured = config["settings"].get("star_align_mem_mb", "")

    if configured not in ("", None):
        base_mem_mb = int(configured)
        if base_mem_mb <= 0:
            raise ValueError("settings.star_align_mem_mb must be > 0")
    else:
        try:
            with open(input.demux_qc, "rb") as handle:
                qc = pickle.load(handle)
            n_reads = int(qc["sample_success"][wildcards.sample_name]["n_pairs_success"])
        except (FileNotFoundError, KeyError, ValueError, TypeError, pickle.PickleError):
            base_mem_mb = min_mem_mb
        else:
            estimated_mem_mb = math.ceil((n_reads * bytes_per_read) / (1024 * 1024))
            base_mem_mb = max(min_mem_mb, estimated_mem_mb)

    return math.ceil(base_mem_mb * (retry_multiplier ** (attempt - 1)))


rule generate_index_STAR:
    input:
        qc=lambda w: get_species_qc_pickles_for_index(w.species),
    output:
        directory(out("resources/index_star/{species}")),
    log:
        out("logs/step3_alignment/generate_index_STAR_{species}.log"),
    threads: 20
    resources:
        mem_mb=1024 * 50,
    benchmark:
        out("benchmarks/generate_index_STAR_{species}.txt")
    params:
        fasta=lambda w: config["species"][w.species]["genome"],
        gtf=lambda w: config["species"][w.species]["genome_gtf"],
        star_index=lambda w: config["species"][w.species]["star_index"],
        star_index_extra=config["settings"]["star_index"],
        star_overhang_script=f"{workflow.basedir}/scripts/demultiplexing/get_star_index_overhang.py",
    conda:
        "../envs/sci-rocket.yaml",
    shell:
        """
        exec > "{log}" 2>&1
        set -euo pipefail
        # Check if STAR_index is given. If not, generate it.
        if [ -n "{params.star_index}" ]; then
            ln -s "$(realpath "{params.star_index}")" "{output}"
        else
            ARGS=$(python "{params.star_overhang_script}" --star-index-extra "{params.star_index_extra}" {input.qc})
            STAR $ARGS --runThreadN {threads} --runMode genomeGenerate --genomeFastaFiles {params.fasta} --genomeDir {output} --sjdbGTFfile {params.gtf}
        fi
        """


rule starSolo_align:
    input:
        R1=out("{experiment_name}/fastp/{sample_name}_R1.fastq.gz"),
        R2=out("{experiment_name}/fastp/{sample_name}_R2.fastq.gz"),
        demux_qc=out("{experiment_name}/demux_reads/{experiment_name}_qc.pickle"),
        index=out("resources/index_star/{species}/"),
        whitelist_p7=out("{experiment_name}/demux_reads/{experiment_name}_whitelist_p7.txt"),
        whitelist_p5=out("{experiment_name}/demux_reads/{experiment_name}_whitelist_p5.txt"),
        whitelist_ligation=out("{experiment_name}/demux_reads/{experiment_name}_whitelist_ligation.txt"),
        whitelist_rt=out("{experiment_name}/demux_reads/{experiment_name}_whitelist_rt.txt"),
    output:
        bam=out("{experiment_name}/alignment/{sample_name}_{species}_Aligned.sortedByCoord.out.bam"),
        sj=out("{experiment_name}/alignment/{sample_name}_{species}_SJ.out.tab"),
        log1=out("{experiment_name}/alignment/{sample_name}_{species}_Log.final.out"),
        log2=out("{experiment_name}/alignment/{sample_name}_{species}_Log.out"),
        log3=out("{experiment_name}/alignment/{sample_name}_{species}_Log.progress.out"),
        dir_tmp=temp(
            directory(out("{experiment_name}/alignment/{sample_name}_{species}__STARtmp"))
        ),
        dir_solo=directory(
            out("{experiment_name}/alignment/{sample_name}_{species}_Solo.out")
        ),
    log:
        out("logs/step3_alignment/star_align_{experiment_name}_{sample_name}_{species}.log"),
    threads: 30
    resources:
        mem_mb=get_star_align_mem_mb,
    benchmark:
        out("benchmarks/{experiment_name}/starSolo_align_{sample_name}_{species}.txt")
    params:
        sampleName=out("{experiment_name}/alignment/{sample_name}_{species}_"),
        extra=config["settings"]["star"],
        path_barcodes=config["path_barcodes"],
        n_expected_cells=lambda w: get_expected_cells(w),
        solo_features=lambda w: get_star_solo_features(),
    conda:
        "../envs/sci-rocket.yaml",
    shell:
        """
        exec > "{log}" 2>&1
        set -euo pipefail

        LOCAL_STAR="{resources.tmpdir}/star_{wildcards.sample_name}_{wildcards.species}"
        mkdir -p "$LOCAL_STAR"
        mkdir -p "$(dirname {params.sampleName})"

        STAR {params.extra} --genomeDir {input.index} --runThreadN {threads} \
        --readFilesIn {input.R2} {input.R1} --readFilesCommand zcat \
        --soloFeatures {params.solo_features} \
        --soloType CB_UMI_Complex \
        --soloCBmatchWLtype Exact \
        --soloCellReadStats Standard \
        --soloCBposition 0_0_0_9 0_10_0_19 0_20_0_29 0_30_0_39 --soloUMIposition 0_40_0_47 \
        --soloCBwhitelist {input.whitelist_p7} {input.whitelist_p5} {input.whitelist_ligation} {input.whitelist_rt} \
        --soloCellFilter CellRanger2.2 {params.n_expected_cells} 0.99 10 \
        --outTmpDir "$LOCAL_STAR/__STARtmp" \
        --outTmpKeep None \
        --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$LOCAL_STAR/"

        mv "$LOCAL_STAR/Aligned.sortedByCoord.out.bam" {output.bam}
        mv "$LOCAL_STAR/SJ.out.tab" {output.sj}
        mv "$LOCAL_STAR/Log.final.out" {output.log1}
        mv "$LOCAL_STAR/Log.out" {output.log2}
        mv "$LOCAL_STAR/Log.progress.out" {output.log3}
        cp -r "$LOCAL_STAR/Solo.out" {output.dir_solo}
        mkdir -p {output.dir_tmp}
        rm -rf "$LOCAL_STAR"

        $CONDA_PREFIX/bin/python3 {workflow.basedir}/scripts/demultiplexing/STARSolo_convertBarcodes.py --solo_out_dir {output.dir_solo} --features {params.solo_features} --barcodes {params.path_barcodes}
        """


rule sambamba_index:
    input:
        out("{experiment_name}/alignment/{sample_name}_{species}_Aligned.sortedByCoord.out.bam"),
    output:
        out("{experiment_name}/alignment/{sample_name}_{species}_Aligned.sortedByCoord.out.bam.bai"),
    log:
        out("logs/step3_alignment/sambamba_index_{experiment_name}_{sample_name}_{species}.log"),
    threads: 8
    resources:
        mem_mb=1024 * 2,
    benchmark:
        out("benchmarks/{experiment_name}/sambamba_index_{sample_name}_{species}.txt")
    conda:
        "../envs/sci-rocket.yaml",
    shell:
        """
        exec > "{log}" 2>&1
        set -euo pipefail
        sambamba index -t {threads} {input} {output}
        """
