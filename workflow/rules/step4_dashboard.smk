#   * Generate QC-dashboard (sci-dash). *
#
#   1. sci-dash:                Generates the QC-dashboard for a given sequencing run.
#############

# Get the samples for a given sequencing run (and sci-dash).
def getsamples_sequencing(wildcards):
    df = samples_unique[samples_unique["experiment_name"] == wildcards.experiment_name]
    files = [out(f"{s.experiment_name}/alignment/{s.sample_name}_{s.species}_Aligned.sortedByCoord.out.bam.bai") for s in df.itertuples(index=False, name='Row')]

    return files


def get_preliminary_umaps(wildcards):
    df = samples_unique[samples_unique["experiment_name"] == wildcards.experiment_name]
    files = [
        out(f"{s.experiment_name}/preliminary_umap/{s.sample_name}_{s.species}_umap.json")
        for s in df.itertuples(index=False, name='Row')
    ]

    return files


rule scanpy_preliminary_umap:
    input:
        dir_solo=out("{experiment_name}/alignment/{sample_name}_{species}_Solo.out"),
    output:
        out("{experiment_name}/preliminary_umap/{sample_name}_{species}_umap.json"),
    log:
        out("logs/step4_dashboard/scanpy_preliminary_umap_{experiment_name}_{sample_name}_{species}.log"),
    threads: 4
    resources:
        mem_mb=1024 * 8,
    benchmark:
        out("benchmarks/{experiment_name}/scanpy_preliminary_umap_{sample_name}_{species}.txt")
    params:
        path_star=out("{experiment_name}/alignment/"),
        solo_features=lambda w: get_star_solo_features(),
        random_seed=0,
    conda:
        "../envs/scanpy-umap.yaml",
    shell:
        """
        exec > "{log}" 2>&1
        set -euo pipefail
        $CONDA_PREFIX/bin/python3 {workflow.basedir}/scripts/demultiplexing/scanpy_preliminary_umap.py \
        --path_star {params.path_star} \
        --sample {wildcards.sample_name} \
        --path_out {output} \
        --random_seed {params.random_seed} \
        --solo_features {params.solo_features:q}
        """

rule sci_dash:
    input:
        lambda w: getsamples_sequencing(w),
        qc=out("{experiment_name}/demux_reads/{experiment_name}_qc.pickle"),
        umaps=lambda w: get_preliminary_umaps(w),
    output:
        dash_folder=directory(out("{experiment_name}/sci-dash/")),
        dash_json=out("{experiment_name}/sci-dash/js/qc_data.js"),
        dash_umap_json=out("{experiment_name}/sci-dash/js/umap_data.js"),
    log:
        out("logs/step4_dashboard/sci_dash_{experiment_name}.log"),
    threads: 1
    resources:
        mem_mb=1024 * 10,
    benchmark:
        out("benchmarks/{experiment_name}/sci_dash.txt")
    params:
        path_star=out("{experiment_name}/alignment/"),
        # Optional hashing output. 
        metrics_hashing=out("{experiment_name}/hashing/{experiment_name}_hashing_metrics.tsv"),
        benchmarks_folder=out("benchmarks/{experiment_name}"),
        solo_features=lambda w: get_star_solo_features(),
    conda:
        "../envs/sci-rocket.yaml",
    shell:
        """
        exec > "{log}" 2>&1
        set -euo pipefail
        # Generate the sci-dashboard report.
        cp -R {workflow.basedir}/scirocket-dash/* {output.dash_folder}

        # Combine the sample-specific QC and STARSolo metrics.
        $CONDA_PREFIX/bin/python3 {workflow.basedir}/scripts/demultiplexing/demux_dash.py \
        --path_out {output.dash_json} \
        --path_pickle {input.qc} \
        --path_star {params.path_star} \
        --path_hashing {params.metrics_hashing} \
        --path_benchmarks {params.benchmarks_folder} \
        --solo_features {params.solo_features:q}

        $CONDA_PREFIX/bin/python3 {workflow.basedir}/scripts/demultiplexing/build_umap_data_js.py \
        --path_out {output.dash_umap_json} \
        --path_inputs {input.umaps}

        # Remove all empty (leftover) folders.
        find {dir_output}/{wildcards.experiment_name}/ -empty -type d -delete
        """
