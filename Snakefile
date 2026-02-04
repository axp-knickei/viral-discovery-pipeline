configfile: "config/default.yaml"

rule all:
    input:
        "outputs/clustering/abundance.tsv"

rule qc:
    input:
        reads=lambda wildcards: ["data/sample_R1.fastq", "data/sample_R2.fastq"]
    output:
        report="outputs/qc/qc_report.json",
        cleaned=["outputs/qc/sample_R1.clean.fastq", "outputs/qc/sample_R2.clean.fastq"]
    conda:
        "workflow/envs/qc.yaml"
    shell:
        "mkdir -p outputs/qc && touch {output.report} && touch {output.cleaned}"

rule host_removal:
    input:
        reads=rules.qc.output.cleaned
    output:
        report="outputs/host/host_removal_report.json",
        filtered=["outputs/host/sample_R1.host_filtered.fastq", "outputs/host/sample_R2.host_filtered.fastq"]
    conda:
        "workflow/envs/host_removal.yaml"
    shell:
        "mkdir -p outputs/host && touch {output.report} && touch {output.filtered}"

rule assembly:
    input:
        reads=rules.host_removal.output.filtered
    output:
        contigs="outputs/assembly/contigs.fasta",
        report="outputs/assembly/assembly_report.json"
    conda:
        "workflow/envs/assembly.yaml"
    shell:
        "mkdir -p outputs/assembly && echo '>contig' > {output.contigs} && echo '{{}}' > {output.report}"

rule viral_detection:
    input:
        contigs=rules.assembly.output.contigs
    output:
        contigs="outputs/viral/viral_contigs.fasta",
        report="outputs/viral/viral_detection_report.json"
    conda:
        "workflow/envs/viral_detection.yaml"
    shell:
        "mkdir -p outputs/viral && echo '>viral' > {output.contigs} && echo '{{}}' > {output.report}"

rule clustering:
    input:
        contigs=rules.viral_detection.output.contigs
    output:
        clusters="outputs/clustering/clusters.tsv",
        abundance="outputs/clustering/abundance.tsv"
    conda:
        "workflow/envs/clustering.yaml"
    shell:
        "mkdir -p outputs/clustering && echo 'cluster\tcontig' > {output.clusters} && echo 'cluster\tabundance' > {output.abundance}"
