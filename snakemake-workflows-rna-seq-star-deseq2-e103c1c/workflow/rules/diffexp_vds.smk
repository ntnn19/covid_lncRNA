def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]


rule limma:
    input:
        counts="results/counts/all.tsv",
    output:
        table=report("results/diffexp/{contrast}.diffexp.tsv", "../report/diffexp.rst"),
    params:
        contrast=get_contrast,
    conda:
        "../envs/limma.yaml"
    log:
        "logs/limma/limma_{contrast}.log",
    script:
        "../scripts/limma_claude.R"

