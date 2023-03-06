# -------------------------------------------------------------------------------------
#
#
#    Script for processing snATAC-seq FASTQ files with Cell Ranger (10X)
#
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"


# -------------  RULES  ---------------

rule CR_cnt_ATAC:
    # Diminishing returns > 128Gs
    output: "../results/01CELLRANGER/{SAMPLE}.stamp" # Need .stamp and touch {output} in shell command as both SM and CR want to mkdir
    params: FASTQ_DIR=config["FASTQ_DIR"],
            REFERENCE=config["REFERENCE_ATAC"]
    log:    "../results/00LOGS/01CELLRANGER/{SAMPLE}.log"
    shell:
            """

            cellranger-atac count --id={wildcards.SAMPLE} \
            --fastqs={params.FASTQ_DIR} \
            --sample={wildcards.SAMPLE} \
            --reference={params.REFERENCE} \
            --localcores=32 \
            --localmem=128 2> {log}
            
            touch {output}

            """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
