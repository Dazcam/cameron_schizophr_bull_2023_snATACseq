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

rule mv_frag_files:
    input: "../results/01CELLRANGER/{SAMPLE}.stamp"
    output: "../results/04FRAGMENT_FILES/{SAMPLE}.stamp"
    log:    "../results/00LOGS/04FRAGMENT_FILES/{SAMPLE}.log"
    params: indir = "../results/01CELLRANGER/{SAMPLE}/outs/",
            outdir = "../results/04FRAGMENT_FILES/"
    shell:
            """ 

            cp {params.indir}fragments.tsv.gz {params.outdir}{wildcards.SAMPLE}.fragments.tsv.gz 
            cp {params.indir}fragments.tsv.gz.tbi {params.outdir}{wildcards.SAMPLE}.fragments.tsv.gz.tbi
            
            touch {output}
            """ 

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
