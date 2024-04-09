# -------------------------------------------------------------------------------------
#
#    QC and processing snATAC-seq data
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARMAS  ----------

configfile: "../config/config.yaml"

# -------------  RULES  ---------------

rule snATAC_seq_QC:
    input:     cellranger = expand("../results/01CELLRANGER/{SAMPLE}.stamp", SAMPLE = config['SAMPLES_ATAC']),
               markdown = "scripts/snATACseq_QC.Rmd"
    output:    "../results/03MARKDOWN/snATACseq_QC_{REGION}.html"
    params:    data_dir = "../results/01CELLRANGER/",
               archR_out_dir = "../results/02ARCHR/{REGION}",
               report_dir = "../results/03MARKDOWN/",
               report_file = "snATACseq_QC_{REGION}.html"
    envmodules: "libgit2/1.1.0", "R/4.2.0", "pandoc/2.7.3"
    log:       "../results/00LOGS/02ARCHR/snATAC_QC_{REGION}.log"
    shell:
               """
               Rscript --vanilla scripts/snATACseq_QC.R {wildcards.REGION} {params.data_dir} \
               {params.archR_out_dir} {input.markdown} {params.report_dir} {params.report_file} 2> {log}
               """

