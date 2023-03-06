# -------------------------------------------------------------------------------------
#
#    QC and processing snATAC-seq data
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARMAS  ----------

configfile: "../config/config.yaml"

# -------------  RULES  ---------------

rule snATAC_seq_QC:
    input:     markdown = "scripts/snATACseq_QC.Rmd"
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

rule snATAC_remove_batch_effects:
    input:  markdown = "scripts/snATACseq_remove_batch_effects.Rmd",
            html = "../results/03MARKDOWN/snATACseq_QC_{REGION}.html", # Needed for rule order
    output: "../results/03MARKDOWN/snATACseq_remove_batch_effects_{REGION}.html"
    params: data_dir = "../results/01CELLRANGER/",
            archR_out_dir = "../results/02ARCHR/{REGION}",
            report_dir = "../results/03MARKDOWN/",
            report_file = "snATACseq_remove_batch_effects_{REGION}.html",
    envmodules: "libgit2/1.1.0", "R/4.2.0", "pandoc/2.7.3"
    log:    "../results/00LOGS/02ARCHR/snATAC_remove_batch_effects_{REGION}.log"
    shell:
            """

            Rscript --vanilla scripts/snATACseq_remove_batch_effects.R {wildcards.REGION} {params.data_dir} \
            {params.archR_out_dir} {input.markdown} {params.report_dir} {params.report_file} 2> {log}
            
            """
