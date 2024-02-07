# -------------------------------------------------------------------------------------
#
#
#    Script for running SLDSC on snATAC-seq data unique GE peaks
#
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"

rule LDSR_unique_remove_MHC_from_peaks:
    input:   bed = "../results/05PEAKS/UNIQUE_PEAKS/{CELL_TYPE}.vs.MSCOFF_unique_peaks.hg19.ext250bp.bed",
             MHC = "../resources/sheets/MHC.hg19.bed"
    output:  "../results/05PEAKS/{CELL_TYPE}.vs.MSCOFF_unique_peaks.hg19.ext250bp.noMHC.bed"
    message: "Removing MHC regions from {wildcards.CELL_TYPE}.vs.MSCOFF_unique_peaks"
    log:     "../results/00LOGS/06LDSR/unique_peaks/{CELL_TYPE}.vs.MSCOFF_unique_peaks.hg19.noMHC.log"
    envmodules: "bedtools"
    shell:
        """
        
	bedtools subtract -a {input.bed} -b {input.MHC} -A > {output}

        """ 

rule LDSR_unique_make_annot:
    # Input can be bed file with gene boundaries or gene set with separate gene coord file
    input:   bed_file = "../results/05PEAKS/{CELL_TYPE}.vs.MSCOFF_unique_peaks.hg19.ext250bp.noMHC.bed",
             bim_file = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}.bim"
    output:  "../results/06LDSR/annotation_files/unique_peaks/snATACseq.{CELL_TYPE}.vs.MSCOFF_unique_peaks.{CHR}.annot.gz"
    conda:   "../envs/ldsr.yml"
    message: "Creating annotation files for snATACseq: {wildcards.CELL_TYPE}, Chr {wildcards.CHR}"
    log:     "../results/00LOGS/06LDSR/unique_peaks/snATACseq.{CELL_TYPE}.vs.MSCOFF_unique_peaks.Chr{CHR}.log"
    shell:
             """
             
             python ../resources/ldsr/make_annot.py \
             --bed-file {input.bed_file} \
             --windowsize 0 \
             --bimfile {input.bim_file} \
             --annot-file {output} 2> {log}             
             
             """


rule LDSR_unique_ld_scores:
    input:   annot = "../results/06LDSR/annotation_files/unique_peaks/snATACseq.{CELL_TYPE}.vs.MSCOFF_unique_peaks.{CHR}.annot.gz",
             bfile_folder = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink",
             snps_folder = "../resources/ldsr/reference_files/hapmap3_snps"
    output:  "../results/06LDSR/annotation_files/unique_peaks/snATACseq.{CELL_TYPE}.vs.MSCOFF_unique_peaks.{CHR}.l2.ldscore.gz"
    conda:   "../envs/ldsr.yml"
    params:  bfile = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}",
             ldscores = "../results/06LDSR/annotation_files/unique_peaks/snATACseq.{CELL_TYPE}.vs.MSCOFF_unique_peaks.{CHR}",
             snps = "../resources/ldsr/reference_files/w_hm3.snplist_rsIds"
    message: "Running LDSR Phase 3 for {wildcards.CELL_TYPE}.vs.MSCOFF_unique_peaks, CHR {wildcards.CHR}" 
    log:     "../results/00LOGS/06LDSR/unique_peaks/snATACseq.{CELL_TYPE}.vs.MSCOFF_unique_peaks.Chr{CHR}_ldsc.log"
    shell:
        "python ../resources/ldsr/ldsc.py --thin-annot --l2 --bfile {params.bfile} --ld-wind-cm 1 "
        "--annot {input.annot} --out {params.ldscores} --print-snps {params.snps} 2> {log}"


rule LDSR_unique_stratified_baseline_v12:
    input:   GWAS = "../results/06LDSR/GWAS_for_LDSR/{GWAS}_hg19_LDSR_ready.sumstats.gz",
             LDSR = expand("../results/06LDSR/annotation_files/unique_peaks/snATACseq.{CELL_TYPE}.vs.MSCOFF_unique_peaks.{CHR}.l2.ldscore.gz", CELL_TYPE = config["UNIQUE_CELL_TYPES"], CHR = range(1,23))
    output:  "../results/06LDSR/unique_peaks/part_herit/baseline_v1.2/snATACseq.{CELL_TYPE}.{GWAS}_baseline.v1.2.results"
    conda:   "../envs/ldsr.yml"
    params:  weights = "../resources/ldsr/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsr/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsr/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/06LDSR/annotation_files/unique_peaks/snATACseq.{CELL_TYPE}.vs.MSCOFF_unique_peaks.",
             union_anns = "../results/06LDSR/annotation_files/snATACseq.union.",
             out_file = "../results/06LDSR/unique_peaks/part_herit/baseline_v1.2/snATACseq.{CELL_TYPE}.{GWAS}_baseline.v1.2"
    message: "Running Prt Hrt with {wildcards.CELL_TYPE}.vs.MSCOFF_unique_peaks and {wildcards.GWAS} GWAS"
    log:     "../results/00LOGS/06LDSR/unique_peaks/snATACseq.{CELL_TYPE}.{GWAS}.baseline.v1.2_partHerit.log"
    shell:
             "python ../resources/ldsr/ldsc.py --h2 {input.GWAS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.union_anns},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

rule LDSR_unique_stratified_summary:
    # Careful here: L2_1 if no covariates, L2_2 if one covariate last number increasing for each additional covariate
    input:   expand("../results/06LDSR/unique_peaks/part_herit/baseline_v1.2/snATACseq.{CELL_TYPE}.{GWAS}_baseline.v1.2.results", CELL_TYPE = config["UNIQUE_CELL_TYPES"], GWAS = config["LDSR_GWAS"])
    output:  "../results/06LDSR/unique_peaks/part_herit/baseline_v1.2/snATACseq_LDSR_{GWAS}_baseline.v1.2_summary.tsv"
    message: "Creating summary file for {wildcards.GWAS} GWAS"
    params:  dir = "../results/06LDSR/unique_peaks/part_herit/baseline_v1.2/",
             cell_types = "../resources/sheets/GE_celltypes_unique.tsv"
    log:     "../results/00LOGS/06LDSR/unique_peaks/snRNAseq.{GWAS}_baseline.v1.2_partHerit.summary.log"
    shell:
             """
             
             head -1 {params.dir}snATACseq.LGE.SCZ_baseline.v1.2.results > {output}
             File={params.cell_types}
             Lines=$(cat $File)
             for Line in $Lines
             do
             grep L2_1 {params.dir}snATACseq."$Line".{wildcards.GWAS}_baseline.v1.2.results | sed "s/L2_1/$Line/g" >> {output} 2> {log}
             done

             """

