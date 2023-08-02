# -------------------------------------------------------------------------------------
#
#
#    Script for running LDSR on snATAC-seq data
#
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"

# -------------  RULES  ---------------
rule ldsr_cond_lift_over:
    input:   mybed = "../results/05PEAKS/ZIFFRA_PEAKS/{ZIF_CELL_TYPE}.bed",
             chain_file = "../resources/liftover/hg38ToHg19.over.chain.gz",
             exe = "../resources/liftover/liftOver"
    output:  "../results/05PEAKS/ZIFFRA_PEAKS/{ZIF_CELL_TYPE}.hg19.bed"
    message: "Lifting {input.mybed} to hg19"
    log:     "../results/00LOGS/06LDSR/ZIFFRA_PEAKS/snATACseq.{ZIF_CELL_TYPE}.hg19_liftOver.log"
    params:  "../results/05PEAKS/ZIFFRA_PEAKS/{ZIF_CELL_TYPE}.hg38_unlifted.bed"
    shell:
             """

             {input.exe} {input.mybed} {input.chain_file} {output} {params} 2> {log}

             """

rule ldsr_cond_remove_MHC_from_peaks:
    input:   bed = "../results/05PEAKS/ZIFFRA_PEAKS/{ZIF_CELL_TYPE}.hg19.bed",
             MHC = "../resources/sheets/MHC.hg19.bed"
    output:  "../results/05PEAKS/ZIFFRA_PEAKS/{ZIF_CELL_TYPE}.hg19.noMHC.bed"
    message: "Removing MHC regions from {wildcards.ZIF_CELL_TYPE}"
    log:     "../results/00LOGS/06LDSR/ZIFFRA_PEAKS/{ZIF_CELL_TYPE}.hg19.noMHC.log"
    envmodules: "bedtools"
    shell:
        """
        
	bedtools subtract -a {input.bed} -b {input.MHC} -A > {output}

        """ 

rule ldsr_cond_make_annot:
    # Input can be bed file with gene boundaries or gene set with separate gene coord file
    input:   bed_file = "../results/05PEAKS/ZIFFRA_PEAKS/{ZIF_CELL_TYPE}.hg19.noMHC.bed",
             bim_file = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}.bim"
    output:  "../results/06LDSR/ZIFFRA_PEAKS/annotation_files/snATACseq.{ZIF_CELL_TYPE}.{CHR}.annot.gz"
    conda:   "../envs/ldsr.yml"
    message: "Creating annotation files for snATACseq: {wildcards.ZIF_CELL_TYPE}, Chr {wildcards.CHR}"
    log:     "../results/00LOGS/06LDSR/ZIFFRA_PEAKS/snATACseq.{ZIF_CELL_TYPE}.Chr{CHR}.log"
    shell:
             """
             
             python ../resources/ldsr/make_annot.py \
             --bed-file {input.bed_file} \
             --bimfile {input.bim_file} \
             --annot-file {output} 2> {log}
             
             """

rule ldsr_cond_ld_scores:
    input:   annot = "../results/06LDSR/ZIFFRA_PEAKS/annotation_files/snATACseq.{ZIF_CELL_TYPE}.{CHR}.annot.gz",
             bfile_folder = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink",
             snps_folder = "../resources/ldsr/reference_files/hapmap3_snps"
    output:  "../results/06LDSR/ZIFFRA_PEAKS/annotation_files/snATACseq.{ZIF_CELL_TYPE}.{CHR}.l2.ldscore.gz"
    conda:   "../envs/ldsr.yml"
    params:  bfile = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}",
             ldscores = "../results/06LDSR/ZIFFRA_PEAKS/annotation_files/snATACseq.{ZIF_CELL_TYPE}.{CHR}",
             snps = "../resources/ldsr/reference_files/w_hm3.snplist_rsIds"
    message: "Running LDSR Phase 3 for {wildcards.ZIF_CELL_TYPE}, CHR {wildcards.CHR}" 
    log:     "../results/00LOGS/06LDSR/ZIFFRA_PEAKS/snATACseq.{ZIF_CELL_TYPE}.Chr{CHR}_ldsc.log"
    shell:
        "python ../resources/ldsr/ldsc.py --thin-annot --l2 --bfile {params.bfile} --ld-wind-cm 1 "
        "--annot {input.annot} --out {params.ldscores} --print-snps {params.snps} 2> {log}"


rule ldsr_cond_stratified_baseline_v12:
    # Needs to be run after standard, non-conditional analyses - note we have union and Ziffra peaks in model
    input:   GWAS = "../results/06LDSR/GWAS_for_LDSR/{GWAS}_hg19_ldsc_ready.sumstats.gz",
             LDSR = expand("../results/06LDSR/annotation_files/snATACseq.{SLDSR_CELL_TYPE}.{CHR}.l2.ldscore.gz", SLDSR_CELL_TYPE = config["SLDSR_CELL_TYPES"], CHR = range(1,23)),
             COND = expand("../results/06LDSR/ZIFFRA_PEAKS/annotation_files/snATACseq.{COND_CELL_TYPE}.{CHR}.l2.ldscore.gz", COND_CELL_TYPE = config["COND_CELL_TYPES"], CHR = range(1,23)),
             UNION = expand("../results/06LDSR/annotation_files/snATACseq.union.{CHR}.l2.ldscore.gz", CHR = range(1,23))   
    output:  "../results/06LDSR/ZIFFRA_PEAKS/part_herit/baseline_v1.2/snATACseq.{SLDSR_CELL_TYPE}_vs_{COND_CELL_TYPE}.{GWAS}_baseline.v1.2.results"
    conda:   "../envs/ldsr.yml"
    params:  weights = "../resources/ldsr/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsr/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsr/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/06LDSR/annotation_files/snATACseq.{SLDSR_CELL_TYPE}.",
             cond_union = "../results/06LDSR/annotation_files/snATACseq.union.",
             cond_anns = "../results/06LDSR/ZIFFRA_PEAKS/annotation_files/snATACseq.{COND_CELL_TYPE}.",
             out_file = "../results/06LDSR/ZIFFRA_PEAKS/part_herit/baseline_v1.2/snATACseq.{SLDSR_CELL_TYPE}_vs_{COND_CELL_TYPE}.{GWAS}_baseline.v1.2"
    message: "Running conditional Prt Hrt with {wildcards.SLDSR_CELL_TYPE} vs. {wildcards.COND_CELL_TYPE} and {wildcards.GWAS} GWAS"
    log:     "../results/00LOGS/06LDSR/ZIFFRA_PEAKS/snATACseq.{SLDSR_CELL_TYPE}_vs_{COND_CELL_TYPE}.{GWAS}.baseline.v1.2_partHerit.log"
    shell:
             "python ../resources/ldsr/ldsc.py --h2 {input.GWAS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.cond_union},{params.cond_anns},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

rule ldsr_cond_stratified_summary:
    # Careful here: L2_1 if no covariates, L2_2 if one covariate last number increasing for each additional covariate
    input:   expand("../results/06LDSR/ZIFFRA_PEAKS/part_herit/baseline_v1.2/snATACseq.{SLDSR_CELL_TYPE}_vs_{COND_CELL_TYPE}.{GWAS}_baseline.v1.2.results", SLDSR_CELL_TYPE = config["SLDSR_CELL_TYPES"], COND_CELL_TYPE = config["COND_CELL_TYPES"], GWAS = config['LDSR_GWAS'])
    output:  "../results/06LDSR/ZIFFRA_PEAKS/part_herit/baseline_v1.2/snATACseq_LDSR_{GWAS}_baseline.v1.2_summary.tsv"
    message: "Creating summary file for {wildcards.GWAS} GWAS"
    params:  dir = "../results/06LDSR/ZIFFRA_PEAKS/part_herit/baseline_v1.2/",
             cell_types = "../resources/sheets/GE_cond_celltypes.tsv"
    log:     "../results/00LOGS/06LDSR/ZIFFRA_PEAKS/snRNAseq.{GWAS}_baseline.v1.2_partHerit.summary.log"
    shell:
             """
             
             head -1 {params.dir}snATACseq.CGE_vs_dlEN_ziffra_macs2peaks.SCZ_baseline.v1.2.results > {output}
             File={params.cell_types}
             Lines=$(cat $File)
             for Line in $Lines
             do
             grep L2_3 {params.dir}snATACseq."$Line".{wildcards.GWAS}_baseline.v1.2.results | sed "s/L2_3/$Line/g" >> {output} 2> {log}
             done

             """
