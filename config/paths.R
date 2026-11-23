# All file locations used by the project, in one place.
# Scripts are run from the project root; every path below is relative to it,
# except the external databases at the bottom (used only by scripts/00_data_generation/).

paths <- list(
    # --- raw inputs (copied or symlinked from the original projects) -----------------------------
    quant_alignment = 'data/raw/quantification/alignment',   # BWA-mem2 (nextflow quantification), TPM/Numreads/gene-characteristics
    quant_mapping   = 'data/raw/quantification/mapping',     # Salmon (nextflow quantification)
    quant_stats     = 'data/raw/quantification/statistics',  # Pelagomonas Tara TPM table, Carradec suppl. tables
    annotation      = 'data/raw/genomic/functional_annotation',
    gene_data       = 'data/raw/gene_data',
    env             = 'data/raw/env',
    placement       = 'data/raw/placement',                  # symlinks to the jplaces / multiplicities

    # --- outputs of the data-generation steps (copied, not recomputed) -------------------------
    intermediate    = 'data/intermediate',

    # --- outputs ---------------------------------------------------------------------------------
    figures         = 'results/figures',
    tables          = 'results/tables'
)

# External databases (only for scripts/00_data_generation/, never read by the figure scripts)
external <- list(
    matou     = '/home/aauladell/datasets/common_databases/MATOU',
    eukprot   = '/home/aauladell/datasets/common_databases/2020_eukprot_v1',
    env_root  = '/home/aauladell/datasets/environmental_data_from_datasets',
    # --- only for scripts/00_data_generation/00_upstream_pipelines/ ---
    matou_translations = '/home/aauladell/datasets/common_databases/MATOU/matou_alltranslations_filt10aa.faa',
    eukprot_proteins   = '/scratch/data4/aauladell/eukprot_v01/proteins',     # EukProt v01 proteomes (hmmsearch db)
    eggnog_db          = '/scratch/data1/aauladell/databases/eggnog',
    hk_taxonomies_richter = '/scratch/data1/drichter/metatranscriptome_read_placement/jurgen_321_genes/sequences/31_gene-specific_taxonomies_23',
    # original projects, kept for provenance
    gene_explorer        = '/scratch/data1/aauladell/projects_big/gene_explorer',
    delving_florenciella = '/scratch/data3/aauladell/project_big/delving_into_florenciella',
    phylotranscriptomics = '/scratch/data3/aauladell/project_big/phylotranscriptomics'
)

# Small helper: build a path inside one of the entries above
p <- function(key, ...) file.path(paths[[key]], ...)
