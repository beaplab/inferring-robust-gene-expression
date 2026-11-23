# Wrappers around gappa (https://github.com/lczech/gappa) used by the data-generation scripts of
# Figs 2 and 7: per-query multiplicities of a jplace, new jplaces with other multiplicities (ratios),
# placement factorization and per-clade sums.
# Only used by scripts/00_data_generation/; needs gappa on the PATH (v0.9.0 was used).

# Run a shell command, stop if it fails
run_command <- function(command) {
    message('$ ', str_trunc(command, 300))
    status <- system(command)
    if (status != 0) stop('command failed: ', command)
    invisible(status)
}

# Low/high category of an environmental variable: "low_<variable>" below the threshold, "high_<variable>"
# at or above it, NA without a value
env_category <- function(values, variable, threshold) {
    levels <- str_c(c('low_', 'high_'), variable)
    factor(case_when(values < threshold ~ levels[1],
                     values >= threshold ~ levels[2]), levels = levels)
}

# Per-query multiplicities of a jplace: tibble(jplace_sample, name, count)
# (gappa edit multiplicity --write-multiplicity-file, in a temporary folder of its own)
gappa_multiplicities <- function(jplace) {
    prefix <- str_remove(basename(jplace), '\\.jplace(\\.gz)?$')
    outdir <- tempfile('multiplicity')
    dir.create(outdir)
    on.exit(unlink(outdir, recursive = TRUE))
    run_command(str_glue('gappa edit multiplicity --jplace-path {jplace} --out-dir {outdir} ',
                         '--write-multiplicity-file --file-prefix {prefix} --threads 1 ',
                         '--allow-file-overwriting > /dev/null 2>&1'))
    read_tsv(file.path(outdir, str_c(prefix, 'multiplicities.csv')),
             col_names = c('jplace_sample', 'name', 'count'), col_types = 'ccd', lazy = FALSE)
}

# Copy of a jplace with new multiplicities, written to out_dir with the same name (gzipped).
# values: two columns (query name, new multiplicity). Queries not listed keep their multiplicity, so
# every query of the jplace has to be in `values` (0 for those without a value).
write_jplace_multiplicities <- function(jplace, values, out_dir) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    tmpfile <- tempfile(fileext = '.tsv')
    on.exit(unlink(tmpfile))
    write_tsv(values, tmpfile, col_names = FALSE)
    run_command(str_glue('gappa edit multiplicity --jplace-path {jplace} --out-dir {out_dir} ',
                         '--multiplicity-file {tmpfile} --threads 1 --compress ',
                         '--allow-file-overwriting > /dev/null 2>&1'))
}

# Pseudo-count for the placement factorization of ratio jplaces.
# gappa adds --pseudo-count-summand-all (default 0.65) to every edge mass before computing the balances.
# Balances are log-ratios, so scaling all masses by k is the same as dividing the pseudo-count by k:
# the ratios sum to much more per sample than the metaT reads, so the default would barely smooth the
# zeros. It is scaled so that it weighs as much as in a factorization of raw metaT counts:
#   0.65 * median(per-sample ratio total) / median(per-sample metaT read total)
factorization_pseudocount <- function(ratio_totals, metat_totals) {
    0.65 * median(ratio_totals) / median(metat_totals, na.rm = TRUE)
}

# gappa analyze placement-factorization of `jplaces` against the column `catcol` of `metadata`
# (a data frame with `Sample Code` and catcol, one row per jplace, in the same order).
# Writes <out_dir>/<prefix>* (factors_tree.newick, factor_taxa.csv, ...) and the log.
gappa_placement_factorization <- function(jplaces, metadata, catcol, out_dir, prefix, pseudocount,
                                          threads = 8) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    stopifnot(nrow(metadata) == length(jplaces))
    write_tsv(select(metadata, `Sample Code`, all_of(catcol)), file.path(out_dir, 'metadata_factorization.tsv'))
    writeLines(jplaces, file.path(out_dir, 'jplace_list.txt'))
    run_command(str_glue('gappa analyze placement-factorization --jplace-path $(cat {out_dir}/jplace_list.txt) ',
                         '--write-newick-tree --write-svg-tree ',
                         '--metadata-table-file {out_dir}/metadata_factorization.tsv ',
                         '--metadata-separator-char tab --metadata-select-columns {catcol} ',
                         '--pseudo-count-summand-all {format(pseudocount, scientific = FALSE, digits = 6)} ',
                         '--threads {threads} --allow-file-overwriting --out-dir {out_dir}/ --file-prefix "{prefix}" ',
                         '> {out_dir}/placement-factorization.log 2>&1'))
}

# Clades of a placement factorization: each reference sequence goes to the highest-numbered factor
# whose non-root side contains it. Written as a gappa clade list (tip label, clade_<factor>).
factorization_clades <- function(factor_taxa_file, out_file) {
    read_tsv(factor_taxa_file, show_col_types = FALSE) %>%
        filter(RootSide == 0) %>%
        # reference sequences only, removing internal nodes
        filter(str_detect(Taxon, 'EP|[A-Z,a-z]')) %>%
        mutate(Factor = str_c('clade_', Factor)) %>%
        group_by(Taxon) %>%
        slice_max(order_by = Factor, n = 1) %>%
        select(Taxon, Factor) %>%
        write_tsv(out_file, col_names = FALSE)
}

# Splits every jplace into one jplace per clade of `clade_file` (gappa edit extract; queries outside
# every clade go to "basal", those spanning several to "uncertain") and writes the per-query
# multiplicities of each into mult_dir. Returns mult_dir.
gappa_clade_multiplicities <- function(jplaces, clade_file, jplace_dir, mult_dir, threads = 8) {
    unlink(c(jplace_dir, mult_dir), recursive = TRUE)
    dir.create(jplace_dir, recursive = TRUE)
    dir.create(mult_dir, recursive = TRUE)
    list_file <- tempfile(fileext = '.txt')
    on.exit(unlink(list_file))
    writeLines(jplaces, list_file)
    run_command(str_glue("cat {list_file} | xargs -P {threads} -I{{}} sh -c ",
                         "'gappa edit extract --threads 1 --jplace-path {{}} --clade-list-file {clade_file} ",
                         "--samples-out-dir {jplace_dir} --samples-file-prefix $(basename {{}} .jplace.gz) ",
                         "--allow-file-overwriting > /dev/null'"))
    run_command(str_glue("ls {jplace_dir}/* | xargs -P {threads} -I{{}} sh -c ",
                         "'gappa edit multiplicity --threads 1 --jplace-path {{}} --write-multiplicity-file ",
                         "--out-dir {mult_dir} --file-prefix $(basename {{}} .jplace) ",
                         "--allow-file-overwriting > /dev/null'"))
    invisible(mult_dir)
}

# Per-unigene multiplicities from the files of gappa_clade_multiplicities(), in one table:
# sample (sample code without "filter<clade>"), clade, geneid, count. The counts are kept as gappa wrote
# them (text), so that a table written from this reads back to exactly the same numbers.
read_clade_multiplicities <- function(mult_dir) {
    tibble(path = list.files(mult_dir, full.names = TRUE)) %>%
        mutate(data = map(path, ~read_tsv(.x, col_names = c('sample_clade', 'geneid', 'count'), col_types = 'ccc'))) %>%
        unnest(data) %>%
        transmute(sample = str_remove(sample_clade, 'filter.*'),
                  clade = str_extract(sample_clade, 'filter(.*)', group = 1),
                  geneid, count)
}

# Summed multiplicity per clade and sample from the files of gappa_clade_multiplicities().
# The files are named <sample code>filter<clade>multiplicities.csv; sample_regex removes everything
# but the sample from the sample code ("filter.*" for metaT; metaG codes carry two more digits)
sum_clade_multiplicities <- function(mult_dir, sample_regex = 'filter.*') {
    tibble(path = list.files(mult_dir, full.names = TRUE)) %>%
        mutate(data = map(path, ~read_tsv(.x, col_names = c('clade', 'query', 'count'), show_col_types = FALSE))) %>%
        unnest(data) %>%
        group_by(clade) %>%
        summarize(count = sum(count)) %>%
        mutate(sample = str_remove(clade, sample_regex),
               clade = str_extract(clade, 'filter(.*)', group = 1)) %>%
        filter(clade != 'uncertain')
}

# Placement factorization of ratio jplaces (metaT/metaG or metaT/HK) against a low/high category.
#   samples:   `Sample Code` + catcol, one row per sample (samples without a category already removed)
#   ratio_dir: folder with <Sample Code>.jplace.gz carrying the ratios as multiplicities
#   metat_dir: the metaT read-count jplaces of the same samples (only for the pseudo-count)
# Samples without a ratio jplace, or whose ratios sum to <= 1 (gappa refuses them), are left out.
# Writes, in out_dir: the gappa outputs (<prefix>*), pseudocount.txt, metadata_factorization.tsv,
# jplace_list.txt. Returns the samples used, with their ratio jplace.
ratio_placement_factorization <- function(samples, catcol, ratio_dir, metat_dir, out_dir, prefix, threads = 8) {
    samples <- samples %>%
        select(`Sample Code`, all_of(catcol)) %>%
        mutate(ratio_jplace = file.path(ratio_dir, str_c(`Sample Code`, '.jplace.gz')),
               metat_jplace = file.path(metat_dir, str_c(`Sample Code`, '.jplace.gz'))) %>%
        filter(file.exists(ratio_jplace))

    message('Summing the multiplicities of ', nrow(samples), ' samples')
    samples <- samples %>%
        mutate(ratio_total = map_dbl(ratio_jplace, jplace_total),
               metat_total = map_dbl(metat_jplace, ~if (file.exists(.x)) jplace_total(.x) else NA_real_)) %>%
        filter(ratio_total > 1)

    pseudocount <- factorization_pseudocount(samples$ratio_total, samples$metat_total)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    writeLines(c(str_glue('samples\t{nrow(samples)}'),
                 str_glue('median_ratio_total\t{median(samples$ratio_total)}'),
                 str_glue('median_metaT_total\t{median(samples$metat_total, na.rm = TRUE)}'),
                 str_glue('pseudo_count_summand_all\t{format(pseudocount, scientific = FALSE, digits = 6)}')),
               file.path(out_dir, 'pseudocount.txt'))
    message(nrow(samples), ' samples (', str_flatten_comma(str_c(names(table(samples[[catcol]])), ' ', table(samples[[catcol]]))),
            '), pseudo-count ', signif(pseudocount, 4))

    unlink(Sys.glob(file.path(out_dir, str_c(prefix, '*'))))
    gappa_placement_factorization(samples$ratio_jplace, samples, catcol, out_dir, prefix, pseudocount, threads)
    samples
}
