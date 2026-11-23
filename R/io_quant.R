# Reading expression tables.

# Gene x sample expression table (first column = gene id) -> samples x genes data.frame
# for the selected genes. `trim_sample_suffix` removes that many characters from the end of the
# sample names (the Carradec Pelagomonas table has a 2-character suffix, e.g. "_1").
expression_to_sample_matrix <- function(expr_table, genes, trim_sample_suffix = 0) {
    gene_col <- names(expr_table)[1]
    mat <- expr_table %>%
        filter(.data[[gene_col]] %in% genes) %>%
        column_to_rownames(var = gene_col) %>%
        t() %>%
        as.data.frame()
    if (trim_sample_suffix > 0) {
        rownames(mat) <- str_sub(rownames(mat), 1, -(trim_sample_suffix + 1))
    }
    mat
}
