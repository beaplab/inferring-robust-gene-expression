# Reading jplace files (phylogenetic placements, EPA-ng / gappa).
#
# A jplace holds the reference tree, with every edge numbered as {n} in the newick, and one entry per
# placed query ("pquery"): its candidate edges (p, with the columns named in "fields", e.g. edge_num,
# like_weight_ratio, pendant_length) and its names with multiplicities (nm) or only names (n).

# Parsed jplace (json list)
read_jplace <- function(file) {
    jsonlite::fromJSON(file, simplifyVector = FALSE)
}

# Multiplicity of a pquery: sum of its name multiplicities, or the number of names without them
pquery_multiplicity <- function(pquery) {
    if (is.null(pquery$nm)) length(pquery$n) else sum(map_dbl(pquery$nm, 2))
}

# Sum of the multiplicities of a jplace (its total "mass")
jplace_total <- function(file) {
    sum(map_dbl(read_jplace(file)$placements, pquery_multiplicity))
}

# Reference tree with the jplace edge numbers {n} moved into the node labels (_E<n>), so that every
# ape node knows the number of the edge above it
read_jplace_tree <- function(newick) {
    newick %>%
        str_replace_all(':([^:{}(),;]+)\\{(\\d+)\\}', '_E\\2:\\1') %>%
        str_replace_all('\\)\\{(\\d+)\\}', ')_E\\1') %>%
        ape::read.tree(text = .)
}

# Tree and one row per candidate placement: pquery id, edge, LWR, multiplicity and whether it is the
# best placement of the pquery
read_placements <- function(file) {
    j <- read_jplace(file)
    fields <- unlist(j$fields)
    edge_i <- which(fields == 'edge_num')
    lwr_i <- which(fields == 'like_weight_ratio')
    pq <- imap_dfr(j$placements, ~{
        tibble(pquery = .y,
               edge_num = map_int(.x$p, ~as.integer(.x[[edge_i]])),
               lwr = map_dbl(.x$p, lwr_i),
               multiplicity = pquery_multiplicity(.x)) %>%
            mutate(best = lwr == max(lwr))
    })
    list(tree = read_jplace_tree(j$tree), pq = pq)
}

# One row per pquery: LWR, pendant length and edge of the best placement, number of candidate
# placements and multiplicity. Also returns the newick of the reference tree.
read_pqueries <- function(file) {
    j <- read_jplace(file)
    fields <- unlist(j$fields)
    lwr_i <- which(fields == 'like_weight_ratio')
    pend_i <- which(fields == 'pendant_length')
    best <- map(j$placements, ~{
        lwrs <- map_dbl(.x$p, lwr_i)
        top <- which.max(lwrs)
        c(lwrs[top], length(lwrs), .x$p[[top]][[pend_i]])
    })
    # as.numeric keeps samples without placements as a 0-row matrix
    best <- matrix(as.numeric(unlist(best)), ncol = 3, byrow = TRUE)
    list(tree = j$tree,
         pq = tibble(lwr1 = best[, 1],
                     n_alt = best[, 2],
                     pendant1 = best[, 3],
                     multiplicity = map_dbl(j$placements, pquery_multiplicity)))
}

# Size and shape of a reference tree (newick of a jplace; edge numbers removed).
# Treeness: fraction of the total branch length on internal branches
tree_stats <- function(newick) {
    tree <- ape::read.tree(text = str_remove_all(newick, '\\{\\d+\\}'))
    internal <- tree$edge[, 2] > length(tree$tip.label)
    tibble(nnodes = tree$Nnode,
           ntips = length(tree$tip.label),
           branch_length = sum(tree$edge.length),
           treeness = sum(tree$edge.length[internal]) / branch_length) %>%
        mutate(branch_length_per_tip = branch_length / ntips)
}
