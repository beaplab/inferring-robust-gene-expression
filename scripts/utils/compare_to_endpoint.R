# Compare a regenerated figure with its endpoint PDF: rasterise both (pdftoppm, 100 dpi) and count the
# pixels that differ. Writes a diff image (differing pixels in red) next to the new figure.
# Usage: /usr/bin/Rscript scripts/utils/compare_to_endpoint.R <endpoint.pdf> <new.pdf>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
tmp <- tempfile()
dir.create(tmp)
raster <- function(pdf, prefix) {
    system2('pdftoppm', c('-r', '100', '-png', '-singlefile', shQuote(pdf), file.path(tmp, prefix)))
    png::readPNG(file.path(tmp, paste0(prefix, '.png')))[, , 1:3]
}
a <- raster(args[1], 'endpoint')
b <- raster(args[2], 'new')
if (!identical(dim(a), dim(b))) {
    cat(sprintf('Different sizes: endpoint %dx%d, new %dx%d\n', dim(a)[2], dim(a)[1], dim(b)[2], dim(b)[1]))
    quit(status = 1)
}
differ <- apply(abs(a - b), c(1, 2), max) > 0.01
cat(sprintf('%s: %d of %d pixels differ (%.3f%%)\n', basename(args[2]), sum(differ), length(differ),
            100 * mean(differ)))
if (any(differ)) {
    out <- b * 0.3 + 0.7
    out[, , 1][differ] <- 1; out[, , 2][differ] <- 0; out[, , 3][differ] <- 0
    png::writePNG(out, sub('\\.pdf$', '_diff.png', args[2]))
}
