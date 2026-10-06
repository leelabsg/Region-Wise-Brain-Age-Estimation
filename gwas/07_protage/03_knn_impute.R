# Stage 03 - KNN imputation of the QC-passed protein matrix (protein columns only).
# impute::impute.knn(t(X), k = 10, rowmax = 0.5, colmax = 0.95, maxp = ncol(X), rng.seed = 362436069):
# rows = proteins, columns = participants; the result is transposed back. No age, sex, clinical or batch
# variables enter the matrix. Observed values must be unchanged.
# Usage: PROTAGE_WORK_DIR=... OMP_NUM_THREADS=1 Rscript 03_knn_impute.R    (impute 1.76.0, arrow, R 4.3)
suppressPackageStartupMessages({library(impute); library(arrow)})
W <- Sys.getenv("PROTAGE_WORK_DIR"); if (W == "") stop("PROTAGE_WORK_DIR is not set")
d <- as.data.frame(read_parquet(file.path(W, "preprocessing", "knn_input.parquet")))
eid <- d$eid; X <- as.matrix(d[, setdiff(names(d), "eid")]); rm(d)
storage.mode(X) <- "double"
cat(sprintf("input %d participants x %d proteins; missing cells %d\n", nrow(X), ncol(X), sum(is.na(X))))
res <- impute.knn(t(X), k = 10, rowmax = 0.5, colmax = 0.95, maxp = ncol(X), rng.seed = 362436069)
Xi <- t(res$data)
stopifnot(identical(dim(Xi), dim(X)), identical(colnames(Xi), colnames(X)))
obs <- !is.na(X)
stopifnot(all(Xi[obs] == X[obs]), !anyNA(Xi), all(is.finite(Xi)))
out <- data.frame(eid = eid, Xi, check.names = FALSE)
write_parquet(out, file.path(W, "preprocessing", "npx_pool_imputed.parquet.tmp"))
file.rename(file.path(W, "preprocessing", "npx_pool_imputed.parquet.tmp"), file.path(W, "preprocessing", "npx_pool_imputed.parquet"))
writeLines(format(Sys.time(), tz = "UTC", usetz = TRUE), file.path(W, "logs", "s02.done"))
cat("done; observed values unchanged; residual NA 0\n")
