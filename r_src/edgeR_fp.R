library(edgeR)

calculate.correlations <- function(geno.df, pheno.df, norm.method="TMM",
    fdr.method="BH") {
    # Calculates gene-response correlations using edgeR.
    #
    # Args:
    #    geno.df: Dataframe of raw read counts oriented such that columns
    #        represent samples and rows represent genes.
    #    pheno.df: Dataframe of values for a single binary phenotype oriented
    #        such that columns represent samples and the one row represents the
    #        one phenotype. The set of sample identifiers must be identical to
    #        geno.df's set of sample identifiers; i.e., both dataframes must
    #        have exactly the same samples in exactly the same order. **This
    #        function does not validate the sample identifiers on the inputs.
    #        For efficiency's sake, that's the caller's responsibility.**
    #    norm.method: one of "TMM", "TMMwsp", "RLE", "upperquartile", and "none"
    #    fdr.method: one of "holm", "hochberg", "hommel", "bonferroni", "BH",
    #        "BY", "fdr", "none"
    #
    # Returns:
    #    a dataframe consisting of a single column named "correlation_array" and
    #    one row for each gene, where the rows are ordered to match geno.df
    
    pheno.group <- factor(pheno.df[1, ])

    # consider removing genes by minimum counts for certain number of samples

    y <- DGEList(counts = geno.df, group = pheno.group, remove.zeros = F)

    y <- calcNormFactors(y, method = norm.method)

    design <- model.matrix(~pheno.group)
    y <- estimateDisp(y, design)

    # to perform quasi-likelihood F-tests: (regular lrt tests worse for bulk)
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit, coef = 2)
    DE_results <- topTags(object = qlf, n = dim(geno.df)[1],
                          sort.by = "PValue", p.value = 1.1,
                          adjust.method = fdr.method )

    correlation_array <- (sign(DE_results$table$logFC) *
                            -1 * log10(DE_results$table$PValue))

    correlation_matrix <- matrix(data = 0,
                                 nrow = dim(geno.df)[1],
                                 ncol = 1,
                                 dimnames = list(rownames(geno.df),
                                                 "correlation_array"))
    correlation_matrix[rownames(DE_results$table), 1] <- correlation_array
    return(correlation_matrix)
}

args = commandArgs(trailingOnly=TRUE)

geno.df <- feather::read_feather(args[1])
rownames(geno.df) <- geno.df$index
geno.df$index <- NULL

pheno.df <- feather::read_feather(args[2])
rownames(pheno.df) <- pheno.df$index
pheno.df$index <- NULL

out.df = calculate.correlations(geno.df, pheno.df)
feather::write_feather(out.df, args[3])
