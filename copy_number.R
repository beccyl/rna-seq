library("VariantAnnotation")
options(scipen=5)

##
get_snp_allelic_frequencies <- function (vcf, seqlevel, dp_cutoff=5) {
    ## Select only those rows matching the seqlevel required
    ad<-geno(vcf)$AD[grepl(paste0("^", seqlevel, ":"),rownames(geno(vcf)$AD)),]
    dp<-geno(vcf)$DP[grepl(paste0("^", seqlevel, ":"),rownames(geno(vcf)$DP)),]
    select_names <- names(dp[dp >= dp_cutoff])
    selection <- mapply(
        function(x,y) x/y, ad[select_names], dp[select_names], SIMPLIFY=FALSE
    )
    allelic_frequencies <- vapply(selection, function(x) max(x), numeric(1))
    allelic_positions <- sapply(names(allelic_frequencies),
        function (x) {
            re=regexpr(":.*_", x)
            as.numeric(substring(x, re+1, re+attr(re, "match.length")-2))
        }
    )
    cbind(position=allelic_positions,frequency=allelic_frequencies, dp=dp[dp >= dp_cutoff])
}

## draw a png graph
## input: vcf
##      : seqlevel - the chromosome to plot
##      : dp_cutoff - the depth cutoff to plot
plot_graph <- function (vcf, seqlevel, dp_cutoff=5, outdir="plots", region="") {
    snp_freq_mat <- get_snp_allelic_frequencies(vcf, seqlevel, dp_cutoff)

    if (length(snp_freq_mat) > 0) {
        if (!file.exists(outdir)) {
            dir.create(outdir)
        }
        png(filename=paste0(outdir, "/chr", seqlevel,region,"_dp", dp_cutoff, "_Rplot.png"),
            width=960, height=480, units="px"
        )
        plot(x=snp_freq_mat[,"position"], y=snp_freq_mat[,"frequency"], ylim=c(0.5,1),
            cex=log(sqrt(snp_freq_mat[,"dp"]))/2, pch=19, 
            xlab="position [bp]", ylab="frequency [allelic depth (AD) / read depth (DP)]",
            main=paste("Chromosome", seqlevel, "SNP frequency")
        )
        dev.off()
    }
}

##
plot_difference_graph <- function (vcf1, vcf2, seqlevel, dp_cutoff=5, outdir="plots") {
    snp_freq_mat1 <- get_snp_allelic_frequencies(vcf1, seqlevel, dp_cutoff)
    snp_freq_mat2 <- get_snp_allelic_frequencies(vcf2, seqlevel, dp_cutoff)
    diff_names <- setdiff(rownames(snp_freq_mat1), rownames(snp_freq_mat2))
    snp_freq_mat3 <- snp_freq_mat1[diff_names,]

    ## number of graphs to plot (0, 1 or 2)
    n_graphs <- (length(snp_freq_mat1) > 0) + (length(snp_freq_mat3) > 0)

    if (n_graphs > 0) {
        if (!file.exists(outdir)) {
            dir.create(outdir)
        }
        png(filename=paste0(outdir, "/chr", seqlevel,"_dp", dp_cutoff, "_Rplot.png"),
            width=960, height=480*n_graphs, units="px"
        )
        par(mfrow = c(n_graphs, 1))
        if (length(snp_freq_mat1) > 0) {
            plot(x=snp_freq_mat1[,"position"], y=snp_freq_mat1[,"frequency"], ylim=c(0.5,1),
                cex=log(sqrt(snp_freq_mat1[,"dp"]))/2, pch=19, 
                xlab="position [bp]", ylab="frequency [allelic depth (AD) / read depth (DP)]",
                main=paste("Chromosome", seqlevel, "SNP frequency")
            )
        } 
        if (length(snp_freq_mat3) > 0) {
            plot(x=snp_freq_mat3[,"position"], y=snp_freq_mat3[,"frequency"], ylim=c(0.5,1),
                cex=log(sqrt(snp_freq_mat3[,"dp"]))/2, pch=19, 
                xlab="position [bp]", ylab="frequency [allelic depth (AD) / read depth (DP)]",
                main=paste("Chromosome", seqlevel, "SNP frequency (potential deleterious omitted)")
            )
        }
        dev.off()
    }
}

## default is HLA-B
graph_region <- function (seqlevel=6, startpos=31353871, endpos=31357211, region="_HLA-B") {
    compressVcf <- bgzip("variants.vcf", tempfile())
    idx <- indexTabix(compressVcf, "vcf")
    tab <- TabixFile(compressVcf, idx)
    gr <- GRanges(seqlevel, IRanges(startpos, endpos, name=region))
    vcf <- readVcf(tab, "hg38", gr)
    plot_graph(vcf, seqlevel, region=region)
}

### 31350000 - 31360000
## Example Usage:
#vcf1<-readVcf("variants.vcf","hg38")
#vcf2<-readVcf("vep_out_sifted.vcf","hg38")

## Which chromosomes do we want to look at? (ALL, skip MT and KI/GL)
#seqlevels <- c(
#    extractSeqlevelsByGroup("Homo_sapiens", style="NCBI", group="auto"),
#    extractSeqlevelsByGroup("Homo_sapiens", style="NCBI", group="sex")
#)

## set dp_cutoff as desired
## set outdir as desired
#sapply(seqlevels, function(x) plot_graph(vcf=vcf, seqlevel=x))
#sapply(seqlevels, function(x) plot_difference_graph(vcf1=vcf1, vcf2=vcf2, seqlevel=x))
 graph_region(6, 31353871, 31357211, "_HLA-B")
