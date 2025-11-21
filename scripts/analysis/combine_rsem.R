library(readr)
library(dplyr)
library(tximport)
library(apeglm)
library(DESeq2)
library(here)
library(fs)

## gene mapping
id_map <- "/lustre1/g/path_my/Reference/Gencode/gencode.hg38.v44/GRCh38.primary_assembly.genecode.v44.id_map.txt"
df_map_isoform <- read_tsv(id_map, col_names = FALSE)
colnames(df_map_isoform) <- c("Gene", "Isoform", "Symbol")
df_map_gene <- dplyr::select(df_map_isoform, Gene, Symbol) %>% 
    dplyr::filter(!duplicated(Gene)) # remove duplicated records at gene level

## Find all the genes.results
rsem_files <- dir_ls(here("data/DFSP/RSEM"), glob = "*genes.results", recurse = TRUE)
samples <- path_file(path_dir(rsem_files))
names(rsem_files) <- samples

txi <- tximport(rsem_files, type = "rsem", txIn = FALSE, txOut = FALSE)

# technically, RSEM produce effect gene length of 0 when gene length < read length
# it can be converted to 1 to allow DESeq2 estimate library size
txi$length[txi$length == 0] <- 1

rownames(txi$counts) <- gsub("\\.\\d+$", "", rownames(txi$counts))
rownames(txi$abundance) <- gsub("\\.\\d+$", "", rownames(txi$abundance))
rownames(txi$length) <- gsub("\\.\\d+$", "", rownames(txi$length))

# convert expression matrix to table and map gene ID
convert_expr_mat <- function(x, type, df_map) {
    if (type == "Gene") {
        df_new <- as_tibble(x, rownames = "Gene")
        df_new <- right_join(df_map, df_new, by = "Gene")
    } else if (type == "Isoform") {
        df_new <- as_tibble(x, rownames = "Isoform")
        df_new <- right_join(df_map, df_new, by = "Isoform")
    }
    return(df_new)
}

calculate_fpkm <- function(counts, lengths) {
    # FPKM = (counts * 10^9) / (effective_length * total_mapped_reads)
    total_counts <- colSums(counts)
    fpkm <- sweep(counts, 2, total_counts, "/") * 1e6  # RPM
    fpkm <- fpkm / lengths * 1000  # FPKM
    return(fpkm)
}

fpkm_matrix <- calculate_fpkm(txi$counts, txi$length)
rownames(fpkm_matrix) <- gsub("\\.\\d+$", "", rownames(fpkm_matrix))

df_rsem_gene_counts <- convert_expr_mat(txi$counts, "Gene", df_map_gene)
df_rsem_gene_tpm <- convert_expr_mat(txi$abundance, "Gene", df_map_gene)
df_rsem_gene_fpkm <- convert_expr_mat(fpkm_matrix, "Gene", df_map_gene)

write_csv(df_rsem_gene_counts, here("data/DFSP/RSEM/gene_expression_matrix_counts.csv"))
write_csv(df_rsem_gene_tpm, here("data/DFSP/RSEM/gene_expression_matrix_tpm.csv"))
write_csv(df_rsem_gene_fpkm, here("data/DFSP/RSEM/gene_expression_matrix_fpkm.csv"))
