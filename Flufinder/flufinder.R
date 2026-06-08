# flufinder.R — flu strain detection pipeline (BIMM143)

# Function 1: upload_fasta() — reads a fasta file, returns list of named protein sequences
upload_fasta <- function(fasta_filename) {
  library(seqinr)
  read.fasta(fasta_filename, seqtype = "AA", as.string = TRUE, set.attributes = FALSE)
}


# Function 2: trypsinize() — splits each protein into peptides after every R and K
trypsinize <- function(proteins) {
  library(stringr)
  lapply(proteins, str_split_1, pattern = "(?<=R|K)")
}


# Function 3: split_peptides() — splits each peptide into individual amino acids
split_peptides <- function(peptides) {
  library(stringr)
  lapply(peptides, str_split, pattern = "")
}


# Function 4: splitpeptides_to_masses() — calculates peptide masses from split amino acids
splitpeptides_to_masses <- function(aa) {
  aa_masses <- c(A=71.037, R=156.101, N=114.042, D=115.026, C=103.009,
                 Q=128.058, E=129.042, G=57.021, H=137.058, I=113.084, L=113.084,
                 K=128.094, M=131.040, F=147.068, P=97.052, S=87.032, T=101.047,
                 W=186.079, Y=163.063, V=99.068)
  peptide_masses <- aa
  for (i in 1:length(aa)) {
    peptide_masses[[i]] <- lapply(aa[[i]], function(x) sum(aa_masses[x]))
  }
  lapply(peptide_masses, unlist)
}


# Function 5: count_matching_masses() — counts how many sample masses match each strain
count_matching_masses <- function(protein_masses, sample) {
  df <- as.data.frame(sapply(protein_masses, function(x)
    sum(as.character(sample) %in% as.character(x))))
  names(df) <- "peptide_counts"
  return(df)
}


# Function 6: ggbarplot() — bar plot of the peptide counts data frame
ggbarplot <- function(peptide_counts_table) {
  library(ggplot2)
  ggplot(peptide_counts_table) +
    aes(rownames(peptide_counts_table), peptide_counts) +
    geom_col(fill = "blue", width = 0.5) +
    theme_bw() +
    labs(x = "Flu Strain", y = "Peptide Counts")
}