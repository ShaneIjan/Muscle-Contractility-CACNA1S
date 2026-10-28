install.packages(c("rentrez", "ggplot2"))

# Load packages whenever you start a new R session.
library(rentrez)
library(ggplot2)

gene_symbol <- "CACNA1S"
transcript_accession <- "NM_000069.3"
clinvar_id <- "17625"
variant_protein <- "p.Arg528His"

for (x in c("data", "results", "figures")) {
  dir.create(x, showWarnings = FALSE)
}

# STEP 2: Retrieve the human CACNA1S gene record


hits <- entrez_search(
  db = "gene",
  term = paste0(gene_symbol, "[Gene Name] AND txid9606[Organism]"),
  retmax = 20
)

if (length(hits$ids) == 0) {
  stop("Gene not found")
}

# Verify the exact gene symbol
gene_id <- NA_character_

for (id in hits$ids) {
  rec <- entrez_summary(db = "gene", id = id)
  
  if (toupper(rec$name) == toupper(gene_symbol)) {
    gene_id <- id
    gene_info <- rec
    break
  }
}

if (is.na(gene_id)) {
  stop("Exact human gene not found")
}

# Display the retrieved gene information
print(gene_info)

# Save the gene summary
writeLines(
  c(
    paste("Gene symbol:", gene_info$name),
    paste("NCBI Gene ID:", gene_id),
    paste("Official name:", gene_info$description),
    paste("Chromosome:", gene_info$chromosome),
    paste(
      "Normal role in muscle:",
      "Encodes the alpha-1S subunit of the voltage-gated",
      "L-type calcium channel involved in skeletal muscle",
      "excitation-contraction coupling."
    ),
    "Source: https://www.ncbi.nlm.nih.gov/gene/779"
  ),
  "results/gene_summary.txt"
)

cat("Gene summary saved successfully.\n")


# STEP 3: Retrieve the matching protein sequence

# Retrieve the assigned transcript's GenBank record
gb <- entrez_fetch(
  db = "nuccore",
  id = transcript_accession,
  rettype = "gb",
  retmode = "text"
)

# Find the RefSeq NP protein accession in the record
pat <- '/protein_id="NP_[0-9]+[.][0-9]+"'
fragment <- regmatches(gb, regexpr(pat, gb))

if (length(fragment) == 0 || fragment == "") {
  stop("NP protein ID not found; inspect the transcript record")
}

protein_accession <- gsub(
  '^/protein_id="|"$',
  '',
  fragment
)

# Retrieve the matching protein sequence in FASTA format
protein_fasta <- entrez_fetch(
  db = "protein",
  id = protein_accession,
  rettype = "fasta",
  retmode = "text"
)

if (!startsWith(protein_fasta, ">")) {
  stop("Invalid FASTA")
}

# Save the protein sequence
fasta_file <- file.path(
  "data",
  paste0(gene_symbol, "_protein.fasta")
)

writeLines(protein_fasta, fasta_file)

# Display and save the protein accession
cat("Matched protein:", protein_accession, "\n")

writeLines(
  paste("NP accession:", protein_accession),
  "results/protein_accession.txt"
)

cat("Protein FASTA and accession saved successfully.\n")


# STEP 4: Analyze protein length and amino-acid composition

# Read the saved protein FASTA file
fasta_lines <- readLines(fasta_file)

# Extract the protein sequence, excluding the FASTA header
sequence <- paste(
  fasta_lines[!grepl("^>", fasta_lines)],
  collapse = ""
)

sequence <- toupper(gsub("\\s+", "", sequence))

# Calculate protein length
aa <- strsplit(sequence, "")[[1]]
protein_length <- length(aa)

# Define the 20 standard amino acids
standard <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

# Count each amino acid
counts <- table(factor(aa, levels = standard))

# Create and save the frequency table
results <- data.frame(
  AminoAcid = standard,
  Count = as.integer(counts)
)

write.csv(
  results,
  "results/amino_acid_composition.csv",
  row.names = FALSE
)

# Display protein length and the three most abundant amino acids
cat("Protein length:", protein_length, "amino acids\n")
cat("Three most abundant amino acids:\n")
print(head(results[order(-results$Count), ], 3))

# Create the amino-acid composition bar chart
fig1 <- ggplot(results, aes(AminoAcid, Count)) +
  geom_col(fill = "grey45") +
  theme_bw() +
  labs(
    title = paste(gene_symbol, "amino-acid composition"),
    x = "Amino acid",
    y = "Number of residues"
  )

print(fig1)

# Save the chart as a PNG image
ggsave(
  "figures/01_amino_acid_composition.png",
  fig1,
  width = 8,
  height = 5,
  dpi = 300
)

cat("Composition CSV and bar chart saved successfully.\n")

