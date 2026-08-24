# Hands-on CRISPR Screen Data Analysis  
## Cutadapt → Bowtie1 → sgRNA Count Matrix → MAGeCK

---

# 1. Overview

In a pooled CRISPR screen, sequencing does **not directly measure gene expression**.

Instead, sequencing measures the abundance of each **sgRNA** in the cell population.

The computational workflow is:

```text
Raw FASTQ
   │
   ▼
Cutadapt
Remove 5' constant sequence
Extract ~20-nt sgRNA
   │
   ▼
Trimmed sgRNA FASTQ
   │
   ▼
Bowtie1
Align reads to sgRNA library
   │
   ▼
SAM alignment
   │
   ▼
Count reads assigned to each sgRNA
   │
   ▼
sgRNA × sample count matrix
   │
   ▼
MAGeCK
Differential sgRNA / gene analysis
   │
   ▼
Enriched / depleted genes
```

The central question is:

> **Which sgRNAs become enriched or depleted after biological selection?**

---

# 2. Example Experimental Design

Suppose we have four samples:

```text
CTRL_R1.fastq.gz
CTRL_R2.fastq.gz
TRT_R1.fastq.gz
TRT_R2.fastq.gz
```

Experimental design:

| Sample | Condition | Replicate |
|---|---|---:|
| CTRL_R1 | Control | 1 |
| CTRL_R2 | Control | 2 |
| TRT_R1 | Treatment | 1 |
| TRT_R2 | Treatment | 2 |

Our final comparison will be:

```text
TRT_R1 + TRT_R2

       vs

CTRL_R1 + CTRL_R2
```

---

# 3. Software

For this tutorial we need:

```text
cutadapt
bowtie1
samtools
MAGeCK
R
```

A Conda environment can be created with:

```bash
conda create -n crispr_screen \
    -c conda-forge \
    -c bioconda \
    cutadapt bowtie samtools mageck r-base
```

Activate it:

```bash
conda activate crispr_screen
```

Check installation:

```bash
cutadapt --version
bowtie --version
samtools --version
mageck --version
```

---

# 4. Understand the Raw FASTQ Read

Before doing anything, inspect several reads.

```bash
zcat CTRL_R1.fastq.gz | head -20
```

A typical CRISPR sequencing read may look conceptually like:

```text
5' constant sequence          sgRNA          downstream sequence
       │                       │                      │
       ▼                       ▼                      ▼

ACACTCTTTCCCTACACGACGCT  GACCTGATCGTGACTGAGTA  GTTTTAGAGCTAGAA
                         └────────20 nt────────┘
```

What we ultimately want is:

```text
GACCTGATCGTGACTGAGTA
```

The exact 5′ constant sequence depends on:

- CRISPR library
- PCR primers
- sequencing primer
- library preparation protocol

Therefore, **inspect the actual FASTQ before choosing the Cutadapt sequence**.

---

# 5. Step 1 — Extract the sgRNA with Cutadapt

In this workflow, we assume the read structure is:

```text
5' constant sequence + 20-nt sgRNA + downstream sequence
```

For example:

```text
ACACTCTTTCCCTACACGACGCT GACCTGATCGTGACTGAGTA GTTTTAGAGCTAGAA
└──── 5' constant ─────┘ └────── sgRNA ──────┘
```

We first remove the known 5′ constant sequence, then retain the next 20 bases.

```bash
cutadapt \
    -g "^ACACTCTTTCCCTACACGACGCT" \
    --discard-untrimmed \
    -l 20 \
    --minimum-length 20 \
    --max-n 0 \
    -o CTRL_R1.sgRNA.fastq.gz \
    CTRL_R1.fastq.gz
```

The logic is:

```text
Raw read

5' constant + sgRNA + downstream sequence
             │
             ▼
Remove 5' constant sequence
             │
             ▼
sgRNA + downstream sequence
             │
             ▼
Keep first 20 nt
             │
             ▼
20-nt sgRNA
```

Example:

```text
Before:

ACACTCTTTCCCTACACGACGCTGACCTGATCGTGACTGAGTAGTTTTAGAGCTAGAA

                          ↓ Cutadapt

After:

GACCTGATCGTGACTGAGTA
```

---

# 6. Trim All Samples

For multiple samples:

```bash
for sample in CTRL_R1 CTRL_R2 TRT_R1 TRT_R2
do

    cutadapt \
        -g "^ACACTCTTTCCCTACACGACGCT" \
        --discard-untrimmed \
        -l 20 \
        --minimum-length 20 \
        --max-n 0 \
        -o ${sample}.sgRNA.fastq.gz \
        ${sample}.fastq.gz

done
```

Replace:

```text
ACACTCTTTCCCTACACGACGCT
```

with the actual 5′ constant sequence from your CRISPR-screen library.

---

# 7. sgRNA Reference Library

Now we need the original CRISPR library.

Example:

```text
sgRNA_ID    Sequence                Gene
sgTP53_1    GACCTGATCGTGACTGAGTA    TP53
sgTP53_2    GTGACTGACTGACTGACTGA    TP53
sgMYC_1     CTGACTGACTGACGTGACTG    MYC
sgMYC_2     ACTGACTGACGTGACTGACT    MYC
sgB2M_1     TGACTGACCTGACTGACTGA    B2M
```

Save it as:

```text
library.tsv
```

The required information is:

```text
sgRNA identifier
       +
sgRNA sequence
       +
target gene
```

---

# 8. Convert the sgRNA Library to FASTA

Bowtie requires a FASTA reference.

Convert:

```bash
awk 'BEGIN{FS="\t"} NR>1 {
    print ">"$1;
    print $2
}' library.tsv > sgRNA_library.fa
```

The resulting file looks like:

```text
>sgTP53_1
GACCTGATCGTGACTGAGTA

>sgTP53_2
GTGACTGACTGACTGACTGA

>sgMYC_1
CTGACTGACTGACGTGACTG

>sgB2M_1
TGACTGACCTGACTGACTGA
```

Each reference sequence represents **one sgRNA**.

---

# 9. Step 2 — Build the Bowtie1 Index

Build the index:

```bash
bowtie-build sgRNA_library.fa sgRNA_index
```

This creates Bowtie1 index files such as:

```text
sgRNA_index.1.ebwt
sgRNA_index.2.ebwt
sgRNA_index.3.ebwt
sgRNA_index.4.ebwt
sgRNA_index.rev.1.ebwt
sgRNA_index.rev.2.ebwt
```

Conceptually:

```text
sgRNA_library.fa
       │
       ▼
 bowtie-build
       │
       ▼
 searchable sgRNA index
```

---

# 10. Why Align Against the sgRNA Library?

We do **not** need to align these reads against the human genome.

Our question is not:

```text
Where in the genome did this read originate?
```

Instead, our question is:

```text
Which sgRNA in the CRISPR library
does this sequencing read represent?
```

Therefore:

```text
20-nt sequencing read
        │
        ▼
      Bowtie
        │
        ▼
20-nt library sgRNA
        │
        ▼
     sgRNA ID
        │
        ▼
    Target gene
```

---

# 11. Step 3 — Bowtie1 Alignment

For a simple teaching workflow, use exact matching:

```bash
zcat CTRL_R1.sgRNA.fastq.gz |
bowtie \
    -q \
    -v 0 \
    -m 1 \
    --best \
    --strata \
    -S \
    -p 4 \
    sgRNA_index \
    - \
    CTRL_R1.sam \
    2> CTRL_R1.bowtie.log
```

---

# 12. Understanding the Bowtie Parameters

### FASTQ input

```bash
-q
```

means the input is FASTQ.

### Exact matching

```bash
-v 0
```

means:

```text
0 mismatches allowed
```

For a short ~20-nt sgRNA, exact matching provides a simple and conservative starting point.

### Unique assignment

```bash
-m 1
```

suppresses reads with more than one reportable alignment.

Therefore:

```text
one read
   ↓
one unique sgRNA
```

### Best alignment

```bash
--best --strata
```

asks Bowtie to report the best available alignment stratum.

### SAM output

```bash
-S
```

produces SAM format.

### Threads

```bash
-p 4
```

uses four CPU threads.

---

# 13. Optional: Allow One Mismatch

A less stringent analysis could use:

```bash
-v 1
```

For example:

```bash
zcat CTRL_R1.sgRNA.fastq.gz |
bowtie \
    -q \
    -v 1 \
    -m 1 \
    --best \
    --strata \
    -S \
    -p 4 \
    sgRNA_index \
    - \
    CTRL_R1.v1.sam \
    2> CTRL_R1.v1.bowtie.log
```

Example:

```text
sgRNA reference:
GACCTGATCGTGACTGAGTA

sequencing read:
GACCTGATCGTGACTGAGTT
                   ↑
              1 mismatch
```

For teaching:

```text
Primary analysis     → -v 0
Sensitivity analysis → -v 1
```

---

# 14. Run Bowtie for All Samples

```bash
for sample in CTRL_R1 CTRL_R2 TRT_R1 TRT_R2
do

    zcat ${sample}.sgRNA.fastq.gz |
    bowtie \
        -q \
        -v 0 \
        -m 1 \
        --best \
        --strata \
        -S \
        -p 4 \
        sgRNA_index \
        - \
        ${sample}.sam \
        2> ${sample}.bowtie.log

done
```

Now we have:

```text
CTRL_R1.sam
CTRL_R2.sam
TRT_R1.sam
TRT_R2.sam
```

---

# 15. Look at the Bowtie Output

Inspect:

```bash
head CTRL_R1.sam
```

A mapped read may look conceptually like:

```text
READ001   0   sgTP53_1   1   255   20M ...
```

The important field is:

```text
sgTP53_1
```

because this tells us:

```text
READ001
   ↓
mapped to
   ↓
sgTP53_1
   ↓
targets
   ↓
TP53
```

---

# 16. Step 4 — Count Reads per sgRNA

Now convert millions of read alignments into one count per sgRNA.

For each sample:

```bash
samtools view -F 4 CTRL_R1.sam |
cut -f3 |
sort |
uniq -c |
awk 'BEGIN{OFS="\t"} {print $2,$1}' \
> CTRL_R1.counts.tsv
```

Example output:

```text
sgTP53_1    1823
sgTP53_2    1672
sgMYC_1     553
sgMYC_2     612
sgB2M_1     2250
```

Interpretation:

```text
sgTP53_1 → 1,823 sequencing reads
sgTP53_2 → 1,672 sequencing reads
sgMYC_1  →   553 sequencing reads
```

---

# 17. Count All Samples

```bash
for sample in CTRL_R1 CTRL_R2 TRT_R1 TRT_R2
do

    samtools view -F 4 ${sample}.sam |
    cut -f3 |
    sort |
    uniq -c |
    awk 'BEGIN{OFS="\t"} {print $2,$1}' \
    > ${sample}.counts.tsv

done
```

We now have:

```text
CTRL_R1.counts.tsv
CTRL_R2.counts.tsv
TRT_R1.counts.tsv
TRT_R2.counts.tsv
```

---

# 18. Why Add Zero Counts?

Suppose a guide exists in the library but is not observed in one sample.

It will be absent from:

```text
CTRL_R1.counts.tsv
```

But MAGeCK needs the complete library:

```text
sgA    100
sgB    50
sgC    0
sgD    75
```

rather than:

```text
sgA    100
sgB    50
sgD    75
```

Therefore, counts should be merged back to the complete sgRNA library and missing guides assigned:

```text
0 reads
```

---

# 19. Step 5 — Build the MAGeCK Count Matrix

Use R:

```r
# ============================================================
# Build MAGeCK count table
# ============================================================

library <- read.delim(
    "library.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

samples <- c(
    "CTRL_R1",
    "CTRL_R2",
    "TRT_R1",
    "TRT_R2"
)

# Start with sgRNA and gene annotation
count_table <- library[, c("sgRNA_ID", "Gene")]

colnames(count_table)[1] <- "sgRNA"

# Add counts from each sample
for (sample in samples) {

    x <- read.delim(
        paste0(sample, ".counts.tsv"),
        header = FALSE,
        stringsAsFactors = FALSE
    )

    colnames(x) <- c("sgRNA", sample)

    count_table <- merge(
        count_table,
        x,
        by = "sgRNA",
        all.x = TRUE,
        sort = FALSE
    )
}

# Missing guides = 0 reads
count_table[samples] <-
    lapply(
        count_table[samples],
        function(x) {
            x[is.na(x)] <- 0
            as.integer(x)
        }
    )

# Rearrange columns
count_table <- count_table[
    ,
    c(
        "sgRNA",
        "Gene",
        samples
    )
]

# Save
write.table(
    count_table,
    "CRISPR_screen.count.txt",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
```

---

# 20. Final MAGeCK Input Table

The result should look like:

```text
sgRNA       Gene     CTRL_R1 CTRL_R2 TRT_R1 TRT_R2
sgTP53_1    TP53       1823    1755    5620   5310
sgTP53_2    TP53       1672    1599    4980   5122
sgMYC_1     MYC         553     590      92     81
sgMYC_2     MYC         612     630     105     98
sgB2M_1     B2M        2250    2180     420    390
```

This is the central input for MAGeCK:

```text
sgRNA
  +
Gene
  +
raw read counts for every sample
```

---

# 21. Step 6 — Run MAGeCK Differential Analysis

Now compare treatment against control:

```bash
mageck test \
    -k CRISPR_screen.count.txt \
    -t TRT_R1,TRT_R2 \
    -c CTRL_R1,CTRL_R2 \
    -n TRT_vs_CTRL \
    --norm-method median
```

The comparison is:

```text
Treatment
TRT_R1 + TRT_R2

       versus

Control
CTRL_R1 + CTRL_R2
```

---

# 22. What MAGeCK Does

Conceptually:

```text
sgRNA count matrix
        │
        ▼
normalization
        │
        ▼
estimate sgRNA abundance changes
        │
        ▼
rank sgRNAs
        │
        ▼
combine sgRNAs targeting same gene
        │
        ▼
gene-level statistics
        │
        ▼
enriched / depleted genes
```

---

# 23. MAGeCK Output

Important files include:

```text
TRT_vs_CTRL.sgrna_summary.txt
TRT_vs_CTRL.gene_summary.txt
TRT_vs_CTRL.log
```

---

# 24. sgRNA-Level Results

Look at:

```bash
head TRT_vs_CTRL.sgrna_summary.txt
```

This file describes how individual sgRNAs behave.

Conceptually:

| sgRNA | Gene | Control | Treatment | Effect |
|---|---|---:|---:|---|
| sgA1 | GeneA | 2,000 | 200 | depleted |
| sgA2 | GeneA | 1,800 | 250 | depleted |
| sgB1 | GeneB | 500 | 3,000 | enriched |

Individual sgRNAs are useful for checking whether a gene-level hit is supported by multiple guides.

---

# 25. Gene-Level Results

Look at:

```bash
head TRT_vs_CTRL.gene_summary.txt
```

MAGeCK reports separate statistics for:

```text
negative selection

and

positive selection
```

Important columns include:

```text
neg|score
neg|p-value
neg|fdr
neg|rank

pos|score
pos|p-value
pos|fdr
pos|rank
```

---

# 26. Negative Selection

Suppose:

```text
Control:

sgGeneA_1   ████████████
sgGeneA_2   ███████████
sgGeneA_3   █████████████


Treatment:

sgGeneA_1   ██
sgGeneA_2   █
sgGeneA_3   ██
```

All Gene A guides became depleted.

Interpretation:

```text
Gene A knockout
      ↓
reduced fitness / survival
      ↓
cells carrying Gene A sgRNAs disappear
      ↓
sgRNAs become depleted
      ↓
negative selection
```

Therefore:

```text
small neg|FDR
```

indicates evidence for negative selection.

---

# 27. Positive Selection

Suppose:

```text
Control:

sgGeneB_1   ██
sgGeneB_2   ███
sgGeneB_3   ██


Treatment:

sgGeneB_1   █████████████
sgGeneB_2   ████████████
sgGeneB_3   ██████████████
```

Interpretation:

```text
Gene B knockout
      ↓
selective advantage
      ↓
cells expand
      ↓
Gene B sgRNAs become enriched
      ↓
positive selection
```

Therefore:

```text
small pos|FDR
```

indicates evidence for positive selection.

---

# 28. Why Gene-Level Analysis Matters

Suppose Gene A has four sgRNAs:

```text
sgA1     strongly depleted
sgA2     strongly depleted
sgA3     strongly depleted
sgA4     strongly depleted
```

That is convincing.

But consider Gene B:

```text
sgB1     extremely depleted
sgB2     unchanged
sgB3     unchanged
sgB4     unchanged
```

The apparent effect may be caused by:

```text
off-target activity
technical artifact
poor guide behavior
```

Therefore, gene-level analysis asks whether **multiple independent guides targeting the same gene show consistent evidence**.

---

# 29. A Useful Hit-Calling Strategy

For initial exploration:

```text
FDR < 0.05
```

can be used as a statistical threshold.

But a strong candidate should ideally satisfy:

```text
Significant MAGeCK FDR
        +
consistent effect direction
        +
multiple supporting sgRNAs
        +
good replicate behavior
        +
biological plausibility
```

---

# 30. Important Interpretation Example

Suppose this is a tumor-cell screen in which:

```text
CRISPR KO tumor cells
        ↓
co-culture with T cells
        ↓
sequence surviving tumor cells
```

If sgRNAs targeting:

```text
Gene X
```

are strongly **depleted** after T-cell exposure:

```text
Gene X knockout
       ↓
tumor cells become more sensitive
to T-cell-mediated killing
       ↓
fewer surviving cells
       ↓
Gene X sgRNAs depleted
```

Therefore:

```text
Gene X may normally protect tumor cells
from immune-mediated killing.
```

The biological meaning of enrichment/depletion always depends on **which cell population was sequenced**.

---

# 31. Complete Demo Pipeline

```bash
# ============================================================
# STEP 1 — CUTADAPT
# ============================================================

for sample in CTRL_R1 CTRL_R2 TRT_R1 TRT_R2
do

    cutadapt \
        -g "^ACACTCTTTCCCTACACGACGCT" \
        --discard-untrimmed \
        -l 20 \
        --minimum-length 20 \
        --max-n 0 \
        -o ${sample}.sgRNA.fastq.gz \
        ${sample}.fastq.gz

done


# ============================================================
# STEP 2 — CREATE sgRNA FASTA
# ============================================================

awk 'BEGIN{FS="\t"} NR>1 {
    print ">"$1;
    print $2
}' library.tsv > sgRNA_library.fa


# ============================================================
# STEP 3 — BUILD BOWTIE1 INDEX
# ============================================================

bowtie-build \
    sgRNA_library.fa \
    sgRNA_index


# ============================================================
# STEP 4 — ALIGN EACH SAMPLE
# ============================================================

for sample in CTRL_R1 CTRL_R2 TRT_R1 TRT_R2
do

    zcat ${sample}.sgRNA.fastq.gz |
    bowtie \
        -q \
        -v 0 \
        -m 1 \
        --best \
        --strata \
        -S \
        -p 4 \
        sgRNA_index \
        - \
        ${sample}.sam \
        2> ${sample}.bowtie.log

done


# ============================================================
# STEP 5 — COUNT READS PER sgRNA
# ============================================================

for sample in CTRL_R1 CTRL_R2 TRT_R1 TRT_R2
do

    samtools view -F 4 ${sample}.sam |
    cut -f3 |
    sort |
    uniq -c |
    awk 'BEGIN{OFS="\t"} {print $2,$1}' \
    > ${sample}.counts.tsv

done


# ============================================================
# STEP 6 — MERGE COUNTS
# ============================================================

Rscript build_count_table.R


# ============================================================
# STEP 7 — MAGeCK
# ============================================================

mageck test \
    -k CRISPR_screen.count.txt \
    -t TRT_R1,TRT_R2 \
    -c CTRL_R1,CTRL_R2 \
    -n TRT_vs_CTRL \
    --norm-method median


# ============================================================
# STEP 8 — INSPECT RESULTS
# ============================================================

head TRT_vs_CTRL.sgrna_summary.txt

head TRT_vs_CTRL.gene_summary.txt
```

---

# 32. The Central Concept

The full analysis can be reduced to four transformations:

```text
SEQUENCE

FASTQ read
   │
   ▼ Cutadapt

20-nt sgRNA
   │
   ▼ Bowtie


IDENTITY

sgRNA ID
   │
   ▼ count reads


ABUNDANCE

sgRNA × sample count matrix
   │
   ▼ MAGeCK


BIOLOGY

Enriched / depleted sgRNAs
   │
   ▼
Gene-level CRISPR hits
```

---

# 33. Take-Home Message

For this workflow:

```text
Cutadapt
```

answers:

> **How do I extract the sgRNA sequence from the sequencing read?**

```text
Bowtie1
```

answers:

> **Which library sgRNA does this read correspond to?**

```text
Count table
```

answers:

> **How abundant is each sgRNA in each sample?**

```text
MAGeCK
```

answers:

> **Which sgRNAs and genes are significantly enriched or depleted between conditions?**

The complete workflow is:

```text
FASTQ
  ↓
Cutadapt
  ↓
20-nt sgRNA reads
  ↓
Bowtie1
  ↓
sgRNA identity
  ↓
read counting
  ↓
sgRNA × sample count matrix
  ↓
MAGeCK
  ↓
sgRNA-level differential signal
  ↓
gene-level statistics
  ↓
CRISPR screen hits
```