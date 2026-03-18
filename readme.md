# Ensembl VEP REST Variant Annotation Pipeline

A Python pipeline for transcript-level human variant annotation using the Ensembl Variant Effect Predictor (VEP) REST API.

## Short description

This script annotates genomic variants through the Ensembl VEP REST API and writes the results to a tab-separated output file. The workflow is designed to preserve **one row per transcript consequence**, allowing transcript-specific consequence, prediction, and ClinVar-related information to be retained.

## Project summary

This project was developed to automate REST-based variant annotation for human genomic variants. The pipeline was designed to:

- annotate curated known variants
- process a reproducible subset of variants from a VCF file
- recover transcript-level consequence information
- extract HGVS, SIFT, PolyPhen, AlphaMissense, and ClinVar-related annotations
- produce structured TSV output for downstream analysis

The workflow was applied to:
- a curated BRCA1 validation set
- the first 200 valid variants from the HG01868 chromosome 22 VCF

## Repository contents

```text
.
├── readme.md
├── vep_annotate.py
├── data/
│   ├── brca1_known.txt
│   └── HG01868.chr22.vcf
├── output/
│   ├── brca1_validation.tsv
│   └── chr22_first200.tsv
└── figures/
    ├── figure1_consequence_terms.png
    ├── figure2_alphamissense_scores.png
    ├── figure3_mean_alphamissense_by_variant.png
    ├── figure4_workflow.png
    ├── table1_brca1_summary.png
    └── table2_chr22_summary.png

```

## Input modes
The script supports two modes.

1. known
Annotates a small curated set of variants from a text file.
Accepted input format:
chrom pos id ref alt
         or
chrom pos ref alt

Example:
17 43057117 . C G
17 43104882 T A
17 43047643 C A

2. vcf_firstN
Reads a VCF file and annotates the first N valid variants, up to a maximum of 200 variants per POST request.

The script:
- ignores VCF header lines
- splits multi-allelic records into separate variants
- removes unsupported or symbolic alleles
- converts valid variants into Ensembl /region format


## Output
The script writes a TSV file containing one row per transcript consequence.

Example output fields include:
- input
- chrom
- pos
- allele_string
- most_severe_consequence
- consequence_terms
- impact
- gene_id
- gene_symbol
- transcript_id
- gene_biotype
- hgvsg
- hgvsc
- hgvsp
- sift_prediction
- sift_score
- polyphen_prediction
- polyphen_score
- alphamissense_score
- alphamissense_prediction
- clinvar_clin_sig
- colocated_variant_ids
- clinvar_phenotypes

If no transcript consequence is returned, the script still writes a minimal row so that non-coding or intergenic variants are retained.

## Requirements
This script uses only the Python standard library.

Recommended:
- Python 3.10+

No third-party packages are required.


## Usage

**Example commands used in this project:**

conda create -n vcf_env python=3.10


conda activate vcf_env

#
To run the curated BRCA1 valdiation set:

python3 vep_annotate.py --mode known --known brca1_known.txt --out brca1_validation.tsv
#
To run the first 200 valid variants from the HG01868 chromosome 22 VCF:

python3 vep_annotate.py --mode vcf_firstN --vcf HG01868.chr22.vcf --limit 200 --out chr22_first200.tsv

****

**General usage**

Annotate a curated known variant set:

python3 vep_annotate.py --mode known --known variants.txt --out known_output.tsv
#
Annotate the first N valid variants from a VCF:

python3 vep_annotate.py --mode vcf_firstN --vcf sample.vcf --limit 200 --out sample_output.tsv


## Command-line arguments

| Argument | Description |
|---|---|
| `--mode` | `known` or `vcf_firstN` |
| `--known` | Path to known variants file |
| `--vcf` | Path to input VCF |
| `--limit` | Number of valid variants to process in `vcf_firstN` mode, max 200 |
| `--species` | Species for Ensembl VEP, default `human` |
| `--out` | Output TSV file |
| `--no-alphamissense` | Disable AlphaMissense output |
| `--no-phenotypes` | Disable phenotype output |


## How the pipeline works
The workflow follows these steps:
1. Read variants from either a curated text file or a VCF
2. Remove unsupported or symbolic alleles
3. Convert valid variants into Ensembl VEP /region format
4. Submit variants to the Ensembl REST API
5. Parse the returend JSON response
6. Extract transcript-level annotations
7. Recover AlphaMissense and ClinVae-related fields where available
8. Write the final transcipt-level output to a TSV file


## AlphaMissense and ClinVar extraction
AlphaMissense is extracted from nested trasncript-level JSON fields when present:

transcript_consequences -> alphamissense -> am_class

transcript_consequences -> alphamissense -> am_pathogenicity
#
ClinVar related information is extracted from:
- colocated_variants
- PHENOTYPES

## Notes
- The script limits batch submissions to 200 variants per POST request.
- Only simple DNA alleles are accepted for /region input.
- Symbolic alleles such as `<DEL>`, `<CN0>`, `*`, and `.` are skipped.
- Multi-allelic VCF entries are split into separate variants before submission.
- One output row is written per transcript consequence, not per genomic variant.
- Non-coding or intergenic variants are still retained in the output.

## Outputs generated in this project
**BRCA1 validation output**

Command: 

python3 vep_annotate.py --mode known --known brca1_known.txt --out brca1_validation.tsv


Output file: 

brca1_validation.tsv

This file contains transcript-level annotation for the curated BRCA1 validation variants.

#

**HG01868 chromosome 22 output**

Command: 

python3 vep_annotate.py --mode vcf_firstN --vcf HG01868.chr22.vcf --limit 200 --out chr22_first200.tsv


Output file: 

chr22_first200.tsv

This file contains transcript-level annotation for the first 200 valid variants from the HG01868 chromosome 22 VCF.


## Reproducibility
The workflow is fully script-based and uses the Ensembl REST API directly. Reproducible use depends on:
- keeping the same input files
- using the same --limit value for VCF mode
- documenting the genome build of the input variants
- storing the output TSV files used for downstream analysis


## Limitations
Current limitations include:
- a maximum of 200 variants per POST request
- no automatic chunking of larger VCF files
- no automatic coordinate conversion between genome builds
- no direct ranking or prioritisation of variants
- designed for simple REST-based annotation rather than full production-scale pipelines

## Future improvements
Possible future extensions include:
- chunked submission for larger VCF files
- automated genome build checking
- additional pathogenicity and conservation predictors
- variant prioritisation and ranking
- broader cohort-scale support
- more flexible filtering for coding or disease-relevant variants

## Author
Developed as part of a transcript-level human variant annotation project using the Ensembl VEP REST API.



