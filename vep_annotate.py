#!/usr/bin/env python3
"""
vep_annotate.py

An Ensembl VEP REST pipeline for two use cases:

1) known
   Annotate a small curated variant set
   Input per line:
     chrom pos id ref alt
   or:
     chrom pos ref alt

2) vcf_firstN
   Annotate the first N valid variants from a VCF (max 200)

Output:
- one TSV row per transcript consequence
- includes consequence, HGVS, SIFT, PolyPhen, AlphaMissense, and ClinVar predictions.
"""

# Allows the use of modern type hints like str | None
from __future__ import annotations

# argparse = read command-line options
# csv = write tabular output
# json = send/receive JSON data to/from the REST API
# time = wait briefly when retrying failed requests
import argparse
import csv
import json
import time

# These are standard-library HTTP tools
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

# Base Ensembl VEP REST URL
VEP_BASE = "https://rest.ensembl.org/vep"

# Ensembl REST allows practical batch submission around this size
MAX_POST_VARIANTS = 200

# Only allow simple DNA bases for region-based submission
VALID_BASES = set("ACGTN")



# HTTP

def add_query_params(url: str, params: dict) -> str:
    """
    Add query parameters to a URL.

    Example:
      base URL: https://rest.ensembl.org/vep/human/region
      params: {"hgvs": 1, "Phenotypes": 1}

    Output:
      https://rest.ensembl.org/vep/human/region?hgvs=1&Phenotypes=1
    """
    # Remove any parameters whose value is None
    params = {k: v for k, v in params.items() if v is not None}

    # If there are no parameters, return the original URL unchanged
    if not params:
        return url

    # Add either ? or & depending on whether the URL already has parameters
    return f"{url}{'&' if '?' in url else '?'}{urlencode(params)}"


def http_post_json(url: str, payload: dict, retries: int = 4, timeout: int = 60):
    """
    Send a JSON POST request and return the JSON response.

    retries:
      Number of attempts for temporary failures

    timeout:
      Max seconds to wait for a response
    """
    # Convert Python dictionary into JSON bytes
    data = json.dumps(payload).encode("utf-8")

    # Build an HTTP POST request
    req = Request(
        url,
        data=data,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )

    # Try the request multiple times if temporary errors occur
    for attempt in range(retries):
        try:
            # Send request and parse response JSON back into Python
            with urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))

        except HTTPError as e:
            # If Ensembl returns 400, show the body to help debugging bad input
            if e.code == 400:
                body = e.read().decode("utf-8", errors="replace")
                raise RuntimeError(
                    f"HTTP 400 for {url}\nResponse body:\n{body}\nPayload preview:\n{str(payload)[:1500]}"
                ) from e

            # Retry on rate-limiting or temporary server errors
            if e.code in (429, 500, 502, 503, 504):
                time.sleep(1.5 ** attempt)
                continue
            raise

        except URLError:
            # Retry on temporary network errors
            time.sleep(1.5 ** attempt)
            continue

    # If all retries fail, stop with an error
    raise RuntimeError(f"POST failed after {retries} retries: {url}")



# Input parsing

def parse_vcf_variants(vcf_path: str):
    """
    Read a VCF file and yield one dictionary per ALT allele.

    Example VCF row:
      1   12345   .   A   G,T

    This becomes two yielded variants:
      {"chrom": "1", "pos": 12345, "ref": "A", "alt": "G"}
      {"chrom": "1", "pos": 12345, "ref": "A", "alt": "T"}
    """
    with open(vcf_path) as f:
        for line in f:
            # Skip metadata and header lines
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue

            # Extract the key VCF columns
            chrom = fields[0].replace("chr", "")
            pos = int(fields[1])
            ref = fields[3].strip()

            # Split multiple ALT alleles into separate variant records
            for alt in fields[4].strip().split(","):
                yield {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt.strip(),
                }


def parse_known_variants(path: str):
    """
    Read a small curated variant file.

    Accept either:
      chrom pos id ref alt
    or:
      chrom pos ref alt

    Returns one dictionary per line.
    """
    with open(path) as f:
        for raw in f:
            line = raw.strip()

            # Skip empty lines or comments
            if not line or line.startswith("#"):
                continue

            parts = line.split()

            # Accept either 5-column or 4-column input
            if len(parts) >= 5:
                chrom, pos, _id, ref, alt = parts[:5]
            elif len(parts) == 4:
                chrom, pos, ref, alt = parts
            else:
                raise ValueError(f"Bad known variant line: {line}")

            yield {
                "chrom": chrom.replace("chr", ""),
                "pos": int(pos),
                "ref": ref.strip(),
                "alt": alt.strip(),
            }



# Variant filtering / formatting

def is_simple_allele(allele: str) -> bool:
    """
    Check whether an allele is simple enough for VEP /region input.

    Allowed:
      A, C, G, T, N
      short strings made from these, e.g. AT, GCA

    Rejected:
      .
      *
      <DEL>
      <CN0>
      symbolic or ambiguous allele strings
    """
    return (
        bool(allele)
        and allele not in {".", "*"}
        and "<" not in allele
        and ">" not in allele
        and all(base in VALID_BASES for base in allele.upper())
    )


def variant_to_region_string(variant: dict) -> str | None:
    """
    Convert one variant into Ensembl VEP /region format:
      chrom start end ref/alt

    Example:
      {"chrom": "17", "pos": 43057117, "ref": "C", "alt": "G"}
    becomes:
      "17 43057117 43057117 C/G"
    """
    ref = variant["ref"]
    alt = variant["alt"]

    # If ref or alt is not a simple allele, skip it
    if not (is_simple_allele(ref) and is_simple_allele(alt)):
        return None

    pos = int(variant["pos"])
    end = pos + len(ref) - 1
    return f"{variant['chrom']} {pos} {end} {ref}/{alt}"


def take_first_n_valid(variants, n: int):
    """
    Keep the first N variants that can be converted into valid VEP /region strings.

    Returns:
      selected = list of valid variants
      skipped = how many invalid variants were ignored
    """
    selected = []
    skipped = 0

    for variant in variants:
        if variant_to_region_string(variant) is None:
            skipped += 1
            continue

        selected.append(variant)
        if len(selected) >= n:
            break

    return selected, skipped



# VEP REST call

def vep_post_region(
    variants: list[dict],
    species: str = "human",
    include_alphamissense: bool = True,
    include_phenotypes: bool = True,
):
    """
    Submit a batch of variants to Ensembl VEP /region.

    Returns:
      records = JSON response from Ensembl
      skipped = how many variants were invalid at submission stage
      len(region_strings) = number of variants actually sent
    """
    # Do not exceed the practical Ensembl REST batch size
    if len(variants) > MAX_POST_VARIANTS:
        raise ValueError(f"Too many variants for one POST ({len(variants)}). Max {MAX_POST_VARIANTS}.")

    # Build endpoint URL with useful annotation options
    url = add_query_params(
        f"{VEP_BASE}/{species}/region",
        {
            "AlphaMissense": 1 if include_alphamissense else 0,
            "Phenotypes": 1 if include_phenotypes else 0,
            "hgvs": 1,
            "vcf_string": 1,
        },
    )

    region_strings = []
    skipped = 0

    # Convert variants into Ensembl /region strings
    for variant in variants:
        region = variant_to_region_string(variant)
        if region is None:
            skipped += 1
        else:
            region_strings.append(region)

    if not region_strings:
        raise ValueError("No valid variants to send.")

    # Send POST request to Ensembl
    records = http_post_json(url, {"variants": region_strings})
    return records, skipped, len(region_strings)



# Response parsing

def first_present(d: dict, keys: list[str], default=""):
    """
    Return the first non-None value found from a list of possible keys.

    Useful when APIs may use slightly different field names.
    """
    for key in keys:
        value = d.get(key)
        if value is not None:
            return value
    return default


def extract_alphamissense(tc: dict):
    """
    Extract AlphaMissense values from one transcript consequence.

    Typical location:
      tc["alphamissense"]["am_class"]
      tc["alphamissense"]["am_pathogenicity"]

    Also checks fallback field names just in case.
    """
    nested = tc.get("alphamissense")
    if isinstance(nested, dict):
        score = nested.get("am_pathogenicity", "")
        pred = nested.get("am_class", "")
        if score != "" or pred != "":
            return score, pred

    # Fallback for alternative names
    score = first_present(tc, ["am_pathogenicity", "alpha_missense_score", "alphamissense_score"], "")
    pred = first_present(tc, ["am_class", "alpha_missense_prediction", "alphamissense_prediction"], "")
    return score, pred


def extract_clinvar_from_colocated(record: dict):
    """
    Extract ClinVar-like information from Ensembl 'colocated_variants'.

    Returns:
      clin_sigs = comma-separated list of clinical significance values
      ids = comma-separated list of colocated variant IDs
    """
    clin_sigs = set()
    ids = set()

    for cv in record.get("colocated_variants") or []:
        if cv.get("id"):
            ids.add(str(cv["id"]))

        clin_sig = cv.get("clin_sig")
        if isinstance(clin_sig, list):
            clin_sigs.update(str(x) for x in clin_sig if x)
        elif clin_sig:
            clin_sigs.add(str(clin_sig))

    return ",".join(sorted(clin_sigs)), ",".join(sorted(ids))


def extract_clinvar_from_phenotypes(record: dict):
    """
    Extract phenotype entries where the source looks like ClinVar.

    Returns one pipe-separated string.
    """
    entries = []

    for p in record.get("PHENOTYPES") or []:
        if not isinstance(p, dict):
            continue

        source = str(p.get("source", ""))
        if "clinvar" in source.lower():
            entries.append(f"{source}:{p.get('phenotype', '')}:{p.get('id', '')}".strip(":"))

    return "|".join(entries)


def parse_vep_record_to_rows(record: dict):
    """
    Convert one Ensembl VEP response record into one or more output rows.

    Important:
      one genomic variant may have multiple transcript consequences,
      so this function returns one TSV row per transcript consequence.
    """
    clin_sig, colocated_ids = extract_clinvar_from_colocated(record)
    clinvar_pheno = extract_clinvar_from_phenotypes(record)

    # Fields shared by every transcript row for this variant
    common = {
        "input": record.get("input", ""),
        "chrom": str(record.get("seq_region_name", "")),
        "pos": str(record.get("start", "")),
        "allele_string": record.get("allele_string", ""),
        "most_severe_consequence": record.get("most_severe_consequence", ""),
        "hgvsg": record.get("hgvsg", ""),
        "clinvar_clin_sig": clin_sig,
        "colocated_variant_ids": colocated_ids,
        "clinvar_phenotypes": clinvar_pheno,
    }

    transcript_consequences = record.get("transcript_consequences") or []

    # If no transcript-level entries exist, still output one minimal row
    if not transcript_consequences:
        return [{
            **common,
            "consequence_terms": "",
            "impact": "",
            "gene_id": "",
            "gene_symbol": "",
            "transcript_id": "",
            "gene_biotype": "",
            "hgvsc": "",
            "hgvsp": "",
            "sift_prediction": "",
            "sift_score": "",
            "polyphen_prediction": "",
            "polyphen_score": "",
            "alphamissense_score": "",
            "alphamissense_prediction": "",
        }]

    rows = []

    # Otherwise create one row per transcript consequence
    for tc in transcript_consequences:
        am_score, am_pred = extract_alphamissense(tc)

        rows.append({
            **common,
            "consequence_terms": ",".join(tc.get("consequence_terms") or []),
            "impact": tc.get("impact", ""),
            "gene_id": tc.get("gene_id", ""),
            "gene_symbol": tc.get("gene_symbol", ""),
            "transcript_id": tc.get("transcript_id", ""),
            "gene_biotype": tc.get("biotype", ""),
            "hgvsc": tc.get("hgvsc", ""),
            "hgvsp": tc.get("hgvsp", ""),
            "sift_prediction": tc.get("sift_prediction", ""),
            "sift_score": tc.get("sift_score", ""),
            "polyphen_prediction": tc.get("polyphen_prediction", ""),
            "polyphen_score": tc.get("polyphen_score", ""),
            "alphamissense_score": str(am_score),
            "alphamissense_prediction": str(am_pred),
        })

    return rows



# Output

def write_tsv(rows: list[dict], out_path: str):
    """
    Write parsed annotation rows to a TSV file.
    """
    if not rows:
        raise ValueError("No rows to write.")

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)



# CLI (Command Line Interface)

def parse_args():
    """
    Read command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Annotate variants with Ensembl VEP REST.")
    parser.add_argument("--mode", choices=["known", "vcf_firstN"], required=True)
    parser.add_argument("--known", help="Path to known variants file")
    parser.add_argument("--vcf", help="Path to input VCF")
    parser.add_argument("--limit", type=int, default=200, help="Number of variants for vcf_firstN (max 200)")
    parser.add_argument("--species", default="human")
    parser.add_argument("--out", required=True)
    parser.add_argument("--no-alphamissense", action="store_true")
    parser.add_argument("--no-phenotypes", action="store_true")
    return parser.parse_args()


def load_variants(args):
    """
    Load variants according to the selected mode.

    known mode:
      reads a small curated list

    vcf_firstN mode:
      reads a VCF and keeps the first N valid variants
    """
    if args.mode == "known":
        if not args.known:
            raise SystemExit("known mode requires --known")

        variants = list(parse_known_variants(args.known))
        if len(variants) > MAX_POST_VARIANTS:
            raise SystemExit(f"known input has >{MAX_POST_VARIANTS} variants.")
        return variants

    if not args.vcf:
        raise SystemExit("vcf_firstN mode requires --vcf")
    if args.limit > MAX_POST_VARIANTS:
        raise SystemExit(f"--limit must be <= {MAX_POST_VARIANTS}.")

    variants, skipped = take_first_n_valid(parse_vcf_variants(args.vcf), args.limit)
    if skipped:
        print(f"[warn] Skipped {skipped} unsupported alleles while collecting first {args.limit} valid variants.")
    return variants


def main():
    """
    Main program order:
      1. read command-line arguments
      2. load variants
      3. send variants to Ensembl VEP
      4. parse JSON into transcript-level rows
      5. write TSV output
    """
    args = parse_args()
    variants = load_variants(args)

    vep_records, skipped, sent = vep_post_region(
        variants,
        species=args.species,
        include_alphamissense=not args.no_alphamissense,
        include_phenotypes=not args.no_phenotypes,
    )

    if skipped:
        print(f"[warn] Skipped {skipped} variants at POST time (symbolic/invalid alleles).")

    rows = []

    # Parse each Ensembl record into output rows
    for record in vep_records:
        rows.extend(parse_vep_record_to_rows(record))

    # Write final TSV
    write_tsv(rows, args.out)

    # Print summary to terminal
    print(f"Done. Wrote TSV: {args.out}")
    print(f"Sent {sent} variants to VEP.")
    print(f"Wrote {len(rows)} transcript-level rows.")


# Run the script only when called directly
if __name__ == "__main__":
    main()