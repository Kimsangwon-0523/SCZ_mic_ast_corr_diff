"""Build the compact references used by the program scoring script."""

import argparse
import gzip
import hashlib
import re
from pathlib import Path

import pandas as pd

from .inputs import read_geo_samples
from .paths import REFERENCE_DIR


def gtf_symbols(path, species):
    """Read gene records and retain only unique symbol to stable ID relationships."""
    rows = []
    with gzip.open(path, "rt") as source:
        for line in source:
            # gene, transcript and exon records share the same annotation table.
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t", 8)
            if fields[2] != "gene":
                continue
            attributes = dict(re.findall(r'(gene_id|gene_name) "([^"]+)"', fields[8]))
            if "gene_id" in attributes and "gene_name" in attributes:
                rows.append((attributes["gene_id"], attributes["gene_name"]))
    keys = [f"{species}_gene_id", f"{species}_symbol"]
    data = pd.DataFrame(rows, columns=keys).drop_duplicates()
    return data.drop_duplicates(subset=keys[0], keep=False).drop_duplicates(subset=keys[1], keep=False)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geo-dir", type=Path, required=True, help="Directory with <GSE>.samples.txt exports")
    parser.add_argument("--ortholog-candidates", type=Path, required=True)
    parser.add_argument("--human-gtf", type=Path, required=True)
    parser.add_argument("--mouse-gtf", type=Path, required=True)
    args = parser.parse_args()
    output = REFERENCE_DIR
    output.mkdir(parents=True, exist_ok=True)
    sources = []
    records = []
    for gse in ("GSE228927", "GSE282665"):
        path = args.geo_dir / f"{gse}.samples.txt"
        sources.append((path, f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=gsm&form=text&view=full"))
        records.append(read_geo_samples(path, gse))
    pd.concat(records, ignore_index=True).to_csv(output / "geo_samples.tsv", sep="\t", index=False, na_rep="NA")

    pairs = pd.read_csv(args.ortholog_candidates)
    selected = pairs.loc[
        pairs["homology_type"].eq("ortholog_one2one") & pairs["is_high_confidence"].eq(1),
        ["human_gene_id", "mouse_gene_id", "homology_type", "is_high_confidence"],
    ].drop_duplicates()
    if selected["human_gene_id"].duplicated().any() or selected["mouse_gene_id"].duplicated().any():
        raise ValueError("High-confidence 1:1 ortholog reference contains conflicting IDs")
    orthology = selected.merge(
        gtf_symbols(args.human_gtf, "human"), validate="one_to_one", on="human_gene_id",
    ).merge(gtf_symbols(args.mouse_gtf, "mouse"), validate="one_to_one", on="mouse_gene_id")
    orthology.to_csv(output / "ortholog_symbols.tsv.gz", sep="\t", index=False,
                     compression={"method": "gzip", "mtime": 0})
    sources.extend([
        (args.ortholog_candidates, "Local datamining_5 Ensembl Compara candidates; see README.md"),
        (args.human_gtf, "Ensembl GRCh38 release 114 GTF from repository data"),
        (args.mouse_gtf, "https://ftp.ensembl.org/pub/release-116/gtf/mus_musculus/Mus_musculus.GRCm39.116.gtf.gz"),
    ])
    pd.DataFrame([
        {"file": path.name, "sha256": sha256(path), "source": source}
        for path, source in sources
    ]).to_csv(output / "provenance.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
