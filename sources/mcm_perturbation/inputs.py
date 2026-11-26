"""Read the downloaded expression matrices and match columns to GEO."""

import re
from pathlib import Path

import pandas as pd


RECIPIENT_FILES = {
    "astrocyte": "GSE282665_astrocyte.GSE282665_microglia.media.salmon.merged.gene_tpm.tsv",
    "motor neuron": "GSE282665_motor.neuron.GSE282665_microglia.media.salmon.merged.gene_tpm.tsv",
}


def read_geo_samples(path, dataset):
    """Extract titles and characteristics from a GEO sample SOFT export."""
    rows = []
    for record in Path(path).read_text().split("^SAMPLE = ")[1:]:
        gsm, *lines = record.splitlines()
        fields = {"dataset": dataset, "gsm_id": gsm}
        for line in lines:
            key, _, value = line.partition(" = ")
            if key == "!Sample_title":
                fields["geo_title"] = value
            elif key == "!Sample_characteristics_ch1":
                field, _, text = value.partition(": ")
                fields[field.lower().replace(" ", "_")] = text
        rows.append(fields)
    if not rows:
        raise ValueError(f"No GEO sample records in {path}")
    return pd.DataFrame(rows)


def read_mouse_fpkm(data_root):
    path = data_root / "GSE228927/GSE228927_FPKM_anno.xls"
    raw = pd.read_csv(path, sep="\t", low_memory=False).set_index("id")
    return raw.loc[:, raw.columns.str.fullmatch(r"(NOT|LPS|Hemin)\d+")]


def read_human_tpm(data_root, recipient):
    path = data_root / "GSE282665" / RECIPIENT_FILES[recipient]
    raw = pd.read_csv(path, sep="\t", keep_default_na=False).set_index("gene_id")
    if raw.index.has_duplicates or (raw.index == "").any():
        raise ValueError(f"Missing or duplicate Ensembl gene ID in {path.name}")
    return raw.drop(columns="gene_name"), raw["gene_name"]


def sample_manifest(geo, dataset, recipient, columns):
    """Assign conditions from GEO with explicit rules for each worksheet."""
    geo = geo.loc[geo["dataset"].eq(dataset)]
    if dataset == "GSE282665":
        cell_type = {
            "astrocyte": "microglia_media_on_ips_astrocytes",
            "motor neuron": "microglia_media_on_ips_motor_neurons",
        }[recipient]
        geo = geo.loc[geo["cell_type"].eq(cell_type)]
    records = []
    for sample in columns:
        pair_id = None
        recipient_line = None
        if dataset == "GSE228927":
            match = re.fullmatch(r"(NOT|LPS|Hemin)(\d+)", sample)
            if match is None:
                raise ValueError(f"Unexpected GSE228927 expression column: {sample}")
            medium, replicate = match.groups()
            condition = {"NOT": "Resting-MCM", "LPS": "LPS-MCM", "Hemin": "Hemin-MCM"}[medium]
            title = f"{condition}, replicate {replicate}"
            planned = True
        else:
            title = sample
        matches = geo.loc[geo["geo_title"].eq(title)]
        if len(matches) != 1:
            raise ValueError(f"{dataset} {sample}: expected one GEO sample, found {len(matches)}")
        record = matches.iloc[0]
        if dataset == "GSE282665":
            recipient_line = record["cell_line"]
            source = record["microglia_media"]
            stimulus = record["lps_treatment"]
            condition = {
                "ctrl": {"untx": "Resting-MCM", "lps": "LPS-MCM"},
                "nil": {"untx": "Fresh", "lps": "Fresh+LPS"},
                "vcp": {"untx": "Resting VCP-MCM", "lps": "LPS VCP-MCM"},
            }[source][stimulus]
            planned = source == "ctrl" and record["genotype"] == "ctrl"
            pair_id = recipient_line if planned else None
        records.append({
            "dataset": dataset, "recipient": recipient, "sample_id": sample,
            "gsm_id": record["gsm_id"], "geo_title": title, "condition": condition,
            "recipient_line": recipient_line, "recipient_genotype": record.get("genotype"),
            "mcm_source_genotype": record.get("microglia_media"), "pair_id": pair_id,
            "planned_sample": planned,
        })
    manifest = pd.DataFrame(records).set_index("sample_id", drop=False)
    if not manifest.index.is_unique or not manifest["gsm_id"].is_unique:
        raise ValueError(f"Duplicate samples for {dataset} {recipient}")
    if dataset == "GSE282665":
        paired = manifest.loc[manifest["planned_sample"]]
        pair_counts = pd.crosstab(paired["pair_id"], paired["condition"])
        if set(pair_counts) != {"Resting-MCM", "LPS-MCM"} or not pair_counts.eq(1).all().all():
            raise ValueError(f"Incomplete or duplicate GEO-confirmed pairs for {recipient}")
    return manifest
