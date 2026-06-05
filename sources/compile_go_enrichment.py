import os
import glob
import re
import csv

def compile_enrichment(cohort_path, cell_type, prefix):
    directory = os.path.join(cohort_path, "WGCNA_GO_enrichr", f"GO_enrichr_{cell_type}")
    output_file = os.path.join(directory, f"{cell_type}_GO_BP_enrichment_summary.csv")
    
    # Find all module GO_BP csv files
    pattern = os.path.join(directory, f"{cell_type}_module_*_GO_BP.csv")
    files = glob.glob(pattern)
    
    all_rows = []
    
    for filepath in files:
        filename = os.path.basename(filepath)
        # Match the module number
        match = re.search(rf"{cell_type}_module_(\d+)_GO_BP\.csv", filename)
        if match:
            module_num = int(match.group(1))
            module_name = f"{prefix}{module_num}"
            
            with open(filepath, "r", encoding="utf-8") as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    # Filter for Adjusted P-value <= 0.05
                    try:
                        adj_p = float(row["Adjusted P-value"])
                    except (ValueError, TypeError, KeyError):
                        continue
                    
                    if adj_p <= 0.05:
                        all_rows.append({
                            "module_num": module_num,
                            "Module": module_name,
                            "GO_term": row.get("Term", ""),
                            "Overlap": row.get("Overlap", ""),
                            "P-value": row.get("P-value", ""),
                            "Adjusted P-value": row.get("Adjusted P-value", ""),
                            "Odds Ratio": row.get("Odds Ratio", ""),
                            "Genes": row.get("Genes", "")
                        })
                        
    # Sort by module number first, then by Adjusted P-value ascending
    all_rows.sort(key=lambda x: (x["module_num"], float(x["Adjusted P-value"])))
    
    # Write output CSV
    headers = ["Module", "GO_term", "Overlap", "P-value", "Adjusted P-value", "Odds Ratio", "Genes"]
    with open(output_file, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        for row in all_rows:
            filtered_row = {k: row[k] for k in headers}
            writer.writerow(filtered_row)
            
    print(f"Successfully compiled {len(all_rows)} rows into {output_file}.")

def main():
    cohorts = [
        "/data/home/swkim0523/research/SCZ_mic_ast_corr_diff/data/primary_cohort",
        "/data/home/swkim0523/research/SCZ_mic_ast_corr_diff/data/validation_cohort/GSE254569"
    ]
    
    configs = [
        ("microglia", "Mic"),
        ("astrocyte", "Ast")
    ]
    
    for cohort_path in cohorts:
        cohort_name = os.path.basename(cohort_path) if "GSE" not in cohort_path else "validation_cohort (GSE254569)"
        print(f"\nProcessing {cohort_name}...")
        for cell_type, prefix in configs:
            compile_enrichment(cohort_path, cell_type, prefix)

if __name__ == "__main__":
    main()
