import os
import glob
import re
import csv

def process_modules(cell_type, prefix):
    base_dir = f"../../data/primary_cohort/WGCNA_modules/modules_{cell_type}"
    output_file = os.path.join(base_dir, f"{cell_type}_module_info.csv")
    
    # Find all module text files
    pattern = os.path.join(base_dir, f"{cell_type}_module_*.txt")
    files = glob.glob(pattern)
    
    module_data = []
    
    for filepath in files:
        filename = os.path.basename(filepath)
        # Match the module number
        match = re.search(rf"{cell_type}_module_(\d+)\.txt", filename)
        if match:
            module_num = int(match.group(1))
            module_name = f"{prefix}{module_num}"
            
            # Read genes from file
            with open(filepath, "r", encoding="utf-8") as f:
                genes = [line.strip() for line in f if line.strip()]
                
            # Format as [GeneA, GeneB, ...]
            genes_str = f"[{', '.join(genes)}]"
            
            module_data.append((module_num, module_name, genes_str))
            
    # Sort modules by number
    module_data.sort(key=lambda x: x[0])
    
    # Write to CSV
    with open(output_file, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Module", "Member_genes"])
        for _, module_name, genes_str in module_data:
            writer.writerow([module_name, genes_str])
            
    print(f"Successfully created {output_file} with {len(module_data)} modules for {cell_type}.")

def main():
    configs = [
        ("microglia", "Mic"),
        ("astrocyte", "Ast")
    ]
    for cell_type, prefix in configs:
        process_modules(cell_type, prefix)

if __name__ == "__main__":
    main()
