import os
import re
import xml.etree.ElementTree as ET

def count_stats():
    data_dir = "../Data"
    
    # 1. Gene Count (GFF)
    gff_path = os.path.join(data_dir, "JCVI_3.0.gff3")
    gene_count = 0
    cds_count = 0
    if os.path.exists(gff_path):
        with open(gff_path, 'r') as f:
            for line in f:
                if "	gene	" in line:
                    gene_count += 1
                if "	CDS	" in line:
                    cds_count += 1
    
    # 2. Metabolic Model Stats (XML)
    xml_path = os.path.join(data_dir, "metabolic_model_iMB155.xml")
    rxn_count = 0
    metab_count = 0
    mb_genes = set()
    if os.path.exists(xml_path):
        import xml.etree.ElementTree as ET
        tree = ET.parse(xml_path)
        root = tree.getroot()
        # Find all tags ignoring namespaces
        all_tags = [elem.tag for elem in root.iter()]
        
        # Count reactions
        rxn_count = sum(1 for tag in all_tags if tag.endswith('reaction'))
        # Count species
        metab_count = sum(1 for tag in all_tags if tag.endswith('species'))
        # Count genes (fbc:geneProduct)
        for elem in root.iter():
            if elem.tag.endswith('geneProduct'):
                label = elem.get('label')
                if label:
                    mb_genes.add(label)

    # 3. Excel Model Stats
    xlsx_path = os.path.join(data_dir, "metabolic_reconstruction.xlsx")
    xl_rxns = 0
    xl_metabs = 0
    if os.path.exists(xlsx_path):
        import pandas as pd
        xls = pd.ExcelFile(xlsx_path)
        if 'Reactions' in xls.sheet_names:
            xl_rxns = len(pd.read_excel(xls, 'Reactions'))
        if 'Metabolites' in xls.sheet_names:
            xl_metabs = len(pd.read_excel(xls, 'Metabolites'))

    # 4. Code Stats
    sloc = 0
    funcs = 0
    classes = 0
    for root_dir, _, files in os.walk("."):
        for file in files:
            if file.endswith(".py"):
                path = os.path.join(root_dir, file)
                with open(path, 'r') as f:
                    content = f.read()
                    sloc += len(content.splitlines())
                    funcs += content.count("def ")
                    classes += content.count("class ")

    print(f"Gene Count (GFF): {gene_count}")
    print(f"CDS Count (GFF): {cds_count}")
    print(f"Reaction Count (iMB155): {rxn_count}")
    print(f"Metabolite Count (iMB155): {metab_count}")
    print(f"Genes in iMB155: {len(mb_genes)}")
    print(f"Code - SLOC: {sloc}")
    print(f"Code - Functions: {funcs}")
    print(f"Code - Classes: {classes}")

if __name__ == "__main__":
    count_stats()
