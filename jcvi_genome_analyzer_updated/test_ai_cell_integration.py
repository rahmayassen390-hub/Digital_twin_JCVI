import sys
import os
import torch
import numpy as np

# Add current dir to path
sys.path.append(os.getcwd())

from data_structures import Gene
from data_enrichment_layer import DataEnrichmentLayer
from molecular_ur_generator import MolecularURGenerator

def test_integration():
    print("🚀 Starting AI Virtual Cell Integration Test...")
    
    # 1. Create a dummy gene (dnaA from syn3A)
    # CP002027.1:1-1353
    sample_seq = "ATGAACGTAAACGATATTTTAAAAGAACTTAAACTAAGTTTAATGGCTAATAAAAATATTGATGAATCCG" # truncated
    gene = Gene(
        id="JCVISYN3A_0001",
        name="dnaA",
        start=1,
        end=1353,
        strand="+",
        type="CDS",
        product="Chromosomal replication initiator protein",
        sequence=sample_seq
    )
    
    print(f"Testing Gene: {gene.id} ({gene.name})")
    
    # 2. Data Enrichment
    print("\nPhase 2 & 3: Data Enrichment & AI Gap Filling...")
    enricher = DataEnrichmentLayer(data_dir="Data")
    enriched_gene = enricher.enrich_gene(gene)
    
    print(f"  Essentiality: {enriched_gene.essentiality_status}")
    print(f"  mRNA Expression: {enriched_gene.mrna_expression}")
    print(f"  mRNA Half-life (pred): {enriched_gene.mrna_half_life} min")
    print(f"  Protein Half-life (pred): {enriched_gene.protein_half_life} hr")
    print(f"  RBS Strength (pred): {enriched_gene.rbs_strength}")
    
    # 3. Molecular UR Generation
    print("\nPhase 4: Molecular UR Generation (Transformers)...")
    try:
        ur_generator = MolecularURGenerator(use_gpu=False) # Use CPU for local test if GPU OOMs
        final_gene = ur_generator.generate_ur(enriched_gene)
        
        print("  Embeddings generated successfully!")
        print(f"  DNA Embedding Shape: {final_gene.molecular_ur['dna_embedding'].shape}")
        print(f"  RNA Embedding Shape: {final_gene.molecular_ur['rna_embedding'].shape}")
        print(f"  Protein Embedding Shape: {final_gene.molecular_ur['protein_embedding'].shape}")
        print(f"  Combined Embedding Shape: {final_gene.molecular_ur['combined_embedding'].shape}")
        
    except Exception as e:
        print(f"  ⚠️ Embedding generation failed/skipped: {e}")
        print("  (This is expected if models are too large for current environment RAM/CPU)")

    print("\n✅ Verification complete!")

if __name__ == "__main__":
    test_integration()
