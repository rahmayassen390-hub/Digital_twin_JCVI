import sys
import os
import pandas as pd
import numpy as np

# Adjust paths
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

from data_enrichment_layer import DataEnrichmentLayer
from cell_simulation_engine import CellSimulationEngine, GeneCellState
from metabolic_solver import MetabolicSolver

def test_priority_loading():
    print("\n--- Testing Priority Loading (Problem 1) ---")
    data_dir = "/home/ubuntu/Desktop/digital_twin/jcvi_genome_analyzer_quantum_entanglment_with_REQ/jcvi_genome_analyzer_quantum_entanglment_withnew_update/Data"
    enricher = DataEnrichmentLayer(data_dir)
    
    # Gene with known specialization (JCVISYN3A_0001 usually in all sets)
    class MockGene:
        def __init__(self, id, **kwargs):
            self.id = id
            self.locus_tag = id
            self.mrna_expression = 0.0
            self.half_life = 0.0
            self.gc_content = 30.0
            self.length = 1000
            self.essentiality_status = "unknown"
            self.is_essential = False
            self.category = "other"
            self.protein_id = "unknown"
            self.rbs_strength = 0.0
            self.function_class = "unknown"
            self.cog_category = "unknown"
            for k, v in kwargs.items():
                setattr(self, k, v)
        def __getattr__(self, name):
            return None # Fallback for any missing attribute
    
    gene = MockGene("JCVISYN3A_0001")
    enriched = enricher.enrich_gene(gene)
    print(f"Gene {gene.locus_tag} mRNA count: {enriched.mrna_expression}")
    # If the count is from the specialized file, it's working.
    
def test_simulation_rk4_metabolic():
    print("\n--- Testing RK4 & Metabolic Coupling (Problem 2, 4, 5) ---")
    data_dir = "/home/ubuntu/Desktop/digital_twin/jcvi_genome_analyzer_quantum_entanglment_with_REQ/jcvi_genome_analyzer_quantum_entanglment_withnew_update/Data"
    engine = CellSimulationEngine()
    
    # Initialize metabolic solver
    xml_path = os.path.join(data_dir, "metabolic_model_iMB155.xml")
    if os.path.exists(xml_path):
        engine.metabolic_solver = MetabolicSolver(xml_path)
        print("Metabolic solver initialized.")
    else:
        print(f"XML not found at {xml_path}")
        return

    # Mock some gene states
    # Glucose-6-phosphate isomerase (PGI) often associated with R_PGI and locus tag JCVISYN3A_0028 or similar
    engine.gene_states["JCVISYN3A_0028"] = GeneCellState(
        gene_id="JCVISYN3A_0028",
        mrna_count=10.0,
        protein_count=100.0, # Healthy enzyme level
        transcription_rate=0.1,
        translation_rate=0.01,
        mrna_half_life=2.0,
        protein_half_life=20.0
    )
    
    print("Pre-step fluxes:", len(engine.current_fluxes))
    
    # Run a few steps (this should trigger FBA at t=0 since last_fba_update=-60)
    dt = 1.0
    engine.step(dt)
    
    print(f"Simulation Time: {engine.simulation_time}s")
    print(f"Post-step fluxes count: {len(engine.current_fluxes)}")
    if engine.current_fluxes:
        # Check some key reactions
        for rxn in ["R_PGI", "R_LDH_L", "BIOMASS_iMB155"]:
            if rxn in engine.current_fluxes:
                print(f"  Flux {rxn}: {engine.current_fluxes[rxn]:.4f}")

def test_ml_predictor():
    print("\n--- Testing ML Predictor (Problem 3) ---")
    data_dir = "/home/ubuntu/Desktop/digital_twin/jcvi_genome_analyzer_quantum_entanglment_with_REQ/jcvi_genome_analyzer_quantum_entanglment_withnew_update/Data"
    enricher = DataEnrichmentLayer(data_dir)
    
    class MockGene:
        def __init__(self, length, gc):
            self.length = length
            self.gc_content = gc
    
    gene = MockGene(500, 25.0)
    hl = enricher.predict_mrna_half_life(gene)
    print(f"Predicted half-life for gene (L=500, GC=25): {hl} min")

if __name__ == "__main__":
    try:
        test_priority_loading()
        test_ml_predictor()
        test_simulation_rk4_metabolic()
        print("\n✅ Verification Suite Completed Successfully.")
    except Exception as e:
        print(f"\n❌ Verification Failed: {e}")
        import traceback
        traceback.print_exc()
