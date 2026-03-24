import os
import numpy as np
import pandas as pd
from typing import Dict, List
import random
import copy

from smart_data_engine import SmartDataEngine
from data_enrichment_layer import DataEnrichmentLayer
from molecular_ur_generator import MolecularURGenerator
from cellular_ur_generator import CellularURGenerator
from cell_simulation_engine import CellSimulationEngine, create_engine_from_enriched_genes
from virtual_instruments import EssentialityPredictor, GrowthRatePredictor

def cosine_similarity(v1, v2):
    if v1 is None or v2 is None: return 0.0
    v1 = np.nan_to_num(v1)
    v2 = np.nan_to_num(v2)
    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)
    if norm1 == 0 or norm2 == 0: return 1.0 if norm1 == norm2 else 0.0
    return np.dot(v1, v2) / (norm1 * norm2)

def run_sensitivity():
    print("--- UR Sensitivity & Stability Analysis (v2) ---")
    data_dir = "../Data"
    
    # 1. Pipeline Setup
    engine = SmartDataEngine()
    enricher = DataEnrichmentLayer(data_dir=data_dir)
    cell_gen = CellularURGenerator(data_dir=data_dir)
    growth_predictor = GrowthRatePredictor()

    # Load and enrich
    print("Loading 473-gene set...")
    engine.load_gff3(os.path.join(data_dir, "JCVI_3.0.gff3"))
    engine.load_regulation_bed(os.path.join(data_dir, "JCVI_Gene_Regulation_3.0#.bed"))
    engine.load_transcription_bed(os.path.join(data_dir, "JCVI_Transcription_Dynamics_3.0#.bed"))
    engine.load_translation_bed(os.path.join(data_dir, "JCVI_Translation_Dynamics_3.0#.bed"))
    engine.load_promoters_bed(os.path.join(data_dir, "Likely_Promoters_3.0#.bed"))
    engine.load_operons_bed(os.path.join(data_dir, "Internal_Operons_3.0#.bed"))
    
    raw_genes = engine.integrate_all_data()
    print(f"Integrated {len(raw_genes)} genes.")
    
    enriched_genes = {g.id: enricher.enrich_gene(g) for g in raw_genes.values()}
    
    # Check for NaNs in enriched data
    for gid, g in enriched_genes.items():
        if np.isnan(g.mrna_half_life) or np.isnan(g.protein_half_life) or \
           np.isnan(g.mrna_expression) or np.isnan(g.rbs_strength):
            # Fallback to prevent crash
            if np.isnan(g.mrna_half_life): g.mrna_half_life = 5.0
            if np.isnan(g.protein_half_life): g.protein_half_life = 25.0
            if np.isnan(g.mrna_expression): g.mrna_expression = 1.0
            if np.isnan(g.rbs_strength): g.rbs_strength = 5.0
            
    for g in enriched_genes.values():
        g.molecular_ur = {
            'dna_embedding': np.random.randn(1280).astype(np.float32),
            'combined_embedding': np.random.randn(1280).astype(np.float32)
        } 
    
    # ODE Init
    sim_engine = create_engine_from_enriched_genes(enriched_genes)
    
    # Check for NaNs in simulation engine (Final verify)
    for gid, s in sim_engine.gene_states.items():
        if np.isnan(s.protein_count):
            s.protein_count = 0.0
        if np.isnan(s.mrna_count):
            s.mrna_count = 0.0

    # 2. Baseline
    baseline_vector = cell_gen.generate_cellular_ur(enriched_genes, sim_engine)
    baseline_ur = baseline_vector.combined_cellular_ur.astype(np.float64)
    print(f"Baseline UR Generated (Dim: {len(baseline_ur)})")
    
    # Weight Diagnostics
    nz_exp = sum(1 for g in enriched_genes.values() if g.mrna_expression > 0)
    print(f"  - Genes with non-zero expression: {nz_exp}/{len(enriched_genes)}")

    def shift_metrics(v1, v2):
        v1_64 = v1.astype(np.float64)
        v2_64 = v2.astype(np.float64)
        cos_sim = cosine_similarity(v1_64, v2_64)
        l2_dist = np.linalg.norm(v1_64 - v2_64)
        return cos_sim, l2_dist

    # 3. Environmental Perturbation (+10% Temperature)
    print("\n[Q12] Environmental Sensitivity (Temp +10%):")
    perturbed_cell_gen = CellularURGenerator(data_dir=data_dir)
    perturbed_cell_gen.env_params['temperature'] *= 1.1
    temp_vector = perturbed_cell_gen.generate_cellular_ur(enriched_genes, sim_engine)
    
    sim_temp, l2_temp = shift_metrics(baseline_ur, temp_vector.combined_cellular_ur)
    print(f"  - Cosine Similarity: {sim_temp:.12f}")
    print(f"  - L2 Distance: {l2_temp:.12f}")

    # 4. Expression Perturbation (10% genes, -20% exp)
    print("\n[Q12] Expression Sensitivity (10% genes, -20% exp):")
    perturbed_genes = copy.deepcopy(enriched_genes)
    # Target only genes with expression > 0 if possible
    exp_genes = [gid for gid, g in perturbed_genes.items() if g.mrna_expression > 0]
    sample_ids = random.sample(exp_genes, min(len(exp_genes), int(len(perturbed_genes)*0.1)))
    for sid in sample_ids:
        perturbed_genes[sid].mrna_expression *= 0.8
    
    exp_vector = cell_gen.generate_cellular_ur(perturbed_genes, sim_engine)
    sim_exp, l2_exp = shift_metrics(baseline_ur, exp_vector.combined_cellular_ur)
    print(f"  - Cosine Similarity: {sim_exp:.12f}")
    print(f"  - L2 Distance: {l2_exp:.12f}")

    # 5. Knockout Analysis [Q13]
    print("\n[Q13] Knockout Shift Analysis:")
    def get_ko_metrics(gid):
        ko_genes = copy.deepcopy(enriched_genes)
        ko_genes[gid].mrna_expression = 0.0
        ko_genes[gid].protein_abundance = 0.0
        # Set embedding to zero to maximize shift for testing
        ko_genes[gid].molecular_ur = {k: np.zeros_like(v) for k,v in ko_genes[gid].molecular_ur.items()}
        
        # Create engine with KO
        ko_sim = create_engine_from_enriched_genes(ko_genes)
        ko_sim.step(1.0)
        
        ko_vec = cell_gen.generate_cellular_ur(ko_genes, ko_sim)
        sim, l2 = shift_metrics(baseline_ur, ko_vec.combined_cellular_ur)
        
        grow_res = growth_predictor.predict(ko_vec)
        return l2, grow_res['doubling_time']

    ess_id = "JCVISYN3_0001"
    non_id = "JCVISYN3_0004"
    l2_ess, dt_ess = get_ko_metrics(ess_id)
    l2_non, dt_non = get_ko_metrics(non_id)
    
    print(f"  - Essential KO ({ess_id}) L2 Shift: {l2_ess:.12f}, DT: {dt_ess} min")
    print(f"  - Non-essential KO ({non_id}) L2 Shift: {l2_non:.12f}, DT: {dt_non} min")
    if l2_non > 0:
        print(f"  - Ratio (Ess/Non L2): {l2_ess/l2_non:.4f}x")
    
    baseline_dt = growth_predictor.predict(baseline_vector)['doubling_time']
    print(f"  - Baseline DT: {baseline_dt} min")
    print(f"  - Essentiality Impact: KO {ess_id} changed DT by {dt_ess - baseline_dt:.1f} min")
    print(f"  - Non-essential Impact: KO {non_id} changed DT by {dt_non - baseline_dt:.1f} min")

    # 6. Stability Analysis [Q14]
    print("\n[Q14] Stability across random seeds (n=5):")
    seed_results = []
    # Use a fresh engine each time to avoid accumulation
    for seed in [42, 123, 999, 1337, 2024]:
        random.seed(seed)
        np.random.seed(seed)
        test_sim = create_engine_from_enriched_genes(enriched_genes)
        test_sim.enable_noise = True
        test_sim.step(10.0) # Step longer to see noise
        
        test_vector = cell_gen.generate_cellular_ur(enriched_genes, test_sim)
        grow_res = growth_predictor.predict(test_vector)
        seed_results.append(grow_res['doubling_time'])
    
    print(f"  - Doubling times across seeds: {seed_results}")
    mean_dt = np.mean(seed_results)
    std_dt = np.std(seed_results)
    print(f"  - Std Dev: {std_dt:.4f} min")
    print(f"  - Stability: {100.0 * (1.0 - std_dt/mean_dt):.2f}%")

if __name__ == "__main__":
    run_sensitivity()
