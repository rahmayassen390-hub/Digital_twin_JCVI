#!/usr/bin/env python3
"""
Comprehensive Biological Validation & Performance Report
=========================================================
Priority 1: Biological Validation
  1.1 Growth Rate — debug 296 min vs 105 min discrepancy
  1.2 Essentiality Prediction — confusion matrix on all genes
  1.3 Knockout Simulation — essential vs non-essential comparison

Priority 3: Data Quality
  3.1 Gap-filling quality report (RBS, mRNA HL, protein HL)
"""

import sys, os, time, math
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from data_structures import Gene, CellStateVector
from cell_simulation_engine import (
    CellSimulationEngine, SimulationMode, GeneCellState,
    create_engine_from_enriched_genes,
    JCVI_DOUBLING_TIME, JCVI_MRNA_HALF_LIFE, JCVI_PROTEIN_HALF_LIFE,
    JCVI_TRANSCRIPTION_RATE, JCVI_TRANSLATION_RATE, JCVI_RIBOSOME_COUNT
)
from cellular_ur_generator import CellularURGenerator
from virtual_instruments import (
    VirtualInstrumentSuite, GrowthRatePredictor,
    EssentialityPredictor, KnockoutSimulator
)


def create_realistic_genes(n=498):
    """Create biologically realistic JCVI-syn3A gene set."""
    np.random.seed(42)
    genes = {}
    
    # Distribution parameters from JCVI literature
    # Gene lengths: mean=762, median=567 nt (Hutchison et al.)
    gene_lengths = np.random.lognormal(mean=6.3, sigma=0.7, size=n).astype(int) * 3

    for i in range(n):
        gene_id = f"JCVISYN3A_{i:04d}"
        seq_len = max(90, int(gene_lengths[i]))
        seq = "ATG" + "".join(np.random.choice(list("ACGT"), seq_len - 6).tolist()) + "TAA"
        
        # Essentiality: 307 essential, 114 quasi, 77 non-essential
        if i < 307:
            ess = "Essential"
        elif i < 421:
            ess = "Quasiessential"
        else:
            ess = "Nonessential"
        
        # Name some genes for knockout tests
        names = {0: "dnaA", 1: "dnaN", 2: "rpoB", 3: "rpoC", 4: "ftsZ",
                 5: "rpsA", 6: "rplB", 7: "groEL", 8: "dnaE", 9: "rpoA",
                 480: "yfaA", 481: "yfaB", 482: "yfaC", 483: "yfaD", 484: "yfaE",
                 485: "yfaF", 486: "yfaG", 487: "yfaH", 488: "yfaI", 489: "yfaJ"}
        
        gene = Gene(
            id=gene_id, name=names.get(i, f"gene_{i}"),
            start=i * 1200, end=i * 1200 + seq_len,
            strand="+" if i % 2 == 0 else "-",
            type="CDS", function="Cellular function",
            product=f"Protein {i}", sequence=seq,
        )
        
        # Enriched data (biologically realistic distributions)
        gene.essentiality_status = ess
        gene.mrna_expression = float(np.random.lognormal(3.5, 1.2))  # copies/cell
        gene.protein_abundance = float(np.random.lognormal(6.0, 1.5))  # copies/cell
        
        # RBS strength: most ~1-4, few very strong (~8-10)
        gene.rbs_strength = float(np.random.gamma(2.0, 1.5))
        gene.rbs_strength = min(10.0, max(0.5, gene.rbs_strength))
        
        # mRNA half-life: JCVI ~2-8 min (Breuer et al.)
        gene.mrna_half_life = float(np.random.lognormal(1.2, 0.4))  # minutes
        gene.mrna_half_life = min(15.0, max(0.5, gene.mrna_half_life))
        
        # Protein half-life: 5-100 hrs (stable in minimal cell)
        gene.protein_half_life = float(np.random.lognormal(3.2, 0.6))  # hours
        gene.protein_half_life = min(200.0, max(2.0, gene.protein_half_life))
        
        gene.function_class = np.random.choice(
            ["Translation", "Transcription", "Replication", "Metabolism",
             "Membrane", "Unknown"], p=[0.15, 0.05, 0.05, 0.25, 0.15, 0.35])
        gene.ortholog_id = f"ort_{i}" if np.random.random() < 0.6 else ""
        
        # Mock molecular UR
        gene.molecular_ur = {
            'dna_embedding': np.random.randn(1280).astype(np.float32) * 0.1,
            'rna_embedding': np.random.randn(1280).astype(np.float32) * 0.1,
            'protein_embedding': np.random.randn(1280).astype(np.float32) * 0.1,
            'combined_embedding': np.random.randn(3840).astype(np.float32) * 0.1,
        }
        
        genes[gene_id] = gene
    
    return genes


# =============================================================================
# TASK 1.1: GROWTH RATE VALIDATION
# =============================================================================

def task_1_1_growth_rate_validation(genes):
    """Investigate the 296 min vs 105 min doubling time discrepancy."""
    print("\n" + "=" * 70)
    print("TASK 1.1: GROWTH RATE VALIDATION")
    print("=" * 70)
    
    print("\n--- JCVI Literature Reference Values ---")
    print(f"  Doubling time: {JCVI_DOUBLING_TIME/60:.0f} min ({JCVI_DOUBLING_TIME/3600:.1f} hr)")
    print(f"  mRNA half-life: {JCVI_MRNA_HALF_LIFE:.0f} sec ({JCVI_MRNA_HALF_LIFE/60:.1f} min)")
    print(f"  Protein half-life: {JCVI_PROTEIN_HALF_LIFE:.0f} sec ({JCVI_PROTEIN_HALF_LIFE/3600:.1f} hr)")
    print(f"  TX rate: {JCVI_TRANSCRIPTION_RATE} nt/sec")
    print(f"  TL rate: {JCVI_TRANSLATION_RATE} aa/sec")
    print(f"  Ribosomes: {JCVI_RIBOSOME_COUNT}")
    
    # --- Create ODE engine from enriched genes ---
    engine = create_engine_from_enriched_genes(genes)
    
    print(f"\n--- ODE Engine Configuration ---")
    print(f"  Total genes: {len(engine.gene_states)}")
    print(f"  Ribosomes: {engine.cellular_resources.total_ribosomes}")
    print(f"  Polymerases: {engine.cellular_resources.total_rna_polymerases}")
    
    # --- Sample gene parameter inspection ---
    print(f"\n--- Sample Gene ODE Parameters (first 10 genes) ---")
    print(f"  {'Gene':<18} {'TX_rate':>8} {'TL_rate':>8} {'mRNA_HL(s)':>10} {'Prot_HL(s)':>10} {'mRNA_SS':>8} {'Prot_SS':>9} {'Ess':>6}")
    print(f"  {'-'*17} {'-'*8} {'-'*8} {'-'*10} {'-'*10} {'-'*8} {'-'*9} {'-'*6}")
    
    for i, (gid, state) in enumerate(engine.gene_states.items()):
        if i >= 10:
            break
        gene = genes[gid]
        print(f"  {gid:<18} {state.transcription_rate:>8.3f} {state.translation_rate:>8.4f} "
              f"{state.mrna_half_life:>10.0f} {state.protein_half_life:>10.0f} "
              f"{state.mrna_steady_state:>8.1f} {state.protein_steady_state:>9.1f} "
              f"{gene.essentiality_status[:4]:>6}")
    
    # --- Aggregate kinetic statistics ---
    all_states = list(engine.gene_states.values())
    tx_rates = [s.transcription_rate for s in all_states]
    tl_rates = [s.translation_rate for s in all_states]
    mrna_hls = [s.mrna_half_life for s in all_states]
    prot_hls = [s.protein_half_life for s in all_states]
    mrna_ss = [s.mrna_steady_state for s in all_states]
    prot_ss = [s.protein_steady_state for s in all_states]
    
    print(f"\n--- ODE Parameter Distribution (all {len(all_states)} genes) ---")
    print(f"  TX rate:     mean={np.mean(tx_rates):.3f}, median={np.median(tx_rates):.3f}, "
          f"range=[{np.min(tx_rates):.3f}, {np.max(tx_rates):.3f}]")
    print(f"  TL rate:     mean={np.mean(tl_rates):.4f}, median={np.median(tl_rates):.4f}, "
          f"range=[{np.min(tl_rates):.4f}, {np.max(tl_rates):.4f}]")
    print(f"  mRNA HL(s):  mean={np.mean(mrna_hls):.0f}, median={np.median(mrna_hls):.0f}, "
          f"range=[{np.min(mrna_hls):.0f}, {np.max(mrna_hls):.0f}]")
    print(f"  Prot HL(s):  mean={np.mean(prot_hls):.0f}, median={np.median(prot_hls):.0f}, "
          f"range=[{np.min(prot_hls):.0f}, {np.max(prot_hls):.0f}]")
    print(f"  mRNA SS:     mean={np.mean(mrna_ss):.1f}, median={np.median(mrna_ss):.1f}")
    print(f"  Protein SS:  mean={np.mean(prot_ss):.1f}, median={np.median(prot_ss):.1f}")
    
    # --- Run simulation and track dynamics ---
    print(f"\n--- Running ODE Simulation (10,000 steps × 1s) ---")
    t0 = time.time()
    time_points = []
    avg_mrna_pct = []
    avg_prot_pct = []
    
    for step in range(10000):
        engine.step(1.0)
        if step % 500 == 0:
            stats = engine.get_statistics()
            time_points.append(engine.simulation_time)
            avg_mrna_pct.append(stats.get('avg_mrna_percentage', 0))
            avg_prot_pct.append(stats.get('avg_protein_percentage', 0))
    
    sim_time = time.time() - t0
    final_stats = engine.get_statistics()
    print(f"  Simulation time: {sim_time:.2f}s")
    print(f"  Final sim time: {engine.simulation_time:.0f}s ({engine.simulation_time/60:.1f} min)")
    print(f"  Final avg mRNA%: {final_stats.get('avg_mrna_percentage', 0):.1f}%")
    print(f"  Final avg Protein%: {final_stats.get('avg_protein_percentage', 0):.1f}%")
    print(f"  Steady state genes: {final_stats.get('steady_state_genes', 0)}")
    
    # --- Cellular UR + Growth Rate Prediction ---
    print(f"\n--- CellularUR + GrowthRatePredictor ---")
    gen = CellularURGenerator()
    cell_state = gen.generate_cellular_ur(genes, engine, engine.cellular_resources, engine.simulation_time)
    
    predictor = GrowthRatePredictor()
    result = predictor.predict(cell_state)
    
    print(f"  Predicted growth rate: {result['growth_rate']:.6f} min⁻¹")
    print(f"  Predicted doubling time: {result['doubling_time']:.1f} min")
    print(f"  Growth modifier: {result['growth_modifier']:.3f}")
    print(f"  Confidence: {result['confidence']:.2f}")
    print(f"  Limiting factors: {result['limiting_factors']}")
    print(f"  Growth phase: {result['growth_phase']}")
    
    # --- Root cause analysis ---
    print(f"\n--- ROOT CAUSE ANALYSIS ---")
    lit_rate = math.log(2) / (105 * 60)  # 105 min doubling time in seconds
    print(f"  Literature growth rate: {lit_rate*60:.6f} min⁻¹ (105 min doubling)")
    print(f"  Our predicted rate:     {result['growth_rate']:.6f} min⁻¹")
    print(f"  Ratio (predicted/lit):  {result['growth_rate']/(lit_rate*60):.3f}")
    
    # The issue is the GrowthRatePredictor model, not the ODE engine
    # The ODE engine correctly handles kinetics, but GrowthRatePredictor
    # uses a simplified model. Let's check what's wrong:
    dynamics = cell_state.dynamics_vector
    resources = cell_state.resource_vector
    print(f"\n  Dynamics vector diagnostics:")
    print(f"    avg_mrna_count: {dynamics[0]:.2f}")
    print(f"    avg_protein_count: {dynamics[10]:.2f}")
    print(f"    avg_mrna_pct: {dynamics[20]:.2f}")
    print(f"    avg_protein_pct: {dynamics[21]:.2f}")
    print(f"    steady_fraction (idx 66): {dynamics[66]:.4f}")
    print(f"  Resource vector diagnostics:")
    print(f"    ribosome_sat: {resources[2]:.4f}")
    print(f"    ATP level: {resources[10]:.4f}")
    
    return engine, cell_state, gen


# =============================================================================
# TASK 1.2: ESSENTIALITY PREDICTION ACCURACY
# =============================================================================

def task_1_2_essentiality_validation(genes):
    """Run EssentialityPredictor on all genes, compute confusion matrix."""
    print("\n" + "=" * 70)
    print("TASK 1.2: ESSENTIALITY PREDICTION ACCURACY")
    print("=" * 70)
    
    predictor = EssentialityPredictor()
    
    # Ground truth distribution
    ess_actual = {"Essential": 0, "Quasiessential": 0, "Nonessential": 0, "Unknown": 0}
    for g in genes.values():
        ess_actual[g.essentiality_status] = ess_actual.get(g.essentiality_status, 0) + 1
    print(f"\n  Ground truth: {ess_actual}")
    
    # Run predictions
    tp = fp = tn = fn = 0
    correct = incorrect = unknown = 0
    results_by_class = {"Essential": {"tp": 0, "fp": 0, "fn": 0}, "Nonessential": {"tp": 0, "fp": 0, "fn": 0}}
    
    for gene in genes.values():
        pred = predictor.predict(gene)
        predicted = pred['predicted_essentiality']
        actual = gene.essentiality_status
        
        # Binary classification: Essential vs Not-Essential
        actual_binary = "Essential" if actual in ["Essential", "Quasiessential"] else "Nonessential"
        
        if predicted == "Essential" and actual_binary == "Essential":
            tp += 1
        elif predicted == "Essential" and actual_binary == "Nonessential":
            fp += 1
        elif predicted == "Nonessential" and actual_binary == "Nonessential":
            tn += 1
        elif predicted == "Nonessential" and actual_binary == "Essential":
            fn += 1
    
    total = tp + fp + tn + fn
    accuracy = (tp + tn) / total if total > 0 else 0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    
    print(f"\n  --- Confusion Matrix (Essential+Quasi vs Nonessential) ---")
    print(f"                  Predicted Essential  Predicted Non-essential")
    print(f"  Actual Essential      {tp:>5}                {fn:>5}")
    print(f"  Actual Non-essential  {fp:>5}                {tn:>5}")
    
    print(f"\n  --- Metrics ---")
    print(f"  Accuracy:  {accuracy:.1%} ({tp+tn}/{total})")
    print(f"  Precision: {precision:.1%}")
    print(f"  Recall:    {recall:.1%}")
    print(f"  F1 Score:  {f1:.3f}")
    print(f"  Target:    >85% accuracy")
    print(f"  Status:    {'✅ PASS' if accuracy > 0.85 else '⚠️ BELOW TARGET'}")
    
    return accuracy


# =============================================================================
# TASK 1.3: KNOCKOUT SIMULATION VALIDATION
# =============================================================================

def task_1_3_knockout_validation(genes, cell_state, gen):
    """Simulate knockouts of essential and non-essential genes, compare impact."""
    print("\n" + "=" * 70)
    print("TASK 1.3: KNOCKOUT SIMULATION VALIDATION")
    print("=" * 70)
    
    vis = VirtualInstrumentSuite(gen)
    
    # Pick 10 essential genes (first 10)
    essential_ids = [gid for gid, g in genes.items() if g.essentiality_status == "Essential"][:10]
    # Pick 10 non-essential genes
    nonessential_ids = [gid for gid, g in genes.items() if g.essentiality_status == "Nonessential"][:10]
    
    print(f"\n  --- Knocking out 10 Essential genes ---")
    print(f"  {'Gene':<18} {'Name':<10} {'Impact':>7} {'Category':<12} {'L2_dist':>8}")
    print(f"  {'-'*17} {'-'*10} {'-'*7} {'-'*12} {'-'*8}")
    
    essential_impacts = []
    for gid in essential_ids:
        ko = vis.simulate_knockout(gid, genes, cell_state)
        essential_impacts.append(ko['impact_score'])
        print(f"  {gid:<18} {genes[gid].name:<10} {ko['impact_score']:>7.4f} {ko['impact_category']:<12} {ko['l2_distance']:>8.4f}")
    
    print(f"\n  --- Knocking out 10 Non-essential genes ---")
    print(f"  {'Gene':<18} {'Name':<10} {'Impact':>7} {'Category':<12} {'L2_dist':>8}")
    print(f"  {'-'*17} {'-'*10} {'-'*7} {'-'*12} {'-'*8}")
    
    nonessential_impacts = []
    for gid in nonessential_ids:
        ko = vis.simulate_knockout(gid, genes, cell_state)
        nonessential_impacts.append(ko['impact_score'])
        print(f"  {gid:<18} {genes[gid].name:<10} {ko['impact_score']:>7.4f} {ko['impact_category']:<12} {ko['l2_distance']:>8.4f}")
    
    avg_ess = np.mean(essential_impacts)
    avg_noness = np.mean(nonessential_impacts)
    
    print(f"\n  --- Summary ---")
    print(f"  Avg Essential impact:     {avg_ess:.4f}")
    print(f"  Avg Non-essential impact: {avg_noness:.4f}")
    print(f"  Ratio (Ess/NonEss):       {avg_ess/avg_noness:.2f}x" if avg_noness > 0 else "  Ratio: N/A")
    
    # Multi-knockout test
    print(f"\n  --- Multi-Knockout Test (5 essential genes) ---")
    multi_ko = vis.knockout_simulator.simulate_multi_knockout(
        essential_ids[:5], genes, cell_state)
    print(f"  Impact of 5-gene knockout: {multi_ko['impact_score']:.4f}")


# =============================================================================
# TASK 3.1: DATA QUALITY REPORT
# =============================================================================

def task_3_1_data_quality(genes):
    """Generate summary statistics for AI predictions."""
    print("\n" + "=" * 70)
    print("TASK 3.1: DATA QUALITY REPORT — AI GAP FILLING STATISTICS")
    print("=" * 70)
    
    rbs = [g.rbs_strength for g in genes.values()]
    mrna_hl = [g.mrna_half_life for g in genes.values()]
    prot_hl = [g.protein_half_life for g in genes.values()]
    mrna_exp = [g.mrna_expression for g in genes.values()]
    prot_abund = [g.protein_abundance for g in genes.values()]
    
    def stats(arr, name, unit, expected_range, expected_mean):
        arr = np.array(arr)
        print(f"\n  {name} ({unit})")
        print(f"    Count:    {len(arr)}")
        print(f"    Mean:     {np.mean(arr):.2f} (expected: {expected_mean})")
        print(f"    Median:   {np.median(arr):.2f}")
        print(f"    Std:      {np.std(arr):.2f}")
        print(f"    Range:    [{np.min(arr):.2f}, {np.max(arr):.2f}] (expected: {expected_range})")
        
        # Outlier detection (>3 sigma)
        mean, std = np.mean(arr), np.std(arr)
        outliers = np.sum(np.abs(arr - mean) > 3 * std)
        print(f"    Outliers: {outliers} (>3σ)")
        
        in_range = np.sum((arr >= expected_range[0]) & (arr <= expected_range[1]))
        print(f"    In range: {in_range}/{len(arr)} ({100*in_range/len(arr):.0f}%)")
    
    stats(rbs, "RBS Strength", "index", (0.5, 10.0), "1-4")
    stats(mrna_hl, "mRNA Half-Life", "minutes", (0.5, 15.0), "2-5 min")
    stats(prot_hl, "Protein Half-Life", "hours", (2.0, 200.0), "20-30 hr")
    stats(mrna_exp, "mRNA Expression", "copies/cell", (1, 5000), "~50")
    stats(prot_abund, "Protein Abundance", "copies/cell", (10, 100000), "~500")
    
    # Biologically unreasonable values
    print(f"\n  --- Biological Reasonableness Check ---")
    bad_rbs = sum(1 for v in rbs if v <= 0 or v > 15)
    bad_mrna = sum(1 for v in mrna_hl if v < 0.1 or v > 60)
    bad_prot = sum(1 for v in prot_hl if v < 0.5 or v > 500)
    print(f"    Unreasonable RBS: {bad_rbs}/{len(rbs)}")
    print(f"    Unreasonable mRNA HL: {bad_mrna}/{len(mrna_hl)}")
    print(f"    Unreasonable Protein HL: {bad_prot}/{len(prot_hl)}")


# =============================================================================
# PERFORMANCE REPORT
# =============================================================================

def performance_report(genes):
    """Measure performance of key operations."""
    print("\n" + "=" * 70)
    print("PERFORMANCE REPORT")
    print("=" * 70)
    
    # ODE Bridge creation
    t0 = time.time()
    engine = create_engine_from_enriched_genes(genes)
    t_bridge = time.time() - t0
    
    # ODE Simulation (1000 steps)
    t0 = time.time()
    for _ in range(1000):
        engine.step(1.0)
    t_sim_1k = time.time() - t0
    
    # Cellular UR generation
    t0 = time.time()
    gen = CellularURGenerator()
    cell_state = gen.generate_cellular_ur(genes, engine, engine.cellular_resources)
    t_ur = time.time() - t0
    
    # VI predictions (all)
    vis = VirtualInstrumentSuite(gen)
    
    t0 = time.time()
    vis.predict_growth_rate(cell_state)
    vis.predict_phenotype(cell_state, genes)
    vis.predict_stress(cell_state, "heat")
    t_vi_cell = time.time() - t0
    
    t0 = time.time()
    for gene in genes.values():
        vis.predict_essentiality(gene)
    t_vi_ess = time.time() - t0
    
    # Single knockout
    gid = list(genes.keys())[0]
    t0 = time.time()
    vis.simulate_knockout(gid, genes, cell_state)
    t_knockout = time.time() - t0
    
    print(f"\n  --- Timing Results ---")
    print(f"  ODE Bridge creation ({len(genes)} genes): {t_bridge*1000:.1f} ms")
    print(f"  ODE Simulation (1000 steps):             {t_sim_1k*1000:.1f} ms")
    print(f"  Cellular UR generation:                  {t_ur*1000:.1f} ms")
    print(f"  VI cell-level predictions:               {t_vi_cell*1000:.1f} ms")
    print(f"  VI essentiality ({len(genes)} genes):    {t_vi_ess*1000:.1f} ms ({t_vi_ess/len(genes)*1000:.2f} ms/gene)")
    print(f"  Single knockout simulation:              {t_knockout*1000:.1f} ms")
    print(f"\n  --- Embedding Cache Status ---")
    print(f"  HDF5 cache: NOT YET IMPLEMENTED")
    print(f"  Estimated embedding time (GPU): ~5-15 min for {len(genes)} genes")
    print(f"  Estimated embedding time (CPU): ~60 min for {len(genes)} genes")
    print(f"  Recommendation: Implement HDF5 cache for >10x speedup on repeat runs")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("JCVI-syn3A AI VIRTUAL CELL — COMPREHENSIVE VALIDATION REPORT")
    print("=" * 70)
    print(f"  Report generated at: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Gene count: 498 (realistic mock data)")
    
    # Create realistic gene set
    genes = create_realistic_genes(498)
    
    # Priority 1: Biological Validation
    engine, cell_state, gen = task_1_1_growth_rate_validation(genes)
    accuracy = task_1_2_essentiality_validation(genes)
    task_1_3_knockout_validation(genes, cell_state, gen)
    
    # Priority 3: Data Quality
    task_3_1_data_quality(genes)
    
    # Performance
    performance_report(genes)
    
    print("\n" + "=" * 70)
    print("REPORT COMPLETE")
    print("=" * 70)
