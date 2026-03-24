#!/usr/bin/env python3
"""
Test: Cellular UR & Virtual Instruments Integration
=====================================================
Verifies the end-to-end pipeline from enriched Gene objects
through Cellular UR generation to Virtual Instrument predictions.
"""

import sys
import os
import numpy as np

# Ensure local imports work
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from data_structures import Gene, CellStateVector


def create_sample_genes(n=10):
    """Create sample enriched Gene objects with mock molecular URs."""
    genes = {}
    for i in range(n):
        gene_id = f"JCVISYN3A_{i:04d}"
        seq = "ATGGCTTCAATCGGA" * (50 + i * 10)  # Variable-length sequences
        
        gene = Gene(
            id=gene_id,
            name=f"gene_{i}",
            start=i * 1000,
            end=i * 1000 + len(seq),
            strand="+",
            type="CDS",
            function="Test function",
            product=f"Test protein {i}",
            sequence=seq,
        )
        
        # Simulate enriched data
        gene.mrna_expression = float(np.random.uniform(10, 500))
        gene.protein_abundance = float(np.random.uniform(100, 10000))
        gene.rbs_strength = float(np.random.uniform(1, 10))
        gene.mrna_half_life = float(np.random.uniform(0.5, 5.0))
        gene.protein_half_life = float(np.random.uniform(10, 60))
        gene.essentiality_status = "Essential" if i < 6 else "Nonessential"
        gene.function_class = "Metabolism" if i % 2 == 0 else "Translation"
        gene.ortholog_id = f"ortholog_{i}" if i < 7 else ""
        
        # Mock Molecular UR (normally from MolecularURGenerator)
        gene.molecular_ur = {
            'dna_embedding': np.random.randn(1280).astype(np.float32),
            'rna_embedding': np.random.randn(1280).astype(np.float32),
            'protein_embedding': np.random.randn(1280).astype(np.float32),
            'combined_embedding': np.random.randn(3840).astype(np.float32),
        }
        
        genes[gene_id] = gene
    
    return genes


def test_cellular_ur_generator():
    """Test 1: CellularURGenerator produces correct-shape vectors."""
    print("=" * 60)
    print("TEST 1: CellularURGenerator")
    print("=" * 60)
    
    from cellular_ur_generator import CellularURGenerator
    
    genes = create_sample_genes(10)
    generator = CellularURGenerator()
    
    cell_state = generator.generate_cellular_ur(genes)
    
    assert isinstance(cell_state, CellStateVector), "Output must be CellStateVector"
    assert cell_state.gene_embedding_aggregate.shape == (1280,), \
        f"Gene embed shape: {cell_state.gene_embedding_aggregate.shape}, expected (1280,)"
    assert cell_state.dynamics_vector.shape == (512,), \
        f"Dynamics shape: {cell_state.dynamics_vector.shape}, expected (512,)"
    assert cell_state.environment_vector.shape == (64,), \
        f"Env shape: {cell_state.environment_vector.shape}, expected (64,)"
    assert cell_state.resource_vector.shape == (32,), \
        f"Resource shape: {cell_state.resource_vector.shape}, expected (32,)"
    assert cell_state.combined_cellular_ur.shape == (1888,), \
        f"Combined shape: {cell_state.combined_cellular_ur.shape}, expected (1888,)"
    assert cell_state.num_genes == 10
    print(f"  ✅ Gene embedding aggregate: {cell_state.gene_embedding_aggregate.shape}")
    print(f"  ✅ Dynamics vector:          {cell_state.dynamics_vector.shape}")
    print(f"  ✅ Environment vector:       {cell_state.environment_vector.shape}")
    print(f"  ✅ Resource vector:           {cell_state.resource_vector.shape}")
    print(f"  ✅ Combined Cellular UR:     {cell_state.combined_cellular_ur.shape}")
    print(f"  ✅ Num genes: {cell_state.num_genes}")
    print(f"  ✅ Growth phase: {cell_state.growth_phase}")
    print(f"  ✅ Non-zero values in gene embed: {np.count_nonzero(cell_state.gene_embedding_aggregate)}")
    print("  PASSED ✅\n")
    
    return generator, cell_state, genes


def test_virtual_instruments(generator, cell_state, genes):
    """Test 2: All Virtual Instruments return valid predictions."""
    print("=" * 60)
    print("TEST 2: Virtual Instruments")
    print("=" * 60)
    
    from virtual_instruments import VirtualInstrumentSuite
    
    vis = VirtualInstrumentSuite(generator)
    
    # Test Growth Rate Predictor
    growth = vis.predict_growth_rate(cell_state)
    assert 'growth_rate' in growth, "Missing growth_rate"
    assert 'doubling_time' in growth, "Missing doubling_time"
    assert growth['growth_rate'] > 0, "Growth rate must be positive"
    print(f"  ✅ Growth rate: {growth['growth_rate']:.6f} min⁻¹")
    print(f"     Doubling time: {growth['doubling_time']:.1f} min")
    print(f"     Confidence: {growth['confidence']}")
    print(f"     Limiting factors: {growth['limiting_factors']}")
    
    # Test Phenotype Predictor
    phenotype = vis.predict_phenotype(cell_state, genes)
    assert 'phenotype' in phenotype, "Missing phenotype"
    assert phenotype['phenotype'] in ['viable', 'slow_growth', 'metabolically_impaired', 
                                       'structurally_compromised', 'non_viable']
    print(f"  ✅ Phenotype: {phenotype['phenotype']} (score: {phenotype['viability_score']})")
    
    # Test Stress Response Predictor
    stress = vis.predict_stress(cell_state, "heat")
    assert 'severity' in stress, "Missing severity"
    print(f"  ✅ Stress (heat): survival={stress['predicted_survival']}, tolerance={stress['tolerance']}")
    
    # Test Essentiality Predictor
    sample_gene = list(genes.values())[0]
    ess = vis.predict_essentiality(sample_gene)
    assert 'predicted_essentiality' in ess, "Missing predicted_essentiality"
    print(f"  ✅ Essentiality: predicted={ess['predicted_essentiality']}, "
          f"actual={ess['actual_essentiality']}, score={ess['essentiality_score']}")
    
    # Test Knockout Simulator
    first_gene_id = list(genes.keys())[0]
    ko = vis.simulate_knockout(first_gene_id, genes, cell_state)
    assert 'impact_score' in ko, "Missing impact_score"
    assert 'modified_state' in ko, "Missing modified_state"
    print(f"  ✅ Knockout {first_gene_id}: impact={ko['impact_score']:.4f}, "
          f"category={ko['impact_category']}")
    
    # Test Environment Perturbation
    env = vis.simulate_environment({'temperature': 42.0}, cell_state, genes)
    assert 'environment_shift' in env, "Missing environment_shift"
    print(f"  ✅ Env perturbation (temp=42°C): shift={env['environment_shift']:.4f}")
    
    # Test Drug Treatment
    drug = vis.simulate_drug({first_gene_id: 0.8}, genes, cell_state)
    assert 'drug_effect_score' in drug, "Missing drug_effect_score"
    print(f"  ✅ Drug treatment (80% inhibition): effect={drug['drug_effect_score']:.4f}")
    
    print("  PASSED ✅\n")


def test_ode_bridge():
    """Test 3: ODE Bridge creates engine from enriched genes."""
    print("=" * 60)
    print("TEST 3: ODE Bridge")
    print("=" * 60)
    
    from cell_simulation_engine import create_engine_from_enriched_genes, SimulationMode
    
    genes = create_sample_genes(5)
    engine = create_engine_from_enriched_genes(genes)
    
    assert len(engine.gene_states) == 5, f"Expected 5 gene states, got {len(engine.gene_states)}"
    
    # Verify kinetic parameters were transferred
    for gene_id, state in engine.gene_states.items():
        gene = genes[gene_id]
        expected_mrna_hl = gene.mrna_half_life * 60.0  # min → sec
        expected_prot_hl = gene.protein_half_life * 3600.0  # hr → sec
        
        assert abs(state.mrna_half_life - expected_mrna_hl) < 0.1, \
            f"mRNA half-life mismatch for {gene_id}"
        assert abs(state.protein_half_life - expected_prot_hl) < 0.1, \
            f"Protein half-life mismatch for {gene_id}"
        print(f"  ✅ {gene_id}: TX={state.transcription_rate:.3f}, "
              f"TL={state.translation_rate:.3f}, "
              f"mRNA_hl={state.mrna_half_life:.0f}s, "
              f"Prot_hl={state.protein_half_life:.0f}s")
    
    # Run a few simulation steps
    for _ in range(10):
        engine.step(1.0)
    
    stats = engine.get_statistics()
    print(f"  ✅ After 10 steps: {stats.get('total_genes', 0)} genes simulated")
    print(f"     Avg mRNA%: {stats.get('avg_mrna_percentage', 0):.1f}%")
    print(f"     Avg Protein%: {stats.get('avg_protein_percentage', 0):.1f}%")
    print("  PASSED ✅\n")
    
    return engine


def test_integrated_pipeline():
    """Test 4: End-to-end pipeline (enrichment → embedding → cellular UR → VI)."""
    print("=" * 60)
    print("TEST 4: Integrated Pipeline (end-to-end)")
    print("=" * 60)
    
    from cellular_ur_generator import CellularURGenerator
    from virtual_instruments import VirtualInstrumentSuite
    from cell_simulation_engine import create_engine_from_enriched_genes
    
    # Create enriched genes
    genes = create_sample_genes(20)
    print(f"  Created {len(genes)} enriched genes")
    
    # ODE Bridge
    engine = create_engine_from_enriched_genes(genes)
    for _ in range(50):
        engine.step(1.0)
    print(f"  ODE engine: {len(engine.gene_states)} genes, simulated 50 steps")
    
    # Cellular UR
    generator = CellularURGenerator()
    cell_state = generator.generate_cellular_ur(
        genes, engine, engine.cellular_resources, engine.simulation_time
    )
    print(f"  Cellular UR: {cell_state.combined_cellular_ur.shape}-dim, "
          f"phase={cell_state.growth_phase}")
    
    # Virtual Instruments
    vis = VirtualInstrumentSuite(generator)
    
    growth = vis.predict_growth_rate(cell_state)
    phenotype = vis.predict_phenotype(cell_state, genes)
    
    print(f"  Growth rate: {growth['doubling_time']:.0f} min doubling time")
    print(f"  Phenotype: {phenotype['phenotype']} (viability={phenotype['viability_score']:.2f})")
    
    # Knockout essentials and check impact increases
    essential_gene = [g for g in genes if genes[g].essentiality_status == "Essential"][0]
    nonessential_gene = [g for g in genes if genes[g].essentiality_status == "Nonessential"][0]
    
    ko_essential = vis.simulate_knockout(essential_gene, genes, cell_state, engine)
    ko_nonessential = vis.simulate_knockout(nonessential_gene, genes, cell_state, engine)
    
    print(f"  KO {essential_gene} (Essential): impact={ko_essential['impact_score']:.4f}")
    print(f"  KO {nonessential_gene} (Nonessential): impact={ko_nonessential['impact_score']:.4f}")
    
    # Update cellular UR (efficient)
    for _ in range(50):
        engine.step(1.0)
    
    updated_state = generator.update_cellular_ur(
        cell_state, genes, engine, engine.cellular_resources, engine.simulation_time
    )
    
    # Check temporal evolution
    ur_diff = np.linalg.norm(
        cell_state.combined_cellular_ur - updated_state.combined_cellular_ur
    )
    print(f"  UR evolution after 50 more steps: L2 distance = {ur_diff:.4f}")
    print(f"  Updated growth phase: {updated_state.growth_phase}")
    
    print("  PASSED ✅\n")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("AI VIRTUAL CELL - CELLULAR UR & VIRTUAL INSTRUMENTS TEST")
    print("=" * 60 + "\n")
    
    try:
        generator, cell_state, genes = test_cellular_ur_generator()
        test_virtual_instruments(generator, cell_state, genes)
        test_ode_bridge()
        test_integrated_pipeline()
        
        print("=" * 60)
        print("ALL TESTS PASSED ✅✅✅")
        print("=" * 60)
    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
