import time
import os
import psutil
import numpy as np
from smart_data_engine import SmartDataEngine
from data_enrichment_layer import DataEnrichmentLayer
from cell_simulation_engine import CellSimulationEngine, create_engine_from_enriched_genes
from cellular_ur_generator import CellularURGenerator
from molecular_ur_generator import MolecularURGenerator

def benchmark_full_pipeline():
    process = psutil.Process(os.getpid())
    results = {}
    
    # Data directory is one level up
    data_dir = "../Data"
    
    # Specific JCVI 3.0 file mapping
    file_map = {
        'gff': os.path.join(data_dir, "JCVI_3.0.gff3"),
        'regulation': os.path.join(data_dir, "JCVI_Gene_Regulation_3.0#.bed"),
        'transcription': os.path.join(data_dir, "JCVI_Transcription_Dynamics_3.0#.bed"),
        'translation': os.path.join(data_dir, "JCVI_Translation_Dynamics_3.0#.bed"),
        'promoters': os.path.join(data_dir, "Likely_Promoters_3.0#.bed"),
        'operons': os.path.join(data_dir, "Internal_Operons_3.0#.bed")
    }
    
    print(f"🚀 Starting Full 498-Gene Pipeline Benchmark...")
    print(f"Using 3.0 Data Files from {data_dir}")
    
    # === 1. Smart Integration ===
    start = time.time()
    engine = SmartDataEngine()
    
    # Explicitly load all files
    engine.load_gff3(file_map['gff'])
    engine.load_regulation_bed(file_map['regulation'])
    engine.load_transcription_bed(file_map['transcription'])
    engine.load_translation_bed(file_map['translation'])
    engine.load_promoters_bed(file_map['promoters'])
    engine.load_operons_bed(file_map['operons'])
    
    integrated_genes = engine.integrate_all_data()
    results['integration_time'] = time.time() - start
    print(f"✅ Integration: {len(integrated_genes)} genes in {results['integration_time']:.3f}s")
    
    if len(integrated_genes) == 0:
        print("❌ Error: No genes integrated. Check data paths and filenames.")
        return
    
    # === 2. Enrichment ===
    start = time.time()
    enricher = DataEnrichmentLayer(data_dir=data_dir)
    enriched_genes = {gid: enricher.enrich_gene(g) for gid, g in integrated_genes.items()}
    results['enrichment_time'] = time.time() - start
    print(f"✅ Enrichment: {results['enrichment_time']:.3f}s")
    
    # === 3. ODE Bridge ===
    start = time.time()
    sim_engine = create_engine_from_enriched_genes(enriched_genes)
    results['ode_init_time'] = time.time() - start
    # No simulate() method, use step() loop
    print(f"✅ ODE Bridge Init: {results['ode_init_time']:.3f}s")
    
    # === 4. ODE Simulation (1000s) ===
    start = time.time()
    dt = 10.0
    steps = 100
    for _ in range(steps):
        sim_engine.step(dt)
    results['ode_sim_time'] = time.time() - start
    print(f"✅ ODE Simulation (1000s total, step-wise): {results['ode_sim_time']:.3f}s")
    
    # === 5. Molecular UR (Cached Lookups) ===
    start = time.time()
    # Test mode: load_models=False, but use cache Lookups
    ur_gen = MolecularURGenerator(load_models=False, cache_dir=data_dir)
    for gene in enriched_genes.values():
        ur_gen.generate_ur(gene)
    results['molecular_ur_time'] = time.time() - start
    print(f"✅ Molecular UR (Cached Lookups): {results['molecular_ur_time']:.3f}s")
    
    # === 6. Cellular UR ===
    start = time.time()
    cell_ur_gen = CellularURGenerator()
    cell_state = cell_ur_gen.generate_cellular_ur(enriched_genes, sim_engine)
    results['cellular_ur_time'] = time.time() - start
    print(f"✅ Cellular UR: {results['cellular_ur_time']:.3f}s")
    
    # === Memory Footprint ===
    mem_mb = process.memory_info().rss / 1024 / 1024
    results['memory_mb'] = mem_mb
    print(f"📊 Final Memory Footprint: {mem_mb:.2f} MB")
    
    return results

if __name__ == "__main__":
    benchmark_full_pipeline()
