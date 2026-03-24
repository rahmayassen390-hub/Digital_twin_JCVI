import os
import sys
import numpy as np
import h5py
import copy
import random
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

# Local imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from smart_data_engine import SmartDataEngine
from data_enrichment_layer import DataEnrichmentLayer
from molecular_ur_generator import MolecularURGenerator
from cellular_ur_generator import CellularURGenerator
from cell_simulation_engine import create_engine_from_enriched_genes
from virtual_instruments import VirtualInstrumentSuite

class PerturbationEngine:
    def __init__(self, data_dir="../Data", output_file="perturbation_dataset_50k.h5"):
        self.data_dir = data_dir
        self.output_file = output_file
        
        # 1. Initialize core engines
        print("Initializing Base Engines...")
        self.smart_engine = SmartDataEngine()
        self.enricher = DataEnrichmentLayer(data_dir=data_dir)
        self.mol_gen = MolecularURGenerator(load_models=False)
        self.cell_gen = CellularURGenerator()
        
        # 2. Load and Enrich Base Genome (Syn3A Set: 493 genes)
        print("Standardizing to 493 protein-coding genes (Syn3A Set)...")
        import pandas as pd
        df_genes = pd.read_csv(os.path.join(data_dir, "JCVI_syn3A_full_genome_dataset.csv"))
        target_locus_tags = df_genes['locus_tag'].head(493).tolist()
        
        self.base_genes = {}
        from data_structures import Gene
        for locus_id in target_locus_tags:
            row = df_genes[df_genes['locus_tag'] == locus_id].iloc[0]
            g = Gene(
                id=locus_id,
                name=str(row.get('gene_name', locus_id)),
                start=0, end=0, strand='+', type='CDS',
                product=str(row.get('product_description', 'Unknown'))
            )
            # Enrich using the local enricher
            enriched_g = self.enricher.enrich_gene(g)
            
            # Mock Molecular UR (1280 dim)
            enriched_g.molecular_ur = {
                'dna_embedding': np.random.randn(1280).astype(np.float32),
                'combined_embedding': np.random.randn(1280).astype(np.float32)
            }
            self.base_genes[locus_id] = enriched_g
        
        self.vi_suite = VirtualInstrumentSuite(self.cell_gen)
        print(f"Base Genome Initialized: {len(self.base_genes)} genes (Syn3A Standard).")

    def generate_single_sample(self, seed_and_worker):
        """Generates a single perturbation sample."""
        seed, worker_idx = seed_and_worker
        np.random.seed(seed)
        random.seed(seed)
        
        # Assign GPU for this worker
        from gpu_orchestrator import orchestrator
        device = orchestrator.get_device_for_worker(worker_idx)
        
        # Local re-initialization of model-heavy components on specific GPU
        # We use a cache-like approach to avoid re-init for every sample
        if not hasattr(self, '_worker_vIs'):
            self._worker_mol_gen = MolecularURGenerator(load_models=False, device=device)
            self._worker_vIs = VirtualInstrumentSuite(self.cell_gen, device=device)
        
        # Deepcopy to avoid state leakage
        sample_genes = copy.deepcopy(self.base_genes)
        
        # 1. Random Perturbation Logic
        perturb_type = random.choices(
            ['knockout', 'expression_shift', 'environment_stress', 'resource_stress'],
            weights=[0.3, 0.4, 0.2, 0.1]
        )[0]
        
        perturb_info = {"type": perturb_type}
        
        if perturb_type == 'knockout':
            # 1-3 random genes knocked out
            num_ko = random.randint(1, 4)
            ko_targets = random.sample(list(sample_genes.keys()), num_ko)
            for gid in ko_targets:
                sample_genes[gid].mrna_expression = 0.0
                sample_genes[gid].protein_abundance = 0.0
                # Zero out embedding to signal absence
                sample_genes[gid].molecular_ur = {k: np.zeros_like(v) for k,v in sample_genes[gid].molecular_ur.items()}
            perturb_info["targets"] = ko_targets
            
        elif perturb_type == 'expression_shift':
            # Gaussian noise on 20% of genes
            target_ids = random.sample(list(sample_genes.keys()), int(len(sample_genes)*0.2))
            for gid in target_ids:
                shift = np.random.normal(1.0, 0.2)
                sample_genes[gid].mrna_expression *= max(0.0, shift)
            perturb_info["shift_scale"] = 0.2
            
        # 2. Simulate
        try:
            sim = create_engine_from_enriched_genes(sample_genes)
            # Perturb environment if needed
            if perturb_type == 'environment_stress':
                temp_shift = np.random.uniform(0.8, 1.2)
                sim.global_params['temperature'] *= temp_shift
                perturb_info["temp_shift"] = temp_shift
            
            # Step to stabilize
            sim.step(10.0)
            
            # 3. Generate UR and Phenotype
            from cellular_ur_generator import CellularURGenerator
            cell_ur_gen = CellularURGenerator() # Lightweight
            cell_state = cell_ur_gen.generate_cellular_ur(sample_genes, sim)
            growth_res = self._worker_vIs.predict_growth_rate(cell_state)
            
            return {
                "ur": cell_state.combined_cellular_ur.astype(np.float32),
                "growth_rate": float(growth_res['growth_rate']),
                "doubling_time": float(growth_res['doubling_time']),
                "viability": 1.0 if growth_res['growth_rate'] > 1e-6 else 0.0
            }
        except Exception as e:
            return None

    def run_batch(self, num_samples=100, output_file="Data/perturbation_results.h5"):
        """Run a batch of perturbations in parallel."""
        seeds = [(random.randint(0, 1000000), i % 4) for i in range(num_samples)]
        
        results = []
        print(f"🚀 Starting Quad-GPU Batch: {num_samples} samples...")
        
        with ProcessPoolExecutor(max_workers=4) as executor:
            # We use a custom mapping to track worker IDs for GPU assignment
            futures = [executor.submit(self.generate_single_sample, s) for s in seeds]
            for f in tqdm(as_completed(futures), total=num_samples, desc="Perturbing"):
                try:
                    res = f.result()
                    if res:
                        results.append(res)
                except Exception as e:
                    print(f"Sample failed: {e}")
        
        # Save results to HDF5
        self._save_results(results, output_file)
        print(f"✅ Batch complete. Saved {len(results)} samples to {output_file}")
        return results

    def _save_results(self, results, output_file):
        """Helper to save collected results to an HDF5 file."""
        if not results:
            print("No results to save.")
            return

        num_samples = len(results)
        ur_dim = results[0]["ur"].shape[0]

        with h5py.File(output_file, 'w') as f:
            f.create_dataset("cell_ur", (num_samples, ur_dim), dtype='f4')
            f.create_dataset("growth_rate", (num_samples,), dtype='f4')
            f.create_dataset("doubling_time", (num_samples,), dtype='f4')
            f.create_dataset("viability", (num_samples,), dtype='f4')
            
            for idx, res in enumerate(results):
                f["cell_ur"][idx] = res["ur"]
                f["growth_rate"][idx] = res["growth_rate"]
                f["doubling_time"][idx] = res["doubling_time"]
                f["viability"][idx] = res["viability"]

if __name__ == "__main__":
    # Paper Alignment: Generate verification batch
    engine = PerturbationEngine()
    engine.run_batch(num_samples=100)
