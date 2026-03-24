import sys
import os
import time

# Add current directory to path
sys.path.append(os.getcwd())

from cell_simulation_engine import create_engine_from_enriched_genes, SimulationMode
from metabolic_solver import MetabolicSolver
from data_structures import Gene

def test_integration():
    print("🧪 Testing Integrated CellSimulationEngine + Excel MetabolicSolver")
    
    # Load solver first to get gene list
    temp_solver = MetabolicSolver("../Data/metabolic_reconstruction.xlsx")
    
    genes = {}
    for gene_id in temp_solver.gene_rxn_map.keys():
        genes[gene_id] = Gene(
            id=gene_id, 
            name=f"gene_{gene_id}", 
            essentiality_status="Essential", 
            mrna_expression=1.0, 
            rbs_strength=5.0, 
            start=1, 
            end=1000, 
            strand="+", 
            type="CDS"
        )
    
    print(f"Mocking {len(genes)} genes for whole-cell simulation")
    
    # 2. Create engine via factory
    # This should automatically load ../Data/metabolic_reconstruction.xlsx
    engine = create_engine_from_enriched_genes(genes, mode=SimulationMode.CONTINUOUS)
    
    if engine.metabolic_solver:
        print(f"✅ MetabolicSolver integrated: {engine.metabolic_solver.model_path}")
        print(f"✅ Reversible reactions: {len([r for r in engine.metabolic_solver.reactions if 'EX_' in r])} exchange reactions detected")
    else:
        print("❌ MetabolicSolver NOT integrated")
        return

    # 3. Simulate and check flux
    print("\nRunning simulation steps...")
    # Force immediate FBA update
    engine.last_fba_update = -1000.0 
    
    engine.step(1.0)
    
    fluxes = engine.current_fluxes
    biomass_flux = fluxes.get('EX_biomass_c', 0.0)
    
    print(f"Biomass Flux (EX_biomass_c): {biomass_flux:.6f}")
    
    if biomass_flux > 0:
        print(f"🎉 SUCCESS: Integrated engine achieved GROWTH flux: {biomass_flux:.6f}")
        # Calculate doubling time
        # mu = biomass_flux (hr-1)
        doubling_time_min = (0.693 / biomass_flux) * 60.0 if biomass_flux > 0 else float('inf')
        print(f"Predicted Doubling Time: {doubling_time_min:.1f} min")
    else:
        print("❌ FAILURE: Zero growth flux in integrated engine")
        print("Checking for possible reasons...")
        # Check objective coeffs
        obj_idx = [i for i, c in enumerate(engine.metabolic_solver.objective_coeffs) if c > 0]
        if obj_idx:
            obj_rxn = engine.metabolic_solver.reactions[obj_idx[0]]
            print(f"Target Objective: {obj_rxn}")
        else:
            print("❌ No objective defined in solver")

if __name__ == "__main__":
    test_integration()
