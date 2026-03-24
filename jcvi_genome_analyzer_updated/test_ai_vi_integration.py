import os
import sys
import numpy as np
import torch

# Ensure local imports work
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cellular_ur_generator import CellularURGenerator
from virtual_instruments import VirtualInstrumentSuite
from data_structures import Gene, CellStateVector

def test_ai_vi():
    print("=" * 60)
    print("TEST: AI VIRTUAL INSTRUMENT 2.0 VERIFICATION")
    print("=" * 60)
    
    # 1. Setup Generator
    generator = CellularURGenerator()
    
    # 2. Setup Suite (should auto-load ai_virtual_instrument_best.pt)
    vis = VirtualInstrumentSuite(generator)
    
    if not vis.ai_active:
        print("❌ AI Model not active! Check if ai_virtual_instrument_best.pt exists.")
        return
    
    # 3. Create a mock CellStateVector
    mock_ur = np.random.randn(1888).astype(np.float32)
    # Ensure no NaNs as we learned that was an issue
    mock_ur = np.nan_to_num(mock_ur)
    
    cell_state = CellStateVector(
        gene_embedding_aggregate=mock_ur[:1280],
        dynamics_vector=mock_ur[1280:1280+512],
        environment_vector=mock_ur[1280+512:1280+512+64],
        resource_vector=mock_ur[1280+512+64:],
        combined_cellular_ur=mock_ur,
        timestamp=0.0,
        num_genes=473,
        num_active_genes=473,
        growth_phase="exponential"
    )
    
    # 4. Predict Growth
    growth = vis.predict_growth_rate(cell_state)
    print(f"✅ AI Growth Rate: {growth['growth_rate']:.6f}")
    print(f"   Doubling Time: {growth['doubling_time']:.1f} min")
    print(f"   Source: {growth.get('prediction_source', 'heuristic')}")
    
    # 5. Predict Phenotype
    # We need some dummy genes for the heuristic part
    genes = {"dummy": Gene(id="dummy", name="dummy", start=0, end=100, strand="+", type="CDS")}
    pheno = vis.predict_phenotype(cell_state, genes)
    print(f"✅ AI Phenotype: {pheno['phenotype']}")
    print(f"   AI Viability Score: {pheno.get('viability_score_ai', 'N/A')}")
    print(f"   Source: {pheno.get('prediction_source', 'heuristic')}")
    
    assert "neural_decoder" in growth.get('prediction_source', ''), "AI was not used for growth prediction"
    assert "neural_decoder" in pheno.get('prediction_source', ''), "AI was not used for phenotype prediction"
    
    print("\n=" * 60)
    print("ALL AI INTEGRATION TESTS PASSED ✅✅✅")
    print("=" * 60)

if __name__ == "__main__":
    test_ai_vi()
