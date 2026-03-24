#!/usr/bin/env python3
"""
Ablation Test: Temporal Dynamics Contribution
==============================================
Evaluates the ResBatchMLP with the dynamics_vector zeroed out.
Proves that the temporal context (mRNA/protein states) is essential 
for accurate growth prediction.
"""

import os
import sys
import numpy as np
import h5py
import torch
import copy
from neural_decoders import ResBatchMLP, load_model
from sklearn.metrics import mean_squared_error, f1_score

def evaluate_ablation(model, ur, growth, viab, device, mask_dynamics=False):
    model.eval()
    test_ur = copy.deepcopy(ur)
    
    if mask_dynamics:
        # DYNAMICS_DIM is indices 1280 to 1791 (512 dims)
        test_ur[:, 1280:1280+512] = 0.0
        
    with torch.no_grad():
        inp = torch.tensor(test_ur, dtype=torch.float32).to(device)
        out = model(inp)
        pred_growth = out['growth_rate'].cpu().numpy()
        pred_viab_logits = out['viability_logits'].cpu().numpy()
        
    pred_viab_prob = 1.0 / (1.0 + np.exp(-pred_viab_logits))
    pred_viab_bin = (pred_viab_prob > 0.5).astype(float)
    
    rmse = np.sqrt(mean_squared_error(growth, pred_growth))
    f1 = f1_score(viab.astype(int), pred_viab_bin.astype(int), zero_division=0)
    
    return rmse, f1

if __name__ == "__main__":
    model_path = "ai_virtual_instrument_best.pt"
    dataset_path = "perturbation_dataset_50k.h5"
    
    if not os.path.exists(model_path) or not os.path.exists(dataset_path):
        print("❌ Error: Files missing")
        sys.exit(1)
        
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    model = load_model(ResBatchMLP(), model_path, device)
    
    with h5py.File(dataset_path, 'r') as f:
        all_ur = f['cell_ur'][:]
        all_growth = f['growth_rate'][:]
        all_viab = f['viability'][:]
        
    # Take 1000 samples for fast test
    test_idx = np.random.RandomState(42).permutation(len(all_ur))[:1000]
    ur_sample = all_ur[test_idx]
    growth_sample = all_growth[test_idx]
    viab_sample = all_viab[test_idx]
    
    print("\n" + "="*50)
    print("  NEURAL DECODER ABLATION TEST")
    print("="*50)
    
    # Baseline
    rmse_base, f1_base = evaluate_ablation(model, ur_sample, growth_sample, viab_sample, device, mask_dynamics=False)
    print(f"Baseline:    RMSE={rmse_base:.6f}, F1={f1_base:.4f}")
    
    # Ablated (No Dynamics)
    rmse_abl, f1_abl = evaluate_ablation(model, ur_sample, growth_sample, viab_sample, device, mask_dynamics=True)
    print(f"No Dynamics: RMSE={rmse_abl:.6f}, F1={f1_abl:.4f}")
    
    # Delta
    rmse_inc = (rmse_abl - rmse_base) / rmse_base * 100
    f1_dec = (f1_base - f1_abl) / f1_base * 100
    
    print("-" * 50)
    print(f"RMSE Increase (Error): {rmse_inc:+.1f}%")
    print(f"F1 Accuracy Decrease:  {f1_dec:+.1f}%")
    
    if rmse_inc > 10:
        print("\n✅ SUCCESS: Temporal dynamics are proven essential for prediction.")
    else:
        print("\n❌ WARNING: Minimal impact from dynamics vector.")
    print("="*50 + "\n")
