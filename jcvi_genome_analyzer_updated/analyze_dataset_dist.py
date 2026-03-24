import h5py
import numpy as np
import os

h5_path = "perturbation_dataset_50k.h5"
if not os.path.exists(h5_path):
    print(f"❌ Error: {h5_path} not found in {os.getcwd()}")
    sys.exit(1)

with h5py.File(h5_path, 'r') as f:
    growth = f['growth_rate'][:]
    viability = f['viability'][:]

print(f"Total samples: {len(growth)}")
print(f"Mean Growth Rate: {np.mean(growth):.6f}")
print(f"Median Growth Rate: {np.median(growth):.6f}")
print(f"Std Growth Rate: {np.std(growth):.6f}")
print(f"Max Growth Rate: {np.max(growth):.6f}")
print(f"Fraction Lethal (0 growth): {np.sum(viability == 0) / len(viability):.1%}")
print(f"Growth Rate Distribution (Quartiles): {np.percentile(growth, [25, 50, 75])}")
