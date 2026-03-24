import os
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader, random_split
import h5py
import numpy as np
from tqdm import tqdm

# Local imports
from neural_decoders import ResBatchMLP, save_model

class URDataset(Dataset):
    def __init__(self, h5_path):
        self.h5_path = h5_path
        with h5py.File(h5_path, 'r') as f:
            self.length = len(f['cell_ur'])
            
    def __len__(self):
        return self.length
    
    def __getitem__(self, idx):
        # We open file in __getitem__ to be thread-safe with workers
        with h5py.File(self.h5_path, 'r') as f:
            ur = f['cell_ur'][idx]
            growth = f['growth_rate'][idx]
            viability = f['viability'][idx]
            
        return torch.tensor(ur, dtype=torch.float32), \
               torch.tensor(growth, dtype=torch.float32), \
               torch.tensor(viability, dtype=torch.float32)

def train_epoch(model, loader, optimizer, criterion_growth, criterion_viability, device):
    model.train()
    total_loss = 0
    
    for ur, growth, viability in tqdm(loader, desc="Training"):
        ur, growth, viability = ur.to(device), growth.to(device), viability.to(device)
        
        optimizer.zero_grad()
        outputs = model(ur)
        
        # Multi-task loss
        loss_growth = criterion_growth(outputs['growth_rate'], growth)
        loss_viability = criterion_viability(outputs['viability_logits'], viability)
        
        # Weighting: Growth is more sensitive, Viability is stable
        loss = loss_growth + 0.1 * loss_viability
        
        loss.backward()
        optimizer.step()
        
        total_loss += loss.item()
        
    return total_loss / len(loader)

def validate(model, loader, criterion_growth, criterion_viability, device):
    model.eval()
    val_loss = 0
    with torch.no_grad():
        for ur, growth, viability in loader:
            ur, growth, viability = ur.to(device), growth.to(device), viability.to(device)
            outputs = model(ur)
            
            loss_growth = criterion_growth(outputs['growth_rate'], growth)
            loss_viability = criterion_viability(outputs['viability_logits'], viability)
            loss = loss_growth + 0.1 * loss_viability
            val_loss += loss.item()
            
    return val_loss / len(loader)

def main():
    # Setup Device
    device = torch.device('cuda:3' if torch.cuda.device_count() > 3 else 'cuda:0' if torch.cuda.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    # Load Dataset
    dataset_path = "perturbation_dataset_50k.h5"
    if not os.path.exists(dataset_path):
        print(f"Error: {dataset_path} not found. Run perturbation_engine.py first.")
        return
        
    full_ds = URDataset(dataset_path)
    train_size = int(0.8 * len(full_ds))
    val_size = len(full_ds) - train_size
    train_ds, val_ds = random_split(full_ds, [train_size, val_size])
    
    train_loader = DataLoader(train_ds, batch_size=64, shuffle=True, num_workers=0)
    val_loader = DataLoader(val_ds, batch_size=64, shuffle=False)
    
    # Initialize Model
    model = ResBatchMLP().to(device)
    optimizer = optim.AdamW(model.parameters(), lr=1e-4, weight_decay=1e-2)
    
    criterion_growth = nn.MSELoss()
    criterion_viability = nn.BCEWithLogitsLoss()
    
    # Training Loop
    epochs = 10
    best_val_loss = float('inf')
    
    for epoch in range(epochs):
        train_loss = train_epoch(model, train_loader, optimizer, criterion_growth, criterion_viability, device)
        val_loss = validate(model, val_loader, criterion_growth, criterion_viability, device)
        
        print(f"Epoch {epoch+1}/{epochs} - Train Loss: {train_loss:.6f} - Val Loss: {val_loss:.6f}")
        
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            save_model(model, "ai_virtual_instrument_best.pt")

if __name__ == "__main__":
    main()
