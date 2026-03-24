import torch
import torch.nn as nn
import torch.nn.functional as F

class ResidualBlock(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(dim, dim),
            nn.BatchNorm1d(dim),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(dim, dim),
            nn.BatchNorm1d(dim)
        )
    
    def forward(self, x):
        return F.relu(x + self.net(x))

class ResBatchMLP(nn.Module):
    """
    Multi-Task Neural Decoder for Cellular UR.
    Maps 1888-dim UR to Growth Rate (Regression) and Viability (Classification).
    """
    def __init__(self, input_dim=1888, hidden_dim=1024, n_blocks=3):
        super().__init__()
        
        # 1. Input Projection
        self.input_proj = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU()
        )
        
        # 2. Shared Residual Trunk
        self.trunk = nn.Sequential(*[
            ResidualBlock(hidden_dim) for _ in range(n_blocks)
        ])
        
        # 3. Growth Regression Head
        self.growth_head = nn.Sequential(
            nn.Linear(hidden_dim, 256),
            nn.ReLU(),
            nn.Linear(256, 64),
            nn.ReLU(),
            nn.Linear(64, 1) # Output: Growth Rate
        )
        
        # 4. Viability Classification Head
        self.viability_head = nn.Sequential(
            nn.Linear(hidden_dim, 256),
            nn.ReLU(),
            nn.Linear(256, 1) # Output: Logit for Sigmoid
        )
        
        print(f"✅ ResBatchMLP Initialized: {input_dim} -> {hidden_dim} (x{n_blocks})")

    def forward(self, x):
        features = self.input_proj(x)
        features = self.trunk(features)
        
        growth = self.growth_head(features)
        viability_logits = self.viability_head(features)
        
        return {
            "growth_rate": growth.squeeze(-1),
            "viability_logits": viability_logits.squeeze(-1)
        }

def save_model(model, path):
    torch.save(model.state_dict(), path)
    print(f"Model saved to {path}")

def load_model(model, path, device='cpu'):
    model.load_state_dict(torch.load(path, map_location=device))
    model.to(device)
    print(f"Model loaded from {path}")
    return model

if __name__ == "__main__":
    # Test architecture
    model = ResBatchMLP()
    test_input = torch.randn(16, 1888)
    output = model(test_input)
    print(f"Output growth shape: {output['growth_rate'].shape}")
    print(f"Output viability shape: {output['viability_logits'].shape}")
