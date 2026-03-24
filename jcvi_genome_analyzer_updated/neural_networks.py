"""
Neural Network Models for JCVI Genome Analyzer
Contains CNN and LSTM models for promoter and transcription analysis
"""

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False
    # Create dummy classes for when PyTorch is not available
    class nn:
        class Module:
            pass


class PromoterCNN(nn.Module):
    """CNN for promoter sequence motif detection"""
    
    def __init__(self, seq_length=100, device='cpu'):
        super(PromoterCNN, self).__init__()
        self.device = torch.device(device) if TORCH_AVAILABLE else None
        
        # One-hot encoding: 4 channels for A, T, G, C
        self.conv1 = nn.Conv1d(4, 32, kernel_size=6, padding=2)
        self.conv2 = nn.Conv1d(32, 64, kernel_size=6, padding=2)
        self.conv3 = nn.Conv1d(64, 128, kernel_size=6, padding=2)
        
        self.pool = nn.MaxPool1d(2)
        self.dropout = nn.Dropout(0.3)
        
        # Dynamically calculate the flattened size
        self.seq_length = seq_length
        self._to_linear = None
        
        self.fc1 = None
        self.fc2 = nn.Linear(256, 64)
        self.fc3 = nn.Linear(64, 3)  # 3 classes: own_promoter, operon_internal, ambiguous
        
        if TORCH_AVAILABLE:
            self.to(self.device)
    
    def _get_conv_output(self, shape):
        """Calculate output size after convolutions"""
        with torch.no_grad():
            dummy_input = torch.zeros(1, *shape).to(self.device)
            x = self.pool(F.relu(self.conv1(dummy_input)))
            x = self.pool(F.relu(self.conv2(x)))
            x = self.pool(F.relu(self.conv3(x)))
            return int(torch.numel(x) / x.size(0))
        
    def forward(self, x):
        # Move input to device
        x = x.to(self.device)
        
        # Initialize fc1 on first forward pass
        if self._to_linear is None:
            self._to_linear = self._get_conv_output((4, self.seq_length))
            self.fc1 = nn.Linear(self._to_linear, 256).to(self.device)
        
        x = self.pool(F.relu(self.conv1(x)))
        x = self.pool(F.relu(self.conv2(x)))
        x = self.pool(F.relu(self.conv3(x)))
        
        x = x.view(x.size(0), -1)
        x = self.dropout(F.relu(self.fc1(x)))
        x = self.dropout(F.relu(self.fc2(x)))
        x = self.fc3(x)
        
        return x


class TranscriptionLSTM(nn.Module):
    """LSTM for transcription state and dynamics prediction"""
    
    def __init__(self, input_size=20, hidden_size=128, num_layers=2, device='cpu'):
        super(TranscriptionLSTM, self).__init__()
        self.device = torch.device(device) if TORCH_AVAILABLE else None
        
        self.lstm = nn.LSTM(input_size, hidden_size, num_layers, 
                           batch_first=True, dropout=0.3)
        
        self.fc_state = nn.Linear(hidden_size, 4)  # 4 transcription states
        self.fc_level = nn.Linear(hidden_size, 1)  # Expression level
        self.fc_timing = nn.Linear(hidden_size, 1)  # Temporal order
        
        if TORCH_AVAILABLE:
            self.to(self.device)
        
    def forward(self, x):
        x = x.to(self.device)
        lstm_out, _ = self.lstm(x)
        last_output = lstm_out[:, -1, :]
        
        state = self.fc_state(last_output)
        level = torch.sigmoid(self.fc_level(last_output))
        timing = torch.relu(self.fc_timing(last_output))
        
        return state, level, timing
