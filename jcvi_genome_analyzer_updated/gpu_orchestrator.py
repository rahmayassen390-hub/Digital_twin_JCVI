import torch
import threading
import time
from typing import Dict, Any, Optional, List

try:
    from neural_networks import PromoterCNN, TranscriptionLSTM
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False

class GPUOrchestrator:
    """
    Manages internal GPU distribution for a single-window session.
    Assigns specific functional roles to dedicated RTX 5000 GPUs.
    """
    _instance = None
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(GPUOrchestrator, cls).__new__(cls)
            cls._instance.initialized = False
        return cls._instance

    def __init__(self):
        if self.initialized:
            return
            
        self.devices = {
            'motif': torch.device('cuda:0') if torch.cuda.device_count() > 0 else torch.device('cpu'),
            'dynamics': torch.device('cuda:1') if torch.cuda.device_count() > 1 else torch.device('cpu'),
            'synthesis': torch.device('cuda:2') if torch.cuda.device_count() > 2 else torch.device('cpu'),
            'master': torch.device('cuda:3') if torch.cuda.device_count() > 3 else torch.device('cpu')
        }
        
        self.models = {}
        self.threads = {}
        self.running = False
        self.initialized = True
        
        print(f"✅ GPU Orchestrator Initialized with {torch.cuda.device_count()} GPUs")

    def get_device_for_worker(self, worker_idx: int) -> torch.device:
        """Assigns a GPU in a round-robin fashion based on worker index"""
        num_gpus = torch.cuda.device_count()
        if num_gpus == 0:
            return torch.device('cpu')
        
        gpu_idx = worker_idx % num_gpus
        return torch.device(f'cuda:{gpu_idx}')

    def get_node_device(self, node_role: str) -> torch.device:
        """Returns the dedicated device for a specific functional node"""
        return self.devices.get(node_role, torch.device('cpu'))

    def start_nodes(self):
        """Pre-initialize models on their respective GPUs"""
        if not TORCH_AVAILABLE:
            return
            
        print("🚀 Initializing Specialized GPU Nodes...")
        
        # 1. Motif Node (GPU 0)
        try:
            self.models['motif'] = PromoterCNN(device=self.devices['motif'])
            print(f"  🧬 Motif Node ready on {self.devices['motif']}")
        except Exception as e:
            print(f"  ❌ Failed to init Motif Node: {e}")

        # 2. Dynamics Node (GPU 1)
        try:
            self.models['dynamics'] = TranscriptionLSTM(device=self.devices['dynamics'])
            print(f"  ⏳ Dynamics Node ready on {self.devices['dynamics']}")
        except Exception as e:
            print(f"  ❌ Failed to init Dynamics Node: {e}")

        # 3. Synthesis Node (GPU 2)
        # (ProteinSynthesisWidget handles its own rendering, but we can offload computation here)
        print(f"  🔬 Synthesis Node affinity set to {self.devices['synthesis']}")
        
        # 4. Master Node (GPU 3)
        print(f"  🪐 Master Node affinity set to {self.devices['master']}")

    def run_parallel_analysis(self, sequence_data: Any):
        """Example of triggering parallel tasks across GPUs"""
        def motif_task():
            # Mock task on GPU 0
            with torch.cuda.device(self.devices['motif']):
                # Perform CNN inference
                pass
                
        def dynamics_task():
            # Mock task on GPU 1
            with torch.cuda.device(self.devices['dynamics']):
                # Perform LSTM prediction
                pass

        threading.Thread(target=motif_task).start()
        threading.Thread(target=dynamics_task).start()

    def get_hpc_metrics(self) -> Dict[str, Any]:
        """
        Unified collection of HPC performance metrics.
        Returns detailed stats for CPU, RAM, and all available GPUs.
        """
        import subprocess
        import psutil
        
        metrics = {
            'cpu_util': psutil.cpu_percent(),
            'ram_util': psutil.virtual_memory().percent,
            'gpus': []
        }
        
        try:
            # Query nvidia-smi for index, temperature, and utilization
            cmd = ["nvidia-smi", "--query-gpu=index,temperature.gpu,utilization.gpu", "--format=csv,noheader,nounits"]
            result = subprocess.check_output(cmd, encoding='utf-8')
            
            for line in result.strip().split('\n'):
                if line.strip():
                    idx, temp, util = map(int, line.split(','))
                    metrics['gpus'].append({
                        'index': idx,
                        'temp': temp,
                        'util': util
                    })
        except Exception as e:
            print(f"⚠️ Warning: Could not fetch GPU metrics: {e}")
            # Fallback for display if no GPUs found
            if not metrics['gpus']:
                for i in range(4):
                    metrics['gpus'].append({'index': i, 'temp': 0, 'util': 0})
            
        return metrics

    def get_gpu_temperatures(self) -> Dict[int, int]:
        """Legacy compatibility wrapper for thermal indicators"""
        metrics = self.get_hpc_metrics()
        return {g['index']: g['temp'] for g in metrics['gpus']}

    def initialize_ddp(self, rank: int, world_size: int):
        """Initialize Distributed Data Parallel for training nodes"""
        if not TORCH_AVAILABLE:
            return
            
        import os
        os.environ['MASTER_ADDR'] = 'localhost'
        os.environ['MASTER_PORT'] = '12355'
        
        # Initialize the process group
        torch.distributed.init_process_group("nccl", rank=rank, world_size=world_size)
        torch.cuda.set_device(rank)
        print(f"📡 DDP Node {rank}/{world_size} initialized on {torch.cuda.get_device_name(rank)}")

    def get_training_subset(self) -> List[torch.device]:
        """Returns GPUs dedicated to the training cluster (GPUs 1-3)"""
        return [self.devices['dynamics'], self.devices['synthesis'], self.devices['master']]

    def cleanup_ddp(self):
        """Cleanup DDP processes"""
        if torch.distributed.is_initialized():
            torch.distributed.destroy_process_group()

orchestrator = GPUOrchestrator()
