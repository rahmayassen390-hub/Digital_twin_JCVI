#!/usr/bin/env python3
import sys
import os

# Check dependencies
try:
    from Bio import SeqIO
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False

try:
    import torch
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False

from PyQt5.QtWidgets import QApplication

# Import the main application class
from gui_main import GenomeAnalyzerApp
from gui_methods import GenomeAnalyzerMethods


# Create combined class with all functionality
class JCVIGenomeAnalyzer(GenomeAnalyzerApp, GenomeAnalyzerMethods):
    """Complete JCVI Genome Analyzer Application"""
    pass


def main():
    """Main entry point"""
    # Ensure all 4 GPUs are visible for the internal multi-GPU distribution
    # PCI_BUS_ID ensures that 0,1,2,3 map correctly to the physical cards
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    # Note: We don't force CUDA_VISIBLE_DEVICES here to allow system-level control,
    # but we assume 4 GPUs are available as per HPC specs.

    if not BIOPYTHON_AVAILABLE:
        print("ERROR: BioPython not installed!")
        print("Install: pip install biopython")
        return
    
    if not NUMPY_AVAILABLE:
        print("WARNING: NumPy not installed - AI features will be limited")
        print("Install: pip install numpy")
    
    if not TORCH_AVAILABLE:
        print("WARNING: PyTorch not installed - Neural networks will not be available")
        print("Install: pip install torch")
    
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    
    window = JCVIGenomeAnalyzer()
    window.setWindowTitle("🧬 JCVI Genome Analyzer v3.0 [Quad-GPU Enabled]")
    window.show()
    
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
