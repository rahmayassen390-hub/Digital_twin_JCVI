<p align="center">
  <h1 align="center">🧬 JCVI Genome Analyzer — Digital Twin of the Minimal Cell</h1>
  <p align="center">
    <em>A GPU-accelerated, quantum-enhanced desktop application for simulating and analyzing the JCVI-syn3A minimal synthetic genome</em>
  </p>
  <p align="center">
    <img src="https://img.shields.io/badge/python-3.10+-blue?logo=python&logoColor=white" alt="Python 3.10+"/>
    <img src="https://img.shields.io/badge/PyQt5-Desktop_GUI-green?logo=qt&logoColor=white" alt="PyQt5"/>
    <img src="https://img.shields.io/badge/PyTorch-GPU_Accelerated-red?logo=pytorch&logoColor=white" alt="PyTorch"/>
    <img src="https://img.shields.io/badge/OQTOPUS-Quantum_Backend-purple" alt="OQTOPUS"/>
    <img src="https://img.shields.io/badge/CUDA-Quad_GPU-76B900?logo=nvidia&logoColor=white" alt="CUDA"/>
    <img src="https://img.shields.io/badge/License-Apache_2.0-orange" alt="License"/>
  </p>
</p>

---

## 👥 Authors

| Author | GitHub |
|---|---|
| **Rahma Yassen** | [@rahmayassen390-hub](https://github.com/rahmayassen390-hub) |
| **Karma** | [@karma119](https://github.com/karma119) |
| **Mohamed S.** | [@mohamedss24](https://github.com/mohamedss24) |
| **Radwa Kamal** | [@RadwaKamalHamdyDeif](https://github.com/RadwaKamalHamdyDeif) |
| **Mohamed H. Lotfy** | [@mhlotfy24-3D](https://github.com/mhlotfy24-3D) |
| **Ahmed E. Gomaa** | [@Ahmed-EGomaa](https://github.com/Ahmed-EGomaa) |

---

## 📖 Overview

**JCVI Genome Analyzer** is a professional-grade **Digital Twin** of the *Mycoplasma mycoides* JCVI-syn3A minimal cell (493 genes, 543 kbp). It combines real-time ODE-based cell simulation, AI-driven gene analysis, quad-GPU distribution, and quantum circuit simulation to provide an unprecedented window into minimal cell biology.

The application enables researchers to:

- **Load & compare** JCVI-syn1.0 (reference) and JCVI-syn3.0 (query) genomes via BLAST
- **Simulate** continuous mRNA/protein dynamics using mass-balance ODEs (RK4 integration)
- **Predict** growth rate, phenotype, essentiality, and stress response via AI Virtual Instruments
- **Visualize** gene knockout cascade effects on the entire cell
- **Run quantum simulations** of gene-cluster entanglement via the OQTOPUS cloud backend
- **Distribute workloads** across 4 NVIDIA GPUs for motif detection, dynamics, synthesis, and orchestration

---

## ✨ Key Features

### 🔬 Cell Simulation Engine
- Continuous-time ODE simulation with **Runge-Kutta 4th-order** integration
- JCVI-syn3A calibrated kinetics (105 min doubling time, 5 min mRNA half-life)
- Genome-derived resource pools (ribosomes, RNA polymerases, ATP synthase)
- Biological noise (stochastic fluctuations) and metabolic coupling via FBA
- Individual gene pause/resume and knockout capabilities

### 🤖 AI Agents & Virtual Instruments
| Agent / Instrument | Function |
|---|---|
| **AI Promoter Agent** | CNN-based promoter motif detection & classification |
| **AI Transcription Agent** | LSTM-based transcription dynamics prediction |
| **AI Translation Agent** | Translation efficiency & protein synthesis modeling |
| **Growth Rate Predictor** | Predicts doubling time from cellular Universal Representation (UR) |
| **Essentiality Predictor** | Multi-feature gene essentiality classification (>85% accuracy) |
| **Phenotype Predictor** | Cell viability & phenotype scoring |
| **Stress Response Predictor** | Environmental stress tolerance estimation |
| **Knockout Simulator** | Expression-weighted knockout impact scoring |
| **Drug Treatment Simulator** | Drug target attenuation modeling |
| **Neural Decoder (v2.0)** | Trained MLP for hybrid AI-heuristic predictions |

### ⚡ Quad-GPU Distribution
| GPU | Role | Function |
|---|---|---|
| GPU 0 | **Motif Node** | PromoterCNN inference |
| GPU 1 | **Dynamics Node** | TranscriptionLSTM prediction |
| GPU 2 | **Synthesis Node** | Protein synthesis computation |
| GPU 3 | **Master Node** | Orchestration & aggregation |

### 🔮 Quantum Simulation (OQTOPUS)
- **Single-gene circuits**: RY rotation encoding transcription rate → measurement probability
- **28-qubit cluster circuits**: Entangled gene clusters with CNOT topology
- **OpenQASM 3.0** circuit generation, submitted to OQTOPUS cloud backend (SVSim device)
- Docker-based local backend with MySQL, MinIO, and mock simulator worker

### 📊 Interactive Dashboards
- **Live Transcription Animation** — Real-time gene expression with colored status indicators
- **Translation Window** — Separate translation dynamics view
- **Gene Dynamics Panel** — 4-tab detail view (mRNA TX, Protein TL, mRNA Degradation, Protein Stability)
- **Metabolic Dashboard** — Real-time FBA visualization
- **Knockout Effects Dashboard** — Gene dependency graph with cascade visualization
- **Resource Competition Dashboard** — Ribosome/polymerase saturation monitoring
- **HPC Monitor** — GPU temperature, utilization, and CPU/RAM metrics
- **Live Data Viewer** — ZeroMQ-based real-time data broadcasting (port 5555)

### 🧠 Smart Data Engine
- Intelligent multi-source data integration (GFF3 + 6 BED file types)
- Automatic category derivation from gene product descriptions
- "Thinking" gap-filling: never leaves "Unknown" if data exists in any source file
- Complete `IntegratedGene` profiles with regulation, transcription, translation, and kinetics

---

## 📁 Project Structure

```
jcvi_genome_analyzer/
│
├── jcvi_genome_analyzer_updated/     # Main application package
│   ├── main.py                       # 🚀 Entry point
│   ├── gui_main.py                   # Main GUI window (PyQt5)
│   ├── gui_methods.py                # GUI event handlers & methods
│   ├── data_structures.py            # Core data classes (Gene, BlastHit, etc.)
│   ├── neural_networks.py            # AI models (PromoterCNN, TranscriptionLSTM)
│   ├── ai_promoter_agent.py          # AI Agent #1: Promoter analysis
│   ├── ai_transcription_agent.py     # AI Agent #2: Transcription dynamics
│   ├── ai_translation_agent.py       # AI Agent #3: Translation dynamics
│   ├── blast_manager.py              # BLAST operations manager
│   ├── genome_manager.py             # Genome data loading (FASTA/GFF3/BED)
│   ├── workers.py                    # Background worker threads
│   ├── cell_simulation_engine.py     # ODE cell simulation (RK4)
│   ├── smart_data_engine.py          # Multi-source data integration
│   ├── data_enrichment_layer.py      # Gene data enrichment (CSV datasets)
│   ├── data_integrator.py            # Data integration pipeline
│   ├── data_auto_loader.py           # Auto-discovery of Data/ files
│   ├── molecular_ur_generator.py     # Molecular Universal Representation
│   ├── cellular_ur_generator.py      # Cellular UR aggregation
│   ├── virtual_instruments.py        # AI Virtual Instruments suite
│   ├── neural_decoders.py            # Trained MLP models (ResBatchMLP)
│   ├── gpu_orchestrator.py           # Quad-GPU distribution (Singleton)
│   ├── perturbation_engine.py        # Batch perturbation dataset generation
│   ├── sensitivity_analysis.py       # UR sensitivity & stability analysis
│   ├── quantum_worker.py             # OQTOPUS quantum job workers
│   ├── cluster_definitions.py        # 28 functional gene clusters
│   ├── knockout_rules_engine.py      # Cascade knockout rules
│   ├── knockout_effects_dashboard.py # Knockout visualization
│   ├── gene_dependency_graph.py      # Gene-gene dependency graph
│   ├── gene_dynamics_panel.py        # 4-tab gene dynamics dialogs
│   ├── metabolic_solver.py           # FBA metabolic solver
│   ├── metabolic_dashboard.py        # Metabolic visualization
│   ├── resource_dashboard.py         # Ribosome/polymerase dashboard
│   ├── protein_synthesis_widget.py   # Amino acid generator integration
│   ├── translation_window.py         # Separate translation window
│   ├── live_data_viewer.py           # ZeroMQ real-time data viewer
│   ├── hpc_monitor.py                # HPC performance monitoring
│   └── Data/                         # Genome data files
│
├── oqtopus-cloud/                    # OQTOPUS quantum backend (Docker)
├── oqtopus-setup/                    # OQTOPUS setup & configuration
├── amino_acid_generator/             # Standalone amino acid generator
├── draft_paper/                      # Research paper draft & figures
│
├── mock_simulator_worker.py          # Mock quantum simulator
├── standalone_quantum_api.py         # Standalone quantum API server
├── start_all_services.sh             # Master startup script
├── requirements.txt                  # Python dependencies
└── Data/                             # Root-level genome datasets
```

---

## 🚀 Getting Started

### Prerequisites

- **Python 3.10+**
- **NVIDIA GPUs** (4× recommended; falls back to CPU)
- **CUDA Toolkit** (for GPU acceleration)
- **NCBI BLAST+** (for genome comparison)
- **Docker & Docker Compose** (optional, for OQTOPUS quantum backend)

### Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/rahmayassen390-hub/Digital_twin_JCVI.git
   cd jcvi-genome-analyzer
   ```

2. **Create a virtual environment**
   ```bash
   python -m venv venv
   source venv/bin/activate  # Linux/macOS
   ```

3. **Install Python dependencies**
   ```bash
   pip install -r requirements.txt
   ```

4. **Install NCBI BLAST+ (required)**
   ```bash
   # Ubuntu/Debian
   sudo apt install ncbi-blast+

   # Or via Conda
   conda install -c bioconda blast

   # Verify
   blastn -version
   ```

5. **Set up OQTOPUS Quantum Backend (optional)**
   ```bash
   # Install Docker
   sudo apt install docker.io docker-compose

   # Start all OQTOPUS services
   bash start_all_services.sh
   ```

### Running the Application

```bash
cd jcvi_genome_analyzer_updated
python main.py
```

The application will:
1. Initialize the GPU Orchestrator (detects 4 GPUs or falls back to CPU)
2. Auto-discover genome data files from the `Data/` directory
3. Automatically load JCVI-syn1.0 (reference) and syn3.0 (query) genomes
4. Launch the interactive PyQt5 GUI

---

## 📋 Dependencies

### Core (Required)
| Package | Version | Purpose |
|---|---|---|
| PyQt5 | ≥ 5.15 | Desktop GUI framework |
| BioPython | ≥ 1.79 | Genome file parsing (FASTA, GFF3) |
| requests | ≥ 2.28 | HTTP client for OQTOPUS API |

### AI & Data Science (Recommended)
| Package | Version | Purpose |
|---|---|---|
| NumPy | ≥ 1.21 | Numerical computation |
| Pandas | ≥ 1.3 | Data analysis & CSV handling |
| PyTorch | ≥ 1.9 | Neural networks (CNN, LSTM, MLP) |
| openpyxl | ≥ 3.0 | Excel file support |

### Optional
| Package | Version | Purpose |
|---|---|---|
| pyzmq | ≥ 22.0 | ZeroMQ real-time data broadcasting |
| Flask | ≥ 2.0 | Standalone quantum API server |
| psutil | ≥ 5.9 | System monitoring & benchmarking |
| h5py | — | HDF5 perturbation dataset I/O |
| tqdm | — | Progress bars for batch processing |
| scikit-learn | — | ML utilities |

### System-Level
| Tool | Required | Purpose |
|---|---|---|
| NCBI BLAST+ | ✅ Yes | Genome sequence comparison |
| Docker | Optional | OQTOPUS quantum backend |
| NVIDIA CUDA | Optional | GPU acceleration |

---

## 🧪 Usage Guide

### 1. Genome Loading
The application automatically discovers and loads genome files from the `Data/` directory:
- **JCVI-syn1.0** (Reference): FASTA + GFF3 + BED annotation files
- **JCVI-syn3.0** (Query): FASTA + optional annotations

### 2. BLAST Comparison
Navigate to the **BLAST Search** tab and click **Run BLAST Analysis** to compare query vs. reference genomes.

### 3. AI Analysis
After BLAST completes, click **Run AI Analysis** to activate:
- **Promoter Analysis** — Identifies promoter motifs and classifies regulation types
- **Transcription Animation** — Real-time visualization of mRNA dynamics
- **Translation Window** — Protein synthesis monitoring

### 4. Live Simulation
The Cell Simulation Engine runs continuously, modeling:
- mRNA transcription and degradation (half-life: ~5 min)
- Protein translation and degradation (half-life: ~1 hr)
- Resource competition (ribosomes, RNA polymerases, ATP)

### 5. Knockout Experiments
Use the **Knockout Effects Dashboard** to:
- Select genes for knockout
- Visualize cascade effects on dependent genes
- Monitor growth rate impact in real-time

### 6. Quantum Simulation
With OQTOPUS running, the quantum features enable:
- Per-gene quantum probability estimation
- 28-qubit entangled cluster simulation
- Measurement-based transcription level validation

---

## 🔬 Scientific Background

This project models the **JCVI-syn3A** genome — the smallest known genome capable of sustaining life (Hutchison et al., 2016). Key biological parameters are calibrated against published literature:

| Parameter | Value | Source |
|---|---|---|
| Genome Size | 543,379 bp | CP016816 (NCBI) |
| Gene Count | ~493 protein-coding | Hutchison et al., 2016 |
| Doubling Time | ~105 min | Breuer et al., 2019 |
| mRNA Half-Life | ~5 min | Minimal cell estimate |
| Protein Half-Life | ~1 hr | Minimal cell estimate |
| Ribosomes | ~4,500 | Calibrated for 105 min doubling |
| RNA Polymerases | ~1,200 | Scaled from rpo gene count |

### References
1. Hutchison, C.A. et al. (2016). Design and synthesis of a minimal bacterial genome. *Science*, 351(6280).
2. Breuer, M. et al. (2019). Essential metabolism for a minimal cell. *eLife*, 8.
3. Thornburg, Z.R. et al. (2022). Fundamental behaviors emerge from simulations of a living minimal cell. *Cell*, 185(2).

---

## 🏗️ Architecture

```
┌────────────────────────────────────────────────────────────────┐
│                    PyQt5 Desktop GUI                           │
│  ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────────────┐  │
│  │Reference │ │  Query   │ │  BLAST   │ │   AI Results     │  │
│  │ Genome   │ │  Genome  │ │  Search  │ │ (Promoter/TX/TL) │  │
│  └──────────┘ └──────────┘ └──────────┘ └──────────────────┘  │
├────────────────────────────────────────────────────────────────┤
│              Smart Data Engine + Data Enrichment               │
│           GFF3 ← Regulation ← TX ← TL ← Promoters            │
├────────────────────────────────────────────────────────────────┤
│         Cell Simulation Engine (RK4 ODE + Noise + FBA)         │
│    mRNA ⇄ Protein dynamics │ Ribosomes │ ATP │ Resources       │
├────────────────────────────────────────────────────────────────┤
│         AI Virtual Instruments (Decoder + Manipulator)         │
│    Growth │ Phenotype │ Essentiality │ Knockout │ Drug         │
├──────────────────────┬─────────────────────────────────────────┤
│  GPU Orchestrator    │        Quantum Worker (OQTOPUS)         │
│  ┌─────┐ ┌─────┐    │   ┌────────────────────────────────┐    │
│  │GPU 0│ │GPU 1│    │   │ OpenQASM 3.0 Circuit Gen       │    │
│  │Motif│ │Dyna │    │   │ 1-qubit (RY) & 28-qubit (CNOT) │    │
│  ├─────┤ ├─────┤    │   │ → OQTOPUS API (SVSim)          │    │
│  │GPU 2│ │GPU 3│    │   └────────────────────────────────┘    │
│  │Synth│ │Mstr │    │                                         │
│  └─────┘ └─────┘    │                                         │
└──────────────────────┴─────────────────────────────────────────┘
```

---

## 📝 Data Files

Place the following files in the `Data/` directory:

| File | Description |
|---|---|
| `JCVI_1.0.fasta` | JCVI-syn1.0 reference genome sequence |
| `JCVI_1.0.gff3` | JCVI-syn1.0 gene annotations |
| `JCVI_3.0.fasta` | JCVI-syn3.0 query genome sequence |
| `JCVI_3.0.gff3` | JCVI-syn3.0 gene annotations |
| `JCVI_Gene_Regulation_*.bed` | Gene regulation types & categories |
| `JCVI_Transcription_Dynamics_*.bed` | Transcription levels, order, timeline |
| `JCVI_Translation_Dynamics_*.bed` | Protein levels, RBS, stability |
| `Likely_Promoters_*.bed` | Promoter locations & motifs |
| `Internal_Operons_*.bed` | Operon structure |
| `JCVI_syn3A_full_genome_dataset.csv` | Comprehensive gene dataset |
| `syn3A_stress_response_data.csv` | Stress response data |

---

## 🧩 Running Tests & Benchmarks

```bash
# Run system verification
cd jcvi_genome_analyzer_updated
python system_verification.py

# Run sensitivity analysis
python sensitivity_analysis.py

# Run production validation
python production_validation.py

# Run perturbation batch (requires GPU)
python perturbation_engine.py

# Benchmark performance
cd ..
python perform_benchmarks.py
```

---

## 🤝 Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

---

## 📄 License

This project is licensed under the **Apache License 2.0** — see the [LICENSE](oqtopus-cloud/LICENSE) file for details.

---

## 🙏 Acknowledgments

- **JCVI** (J. Craig Venter Institute) — for the syn3A minimal cell genome
- **OQTOPUS** — open-source quantum cloud platform for quantum simulation backend
- **NCBI** — for BLAST+ sequence comparison tools
- **PyTorch** — for GPU-accelerated neural network inference

---

<p align="center">
  <em>Built with ❤️ for synthetic biology research</em>
</p>
