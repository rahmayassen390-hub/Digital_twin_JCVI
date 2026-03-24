"""
Main GUI Application for JCVI Genome Analyzer - Part 1
Contains the main window class and UI setup
"""

import os
from collections import defaultdict

from PyQt5.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QTextEdit, QFileDialog, QProgressBar,
    QMessageBox, QTabWidget, QGroupBox,
    QTableWidget, QTableWidgetItem, QSpinBox, QDoubleSpinBox,
    QFormLayout, QSplitter, QScrollArea
)
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QFont, QColor

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False

try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False

from genome_manager import GenomeManager
from blast_manager import BlastManager
from workers import LoadGenomeWorker, BlastSearchWorker, AIAnalysisWorker
from live_data_viewer import LiveDataViewerDialog
from protein_synthesis_widget import ProteinSynthesisWidget
from gpu_orchestrator import orchestrator
from data_auto_loader import auto_discover_data

# Knockout Effects imports
try:
    from gene_dependency_graph import GeneDependencyGraph, get_effect_color, get_effect_icon
    from knockout_effects_dashboard import KnockoutEffectsDashboard
    KNOCKOUT_EFFECTS_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Knockout effects modules not available: {e}")
    KNOCKOUT_EFFECTS_AVAILABLE = False

# Metabolic Dashboard
try:
    from metabolic_dashboard import MetabolicDashboard
    METABOLIC_DASHBOARD_AVAILABLE = True
except ImportError:
    METABOLIC_DASHBOARD_AVAILABLE = False
    print("Warning: metabolic_dashboard not available")

# Gene Dynamics Panel (4-tab gene detail view)
try:
    from gene_dynamics_panel import GeneDynamicsPanel
    GENE_DYNAMICS_AVAILABLE = True
except ImportError:
    GENE_DYNAMICS_AVAILABLE = False
    print("Warning: gene_dynamics_panel not available")
    
# HPC Monitor Component
try:
    from hpc_monitor import HPCMonitorWidget
    HPC_MONITOR_AVAILABLE = True
except ImportError:
    HPC_MONITOR_AVAILABLE = False
    print("Warning: hpc_monitor not available")

# Resource Competition Dashboard
try:
    from resource_dashboard import ResourceDashboard
    RESOURCE_DASHBOARD_AVAILABLE = True
except ImportError:
    RESOURCE_DASHBOARD_AVAILABLE = False
    print("Warning: resource_dashboard not available")

# Enhanced Knockout Dialog
try:
    from enhanced_knockout_dialog import EnhancedKnockoutDialog
    ENHANCED_KNOCKOUT_AVAILABLE = True
except ImportError:
    ENHANCED_KNOCKOUT_AVAILABLE = False
    print("Warning: enhanced_knockout_dialog not available")

# ZeroMQ for real-time data broadcasting
try:
    import zmq
    import json
    ZMQ_AVAILABLE = True
except ImportError:
    ZMQ_AVAILABLE = False
    print("Warning: pyzmq not installed. Real-time broadcasting disabled.")
    print("Install with: pip install pyzmq")


class GenomeAnalyzerApp(QMainWindow):
    """Main application window"""
    
    def __init__(self):
        super().__init__()
        
        # Initialize GPU Orchestrator for 4-GPU Distribution
        self.gpu_orchestrator = orchestrator
        self.gpu_orchestrator.start_nodes()
        
        self.genome_manager = GenomeManager()
        self.query_genome_manager = GenomeManager()
        self.blast_manager = BlastManager()
        
        self.fasta_path = None
        self.gff_path = None
        self.likely_promoters_path = None
        self.internal_operons_path = None
        self.ref_complex_cases_path = None
        self.ref_regulation_path = None
        self.ref_transcription_path = None
        self.ref_translation_path = None
        
        self.query_fasta_path = None
        self.query_gff_path = None
        self.query_likely_promoters_path = None
        self.query_internal_operons_path = None
        self.query_complex_cases_path = None
        self.query_regulation_path = None
        self.query_transcription_path = None
        self.query_translation_path = None
        
        self.search_results = []
        self.ai_results = {}
        
        # Animation state - Transcription
        self.animation_running = False
        self.animation_time = 0.0
        self.animation_timer = None
        self.gene_transcription_data = []
        
        # Animation state - Translation
        self.translation_animation_running = False
        self.translation_animation_time = 0.0
        self.translation_animation_timer = None
        self.gene_translation_data = []
        
        # Paused genes tracking for individual gene control
        self.paused_genes = set()  # Set of gene indices that are paused
        self.gene_pause_buttons = {}  # Map gene index to pause button
        
        # Quantum Simulation state
        self.quantum_jobs = {}  # Map gene_id to QuantumJobWorker
        self.quantum_genes_count = 0  # Limit to 10 genes
        
        # === KNOCKOUT EFFECTS STATE ===
        self.dependency_graph = None  # GeneDependencyGraph instance
        self.knockout_dashboard = None  # KnockoutEffectsDashboard instance
        self.affected_genes_map = {}  # knocked_gene_id -> list of affected gene_ids
        self.gene_row_highlights = {}  # gene_id -> original background color (for restoration)
        
        # Create separate translation window (initially hidden)
        from translation_window import TranslationWindow
        self.translation_window = TranslationWindow(self)
        
        # Connect translation window controls to methods
        self.translation_window.start_btn.clicked.connect(self.start_translation_animation)
        self.translation_window.stop_btn.clicked.connect(self.stop_translation_animation)
        self.translation_window.reset_btn.clicked.connect(self.reset_translation_animation)
        
        # Use translation window's widgets
        self.translation_table = self.translation_window.translation_table
        self.transl_time_label = self.translation_window.time_label
        self.transl_speed_spinbox = self.translation_window.speed_spinbox
        self.start_translation_btn = self.translation_window.start_btn
        self.stop_translation_btn = self.translation_window.stop_btn
        self.reset_translation_btn = self.translation_window.reset_btn
        
        # Connect translation table click to gene dynamics panel
        self.translation_table.cellClicked.connect(self._on_translation_table_clicked)
        
        # Create Live Data Viewer dialog (shows translation data in real-time)
        self.live_data_viewer = LiveDataViewerDialog(self)
        
        # Create Metabolic Dashboard (real-time FBA visualization)
        self.metabolic_dashboard = None
        if METABOLIC_DASHBOARD_AVAILABLE:
            self.metabolic_dashboard = MetabolicDashboard(self)
        
        # Create Gene Dynamics Panel (4-tab gene detail dialog)
        self.gene_dynamics_panel = None
        if GENE_DYNAMICS_AVAILABLE:
            self.gene_dynamics_panel = GeneDynamicsPanel(self)
        
        # Create Resource Dashboard (ribosome/polymerase competition)
        self.resource_dashboard = None
        if RESOURCE_DASHBOARD_AVAILABLE:
            self.resource_dashboard = ResourceDashboard()
        
        # Create Protein Synthesis widget (integrated version of Amino Acid Generator)
        self.protein_synthesis_widget = ProteinSynthesisWidget(self)
        
        # Create HPC Monitor Sidebar
        self.hpc_monitor = None
        if HPC_MONITOR_AVAILABLE:
            self.hpc_monitor = HPCMonitorWidget(self)
        
        # Initialize ZeroMQ publisher for real-time broadcasting
        self.zmq_context = None
        self.zmq_socket = None
        self._init_zmq_publisher()
        
        # Thermal monitoring state
        self.gpu_temp_labels = {}
        self.thermal_timer = QTimer()
        self.thermal_timer.timeout.connect(self.update_gpu_thermals)
        self.thermal_timer.start(5000)  # Update every 5 seconds
        
        self.init_ui()
        
        # Auto-discover data files from /Data directory
        self._auto_load_from_data_dir()
        
        # Check BLAST
        blast_ok, blast_msg = self.blast_manager.check_blast()
        if not blast_ok:
            QMessageBox.warning(self, "BLAST Not Found", 
                              f"{blast_msg}\n\nPlease install BLAST to use this application.")
    
    def _auto_load_from_data_dir(self):
        """Auto-discover and populate file paths from the Data directory."""
        print("\n🚀 Starting Dual-Genome Auto-Discovery (1.0 first, then 3.0)...")
        
        # 1. Discover JCVI 1.0 (Reference)
        try:
            manifest_1_0 = auto_discover_data(genome_version="1.0")
            if manifest_1_0:
                ref_mapping = [
                    ('fasta_path',              'fasta_path',              'ref_fasta_label'),
                    ('gff3_path',               'gff_path',                'ref_gff_label'),
                    ('likely_promoters_path',   'likely_promoters_path',   'ref_likely_label'),
                    ('internal_operons_path',   'internal_operons_path',   'ref_operon_label'),
                    ('complex_cases_path',      'ref_complex_cases_path',  'ref_complex_label'),
                    ('regulation_bed_path',     'ref_regulation_path',     'ref_regulation_label'),
                    ('transcription_bed_path',  'ref_transcription_path',  'ref_transcription_label'),
                    ('translation_bed_path',    'ref_translation_path',    'ref_translation_label'),
                ]
                self._apply_manifest_mapping(manifest_1_0, ref_mapping, "Reference (1.0)")
        except Exception as e:
            print(f"Auto-discovery for 1.0 failed: {e}")
            
        # 2. Discover JCVI 3.0 (Query)
        try:
            manifest_3_0 = auto_discover_data(genome_version="3.0")
            if manifest_3_0:
                query_mapping = [
                    ('fasta_path',              'query_fasta_path',              'query_fasta_label'),
                    ('gff3_path',               'query_gff_path',                'query_gff_label'),
                    ('likely_promoters_path',   'query_likely_promoters_path',   'query_likely_label'),
                    ('internal_operons_path',   'query_internal_operons_path',   'query_operon_label'),
                    ('complex_cases_path',      'query_complex_cases_path',      'query_complex_label'),
                    ('regulation_bed_path',     'query_regulation_path',         'query_regulation_label'),
                    ('transcription_bed_path',  'query_transcription_path',      'query_transcription_label'),
                    ('translation_bed_path',    'query_translation_path',      'query_translation_label'),
                ]
                self._apply_manifest_mapping(manifest_3_0, query_mapping, "Query (3.0)")
        except Exception as e:
            print(f"Auto-discovery for 3.0 failed: {e}")

        # 3. AUTOMATION: Kick off the sequential load (1.0 first, which then triggers 3.0)
        if self.fasta_path and self.gff_path:
            print("🚀 Paths discovered. Automatically starting JCVI 1.0 (Reference) load...")
            QTimer.singleShot(500, self.load_reference_genome)

    def _apply_manifest_mapping(self, manifest, mapping, version_name):
        """Helper to apply discovered manifest to GUI attributes and labels."""
        loaded_count = 0
        for manifest_attr, self_attr, label_attr in mapping:
            path = getattr(manifest, manifest_attr, None)
            if path and os.path.exists(path):
                setattr(self, self_attr, path)
                label = getattr(self, label_attr, None)
                if label:
                    label.setText(f"✅ {os.path.basename(path)}")
                loaded_count += 1
        
        if loaded_count > 0:
            print(f"✅ Auto-loaded {loaded_count} files for {version_name}")
            if "Reference" in version_name:
                self.ref_stats_text.setPlainText(
                    f"Auto-Discovery ({version_name}): {loaded_count} files loaded.\n"
                    f"Ready for analysis."
                )
            else:
                self.query_stats_text.setPlainText(
                    f"Auto-Discovery ({version_name}): {loaded_count} files loaded.\n"
                    f"Ready for comparison."
                )
    
    def _init_zmq_publisher(self):
        """Initialize ZeroMQ publisher socket for real-time broadcasting"""
        if not ZMQ_AVAILABLE:
            self.live_data_viewer.set_zmq_status("🔌 ZMQ: Not available (install pyzmq)")
            return
        
        try:
            self.zmq_context = zmq.Context()
            self.zmq_socket = self.zmq_context.socket(zmq.PUB)
            self.zmq_socket.bind("tcp://*:5555")
            self.live_data_viewer.set_zmq_status("🟢 ZMQ Publisher: Active on port 5555")
            print("DEBUG: ZMQ publisher initialized on port 5555")
        except Exception as e:
            print(f"Warning: Could not initialize ZMQ publisher: {e}")
            self.live_data_viewer.set_zmq_status(f"🔴 ZMQ Error: {e}")
            self.zmq_socket = None
    
    def closeEvent(self, event):
        """Clean up resources on app close"""
        # Stop any running animations first to prevent ZMQ errors
        if hasattr(self, 'translation_animation_timer') and self.translation_animation_timer:
            self.translation_animation_timer.stop()
        if hasattr(self, 'animation_timer') and self.animation_timer:
            self.animation_timer.stop()
        
        # Now safely close ZMQ
        if self.zmq_socket:
            self.zmq_socket.close()
        if self.zmq_context:
            self.zmq_context.term()
        print("DEBUG: ZMQ resources cleaned up")
        event.accept()
    
    def init_ui(self):
        self.setWindowTitle("JCVI Genome Analyzer with AI Agents")
        self.setGeometry(50, 50, 1200, 800)
        self.setMinimumSize(800, 600)
        
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        
        # Main Horizontal Layout (Tabs | HPC Monitor)
        content_layout = QHBoxLayout()
        main_layout.addLayout(content_layout)
        
        self.tabs = QTabWidget()
        content_layout.addWidget(self.tabs, 1) # Tabs take priority space
        
        if self.hpc_monitor:
            content_layout.addWidget(self.hpc_monitor)
        
        self.create_reference_tab()
        self.create_query_tab()
        self.create_search_tab()
        self.create_results_tab()
        self.create_ai_results_tab()
        
        # Initialize Status Bar with Thermal Monitoring
        self.init_status_bar()
        
        self.tabs.setTabEnabled(1, False)
        self.tabs.setTabEnabled(2, False)
        self.tabs.setTabEnabled(3, False)
        self.tabs.setTabEnabled(4, False)
    
    def create_reference_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        title = QLabel("📊 Reference Genome (JCVI-syn1.0)")
        title.setFont(QFont("Arial", 13, QFont.Bold))
        layout.addWidget(title)
        
        file_group = QGroupBox("Select Reference Files")
        file_layout = QFormLayout()
        
        fasta_layout = QHBoxLayout()
        self.ref_fasta_label = QLabel("No file selected")
        fasta_btn = QPushButton("Browse FASTA")
        fasta_btn.clicked.connect(self.select_ref_fasta)
        fasta_layout.addWidget(self.ref_fasta_label, 1)
        fasta_layout.addWidget(fasta_btn)
        file_layout.addRow("FASTA File:", fasta_layout)
        
        gff_layout = QHBoxLayout()
        self.ref_gff_label = QLabel("No file selected")
        gff_btn = QPushButton("Browse GFF3")
        gff_btn.clicked.connect(self.select_ref_gff)
        gff_layout.addWidget(self.ref_gff_label, 1)
        gff_layout.addWidget(gff_btn)
        file_layout.addRow("GFF3 File:", gff_layout)
        
        likely_layout = QHBoxLayout()
        self.ref_likely_label = QLabel("No file selected")
        likely_btn = QPushButton("Browse BED")
        likely_btn.clicked.connect(self.select_ref_likely_promoters)
        likely_layout.addWidget(self.ref_likely_label, 1)
        likely_layout.addWidget(likely_btn)
        file_layout.addRow("Likely Promoters (BED):", likely_layout)
        
        operon_layout = QHBoxLayout()
        self.ref_operon_label = QLabel("No file selected")
        operon_btn = QPushButton("Browse BED")
        operon_btn.clicked.connect(self.select_ref_internal_operons)
        operon_layout.addWidget(self.ref_operon_label, 1)
        operon_layout.addWidget(operon_btn)
        file_layout.addRow("Internal Operons (BED):", operon_layout)
        
        # === NEW: Additional BED files for Reference ===
        
        # Complex Cases BED
        ref_complex_layout = QHBoxLayout()
        self.ref_complex_label = QLabel("No file selected (optional)")
        ref_complex_btn = QPushButton("Browse BED")
        ref_complex_btn.clicked.connect(self.select_ref_complex_cases)
        ref_complex_layout.addWidget(self.ref_complex_label, 1)
        ref_complex_layout.addWidget(ref_complex_btn)
        file_layout.addRow("Complex Cases (BED):", ref_complex_layout)
        
        # Regulation Analysis BED
        ref_regulation_layout = QHBoxLayout()
        self.ref_regulation_label = QLabel("No file selected")
        ref_regulation_btn = QPushButton("Browse BED")
        ref_regulation_btn.clicked.connect(self.select_ref_regulation)
        ref_regulation_layout.addWidget(self.ref_regulation_label, 1)
        ref_regulation_layout.addWidget(ref_regulation_btn)
        file_layout.addRow("Regulation Analysis (BED):", ref_regulation_layout)
        
        # Transcription Dynamics BED
        ref_transcription_layout = QHBoxLayout()
        self.ref_transcription_label = QLabel("No file selected")
        ref_transcription_btn = QPushButton("Browse BED")
        ref_transcription_btn.clicked.connect(self.select_ref_transcription)
        ref_transcription_layout.addWidget(self.ref_transcription_label, 1)
        ref_transcription_layout.addWidget(ref_transcription_btn)
        file_layout.addRow("Transcription Dynamics (BED):", ref_transcription_layout)
        
        # Translation Dynamics BED
        ref_translation_layout = QHBoxLayout()
        self.ref_translation_label = QLabel("No file selected")
        ref_translation_btn = QPushButton("Browse BED")
        ref_translation_btn.clicked.connect(self.select_ref_translation)
        ref_translation_layout.addWidget(self.ref_translation_label, 1)
        ref_translation_layout.addWidget(ref_translation_btn)
        file_layout.addRow("Translation Dynamics (BED):", ref_translation_layout)
        
        file_group.setLayout(file_layout)
        layout.addWidget(file_group)
        
        self.load_ref_btn = QPushButton("🚀 Load Reference Genome")
        self.load_ref_btn.setFont(QFont("Arial", 10, QFont.Bold))
        self.load_ref_btn.setMinimumHeight(35)
        self.load_ref_btn.clicked.connect(self.load_reference_genome)
        layout.addWidget(self.load_ref_btn)
        
        self.ref_progress_bar = QProgressBar()
        self.ref_progress_bar.setVisible(False)
        layout.addWidget(self.ref_progress_bar)
        
        self.ref_status_label = QLabel("")
        layout.addWidget(self.ref_status_label)
        
        self.ref_stats_text = QTextEdit()
        self.ref_stats_text.setReadOnly(True)
        self.ref_stats_text.setMaximumHeight(120)
        layout.addWidget(self.ref_stats_text)
        
        layout.addStretch()
        self.tabs.addTab(tab, "1. Reference Genome")

    def init_status_bar(self):
        """Initialize the status bar with GPU thermal indicators"""
        status_bar = self.statusBar()
        status_bar.setStyleSheet("QStatusBar { color: #555; background-color: #f8f8f8; }")
        
        # Add a permanent widget for GPU thermals
        thermal_widget = QWidget()
        thermal_layout = QHBoxLayout(thermal_widget)
        thermal_layout.setContentsMargins(5, 0, 15, 0)
        thermal_layout.setSpacing(15)
        
        for i in range(4):
            label = QLabel(f"GPU {i}: --°C")
            label.setFont(QFont("Monospace", 9))
            label.setStyleSheet("color: #666;")
            thermal_layout.addWidget(label)
            self.gpu_temp_labels[i] = label
            
        status_bar.addPermanentWidget(thermal_widget)
        status_bar.showMessage("Ready. Quad-GPU Orchestration Active.")

    def create_query_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        title = QLabel("🔬 Query Genome (JCVI-syn3.0)")
        title.setFont(QFont("Arial", 13, QFont.Bold))
        layout.addWidget(title)
        
        file_group = QGroupBox("Select Query Files")
        file_layout = QFormLayout()
        
        fasta_layout = QHBoxLayout()
        self.query_fasta_label = QLabel("No file selected")
        fasta_btn = QPushButton("Browse FASTA")
        fasta_btn.clicked.connect(self.select_query_fasta)
        fasta_layout.addWidget(self.query_fasta_label, 1)
        fasta_layout.addWidget(fasta_btn)
        file_layout.addRow("FASTA File:", fasta_layout)
        
        gff_layout = QHBoxLayout()
        self.query_gff_label = QLabel("No file selected (optional)")
        gff_btn = QPushButton("Browse GFF3")
        gff_btn.clicked.connect(self.select_query_gff)
        gff_layout.addWidget(self.query_gff_label, 1)
        gff_layout.addWidget(gff_btn)
        file_layout.addRow("GFF3 File (Optional):", gff_layout)
        
        likely_layout = QHBoxLayout()
        self.query_likely_label = QLabel("No file selected")
        likely_btn = QPushButton("Browse BED")
        likely_btn.clicked.connect(self.select_query_likely_promoters)
        likely_layout.addWidget(self.query_likely_label, 1)
        likely_layout.addWidget(likely_btn)
        file_layout.addRow("Likely Promoters (BED):", likely_layout)
        
        operon_layout = QHBoxLayout()
        self.query_operon_label = QLabel("No file selected")
        operon_btn = QPushButton("Browse BED")
        operon_btn.clicked.connect(self.select_query_internal_operons)
        operon_layout.addWidget(self.query_operon_label, 1)
        operon_layout.addWidget(operon_btn)
        file_layout.addRow("Internal Operons (BED):", operon_layout)
        
        # === New BED file inputs ===
        
        # Complex Cases BED
        complex_layout = QHBoxLayout()
        self.query_complex_label = QLabel("No file selected (optional)")
        complex_btn = QPushButton("Browse BED")
        complex_btn.clicked.connect(self.select_query_complex_cases)
        complex_layout.addWidget(self.query_complex_label, 1)
        complex_layout.addWidget(complex_btn)
        file_layout.addRow("Complex Cases (BED):", complex_layout)
        
        # Regulation Analysis BED
        regulation_layout = QHBoxLayout()
        self.query_regulation_label = QLabel("No file selected → Promoter Analysis")
        regulation_btn = QPushButton("Browse BED")
        regulation_btn.clicked.connect(self.select_query_regulation)
        regulation_layout.addWidget(self.query_regulation_label, 1)
        regulation_layout.addWidget(regulation_btn)
        file_layout.addRow("Regulation Analysis (BED):", regulation_layout)
        
        # Transcription Dynamics BED
        transcription_layout = QHBoxLayout()
        self.query_transcription_label = QLabel("No file selected → Transcription Animation")
        transcription_btn = QPushButton("Browse BED")
        transcription_btn.clicked.connect(self.select_query_transcription)
        transcription_layout.addWidget(self.query_transcription_label, 1)
        transcription_layout.addWidget(transcription_btn)
        file_layout.addRow("Transcription Dynamics (BED):", transcription_layout)
        
        # Translation Dynamics BED
        translation_layout = QHBoxLayout()
        self.query_translation_label = QLabel("No file selected → Translation Animation")
        translation_btn = QPushButton("Browse BED")
        translation_btn.clicked.connect(self.select_query_translation)
        translation_layout.addWidget(self.query_translation_label, 1)
        translation_layout.addWidget(translation_btn)
        file_layout.addRow("Translation Dynamics (BED):", translation_layout)
        
        file_group.setLayout(file_layout)
        layout.addWidget(file_group)
        
        self.load_query_btn = QPushButton("🚀 Load Query Genome")
        self.load_query_btn.setFont(QFont("Arial", 10, QFont.Bold))
        self.load_query_btn.setMinimumHeight(35)
        self.load_query_btn.clicked.connect(self.load_query_genome)
        layout.addWidget(self.load_query_btn)
        
        self.query_progress_bar = QProgressBar()
        self.query_progress_bar.setVisible(False)
        layout.addWidget(self.query_progress_bar)
        
        self.query_status_label = QLabel("")
        layout.addWidget(self.query_status_label)
        
        self.query_stats_text = QTextEdit()
        self.query_stats_text.setReadOnly(True)
        self.query_stats_text.setMaximumHeight(120)
        layout.addWidget(self.query_stats_text)
        
        layout.addStretch()
        self.tabs.addTab(tab, "2. Query Genome")
    
    def create_search_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        title = QLabel("🔍 BLAST Comparison")
        title.setFont(QFont("Arial", 13, QFont.Bold))
        layout.addWidget(title)
        
        info = QLabel("Query genome will be compared against reference genome")
        layout.addWidget(info)
        
        self.search_btn = QPushButton("🚀 Run BLAST Analysis")
        self.search_btn.setFont(QFont("Arial", 10, QFont.Bold))
        self.search_btn.setMinimumHeight(35)
        self.search_btn.clicked.connect(self.run_blast_search)
        layout.addWidget(self.search_btn)
        
        self.search_progress_bar = QProgressBar()
        self.search_progress_bar.setVisible(False)
        layout.addWidget(self.search_progress_bar)
        
        self.search_status_label = QLabel("")
        layout.addWidget(self.search_status_label)
        
        layout.addStretch()
        self.tabs.addTab(tab, "3. BLAST Search")
    
    def create_results_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        title = QLabel("📊 BLAST Results")
        title.setFont(QFont("Arial", 13, QFont.Bold))
        layout.addWidget(title)
        
        self.summary_label = QLabel("")
        self.summary_label.setFont(QFont("Arial", 10))
        layout.addWidget(self.summary_label)
        
        results_tabs = QTabWidget()
        
        hits_tab = QWidget()
        hits_layout = QVBoxLayout(hits_tab)
        
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(7)
        self.results_table.setHorizontalHeaderLabels([
            "Identity %", "Position (Query | Ref)", "E-value", 
            "Query Genes (syn3.0)", "Function", "Protein", "Details"
        ])
        self.results_table.horizontalHeader().setStretchLastSection(True)
        self.results_table.cellClicked.connect(self.show_result_details)
        hits_layout.addWidget(self.results_table)
        
        self.details_text = QTextEdit()
        self.details_text.setReadOnly(True)
        self.details_text.setMaximumHeight(120)
        hits_layout.addWidget(self.details_text)
        
        results_tabs.addTab(hits_tab, "All Hits")
        
        gene_summary_tab = QWidget()
        gene_summary_layout = QVBoxLayout(gene_summary_tab)
        
        gene_info = QLabel("📋 Gene-level Summary (Query Genome)")
        gene_info.setFont(QFont("Arial", 11, QFont.Bold))
        gene_summary_layout.addWidget(gene_info)
        
        self.gene_summary_table = QTableWidget()
        self.gene_summary_table.setColumnCount(6)
        self.gene_summary_table.setHorizontalHeaderLabels([
            "Gene ID", "Gene Name", "Function", "Best Identity %", "Coverage", "Status"
        ])
        self.gene_summary_table.horizontalHeader().setStretchLastSection(True)
        gene_summary_layout.addWidget(self.gene_summary_table)
        
        results_tabs.addTab(gene_summary_tab, "Gene Summary")
        
        layout.addWidget(results_tabs)
        
        export_btn = QPushButton("📥 Export Results to CSV")
        export_btn.clicked.connect(self.export_results)
        layout.addWidget(export_btn)
        
        ai_btn = QPushButton("🤖 Run AI Analysis")
        ai_btn.setFont(QFont("Arial", 10, QFont.Bold))
        ai_btn.setMinimumHeight(35)
        ai_btn.clicked.connect(self.run_ai_analysis)
        layout.addWidget(ai_btn)
        
        self.tabs.addTab(tab, "4. BLAST Results")
    
    def create_ai_results_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        title = QLabel("🤖 AI Analysis Results")
        title.setFont(QFont("Arial", 13, QFont.Bold))
        layout.addWidget(title)
        
        ai_tabs = QTabWidget()
        
        # Promoter analysis tab
        promoter_tab = QWidget()
        promoter_layout = QVBoxLayout(promoter_tab)
        
        promoter_title = QLabel("Agent #1: Promoter Analysis - QUERY GENOME (syn3.0)")
        promoter_title.setFont(QFont("Arial", 12, QFont.Bold))
        promoter_title.setStyleSheet("color: #1565C0;")
        promoter_layout.addWidget(promoter_title)
        
        self.promoter_table = QTableWidget()
        self.promoter_table.setColumnCount(8)
        self.promoter_table.setHorizontalHeaderLabels([
            "Gene ID", "Regulation Type", "Category", "Function", 
            "Strength", "TX Level", "Protein", "Sources"
        ])
        promoter_layout.addWidget(self.promoter_table)
        
        ai_tabs.addTab(promoter_tab, "Promoter Analysis (Query)")  # Will be pushed to index 1 after insertTab(0,...)
        
        # Transcription animation tab (full width) — wrapped in scroll area
        transcription_scroll = QScrollArea()
        transcription_scroll.setWidgetResizable(True)
        transcription_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        transcription_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        
        transcription_tab = QWidget()
        transcription_tab.setMinimumHeight(1000) # Ensure it's large enough to scroll
        transcription_layout = QVBoxLayout(transcription_tab)
        
        trans_title = QLabel("Agent #2: Live Transcription Animation 🧬 - QUERY GENOME (syn3.0)")
        trans_title.setFont(QFont("Arial", 12, QFont.Bold))
        trans_title.setStyleSheet("color: #2E7D32;")
        transcription_layout.addWidget(trans_title)
        
        # MOVED TO TOP: Button to open Resource Dashboard dialog
        if RESOURCE_DASHBOARD_AVAILABLE and self.resource_dashboard:
            open_resource_btn = QPushButton("⚡ Open Resource Competition Dashboard")
            open_resource_btn.setFont(QFont("Arial", 9, QFont.Bold))
            open_resource_btn.setMinimumHeight(30)
            open_resource_btn.setStyleSheet("background-color: #FF8F00; color: white;")
            open_resource_btn.clicked.connect(lambda: self.resource_dashboard.show())
            transcription_layout.addWidget(open_resource_btn)
        
        # Button to open translation window
        open_translation_btn = QPushButton("🔬 Open Translation Window (Separate)")
        open_translation_btn.setFont(QFont("Arial", 9, QFont.Bold))
        open_translation_btn.setMinimumHeight(30)
        open_translation_btn.setStyleSheet("background-color: #FF5722; color: white;")
        open_translation_btn.clicked.connect(lambda: self.translation_window.show())
        transcription_layout.addWidget(open_translation_btn)
        
        # Button to open Live Data Viewer
        open_data_viewer_btn = QPushButton("📊 Open Live Data Viewer (Real-Time)")
        open_data_viewer_btn.setFont(QFont("Arial", 9, QFont.Bold))
        open_data_viewer_btn.setMinimumHeight(30)
        open_data_viewer_btn.setStyleSheet("background-color: #2196F3; color: white;")
        open_data_viewer_btn.clicked.connect(lambda: self.live_data_viewer.show())
        transcription_layout.addWidget(open_data_viewer_btn)
        
        # Button to open Metabolic Dashboard
        if METABOLIC_DASHBOARD_AVAILABLE:
            open_metabolic_btn = QPushButton("🧬 Open Metabolic Dashboard (Real-Time FBA)")
            open_metabolic_btn.setFont(QFont("Arial", 9, QFont.Bold))
            open_metabolic_btn.setMinimumHeight(30)
            open_metabolic_btn.setStyleSheet("background-color: #00897B; color: white;")
            open_metabolic_btn.clicked.connect(lambda: self.metabolic_dashboard.show())
            transcription_layout.addWidget(open_metabolic_btn)
        
        # Buttons to open individual Gene Dynamics dialogs
        if GENE_DYNAMICS_AVAILABLE:
            dynamics_btn_layout = QHBoxLayout()

            btn_mrna_tx = QPushButton("📝 mRNA Transcription")
            btn_mrna_tx.setFont(QFont("Arial", 8, QFont.Bold))
            btn_mrna_tx.setMinimumHeight(28)
            btn_mrna_tx.setStyleSheet("background-color: #1565C0; color: white;")
            btn_mrna_tx.clicked.connect(lambda: self._open_gene_dynamics_panel('mrna_tx'))
            dynamics_btn_layout.addWidget(btn_mrna_tx)

            btn_protein_tl = QPushButton("🧬 Protein Translation")
            btn_protein_tl.setFont(QFont("Arial", 8, QFont.Bold))
            btn_protein_tl.setMinimumHeight(28)
            btn_protein_tl.setStyleSheet("background-color: #2E7D32; color: white;")
            btn_protein_tl.clicked.connect(lambda: self._open_gene_dynamics_panel('protein_tl'))
            dynamics_btn_layout.addWidget(btn_protein_tl)

            btn_mrna_deg = QPushButton("💀 mRNA Degradation")
            btn_mrna_deg.setFont(QFont("Arial", 8, QFont.Bold))
            btn_mrna_deg.setMinimumHeight(28)
            btn_mrna_deg.setStyleSheet("background-color: #E65100; color: white;")
            btn_mrna_deg.clicked.connect(lambda: self._open_gene_dynamics_panel('mrna_deg'))
            dynamics_btn_layout.addWidget(btn_mrna_deg)

            btn_prot_deg = QPushButton("🛡️ Protein Stability")
            btn_prot_deg.setFont(QFont("Arial", 8, QFont.Bold))
            btn_prot_deg.setMinimumHeight(28)
            btn_prot_deg.setStyleSheet("background-color: #6A1B9A; color: white;")
            btn_prot_deg.clicked.connect(lambda: self._open_gene_dynamics_panel('protein_deg'))
            dynamics_btn_layout.addWidget(btn_prot_deg)

            transcription_layout.addLayout(dynamics_btn_layout)
        
        control_layout = QHBoxLayout()
        
        self.start_animation_btn = QPushButton("▶️ Start Transcription")
        self.start_animation_btn.setFont(QFont("Arial", 9, QFont.Bold))
        self.start_animation_btn.setMinimumHeight(30)
        self.start_animation_btn.clicked.connect(self.start_transcription_animation)
        self.start_animation_btn.setStyleSheet("background-color: #4CAF50; color: white;")
        control_layout.addWidget(self.start_animation_btn)
        
        self.stop_animation_btn = QPushButton("⏸️ Pause All")
        self.stop_animation_btn.setEnabled(False)
        self.stop_animation_btn.clicked.connect(self.stop_transcription_animation)
        control_layout.addWidget(self.stop_animation_btn)
        
        self.reset_animation_btn = QPushButton("🔄 Reset")
        self.reset_animation_btn.clicked.connect(self.reset_transcription_animation)
        control_layout.addWidget(self.reset_animation_btn)
        
        transcription_layout.addLayout(control_layout)
        
        gene_limit_layout = QHBoxLayout()
        self.gene_limit_spinbox = QSpinBox()
        self.gene_limit_spinbox.setRange(10, 1000)
        self.gene_limit_spinbox.setValue(473)
        self.gene_limit_spinbox.setSuffix(" genes")
        self.gene_limit_spinbox.valueChanged.connect(self.on_gene_limit_changed)
        gene_limit_layout.addWidget(QLabel("Number of genes:"))
        gene_limit_layout.addWidget(self.gene_limit_spinbox)
        gene_limit_layout.addStretch()
        transcription_layout.addLayout(gene_limit_layout)
        
        time_layout = QHBoxLayout()
        self.time_label = QLabel("⏱️ Time: 0.0s")
        self.time_label.setFont(QFont("Arial", 11, QFont.Bold))
        time_layout.addWidget(self.time_label)
        
        self.speed_spinbox = QDoubleSpinBox()
        self.speed_spinbox.setRange(0.1, 10.0)
        self.speed_spinbox.setValue(1.0)
        self.speed_spinbox.setSuffix("x")
        time_layout.addWidget(QLabel("Speed:"))
        time_layout.addWidget(self.speed_spinbox)
        time_layout.addStretch()
        transcription_layout.addLayout(time_layout)
        
        # === CONTINUOUS MODE CONTROLS ===
        from PyQt5.QtWidgets import QCheckBox
        
        continuous_group = QGroupBox("🔄 Continuous Cell Simulation Mode")
        continuous_group.setStyleSheet("QGroupBox { font-weight: bold; color: #1565C0; }")
        continuous_layout = QHBoxLayout()
        
        # Mode toggle
        self.continuous_mode_checkbox = QCheckBox("Enable Continuous Mode")
        self.continuous_mode_checkbox.setChecked(True)  # Default to continuous
        self.continuous_mode_checkbox.setToolTip(
            "Linear: Genes reach 100% and stop\n"
            "Continuous: Genes fluctuate with synthesis + degradation"
        )
        continuous_layout.addWidget(self.continuous_mode_checkbox)
        
        # Quantum Simulation Mode
        self.enable_quantum_checkbox = QCheckBox("Enable Quantum Simulation (10 Genes Max)")
        self.enable_quantum_checkbox.setChecked(True)
        self.enable_quantum_checkbox.setToolTip("Submit transcription rates to OQTOPUS Cloud for probability calculation")
        self.enable_quantum_checkbox.setStyleSheet("color: #673AB7; font-weight: bold;")
        continuous_layout.addWidget(self.enable_quantum_checkbox)
        
        # mRNA half-life
        continuous_layout.addWidget(QLabel("mRNA Half-Life:"))
        self.mrna_halflife_spinbox = QSpinBox()
        self.mrna_halflife_spinbox.setRange(30, 1800)  # 30 sec to 30 min
        self.mrna_halflife_spinbox.setValue(300)  # 5 min default
        self.mrna_halflife_spinbox.setSuffix(" sec")
        self.mrna_halflife_spinbox.setToolTip("Time for mRNA to degrade by 50% (typical: 2-10 min)")
        continuous_layout.addWidget(self.mrna_halflife_spinbox)
        
        # Protein half-life
        continuous_layout.addWidget(QLabel("Protein Half-Life:"))
        self.protein_halflife_spinbox = QSpinBox()
        self.protein_halflife_spinbox.setRange(300, 36000)  # 5 min to 10 hours
        self.protein_halflife_spinbox.setValue(3600)  # 1 hr default
        self.protein_halflife_spinbox.setSuffix(" sec")
        self.protein_halflife_spinbox.setToolTip("Time for protein to degrade by 50% (typical: 30min - 24hr)")
        continuous_layout.addWidget(self.protein_halflife_spinbox)
        
        # Steady-state indicator
        self.steady_state_label = QLabel("🟡 Initializing...")
        self.steady_state_label.setStyleSheet("font-weight: bold; padding: 5px;")
        continuous_layout.addWidget(self.steady_state_label)
        
        continuous_layout.addStretch()
        continuous_group.setLayout(continuous_layout)
        transcription_layout.addWidget(continuous_group)
        
        self.transcription_table = QTableWidget()
        self.transcription_table.setColumnCount(12)  # Added Quantum P column
        self.transcription_table.setHorizontalHeaderLabels([
            "Gene ID", "Level", "Status", "Current", "Phase", "Category", 
            "Regulation", "Order", "Control", "Status Icon", "Paused", "Quantum P"
        ])
        # Set column widths to make Control buttons visible
        self.transcription_table.setColumnWidth(0, 85)   # Gene ID
        self.transcription_table.setColumnWidth(1, 40)   # Level
        self.transcription_table.setColumnWidth(2, 70)   # Status
        self.transcription_table.setColumnWidth(3, 50)   # Current
        self.transcription_table.setColumnWidth(4, 105)  # Phase
        self.transcription_table.setColumnWidth(5, 70)   # Category
        self.transcription_table.setColumnWidth(6, 70)   # Regulation
        self.transcription_table.setColumnWidth(7, 35)   # Order
        self.transcription_table.setColumnWidth(8, 65)   # Control (BUTTONS HERE!)
        self.transcription_table.setColumnWidth(9, 85)   # Status Icon
        self.transcription_table.setColumnWidth(10, 50)  # Paused
        self.transcription_table.horizontalHeader().setStretchLastSection(True)
        transcription_layout.addWidget(self.transcription_table)
        # Removed stretch to allow table to fill the space

        
        # Connect table click to gene dynamics panel
        self.transcription_table.cellClicked.connect(self._on_transcription_table_clicked)
        
        # Button to open Resource Dashboard dialog
        transcription_scroll.setWidget(transcription_tab)
        
        # Live Transcription tab added FIRST for top visibility
        ai_tabs.insertTab(0, transcription_scroll, "Live Transcription (Query)")
        ai_tabs.setCurrentIndex(0)
        
        # Comparison tab
        comparison_tab = QWidget()
        comparison_layout = QVBoxLayout(comparison_tab)
        self.comparison_text = QTextEdit()
        self.comparison_text.setReadOnly(True)
        comparison_layout.addWidget(self.comparison_text)
        ai_tabs.addTab(comparison_tab, "Comparison")
        
        # Integrated Protein Synthesis tab
        ai_tabs.addTab(self.protein_synthesis_widget, "🧬 Protein Synthesis")
        
        layout.addWidget(ai_tabs)
        
        export_ai_btn = QPushButton("📥 Export AI Results to Excel")
        export_ai_btn.clicked.connect(self.export_ai_results)
        layout.addWidget(export_ai_btn)
        
        self.tabs.addTab(tab, "5. AI Results")
    
    def _open_gene_dynamics_panel(self, dialog_name: str = 'mrna_tx'):
        """Open a specific Gene Dynamics dialog by name."""
        if not hasattr(self, 'gene_dynamics_panel') or not self.gene_dynamics_panel:
            return
        
        # Try to determine current selected row from the active table
        current_row = -1
        if self.tabs.currentIndex() == 0: # Transcription tab
            current_row = self.transcription_table.currentRow()
        elif self.tabs.currentIndex() == 1: # AI Results? (Actually we want Translation window)
            current_row = self.translation_table.currentRow()
            
        panel = self.gene_dynamics_panel
        if current_row >= 0:
            panel.set_gene(current_row)
            
        dialog_map = {
            'mrna_tx': panel.mrna_tx,
            'protein_tl': panel.protein_tl,
            'mrna_deg': panel.mrna_deg,
            'protein_deg': panel.protein_deg,
        }
        dlg = dialog_map.get(dialog_name)
        if dlg:
            if self.gene_transcription_data:
                dlg.populate_genes(self.gene_transcription_data)
            dlg.show()
            dlg.raise_()
    
    def _on_transcription_table_clicked(self, row, col):
        """When user clicks a gene row, open mRNA Transcription dialog for that gene."""
        if hasattr(self, 'gene_dynamics_panel') and self.gene_dynamics_panel:
            panel = self.gene_dynamics_panel
            if self.gene_transcription_data:
                panel.populate_genes(self.gene_transcription_data)
            panel.set_gene(row)
            panel.mrna_tx.show()
            panel.mrna_tx.raise_()

    def _on_translation_table_clicked(self, row, col):
        """When user clicks a gene row in translation window, open Protein Translation dialog."""
        if hasattr(self, 'gene_dynamics_panel') and self.gene_dynamics_panel:
            panel = self.gene_dynamics_panel
            if self.gene_transcription_data:
                panel.populate_genes(self.gene_transcription_data)
            panel.set_gene(row)
            panel.protein_tl.show()
            panel.protein_tl.raise_()
    
    # File selection methods
    def select_ref_fasta(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Reference FASTA", "", "FASTA Files (*.fasta *.fa *.fna)")
        if path:
            self.fasta_path = path
            self.ref_fasta_label.setText(os.path.basename(path))
    
    def select_ref_gff(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Reference GFF3", "", "GFF Files (*.gff *.gff3)")
        if path:
            self.gff_path = path
            self.ref_gff_label.setText(os.path.basename(path))
    
    def select_ref_likely_promoters(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Likely Promoters BED", "", "BED Files (*.bed)")
        if path:
            self.likely_promoters_path = path
            self.ref_likely_label.setText(os.path.basename(path))
    
    def select_ref_internal_operons(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Internal Operons BED", "", "BED Files (*.bed)")
        if path:
            self.internal_operons_path = path
            self.ref_operon_label.setText(os.path.basename(path))
    
    # === Reference: Additional BED file selection methods ===
    def select_ref_complex_cases(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Reference Complex Cases BED", "", "BED Files (*.bed)")
        if path:
            self.ref_complex_cases_path = path
            self.ref_complex_label.setText(os.path.basename(path))
    
    def select_ref_regulation(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Reference Regulation Analysis BED", "", "BED Files (*.bed)")
        if path:
            self.ref_regulation_path = path
            self.ref_regulation_label.setText(f"✅ {os.path.basename(path)}")
    
    def select_ref_transcription(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Reference Transcription Dynamics BED", "", "BED Files (*.bed)")
        if path:
            self.ref_transcription_path = path
            self.ref_transcription_label.setText(f"✅ {os.path.basename(path)}")
    
    def select_ref_translation(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Reference Translation Dynamics BED", "", "BED Files (*.bed)")
        if path:
            self.ref_translation_path = path
            self.ref_translation_label.setText(f"✅ {os.path.basename(path)}")
    
    def select_query_fasta(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Query FASTA", "", "FASTA Files (*.fasta *.fa *.fna)")
        if path:
            self.query_fasta_path = path
            self.query_fasta_label.setText(os.path.basename(path))
    
    def select_query_gff(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Query GFF3", "", "GFF Files (*.gff *.gff3)")
        if path:
            self.query_gff_path = path
            self.query_gff_label.setText(os.path.basename(path))
    
    def select_query_likely_promoters(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Query Likely Promoters BED", "", "BED Files (*.bed)")
        if path:
            self.query_likely_promoters_path = path
            self.query_likely_label.setText(os.path.basename(path))
    
    def select_query_internal_operons(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Query Internal Operons BED", "", "BED Files (*.bed)")
        if path:
            self.query_internal_operons_path = path
            self.query_operon_label.setText(os.path.basename(path))
    
    # === NEW: Additional file selection methods ===
    def select_query_complex_cases(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Complex Cases BED", "", "BED Files (*.bed)")
        if path:
            self.query_complex_cases_path = path
            self.query_complex_label.setText(os.path.basename(path))
    
    def select_query_regulation(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Regulation Analysis BED", "", "BED Files (*.bed)")
        if path:
            self.query_regulation_path = path
            self.query_regulation_label.setText(f"✅ {os.path.basename(path)}")
    
    def select_query_transcription(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Transcription Dynamics BED", "", "BED Files (*.bed)")
        if path:
            self.query_transcription_path = path
            self.query_transcription_label.setText(f"✅ {os.path.basename(path)}")
    
    def select_query_translation(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Translation Dynamics BED", "", "BED Files (*.bed)")
        if path:
            self.query_translation_path = path
            self.query_translation_label.setText(f"✅ {os.path.basename(path)}")
