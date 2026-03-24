#!/usr/bin/env python3
"""
Amino Acid Generator Application
Receives translation data via ZeroMQ and displays protein sequences in real-time

Features:
- ZeroMQ subscriber (connects to JCVI Genome Analyzer)
- Real-time CDS to Amino Acid conversion
- Progressive protein display based on translation %
- Active and Completed sections
- Three-letter amino acid format
- Secondary Structure Prediction (Chou-Fasman method)
"""

import sys
import json
from typing import Dict, Optional, List, Tuple

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QTableWidget, QTableWidgetItem, QPushButton, QGroupBox,
    QHeaderView, QSplitter, QFrame, QMessageBox, QTabWidget,
    QTextEdit, QScrollArea
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt5.QtGui import QFont, QColor, QPalette

try:
    import zmq
    ZMQ_AVAILABLE = True
except ImportError:
    ZMQ_AVAILABLE = False
    print("WARNING: PyZMQ not installed. Install with: pip install pyzmq")


# =============================================================================
# Genetic Code - Codon to Amino Acid mapping
# =============================================================================

CODON_TABLE = {
    # Phenylalanine
    'TTT': 'Phe', 'TTC': 'Phe',
    # Leucine
    'TTA': 'Leu', 'TTG': 'Leu', 'CTT': 'Leu', 'CTC': 'Leu', 'CTA': 'Leu', 'CTG': 'Leu',
    # Isoleucine
    'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile',
    # Methionine (Start)
    'ATG': 'Met',
    # Valine
    'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val',
    # Serine
    'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser', 'AGT': 'Ser', 'AGC': 'Ser',
    # Proline
    'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    # Threonine
    'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    # Alanine
    'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    # Tyrosine
    'TAT': 'Tyr', 'TAC': 'Tyr',
    # Stop codons
    'TAA': 'Stop', 'TAG': 'Stop', 'TGA': 'Stop',
    # Histidine
    'CAT': 'His', 'CAC': 'His',
    # Glutamine
    'CAA': 'Gln', 'CAG': 'Gln',
    # Asparagine
    'AAT': 'Asn', 'AAC': 'Asn',
    # Lysine
    'AAA': 'Lys', 'AAG': 'Lys',
    # Aspartic Acid
    'GAT': 'Asp', 'GAC': 'Asp',
    # Glutamic Acid
    'GAA': 'Glu', 'GAG': 'Glu',
    # Cysteine
    'TGT': 'Cys', 'TGC': 'Cys',
    # Tryptophan
    'TGG': 'Trp',
    # Arginine
    'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'AGA': 'Arg', 'AGG': 'Arg',
    # Glycine
    'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
}

# One-letter to Three-letter mapping
ONE_TO_THREE = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
}

THREE_TO_ONE = {v: k for k, v in ONE_TO_THREE.items()}


# =============================================================================
# Chou-Fasman Parameters for Secondary Structure Prediction
# =============================================================================

# Propensity values for Alpha Helix (Pα)
HELIX_PROPENSITY = {
    'Ala': 1.42, 'Arg': 0.98, 'Asn': 0.67, 'Asp': 1.01, 'Cys': 0.70,
    'Gln': 1.11, 'Glu': 1.51, 'Gly': 0.57, 'His': 1.00, 'Ile': 1.08,
    'Leu': 1.21, 'Lys': 1.16, 'Met': 1.45, 'Phe': 1.13, 'Pro': 0.57,
    'Ser': 0.77, 'Thr': 0.83, 'Trp': 1.08, 'Tyr': 0.69, 'Val': 1.06
}

# Propensity values for Beta Sheet (Pβ)
SHEET_PROPENSITY = {
    'Ala': 0.83, 'Arg': 0.93, 'Asn': 0.89, 'Asp': 0.54, 'Cys': 1.19,
    'Gln': 1.10, 'Glu': 0.37, 'Gly': 0.75, 'His': 0.87, 'Ile': 1.60,
    'Leu': 1.30, 'Lys': 0.74, 'Met': 1.05, 'Phe': 1.38, 'Pro': 0.55,
    'Ser': 0.75, 'Thr': 1.19, 'Trp': 1.37, 'Tyr': 1.47, 'Val': 1.70
}

# Propensity values for Turn (Pt)
TURN_PROPENSITY = {
    'Ala': 0.66, 'Arg': 0.95, 'Asn': 1.56, 'Asp': 1.46, 'Cys': 1.19,
    'Gln': 0.98, 'Glu': 0.74, 'Gly': 1.56, 'His': 0.95, 'Ile': 0.47,
    'Leu': 0.59, 'Lys': 1.01, 'Met': 0.60, 'Phe': 0.60, 'Pro': 1.52,
    'Ser': 1.43, 'Thr': 0.96, 'Trp': 0.96, 'Tyr': 1.14, 'Val': 0.50
}

# Thresholds for structure assignment
HELIX_THRESHOLD = 1.03
SHEET_THRESHOLD = 1.05
TURN_THRESHOLD = 1.00


def predict_secondary_structure(amino_acids: List[str]) -> Tuple[str, Dict]:
    """
    Predict secondary structure using Chou-Fasman method
    
    Args:
        amino_acids: List of three-letter amino acid codes
    
    Returns:
        Tuple of (structure_string, statistics_dict)
        structure_string: H=Helix, E=Sheet, T=Turn, C=Coil
    """
    if not amino_acids:
        return "", {}
    
    n = len(amino_acids)
    structure = ['C'] * n  # Default to coil
    
    # Calculate propensities for each position
    helix_scores = []
    sheet_scores = []
    turn_scores = []
    
    for aa in amino_acids:
        helix_scores.append(HELIX_PROPENSITY.get(aa, 1.0))
        sheet_scores.append(SHEET_PROPENSITY.get(aa, 1.0))
        turn_scores.append(TURN_PROPENSITY.get(aa, 1.0))
    
    # Identify helix nucleation sites (4 out of 6 residues with high helix propensity)
    window_size = 6
    for i in range(n - window_size + 1):
        window = helix_scores[i:i + window_size]
        high_propensity_count = sum(1 for p in window if p > HELIX_THRESHOLD)
        if high_propensity_count >= 4:
            # Extend helix in both directions
            for j in range(i, min(i + window_size, n)):
                if helix_scores[j] > 0.9:  # Allow some flexibility
                    structure[j] = 'H'
    
    # Identify sheet nucleation sites (3 out of 5 residues with high sheet propensity)
    window_size = 5
    for i in range(n - window_size + 1):
        window = sheet_scores[i:i + window_size]
        high_propensity_count = sum(1 for p in window if p > SHEET_THRESHOLD)
        if high_propensity_count >= 3:
            for j in range(i, min(i + window_size, n)):
                if structure[j] == 'C' and sheet_scores[j] > 0.9:
                    structure[j] = 'E'
    
    # Identify turns (4-residue segments)
    for i in range(n - 3):
        turn_score = sum(turn_scores[i:i + 4]) / 4
        if turn_score > TURN_THRESHOLD:
            # Check if not already assigned to helix or sheet
            if all(structure[j] == 'C' for j in range(i, i + 4)):
                for j in range(i, i + 4):
                    structure[j] = 'T'
    
    structure_string = ''.join(structure)
    
    # Calculate statistics
    stats = {
        'helix_count': structure_string.count('H'),
        'sheet_count': structure_string.count('E'),
        'turn_count': structure_string.count('T'),
        'coil_count': structure_string.count('C'),
        'total': n,
        'helix_percent': (structure_string.count('H') / n * 100) if n > 0 else 0,
        'sheet_percent': (structure_string.count('E') / n * 100) if n > 0 else 0,
        'turn_percent': (structure_string.count('T') / n * 100) if n > 0 else 0,
        'coil_percent': (structure_string.count('C') / n * 100) if n > 0 else 0,
    }
    
    return structure_string, stats


def translate_cds_to_amino_acids(cds_sequence: str, percentage: float = 100.0) -> str:
    """
    Translate CDS sequence to amino acid sequence (three-letter format)
    
    Args:
        cds_sequence: DNA sequence (must be divisible by 3)
        percentage: Translation percentage (0-100), determines how much to translate
    
    Returns:
        Amino acid sequence in three-letter format (e.g., "Met-Lys-Arg-Phe")
    """
    if not cds_sequence:
        return ""
    
    # Clean sequence
    cds_sequence = cds_sequence.upper().replace(' ', '').replace('\n', '')
    
    # Calculate how many codons to translate based on percentage
    total_codons = len(cds_sequence) // 3
    codons_to_translate = int(total_codons * (percentage / 100.0))
    
    if codons_to_translate == 0:
        return ""
    
    # Translate codons
    amino_acids = []
    for i in range(codons_to_translate):
        codon = cds_sequence[i*3:(i+1)*3]
        if len(codon) == 3:
            aa = CODON_TABLE.get(codon, 'Unk')
            if aa == 'Stop':
                break
            amino_acids.append(aa)
    
    return '-'.join(amino_acids)


def get_amino_acid_list(aa_string: str) -> List[str]:
    """Convert amino acid string to list"""
    if not aa_string:
        return []
    return aa_string.split('-')


# =============================================================================
# ZeroMQ Receiver Thread
# =============================================================================

class ZMQReceiverThread(QThread):
    """Background thread to receive ZeroMQ messages"""
    
    data_received = pyqtSignal(dict)
    connection_status = pyqtSignal(str)
    error_occurred = pyqtSignal(str)
    
    def __init__(self, port: int = 5555):
        super().__init__()
        self.port = port
        self.running = False
        self.context = None
        self.socket = None
    
    def run(self):
        """Main thread loop - receive messages"""
        if not ZMQ_AVAILABLE:
            self.error_occurred.emit("PyZMQ not installed!")
            return
        
        try:
            self.context = zmq.Context()
            self.socket = self.context.socket(zmq.SUB)
            self.socket.connect(f"tcp://localhost:{self.port}")
            self.socket.setsockopt_string(zmq.SUBSCRIBE, "")
            self.socket.setsockopt(zmq.RCVTIMEO, 1000)  # 1 second timeout
            
            self.running = True
            self.connection_status.emit(f"🟢 Connected to port {self.port}")
            
            while self.running:
                try:
                    message = self.socket.recv_string()
                    data = json.loads(message)
                    self.data_received.emit(data)
                except zmq.Again:
                    # Timeout - no message received, continue loop
                    continue
                except json.JSONDecodeError as e:
                    self.error_occurred.emit(f"JSON Error: {e}")
                except Exception as e:
                    if self.running:
                        self.error_occurred.emit(f"Receive Error: {e}")
        
        except Exception as e:
            self.error_occurred.emit(f"Connection Error: {e}")
        
        finally:
            self.cleanup()
    
    def stop(self):
        """Stop the receiver thread"""
        self.running = False
        self.connection_status.emit("🔴 Disconnected")
    
    def cleanup(self):
        """Clean up ZeroMQ resources"""
        if self.socket:
            self.socket.close()
        if self.context:
            self.context.term()


# =============================================================================
# Secondary Structure Visualization Widget
# =============================================================================

class SecondaryStructureWidget(QWidget):
    """Widget to display secondary structure prediction"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout(self)
        
        # Gene selector info
        self.gene_label = QLabel("Select a gene to view secondary structure")
        self.gene_label.setFont(QFont("Arial", 12, QFont.Bold))
        self.gene_label.setStyleSheet("color: #2c3e50; padding: 10px;")
        layout.addWidget(self.gene_label)
        
        # Structure visualization
        self.structure_display = QTextEdit()
        self.structure_display.setReadOnly(True)
        self.structure_display.setFont(QFont("Consolas", 10))
        self.structure_display.setMinimumHeight(150)
        layout.addWidget(self.structure_display)
        
        # Legend
        legend_layout = QHBoxLayout()
        
        legend_items = [
            ("H - α Helix", "#e74c3c"),
            ("E - β Sheet", "#3498db"),
            ("T - Turn", "#2ecc71"),
            ("C - Coil", "#95a5a6")
        ]
        
        for text, color in legend_items:
            label = QLabel(f"  ■ {text}  ")
            label.setStyleSheet(f"color: {color}; font-weight: bold; font-size: 11px;")
            legend_layout.addWidget(label)
        
        legend_layout.addStretch()
        layout.addLayout(legend_layout)
        
        # Statistics
        self.stats_label = QLabel("")
        self.stats_label.setStyleSheet("padding: 10px; background-color: #ecf0f1; border-radius: 5px;")
        layout.addWidget(self.stats_label)
    
    def update_display(self, gene_id: str, amino_acids: List[str], structure: str, stats: Dict):
        """Update the display with new data"""
        self.gene_label.setText(f"🧬 Gene: {gene_id}")
        
        # Format the display with colors
        html_parts = ["<pre style='font-family: Consolas; font-size: 11px;'>"]
        
        # Show amino acids (abbreviated to single letter for display)
        html_parts.append("<b>Sequence:</b>\n")
        
        # Split into lines of 60 characters
        line_length = 60
        
        for i in range(0, len(amino_acids), line_length):
            # Position numbers
            html_parts.append(f"<span style='color: #888;'>{i+1:5d}</span>  ")
            
            # Amino acids for this line
            line_aa = amino_acids[i:i + line_length]
            line_struct = structure[i:i + line_length] if i < len(structure) else ""
            
            for j, aa in enumerate(line_aa):
                one_letter = THREE_TO_ONE.get(aa, '?')
                struct_char = line_struct[j] if j < len(line_struct) else 'C'
                
                # Color based on structure
                if struct_char == 'H':
                    color = "#e74c3c"  # Red for helix
                elif struct_char == 'E':
                    color = "#3498db"  # Blue for sheet
                elif struct_char == 'T':
                    color = "#2ecc71"  # Green for turn
                else:
                    color = "#666666"  # Gray for coil
                
                html_parts.append(f"<span style='color: {color};'>{one_letter}</span>")
            
            html_parts.append("\n")
            
            # Structure line
            html_parts.append("        ")
            for char in line_struct:
                if char == 'H':
                    color = "#e74c3c"
                elif char == 'E':
                    color = "#3498db"
                elif char == 'T':
                    color = "#2ecc71"
                else:
                    color = "#666666"
                html_parts.append(f"<span style='color: {color};'>{char}</span>")
            html_parts.append("\n\n")
        
        html_parts.append("</pre>")
        
        self.structure_display.setHtml(''.join(html_parts))
        
        # Update statistics
        stats_text = (
            f"📊 <b>Secondary Structure Statistics:</b><br>"
            f"&nbsp;&nbsp;• α Helix: {stats.get('helix_count', 0)} residues ({stats.get('helix_percent', 0):.1f}%)<br>"
            f"&nbsp;&nbsp;• β Sheet: {stats.get('sheet_count', 0)} residues ({stats.get('sheet_percent', 0):.1f}%)<br>"
            f"&nbsp;&nbsp;• Turn: {stats.get('turn_count', 0)} residues ({stats.get('turn_percent', 0):.1f}%)<br>"
            f"&nbsp;&nbsp;• Coil: {stats.get('coil_count', 0)} residues ({stats.get('coil_percent', 0):.1f}%)<br>"
            f"&nbsp;&nbsp;• <b>Total: {stats.get('total', 0)} amino acids</b>"
        )
        self.stats_label.setText(stats_text)


# =============================================================================
# Main Application Window
# =============================================================================

class AminoAcidGeneratorApp(QMainWindow):
    """Main application window for Amino Acid Generator"""
    
    def __init__(self):
        super().__init__()
        
        # Data storage
        self.active_genes: Dict[str, dict] = {}
        self.completed_genes: Dict[str, dict] = {}
        self.message_count = 0
        
        # ZeroMQ receiver
        self.zmq_thread: Optional[ZMQReceiverThread] = None
        
        # Initialize UI
        self.init_ui()
        
        # Auto-connect on startup
        QTimer.singleShot(500, self.connect_zmq)
    
    def init_ui(self):
        """Initialize the user interface"""
        self.setWindowTitle("🧬 Amino Acid Generator - Live Protein Synthesis & Structure Prediction")
        self.setGeometry(100, 100, 1400, 900)
        self.setStyleSheet("""
            QMainWindow {
                background-color: #f5f5f5;
            }
            QGroupBox {
                font-weight: bold;
                border: 2px solid #3498db;
                border-radius: 8px;
                margin-top: 10px;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px;
            }
            QTableWidget {
                background-color: white;
                gridline-color: #ddd;
                border-radius: 4px;
            }
            QTableWidget::item {
                padding: 5px;
            }
            QPushButton {
                padding: 8px 16px;
                border-radius: 4px;
                font-weight: bold;
            }
            QTabWidget::pane {
                border: 1px solid #ccc;
                border-radius: 4px;
            }
            QTabBar::tab {
                padding: 8px 16px;
                margin-right: 2px;
            }
            QTabBar::tab:selected {
                background-color: #3498db;
                color: white;
            }
        """)
        
        # Central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        
        # Header
        header = QLabel("🧬 Amino Acid Generator - Real-Time Protein Synthesis & Structure Prediction")
        header.setFont(QFont("Arial", 16, QFont.Bold))
        header.setStyleSheet("""
            background-color: #2c3e50;
            color: white;
            padding: 15px;
            border-radius: 8px;
        """)
        header.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(header)
        
        # Connection status bar
        status_layout = QHBoxLayout()
        
        self.status_label = QLabel("🔴 Not Connected")
        self.status_label.setFont(QFont("Arial", 11))
        self.status_label.setStyleSheet("padding: 5px;")
        status_layout.addWidget(self.status_label)
        
        self.message_count_label = QLabel("Messages: 0")
        self.message_count_label.setFont(QFont("Arial", 11))
        status_layout.addWidget(self.message_count_label)
        
        status_layout.addStretch()
        
        self.connect_btn = QPushButton("🔌 Connect")
        self.connect_btn.setStyleSheet("background-color: #27ae60; color: white;")
        self.connect_btn.clicked.connect(self.connect_zmq)
        status_layout.addWidget(self.connect_btn)
        
        self.disconnect_btn = QPushButton("⏹ Disconnect")
        self.disconnect_btn.setStyleSheet("background-color: #e74c3c; color: white;")
        self.disconnect_btn.clicked.connect(self.disconnect_zmq)
        self.disconnect_btn.setEnabled(False)
        status_layout.addWidget(self.disconnect_btn)
        
        self.clear_btn = QPushButton("🗑 Clear All")
        self.clear_btn.setStyleSheet("background-color: #95a5a6; color: white;")
        self.clear_btn.clicked.connect(self.clear_all_data)
        status_layout.addWidget(self.clear_btn)
        
        main_layout.addLayout(status_layout)
        
        # Tab widget for Primary and Secondary structure
        self.tab_widget = QTabWidget()
        
        # Tab 1: Primary Structure (Amino Acid Sequence)
        primary_tab = QWidget()
        primary_layout = QVBoxLayout(primary_tab)
        
        # Splitter for Active and Completed sections
        splitter = QSplitter(Qt.Vertical)
        
        # Active Genes Section
        active_group = QGroupBox("🔄 Active Translation - Primary Structure")
        active_group.setStyleSheet("""
            QGroupBox {
                border-color: #e67e22;
            }
            QGroupBox::title {
                color: #e67e22;
            }
        """)
        active_layout = QVBoxLayout(active_group)
        
        self.active_count_label = QLabel("Active genes: 0")
        self.active_count_label.setFont(QFont("Arial", 10))
        active_layout.addWidget(self.active_count_label)
        
        self.active_table = QTableWidget()
        self.active_table.setColumnCount(4)
        self.active_table.setHorizontalHeaderLabels([
            "Gene ID", "Progress", "Amino Acids (Primary Structure)", "Full Length"
        ])
        self.active_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)
        self.active_table.setAlternatingRowColors(True)
        self.active_table.cellClicked.connect(self.on_active_gene_clicked)
        active_layout.addWidget(self.active_table)
        
        splitter.addWidget(active_group)
        
        # Completed Genes Section
        completed_group = QGroupBox("✅ Completed Proteins - Primary + Secondary Structure (LIVE)")
        completed_group.setStyleSheet("""
            QGroupBox {
                border-color: #27ae60;
            }
            QGroupBox::title {
                color: #27ae60;
            }
        """)
        completed_layout = QVBoxLayout(completed_group)
        
        self.completed_count_label = QLabel("Completed proteins: 0")
        self.completed_count_label.setFont(QFont("Arial", 10))
        completed_layout.addWidget(self.completed_count_label)
        
        self.completed_table = QTableWidget()
        self.completed_table.setColumnCount(5)
        self.completed_table.setHorizontalHeaderLabels([
            "Gene ID", "Amino Acids", "Secondary Structure", "Structure Stats", "Sequence"
        ])
        self.completed_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)
        self.completed_table.horizontalHeader().setSectionResizeMode(4, QHeaderView.Stretch)
        self.completed_table.setAlternatingRowColors(True)
        self.completed_table.cellClicked.connect(self.on_completed_gene_clicked)
        completed_layout.addWidget(self.completed_table)
        
        splitter.addWidget(completed_group)
        splitter.setSizes([400, 300])
        
        primary_layout.addWidget(splitter)
        
        self.tab_widget.addTab(primary_tab, "🔤 Primary Structure")
        
        # Tab 2: Secondary Structure
        secondary_tab = QWidget()
        secondary_layout = QVBoxLayout(secondary_tab)
        
        # Info label
        info = QLabel(
            "🔴 LIVE: Secondary structure updates automatically when a gene completes translation!\n"
            "Using Chou-Fasman method (~65% accuracy)"
        )
        info.setStyleSheet("color: #c0392b; font-weight: bold; padding: 10px; background-color: #fadbd8; border-radius: 5px;")
        secondary_layout.addWidget(info)
        
        # Secondary structure widget
        self.secondary_widget = SecondaryStructureWidget()
        secondary_layout.addWidget(self.secondary_widget)
        
        self.tab_widget.addTab(secondary_tab, "🌀 Secondary Structure")
        
        main_layout.addWidget(self.tab_widget)
        
        # Info footer
        info_label = QLabel(
            "💡 This app receives translation data from JCVI Genome Analyzer via ZeroMQ (port 5555), "
            "converts CDS sequences to amino acid sequences, and predicts secondary structure in real-time."
        )
        info_label.setStyleSheet("color: #666; font-style: italic; padding: 10px;")
        info_label.setWordWrap(True)
        main_layout.addWidget(info_label)
    
    def on_active_gene_clicked(self, row, column):
        """Handle click on active gene - show secondary structure"""
        gene_ids = list(self.active_genes.keys())
        if row < len(gene_ids):
            gene_id = gene_ids[row]
            self.show_secondary_structure(gene_id, self.active_genes[gene_id])
    
    def on_completed_gene_clicked(self, row, column):
        """Handle click on completed gene - show secondary structure"""
        gene_ids = list(self.completed_genes.keys())
        if row < len(gene_ids):
            gene_id = gene_ids[row]
            self.show_secondary_structure(gene_id, self.completed_genes[gene_id])
    
    def show_secondary_structure(self, gene_id: str, data: dict):
        """Show secondary structure for selected gene"""
        aa_list = get_amino_acid_list(data.get('full_amino_acids', data.get('amino_acids', '')))
        
        if not aa_list:
            return
        
        # Predict secondary structure
        structure, stats = predict_secondary_structure(aa_list)
        
        # Update display
        self.secondary_widget.update_display(gene_id, aa_list, structure, stats)
        
        # Switch to secondary structure tab
        self.tab_widget.setCurrentIndex(1)
    
    def connect_zmq(self):
        """Connect to ZeroMQ publisher"""
        if self.zmq_thread and self.zmq_thread.isRunning():
            return
        
        self.zmq_thread = ZMQReceiverThread(port=5555)
        self.zmq_thread.data_received.connect(self.on_data_received)
        self.zmq_thread.connection_status.connect(self.on_connection_status)
        self.zmq_thread.error_occurred.connect(self.on_error)
        self.zmq_thread.start()
        
        self.connect_btn.setEnabled(False)
        self.disconnect_btn.setEnabled(True)
    
    def disconnect_zmq(self):
        """Disconnect from ZeroMQ"""
        if self.zmq_thread:
            self.zmq_thread.stop()
            self.zmq_thread.wait(2000)
            self.zmq_thread = None
        
        self.connect_btn.setEnabled(True)
        self.disconnect_btn.setEnabled(False)
        self.status_label.setText("🔴 Disconnected")
    
    def on_connection_status(self, status: str):
        """Handle connection status updates"""
        self.status_label.setText(status)
    
    def on_error(self, error: str):
        """Handle errors"""
        print(f"Error: {error}")
        self.status_label.setText(f"⚠️ {error}")
    
    def on_data_received(self, data: dict):
        """Handle received translation data"""
        self.message_count += 1
        self.message_count_label.setText(f"Messages: {self.message_count}")
        
        gene_id = data.get('gene_id', '')
        protein_level = data.get('protein_level', 0)
        status = data.get('status', '')
        cds_sequence = data.get('cds_sequence', '')
        
        if not gene_id or not cds_sequence:
            return
        
        # Translate CDS to amino acids based on current progress
        amino_acids = translate_cds_to_amino_acids(cds_sequence, protein_level)
        full_amino_acids = translate_cds_to_amino_acids(cds_sequence, 100)
        
        # Predict secondary structure for full protein
        aa_list = get_amino_acid_list(full_amino_acids)
        structure, struct_stats = predict_secondary_structure(aa_list)
        
        # Store data
        gene_data = {
            'gene_id': gene_id,
            'protein_level': protein_level,
            'status': status,
            'cds_sequence': cds_sequence,
            'amino_acids': amino_acids,
            'full_amino_acids': full_amino_acids,
            'total_aa_count': len(aa_list),
            'secondary_structure': structure,
            'structure_stats': struct_stats
        }
        
        # Check if completed
        if protein_level >= 100 or status == 'complete':
            # Move to completed
            if gene_id in self.active_genes:
                del self.active_genes[gene_id]
            gene_data['amino_acids'] = full_amino_acids
            self.completed_genes[gene_id] = gene_data
        else:
            # Update active
            self.active_genes[gene_id] = gene_data
        
        # Update UI
        self.update_tables()
    
    def update_tables(self):
        """Update both tables with current data"""
        # Update Active Table
        self.active_table.setRowCount(len(self.active_genes))
        for row, (gene_id, data) in enumerate(self.active_genes.items()):
            # Gene ID
            id_item = QTableWidgetItem(gene_id)
            id_item.setFont(QFont("Consolas", 10, QFont.Bold))
            self.active_table.setItem(row, 0, id_item)
            
            # Progress
            progress = f"{data['protein_level']:.1f}%"
            progress_item = QTableWidgetItem(progress)
            progress_item.setTextAlignment(Qt.AlignCenter)
            # Color based on progress
            if data['protein_level'] < 30:
                progress_item.setBackground(QColor(255, 235, 235))
            elif data['protein_level'] < 70:
                progress_item.setBackground(QColor(255, 250, 230))
            else:
                progress_item.setBackground(QColor(235, 255, 235))
            self.active_table.setItem(row, 1, progress_item)
            
            # Amino Acids (truncate if too long)
            aa = data['amino_acids']
            if len(aa) > 100:
                aa = aa[:100] + "..."
            aa_item = QTableWidgetItem(aa)
            aa_item.setFont(QFont("Consolas", 9))
            self.active_table.setItem(row, 2, aa_item)
            
            # Full length
            length_item = QTableWidgetItem(f"{data['total_aa_count']} aa")
            length_item.setTextAlignment(Qt.AlignCenter)
            self.active_table.setItem(row, 3, length_item)
        
        self.active_count_label.setText(f"Active genes: {len(self.active_genes)}")
        
        # Update Completed Table with Secondary Structure (LIVE)
        self.completed_table.setRowCount(len(self.completed_genes))
        for row, (gene_id, data) in enumerate(self.completed_genes.items()):
            # Gene ID
            id_text = f"✅ {gene_id}"
            id_item = QTableWidgetItem(id_text)
            id_item.setFont(QFont("Consolas", 10, QFont.Bold))
            id_item.setForeground(QColor(39, 174, 96))
            self.completed_table.setItem(row, 0, id_item)
            
            # Total Amino Acids
            count_item = QTableWidgetItem(f"{data['total_aa_count']} aa")
            count_item.setTextAlignment(Qt.AlignCenter)
            self.completed_table.setItem(row, 1, count_item)
            
            # Secondary Structure (live display)
            structure = data.get('secondary_structure', '')
            if len(structure) > 50:
                structure_display = structure[:50] + "..."
            else:
                structure_display = structure
            struct_item = QTableWidgetItem(structure_display)
            struct_item.setFont(QFont("Consolas", 9))
            struct_item.setToolTip(f"Full structure: {structure}")
            self.completed_table.setItem(row, 2, struct_item)
            
            # Structure Stats (live display)
            stats = data.get('structure_stats', {})
            helix_pct = stats.get('helix_percent', 0)
            sheet_pct = stats.get('sheet_percent', 0)
            turn_pct = stats.get('turn_percent', 0)
            coil_pct = stats.get('coil_percent', 0)
            
            stats_text = f"H:{helix_pct:.0f}% E:{sheet_pct:.0f}% T:{turn_pct:.0f}% C:{coil_pct:.0f}%"
            stats_item = QTableWidgetItem(stats_text)
            stats_item.setTextAlignment(Qt.AlignCenter)
            stats_item.setFont(QFont("Consolas", 9))
            # Color code based on dominant structure
            if helix_pct > sheet_pct and helix_pct > coil_pct:
                stats_item.setBackground(QColor(255, 230, 230))  # Red tint for helix
            elif sheet_pct > helix_pct and sheet_pct > coil_pct:
                stats_item.setBackground(QColor(230, 240, 255))  # Blue tint for sheet
            self.completed_table.setItem(row, 3, stats_item)
            
            # Protein Sequence
            aa = data['full_amino_acids']
            if len(aa) > 60:
                aa = aa[:60] + "..."
            aa_item = QTableWidgetItem(aa)
            aa_item.setFont(QFont("Consolas", 8))
            self.completed_table.setItem(row, 4, aa_item)
        
        self.completed_count_label.setText(f"Completed proteins: {len(self.completed_genes)}")
        
        # Auto-update secondary structure visualization for the LATEST completed gene
        if self.completed_genes:
            latest_gene_id = list(self.completed_genes.keys())[-1]
            latest_data = self.completed_genes[latest_gene_id]
            aa_list = get_amino_acid_list(latest_data.get('full_amino_acids', ''))
            structure = latest_data.get('secondary_structure', '')
            stats = latest_data.get('structure_stats', {})
            if aa_list and structure:
                self.secondary_widget.update_display(latest_gene_id, aa_list, structure, stats)
    
    def clear_all_data(self):
        """Clear all data"""
        self.active_genes.clear()
        self.completed_genes.clear()
        self.message_count = 0
        self.message_count_label.setText("Messages: 0")
        self.update_tables()
    
    def closeEvent(self, event):
        """Handle window close"""
        self.disconnect_zmq()
        event.accept()


# =============================================================================
# Main Entry Point
# =============================================================================

def main():
    """Main entry point"""
    if not ZMQ_AVAILABLE:
        print("=" * 50)
        print("ERROR: PyZMQ is required!")
        print("Install with: pip install pyzmq")
        print("=" * 50)
        return
    
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    
    window = AminoAcidGeneratorApp()
    window.show()
    
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
