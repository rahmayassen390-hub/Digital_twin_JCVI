"""
Protein Synthesis Widget for JCVI Genome Analyzer
Integrated version of the Amino Acid Generator
"""

import json
from typing import Dict, Optional, List, Tuple

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
    QTableWidget, QTableWidgetItem, QPushButton, QGroupBox,
    QHeaderView, QSplitter, QTabWidget, QTextEdit
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont, QColor

# =============================================================================
# Genetic Code - Codon to Amino Acid mapping
# =============================================================================

CODON_TABLE = {
    'TTT': 'Phe', 'TTC': 'Phe',
    'TTA': 'Leu', 'TTG': 'Leu', 'CTT': 'Leu', 'CTC': 'Leu', 'CTA': 'Leu', 'CTG': 'Leu',
    'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile',
    'ATG': 'Met',
    'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val',
    'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser', 'AGT': 'Ser', 'AGC': 'Ser',
    'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'TAT': 'Tyr', 'TAC': 'Tyr',
    'TAA': 'Stop', 'TAG': 'Stop', 'TGA': 'Stop',
    'CAT': 'His', 'CAC': 'His',
    'CAA': 'Gln', 'CAG': 'Gln',
    'AAT': 'Asn', 'AAC': 'Asn',
    'AAA': 'Lys', 'AAG': 'Lys',
    'GAT': 'Asp', 'GAC': 'Asp',
    'GAA': 'Glu', 'GAG': 'Glu',
    'TGT': 'Cys', 'TGC': 'Cys',
    'TGG': 'Trp',
    'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'AGA': 'Arg', 'AGG': 'Arg',
    'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
}

ONE_TO_THREE = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
}

THREE_TO_ONE = {v: k for k, v in ONE_TO_THREE.items()}

# =============================================================================
# Chou-Fasman Parameters
# =============================================================================

HELIX_PROPENSITY = {
    'Ala': 1.42, 'Arg': 0.98, 'Asn': 0.67, 'Asp': 1.01, 'Cys': 0.70,
    'Gln': 1.11, 'Glu': 1.51, 'Gly': 0.57, 'His': 1.00, 'Ile': 1.08,
    'Leu': 1.21, 'Lys': 1.16, 'Met': 1.45, 'Phe': 1.13, 'Pro': 0.57,
    'Ser': 0.77, 'Thr': 0.83, 'Trp': 1.08, 'Tyr': 0.69, 'Val': 1.06
}

SHEET_PROPENSITY = {
    'Ala': 0.83, 'Arg': 0.93, 'Asn': 0.89, 'Asp': 0.54, 'Cys': 1.19,
    'Gln': 1.10, 'Glu': 0.37, 'Gly': 0.75, 'His': 0.87, 'Ile': 1.60,
    'Leu': 1.30, 'Lys': 0.74, 'Met': 1.05, 'Phe': 1.38, 'Pro': 0.55,
    'Ser': 0.75, 'Thr': 1.19, 'Trp': 1.37, 'Tyr': 1.47, 'Val': 1.70
}

TURN_PROPENSITY = {
    'Ala': 0.66, 'Arg': 0.95, 'Asn': 1.56, 'Asp': 1.46, 'Cys': 1.19,
    'Gln': 0.98, 'Glu': 0.74, 'Gly': 1.56, 'His': 0.95, 'Ile': 0.47,
    'Leu': 0.59, 'Lys': 1.01, 'Met': 0.60, 'Phe': 0.60, 'Pro': 1.52,
    'Ser': 1.43, 'Thr': 0.96, 'Trp': 0.96, 'Tyr': 1.14, 'Val': 0.50
}

HELIX_THRESHOLD = 1.03
SHEET_THRESHOLD = 1.05
TURN_THRESHOLD = 1.00

def predict_secondary_structure(amino_acids: List[str]) -> Tuple[str, Dict]:
    if not amino_acids: return "", {}
    n = len(amino_acids)
    structure = ['C'] * n
    helix_scores = [HELIX_PROPENSITY.get(aa, 1.0) for aa in amino_acids]
    sheet_scores = [SHEET_PROPENSITY.get(aa, 1.0) for aa in amino_acids]
    turn_scores = [TURN_PROPENSITY.get(aa, 1.0) for aa in amino_acids]
    
    for i in range(n - 5):
        window = helix_scores[i:i+6]
        if sum(1 for p in window if p > HELIX_THRESHOLD) >= 4:
            for j in range(i, min(i+6, n)):
                if helix_scores[j] > 0.9: structure[j] = 'H'
    
    for i in range(n - 4):
        window = sheet_scores[i:i+5]
        if sum(1 for p in window if p > SHEET_THRESHOLD) >= 3:
            for j in range(i, min(i+5, n)):
                if structure[j] == 'C' and sheet_scores[j] > 0.9: structure[j] = 'E'
                
    for i in range(n - 3):
        window = turn_scores[i:i+4]
        if sum(window)/4 > TURN_THRESHOLD:
            if all(structure[j] == 'C' for j in range(i, i+4)):
                for j in range(i, i+4): structure[j] = 'T'
                
    ss = ''.join(structure)
    stats = {
        'helix_count': ss.count('H'), 'sheet_count': ss.count('E'),
        'turn_count': ss.count('T'), 'coil_count': ss.count('C'),
        'total': n,
        'helix_percent': (ss.count('H')/n*100) if n > 0 else 0,
        'sheet_percent': (ss.count('E')/n*100) if n > 0 else 0,
        'turn_percent': (ss.count('T')/n*100) if n > 0 else 0,
        'coil_percent': (ss.count('C')/n*100) if n > 0 else 0
    }
    return ss, stats

def translate_cds(cds: str, percentage: float = 100.0) -> str:
    if not cds: return ""
    cds = cds.upper().replace(' ', '').replace('\n', '')
    total = len(cds) // 3
    to_trans = int(total * (percentage / 100.0))
    if to_trans == 0: return ""
    aas = []
    for i in range(to_trans):
        codon = cds[i*3:(i+1)*3]
        if len(codon) < 3: break
        aa = CODON_TABLE.get(codon, 'Unk')
        if aa == 'Stop': break
        aas.append(aa)
    
    # Ensure at least one AA is shown if synthesis has started
    if not aas and to_trans > 0:
        codon = cds[:3]
        if len(codon) == 3:
            aa = CODON_TABLE.get(codon, 'Unk')
            if aa != 'Stop':
                aas.append(aa)
                
    return '-'.join(aas)

# =============================================================================
# UI Widgets
# =============================================================================

class SecondaryStructureWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        self.gene_label = QLabel("Select a gene to view secondary structure")
        self.gene_label.setFont(QFont("Arial", 12, QFont.Bold))
        layout.addWidget(self.gene_label)
        self.display = QTextEdit()
        self.display.setReadOnly(True)
        self.display.setFont(QFont("Consolas", 10))
        layout.addWidget(self.display)
        self.stats_label = QLabel("")
        layout.addWidget(self.stats_label)

    def update_display(self, gene_id: str, amino_acids: List[str], structure: str, stats: Dict):
        self.gene_label.setText(f"🧬 Gene: {gene_id}")
        html = ["<pre style='font-family: Consolas; font-size: 11px;'><b>Sequence:</b>\n"]
        for i in range(0, len(amino_acids), 60):
            html.append(f"<span style='color: #888;'>{i+1:5d}</span>  ")
            line_aa = amino_acids[i:i+60]
            line_struct = structure[i:i+60]
            for j, aa in enumerate(line_aa):
                char = THREE_TO_ONE.get(aa, '?')
                s = line_struct[j] if j < len(line_struct) else 'C'
                color = {"H": "#e74c3c", "E": "#3498db", "T": "#2ecc71"}.get(s, "#666666")
                html.append(f"<span style='color: {color};'>{char}</span>")
            html.append("\n        ")
            for char in line_struct:
                color = {"H": "#e74c3c", "E": "#3498db", "T": "#2ecc71"}.get(char, "#666666")
                html.append(f"<span style='color: {color};'>{char}</span>")
            html.append("\n\n")
        html.append("</pre>")
        self.display.setHtml(''.join(html))
        self.stats_label.setText(f"📊 <b>Stats:</b> H:{stats['helix_percent']:.1f}% E:{stats['sheet_percent']:.1f}% T:{stats['turn_percent']:.1f}% C:{stats['coil_percent']:.1f}%")

class ProteinSynthesisWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.active_genes = {}
        self.completed_genes = {}
        self._gene_id_to_row = {'active': {}, 'completed': {}}
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout(self)
        self.tabs = QTabWidget()
        layout.addWidget(self.tabs)
        
        self.selected_gene_id = None
        
        # Primary Tab
        primary = QWidget()
        p_layout = QVBoxLayout(primary)
        splitter = QSplitter(Qt.Vertical)
        
        self.active_table = QTableWidget(0, 4)
        self.active_table.setHorizontalHeaderLabels(["Gene ID", "Progress", "Sequence", "Length"])
        self.active_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)
        self.active_table.cellClicked.connect(self.on_active_clicked)
        splitter.addWidget(self.active_table)
        
        self.completed_table = QTableWidget(0, 5)
        self.completed_table.setHorizontalHeaderLabels(["Gene ID", "AA Count", "Helix%", "Sheet%", "Sequence"])
        self.completed_table.horizontalHeader().setSectionResizeMode(4, QHeaderView.Stretch)
        self.completed_table.cellClicked.connect(self.on_completed_clicked)
        self.completed_table.cellDoubleClicked.connect(self.on_completed_double_clicked)
        splitter.addWidget(self.completed_table)
        
        # Add a full sequence viewer at the bottom of the primary tab
        viewer_group = QGroupBox("🧬 Detailed Structural Analysis (Primary + Secondary)")
        v_layout = QVBoxLayout(viewer_group)
        
        h_ctrl = QHBoxLayout()
        self.viewer_info = QLabel("Click a gene to analyze...")
        h_ctrl.addWidget(self.viewer_info)
        h_ctrl.addStretch()
        self.view_detail_btn = QPushButton("🔍 Open Full Detail Tab")
        self.view_detail_btn.setEnabled(False)
        self.view_detail_btn.clicked.connect(self.open_current_in_detail_tab)
        h_ctrl.addWidget(self.view_detail_btn)
        v_layout.addLayout(h_ctrl)
        
        self.full_sequence_viewer = QTextEdit()
        self.full_sequence_viewer.setReadOnly(True)
        self.full_sequence_viewer.setFont(QFont("Consolas", 10))
        v_layout.addWidget(self.full_sequence_viewer)
        
        splitter.addWidget(viewer_group)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 1)
        splitter.setStretchFactor(2, 2)
        
        p_layout.addWidget(splitter)
        self.tabs.addTab(primary, "🔤 Primary Structure")
        
        # Secondary Tab
        self.ss_widget = SecondaryStructureWidget()
        self.tabs.addTab(self.ss_widget, "🌀 Secondary Structure")

    def add_row(self, data: dict):
        gene_id = data.get('gene_id')
        level = data.get('protein_current', 0) / data.get('protein_total', 1) * 100 if data.get('protein_total') else 0
        cds = data.get('cds_sequence')
        if not gene_id or not cds: return
        
        aas = translate_cds(cds, level)
        full_aas = translate_cds(cds, 100)
        aa_list = full_aas.split('-') if full_aas else []
        ss, stats = predict_secondary_structure(aa_list)
        
        entry = {
            'gene_id': gene_id, 'level': level, 'aas': aas, 'full_aas': full_aas,
            'aa_count': len(aa_list), 'ss': ss, 'stats': stats
        }
        
        if level >= 100:
            if gene_id in self.active_genes:
                self.active_genes.pop(gene_id)
                self.rebuild_active_table_map()
            
            is_new = gene_id not in self.completed_genes
            self.completed_genes[gene_id] = entry
            self.update_completed_row(gene_id, entry)
            
            # AUTOMATIC SELECTION: If it just finished, show it in the viewer
            if is_new:
                self.selected_gene_id = gene_id
                self.update_viewer(entry)
        else:
            self.active_genes[gene_id] = entry
            self.update_active_row(gene_id, entry)
        
        # Update viewer if this is the currently selected gene
        if self.selected_gene_id == gene_id:
            self.update_viewer(entry)

    def rebuild_active_table_map(self):
        """Rebuild mapping when an active gene is removed"""
        self.active_table.setRowCount(len(self.active_genes))
        self._gene_id_to_row['active'] = {}
        for i, gid in enumerate(self.active_genes.keys()):
            self._gene_id_to_row['active'][gid] = i
        self.refresh_tables()

    def update_active_row(self, gene_id, d):
        """Update or add a single row in the active table"""
        if gene_id not in self._gene_id_to_row['active']:
            row = self.active_table.rowCount()
            self.active_table.insertRow(row)
            self._gene_id_to_row['active'][gene_id] = row
            self.active_table.setItem(row, 0, QTableWidgetItem(gene_id))
        
        row = self._gene_id_to_row['active'][gene_id]
        self.active_table.setItem(row, 1, QTableWidgetItem(f"{d['level']:.1f}%"))
        
        seq_item = QTableWidgetItem(d['aas'])
        seq_item.setToolTip(d['aas'])
        self.active_table.setItem(row, 2, seq_item)
        self.active_table.setItem(row, 3, QTableWidgetItem(f"{d['aa_count']}"))

    def update_completed_row(self, gene_id, d):
        """Update or add a single row in the completed table"""
        if gene_id not in self._gene_id_to_row['completed']:
            row = self.completed_table.rowCount()
            self.completed_table.insertRow(row)
            self._gene_id_to_row['completed'][gene_id] = row
            self.completed_table.setItem(row, 0, QTableWidgetItem(gene_id))
        
        row = self._gene_id_to_row['completed'][gene_id]
        self.completed_table.setItem(row, 1, QTableWidgetItem(f"{d['aa_count']}"))
        self.completed_table.setItem(row, 2, QTableWidgetItem(f"{d['stats']['helix_percent']:.1f}%"))
        self.completed_table.setItem(row, 3, QTableWidgetItem(f"{d['stats']['sheet_percent']:.1f}%"))
        
        full_seq_item = QTableWidgetItem(d['full_aas'])
        full_seq_item.setToolTip(d['full_aas'])
        self.completed_table.setItem(row, 4, full_seq_item)

    def refresh_tables(self):
        self.active_table.setRowCount(len(self.active_genes))
        for i, (gid, d) in enumerate(self.active_genes.items()):
            self.active_table.setItem(i, 0, QTableWidgetItem(gid))
            self.active_table.setItem(i, 1, QTableWidgetItem(f"{d['level']:.1f}%"))
            # Show actual growing sequence in the table row with a tooltip for the full string
            seq_item = QTableWidgetItem(d['aas'])
            seq_item.setToolTip(d['aas'])
            self.active_table.setItem(i, 2, seq_item)
            self.active_table.setItem(i, 3, QTableWidgetItem(f"{d['aa_count']}"))
            
        self.completed_table.setRowCount(len(self.completed_genes))
        for i, (gid, d) in enumerate(self.completed_genes.items()):
            self.completed_table.setItem(i, 0, QTableWidgetItem(gid))
            self.completed_table.setItem(i, 1, QTableWidgetItem(f"{d['aa_count']}"))
            self.completed_table.setItem(i, 2, QTableWidgetItem(f"{d['stats']['helix_percent']:.1f}%"))
            self.completed_table.setItem(i, 3, QTableWidgetItem(f"{d['stats']['sheet_percent']:.1f}%"))
            # Show full sequence in completed table with tooltip
            full_seq_item = QTableWidgetItem(d['full_aas'])
            full_seq_item.setToolTip(d['full_aas'])
            self.completed_table.setItem(i, 4, full_seq_item)

    def on_active_clicked(self, r, c):
        gid = list(self.active_genes.keys())[r]
        self.selected_gene_id = gid
        self.update_viewer(self.active_genes[gid])

    def on_completed_clicked(self, r, c):
        gid = list(self.completed_genes.keys())[r]
        self.selected_gene_id = gid
        self.update_viewer(self.completed_genes[gid])
        
    def on_completed_double_clicked(self, r, c):
        self.open_current_in_detail_tab()

    def open_current_in_detail_tab(self):
        if not self.selected_gene_id: return
        d = self.active_genes.get(self.selected_gene_id) or self.completed_genes.get(self.selected_gene_id)
        if d:
            self.show_ss(self.selected_gene_id, d)

    def update_viewer(self, d):
        gid = d['gene_id']
        level = d['level']
        aas = d['aas']
        full_aas = d['full_aas']
        
        self.viewer_info.setText(f"🧬 Analyzing <b>{gid}</b> ({level:.1f}% complete)")
        self.view_detail_btn.setEnabled(True)
        
        # Perform "live" prediction for the current synthesized part
        current_aas = aas.split('-') if aas else []
        ss, stats = predict_secondary_structure(current_aas)
        
        # Build Aligned View
        html = ["<div style='font-family: Consolas; font-size: 11px;'>"]
        
        full_list = full_aas.split('-') if full_aas else []
        trans_len = len(current_aas)
        
        for i in range(0, len(full_list), 60):
            # Line 1: Amino Acids
            line_aas = full_list[i:i+60]
            aa_row = [f"<span style='color: #888;'>{i+1:4d} AA: </span>"]
            for j, aa in enumerate(line_aas):
                if aa == 'Stop':
                    char = '*'
                else:
                    char = THREE_TO_ONE.get(aa, '?')
                color = "#2ecc71" if (i+j) < trans_len else "#888888"
                aa_row.append(f"<span style='color: {color};'>{char}</span>")
                if (j+1) % 10 == 0: aa_row.append(" ")
            html.append("".join(aa_row) + "<br>")
            
            # Line 2: Secondary Structure (only for translated part)
            ss_row = [f"<span style='color: #888;'>     SS: </span>"]
            for j in range(len(line_aas)):
                idx = i + j
                if idx < trans_len:
                    char = ss[idx]
                    color = {"H": "#e74c3c", "E": "#3498db", "T": "#2ecc71"}.get(char, "#666666")
                    ss_row.append(f"<span style='color: {color};'>{char}</span>")
                else:
                    ss_row.append("<span style='color: #444;'>.</span>")
                if (j+1) % 10 == 0: ss_row.append(" ")
            html.append("".join(ss_row) + "<br><br>")
            
        html.append("</div>")
        self.full_sequence_viewer.setHtml("".join(html))

    def show_ss(self, gid, d):
        aa_list = d['full_aas'].split('-')
        self.ss_widget.update_display(gid, aa_list, d['ss'], d['stats'])
        self.tabs.setCurrentIndex(1)

    def clear(self):
        self.active_genes.clear()
        self.completed_genes.clear()
        self.refresh_tables()
