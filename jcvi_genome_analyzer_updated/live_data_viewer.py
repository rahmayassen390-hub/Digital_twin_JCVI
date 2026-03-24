"""
Live Data Viewer Dialog for JCVI Genome Analyzer
Shows translation data in real-time as genes are translated
"""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QPushButton, 
    QLabel, QTableWidget, QTableWidgetItem, QHeaderView
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont, QColor


class LiveDataViewerDialog(QDialog):
    """Live data viewer showing translation data in real-time"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent_app = parent
        self.init_ui()
        
    def init_ui(self):
        """Initialize the UI"""
        self.setWindowTitle("📊 Live Translation Data - JCVI Genome Analyzer")
        self.setGeometry(120, 120, 800, 550)
        
        layout = QVBoxLayout(self)
        
        # Title
        title = QLabel("📊 Live Translation Data Viewer")
        title.setFont(QFont("Arial", 14, QFont.Bold))
        title.setStyleSheet("color: #2196F3; background-color: #E3F2FD; padding: 10px;")
        layout.addWidget(title)
        
        # Info bar
        info_layout = QHBoxLayout()
        
        self.row_count_label = QLabel("Rows: 0")
        self.row_count_label.setFont(QFont("Arial", 11))
        info_layout.addWidget(self.row_count_label)
        
        self.status_label = QLabel("⏸️ Waiting for translation to start...")
        self.status_label.setStyleSheet("color: #666;")
        info_layout.addWidget(self.status_label)
        
        info_layout.addStretch()
        
        # Clear button
        self.clear_btn = QPushButton("🗑️ Clear Data")
        self.clear_btn.clicked.connect(self.clear_data)
        info_layout.addWidget(self.clear_btn)
        
        # Auto-scroll toggle
        self.auto_scroll_btn = QPushButton("📜 Auto-Scroll: ON")
        self.auto_scroll_btn.setCheckable(True)
        self.auto_scroll_btn.setChecked(True)
        self.auto_scroll_btn.clicked.connect(self.toggle_auto_scroll)
        info_layout.addWidget(self.auto_scroll_btn)
        
        layout.addLayout(info_layout)
        
        # Data table with all biological columns
        self.data_table = QTableWidget()
        self.data_table.setColumnCount(14)  # Added Protein State column
        self.data_table.setHorizontalHeaderLabels([
            "Time", "Gene ID", "Status", "Prot State",
            "mRNA (Cur/Tot)", "TX Rate",
            "Protein (Cur/Tot)", "TL Rate",
            "Gene Len", "Prot Len", "Promoter",
            "Strand", "Start", "End"
        ])
        
        # Set column widths
        header = self.data_table.horizontalHeader()
        for i in range(14):  # Updated count
            header.setSectionResizeMode(i, QHeaderView.ResizeToContents)
        
        self.data_table.setAlternatingRowColors(True)
        self.data_table.setStyleSheet("""
            QTableWidget {
                alternate-background-color: #f5f5f5;
                font-size: 11px;
            }
        """)
        layout.addWidget(self.data_table)
        
        # ZMQ status
        self.zmq_status_label = QLabel("🔌 ZMQ Publisher: Not initialized")
        self.zmq_status_label.setStyleSheet("color: #999; font-size: 10px;")
        layout.addWidget(self.zmq_status_label)
        
        # Auto-scroll flag
        self.auto_scroll = True
    
    def toggle_auto_scroll(self):
        """Toggle auto-scroll feature"""
        self.auto_scroll = self.auto_scroll_btn.isChecked()
        if self.auto_scroll:
            self.auto_scroll_btn.setText("📜 Auto-Scroll: ON")
        else:
            self.auto_scroll_btn.setText("📜 Auto-Scroll: OFF")
    
    def add_row(self, data: dict):
        """Add a new row to the table (log-style, always inserts new row)"""
        row = self.data_table.rowCount()
        self.data_table.insertRow(row)
        
        # Time
        self.data_table.setItem(row, 0, QTableWidgetItem(f"{data.get('timestamp', 0):.1f}s"))
        
        # Gene ID
        self.data_table.setItem(row, 1, QTableWidgetItem(data.get('gene_id', '')))
        
        # Status with color
        status = data.get('status', 'waiting')
        status_item = QTableWidgetItem(status)
        if status == 'translating':
            status_item.setBackground(QColor(255, 255, 200))  # Yellow
        elif status == 'complete':
            status_item.setBackground(QColor(200, 255, 200))  # Green
        else:
            status_item.setBackground(QColor(240, 240, 240))  # Gray
        self.data_table.setItem(row, 2, status_item)
        
        # Protein State (NEW - for downstream modules)
        protein_state = data.get('protein_state', 'unknown')
        state_icons = {
            'depleted': '❌ Depleted',
            'declining_fast': '📉 Declining Fast',
            'declining': '📉 Declining',
            'initializing': '⏳ Initializing',
            'rising': '📈 Rising',
            'active': '✅ Active'
        }
        state_item = QTableWidgetItem(state_icons.get(protein_state, protein_state))
        # Color coding for protein state
        if protein_state == 'depleted':
            state_item.setBackground(QColor(255, 200, 200))  # Red - no protein
        elif protein_state in ['declining', 'declining_fast']:
            state_item.setBackground(QColor(255, 230, 200))  # Orange - decaying
        elif protein_state == 'active':
            state_item.setBackground(QColor(200, 255, 200))  # Green - steady state
        elif protein_state == 'rising':
            state_item.setBackground(QColor(255, 255, 200))  # Yellow - accumulating
        else:
            state_item.setBackground(QColor(240, 240, 240))  # Gray - initializing
        self.data_table.setItem(row, 3, state_item)
        
        # mRNA (Current/Total)
        mrna_cur = data.get('mrna_current', 0)
        mrna_tot = data.get('mrna_copies', 0)
        mrna_item = QTableWidgetItem(f"{mrna_cur:.0f}/{mrna_tot}")
        if mrna_cur >= mrna_tot * 0.9:
            mrna_item.setBackground(QColor(200, 255, 200))  # Complete
        elif mrna_cur > 0:
            mrna_item.setBackground(QColor(255, 255, 200))  # In progress
        self.data_table.setItem(row, 4, mrna_item)
        
        # TX Rate (copies/sec)
        tx_rate = data.get('transcription_rate', 0)
        self.data_table.setItem(row, 5, QTableWidgetItem(f"{tx_rate:.2f}/s"))
        
        # Protein (Current/Total)
        prot_cur = data.get('protein_current', 0)
        prot_tot = data.get('protein_total', 0)
        prot_item = QTableWidgetItem(f"{prot_cur:.0f}/{prot_tot}")
        if prot_cur >= prot_tot * 0.9:
            prot_item.setBackground(QColor(200, 255, 200))  # Complete
        elif prot_cur > 0:
            prot_item.setBackground(QColor(255, 255, 200))  # In progress
        self.data_table.setItem(row, 6, prot_item)
        
        # TL Rate (proteins/sec)
        tl_rate = data.get('translation_rate', 0)
        self.data_table.setItem(row, 7, QTableWidgetItem(f"{tl_rate:.1f}/s"))
        
        # Gene Length
        self.data_table.setItem(row, 8, QTableWidgetItem(f"{data.get('gene_length', 0)} bp"))
        
        # Protein Length
        self.data_table.setItem(row, 9, QTableWidgetItem(f"{data.get('protein_length', 0)} aa"))
        
        # Promoter Type
        promoter = data.get('promoter_type', 'None')
        prom_item = QTableWidgetItem(promoter)
        if promoter == 'BOTH':
            prom_item.setBackground(QColor(200, 255, 200))  # Strong
        elif promoter == 'ONE':
            prom_item.setBackground(QColor(255, 255, 200))  # Weak
        self.data_table.setItem(row, 10, prom_item)
        
        # Strand
        self.data_table.setItem(row, 11, QTableWidgetItem(data.get('strand', '+')))
        
        # Start
        self.data_table.setItem(row, 12, QTableWidgetItem(str(data.get('start', 0))))
        
        # End
        self.data_table.setItem(row, 13, QTableWidgetItem(str(data.get('end', 0))))
        
        # Update row count
        self.row_count_label.setText(f"Rows: {row + 1}")
        
        # Auto-scroll to bottom
        if self.auto_scroll:
            self.data_table.scrollToBottom()
    
    def clear_data(self):
        """Clear all rows from the table"""
        self.data_table.setRowCount(0)
        self.row_count_label.setText("Rows: 0")
        self.status_label.setText("⏸️ Data cleared")
    
    def set_status(self, status: str):
        """Update the status label"""
        self.status_label.setText(status)
    
    def set_zmq_status(self, status: str):
        """Update ZMQ status label"""
        self.zmq_status_label.setText(status)
    
    def closeEvent(self, event):
        """Handle window close event"""
        # Just hide instead of close so it can be reopened
        event.ignore()
        self.hide()
