"""
Separate Translation Animation Window
Displays translation dynamics in a standalone window for better visibility
"""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QPushButton, 
    QLabel, QTableWidget, QDoubleSpinBox
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont


class TranslationWindow(QDialog):
    """Standalone window for translation animation"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent_app = parent
        self.init_ui()
        
    def init_ui(self):
        """Initialize the UI"""
        self.setWindowTitle("🔬 Translation Animation - JCVI Genome Analyzer")
        self.setGeometry(100, 100, 850, 600)
        
        layout = QVBoxLayout(self)
        
        # Title
        title = QLabel("Agent #3: Live Translation Animation 🔬 - QUERY GENOME (syn3.0)")
        title.setFont(QFont("Arial", 14, QFont.Bold))
        title.setStyleSheet("color: #FF5722; background-color: #FBE9E7; padding: 10px;")
        layout.addWidget(title)
        
        # Controls
        control_layout = QHBoxLayout()
        
        self.start_btn = QPushButton("▶️ Start Translation")
        self.start_btn.setFont(QFont("Arial", 11, QFont.Bold))
        self.start_btn.setMinimumHeight(40)
        self.start_btn.setStyleSheet("background-color: #FF5722; color: white;")
        control_layout.addWidget(self.start_btn)
        
        self.stop_btn = QPushButton("⏸️ Pause")
        self.stop_btn.setEnabled(False)
        control_layout.addWidget(self.stop_btn)
        
        self.reset_btn = QPushButton("🔄 Reset")
        control_layout.addWidget(self.reset_btn)
        
        layout.addLayout(control_layout)
        
        # Time and speed controls
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
        layout.addLayout(time_layout)
        
        # Translation table
        self.translation_table = QTableWidget()
        self.translation_table.setColumnCount(10)
        self.translation_table.setHorizontalHeaderLabels([
            "Gene ID", "Protein %", "Phase", "Category", "RBS", 
            "Expression Status", "Reason", "Order", "Status", "Current"
        ])
        self.translation_table.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(self.translation_table)
        
        # Info label
        info = QLabel("💡 Tip: When you pause a gene in the Transcription window, it will be highlighted here in yellow")
        info.setStyleSheet("color: #666; font-style: italic; padding: 5px;")
        layout.addWidget(info)
    
    def closeEvent(self, event):
        """Handle window close event"""
        # Just hide instead of close so it can be reopened
        event.ignore()
        self.hide()
