"""
Paused Gene Real-Time Dashboard for JCVI Genome Analyzer
=========================================================

Shows real-time decay dynamics when a gene's transcription is paused.
Displays mRNA and protein decay curves, state indicators, and rate values.
"""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QFrame,
    QGridLayout, QGroupBox, QSizePolicy
)
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QFont, QPainter, QPen, QColor, QPainterPath
from collections import deque
import math


class SimpleLineChart(QFrame):
    """Simple line chart widget using QPainter"""
    
    def __init__(self, title: str, color: QColor, max_points: int = 100, parent=None):
        super().__init__(parent)
        self.title = title
        self.line_color = color
        self.max_points = max_points
        self.data_points = deque(maxlen=max_points)
        self.time_points = deque(maxlen=max_points)
        self.current_value = 0.0
        self.setMinimumSize(250, 120)
        self.setFrameStyle(QFrame.Box | QFrame.Sunken)
        self.setStyleSheet("background-color: #1a1a2e;")
    
    def add_point(self, time_val: float, value: float):
        """Add a data point"""
        self.time_points.append(time_val)
        self.data_points.append(value)
        self.current_value = value
        self.update()
    
    def clear_data(self):
        """Clear all data points"""
        self.data_points.clear()
        self.time_points.clear()
        self.current_value = 0.0
        self.update()
    
    def paintEvent(self, event):
        super().paintEvent(event)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        
        # Margins
        margin = 10
        title_height = 20
        label_width = 40
        
        # Draw title
        painter.setPen(Qt.white)
        painter.setFont(QFont("Arial", 9, QFont.Bold))
        painter.drawText(margin, margin + 12, self.title)
        
        # Draw current value
        value_text = f"{self.current_value:.1f}"
        painter.setFont(QFont("Arial", 10))
        painter.setPen(self.line_color)
        painter.drawText(self.width() - 60, margin + 12, value_text)
        
        # Chart area
        chart_x = margin + label_width
        chart_y = margin + title_height
        chart_width = self.width() - 2 * margin - label_width
        chart_height = self.height() - 2 * margin - title_height - 10
        
        # Draw axis
        painter.setPen(QPen(QColor(80, 80, 80), 1))
        painter.drawLine(chart_x, chart_y + chart_height, 
                        chart_x + chart_width, chart_y + chart_height)
        painter.drawLine(chart_x, chart_y, chart_x, chart_y + chart_height)
        
        if len(self.data_points) < 2:
            # Draw "No data" text
            painter.setPen(Qt.gray)
            painter.drawText(chart_x + chart_width // 3, chart_y + chart_height // 2, 
                           "Waiting for data...")
            return
        
        # Calculate scaling
        max_val = max(self.data_points) if self.data_points else 1
        min_val = min(self.data_points) if self.data_points else 0
        if max_val == min_val:
            max_val = min_val + 1
        
        # Draw Y axis labels
        painter.setPen(Qt.gray)
        painter.setFont(QFont("Arial", 7))
        painter.drawText(margin, chart_y + 10, f"{max_val:.0f}")
        painter.drawText(margin, chart_y + chart_height, f"{min_val:.0f}")
        
        # Draw line
        painter.setPen(QPen(self.line_color, 2))
        path = QPainterPath()
        
        for i, (t, val) in enumerate(zip(self.time_points, self.data_points)):
            x = chart_x + (i / max(1, len(self.data_points) - 1)) * chart_width
            y = chart_y + chart_height - ((val - min_val) / (max_val - min_val)) * chart_height
            
            if i == 0:
                path.moveTo(x, y)
            else:
                path.lineTo(x, y)
        
        painter.drawPath(path)
        
        # Draw decay indicator arrow
        if len(self.data_points) >= 2:
            if self.data_points[-1] < self.data_points[-2]:
                painter.setPen(QPen(QColor(255, 100, 100), 2))
                painter.drawText(self.width() - 25, margin + 12, "↘")


class StatusIndicator(QFrame):
    """Color-coded status indicator"""
    
    STATUS_COLORS = {
        'active': QColor(46, 204, 113),      # Green
        'rising': QColor(52, 152, 219),       # Blue
        'declining': QColor(241, 196, 15),    # Yellow
        'declining_fast': QColor(230, 126, 34),  # Orange
        'depleted': QColor(231, 76, 60),      # Red
        'initializing': QColor(149, 165, 166),  # Gray
        'paused': QColor(155, 89, 182),       # Purple
    }
    
    def __init__(self, label: str, parent=None):
        super().__init__(parent)
        self.label = label
        self.status = "initializing"
        self.setFixedSize(150, 30)
    
    def set_status(self, status: str):
        self.status = status.lower()
        self.update()
    
    def paintEvent(self, event):
        super().paintEvent(event)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        
        # Draw label
        painter.setPen(Qt.white)
        painter.setFont(QFont("Arial", 9))
        painter.drawText(5, 20, self.label)
        
        # Draw status circle
        color = self.STATUS_COLORS.get(self.status, QColor(128, 128, 128))
        painter.setBrush(color)
        painter.setPen(QPen(color.darker(120), 1))
        painter.drawEllipse(100, 8, 14, 14)
        
        # Draw status text
        painter.setPen(color)
        painter.setFont(QFont("Arial", 8))
        status_text = self.status.replace('_', ' ').title()
        painter.drawText(120, 20, status_text[:10])


class PausedGeneDashboard(QDialog):
    """
    Real-time dashboard showing decay dynamics for a paused gene.
    
    Displays:
    - mRNA count decay curve
    - Protein count decay curve  
    - Protein state indicator
    - Structure prediction status
    - Rate values
    """
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent_app = parent
        self.gene_id = None
        self.pause_time = 0.0
        self.init_ui()
    
    def init_ui(self):
        """Initialize the dashboard UI"""
        self.setWindowTitle("📊 Paused Gene Decay Dashboard")
        self.setMinimumSize(500, 400)
        self.setStyleSheet("""
            QDialog {
                background-color: #0f0f23;
            }
            QLabel {
                color: white;
            }
            QGroupBox {
                color: #aaa;
                border: 1px solid #333;
                border-radius: 5px;
                margin-top: 10px;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px;
            }
        """)
        
        layout = QVBoxLayout(self)
        
        # Header
        self.header_label = QLabel("⏸️ No gene paused")
        self.header_label.setFont(QFont("Arial", 14, QFont.Bold))
        self.header_label.setStyleSheet("color: #00d4ff;")
        layout.addWidget(self.header_label)
        
        self.time_label = QLabel("Time since pause: 0.0s")
        self.time_label.setStyleSheet("color: #888;")
        layout.addWidget(self.time_label)
        
        # Charts section
        charts_layout = QHBoxLayout()
        
        # mRNA chart
        self.mrna_chart = SimpleLineChart("mRNA Count", QColor(52, 152, 219))
        charts_layout.addWidget(self.mrna_chart)
        
        # Protein chart
        self.protein_chart = SimpleLineChart("Protein Count", QColor(46, 204, 113))
        charts_layout.addWidget(self.protein_chart)
        
        layout.addLayout(charts_layout)
        
        # Status indicators section
        status_group = QGroupBox("Status Indicators")
        status_layout = QGridLayout(status_group)
        
        self.protein_state_indicator = StatusIndicator("Protein State:")
        status_layout.addWidget(self.protein_state_indicator, 0, 0)
        
        self.structure_pred_label = QLabel("Structure Prediction: ❓")
        status_layout.addWidget(self.structure_pred_label, 0, 1)
        
        self.sequence_valid_label = QLabel("Sequence Valid: ❓")
        status_layout.addWidget(self.sequence_valid_label, 1, 0)
        
        self.size_category_label = QLabel("Size Category: —")
        status_layout.addWidget(self.size_category_label, 1, 1)
        
        layout.addWidget(status_group)
        
        # Rates section
        rates_group = QGroupBox("Rates & Parameters")
        rates_layout = QGridLayout(rates_group)
        
        self.synthesis_rate_label = QLabel("Synthesis Rate: —")
        rates_layout.addWidget(self.synthesis_rate_label, 0, 0)
        
        self.degradation_rate_label = QLabel("Degradation Rate: —")
        rates_layout.addWidget(self.degradation_rate_label, 0, 1)
        
        self.mrna_halflife_label = QLabel("mRNA Half-life: —")
        rates_layout.addWidget(self.mrna_halflife_label, 1, 0)
        
        self.protein_halflife_label = QLabel("Protein Half-life: —")
        rates_layout.addWidget(self.protein_halflife_label, 1, 1)
        
        layout.addWidget(rates_group)
        
        # Numerical values section
        values_group = QGroupBox("Current Values")
        values_layout = QHBoxLayout(values_group)
        
        self.mrna_value_label = QLabel("mRNA: 0.0")
        self.mrna_value_label.setFont(QFont("Arial", 12, QFont.Bold))
        self.mrna_value_label.setStyleSheet("color: #3498db;")
        values_layout.addWidget(self.mrna_value_label)
        
        self.protein_value_label = QLabel("Protein: 0.0")
        self.protein_value_label.setFont(QFont("Arial", 12, QFont.Bold))
        self.protein_value_label.setStyleSheet("color: #2ecc71;")
        values_layout.addWidget(self.protein_value_label)
        
        self.protein_fraction_label = QLabel("Fraction: 0%")
        self.protein_fraction_label.setFont(QFont("Arial", 12))
        self.protein_fraction_label.setStyleSheet("color: #f1c40f;")
        values_layout.addWidget(self.protein_fraction_label)
        
        layout.addWidget(values_group)
    
    def set_gene(self, gene_id: str, pause_time: float):
        """Set the gene to monitor"""
        self.gene_id = gene_id
        self.pause_time = pause_time
        self.header_label.setText(f"⏸️ Gene: {gene_id}")
        self.time_label.setText(f"Paused at t={pause_time:.1f}s")
        
        # Clear previous data
        self.mrna_chart.clear_data()
        self.protein_chart.clear_data()
        
        self.show()
        self.raise_()
    
    def update_data(self, data: dict):
        """Update dashboard with new data point"""
        incoming_gene = data.get('gene_id')
        if incoming_gene != self.gene_id:
            # Ignore data for other genes
            return
        
        current_time = data.get('timestamp', 0.0)
        time_since_pause = current_time - self.pause_time
        
        # Update time label
        self.time_label.setText(f"Paused at t={self.pause_time:.1f}s | Elapsed: {time_since_pause:.1f}s")
        
        # Update charts
        mrna_count = data.get('mrna_count', 0.0)
        protein_count = data.get('protein_count', 0.0)
        
        self.mrna_chart.add_point(current_time, mrna_count)
        self.protein_chart.add_point(current_time, protein_count)
        
        # Update numerical values
        self.mrna_value_label.setText(f"mRNA: {mrna_count:.1f}")
        self.protein_value_label.setText(f"Protein: {protein_count:.1f}")
        
        protein_fraction = data.get('protein_fraction', 0.0)
        self.protein_fraction_label.setText(f"Fraction: {protein_fraction*100:.1f}%")
        
        # Update status indicators
        protein_state = data.get('protein_state', 'initializing')
        self.protein_state_indicator.set_status(protein_state)
        
        # Structure prediction
        should_predict = data.get('should_predict_structure', False)
        if should_predict:
            self.structure_pred_label.setText("Structure Prediction: ✅ Enabled")
            self.structure_pred_label.setStyleSheet("color: #2ecc71;")
        else:
            self.structure_pred_label.setText("Structure Prediction: ❌ Disabled")
            self.structure_pred_label.setStyleSheet("color: #e74c3c;")
        
        # Sequence validity
        sequence_valid = data.get('sequence_valid', False)
        if sequence_valid:
            self.sequence_valid_label.setText("Sequence Valid: ✅ Yes")
            self.sequence_valid_label.setStyleSheet("color: #2ecc71;")
        else:
            self.sequence_valid_label.setText("Sequence Valid: ❌ No")
            self.sequence_valid_label.setStyleSheet("color: #e74c3c;")
        
        # Size category
        size_cat = data.get('size_category', '—')
        self.size_category_label.setText(f"Size Category: {size_cat}")
        
        # Rates
        synthesis_rate = data.get('synthesis_rate', 0.0)
        degradation_rate = data.get('degradation_rate', 0.0)
        
        self.synthesis_rate_label.setText(f"Synthesis Rate: {synthesis_rate:.4f}")
        if synthesis_rate == 0:
            self.synthesis_rate_label.setStyleSheet("color: #e74c3c;")
        else:
            self.synthesis_rate_label.setStyleSheet("color: white;")
        
        self.degradation_rate_label.setText(f"Degradation Rate: {degradation_rate:.6f}/s")
        
        # Half-lives
        mrna_hl = data.get('mrna_half_life', 300.0)
        protein_hl = data.get('protein_half_life', 3600.0)
        
        self.mrna_halflife_label.setText(f"mRNA Half-life: {mrna_hl:.0f}s")
        self.protein_halflife_label.setText(f"Protein Half-life: {protein_hl:.0f}s")
    
    def clear_gene(self):
        """Clear the dashboard when gene is resumed"""
        self.gene_id = None
        self.header_label.setText("⏸️ No gene paused")
        self.mrna_chart.clear_data()
        self.protein_chart.clear_data()
