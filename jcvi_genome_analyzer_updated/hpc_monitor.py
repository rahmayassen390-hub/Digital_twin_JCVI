from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
    QProgressBar, QFrame, QGridLayout
)
from PyQt5.QtCore import Qt, QSize
from PyQt5.QtGui import QFont, QColor

class MetricBar(QWidget):
    """Compact progress bar with label and value"""
    def __init__(self, label_text, unit="%", parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 5)
        layout.setSpacing(2)
        
        header = QHBoxLayout()
        self.label = QLabel(label_text)
        self.label.setFont(QFont("Segoe UI", 8, QFont.Bold))
        self.label.setStyleSheet("color: #aaa;")
        
        self.value_label = QLabel(f"0{unit}")
        self.value_label.setFont(QFont("Consolas", 8))
        self.value_label.setStyleSheet("color: #fff;")
        self.value_label.setAlignment(Qt.AlignRight)
        
        header.addWidget(self.label)
        header.addStretch()
        header.addWidget(self.value_label)
        layout.addLayout(header)
        
        self.bar = QProgressBar()
        self.bar.setFixedHeight(6)
        self.bar.setTextVisible(False)
        self.bar.setRange(0, 100)
        self.bar.setStyleSheet("""
            QProgressBar {
                background-color: #1a1a1a;
                border: none;
                border-radius: 3px;
            }
            QProgressBar::chunk {
                background-color: #4CAF50;
                border-radius: 3px;
            }
        """)
        layout.addWidget(self.bar)
        self.unit = unit

    def update_value(self, value):
        self.bar.setValue(int(value))
        self.value_label.setText(f"{value:.1f}{self.unit}")
        
        # Color thresholds
        if value < 70:
            color = "#4CAF50" # Green
        elif value < 90:
            color = "#FF9800" # Orange
        else:
            color = "#F44336" # Red
            
        self.bar.setStyleSheet(f"""
            QProgressBar {{ background-color: #1a1a1a; border: none; border-radius: 3px; }}
            QProgressBar::chunk {{ background-color: {color}; border-radius: 3px; }}
        """)

class GPUCard(QFrame):
    """Compact card for a single GPU"""
    def __init__(self, index, parent=None):
        super().__init__(parent)
        self.setFrameShape(QFrame.StyledPanel)
        self.setFixedWidth(140)
        self.setStyleSheet("""
            GPUCard {
                background-color: #1e1e1e;
                border: 1px solid #333;
                border-radius: 6px;
                padding: 4px;
            }
        """)
        
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)
        
        header = QLabel(f"GPU {index}")
        header.setFont(QFont("Segoe UI", 9, QFont.Bold))
        header.setStyleSheet("color: #4fc3f7;")
        layout.addWidget(header)
        
        self.util_bar = MetricBar("Load", "%")
        self.temp_bar = MetricBar("Temp", "°C")
        layout.addWidget(self.util_bar)
        layout.addWidget(self.temp_bar)

    def update_metrics(self, util, temp):
        self.util_bar.update_value(util)
        self.temp_bar.update_value(temp)

class HPCMonitorWidget(QWidget):
    """Unified HPC Monitoring Panel Sidebar"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFixedWidth(160)
        self.setStyleSheet("background-color: #121212; border-left: 1px solid #222;")
        
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 15, 10, 10)
        layout.setSpacing(15)
        
        title = QLabel("HPC MONITOR")
        title.setFont(QFont("Segoe UI Black", 10))
        title.setStyleSheet("color: #666; letter-spacing: 1px;")
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)
        
        # System Stats
        sys_group = QVBoxLayout()
        sys_group.setSpacing(8)
        self.cpu_bar = MetricBar("CPU CLUSTER")
        self.ram_bar = MetricBar("SYSTEM RAM")
        sys_group.addWidget(self.cpu_bar)
        sys_group.addWidget(self.ram_bar)
        layout.addLayout(sys_group)
        
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setStyleSheet("background-color: #222;")
        layout.addWidget(line)
        
        # GPU Grid
        gpu_title = QLabel("QUAD-GPU CLUSTER")
        gpu_title.setFont(QFont("Segoe UI", 8, QFont.Bold))
        gpu_title.setStyleSheet("color: #555;")
        layout.addWidget(gpu_title)
        
        self.gpu_cards = []
        for i in range(4):
            card = GPUCard(i)
            self.gpu_cards.append(card)
            layout.addWidget(card)
            
        layout.addStretch()
        
        # Status footer
        self.status = QLabel("● SYSTEM ACTIVE")
        self.status.setFont(QFont("Segoe UI", 7, QFont.Bold))
        self.status.setStyleSheet("color: #4CAF50;")
        self.status.setAlignment(Qt.AlignCenter)
        layout.addWidget(self.status)

    def update_all(self, metrics):
        """Update all UI elements with new metrics dict"""
        self.cpu_bar.update_value(metrics.get('cpu_util', 0))
        self.ram_bar.update_value(metrics.get('ram_util', 0))
        
        gpu_data = metrics.get('gpus', [])
        for i, card in enumerate(self.gpu_cards):
            if i < len(gpu_data):
                d = gpu_data[i]
                card.update_metrics(d.get('util', 0), d.get('temp', 0))
