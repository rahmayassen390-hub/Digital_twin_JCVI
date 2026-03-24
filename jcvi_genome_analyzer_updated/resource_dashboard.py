"""
Resource Dashboard — Real-Time Cellular Resource Competition Visualizer
=======================================================================

Embeddable QWidget showing live ribosome & polymerase allocation,
saturation gauges, top consumers, and competition dynamics.

Derived from JCVI-syn3.0 genome content:
  - Ribosomal protein genes → ribosome pool (~187)
  - RNA polymerase genes (rpoA/B/C/D) → polymerase pool (~50)
"""

import math
from collections import deque

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QGroupBox,
    QTableWidget, QTableWidgetItem, QFrame, QHeaderView,
    QSizePolicy, QPushButton, QDialog
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont, QColor

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


# ─────────────────────────────────────────────────────────
#  SATURATION GAUGE WIDGET
# ─────────────────────────────────────────────────────────
class SaturationGauge(QWidget):
    """Horizontal bar showing resource saturation with color transitions."""

    def __init__(self, title: str, total: int, color_low: str = "#4CAF50",
                 color_mid: str = "#FF9800", color_high: str = "#F44336", parent=None):
        super().__init__(parent)
        self._total = total
        self._used = 0
        self._color_low = color_low
        self._color_mid = color_mid
        self._color_high = color_high

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)

        # Header
        self.header = QLabel(f"{title}: 0 / {total}")
        self.header.setFont(QFont("Arial", 10, QFont.Bold))
        self.header.setStyleSheet("color: white;")
        layout.addWidget(self.header)

        # Bar
        self.bar = QLabel()
        self.bar.setMinimumHeight(22)
        self.bar.setStyleSheet("background: #333; border-radius: 4px;")
        layout.addWidget(self.bar)

        # Status text
        self.status = QLabel("🟢 Adequate capacity")
        self.status.setFont(QFont("Arial", 8))
        self.status.setStyleSheet("color: #4CAF50;")
        layout.addWidget(self.status)

    def set_values(self, used: int, total: int = None):
        if total is not None:
            self._total = total
        self._used = used
        pct = (used / self._total * 100) if self._total > 0 else 0
        pct = min(pct, 100)

        # Color
        if pct < 60:
            bar_color = self._color_low
            status_icon = "✅"
            status_text = "Adequate capacity"
            status_color = self._color_low
        elif pct < 80:
            bar_color = self._color_mid
            status_icon = "⚠️"
            status_text = "MODERATE load"
            status_color = self._color_mid
        else:
            bar_color = self._color_high
            status_icon = "🔴"
            status_text = "HIGH — approaching saturation"
            status_color = self._color_high

        self.header.setText(f"{used} / {self._total}  ({pct:.1f}%)")
        self.header.setStyleSheet(f"color: {bar_color};")

        # Gradient bar
        stop = max(0.01, min(pct / 100, 0.99))
        self.bar.setStyleSheet(
            f"background: qlineargradient(x1:0,y1:0,x2:1,y2:0,"
            f"stop:0 {bar_color}, stop:{stop:.3f} {bar_color},"
            f"stop:{min(stop+0.01, 1.0):.3f} #222, stop:1 #222);"
            f"border-radius: 4px; border: 1px solid #444;"
        )

        self.status.setText(f"  {status_icon} {status_text}")
        self.status.setStyleSheet(f"color: {status_color};")

        return pct


# ─────────────────────────────────────────────────────────
#  RESOURCE DASHBOARD
# ─────────────────────────────────────────────────────────
class ResourceDashboard(QDialog):
    """
    Standalone dialog showing live ribosome and polymerase competition.
    Opened via a button in the Live Transcription tab.
    """

    TOTAL_RIBOSOMES = 187   # Genome-derived estimate for JCVI-syn3.0
    TOTAL_POLYMERASES = 50   # Genome-derived estimate

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("⚡ Cellular Resources — Live Competition")
        self.setMinimumSize(700, 600)
        self.resize(740, 650)
        self.setStyleSheet("background: #0f0f23; color: #e0e0e0;")

        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(6, 4, 6, 4)
        main_layout.setSpacing(4)

        # Title
        title = QLabel("⚡ Cellular Resources — Live Competition")
        title.setFont(QFont("Arial", 11, QFont.Bold))
        title.setStyleSheet("color: #FFD54F; padding: 2px;")
        main_layout.addWidget(title)

        # Two-column layout for gauges
        gauges_layout = QHBoxLayout()

        # === LEFT: Ribosomes ===
        ribo_group = QGroupBox("🔩 RIBOSOMES (Translation)")
        ribo_group.setStyleSheet(
            "QGroupBox { color: #4FC3F7; border: 1px solid #333; border-radius: 4px; margin-top: 8px; padding-top: 16px; }"
            "QGroupBox::title { subcontrol-position: top left; padding: 0 8px; }"
        )
        ribo_layout = QVBoxLayout(ribo_group)

        self.ribo_gauge = SaturationGauge("Ribosomes", self.TOTAL_RIBOSOMES)
        ribo_layout.addWidget(self.ribo_gauge)

        # Top consumers
        self.ribo_table = QTableWidget(5, 3)
        self.ribo_table.setHorizontalHeaderLabels(["Gene", "Ribosomes", "%"])
        self.ribo_table.setStyleSheet(
            "QTableWidget { background: #16213e; color: #ccc; gridline-color: #333; border: none; }"
            "QHeaderView::section { background: #1a1a2e; color: #888; border: 1px solid #333; padding: 2px; }"
        )
        self.ribo_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.ribo_table.verticalHeader().setVisible(False)
        self.ribo_table.setMaximumHeight(140)
        ribo_layout.addWidget(self.ribo_table)

        gauges_layout.addWidget(ribo_group)

        # === RIGHT: Polymerases ===
        poly_group = QGroupBox("🧪 POLYMERASES (Transcription)")
        poly_group.setStyleSheet(
            "QGroupBox { color: #81C784; border: 1px solid #333; border-radius: 4px; margin-top: 8px; padding-top: 16px; }"
            "QGroupBox::title { subcontrol-position: top left; padding: 0 8px; }"
        )
        poly_layout = QVBoxLayout(poly_group)

        self.poly_gauge = SaturationGauge("Polymerases", self.TOTAL_POLYMERASES,
                                          "#81C784", "#FF9800", "#F44336")
        poly_layout.addWidget(self.poly_gauge)

        self.poly_table = QTableWidget(5, 3)
        self.poly_table.setHorizontalHeaderLabels(["Gene", "Polymerases", "%"])
        self.poly_table.setStyleSheet(
            "QTableWidget { background: #16213e; color: #ccc; gridline-color: #333; border: none; }"
            "QHeaderView::section { background: #1a1a2e; color: #888; border: 1px solid #333; padding: 2px; }"
        )
        self.poly_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.poly_table.verticalHeader().setVisible(False)
        self.poly_table.setMaximumHeight(140)
        poly_layout.addWidget(self.poly_table)

        gauges_layout.addWidget(poly_group)
        main_layout.addLayout(gauges_layout)

        # Competition dynamics summary
        dynamics_frame = QFrame()
        dynamics_frame.setStyleSheet(
            "background: #16213e; border: 1px solid #333; border-radius: 4px; padding: 6px;"
        )
        dyn_layout = QVBoxLayout(dynamics_frame)
        dyn_layout.setSpacing(2)

        dyn_header = QLabel("📊 Competition Dynamics")
        dyn_header.setFont(QFont("Arial", 9, QFont.Bold))
        dyn_header.setStyleSheet("color: #FFD54F;")
        dyn_layout.addWidget(dyn_header)

        self.limiting_label = QLabel("Limiting Factor: —")
        self.limiting_label.setFont(QFont("Consolas", 9))
        self.limiting_label.setStyleSheet("color: #B0BEC5;")
        dyn_layout.addWidget(self.limiting_label)

        self.bottleneck_label = QLabel("Bottleneck Impact: —")
        self.bottleneck_label.setFont(QFont("Consolas", 9))
        self.bottleneck_label.setStyleSheet("color: #B0BEC5;")
        dyn_layout.addWidget(self.bottleneck_label)

        self.competitive_label = QLabel("Most Competitive: —")
        self.competitive_label.setFont(QFont("Consolas", 9))
        self.competitive_label.setStyleSheet("color: #B0BEC5;")
        dyn_layout.addWidget(self.competitive_label)

        main_layout.addWidget(dynamics_frame)

        # Saturation history plot
        self.fig = Figure(figsize=(5, 1.5), dpi=100)
        self.fig.patch.set_facecolor('#1a1a2e')
        self.ax = self.fig.add_subplot(111)
        self.ax.set_facecolor('#16213e')
        self.ax.set_title("Saturation Over Time", fontsize=8, color='white')
        self.ax.set_ylabel("%", fontsize=7, color='#888')
        self.ax.tick_params(colors='#666', labelsize=6)
        for spine in self.ax.spines.values():
            spine.set_color('#333')
        self.ribo_line, = self.ax.plot([], [], color='#4FC3F7', linewidth=1.2, label='Ribosomes')
        self.poly_line, = self.ax.plot([], [], color='#81C784', linewidth=1.2, label='Polymerases')
        self.ax.legend(loc='upper left', fontsize=6, frameon=False, labelcolor='#888')
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setMaximumHeight(110)
        self.fig.tight_layout(pad=1.0)
        main_layout.addWidget(self.canvas)

        # History
        self.times = deque(maxlen=200)
        self.ribo_history = deque(maxlen=200)
        self.poly_history = deque(maxlen=200)
        self._frame_count = 0

    def update_from_simulation(self, gene_data_list: list, sim_time: float):
        """Called every frame. Computes resource allocation from gene data."""
        if not self.isVisible():
            return

        self._frame_count += 1
        if self._frame_count % 3 != 0:  # Update every 3rd frame
            return

        # Compute per-gene resource allocation
        gene_ribo = []
        gene_poly = []
        total_ribo_used = 0
        total_poly_used = 0

        for d in gene_data_list:
            gid = d.get('gene_id', '?')
            gname = d.get('gene_name', gid)
            mrna = d.get('mrna_count', 0.0)
            protein = d.get('protein_count', 0.0)
            tx_rate = d.get('transcription_rate', 0.0)

            # Ribosome allocation: proportional to mRNA count (simplified model)
            ribo = max(0, int(mrna * 0.02))  # ~2% of mRNA count as ribosomes
            gene_ribo.append((gid, gname, ribo))
            total_ribo_used += ribo

            # Polymerase allocation: proportional to transcription rate
            poly = max(0, tx_rate * 200)  # scaled
            gene_poly.append((gid, gname, poly))
            total_poly_used += poly

        # Cap at totals
        total_ribo_used = min(total_ribo_used, self.TOTAL_RIBOSOMES)
        total_poly_used = min(total_poly_used, self.TOTAL_POLYMERASES)

        # Update gauges
        ribo_pct = self.ribo_gauge.set_values(int(total_ribo_used))
        poly_pct = self.poly_gauge.set_values(int(total_poly_used))

        # Sort and show top 5 consumers
        gene_ribo.sort(key=lambda x: x[2], reverse=True)
        gene_poly.sort(key=lambda x: x[2], reverse=True)

        self._fill_table(self.ribo_table, gene_ribo[:5], total_ribo_used)
        self._fill_table(self.poly_table, gene_poly[:5], total_poly_used)

        # Competition dynamics
        if ribo_pct > poly_pct:
            self.limiting_label.setText(f"🔩 Limiting Factor: RIBOSOMES  ({ribo_pct:.0f}% sat)")
            self.limiting_label.setStyleSheet("color: #F44336;" if ribo_pct > 80 else "color: #FF9800;")
            impact = min(50, max(0, (ribo_pct - 50) * 0.8))
            self.bottleneck_label.setText(f"📉 Bottleneck Impact: -{impact:.0f}% translation rate")
        else:
            self.limiting_label.setText(f"🧪 Limiting Factor: POLYMERASES  ({poly_pct:.0f}% sat)")
            self.limiting_label.setStyleSheet("color: #F44336;" if poly_pct > 80 else "color: #FF9800;")
            impact = min(50, max(0, (poly_pct - 50) * 0.8))
            self.bottleneck_label.setText(f"📉 Bottleneck Impact: -{impact:.0f}% transcription rate")

        if gene_ribo:
            top3 = ", ".join([f"{g[1] or g[0]}" for g in gene_ribo[:3]])
            self.competitive_label.setText(f"🏆 Most Competitive: {top3}")

        # History
        self.times.append(sim_time)
        self.ribo_history.append(ribo_pct)
        self.poly_history.append(poly_pct)

        # Refresh plot
        if len(self.times) > 2:
            t = list(self.times)
            self.ribo_line.set_data(t, list(self.ribo_history))
            self.poly_line.set_data(t, list(self.poly_history))
            self.ax.set_xlim(t[0], t[-1])
            self.ax.set_ylim(0, 105)
            self.canvas.draw_idle()

    def _fill_table(self, table: QTableWidget, entries: list, total: float):
        table.setRowCount(len(entries))
        for i, (gid, gname, count) in enumerate(entries):
            label = gname[:12] if gname else gid[:12]
            pct = (count / total * 100) if total > 0 else 0

            item_name = QTableWidgetItem(label)
            item_name.setForeground(QColor("#ccc"))
            item_count = QTableWidgetItem(f"{count:.0f}")
            item_count.setForeground(QColor("#4FC3F7"))
            item_pct = QTableWidgetItem(f"{pct:.1f}%")
            pct_color = "#4CAF50" if pct < 10 else "#FF9800" if pct < 20 else "#F44336"
            item_pct.setForeground(QColor(pct_color))

            table.setItem(i, 0, item_name)
            table.setItem(i, 1, item_count)
            table.setItem(i, 2, item_pct)
