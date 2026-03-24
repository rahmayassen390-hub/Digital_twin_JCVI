"""
Gene Dynamics Panel — Real-Time Visualization of All Gene-Level Dynamics
=========================================================================

4-tab QDialog showing mRNA Transcription, Protein Translation,
mRNA Degradation, and Protein Degradation monitors for a selected gene.

Each tab provides:
  - Live time-series plot (matplotlib)
  - Parameter table with live values
  - Dynamic equation with substituted values
  - Color-coded status indicator
  - Predicted steady-state marker
"""

import math
import numpy as np
from collections import deque

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QTabWidget, QWidget,
    QLabel, QGroupBox, QGridLayout, QFrame, QScrollArea,
    QPushButton, QComboBox, QSizePolicy
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont, QColor

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt


# ─────────────────────────────────────────────────────────
#  PARAMETER ROW WIDGET
# ─────────────────────────────────────────────────────────
class ParamRow(QLabel):
    """Single parameter display row with icon, name, and live value."""

    def __init__(self, icon: str, name: str, initial: str = "—", parent=None):
        super().__init__(parent)
        self._icon = icon
        self._name = name
        self.set_value(initial)
        self.setFont(QFont("Consolas", 9))
        self.setStyleSheet("color: #B0BEC5; padding: 1px 4px;")

    def set_value(self, value: str, color: str = "#B0BEC5"):
        self.setText(f"{self._icon} {self._name}: {value}")
        self.setStyleSheet(f"color: {color}; padding: 1px 4px;")


# ─────────────────────────────────────────────────────────
#  MINI PLOT WIDGET (Shared Component)
# ─────────────────────────────────────────────────────────
class DynamicsPlotWidget(QWidget):
    """Small matplotlib time-series plot for a single metric."""

    MAX_HISTORY = 400  # data points kept

    def __init__(self, title: str, ylabel: str, color: str = "#4FC3F7", parent=None):
        super().__init__(parent)
        self._title = title
        self._color = color

        self.times = deque(maxlen=self.MAX_HISTORY)
        self.values = deque(maxlen=self.MAX_HISTORY)
        self.steady_state = None

        # Dark figure
        self.fig = Figure(figsize=(4.5, 2.0), dpi=100)
        self.fig.patch.set_facecolor('#1a1a2e')
        self.ax = self.fig.add_subplot(111)
        self.ax.set_facecolor('#16213e')
        self.ax.set_title(title, fontsize=9, color='white', fontweight='bold')
        self.ax.set_xlabel("Time (s)", fontsize=7, color='#888')
        self.ax.set_ylabel(ylabel, fontsize=7, color='#888')
        self.ax.tick_params(colors='#666', labelsize=7)
        for spine in self.ax.spines.values():
            spine.set_color('#333')

        self.line, = self.ax.plot([], [], color=color, linewidth=1.5)
        self.ss_line = None

        self.canvas = FigureCanvas(self.fig)
        self.canvas.setMinimumHeight(120)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.canvas)
        self.fig.tight_layout(pad=1.5)

    def add_point(self, t: float, v: float):
        self.times.append(t)
        self.values.append(v)

    def set_steady_state(self, ss: float):
        self.steady_state = ss

    def refresh(self):
        if len(self.times) < 2:
            return
        t = list(self.times)
        v = list(self.values)
        self.line.set_data(t, v)
        self.ax.set_xlim(t[0], t[-1])
        vmin = min(v)
        vmax = max(v)
        margin = max((vmax - vmin) * 0.15, 1.0)
        self.ax.set_ylim(max(0, vmin - margin), vmax + margin)

        # Steady-state line
        if self.ss_line:
            self.ss_line.remove()
            self.ss_line = None
        if self.steady_state is not None and self.steady_state > 0:
            self.ss_line = self.ax.axhline(
                self.steady_state, color='#FFD54F', linestyle='--',
                linewidth=0.8, alpha=0.7
            )

        self.canvas.draw_idle()

    def clear(self):
        self.times.clear()
        self.values.clear()
        self.steady_state = None


# ─────────────────────────────────────────────────────────
#  DECAY PROJECTION PLOT
# ─────────────────────────────────────────────────────────
class DecayCurvePlot(QWidget):
    """Shows projected exponential decay curve from current count."""

    def __init__(self, title: str, time_unit: str = "min", color: str = "#FF7043", parent=None):
        super().__init__(parent)
        self.fig = Figure(figsize=(4.5, 2.0), dpi=100)
        self.fig.patch.set_facecolor('#1a1a2e')
        self.ax = self.fig.add_subplot(111)
        self.ax.set_facecolor('#16213e')
        self.ax.set_title(title, fontsize=9, color='white', fontweight='bold')
        self.ax.set_xlabel(f"Time ({time_unit})", fontsize=7, color='#888')
        self.ax.set_ylabel("Molecules", fontsize=7, color='#888')
        self.ax.tick_params(colors='#666', labelsize=7)
        for spine in self.ax.spines.values():
            spine.set_color('#333')

        self._color = color
        self._time_unit = time_unit
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setMinimumHeight(120)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.canvas)
        self.fig.tight_layout(pad=1.5)

    def update_curve(self, current_count: float, half_life_seconds: float):
        """Draw projected decay from current count."""
        self.ax.clear()
        self.ax.set_facecolor('#16213e')
        self.ax.tick_params(colors='#666', labelsize=7)
        for spine in self.ax.spines.values():
            spine.set_color('#333')

        if current_count <= 0 or half_life_seconds <= 0:
            self.ax.text(0.5, 0.5, "No protein present\nor zero half-life",
                        ha='center', va='center', transform=self.ax.transAxes,
                        color='gray', fontsize=10)
            self.ax.set_title("Decay Projection (No Data)", fontsize=9, color='white', fontweight='bold')
            self.ax.set_xlabel("") # Clear labels for no data
            self.ax.set_ylabel("")
            self.canvas.draw_idle()
            return

        divisor = 60.0 if self._time_unit == "min" else 3600.0
        k = math.log(2) / half_life_seconds
        hl_display = half_life_seconds / divisor

        # Generate decay curve (0 to 5 half-lives)
        t_max = 5 * hl_display
        t = np.linspace(0, t_max, 100)
        counts = current_count * np.exp(-k * t * divisor)

        self.ax.plot(t, counts, color=self._color, linewidth=2)
        self.ax.fill_between(t, 0, counts, color=self._color, alpha=0.15)

        # Mark half-life points
        for i in range(1, 4):
            t_hl = i * hl_display
            c_hl = current_count * (0.5 ** i)
            self.ax.plot(t_hl, c_hl, 'o', color='#FFD54F', markersize=5, zorder=5)
            pct = 100 * (0.5 ** i)
            self.ax.annotate(
                f"{pct:.0f}%", (t_hl, c_hl), textcoords="offset points",
                xytext=(8, 5), fontsize=7, color='#FFD54F'
            )

        self.ax.set_title(
            f"Decay Projection  (t½ = {hl_display:.1f} {self._time_unit})",
            fontsize=9, color='white', fontweight='bold'
        )
        self.ax.set_xlabel(f"Time ({self._time_unit})", fontsize=7, color='#888')
        self.ax.set_ylabel("Molecules", fontsize=7, color='#888')
        self.fig.tight_layout(pad=1.5)
        self.canvas.draw_idle()


# ─────────────────────────────────────────────────────────
#  STATUS BADGE
# ─────────────────────────────────────────────────────────
class StatusBadge(QLabel):
    COLORS = {
        "rising": ("#4CAF50", "🟢"),
        "steady": ("#2196F3", "🔵"),
        "declining": ("#FF9800", "🟠"),
        "depleted": ("#F44336", "🔴"),
        "paused": ("#9E9E9E", "⏸️"),
    }

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFont(QFont("Arial", 10, QFont.Bold))
        self.set_status("rising")

    def set_status(self, status: str):
        color, icon = self.COLORS.get(status, ("#9E9E9E", "❓"))
        label_map = {
            "rising": "Approaching Steady State",
            "steady": "At Steady State",
            "steady_state": "At Steady State",
            "declining": "Declining",
            "declining_fast": "Declining (Fast)",
            "depleted": "Depleted",
            "paused": "Transcription Paused",
        }
        text = label_map.get(status, status.replace("_", " ").title())
        self.setText(f"  {icon} {text}")
        self.setStyleSheet(
            f"color: {color}; background: #1a1a2e; border-radius: 4px; padding: 3px 8px;"
        )


# ═══════════════════════════════════════════════════════════
#  TAB 1: mRNA TRANSCRIPTION MONITOR
# ═══════════════════════════════════════════════════════════
class MRNATranscriptionTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setStyleSheet("background: #0f0f23; color: #e0e0e0;")
        layout = QVBoxLayout(self)
        layout.setSpacing(4)

        # Header
        self.gene_label = QLabel("Select a gene to monitor")
        self.gene_label.setFont(QFont("Arial", 11, QFont.Bold))
        self.gene_label.setStyleSheet("color: #4FC3F7;")
        layout.addWidget(self.gene_label)

        # Status badge
        self.status = StatusBadge()
        layout.addWidget(self.status)

        # Plot
        self.plot = DynamicsPlotWidget("mRNA Count Over Time", "Molecules", "#4FC3F7")
        layout.addWidget(self.plot)

        # Parameters group
        params = QGroupBox("Live Parameters")
        params.setStyleSheet(
            "QGroupBox { color: #4FC3F7; border: 1px solid #333; border-radius: 4px; margin-top: 6px; padding-top: 14px; }"
            "QGroupBox::title { subcontrol-position: top left; padding: 0 6px; }"
        )
        pg = QVBoxLayout(params)
        pg.setSpacing(1)

        self.p_promoter_type = ParamRow("🎯", "Promoter Type")
        self.p_tx_rate = ParamRow("⚡", "Transcription Rate")
        self.p_tx_level = ParamRow("📊", "TX Level (BED)")
        self.p_regulation = ParamRow("🔗", "Regulation")
        self.p_polymerase = ParamRow("🧪", "RNA Polymerase Alloc")
        self.p_knockout = ParamRow("🔧", "Knockout Modifier")
        self.p_current_rate = ParamRow("📈", "Current Net Rate")
        self.p_steady_state = ParamRow("🎯", "Predicted Steady State")

        for p in [self.p_promoter_type, self.p_tx_rate, self.p_tx_level,
                  self.p_regulation, self.p_polymerase, self.p_knockout,
                  self.p_current_rate, self.p_steady_state]:
            pg.addWidget(p)

        layout.addWidget(params)

        # Equation
        self.equation_label = QLabel("dm/dt = TX_rate × modifier − k_deg × m")
        self.equation_label.setFont(QFont("Consolas", 9))
        self.equation_label.setStyleSheet(
            "color: #FFD54F; background: #1a1a2e; padding: 6px; border: 1px solid #333; border-radius: 4px;"
        )
        self.equation_label.setAlignment(Qt.AlignCenter)
        self.equation_label.setWordWrap(True)
        layout.addWidget(self.equation_label)

    def update_data(self, d: dict, sim_time: float):
        gene_id = d.get('gene_id', '?')
        gene_name = d.get('gene_name', '')
        self.gene_label.setText(f"Gene: {gene_name} ({gene_id})")

        mrna = d.get('mrna_count', 0.0)
        tx_rate = d.get('transcription_rate', 0.0)
        status = d.get('status', 'rising')
        regulation = d.get('regulation_type', 'Unknown')
        phase = d.get('phase', '')
        ko_mod = d.get('knockout_modifier', 1.0)
        mrna_hl = d.get('mrna_half_life', 300.0)
        k_deg = math.log(2) / mrna_hl if mrna_hl > 0 else 0.01

        # Compute steady state
        effective_rate = tx_rate * ko_mod
        ss = effective_rate / k_deg if k_deg > 0 else 0
        net_rate = effective_rate - k_deg * mrna

        # Update plot
        self.plot.add_point(sim_time, mrna)
        self.plot.set_steady_state(ss)

        # Update status
        st = status.replace("_fast", "")
        if d.get('is_paused', False):
            st = "paused"
        self.status.set_status(st)

        # Update params
        self.p_promoter_type.set_value(d.get('has_own_promoter', 'Yes'))
        self.p_tx_rate.set_value(f"{tx_rate:.4f} mol/s", "#4CAF50")
        level_pct = d.get('transcription_level', 0.0)
        level_stars = "★" * max(1, int(level_pct / 33) + 1)
        self.p_tx_level.set_value(f"{level_pct:.0f}% ({level_stars})")
        self.p_regulation.set_value(regulation)
        self.p_polymerase.set_value(f"{max(1, mrna * 0.003):.1f} polymerases")
        ko_color = "#4CAF50" if ko_mod >= 0.9 else "#FF9800" if ko_mod > 0.1 else "#F44336"
        self.p_knockout.set_value(f"{ko_mod:.2f}", ko_color)
        rate_color = "#4CAF50" if net_rate > 0.001 else "#2196F3" if abs(net_rate) < 0.001 else "#FF9800"
        self.p_current_rate.set_value(f"{net_rate:+.4f} mol/s", rate_color)
        self.p_steady_state.set_value(f"{ss:.0f} molecules", "#FFD54F")

        # Equation with live values
        self.equation_label.setText(
            f"dm/dt = {tx_rate:.4f} × {ko_mod:.2f} − {k_deg:.5f} × {mrna:.0f} = {net_rate:+.4f} mol/s"
        )


# ═══════════════════════════════════════════════════════════
#  TAB 2: PROTEIN TRANSLATION MONITOR
# ═══════════════════════════════════════════════════════════
class ProteinTranslationTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setStyleSheet("background: #0f0f23; color: #e0e0e0;")
        layout = QVBoxLayout(self)
        layout.setSpacing(4)

        self.gene_label = QLabel("Select a gene to monitor")
        self.gene_label.setFont(QFont("Arial", 11, QFont.Bold))
        self.gene_label.setStyleSheet("color: #81C784;")
        layout.addWidget(self.gene_label)

        self.status = StatusBadge()
        layout.addWidget(self.status)

        self.plot = DynamicsPlotWidget("Protein Count Over Time", "Molecules", "#81C784")
        layout.addWidget(self.plot)

        # Ribosome allocation bar
        ribo_frame = QFrame()
        ribo_frame.setStyleSheet("background: #16213e; border: 1px solid #333; border-radius: 4px; padding: 4px;")
        ribo_layout = QVBoxLayout(ribo_frame)
        ribo_layout.setContentsMargins(6, 4, 6, 4)
        self.ribo_label = QLabel("Ribosome Allocation: — / —")
        self.ribo_label.setFont(QFont("Consolas", 9))
        self.ribo_label.setStyleSheet("color: #81C784;")
        ribo_layout.addWidget(self.ribo_label)
        self.ribo_bar = QLabel()
        self.ribo_bar.setMinimumHeight(16)
        self.ribo_bar.setStyleSheet("background: #333; border-radius: 3px;")
        ribo_layout.addWidget(self.ribo_bar)
        layout.addWidget(ribo_frame)

        # Parameters
        params = QGroupBox("Live Parameters")
        params.setStyleSheet(
            "QGroupBox { color: #81C784; border: 1px solid #333; border-radius: 4px; margin-top: 6px; padding-top: 14px; }"
            "QGroupBox::title { subcontrol-position: top left; padding: 0 6px; }"
        )
        pg = QVBoxLayout(params)
        pg.setSpacing(1)

        self.p_rbs = ParamRow("🧬", "RBS Strength")
        self.p_tl_rate = ParamRow("⚡", "Translation Rate")
        self.p_tl_level = ParamRow("📊", "TL Level (BED)")
        self.p_mrna_template = ParamRow("📋", "mRNA Template Count")
        self.p_ribosome_alloc = ParamRow("🔩", "Ribosome Allocation")
        self.p_protein_len = ParamRow("📏", "Protein Length")
        self.p_steady_state = ParamRow("🎯", "Predicted Steady State")

        for p in [self.p_rbs, self.p_tl_rate, self.p_tl_level, self.p_mrna_template,
                  self.p_ribosome_alloc, self.p_protein_len, self.p_steady_state]:
            pg.addWidget(p)

        layout.addWidget(params)

        self.equation_label = QLabel("dp/dt = TL_rate × mRNA − k_deg_p × P")
        self.equation_label.setFont(QFont("Consolas", 9))
        self.equation_label.setStyleSheet(
            "color: #FFD54F; background: #1a1a2e; padding: 6px; border: 1px solid #333; border-radius: 4px;"
        )
        self.equation_label.setAlignment(Qt.AlignCenter)
        self.equation_label.setWordWrap(True)
        layout.addWidget(self.equation_label)

    def update_data(self, d: dict, sim_time: float):
        gene_id = d.get('gene_id', '?')
        gene_name = d.get('gene_name', '')
        self.gene_label.setText(f"Gene: {gene_name} ({gene_id})")

        protein = d.get('protein_count', 0.0)
        mrna = d.get('mrna_count', 0.0)
        tl_rate = d.get('translation_rate', 0.1)
        status = d.get('status', 'rising')
        prot_hl = d.get('protein_half_life', 90000.0)
        k_deg_p = math.log(2) / prot_hl if prot_hl > 0 else 1e-6
        rbs = d.get('rbs_strength', 1.0)

        effective_tl = tl_rate * mrna
        ss = effective_tl / k_deg_p if k_deg_p > 0 else 0
        net_rate = effective_tl - k_deg_p * protein

        # Plot update
        self.plot.add_point(sim_time, protein)
        self.plot.set_steady_state(ss)

        # Status
        st = status.replace("_fast", "")
        if d.get('is_paused', False):
            st = "paused"
        self.status.set_status(st)

        # Ribosome bar
        ribo_per_gene = max(1, int(mrna * 0.018))
        total_ribo = 187
        pct = min(100, ribo_per_gene / total_ribo * 100)
        bar_width = int(pct * 3)
        bar_color = "#4CAF50" if pct < 5 else "#FF9800" if pct < 15 else "#F44336"
        self.ribo_label.setText(f"Ribosome Allocation: {ribo_per_gene} / {total_ribo} ({pct:.1f}%)")
        self.ribo_bar.setStyleSheet(
            f"background: qlineargradient(x1:0,y1:0,x2:1,y2:0,"
            f"stop:0 {bar_color}, stop:{min(pct/100, 0.99):.3f} {bar_color},"
            f"stop:{min(pct/100 + 0.01, 1.0):.3f} #333, stop:1 #333);"
            f"border-radius: 3px;"
        )

        # Parameters
        rbs_label = "Strong" if rbs > 0.7 else "Medium" if rbs > 0.3 else "Weak"
        self.p_rbs.set_value(f"{rbs:.2f} ({rbs_label})")
        self.p_tl_rate.set_value(f"{tl_rate:.4f} prot/mRNA/s", "#81C784")
        tl_level = d.get('protein_level', 0.0)
        self.p_tl_level.set_value(f"{tl_level:.0f}%")
        self.p_mrna_template.set_value(f"{mrna:.0f} molecules", "#4FC3F7")
        self.p_ribosome_alloc.set_value(f"{ribo_per_gene} ribosomes/gene")
        gene_len = d.get('gene_length', 0)
        aa_len = gene_len // 3 if gene_len > 0 else 0
        time_per = aa_len / 19 if aa_len > 0 else 0  # 19 aa/s elongation
        self.p_protein_len.set_value(f"{aa_len} aa → {time_per:.1f}s per protein")
        self.p_steady_state.set_value(f"{ss:,.0f} molecules", "#FFD54F")

        self.equation_label.setText(
            f"dp/dt = {tl_rate:.4f} × {mrna:.0f} − {k_deg_p:.7f} × {protein:.0f} = {net_rate:+.2f} mol/s"
        )


# ═══════════════════════════════════════════════════════════
#  TAB 3: mRNA DEGRADATION MONITOR
# ═══════════════════════════════════════════════════════════
class MRNADegradationTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setStyleSheet("background: #0f0f23; color: #e0e0e0;")
        layout = QVBoxLayout(self)
        layout.setSpacing(4)

        self.gene_label = QLabel("Select a gene to monitor")
        self.gene_label.setFont(QFont("Arial", 11, QFont.Bold))
        self.gene_label.setStyleSheet("color: #FF7043;")
        layout.addWidget(self.gene_label)

        # ── LIVE time-series plot (actual mRNA count over time) ──
        self.plot = DynamicsPlotWidget(
            "mRNA Count (Live)", "Molecules", "#FF7043"
        )
        layout.addWidget(self.plot)

        # Decay curve plot (theoretical projection)
        self.decay_plot = DecayCurvePlot(
            "Decay Projection (if TX stopped now)", "min", "#FF7043"
        )
        layout.addWidget(self.decay_plot)

        # Explanation label for the decay projection
        decay_info = QLabel(
            "📖 This curve predicts how fast mRNA molecules would degrade "
            "if transcription were stopped at this moment. It follows "
            "first-order kinetics: N(t) = N₀ × e⁻ᵏᵗ, where k = ln(2)/t½. "
            "The 50%→25%→12% markers show molecule count at 1×, 2×, and 3× "
            "half-lives. Unlike protein, mRNA degrades within minutes, "
            "meaning gene expression responds rapidly to transcription changes."
        )
        decay_info.setWordWrap(True)
        decay_info.setFont(QFont("Arial", 8))
        decay_info.setStyleSheet(
            "color: #8899AA; background: #12152a; border: 1px solid #2a2a4a; "
            "border-radius: 4px; padding: 6px; margin: 2px 0;"
        )
        layout.addWidget(decay_info)

        # Parameters
        params = QGroupBox("Degradation Parameters")
        params.setStyleSheet(
            "QGroupBox { color: #FF7043; border: 1px solid #333; border-radius: 4px; margin-top: 6px; padding-top: 14px; }"
            "QGroupBox::title { subcontrol-position: top left; padding: 0 6px; }"
        )
        pg = QVBoxLayout(params)
        pg.setSpacing(1)

        self.p_half_life = ParamRow("⏱️", "Half-Life")
        self.p_source = ParamRow("📋", "Data Source")
        self.p_confidence = ParamRow("📊", "Confidence")
        self.p_k_deg = ParamRow("📉", "Degradation Constant")
        self.p_current_decay = ParamRow("💀", "Current Decay Rate")
        self.p_time_50 = ParamRow("⏳", "Time to 50%")
        self.p_time_10 = ParamRow("⏳", "Time to 10%")
        self.p_time_1 = ParamRow("⏳", "Time to 1%")

        for p in [self.p_half_life, self.p_source, self.p_confidence,
                  self.p_k_deg, self.p_current_decay,
                  self.p_time_50, self.p_time_10, self.p_time_1]:
            pg.addWidget(p)

        layout.addWidget(params)

        # Real-time degradation stats
        rt_frame = QFrame()
        rt_frame.setStyleSheet("background: #16213e; border: 1px solid #333; border-radius: 4px; padding: 6px;")
        rt_layout = QVBoxLayout(rt_frame)
        rt_layout.setSpacing(2)
        rt_header = QLabel("Real-Time Degradation")
        rt_header.setFont(QFont("Arial", 9, QFont.Bold))
        rt_header.setStyleSheet("color: #FF7043;")
        rt_layout.addWidget(rt_header)
        self.rt_per_sec = ParamRow("💥", "Molecules degraded/s")
        self.rt_cumulative = ParamRow("📊", "Cumulative degraded")
        self.rt_position = ParamRow("📍", "Half-life position")
        rt_layout.addWidget(self.rt_per_sec)
        rt_layout.addWidget(self.rt_cumulative)
        rt_layout.addWidget(self.rt_position)
        layout.addWidget(rt_frame)

        # Compare label
        self.compare_label = QLabel("")
        self.compare_label.setFont(QFont("Consolas", 8))
        self.compare_label.setStyleSheet("color: #888; padding: 4px;")
        self.compare_label.setWordWrap(True)
        layout.addWidget(self.compare_label)

        self.last_sync_label = QLabel("Last Sync: —")
        self.last_sync_label.setFont(QFont("Consolas", 8))
        self.last_sync_label.setStyleSheet("color: #4fc3f7; padding: 2px;")
        self.last_sync_label.setAlignment(Qt.AlignRight)
        layout.addWidget(self.last_sync_label)

        self._cumulative_deg = 0.0
        self._last_mrna = 0.0

    def update_data(self, d: dict, sim_time: float):
        gene_id = d.get('gene_id', '?')
        gene_name = d.get('gene_name', '')
        self.gene_label.setText(f"mRNA Degradation: {gene_name} ({gene_id})")

        mrna = d.get('mrna_count', 0.0)
        mrna_hl = d.get('mrna_half_life', 300.0)
        k_deg = math.log(2) / mrna_hl if mrna_hl > 0 else 0.01
        current_decay = k_deg * mrna
        hl_used_default = d.get('uses_default_params', True)

        # ── Update LIVE time-series plot ──
        self.plot.add_point(sim_time, mrna)
        tx_rate = d.get('transcription_rate', 0.01)
        ko = d.get('knockout_modifier', 1.0)
        ss_mrna = tx_rate * ko / k_deg if k_deg > 0 else 0
        self.plot.set_steady_state(ss_mrna)

        # Update decay projection
        self.decay_plot.update_curve(mrna, mrna_hl)

        # Parameters
        hl_min = mrna_hl / 60.0
        self.p_half_life.set_value(f"{hl_min:.1f} minutes ({mrna_hl:.0f} s)")

        if hl_used_default:
            self.p_source.set_value("ML-Predicted ⚠️ (Linear Regression)", "#FF9800")
            self.p_confidence.set_value("Medium (11 training genes)", "#FF9800")
        else:
            self.p_source.set_value("Measured ✅ (syn3A_degradation_rates.csv)", "#4CAF50")
            self.p_confidence.set_value("High (experimental)", "#4CAF50")

        self.p_k_deg.set_value(f"{k_deg:.5f} /s  [ln(2) / {mrna_hl:.0f}]")
        self.p_current_decay.set_value(f"{current_decay:.2f} mol/s", "#FF7043")

        t50 = mrna_hl / 60.0
        t10 = mrna_hl * 3.32 / 60.0  # log2(10) half-lives
        t1 = mrna_hl * 6.64 / 60.0   # log2(100) half-lives
        self.p_time_50.set_value(f"{t50:.1f} min")
        self.p_time_10.set_value(f"{t10:.1f} min")
        self.p_time_1.set_value(f"{t1:.1f} min")

        # Real-time stats
        deg_per_sec = k_deg * mrna
        self.rt_per_sec.set_value(f"{deg_per_sec:.1f}", "#FF7043")
        self._cumulative_deg += deg_per_sec * 5.0  # ~5s per frame
        self.rt_cumulative.set_value(f"{self._cumulative_deg:,.0f}")
        # Position in half-life cycle
        tx_rate = d.get('transcription_rate', 0.01)
        ko = d.get('knockout_modifier', 1.0)
        ss = tx_rate * ko / k_deg if k_deg > 0 else 1
        self.last_sync_label.setText(f"Last Sync: {sim_time:.1f}s")
        hl_pos = mrna / ss if ss > 0 else 0
        self.rt_position.set_value(f"{hl_pos:.2f} ({hl_pos*100:.0f}% of steady state)")

        # Compare
        self.compare_label.setText(
            f"Predicted (ML): {hl_min:.1f} min | "
            f"{'Measured available ✅' if not hl_used_default else 'No experimental data available ⚠️'}"
        )


# ═══════════════════════════════════════════════════════════
#  TAB 4: PROTEIN DEGRADATION MONITOR
# ═══════════════════════════════════════════════════════════
class ProteinDegradationTab(QWidget):
    N_END_RULE = {
        'M': ('Stabilizing', 'Long', '#4CAF50'), 'A': ('Stabilizing', 'Long', '#4CAF50'),
        'G': ('Stabilizing', 'Long', '#4CAF50'), 'P': ('Stabilizing', 'Long', '#4CAF50'),
        'S': ('Stabilizing', 'Long', '#4CAF50'), 'T': ('Stabilizing', 'Long', '#4CAF50'),
        'V': ('Stabilizing', 'Long', '#4CAF50'),
        'D': ('Secondary destabilizing', 'Medium', '#FF9800'),
        'E': ('Secondary destabilizing', 'Medium', '#FF9800'),
        'N': ('Tertiary destabilizing', 'Medium', '#FF9800'),
        'Q': ('Tertiary destabilizing', 'Medium', '#FF9800'),
        'C': ('Secondary destabilizing', 'Medium', '#FF9800'),
        'R': ('Primary destabilizing', 'Short', '#F44336'),
        'K': ('Primary destabilizing', 'Short', '#F44336'),
        'H': ('Primary destabilizing', 'Short', '#F44336'),
        'F': ('Primary destabilizing', 'Short', '#F44336'),
        'L': ('Primary destabilizing', 'Short', '#F44336'),
        'W': ('Primary destabilizing', 'Short', '#F44336'),
        'Y': ('Primary destabilizing', 'Short', '#F44336'),
        'I': ('Primary destabilizing', 'Short', '#F44336'),
    }

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setStyleSheet("background: #0f0f23; color: #e0e0e0;")
        layout = QVBoxLayout(self)
        layout.setSpacing(4)

        self.gene_label = QLabel("Select a gene to monitor")
        self.gene_label.setFont(QFont("Arial", 11, QFont.Bold))
        self.gene_label.setStyleSheet("color: #CE93D8;")
        layout.addWidget(self.gene_label)

        # ── LIVE time-series plot (actual protein count over time) ──
        self.plot = DynamicsPlotWidget(
            "Protein Count (Live)", "Molecules", "#CE93D8"
        )
        layout.addWidget(self.plot)

        # Decay curve (theoretical projection)
        self.decay_plot = DecayCurvePlot(
            "Protein Decay Projection", "hr", "#CE93D8"
        )
        layout.addWidget(self.decay_plot)

        # Explanation label for the decay projection
        decay_info = QLabel(
            "📖 This curve predicts how fast protein molecules would degrade "
            "if translation were stopped at this moment. Protein half-life is "
            "determined by the N-End Rule — the identity of the first amino acid "
            "controls ubiquitin-dependent degradation. Stabilizing residues (M, A, G…) "
            "give half-lives of 20–100 hr; destabilizing residues (R, K, F…) give "
            "0.5–5 hr. Both curves use the same exponential decay formula "
            "N(t) = N₀ × e⁻ᵏᵗ, but proteins persist far longer than mRNA (hours vs minutes)."
        )
        decay_info.setWordWrap(True)
        decay_info.setFont(QFont("Arial", 8))
        decay_info.setStyleSheet(
            "color: #8899AA; background: #12152a; border: 1px solid #2a2a4a; "
            "border-radius: 4px; padding: 6px; margin: 2px 0;"
        )
        layout.addWidget(decay_info)

        # Parameters
        params = QGroupBox("Degradation Parameters")
        params.setStyleSheet(
            "QGroupBox { color: #CE93D8; border: 1px solid #333; border-radius: 4px; margin-top: 6px; padding-top: 14px; }"
            "QGroupBox::title { subcontrol-position: top left; padding: 0 6px; }"
        )
        pg = QVBoxLayout(params)
        pg.setSpacing(1)

        self.p_half_life = ParamRow("⏱️", "Half-Life")
        self.p_source = ParamRow("📋", "Prediction Source")
        self.p_n_terminal = ParamRow("🧬", "N-Terminal Residue")
        self.p_stability = ParamRow("🛡️", "Stability Class")
        self.p_k_deg = ParamRow("📉", "Degradation Constant")
        self.p_current_decay = ParamRow("💀", "Current Decay Rate")
        self.p_time_50 = ParamRow("⏳", "Time to 50%")
        self.p_time_10 = ParamRow("⏳", "Time to 10%")

        for p in [self.p_half_life, self.p_source, self.p_n_terminal,
                  self.p_stability, self.p_k_deg, self.p_current_decay,
                  self.p_time_50, self.p_time_10]:
            pg.addWidget(p)

        layout.addWidget(params)

        # N-End Rule detail box
        self.nend_frame = QFrame()
        self.nend_frame.setStyleSheet("background: #16213e; border: 1px solid #333; border-radius: 4px; padding: 6px;")
        nend_layout = QVBoxLayout(self.nend_frame)
        nend_layout.setSpacing(2)
        nend_header = QLabel("🧬 N-End Rule Details")
        nend_header.setFont(QFont("Arial", 9, QFont.Bold))
        nend_header.setStyleSheet("color: #CE93D8;")
        nend_layout.addWidget(nend_header)
        self.nend_seq = ParamRow("🔤", "N-Terminal Sequence")
        self.nend_class = ParamRow("📋", "Rule Classification")
        self.nend_range = ParamRow("📊", "Expected HL Range")
        nend_layout.addWidget(self.nend_seq)
        nend_layout.addWidget(self.nend_class)
        nend_layout.addWidget(self.nend_range)
        layout.addWidget(self.nend_frame)

        # Real-time
        rt_frame = QFrame()
        rt_frame.setStyleSheet("background: #16213e; border: 1px solid #333; border-radius: 4px; padding: 6px;")
        rt_layout = QVBoxLayout(rt_frame)
        rt_layout.setSpacing(2)
        self.rt_stability = ParamRow("🛡️", "Protein is")
        self.rt_turnover = ParamRow("🔄", "Turnover rate")
        rt_layout.addWidget(self.rt_stability)
        rt_layout.addWidget(self.rt_turnover)
        layout.addWidget(rt_frame)

        self.last_sync_label = QLabel("Last Sync: —")
        self.last_sync_label.setFont(QFont("Consolas", 8))
        self.last_sync_label.setStyleSheet("color: #CE93D8; padding: 2px;")
        self.last_sync_label.setAlignment(Qt.AlignRight)
        layout.addWidget(self.last_sync_label)

    def update_data(self, d: dict, sim_time: float):
        gene_id = d.get('gene_id', '?')
        gene_name = d.get('gene_name', '')
        self.gene_label.setText(f"Protein Degradation: {gene_name} ({gene_id})")

        protein = d.get('protein_count', 0.0)
        prot_hl = d.get('protein_half_life', 90000.0)
        k_deg = math.log(2) / prot_hl if prot_hl > 0 else 1e-6
        current_decay = k_deg * protein

        # ── Update LIVE time-series plot ──
        self.plot.add_point(sim_time, protein)
        tl_rate = d.get('translation_rate', 0.1)
        mrna = d.get('mrna_count', 0.0)
        ss_prot = (tl_rate * mrna) / k_deg if k_deg > 0 else 0
        self.plot.set_steady_state(ss_prot)

        # Decay plot
        self.decay_plot.update_curve(protein, prot_hl)

        # Parameters
        hl_hr = prot_hl / 3600.0
        self.p_half_life.set_value(f"{hl_hr:.1f} hours ({prot_hl:.0f} s)")
        self.p_source.set_value("N-End Rule + Gene Features 🧬", "#CE93D8")

        # N-End Rule
        sequence = d.get('protein_sequence', 'M')
        n_term = sequence[0] if sequence else 'M'
        rule_class, life_cat, color = self.N_END_RULE.get(n_term, ('Unknown', 'Unknown', '#888'))
        self.p_n_terminal.set_value(f"{n_term} ({self._aa_name(n_term)})", color)
        self.p_stability.set_value(f"{rule_class.upper()} {'✅' if 'Stab' in rule_class else '⚠️'}", color)
        self.p_k_deg.set_value(f"{k_deg:.7f} /s  [ln(2) / {prot_hl:.0f}]")
        self.p_current_decay.set_value(f"{current_decay:.3f} mol/s", "#CE93D8")

        t50 = prot_hl / 3600.0
        t10 = prot_hl * 3.32 / 3600.0
        self.p_time_50.set_value(f"{t50:.1f} hours")
        self.p_time_10.set_value(f"{t10:.1f} hours ({t10/24:.1f} days)")

        # N-End detail
        n_seq = sequence[:20] if len(sequence) > 20 else sequence
        self.nend_seq.set_value(n_seq if n_seq else "Unknown")
        self.nend_class.set_value(f"{rule_class} ({life_cat} half-life)", color)
        hl_range = "20–100 hr" if "Stab" in rule_class else "2–20 hr" if "Second" in rule_class else "0.5–5 hr"
        self.nend_range.set_value(hl_range)

        # Stability summary - align with N-End Rule if data is ambiguous
        n_end_stab = "Stabilizing" if "Stab" in rule_class else "Destabilizing"
        
        if hl_hr > 10 or (hl_hr == 0 and n_end_stab == "Stabilizing"):
            self.rt_stability.set_value("VERY STABLE ✅", "#4CAF50")
            self.rt_turnover.set_value("SLOW (essential gene characteristic)", "#4CAF50")
        elif hl_hr > 2:
            self.rt_stability.set_value("MODERATELY STABLE ⚠️", "#FF9800")
            self.rt_turnover.set_value("MEDIUM (regulatory protein range)", "#FF9800")
        else:
            # If N-End rule says stabilizing but half-life is low (e.g. 0), 
            # we respect the rule for the warning level.
            if n_end_stab == "Stabilizing":
                self.rt_stability.set_value("STABLE (N-End Rule) ℹ️", "#81C784")
                self.rt_turnover.set_value("NORMAL (predicted from sequence)", "#81C784")
            else:
                self.rt_stability.set_value("UNSTABLE 🔴", "#F44336")
                self.rt_turnover.set_value("FAST (rapid response protein)", "#F44336")
        
        self.last_sync_label.setText(f"Last Sync: {sim_time:.1f}s")

    @staticmethod
    def _aa_name(aa: str) -> str:
        names = {
            'M': 'Methionine', 'A': 'Alanine', 'G': 'Glycine', 'V': 'Valine',
            'L': 'Leucine', 'I': 'Isoleucine', 'P': 'Proline', 'F': 'Phenylalanine',
            'W': 'Tryptophan', 'S': 'Serine', 'T': 'Threonine', 'C': 'Cysteine',
            'Y': 'Tyrosine', 'N': 'Asparagine', 'Q': 'Glutamine', 'D': 'Aspartate',
            'E': 'Glutamate', 'K': 'Lysine', 'R': 'Arginine', 'H': 'Histidine',
        }
        return names.get(aa, 'Unknown')


# ═══════════════════════════════════════════════════════════
#  BASE DIALOG MIXIN (shared gene selector logic)
# ═══════════════════════════════════════════════════════════
class _DynamicsDialogBase(QDialog):
    """Shared logic for all 4 standalone dynamics dialogs."""

    _TITLE = "Dynamics Monitor"
    _MIN_W, _MIN_H = 640, 920

    def __init__(self, inner_widget_cls, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self._TITLE)
        self.setMinimumSize(self._MIN_W, self._MIN_H)
        self.resize(self._MIN_W, self._MIN_H + 30)
        self.setStyleSheet("background: #0f0f23; color: #e0e0e0;")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)

        # Gene selector
        selector_layout = QHBoxLayout()
        selector_layout.addWidget(QLabel("Gene:"))
        self.gene_combo = QComboBox()
        self.gene_combo.setStyleSheet(
            "QComboBox { background: #16213e; color: white; border: 1px solid #444; padding: 4px; }"
        )
        self.gene_combo.currentIndexChanged.connect(self._on_gene_changed)
        selector_layout.addWidget(self.gene_combo, 1)
        layout.addLayout(selector_layout)

        # Inner widget wrapped in scroll area
        self._inner = inner_widget_cls()
        scroll = QScrollArea()
        scroll.setWidget(self._inner)
        scroll.setWidgetResizable(True)
        scroll.setStyleSheet("QScrollArea { border: none; background: #0f0f23; }")
        layout.addWidget(scroll)

        # State
        self._gene_data_list = []
        self._selected_index = 0
        self._frame_count = 0

    def populate_genes(self, gene_data_list: list):
        self._gene_data_list = gene_data_list
        self.gene_combo.blockSignals(True)
        self.gene_combo.clear()
        for i, d in enumerate(gene_data_list):
            gid = d.get('gene_id', f'Gene_{i}')
            gname = d.get('gene_name', '')
            label = f"{gid} ({gname})" if gname else gid
            self.gene_combo.addItem(label)
        self.gene_combo.blockSignals(False)

    def set_gene(self, index: int):
        if 0 <= index < self.gene_combo.count():
            self.gene_combo.setCurrentIndex(index)
            self._selected_index = index
            if not self.isVisible():
                self.show()
                self.raise_()

    def _on_gene_changed(self, index: int):
        self._selected_index = index
        # Clear plot histories if the inner widget has them
        if hasattr(self._inner, 'plot'):
            self._inner.plot.clear()
        if hasattr(self._inner, '_cumulative_deg'):
            self._inner._cumulative_deg = 0.0

    def _do_update(self, gene_data_list: list, sim_time: float):
        """Override-friendly update hook."""
        raise NotImplementedError

    def update_from_simulation(self, gene_data_list: list, sim_time: float):
        if not self.isVisible():
            return
        self._gene_data_list = gene_data_list
        self._frame_count += 1
        if self._frame_count % 3 != 0:
            return
        idx = self._selected_index
        if idx < 0 or idx >= len(gene_data_list):
            return
        d = gene_data_list[idx]
        self._inner.update_data(d, sim_time)
        self._do_update(gene_data_list, sim_time)


# ═══════════════════════════════════════════════════════════
#  4 STANDALONE DIALOGS
# ═══════════════════════════════════════════════════════════

class MRNATranscriptionDialog(_DynamicsDialogBase):
    _TITLE = "📝 mRNA Transcription — Live Monitor"

    def __init__(self, parent=None):
        super().__init__(MRNATranscriptionTab, parent)

    def _do_update(self, gene_data_list, sim_time):
        self._inner.plot.refresh()


class ProteinTranslationDialog(_DynamicsDialogBase):
    _TITLE = "🧬 Protein Translation — Live Monitor"

    def __init__(self, parent=None):
        super().__init__(ProteinTranslationTab, parent)

    def _do_update(self, gene_data_list, sim_time):
        self._inner.plot.refresh()


class MRNADegradationDialog(_DynamicsDialogBase):
    _TITLE = "💀 mRNA Degradation — Live Monitor"
    _MIN_W, _MIN_H = 640, 1050

    def __init__(self, parent=None):
        super().__init__(MRNADegradationTab, parent)

    def _do_update(self, gene_data_list, sim_time):
        self._inner.plot.refresh()  # refresh the live time-series plot


class ProteinDegradationDialog(_DynamicsDialogBase):
    _TITLE = "🛡️ Protein Stability — Live Monitor"
    _MIN_W, _MIN_H = 640, 1050

    def __init__(self, parent=None):
        super().__init__(ProteinDegradationTab, parent)

    def _do_update(self, gene_data_list, sim_time):
        self._inner.plot.refresh()  # refresh the live time-series plot


# ═══════════════════════════════════════════════════════════
#  BACKWARD-COMPATIBLE WRAPPER (kept for gui_main references)
# ═══════════════════════════════════════════════════════════
class GeneDynamicsPanel:
    """
    Facade that manages the 4 individual dialogs.
    gui_main.py still calls populate_genes / set_gene / update_from_simulation
    on this object — it fans out to each dialog.
    """

    def __init__(self, parent=None):
        self.mrna_tx = MRNATranscriptionDialog(parent)
        self.protein_tl = ProteinTranslationDialog(parent)
        self.mrna_deg = MRNADegradationDialog(parent)
        self.protein_deg = ProteinDegradationDialog(parent)
        self._dialogs = [self.mrna_tx, self.protein_tl, self.mrna_deg, self.protein_deg]

    def populate_genes(self, gene_data_list: list):
        for d in self._dialogs:
            d.populate_genes(gene_data_list)

    def set_gene(self, index: int):
        for d in self._dialogs:
            d.set_gene(index)

    def show(self):
        for d in self._dialogs:
            d.show()

    def raise_(self):
        for d in self._dialogs:
            d.raise_()

    def isVisible(self):
        return any(d.isVisible() for d in self._dialogs)

    def update_from_simulation(self, gene_data_list: list, sim_time: float):
        for d in self._dialogs:
            d.update_from_simulation(gene_data_list, sim_time)

