"""
Metabolic Dashboard — Real-time visualization of metabolic solver state.
Features:
  1. Flux Panel: Top reactions by magnitude with color-coded bars
  2. Metabolite Plots: Time-series of key metabolite proxies
  3. FBA Status Indicator: Solver state, call count, doubling time
"""

import math
from collections import deque

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QWidget,
    QTableWidget, QTableWidgetItem, QHeaderView,
    QLabel, QGroupBox, QSplitter, QFrame,
    QProgressBar, QPushButton, QSizePolicy
)
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QFont, QColor, QPainter, QLinearGradient

# Matplotlib embedding
try:
    import matplotlib
    matplotlib.use('Qt5Agg')
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


# ═══════════════════════════════════════════════════════════════════════
# Styled ProgressBar for flux magnitude
# ═══════════════════════════════════════════════════════════════════════
class FluxBar(QProgressBar):
    """A styled progress bar representing flux magnitude."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setTextVisible(False)
        self.setMinimum(0)
        self.setMaximum(1000)
        self.setFixedHeight(16)
        self._color = "#4CAF50"
    
    def set_flux(self, value, max_val):
        """Set bar value scaled to max_val."""
        if max_val <= 0:
            self.setValue(0)
            return
        pct = min(1.0, abs(value) / max_val)
        self.setValue(int(pct * 1000))
        
        # Color by magnitude
        if value < 0:
            self._color = "#E53935"  # Red for reverse
        elif abs(value) > 10:
            self._color = "#43A047"  # Green for high
        elif abs(value) > 1:
            self._color = "#FFB300"  # Amber for medium
        else:
            self._color = "#78909C"  # Blue-gray for low
        
        self.setStyleSheet(f"""
            QProgressBar {{
                border: 1px solid #37474F;
                border-radius: 4px;
                background: #263238;
            }}
            QProgressBar::chunk {{
                background: {self._color};
                border-radius: 3px;
            }}
        """)


# ═══════════════════════════════════════════════════════════════════════
# Feature 1: Flux Table Widget
# ═══════════════════════════════════════════════════════════════════════
class FluxTableWidget(QWidget):
    """Displays top N metabolic reactions sorted by flux magnitude."""
    
    # Common reaction name mapping for readability
    RXN_NAMES = {
        'EX_biomass_c': 'Biomass Export',
        'EX_glc__D_e': 'Glucose Uptake',
        'EX_o2_e': 'O₂ Uptake',
        'EX_co2_e': 'CO₂ Export',
        'EX_h2o_e': 'H₂O Exchange',
        'EX_pyr_e': 'Pyruvate Uptake',
        'EX_ac_e': 'Acetate Export',
        'EX_lac__L_e': 'Lactate Export',
        'EX_h_e': 'H⁺ Export',
        'PFK': 'Phosphofructokinase',
        'PGK': 'Phosphoglycerate Kinase',
        'PYK': 'Pyruvate Kinase',
        'GAPD': 'GAPDH',
        'TPI': 'Triose-P Isomerase',
        'ENO': 'Enolase',
        'PGI': 'Glucose-6P Isomerase',
        'LDH_L': 'Lactate Dehydrogenase',
        'ATPS4r': 'ATP Synthase',
        'GLCpts': 'Glucose PTS Transport',
    }
    
    def __init__(self, parent=None, n_rows=10):
        super().__init__(parent)
        self.n_rows = n_rows
        self.flux_bars = []
        self._init_ui()
    
    def _init_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Title
        title = QLabel("⚡ Active Metabolic Reactions (Top 10)")
        title.setFont(QFont("Segoe UI", 11, QFont.Bold))
        title.setStyleSheet("color: #E0E0E0; padding: 4px;")
        layout.addWidget(title)
        
        # Table
        self.table = QTableWidget(self.n_rows, 4)
        self.table.setHorizontalHeaderLabels(["Reaction", "Flux (mmol/gDW/hr)", "Direction", "Activity"])
        self.table.horizontalHeader().setStyleSheet(
            "QHeaderView::section { background: #1A237E; color: white; font-weight: bold; padding: 4px; }"
        )
        self.table.setColumnWidth(0, 200)
        self.table.setColumnWidth(1, 140)
        self.table.setColumnWidth(2, 70)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setAlternatingRowColors(True)
        self.table.setStyleSheet("""
            QTableWidget {
                background-color: #1B2631;
                alternate-background-color: #212F3D;
                color: #E0E0E0;
                gridline-color: #37474F;
                border: 1px solid #37474F;
                border-radius: 6px;
                font-size: 11px;
            }
            QTableWidget::item:selected {
                background-color: #1565C0;
            }
        """)
        
        # Pre-create flux bars
        for i in range(self.n_rows):
            bar = FluxBar()
            self.flux_bars.append(bar)
            self.table.setCellWidget(i, 3, bar)
            self.table.setRowHeight(i, 24)
        
        layout.addWidget(self.table)
    
    def update_fluxes(self, fluxes: dict):
        """Update the table with new flux data."""
        if not fluxes:
            return
        
        # Sort by absolute magnitude
        sorted_fluxes = sorted(fluxes.items(), key=lambda x: abs(x[1]), reverse=True)
        top = sorted_fluxes[:self.n_rows]
        
        max_flux = abs(top[0][1]) if top else 1.0
        
        for i in range(self.n_rows):
            if i < len(top):
                rxn_id, flux_val = top[i]
                display_name = self.RXN_NAMES.get(rxn_id, rxn_id)
                direction = "→" if flux_val >= 0 else "←"
                
                name_item = QTableWidgetItem(display_name)
                name_item.setToolTip(f"ID: {rxn_id}")
                self.table.setItem(i, 0, name_item)
                
                flux_item = QTableWidgetItem(f"{flux_val:+.4f}")
                flux_item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                if flux_val < 0:
                    flux_item.setForeground(QColor("#EF5350"))
                elif flux_val > 10:
                    flux_item.setForeground(QColor("#66BB6A"))
                else:
                    flux_item.setForeground(QColor("#FFB74D"))
                self.table.setItem(i, 1, flux_item)
                
                dir_item = QTableWidgetItem(direction)
                dir_item.setTextAlignment(Qt.AlignCenter)
                self.table.setItem(i, 2, dir_item)
                
                self.flux_bars[i].set_flux(flux_val, max_flux)
            else:
                for col in range(3):
                    self.table.setItem(i, col, QTableWidgetItem(""))
                self.flux_bars[i].setValue(0)


# ═══════════════════════════════════════════════════════════════════════
# Feature 2: Metabolite Time-Series Plots
# ═══════════════════════════════════════════════════════════════════════
class MetabolitePlotWidget(QWidget):
    """Matplotlib-embedded time-series plots for key metabolic indicators."""
    
    MAX_HISTORY = 200  # Max data points to display
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.time_data = deque(maxlen=self.MAX_HISTORY)
        self.biomass_data = deque(maxlen=self.MAX_HISTORY)
        self.carbon_data = deque(maxlen=self.MAX_HISTORY)
        self.energy_data = deque(maxlen=self.MAX_HISTORY)
        self.redox_data = deque(maxlen=self.MAX_HISTORY)
        self._init_ui()
    
    def _init_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        title = QLabel("📊 Metabolic Indicators Over Time")
        title.setFont(QFont("Segoe UI", 11, QFont.Bold))
        title.setStyleSheet("color: #E0E0E0; padding: 4px;")
        layout.addWidget(title)
        
        if not MATPLOTLIB_AVAILABLE:
            fallback = QLabel("⚠️ matplotlib not installed. Install with: pip install matplotlib")
            fallback.setStyleSheet("color: #FF9800; padding: 20px;")
            layout.addWidget(fallback)
            return
        
        # Dark-themed figure
        self.fig = Figure(figsize=(8, 5), dpi=90, facecolor='#1B2631')
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        # Create 4 subplots
        gs = self.fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)
        
        self.ax_biomass = self.fig.add_subplot(gs[0, 0])
        self.ax_carbon = self.fig.add_subplot(gs[0, 1])
        self.ax_energy = self.fig.add_subplot(gs[1, 0])
        self.ax_redox = self.fig.add_subplot(gs[1, 1])
        
        for ax, title_text, color in [
            (self.ax_biomass, "Biomass Growth Rate", "#66BB6A"),
            (self.ax_carbon, "Carbon Flux (PFK+PYK)", "#42A5F5"),
            (self.ax_energy, "Energy (PGK Flux)", "#FFA726"),
            (self.ax_redox, "Total Protein", "#AB47BC")
        ]:
            ax.set_facecolor('#263238')
            ax.set_title(title_text, fontsize=9, color=color, fontweight='bold')
            ax.tick_params(colors='#90A4AE', labelsize=7)
            ax.spines['bottom'].set_color('#455A64')
            ax.spines['left'].set_color('#455A64')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.grid(True, alpha=0.15, color='#546E7A')
            ax.set_xlabel("Time (s)", fontsize=7, color='#90A4AE')
        
        self.fig.tight_layout(pad=1.5)
        layout.addWidget(self.canvas)
    
    def append_data(self, sim_time, fluxes, total_protein=0):
        """Append a new data point from FBA results."""
        self.time_data.append(sim_time)
        self.biomass_data.append(fluxes.get('EX_biomass_c', 0.0))
        
        # Carbon flux = sum of glycolytic fluxes
        carbon = abs(fluxes.get('PFK', 0.0)) + abs(fluxes.get('PYK', 0.0))
        self.carbon_data.append(carbon)
        
        # Energy proxy = PGK flux (ATP-generating step)
        self.energy_data.append(abs(fluxes.get('PGK', 0.0)))
        
        # Total protein
        self.redox_data.append(total_protein)
    
    def refresh_plots(self):
        """Redraw all 4 plots with current data."""
        if not MATPLOTLIB_AVAILABLE or len(self.time_data) < 1:
            return
        
        t = list(self.time_data)
        
        for ax, data, color in [
            (self.ax_biomass, list(self.biomass_data), '#66BB6A'),
            (self.ax_carbon, list(self.carbon_data), '#42A5F5'),
            (self.ax_energy, list(self.energy_data), '#FFA726'),
            (self.ax_redox, list(self.redox_data), '#AB47BC'),
        ]:
            ax.clear()
            ax.set_facecolor('#263238')
            ax.plot(t, data, color=color, linewidth=1.5, alpha=0.9)
            ax.fill_between(t, data, alpha=0.15, color=color)
            ax.tick_params(colors='#90A4AE', labelsize=7)
            ax.spines['bottom'].set_color('#455A64')
            ax.spines['left'].set_color('#455A64')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.grid(True, alpha=0.15, color='#546E7A')
            ax.set_xlabel("Time (s)", fontsize=7, color='#90A4AE')
        
        # Restore titles
        self.ax_biomass.set_title("Biomass Growth (hr⁻¹)", fontsize=9, color='#66BB6A', fontweight='bold')
        self.ax_carbon.set_title("Carbon Flux (PFK+PYK)", fontsize=9, color='#42A5F5', fontweight='bold')
        self.ax_energy.set_title("Energy Flux (PGK)", fontsize=9, color='#FFA726', fontweight='bold')
        self.ax_redox.set_title("Total Protein Count", fontsize=9, color='#AB47BC', fontweight='bold')
        
        self.fig.tight_layout(pad=1.5)
        self.canvas.draw_idle()
    
    def clear(self):
        """Reset all data."""
        self.time_data.clear()
        self.biomass_data.clear()
        self.carbon_data.clear()
        self.energy_data.clear()
        self.redox_data.clear()


# ═══════════════════════════════════════════════════════════════════════
# Feature 3: FBA Status Indicator
# ═══════════════════════════════════════════════════════════════════════
class FBAStatusWidget(QWidget):
    """Shows current solver state at a glance."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self._init_ui()
        self.fba_call_count = 0
    
    def _init_ui(self):
        layout = QHBoxLayout(self)
        layout.setContentsMargins(8, 4, 8, 4)
        
        self.setStyleSheet("""
            QWidget {
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                    stop:0 #1A237E, stop:1 #0D47A1);
                border-radius: 8px;
                border: 1px solid #1565C0;
            }
        """)
        self.setFixedHeight(50)
        
        # Status dot + label
        self.status_dot = QLabel("●")
        self.status_dot.setFont(QFont("Arial", 16))
        self.status_dot.setStyleSheet("color: #78909C; border: none; background: transparent;")
        layout.addWidget(self.status_dot)
        
        self.status_label = QLabel("Metabolic Solver: INACTIVE")
        self.status_label.setFont(QFont("Segoe UI", 10, QFont.Bold))
        self.status_label.setStyleSheet("color: #E0E0E0; border: none; background: transparent;")
        layout.addWidget(self.status_label)
        
        layout.addSpacing(20)
        
        # FBA info
        info_frame = QWidget()
        info_frame.setStyleSheet("background: transparent; border: none;")
        info_layout = QVBoxLayout(info_frame)
        info_layout.setContentsMargins(0, 0, 0, 0)
        info_layout.setSpacing(0)
        
        self.fba_time_label = QLabel("Last FBA: —")
        self.fba_time_label.setFont(QFont("Segoe UI", 8))
        self.fba_time_label.setStyleSheet("color: #90A4AE; border: none; background: transparent;")
        info_layout.addWidget(self.fba_time_label)
        
        self.fba_count_label = QLabel("FBA Calls: 0")
        self.fba_count_label.setFont(QFont("Segoe UI", 8))
        self.fba_count_label.setStyleSheet("color: #90A4AE; border: none; background: transparent;")
        info_layout.addWidget(self.fba_count_label)
        
        layout.addWidget(info_frame)
        layout.addSpacing(20)
        
        # Growth metrics
        metrics_frame = QWidget()
        metrics_frame.setStyleSheet("background: transparent; border: none;")
        metrics_layout = QVBoxLayout(metrics_frame)
        metrics_layout.setContentsMargins(0, 0, 0, 0)
        metrics_layout.setSpacing(0)
        
        self.biomass_label = QLabel("Biomass Flux: —")
        self.biomass_label.setFont(QFont("Segoe UI", 9, QFont.Bold))
        self.biomass_label.setStyleSheet("color: #66BB6A; border: none; background: transparent;")
        metrics_layout.addWidget(self.biomass_label)
        
        self.doubling_label = QLabel("Doubling Time: —")
        self.doubling_label.setFont(QFont("Segoe UI", 9, QFont.Bold))
        self.doubling_label.setStyleSheet("color: #42A5F5; border: none; background: transparent;")
        metrics_layout.addWidget(self.doubling_label)
        
        layout.addWidget(metrics_frame)
        layout.addStretch()
    
    def update_status(self, active, sim_time=0, biomass_flux=0, fba_calls=0):
        """Update all status indicators."""
        self.fba_call_count = fba_calls
        
        if active:
            self.status_dot.setStyleSheet("color: #66BB6A; border: none; background: transparent;")
            self.status_label.setText("Metabolic Solver: ACTIVE")
            self.status_label.setStyleSheet("color: #66BB6A; border: none; background: transparent;")
        else:
            self.status_dot.setStyleSheet("color: #78909C; border: none; background: transparent;")
            self.status_label.setText("Metabolic Solver: INACTIVE")
            self.status_label.setStyleSheet("color: #78909C; border: none; background: transparent;")
        
        self.fba_time_label.setText(f"Last FBA: t={sim_time:.0f}s")
        self.fba_count_label.setText(f"FBA Calls: {fba_calls}")
        
        if biomass_flux > 0:
            dt_min = (0.693 / biomass_flux) * 60.0
            self.biomass_label.setText(f"Biomass Flux: {biomass_flux:.4f} hr⁻¹")
            self.doubling_label.setText(f"Doubling Time: {dt_min:.1f} min")
        else:
            self.biomass_label.setText("Biomass Flux: 0.0")
            self.doubling_label.setText("Doubling Time: ∞")


# ═══════════════════════════════════════════════════════════════════════
# Main Dashboard Window
# ═══════════════════════════════════════════════════════════════════════
class MetabolicDashboard(QDialog):
    """
    Unified Metabolic Dashboard window combining all 3 features.
    Opened from AI Results tab via button click.
    """
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("🧬 Metabolic Dashboard — Real-Time FBA Visualization")
        self.setMinimumSize(900, 700)
        self.resize(1000, 750)
        self._init_ui()
        self._fba_call_count = 0
    
    def _init_ui(self):
        self.setStyleSheet("""
            QDialog {
                background-color: #0D1117;
            }
            QLabel {
                color: #E0E0E0;
            }
            QGroupBox {
                color: #E0E0E0;
                border: 1px solid #30363D;
                border-radius: 6px;
                margin-top: 8px;
                font-weight: bold;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 12px;
                padding: 0 4px;
            }
        """)
        
        main_layout = QVBoxLayout(self)
        main_layout.setSpacing(8)
        
        # ── Feature 3: FBA Status Bar (top) ──
        self.fba_status = FBAStatusWidget()
        main_layout.addWidget(self.fba_status)
        
        # ── Splitter: Flux Table (left) + Plots (right) ──
        splitter = QSplitter(Qt.Horizontal)
        splitter.setStyleSheet("QSplitter::handle { background: #30363D; width: 3px; }")
        
        # Left: Flux Table (Feature 1)
        self.flux_table = FluxTableWidget(n_rows=10)
        splitter.addWidget(self.flux_table)
        
        # Right: Metabolite Plots (Feature 2)
        self.metabolite_plots = MetabolitePlotWidget()
        splitter.addWidget(self.metabolite_plots)
        
        splitter.setSizes([400, 600])
        main_layout.addWidget(splitter)
        
        # ── Bottom bar with exchange reactions ──
        exchange_group = QGroupBox("🔄 Key Exchange Fluxes")
        exchange_layout = QHBoxLayout()
        
        self.exchange_labels = {}
        for name, default in [
            ("Glucose", "—"), ("O₂", "—"), ("CO₂", "—"), 
            ("Pyruvate", "—"), ("Biomass", "—")
        ]:
            lbl = QLabel(f"{name}: {default}")
            lbl.setFont(QFont("Segoe UI", 9))
            lbl.setStyleSheet("color: #B0BEC5; padding: 2px 8px;")
            exchange_layout.addWidget(lbl)
            self.exchange_labels[name] = lbl
        
        exchange_group.setLayout(exchange_layout)
        main_layout.addWidget(exchange_group)
    
    def update_from_fba(self, sim_time, fluxes, total_protein=0, fba_calls=0):
        """
        Called after each FBA solve to refresh the entire dashboard.
        
        Args:
            sim_time: Current simulation time in seconds
            fluxes: Dict of reaction_id -> flux_value
            total_protein: Total protein count from ODE engine
            fba_calls: Total number of FBA calls made
        """
        if not fluxes:
            return
        
        self._fba_call_count = fba_calls
        biomass_flux = fluxes.get('EX_biomass_c', 0.0)
        
        # Feature 1: Update flux table
        self.flux_table.update_fluxes(fluxes)
        
        # Feature 2: Append data and redraw plots
        self.metabolite_plots.append_data(sim_time, fluxes, total_protein)
        self.metabolite_plots.refresh_plots()
        
        # Feature 3: Update status bar
        self.fba_status.update_status(
            active=True,
            sim_time=sim_time,
            biomass_flux=biomass_flux,
            fba_calls=fba_calls
        )
        
        # Update exchange labels
        exchange_map = {
            "Glucose": ("EX_glc__D_e", "mmol/gDW/hr"),
            "O₂": ("EX_o2_e", "mmol/gDW/hr"),
            "CO₂": ("EX_co2_e", "mmol/gDW/hr"),
            "Pyruvate": ("EX_pyr_e", "mmol/gDW/hr"),
            "Biomass": ("EX_biomass_c", "hr⁻¹"),
        }
        
        for name, (rxn_id, unit) in exchange_map.items():
            val = fluxes.get(rxn_id, 0.0)
            direction = "↓" if val < 0 else "↑" if val > 0 else "—"
            color = "#EF5350" if val < -0.01 else "#66BB6A" if val > 0.01 else "#78909C"
            self.exchange_labels[name].setText(f"{name}: {val:+.2f} {unit} {direction}")
            self.exchange_labels[name].setStyleSheet(f"color: {color}; padding: 2px 8px;")
    
    def reset(self):
        """Reset all dashboard state."""
        self._fba_call_count = 0
        self.metabolite_plots.clear()
        self.fba_status.update_status(active=False)
