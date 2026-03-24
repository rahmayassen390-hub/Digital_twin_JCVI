"""
Knockout Effects Dashboard for JCVI Genome Analyzer
====================================================

Visual dashboard showing real-time cascade effects when genes are knocked out.

Features:
- Knocked out gene info panel
- Affected genes list with color-coded effect strength
- Cluster impact visualization
- Real-time effect propagation display
"""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QFrame,
    QGridLayout, QGroupBox, QScrollArea, QWidget, QProgressBar,
    QPushButton, QTableWidget, QTableWidgetItem, QHeaderView,
    QSplitter, QSizePolicy, QTabWidget
)
from PyQt5.QtCore import Qt, QTimer, pyqtSignal
from PyQt5.QtGui import QFont, QColor, QPalette, QBrush
from typing import Dict, List, Tuple, Optional


# Effect color gradient (matches gene_dependency_graph.py)
def get_effect_color(effect_strength: float) -> QColor:
    """Get QColor based on effect strength."""
    if effect_strength >= 0.7:
        return QColor(255, 50, 50)      # Red - high impact
    elif effect_strength >= 0.4:
        return QColor(255, 165, 0)      # Orange - medium impact
    elif effect_strength >= 0.1:
        return QColor(255, 255, 0)      # Yellow - low impact
    else:
        return QColor(255, 255, 180)    # Light yellow - minimal


def get_effect_icon(effect_strength: float) -> str:
    """Get emoji icon based on effect strength."""
    if effect_strength >= 0.7:
        return "🔴"
    elif effect_strength >= 0.4:
        return "🟠"
    elif effect_strength >= 0.1:
        return "🟡"
    else:
        return "⚪"


class EffectProgressBar(QProgressBar):
    """Custom progress bar showing effect strength."""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setTextVisible(True)
        self.setMinimum(0)
        self.setMaximum(100)
        self.setFixedHeight(20)
    
    def set_effect(self, effect_strength: float):
        """Set effect strength (0.0 to 1.0)."""
        self.setValue(int(effect_strength * 100))
        
        # Color based on effect
        color = get_effect_color(effect_strength)
        self.setStyleSheet(f"""
            QProgressBar {{
                border: 1px solid #444;
                border-radius: 3px;
                text-align: center;
                background-color: #2a2a3a;
                color: white;
            }}
            QProgressBar::chunk {{
                background-color: {color.name()};
                border-radius: 2px;
            }}
        """)


class ClusterImpactWidget(QFrame):
    """Widget showing impact on a single cluster."""
    
    def __init__(self, cluster_name: str, effect: float, parent=None):
        super().__init__(parent)
        self.setFrameStyle(QFrame.StyledPanel)
        self.setStyleSheet("background-color: #1a1a2e; border-radius: 5px; padding: 5px;")
        
        layout = QHBoxLayout(self)
        layout.setContentsMargins(8, 4, 8, 4)
        
        # Cluster name
        name_label = QLabel(f"{get_effect_icon(effect)} {cluster_name}")
        name_label.setStyleSheet("color: white; font-weight: bold;")
        name_label.setFixedWidth(200)
        layout.addWidget(name_label)
        
        # Effect bar
        self.progress = EffectProgressBar()
        self.progress.set_effect(effect)
        layout.addWidget(self.progress)
        
        # Percentage
        pct_label = QLabel(f"{int(effect * 100)}%")
        pct_label.setStyleSheet("color: white; font-weight: bold;")
        pct_label.setFixedWidth(50)
        pct_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        layout.addWidget(pct_label)


class KnockoutEffectsDashboard(QDialog):
    """
    Dashboard showing real-time cascade effects when genes are knocked out.
    
    Supports MULTIPLE knocked genes using tabs.
    
    Displays:
    - Knocked out gene info
    - List of affected genes with effect strength
    - Cluster impact summary
    - Real-time decay visualization
    """
    
    # Signal emitted when dashboard is closed
    dashboard_closed = pyqtSignal()
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("🛑 Knockout Effects Dashboard")
        self.setMinimumSize(600, 400)
        self.resize(800, 500)
        self.setStyleSheet("""
            QDialog {
                background-color: #0f0f23;
            }
            QLabel {
                color: #e0e0e0;
            }
            QGroupBox {
                color: #00d4ff;
                font-weight: bold;
                border: 1px solid #00d4ff;
                border-radius: 5px;
                margin-top: 10px;
                padding-top: 10px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px;
            }
            QTabWidget::pane {
                border: 1px solid #444;
                background-color: #1a1a2e;
            }
            QTabBar::tab {
                background-color: #2a2a4a;
                color: white;
                padding: 8px 15px;
                border-top-left-radius: 5px;
                border-top-right-radius: 5px;
            }
            QTabBar::tab:selected {
                background-color: #ff4444;
                color: white;
            }
        """)
        
        # Track multiple knocked genes: gene_id -> {cascade_effects, gene_info, tab_index}
        self.knocked_genes: Dict[str, Dict] = {}
        
        self.init_ui()
    
    def init_ui(self):
        """Initialize the dashboard UI with tab support for multiple knocked genes."""
        layout = QVBoxLayout(self)
        layout.setSpacing(10)
        
        # === HEADER ===
        header = QLabel("🛑 KNOCKOUT CASCADE EFFECTS")
        header.setStyleSheet("""
            font-size: 18px;
            font-weight: bold;
            color: #ff4444;
            padding: 10px;
            background-color: #1a1a2e;
            border-radius: 5px;
        """)
        header.setAlignment(Qt.AlignCenter)
        layout.addWidget(header)
        
        # === TAB WIDGET FOR MULTIPLE KNOCKED GENES ===
        self.tab_widget = QTabWidget()
        self.tab_widget.setTabsClosable(True)
        self.tab_widget.tabCloseRequested.connect(self._close_tab)
        layout.addWidget(self.tab_widget)
        
        # === BUTTONS ===
        btn_layout = QHBoxLayout()
        
        self.close_btn = QPushButton("✕ Close Dashboard")
        self.close_btn.setStyleSheet("""
            QPushButton {
                background-color: #444;
                color: white;
                padding: 10px 20px;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #666;
            }
        """)
        self.close_btn.clicked.connect(self.close)
        btn_layout.addStretch()
        btn_layout.addWidget(self.close_btn)
        
        layout.addLayout(btn_layout)
    
    def _close_tab(self, index: int):
        """Handle tab close request."""
        widget = self.tab_widget.widget(index)
        if widget:
            gene_id = widget.property("gene_id")
            if gene_id and gene_id in self.knocked_genes:
                del self.knocked_genes[gene_id]
        
        self.tab_widget.removeTab(index)
        
        # Close dashboard if no tabs left
        if self.tab_widget.count() == 0:
            self.close()
    
    def set_knockout(self, gene_id: str, cascade_effects: Dict, gene_info: Dict = None):
        """
        Display knockout effects for a gene (adds as a tab if multiple genes).
        
        Args:
            gene_id: The knocked out gene ID
            cascade_effects: Result from GeneDependencyGraph.calculate_cascade()
            gene_info: Optional additional gene info (name, product, etc.)
        """
        # If gene already has a tab, update it instead of creating new
        if gene_id in self.knocked_genes:
            tab_index = self.knocked_genes[gene_id].get('tab_index')
            if tab_index is not None and tab_index < self.tab_widget.count():
                self.tab_widget.setCurrentIndex(tab_index)
                self.show()
                return
        
        # Create new tab for this gene
        tab_widget = self._create_gene_tab(gene_id, cascade_effects, gene_info)
        tab_widget.setProperty("gene_id", gene_id)
        
        # Determine tab name
        gene_name = gene_info.get('gene_name', gene_info.get('name', '')) if gene_info else ''
        tab_name = f"{gene_id}" if not gene_name else f"{gene_name}"
        
        tab_index = self.tab_widget.addTab(tab_widget, f"🛑 {tab_name}")
        self.tab_widget.setCurrentIndex(tab_index)
        
        # Store reference
        self.knocked_genes[gene_id] = {
            'cascade_effects': cascade_effects,
            'gene_info': gene_info,
            'tab_index': tab_index
        }
        
        self.show()
    
    def remove_knockout(self, gene_id: str):
        """Remove a knocked gene tab."""
        if gene_id in self.knocked_genes:
            tab_info = self.knocked_genes[gene_id]
            tab_index = tab_info.get('tab_index')
            if tab_index is not None:
                # Find actual index (may have shifted)
                for i in range(self.tab_widget.count()):
                    widget = self.tab_widget.widget(i)
                    if widget.property("gene_id") == gene_id:
                        self._close_tab(i)
                        break
    
    def _create_gene_tab(self, gene_id: str, cascade_effects: Dict, gene_info: Dict = None) -> QWidget:
        """Create a tab widget for a single knocked gene."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        # === GENE INFO HEADER ===
        gene_info_group = QGroupBox("Knocked Out Gene")
        gene_info_layout = QGridLayout(gene_info_group)
        
        gene_name = gene_info.get('gene_name', gene_info.get('name', '')) if gene_info else ''
        display_name = f"{gene_id}" + (f" ({gene_name})" if gene_name else "")
        
        gene_id_label = QLabel(f"Gene ID: {display_name}")
        gene_id_label.setStyleSheet("font-size: 14px; font-weight: bold; color: #ff6666;")
        gene_info_layout.addWidget(gene_id_label, 0, 0)
        
        cluster_name = cascade_effects.get('cluster_name', 'Unknown')
        cluster_label = QLabel(f"Cluster: {cluster_name}")
        cluster_label.setStyleSheet("font-size: 12px;")
        gene_info_layout.addWidget(cluster_label, 0, 1)
        
        # Impact info
        is_high_impact = cascade_effects.get('is_high_impact', False)
        hub_info = cascade_effects.get('hub_info')
        uses_scientific_rules = cascade_effects.get('uses_scientific_rules', False)
        
        if is_high_impact and hub_info:
            hub_type, hub_desc = hub_info
            impact_text = f"🌐 REGULATORY HUB: {hub_type}"
        elif is_high_impact:
            impact_text = "⚠️ HIGH IMPACT"
        else:
            impact_text = "Normal Impact"
        
        if uses_scientific_rules:
            impact_text += " • 🔬 Scientific Rules"
        
        impact_label = QLabel(f"Impact: {impact_text}")
        impact_label.setStyleSheet("font-size: 12px;")
        gene_info_layout.addWidget(impact_label, 1, 0)
        
        total_affected = cascade_effects.get('total_affected', 0)
        total_label = QLabel(f"Affected Genes: {total_affected}")
        total_label.setStyleSheet("font-size: 12px;")
        gene_info_layout.addWidget(total_label, 1, 1)
        
        layout.addWidget(gene_info_group)
        
        # === SPLITTER FOR TABLE AND CLUSTER ===
        splitter = QSplitter(Qt.Horizontal)
        
        # === AFFECTED GENES TABLE ===
        affected_group = QGroupBox("Affected Genes")
        affected_layout = QVBoxLayout(affected_group)
        
        affected_table = QTableWidget()
        affected_table.setColumnCount(6)
        affected_table.setHorizontalHeaderLabels(["🎯", "Gene ID", "Effect", "Confidence", "Reason", "Cluster"])
        affected_table.horizontalHeader().setStretchLastSection(True)
        affected_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Fixed)
        affected_table.setColumnWidth(0, 30)
        affected_table.setColumnWidth(1, 120)
        affected_table.setColumnWidth(2, 60)
        affected_table.setColumnWidth(3, 80)
        affected_table.setColumnWidth(4, 220)
        affected_table.verticalHeader().setVisible(False)
        affected_table.setStyleSheet("""
            QTableWidget {
                background-color: #1a1a2e;
                color: white;
                gridline-color: #333;
            }
            QHeaderView::section {
                background-color: #2a2a4a;
                color: #00d4ff;
                padding: 5px;
                border: none;
            }
        """)
        
        # Populate table
        self._populate_table(affected_table, cascade_effects)
        
        affected_layout.addWidget(affected_table)
        splitter.addWidget(affected_group)
        
        # === CLUSTER IMPACT ===
        cluster_group = QGroupBox("Cluster Impact")
        cluster_layout = QVBoxLayout(cluster_group)
        
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setStyleSheet("border: none;")
        
        cluster_container = QWidget()
        cluster_inner_layout = QVBoxLayout(cluster_container)
        cluster_inner_layout.setSpacing(5)
        
        # Populate cluster impacts
        cluster_impact = cascade_effects.get('cluster_impact', {})
        for cluster_name, effect in sorted(cluster_impact.items(), key=lambda x: -x[1]):
            widget = ClusterImpactWidget(cluster_name, effect)
            cluster_inner_layout.addWidget(widget)
        
        cluster_inner_layout.addStretch()
        scroll.setWidget(cluster_container)
        cluster_layout.addWidget(scroll)
        splitter.addWidget(cluster_group)
        
        splitter.setSizes([500, 300])
        layout.addWidget(splitter)
        
        return tab
    
    def _populate_table(self, table: QTableWidget, cascade_effects: Dict):
        """Populate the affected genes table with scientific rule information."""
        table.setRowCount(0)
        
        # Combine same cluster and connected genes
        all_affected = []
        knocked_cluster_name = cascade_effects.get('cluster_name', '')
        
        # Handle new format with reason/confidence or legacy format
        for entry in cascade_effects.get('same_cluster_genes', []):
            if len(entry) >= 5:
                # New format: (gene_id, effect, reason, confidence, rule_type)
                gene_id, effect, reason, confidence, rule_type = entry[:5]
                all_affected.append((gene_id, effect, reason, confidence, rule_type, knocked_cluster_name))
            elif len(entry) == 2:
                # Legacy format: (gene_id, effect)
                gene_id, effect = entry
                all_affected.append((gene_id, effect, "Same cluster", "medium", "core_process", knocked_cluster_name))
        
        for entry in cascade_effects.get('connected_genes', []):
            if len(entry) >= 6:
                # New format: (gene_id, effect, reason, confidence, rule_type, cluster_name)
                gene_id, effect, reason, confidence, rule_type, cluster_name = entry[:6]
                all_affected.append((gene_id, effect, reason, confidence, rule_type, cluster_name))
            elif len(entry) == 4:
                # Legacy format: (gene_id, effect, distance, cluster_name)
                gene_id, effect, distance, cluster_name = entry
                all_affected.append((gene_id, effect, f"Connected cluster (distance {distance})", "low", "proximity", cluster_name))
        
        # Sort by effect strength (highest first)
        all_affected.sort(key=lambda x: -x[1])
        
        # Limit to top 100 for performance
        all_affected = all_affected[:100]
        
        table.setRowCount(len(all_affected))
        
        # Confidence icons
        confidence_icons = {"high": "🔴", "medium": "🟠", "low": "🟡"}
        rule_icons = {
            "operon": "📜",
            "core_process": "⚙️",
            "pathway": "🔬",
            "hub": "🌐",
            "proximity": "📍",
            "stress": "⚠️"
        }
        
        for row, (gene_id, effect, reason, confidence, rule_type, cluster_name) in enumerate(all_affected):
            # Rule type icon
            rule_icon = rule_icons.get(rule_type, "❓")
            icon_item = QTableWidgetItem(rule_icon)
            icon_item.setTextAlignment(Qt.AlignCenter)
            icon_item.setToolTip(f"Rule: {rule_type}")
            table.setItem(row, 0, icon_item)
            
            # Gene ID with color background
            gene_item = QTableWidgetItem(gene_id)
            gene_item.setBackground(QBrush(get_effect_color(effect)))
            gene_item.setForeground(QBrush(QColor("black")))
            table.setItem(row, 1, gene_item)
            
            # Effect percentage
            effect_item = QTableWidgetItem(f"{int(effect * 100)}%")
            effect_item.setTextAlignment(Qt.AlignCenter)
            table.setItem(row, 2, effect_item)
            
            # Confidence with icon
            conf_icon = confidence_icons.get(confidence, "⚪")
            conf_item = QTableWidgetItem(f"{conf_icon} {confidence.capitalize()}")
            table.setItem(row, 3, conf_item)
            
            # Reason
            reason_item = QTableWidgetItem(reason)
            reason_item.setToolTip(reason)  # Full text on hover
            table.setItem(row, 4, reason_item)
            
            # Cluster
            cluster_item = QTableWidgetItem(cluster_name)
            table.setItem(row, 5, cluster_item)
    
    def _populate_cluster_impact(self, cluster_impact: Dict):
        """Populate the cluster impact section."""
        # Clear existing
        while self.cluster_layout.count():
            child = self.cluster_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
        
        # Import cluster names
        try:
            from cluster_definitions import get_cluster_name
        except ImportError:
            def get_cluster_name(cid):
                return f"Cluster {cid}"
        
        # Sort clusters by impact (highest first)
        sorted_clusters = sorted(cluster_impact.items(), key=lambda x: -x[1])
        
        for cluster_id, effect in sorted_clusters:
            if effect > 0.01:  # Only show non-trivial effects
                cluster_name = get_cluster_name(cluster_id)
                widget = ClusterImpactWidget(cluster_name, effect)
                self.cluster_layout.addWidget(widget)
        
        self.cluster_layout.addStretch()
    
    def update_realtime(self, current_levels: Dict):
        """
        Update display with current transcription/translation levels.
        
        Args:
            current_levels: Dict of gene_id -> current_level (0-100)
        """
        # Update visual effects based on current levels if needed
        pass
    
    def clear_knockout(self):
        """Clear the dashboard when gene is resumed."""
        self.knocked_gene = None
        self.cascade_effects = None
        self.gene_id_label.setText("Gene ID: --")
        self.cluster_label.setText("Cluster: --")
        self.impact_label.setText("Impact: --")
        self.total_affected_label.setText("Affected Genes: --")
        self.affected_table.setRowCount(0)
        
        # Clear cluster widgets
        while self.cluster_layout.count():
            child = self.cluster_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
    
    def closeEvent(self, event):
        """Handle close event."""
        self.dashboard_closed.emit()
        super().closeEvent(event)


# =============================================================================
# STANDALONE TESTING
# =============================================================================

if __name__ == "__main__":
    import sys
    from PyQt5.QtWidgets import QApplication
    
    app = QApplication(sys.argv)
    
    # Create test dashboard
    dashboard = KnockoutEffectsDashboard()
    
    # Test data
    test_cascade = {
        'knocked_gene': 'JCVISYN3_0001',
        'cluster_name': 'DNA Replication',
        'cluster_id': 3,
        'is_high_impact': True,
        'total_affected': 25,
        'same_cluster_genes': [
            ('JCVISYN3_0002', 0.85),
            ('JCVISYN3_0003', 0.85),
            ('JCVISYN3_0004', 0.85),
        ],
        'connected_genes': [
            ('JCVISYN3_0010', 0.80, 1, 'Transcription'),
            ('JCVISYN3_0011', 0.80, 1, 'Transcription'),
            ('JCVISYN3_0020', 0.50, 2, 'Cell Division'),
            ('JCVISYN3_0021', 0.50, 2, 'Cell Division'),
            ('JCVISYN3_0030', 0.30, 3, 'Protein Folding'),
        ],
        'cluster_impact': {
            3: 1.0,   # DNA Replication
            4: 0.80,  # Transcription
            7: 0.50,  # Cell Division
            11: 0.30, # Protein Folding
        }
    }
    
    dashboard.set_knockout('JCVISYN3_0001', test_cascade, {'gene_name': 'dnaA'})
    
    sys.exit(app.exec_())
