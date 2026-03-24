"""
Enhanced Knockout Dialog — Before/After Cascade Analysis
=========================================================

Shows detailed 6-rule cascade analysis with:
  - Primary effect (mRNA decay countdown)
  - Each biological rule result with affected genes
  - Predicted outcome (growth rate, viability, lethality)
  - Essentiality warning banner
  - Simulation preview
"""

import math
from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QGroupBox,
    QFrame, QPushButton, QScrollArea, QWidget, QSizePolicy
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont, QColor


class RuleBadge(QFrame):
    """Compact display for one knockout rule result."""

    def __init__(self, rule_num: int, rule_name: str, parent=None):
        super().__init__(parent)
        self.setStyleSheet(
            "background: #16213e; border: 1px solid #333; border-radius: 4px; padding: 6px;"
        )
        self._layout = QVBoxLayout(self)
        self._layout.setContentsMargins(6, 4, 6, 4)
        self._layout.setSpacing(2)

        self.header = QLabel(f"{rule_num}️⃣ {rule_name}")
        self.header.setFont(QFont("Arial", 9, QFont.Bold))
        self.header.setStyleSheet("color: #B0BEC5;")
        self._layout.addWidget(self.header)

        self.result_label = QLabel("Evaluating...")
        self.result_label.setFont(QFont("Consolas", 8))
        self.result_label.setStyleSheet("color: #888;")
        self.result_label.setWordWrap(True)
        self._layout.addWidget(self.result_label)

        self.detail_label = QLabel("")
        self.detail_label.setFont(QFont("Consolas", 8))
        self.detail_label.setStyleSheet("color: #666;")
        self.detail_label.setWordWrap(True)
        self._layout.addWidget(self.detail_label)
        self.detail_label.hide()

    def set_result(self, icon: str, text: str, color: str = "#B0BEC5",
                   detail: str = ""):
        self.result_label.setText(f"  {icon} {text}")
        self.result_label.setStyleSheet(f"color: {color};")
        if detail:
            self.detail_label.setText(f"    {detail}")
            self.detail_label.show()
        else:
            self.detail_label.hide()


class EnhancedKnockoutDialog(QDialog):
    """
    Pre-pause analysis dialog showing predicted cascade effects.

    Opened BEFORE actually pausing a gene. User can review analysis
    and choose to proceed or cancel.

    Returns:
        QDialog.Accepted if user clicks "Pause Anyway"
        QDialog.Rejected if user clicks "Cancel"
    """

    def __init__(self, gene_data: dict, cascade_effects: list = None,
                 total_genes: int = 498, parent=None):
        super().__init__(parent)
        self.setWindowTitle("🧬 Knockout Cascade Analysis")
        self.setMinimumSize(550, 600)
        self.resize(580, 700)
        self.setStyleSheet("background: #0f0f23; color: #e0e0e0;")

        self._gene_data = gene_data
        self._cascade = cascade_effects or []
        self._total_genes = total_genes
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setSpacing(4)

        gid = self._gene_data.get('gene_id', '?')
        gname = self._gene_data.get('gene_name', '')

        # Title
        title = QLabel(f"Knockout Cascade Analysis: {gname or gid}")
        title.setFont(QFont("Arial", 12, QFont.Bold))
        title.setStyleSheet("color: #FF7043;")
        layout.addWidget(title)

        # Scroll area for content
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setStyleSheet("QScrollArea { border: none; }")
        content = QWidget()
        content_layout = QVBoxLayout(content)
        content_layout.setSpacing(6)

        # ── PRIMARY EFFECT ──
        primary = QGroupBox("PRIMARY EFFECT")
        primary.setStyleSheet(
            "QGroupBox { color: #FF7043; border: 1px solid #444; border-radius: 4px; margin-top: 8px; padding-top: 16px; }"
            "QGroupBox::title { subcontrol-position: top left; padding: 0 8px; }"
        )
        p_layout = QVBoxLayout(primary)
        p_layout.setSpacing(2)

        tx_rate = self._gene_data.get('transcription_rate', 0.01)
        mrna = self._gene_data.get('mrna_count', 0.0)
        mrna_hl = self._gene_data.get('mrna_half_life', 300.0)
        hl_min = mrna_hl / 60.0

        p1 = QLabel(f"├─ {gname or gid} → Transcription STOPPED")
        p1.setFont(QFont("Consolas", 9))
        p1.setStyleSheet("color: #F44336;")
        p_layout.addWidget(p1)

        p2 = QLabel(f"└─ mRNA decay: {hl_min:.0f} seconds to 50%   (started immediately)")
        p2.setFont(QFont("Consolas", 9))
        p2.setStyleSheet("color: #FF9800;")
        p_layout.addWidget(p2)

        content_layout.addWidget(primary)

        # ── CASCADE RULES ──
        rules_group = QGroupBox("CASCADE EFFECTS (6 rules applied)")
        rules_group.setStyleSheet(
            "QGroupBox { color: #FFD54F; border: 1px solid #444; border-radius: 4px; margin-top: 8px; padding-top: 16px; }"
            "QGroupBox::title { subcontrol-position: top left; padding: 0 8px; }"
        )
        rules_layout = QVBoxLayout(rules_group)
        rules_layout.setSpacing(4)

        # Analyze cascade effects by type
        operon_effects = [e for e in self._cascade if 'operon' in e.get('rule_type', '').lower()
                          or 'polar' in e.get('reason', '').lower()]
        pathway_effects = [e for e in self._cascade if 'pathway' in e.get('rule_type', '').lower()
                           or 'pathway' in e.get('reason', '').lower()]
        hub_effects = [e for e in self._cascade if 'hub' in e.get('rule_type', '').lower()
                       or 'regulatory' in e.get('reason', '').lower()]
        cluster_effects = [e for e in self._cascade if 'cluster' in e.get('rule_type', '').lower()
                           or 'entangle' in e.get('reason', '').lower()]
        resource_effects = [e for e in self._cascade if 'resource' in e.get('rule_type', '').lower()
                            or 'ribosome' in e.get('reason', '').lower()]

        # Rule 1: Operon Polarity
        r1 = RuleBadge(1, "Operon Polarity")
        if operon_effects:
            genes = ", ".join([e.get('gene_id', '?')[:12] for e in operon_effects[:3]])
            r1.set_result("⚠️", f"{len(operon_effects)} downstream genes affected", "#FF9800",
                          f"Polar effect on: {genes}")
        else:
            r1.set_result("❌", "None (no downstream operon members)", "#4CAF50")
        rules_layout.addWidget(r1)

        # Rule 2: Pathway Bottleneck
        r2 = RuleBadge(2, "Pathway Bottleneck")
        if pathway_effects:
            genes = ", ".join([e.get('gene_id', '?')[:12] for e in pathway_effects[:3]])
            r2.set_result("⚠️", f"{len(pathway_effects)} pathway genes affected", "#FF9800",
                          f"Disrupted: {genes}")
        else:
            r2.set_result("❌", "None (no shared pathway)", "#4CAF50")
        rules_layout.addWidget(r2)

        # Rule 3: Regulatory Hub
        r3 = RuleBadge(3, "Regulatory Hub")
        if hub_effects:
            r3.set_result("🔴", f"{len(hub_effects)} genes affected (hub gene!)", "#F44336",
                          "RNA Polymerase / Sigma factor — affects ALL transcription")
        else:
            r3.set_result("❌", "None (not a transcription factor)", "#4CAF50")
        rules_layout.addWidget(r3)

        # Rule 4: Cluster Entanglement
        r4 = RuleBadge(4, "Cluster Entanglement")
        if cluster_effects:
            genes = ", ".join([e.get('gene_id', '?')[:12] for e in cluster_effects[:3]])
            r4.set_result("⚠️", f"{len(cluster_effects)} quantum-linked genes", "#FF9800",
                          f"Entangled with: {genes}")
        else:
            r4.set_result("❌", "None (no cluster entanglement)", "#4CAF50")
        rules_layout.addWidget(r4)

        # Rule 5: Ribosome/Polymerase Release
        r5 = RuleBadge(5, "Resource Release")
        freed_ribo = max(1, int(mrna * 0.02))
        r5.set_result("✅", f"Minor benefit — freed {freed_ribo} ribosomes", "#81C784",
                       f"Global translation speedup: +{freed_ribo * 0.05:.1f}%")
        rules_layout.addWidget(r5)

        # Rule 6: Essentiality
        r6 = RuleBadge(6, "Essentiality Check")
        category = self._gene_data.get('category', '').lower()
        func = self._gene_data.get('function', '').lower()
        is_essential = ('essential' in category or 'dna' in func or
                        'ribosom' in func or 'rrna' in func or 'rpo' in func or
                        'dnaa' in func.lower() or 'lipid' in func)
        if is_essential:
            r6.set_result("🔴", "ESSENTIAL GENE — LETHAL KNOCKOUT", "#F44336",
                          "Predicted: Cell death within 2–3 divisions")
        else:
            r6.set_result("✅", "Non-essential — cell should survive", "#4CAF50")
        rules_layout.addWidget(r6)

        content_layout.addWidget(rules_group)

        # ── PREDICTED OUTCOME ──
        outcome = QGroupBox("PREDICTED OUTCOME")
        outcome.setStyleSheet(
            "QGroupBox { color: #CE93D8; border: 1px solid #444; border-radius: 4px; margin-top: 8px; padding-top: 16px; }"
            "QGroupBox::title { subcontrol-position: top left; padding: 0 8px; }"
        )
        o_layout = QVBoxLayout(outcome)
        o_layout.setSpacing(2)

        total_affected = len(self._cascade)
        pct_affected = total_affected / self._total_genes * 100 if self._total_genes > 0 else 0

        growth_before = "105 min"
        viability_before = "100%"

        if is_essential:
            growth_after = "∞ (no division)"
            viability_after = "0% (LETHAL)"
            viab_color = "#F44336"
            time_to_death = "~210 minutes (2 cell cycles)"
        elif total_affected > 20:
            growth_after = "180–300 min (severely impaired)"
            viability_after = "30–60% (impaired)"
            viab_color = "#FF9800"
            time_to_death = "Cell may survive with reduced fitness"
        else:
            growth_after = "110–130 min (mildly impaired)"
            viability_after = "80–100% (viable)"
            viab_color = "#4CAF50"
            time_to_death = "Cell should survive normally"

        for text, color in [
            (f"├─ Growth rate: {growth_before} → {growth_after}", "#B0BEC5"),
            (f"├─ Viability: {viability_before} → {viability_after}", viab_color),
            (f"├─ Prediction: {time_to_death}", viab_color),
            (f"└─ Affected genes: {total_affected} total ({pct_affected:.1f}% of genome)", "#B0BEC5"),
        ]:
            lbl = QLabel(text)
            lbl.setFont(QFont("Consolas", 9))
            lbl.setStyleSheet(f"color: {color};")
            o_layout.addWidget(lbl)

        content_layout.addWidget(outcome)

        # ── WARNING BANNER ──
        if is_essential:
            warning = QLabel("  ⚠️ WARNING: Pausing this gene will likely KILL the cell!  ")
            warning.setFont(QFont("Arial", 11, QFont.Bold))
            warning.setAlignment(Qt.AlignCenter)
            warning.setStyleSheet(
                "color: white; background: #D32F2F; border-radius: 6px; padding: 10px; margin: 4px;"
            )
            content_layout.addWidget(warning)

        scroll.setWidget(content)
        layout.addWidget(scroll)

        # ── BUTTONS ──
        btn_layout = QHBoxLayout()

        cancel_btn = QPushButton("Cancel")
        cancel_btn.setFont(QFont("Arial", 10))
        cancel_btn.setMinimumHeight(32)
        cancel_btn.setStyleSheet(
            "background: #455A64; color: white; border-radius: 4px; padding: 6px 20px;"
        )
        cancel_btn.clicked.connect(self.reject)
        btn_layout.addWidget(cancel_btn)

        btn_layout.addStretch()

        pause_btn = QPushButton("⏸️ Pause Anyway")
        pause_btn.setFont(QFont("Arial", 10, QFont.Bold))
        pause_btn.setMinimumHeight(32)
        color = "#D32F2F" if is_essential else "#FF9800" if total_affected > 5 else "#4CAF50"
        pause_btn.setStyleSheet(
            f"background: {color}; color: white; border-radius: 4px; padding: 6px 20px;"
        )
        pause_btn.clicked.connect(self.accept)
        btn_layout.addWidget(pause_btn)

        layout.addLayout(btn_layout)
