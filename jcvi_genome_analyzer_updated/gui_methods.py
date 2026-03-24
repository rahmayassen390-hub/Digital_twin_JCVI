"""
GUI Methods for JCVI Genome Analyzer - Part 2
Contains the event handlers, data processing and animation methods
"""

import os
import json
from collections import defaultdict

from PyQt5.QtWidgets import QMessageBox, QFileDialog, QTableWidgetItem, QPushButton
from PyQt5.QtCore import QTimer, Qt
from PyQt5.QtGui import QColor

# Import paused gene dashboard
try:
    from paused_gene_dashboard import PausedGeneDashboard
    DASHBOARD_AVAILABLE = True
except ImportError:
    DASHBOARD_AVAILABLE = False
    print("Warning: paused_gene_dashboard not found")

# Import knockout effects modules
try:
    from gene_dependency_graph import GeneDependencyGraph, get_effect_color, get_effect_icon
    from knockout_effects_dashboard import KnockoutEffectsDashboard
    KNOCKOUT_EFFECTS_AVAILABLE = True
except ImportError as e:
    KNOCKOUT_EFFECTS_AVAILABLE = False
    print(f"Warning: Knockout effects modules not available: {e}")

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False

try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False

from workers import LoadGenomeWorker, BlastSearchWorker, AIAnalysisWorker, QuadGPUOrchestrationWorker
import math

# Import continuous simulation engine
try:
    from cell_simulation_engine import (
        CellSimulationEngine, SimulationMode, 
        create_engine_from_transcription_data,
        derive_resources_from_genome,
        validate_jcvi_organism,
        JCVI_MRNA_HALF_LIFE, JCVI_PROTEIN_HALF_LIFE,
        JCVI_DEFAULT_NOISE, ParameterSource
    )
    CONTINUOUS_ENGINE_AVAILABLE = True
except ImportError:
    CONTINUOUS_ENGINE_AVAILABLE = False
    print("Warning: cell_simulation_engine not found - continuous mode unavailable")

# Import MetabolicSolver for GUI FBA coupling
try:
    from metabolic_solver import MetabolicSolver
    METABOLIC_SOLVER_AVAILABLE = True
except ImportError:
    METABOLIC_SOLVER_AVAILABLE = False
    print("Warning: metabolic_solver not found - metabolic visualization unavailable")

# Import Enhanced Knockout Dialog for pre-pause analysis
try:
    from enhanced_knockout_dialog import EnhancedKnockoutDialog
    ENHANCED_KNOCKOUT_DIALOG_AVAILABLE = True
except ImportError:
    ENHANCED_KNOCKOUT_DIALOG_AVAILABLE = False
    print("Warning: enhanced_knockout_dialog not found")

# Import quantum worker
try:
    from quantum_worker import QuantumJobWorker, QuantumClusterWorker
    QUANTUM_AVAILABLE = True
except ImportError:
    QUANTUM_AVAILABLE = False
    print("Warning: quantum_worker not found - quantum simulation unavailable")

# Import cluster definitions for 28-qubit quantum model
try:
    from cluster_definitions import (
        CLUSTER_DEFINITIONS,
        ENTANGLEMENT_PAIRS,
        assign_gene_to_cluster,
        build_gene_cluster_mapping,
        aggregate_cluster_levels,
        get_cluster_gene_counts,
        get_cluster_name
    )
    CLUSTER_DEFINITIONS_AVAILABLE = True
except ImportError:
    CLUSTER_DEFINITIONS_AVAILABLE = False
    print("Warning: cluster_definitions not found - 28-qubit mode unavailable")


# =============================================================================
# JCVI syn3.0 PROTEIN VALIDATION CONSTANTS
# =============================================================================
# These thresholds are derived from JCVI syn3.0 minimal genome analysis
# Source: Hutchison et al. 2016, NCBI CP016816

# JCVI syn3.0 Protein Size Statistics (from 473 genes)
JCVI_SHORTEST_PROTEIN = 29       # Shortest annotated protein in JCVI (rpmJ-like)
JCVI_MEDIAN_PROTEIN = 189        # Median protein length in amino acids
JCVI_MEAN_PROTEIN = 254          # Mean protein length
JCVI_LONGEST_PROTEIN = 1524      # Longest protein (rpoC - RNA polymerase)

# Minimum thresholds based on JCVI biology
MIN_PROTEIN_LENGTH = 30          # Minimum for JCVI (based on shortest real protein)
JCVI_SUBUNIT_THRESHOLD = 60      # Proteins ≤60aa likely subunits/fragments
JCVI_ESSENTIAL_MIN = 50          # Most essential genes encode ≥50aa proteins

# Known ultra-short JCVI proteins (legitimate exceptions)
JCVI_SHORT_PROTEINS = {
    'rpmJ': 29,   # Ribosomal protein L36
    'rpmI': 33,   # Ribosomal protein L35  
    'rpsU': 35,   # Ribosomal protein S21
    'rpsT': 35,   # Ribosomal protein S20
    'csaA': 42,   # Chaperone-like
}

# Start/Stop codons (JCVI uses standard bacterial codons)
VALID_START_CODONS = {'ATG', 'GTG', 'TTG'}  # ATG is most common in JCVI
VALID_STOP_CODONS = {'TAA', 'TAG', 'TGA'}   # TAA is most common in JCVI

# JCVI-specific codon usage preference
JCVI_PREFERRED_START = 'ATG'     # >90% of JCVI genes
JCVI_PREFERRED_STOP = 'TAA'      # Most common stop in minimal genome

def validate_protein_sequence(cds_sequence: str, sequence_source: str = "unknown") -> dict:
    """
    Validate a CDS nucleotide sequence for biological validity.
    
    JCVI syn3.0 specific validation:
    - Minimum 30 amino acids (90 nucleotides) for valid protein
    - Ultra-short sequences (<10aa) are biologically impossible
    - Short proteins (30-60aa) flagged as potential subunits
    
    Checks:
    1. Minimum length (≥30 amino acids / 90 nucleotides)
    2. Valid start codon (ATG, GTG, or TTG for bacteria)
    3. Proper stop codon (TAA, TAG, or TGA)
    4. Length is multiple of 3 (complete codons)
    5. No internal stop codons (premature termination)
    
    Args:
        cds_sequence: CDS nucleotide sequence
        sequence_source: Origin of sequence ("bed", "gff", "cds", "computed", "unknown")
    
    Returns dict with validation results and classification
    """
    issues = []
    warnings = []
    cds_upper = cds_sequence.upper().strip() if cds_sequence else ""
    
    # === ULTRA-SHORT REJECTION ===
    # Sequences like "YC" (2aa) are biologically impossible
    if not cds_upper or len(cds_upper) < 9:  # Less than 3 codons
        return {
            'is_valid': False,
            'protein_length': len(cds_upper) // 3 if cds_upper else 0,
            'has_start_codon': False,
            'has_stop_codon': False,
            'issues': ['Ultra-short sequence - biologically impossible'],
            'warnings': [],
            'validity_score': 0.0,
            'size_category': 'invalid',
            'prediction_mode': 'none',
            'sequence_source': sequence_source,
            'is_subunit_candidate': False,
            'needs_review': False
        }
    
    # Length checks
    nucleotide_length = len(cds_upper)
    protein_length = nucleotide_length // 3
    
    # Check for complete codons
    has_complete_codons = nucleotide_length % 3 == 0
    if not has_complete_codons:
        issues.append(f'Incomplete codon (length {nucleotide_length} not divisible by 3)')
    
    # === JCVI-SPECIFIC SIZE CATEGORY CLASSIFICATION ===
    # Based on JCVI syn3.0 protein size distribution
    # Shortest real JCVI protein is rpmJ at 29aa
    if protein_length < 10:
        size_category = 'invalid'  # Biologically impossible for any organism
        prediction_mode = 'none'
        is_subunit = False
        issues.append(f'Impossibly short ({protein_length}aa < 10aa minimum)')
    elif protein_length < JCVI_SHORTEST_PROTEIN:  # <29aa (shorter than any JCVI protein)
        size_category = 'marginal'
        prediction_mode = 'flagged'
        is_subunit = True
        issues.append(f'Below JCVI minimum ({protein_length}aa < {JCVI_SHORTEST_PROTEIN}aa)')
        warnings.append('Shorter than any known JCVI syn3.0 protein')
    elif protein_length < MIN_PROTEIN_LENGTH:  # 29aa exactly (like rpmJ)
        size_category = 'minimal'
        prediction_mode = 'limited'
        is_subunit = True
        warnings.append(f'At JCVI minimum ({protein_length}aa) - like ribosomal protein L36')
    elif protein_length <= JCVI_SUBUNIT_THRESHOLD:  # 30-60aa
        size_category = 'short'
        prediction_mode = 'limited'  # Monomer-only prediction
        is_subunit = True
        warnings.append(f'Short JCVI protein ({protein_length}aa) - likely ribosomal subunit')
    elif protein_length <= 150:  # 60-150aa (small proteins in JCVI)
        size_category = 'small'
        prediction_mode = 'standard'
        is_subunit = False
    elif protein_length <= JCVI_MEAN_PROTEIN:  # 150-254aa (below average)
        size_category = 'below_median'
        prediction_mode = 'standard'
        is_subunit = False
    elif protein_length <= 500:  # 254-500aa (above average)
        size_category = 'medium'
        prediction_mode = 'standard'
        is_subunit = False
    elif protein_length <= JCVI_LONGEST_PROTEIN:  # 500-1524aa
        size_category = 'large'
        prediction_mode = 'standard'
        is_subunit = False
    else:  # >1524aa (longer than any JCVI protein!)
        size_category = 'enormous'
        prediction_mode = 'flagged'
        is_subunit = False
        warnings.append(f'Longer than any JCVI protein ({protein_length}aa > {JCVI_LONGEST_PROTEIN}aa)')
    
    # Start codon check with JCVI preference
    start_codon = cds_upper[:3] if len(cds_upper) >= 3 else ""
    has_start_codon = start_codon in VALID_START_CODONS
    if not has_start_codon:
        issues.append(f'Invalid start codon ({start_codon})')
    elif start_codon != JCVI_PREFERRED_START:
        warnings.append(f'Non-preferred start codon ({start_codon}) - JCVI typically uses {JCVI_PREFERRED_START}')
    
    # Stop codon check with JCVI preference
    stop_codon = cds_upper[-3:] if len(cds_upper) >= 3 else ""
    has_stop_codon = stop_codon in VALID_STOP_CODONS
    if not has_stop_codon:
        issues.append(f'Missing stop codon ({stop_codon})')
    
    # === INTERNAL STOP CODON CHECK ===
    # Check for premature termination
    internal_stops = 0
    if has_complete_codons and len(cds_upper) >= 6:
        for i in range(3, len(cds_upper) - 3, 3):  # Skip first and last codons
            codon = cds_upper[i:i+3]
            if codon in VALID_STOP_CODONS:
                internal_stops += 1
        if internal_stops > 0:
            issues.append(f'Internal stop codon(s) found ({internal_stops})')
    
    # Calculate validity score
    score = 0.0
    if protein_length >= MIN_PROTEIN_LENGTH:
        score += 0.4
    elif protein_length >= 10:
        score += 0.1  # Minimal credit for very short
    
    if has_start_codon:
        score += 0.25
    if has_stop_codon:
        score += 0.15
    if has_complete_codons:
        score += 0.1
    if internal_stops == 0:
        score += 0.1
    
    # Is fully valid?
    is_valid = (
        protein_length >= MIN_PROTEIN_LENGTH and
        has_start_codon and
        has_stop_codon and
        has_complete_codons and
        internal_stops == 0
    )
    
    # Needs review flag
    needs_review = (
        size_category in ['marginal', 'short'] or
        internal_stops > 0 or
        not has_start_codon or
        sequence_source == 'unknown'
    )
    
    return {
        'is_valid': is_valid,
        'protein_length': protein_length,
        'has_start_codon': has_start_codon,
        'has_stop_codon': has_stop_codon,
        'has_complete_codons': has_complete_codons,
        'internal_stop_count': internal_stops,
        'issues': issues,
        'warnings': warnings,
        'validity_score': round(score, 2),
        # === ENHANCED CLASSIFICATION ===
        'size_category': size_category,  # invalid/marginal/short/small/medium/large
        'prediction_mode': prediction_mode,  # none/flagged/limited/standard
        'sequence_source': sequence_source,
        'is_subunit_candidate': is_subunit,
        'needs_review': needs_review
    }


# === QUATERNARY STRUCTURE DATABASE (Issue 6) ===
# Known multimeric protein complexes in JCVI syn3.0 minimal cell
# Format: complex_name -> {subunits: [(gene_pattern, copy_number)], stoichiometry}

PROTEIN_COMPLEXES = {
    # Ribosome (70S = 30S + 50S)
    'ribosome_30S': {
        'name': 'Small ribosomal subunit (30S)',
        'subunit_patterns': ['rpsA', 'rpsB', 'rpsC', 'rpsD', 'rpsE', 'rpsF', 'rpsG', 'rpsH',
                            'rpsI', 'rpsJ', 'rpsK', 'rpsL', 'rpsM', 'rpsN', 'rpsO', 'rpsP',
                            'rpsQ', 'rpsR', 'rpsS', 'rpsT', 'rpsU'],
        'stoichiometry': '1:1:1:1... (21 proteins)',
        'assembly_order': 1,
        'parent_complex': 'ribosome_70S'
    },
    'ribosome_50S': {
        'name': 'Large ribosomal subunit (50S)',
        'subunit_patterns': ['rplA', 'rplB', 'rplC', 'rplD', 'rplE', 'rplF', 'rplI',
                            'rplJ', 'rplK', 'rplL', 'rplM', 'rplN', 'rplO', 'rplP',
                            'rplQ', 'rplR', 'rplS', 'rplT', 'rplU', 'rplV', 'rplW', 'rplX'],
        'stoichiometry': '1:1:1:1... (29 proteins)',
        'assembly_order': 2,
        'parent_complex': 'ribosome_70S'
    },
    'ribosome_70S': {
        'name': 'Complete ribosome (70S)',
        'subunit_patterns': [],  # Composed of 30S + 50S
        'stoichiometry': '30S:50S = 1:1',
        'assembly_order': 3,
        'subassemblies': ['ribosome_30S', 'ribosome_50S']
    },
    # ATP Synthase (F0F1)
    'atp_synthase_F0': {
        'name': 'ATP synthase F0 sector',
        'subunit_patterns': ['atpB', 'atpE', 'atpF'],
        'stoichiometry': 'a:b:c = 1:2:9-12',
        'assembly_order': 1,
        'parent_complex': 'atp_synthase'
    },
    'atp_synthase_F1': {
        'name': 'ATP synthase F1 sector',
        'subunit_patterns': ['atpA', 'atpD', 'atpG', 'atpH', 'atpC'],
        'stoichiometry': 'α:β:γ:δ:ε = 3:3:1:1:1',
        'assembly_order': 2,
        'parent_complex': 'atp_synthase'
    },
    'atp_synthase': {
        'name': 'Complete ATP synthase (F0F1)',
        'subunit_patterns': [],
        'stoichiometry': 'F0:F1 = 1:1',
        'assembly_order': 3,
        'subassemblies': ['atp_synthase_F0', 'atp_synthase_F1']
    },
    # DNA Polymerase III
    'dna_pol_III': {
        'name': 'DNA Polymerase III holoenzyme',
        'subunit_patterns': ['dnaE', 'dnaN', 'dnaQ', 'holA', 'holB'],
        'stoichiometry': 'α:β:ε... (10+ subunits)',
        'assembly_order': 1,
        'is_holoenzyme': True
    },
    # RNA Polymerase
    'rna_polymerase': {
        'name': 'RNA Polymerase core enzyme',
        'subunit_patterns': ['rpoA', 'rpoB', 'rpoC', 'rpoD', 'rpoE'],
        'stoichiometry': 'α2:β:β\':ω = 2:1:1:1',
        'assembly_order': 1,
        'is_holoenzyme': True
    },
    # Chaperonins
    'groel_groes': {
        'name': 'GroEL/GroES chaperonin',
        'subunit_patterns': ['groEL', 'groES', 'groL', 'groS'],
        'stoichiometry': 'GroEL14:GroES7 = 1:1-2',
        'assembly_order': 1,
    },
    # Membrane protein complexes
    'sec_translocase': {
        'name': 'Sec translocase',
        'subunit_patterns': ['secA', 'secY', 'secE', 'secG', 'secD', 'secF'],
        'stoichiometry': 'SecYEG:SecA = 1:1',
        'assembly_order': 1,
    }
}

def classify_protein_quaternary(gene_id: str, product: str = "") -> dict:
    """
    Classify a protein's quaternary structure based on gene ID and product.
    
    Returns:
        dict with:
        - quaternary_type: 'monomer', 'subunit', 'complex'
        - complex_name: Name of parent complex if subunit
        - subunit_role: Role description if subunit
        - assembly_partners: List of partner gene patterns
        - stoichiometry: Expected copy numbers
        - is_assembled: Whether this needs other subunits to function
    """
    gene_lower = gene_id.lower()
    product_lower = product.lower() if product else ""
    
    # Check each complex
    for complex_id, complex_info in PROTEIN_COMPLEXES.items():
        for pattern in complex_info.get('subunit_patterns', []):
            if pattern.lower() in gene_lower or pattern.lower() in product_lower:
                return {
                    'quaternary_type': 'subunit',
                    'complex_name': complex_info['name'],
                    'complex_id': complex_id,
                    'subunit_role': f"Subunit of {complex_info['name']}",
                    'assembly_partners': complex_info.get('subunit_patterns', []),
                    'stoichiometry': complex_info.get('stoichiometry', ''),
                    'assembly_order': complex_info.get('assembly_order', 0),
                    'parent_complex': complex_info.get('parent_complex', ''),
                    'is_assembled': False,  # Needs partners to function
                    'requires_partners': True
                }
    
    # Check for common subunit indicators in product description
    subunit_keywords = ['subunit', 'component', 'chain', 'domain']
    if any(kw in product_lower for kw in subunit_keywords):
        return {
            'quaternary_type': 'subunit',
            'complex_name': 'Unknown complex',
            'complex_id': '',
            'subunit_role': 'Subunit (complex unknown)',
            'assembly_partners': [],
            'stoichiometry': '',
            'assembly_order': 0,
            'parent_complex': '',
            'is_assembled': False,
            'requires_partners': True
        }
    
    # Default: monomeric protein (functions alone)
    return {
        'quaternary_type': 'monomer',
        'complex_name': '',
        'complex_id': '',
        'subunit_role': 'Functions independently',
        'assembly_partners': [],
        'stoichiometry': '1',
        'assembly_order': 0,
        'parent_complex': '',
        'is_assembled': True,  # Monomers are "assembled" by default
        'requires_partners': False
    }


# === GENE-SPECIFIC NOISE DERIVATION (Issue: Regulation-Driven Variability) ===
# Derives noise levels from regulatory BED data

# Noise mapping based on regulation type (JCVI minimal cell = low regulation)
REGULATION_NOISE_MAP = {
    'constitutive': 0.03,        # Minimal noise - always on, stable
    'housekeeping': 0.03,        # Essential genes - very stable
    'operonMember': 0.05,        # Part of operon - coordinated expression
    'unknown': 0.08,             # Unknown - moderate noise
    'potentiallyinducible': 0.15,  # Inducible - higher variability
    'potentiallyrepressible': 0.12,  # Repressible - moderate variability
    'inducible': 0.20,           # Strongly inducible - high noise
    'repressible': 0.15,         # Repressible - moderate-high noise
}

def derive_gene_noise_level(
    gene_id: str,
    regulation_type: str = "",
    regulation_strength: str = "",
    has_promoter: bool = True,
    is_operon_internal: bool = False,
    promoter_type: str = ""
) -> dict:
    """
    Derive gene-specific noise level from regulatory complexity.
    
    JCVI syn3.0 minimal cell has reduced regulation, so base noise is low.
    
    Args:
        gene_id: Gene identifier
        regulation_type: From regulation BED (constitutive, inducible, etc.)
        regulation_strength: high, medium, low
        has_promoter: Whether gene has its own promoter
        is_operon_internal: True if gene is internal to an operon
        promoter_type: Type of promoter (BOTH, ONE, etc.)
    
    Returns:
        dict with noise_level (0-1), noise_source, and explanation
    """
    # Base noise from regulation type
    reg_lower = regulation_type.lower() if regulation_type else 'unknown'
    base_noise = REGULATION_NOISE_MAP.get(reg_lower, 0.08)
    
    explanation = []
    
    # Modify based on regulation strength
    strength_modifier = 1.0
    if regulation_strength:
        strength_lower = regulation_strength.lower()
        if strength_lower == 'high':
            strength_modifier = 0.8   # Strong regulation = lower noise
            explanation.append("strong regulation (-20%)")
        elif strength_lower == 'low':
            strength_modifier = 1.3   # Weak regulation = higher noise
            explanation.append("weak regulation (+30%)")
    
    # Operon-internal genes have coordinated expression = lower noise
    if is_operon_internal or not has_promoter:
        base_noise *= 0.7  # 30% reduction
        explanation.append("operon-internal (-30%)")
    
    # Promoter complexity
    if promoter_type:
        if promoter_type.upper() == 'BOTH':
            # Both -10 and -35 boxes = well-regulated = lower noise
            base_noise *= 0.85
            explanation.append("complete promoter (-15%)")
        elif promoter_type.upper() == 'NONE':
            # Missing promoter elements = higher noise
            base_noise *= 1.2
            explanation.append("weak promoter (+20%)")
    
    final_noise = base_noise * strength_modifier
    
    # Clamp to valid range
    final_noise = max(0.01, min(0.5, final_noise))
    
    return {
        'noise_level': round(final_noise, 3),
        'regulation_type': reg_lower,
        'noise_source': 'regulation_bed',
        'explanation': ', '.join(explanation) if explanation else 'base regulation type'
    }


# === DYNAMIC STRUCTURE PREDICTION THRESHOLDS (Issue 4) ===
# Compute thresholds from protein abundance distributions in BED data

def compute_structure_prediction_thresholds(translation_dynamics: dict) -> dict:
    """
    Compute dynamic thresholds for structure prediction eligibility
    based on protein abundance distributions in translation BED data.
    
    Instead of fixed 10% threshold, derives thresholds from:
    - Protein level distribution (percentiles)
    - Expression efficiency distribution
    - RBS strength distribution
    
    Args:
        translation_dynamics: Dict of gene_id -> TranslationDynamics from BED
    
    Returns:
        dict with:
        - min_protein_fraction: Fraction threshold for prediction (derived)
        - min_protein_level: Absolute level threshold
        - confidence_percentile: What percentile this represents
        - threshold_source: 'expression_data' or 'default'
    """
    if not translation_dynamics:
        # No data - use defaults
        return {
            'min_protein_fraction': 0.1,  # Default 10%
            'min_protein_level': 10.0,
            'confidence_percentile': 10,
            'threshold_source': 'default',
            'distribution_stats': {}
        }
    
    # Collect protein levels from all genes
    protein_levels = []
    efficiencies = []
    
    for gene_id, data in translation_dynamics.items():
        if hasattr(data, 'protein_level') and data.protein_level > 0:
            protein_levels.append(data.protein_level)
        if hasattr(data, 'efficiency') and data.efficiency > 0:
            efficiencies.append(data.efficiency)
    
    if not protein_levels:
        return {
            'min_protein_fraction': 0.1,
            'min_protein_level': 10.0,
            'confidence_percentile': 10,
            'threshold_source': 'default',
            'distribution_stats': {}
        }
    
    # Sort for percentile calculations
    protein_levels.sort()
    n = len(protein_levels)
    
    # Compute distribution statistics
    mean_level = sum(protein_levels) / n
    median_level = protein_levels[n // 2]
    min_level = protein_levels[0]
    max_level = protein_levels[-1]
    
    # Compute percentiles
    p10_idx = max(0, int(n * 0.10) - 1)
    p25_idx = max(0, int(n * 0.25) - 1)
    p75_idx = min(n - 1, int(n * 0.75))
    
    p10_level = protein_levels[p10_idx]
    p25_level = protein_levels[p25_idx]
    p75_level = protein_levels[p75_idx]
    
    # === DYNAMIC THRESHOLD LOGIC ===
    # Use the 10th percentile as the minimum for structure prediction
    # This means only genes with protein levels in top 90% are eligible
    # Adjusts automatically based on the actual expression distribution
    
    # Convert percentile level to fraction of max
    dynamic_fraction = p10_level / max_level if max_level > 0 else 0.1
    
    # Ensure reasonable bounds (at least 5%, no more than 30%)
    dynamic_fraction = max(0.05, min(0.30, dynamic_fraction))
    
    distribution_stats = {
        'n_genes': n,
        'mean_level': round(mean_level, 1),
        'median_level': round(median_level, 1),
        'min_level': round(min_level, 1),
        'max_level': round(max_level, 1),
        'p10_level': round(p10_level, 1),
        'p25_level': round(p25_level, 1),
        'p75_level': round(p75_level, 1)
    }
    
    print(f"DEBUG: Computed dynamic structure thresholds from {n} genes:")
    print(f"  - Distribution: min={min_level:.0f}%, median={median_level:.0f}%, max={max_level:.0f}%")
    print(f"  - 10th percentile: {p10_level:.1f}% → threshold fraction: {dynamic_fraction:.2f}")
    
    return {
        'min_protein_fraction': round(dynamic_fraction, 3),
        'min_protein_level': p10_level,
        'confidence_percentile': 10,
        'threshold_source': 'expression_data',
        'distribution_stats': distribution_stats
    }



class GenomeAnalyzerMethods:
    def verify_all_paths(self, is_query=False):
        """Verify that all selected paths exist and are files, not directories"""
        missing = []
        invalid_dirs = []
        
        if not is_query:
            paths = {
                "FASTA": self.fasta_path,
                "GFF3": self.gff_path,
                "Promoters": self.likely_promoters_path,
                "Operons": self.internal_operons_path,
                "Complex Cases": getattr(self, 'ref_complex_cases_path', None),
                "Regulation": getattr(self, 'ref_regulation_path', None),
                "Transcription": getattr(self, 'ref_transcription_path', None),
                "Translation": getattr(self, 'ref_translation_path', None)
            }
        else:
            paths = {
                "FASTA": self.query_fasta_path,
                "GFF3": self.query_gff_path,
                "Promoters": self.query_likely_promoters_path,
                "Operons": self.query_internal_operons_path,
                "Complex Cases": getattr(self, 'query_complex_cases_path', None),
                "Regulation": getattr(self, 'query_regulation_path', None),
                "Transcription": getattr(self, 'query_transcription_path', None),
                "Translation": getattr(self, 'query_translation_path', None)
            }
            
        for name, path in paths.items():
            if path:
                if not os.path.exists(path):
                    missing.append(f"{name}: {path}")
                elif os.path.isdir(path):
                    invalid_dirs.append(f"{name}: {path}")
                    
        return missing, invalid_dirs
    """Mixin class with all the methods for GenomeAnalyzerApp"""
    
    # Loading methods
    def load_reference_genome(self):
        if not self.fasta_path or not self.gff_path:
            QMessageBox.warning(self, "Missing Files", "Please select both FASTA and GFF3 files")
            return
        self.ref_progress_bar.setVisible(True)
        self.load_ref_btn.setEnabled(False)
        
        self.ref_worker = LoadGenomeWorker(
            self.fasta_path, self.gff_path,
            self.genome_manager, self.blast_manager,
            self.likely_promoters_path, self.internal_operons_path,
            self.ref_complex_cases_path if hasattr(self, 'ref_complex_cases_path') else None,
            self.ref_regulation_path if hasattr(self, 'ref_regulation_path') else None,
            self.ref_transcription_path if hasattr(self, 'ref_transcription_path') else None,
            self.ref_translation_path if hasattr(self, 'ref_translation_path') else None
        )
        
        # Verify paths before starting
        missing, invalid_dirs = self.verify_all_paths(is_query=False)
        if missing or invalid_dirs:
            error_msg = "Cannot load reference genome due to path errors:\n\n"
            if missing:
                error_msg += "Missing Files:\n- " + "\n- ".join(missing) + "\n\n"
            if invalid_dirs:
                error_msg += "Directories selected instead of files:\n- " + "\n- ".join(invalid_dirs) + "\n\n"
            error_msg += "Please use the Browse buttons to select the correct files."
            
            QMessageBox.critical(self, "Path Error", error_msg)
            self.ref_progress_bar.setVisible(False)
            self.load_ref_btn.setEnabled(True)
            return

        self.ref_worker.progress.connect(self.update_ref_progress)
        self.ref_worker.finished.connect(self.ref_load_finished)
        self.ref_worker.start()
    
    def load_query_genome(self):
        if not self.query_fasta_path:
            QMessageBox.warning(self, "Missing Files", "Please select query FASTA file")
            return
        
        self.query_progress_bar.setVisible(True)
        self.load_query_btn.setEnabled(False)
        
        gff_path = self.query_gff_path if self.query_gff_path else self.gff_path
        
        # Verify paths before starting
        missing, invalid_dirs = self.verify_all_paths(is_query=True)
        if missing or invalid_dirs:
            error_msg = "Cannot load query genome due to path errors:\n\n"
            if missing:
                error_msg += "Missing Files:\n- " + "\n- ".join(missing) + "\n\n"
            if invalid_dirs:
                error_msg += "Directories selected instead of files:\n- " + "\n- ".join(invalid_dirs) + "\n\n"
            error_msg += "Please use the Browse buttons to select the correct files."
            
            QMessageBox.critical(self, "Path Error", error_msg)
            self.query_progress_bar.setVisible(False)
            self.load_query_btn.setEnabled(True)
            return

        self.query_worker = LoadGenomeWorker(
            self.query_fasta_path, gff_path,
            self.query_genome_manager, None,  # No BLAST for query
            self.query_likely_promoters_path, self.query_internal_operons_path,
            self.query_complex_cases_path if hasattr(self, 'query_complex_cases_path') else None,
            self.query_regulation_path if hasattr(self, 'query_regulation_path') else None,
            self.query_transcription_path if hasattr(self, 'query_transcription_path') else None,
            self.query_translation_path if hasattr(self, 'query_translation_path') else None
        )
        self.query_worker.progress.connect(self.update_query_progress)
        self.query_worker.finished.connect(self.query_load_finished)
        self.query_worker.start()
    
    def update_ref_progress(self, value, message):
        self.ref_progress_bar.setValue(value)
        self.ref_status_label.setText(message)
    
    def update_query_progress(self, value, message):
        self.query_progress_bar.setValue(value)
        self.query_status_label.setText(message)
    
    def ref_load_finished(self, success, message, stats):
        self.ref_progress_bar.setVisible(False)
        self.load_ref_btn.setEnabled(True)
        
        if success:
            self.ref_status_label.setText("✅ Reference genome loaded!")
            self.ref_status_label.setStyleSheet("color: green;")
            
            stats_text = f"""Reference Genome Statistics:

Sequence Name: {stats['seq_name']}
Sequence Length: {stats['seq_len']:,} bp
Total Genes: {stats['gene_count']:,}

Annotations Loaded:
  • Likely Promoters: {stats.get('likely_promoters_count', 0):,}
  • Internal Operons: {stats.get('internal_operons_count', 0):,}
  • Regulation Data: {stats.get('regulation_count', 0):,} genes
  • Transcription Data: {stats.get('transcription_count', 0):,} genes
  • Translation Data: {stats.get('translation_count', 0):,} genes

✅ BLAST database created successfully!
"""
            self.ref_stats_text.setText(stats_text)
            self.tabs.setTabEnabled(1, True)
            
            # AUTOMATION: Trigger JCVI 3.0 (Query) load after 1.0 (Reference) is done
            if hasattr(self, 'query_fasta_path') and self.query_fasta_path:
                print(f"🔄 JCVI 1.0 loaded. Automatically starting JCVI 3.0 load...")
                self.load_query_genome()
        else:
            QMessageBox.critical(self, "Error", f"Failed to load reference:\n{message}")
    
    def query_load_finished(self, success, message, stats):
        self.query_progress_bar.setVisible(False)
        self.load_query_btn.setEnabled(True)
        
        if success:
            self.query_status_label.setText("✅ Query genome loaded!")
            self.query_status_label.setStyleSheet("color: green;")
            
            stats_text = f"""Query Genome Statistics:

Sequence Name: {stats['seq_name']}
Sequence Length: {stats['seq_len']:,} bp
Total Genes: {stats.get('gene_count', 0):,}

Annotations Loaded:
  • Likely Promoters: {stats.get('likely_promoters_count', 0):,}
  • Internal Operons: {stats.get('internal_operons_count', 0):,}
  • Regulation Data: {stats.get('regulation_count', 0):,} genes
  • Transcription Data: {stats.get('transcription_count', 0):,} genes
  • Translation Data: {stats.get('translation_count', 0):,} genes

✅ Ready for analysis!
"""
            self.query_stats_text.setText(stats_text)
            self.tabs.setTabEnabled(2, True)
            
            # Enable AI Results tab if BED data is loaded
            has_bed_data = (
                stats.get('regulation_count', 0) > 0 or 
                stats.get('transcription_count', 0) > 0 or 
                stats.get('translation_count', 0) > 0
            )
            
            if has_bed_data:
                self.tabs.setTabEnabled(4, True)
                
                # Prepare translation animation data if available
                if stats.get('translation_count', 0) > 0:
                    self.prepare_translation_animation_data()
        else:
            QMessageBox.critical(self, "Error", f"Failed to load query:\n{message}")
    
    def update_gpu_thermals(self):
        """Update GPU thermal indicators and HPC monitor metrics"""
        try:
            metrics = self.gpu_orchestrator.get_hpc_metrics()
            
            # 1. Update Sidebar Monitor (if available)
            if hasattr(self, 'hpc_monitor') and self.hpc_monitor:
                self.hpc_monitor.update_all(metrics)
            
            # 2. Update Status Bar (legacy compatibility)
            for gpu in metrics.get('gpus', []):
                idx = gpu['index']
                temp = gpu['temp']
                if idx in self.gpu_temp_labels:
                    label = self.gpu_temp_labels[idx]
                    label.setText(f"GPU {idx}: {temp}°C")
                    
                    # Update colors based on temperature thresholds
                    if temp < 60:
                        label.setStyleSheet("color: #2E7D32; font-weight: bold;") # Green
                    elif temp < 80:
                        label.setStyleSheet("color: #F57C00; font-weight: bold;") # Orange
                    else:
                        label.setStyleSheet("color: #D32F2F; font-weight: bold;") # Red
                        
        except Exception as e:
            # Silently fail to avoid interrupting the user if something goes wrong with polling
            print(f"DEBUG: HPC metrics poll failed: {e}")

    def run_blast_search(self):
        if not self.query_genome_manager.sequence:
            QMessageBox.warning(self, "No Query", "Please load query genome first")
            return
        
        self.search_progress_bar.setVisible(True)
        self.search_btn.setEnabled(False)
        
        self.blast_worker = BlastSearchWorker(
            self.query_genome_manager.sequence,
            self.genome_manager,
            self.blast_manager,
            self.query_genome_manager
        )
        self.blast_worker.progress.connect(self.update_search_progress)
        self.blast_worker.finished.connect(self.search_finished)
        self.blast_worker.start()
    
    def update_search_progress(self, value, message):
        self.search_progress_bar.setValue(value)
        self.search_status_label.setText(message)
    
    def search_finished(self, success, message, results):
        self.search_progress_bar.setVisible(False)
        self.search_btn.setEnabled(True)
        
        if success:
            self.search_results = results
            self.display_results(results)
            self.tabs.setTabEnabled(3, True)
            self.tabs.setCurrentIndex(3)
        else:
            QMessageBox.critical(self, "Error", f"BLAST search failed:\n{message}")
    
    def display_results(self, results):
        identical = sum(1 for r in results if r.identity >= 99.9)
        high = sum(1 for r in results if 90 <= r.identity < 99.9)
        
        summary = f"🔬 Query Genome (syn3.0) Analysis | Total Hits: {len(results)} | Perfect Match: {identical} | High Similarity (>90%): {high}"
        self.summary_label.setText(summary)
        
        self.results_table.setRowCount(len(results))
        
        for i, result in enumerate(results):
            self.results_table.setItem(i, 0, QTableWidgetItem(f"{result.identity:.2f}"))
            
            query_pos = f"Q:{result.query_start}-{result.query_end}"
            ref_pos = f"R:{result.target_start}-{result.target_end}"
            self.results_table.setItem(i, 1, QTableWidgetItem(f"{query_pos} | {ref_pos}"))
            
            self.results_table.setItem(i, 2, QTableWidgetItem(f"{result.evalue:.2e}"))
            self.results_table.setItem(i, 3, QTableWidgetItem(", ".join(result.genes[:3]) if result.genes else "No genes"))
            self.results_table.setItem(i, 4, QTableWidgetItem(result.gene_functions[0] if result.gene_functions else ""))
            self.results_table.setItem(i, 5, QTableWidgetItem(f"{len(result.protein_sequence)} aa" if result.protein_sequence else ""))
            self.results_table.setItem(i, 6, QTableWidgetItem("Click for details"))
            
            if result.identity >= 99.9:
                color = QColor(200, 255, 200)
            elif result.identity >= 90:
                color = QColor(255, 255, 200)
            else:
                color = QColor(255, 230, 230)
            
            for j in range(7):
                self.results_table.item(i, j).setBackground(color)
        
        self.display_gene_summary(results)
    
    def display_gene_summary(self, results):
        if not self.query_genome_manager or not self.query_genome_manager.genes:
            return
        
        gene_stats = {}
        
        for gene in self.query_genome_manager.genes:
            gene_stats[gene.id] = {
                'name': gene.name,
                'function': gene.function,
                'best_identity': 0.0,
                'total_coverage': 0,
                'hit_count': 0,
                'covered_positions': set()
            }
        
        for result in results:
            for gene in self.query_genome_manager.genes:
                overlap_start = max(result.query_start, gene.start)
                overlap_end = min(result.query_end, gene.end)
                overlap_len = max(0, overlap_end - overlap_start + 1)
                
                if overlap_len > 0:
                    gene_stats[gene.id]['hit_count'] += 1
                    gene_stats[gene.id]['best_identity'] = max(
                        gene_stats[gene.id]['best_identity'], 
                        result.identity
                    )
                    for pos in range(overlap_start, overlap_end + 1):
                        gene_stats[gene.id]['covered_positions'].add(pos)
        
        gene_rows = []
        for gene_id, stats in gene_stats.items():
            gene = next((g for g in self.query_genome_manager.genes if g.id == gene_id), None)
            if gene:
                gene_length = abs(gene.end - gene.start) + 1
                coverage = len(stats['covered_positions']) / gene_length * 100 if gene_length > 0 else 0
                
                if stats['best_identity'] >= 99.9:
                    status = "✅ Conserved"
                elif stats['best_identity'] >= 90:
                    status = "⚠️ Modified"
                elif stats['best_identity'] > 0:
                    status = "🔄 Diverged"
                else:
                    status = "❌ Not found"
                
                gene_rows.append({
                    'id': gene_id,
                    'name': stats['name'],
                    'function': stats['function'],
                    'identity': stats['best_identity'],
                    'coverage': coverage,
                    'status': status
                })
        
        gene_rows.sort(key=lambda x: x['identity'], reverse=True)
        
        self.gene_summary_table.setRowCount(len(gene_rows))
        
        for i, row in enumerate(gene_rows):
            self.gene_summary_table.setItem(i, 0, QTableWidgetItem(row['id']))
            self.gene_summary_table.setItem(i, 1, QTableWidgetItem(row['name']))
            self.gene_summary_table.setItem(i, 2, QTableWidgetItem(row['function']))
            self.gene_summary_table.setItem(i, 3, QTableWidgetItem(f"{row['identity']:.2f}"))
            self.gene_summary_table.setItem(i, 4, QTableWidgetItem(f"{row['coverage']:.1f}%"))
            self.gene_summary_table.setItem(i, 5, QTableWidgetItem(row['status']))
            
            if row['identity'] >= 99.9:
                color = QColor(200, 255, 200)
            elif row['identity'] >= 90:
                color = QColor(255, 255, 200)
            elif row['identity'] > 0:
                color = QColor(255, 230, 230)
            else:
                color = QColor(220, 220, 220)
            
            for j in range(6):
                self.gene_summary_table.item(i, j).setBackground(color)
    
    def show_result_details(self, row, column):
        if row >= len(self.search_results):
            return
        
        result = self.search_results[row]
        
        details = f"""
=== HIT #{row + 1} - Query Genome (syn3.0) Analysis ===

🎯 Sequence Identity: {result.identity:.2f}%
   - Alignment Length: {result.alignment_length} bp
   - Mismatches: {result.mismatches}
   - Gap Opens: {result.gap_opens}

📍 Query Position (syn3.0): {result.query_start:,} - {result.query_end:,} bp
📍 Reference Position (syn1.0): {result.target_start:,} - {result.target_end:,} bp

🧬 Query Genes (syn3.0): {', '.join(result.genes) if result.genes else 'No overlapping genes'}

⚙️ Gene Functions:
{chr(10).join('  • ' + f for f in result.gene_functions) if result.gene_functions else '  • No function annotations'}

📊 Statistical Metrics:
   - E-value: {result.evalue:.2e}
   - Bit Score: {result.bit_score:.1f}

🧪 Protein Sequence ({len(result.protein_sequence)} amino acids):
{result.protein_sequence[:200]}{'...' if len(result.protein_sequence) > 200 else ''}

🔬 Sequence Alignment (first 80 bp):
Query (syn3.0):  {result.query_seq[:80]}
Target (syn1.0): {result.target_seq[:80]}
        """
        
        self.details_text.setText(details)
    
    def export_results(self):
        if not self.search_results:
            QMessageBox.warning(self, "No Results", "No search results to export")
            return
        
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Results", "blast_results.csv", "CSV Files (*.csv)"
        )
        
        if file_path and PANDAS_AVAILABLE:
            df_data = []
            for r in self.search_results:
                df_data.append({
                    'Identity_%': round(r.identity, 2),
                    'Query_Start': r.query_start,
                    'Query_End': r.query_end,
                    'Target_Start': r.target_start,
                    'Target_End': r.target_end,
                    'E-value': r.evalue,
                    'Bit_Score': r.bit_score,
                    'Genes': ', '.join(r.genes),
                    'Functions': ' | '.join(r.gene_functions[:2]) if r.gene_functions else ""
                })
            
            df = pd.DataFrame(df_data)
            df.to_csv(file_path, index=False)
            
            QMessageBox.information(self, "Export Complete", f"Results exported to:\n{file_path}")
    
    def run_ai_analysis(self):
        # Check if we have any data loaded
        has_data = (
            (self.query_genome_manager and self.query_genome_manager.genes) or
            (self.query_genome_manager and self.query_genome_manager.regulation_data) or
            (self.query_genome_manager and self.query_genome_manager.transcription_dynamics) or
            (self.query_genome_manager and self.query_genome_manager.translation_dynamics)
        )
        
        if not has_data:
            QMessageBox.warning(self, "No Data", 
                "No data loaded!\n\n"
                "Please load Query genome files first:\n"
                "• GFF3 (for gene functions)\n"
                "• BED files (Regulation, Transcription, Translation)")
            return
        
        # Collect ALL file paths for smart engine
        all_file_paths = {
            'ref_gff': self.gff_path,
            'query_gff': self.query_gff_path,
            'query_regulation': getattr(self, 'query_regulation_path', None),
            'query_transcription': getattr(self, 'query_transcription_path', None),
            'query_translation': getattr(self, 'query_translation_path', None),
            'query_promoters': getattr(self, 'query_promoters_path', None),
            'query_operons': getattr(self, 'query_operons_path', None),
        }
        
        self.ai_worker = AIAnalysisWorker(
            self.genome_manager,
            self.search_results if self.search_results else [],
            self.query_genome_manager if self.query_genome_manager.genes else None,
            gff_path=self.gff_path,
            query_gff_path=self.query_gff_path,
            all_file_paths=all_file_paths
        )
        self.ai_worker.progress.connect(self.update_search_progress)
        self.ai_worker.finished.connect(self.ai_analysis_finished)
        self.ai_worker.start()
        
        # === START QUAD-GPU ORCHESTRATION ===
        # This will engage all 4 specialized GPU nodes simultaneously
        try:
            self.gpu_worker = QuadGPUOrchestrationWorker()
            # Link node status to a console or specific UI element if needed
            self.gpu_worker.node_status.connect(lambda name, status, dev: 
                print(f"🚀 GPU NODE [{name.upper()}]: {status} on {dev}"))
            self.gpu_worker.start()
        except Exception as e:
            print(f"Warning: Could not start Quad-GPU worker: {e}")
        
        self.search_progress_bar.setVisible(True)
    
    def ai_analysis_finished(self, success, message, results):
        self.search_progress_bar.setVisible(False)
        
        if success:
            self.ai_results = results
            self.display_ai_results(results)
            self.tabs.setTabEnabled(4, True)
            self.tabs.setCurrentIndex(4)
        else:
            QMessageBox.critical(self, "Error", f"AI analysis failed:\n{message}")
    
    def display_ai_results(self, results):
        """Display AI results - uses integrated genes from Smart Data Engine"""
        
        # Get integrated genes from Smart Engine
        integrated_genes = results.get('integrated_genes', {})
        stats = results.get('statistics', {})
        
        # Store for later use
        self.integrated_genes = integrated_genes
        
        # Get data from genome_manager (still needed for some operations)
        regulation_data = self.query_genome_manager.regulation_data if self.query_genome_manager else {}
        transcription_data = self.query_genome_manager.transcription_dynamics if self.query_genome_manager else {}
        translation_data = self.query_genome_manager.translation_dynamics if self.query_genome_manager else {}
        
        # === Display Regulation Data with SMART categories ===
        if integrated_genes or regulation_data:
            # Use integrated genes if available, otherwise fall back to regulation_data
            display_data = integrated_genes if integrated_genes else regulation_data
            
            # Limit to first 1000 for display
            limit = min(1000, len(display_data) if display_data else 0)
            display_data = dict(list(display_data.items())[:limit]) if display_data else {}
            
            self.promoter_table.setRowCount(len(display_data))
            self.promoter_table.setColumnCount(8)
            self.promoter_table.setHorizontalHeaderLabels([
                "Gene ID", "Regulation Type", "Category", "Function", 
                "Strength", "TX Level", "Protein", "Sources"
            ])
            
            for i, (gene_id, gene) in enumerate(display_data.items()):
                self.promoter_table.setItem(i, 0, QTableWidgetItem(gene_id))
                
                # Regulation type
                reg_type = getattr(gene, 'regulation_type', '')
                self.promoter_table.setItem(i, 1, QTableWidgetItem(reg_type))
                
                # Category - use function_category if regulation_category is Unknown
                category = getattr(gene, 'regulation_category', '')
                if category in ['', 'Unknown']:
                    category = getattr(gene, 'function_category', '') or 'Unknown'
                self.promoter_table.setItem(i, 2, QTableWidgetItem(category))
                
                # Function (product)
                product = getattr(gene, 'product', '')
                func_short = product[:45] + "..." if len(product) > 45 else product
                self.promoter_table.setItem(i, 3, QTableWidgetItem(func_short))
                
                # Strength
                strength = getattr(gene, 'regulation_strength', '')
                self.promoter_table.setItem(i, 4, QTableWidgetItem(strength))
                
                # TX Level
                tx_level = getattr(gene, 'transcription_level', 0)
                self.promoter_table.setItem(i, 5, QTableWidgetItem(f"{tx_level:.0f}%" if tx_level else ""))
                
                # Protein Level
                protein = getattr(gene, 'protein_level', 0)
                self.promoter_table.setItem(i, 6, QTableWidgetItem(f"{protein:.0f}%" if protein else ""))
                
                # Data Sources
                sources = getattr(gene, 'data_sources', [])
                self.promoter_table.setItem(i, 7, QTableWidgetItem(", ".join(sources) if sources else ""))
                
                # Color by regulation type
                if reg_type == 'constitutive':
                    color = QColor(200, 255, 200)  # Green
                elif reg_type == 'operonMember':
                    color = QColor(255, 220, 200)  # Orange
                elif reg_type == 'potentiallyInducible':
                    color = QColor(200, 220, 255)  # Blue
                elif reg_type == 'potentiallyRepressible':
                    color = QColor(255, 200, 200)  # Red
                else:
                    color = QColor(240, 240, 240)  # Gray
                
                for j in range(8):
                    if self.promoter_table.item(i, j):
                        self.promoter_table.item(i, j).setBackground(color)
            
            print(f"DEBUG: Displayed {len(display_data)} genes with SMART integration")
        
        # === Prepare Transcription Animation Data ===
        self.prepare_transcription_animation_data_from_bed()
        
        # === Prepare Translation Animation Data ===
        if translation_data:
            self.prepare_translation_animation_data()
        
        # === Generate Comparison Report ===
        comparison_text = self.generate_comparison_report_from_bed()
        
        # Add engine report if available
        engine_report = results.get('engine_report', '')
        if engine_report:
            comparison_text += "\n\n" + engine_report
        else:
            comparison_text += self._generate_files_used_report()
        
        self.comparison_text.setText(comparison_text)
        
        # Show statistics
        unknown_count = stats.get('unknown_count', 0)
        total = stats.get('total_genes', len(integrated_genes))
        known = total - unknown_count - stats.get('hypothetical_count', 0)
        
        QMessageBox.information(self, "🧠 SMART Integration Complete", 
            f"✅ Data Integrated from All Sources!\n\n"
            f"📊 Statistics:\n"
            f"• Total Genes: {total}\n"
            f"• Known Function: {known}\n"
            f"• Hypothetical: {stats.get('hypothetical_count', 0)}\n"
            f"• Unknown: {unknown_count}\n"
            f"• ncRNA: {stats.get('ncrna_count', 0)}\n\n"
            f"📁 Files Used: {len(stats.get('file_counts', {}))}\n\n"
            f"🎬 Transcription window is ready!\n"
            f"🔬 Translation window will open automatically!"
        )
        
        # Auto-open translation window if we have translation data
        if translation_data and self.gene_translation_data:
            self.translation_window.show()
            print("DEBUG: Translation window opened automatically")
    
    def _generate_files_used_report(self) -> str:
        """Generate report of which files are used"""
        return """

=== FILES USED FOR EACH DATA TYPE ===

📁 Promoter/Regulation Analysis:
   • Source: Regulation Analysis BED (JCVI_syn3_Gene_Regulation_3_0_.bed)
   • Fields: regulation_type, category, strength
   • Category "Unknown" → Replaced with function from GFF3

📁 Gene Functions:
   • Source: GFF3 file (JCVI_3_0.gff3)
   • Used to: Get gene product/function descriptions
   • Replaces "Unknown" categories

📁 Transcription Dynamics:
   • Source: Transcription Dynamics BED (JCVI_Transcription_Dynamics_3_0_.bed)
   • Track 1 - Transcription_Level: Level percentage
   • Track 2 - Transcription_Order: Phase, Category, Order
   • Track 3 - Transcription_Timeline: Time range
   
   PHASE COLORS:
   • 🔴 Phase 1-2 (Early) - Essential genes first
   • 🟡 Phase 3-4 (Middle) 
   • 🟢 Phase 5 (Late)
   • 🔵 Phase 6 (Final)

📁 Translation Dynamics:
   • Source: Translation Dynamics BED (JCVI_Translation_Dynamics_3_0_.bed)
   • Track 1 - Translation_Level: Protein percentage
   • Track 2 - Translation_Order: Phase, Category
   • Track 3 - RBS_Strength: RBS strength, Stability
   • Track 4 - Protein_Stability: Half-life
   • Track 5 - ncRNA_No_Translation: Non-coding RNAs
   • Track 6 - TX_vs_TL: Transcription vs Translation efficiency
   
   EXPRESSION STATUS:
   • 🟢 High - Protein level >= 80%
   • 🟡 Low Translation - TX high but TL low (Red in BED = 255,0,0)
   • 🔴 No Translation (ncRNA) - Orange in BED (255,69,0)
   • 🟠 Low - Protein level < 50%
"""
    
    def prepare_transcription_animation_data_from_bed(self):
        """Prepare transcription animation data directly from BED file"""
        self.gene_transcription_data = []
        self.paused_genes.clear()  # Clear paused genes on data reload
        self.gene_pause_buttons.clear()
        
        transcription_data = self.query_genome_manager.transcription_dynamics
        regulation_data = self.query_genome_manager.regulation_data
        
        if not transcription_data:
            print("DEBUG: No transcription data available")
            return
        
        # Get gene functions from stored data
        gene_functions = getattr(self, 'gene_functions', {})
        
        # Sort by order (if available) then by level
        sorted_data = sorted(
            transcription_data.items(),
            key=lambda x: (x[1].order if x[1].order > 0 else 9999, -x[1].transcription_level)
        )
        
        for idx, (gene_id, data) in enumerate(sorted_data):
            order_num = data.order if data.order > 0 else idx + 1
            start_time = idx * 0.3
            rate = 0.1 + (data.transcription_level / 100.0) * 0.2
            
            # Get regulation type if available
            reg_type = ""
            if regulation_data and gene_id in regulation_data:
                reg_type = regulation_data[gene_id].regulation_type
            
            # Get gene function from GFF3
            gene_function = gene_functions.get(gene_id, "")
            func_short = gene_function[:40] + "..." if len(gene_function) > 40 else gene_function
            
            # Determine phase color and status based on Phase number
            phase = data.phase or ""
            phase_num = 0
            if "Phase1" in phase:
                phase_num = 1
            elif "Phase2" in phase:
                phase_num = 2
            elif "Phase3" in phase:
                phase_num = 3
            elif "Phase4" in phase:
                phase_num = 4
            elif "Phase5" in phase:
                phase_num = 5
            elif "Phase6" in phase:
                phase_num = 6
            
            # Phase status with emoji
            if phase_num in [1, 2]:
                phase_status = f"🔴 {phase} (Early)"
            elif phase_num in [3, 4]:
                phase_status = f"🟡 {phase} (Middle)"
            elif phase_num == 5:
                phase_status = f"🟢 {phase} (Late)"
            elif phase_num == 6:
                phase_status = f"🔵 {phase} (Final)"
            else:
                phase_status = phase if phase else "Unknown"
            
            # Expression level status
            level = data.transcription_level
            if level >= 80:
                level_status = "🟢 Very High"
            elif level >= 60:
                level_status = "🟢 High"
            elif level >= 40:
                level_status = "🟡 Medium"
            elif level >= 20:
                level_status = "🟠 Low"
            else:
                level_status = "🔴 Very Low"
            
            self.gene_transcription_data.append({
                'gene_id': gene_id,
                'gene_name': gene_id,
                'function': func_short,
                'start_time': start_time,
                'max_level': data.transcription_level / 100.0,
                'level_percent': data.transcription_level,
                'level_status': level_status,
                'rate': rate,
                'current_level': 0.0,
                'status': 'Waiting',
                'phase': phase,
                'phase_num': phase_num,
                'phase_status': phase_status,
                'category': data.category or "",
                'time_range': data.time_range or "",
                'regulation': reg_type,
                'temporal_order': order_num,
                'strand': data.strand,
                'color': data.color,
                'paused': False  # Track if gene is paused
            })
        
        # Update table with individual control buttons
        self.transcription_table.setRowCount(len(self.gene_transcription_data))
        self.transcription_table.setColumnCount(12)
        self.transcription_table.setHorizontalHeaderLabels([
            "Gene ID", "Level", "Status", "Current", "Phase", "Category", 
            "Regulation", "Order", "Control", "Status Icon", "Paused", "Quantum P"
        ])
        
        for i, data in enumerate(self.gene_transcription_data):
            self.transcription_table.setItem(i, 0, QTableWidgetItem(data['gene_id']))
            self.transcription_table.setItem(i, 1, QTableWidgetItem(f"{data['level_percent']:.0f}%"))
            self.transcription_table.setItem(i, 2, QTableWidgetItem(data['level_status']))
            self.transcription_table.setItem(i, 3, QTableWidgetItem("0.00"))
            self.transcription_table.setItem(i, 4, QTableWidgetItem(data['phase_status']))
            self.transcription_table.setItem(i, 5, QTableWidgetItem(data['category']))
            self.transcription_table.setItem(i, 6, QTableWidgetItem(data['regulation']))
            self.transcription_table.setItem(i, 7, QTableWidgetItem(str(data['temporal_order'])))
            
            # Create pause/resume button for this gene
            pause_btn = QPushButton("⏸️ Pause")
            pause_btn.setMinimumWidth(85)
            pause_btn.setMaximumHeight(30)
            pause_btn.setProperty("gene_index", i)  # Store gene index
            pause_btn.clicked.connect(lambda checked, idx=i: self.toggle_gene_pause(idx))
            self.transcription_table.setCellWidget(i, 8, pause_btn)
            self.gene_pause_buttons[i] = pause_btn
            
            self.transcription_table.setItem(i, 9, QTableWidgetItem("⏸️ Waiting"))
            self.transcription_table.setItem(i, 10, QTableWidgetItem("No"))
            self.transcription_table.setItem(i, 11, QTableWidgetItem("-"))
            
            # Color by Phase (as specified in BED file meaning)
            phase_num = data['phase_num']
            if phase_num in [1, 2]:
                color = QColor(255, 150, 150)  # Red - Early (Phase 1-2)
            elif phase_num in [3, 4]:
                color = QColor(255, 255, 150)  # Yellow - Middle (Phase 3-4)
            elif phase_num == 5:
                color = QColor(150, 255, 150)  # Green - Late (Phase 5)
            elif phase_num == 6:
                color = QColor(150, 150, 255)  # Blue - Final (Phase 6)
            else:
                # Fallback: color by level
                level = data['level_percent']
                if level >= 70:
                    color = QColor(200, 255, 200)
                elif level >= 40:
                    color = QColor(255, 255, 200)
                else:
                    color = QColor(255, 230, 230)
            
            for j in range(12):
                if self.transcription_table.item(i, j):
                    self.transcription_table.item(i, j).setBackground(color)
        
        print(f"DEBUG: Prepared {len(self.gene_transcription_data)} genes for transcription animation from BED")
        
        # === CLUSTER INITIALIZATION FOR 28-QUBIT QUANTUM MODEL ===
        if CLUSTER_DEFINITIONS_AVAILABLE:
            self._init_cluster_mapping()

    def _init_cluster_mapping(self):
        """Initialize gene-to-cluster mapping for 28-qubit quantum model."""
        if not CLUSTER_DEFINITIONS_AVAILABLE:
            return
        
        # Build gene-to-cluster mapping
        self.gene_cluster_mapping = build_gene_cluster_mapping(self.gene_transcription_data)
        self.cluster_quantum_p = {i: 0.0 for i in range(28)}  # 28 cluster Quantum P values
        self.cluster_worker = None  # Active cluster worker
        self.cluster_job_pending = False
        self.last_cluster_trigger_time = 0.0
        
        # Get cluster gene counts for display
        self.cluster_gene_counts = get_cluster_gene_counts(self.gene_cluster_mapping)
        
        # Update gene data with cluster IDs
        for data in self.gene_transcription_data:
            gene_id = data['gene_id']
            data['cluster_id'] = self.gene_cluster_mapping.get(gene_id, 27)
        
        print(f"DEBUG: [Cluster] Mapped {len(self.gene_cluster_mapping)} genes to 28 clusters")
        for cluster_id, count in sorted(self.cluster_gene_counts.items()):
            if count > 0:
                print(f"  Cluster {cluster_id} ({get_cluster_name(cluster_id)}): {count} genes")
    
    def _trigger_cluster_quantum_job(self):
        """Trigger a 28-qubit cluster quantum circuit job."""
        if not QUANTUM_AVAILABLE or not CLUSTER_DEFINITIONS_AVAILABLE:
            return
        
        if self.cluster_job_pending:
            return  # Already have a pending job
        
        # Aggregate current transcription levels by cluster
        cluster_levels = aggregate_cluster_levels(
            self.gene_transcription_data,
            self.gene_cluster_mapping,
            level_key='current_level'
        )
        
        # Create and start cluster worker
        self.cluster_worker = QuantumClusterWorker(cluster_levels, ENTANGLEMENT_PAIRS)
        self.cluster_worker.cluster_finished.connect(self._on_cluster_job_finished)
        self.cluster_worker.cluster_failed.connect(self._on_cluster_job_failed)
        self.cluster_job_pending = True
        self.cluster_worker.start()
        
        print(f"DEBUG: [Quantum] Triggered 28-qubit cluster job at t={self.animation_time:.1f}s")
    
    def _on_cluster_job_finished(self, cluster_probs: dict, counts: dict):
        """Handle 28-qubit cluster job completion."""
        self.cluster_job_pending = False
        self.cluster_quantum_p = cluster_probs
        
        print(f"DEBUG: [Quantum] 28-qubit cluster results received:")
        for cluster_id, prob in sorted(cluster_probs.items()):
            if self.cluster_gene_counts.get(cluster_id, 0) > 0:
                print(f"  Cluster {cluster_id} ({get_cluster_name(cluster_id)}): {prob:.1f}%")
        
        # Update each gene's Quantum P from its cluster
        for i, data in enumerate(self.gene_transcription_data):
            cluster_id = data.get('cluster_id', 27)
            gene_quantum_p = cluster_probs.get(cluster_id, 0.0)
            
            # Update gene data
            data['quantum_activation'] = gene_quantum_p
            
            # Update table display
            display_text = f"✅ {gene_quantum_p:.1f}%"
            item = QTableWidgetItem(display_text)
            item.setBackground(QColor(200, 255, 200))
            item.setToolTip(f"Cluster: {get_cluster_name(cluster_id)} (ID: {cluster_id})")
            self.transcription_table.setItem(i, 11, item)
    
    def _on_cluster_job_failed(self, error_msg: str):
        """Handle cluster job failure."""
        self.cluster_job_pending = False
        print(f"ERROR: [Quantum] Cluster job failed: {error_msg}")
        
        # Update all genes to show error
        for i in range(len(self.gene_transcription_data)):
            self.transcription_table.setItem(i, 11, QTableWidgetItem("❌ Cluster Error"))

    
    def generate_comparison_report_from_bed(self):
        """Generate comparison report from BED data"""
        report = "=== Query Genome Analysis (from BED files) ===\n\n"
        
        # Regulation statistics
        regulation_data = self.query_genome_manager.regulation_data
        if regulation_data:
            type_counts = {}
            for gene_id, data in regulation_data.items():
                reg_type = data.regulation_type
                type_counts[reg_type] = type_counts.get(reg_type, 0) + 1
            
            report += f"📊 Regulation Analysis: {len(regulation_data)} genes\n"
            report += "Regulation Types:\n"
            for reg_type, count in sorted(type_counts.items(), key=lambda x: -x[1]):
                report += f"  • {reg_type}: {count}\n"
            report += "\n"
        
        # Transcription statistics
        transcription_data = self.query_genome_manager.transcription_dynamics
        if transcription_data:
            levels = [d.transcription_level for d in transcription_data.values()]
            high = sum(1 for l in levels if l >= 70)
            medium = sum(1 for l in levels if 40 <= l < 70)
            low = sum(1 for l in levels if l < 40)
            
            report += f"🧬 Transcription Dynamics: {len(transcription_data)} genes\n"
            report += f"  • High (≥70%): {high}\n"
            report += f"  • Medium (40-70%): {medium}\n"
            report += f"  • Low (<40%): {low}\n"
            report += f"  • Max Level: {max(levels):.0f}%\n"
            report += f"  • Average: {sum(levels)/len(levels):.1f}%\n\n"
        
        # Translation statistics
        translation_data = self.query_genome_manager.translation_dynamics
        if translation_data:
            levels = [d.protein_level for d in translation_data.values()]
            high = sum(1 for l in levels if l >= 80)
            medium = sum(1 for l in levels if 50 <= l < 80)
            low = sum(1 for l in levels if l < 50)
            
            # RBS counts
            rbs_counts = {}
            for d in translation_data.values():
                rbs = d.rbs_strength or 'unknown'
                rbs_counts[rbs] = rbs_counts.get(rbs, 0) + 1
            
            report += f"🔬 Translation Dynamics: {len(translation_data)} genes\n"
            report += f"  • High Protein (≥80%): {high}\n"
            report += f"  • Medium (50-80%): {medium}\n"
            report += f"  • Low (<50%): {low}\n"
            report += f"  • Max Level: {max(levels):.0f}%\n"
            report += f"  • Average: {sum(levels)/len(levels):.1f}%\n"
            report += f"RBS Strength Distribution:\n"
            for rbs, count in sorted(rbs_counts.items(), key=lambda x: -x[1]):
                report += f"  • {rbs}: {count}\n"
        
        return report
    
    def prepare_transcription_animation_data(self, results, gene_limit=None):
        self.gene_transcription_data = []
        
        if gene_limit is None:
            gene_limit = getattr(self, 'gene_limit_spinbox', None)
            gene_limit = gene_limit.value() if gene_limit else 473
        
        trans_states = results.get('query_transcription_states', {})
        promoter_preds = results.get('query_promoter_predictions', {})
        
        if not trans_states:
            trans_states = results.get('transcription_states', {})
            promoter_preds = results.get('promoter_predictions', {})
            genes_to_use = self.genome_manager.genes
        else:
            genes_to_use = self.query_genome_manager.genes
        
        gene_completion_times = {}
        
        sorted_states = sorted(trans_states.items(), key=lambda x: x[1].temporal_order)
        genes_to_animate_ids = set(gene_id for gene_id, _ in sorted_states[:gene_limit])
        
        for gene_id, state in sorted_states:
            if gene_id not in genes_to_animate_ids:
                continue
            
            if state.state == 'silent':
                start_time = 999999
                duration = 0
            else:
                if hasattr(state, 'dependency_genes') and state.dependency_genes:
                    dep_completion_times = []
                    for dep_id in state.dependency_genes:
                        if dep_id in gene_completion_times:
                            dep_completion_times.append(gene_completion_times[dep_id])
                    
                    if dep_completion_times:
                        start_time = max(dep_completion_times) + 2.0
                    else:
                        start_time = state.temporal_order * 5.0
                else:
                    if state.state == 'constitutive':
                        start_time = state.temporal_order * 3.0
                    elif state.state == 'induced':
                        start_time = state.temporal_order * 5.0
                    else:
                        start_time = state.temporal_order * 8.0
                
                base_duration = 8.0
                if state.state == 'constitutive':
                    duration = base_duration * (1.2 - state.expression_level * 0.4)
                elif state.state == 'induced':
                    duration = base_duration * (1.0 - state.expression_level * 0.3)
                else:
                    duration = base_duration * (1.5 - state.expression_level * 0.5)
            
            if state.state != 'silent':
                gene_completion_times[gene_id] = start_time + duration
            
            rate = 0.2 / duration if duration > 0 and state.state != 'silent' else 0.0
            
            deps_display = "None"
            if hasattr(state, 'dependency_genes') and state.dependency_genes:
                deps_display = ', '.join(state.dependency_genes[:3])
            
            self.gene_transcription_data.append({
                'gene_id': gene_id,
                'gene_name': next((g.name for g in genes_to_use if g.id == gene_id), gene_id),
                'start_time': start_time,
                'max_level': state.expression_level,
                'rate': rate,
                'current_level': 0.0,
                'status': 'Waiting',
                'state': state.state,
                'temporal_order': state.temporal_order if state.temporal_order < 999999 else 'N/A',
                'confidence': state.confidence,
                'dependencies': deps_display,
                'quantum_milestones_triggered': set()
            })
        
        self.transcription_table.setRowCount(len(self.gene_transcription_data))
        self.transcription_table.setColumnCount(12)
        self.transcription_table.setHorizontalHeaderLabels([
            "Gene ID", "Level", "Status", "Current", "Phase", "Category", 
            "Regulation", "Order", "Control", "Status Icon", "Paused", "Quantum P"
        ])
        
        state_icons = {
            'constitutive': '🔵 Constitutive',
            'induced': '🟢 Induced',
            'conditional': '🟡 Conditional',
            'silent': '🔴 Silent'
        }
        
        for i, data in enumerate(self.gene_transcription_data):
            self.transcription_table.setItem(i, 0, QTableWidgetItem(data['gene_id']))
            self.transcription_table.setItem(i, 1, QTableWidgetItem(state_icons.get(data['state'], data['state'])))
            self.transcription_table.setItem(i, 2, QTableWidgetItem("⏸️ Waiting"))
            self.transcription_table.setItem(i, 3, QTableWidgetItem("0.00"))
            self.transcription_table.setItem(i, 4, QTableWidgetItem(str(data['temporal_order'])))
            self.transcription_table.setItem(i, 5, QTableWidgetItem(f"{data['confidence']:.2f}"))
            self.transcription_table.setItem(i, 6, QTableWidgetItem(data['dependencies']))
            self.transcription_table.setItem(i, 11, QTableWidgetItem("-"))
    
    def on_gene_limit_changed(self, value):
        if hasattr(self, 'ai_results') and self.ai_results:
            reply = QMessageBox.question(
                self, "Update Animation",
                f"Reload animation with {value} genes?",
                QMessageBox.Yes | QMessageBox.No
            )
            if reply == QMessageBox.Yes:
                self.prepare_transcription_animation_data(self.ai_results, value)
                self.reset_transcription_animation()
    
    def on_quantum_job_finished(self, gene_id, probability, counts):
        """Callback for successful quantum job completion"""
        # Find the row for this gene
        for i, data in enumerate(self.gene_transcription_data):
            if data['gene_id'] == gene_id:
                # Update the result
                # Update the result
                data['quantum_activation'] = probability
                
                # Check progress (how many milestones triggered)
                milestones_done = len(data.get('quantum_milestones_triggered', set()))
                total_milestones = 4
                
                # Update table column 11 (Quantum P)
                if milestones_done >= total_milestones:
                     display_text = f"✅ {probability:.1f}%"
                     bg_color = QColor(200, 255, 200) # Greenish for complete
                else:
                     display_text = f"{probability:.1f}% ({milestones_done}/{total_milestones})"
                     bg_color = QColor(220, 200, 255) # Light purple for ongoing
                
                item = QTableWidgetItem(display_text)
                item.setBackground(bg_color)
                item.setToolTip(f"Counts: {json.dumps(counts)}")
                self.transcription_table.setItem(i, 11, item)
                
                print(f"DEBUG: [Quantum] Result for {gene_id}: {probability}% (Milestone {milestones_done}/{total_milestones})")
                break
        
        # Cleanup worker
        if gene_id in self.quantum_jobs:
            del self.quantum_jobs[gene_id]

    def on_quantum_job_failed(self, gene_id, error_message):
        """Callback for failed quantum jobs"""
        print(f"ERROR: [Quantum] Job for {gene_id} failed: {error_message}")
        
        # Update table to show error
        for i, data in enumerate(self.gene_transcription_data):
            if data['gene_id'] == gene_id:
                self.transcription_table.setItem(i, 11, QTableWidgetItem("❌ Error"))
                break
                
        # Cleanup worker
        if gene_id in self.quantum_jobs:
            del self.quantum_jobs[gene_id]

    def start_transcription_animation(self):
        if not self.gene_transcription_data:
            QMessageBox.warning(self, "No Data", "Please run AI analysis first")
            return
        
        self.animation_running = True
        self.start_animation_btn.setEnabled(False)
        self.stop_animation_btn.setEnabled(True)
        
        self.animation_timer = QTimer()
        self.animation_timer.timeout.connect(self.update_transcription_animation)
        self.animation_timer.start(100)
        
        # BIOLOGICAL SIMULATION: Auto-start translation when transcription starts
        if self.gene_translation_data and not self.translation_animation_running:
            self.start_translation_animation()
            print("DEBUG: Translation auto-started with transcription (biological coupling)")
    
    def stop_transcription_animation(self):
        self.animation_running = False
        self.start_animation_btn.setEnabled(True)
        self.stop_animation_btn.setEnabled(False)
        if self.animation_timer:
            self.animation_timer.stop()
        
        # BIOLOGICAL SYNC: Stop translation animation too
        self.stop_translation_animation()
        
        # Broadcast 'stopped' status to all layers
        self._broadcast_status('stopped')
    
    def reset_transcription_animation(self):
        self.stop_transcription_animation()
        self.animation_time = 0.0
        self.time_label.setText("⏱️ Time: 0.0s")
        
        # Clear all paused genes and quantum jobs
        self.paused_genes.clear()
        self.clear_all_translation_highlights()
        
        # Stop all pending quantum workers
        for worker in self.quantum_jobs.values():
            worker.stop()
            worker.wait()
        self.quantum_jobs.clear()
        self.quantum_genes_count = 0
        
        for data in self.gene_transcription_data:
            data['current_level'] = 0.0
            data['status'] = 'Waiting'
            data['paused'] = False
            if 'quantum_milestones_triggered' in data:
                data['quantum_milestones_triggered'].clear()
            if 'quantum_activation' in data:
                del data['quantum_activation']
        
        # BIOLOGICAL SYNC: Reset translation animation too
        self.reset_translation_animation()
        
        # Clear Live Data Viewer
        if hasattr(self, 'live_data_viewer') and self.live_data_viewer:
            self.live_data_viewer.clear_data()
            self.live_data_viewer.set_status("⏸️ System Reset")
        
        # Broadcast 'reset' status to all layers
        self._broadcast_status('reset')
        
        for i, data in enumerate(self.gene_transcription_data):
            # Reset all buttons
            if i in self.gene_pause_buttons:
                self.gene_pause_buttons[i].setText("⏸️ Pause")
                self.gene_pause_buttons[i].setStyleSheet("")
            
            # Update columns: 3=Current, 9=Status Icon, 10=Paused, 11=Quantum P
            self.transcription_table.setItem(i, 3, QTableWidgetItem("0.00"))
            self.transcription_table.setItem(i, 9, QTableWidgetItem("⏸️ Waiting"))
            self.transcription_table.setItem(i, 10, QTableWidgetItem("No"))
            self.transcription_table.setItem(i, 11, QTableWidgetItem("-"))
            
            # Reset color based on Phase
            phase_num = data.get('phase_num', 0)
            if phase_num in [1, 2]:
                color = QColor(255, 150, 150)  # Red - Early
            elif phase_num in [3, 4]:
                color = QColor(255, 255, 150)  # Yellow - Middle
            elif phase_num == 5:
                color = QColor(150, 255, 150)  # Green - Late
            elif phase_num == 6:
                color = QColor(150, 150, 255)  # Blue - Final
            else:
                level = data.get('level_percent', 0)
                if level >= 70:
                    color = QColor(200, 255, 200)
                elif level >= 40:
                    color = QColor(255, 255, 200)
                else:
                    color = QColor(255, 230, 230)
            
            for j in range(12):
                if self.transcription_table.item(i, j):
                    self.transcription_table.item(i, j).setBackground(color)
    
    def update_transcription_animation(self):
        if not self.animation_running:
            return
        
        speed = self.speed_spinbox.value()
        dt = 0.1 * speed  # Time step in seconds
        self.animation_time += dt
        self.time_label.setText(f"⏱️ Time: {self.animation_time:.1f}s")
        
        # Check if continuous mode is enabled
        continuous_mode = getattr(self, 'continuous_mode_checkbox', None)
        is_continuous = continuous_mode.isChecked() if continuous_mode else False
        
        # Get half-lives from UI (convert to decay constants)
        mrna_halflife = getattr(self, 'mrna_halflife_spinbox', None)
        mrna_half_life = mrna_halflife.value() if mrna_halflife else 300.0
        mrna_decay_constant = math.log(2) / mrna_half_life if mrna_half_life > 0 else 0
        
        # Track steady-state genes
        steady_state_count = 0
        active_count = 0
        
        # === 28-QUBIT CLUSTER QUANTUM TRIGGER (Every 5 seconds) ===
        # Triggered once per frame, not per-gene (Performance optimization)
        if (QUANTUM_AVAILABLE and CLUSTER_DEFINITIONS_AVAILABLE and 
            getattr(self, 'enable_quantum_checkbox', None) and 
            self.enable_quantum_checkbox.isChecked()):
            
            trigger_interval = 5.0
            last_trigger = getattr(self, 'last_cluster_trigger_time', 0.0)
            
            if self.animation_time - last_trigger >= trigger_interval:
                if not getattr(self, 'cluster_job_pending', False):
                    self._trigger_cluster_quantum_job()
                    self.last_cluster_trigger_time = self.animation_time
                    
                    # Update all genes in table once (O(n), not O(n^2))
                    for j in range(len(self.gene_transcription_data)):
                        self.transcription_table.setItem(j, 11, QTableWidgetItem("🛰️ Cluster..."))

        for i, data in enumerate(self.gene_transcription_data):
            # Skip animation update if this gene is paused (but degradation still happens in continuous mode)
            is_paused = i in self.paused_genes
            
            # Initialize mRNA count if not present
            if 'mrna_count' not in data:
                data['mrna_count'] = 0.0
            
            if self.animation_time >= data['start_time']:
                active_count += 1
                
                if data['status'] == 'Waiting':
                    data['status'] = 'Transcribing'
                
                time_since_start = self.animation_time - data['start_time']
                target_level = data.get('level_percent', data['max_level'] * 100)
                
                if is_continuous and target_level > 0:
                    # === CONTINUOUS MODE: Synthesis + Degradation ===
                    # Transcription rate derived from level (higher level = higher rate)
                    transcription_rate = (target_level / 100.0) * 2.0  # Scale to reasonable rate

                    # Calculate steady-state mRNA count
                    mrna_steady_state = transcription_rate / mrna_decay_constant if mrna_decay_constant > 0 else target_level
                    
                    # Synthesis (only if not paused)
                    if not is_paused:
                        mrna_synthesis = transcription_rate * dt
                    else:
                        mrna_synthesis = 0.0
                    
                    # Degradation (always happens)
                    mrna_decay = mrna_decay_constant * data['mrna_count'] * dt
                    
                    # Update mRNA count
                    data['mrna_count'] = max(0.0, data['mrna_count'] + mrna_synthesis - mrna_decay)
                    
                    # Calculate percentage of steady state for display
                    if mrna_steady_state > 0:
                        data['current_level'] = (data['mrna_count'] / mrna_steady_state) * target_level
                    else:
                        data['current_level'] = 0.0
                    
                    # Status based on how close to steady state
                    pct_of_steady = (data['mrna_count'] / mrna_steady_state * 100) if mrna_steady_state > 0 else 0
                    
                    if is_paused:
                        status_icon = "⏸️ Declining"
                        # Color gradient showing decline
                        decline_pct = min(1.0, data['current_level'] / target_level) if target_level > 0 else 0
                        color = QColor(255, int(180 + decline_pct * 50), int(150 + decline_pct * 50))
                    elif pct_of_steady >= 90:
                        status_icon = "🔄 Steady State"
                        steady_state_count += 1
                        # Pulsing green to show active equilibrium
                        pulse = int(10 * math.sin(self.animation_time * 2))  # Subtle pulse
                        color = QColor(160 + pulse, 255, 160 + pulse)
                    else:
                        status_icon = "🧬 Rising"
                        progress = data['current_level'] / target_level if target_level > 0 else 0
                        # Color based on phase during animation
                        phase_num = data.get('phase_num', 0)
                        if phase_num in [1, 2]:
                            color = QColor(255, int(150 + progress * 50), int(150 + progress * 50))
                        elif phase_num in [3, 4]:
                            color = QColor(255, 255, int(150 + progress * 50))
                        elif phase_num == 5:
                            color = QColor(int(150 + progress * 50), 255, int(150 + progress * 50))
                        elif phase_num == 6:
                            color = QColor(int(150 + progress * 50), int(150 + progress * 50), 255)
                        else:
                            intensity = int(200 + progress * 55)
                            color = QColor(200, intensity, 200)
                else:
                    # === LINEAR MODE: Original behavior ===
                    if is_paused:
                        # In linear mode, pausing just freezes the level
                        pass
                    elif target_level > 0:
                        if NUMPY_AVAILABLE:
                            growth = 1 - np.exp(-data['rate'] * time_since_start)
                        else:
                            growth = min(1.0, data['rate'] * time_since_start)
                        data['current_level'] = min(target_level, target_level * growth)
                    else:
                        data['current_level'] = 0.0
                    
                    if data['current_level'] >= target_level * 0.99:
                        status_icon = "✅ Complete"
                        steady_state_count += 1
                        color = QColor(180, 255, 180)
                    else:
                        status_icon = "🧬 Transcribing"
                        progress = data['current_level'] / target_level if target_level > 0 else 0
                        phase_num = data.get('phase_num', 0)
                        if phase_num in [1, 2]:
                            color = QColor(255, int(150 + progress * 50), int(150 + progress * 50))
                        elif phase_num in [3, 4]:
                            color = QColor(255, 255, int(150 + progress * 50))
                        elif phase_num == 5:
                            color = QColor(int(150 + progress * 50), 255, int(150 + progress * 50))
                        elif phase_num == 6:
                            color = QColor(int(150 + progress * 50), int(150 + progress * 50), 255)
                        else:
                            intensity = int(200 + progress * 55)
                            color = QColor(200, intensity, 200)
                
                # Update columns: 3=Current, 9=Status Icon
                self.transcription_table.setItem(i, 3, QTableWidgetItem(f"{data['current_level']:.1f}%"))
                self.transcription_table.setItem(i, 9, QTableWidgetItem(status_icon))
                
                for j in range(12):
                    if self.transcription_table.item(i, j):
                        self.transcription_table.item(i, j).setBackground(color)
            else:
                self.transcription_table.setItem(i, 9, QTableWidgetItem("⏸️ Waiting"))
                for j in range(12):
                    if self.transcription_table.item(i, j):
                        self.transcription_table.item(i, j).setBackground(QColor(240, 240, 240))
        
        # Update steady-state indicator
        if hasattr(self, 'steady_state_label') and active_count > 0:
            pct_steady = (steady_state_count / active_count) * 100
            if is_continuous:
                if pct_steady >= 90:
                    self.steady_state_label.setText(f"🟢 Steady State ({steady_state_count}/{active_count} genes)")
                    self.steady_state_label.setStyleSheet("font-weight: bold; padding: 5px; color: #2E7D32;")
                elif pct_steady >= 50:
                    self.steady_state_label.setText(f"🟡 Approaching ({steady_state_count}/{active_count} genes)")
                    self.steady_state_label.setStyleSheet("font-weight: bold; padding: 5px; color: #F57C00;")
                else:
                    self.steady_state_label.setText(f"🔴 Rising ({steady_state_count}/{active_count} genes)")
                    self.steady_state_label.setStyleSheet("font-weight: bold; padding: 5px; color: #C62828;")
            else:
                self.steady_state_label.setText(f"📊 Linear Mode ({steady_state_count}/{active_count} complete)")
                self.steady_state_label.setStyleSheet("font-weight: bold; padding: 5px; color: #1565C0;")
        
        # ═══════ METABOLIC DASHBOARD FBA UPDATE ═══════
        # Run FBA every 250 sim-seconds and push results to the dashboard
        self._update_metabolic_dashboard_if_due()

        # ═══════ GENE DYNAMICS PANEL UPDATE ═══════
        # Push per-gene data to the 4-tab gene detail dialog
        gene_dynamics = getattr(self, 'gene_dynamics_panel', None)
        if gene_dynamics is not None and gene_dynamics.isVisible():
            gene_dynamics.update_from_simulation(
                self.gene_transcription_data, self.animation_time
            )

        # ═══════ RESOURCE DASHBOARD UPDATE ═══════
        # Push per-gene data to the resource competition gauges
        res_dash = getattr(self, 'resource_dashboard', None)
        if res_dash is not None:
            res_dash.update_from_simulation(
                self.gene_transcription_data, self.animation_time
            )

    
    # ═══════════════════════════════════════════════════════════════════════
    # METABOLIC DASHBOARD INTEGRATION
    # ═══════════════════════════════════════════════════════════════════════
    def _update_metabolic_dashboard_if_due(self):
        """
        Run FBA at defined intervals and push results to the MetabolicDashboard.
        Lazily initializes the MetabolicSolver on first call.
        """
        # Skip if no dashboard
        dashboard = getattr(self, 'metabolic_dashboard', None)
        if dashboard is None:
            return
        
        # Lazy-init solver and state
        if not hasattr(self, '_metabolic_solver'):
            self._metabolic_solver = None
            self._fba_update_interval = 250.0  # seconds of sim time
            self._last_fba_time = -999.0  # Force first update
            self._fba_call_count = 0
            
            if METABOLIC_SOLVER_AVAILABLE:
                import os
                defaults = [
                    "../Data/metabolic_reconstruction.xlsx",
                    "Data/metabolic_reconstruction.xlsx",
                ]
                for d in defaults:
                    if os.path.exists(d):
                        try:
                            self._metabolic_solver = MetabolicSolver(d)
                            print(f"[MetabolicDashboard] Solver loaded from {d}")
                        except Exception as e:
                            print(f"[MetabolicDashboard] Could not load solver: {e}")
                        break
        
        # Skip if solver failed to load
        if self._metabolic_solver is None:
            return
        
        # Check if it's time to run FBA
        if self.animation_time - self._last_fba_time < self._fba_update_interval:
            return
        
        # Gather protein levels from gene transcription data
        protein_levels = {}
        total_protein = 0.0
        for data in getattr(self, 'gene_transcription_data', []):
            gene_id = data.get('gene_id', '')
            # mrna_count acts as a proxy for protein abundance
            protein_count = data.get('mrna_count', 0.0) * 10.0  # Amplify for FBA
            protein_levels[gene_id] = protein_count
            total_protein += protein_count
        
        # Build FBA constraints from protein levels
        kcat_val = 0.001  # Calibrated kcat
        constraints = {}
        solver = self._metabolic_solver
        
        for rxn_id in solver.reactions:
            genes_for_rxn = solver.gene_rxn_map.get(rxn_id, [])
            vmax_total = 0.0
            for g in genes_for_rxn:
                p = protein_levels.get(g, 0.0)
                vmax_total += p * kcat_val
            
            if not genes_for_rxn:
                continue  # Unconstrained
            
            vmax_total = max(vmax_total, 1e-6)
            idx = solver.reactions.index(rxn_id)
            orig_lb = solver.lb[idx]
            orig_ub = solver.ub[idx]
            
            if orig_lb >= 0:
                constraints[rxn_id] = (orig_lb, min(orig_ub, vmax_total))
            else:
                constraints[rxn_id] = (max(orig_lb, -float(vmax_total)), min(orig_ub, float(vmax_total)))
        
        # Run FBA
        try:
            result = solver.solve(constraints)
            if result["status"] == "success":
                self._fba_call_count += 1
                self._last_fba_time = self.animation_time
                fluxes = result["fluxes"]
                
                # Push to dashboard
                dashboard.update_from_fba(
                    sim_time=self.animation_time,
                    fluxes=fluxes,
                    total_protein=total_protein,
                    fba_calls=self._fba_call_count
                )
        except Exception as e:
            print(f"[MetabolicDashboard] FBA error: {e}")

    def toggle_gene_pause(self, gene_index):
        """Toggle pause state for a specific gene with knockout cascade effects"""
        if gene_index >= len(self.gene_transcription_data):
            return
        
        gene_data = self.gene_transcription_data[gene_index]
        gene_id = gene_data['gene_id']
        
        if gene_index in self.paused_genes:
            # === RESUME THE GENE ===
            self.paused_genes.remove(gene_index)
            gene_data['paused'] = False
            
            # Update button
            if gene_index in self.gene_pause_buttons:
                self.gene_pause_buttons[gene_index].setText("⏸️ Pause")
                self.gene_pause_buttons[gene_index].setStyleSheet("")
            
            # Update paused column
            self.transcription_table.setItem(gene_index, 10, QTableWidgetItem("No"))
            
            # BIOLOGICAL SIMULATION: Resuming transcription resumes translation
            self.resume_translation_for_gene(gene_id)
            
            # Clear translation highlight
            self.clear_translation_highlight_for_gene(gene_id)
            
            # === CLEAR KNOCKOUT CASCADE EFFECTS ===
            if KNOCKOUT_EFFECTS_AVAILABLE and hasattr(self, 'dependency_graph') and self.dependency_graph:
                # Clear affected gene highlights
                self._clear_knockout_affected_highlights(gene_id)
                
                # Close knockout dashboard if this was the gene it was showing
                if hasattr(self, 'knockout_dashboard') and self.knockout_dashboard:
                    if self.knockout_dashboard.knocked_gene == gene_id:
                        self.knockout_dashboard.clear_knockout()
                        self.knockout_dashboard.hide()
                
                # Clear from affected genes map
                if gene_id in self.affected_genes_map:
                    del self.affected_genes_map[gene_id]
        else:
            # === PAUSE THE GENE (KNOCKOUT) ===
            self.paused_genes.add(gene_index)
            gene_data['paused'] = True
            
            # Update button
            if gene_index in self.gene_pause_buttons:
                self.gene_pause_buttons[gene_index].setText("▶️ Resume")
                self.gene_pause_buttons[gene_index].setStyleSheet("background-color: #FFC107; color: black;")
            
            # Update paused column
            self.transcription_table.setItem(gene_index, 10, QTableWidgetItem("⏸️ PAUSED"))
            
            # BIOLOGICAL SIMULATION: Pausing transcription affects translation
            self.pause_translation_for_gene(gene_id)
            
            # Highlight corresponding translation row
            self.highlight_translation_for_gene(gene_id)
            
            # === SHOW ENHANCED KNOCKOUT ANALYSIS DIALOG ===
            if ENHANCED_KNOCKOUT_DIALOG_AVAILABLE:
                try:
                    # Gather cascade effects from dependency graph if available
                    cascade_effects = []
                    if KNOCKOUT_EFFECTS_AVAILABLE and hasattr(self, 'dependency_graph') and self.dependency_graph:
                        try:
                            affected = self.dependency_graph.get_downstream_effects(gene_id)
                            for aff in affected:
                                cascade_effects.append({
                                    'gene_id': aff.get('gene_id', '?'),
                                    'rule_type': aff.get('rule_type', 'unknown'),
                                    'reason': aff.get('reason', ''),
                                    'tx_modifier': aff.get('tx_modifier', 1.0),
                                })
                        except Exception:
                            pass
                    
                    dlg = EnhancedKnockoutDialog(
                        gene_data=gene_data,
                        cascade_effects=cascade_effects,
                        total_genes=len(self.gene_transcription_data),
                        parent=self
                    )
                    dlg.exec_()  # Show but don't block — gene is already paused
                except Exception as e:
                    print(f"[EnhancedKnockout] Dialog error: {e}")
            
            # === CALCULATE AND DISPLAY KNOCKOUT CASCADE EFFECTS ===
            if KNOCKOUT_EFFECTS_AVAILABLE:
                self._apply_knockout_cascade_effects(gene_index, gene_id)
    
    def highlight_translation_for_gene(self, gene_id):
        """Highlight the corresponding row in translation table for paused gene"""
        if not self.gene_translation_data:
            return
        
        # Find the gene in translation table
        for i, trans_data in enumerate(self.gene_translation_data):
            if trans_data['gene_id'] == gene_id:
                # Apply bright highlight color
                highlight_color = QColor(255, 255, 0, 180)  # Bright yellow with transparency
                
                for j in range(10):
                    if self.translation_table.item(i, j):
                        item = self.translation_table.item(i, j)
                        # Store original background if not already stored
                        if not item.data(Qt.UserRole):
                            item.setData(Qt.UserRole, item.background())
                        item.setBackground(highlight_color)
                        
                        # Add bold font to emphasize
                        font = item.font()
                        font.setBold(True)
                        item.setFont(font)
                break
    
    def clear_translation_highlight_for_gene(self, gene_id):
        """Clear highlight from translation table for a specific gene"""
        if not self.gene_translation_data:
            return
        
        # Find the gene in translation table
        for i, trans_data in enumerate(self.gene_translation_data):
            if trans_data['gene_id'] == gene_id:
                # Restore original colors
                for j in range(10):
                    if self.translation_table.item(i, j):
                        item = self.translation_table.item(i, j)
                        # Restore original background if stored
                        orig_bg = item.data(Qt.UserRole)
                        if orig_bg:
                            item.setBackground(orig_bg)
                            item.setData(Qt.UserRole, None)
                        
                        # Remove bold font
                        font = item.font()
                        font.setBold(False)
                        item.setFont(font)
                break
    
    def clear_all_translation_highlights(self):
        """Clear all highlights from translation table"""
        for trans_data in self.gene_translation_data:
            self.clear_translation_highlight_for_gene(trans_data['gene_id'])
    
    def pause_translation_for_gene(self, gene_id):
        """Pause translation animation for a specific gene (biological coupling)"""
        for i, trans_data in enumerate(self.gene_translation_data):
            if trans_data['gene_id'] == gene_id:
                # Mark this translation as paused
                trans_data['tx_paused'] = True
                trans_data['pause_time'] = self.translation_animation_time
                # Update visual indicator
                status_item = self.translation_table.item(i, 8)
                if status_item:
                    status_item.setText("⏸️ TX Paused")
                print(f"DEBUG: Translation paused for {gene_id} (transcription paused)")
                
                # === OPEN PAUSED GENE DASHBOARD ===
                if hasattr(self, 'paused_gene_dashboard') and self.paused_gene_dashboard:
                    self.paused_gene_dashboard.set_gene(gene_id, self.translation_animation_time)
                    self.paused_gene_dashboard.show()  # Make dashboard visible
                    self.paused_gene_dashboard.raise_()  # Bring to front
                break
    
    def resume_translation_for_gene(self, gene_id):
        """Resume translation animation for a specific gene (biological coupling)"""
        for i, trans_data in enumerate(self.gene_translation_data):
            if trans_data['gene_id'] == gene_id:
                # Unmark this translation as paused
                trans_data['tx_paused'] = False
                print(f"DEBUG: Translation resumed for {gene_id} (transcription resumed)")
                break
    
    # =========================================================================
    # KNOCKOUT CASCADE EFFECTS METHODS
    # =========================================================================
    
    def _initialize_dependency_graph(self):
        """Initialize the gene dependency graph for knockout cascade calculations."""
        if not KNOCKOUT_EFFECTS_AVAILABLE:
            return
        
        if hasattr(self, 'dependency_graph') and self.dependency_graph:
            return  # Already initialized
        
        try:
            # Create dependency graph
            self.dependency_graph = GeneDependencyGraph()
            
            # Build from transcription data
            gene_data_list = []
            for data in self.gene_transcription_data:
                gene_data_list.append({
                    'gene_id': data['gene_id'],
                    'gene_name': data.get('gene_name', data['gene_id']),
                    'category': data.get('category', ''),
                    'function': data.get('phase', ''),
                })
            
            self.dependency_graph.build_from_gene_data(gene_data_list)
            print(f"DEBUG: Dependency graph initialized with {len(gene_data_list)} genes")
            
            # Initialize tracking structures
            self.affected_genes_map = {}
            self.gene_row_highlights = {}
            
        except Exception as e:
            print(f"ERROR: Failed to initialize dependency graph: {e}")
            self.dependency_graph = None
    
    def _apply_knockout_cascade_effects(self, gene_index: int, gene_id: str):
        """
        Calculate and apply knockout cascade effects when a gene is paused.
        
        This method:
        1. Initializes the dependency graph if needed
        2. Calculates cascade effects using the graph
        3. Highlights affected genes in transcription table
        4. Shows the knockout effects dashboard
        """
        # Initialize graph if needed
        self._initialize_dependency_graph()
        
        if not self.dependency_graph:
            print("DEBUG: Dependency graph not available, skipping cascade effects")
            return
        
        # Calculate cascade effects
        cascade_effects = self.dependency_graph.calculate_cascade(gene_id)
        
        if cascade_effects.get('total_affected', 0) == 0:
            print(f"DEBUG: No cascade effects for {gene_id}")
            return
        
        print(f"DEBUG: [KNOCKOUT] {gene_id} affects {cascade_effects['total_affected']} genes")
        
        # Store affected genes for later cleanup
        affected_ids = []
        # Handle both new format (5 values) and legacy format (2 values)
        for entry in cascade_effects.get('same_cluster_genes', []):
            gid = entry[0]  # First element is always gene_id
            affected_ids.append(gid)
        # Handle both new format (6 values) and legacy format (4 values)
        for entry in cascade_effects.get('connected_genes', []):
            gid = entry[0]  # First element is always gene_id
            affected_ids.append(gid)
        self.affected_genes_map[gene_id] = affected_ids
        
        # Highlight affected genes in transcription table
        self._highlight_affected_genes(cascade_effects)
        
        # Show knockout effects dashboard
        self._show_knockout_dashboard(gene_id, cascade_effects)
    
    def _highlight_affected_genes(self, cascade_effects: dict):
        """
        Highlight affected genes in the transcription table based on effect strength.
        
        Colors:
        - Red (high effect >= 0.7)  
        - Orange (medium effect 0.4-0.7)
        - Yellow (low effect 0.1-0.4)
        """
        # Build lookup of gene_id -> effect strength
        effect_lookup = {}
        
        # Handle both new format (5 values: gene_id, effect, reason, confidence, rule_type)
        # and legacy format (2 values: gene_id, effect)
        for entry in cascade_effects.get('same_cluster_genes', []):
            gene_id = entry[0]
            effect = entry[1]
            effect_lookup[gene_id] = effect
        
        # Handle both new format (6 values) and legacy format (4 values)
        for entry in cascade_effects.get('connected_genes', []):
            gene_id = entry[0]
            effect = entry[1]
            effect_lookup[gene_id] = effect
        
        # Apply highlights to transcription table
        for row_idx, data in enumerate(self.gene_transcription_data):
            gene_id = data['gene_id']
            
            if gene_id in effect_lookup:
                effect = effect_lookup[gene_id]
                
                # Store original color for restoration
                if gene_id not in self.gene_row_highlights:
                    first_item = self.transcription_table.item(row_idx, 0)
                    if first_item:
                        self.gene_row_highlights[gene_id] = first_item.background()
                
                # Get highlight color based on effect strength
                highlight_color = get_effect_color(effect)
                qt_color = QColor(highlight_color[0], highlight_color[1], highlight_color[2])
                
                # Apply to all columns
                for col_idx in range(self.transcription_table.columnCount()):
                    item = self.transcription_table.item(row_idx, col_idx)
                    if item:
                        item.setBackground(qt_color)
                
                # Update status column to show affected status
                icon = get_effect_icon(effect)
                effect_pct = int(effect * 100)
                self.transcription_table.setItem(row_idx, 9, 
                    QTableWidgetItem(f"{icon} Affected ({effect_pct}%)"))
    
    def _clear_knockout_affected_highlights(self, knocked_gene_id: str):
        """
        Clear highlights from genes that were affected by a specific knockout.
        """
        if knocked_gene_id not in self.affected_genes_map:
            return
        
        affected_gene_ids = self.affected_genes_map.get(knocked_gene_id, [])
        
        for row_idx, data in enumerate(self.gene_transcription_data):
            gene_id = data['gene_id']
            
            if gene_id in affected_gene_ids:
                # Check if still affected by another knockout
                still_affected = False
                for other_knocked, other_affected in self.affected_genes_map.items():
                    if other_knocked != knocked_gene_id and gene_id in other_affected:
                        still_affected = True
                        break
                
                if not still_affected:
                    # Restore original background color
                    if gene_id in self.gene_row_highlights:
                        orig_color = self.gene_row_highlights[gene_id]
                        for col_idx in range(self.transcription_table.columnCount()):
                            item = self.transcription_table.item(row_idx, col_idx)
                            if item:
                                item.setBackground(orig_color)
                        del self.gene_row_highlights[gene_id]
                    
                    # Reset status column
                    self.transcription_table.setItem(row_idx, 9, 
                        QTableWidgetItem("🧬 Recovering"))
        
        print(f"DEBUG: Cleared highlights for {len(affected_gene_ids)} genes affected by {knocked_gene_id}")
    
    def _show_knockout_dashboard(self, gene_id: str, cascade_effects: dict):
        """Show the knockout effects dashboard with cascade information."""
        if not KNOCKOUT_EFFECTS_AVAILABLE:
            return
        
        # Create dashboard if not exists
        if not hasattr(self, 'knockout_dashboard') or not self.knockout_dashboard:
            self.knockout_dashboard = KnockoutEffectsDashboard(self)
        
        # Get gene info for display
        gene_info = {}
        for data in self.gene_transcription_data:
            if data['gene_id'] == gene_id:
                gene_info = {
                    'gene_name': data.get('gene_name', ''),
                    'category': data.get('category', ''),
                    'phase': data.get('phase', ''),
                }
                break
        
        # Show dashboard with cascade effects
        self.knockout_dashboard.set_knockout(gene_id, cascade_effects, gene_info)
        self.knockout_dashboard.show()
        self.knockout_dashboard.raise_()
        
        print(f"DEBUG: Showing knockout dashboard for {gene_id}")

    
    def generate_comparison_report(self, results):
        report = "=== Reference vs Query Genome Comparison ===\n\n"
        
        ref_promoters = results.get('promoter_predictions', {})
        ref_trans = results.get('transcription_states', {})
        
        report += f"Reference Genome (JCVI-syn1.0):\n"
        report += f"  Total genes analyzed: {len(ref_promoters)}\n"
        report += f"  Genes with own promoters: {sum(1 for p in ref_promoters.values() if p.has_promoter)}\n"
        report += f"  Operon internal genes: {sum(1 for p in ref_promoters.values() if p.classification == 'operon_internal')}\n\n"
        
        if ref_trans:
            state_counts = defaultdict(int)
            for state in ref_trans.values():
                state_counts[state.state] += 1
            
            report += "Transcription States:\n"
            for state_name, count in state_counts.items():
                report += f"  {state_name.capitalize()}: {count}\n"
        
        query_promoters = results.get('query_promoter_predictions', {})
        query_trans = results.get('query_transcription_states', {})
        
        if query_promoters:
            report += f"\n\nQuery Genome (JCVI-syn3.0):\n"
            report += f"  Total genes analyzed: {len(query_promoters)}\n"
            report += f"  Genes with own promoters: {sum(1 for p in query_promoters.values() if p.has_promoter)}\n"
            report += f"  Operon internal genes: {sum(1 for p in query_promoters.values() if p.classification == 'operon_internal')}\n\n"
            
            if query_trans:
                state_counts = defaultdict(int)
                for state in query_trans.values():
                    state_counts[state.state] += 1
                
                report += "Transcription States:\n"
                for state_name, count in state_counts.items():
                    report += f"  {state_name.capitalize()}: {count}\n"
        
        return report
    
    def export_ai_results(self):
        if not self.ai_results:
            QMessageBox.warning(self, "No Results", "No AI results to export")
            return
        
        try:
            import openpyxl
        except ImportError:
            QMessageBox.warning(self, "Missing Library", "Please install openpyxl: pip install openpyxl")
            return
        
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save AI Results", "ai_analysis_results.xlsx", "Excel Files (*.xlsx)"
        )
        
        if file_path and PANDAS_AVAILABLE:
            with pd.ExcelWriter(file_path, engine='openpyxl') as writer:
                promoter_data = []
                query_promoters = self.ai_results.get('query_promoter_predictions', {})
                if not query_promoters:
                    query_promoters = self.ai_results.get('promoter_predictions', {})
                
                for gene_id, pred in query_promoters.items():
                    promoter_data.append({
                        'Gene_ID': gene_id,
                        'Has_Promoter': pred.has_promoter,
                        'Confidence': pred.confidence,
                        'Promoter_Strength': pred.promoter_strength,
                        'Minus10_Score': pred.motif_scores.get('minus10', 0),
                        'Minus35_Score': pred.motif_scores.get('minus35', 0),
                        'Classification': pred.classification
                    })
                
                if promoter_data:
                    pd.DataFrame(promoter_data).to_excel(writer, sheet_name='Promoters', index=False)
                
                trans_data = []
                query_trans = self.ai_results.get('query_transcription_states', {})
                if not query_trans:
                    query_trans = self.ai_results.get('transcription_states', {})
                
                for gene_id, state in query_trans.items():
                    trans_data.append({
                        'Gene_ID': gene_id,
                        'State': state.state,
                        'Expression_Level': state.expression_level,
                        'Temporal_Order': state.temporal_order,
                        'Confidence': state.confidence,
                        'Regulatory_Factors': ', '.join(state.regulatory_factors)
                    })
                
                if trans_data:
                    pd.DataFrame(trans_data).to_excel(writer, sheet_name='Transcription', index=False)
                
                # Export Translation data
                transl_data = []
                for data in self.gene_translation_data:
                    transl_data.append({
                        'Gene_ID': data['gene_id'],
                        'Protein_Level': data['protein_level'],
                        'RBS_Strength': data['rbs_strength'],
                        'TX_Level': data.get('tx_level', 0),
                        'TL_Level': data.get('tl_level', 0),
                        'Efficiency': data.get('efficiency', 0),
                        'Temporal_Order': data['temporal_order']
                    })
                
                if transl_data:
                    pd.DataFrame(transl_data).to_excel(writer, sheet_name='Translation', index=False)
            
            QMessageBox.information(self, "Export Complete", f"AI results exported to:\n{file_path}")
    
    # ======== Translation Animation Methods ========
    
    def prepare_translation_animation_data(self):
        """Prepare translation data for animation from BED file"""
        self.gene_translation_data = []
        
        # Get translation data from genome manager
        translation_dynamics = self.query_genome_manager.translation_dynamics
        
        if not translation_dynamics:
            QMessageBox.warning(self, "No Data", 
                "No Translation Dynamics data loaded.\nPlease load a Translation Dynamics BED file.")
            return
        
        # Sort by order (if available) or by protein level
        sorted_data = sorted(
            translation_dynamics.items(),
            key=lambda x: (x[1].order if x[1].order > 0 else 9999, -x[1].protein_level)
        )
        
        for idx, (gene_id, data) in enumerate(sorted_data):
            main_level = data.protein_level if data.protein_level > 0 else data.tl_level
            order_num = data.order if data.order > 0 else idx + 1
            
            # Calculate start time based on order
            start_time = idx * 0.3
            
            # Rate based on RBS strength
            rate_map = {'strong': 0.3, 'medium': 0.2, 'weak': 0.1}
            rate = rate_map.get(data.rbs_strength.lower() if data.rbs_strength else '', 0.15)
            
            # Determine expression status based on BED data
            expression_status = ""
            status_reason = ""
            
            if data.is_ncrna:
                expression_status = "🔴 No Translation (ncRNA)"
                status_reason = "Non-coding RNA"
            elif data.color == (255, 0, 0):
                expression_status = "🟡 Low Translation"
                status_reason = f"TX:{data.tx_level}% TL:{data.tl_level}%"
            elif main_level >= 80:
                expression_status = "🟢 High"
                status_reason = f"Protein: {main_level}%"
            elif main_level >= 50:
                expression_status = "🟡 Medium"
                status_reason = f"Protein: {main_level}%"
            else:
                expression_status = "🟠 Low"
                status_reason = f"Protein: {main_level}%"
            
            self.gene_translation_data.append({
                'gene_id': gene_id,
                'protein_level': main_level,
                'rbs_strength': data.rbs_strength or '',
                'stability': data.stability or '',
                'phase': data.phase or '',
                'category': data.category or '',
                'tx_level': data.tx_level,
                'tl_level': data.tl_level,
                'efficiency': data.efficiency,
                'is_ncrna': data.is_ncrna,
                'expression_status': expression_status,
                'status_reason': status_reason,
                'temporal_order': order_num,
                'start_time': start_time,
                'rate': rate,
                'current_level': 0.0,
                'status': 'waiting',
                'color': data.color,
                'tx_paused': False,  # Track if transcription is paused for this gene
                # CSV export fields
                'strand': data.strand,
                'start': data.start,
                'end': data.end,
                'cds_sequence': data.cds_sequence,
                # mRNA copies and transcription rate will be calculated after this loop
                'mrna_copies': 0,
                'mrna_current': 0.0,
                'transcription_rate': 0.0,
                'time_to_complete': 0.0,
                # Translation fields - will be calculated after this loop
                'gene_length': data.end - data.start,
                'promoter_type': '',
                'protein_length': 0,
                'protein_total': 0,
                'protein_current': 0.0,
                'translation_rate': 0.0,
                # === DATA-DRIVEN KINETICS ===
                'protein_half_life': data.protein_half_life,
                'kinetic_source': data.kinetic_params_source,
            })
            
            # === DERIVE GENE-SPECIFIC NOISE from regulation BED ===
            # Look up regulation data for this gene
            reg_data = self.query_genome_manager.regulation_data.get(gene_id, None)
            if reg_data:
                noise_info = derive_gene_noise_level(
                    gene_id=gene_id,
                    regulation_type=reg_data.regulation_type,
                    regulation_strength=reg_data.strength,
                    has_promoter=True,  # Will be updated if promoter data available
                    is_operon_internal=(reg_data.regulation_type.lower() == 'operonmember'),
                    promoter_type=''  # Will be updated from promoter BED
                )
            else:
                # No regulation data - use default minimal noise (JCVI minimal cell)
                noise_info = {
                    'noise_level': 0.05,
                    'regulation_type': 'unknown',
                    'noise_source': 'default',
                    'explanation': 'no regulation data'
                }
            
            # Add noise info to last appended gene data
            self.gene_translation_data[-1]['gene_noise_level'] = noise_info['noise_level']
            self.gene_translation_data[-1]['noise_source'] = noise_info['noise_source']
            self.gene_translation_data[-1]['regulation_type'] = noise_info['regulation_type']
        
        # Populate table with expression status
        self.translation_table.setRowCount(len(self.gene_translation_data))
        self.translation_table.setColumnCount(12)
        self.translation_table.setHorizontalHeaderLabels([
            "Gene ID", "Protein %", "mRNA (Cur/Tot)", "Proteins (Cur/Tot)", "Phase", "Category", "RBS", 
            "Expression Status", "Reason", "Order", "Status", "Current"
        ])
        
        for i, data in enumerate(self.gene_translation_data):
            self.translation_table.setItem(i, 0, QTableWidgetItem(data['gene_id']))
            self.translation_table.setItem(i, 1, QTableWidgetItem(f"{data['protein_level']:.0f}%"))
            self.translation_table.setItem(i, 2, QTableWidgetItem("0/0"))  # mRNA placeholder
            self.translation_table.setItem(i, 3, QTableWidgetItem("0/0"))  # Proteins placeholder
            self.translation_table.setItem(i, 4, QTableWidgetItem(data['phase']))
            self.translation_table.setItem(i, 5, QTableWidgetItem(data['category']))
            self.translation_table.setItem(i, 6, QTableWidgetItem(data['rbs_strength']))
            self.translation_table.setItem(i, 7, QTableWidgetItem(data['expression_status']))
            self.translation_table.setItem(i, 8, QTableWidgetItem(data['status_reason']))
            self.translation_table.setItem(i, 9, QTableWidgetItem(str(data['temporal_order'])))
            self.translation_table.setItem(i, 10, QTableWidgetItem("⏸️ Waiting"))
            self.translation_table.setItem(i, 11, QTableWidgetItem("0.00"))
            
            # Color based on expression status (interpret BED colors)
            if data['is_ncrna'] or data.get('color') == (255, 69, 0):
                color = QColor(255, 100, 100)  # Red - No translation
            elif data.get('color') == (255, 0, 0):
                color = QColor(255, 200, 100)  # Orange - Low efficiency
            elif "🟢" in data['expression_status']:
                color = QColor(200, 255, 200)  # Green - High
            elif "🟡" in data['expression_status']:
                color = QColor(255, 255, 200)  # Yellow - Medium/Low translation
            elif "🟠" in data['expression_status']:
                color = QColor(255, 230, 200)  # Light orange - Low
            else:
                color = QColor(240, 240, 240)
            
            for j in range(10):
                self.translation_table.item(i, j).setBackground(color)
        
        print(f"DEBUG: Prepared {len(self.gene_translation_data)} genes for translation animation")
        
        # Calculate all biological rates for each gene (after all data is loaded)
        for data in self.gene_translation_data:
            gene_id = data['gene_id']
            gene_length = data['end'] - data['start']
            data['gene_length'] = gene_length
            
            # Calculate mRNA copies
            data['mrna_copies'] = self._calculate_mrna_copies(gene_id)
            
            # Get promoter type
            data['promoter_type'] = self._get_promoter_type(gene_id)
            
            # Calculate transcription rate
            tx_rate, time_complete = self._calculate_transcription_rate(
                gene_id, gene_length, data['mrna_copies']
            )
            data['transcription_rate'] = tx_rate
            data['time_to_complete'] = time_complete
            
            # Calculate translation rate and protein counts
            cds_length = len(data.get('cds_sequence', ''))
            protein_length, protein_total, tl_rate = self._calculate_translation_rate(
                gene_id, cds_length, data['mrna_copies'], data['rbs_strength']
            )
            data['protein_length'] = protein_length
            data['protein_total'] = protein_total
            data['translation_rate'] = tl_rate
        
        print(f"DEBUG: Calculated biological rates for {len(self.gene_translation_data)} genes")
    
    def _calculate_mrna_copies(self, gene_id: str) -> int:
        """Calculate mRNA copy number for a gene
        
        Formula: Base Copies × (Tx Level / 100) × Promoter Multiplier
        
        Base Copies by Regulation Strength:
        - high: 100
        - medium: 50
        - low: 20
        
        Promoter Multiplier:
        - BOTH (strong): 1.5
        - ONE (weak): 1.0
        - No promoter (operon): 0.8
        """
        # Get regulation strength (default to medium)
        base_copies = 50  # Default medium
        if hasattr(self, 'query_genome_manager') and self.query_genome_manager:
            regulation_data = self.query_genome_manager.regulation_data
            if gene_id in regulation_data:
                strength = regulation_data[gene_id].strength.lower()
                if strength == 'high':
                    base_copies = 100
                elif strength == 'medium':
                    base_copies = 50
                elif strength == 'low':
                    base_copies = 20
        
        # Get transcription level (default to 50%)
        tx_level = 50.0
        if hasattr(self, 'query_genome_manager') and self.query_genome_manager:
            tx_dynamics = self.query_genome_manager.transcription_dynamics
            if gene_id in tx_dynamics:
                tx_level = tx_dynamics[gene_id].transcription_level
        
        # Get promoter multiplier (default to 1.0)
        promoter_multiplier = 1.0
        if hasattr(self, 'query_genome_manager') and self.query_genome_manager:
            # Check likely_promoter annotations for this gene
            promoter_annotations = self.query_genome_manager.promoter_annotations.get('likely_promoters', {})
            for key, annotation in promoter_annotations.items():
                if annotation.gene_id == gene_id:
                    strand_info = annotation.strand_info.upper()
                    if strand_info == 'BOTH':
                        promoter_multiplier = 1.5
                    elif strand_info == 'ONE':
                        promoter_multiplier = 1.0
                    break
            else:
                # No promoter found - might be operon internal
                promoter_multiplier = 0.8
        
        # Calculate final mRNA copies
        mrna_copies = int(base_copies * (tx_level / 100.0) * promoter_multiplier)
        return max(1, mrna_copies)  # At least 1 copy
    
    def _calculate_transcription_rate(self, gene_id: str, gene_length: int, mrna_copies: int) -> tuple:
        """Calculate transcription rate for a gene
        
        Formula: Rate = (1 ÷ (Gene Length ÷ RNAP_SPEED)) × Promoter Multiplier
        
        Where:
        - RNAP_SPEED = 50 nucleotides per second (RNA Polymerase speed)
        - Promoter Multiplier affects number of simultaneous polymerases:
          - BOTH (strong): 3× (multiple polymerases)
          - ONE (weak): 1.5×
          - No promoter: 1×
        
        Returns:
            tuple: (transcription_rate, time_to_complete)
                - transcription_rate: copies per second
                - time_to_complete: seconds to reach total mrna_copies
        """
        RNAP_SPEED = 50  # nucleotides per second
        
        # Calculate time to make one copy
        if gene_length <= 0:
            gene_length = 500  # Default gene length
        time_per_copy = gene_length / RNAP_SPEED
        
        # Base rate (copies per second)
        base_rate = 1.0 / time_per_copy
        
        # Get promoter polymerase multiplier (different from mRNA copies multiplier)
        polymerase_multiplier = 1.0
        if hasattr(self, 'query_genome_manager') and self.query_genome_manager:
            promoter_annotations = self.query_genome_manager.promoter_annotations.get('likely_promoters', {})
            for key, annotation in promoter_annotations.items():
                if annotation.gene_id == gene_id:
                    strand_info = annotation.strand_info.upper()
                    if strand_info == 'BOTH':
                        polymerase_multiplier = 3.0  # Multiple polymerases
                    elif strand_info == 'ONE':
                        polymerase_multiplier = 1.5
                    break
        
        # Final transcription rate
        transcription_rate = base_rate * polymerase_multiplier
        
        # Time to reach total copies
        time_to_complete = mrna_copies / transcription_rate if transcription_rate > 0 else 0
        
        return round(transcription_rate, 4), round(time_to_complete, 1)
    
    def _get_promoter_type(self, gene_id: str) -> str:
        """Get promoter type for a gene from multiple data sources
        
        Searches in order:
        1. integrated_genes from smart data engine
        2. transcription_dynamics from genome_manager
        3. regulation_data from genome_manager
        4. promoter_annotations from genome_manager
        
        Returns:
            - Promoter type string (e.g., 'BOTH', 'FWD', 'REV', 'ONE')
            - '-35/-10' or similar if promoter elements found
            - 'Yes' if promoter exists but type unknown
            - 'None' if no promoter found
        """
        # 1. Check integrated_genes from smart data engine
        if hasattr(self, 'integrated_genes') and self.integrated_genes:
            gene = self.integrated_genes.get(gene_id)
            if gene:
                if gene.has_promoter:
                    if gene.promoter_type:
                        ptype = gene.promoter_type
                        if 'BOTH' in ptype.upper():
                            return 'BOTH'
                        elif 'FWD' in ptype.upper() or 'ONE' in ptype.upper():
                            return 'ONE'
                        else:
                            return ptype[:15] if len(ptype) > 15 else ptype
                    return 'Yes'
        
        # 2. Check transcription_dynamics for promoter info
        if hasattr(self, 'query_genome_manager') and self.query_genome_manager:
            trans_data = self.query_genome_manager.transcription_dynamics.get(gene_id)
            if trans_data:
                # Check for promoter elements in transcription data
                if hasattr(trans_data, 'has_minus35') and hasattr(trans_data, 'has_minus10'):
                    if trans_data.has_minus35 and trans_data.has_minus10:
                        return 'BOTH'
                    elif trans_data.has_minus35 or trans_data.has_minus10:
                        return 'ONE'
                # Check regulation attribute
                if hasattr(trans_data, 'regulation') and trans_data.regulation:
                    reg = trans_data.regulation
                    if reg and reg.lower() not in ['unknown', 'none', '']:
                        return reg[:15] if len(reg) > 15 else reg
        
        # 3. Check regulation_data for regulation type (can indicate promoter)
        if hasattr(self, 'query_genome_manager') and self.query_genome_manager:
            reg_data = self.query_genome_manager.regulation_data.get(gene_id)
            if reg_data and hasattr(reg_data, 'regulation_type'):
                reg_type = reg_data.regulation_type
                if reg_type and reg_type.lower() not in ['unknown', 'none', '']:
                    # Constitutive genes have promoters
                    if 'constitutive' in reg_type.lower():
                        return 'Yes'
                    return reg_type[:12]
        
        # 4. Fallback to promoter_annotations from genome_manager
        if hasattr(self, 'query_genome_manager') and self.query_genome_manager:
            promoter_annotations = self.query_genome_manager.promoter_annotations.get('likely_promoters', {})
            for key, annotation in promoter_annotations.items():
                if annotation.gene_id == gene_id:
                    return annotation.strand_info.upper() if hasattr(annotation, 'strand_info') else 'Yes'
        
        return "—"  # Use dash instead of None for clearer display
    
    def _calculate_translation_rate(self, gene_id: str, cds_length: int, mrna_copies: int, rbs_strength: str) -> tuple:
        """Calculate translation rate for a gene
        
        Formula: Rate = (1 ÷ (Protein Length ÷ RIBOSOME_SPEED)) × RBS Multiplier × mRNA Copies
        
        Where:
        - Protein Length = CDS Length ÷ 3 (amino acids)
        - RIBOSOME_SPEED = 15 amino acids per second
        - RBS Multiplier affects ribosome binding:
          - strong: 2×
          - medium: 1×
          - weak: 0.5×
        - Proteins per mRNA = 10 (average)
        
        Returns:
            tuple: (protein_length, protein_total, translation_rate)
                - protein_length: amino acids
                - protein_total: total proteins from all mRNA copies
                - translation_rate: proteins per second
        """
        RIBOSOME_SPEED = 15  # amino acids per second
        PROTEINS_PER_MRNA = 10  # average protein copies per mRNA
        
        # Calculate protein length
        if cds_length <= 0:
            cds_length = 300  # Default ~100 amino acids
        protein_length = cds_length // 3
        
        # Calculate total proteins
        protein_total = mrna_copies * PROTEINS_PER_MRNA
        
        # Calculate time per protein
        if protein_length <= 0:
            protein_length = 100
        time_per_protein = protein_length / RIBOSOME_SPEED
        
        # Base rate (one protein per time_per_protein)
        base_rate = 1.0 / time_per_protein
        
        # RBS multiplier
        rbs_multiplier = 1.0
        rbs = rbs_strength.lower() if rbs_strength else ''
        if rbs == 'strong':
            rbs_multiplier = 2.0
        elif rbs == 'medium':
            rbs_multiplier = 1.0
        elif rbs == 'weak':
            rbs_multiplier = 0.5
        
        # Translation rate = base rate × rbs multiplier × mRNA copies
        translation_rate = base_rate * rbs_multiplier * mrna_copies
        
        return protein_length, protein_total, round(translation_rate, 2)
    
    def _broadcast_status(self, status: str):
        """Broadcast a status message to Live Data Viewer and ZMQ"""
        status_data = {
            'type': 'status',
            'status': status,
            'timestamp': getattr(self, 'translation_animation_time', 0.0)
        }
        
        # Update Live Data Viewer
        if hasattr(self, 'live_data_viewer') and self.live_data_viewer:
            status_icons = {
                'started': '🔬 Translation running...',
                'completed': '✅ Translation completed - Data stable',
                'updated': '🔄 Data updated',
                'reset': '🔄 Reset - Data cleared'
            }
            self.live_data_viewer.set_status(status_icons.get(status, status))
        
        # Broadcast via ZMQ
        if hasattr(self, 'zmq_socket') and self.zmq_socket:
            try:
                self.zmq_socket.send_string(json.dumps(status_data))
            except Exception as e:
                print(f"Warning: ZMQ status broadcast error: {e}")
        
        print(f"DEBUG: Broadcast status: {status}")
    
    def _broadcast_data(self, data, force=False):
        """Broadcast gene data to Live Data Viewer and ZMQ socket
        
        Args:
            data: Gene data dictionary
            force: If True, broadcast even if data hasn't changed
            
        The broadcast includes protein availability state so downstream modules
        (e.g., protein structure predictors) can respond to gene deactivation.
        """
        # Track previous levels to detect changes
        if not hasattr(self, '_last_broadcast_levels'):
            self._last_broadcast_levels = {}
        
        gene_id = data['gene_id']
        current_level = round(data['current_level'], 1)  # Round to 1 decimal for comparison
        
        # Check if gene is paused - paused genes ALWAYS need dashboard updates
        is_tx_paused = data.get('tx_paused', False)
        
        # === PAUSED GENE DASHBOARD - Always update BEFORE any early return ===
        # This ensures real-time decay visibility even when levels change slowly
        if is_tx_paused and hasattr(self, 'paused_gene_dashboard') and self.paused_gene_dashboard:
            # === Calculate REAL values from simulation state ===
            
            # mRNA count: Get from transcription data for this gene
            mrna_count = 0.0
            mrna_total_copies = data.get('mrna_copies', 50)  # Default 50 copies
            mrna_current_pct = data.get('mrna_current', 0)
            
            # If mrna_current is stored as actual count
            if isinstance(mrna_current_pct, (int, float)) and mrna_current_pct > 0:
                mrna_count = float(mrna_current_pct)
            else:
                # Try to get from transcription table
                for tx_data in getattr(self, 'gene_transcription_data', []):
                    if tx_data.get('gene_id') == gene_id:
                        tx_level = tx_data.get('current_level', 0)
                        tx_max = tx_data.get('level_percent', 100)
                        if tx_max > 0:
                            mrna_fraction = tx_level / tx_max
                            mrna_count = mrna_fraction * mrna_total_copies
                        break
            
            # Protein count
            protein_count = data.get('protein_count', data.get('protein_current', 0))
            if isinstance(protein_count, (int,)):
                protein_count = float(protein_count)
            protein_total = data.get('protein_total', 100)
            
            # Get half-lives (use defaults if not set)
            mrna_halflife_spinbox = getattr(self, 'mrna_halflife_spinbox', None)
            mrna_half_life = mrna_halflife_spinbox.value() if mrna_halflife_spinbox else 300.0
            
            protein_halflife_spinbox = getattr(self, 'protein_halflife_spinbox', None)
            default_protein_half_life = protein_halflife_spinbox.value() if protein_halflife_spinbox else 3600.0
            
            # Use per-gene half-life from BED if available
            gene_protein_half_life = data.get('protein_half_life', 0.0)
            if gene_protein_half_life > 0:
                protein_half_life = gene_protein_half_life
            else:
                protein_half_life = default_protein_half_life
            
            # Calculate degradation rates from half-lives
            import math
            mrna_decay_rate = math.log(2) / mrna_half_life if mrna_half_life > 0 else 0.0
            protein_decay_rate = math.log(2) / protein_half_life if protein_half_life > 0 else 0.0
            
            # Determine protein state
            if protein_count <= 0:
                protein_state = 'depleted'
            elif mrna_count <= 0:
                protein_state = 'declining_fast'
            else:
                protein_state = 'declining'
            
            dashboard_data = {
                'gene_id': gene_id,
                'timestamp': self.translation_animation_time,
                'mrna_count': mrna_count,
                'protein_count': protein_count,
                'protein_fraction': protein_count / protein_total if protein_total > 0 else 0,
                'protein_state': protein_state,
                'should_predict_structure': data.get('should_predict_structure', protein_count > 10),
                'sequence_valid': data.get('sequence_valid', True),
                'size_category': data.get('size_category', 'medium'),
                'synthesis_rate': 0.0,  # Always 0 when paused
                'degradation_rate': protein_decay_rate,
                'mrna_half_life': mrna_half_life,
                'protein_half_life': protein_half_life,
            }
            self.paused_gene_dashboard.update_data(dashboard_data)
        
        # Skip full broadcast if no significant change (unless forced or paused)
        if not force and not is_tx_paused:
            last_level = self._last_broadcast_levels.get(gene_id, -1)
            # Only broadcast if level changed by at least 0.2% (smoother for synthesis tab)
            if abs(current_level - last_level) < 0.2:
                return
        
        # Update tracking
        self._last_broadcast_levels[gene_id] = current_level
        
        # === PROTEIN AVAILABILITY STATE ===
        # Determine if gene is paused and protein availability
        mrna_count = data.get('mrna_count', data.get('mrna_current', 0))
        protein_count = data.get('protein_count', data.get('protein_current', 0))
        
        # Protein is available if there are translated proteins
        protein_available = protein_count > 0
        
        # Determine protein state for downstream processing
        if is_tx_paused:
            if protein_count <= 0:
                protein_state = 'depleted'  # No protein left
            elif mrna_count <= 0:
                protein_state = 'declining_fast'  # mRNA gone, protein decaying
            else:
                protein_state = 'declining'  # Both still present but decaying
        else:
            if protein_count <= 0:
                protein_state = 'initializing'  # Not yet translated
            elif current_level >= data.get('protein_level', 100) * 0.9:
                protein_state = 'active'  # At or near steady state
            else:
                protein_state = 'rising'  # Still accumulating
        
        # Clear CDS sequence when protein is depleted (downstream should stop processing)
        cds_sequence = data.get('cds_sequence', '')
        if not protein_available:
            cds_sequence = ''  # Signal to downstream: no protein to process
        
        # === PROTEIN STRUCTURE COUPLING (Issue 4) ===
        # Calculate metrics for downstream structure prediction modules
        
        # Protein fraction of expected steady-state (0.0 to 1.0+)
        protein_total = data.get('protein_total', 1)
        protein_fraction = protein_count / protein_total if protein_total > 0 else 0.0
        
        # Structure confidence: higher when protein pool is stable and abundant
        # Low when depleted, declining, or just starting
        if protein_state == 'active':
            structure_confidence = min(1.0, protein_fraction)  # Max when at steady-state
        elif protein_state == 'rising':
            structure_confidence = min(0.8, protein_fraction * 0.8)  # Partial confidence
        elif protein_state in ['declining', 'declining_fast']:
            structure_confidence = protein_fraction * 0.5  # Reduce for declining
        else:  # depleted or initializing
            structure_confidence = 0.0
        
        # === DYNAMIC THRESHOLD (Issue 4) ===
        # Use expression-derived threshold instead of fixed 0.1
        thresholds = getattr(self, 'structure_thresholds', None)
        if thresholds and thresholds.get('threshold_source') == 'expression_data':
            min_fraction = thresholds['min_protein_fraction']
            threshold_source = 'expression_data'
        else:
            min_fraction = 0.1  # Fallback default
            threshold_source = 'default'
        
        # Explicit flag: should downstream process this for structure prediction?
        # True only when protein exceeds dynamically-derived threshold
        should_predict_structure = protein_available and protein_fraction >= min_fraction
        
        # Get half-lives from UI for rate information
        mrna_halflife_spinbox = getattr(self, 'mrna_halflife_spinbox', None)
        protein_halflife_spinbox = getattr(self, 'protein_halflife_spinbox', None)
        mrna_half_life = mrna_halflife_spinbox.value() if mrna_halflife_spinbox else 300.0
        protein_half_life = protein_halflife_spinbox.value() if protein_halflife_spinbox else 3600.0
        
        # === SEQUENCE VALIDATION (Enhanced for ultra-short rejection) ===
        # Validate CDS sequence for biological validity before sending to structure modules
        # Pass source as 'bed' since CDS comes from BED/GFF coordinate extraction
        sequence_validation = validate_protein_sequence(cds_sequence, sequence_source='bed')
        
        # Override should_predict_structure based on validation
        if not sequence_validation['is_valid']:
            should_predict_structure = False
        
        # Adjust prediction mode for short proteins
        prediction_mode = sequence_validation.get('prediction_mode', 'standard')
        if prediction_mode == 'none':
            should_predict_structure = False
        elif prediction_mode == 'limited':
            # Short proteins: only monomer prediction, no complex assembly
            pass  # Structure prediction allowed but limited
        
        # === STRICT FILTERING: Block biologically impossible sequences ===
        # Clear sequence for invalid/marginal categories to prevent "YC" type outputs
        size_category = sequence_validation.get('size_category', 'unknown')
        if size_category in ['invalid', 'marginal'] or sequence_validation['validity_score'] < 0.3:
            cds_sequence = ''  # Absolutely don't send to structure modules
            should_predict_structure = False
        
        # === QUATERNARY STRUCTURE (Issue 6) ===
        # Classify protein as monomer vs subunit of multimeric complex
        product = data.get('product', data.get('gene_name', ''))
        quaternary_info = classify_protein_quaternary(gene_id, product)
        
        # Prepare data dict with timestamp and protein state
        broadcast_data = {
            'type': 'data',
            'timestamp': self.translation_animation_time,
            'gene_id': gene_id,
            'status': data['status'],
            'protein_level': current_level,
            # === Protein availability for downstream modules ===
            'protein_available': protein_available,
            'protein_state': protein_state,
            'mrna_count': mrna_count,
            'protein_count': protein_count,
            'is_paused': is_tx_paused,
            # === STRUCTURE COUPLING (Issue 4) ===
            'protein_fraction': round(protein_fraction, 3),
            'structure_confidence': round(structure_confidence, 3),
            'should_predict_structure': should_predict_structure,
            'synthesis_rate': data.get('translation_rate', 0.0),
            'degradation_rate': round(0.693 / protein_half_life, 6) if protein_half_life > 0 else 0,
            'mrna_half_life': mrna_half_life,
            'protein_half_life': protein_half_life,
            # === ENHANCED SEQUENCE VALIDATION ===
            'sequence_valid': sequence_validation['is_valid'],
            'sequence_issues': sequence_validation['issues'],
            'sequence_warnings': sequence_validation.get('warnings', []),
            'validity_score': sequence_validation['validity_score'],
            'has_start_codon': sequence_validation['has_start_codon'],
            'has_stop_codon': sequence_validation['has_stop_codon'],
            'internal_stop_count': sequence_validation.get('internal_stop_count', 0),
            # === PROTEIN SIZE CLASSIFICATION ===
            'size_category': size_category,  # invalid/marginal/short/small/medium/large
            'prediction_mode': prediction_mode,  # none/flagged/limited/standard
            'is_subunit_candidate': sequence_validation.get('is_subunit_candidate', False),
            'needs_review': sequence_validation.get('needs_review', False),
            'sequence_source': sequence_validation.get('sequence_source', 'unknown'),
            # === QUATERNARY STRUCTURE (Issue 6) ===
            'quaternary_type': quaternary_info['quaternary_type'],
            'complex_name': quaternary_info['complex_name'],
            'complex_id': quaternary_info['complex_id'],
            'subunit_role': quaternary_info['subunit_role'],
            'requires_partners': quaternary_info['requires_partners'],
            'stoichiometry': quaternary_info['stoichiometry'],
            'is_assembled': quaternary_info['is_assembled'],
            'assembly_partners': quaternary_info['assembly_partners'],
            # Transcription data
            'gene_length': data.get('gene_length', 0),
            'promoter_type': data.get('promoter_type', 'None'),
            'mrna_copies': data.get('mrna_copies', 0),
            'mrna_current': data.get('mrna_current', 0.0),
            'transcription_rate': data.get('transcription_rate', 0.0),
            'time_to_complete': data.get('time_to_complete', 0.0),
            # Translation data
            'protein_length': data.get('protein_length', 0),
            'protein_total': data.get('protein_total', 0),
            'protein_current': data.get('protein_current', 0.0),
            'translation_rate': data.get('translation_rate', 0.0),
            # Gene location
            'strand': data.get('strand', '+'),
            'start': data.get('start', 0),
            'end': data.get('end', 0),
            'cds_sequence': cds_sequence,  # Empty when invalid or protein depleted
            # === JCVI SPECIFICITY (Issue 5) ===
            'is_jcvi_genome': getattr(self, 'jcvi_validation', {}).get('is_jcvi', False),
            'jcvi_confidence': getattr(self, 'jcvi_validation', {}).get('confidence', 0),
            'kinetic_source': data.get('kinetic_source', 'default'),
            'noise_source': data.get('noise_source', 'default'),
            'structure_threshold_source': threshold_source
        }
        
        # Add row to Live Data Viewer dialog
        if hasattr(self, 'live_data_viewer') and self.live_data_viewer:
            self.live_data_viewer.add_row(broadcast_data)
        
        # === INTEGRATED AMINO ACID GENERATOR ===
        # Update the integrated Protein Synthesis widget
        if hasattr(self, 'protein_synthesis_widget') and self.protein_synthesis_widget:
            self.protein_synthesis_widget.add_row(broadcast_data)
        
        # Note: Paused gene dashboard is updated earlier in this function,
        # BEFORE the early return check, to ensure real-time decay visibility
        
        # Broadcast via ZMQ if available
        if hasattr(self, 'zmq_socket') and self.zmq_socket:
            try:
                self.zmq_socket.send_string(json.dumps(broadcast_data))
            except Exception as e:
                print(f"Warning: ZMQ broadcast error: {e}")

    def start_translation_animation(self):
        """Start translation animation"""
        if not self.gene_translation_data:
            self.prepare_translation_animation_data()
        
        if not self.gene_translation_data:
            return
        
        # Populate CDS sequences if not already done
        if self.query_genome_manager:
            self.query_genome_manager.populate_translation_cds_sequences()
            
            # === GENOME-DERIVED RESOURCES (Issue 3) ===
            # Derive cellular resources from actual gene content
            if CONTINUOUS_ENGINE_AVAILABLE and hasattr(self.query_genome_manager, 'genes'):
                self.cellular_resources = derive_resources_from_genome(
                    self.query_genome_manager.genes,
                    self.query_genome_manager.translation_dynamics
                )
            else:
                self.cellular_resources = None
            
            # === JCVI ORGANISM VALIDATION (Issue 5) ===
            # Validate that loaded data is from JCVI syn3.0 genome
            if CONTINUOUS_ENGINE_AVAILABLE:
                self.jcvi_validation = validate_jcvi_organism(
                    genome_sequence=self.query_genome_manager.sequence,
                    genes=self.query_genome_manager.genes
                )
                if self.jcvi_validation['warnings']:
                    print("⚠️ JCVI Specificity Warnings:")
                    for warning in self.jcvi_validation['warnings']:
                        print(f"   - {warning}")
                if self.jcvi_validation['is_jcvi']:
                    print(f"✅ JCVI syn3.0 genome confirmed (confidence: {self.jcvi_validation['confidence']:.0%})")
            else:
                self.jcvi_validation = {'is_jcvi': False, 'warnings': ['Engine not available']}
            
            # === DYNAMIC STRUCTURE THRESHOLDS (Issue 4) ===
            # Compute thresholds from protein expression distribution
            self.structure_thresholds = compute_structure_prediction_thresholds(
                self.query_genome_manager.translation_dynamics
            )
            
            # Update gene_translation_data with CDS sequences and recalculate protein values
            for data in self.gene_translation_data:
                gene_id = data['gene_id']
                trans_data = self.query_genome_manager.translation_dynamics.get(gene_id)
                if trans_data:
                    data['cds_sequence'] = trans_data.cds_sequence
                    
                    # Recalculate protein length and translation rate with actual CDS
                    cds_length = len(trans_data.cds_sequence)
                    protein_length, protein_total, tl_rate = self._calculate_translation_rate(
                        gene_id, cds_length, data['mrna_copies'], data['rbs_strength']
                    )
                    data['protein_length'] = protein_length
                    data['protein_total'] = protein_total
                    data['translation_rate'] = tl_rate
        
        # Clear and show Live Data Viewer
        if hasattr(self, 'live_data_viewer') and self.live_data_viewer:
            self.live_data_viewer.clear_data()
            self.live_data_viewer.show()
        
        # === INITIALIZE PAUSED GENE DASHBOARD ===
        if DASHBOARD_AVAILABLE:
            if not hasattr(self, 'paused_gene_dashboard') or not self.paused_gene_dashboard:
                self.paused_gene_dashboard = PausedGeneDashboard(self)
        else:
            self.paused_gene_dashboard = None
        
        # Reset change tracking
        self._last_broadcast_levels = {}
        self._all_complete_sent = False
        
        # Broadcast 'started' status
        self._broadcast_status('started')
        
        self.translation_animation_running = True
        self.start_translation_btn.setEnabled(False)
        self.stop_translation_btn.setEnabled(True)
        
        self.translation_animation_timer = QTimer()
        self.translation_animation_timer.timeout.connect(self.update_translation_animation)
        self.translation_animation_timer.start(100)
    
    def stop_translation_animation(self):
        """Stop translation animation"""
        self.translation_animation_running = False
        self.start_translation_btn.setEnabled(True)
        self.stop_translation_btn.setEnabled(False)
        if self.translation_animation_timer:
            self.translation_animation_timer.stop()
        
        # Broadcast 'completed' status
        self._broadcast_status('completed')
    
    def reset_translation_animation(self):
        """Reset translation animation"""
        self.stop_translation_animation()
        self.translation_animation_time = 0.0
        self.transl_time_label.setText("⏱️ Time: 0.0 seconds")
        
        for data in self.gene_translation_data:
            data['current_level'] = 0.0
            data['status'] = 'waiting'
        
        for i, data in enumerate(self.gene_translation_data):
            self.translation_table.setItem(i, 6, QTableWidgetItem("⏸️ Waiting"))
            self.translation_table.setItem(i, 7, QTableWidgetItem("0.00"))
            
            # Reset color
            level = data['protein_level']
            if level >= 80:
                color = QColor(200, 255, 200)
            elif level >= 50:
                color = QColor(255, 255, 200)
            else:
                color = QColor(255, 230, 230)
            
            for j in range(8):
                self.translation_table.item(i, j).setBackground(color)
        
        # Reset tracking and broadcast reset status
        self._last_broadcast_levels = {}
        self._all_complete_sent = False
        self._broadcast_status('reset')
        
        # Clear Live Data Viewer
        if hasattr(self, 'live_data_viewer') and self.live_data_viewer:
            self.live_data_viewer.clear_data()
    
    def update_translation_animation(self):
        """Update translation animation frame - supports continuous mode"""
        if not self.translation_animation_running:
            return
        
        speed = self.transl_speed_spinbox.value()
        dt = 0.1 * speed  # Time step in seconds
        self.translation_animation_time += dt
        self.transl_time_label.setText(f"⏱️ Time: {self.translation_animation_time:.1f}s")
        
        # Check if continuous mode is enabled
        continuous_mode = getattr(self, 'continuous_mode_checkbox', None)
        is_continuous = continuous_mode.isChecked() if continuous_mode else False
        
        # Get default half-lives from UI (fallback values)
        protein_halflife_spinbox = getattr(self, 'protein_halflife_spinbox', None)
        default_protein_half_life = protein_halflife_spinbox.value() if protein_halflife_spinbox else 3600.0
        
        all_complete = True  # Track for linear mode only
        
        for i, data in enumerate(self.gene_translation_data):
            # BIOLOGICAL SIMULATION: Check if transcription is paused for this gene
            is_tx_paused = data.get('tx_paused', False)
            
            # === DATA-DRIVEN: Use per-gene half-life from BED if available ===
            # Check if this gene has BED-derived half-life
            gene_protein_half_life = data.get('protein_half_life', 0.0)
            if gene_protein_half_life > 0:
                protein_half_life = gene_protein_half_life
                data['kinetic_source'] = 'bed_file'
            else:
                # Use UI default if no BED data, but ensure it's at least 300s (5m)
                protein_half_life = max(300.0, default_protein_half_life)
                data['kinetic_source'] = 'ui_default'
            
            # Safety floor to prevent division by zero or UI freezes
            protein_half_life = max(1.0, protein_half_life)
            data['protein_half_life'] = protein_half_life  # Update record for consistent indexing
            
            protein_decay_constant = math.log(2) / protein_half_life if protein_half_life > 0 else 0
            
            # === GENE-SPECIFIC NOISE from regulation BED ===
            # Use per-gene noise level instead of global noise
            gene_noise_level = data.get('gene_noise_level', 0.05)
            
            # Initialize protein count if not present
            if 'protein_count' not in data:
                data['protein_count'] = 0.0
                
            if self.translation_animation_time >= data['start_time']:
                if data['status'] == 'waiting':
                    data['status'] = 'translating'
                
                time_since_start = self.translation_animation_time - data['start_time']
                target_level = data['protein_level']
                
                if is_continuous and target_level > 0:
                    # === CONTINUOUS MODE: Protein synthesis + degradation ===
                    
                    # Get current mRNA from transcription (affects protein synthesis rate)
                    # Look up corresponding transcription data
                    gene_id = data.get('gene_id', '')
                    mrna_fraction = 1.0
                    for tx_data in self.gene_transcription_data:
                        if tx_data.get('gene_id') == gene_id:
                            tx_target = tx_data.get('level_percent', 100)
                            tx_current = tx_data.get('current_level', 0)
                            mrna_fraction = tx_current / tx_target if tx_target > 0 else 0
                            break
                    
                    # === BIOLOGICAL FIX: When transcription is paused, mRNA decays ===
                    # This means protein synthesis should STOP, not just reduce
                    if is_tx_paused:
                        # Calculate time since pause started
                        pause_time = data.get('pause_time', self.translation_animation_time)
                        time_since_pause = self.translation_animation_time - pause_time
                        
                        # mRNA decays with half-life (default 300s = 5 minutes for bacteria)
                        mrna_halflife_spinbox = getattr(self, 'mrna_halflife_spinbox', None)
                        mrna_half_life = mrna_halflife_spinbox.value() if mrna_halflife_spinbox else 300.0
                        mrna_decay_constant = math.log(2) / mrna_half_life if mrna_half_life > 0 else 0.002
                        
                        # mRNA level decays exponentially: N(t) = N0 * e^(-kt)
                        mrna_fraction = mrna_fraction * math.exp(-mrna_decay_constant * time_since_pause)
                        
                        # Store for dashboard display
                        data['mrna_decay_fraction'] = mrna_fraction
                    
                    # Translation rate (scale based on target protein level)
                    translation_rate = (target_level / 100.0) * 1.0  # Scale to reasonable rate
                    
                    # Protein steady-state (at full mRNA)
                    protein_steady_state = translation_rate / protein_decay_constant if protein_decay_constant > 0 else target_level
                    
                    # Protein synthesis - ONLY if mRNA is present (proportional to mRNA level)
                    # When paused and mRNA depletes, synthesis approaches zero
                    protein_synthesis = translation_rate * mrna_fraction * dt
                    
                    # Protein degradation (always happens)
                    protein_decay = protein_decay_constant * data['protein_count'] * dt
                    
                    # Update protein count
                    data['protein_count'] = max(0.0, data['protein_count'] + protein_synthesis - protein_decay)
                    
                    # Calculate level percentage for display
                    if protein_steady_state > 0:
                        data['current_level'] = (data['protein_count'] / protein_steady_state) * target_level
                    else:
                        data['current_level'] = 0.0
                    
                    # Status based on state
                    pct_of_steady = (data['protein_count'] / protein_steady_state * 100) if protein_steady_state > 0 else 0
                    
                    if is_tx_paused:
                        status_icon = "📉 Declining"
                        data['status'] = 'declining'
                        decline_pct = min(1.0, data['current_level'] / target_level) if target_level > 0 else 0
                        color = QColor(255, int(200 - decline_pct * 50), int(150 + decline_pct * 50))
                    elif pct_of_steady >= 90:
                        status_icon = "🔄 Steady"
                        data['status'] = 'steady_state'
                        pulse = int(10 * math.sin(self.translation_animation_time * 2))
                        color = QColor(160 + pulse, 255, 160 + pulse)
                    else:
                        status_icon = "🔬 Rising"
                        data['status'] = 'rising'
                        progress = data['current_level'] / target_level if target_level > 0 else 0
                        intensity = int(200 + progress * 55)
                        color = QColor(255, intensity, 180)
                    
                    # Update displays
                    mrna_cur = int(mrna_fraction * data.get('mrna_copies', 10))
                    mrna_tot = data.get('mrna_copies', 10)
                    data['mrna_current'] = mrna_cur
                    data['protein_current'] = int(data['protein_count'])
                    
                    # === SYNC: Push protein data to transcription list for Gene Dynamics Panel ===
                    for tx_data in self.gene_transcription_data:
                        if tx_data.get('gene_id') == gene_id:
                            tx_data['protein_count'] = data['protein_count']
                            tx_data['protein_level'] = data['current_level']
                            tx_data['translation_rate'] = translation_rate
                            tx_data['protein_half_life'] = protein_half_life
                            tx_data['protein_sequence'] = data.get('cds_sequence', 'M') # Simplified
                            break
                    
                else:
                    # === LINEAR MODE: Original behavior ===
                    if is_tx_paused:
                        all_complete = False
                        continue
                    
                    tx_rate = data.get('transcription_rate', 0)
                    mrna_total = data.get('mrna_copies', 0)
                    mrna_current = min(mrna_total, tx_rate * time_since_start)
                    data['mrna_current'] = round(mrna_current, 1)
                    
                    tl_rate = data.get('translation_rate', 0)
                    protein_total = data.get('protein_total', 0)
                    mrna_fraction = mrna_current / mrna_total if mrna_total > 0 else 0
                    effective_tl_rate = tl_rate * mrna_fraction
                    protein_current = min(protein_total, effective_tl_rate * time_since_start)
                    data['protein_current'] = round(protein_current, 0)
                    
                    if target_level > 0:
                        if NUMPY_AVAILABLE:
                            growth = 1 - np.exp(-data['rate'] * time_since_start)
                        else:
                            growth = min(1.0, data['rate'] * time_since_start)
                        data['current_level'] = min(target_level, target_level * growth)
                    else:
                        data['current_level'] = 0.0
                    
                    is_complete = data['current_level'] >= target_level * 0.99
                    if is_complete:
                        data['status'] = 'complete'
                        data['mrna_current'] = data.get('mrna_copies', 0)
                        data['protein_current'] = data.get('protein_total', 0)
                        data['protein_count'] = data.get('protein_total', 0)  # For broadcast
                        status_icon = "✅ Complete"
                        color = QColor(180, 255, 180)
                    else:
                        all_complete = False
                        status_icon = "🔬 Translating"
                        progress = data['current_level'] / target_level if target_level > 0 else 0
                        intensity = int(200 + progress * 55)
                        color = QColor(255, intensity, 180)
                        # Set protein_count for broadcast based on progress
                        data['protein_count'] = data.get('protein_current', 0)
                
                # === SYNC: Push protein data to transcription list for Gene Dynamics Panel (Linear) ===
                gene_id = data.get('gene_id', '')
                for tx_data in self.gene_transcription_data:
                    if tx_data.get('gene_id') == gene_id:
                        tx_data['protein_count'] = data.get('protein_count', 0.0)
                        tx_data['protein_level'] = data['current_level']
                        tx_data['translation_rate'] = data.get('rate', 0.1)
                        # Ensure non-zero sync
                        tx_data['protein_half_life'] = data.get('protein_half_life') or default_protein_half_life or 3600.0
                        tx_data['protein_sequence'] = data.get('cds_sequence', 'M')
                        break

                # Update mRNA column
                mrna_cur = int(data.get('mrna_current', 0))
                mrna_tot = data.get('mrna_copies', 0)
                self.translation_table.setItem(i, 2, QTableWidgetItem(f"{mrna_cur}/{mrna_tot}"))
                
                # Update Proteins column
                prot_cur = int(data.get('protein_current', 0))
                prot_tot = data.get('protein_total', 0)
                self.translation_table.setItem(i, 3, QTableWidgetItem(f"{prot_cur}/{prot_tot}"))
                
                # Update status and current level
                self.translation_table.setItem(i, 10, QTableWidgetItem(status_icon))
                self.translation_table.setItem(i, 11, QTableWidgetItem(f"{data['current_level']:.1f}%"))
                
                # Broadcast to Live Data Viewer and ZMQ
                self._broadcast_data(data)
                
                for j in range(12):
                    if self.translation_table.item(i, j):
                        self.translation_table.item(i, j).setBackground(color)
            else:
                all_complete = False
                self.translation_table.setItem(i, 10, QTableWidgetItem("⏸️ Waiting"))
                for j in range(12):
                    if self.translation_table.item(i, j):
                        self.translation_table.item(i, j).setBackground(QColor(240, 240, 240))
        
        # ═══════ GENE DYNAMICS PANEL UPDATE ═══════
        # Push synced data to the dynamics monitors
        gene_dynamics = getattr(self, 'gene_dynamics_panel', None)
        if gene_dynamics is not None and gene_dynamics.isVisible():
            # Pass transcription data (which now has synced protein counts)
            gene_dynamics.update_from_simulation(
                self.gene_transcription_data, self.translation_animation_time
            )

        # In LINEAR mode only, check if all complete
        if not is_continuous and all_complete and not getattr(self, '_all_complete_sent', False):
            self._all_complete_sent = True
            self._broadcast_status('completed')
            print("DEBUG: All genes complete - linear mode finished")

