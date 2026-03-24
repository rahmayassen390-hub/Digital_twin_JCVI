"""
Cell Simulation Engine for JCVI Genome Analyzer
================================================

Continuous time-dependent simulation of bacterial cell processes.
Models mRNA and protein dynamics with synthesis AND degradation.

Key concepts:
- mRNA is continuously synthesized (transcription) and degraded
- Proteins are continuously synthesized (translation) from mRNA and degraded  
- System reaches steady-state where synthesis = degradation
- Pausing transcription causes mRNA to decay, which affects translation
- Biological noise adds stochastic fluctuations to gene expression
- Cellular resources (ribosomes, ATP) limit concurrent translation
"""

import math
import random
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Callable, Tuple, Set
from enum import Enum
import sys
import os
import re

# Ensure local directory is in path for metabolic_solver
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
try:
    from metabolic_solver import MetabolicSolver
    METABOLIC_AVAILABLE = True
except ImportError:
    MetabolicSolver = None
    METABOLIC_AVAILABLE = False


class SimulationMode(Enum):
    """Simulation mode selection"""
    LINEAR = "linear"          # Original one-time pipeline (genes reach 100% and stop)
    CONTINUOUS = "continuous"  # Living cell (continuous synthesis + degradation)


# =============================================================================
# JCVI syn3.0 SPECIFIC CONSTANTS (Issue 5: JCVI Specificity Enforcement)
# =============================================================================
# These constants are derived from JCVI syn3.0 minimal cell literature and data
# Do NOT use E. coli or other bacterial assumptions without explicit flagging

# JCVI syn3.0 Genome Characteristics
JCVI_GENOME_ID = "JCVISYN3A"           # Official NCBI accession prefix
JCVI_CHROMOSOME_ID = "CP016816"        # NCBI chromosome ID
JCVI_GENOME_SIZE = 543379              # bp (minimal genome)
JCVI_GENE_COUNT_EXPECTED = 493         # Protein-coding genes
JCVI_DOUBLING_TIME = 1.75 * 3600       # ~105 min (seconds)

# JCVI-Specific Kinetic Parameters (from literature)
# Source: Hutchison et al., 2016 Science; Breuer et al., 2019
JCVI_MRNA_HALF_LIFE = 300.0            # 5 minutes - typical for minimal cell
JCVI_PROTEIN_HALF_LIFE = 3600.0        # 1 hour - stable proteins
JCVI_TRANSCRIPTION_RATE = 30.0         # nucleotides/second (slower than E. coli)
JCVI_TRANSLATION_RATE = 12.0           # amino acids/second (slower than E. coli)
JCVI_RIBOSOME_COUNT = 4500             # Scaled to achieve 105 min doubling time
JCVI_POLYMERASE_COUNT = 1200           # Scaled to maintain TX-TL ratio

# JCVI-Specific Regulatory Characteristics
# JCVI syn3.0 has MINIMAL regulation - mostly constitutive expression
JCVI_REGULATORY_COMPLEXITY = "minimal"  # Most genes are constitutive
JCVI_OPERON_COUNT = 50                  # Estimated operons in minimal cell
JCVI_DEFAULT_NOISE = 0.05               # Low noise for minimal cell
JCVI_CELL_DRY_WEIGHT = 25.0             # femtograms (fg) - minimal cell estimate
JCVI_AVOGADRO = 6.022e23

# Flags for non-JCVI assumptions
class ParameterSource:
    """Track where parameters come from for transparency"""
    JCVI_BED = "jcvi_bed_file"       # Derived from provided JCVI BED data
    JCVI_LITERATURE = "jcvi_lit"     # From JCVI syn3.0 literature
    ECOLI_ASSUMPTION = "ecoli_assumption"  # E. coli-like (flagged!)
    GENERIC_BACTERIAL = "generic"     # Generic bacterial (flagged!)
    USER_OVERRIDE = "user_override"   # User explicitly set
    COMPUTED = "computed"             # Computed from other JCVI data


def validate_jcvi_organism(genome_sequence: str = None, genes: list = None) -> dict:
    """
    Validate that loaded data is from JCVI syn3.0 genome.
    
    Returns dict with:
    - is_jcvi: bool - True if data appears to be from JCVI
    - confidence: float 0-1
    - warnings: list of non-JCVI assumptions detected
    """
    warnings = []
    confidence = 0.0
    
    if genome_sequence:
        # Check genome size
        size = len(genome_sequence)
        if abs(size - JCVI_GENOME_SIZE) < 1000:
            confidence += 0.4
        else:
            warnings.append(f"Genome size {size} differs from JCVI ({JCVI_GENOME_SIZE})")
    
    if genes:
        # Check gene count
        gene_count = len(genes)
        if 400 <= gene_count <= 550:
            confidence += 0.3
        else:
            warnings.append(f"Gene count {gene_count} differs from expected JCVI (~473)")
        
        # Check for JCVI gene ID patterns
        jcvi_pattern_count = sum(1 for g in genes if 'JCVISYN3' in (g.id or '').upper())
        if jcvi_pattern_count > len(genes) * 0.5:
            confidence += 0.3
        elif jcvi_pattern_count == 0:
            warnings.append("No JCVISYN3* gene IDs found - may not be JCVI genome")
    
    is_jcvi = confidence >= 0.5
    
    if not is_jcvi and not warnings:
        warnings.append("Unable to confirm JCVI syn3.0 genome - using generic parameters")
    
    return {
        'is_jcvi': is_jcvi,
        'confidence': confidence,
        'warnings': warnings,
        'genome_size': len(genome_sequence) if genome_sequence else 0,
        'gene_count': len(genes) if genes else 0
    }


# === CELLULAR RESOURCES (Issue 7 + Issue 3: Genome-Derived) ===
# Models resource constraints derived from JCVI genome content

# Gene patterns for resource inference
RIBOSOMAL_PROTEIN_PATTERNS = ['rps', 'rpl', 'rpm']  # Small (S), Large (L), other ribosomal
RNA_POLYMERASE_PATTERNS = ['rpoa', 'rpob', 'rpoc', 'rpod', 'rpoe', 'rpoz']
ATP_SYNTHASE_PATTERNS = ['atpa', 'atpb', 'atpc', 'atpd', 'atpe', 'atpf', 'atpg', 'atph']

@dataclass
class CellularResources:
    """
    Cellular resource pool for resource-limited simulation.
    
    JCVI syn3.0 minimal cell resources are DERIVED from genome content:
    - Ribosomes: scaled from ribosomal protein gene count
    - RNA Polymerases: scaled from rpo gene count
    - ATP: scaled from ATP synthase component count
    """
    # Total available resources (derived from genome or defaults)
    total_ribosomes: int = 2000         # Free ribosomes in cell
    total_rna_polymerases: int = 200    # RNA Pol for transcription
    total_atp: float = 1e6              # ATP molecules (regenerates)
    total_amino_acids: float = 1e7      # Amino acid pool
    
    # Currently in use (dynamically updated)
    ribosomes_in_use: int = 0
    polymerases_in_use: int = 0
    
    # Regeneration rates (scaled by ATP synthase capacity)
    atp_regeneration_rate: float = 1e4  # ATP/second
    aa_regeneration_rate: float = 1e3   # amino acids/second
    
    # === GENOME DERIVATION TRACKING (Issue 3) ===
    resource_source: str = "default"    # "genome", "default"
    ribosomal_gene_count: int = 0       # Number of rps*/rpl* genes found
    polymerase_gene_count: int = 0      # Number of rpo* genes found
    atp_synthase_gene_count: int = 0    # Number of atp* genes found
    
    @property
    def available_ribosomes(self) -> int:
        return max(0, self.total_ribosomes - self.ribosomes_in_use)
    
    @property
    def available_polymerases(self) -> int:
        return max(0, self.total_rna_polymerases - self.polymerases_in_use)
    
    @property
    def ribosome_saturation(self) -> float:
        """Fraction of ribosomes in use (0-1)"""
        return self.ribosomes_in_use / self.total_ribosomes if self.total_ribosomes > 0 else 0
    
    @property
    def polymerase_saturation(self) -> float:
        """Fraction of polymerases in use (0-1)"""
        return self.polymerases_in_use / self.total_rna_polymerases if self.total_rna_polymerases > 0 else 0
    
    def allocate_ribosomes(self, count: int) -> int:
        """Try to allocate ribosomes, return actual allocated"""
        allocated = min(count, self.available_ribosomes)
        self.ribosomes_in_use += allocated
        return allocated
    
    def release_ribosomes(self, count: int):
        """Release ribosomes back to pool"""
        self.ribosomes_in_use = max(0, self.ribosomes_in_use - count)
    
    def step(self, dt: float):
        """Regenerate resources each time step"""
        self.total_atp += self.atp_regeneration_rate * dt
        self.total_amino_acids += self.aa_regeneration_rate * dt


def derive_resources_from_genome(genes: list, translation_dynamics: dict = None) -> CellularResources:
    """
    Derive cellular resource pools from actual genome content.
    
    Uses gene annotations to count:
    - Ribosomal proteins (rpsA-rpsU, rplA-rplX) → ribosomes
    - RNA polymerase subunits (rpoA, rpoB, rpoC, rpoD) → polymerases  
    - ATP synthase components (atpA-atpH) → ATP regeneration capacity
    
    JCVI syn3.0 has ~50 ribosomal proteins, 4-5 rpo genes, 8 atp genes.
    Scale resources accordingly (minimal cell = reduced capacity).
    
    Args:
        genes: List of Gene objects from GFF
        translation_dynamics: Dict of TranslationDynamics (optional, for protein levels)
    
    Returns:
        CellularResources with genome-derived values
    """
    ribosomal_genes = []
    polymerase_genes = []
    atp_genes = []
    
    for gene in genes:
        gene_name = (gene.name or gene.id).lower()
        gene_product = (getattr(gene, 'product', '') or '').lower()
        
        # Check for ribosomal proteins
        is_ribosomal = False
        for pattern in RIBOSOMAL_PROTEIN_PATTERNS:
            if gene_name.startswith(pattern) or pattern in gene_product:
                is_ribosomal = True
                break
        
        # JCVI-specific: S-prefix or L-prefix for ribosomal proteins
        if not is_ribosomal:
            if 'ribosomal protein' in gene_product:
                is_ribosomal = True
            elif gene_product.startswith('s') and any(f's{i}:' in gene_product for i in range(1, 32)):
                is_ribosomal = True
            elif gene_product.startswith('l') and any(f'l{i}:' in gene_product for i in range(1, 40)):
                is_ribosomal = True
                
        if is_ribosomal:
            ribosomal_genes.append(gene.id)
            continue # Move to next gene
        
        # Check for RNA polymerase subunits
        is_polymerase = False
        for pattern in RNA_POLYMERASE_PATTERNS:
            if pattern in gene_name or 'rna polymerase' in gene_product:
                is_polymerase = True
                break
        
        if is_polymerase:
            polymerase_genes.append(gene.id)
            continue
        
        # Check for ATP synthase
        is_atp = False
        for pattern in ATP_SYNTHASE_PATTERNS:
            if gene_name.startswith(pattern) or 'atp synthase' in gene_product:
                is_atp = True
                break
                
        if is_atp:
            atp_genes.append(gene.id)
            continue
    
    # Calculate resource pools based on gene counts
    # JCVI syn3.0 scaling factors (reduced from E. coli)
    
    # Ribosomes: ~50 ribosomal protein genes in JCVI → ~4500 ribosomes
    # Paper target: 105.2 min doubling time requires ~4500 functional ribosomes
    ribosomal_count = len(set(ribosomal_genes))
    if ribosomal_count > 0:
        # Scale: ~88 ribosomes per ribosomal protein gene for minimal cell
        # (accounts for multi-copy ribosomal genes and high translation demand)
        total_ribosomes = ribosomal_count * 88
    else:
        total_ribosomes = JCVI_RIBOSOME_COUNT  # Use calibrated constant
    
    # RNA Polymerases: ~4-5 rpo genes in JCVI → ~50-100 polymerases
    # (E. coli has ~4 rpo genes and ~2000 polymerases)
    polymerase_count = len(set(polymerase_genes))
    if polymerase_count > 0:
        # Scale: ~200 polymerases per rpo gene for minimal cell
        total_polymerases = polymerase_count * 200
    else:
        total_polymerases = JCVI_POLYMERASE_COUNT  # Use calibrated constant
    
    # ATP regeneration: scaled by ATP synthase genes
    atp_gene_count = len(set(atp_genes))
    if atp_gene_count > 0:
        # Scale ATP regeneration by synthase component completeness
        # Full synthase (8 components) = full rate, partial = reduced
        atp_completeness = min(1.0, atp_gene_count / 8.0)
        atp_regen_rate = 5e4 * atp_completeness  # Higher ATP flux for fast growth
    else:
        atp_regen_rate = 1e4  # Default
    
    resources = CellularResources(
        total_ribosomes=total_ribosomes,
        total_rna_polymerases=total_polymerases,
        total_atp=5e5,  # Minimal cell ATP pool
        total_amino_acids=5e6,  # Minimal cell AA pool
        atp_regeneration_rate=atp_regen_rate,
        aa_regeneration_rate=5e2, # Reduced for minimal cell
        resource_source="genome",
        ribosomal_gene_count=ribosomal_count,
        polymerase_gene_count=polymerase_count,
        atp_synthase_gene_count=atp_gene_count
    )
    
    print(f"DEBUG: Derived resources from genome:")
    print(f"  - Ribosomal genes: {ribosomal_count} → {total_ribosomes} ribosomes")
    print(f"  - Polymerase genes: {polymerase_count} → {total_polymerases} polymerases")
    print(f"  - ATP synthase genes: {atp_gene_count} → {atp_regen_rate:.0f} ATP/sec")
    
    return resources


@dataclass
class GeneCellState:
    """
    Dynamic state of a gene's molecular products (Digital Twin).
    
    Tracks current mRNA and protein counts as they fluctuate over time.
    Uses first-order kinetics for both synthesis and degradation.
    """
    gene_id: str
    mrna_count: float = 0.0          # Number of mRNA molecules
    protein_count: float = 0.0       # Number of protein molecules
    transcription_rate: float = 0.0   # Current synthesis rate
    translation_rate: float = 0.0     # Current synthesis rate
    mrna_half_life: float = 1.0       # Minutes
    protein_half_life: float = 24.0   # Hours
    status: str = "stable"           # "stable", "transcribing", "translation_paused"
    
    # Steady-state targets (calculated from rates and half-lives)
    mrna_steady_state: float = 0.0
    protein_steady_state: float = 0.0
    
    # ... cascade dynamics ...
    
    # Cascade dynamics
    is_transcription_paused: bool = False
    is_translation_paused: bool = False
    
    # Internal tracking for steady-state and knockouts
    _original_transcription_rate: float = 0.0
    _knockout_modifier: float = 1.0
    _affected_by_genes: Set[str] = field(default_factory=set)
    
    # === GENE-SPECIFIC NOISE (from regulation BED) ===
    gene_noise_level: float = 0.05    # Per-gene noise (0-1), derived from regulation complexity

    # === MASS BALANCE TRACKING (for Digital Twin transparency) ===
    # Last step production and degradation values
    last_mrna_production: float = 0.0
    last_mrna_degradation: float = 0.0
    last_protein_production: float = 0.0
    last_protein_degradation: float = 0.0
    
    def __post_init__(self):
        """Calculate steady-state values after initialization"""
        self._calculate_steady_states()
    
    def _calculate_steady_states(self):
        """
        Calculate steady-state molecule counts.
        
        At steady-state: synthesis rate = degradation rate
        For first-order decay: degradation_rate = decay_constant * count
        Therefore: steady_state = synthesis_rate / decay_constant
        
        decay_constant = ln(2) / half_life
        """
        ln2 = math.log(2)
        
        # mRNA steady state
        if self.mrna_half_life > 0:
            mrna_decay_constant = ln2 / self.mrna_half_life
            self.mrna_steady_state = self.transcription_rate / mrna_decay_constant
        else:
            self.mrna_steady_state = 0.0
        
        if self.protein_half_life > 0:
            protein_decay_constant = ln2 / self.protein_half_life
            # Protein synthesis rate = translation_rate * mrna_count
            # At mRNA steady state: protein_synthesis = translation_rate * mrna_steady_state
            protein_synthesis_at_steady = self.translation_rate * self.mrna_steady_state
            self.protein_steady_state = protein_synthesis_at_steady / protein_decay_constant
        else:
            self.protein_steady_state = 0.0
            
        # DEBUG: Print for the first few genes added
        if self.gene_id == 'JCVISYN3A_0001':
            print(f"DEBUG: {self.gene_id} SS: mRNA={self.mrna_steady_state:.2f}, Prot={self.protein_steady_state:.1f} (TX={self.transcription_rate}, TL={self.translation_rate})")
    
    @property
    def mrna_decay_constant(self) -> float:
        """Decay constant for mRNA (1/seconds)"""
        if self.mrna_half_life > 0:
            return math.log(2) / self.mrna_half_life
        return 0.0
    
    @property
    def protein_decay_constant(self) -> float:
        """Decay constant for protein (1/seconds)"""
        if self.protein_half_life > 0:
            return math.log(2) / self.protein_half_life
        return 0.0
    
    @property
    def mrna_percentage(self) -> float:
        """Current mRNA as percentage of steady state"""
        if self.mrna_steady_state > 0:
            return min(100.0, (self.mrna_count / self.mrna_steady_state) * 100)
        return 0.0
    
    @property
    def protein_percentage(self) -> float:
        """Current protein as percentage of steady state"""
        if self.protein_steady_state > 0:
            return min(100.0, (self.protein_count / self.protein_steady_state) * 100)
        return 0.0


class CellSimulationEngine:
    """
    Continuous cell simulation engine.
    
    Manages the time-evolution of all genes in the cell,
    applying synthesis and degradation at each time step.
    
    Usage:
        engine = CellSimulationEngine()
        engine.add_gene("GENE1", transcription_rate=2.0, translation_rate=0.5)
        
        # In animation loop:
        engine.step(dt=0.1)  # Advance by 0.1 seconds
        state = engine.get_gene_state("GENE1")
        print(f"mRNA: {state.mrna_count}, Protein: {state.protein_count}")
    """
    
    def __init__(self, mode: SimulationMode = SimulationMode.CONTINUOUS):
        self.mode = mode
        self.gene_states: Dict[str, GeneCellState] = {}
        self.simulation_time: float = 0.0
        self.is_running: bool = False
        
        # Global half-life defaults (can be overridden per gene)
        self.default_mrna_half_life: float = 300.0     # 5 minutes
        self.default_protein_half_life: float = 3600.0  # 1 hour
        
        # === BIOLOGICAL NOISE PARAMETERS (Issue 7) ===
        self.noise_level: float = 0.1         # 0.0 = deterministic, 1.0 = high noise
        self.enable_noise: bool = True        # Toggle stochastic simulation
        
        # === CELLULAR RESOURCE CONSTRAINTS (Issue 7) ===
        self.enable_resource_limits: bool = True
        self.cellular_resources = CellularResources()
        
        # === METABOLIC COUPLING (Issue 2 + 5) ===
        self.enable_metabolic_coupling: bool = True
        self.metabolic_solver: Optional[MetabolicSolver] = None
        self.fba_update_interval: float = 60.0  # Run FBA every 60 simulation seconds
        self.last_fba_update: float = -1e9  # Force immediate update on first step
        self.current_fluxes: Dict[str, float] = {}
        
        # Callbacks for state changes
        self._on_state_change: Optional[Callable[[str, GeneCellState], None]] = None
    
    def set_noise_level(self, level: float):
        """Set noise level (0.0 = deterministic, 1.0 = high noise)"""
        self.noise_level = max(0.0, min(1.0, level))
    
    def set_mode(self, mode: SimulationMode):
        """Switch simulation mode"""
        self.mode = mode
    
    def add_gene(
        self,
        gene_id: str,
        transcription_rate: float = 1.0,
        translation_rate: float = 0.1,
        mrna_half_life: Optional[float] = None,
        protein_half_life: Optional[float] = None,
        initial_mrna: float = 0.0,
        initial_protein: float = 0.0
    ):
        """
        Add a gene to the simulation.
        
        Args:
            gene_id: Unique identifier for the gene
            transcription_rate: mRNA molecules produced per second
            translation_rate: Proteins per mRNA per second
            mrna_half_life: Time for mRNA to decay by 50% (uses default if None)
            protein_half_life: Time for protein to decay by 50% (uses default if None)
            initial_mrna: Starting mRNA count
            initial_protein: Starting protein count
        """
        state = GeneCellState(
            gene_id=gene_id,
            mrna_count=initial_mrna,
            protein_count=initial_protein,
            transcription_rate=transcription_rate,
            translation_rate=translation_rate,
            mrna_half_life=mrna_half_life if mrna_half_life is not None else self.default_mrna_half_life,
            protein_half_life=protein_half_life if protein_half_life is not None else self.default_protein_half_life,
            status="initializing"
        )
        self.gene_states[gene_id] = state
    
    def get_gene_state(self, gene_id: str) -> Optional[GeneCellState]:
        """Get current state of a gene"""
        return self.gene_states.get(gene_id)
    
    def pause_transcription(self, gene_id: str):
        """Pause transcription for a gene (mRNA will still degrade)"""
        if gene_id in self.gene_states:
            self.gene_states[gene_id].is_transcription_paused = True
            self.gene_states[gene_id].status = "paused"
    
    def resume_transcription(self, gene_id: str):
        """Resume transcription for a gene"""
        if gene_id in self.gene_states:
            self.gene_states[gene_id].is_transcription_paused = False
            self.gene_states[gene_id].status = "rising"
    
    def _apply_noise(self, value: float, noise_scale: float = 1.0, gene_noise_level: float = None) -> float:
        """
        Apply Gaussian noise to a value.
        
        Uses gene-specific noise level if provided, otherwise falls back to global noise.
        """
        if not self.enable_noise:
            return value
        
        # Use gene-specific noise if provided, else global
        effective_noise = gene_noise_level if gene_noise_level is not None else self.noise_level
        
        if effective_noise <= 0:
            return value
        
        # Gaussian noise with coefficient of variation = noise_level
        # CV = std_dev / mean, so std_dev = CV * mean
        std_dev = effective_noise * value * noise_scale
        noise = random.gauss(0, std_dev) if std_dev > 0 else 0
        return max(0.0, value + noise)
    
    def step(self, dt: float):
        """
        Advance simulation by dt seconds.
        """
        self.simulation_time += dt
        
        # Update cellular resources
        if self.enable_resource_limits:
            self.cellular_resources.step(dt)
            
        # 1. Update Metabolic Fluxes (Dynamic Coupling)
        if self.enable_metabolic_coupling and self.metabolic_solver:
            if self.simulation_time - self.last_fba_update >= self.fba_update_interval:
                self._update_fluxes()
                self.last_fba_update = self.simulation_time
        
        # 2. Step each gene dynamics
        for gene_id, state in self.gene_states.items():
            self._step_continuous(state, dt)
            
            # Notify listeners of state change
            if self._on_state_change:
                self._on_state_change(gene_id, state)

    def _update_fluxes(self):
        """
        Run FBA to update metabolic fluxes based on current enzyme (protein) levels.
        Implements: Vmax_j = kcat_j * [Enzyme_j]
        """
        if not self.metabolic_solver:
            return
            
        constraints = {}
        # Iterate over all reactions in the metabolic model
        for idx, rxn_id in enumerate(self.metabolic_solver.reactions):
            # 1. Start with default bounds (thermodynamics + medium constraints)
            orig_lb = self.metabolic_solver.lb[idx]
            orig_ub = self.metabolic_solver.ub[idx]
            
            # 2. Check if this reaction has an associated GPR (Gene-Protein-Reaction rule)
            if rxn_id in self.metabolic_solver.reaction_gprs:
                gpr = self.metabolic_solver.reaction_gprs[rxn_id]
                vmax_gene = 0.0
                
                # Extract all gene IDs referenced in this GPR
                # Supports MMSYN1_xxxx, JCVISYN3_xxxx, JCVISYN3A_xxxx
                referenced_genes = re.findall(r'(?:MMSYN1_|JCVISYN3A_|JCVISYN3)_?\d+', gpr)
                
                for r_gid in referenced_genes:
                    # Normalize to engine's internal ID format (JCVISYN3A_xxxx)
                    # Extract only the numeric part (usually 4 digits) to avoid '1' from MMSYN1
                    m = re.search(r'(\d{4})', r_gid)
                    if m:
                        norm_id = f"JCVISYN3A_{m.group(1)}"
                    else:
                        norm_id = r_gid
                    
                    state = self.get_gene_state(norm_id)
                    if not state:
                        # Fallback: check with possible direct matches
                        state = self.get_gene_state(r_gid)
                    
                    if state:
                        # Catalytic capacity calculation
                        # kcat: calibrated to reach 105 min doubling time
                        # Glucose uptake and glycolysis are historically fast
                        kcat_val = 450.0 if rxn_id in ['PGI', 'GLCpts', 'GAPD'] else 100.0
                        
                        # Units: molecules / cell / sec
                        cap_molecules_sec = kcat_val * state.protein_count
                        
                        # Unit Conversion: molecules/sec -> mmol / gDW / hr
                        # (molecules/sec * 3600 sec/hr) / (N_A * CDW_g) / 1000 [mol->mmol cancels 10^3]
                        cdw_g = JCVI_CELL_DRY_WEIGHT * 1e-15
                        conv = 3600.0 / (JCVI_AVOGADRO * cdw_g)
                        # Avogadro is molecules/mol. To get mmol, we divide mol by 1000.
                        # So: molec/sec * (1/N_A) * 1000 [mmol/mol] * 3600 [sec/hr] / CDW
                        vmax_gene += (cap_molecules_sec * conv) * 1000.0
                
                # Minimum leak for modeled reactions to prevent total blockage due to noise
                vmax_gene = max(vmax_gene, 0.1)
                
                # Constrain the original bounds by the calculated catalytic capacity
                # Ensures enzymes can only slow down a reaction, not violate thermodynamics (LB/UB)
                if orig_lb >= 0:
                    new_ub = min(orig_ub, vmax_gene)
                    constraints[rxn_id] = (orig_lb, new_ub)
                else:
                    # Reversible: cap both directions
                    new_lb = max(orig_lb, -vmax_gene)
                    new_ub = min(orig_ub, vmax_gene)
                    constraints[rxn_id] = (new_lb, new_ub)
                
                # No sampling for speed
                pass
            else:
                # 3. No GPR (e.g., Exchange reactions, non-enzymatic diffusion)
                # Keep original bounds (e.g., global medium constraints for Glucose)
                constraints[rxn_id] = (orig_lb, orig_ub)
                
        # solve FBA with updated capacity constraints
        result = self.metabolic_solver.solve(constraints)
        if result["status"] == "success":
            self.current_fluxes = result["fluxes"]
            # Trigger potential downstream events (cell growth, etc.)
    
    def _step_linear(self, state: GeneCellState, dt: float):
        """
        Linear mode: genes accumulate to 100% and stop.
        This preserves the original behavior for compatibility.
        """
        if state.is_transcription_paused:
            return
        
        # Simple exponential approach to steady state
        if state.mrna_count < state.mrna_steady_state * 0.99:
            # Growth phase
            growth_rate = state.transcription_rate * 0.01  # Slower for animation
            state.mrna_count += growth_rate * dt
            state.mrna_count = min(state.mrna_count, state.mrna_steady_state)
            state.status = "rising"
        else:
            state.mrna_count = state.mrna_steady_state
            state.status = "complete"
        
        # Protein follows mRNA with delay
        if state.protein_count < state.protein_steady_state * 0.99:
            protein_rate = state.translation_rate * state.mrna_count * 0.01
            state.protein_count += protein_rate * dt
            state.protein_count = min(state.protein_count, state.protein_steady_state)
        else:
            state.protein_count = state.protein_steady_state
    
    def _get_derivatives(self, state: GeneCellState, mrna: float, protein: float) -> Tuple[float, float, float, float]:
        """
        Calculate production and degradation rates for current state.
        Returns: (mrna_prod, mrna_deg, prot_prod, prot_deg)
        """
        # mRNA Synthesis (only if not paused)
        if not state.is_transcription_paused:
            mrna_prod = state.transcription_rate
        else:
            mrna_prod = 0.0
        
        # mRNA Degradation
        mrna_deg = state.mrna_decay_constant * mrna
        
        # Protein Synthesis (constrained by mRNA and resources)
        if not state.is_translation_paused and mrna > 0:
            prot_prod = state.translation_rate * mrna
            if self.enable_resource_limits:
                ribosome_fraction = 1.0 - self.cellular_resources.ribosome_saturation
                prot_prod *= max(0.1, ribosome_fraction)
        else:
            prot_prod = 0.0
            
        # Protein Degradation
        prot_deg = state.protein_decay_constant * protein
        
        return mrna_prod, mrna_deg, prot_prod, prot_deg

    def _step_continuous(self, state: GeneCellState, dt: float):
        """
        Digital Twin mass-balance simulation using Runge-Kutta 4th Order (RK4).
        
        Explicit ODEs:
            dm/dt = f(m, p) = TX_rate - k_deg_m * m
            dp/dt = g(m, p) = TL_rate * m - k_deg_p * p
        """
        m0, p0 = state.mrna_count, state.protein_count
        
        def f(m, p):
            mp, md, pp, pd = self._get_derivatives(state, m, p)
            return mp - md, pp - pd
            
        # RK4 Coefficients
        k1m, k1p = f(m0, p0)
        k2m, k2p = f(m0 + 0.5 * dt * k1m, p0 + 0.5 * dt * k1p)
        k3m, k3p = f(m0 + 0.5 * dt * k2m, p0 + 0.5 * dt * k2p)
        k4m, k4p = f(m0 + dt * k3m, p0 + dt * k3p)
        
        # Compute combined deterministic delta
        dm_det = (dt / 6.0) * (k1m + 2*k2m + 2*k3m + k4m)
        dp_det = (dt / 6.0) * (k1p + 2*k2p + 2*k3p + k4p)
        
        # Apply biological noise to the increments
        if self.enable_noise:
            dm = self._apply_noise(dm_det, noise_scale=0.5, gene_noise_level=state.gene_noise_level)
            dp = self._apply_noise(dp_det, noise_scale=0.5, gene_noise_level=state.gene_noise_level)
        else:
            dm, dp = dm_det, dp_det
            
        # Update counts
        state.mrna_count = max(0.0, m0 + dm)
        p_new = p0 + dp
        if math.isnan(p_new):
            p_new = p0 # Fallback if noise/derivs fail
        state.protein_count = max(0.0, p_new)
        
        # Store last rates for transparency (use k1 as instantaneous proxy)
        state.last_mrna_production, state.last_mrna_degradation, \
        state.last_protein_production, state.last_protein_degradation = self._get_derivatives(state, m0, p0)
        
        # Adjust for display by multiplying by dt to show "per step" change
        state.last_mrna_production *= dt
        state.last_mrna_degradation *= dt
        state.last_protein_production *= dt
        state.last_protein_degradation *= dt
        
        # === DETERMINISTIC STATE FROM NET BALANCE ===
        state.status = self._determine_state_from_balance(
            mrna_prod=state.last_mrna_production,
            mrna_deg=state.last_mrna_degradation,
            prot_prod=state.last_protein_production,
            prot_deg=state.last_protein_degradation,
            mrna_count=state.mrna_count,
            protein_count=state.protein_count,
            is_paused=state.is_transcription_paused
        )
    
    def _determine_state_from_balance(self, mrna_prod, mrna_deg, prot_prod, prot_deg,
                                       mrna_count, protein_count, is_paused):
        """
        Deterministic state based on net mass balance.
        
        Rising: production > degradation
        Declining: degradation > production  
        Steady: production ≈ degradation (within 5%)
        Depleted: protein = 0
        
        Biological constraints enforced:
        - Protein = 0 cannot be "declining" (must be "depleted")
        - mRNA = 0 with protein > 0 → "declining_fast"
        """
        # CONSTRAINT: Protein = 0 → depleted, never declining
        if protein_count <= 0.01:
            return "depleted"
        
        # CONSTRAINT: No mRNA but still have protein → declining fast
        if mrna_count <= 0.01 and protein_count > 0.01:
            return "declining_fast"
        
        # Paused transcription with remaining molecules → declining
        if is_paused:
            return "declining"
        
        # Calculate net protein balance
        net_protein = prot_prod - prot_deg
        
        # Threshold for "steady" = within 1% of current protein count
        steady_threshold = 0.01 * protein_count + 0.001
        
        if abs(net_protein) < steady_threshold:
            return "steady_state"
        elif net_protein > 0:
            return "rising"
        else:
            return "declining"
    
    def reset(self):
        """Reset all genes to initial state"""
        self.simulation_time = 0.0
        for state in self.gene_states.values():
            state.mrna_count = 0.0
            state.protein_count = 0.0
            state.is_transcription_paused = False
            state.is_translation_paused = False
            state.status = "initializing"
    
    def get_statistics(self) -> Dict:
        """Get simulation statistics"""
        if not self.gene_states:
            return {}
        
        mrna_counts = [s.mrna_count for s in self.gene_states.values()]
        protein_counts = [s.protein_count for s in self.gene_states.values()]
        
        steady_count = sum(1 for s in self.gene_states.values() if s.status == "steady_state")
        rising_count = sum(1 for s in self.gene_states.values() if s.status == "rising")
        paused_count = sum(1 for s in self.gene_states.values() if s.is_transcription_paused)
        
        return {
            'simulation_time': self.simulation_time,
            'total_genes': len(self.gene_states),
            'steady_state_genes': steady_count,
            'rising_genes': rising_count,
            'paused_genes': paused_count,
            'total_mrna': sum(mrna_counts),
            'total_protein': sum(protein_counts),
            'avg_mrna_percentage': sum(s.mrna_percentage for s in self.gene_states.values()) / len(self.gene_states),
            'avg_protein_percentage': sum(s.protein_percentage for s in self.gene_states.values()) / len(self.gene_states)
        }
    
    def set_on_state_change(self, callback: Callable[[str, GeneCellState], None]):
        """Set callback for state changes"""
        self._on_state_change = callback
    
    # =========================================================================
    # KNOCKOUT CASCADE EFFECTS
    # =========================================================================
    
    def apply_knockout_cascade(self, knocked_gene: str, cascade_effects: Dict):
        """
        Apply knockout cascade effects to dependent genes.
        
        When a gene is knocked out, this modifies the transcription rates
        of all dependent genes based on the cascade effect strength.
        
        Args:
            knocked_gene: The gene that was knocked out
            cascade_effects: Result from GeneDependencyGraph.calculate_cascade()
                Expected keys:
                - same_cluster_genes: List[(gene_id, effect_strength)]
                - connected_genes: List[(gene_id, effect_strength, distance, cluster_name)]
        """
        affected_genes = []
        
        # Process same cluster genes (high effect)  
        for gene_id, effect in cascade_effects.get('same_cluster_genes', []):
            if gene_id in self.gene_states:
                state = self.gene_states[gene_id]
                # Reduce transcription rate by effect strength
                # Effect of 0.8 means 80% reduction (20% of normal)
                modifier = max(0.0, 1.0 - effect)
                self._apply_knockout_modifier(gene_id, modifier, knocked_gene, 
                                               cascade_effects.get('cluster_name', ''))
                affected_genes.append(gene_id)
        
        # Process connected cluster genes (distance-based effect)
        for entry in cascade_effects.get('connected_genes', []):
            gene_id, effect, distance, cluster_name = entry
            if gene_id in self.gene_states:
                modifier = max(0.0, 1.0 - effect)
                self._apply_knockout_modifier(gene_id, modifier, knocked_gene, cluster_name)
                affected_genes.append(gene_id)
        
        return affected_genes
    
    def _apply_knockout_modifier(self, gene_id: str, modifier: float, 
                                   source_gene: str, source_cluster: str):
        """Apply knockout modifier to a single gene."""
        if gene_id not in self.gene_states:
            return
        
        state = self.gene_states[gene_id]
        
        # Store original rate if not already affected
        if not hasattr(state, '_original_transcription_rate'):
            state._original_transcription_rate = state.transcription_rate
        
        # Combine with existing modifier (multiple knockouts can stack)
        current_modifier = getattr(state, '_knockout_modifier', 1.0)
        new_modifier = min(current_modifier, modifier)  # Take the strongest effect
        
        state._knockout_modifier = new_modifier
        state.transcription_rate = state._original_transcription_rate * new_modifier
        
        # Track which genes caused the effect
        if not hasattr(state, '_affected_by_genes'):
            state._affected_by_genes = set()
        state._affected_by_genes.add(source_gene)
        
        # Update status to show effect
        if new_modifier < 1.0:
            state.status = "affected_by_knockout"
    
    def clear_knockout_effects(self, resumed_gene: str):
        """
        Clear knockout effects when a gene is resumed.
        
        Restores transcription rates for genes that were affected
        by the resumed gene's knockout.
        
        Args:
            resumed_gene: The gene that was resumed
            
        Returns:
            List of gene IDs that were cleared
        """
        cleared_genes = []
        
        for gene_id, state in self.gene_states.items():
            if hasattr(state, '_affected_by_genes') and resumed_gene in state._affected_by_genes:
                state._affected_by_genes.remove(resumed_gene)
                
                # If no more knockout sources, restore original rate
                if len(state._affected_by_genes) == 0:
                    if hasattr(state, '_original_transcription_rate'):
                        state.transcription_rate = state._original_transcription_rate
                        state._knockout_modifier = 1.0
                    if state.status == "affected_by_knockout":
                        state.status = "rising"  # Start recovering
                
                cleared_genes.append(gene_id)
        
        return cleared_genes
    
    def get_affected_genes(self) -> List[str]:
        """Get list of genes currently affected by knockouts."""
        affected = []
        for gene_id, state in self.gene_states.items():
            if hasattr(state, '_affected_by_genes') and len(state._affected_by_genes) > 0:
                affected.append(gene_id)
        return affected
    
    def get_knockout_info(self, gene_id: str) -> Dict:
        """Get knockout effect information for a gene."""
        if gene_id not in self.gene_states:
            return {}
        
        state = self.gene_states[gene_id]
        return {
            'gene_id': gene_id,
            'is_affected': hasattr(state, '_affected_by_genes') and len(state._affected_by_genes) > 0,
            'affected_by': list(getattr(state, '_affected_by_genes', set())),
            'modifier': getattr(state, '_knockout_modifier', 1.0),
            'original_rate': getattr(state, '_original_transcription_rate', state.transcription_rate),
            'current_rate': state.transcription_rate,
            'reduction_pct': (1.0 - getattr(state, '_knockout_modifier', 1.0)) * 100
        }



def create_engine_from_transcription_data(
    transcription_data: Dict,
    translation_data: Dict = None,
    mode: SimulationMode = SimulationMode.CONTINUOUS,
    mrna_half_life: float = 300.0,
    protein_half_life: float = 3600.0
) -> CellSimulationEngine:
    """
    Create a simulation engine from BED file data.
    
    Args:
        transcription_data: Dict of gene_id -> TranscriptionDynamics
        translation_data: Dict of gene_id -> TranslationDynamics (optional)
        mode: Simulation mode (LINEAR or CONTINUOUS)
        mrna_half_life: Default mRNA half-life in seconds
        protein_half_life: Default protein half-life in seconds
    
    Returns:
        Configured CellSimulationEngine
    """
    engine = CellSimulationEngine(mode=mode)
    engine.default_mrna_half_life = mrna_half_life
    engine.default_protein_half_life = protein_half_life
    
    for gene_id, tx_data in transcription_data.items():
        # Derive transcription rate from level (higher level = higher rate)
        tx_level = getattr(tx_data, 'transcription_level', 50.0)
        transcription_rate = (tx_level / 100.0) * 2.0  # Scale to 0-2 mRNA/sec
        
        # Get translation data if available
        tl_level = 50.0
        rbs_strength = "medium"
        if translation_data and gene_id in translation_data:
            tl_data = translation_data[gene_id]
            tl_level = getattr(tl_data, 'protein_level', 50.0)
            rbs_strength = getattr(tl_data, 'rbs_strength', 'medium')
        
        # Derive translation rate from RBS strength
        rbs_multiplier = {'strong': 2.0, 'medium': 1.0, 'weak': 0.5}.get(rbs_strength.lower(), 1.0)
        translation_rate = (tl_level / 100.0) * 0.2 * rbs_multiplier  # Scale to proteins/mRNA/sec
        
        engine.add_gene(
            gene_id=gene_id,
            transcription_rate=transcription_rate,
            translation_rate=translation_rate
        )
    
    return engine


def create_engine_from_enriched_genes(
    enriched_genes: dict,
    mode: SimulationMode = SimulationMode.CONTINUOUS,
    metabolic_model_path: Optional[str] = None
) -> CellSimulationEngine:
    """
    Create simulation engine from enriched Gene objects (AI Virtual Cell bridge).
    
    Uses the AI-predicted or data-derived kinetic parameters from the
    DataEnrichmentLayer to configure mass-balance ODEs:
    
        d(mRNA)/dt  = TX_rate × RBS_modifier - ln(2)/mrna_half_life × mRNA
        d(Protein)/dt = TL_rate × mRNA       - ln(2)/protein_half_life × Protein
    
    Args:
        enriched_genes: Dict[gene_id, Gene] with enriched fields
        mode: Simulation mode
        metabolic_model_path: Path to SBML or Excel reconstruction
    
    Returns:
        CellSimulationEngine configured with enriched kinetics and metabolism
    """
    engine = CellSimulationEngine(mode=mode)
    
    # 0. Initialize Metabolic Solver if available
    if METABOLIC_AVAILABLE:
        model_path = metabolic_model_path
        if not model_path:
            # Try default paths
            defaults = [
                "../Data/metabolic_reconstruction.xlsx",
                "Data/metabolic_reconstruction.xlsx",
                "../Data/metabolic_model_iMB155.xml",
                "Data/metabolic_model_iMB155.xml"
            ]
            for d in defaults:
                if os.path.exists(d):
                    model_path = d
                    break
        
        if model_path and os.path.exists(model_path):
            try:
                engine.metabolic_solver = MetabolicSolver(model_path)
                engine.enable_metabolic_coupling = True
                engine.fba_update_interval = 250.0 # Standard interval for stability
                print(f"Metabolic Bridge: Attached solver using {model_path}")
            except Exception as e:
                print(f"Warning: Could not initialize MetabolicSolver: {e}")
    
    # Derive resources from gene list
    gene_list = list(enriched_genes.values())
    engine.cellular_resources = derive_resources_from_genome(gene_list)
    
    for gene_id, gene in enriched_genes.items():
        # === Transcription rate from mRNA expression or default ===
        if gene.mrna_expression > 0:
            # Scale to synthesis rate
            # Target average mass ~10^7 molecules total
            transcription_rate = (gene.mrna_expression / 50.0) * 0.02
        else:
            # Default constitutive expression for JCVI minimal cell
            transcription_rate = 0.02
        
        # === Translation rate modulated by RBS strength ===
        # RBS strength 1-10 maps to translation efficiency
        rbs_raw = getattr(gene, 'rbs_strength', 5.0)
        if math.isnan(rbs_raw) or rbs_raw <= 0:
            rbs_raw = 5.0
        rbs_modifier = rbs_raw / 5.0  # Normalize around 1.0
        # Scaling: around 0.5 prot/mRNA/sec to reach biological counts
        translation_rate = 0.5 * rbs_modifier
        if math.isnan(translation_rate):
            translation_rate = 0.5
        
        # === Half-lives from enriched data (AI-predicted or experimental) ===
        # DataEnrichmentLayer stores mrna_half_life in MINUTES → convert to SECONDS
        mrna_hl_seconds = gene.mrna_half_life * 60.0
        # protein_half_life in HOURS → convert to SECONDS
        protein_hl_seconds = gene.protein_half_life * 3600.0
        
        # === Gene-specific noise based on essentiality ===
        noise = 0.03 if gene.essentiality_status == "Essential" else 0.08
        
        # === Initialize with biologically realistic starting conditions ===
        state = GeneCellState(
            gene_id=gene_id,
            transcription_rate=transcription_rate,
            translation_rate=translation_rate,
            mrna_half_life=mrna_hl_seconds,
            protein_half_life=protein_hl_seconds,
            gene_noise_level=noise,
            status="initializing"
        )
        
        # Start at steady state to skip initialization lag in validation
        state.mrna_count = state.mrna_steady_state
        state.protein_count = state.protein_steady_state
        
        engine.gene_states[gene_id] = state
    
    total_iv_p = sum(s.protein_count for s in engine.gene_states.values())
    print(f"ODE Bridge: Created engine with {len(engine.gene_states)} genes. Total initial mass: {total_iv_p:.1f}")
    print(f"  Ribosome pool: {engine.cellular_resources.total_ribosomes}")
    print(f"  Polymerase pool: {engine.cellular_resources.total_rna_polymerases}")
    
    return engine

