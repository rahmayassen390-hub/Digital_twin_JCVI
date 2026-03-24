"""
Data Structures for JCVI Genome Analyzer
Contains all dataclass definitions used throughout the application
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional


@dataclass
class Gene:
    """Represents a gene with its properties and annotations"""
    id: str
    name: str
    start: int
    end: int
    strand: str
    type: str
    function: str = "Unknown"
    product: str = "Unknown"
    sequence: str = ""
    has_promoter: bool = False
    is_operon_internal: bool = False
    operon_id: Optional[str] = None
    
    # === NEW FIELDS FOR AI VIRTUAL CELL INTEGRATION ===
    mrna_expression: float = 0.0
    protein_abundance: float = 0.0
    kinetics: Dict = field(default_factory=dict)
    essentiality_status: str = "Unknown"  # Essential, Quasiessential, Nonessential
    function_class: str = "Unknown"
    ortholog_id: str = ""
    kegg_id: str = ""
    rbs_strength: float = 1.0
    mrna_half_life: float = 1.0      # in minutes
    protein_half_life: float = 25.0  # in hours
    molecular_ur: Dict = field(default_factory=dict)  # DNA, RNA, Protein embeddings


@dataclass
class CellStateVector:
    """
    Cellular Universal Representation (Cellular UR).
    Aggregates all gene-level Molecular URs with temporal dynamics,
    environmental context, and resource state into a unified cell-state vector.
    
    Dimensions:
        gene_embedding_aggregate: 1280-dim (attention-weighted mean of gene URs)
        dynamics_vector: 512-dim (compressed temporal state from ODE engine)
        environment_vector: 64-dim (environmental context)
        resource_vector: 32-dim (ribosome/polymerase/ATP saturation)
        combined_cellular_ur: 1888-dim (full concatenated cell state)
    """
    gene_embedding_aggregate: object = None    # np.ndarray (1280,)
    dynamics_vector: object = None             # np.ndarray (512,)
    environment_vector: object = None          # np.ndarray (64,)
    resource_vector: object = None             # np.ndarray (32,)
    combined_cellular_ur: object = None        # np.ndarray (1888,)
    
    # Metadata
    timestamp: float = 0.0                     # Simulation time when snapshot was taken
    num_genes: int = 0                         # Number of genes aggregated
    num_active_genes: int = 0                  # Genes with active transcription
    growth_phase: str = "unknown"              # lag, exponential, stationary, death


@dataclass
class PromoterAnnotation:
    """Represents a promoter annotation from BED file"""
    chrom: str
    start: int
    end: int
    name: str
    score: float
    strand: str
    annotation_type: str  # 'likely_promoter' or 'internal_operon'
    # New fields for updated likely promoter BED format
    gene_id: str = ""
    strand_info: str = ""  # BOTH, FWD, or REV
    motif_10_seq: str = ""  # -10 box motif sequence
    motif_10_score: int = 0  # -10 box score
    motif_35_seq: str = ""  # -35 box motif sequence
    motif_35_score: int = 0  # -35 box score
    gc_content: float = 0.0  # GC percentage


@dataclass
class BlastHit:
    """Represents a BLAST alignment hit"""
    query_id: str
    target_id: str
    identity: float
    alignment_length: int
    mismatches: int
    gap_opens: int
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    evalue: float
    bit_score: float
    query_seq: str
    target_seq: str
    genes: List[str]
    gene_functions: List[str]
    protein_sequence: str


@dataclass
class PromoterPrediction:
    """Represents AI prediction for promoter classification"""
    gene_id: str
    has_promoter: bool
    confidence: float
    promoter_strength: float
    motif_scores: Dict[str, float]
    classification: str  # 'own_promoter', 'operon_internal', 'ambiguous'


@dataclass
class TranscriptionState:
    """Represents predicted transcription state for a gene"""
    gene_id: str
    state: str  # 'silent', 'constitutive', 'conditional', 'induced'
    expression_level: float  # 0.0 to 1.0
    temporal_order: int  # Order of transcription
    confidence: float
    regulatory_factors: List[str]
    dependency_genes: List[str] = field(default_factory=list)  # Genes that must be expressed first


@dataclass
class TranslationState:
    """Represents predicted translation state for a gene"""
    gene_id: str
    protein_level: float  # 0-100%
    rbs_strength: str  # 'strong', 'medium', 'weak'
    tx_level: float  # Transcription level (if available)
    tl_level: float  # Translation level (if available)
    efficiency: float  # Translation efficiency (%)
    temporal_order: int  # Order of translation
    confidence: float
    status: str = "waiting"  # 'waiting', 'translating', 'complete'


@dataclass
class RegulationAnnotation:
    """Gene regulation annotation from Regulation Analysis BED"""
    gene_id: str
    regulation_type: str  # constitutive, operonMember, potentiallyInducible, potentiallyRepressible, unknown
    category: str  # Housekeeping, Inducible, Repressible, Unknown
    strength: str  # high, medium, low
    start: int
    end: int
    strand: str
    score: float


@dataclass
class TranscriptionDynamics:
    """Transcription dynamics from BED file - multiple tracks"""
    gene_id: str
    transcription_level: float  # Level from Track 1
    time_range: str = ""  # From Track 3 (e.g., "0-5min")
    phase: str = ""  # From Track 2 (e.g., "Phase1")
    category: str = ""  # From Track 2 (e.g., "DNA/Transcription")
    order: int = 0  # From Track 2 (transcription order)
    start: int = 0
    end: int = 0
    strand: str = "+"
    score: float = 0.0
    color: tuple = (128, 128, 128)


@dataclass
class TranslationDynamics:
    """Translation dynamics from BED file - multiple tracks"""
    gene_id: str
    protein_level: float  # From Track 1
    rbs_strength: str = ""  # From Track 3
    stability: str = ""  # From Track 3 & 4 (e.g., "stable", "unstable")
    phase: str = ""  # From Track 2 (e.g., "Phase2")
    category: str = ""  # From Track 2 (e.g., "Ribosome Assembly")
    order: int = 0  # From Track 2 (translation order)
    tx_level: float = 0.0  # From Track 6
    tl_level: float = 0.0  # From Track 6
    efficiency: float = 0.0  # From Track 6
    is_ncrna: bool = False  # From Track 5
    start: int = 0
    end: int = 0
    strand: str = "+"
    score: float = 0.0
    color: tuple = (128, 128, 128)
    cds_sequence: str = ""  # CDS nucleotide sequence from genome
    
    # === DATA-DRIVEN KINETIC PARAMETERS (Issue: Data-Driven Enforcement) ===
    # These should be parsed from BED Track 4: Protein_Stability
    protein_half_life: float = 0.0      # Parsed from HalfLife:XXmin in BED (0 = use default)
    mrna_half_life: float = 0.0         # If provided in data (0 = use default)
    
    # Parameter source tracking - flag when defaults are used
    kinetic_params_source: str = "default"  # "bed_file", "default", "computed"
    uses_default_params: bool = True        # True if any kinetic params use defaults
    
    # === GENE-SPECIFIC NOISE (Issue: Regulation-Driven Variability) ===
    # Derived from regulatory complexity in BED files
    gene_noise_level: float = 0.05      # 0.0-1.0, derived from regulation data (JCVI minimal = low)
    regulation_complexity: str = ""     # "constitutive", "operonMember", "inducible", etc.
    has_own_promoter: bool = True       # False if operon-internal (lower noise)
    num_regulators: int = 0             # Number of regulatory factors (JCVI = minimal)


@dataclass
class ContinuousCellState:
    """
    State for continuous cell simulation.
    
    Tracks dynamic mRNA and protein counts that fluctuate over time
    due to continuous synthesis and degradation.
    """
    gene_id: str
    
    # Current counts (dynamic, updated each time step)
    mrna_count: float = 0.0
    protein_count: float = 0.0
    
    # Synthesis rates (derived from BED file data)
    transcription_rate: float = 1.0      # mRNA molecules per second
    translation_rate: float = 0.1        # Proteins per mRNA per second
    
    # Degradation half-lives (biological constants)
    mrna_half_life: float = 300.0        # 5 minutes (typical bacterial mRNA)
    protein_half_life: float = 3600.0    # 1 hour (stable protein)
    
    # Steady-state targets (for display as percentages)
    mrna_steady_state: float = 0.0
    protein_steady_state: float = 0.0
    
    # Control flags
    is_paused: bool = False
    
    # Status for display
    status: str = "initializing"  # initializing, rising, steady_state, declining, paused
    
    # === KNOCKOUT CASCADE EFFECT FIELDS ===
    # Track how upstream gene knockouts affect this gene
    knockout_modifier: float = 1.0           # Multiplier from upstream knockouts (1.0 = normal, 0.0 = fully blocked)
    is_affected_by_knockout: bool = False    # True if this gene is affected by a knocked out gene
    affected_by_genes: List[str] = field(default_factory=list)  # List of knocked genes affecting this one
    knockout_effect_strength: float = 0.0    # Combined effect strength from all knockouts (0.0-1.0)
    knockout_effect_source_cluster: str = "" # Name of the cluster causing the knockout effect
