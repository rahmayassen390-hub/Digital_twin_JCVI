"""
Knockout Rules Engine for JCVI Genome Analyzer
===============================================

Scientific, rule-based gene knockout impact system for minimal bacterial genomes.

Implements 6 biological rules with confidence levels:
1. Operon/Polar Effect (HIGH) - Polycistronic transcription units
2. Core Cellular Process (MEDIUM) - Central dogma machinery
3. Metabolic Pathway (MEDIUM) - Short pathway bottlenecks
4. Regulatory Hub (MEDIUM) - RNA polymerase, ribosomes, chaperones
5. Genome Proximity (LOW) - Compact genome co-regulation
6. Essentiality/Stress (LOW-MEDIUM) - Global stress response

Each affected gene includes:
- Effect strength (0.0-1.0)
- Confidence level (high/medium/low)
- Reason explaining WHY it's affected
- Rule type (operon/core_process/pathway/hub/proximity/stress)
"""

from typing import Dict, List, Set, Tuple, Optional
from dataclasses import dataclass, field
from collections import defaultdict
import re

# Import cluster definitions for core process identification
try:
    from cluster_definitions import (
        CLUSTER_DEFINITIONS, ENTANGLEMENT_PAIRS,
        assign_gene_to_cluster, get_cluster_name
    )
    CLUSTERS_AVAILABLE = True
except ImportError:
    CLUSTERS_AVAILABLE = False
    CLUSTER_DEFINITIONS = {}
    ENTANGLEMENT_PAIRS = []


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class KnockoutEffect:
    """
    Represents the effect of a gene knockout on another gene.
    
    Attributes:
        gene_id: The affected gene identifier
        effect_strength: Impact severity (0.0-1.0, higher = stronger effect)
        confidence: Confidence level ("high", "medium", "low")
        reason: Human-readable explanation of why this gene is affected
        rule_type: Which rule triggered this effect
        distance: Graph distance (0=same operon, 1=same cluster, 2+=connected)
        cluster_name: Name of the gene's functional cluster
    """
    gene_id: str
    effect_strength: float
    confidence: str  # "high", "medium", "low"
    reason: str
    rule_type: str  # "operon", "core_process", "pathway", "hub", "proximity", "stress"
    distance: int = 0
    cluster_name: str = ""
    
    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization."""
        return {
            'gene_id': self.gene_id,
            'effect_strength': self.effect_strength,
            'confidence': self.confidence,
            'reason': self.reason,
            'rule_type': self.rule_type,
            'distance': self.distance,
            'cluster_name': self.cluster_name
        }


@dataclass
class GeneInfo:
    """Minimal gene information for knockout calculations."""
    gene_id: str
    gene_name: str = ""
    start: int = 0
    end: int = 0
    strand: str = "+"
    category: str = ""
    function: str = ""
    is_operon_internal: bool = False
    operon_id: Optional[str] = None
    cluster_id: int = -1


# =============================================================================
# CONSTANTS
# =============================================================================

# Core cellular process clusters (central dogma - highest priority)
CORE_PROCESS_CLUSTERS = {
    0: "Ribosome Assembly",
    1: "Translation",
    2: "tRNA Synthetases",
    3: "DNA Replication",
    4: "Transcription",
    5: "RNA Processing",
}

# Cluster importance multipliers
CLUSTER_IMPORTANCE = {
    0: 1.5,   # Ribosome Assembly - critical
    1: 1.5,   # Translation - critical
    2: 1.3,   # tRNA Synthetases
    3: 1.4,   # DNA Replication
    4: 1.5,   # Transcription - critical
    5: 1.2,   # RNA Processing
    7: 1.2,   # Cell Division
    11: 1.3,  # Protein Folding
    14: 1.3,  # Energy/ATP
}

# Regulatory hub gene patterns (affect many downstream genes)
HUB_GENE_PATTERNS = [
    (r'^rpo[A-Z]', "RNA polymerase", "Affects ALL transcription"),
    (r'^rps[A-Z]', "30S ribosomal protein", "Affects ALL translation"),
    (r'^rpl[A-Z]', "50S ribosomal protein", "Affects ALL translation"),
    (r'^rpm[A-Z]', "Ribosomal protein", "Affects translation"),
    (r'^gro[ELS]', "Chaperonin", "Affects protein folding globally"),
    (r'^dna[KJ]', "Chaperone", "Affects protein folding"),
    (r'^tig', "Trigger factor", "Affects nascent protein folding"),
    (r'^fts[ZAI]', "Cell division", "Affects cell division globally"),
    (r'^atp[A-H]', "ATP synthase", "Affects cellular energy"),
]

# Essential gene categories (knockout triggers stress response)
ESSENTIAL_CATEGORIES = [
    "Ribosome", "ribosomal", "Translation", "Transcription",
    "DNA replication", "tRNA", "aminoacyl", "ATP synthase"
]

# Stress response cluster
STRESS_CLUSTER = 18

# Short metabolic pathways in JCVI-syn (for pathway rule)
METABOLIC_PATHWAYS = {
    "glycolysis": {
        "genes": ["pgi", "pfkA", "fbaA", "gapA", "pgk", "gpmA", "eno", "pykA"],
        "description": "Glycolysis pathway"
    },
    "lipid_synthesis": {
        "genes": ["accA", "accB", "accC", "accD", "fabD", "fabH", "fabG", "fabZ", "fabI"],
        "description": "Fatty acid biosynthesis"
    },
    "nucleotide_synthesis": {
        "genes": ["pyrB", "pyrC", "pyrD", "pyrE", "pyrF"],
        "description": "Pyrimidine biosynthesis"
    }
}

# Maximum intergenic distance for operon detection (bp)
OPERON_MAX_DISTANCE = 200


# =============================================================================
# OPERON DETECTOR
# =============================================================================

class OperonDetector:
    """
    Detects operon membership and polar effects.
    
    In bacterial genomes, genes are often co-transcribed in polycistronic mRNA.
    When an upstream gene in an operon is knocked out, downstream genes
    experience "polar effects" - reduced or absent transcription.
    
    Detection logic:
    1. Use Gene.is_operon_internal and Gene.operon_id if available
    2. Otherwise, infer operons from strand + proximity (genes on same strand
       within ~200bp are likely co-transcribed)
    """
    
    def __init__(self):
        self.genes: Dict[str, GeneInfo] = {}
        self.operons: Dict[str, List[str]] = {}  # operon_id -> sorted gene list
        self.gene_to_operon: Dict[str, str] = {}  # gene_id -> operon_id
    
    def build_from_genes(self, genes: List[GeneInfo]):
        """Build operon map from gene list."""
        self.genes = {g.gene_id: g for g in genes}
        self.operons.clear()
        self.gene_to_operon.clear()
        
        # First, use explicit operon assignments if available
        for gene in genes:
            if gene.operon_id:
                if gene.operon_id not in self.operons:
                    self.operons[gene.operon_id] = []
                self.operons[gene.operon_id].append(gene.gene_id)
                self.gene_to_operon[gene.gene_id] = gene.operon_id
        
        # Sort each operon by genomic position
        for operon_id in self.operons:
            self.operons[operon_id].sort(
                key=lambda gid: self.genes[gid].start if gid in self.genes else 0
            )
        
        # Infer additional operons from proximity + strand
        self._infer_operons_from_proximity(genes)
    
    def _infer_operons_from_proximity(self, genes: List[GeneInfo]):
        """Infer operon membership from genomic proximity."""
        # Sort genes by position
        sorted_genes = sorted(genes, key=lambda g: g.start)
        
        inferred_operon_id = 1000  # Start from 1000 to avoid conflicts
        
        for i in range(len(sorted_genes) - 1):
            gene1 = sorted_genes[i]
            gene2 = sorted_genes[i + 1]
            
            # Skip if either gene already has an operon
            if gene1.gene_id in self.gene_to_operon or gene2.gene_id in self.gene_to_operon:
                continue
            
            # Check if on same strand and close together
            if gene1.strand == gene2.strand:
                distance = gene2.start - gene1.end
                if 0 <= distance <= OPERON_MAX_DISTANCE:
                    # These genes are likely co-transcribed
                    operon_id = f"inferred_operon_{inferred_operon_id}"
                    
                    if gene1.gene_id not in self.gene_to_operon:
                        self.gene_to_operon[gene1.gene_id] = operon_id
                        if operon_id not in self.operons:
                            self.operons[operon_id] = []
                        self.operons[operon_id].append(gene1.gene_id)
                    
                    # Use the same operon for gene2
                    current_operon = self.gene_to_operon.get(gene1.gene_id, operon_id)
                    self.gene_to_operon[gene2.gene_id] = current_operon
                    if gene2.gene_id not in self.operons.get(current_operon, []):
                        self.operons[current_operon].append(gene2.gene_id)
                    
                    inferred_operon_id += 1
    
    def get_polar_effects(self, knocked_gene_id: str) -> List[Tuple[str, float, str]]:
        """
        Get downstream genes affected by polar effects.
        
        Returns:
            List of (gene_id, effect_strength, reason) tuples
        """
        effects = []
        
        if knocked_gene_id not in self.gene_to_operon:
            return effects
        
        operon_id = self.gene_to_operon[knocked_gene_id]
        operon_genes = self.operons.get(operon_id, [])
        
        if len(operon_genes) <= 1:
            return effects
        
        knocked_gene = self.genes.get(knocked_gene_id)
        if not knocked_gene:
            return effects
        
        knocked_position = knocked_gene.start
        knocked_strand = knocked_gene.strand
        
        # Find downstream genes (based on strand direction)
        for i, gene_id in enumerate(operon_genes):
            if gene_id == knocked_gene_id:
                continue
            
            gene = self.genes.get(gene_id)
            if not gene:
                continue
            
            # Determine if this gene is downstream
            is_downstream = False
            if knocked_strand == '+':
                is_downstream = gene.start > knocked_position
            else:
                is_downstream = gene.start < knocked_position
            
            if is_downstream:
                # Calculate position offset in operon
                genes_in_operon = len(operon_genes)
                position_in_operon = operon_genes.index(gene_id)
                knocked_position_in_operon = operon_genes.index(knocked_gene_id)
                offset = abs(position_in_operon - knocked_position_in_operon)
                
                # Polar effect decreases with distance
                # First downstream gene: 95%, second: 90%, etc.
                effect = max(0.5, 0.95 - (0.05 * (offset - 1)))
                
                reason = f"Downstream in operon (polar effect, position +{offset})"
                effects.append((gene_id, effect, reason))
        
        return effects


# =============================================================================
# PATHWAY MAPPER
# =============================================================================

class PathwayMapper:
    """
    Maps genes to metabolic pathways for bottleneck detection.
    
    In JCVI-syn minimal genome, metabolic pathways are short (2-5 genes).
    Knocking out any enzyme can disrupt the entire pathway.
    """
    
    def __init__(self):
        self.gene_to_pathway: Dict[str, str] = {}
        self.pathway_genes: Dict[str, List[str]] = {}
    
    def build_from_genes(self, genes: List[GeneInfo]):
        """Build pathway map from gene list."""
        self.gene_to_pathway.clear()
        self.pathway_genes.clear()
        
        for gene in genes:
            gene_name = gene.gene_name.lower() if gene.gene_name else gene.gene_id.lower()
            
            for pathway_id, pathway_info in METABOLIC_PATHWAYS.items():
                for pattern_gene in pathway_info["genes"]:
                    if pattern_gene.lower() in gene_name:
                        self.gene_to_pathway[gene.gene_id] = pathway_id
                        if pathway_id not in self.pathway_genes:
                            self.pathway_genes[pathway_id] = []
                        self.pathway_genes[pathway_id].append(gene.gene_id)
                        break
    
    def get_pathway_effects(self, knocked_gene_id: str) -> List[Tuple[str, float, str]]:
        """
        Get genes affected by pathway disruption.
        
        Returns:
            List of (gene_id, effect_strength, reason) tuples
        """
        effects = []
        
        if knocked_gene_id not in self.gene_to_pathway:
            return effects
        
        pathway_id = self.gene_to_pathway[knocked_gene_id]
        pathway_genes = self.pathway_genes.get(pathway_id, [])
        pathway_info = METABOLIC_PATHWAYS.get(pathway_id, {})
        pathway_name = pathway_info.get("description", pathway_id)
        
        if len(pathway_genes) <= 1:
            return effects
        
        # Shorter pathways = higher confidence in disruption
        pathway_length = len(pathway_genes)
        base_effect = min(0.8, 0.9 - (0.05 * pathway_length))
        
        for i, gene_id in enumerate(pathway_genes):
            if gene_id == knocked_gene_id:
                continue
            
            # Effect based on position in pathway
            distance = abs(i - pathway_genes.index(knocked_gene_id))
            effect = max(0.3, base_effect - (0.1 * distance))
            
            reason = f"Same {pathway_name} (pathway disruption)"
            effects.append((gene_id, effect, reason))
        
        return effects


# =============================================================================
# REGULATORY HUB DETECTOR
# =============================================================================

class RegulatoryHubDetector:
    """
    Identifies regulatory hub genes that affect many downstream targets.
    
    Hub genes include:
    - RNA polymerase subunits (rpoA/B/C/D) - affect ALL transcription
    - Ribosomal proteins (rps/rpl) - affect ALL translation
    - Chaperones (groEL/ES, dnaK/J) - affect protein folding
    """
    
    def __init__(self):
        self.hub_genes: Dict[str, Tuple[str, str]] = {}  # gene_id -> (hub_type, description)
    
    def build_from_genes(self, genes: List[GeneInfo]):
        """Identify hub genes from gene list."""
        self.hub_genes.clear()
        
        for gene in genes:
            gene_name = gene.gene_name if gene.gene_name else gene.gene_id
            
            for pattern, hub_type, description in HUB_GENE_PATTERNS:
                if re.search(pattern, gene_name, re.IGNORECASE):
                    self.hub_genes[gene.gene_id] = (hub_type, description)
                    break
    
    def is_hub_gene(self, gene_id: str) -> bool:
        """Check if gene is a regulatory hub."""
        return gene_id in self.hub_genes
    
    def get_hub_info(self, gene_id: str) -> Optional[Tuple[str, str]]:
        """Get hub type and description for a gene."""
        return self.hub_genes.get(gene_id)
    
    def get_hub_targets(self, hub_gene_id: str, all_genes: List[GeneInfo]) -> List[Tuple[str, float, str]]:
        """
        Get genes affected by knocking out a hub gene.
        
        Returns:
            List of (gene_id, effect_strength, reason) tuples
        """
        effects = []
        
        if hub_gene_id not in self.hub_genes:
            return effects
        
        hub_type, description = self.hub_genes[hub_gene_id]
        hub_name = all_genes[0].gene_name if all_genes else hub_gene_id  # Get actual name
        
        for gene in all_genes:
            if gene.gene_id == hub_gene_id:
                continue
            
            # Determine affected targets based on hub type
            is_affected = False
            effect = 0.0
            reason = ""
            
            if "polymerase" in hub_type.lower():
                # RNA polymerase affects all genes
                is_affected = True
                effect = 0.65
                reason = f"Regulatory cascade from {hub_type} knockout"
                
            elif "ribosomal" in hub_type.lower():
                # Ribosomal proteins affect translation
                if gene.category and any(cat in gene.category.lower() 
                                        for cat in ["translation", "protein", "ribosome"]):
                    is_affected = True
                    effect = 0.70
                    reason = f"Translation machinery affected by {hub_type} knockout"
                else:
                    is_affected = True
                    effect = 0.40
                    reason = f"Protein synthesis reduced by {hub_type} knockout"
                    
            elif "chaperonin" in hub_type.lower() or "chaperone" in hub_type.lower():
                # Chaperones affect protein folding
                is_affected = True
                effect = 0.55
                reason = f"Protein folding affected by {hub_type} knockout"
                
            elif "factor" in hub_type.lower():
                is_affected = True
                effect = 0.50
                reason = f"Translation affected by {hub_type} knockout"
                
            elif "synthase" in hub_type.lower():
                is_affected = True
                effect = 0.45
                reason = f"Energy metabolism affected by {hub_type} knockout"
                
            elif "division" in hub_type.lower():
                if gene.category and "division" in gene.category.lower():
                    is_affected = True
                    effect = 0.75
                    reason = f"Cell division cascade from {hub_type} knockout"
            
            if is_affected:
                effects.append((gene.gene_id, effect, reason))
        
        return effects


# =============================================================================
# MAIN KNOCKOUT RULES ENGINE
# =============================================================================

class KnockoutRulesEngine:
    """
    Main engine that applies all 6 scientific rules for knockout impact.
    
    Usage:
        engine = KnockoutRulesEngine()
        engine.build_from_gene_data(genes)
        effects = engine.calculate_knockout_effects("JCVISYN3_0001")
    """
    
    def __init__(self):
        self.genes: Dict[str, GeneInfo] = {}
        self.gene_list: List[GeneInfo] = []
        self.gene_to_cluster: Dict[str, int] = {}
        
        # Sub-engines for each rule
        self.operon_detector = OperonDetector()
        self.pathway_mapper = PathwayMapper()
        self.hub_detector = RegulatoryHubDetector()
    
    def build_from_gene_data(self, gene_data: List[Dict]):
        """
        Build the rules engine from gene data.
        
        Args:
            gene_data: List of gene dictionaries with keys:
                - gene_id (required)
                - gene_name
                - start, end, strand
                - category, function
                - is_operon_internal, operon_id
        """
        self.gene_list = []
        self.genes = {}
        self.gene_to_cluster = {}
        
        for data in gene_data:
            gene_id = data.get('gene_id', '')
            if not gene_id:
                continue
            
            gene = GeneInfo(
                gene_id=gene_id,
                gene_name=data.get('gene_name', data.get('name', gene_id)),
                start=data.get('start', 0),
                end=data.get('end', 0),
                strand=data.get('strand', '+'),
                category=data.get('category', ''),
                function=data.get('function', data.get('product', '')),
                is_operon_internal=data.get('is_operon_internal', False),
                operon_id=data.get('operon_id', None)
            )
            
            self.gene_list.append(gene)
            self.genes[gene_id] = gene
            
            # Assign to cluster
            if CLUSTERS_AVAILABLE:
                cluster_id = assign_gene_to_cluster(
                    gene_id, gene.gene_name, gene.category, gene.function
                )
                self.gene_to_cluster[gene_id] = cluster_id
                gene.cluster_id = cluster_id
        
        # Build sub-engines
        self.operon_detector.build_from_genes(self.gene_list)
        self.pathway_mapper.build_from_genes(self.gene_list)
        self.hub_detector.build_from_genes(self.gene_list)
        
        print(f"DEBUG: KnockoutRulesEngine built with {len(self.gene_list)} genes")
        print(f"DEBUG: Found {len(self.operon_detector.operons)} operons")
        print(f"DEBUG: Found {len(self.hub_detector.hub_genes)} hub genes")
    
    def calculate_knockout_effects(self, knocked_gene_id: str) -> Dict:
        """
        Calculate all knockout effects using the 6 scientific rules.
        
        Returns:
            Dict containing:
            - knocked_gene: The knocked out gene ID
            - knocked_gene_info: Information about the knocked gene
            - is_hub_gene: Whether this is a regulatory hub
            - effects: List of KnockoutEffect objects
            - effects_by_rule: Dict grouping effects by rule type
            - total_affected: Total number of affected genes
        """
        if knocked_gene_id not in self.genes:
            return {
                'knocked_gene': knocked_gene_id,
                'knocked_gene_info': {},
                'is_hub_gene': False,
                'effects': [],
                'effects_by_rule': {},
                'total_affected': 0
            }
        
        knocked_gene = self.genes[knocked_gene_id]
        is_hub = self.hub_detector.is_hub_gene(knocked_gene_id)
        
        all_effects: Dict[str, KnockoutEffect] = {}
        
        # === RULE 1: Operon/Polar Effects (HIGH confidence) ===
        polar_effects = self.operon_detector.get_polar_effects(knocked_gene_id)
        for gene_id, effect, reason in polar_effects:
            cluster_name = self._get_cluster_name(gene_id)
            all_effects[gene_id] = KnockoutEffect(
                gene_id=gene_id,
                effect_strength=effect,
                confidence="high",
                reason=reason,
                rule_type="operon",
                distance=0,
                cluster_name=cluster_name
            )
        
        # === RULE 2: Core Cellular Process (MEDIUM confidence) ===
        # Genes in same cluster AND connected clusters get effects with INDIVIDUAL VARIATION
        knocked_cluster = self.gene_to_cluster.get(knocked_gene_id, -1)
        knocked_gene_pos = knocked_gene.start if knocked_gene.start > 0 else 1
        
        # Same cluster effects - with variation based on genomic distance
        cluster_name = CORE_PROCESS_CLUSTERS.get(knocked_cluster, self._get_cluster_name(knocked_gene_id))
        importance = CLUSTER_IMPORTANCE.get(knocked_cluster, 1.0)
        
        for gene_id, gene in self.genes.items():
            if gene_id == knocked_gene_id or gene_id in all_effects:
                continue
            
            gene_cluster = self.gene_to_cluster.get(gene_id, -1)
            
            if gene_cluster == knocked_cluster:
                # Same cluster: Effect varies by genomic distance
                gene_pos = gene.start if gene.start > 0 else 1
                genomic_distance = abs(gene_pos - knocked_gene_pos)
                
                # Closer genes = higher effect (range 45%-85%)
                # Every 10kb reduces effect by ~5%
                distance_decay = min(0.4, genomic_distance / 100000)  # Max 40% reduction
                base_effect = 0.85 * importance
                effect = max(0.45, base_effect - distance_decay)
                
                # Add small random variation for biological realism (±3%)
                import random
                noise = random.uniform(-0.03, 0.03)
                effect = max(0.40, min(0.90, effect + noise))
                
                # Create reason with distance info
                if genomic_distance < 5000:
                    proximity_desc = "nearby"
                elif genomic_distance < 20000:
                    proximity_desc = "moderately distant"
                else:
                    proximity_desc = "distant"
                
                reason = f"Same core {cluster_name} process ({proximity_desc})"
                
                all_effects[gene_id] = KnockoutEffect(
                    gene_id=gene_id,
                    effect_strength=effect,
                    confidence="medium",
                    reason=reason,
                    rule_type="core_process",
                    distance=1,
                    cluster_name=cluster_name
                )
        
        # Connected clusters (from ENTANGLEMENT_PAIRS) - lower effect
        if CLUSTERS_AVAILABLE:
            from cluster_definitions import ENTANGLEMENT_PAIRS
            connected_clusters = set()
            for c1, c2 in ENTANGLEMENT_PAIRS:
                if c1 == knocked_cluster:
                    connected_clusters.add(c2)
                elif c2 == knocked_cluster:
                    connected_clusters.add(c1)
            
            for target_cluster in connected_clusters:
                target_cluster_name = get_cluster_name(target_cluster)
                for gene_id, gene in self.genes.items():
                    if gene_id == knocked_gene_id or gene_id in all_effects:
                        continue
                    
                    if self.gene_to_cluster.get(gene_id, -1) == target_cluster:
                        # Connected cluster effect (30%-55%)
                        gene_pos = gene.start if gene.start > 0 else 1
                        genomic_distance = abs(gene_pos - knocked_gene_pos)
                        
                        # Further decay for connected clusters
                        distance_decay = min(0.25, genomic_distance / 200000)
                        effect = max(0.30, 0.55 - distance_decay)
                        
                        # Add noise
                        import random
                        noise = random.uniform(-0.03, 0.03)
                        effect = max(0.25, min(0.60, effect + noise))
                        
                        reason = f"Connected to {cluster_name} via {target_cluster_name}"
                        
                        all_effects[gene_id] = KnockoutEffect(
                            gene_id=gene_id,
                            effect_strength=effect,
                            confidence="low",
                            reason=reason,
                            rule_type="core_process",
                            distance=2,
                            cluster_name=target_cluster_name
                        )
        
        # === RULE 3: Metabolic Pathway (MEDIUM confidence) ===
        pathway_effects = self.pathway_mapper.get_pathway_effects(knocked_gene_id)
        for gene_id, effect, reason in pathway_effects:
            if gene_id not in all_effects:
                cluster_name = self._get_cluster_name(gene_id)
                all_effects[gene_id] = KnockoutEffect(
                    gene_id=gene_id,
                    effect_strength=effect,
                    confidence="medium",
                    reason=reason,
                    rule_type="pathway",
                    distance=1,
                    cluster_name=cluster_name
                )
        
        # === RULE 4: Regulatory Hub (MEDIUM confidence) ===
        # Hub effects also have individual variation
        if is_hub:
            hub_type, hub_desc = self.hub_detector.get_hub_info(knocked_gene_id)
            
            for gene in self.gene_list:
                if gene.gene_id == knocked_gene_id or gene.gene_id in all_effects:
                    continue
                
                # Calculate base effect based on hub type
                base_effect = 0.0
                reason = ""
                
                if "polymerase" in hub_type.lower():
                    base_effect = 0.65
                    reason = f"Transcription reduced by {hub_type} knockout"
                elif "ribosomal" in hub_type.lower():
                    if gene.category and any(cat in gene.category.lower() 
                                            for cat in ["translation", "protein", "ribosome"]):
                        base_effect = 0.70
                    else:
                        base_effect = 0.40
                    reason = f"Translation affected by {hub_type} knockout"
                elif "chaperonin" in hub_type.lower() or "chaperone" in hub_type.lower():
                    base_effect = 0.55
                    reason = f"Protein folding affected by {hub_type} knockout"
                elif "synthase" in hub_type.lower():
                    base_effect = 0.45
                    reason = f"Energy metabolism affected by {hub_type} knockout"
                else:
                    base_effect = 0.35
                    reason = f"Regulatory cascade from {hub_type}"
                
                if base_effect > 0:
                    # Add variation based on genomic distance
                    gene_pos = gene.start if gene.start > 0 else 1
                    genomic_distance = abs(gene_pos - knocked_gene_pos)
                    distance_decay = min(0.20, genomic_distance / 150000)
                    effect = max(0.20, base_effect - distance_decay)
                    
                    # Add noise
                    import random
                    noise = random.uniform(-0.04, 0.04)
                    effect = max(0.15, min(0.80, effect + noise))
                    
                    all_effects[gene.gene_id] = KnockoutEffect(
                        gene_id=gene.gene_id,
                        effect_strength=effect,
                        confidence="medium",
                        reason=reason,
                        rule_type="hub",
                        distance=2,
                        cluster_name=self._get_cluster_name(gene.gene_id)
                    )
        
        # === RULE 5: Genome Proximity (LOW confidence) ===
        # Only add if not already affected by higher-confidence rules
        if knocked_gene.start > 0:
            for gene in self.gene_list:
                if gene.gene_id == knocked_gene_id or gene.gene_id in all_effects:
                    continue
                
                if gene.start > 0:
                    distance = abs(gene.start - knocked_gene.start)
                    if distance < 5000:  # Within 5kb
                        effect = 0.15 * (1 - distance / 5000)
                        if effect > 0.05:  # Only include if meaningful
                            cluster_name = self._get_cluster_name(gene.gene_id)
                            all_effects[gene.gene_id] = KnockoutEffect(
                                gene_id=gene.gene_id,
                                effect_strength=effect,
                                confidence="low",
                                reason="Genome proximity (possible co-regulation)",
                                rule_type="proximity",
                                distance=3,
                                cluster_name=cluster_name
                            )
        
        # === RULE 6: Essentiality/Stress Response (LOW-MEDIUM confidence) ===
        is_essential = self._is_essential_gene(knocked_gene)
        if is_essential:
            for gene in self.gene_list:
                if gene.gene_id == knocked_gene_id or gene.gene_id in all_effects:
                    continue
                
                gene_cluster = self.gene_to_cluster.get(gene.gene_id, -1)
                
                # Stress response genes
                if gene_cluster == STRESS_CLUSTER:
                    all_effects[gene.gene_id] = KnockoutEffect(
                        gene_id=gene.gene_id,
                        effect_strength=0.40,
                        confidence="low",
                        reason="Stress response to essential gene knockout",
                        rule_type="stress",
                        distance=4,
                        cluster_name=self._get_cluster_name(gene.gene_id)
                    )
                
                # Energy cluster
                elif gene_cluster == 14:  # Energy/ATP
                    all_effects[gene.gene_id] = KnockoutEffect(
                        gene_id=gene.gene_id,
                        effect_strength=0.35,
                        confidence="low",
                        reason="Energy rebalance due to essential gene knockout",
                        rule_type="stress",
                        distance=4,
                        cluster_name=self._get_cluster_name(gene.gene_id)
                    )
        
        # Convert to list and sort by effect strength
        effects_list = list(all_effects.values())
        effects_list.sort(key=lambda e: (-e.effect_strength, e.distance))
        
        # Group by rule type
        effects_by_rule = defaultdict(list)
        for effect in effects_list:
            effects_by_rule[effect.rule_type].append(effect)
        
        return {
            'knocked_gene': knocked_gene_id,
            'knocked_gene_info': {
                'gene_name': knocked_gene.gene_name,
                'category': knocked_gene.category,
                'function': knocked_gene.function,
                'cluster_id': knocked_gene.cluster_id,
                'cluster_name': self._get_cluster_name(knocked_gene_id)
            },
            'is_hub_gene': is_hub,
            'hub_info': self.hub_detector.get_hub_info(knocked_gene_id) if is_hub else None,
            'effects': effects_list,
            'effects_by_rule': dict(effects_by_rule),
            'total_affected': len(effects_list)
        }
    
    def _get_cluster_name(self, gene_id: str) -> str:
        """Get cluster name for a gene."""
        cluster_id = self.gene_to_cluster.get(gene_id, -1)
        if CLUSTERS_AVAILABLE and cluster_id >= 0:
            return get_cluster_name(cluster_id)
        return "Unknown"
    
    def _is_essential_gene(self, gene: GeneInfo) -> bool:
        """Check if gene is essential or near-essential."""
        # Check category
        if gene.category:
            for cat in ESSENTIAL_CATEGORIES:
                if cat.lower() in gene.category.lower():
                    return True
        
        # Check function
        if gene.function:
            for cat in ESSENTIAL_CATEGORIES:
                if cat.lower() in gene.function.lower():
                    return True
        
        # Check if it's a hub gene
        return self.hub_detector.is_hub_gene(gene.gene_id)


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def get_confidence_icon(confidence: str) -> str:
    """Get emoji icon for confidence level."""
    icons = {
        "high": "🔴",
        "medium": "🟠",
        "low": "🟡"
    }
    return icons.get(confidence, "⚪")


def get_rule_icon(rule_type: str) -> str:
    """Get emoji icon for rule type."""
    icons = {
        "operon": "📜",      # Scroll - transcription unit
        "core_process": "⚙️", # Gear - machinery
        "pathway": "🔬",      # Microscope - metabolic
        "hub": "🌐",          # Globe - network
        "proximity": "📍",    # Pin - location
        "stress": "⚠️"        # Warning - stress response
    }
    return icons.get(rule_type, "❓")


def format_effect_for_display(effect: KnockoutEffect) -> Dict:
    """Format effect for UI display."""
    return {
        'gene_id': effect.gene_id,
        'effect_percent': f"{int(effect.effect_strength * 100)}%",
        'confidence_icon': get_confidence_icon(effect.confidence),
        'confidence': effect.confidence,
        'rule_icon': get_rule_icon(effect.rule_type),
        'rule_type': effect.rule_type,
        'reason': effect.reason,
        'cluster_name': effect.cluster_name
    }
