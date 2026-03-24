"""
Gene Dependency Graph for JCVI Genome Analyzer
===============================================

Builds and manages gene-to-gene dependency networks for knockout cascade effects.

Uses cluster entanglement pairs to model biological dependencies:
- Genes in the same cluster are strongly coupled
- Genes in entangled clusters affect each other
- Effect strength decays with graph distance

Example Usage:
    graph = GeneDependencyGraph()
    graph.build_from_gene_data(gene_data)
    
    # When gene is knocked out:
    effects = graph.calculate_cascade("JCVISYN3_0001")
    # effects = {
    #     'knocked_gene': 'JCVISYN3_0001',
    #     'directly_affected': [('JCVISYN3_0002', 0.8), ...],
    #     'indirectly_affected': [('JCVISYN3_0010', 0.5, 2), ...],
    #     'cluster_impact': {0: 1.0, 1: 0.8, ...}
    # }
"""

from typing import Dict, List, Set, Tuple, Optional
from collections import defaultdict, deque
import re

# Import the new scientific rules engine
try:
    from knockout_rules_engine import (
        KnockoutRulesEngine, KnockoutEffect,
        get_confidence_icon, get_rule_icon, format_effect_for_display
    )
    RULES_ENGINE_AVAILABLE = True
except ImportError:
    RULES_ENGINE_AVAILABLE = False
    print("Warning: knockout_rules_engine not available, using legacy mode")

# Import cluster definitions
try:
    from cluster_definitions import (
        CLUSTER_DEFINITIONS, ENTANGLEMENT_PAIRS,
        assign_gene_to_cluster, get_cluster_name
    )
except ImportError:
    # Fallback definitions if import fails
    CLUSTER_DEFINITIONS = {}
    ENTANGLEMENT_PAIRS = []
    def assign_gene_to_cluster(gene_id, gene_name="", category="", function=""):
        return 27
    def get_cluster_name(cluster_id):
        return f"Cluster {cluster_id}"


# =============================================================================
# EFFECT STRENGTH CONSTANTS
# =============================================================================

# Effect strength by distance (biological plausibility)
EFFECT_STRENGTH_BY_DISTANCE = {
    0: 1.0,    # Same gene = 100% (the knocked out gene itself)
    1: 0.80,   # Direct dependency = 80%
    2: 0.50,   # 2nd order = 50%
    3: 0.30,   # 3rd order = 30%
    4: 0.15,   # 4th order = 15%
    5: 0.05,   # 5th order = 5%
}

# Genes with high impact when knocked out (essential machinery)
HIGH_IMPACT_GENE_PATTERNS = [
    r'^rpo[A-Z]',   # RNA polymerase - affects ALL transcription
    r'^rps[A-Z]',   # Ribosomal proteins - affects ALL translation
    r'^rpl[A-Z]',   # Ribosomal proteins
    r'^gro[ELS]',   # Chaperones - affects protein folding
    r'^dna[A-Z]',   # DNA replication - affects cell division
    r'^fts[A-Z]',   # Cell division
    r'^atp[A-Z]',   # ATP synthase - affects energy
]

# Cluster importance multiplier (central dogma clusters are critical)
CLUSTER_IMPORTANCE = {
    0: 1.5,   # Ribosome Assembly - critical
    1: 1.5,   # Translation - critical
    2: 1.3,   # tRNA Synthetases - very important
    3: 1.4,   # DNA Replication - critical
    4: 1.5,   # Transcription - critical
    5: 1.2,   # RNA Processing
    6: 1.0,   # DNA Repair
    7: 1.2,   # Cell Division
    11: 1.3,  # Protein Folding - important
    14: 1.3,  # Energy/ATP - important
}


class GeneDependencyGraph:
    """
    Gene-to-gene dependency graph for modeling knockout cascade effects.
    
    The graph has two levels:
    1. CLUSTER LEVEL: Uses ENTANGLEMENT_PAIRS to connect functional clusters
    2. GENE LEVEL: Genes within same/connected clusters affect each other
    
    When a gene is knocked out:
    1. All genes in the same cluster are directly affected
    2. Genes in entangled clusters are indirectly affected
    3. Effect propagates through the cluster network with decay
    """
    
    def __init__(self):
        # Gene data
        self.genes: Dict[str, Dict] = {}              # gene_id -> gene info
        self.gene_to_cluster: Dict[str, int] = {}     # gene_id -> cluster_id
        self.cluster_to_genes: Dict[int, Set[str]] = defaultdict(set)  # cluster_id -> set of gene_ids
        
        # Cluster graph (adjacency list from ENTANGLEMENT_PAIRS)
        self.cluster_adjacency: Dict[int, Set[int]] = defaultdict(set)
        
        # Precomputed cluster distances
        self.cluster_distances: Dict[Tuple[int, int], int] = {}
        
        # Scientific rules engine for biologically accurate knockout effects
        self.rules_engine: Optional[KnockoutRulesEngine] = None
        if RULES_ENGINE_AVAILABLE:
            self.rules_engine = KnockoutRulesEngine()
        
        # Build initial cluster graph from ENTANGLEMENT_PAIRS
        self._build_cluster_graph()
    
    def _build_cluster_graph(self):
        """Build cluster adjacency graph from ENTANGLEMENT_PAIRS."""
        for cluster_a, cluster_b in ENTANGLEMENT_PAIRS:
            self.cluster_adjacency[cluster_a].add(cluster_b)
            self.cluster_adjacency[cluster_b].add(cluster_a)
        
        # Precompute distances between all cluster pairs (BFS)
        self._compute_all_cluster_distances()
    
    def _compute_all_cluster_distances(self):
        """Precompute shortest distances between all cluster pairs using BFS."""
        for start_cluster in range(28):
            distances = self._bfs_cluster_distances(start_cluster)
            for end_cluster, distance in distances.items():
                self.cluster_distances[(start_cluster, end_cluster)] = distance
    
    def _bfs_cluster_distances(self, start: int) -> Dict[int, int]:
        """BFS to find shortest paths from start cluster to all others."""
        distances = {start: 0}
        queue = deque([start])
        
        while queue:
            current = queue.popleft()
            current_dist = distances[current]
            
            for neighbor in self.cluster_adjacency[current]:
                if neighbor not in distances:
                    distances[neighbor] = current_dist + 1
                    queue.append(neighbor)
        
        return distances
    
    def build_from_gene_data(self, gene_data: List[Dict]):
        """
        Build the dependency graph from gene data.
        
        Args:
            gene_data: List of gene dictionaries with keys:
                - gene_id: Unique identifier
                - gene_name: Optional gene name
                - category: Optional category
                - function: Optional function description
                - product: Optional product description
                - start, end, strand: Genomic coordinates
        """
        self.genes.clear()
        self.gene_to_cluster.clear()
        self.cluster_to_genes.clear()
        
        for gene in gene_data:
            gene_id = gene.get('gene_id', '')
            if not gene_id:
                continue
            
            # Store gene info
            self.genes[gene_id] = gene
            
            # Assign to cluster
            gene_name = gene.get('gene_name', gene.get('name', gene_id))
            category = gene.get('category', '')
            function = gene.get('function', gene.get('product', ''))
            
            cluster_id = assign_gene_to_cluster(gene_id, gene_name, category, function)
            self.gene_to_cluster[gene_id] = cluster_id
            self.cluster_to_genes[cluster_id].add(gene_id)
        
        # Build the scientific rules engine with complete gene data
        if self.rules_engine:
            self.rules_engine.build_from_gene_data(gene_data)
        
        print(f"DEBUG: GeneDependencyGraph built with {len(self.genes)} genes in {len(self.cluster_to_genes)} clusters")
    
    def get_cluster_distance(self, cluster_a: int, cluster_b: int) -> int:
        """Get shortest path distance between two clusters."""
        if cluster_a == cluster_b:
            return 0
        return self.cluster_distances.get((cluster_a, cluster_b), 999)
    
    def _is_high_impact_gene(self, gene_id: str) -> bool:
        """Check if gene is in the high-impact list (essential machinery)."""
        gene_name = self.genes.get(gene_id, {}).get('gene_name', gene_id).lower()
        
        for pattern in HIGH_IMPACT_GENE_PATTERNS:
            if re.search(pattern, gene_name, re.IGNORECASE):
                return True
        return False
    
    def calculate_cascade(self, knocked_out_gene: str) -> Dict:
        """
        Calculate all genes affected by knocking out a gene.
        
        Uses the scientific rules engine for biologically accurate effects:
        1. Operon/Polar Effect (HIGH confidence)
        2. Core Cellular Process (MEDIUM confidence)
        3. Metabolic Pathway (MEDIUM confidence)
        4. Regulatory Hub (MEDIUM confidence)
        5. Genome Proximity (LOW confidence)
        6. Essentiality/Stress Response (LOW-MEDIUM confidence)
        
        Returns:
            Dict with:
            - knocked_gene: The knocked out gene ID
            - knocked_gene_info: Info about the knocked gene
            - cluster_id: Cluster of the knocked gene
            - cluster_name: Name of the cluster
            - is_high_impact: Whether this is a regulatory hub
            - effects: List of KnockoutEffect objects (if rules engine available)
            - effects_by_rule: Effects grouped by rule type
            - same_cluster_genes: Legacy format - List of (gene_id, effect, reason, confidence)
            - connected_genes: Legacy format - included for backward compatibility
            - cluster_impact: Dict of cluster_id -> effect percentage
            - total_affected: Total number of affected genes
        """
        # === USE SCIENTIFIC RULES ENGINE (preferred) ===
        if self.rules_engine and RULES_ENGINE_AVAILABLE:
            result = self.rules_engine.calculate_knockout_effects(knocked_out_gene)
            
            # Convert to combined format for backward compatibility
            knocked_cluster = self.gene_to_cluster.get(knocked_out_gene, -1)
            
            # Build same_cluster_genes and connected_genes from effects
            same_cluster_genes = []
            connected_genes = []
            cluster_impact = {knocked_cluster: 1.0} if knocked_cluster >= 0 else {}
            
            for effect in result.get('effects', []):
                gene_cluster = self.gene_to_cluster.get(effect.gene_id, -1)
                
                # Tuple format: (gene_id, effect_strength, reason, confidence, rule_type)
                effect_tuple = (
                    effect.gene_id,
                    effect.effect_strength,
                    effect.reason,
                    effect.confidence,
                    effect.rule_type
                )
                
                if gene_cluster == knocked_cluster:
                    same_cluster_genes.append(effect_tuple)
                else:
                    connected_genes.append(effect_tuple + (effect.cluster_name,))
                    
                    # Update cluster impact
                    if gene_cluster >= 0 and gene_cluster not in cluster_impact:
                        cluster_impact[gene_cluster] = effect.effect_strength
            
            return {
                'knocked_gene': knocked_out_gene,
                'knocked_gene_info': result.get('knocked_gene_info', {}),
                'cluster_id': knocked_cluster,
                'cluster_name': get_cluster_name(knocked_cluster) if knocked_cluster >= 0 else 'Unknown',
                'is_high_impact': result.get('is_hub_gene', False),
                'hub_info': result.get('hub_info'),
                'effects': result.get('effects', []),
                'effects_by_rule': result.get('effects_by_rule', {}),
                'same_cluster_genes': same_cluster_genes,
                'connected_genes': connected_genes,
                'cluster_impact': cluster_impact,
                'total_affected': result.get('total_affected', 0),
                'uses_scientific_rules': True
            }
        
        # === FALLBACK: Legacy cluster-based calculation ===
        return self._calculate_cascade_legacy(knocked_out_gene)
    
    def _calculate_cascade_legacy(self, knocked_out_gene: str) -> Dict:
        """Legacy cascade calculation (cluster-based, no scientific rules)."""
        if knocked_out_gene not in self.gene_to_cluster:
            return {
                'knocked_gene': knocked_out_gene,
                'knocked_gene_info': {},
                'cluster_id': -1,
                'cluster_name': 'Unknown',
                'is_high_impact': False,
                'same_cluster_genes': [],
                'connected_genes': [],
                'cluster_impact': {},
                'total_affected': 0,
                'uses_scientific_rules': False
            }
        
        knocked_cluster = self.gene_to_cluster[knocked_out_gene]
        is_high_impact = self._is_high_impact_gene(knocked_out_gene)
        
        # Impact multiplier for high-impact genes
        impact_multiplier = 1.5 if is_high_impact else 1.0
        cluster_importance = CLUSTER_IMPORTANCE.get(knocked_cluster, 1.0)
        combined_multiplier = min(impact_multiplier * cluster_importance, 2.0)
        
        # Genes in same cluster
        same_cluster_genes = []
        for gene_id in self.cluster_to_genes[knocked_cluster]:
            if gene_id != knocked_out_gene:
                effect = min(0.70 * combined_multiplier, 1.0)
                reason = "Same functional cluster"
                same_cluster_genes.append((gene_id, effect, reason, "medium", "core_process"))
        
        # Genes in connected clusters
        connected_genes = []
        cluster_impact = {knocked_cluster: 1.0}
        
        for cluster_id in range(28):
            if cluster_id == knocked_cluster:
                continue
            
            distance = self.get_cluster_distance(knocked_cluster, cluster_id)
            if distance >= 999 or distance > 3:
                continue
            
            base_effect = EFFECT_STRENGTH_BY_DISTANCE.get(distance, 0.0)
            if base_effect <= 0:
                continue
            
            cluster_effect = min(base_effect * combined_multiplier, 1.0)
            cluster_impact[cluster_id] = cluster_effect
            cluster_name = get_cluster_name(cluster_id)
            
            for gene_id in self.cluster_to_genes[cluster_id]:
                reason = f"Connected cluster ({cluster_name})"
                connected_genes.append((gene_id, cluster_effect, reason, "low", "proximity", cluster_name))
        
        return {
            'knocked_gene': knocked_out_gene,
            'knocked_gene_info': self.genes.get(knocked_out_gene, {}),
            'cluster_id': knocked_cluster,
            'cluster_name': get_cluster_name(knocked_cluster),
            'is_high_impact': is_high_impact,
            'same_cluster_genes': same_cluster_genes,
            'connected_genes': connected_genes,
            'cluster_impact': cluster_impact,
            'total_affected': len(same_cluster_genes) + len(connected_genes),
            'uses_scientific_rules': False
        }
    
    def _extract_gene_number(self, gene_id: str) -> int:
        """Extract numeric suffix from gene ID (e.g., JCVISYN3_0654 -> 654)."""
        match = re.search(r'_(\d+)$', gene_id)
        if match:
            return int(match.group(1))
        return 0
    
    def get_effect_on_gene(self, knocked_gene: str, target_gene: str) -> float:
        """
        Get the effect strength on a target gene when knocked_gene is knocked out.
        
        Returns:
            Effect strength (0.0 to 1.0), where:
            - 0.0 = no effect
            - 1.0 = fully affected (blocked)
        """
        if knocked_gene not in self.gene_to_cluster or target_gene not in self.gene_to_cluster:
            return 0.0
        
        if knocked_gene == target_gene:
            return 1.0  # Knocked out gene is fully affected
        
        knocked_cluster = self.gene_to_cluster[knocked_gene]
        target_cluster = self.gene_to_cluster[target_gene]
        
        # Same cluster
        if knocked_cluster == target_cluster:
            is_high_impact = self._is_high_impact_gene(knocked_gene)
            multiplier = 1.5 if is_high_impact else 1.0
            return min(0.85 * multiplier, 1.0)
        
        # Different clusters - use distance
        distance = self.get_cluster_distance(knocked_cluster, target_cluster)
        base_effect = EFFECT_STRENGTH_BY_DISTANCE.get(distance, 0.0)
        
        is_high_impact = self._is_high_impact_gene(knocked_gene)
        cluster_importance = CLUSTER_IMPORTANCE.get(knocked_cluster, 1.0)
        multiplier = min(1.5 if is_high_impact else 1.0 * cluster_importance, 2.0)
        
        return min(base_effect * multiplier, 1.0)
    
    def get_genes_affecting(self, gene_id: str) -> List[Tuple[str, float]]:
        """
        Get all genes that would affect this gene if knocked out.
        
        Returns:
            List of (affecting_gene_id, effect_strength) tuples
        """
        affecting = []
        
        for other_gene in self.genes:
            if other_gene == gene_id:
                continue
            effect = self.get_effect_on_gene(other_gene, gene_id)
            if effect > 0:
                affecting.append((other_gene, effect))
        
        # Sort by effect strength
        affecting.sort(key=lambda x: -x[1])
        return affecting
    
    def get_statistics(self) -> Dict:
        """Get graph statistics."""
        return {
            'total_genes': len(self.genes),
            'total_clusters_used': len(self.cluster_to_genes),
            'genes_per_cluster': {
                get_cluster_name(cid): len(genes) 
                for cid, genes in self.cluster_to_genes.items()
            },
            'cluster_connections': len(ENTANGLEMENT_PAIRS),
            'high_impact_genes': sum(
                1 for gene_id in self.genes 
                if self._is_high_impact_gene(gene_id)
            )
        }


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def get_effect_color(effect_strength: float) -> Tuple[int, int, int]:
    """
    Get RGB color based on effect strength for visualization.
    
    Returns color gradient:
    - Red (255, 50, 50) for high effect (>0.7)
    - Orange (255, 165, 0) for medium effect (0.4-0.7)
    - Yellow (255, 255, 0) for low effect (0.1-0.4)
    - Light yellow (255, 255, 180) for very low effect (<0.1)
    """
    if effect_strength >= 0.7:
        return (255, 50, 50)      # Red
    elif effect_strength >= 0.4:
        return (255, 165, 0)      # Orange
    elif effect_strength >= 0.1:
        return (255, 255, 0)      # Yellow
    else:
        return (255, 255, 180)    # Light yellow


def get_effect_icon(effect_strength: float) -> str:
    """Get emoji icon based on effect strength."""
    if effect_strength >= 0.7:
        return "🔴"  # High impact
    elif effect_strength >= 0.4:
        return "🟠"  # Medium impact
    elif effect_strength >= 0.1:
        return "🟡"  # Low impact
    else:
        return "⚪"  # Minimal impact


def format_effect_percentage(effect_strength: float) -> str:
    """Format effect strength as percentage string."""
    return f"{int(effect_strength * 100)}%"
